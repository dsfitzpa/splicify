"""External plasmid-repository search and fetch.

Helpers for finding plasmids by free-text description in public
repositories (currently Addgene) and pulling back enough metadata
+ the actual GenBank file to feed into the standard annotation
pipeline. Plus PubMed / CrossRef lookups for the depositor paper.

Pure HTTP — no Anthropic / agent dependencies. Higher layers
(PartScout tools, the interpreter's escalate tool) call into here.
"""
from __future__ import annotations

import asyncio
import logging
import re
from dataclasses import dataclass
from typing import Any, Optional

import httpx
from bs4 import BeautifulSoup


logger = logging.getLogger(__name__)


_USER_AGENT = "AIPlasmidDesign-agent_v2/1.0 (research tool)"

# Process-wide rate gate for NCBI E-utilities. The cap is 10 req/sec with
# an api_key (3/sec without). The asyncio.Lock serialises only the gate-
# keeping — actual HTTP requests fire in parallel AFTER passing through.
# Each pass enforces >= _NCBI_PACING_SECONDS since the previous gate exit,
# capping the GLOBAL outgoing rate regardless of how many search /fetch
# coroutines run at once. Without this, 3 parallel CGAS+KEAP1+TP53
# lookups blast 12 calls in under 1 second and trip 429.
_NCBI_PACING_SECONDS = 0.12   # -> ~8 req/sec global ceiling
import time as _time
_NCBI_RATE_LOCK = asyncio.Lock()
_NCBI_LAST_CALL_TS = 0.0


async def _pace_ncbi():
    """Block until the next NCBI request slot is available."""
    global _NCBI_LAST_CALL_TS
    async with _NCBI_RATE_LOCK:
        now = _time.monotonic()
        elapsed = now - _NCBI_LAST_CALL_TS
        if elapsed < _NCBI_PACING_SECONDS:
            await asyncio.sleep(_NCBI_PACING_SECONDS - elapsed)
        _NCBI_LAST_CALL_TS = _time.monotonic()
_HTTP_TIMEOUT = 15.0
_ADDGENE_SEARCH = "https://www.addgene.org/search/catalog/plasmids/"
_ADDGENE_ENTRY = "https://www.addgene.org/{aid}/"
_ADDGENE_SEQS = "https://www.addgene.org/{aid}/sequences/"
_PUBMED_ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
_PUBMED_ABSTRACT = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
_CROSSREF_WORKS = "https://api.crossref.org/works/{doi}"


@dataclass
class AddgeneCandidate:
    addgene_id: str
    name: str
    url: str
    snippet: Optional[str] = None

    def to_dict(self) -> dict[str, Any]:
        return {"addgene_id": self.addgene_id, "name": self.name,
                "url": self.url, "snippet": self.snippet}


@dataclass
class AddgeneEntry:
    addgene_id: str
    name: str
    url: str
    description: Optional[str]
    depositor: Optional[str]
    pmid: Optional[str]
    doi: Optional[str]
    paper_title: Optional[str]
    gb_url: Optional[str]

    def to_dict(self) -> dict[str, Any]:
        return {
            "addgene_id": self.addgene_id, "name": self.name, "url": self.url,
            "description": self.description, "depositor": self.depositor,
            "pmid": self.pmid, "doi": self.doi, "paper_title": self.paper_title,
            "gb_url": self.gb_url,
        }


# ─────────────────────────────────────────────────────────────────────
# Public helpers
# ─────────────────────────────────────────────────────────────────────
async def search_addgene(query: str, *, max_results: int = 5,
                          client: Optional[httpx.AsyncClient] = None) -> list[AddgeneCandidate]:
    """Free-text search of Addgene's plasmid catalog. Returns
    AddgeneCandidate objects sourced from the search-results HTML."""
    if not (query or "").strip():
        return []
    own_client = False
    if client is None:
        client = httpx.AsyncClient(timeout=_HTTP_TIMEOUT, headers={"User-Agent": _USER_AGENT},
                                    follow_redirects=True)
        own_client = True
    try:
        params = {"q": query.strip()}
        r = await client.get(_ADDGENE_SEARCH, params=params)
        r.raise_for_status()
        return _parse_addgene_search(r.text, max_results=max_results)
    except Exception as e:
        logger.warning("addgene search failed for %r: %s", query, e)
        return []
    finally:
        if own_client:
            await client.aclose()


async def fetch_addgene_entry(addgene_id: str, *,
                               client: Optional[httpx.AsyncClient] = None) -> Optional[AddgeneEntry]:
    """Pull metadata + .gb URL for one Addgene entry. Returns None
    when the page can't be fetched or parsed."""
    aid = (addgene_id or "").strip().lstrip("#")
    if not aid.isdigit():
        return None
    own_client = False
    if client is None:
        client = httpx.AsyncClient(timeout=_HTTP_TIMEOUT, headers={"User-Agent": _USER_AGENT},
                                    follow_redirects=True)
        own_client = True
    try:
        entry_url = _ADDGENE_ENTRY.format(aid=aid)
        seqs_url = _ADDGENE_SEQS.format(aid=aid)
        try:
            ep = await client.get(entry_url)
            ep.raise_for_status()
        except Exception as e:
            logger.warning("addgene entry %s fetch failed: %s", aid, e)
            return None
        meta = _parse_addgene_entry(ep.text, aid=aid, entry_url=entry_url)
        # Two-hop GenBank URL discovery: the /<aid>/sequences/ page lists
        # available sequences (with section headers like "Addgene-Verified
        # Full Sequences") that link to /browse/sequence/<seq_id>/ pages.
        # The actual .gbk URL lives on the per-sequence browse page,
        # embedded in a JS payload. We pick the most-canonical sequence
        # (verified full > depositor full > partials) and scrape its .gbk.
        try:
            sp = await client.get(seqs_url)
            if sp.status_code == 200:
                # Fallback first: maybe a direct .gbk link is on the page.
                gb_url = _parse_addgene_gb_url(sp.text, aid=aid)
                if not gb_url:
                    seq_ids = _parse_addgene_seq_ids(sp.text)
                    for seq_id in seq_ids[:3]:
                        try:
                            bp = await client.get(
                                f"https://www.addgene.org/browse/sequence/{seq_id}/"
                            )
                            if bp.status_code == 200:
                                gb_url = _parse_addgene_gb_url(bp.text, aid=aid)
                                if gb_url:
                                    break
                        except Exception:
                            continue
                meta.gb_url = gb_url
        except Exception:
            pass
        if meta.pmid:
            try:
                paper = await fetch_pubmed_summary(meta.pmid, client=client)
                if paper and paper.get("title"):
                    meta.paper_title = paper["title"]
            except Exception:
                pass
        elif meta.doi:
            try:
                paper = await fetch_crossref_summary(meta.doi, client=client)
                if paper and paper.get("title"):
                    meta.paper_title = paper["title"]
            except Exception:
                pass
        return meta
    finally:
        if own_client:
            await client.aclose()


async def download_addgene_gb(addgene_id: str, *, gb_url: Optional[str] = None,
                               client: Optional[httpx.AsyncClient] = None) -> Optional[str]:
    """Download a plasmid's GenBank file. If `gb_url` isn't supplied,
    fetches the entry first to discover it. Returns the raw .gb text.

    The .gbk files live on `media.addgene.org`, which requires the
    `__Secure_media_edge_auth` cookie set by an `addgene.org` page
    visit. When `client` is None we open our own client and warm it
    against `addgene.org/<aid>/` so the cookie is in the jar before
    the media.addgene.org request.
    """
    own_client = False
    if client is None:
        client = httpx.AsyncClient(timeout=_HTTP_TIMEOUT,
                                    headers={"User-Agent": _USER_AGENT},
                                    follow_redirects=True)
        own_client = True
    try:
        if not gb_url:
            entry = await fetch_addgene_entry(addgene_id, client=client)
            if not entry or not entry.gb_url:
                return None
            gb_url = entry.gb_url
        # Seed the auth cookie before hitting media.addgene.org. The
        # cookie's Domain is .addgene.org, so a single warmup on the
        # plasmid page is enough — only needed when the caller didn't
        # already visit an Addgene page through this client.
        if "media.addgene.org" in (gb_url or "") and own_client:
            try:
                await client.get(_ADDGENE_ENTRY.format(aid=addgene_id))
            except Exception:
                pass
        try:
            r = await client.get(gb_url)
            r.raise_for_status()
            text = r.text
            if "ORIGIN" in text and "LOCUS" in text:
                return text
            logger.warning("addgene gb download for %s did not look like GenBank", addgene_id)
            return None
        except Exception as e:
            logger.warning("addgene gb download %s failed: %s", addgene_id, e)
            return None
    finally:
        if own_client:
            await client.aclose()


async def fetch_pubmed_summary(pmid: str, *,
                                client: Optional[httpx.AsyncClient] = None) -> Optional[dict[str, Any]]:
    """Get title + abstract for a PubMed ID via NCBI E-utilities."""
    pid = (pmid or "").strip()
    if not pid.isdigit():
        return None
    own_client = False
    if client is None:
        client = httpx.AsyncClient(timeout=_HTTP_TIMEOUT, headers={"User-Agent": _USER_AGENT})
        own_client = True
    try:
        try:
            sr = await client.get(_PUBMED_ESUMMARY, params={"db": "pubmed", "id": pid, "retmode": "json"})
            sr.raise_for_status()
            sj = sr.json()
            doc = (sj.get("result") or {}).get(pid) or {}
            title = doc.get("title")
        except Exception:
            title = None
        try:
            ar = await client.get(_PUBMED_ABSTRACT, params={"db": "pubmed", "id": pid, "rettype": "abstract",
                                                             "retmode": "text"})
            ar.raise_for_status()
            abstract = ar.text.strip()
        except Exception:
            abstract = None
        return {"pmid": pid, "title": title, "abstract": abstract}
    finally:
        if own_client:
            await client.aclose()


async def fetch_crossref_summary(doi: str, *,
                                  client: Optional[httpx.AsyncClient] = None) -> Optional[dict[str, Any]]:
    """Get title + authors for a DOI via CrossRef."""
    d = (doi or "").strip()
    if not d:
        return None
    own_client = False
    if client is None:
        client = httpx.AsyncClient(timeout=_HTTP_TIMEOUT, headers={"User-Agent": _USER_AGENT})
        own_client = True
    try:
        try:
            r = await client.get(_CROSSREF_WORKS.format(doi=d))
            r.raise_for_status()
            msg = (r.json() or {}).get("message") or {}
            title = " ".join(msg.get("title", [])).strip() or None
            authors = [
                f"{(a.get('given') or '').strip()} {(a.get('family') or '').strip()}".strip()
                for a in (msg.get("author") or [])
            ]
            return {"doi": d, "title": title, "authors": authors[:8]}
        except Exception as e:
            logger.warning("crossref lookup failed for %s: %s", d, e)
            return None
    finally:
        if own_client:
            await client.aclose()


# ─────────────────────────────────────────────────────────────────────
# HTML parsers
# ─────────────────────────────────────────────────────────────────────
_ADDGENE_ID_IN_URL = re.compile(r"/(\d{2,7})/?$|/(\d{2,7})/[^/]")


def _parse_addgene_search(html: str, *, max_results: int) -> list[AddgeneCandidate]:
    """Pull plasmid candidates out of the Addgene search page HTML.

    The page renders each hit inside a `.search-result-item` block (or
    equivalent class — Addgene's HTML has been stable but selectors
    can drift, so we fall back to scanning every anchor whose href
    matches the canonical plasmid-id pattern).
    """
    soup = BeautifulSoup(html, "html.parser")
    seen: set[str] = set()
    out: list[AddgeneCandidate] = []
    # First pass: typed search-result containers if Addgene's HTML still
    # has them. Each container usually has an anchor with the name.
    for block in soup.select(".search-result-item, .search-result, .results .item"):
        a = block.find("a", href=re.compile(r"/\d+/"))
        if not a:
            continue
        aid = _extract_addgene_id(a.get("href") or "")
        if not aid or aid in seen:
            continue
        name = (a.get_text() or "").strip() or f"Addgene #{aid}"
        snippet_el = block.find("p") or block.find(".description")
        snippet = (snippet_el.get_text(" ", strip=True) if snippet_el else None) or None
        out.append(AddgeneCandidate(addgene_id=aid, name=name,
                                     url=_ADDGENE_ENTRY.format(aid=aid), snippet=snippet))
        seen.add(aid)
        if len(out) >= max_results:
            return out
    # Fallback: any anchor pointing at a plasmid-id-shaped path.
    for a in soup.find_all("a", href=re.compile(r"^/\d+/?$|^/\d+/$|/\d+/$")):
        aid = _extract_addgene_id(a.get("href") or "")
        if not aid or aid in seen:
            continue
        name = (a.get_text() or "").strip() or f"Addgene #{aid}"
        out.append(AddgeneCandidate(addgene_id=aid, name=name,
                                     url=_ADDGENE_ENTRY.format(aid=aid), snippet=None))
        seen.add(aid)
        if len(out) >= max_results:
            return out
    return out


def _parse_addgene_entry(html: str, *, aid: str, entry_url: str) -> AddgeneEntry:
    soup = BeautifulSoup(html, "html.parser")
    title_el = soup.select_one("h1") or soup.title
    name = (title_el.get_text(" ", strip=True) if title_el else f"Addgene #{aid}")
    name = re.sub(r"\s+", " ", name).strip()[:200]

    description = None
    for sel in ("#purpose", ".description", ".plasmid-description", "[data-purpose]"):
        el = soup.select_one(sel)
        if el and el.get_text(strip=True):
            description = el.get_text(" ", strip=True)
            break

    depositor = None
    for label in soup.find_all(["strong", "th", "dt"]):
        txt = (label.get_text() or "").strip().rstrip(":").lower()
        if txt in {"depositing lab", "depositor", "principal investigator"}:
            sib = label.find_next_sibling()
            if sib:
                depositor = (sib.get_text(" ", strip=True) or None)
                break
    pmid = None
    doi = None
    pmid_match = re.search(r"PubMed[^0-9]*(\d{5,9})", html)
    if pmid_match:
        pmid = pmid_match.group(1)
    doi_match = re.search(r"\b(10\.\d{4,9}/[\w\.\-_/\(\)]+?)(?:[\"\'\s<]|$)", html)
    if doi_match:
        doi = doi_match.group(1).rstrip(".)\"'")

    return AddgeneEntry(
        addgene_id=aid, name=name, url=entry_url,
        description=description, depositor=depositor,
        pmid=pmid, doi=doi, paper_title=None, gb_url=None,
    )


def _parse_addgene_gb_url(html: str, *, aid: str) -> Optional[str]:
    """Look for a direct `.gbk` / `.gb` link in the supplied page HTML.

    The `/browse/sequence/<seq_id>/` page embeds the URL inside a JS
    payload with `\\u002D` (escaped hyphens) and `\\u002F` (escaped
    slashes). We unescape those first, then run a plain URL regex.
    """
    decoded = _unescape_js(html)
    m = re.search(
        r"https://media\.addgene\.org/[^\"' >]+?\.gbk",
        decoded,
    )
    if m:
        return m.group(0)
    # Older anchor-based fallback layout.
    soup = BeautifulSoup(html, "html.parser")
    for a in soup.find_all("a", href=True):
        href = a["href"]
        if href.lower().endswith((".gbk", ".gb", ".genbank")):
            return _absolutize(_unescape_js(href))
        text = (a.get_text() or "").lower()
        if "genbank" in text and href:
            return _absolutize(_unescape_js(href))
    return None


# Addgene's /<aid>/sequences/ page groups sequences under <h2> headings
# like "Addgene-Verified Full Sequences", "Depositor Full Sequences",
# etc. We prefer the most-canonical full-length sequence available.
_ADDGENE_SEQ_SECTION_PREFS = (
    "addgene-verified full",
    "addgene full",
    "depositor full",
    "addgene-verified partial",
    "depositor partial",
)


def _parse_addgene_seq_ids(html: str) -> list[str]:
    """Pull sequence IDs from the /<aid>/sequences/ page, ordered by
    preference (Addgene-verified full first, partials last)."""
    soup = BeautifulSoup(html, "html.parser")
    # Group seq_ids by the <h2> section they live under.
    groups: dict[str, list[str]] = {}
    current_section = "unknown"
    for el in soup.find_all(["h2", "a"]):
        if el.name == "h2":
            current_section = re.sub(r"\(\d+\)", "", el.get_text(" ", strip=True)).strip().lower()
            current_section = re.sub(r"\s+sequences?$", "", current_section).strip()
        else:
            href = el.get("href", "")
            m = re.match(r"^/browse/sequence/(\d+)/?$", href)
            if m:
                groups.setdefault(current_section, [])
                if m.group(1) not in groups[current_section]:
                    groups[current_section].append(m.group(1))
    # Order sections by preference, then append leftovers.
    ordered: list[str] = []
    for pref in _ADDGENE_SEQ_SECTION_PREFS:
        for sect, ids in list(groups.items()):
            if sect.startswith(pref):
                ordered.extend(i for i in ids if i not in ordered)
    for sect, ids in groups.items():
        ordered.extend(i for i in ids if i not in ordered)
    return ordered


def _unescape_js(s: str) -> str:
    """Convert JS unicode-escapes (`\\u002D`, etc.) to plain chars.

    A full unicode_escape decode would also mangle real UTF-8 in the
    page (e.g. accented characters in depositor names), so we only
    replace the two escapes we actually care about for URL extraction.
    """
    if "\\u" not in s:
        return s
    return (s.replace("\\u002D", "-")
              .replace("\\u002F", "/")
              .replace("\\u003A", ":"))


def _extract_addgene_id(href: str) -> Optional[str]:
    m = re.search(r"^/?(\d{2,7})/?", href)
    if m:
        return m.group(1)
    m = _ADDGENE_ID_IN_URL.search(href)
    if m:
        return m.group(1) or m.group(2)
    return None


def _absolutize(href: str) -> str:
    if href.startswith("http://") or href.startswith("https://"):
        return href
    if href.startswith("/"):
        return "https://www.addgene.org" + href
    return "https://www.addgene.org/" + href


# ─────────────────────────────────────────────────────────────────────
# NCBI Gene / RefSeqGene lookup (for the CRISPR pipeline)
# ─────────────────────────────────────────────────────────────────────
import os as _os

_NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
_NCBI_ESEARCH = f"{_NCBI_BASE}/esearch.fcgi"
_NCBI_EFETCH = f"{_NCBI_BASE}/efetch.fcgi"
_NCBI_ESUMMARY = f"{_NCBI_BASE}/esummary.fcgi"


def _ncbi_api_key() -> Optional[str]:
    """Resolve the NCBI API key. Env wins so production can rotate
    without touching code. Returns None when no key is configured —
    NCBI then enforces 3 req/sec (vs 10 req/sec with a key)."""
    return _os.environ.get("NCBI_API_KEY") or None


def _ncbi_params(**extra) -> dict[str, Any]:
    """Build the base param dict + the api_key when set, + retmode=json."""
    params: dict[str, Any] = {"retmode": "json"}
    key = _ncbi_api_key()
    if key:
        params["api_key"] = key
    params.update(extra)
    return params


@dataclass
class NCBIGeneHit:
    """One hit from an NCBI gene-symbol search.

    For db_source='chromosomal_slice': accession is a chromosome
    record (NC_*), and seq_start/seq_stop define the 1-indexed
    inclusive slice (gene span + flanking_bp on each side). Pass them
    through to fetch_ncbi_genbank to download just the slice.

    For db_source='refseqgene' or 'mrna': accession is the standalone
    record's accessionversion (NG_* or NM_*); slice fields are None.
    """
    gene_id: str           # NCBI Gene ID (numeric, as string)
    accession: str          # nuccore accession (NG_* / NM_* / chromosome NC_*)
    db_source: str          # "refseqgene" | "chromosomal_slice" | "mrna"
    organism: str
    title: str | None = None
    seq_start: int | None = None     # 1-indexed inclusive
    seq_stop: int | None = None      # 1-indexed inclusive
    flanking_bp: int | None = None   # bp added to each side of the gene span
    gene_chr_start: int | None = None  # un-flanked gene span on chromosome (for downstream context)
    gene_chr_stop: int | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "gene_id": self.gene_id, "accession": self.accession,
            "db_source": self.db_source, "organism": self.organism,
            "title": self.title,
            "seq_start": self.seq_start, "seq_stop": self.seq_stop,
            "flanking_bp": self.flanking_bp,
            "gene_chr_start": self.gene_chr_start,
            "gene_chr_stop": self.gene_chr_stop,
        }


async def search_ncbi_gene(
    gene_symbol: str,
    organism: str = "Homo sapiens",
    *,
    client: Optional[httpx.AsyncClient] = None,
) -> Optional[NCBIGeneHit]:
    """Resolve a gene symbol + organism to a single nuccore accession.

    Prefers RefSeqGene (NG_*) records — genomic with introns + flanking,
    what CRISPR / pegRNA design needs. Falls back to RefSeq mRNA (NM_*)
    when no RefSeqGene exists. Returns None on no-hit.
    """
    gene = (gene_symbol or "").strip()
    org = (organism or "Homo sapiens").strip()
    if not gene:
        return None

    own_client = False
    if client is None:
        client = httpx.AsyncClient(timeout=_HTTP_TIMEOUT,
                                    headers={"User-Agent": _USER_AGENT})
        own_client = True
    try:
        # Step 1 — find the NCBI Gene ID for symbol + organism.
        try:
            await _pace_ncbi()
            r = await client.get(
                _NCBI_ESEARCH,
                params=_ncbi_params(db="gene", term=f"{gene}[Gene Name] AND {org}[orgn]",
                                     retmax="3", sort="relevance"),
            )
            r.raise_for_status()
        except httpx.HTTPError as e:
            logger.warning("NCBI gene esearch failed for %r %r: %s", gene, org, e)
            return None
        ids = ((r.json().get("esearchresult") or {}).get("idlist") or [])
        if not ids:
            logger.info("NCBI gene search returned no hits for %r %r", gene, org)
            return None
        gene_id = ids[0]

        # Step 2 — prefer a RefSeqGene (NG_*) record.
        try:
            await _pace_ncbi()
            r = await client.get(
                _NCBI_ESEARCH,
                params=_ncbi_params(
                    db="nuccore",
                    term=(f"{gene}[Gene Name] AND {org}[orgn] AND "
                           "refseqgene[filter]"),
                    retmax="1",
                ),
            )
            r.raise_for_status()
        except httpx.HTTPError as e:
            logger.warning("NCBI refseqgene esearch failed for %r: %s", gene, e)
            return None
        nuc_ids = ((r.json().get("esearchresult") or {}).get("idlist") or [])
        db_source = "refseqgene"

        if not nuc_ids:
            # Step 2b — try a chromosomal slice from the gene's primary
            # assembly placement. esummary against the gene db returns
            # `genomicinfo`: a list of dicts with chraccver / chrstart /
            # chrstop. Guides + pegRNA design need genomic context to
            # handle intron/exon boundaries correctly, so this beats the
            # mRNA fallback below.
            try:
                await _pace_ncbi()
                r = await client.get(
                    _NCBI_ESUMMARY,
                    params=_ncbi_params(db="gene", id=gene_id),
                )
                r.raise_for_status()
            except httpx.HTTPError as e:
                logger.warning("NCBI gene esummary failed for gene_id=%s: %s", gene_id, e)
                return None
            gene_rec = ((r.json().get("result") or {}).get(gene_id) or {})
            genomic_info = gene_rec.get("genomicinfo") or []
            primary = genomic_info[0] if genomic_info else None
            if primary and primary.get("chraccver"):
                chraccver = primary["chraccver"]
                # chrstart / chrstop are 0-indexed in the JSON but can be
                # in either orientation when the gene is on the minus
                # strand (chrstart > chrstop). Normalise to a 1-indexed
                # inclusive [lo, hi] window then flank.
                a = int(primary.get("chrstart") or 0)
                b = int(primary.get("chrstop") or 0)
                lo, hi = (min(a, b), max(a, b))
                # +1 to convert from NCBI's 0-indexed coords to 1-indexed
                # inclusive coords efetch expects.
                gene_start_1based = lo + 1
                gene_stop_1based = hi + 1
                flank = int(_os.environ.get("NCBI_FLANKING_BP", "2000"))
                slice_start = max(1, gene_start_1based - flank)
                slice_stop = gene_stop_1based + flank

                # esummary on the chromosome record to get a readable title.
                # Some chromosomes are millions of bp so we don't want to
                # download the whole record just for the title — esummary
                # returns it cheap.
                title = None
                try:
                    await _pace_ncbi()
                    r = await client.get(
                        _NCBI_ESUMMARY,
                        params=_ncbi_params(db="nuccore", id=chraccver),
                    )
                    r.raise_for_status()
                    res = (r.json().get("result") or {})
                    for k, v in res.items():
                        if k == "uids":
                            continue
                        title = v.get("title") or None
                        break
                except httpx.HTTPError as e:
                    logger.warning("NCBI chromosome esummary failed for %r: %s",
                                     chraccver, e)
                except Exception:
                    pass

                return NCBIGeneHit(
                    gene_id=gene_id, accession=chraccver,
                    db_source="chromosomal_slice", organism=org,
                    title=title or f"{gene} chromosomal locus on {chraccver}",
                    seq_start=slice_start, seq_stop=slice_stop,
                    flanking_bp=flank,
                    gene_chr_start=gene_start_1based,
                    gene_chr_stop=gene_stop_1based,
                )

        if not nuc_ids:
            # Step 2c (final fallback) — RefSeq mRNA (NM_*). Useful only
            # for CDS-internal substitutions; introns + UTRs are absent.
            try:
                await _pace_ncbi()
                r = await client.get(
                    _NCBI_ESEARCH,
                    params=_ncbi_params(
                        db="nuccore",
                        term=(f"{gene}[Gene Name] AND {org}[orgn] AND "
                               "refseq[filter] AND biomol_mrna[prop]"),
                        retmax="1", sort="relevance",
                    ),
                )
                r.raise_for_status()
            except httpx.HTTPError as e:
                logger.warning("NCBI mRNA esearch failed for %r: %s", gene, e)
                return None
            nuc_ids = ((r.json().get("esearchresult") or {}).get("idlist") or [])
            db_source = "mrna"

        if not nuc_ids:
            logger.info("NCBI nuccore search returned no hits for %r %r", gene, org)
            return None

        nuc_id = nuc_ids[0]

        # Step 3 — esummary to get the accession + title.
        try:
            await _pace_ncbi()
            r = await client.get(
                _NCBI_ESUMMARY,
                params=_ncbi_params(db="nuccore", id=nuc_id),
            )
            r.raise_for_status()
        except httpx.HTTPError as e:
            logger.warning("NCBI final esummary failed for %r: %s", nuc_id, e)
            return None
        rec = ((r.json().get("result") or {}).get(nuc_id) or {})
        accession = rec.get("accessionversion") or rec.get("caption") or nuc_id
        title = rec.get("title") or None

        return NCBIGeneHit(
            gene_id=gene_id, accession=accession,
            db_source=db_source, organism=org, title=title,
        )
    finally:
        if own_client:
            await client.aclose()


async def fetch_ncbi_genbank(
    accession: str,
    *,
    client: Optional[httpx.AsyncClient] = None,
    seq_start: int | None = None,
    seq_stop: int | None = None,
) -> Optional[str]:
    """Download a .gb (GenBank flat-file) record by accession.

    When `accession` is a chromosome record (NC_*), pass `seq_start`
    and `seq_stop` (1-indexed inclusive) to fetch just that slice
    rather than the whole chromosome. NCBI returns the slice's features
    in the FEATURES section, intact for the genomic_annotator.

    Returns None on HTTP error or empty body.
    """
    acc = (accession or "").strip()
    if not acc:
        return None
    own_client = False
    if client is None:
        client = httpx.AsyncClient(timeout=60.0,
                                    headers={"User-Agent": _USER_AGENT})
        own_client = True
    try:
        params = {**_ncbi_params(db="nuccore", id=acc,
                                   rettype="gb", retmode="text"),
                   "retmode": "text"}
        if seq_start is not None:
            params["seq_start"] = str(int(seq_start))
        if seq_stop is not None:
            params["seq_stop"] = str(int(seq_stop))
        await _pace_ncbi()
        r = await client.get(_NCBI_EFETCH, params=params)
        r.raise_for_status()
        body = r.text
        if not body or "LOCUS" not in body.splitlines()[0]:
            return None
        return body
    except Exception as e:
        logger.warning("NCBI efetch failed for %r: %s", acc, e)
        return None
    finally:
        if own_client:
            await client.aclose()
