"""Unit tests for the NCBI gene-lookup helpers.

HTTP is fully mocked via httpx.MockTransport — no real network in CI.
The mock seeds NCBI E-utilities JSON payloads matching what eutils
returns for a real gene-symbol search.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))

import json
import pytest
import httpx

from splicify_api.external_search import (
    search_ncbi_gene,
    fetch_ncbi_genbank,
    NCBIGeneHit,
)


# ---------------------------------------------------------------------------
# Mock transport
# ---------------------------------------------------------------------------
def _make_handler(routes):
    """routes: list of (url-substring, response_factory) pairs.
    Each response_factory takes the request and returns an httpx.Response."""
    def handler(request: httpx.Request) -> httpx.Response:
        url = str(request.url)
        for needle, factory in routes:
            if needle in url:
                return factory(request)
        return httpx.Response(404, text=f"no mock for {url}")
    return handler


def test_search_ncbi_gene_returns_refseqgene_when_available():
    import asyncio
    async def _coro():
        """Happy path: CGAS resolves to a Gene ID, then to an NG_* RefSeqGene
        accession, then esummary returns the accessionversion + title."""
        routes = [
            ("/esearch.fcgi", lambda r: _esearch_response(r)),
            ("/esummary.fcgi", lambda r: httpx.Response(200, json={
                "result": {"123456": {"accessionversion": "NG_046896.1",
                                        "caption": "NG_046896",
                                        "title": "Homo sapiens cyclic GMP-AMP synthase (CGAS), RefSeqGene"}},
            })),
        ]
        client = httpx.AsyncClient(transport=httpx.MockTransport(_make_handler(routes)))
        hit = await search_ncbi_gene("CGAS", "Homo sapiens", client=client)
        await client.aclose()

        assert isinstance(hit, NCBIGeneHit)
        assert hit.accession == "NG_046896.1"
        assert hit.db_source == "refseqgene"
        assert hit.organism == "Homo sapiens"
        assert hit.gene_id == "115004"   # from esearch mock
        assert hit.title and "CGAS" in hit.title


    asyncio.run(_coro())

def _esearch_response(request):
    url = str(request.url)
    # Stage A: db=gene -> NCBI Gene ID for the symbol
    if "db=gene" in url and "db%3Dgene" not in url and "db=gene" in url:
        return httpx.Response(200, json={
            "esearchresult": {"idlist": ["115004"], "count": "1"},
        })
    # Stage B: db=nuccore + refseqgene[filter] -> nuccore record ID
    if "db=nuccore" in url and "refseqgene" in url:
        return httpx.Response(200, json={
            "esearchresult": {"idlist": ["123456"], "count": "1"},
        })
    # Stage B fallback: db=nuccore + biomol_mrna -> NM_* ID
    if "db=nuccore" in url and "biomol_mrna" in url:
        return httpx.Response(200, json={
            "esearchresult": {"idlist": ["789012"], "count": "1"},
        })
    return httpx.Response(200, json={"esearchresult": {"idlist": [], "count": "0"}})


def test_search_ncbi_gene_falls_back_to_mrna_when_no_refseqgene():
    import asyncio
    async def _coro():
        """When refseqgene query returns nothing, the search retries with
        biomol_mrna and labels the hit db_source='mrna'."""
        def handler(request: httpx.Request) -> httpx.Response:
            url = str(request.url)
            if "db=gene" in url and "db%3Dgene" not in url and "db=gene" in url:
                return httpx.Response(200, json={"esearchresult": {"idlist": ["1"]}})
            # RefSeqGene query returns empty.
            if "refseqgene" in url:
                return httpx.Response(200, json={"esearchresult": {"idlist": []}})
            # mRNA fallback returns a hit.
            if "biomol_mrna" in url:
                return httpx.Response(200, json={"esearchresult": {"idlist": ["789012"]}})
            if "esummary.fcgi" in url:
                return httpx.Response(200, json={
                    "result": {"789012": {"accessionversion": "NM_138441.3",
                                           "title": "Homo sapiens FAKE mRNA"}},
                })
            return httpx.Response(200, json={"esearchresult": {"idlist": []}})

        client = httpx.AsyncClient(transport=httpx.MockTransport(handler))
        hit = await search_ncbi_gene("FAKE", "Homo sapiens", client=client)
        await client.aclose()

        assert hit is not None
        assert hit.accession == "NM_138441.3"
        assert hit.db_source == "mrna"


    asyncio.run(_coro())

def test_search_ncbi_gene_returns_none_when_no_hits():
    import asyncio
    async def _coro():
        def handler(request: httpx.Request) -> httpx.Response:
            return httpx.Response(200, json={"esearchresult": {"idlist": []}})
        client = httpx.AsyncClient(transport=httpx.MockTransport(handler))
        hit = await search_ncbi_gene("UNOBTANIUM", "Homo sapiens", client=client)
        await client.aclose()
        assert hit is None


    asyncio.run(_coro())

def test_search_ncbi_gene_empty_symbol_returns_none():
    import asyncio
    async def _coro():
        hit = await search_ncbi_gene("", "Homo sapiens")
        assert hit is None


    asyncio.run(_coro())

def test_fetch_ncbi_genbank_returns_text_when_locus_present():
    import asyncio
    async def _coro():
        fake_gb = "LOCUS       NG_046896           50000 bp    DNA     linear   PRI 16-MAY-2026\nDEFINITION  Homo sapiens cyclic GMP-AMP synthase (CGAS).\nFEATURES\n     gene            1..50000\nORIGIN\n        1 acgtacgtac\n//\n"
        def handler(request):
            assert "efetch.fcgi" in str(request.url)
            return httpx.Response(200, text=fake_gb)
        client = httpx.AsyncClient(transport=httpx.MockTransport(handler))
        body = await fetch_ncbi_genbank("NG_046896.1", client=client)
        await client.aclose()
        assert body is not None
        assert body.startswith("LOCUS")
        assert "CGAS" in body


    asyncio.run(_coro())

def test_fetch_ncbi_genbank_returns_none_on_empty_body():
    import asyncio
    async def _coro():
        def handler(request):
            return httpx.Response(200, text="")
        client = httpx.AsyncClient(transport=httpx.MockTransport(handler))
        body = await fetch_ncbi_genbank("NG_NOTHING", client=client)
        await client.aclose()
        assert body is None


    asyncio.run(_coro())

def test_fetch_ncbi_genbank_returns_none_on_http_error():
    import asyncio
    async def _coro():
        def handler(request):
            return httpx.Response(500, text="upstream error")
        client = httpx.AsyncClient(transport=httpx.MockTransport(handler))
        body = await fetch_ncbi_genbank("NG_X", client=client)
        await client.aclose()
        assert body is None
    asyncio.run(_coro())


# ---------------------------------------------------------------------------
# Chromosomal-slice fallback (NEW): preferred over mRNA, falls between
# RefSeqGene and mRNA in the search ladder.
# ---------------------------------------------------------------------------
def test_search_ncbi_gene_chromosomal_slice_fallback():
    """When RefSeqGene is absent but the gene esummary returns genomicinfo,
    the search returns a chromosomal_slice hit with the chromosome
    accession + flanked seq_start/seq_stop coordinates."""
    import asyncio
    async def _coro():
        def handler(request):
            url = str(request.url)
            # Stage A: gene db search -> Gene ID
            if "db=gene" in url and "Gene+Name" in url:
                return httpx.Response(200, json={
                    "esearchresult": {"idlist": ["115004"]},
                })
            # Stage B: RefSeqGene search -> empty
            if "db=nuccore" in url and "refseqgene" in url:
                return httpx.Response(200, json={"esearchresult": {"idlist": []}})
            # Stage 2b: esummary against gene db for genomicinfo
            if "db=gene" in url and "esummary.fcgi" in url:
                return httpx.Response(200, json={
                    "result": {"115004": {
                        "name": "CGAS",
                        "genomicinfo": [{
                            "chraccver": "NC_000006.12",
                            "chrloc": "6",
                            "chrstart": 73852205,   # 0-indexed
                            "chrstop":  73887166,   # 0-indexed
                            "exoncount": 6,
                        }],
                    }},
                })
            # esummary on chromosome record -> title
            if "esummary.fcgi" in url and "NC_000006" in url:
                return httpx.Response(200, json={
                    "result": {"uids": ["568815592"],
                                "568815592": {"title": "Homo sapiens chromosome 6, GRCh38.p14"}},
                })
            return httpx.Response(200, json={"esearchresult": {"idlist": []}})

        client = httpx.AsyncClient(transport=httpx.MockTransport(handler))
        hit = await search_ncbi_gene("CGAS", "Homo sapiens", client=client)
        await client.aclose()

        assert hit is not None
        assert hit.db_source == "chromosomal_slice"
        assert hit.accession == "NC_000006.12"
        # 1-indexed inclusive coords with 2 kb default flanking on each side.
        assert hit.gene_chr_start == 73852206
        assert hit.gene_chr_stop == 73887167
        assert hit.flanking_bp == 2000
        assert hit.seq_start == 73852206 - 2000
        assert hit.seq_stop == 73887167 + 2000
        assert hit.title and ("CGAS" in hit.title or "chromosome" in hit.title.lower())

    asyncio.run(_coro())


def test_chromosomal_slice_handles_minus_strand_swapped_coords():
    """Some gene records present chrstart > chrstop when the gene is on the
    minus strand. Slice math must normalise to lo <= hi before flanking."""
    import asyncio
    async def _coro():
        def handler(request):
            url = str(request.url)
            if "db=gene" in url and "Gene+Name" in url:
                return httpx.Response(200, json={
                    "esearchresult": {"idlist": ["999"]},
                })
            if "db=nuccore" in url and "refseqgene" in url:
                return httpx.Response(200, json={"esearchresult": {"idlist": []}})
            if "db=gene" in url and "esummary.fcgi" in url:
                # minus-strand gene: chrstart > chrstop in the JSON
                return httpx.Response(200, json={
                    "result": {"999": {
                        "genomicinfo": [{
                            "chraccver": "NC_000019.10",
                            "chrstart": 14741,   # higher first -> minus strand
                            "chrstop":  1724,
                        }],
                    }},
                })
            if "esummary.fcgi" in url:
                return httpx.Response(200, json={"result": {"uids": []}})
            return httpx.Response(200, json={"esearchresult": {"idlist": []}})

        client = httpx.AsyncClient(transport=httpx.MockTransport(handler))
        hit = await search_ncbi_gene("MINUS_STRAND_GENE", "Homo sapiens", client=client)
        await client.aclose()
        assert hit is not None
        assert hit.db_source == "chromosomal_slice"
        # Normalised to lo <= hi, then +1 for 1-indexed; flanked by 2000.
        assert hit.gene_chr_start == 1725
        assert hit.gene_chr_stop == 14742
        assert hit.seq_start == max(1, 1725 - 2000)
        assert hit.seq_stop == 14742 + 2000


def test_chromosomal_slice_clamps_seq_start_to_1():
    """A gene at the very start of a chromosome (chrstart near 0) must
    clamp seq_start to >=1 — efetch rejects 0 or negative values."""
    import asyncio
    async def _coro():
        def handler(request):
            url = str(request.url)
            if "db=gene" in url and "Gene+Name" in url:
                return httpx.Response(200, json={"esearchresult": {"idlist": ["111"]}})
            if "db=nuccore" in url and "refseqgene" in url:
                return httpx.Response(200, json={"esearchresult": {"idlist": []}})
            if "db=gene" in url and "esummary.fcgi" in url:
                return httpx.Response(200, json={
                    "result": {"111": {"genomicinfo": [{
                        "chraccver": "NC_TEST.1",
                        "chrstart": 50, "chrstop": 5000,
                    }]}},
                })
            if "esummary.fcgi" in url:
                return httpx.Response(200, json={"result": {"uids": []}})
            return httpx.Response(200, json={"esearchresult": {"idlist": []}})

        client = httpx.AsyncClient(transport=httpx.MockTransport(handler))
        hit = await search_ncbi_gene("TINY", "Homo sapiens", client=client)
        await client.aclose()
        assert hit is not None
        assert hit.seq_start == 1          # max(1, 51 - 2000)
        assert hit.seq_stop == 5001 + 2000


def test_fetch_ncbi_genbank_forwards_seq_start_seq_stop():
    """When the helper is called with slice coords, efetch URL carries them."""
    import asyncio
    async def _coro():
        captured = {}
        def handler(request):
            captured["url"] = str(request.url)
            return httpx.Response(200, text="LOCUS NC_TEST 8000 bp DNA linear PRI 19-MAY-2026\nFEATURES\nORIGIN\n//\n")
        client = httpx.AsyncClient(transport=httpx.MockTransport(handler))
        body = await fetch_ncbi_genbank("NC_000006.12", client=client,
                                          seq_start=73850206, seq_stop=73889167)
        await client.aclose()
        assert body is not None
        assert "seq_start=73850206" in captured["url"]
        assert "seq_stop=73889167" in captured["url"]

    asyncio.run(_coro())
