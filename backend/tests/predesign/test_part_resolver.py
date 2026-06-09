"""
Unit tests for part_resolver module.

Tests the PartResolver class and ResolvedPart dataclass.
"""

import pytest
from splicify_api.predesign.design_request import InputSource, PartSpecification
from splicify_api.predesign.part_resolver import PartResolver, ResolvedPart


class TestResolvedPart:
    """Test ResolvedPart dataclass."""

    def test_basic_creation(self):
        """Test creating a basic resolved part."""
        part = ResolvedPart(
            name="Fragment1",
            sequence="ATGCATGC",
            length=8,
            source=InputSource.DIRECT_SEQUENCE,
        )
        assert part.name == "Fragment1"
        assert part.sequence == "ATGCATGC"
        assert part.length == 8
        assert part.source == InputSource.DIRECT_SEQUENCE

    def test_sequence_normalized_to_uppercase(self):
        """Test that sequence is normalized to uppercase."""
        part = ResolvedPart(
            name="test",
            sequence="atgcatgc",
            length=8,
            source=InputSource.DIRECT_SEQUENCE,
        )
        assert part.sequence == "ATGCATGC"

    def test_length_updated_from_sequence(self):
        """Test that length is updated from sequence."""
        part = ResolvedPart(
            name="test",
            sequence="ATGCATGC",
            length=0,  # Will be overridden
            source=InputSource.DIRECT_SEQUENCE,
        )
        assert part.length == 8

    def test_canonical_id_defaults_to_name(self):
        """Test that canonical_id defaults to name if not provided."""
        part = ResolvedPart(
            name="Fragment1",
            sequence="ATGC",
            length=4,
            source=InputSource.DIRECT_SEQUENCE,
        )
        assert part.canonical_id == "Fragment1"

    def test_canonical_id_preserved(self):
        """Test that explicit canonical_id is preserved."""
        part = ResolvedPart(
            name="CMV promoter",
            sequence="ATGC",
            length=4,
            source=InputSource.FEATURE_NAME,
            canonical_id="CMV",
        )
        assert part.canonical_id == "CMV"

    def test_invalid_bases_warning(self):
        """Test that invalid DNA bases generate warnings."""
        part = ResolvedPart(
            name="test",
            sequence="ATGCXYZ",
            length=7,
            source=InputSource.DIRECT_SEQUENCE,
        )
        assert len(part.warnings) > 0
        assert "non-DNA characters" in part.warnings[0]

    def test_to_module_dict(self):
        """Test conversion to module dict format."""
        part = ResolvedPart(
            name="Fragment1",
            sequence="ATGC",
            length=4,
            source=InputSource.DIRECT_SEQUENCE,
            canonical_id="frag1",
            description="Test fragment",
            role="promoter",
            origin="user_input",
            source_detail="direct",
        )
        d = part.to_module_dict()
        assert d["canonical_id"] == "frag1"
        assert d["description"] == "Test fragment"
        assert d["role"] == "promoter"
        assert d["sequence"] == "ATGC"
        assert d["length"] == 4
        assert d["origin"] == "user_input"
        assert d["source"] == "direct"

    def test_to_module_dict_defaults(self):
        """Test module dict with default values."""
        part = ResolvedPart(
            name="Fragment1",
            sequence="ATGC",
            length=4,
            source=InputSource.DIRECT_SEQUENCE,
        )
        d = part.to_module_dict()
        assert d["canonical_id"] == "Fragment1"
        assert d["role"] == "unknown"
        assert d["origin"] == "user_provided"
        assert d["source"] == "user_provided"


class TestPartResolverDirect:
    """Test PartResolver direct sequence resolution."""

    @pytest.fixture
    def resolver(self):
        """Create a PartResolver instance."""
        return PartResolver()

    @pytest.mark.asyncio
    async def test_resolve_direct_sequence(self, resolver):
        """Test resolving a direct sequence."""
        part_spec = PartSpecification(
            name="Fragment1",
            source=InputSource.DIRECT_SEQUENCE,
            sequence="ATGCATGC",
        )
        resolved = resolver._resolve_direct(part_spec)

        assert resolved.name == "Fragment1"
        assert resolved.sequence == "ATGCATGC"
        assert resolved.length == 8
        assert resolved.source == InputSource.DIRECT_SEQUENCE
        assert resolved.confidence == 1.0
        assert resolved.origin == "user_input"

    @pytest.mark.asyncio
    async def test_resolve_direct_with_whitespace(self, resolver):
        """Test that whitespace is removed from sequence."""
        part_spec = PartSpecification(
            name="Fragment1",
            source=InputSource.DIRECT_SEQUENCE,
            sequence="ATGC ATGC\nATGC",
        )
        resolved = resolver._resolve_direct(part_spec)

        assert resolved.sequence == "ATGCATGCATGC"
        assert resolved.length == 12

    @pytest.mark.asyncio
    async def test_resolve_direct_with_role(self, resolver):
        """Test that role is preserved."""
        part_spec = PartSpecification(
            name="CMV",
            source=InputSource.DIRECT_SEQUENCE,
            sequence="ATGC",
            role="promoter",
        )
        resolved = resolver._resolve_direct(part_spec)

        assert resolved.role == "promoter"

    @pytest.mark.asyncio
    async def test_resolve_direct_empty_sequence(self, resolver):
        """Test that empty sequence raises error."""
        part_spec = PartSpecification(
            name="Fragment1",
            source=InputSource.DIRECT_SEQUENCE,
            sequence="",
        )
        with pytest.raises(ValueError, match="sequence is empty"):
            resolver._resolve_direct(part_spec)

    @pytest.mark.asyncio
    async def test_resolve_direct_whitespace_only(self, resolver):
        """Test that whitespace-only sequence raises error."""
        part_spec = PartSpecification(
            name="Fragment1",
            source=InputSource.DIRECT_SEQUENCE,
            sequence="   \n  ",
        )
        with pytest.raises(ValueError, match="sequence is empty"):
            resolver._resolve_direct(part_spec)


class TestPartResolverFile:
    """Test PartResolver file resolution."""

    @pytest.fixture
    def resolver(self):
        """Create a PartResolver instance."""
        return PartResolver()

    def test_detect_format_genbank_extension(self, resolver):
        """Test format detection from .gb extension."""
        assert resolver._detect_format("plasmid.gb", "") == "genbank"
        assert resolver._detect_format("plasmid.gbk", "") == "genbank"
        assert resolver._detect_format("plasmid.genbank", "") == "genbank"

    def test_detect_format_fasta_extension(self, resolver):
        """Test format detection from .fasta extension."""
        assert resolver._detect_format("sequence.fasta", "") == "fasta"
        assert resolver._detect_format("sequence.fa", "") == "fasta"
        assert resolver._detect_format("sequence.fna", "") == "fasta"

    def test_detect_format_genbank_content(self, resolver):
        """Test format detection from LOCUS header."""
        content = "LOCUS       pUC19                   2686 bp    DNA"
        assert resolver._detect_format("unknown.txt", content) == "genbank"

    def test_detect_format_fasta_content(self, resolver):
        """Test format detection from > header."""
        content = ">sequence1\nATGCATGC"
        assert resolver._detect_format("unknown.txt", content) == "fasta"

    def test_detect_format_default(self, resolver):
        """Test default format is genbank."""
        assert resolver._detect_format("unknown.txt", "ATGC") == "genbank"

    @pytest.mark.asyncio
    async def test_resolve_file_genbank(self, resolver):
        """Test resolving a GenBank file."""
        # Simple GenBank format
        genbank_content = """LOCUS       test                      12 bp    DNA     circular
DEFINITION  Test plasmid
FEATURES             Location/Qualifiers
     promoter        1..4
                     /label="test_promoter"
ORIGIN
        1 atgcatgcat gc
//
"""
        part_spec = PartSpecification(
            name="test_plasmid",
            source=InputSource.UPLOADED_FILE,
            file_id="file_123",
        )
        context = {
            "file_cache": {
                "file_123": {
                    "name": "test.gb",
                    "content": genbank_content,
                }
            }
        }

        resolved = await resolver._resolve_file(part_spec, context)

        assert resolved.name == "test_plasmid"
        assert resolved.sequence == "ATGCATGCATGC"
        assert resolved.length == 12
        assert resolved.source == InputSource.UPLOADED_FILE
        assert resolved.origin == "uploaded_file"
        assert resolved.source_detail == "test.gb"
        assert len(resolved.features) > 0

    @pytest.mark.asyncio
    async def test_resolve_file_fasta(self, resolver):
        """Test resolving a FASTA file."""
        fasta_content = """>test_sequence
ATGCATGCATGC
"""
        part_spec = PartSpecification(
            name="test_seq",
            source=InputSource.UPLOADED_FILE,
            file_id="file_456",
        )
        context = {
            "file_cache": {
                "file_456": {
                    "name": "test.fasta",
                    "content": fasta_content,
                }
            }
        }

        resolved = await resolver._resolve_file(part_spec, context)

        assert resolved.name == "test_seq"
        assert resolved.sequence == "ATGCATGCATGC"
        assert resolved.length == 12

    @pytest.mark.asyncio
    async def test_resolve_file_missing_cache(self, resolver):
        """Test that missing file_id raises error."""
        part_spec = PartSpecification(
            name="test",
            source=InputSource.UPLOADED_FILE,
            file_id="missing_file",
        )
        context = {"file_cache": {}}

        with pytest.raises(ValueError, match="not found in cache"):
            await resolver._resolve_file(part_spec, context)

    @pytest.mark.asyncio
    async def test_resolve_file_empty(self, resolver):
        """Test that empty file raises error."""
        part_spec = PartSpecification(
            name="test",
            source=InputSource.UPLOADED_FILE,
            file_id="file_empty",
        )
        context = {
            "file_cache": {
                "file_empty": {
                    "name": "empty.gb",
                    "content": "",
                }
            }
        }

        with pytest.raises(ValueError, match="No sequences found"):
            await resolver._resolve_file(part_spec, context)


class TestPartResolverAll:
    """Test PartResolver.resolve_all() method."""

    @pytest.fixture
    def resolver(self):
        """Create a PartResolver instance."""
        return PartResolver()

    @pytest.mark.asyncio
    async def test_resolve_all_direct_sequences(self, resolver):
        """Test resolving multiple direct sequences."""
        parts = [
            PartSpecification(
                name="Part1",
                source=InputSource.DIRECT_SEQUENCE,
                sequence="ATGC",
            ),
            PartSpecification(
                name="Part2",
                source=InputSource.DIRECT_SEQUENCE,
                sequence="GCTA",
            ),
        ]

        resolved = await resolver.resolve_all(parts)

        assert len(resolved) == 2
        assert resolved[0].name == "Part1"
        assert resolved[0].sequence == "ATGC"
        assert resolved[1].name == "Part2"
        assert resolved[1].sequence == "GCTA"

    @pytest.mark.asyncio
    async def test_resolve_all_mixed_sources(self, resolver):
        """Test resolving parts from multiple sources."""
        genbank_content = """LOCUS       test                       4 bp    DNA
ORIGIN
        1 atgc
//
"""
        parts = [
            PartSpecification(
                name="Direct",
                source=InputSource.DIRECT_SEQUENCE,
                sequence="ATGC",
            ),
            PartSpecification(
                name="File",
                source=InputSource.UPLOADED_FILE,
                file_id="file_1",
            ),
        ]
        context = {
            "file_cache": {
                "file_1": {
                    "name": "test.gb",
                    "content": genbank_content,
                }
            }
        }

        resolved = await resolver.resolve_all(parts, context)

        assert len(resolved) == 2
        assert resolved[0].source == InputSource.DIRECT_SEQUENCE
        assert resolved[1].source == InputSource.UPLOADED_FILE

    @pytest.mark.asyncio
    async def test_resolve_all_preserves_order(self, resolver):
        """Test that resolve_all preserves part order."""
        parts = [
            PartSpecification(name="Part3", source=InputSource.DIRECT_SEQUENCE, sequence="AAAA"),
            PartSpecification(name="Part1", source=InputSource.DIRECT_SEQUENCE, sequence="CCCC"),
            PartSpecification(name="Part2", source=InputSource.DIRECT_SEQUENCE, sequence="GGGG"),
        ]

        resolved = await resolver.resolve_all(parts)

        assert [p.name for p in resolved] == ["Part3", "Part1", "Part2"]

    @pytest.mark.asyncio
    async def test_resolve_all_failure_raises(self, resolver):
        """Test that resolution failure raises error."""
        parts = [
            PartSpecification(
                name="Invalid",
                source=InputSource.DIRECT_SEQUENCE,
                sequence="",  # Empty sequence
            ),
        ]

        with pytest.raises(ValueError, match="Failed to resolve part"):
            await resolver.resolve_all(parts)
