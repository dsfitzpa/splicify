"""
Unit tests for CDSProcessingPipeline

Tests each stage and pipeline integration.
"""

import pytest
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from annotation_pipeline.cds_processing import (
    CDSProcessingPipeline,
    ProcessingContext,
    IntronSplicingStage,
    CodonAlignmentStage,
    SubmoduleDetectionStage,
    GapAnalysisStage,
    BoundaryCorrectionStage,
)


class TestIntronSplicingStage:
    """Test IntronSplicingStage"""

    def test_no_introns(self):
        """Test CDS without introns"""
        stage = IntronSplicingStage()

        # Simple CDS: ATG + 30 codons + TAA = 93 bp
        sequence = "ATG" + "GCT" * 30 + "TAA"

        cds_module = {
            'start': 0,
            'end': len(sequence),
            'strand': 1
        }

        ctx = ProcessingContext(
            cds_module=cds_module,
            features=[],
            sequence=sequence
        )

        ctx = stage.execute(ctx)

        # Should have spliced sequence equal to original
        assert ctx.spliced_sequence == sequence
        assert len(ctx.position_map) == len(sequence)
        assert ctx.introns == []

    def test_with_introns(self):
        """Test CDS with introns"""
        stage = IntronSplicingStage()

        # Exon1 (30bp) + Intron (50bp) + Exon2 (30bp) = 110bp total
        exon1 = "ATG" + "GCT" * 9  # 30 bp
        intron = "G" * 50
        exon2 = "CAG" * 9 + "TAA"  # 30 bp
        sequence = exon1 + intron + exon2

        cds_module = {
            'start': 0,
            'end': len(sequence),
            'strand': 1
        }

        intron_feature = {
            'canonical_type': 'intron',
            'start': len(exon1),
            'end': len(exon1) + len(intron),
            'strand': 1
        }

        ctx = ProcessingContext(
            cds_module=cds_module,
            features=[intron_feature],
            sequence=sequence
        )

        ctx = stage.execute(ctx)

        # Should have spliced out intron
        expected_spliced = exon1 + exon2
        assert ctx.spliced_sequence == expected_spliced
        assert len(ctx.introns) == 1


class TestCodonAlignmentStage:
    """Test CodonAlignmentStage"""

    def test_feature_alignment(self):
        """Test alignment of features to codon boundaries"""
        stage = CodonAlignmentStage()

        sequence = "ATG" + "GCT" * 30 + "TAA"  # 93 bp

        cds_module = {
            'start': 0,
            'end': len(sequence),
            'strand': 1
        }

        # Feature not aligned to codon (starts at position 5, should snap to 3 or 6)
        feature = {
            'id': 'feat1',
            'name': 'protein',
            'canonical_type': 'cds',
            'start': 5,
            'end': 35,
            'strand': 1
        }

        ctx = ProcessingContext(
            cds_module=cds_module,
            features=[feature],
            sequence=sequence,
            spliced_sequence=sequence
        )

        ctx = stage.execute(ctx)

        # Should have boundary correction
        assert len(ctx.boundary_corrections) > 0
        correction = ctx.boundary_corrections[0]
        assert correction['original_start'] == 5
        assert correction['reason'] == 'codon_alignment'

        # Aligned features should exist
        assert len(ctx.aligned_features) > 0


class TestSubmoduleDetectionStage:
    """Test SubmoduleDetectionStage"""

    def test_protein_categorization(self):
        """Test categorization of protein features"""
        stage = SubmoduleDetectionStage()

        sequence = "ATG" + "GCT" * 30 + "TAA"

        cds_module = {
            'start': 0,
            'end': len(sequence),
            'strand': 1
        }

        features = [
            {
                'id': 'prot1',
                'name': 'GFP',
                'canonical_type': 'cds',
                'start': 0,
                'end': 90,
                'strand': 1
            }
        ]

        ctx = ProcessingContext(
            cds_module=cds_module,
            features=features,
            sequence=sequence,
            aligned_features=features
        )

        ctx = stage.execute(ctx)

        # Should detect protein submodule
        assert len(ctx.submodules) > 0
        assert ctx.submodules[0]['type'] == 'protein_module'

    def test_linker_categorization(self):
        """Test categorization of linker features"""
        stage = SubmoduleDetectionStage()

        sequence = "ATG" + "GCT" * 30 + "TAA"

        cds_module = {
            'start': 0,
            'end': len(sequence),
            'strand': 1
        }

        features = [
            {
                'id': 'link1',
                'name': 'GSG linker',
                'canonical_type': 'linker',
                'start': 30,
                'end': 45,
                'strand': 1
            }
        ]

        ctx = ProcessingContext(
            cds_module=cds_module,
            features=features,
            sequence=sequence,
            aligned_features=features
        )

        ctx = stage.execute(ctx)

        # Should detect linker submodule
        assert len(ctx.submodules) > 0
        assert ctx.submodules[0]['type'] == 'linker_module'

    def test_coverage_filtering(self):
        """Test >90% coverage filtering"""
        stage = SubmoduleDetectionStage()

        sequence = "ATG" + "GCT" * 30 + "TAA"

        cds_module = {
            'start': 0,
            'end': len(sequence),
            'strand': 1
        }

        # Two overlapping proteins: one large, one small (95% covered)
        features = [
            {
                'id': 'prot1',
                'name': 'Large protein',
                'canonical_type': 'cds',
                'start': 0,
                'end': 90,
                'strand': 1
            },
            {
                'id': 'prot2',
                'name': 'Small protein',
                'canonical_type': 'cds',
                'start': 0,
                'end': 60,  # 60/90 = 66% of large, but large covers 100% of small
                'strand': 1
            }
        ]

        ctx = ProcessingContext(
            cds_module=cds_module,
            features=features,
            sequence=sequence,
            aligned_features=features
        )

        ctx = stage.execute(ctx)

        # Should keep only one protein (larger one)
        protein_submodules = [s for s in ctx.submodules if s['type'] == 'protein_module']
        assert len(protein_submodules) == 1

        # Smaller protein should be in filtered list
        assert 'Small protein' in ctx.filtered_features


class TestGapAnalysisStage:
    """Test GapAnalysisStage"""

    def test_gap_detection(self):
        """Test detection of gaps between submodules"""
        stage = GapAnalysisStage()

        sequence = "ATG" + "GCT" * 50 + "TAA"

        cds_module = {
            'start': 0,
            'end': len(sequence),
            'strand': 1
        }

        # Create context with two submodules with gap between them
        ctx = ProcessingContext(
            cds_module=cds_module,
            features=[],
            sequence=sequence,
            submodules=[
                {
                    'type': 'protein_module',
                    'name': 'Protein1',
                    'start': 0,
                    'end': 60,
                    'strand': 1
                },
                {
                    'type': 'protein_module',
                    'name': 'Protein2',
                    'start': 100,
                    'end': 150,
                    'strand': 1
                }
            ]
        )

        ctx = stage.execute(ctx)

        # Should detect gap
        assert len(ctx.gaps) > 0
        assert ctx.gaps[0]['start'] == 60
        assert ctx.gaps[0]['end'] == 100
        assert ctx.gaps[0]['length'] == 40


class TestBoundaryCorrectionStage:
    """Test BoundaryCorrectionStage"""

    def test_boundary_correction_outside_cds(self):
        """Test correction of submodules outside CDS bounds"""
        stage = BoundaryCorrectionStage()

        cds_module = {
            'start': 100,
            'end': 500,
            'strand': 1
        }

        # Submodule extending outside CDS
        ctx = ProcessingContext(
            cds_module=cds_module,
            features=[],
            sequence="A" * 1000,
            submodules=[
                {
                    'type': 'protein_module',
                    'name': 'Protein1',
                    'start': 50,  # Before CDS start
                    'end': 200,
                    'strand': 1
                }
            ]
        )

        ctx = stage.execute(ctx)

        # Should correct start to CDS start
        assert ctx.submodules[0]['start'] == 100
        assert len(ctx.boundary_corrections) > 0


class TestCDSProcessingPipeline:
    """Test CDSProcessingPipeline integration"""

    def test_pipeline_execution(self):
        """Test full pipeline execution"""
        pipeline = CDSProcessingPipeline()

        sequence = "ATG" + "GCT" * 30 + "TAA"

        cds_module = {
            'start': 0,
            'end': len(sequence),
            'strand': 1
        }

        features = [
            {
                'id': 'prot1',
                'name': 'GFP',
                'canonical_type': 'cds',
                'start': 0,
                'end': 90,
                'strand': 1
            }
        ]

        result = pipeline.process(cds_module, features, sequence)

        # Should return CDSResult
        assert result.submodules is not None
        assert result.filtered_features is not None
        assert result.boundary_corrections is not None
        assert result.processing_log is not None

        # Should have processing log entries
        assert len(result.processing_log) > 0

    def test_pipeline_stages_order(self):
        """Test that pipeline has correct stages in order"""
        pipeline = CDSProcessingPipeline()

        # Should have 5 stages
        assert len(pipeline.stages) == 5

        # Check order
        assert isinstance(pipeline.stages[0], IntronSplicingStage)
        assert isinstance(pipeline.stages[1], CodonAlignmentStage)
        assert isinstance(pipeline.stages[2], SubmoduleDetectionStage)
        assert isinstance(pipeline.stages[3], GapAnalysisStage)
        assert isinstance(pipeline.stages[4], BoundaryCorrectionStage)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
