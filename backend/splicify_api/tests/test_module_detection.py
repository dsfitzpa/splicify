"""
Unit tests for ModuleDetectionPipeline

Tests detector registration, priority-based deduplication, and individual detectors.
"""

import pytest
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from annotation_pipeline.module_detection import (
    ModuleDetectionPipeline,
    ModuleCandidate,
    ModuleDetector,
    OriginDetector,
    MarkerDetector,
)


class MockDetector(ModuleDetector):
    """Mock detector for testing"""

    def __init__(self, name, priority, candidates):
        self._name = name
        self._priority = priority
        self._candidates = candidates

    @property
    def name(self):
        return self._name

    @property
    def priority(self):
        return self._priority

    def detect(self, features, sequence):
        return self._candidates


class TestModuleCandidate:
    """Test ModuleCandidate data class"""

    def test_length_property(self):
        """Test length calculation"""
        candidate = ModuleCandidate(
            module_type='test',
            start=100,
            end=200,
            strand=1,
            priority=50,
            detector_name='test'
        )
        assert candidate.length == 100

    def test_to_module(self):
        """Test conversion to module dict"""
        candidate = ModuleCandidate(
            module_type='pol2_expression',
            start=100,
            end=500,
            strand=1,
            priority=80,
            detector_name='grammar',
            features=['f1', 'f2'],
            metadata={'test': 'value'}
        )

        module = candidate.to_module()

        assert module['module_type'] == 'pol2_expression'
        assert module['start'] == 100
        assert module['end'] == 500
        assert module['strand'] == 1
        assert module['features'] == ['f1', 'f2']
        assert module['metadata']['detector'] == 'grammar'
        assert module['metadata']['priority'] == 80
        assert module['metadata']['test'] == 'value'


class TestModuleDetectionPipeline:
    """Test ModuleDetectionPipeline"""

    def test_detector_registration(self):
        """Test detector registration and sorting by priority"""
        pipeline = ModuleDetectionPipeline()

        detector1 = MockDetector('low', 10, [])
        detector2 = MockDetector('high', 100, [])
        detector3 = MockDetector('med', 50, [])

        pipeline.register_detector(detector1)
        pipeline.register_detector(detector2)
        pipeline.register_detector(detector3)

        # Should be sorted by priority (high to low)
        assert pipeline.detectors[0].name == 'high'
        assert pipeline.detectors[1].name == 'med'
        assert pipeline.detectors[2].name == 'low'

    def test_priority_based_deduplication(self):
        """Test that higher priority modules override lower priority ones"""
        pipeline = ModuleDetectionPipeline()

        # High priority detector finds module at 100-500
        high_priority_candidate = ModuleCandidate(
            module_type='pol2_expression',
            start=100,
            end=500,
            strand=1,
            priority=100,
            detector_name='grammar'
        )

        # Low priority detector finds overlapping module at 150-600
        low_priority_candidate = ModuleCandidate(
            module_type='expression',
            start=150,
            end=600,
            strand=1,
            priority=50,
            detector_name='heuristic'
        )

        pipeline.register_detector(MockDetector('high', 100, [high_priority_candidate]))
        pipeline.register_detector(MockDetector('low', 50, [low_priority_candidate]))

        modules = pipeline.detect_all([], "")

        # Should keep only the high priority module
        assert len(modules) == 1
        assert modules[0]['module_type'] == 'pol2_expression'
        assert modules[0]['metadata']['detector'] == 'grammar'

    def test_non_overlapping_modules_kept(self):
        """Test that non-overlapping modules are all kept"""
        pipeline = ModuleDetectionPipeline()

        candidates = [
            ModuleCandidate(
                module_type='module1',
                start=0,
                end=100,
                strand=1,
                priority=50,
                detector_name='test'
            ),
            ModuleCandidate(
                module_type='module2',
                start=200,
                end=300,
                strand=1,
                priority=50,
                detector_name='test'
            ),
            ModuleCandidate(
                module_type='module3',
                start=400,
                end=500,
                strand=1,
                priority=50,
                detector_name='test'
            ),
        ]

        pipeline.register_detector(MockDetector('test', 50, candidates))
        modules = pipeline.detect_all([], "")

        # All three should be kept (no overlap)
        assert len(modules) == 3

    def test_opposite_strand_not_considered_overlap(self):
        """Test that modules on opposite strands don't conflict"""
        pipeline = ModuleDetectionPipeline()

        candidates = [
            ModuleCandidate(
                module_type='forward',
                start=100,
                end=500,
                strand=1,
                priority=50,
                detector_name='test'
            ),
            ModuleCandidate(
                module_type='reverse',
                start=100,
                end=500,
                strand=-1,
                priority=50,
                detector_name='test'
            ),
        ]

        pipeline.register_detector(MockDetector('test', 50, candidates))
        modules = pipeline.detect_all([], "")

        # Both should be kept (different strands)
        assert len(modules) == 2


class TestOriginDetector:
    """Test OriginDetector"""

    def test_detect_origin_feature(self):
        """Test detection of origin features"""
        detector = OriginDetector()

        features = [
            {
                'id': 'ori1',
                'canonical_type': 'origin',
                'name': 'pUC ori',
                'start': 1000,
                'end': 1500,
                'strand': 1
            }
        ]

        sequence = "A" * 5000

        candidates = detector.detect(features, sequence)

        # Should find origin module
        assert len(candidates) > 0
        assert candidates[0].module_type == 'origin_module'

    def test_origin_type_classification(self):
        """Test classification of origin types"""
        detector = OriginDetector()

        # Test bacterial origin
        features = [{
            'id': 'ori1',
            'canonical_type': 'origin',
            'name': 'pUC ori',
            'start': 1000,
            'end': 1500,
            'strand': 1
        }]

        candidates = detector.detect(features, "A" * 5000)
        assert candidates[0].metadata['origin_type'] == 'bacterial'


class TestMarkerDetector:
    """Test MarkerDetector"""

    def test_detect_marker_with_promoter_and_terminator(self):
        """Test detection of full marker cassette"""
        detector = MarkerDetector()

        features = [
            {
                'id': 'prom1',
                'canonical_type': 'promoter',
                'name': 'amp promoter',
                'start': 1000,
                'end': 1100,
                'strand': 1
            },
            {
                'id': 'amp1',
                'canonical_type': 'cds',
                'name': 'AmpR',
                'start': 1150,
                'end': 2000,
                'strand': 1
            },
            {
                'id': 'term1',
                'canonical_type': 'terminator',
                'name': 'amp terminator',
                'start': 2050,
                'end': 2200,
                'strand': 1
            }
        ]

        sequence = "A" * 5000

        candidates = detector.detect(features, sequence)

        # Should find full cassette
        assert len(candidates) > 0
        full_cassette = [c for c in candidates if c.module_type == 'bacterial_marker_cassette']
        assert len(full_cassette) > 0
        assert full_cassette[0].metadata['has_promoter'] is True
        assert full_cassette[0].metadata['has_terminator'] is True

    def test_detect_standalone_marker(self):
        """Test detection of marker without promoter/terminator"""
        detector = MarkerDetector()

        features = [
            {
                'id': 'amp1',
                'canonical_type': 'cds',
                'name': 'AmpR',
                'start': 1150,
                'end': 2000,
                'strand': 1
            }
        ]

        sequence = "A" * 5000

        candidates = detector.detect(features, sequence)

        # Should find standalone marker
        assert len(candidates) > 0
        assert candidates[0].module_type == 'marker_gene'
        assert candidates[0].metadata['has_promoter'] is False


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
