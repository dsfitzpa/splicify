"""
Unit tests for ValidationPipeline

Tests each checker and the overall pipeline integration.
"""

import pytest
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from annotation_pipeline.validation_pipeline import (
    ValidationPipeline,
    ValidationIssue,
    ComponentChecker,
    DependencyChecker,
    OrientationChecker,
    FrameChecker,
)


class TestComponentChecker:
    """Test ComponentChecker"""

    def test_missing_required_component(self):
        """Test detection of missing required component"""
        schema = {
            'pol2_expression': {
                'promoter': {'required': True, 'multiple': True},
                'orf': {'required': True, 'multiple': False},
            }
        }

        checker = ComponentChecker(schema)

        # Module missing ORF
        modules = [{
            'id': 'mod1',
            'module_type': 'pol2_expression',
            'start': 0,
            'end': 1000,
            'strand': 1
        }]

        # Only has promoter
        features = [{
            'id': 'feat1',
            'canonical_type': 'promoter',
            'start': 0,
            'end': 100,
            'strand': 1
        }]

        issues = checker.check(features, modules)

        # Should have issue for missing ORF
        assert len(issues) > 0
        assert any('orf' in issue.message.lower() for issue in issues)

    def test_multiple_not_allowed(self):
        """Test detection of multiple components when only one allowed"""
        schema = {
            'pol2_expression': {
                'orf': {'required': True, 'multiple': False},
            }
        }

        checker = ComponentChecker(schema)

        modules = [{
            'id': 'mod1',
            'module_type': 'pol2_expression',
            'start': 0,
            'end': 1000,
            'strand': 1
        }]

        # Two ORFs
        features = [
            {
                'id': 'feat1',
                'canonical_type': 'orf',
                'start': 0,
                'end': 300,
                'strand': 1
            },
            {
                'id': 'feat2',
                'canonical_type': 'orf',
                'start': 400,
                'end': 700,
                'strand': 1
            }
        ]

        issues = checker.check(features, modules)

        # Should have warning for multiple ORFs
        assert len(issues) > 0
        assert any(issue.severity == 'warning' for issue in issues)


class TestDependencyChecker:
    """Test DependencyChecker"""

    def test_required_dependency_missing(self):
        """Test detection of missing required dependency"""
        dependencies = [
            {
                'source': 'lentiviral_payload',
                'target': 'lentiviral_cis_module',
                'relationship': 'requires',
                'severity': 'error',
                'description': 'Test dependency'
            }
        ]

        checker = DependencyChecker(dependencies)

        # Has payload but no cis module
        modules = [{
            'id': 'mod1',
            'module_type': 'lentiviral_payload',
            'start': 0,
            'end': 1000,
            'strand': 1
        }]

        issues = checker.check([], modules)

        # Should have error for missing cis module
        assert len(issues) > 0
        assert issues[0].severity == 'error'
        assert 'lentiviral_cis_module' in issues[0].message

    def test_prefers_dependency(self):
        """Test prefers relationship (info level)"""
        dependencies = [
            {
                'source': 'shuttle_vector_backbone',
                'target': 'yeast_pol2_expression',
                'relationship': 'prefers',
                'severity': 'info',
                'description': 'Test preference'
            }
        ]

        checker = DependencyChecker(dependencies)

        modules = [{
            'id': 'mod1',
            'module_type': 'shuttle_vector_backbone',
            'start': 0,
            'end': 1000,
            'strand': 1
        }]

        issues = checker.check([], modules)

        # Should have info-level suggestion
        assert len(issues) > 0
        assert issues[0].severity == 'info'


class TestOrientationChecker:
    """Test OrientationChecker"""

    def test_orientation_mismatch(self):
        """Test detection of orientation mismatches"""
        checker = OrientationChecker()

        modules = [{
            'id': 'mod1',
            'module_type': 'pol2_expression',
            'start': 0,
            'end': 1000,
            'strand': 1  # Forward strand
        }]

        # Features with mixed orientations
        features = [
            {
                'id': 'feat1',
                'start': 0,
                'end': 300,
                'strand': 1  # Correct
            },
            {
                'id': 'feat2',
                'start': 400,
                'end': 700,
                'strand': -1  # Wrong!
            }
        ]

        issues = checker.check(features, modules)

        # Should detect misoriented feature
        assert len(issues) > 0
        assert 'feat2' in issues[0].feature_ids


class TestFrameChecker:
    """Test FrameChecker"""

    def test_frame_not_multiple_of_3(self):
        """Test detection of CDS not multiple of 3"""
        checker = FrameChecker()

        # CDS with length 100 (not divisible by 3)
        features = [{
            'id': 'cds1',
            'canonical_type': 'cds',
            'start': 0,
            'end': 100,  # Length = 100
            'strand': 1
        }]

        issues = checker.check(features, [])

        # Should have warning
        assert len(issues) > 0
        assert '100' in issues[0].message

    def test_frame_correct(self):
        """Test that correct CDS length passes"""
        checker = FrameChecker()

        # CDS with length 99 (divisible by 3)
        features = [{
            'id': 'cds1',
            'canonical_type': 'cds',
            'start': 0,
            'end': 99,  # Length = 99
            'strand': 1
        }]

        issues = checker.check(features, [])

        # Should have no issues
        assert len(issues) == 0


class TestValidationPipeline:
    """Test ValidationPipeline integration"""

    def test_pipeline_runs_all_checkers(self):
        """Test that pipeline runs all checkers"""
        # Create pipeline without YAML files (uses empty defaults)
        pipeline = ValidationPipeline(rules_dir=None)

        # Should have all 4 checkers registered
        assert len(pipeline.checkers) == 4

        # Run validation (should not crash)
        features = []
        modules = []
        issues = pipeline.validate(features, modules)

        # Should return list (may be empty)
        assert isinstance(issues, list)

    def test_pipeline_to_dict(self):
        """Test validation output as dictionaries"""
        pipeline = ValidationPipeline(rules_dir=None)

        features = []
        modules = []
        issues_dict = pipeline.validate_to_dict(features, modules)

        # Should return list of dicts
        assert isinstance(issues_dict, list)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
