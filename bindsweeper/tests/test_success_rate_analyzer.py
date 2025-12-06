#!/usr/bin/env python3
"""Tests for success rate analyzer."""

import json
import os
import tempfile
from pathlib import Path

import pytest

from bindsweeper.success_rate_analyzer import SuccessMetrics, SuccessRateAnalyzer


class TestSuccessMetrics:
    """Test SuccessMetrics dataclass."""

    def test_from_dict_full_pipeline(self):
        """Test creating SuccessMetrics from dict with all stages run."""
        data = {
            'parameter_combination': 'test_full',
            'total_designs': 100,
            'successful_designs': 15,
            'success_rate': 0.15,
            'fold_generated': 100,
            'fold_filtered': 95,
            'seq_generated': 100,
            'seq_filtered': 85,
            'pred_filtered': 75,
            'analysis_filtered': 15,
            'timestamp': '2024-01-01T00:00:00Z',
            'pipeline_metrics': {
                'fold_retention_rate': 0.95,
                'seq_retention_rate': 0.85,
                'pred_retention_rate': 0.88,
                'analysis_retention_rate': 0.20,
                'overall_retention_rate': 0.15
            }
        }

        metrics = SuccessMetrics.from_dict(data)

        assert metrics.parameter_combination == 'test_full'
        assert metrics.total_designs == 100
        assert metrics.successful_designs == 15
        assert metrics.success_rate == 0.15
        assert metrics.fold_generated == 100
        assert metrics.fold_filtered == 95
        assert metrics.seq_generated == 100
        assert metrics.seq_filtered == 85
        assert metrics.pred_filtered == 75
        assert metrics.analysis_filtered == 15
        assert metrics.pipeline_metrics['fold_retention_rate'] == 0.95
        assert metrics.pipeline_metrics['seq_retention_rate'] == 0.85
        assert metrics.pipeline_metrics['pred_retention_rate'] == 0.88
        assert metrics.pipeline_metrics['analysis_retention_rate'] == 0.20
        assert metrics.pipeline_metrics['overall_retention_rate'] == 0.15

    def test_from_dict_skip_fold(self):
        """Test creating SuccessMetrics with fold stage skipped (None values)."""
        data = {
            'parameter_combination': 'test_skip_fold',
            'total_designs': 100,
            'successful_designs': 15,
            'success_rate': 0.15,
            'fold_generated': 0,
            'fold_filtered': 0,
            'seq_generated': 100,
            'seq_filtered': 85,
            'pred_filtered': 75,
            'analysis_filtered': 15,
            'timestamp': '2024-01-01T00:00:00Z',
            'pipeline_metrics': {
                'fold_retention_rate': None,
                'seq_retention_rate': 0.85,
                'pred_retention_rate': 0.88,
                'analysis_retention_rate': 0.20,
                'overall_retention_rate': 0.15
            }
        }

        metrics = SuccessMetrics.from_dict(data)

        assert metrics.parameter_combination == 'test_skip_fold'
        assert metrics.total_designs == 100
        assert metrics.successful_designs == 15
        assert metrics.success_rate == 0.15
        assert metrics.fold_generated == 0
        assert metrics.fold_filtered == 0
        assert metrics.seq_generated == 100
        # Fold stage was skipped
        assert metrics.pipeline_metrics['fold_retention_rate'] is None
        # Other stages ran normally
        assert metrics.pipeline_metrics['seq_retention_rate'] == 0.85
        assert metrics.pipeline_metrics['overall_retention_rate'] == 0.15

    def test_from_dict_skip_fold_seq(self):
        """Test creating SuccessMetrics with fold and seq stages skipped."""
        data = {
            'parameter_combination': 'test_skip_fold_seq',
            'total_designs': 75,
            'successful_designs': 15,
            'success_rate': 0.20,
            'fold_generated': 0,
            'fold_filtered': 0,
            'seq_generated': 0,
            'seq_filtered': 75,
            'pred_filtered': 75,
            'analysis_filtered': 15,
            'timestamp': '2024-01-01T00:00:00Z',
            'pipeline_metrics': {
                'fold_retention_rate': None,
                'seq_retention_rate': None,
                'pred_retention_rate': 1.0,
                'analysis_retention_rate': 0.20,
                'overall_retention_rate': 0.20
            }
        }

        metrics = SuccessMetrics.from_dict(data)

        assert metrics.total_designs == 75  # Entry point is filter_seq_count
        assert metrics.successful_designs == 15
        assert metrics.success_rate == 0.20
        # First two stages skipped
        assert metrics.pipeline_metrics['fold_retention_rate'] is None
        assert metrics.pipeline_metrics['seq_retention_rate'] is None
        # Prediction and analysis ran
        assert metrics.pipeline_metrics['pred_retention_rate'] == 1.0
        assert metrics.pipeline_metrics['analysis_retention_rate'] == 0.20

    def test_from_dict_skip_fold_seq_pred(self):
        """Test creating SuccessMetrics with only analysis stage run."""
        data = {
            'parameter_combination': 'test_skip_fold_seq_pred',
            'total_designs': 75,
            'successful_designs': 15,
            'success_rate': 0.20,
            'fold_generated': 0,
            'fold_filtered': 0,
            'seq_generated': 0,
            'seq_filtered': 0,
            'pred_filtered': 75,
            'analysis_filtered': 15,
            'timestamp': '2024-01-01T00:00:00Z',
            'pipeline_metrics': {
                'fold_retention_rate': None,
                'seq_retention_rate': None,
                'pred_retention_rate': None,
                'analysis_retention_rate': 0.20,
                'overall_retention_rate': 0.20
            }
        }

        metrics = SuccessMetrics.from_dict(data)

        assert metrics.total_designs == 75  # Entry point is filter_pred_count
        assert metrics.successful_designs == 15
        assert metrics.success_rate == 0.20
        # Only analysis stage ran
        assert metrics.pipeline_metrics['fold_retention_rate'] is None
        assert metrics.pipeline_metrics['seq_retention_rate'] is None
        assert metrics.pipeline_metrics['pred_retention_rate'] is None
        assert metrics.pipeline_metrics['analysis_retention_rate'] == 0.20

    def test_from_dict_run_fold_only(self):
        """Test creating SuccessMetrics with only fold stage run."""
        data = {
            'parameter_combination': 'test_fold_only',
            'total_designs': 100,
            'successful_designs': 95,
            'success_rate': 0.95,
            'fold_generated': 100,
            'fold_filtered': 95,
            'seq_generated': 0,
            'seq_filtered': 0,
            'pred_filtered': 0,
            'analysis_filtered': 0,
            'timestamp': '2024-01-01T00:00:00Z',
            'pipeline_metrics': {
                'fold_retention_rate': 0.95,
                'seq_retention_rate': None,
                'pred_retention_rate': None,
                'analysis_retention_rate': None,
                'overall_retention_rate': 0.95
            }
        }

        metrics = SuccessMetrics.from_dict(data)

        assert metrics.total_designs == 100
        assert metrics.successful_designs == 95
        assert metrics.success_rate == 0.95
        # Only fold stage ran
        assert metrics.pipeline_metrics['fold_retention_rate'] == 0.95
        assert metrics.pipeline_metrics['seq_retention_rate'] is None
        assert metrics.pipeline_metrics['pred_retention_rate'] is None
        assert metrics.pipeline_metrics['analysis_retention_rate'] is None
        assert metrics.pipeline_metrics['overall_retention_rate'] == 0.95

    def test_to_dict_preserves_none_values(self):
        """Test that to_dict() preserves None values for skipped stages."""
        data = {
            'parameter_combination': 'test',
            'total_designs': 100,
            'successful_designs': 15,
            'success_rate': 0.15,
            'fold_generated': 0,
            'fold_filtered': 0,
            'seq_generated': 100,
            'seq_filtered': 85,
            'pred_filtered': 75,
            'analysis_filtered': 15,
            'timestamp': '2024-01-01T00:00:00Z',
            'pipeline_metrics': {
                'fold_retention_rate': None,
                'seq_retention_rate': 0.85,
                'pred_retention_rate': 0.88,
                'analysis_retention_rate': 0.20,
                'overall_retention_rate': 0.15
            }
        }

        metrics = SuccessMetrics.from_dict(data)
        result_dict = metrics.to_dict()

        assert result_dict['pipeline_metrics']['fold_retention_rate'] is None
        assert result_dict['pipeline_metrics']['seq_retention_rate'] == 0.85
        assert result_dict['pipeline_metrics']['overall_retention_rate'] == 0.15


class TestSuccessRateAnalyzer:
    """Test SuccessRateAnalyzer class."""

    def test_collect_metrics_with_mixed_stages(self, tmp_path):
        """Test collecting metrics from multiple combinations with different skip scenarios."""
        # Create temporary directory structure
        base_dir = tmp_path / "sweep_results"
        base_dir.mkdir()

        # Combination 1: Full pipeline
        combo1_dir = base_dir / "combo1" / "results"
        combo1_dir.mkdir(parents=True)
        metrics1 = {
            'parameter_combination': 'full_pipeline',
            'total_designs': 100,
            'successful_designs': 20,
            'success_rate': 0.20,
            'fold_generated': 100,
            'fold_filtered': 95,
            'seq_generated': 100,
            'seq_filtered': 85,
            'pred_filtered': 75,
            'analysis_filtered': 20,
            'timestamp': '2024-01-01T00:00:00Z',
            'pipeline_metrics': {
                'fold_retention_rate': 0.95,
                'seq_retention_rate': 0.85,
                'pred_retention_rate': 0.88,
                'analysis_retention_rate': 0.27,
                'overall_retention_rate': 0.20
            }
        }
        with open(combo1_dir / "success_metrics.json", 'w') as f:
            json.dump(metrics1, f)

        # Combination 2: Skip fold
        combo2_dir = base_dir / "combo2" / "results"
        combo2_dir.mkdir(parents=True)
        metrics2 = {
            'parameter_combination': 'skip_fold',
            'total_designs': 100,
            'successful_designs': 25,
            'success_rate': 0.25,
            'fold_generated': 0,
            'fold_filtered': 0,
            'seq_generated': 100,
            'seq_filtered': 90,
            'pred_filtered': 80,
            'analysis_filtered': 25,
            'timestamp': '2024-01-01T00:00:00Z',
            'pipeline_metrics': {
                'fold_retention_rate': None,
                'seq_retention_rate': 0.90,
                'pred_retention_rate': 0.89,
                'analysis_retention_rate': 0.31,
                'overall_retention_rate': 0.25
            }
        }
        with open(combo2_dir / "success_metrics.json", 'w') as f:
            json.dump(metrics2, f)

        # Combination 3: Fold only
        combo3_dir = base_dir / "combo3" / "results"
        combo3_dir.mkdir(parents=True)
        metrics3 = {
            'parameter_combination': 'fold_only',
            'total_designs': 100,
            'successful_designs': 98,
            'success_rate': 0.98,
            'fold_generated': 100,
            'fold_filtered': 98,
            'seq_generated': 0,
            'seq_filtered': 0,
            'pred_filtered': 0,
            'analysis_filtered': 0,
            'timestamp': '2024-01-01T00:00:00Z',
            'pipeline_metrics': {
                'fold_retention_rate': 0.98,
                'seq_retention_rate': None,
                'pred_retention_rate': None,
                'analysis_retention_rate': None,
                'overall_retention_rate': 0.98
            }
        }
        with open(combo3_dir / "success_metrics.json", 'w') as f:
            json.dump(metrics3, f)

        # Analyze results
        analyzer = SuccessRateAnalyzer(str(base_dir))
        analyzer.collect_success_metrics([])

        assert len(analyzer.success_metrics) == 3

        # Verify metrics were loaded correctly
        combo_names = [m.parameter_combination for m in analyzer.success_metrics]
        assert 'full_pipeline' in combo_names
        assert 'skip_fold' in combo_names
        assert 'fold_only' in combo_names

        # Check that best combination is fold_only (highest success rate)
        best = max(analyzer.success_metrics, key=lambda x: x.success_rate)
        assert best.parameter_combination == 'fold_only'
        assert best.success_rate == 0.98

    def test_generate_summary_csv_with_none_values(self, tmp_path):
        """Test CSV generation handles None values correctly."""
        base_dir = tmp_path / "sweep_results"
        base_dir.mkdir()

        # Create a combination with skipped stages
        combo_dir = base_dir / "test_combo" / "results"
        combo_dir.mkdir(parents=True)
        metrics = {
            'parameter_combination': 'skip_fold',
            'total_designs': 100,
            'successful_designs': 25,
            'success_rate': 0.25,
            'fold_generated': 0,
            'fold_filtered': 0,
            'seq_generated': 100,
            'seq_filtered': 90,
            'pred_filtered': 80,
            'analysis_filtered': 25,
            'timestamp': '2024-01-01T00:00:00Z',
            'pipeline_metrics': {
                'fold_retention_rate': None,
                'seq_retention_rate': 0.90,
                'pred_retention_rate': 0.89,
                'analysis_retention_rate': 0.31,
                'overall_retention_rate': 0.25
            }
        }
        with open(combo_dir / "success_metrics.json", 'w') as f:
            json.dump(metrics, f)

        # Analyze and generate CSV
        analyzer = SuccessRateAnalyzer(str(base_dir))
        analyzer.collect_success_metrics([])

        csv_path = tmp_path / "summary.csv"
        analyzer.generate_success_summary_csv(str(csv_path))

        # Read CSV and verify
        assert csv_path.exists()
        with open(csv_path, 'r') as f:
            lines = f.readlines()
            assert len(lines) == 2  # Header + 1 data row
            header = lines[0].strip()
            assert 'fold_retention_rate' in header
            assert 'seq_retention_rate' in header
            
            # Check data row - fold_retention_rate should be empty string for None
            data = lines[1].strip()
            # The CSV should have empty string for None, not "nan"
            assert ',,' in data or ',0.9' in data  # Empty field before seq_retention_rate

    def test_zero_success_rate_handling(self):
        """Test that zero success rate is handled correctly (no designs survived)."""
        data = {
            'parameter_combination': 'failed_run',
            'total_designs': 100,
            'successful_designs': 0,
            'success_rate': 0.0,
            'fold_generated': 100,
            'fold_filtered': 95,
            'seq_generated': 100,
            'seq_filtered': 50,
            'pred_filtered': 10,
            'analysis_filtered': 0,
            'timestamp': '2024-01-01T00:00:00Z',
            'pipeline_metrics': {
                'fold_retention_rate': 0.95,
                'seq_retention_rate': 0.50,
                'pred_retention_rate': 0.20,
                'analysis_retention_rate': 0.0,
                'overall_retention_rate': 0.0
            }
        }

        metrics = SuccessMetrics.from_dict(data)

        assert metrics.success_rate == 0.0
        assert metrics.successful_designs == 0
        assert metrics.pipeline_metrics['analysis_retention_rate'] == 0.0
        assert metrics.pipeline_metrics['overall_retention_rate'] == 0.0
