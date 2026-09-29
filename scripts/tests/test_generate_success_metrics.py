"""Tests for generate_success_metrics.py's pure calculation functions."""
import argparse

import pytest

from generate_success_metrics import (
    calculate_success_rate,
    calculate_overall_success_rate,
    generate_success_metrics,
)


class TestCalculateSuccessRate:
    def test_normal_ratio(self):
        assert calculate_success_rate(5, 10) == 0.5

    def test_zero_total_returns_zero(self):
        assert calculate_success_rate(0, 0) == 0.0

    def test_zero_successful(self):
        assert calculate_success_rate(0, 10) == 0.0

    def test_all_successful(self):
        assert calculate_success_rate(10, 10) == 1.0


def _base_args(**overrides):
    """Full set of args generate_success_metrics()/calculate_overall_success_rate() reads,
    defaulted to 0/None so tests only need to override the fields relevant to the scenario."""
    defaults = dict(
        fold_count=0,
        filter_fold_count=0,
        seq_count=0,
        filter_seq_count=0,
        pred_count=0,
        filter_pred_count=0,
        filter_analysis_count=0,
        final_designs_count=0,
        af2_count=0,
        af2_filter_count=0,
        boltz_count=0,
        boltz_filter_count=0,
        parameter_combination="unknown",
    )
    defaults.update(overrides)
    return argparse.Namespace(**defaults)


class TestCalculateOverallSuccessRate:
    def test_full_pipeline_uses_seq_count_as_denominator(self):
        args = _base_args(fold_count=100, seq_count=50, pred_count=50, final_designs_count=5)
        rate, total = calculate_overall_success_rate(args)
        assert total == 50
        assert rate == pytest.approx(0.1)

    def test_skip_fold_seq_uses_pred_count(self):
        args = _base_args(pred_count=20, final_designs_count=4)
        rate, total = calculate_overall_success_rate(args)
        assert total == 20
        assert rate == pytest.approx(0.2)

    def test_skip_fold_seq_pred_uses_filter_pred_count(self):
        args = _base_args(filter_pred_count=8, final_designs_count=2)
        rate, total = calculate_overall_success_rate(args)
        assert total == 8
        assert rate == pytest.approx(0.25)

    def test_fold_only_uses_fold_count(self):
        args = _base_args(fold_count=40, final_designs_count=10)
        rate, total = calculate_overall_success_rate(args)
        assert total == 40
        assert rate == pytest.approx(0.25)

    def test_no_valid_entry_point_raises(self):
        args = _base_args()
        with pytest.raises(ValueError, match="No valid pipeline entry point"):
            calculate_overall_success_rate(args)

    def test_seq_count_takes_precedence_over_pred_count(self):
        # seq_count > 0 should win even if pred_count is also set (full pipeline case)
        args = _base_args(seq_count=30, pred_count=99, final_designs_count=3)
        _, total = calculate_overall_success_rate(args)
        assert total == 30

    def test_af2_boltz_all_filtered_uses_af2_count(self):
        # Regression: all designs failing AF2 filtering left pred_count=0, previously raising
        args = _base_args(af2_count=10, af2_filter_count=0, boltz_count=0,
                           boltz_filter_count=0, final_designs_count=0)
        rate, total = calculate_overall_success_rate(args)
        assert total == 10
        assert rate == 0.0

    def test_af2_boltz_partial_filter_uses_af2_count_not_pred_count(self):
        # af2_count is the true entry point; pred_count only reflects the Boltz sub-stage
        # and must not be used as the denominator when af2_count is available
        args = _base_args(af2_count=10, af2_filter_count=6, pred_count=6, final_designs_count=2)
        rate, total = calculate_overall_success_rate(args)
        assert total == 10
        assert rate == pytest.approx(0.2)

    def test_seq_count_still_takes_precedence_over_af2_count(self):
        args = _base_args(seq_count=10, af2_count=10, final_designs_count=1)
        _, total = calculate_overall_success_rate(args)
        assert total == 10


class TestGenerateSuccessMetrics:
    def test_stages_not_run_report_none_retention_rate(self):
        # Entry point via pred_count only; fold/seq stages were skipped
        args = _base_args(pred_count=10, filter_pred_count=8, filter_analysis_count=6,
                           final_designs_count=4)
        metrics = generate_success_metrics(args)
        assert metrics["pipeline_metrics"]["fold_retention_rate"] is None
        assert metrics["pipeline_metrics"]["seq_retention_rate"] is None
        assert metrics["pipeline_metrics"]["pred_retention_rate"] == pytest.approx(0.8)
        assert metrics["pipeline_metrics"]["analysis_retention_rate"] == pytest.approx(0.75)

    def test_af2_boltz_cascade_retention_rates_present(self):
        args = _base_args(seq_count=10, pred_count=10, af2_count=10, af2_filter_count=5,
                           boltz_count=5, boltz_filter_count=2, final_designs_count=2)
        metrics = generate_success_metrics(args)
        assert metrics["pipeline_metrics"]["af2_retention_rate"] == pytest.approx(0.5)
        assert metrics["pipeline_metrics"]["boltz_retention_rate"] == pytest.approx(0.4)

    def test_af2_boltz_retention_none_when_not_used(self):
        args = _base_args(seq_count=10, pred_count=10, final_designs_count=1)
        metrics = generate_success_metrics(args)
        assert metrics["pipeline_metrics"]["af2_retention_rate"] is None
        assert metrics["pipeline_metrics"]["boltz_retention_rate"] is None

    def test_overall_retention_rate_rounded_to_4dp(self):
        args = _base_args(seq_count=3, final_designs_count=1)
        metrics = generate_success_metrics(args)
        assert metrics["pipeline_metrics"]["overall_retention_rate"] == round(1 / 3, 4)
        assert metrics["success_rate"] == round(1 / 3, 4)

    def test_parameter_combination_passthrough(self):
        args = _base_args(seq_count=1, final_designs_count=1, parameter_combination="combo_7")
        metrics = generate_success_metrics(args)
        assert metrics["parameter_combination"] == "combo_7"

    def test_timestamp_is_iso_format_string(self):
        args = _base_args(seq_count=1, final_designs_count=1)
        metrics = generate_success_metrics(args)
        # Just confirm it round-trips as an ISO 8601 string, not exact value (uses real clock)
        from datetime import datetime
        datetime.fromisoformat(metrics["timestamp"])

    def test_af2_boltz_all_filtered_does_not_raise(self):
        # End-to-end reproduction of the reported crash: skip_fold_seq=true, af2_boltz,
        # every design fails AF2 filtering so the Boltz stage never runs
        args = _base_args(af2_count=10, af2_filter_count=0, boltz_count=0,
                           boltz_filter_count=0, final_designs_count=0)
        metrics = generate_success_metrics(args)
        assert metrics["total_designs"] == 10
        assert metrics["pipeline_metrics"]["af2_retention_rate"] == 0.0
        assert metrics["pipeline_metrics"]["boltz_retention_rate"] is None
