"""核对合成对照协议、KKT、lambda 留出选择与空间 placebo 诊断。"""
import json
import math
import subprocess
import sys
import unittest

from analysis_data.io import load
from analysis_data.synthetic import analyze, fit_weights, validate


def _analyze(path):
    return analyze(load(path)[0], source_sha256=load(path)[1])


class SyntheticControlTests(unittest.TestCase):
    def test_invalid_panel_rejected(self):
        with self.assertRaises(ValueError):
            analyze(load("data/invalid-unbalanced.json")[0])

    def test_legal_weights_sum_and_kkt(self):
        result = _analyze("data/normal.json")
        weights = [entry["weight"] for entry in result["fit"]["weights"]]
        self.assertTrue(all(weight >= 0.0 for weight in weights))
        self.assertAlmostEqual(sum(weights), 1.0, places=12)
        self.assertEqual(result["fit"]["solver"]["optimal"], True)
        self.assertLessEqual(result["fit"]["kkt"]["max_residual"], 1e-8)
        self.assertIn("objective", result["lambda_selection"]["candidates"][0])
        # donor-a 与 donor-c 故意共线；仍须给出合法解而不是逆矩阵失败。
        self.assertTrue(all(math.isfinite(weight) for weight in weights))

    def test_objective_matches_definition(self):
        result = _analyze("data/hand-example.json")
        fit = result["fit"]
        self.assertAlmostEqual(
            fit["objective"],
            fit["training_mse"] + fit["lambda"] * sum(entry["weight"] ** 2 for entry in fit["weights"]),
            places=12,
        )

    def test_kkt_gradient_is_independently_checkable(self):
        value = load("data/normal.json")[0]
        result = analyze(value)
        fit = result["fit"]
        donors = [site for site in value["sites"] if site != value["treated"]]
        pre = [week for week in value["weeks"] if week < value["treatment_week"]]
        weights = fit["donor_weights"]
        lam = fit["lambda"]
        treated = {row["week"]: row["value"] for row in value["rows"] if row["site"] == "treated"}
        donor_values = {
            site: {row["week"]: row["value"] for row in value["rows"] if row["site"] == site}
            for site in donors
        }
        length = len(pre)
        gradient = []
        for donor in donors:
            total = 0.0
            for week in pre:
                synthetic = sum(weights[name] * donor_values[name][week] for name in donors)
                total += donor_values[donor][week] * (synthetic - treated[week])
            gradient.append(2.0 / length * total + 2.0 * lam * weights[donor])
        for reported, independent in zip(fit["gradient"], gradient):
            self.assertAlmostEqual(reported, independent, places=10)
        # mu = gradient + nu 是输出的非负约束乘子。
        for reported, value_ in zip(fit["mu"], gradient):
            self.assertAlmostEqual(reported, value_ + fit["equality_multiplier"], places=10)

    def test_duplicate_donor_columns(self):
        fitted = fit_weights([[1.0, 2.0, 3.0], [1.0, 2.0, 3.0]], [2.0, 4.0, 6.0], 0.0)
        self.assertTrue(all(weight >= 0.0 for weight in fitted["weights"]))
        self.assertAlmostEqual(sum(fitted["weights"]), 1.0, places=12)
        self.assertLessEqual(fitted["kkt"]["max_residual"], 1e-8)
        self.assertEqual(fitted["converged"], True)

    def test_lambda_holdout_and_tie_rule(self):
        result = _analyze("data/hand-example.json")
        selection = result["lambda_selection"]
        self.assertEqual(selection["selected"], 0.01)
        self.assertEqual([candidate["lambda"] for candidate in selection["candidates"]], [0.0, 0.01])
        for candidate in selection["candidates"]:
            self.assertAlmostEqual(candidate["holdout_mse"], 0.0, places=12)
        self.assertEqual(result["fit"]["lambda"], 0.01)
        self.assertEqual(result["protocol"]["training_weeks"], [0, 1, 2])
        self.assertEqual(result["protocol"]["holdout_weeks"], [3])

    def test_hand_example_weights(self):
        result = _analyze("data/hand-example.json")
        weights = {entry["donor"]: entry["weight"] for entry in result["fit"]["weights"]}
        self.assertAlmostEqual(weights["donor-1"], 0.5, places=12)
        self.assertAlmostEqual(weights["donor-2"], 0.5, places=12)
        self.assertEqual(result["fit"]["ratio_kind"], "zero_over_zero")
        self.assertEqual(result["fit"]["same_period_ratio"], 0.0)

    def test_post_treatment_values_do_not_change_weights(self):
        normal = _analyze("data/normal.json")
        changed = _analyze("data/changed-post-treatment.json")
        normal_weights = [entry["weight"] for entry in normal["fit"]["weights"]]
        changed_weights = [entry["weight"] for entry in changed["fit"]["weights"]]
        self.assertEqual(normal_weights, changed_weights)
        self.assertEqual(normal["lambda_selection"]["selected"], changed["lambda_selection"]["selected"])
        self.assertEqual(
            [candidate["holdout_mse"] for candidate in normal["lambda_selection"]["candidates"]],
            [candidate["holdout_mse"] for candidate in changed["lambda_selection"]["candidates"]],
        )
        self.assertEqual(
            [entry for entry in normal["weekly"] if entry["period"] == "pre"],
            [entry for entry in changed["weekly"] if entry["period"] == "pre"],
        )
        self.assertNotEqual(
            [entry["gap"] for entry in normal["weekly"] if entry["period"] == "post"],
            [entry["gap"] for entry in changed["weekly"] if entry["period"] == "post"],
        )

    def test_zero_pre_error_conventions(self):
        result = _analyze("data/zero-pre-error.json")
        self.assertEqual(result["fit"]["ratio_kind"], "zero_over_zero")
        self.assertEqual(result["fit"]["same_period_ratio"], 0.0)
        self.assertEqual(result["fit"]["post_msp"], 0.0)
        json.dumps(result, allow_nan=False)

        value = load("data/zero-pre-error.json")[0]
        for row in value["rows"]:
            if row["site"] == "treated" and row["week"] >= value["treatment_week"]:
                row["value"] = 2.0
        positive = analyze(value)
        self.assertEqual(positive["fit"]["ratio_kind"], "positive_over_zero")
        self.assertEqual(positive["fit"]["same_period_ratio"], "+infinity")
        json.dumps(positive, allow_nan=False)

    def test_placebo_protocol_and_rank(self):
        result = _analyze("data/normal.json")
        placebo = result["spatial_placebo"]
        self.assertTrue(placebo["not_a_randomization_p_value"])
        donors = ["donor-a", "donor-b", "donor-c"]
        self.assertEqual([entry["donor"] for entry in placebo["placebos"]], donors)
        for entry in placebo["placebos"]:
            self.assertNotIn("treated", entry["donor_pool"])
            self.assertNotIn(entry["donor"], entry["donor_pool"])
        eligible = [entry for entry in placebo["placebos"] if entry["eligible"]]
        self.assertEqual(placebo["included_denominator"], 1 + len(eligible))
        # 用公开输出重新核对秩：+infinity 需要按约定的正/零处理。
        treated = result["fit"]["same_period_ratio"]
        count = 0
        for entry in eligible:
            ratio = entry["same_period_ratio"]
            if isinstance(ratio, str) and ratio == "+infinity":
                count += 1 if treated == "+infinity" else 0
            elif isinstance(treated, str):
                count += 0
            elif ratio >= treated:
                count += 1
        self.assertEqual(placebo["included_numerator"], 1 + count)
        self.assertAlmostEqual(
            placebo["tail_rank_proportion"],
            (1 + count) / (1 + len(eligible)),
            places=12,
        )
        self.assertTrue(placebo["ties_included_in_numerator"])

    def test_placebo_exclusion_reasons(self):
        result = _analyze("data/normal.json")
        excluded = result["spatial_placebo"]["excluded"]
        self.assertTrue(excluded)
        for entry in excluded:
            self.assertEqual(entry["reason"], "pre_rmspe_exceeds_threshold")
            self.assertGreater(entry["pre_rmspe"], entry["threshold"])

    def test_treated_mismatch_still_reports_results(self):
        value = load("data/normal.json")[0]
        value["pre_rmspe_max"] = 0.0
        result = analyze(value)
        self.assertIsNotNone(result["diagnostics"]["mismatch_warning"])
        self.assertFalse(result["diagnostics"]["pre_fit_within_threshold"])
        self.assertEqual(result["spatial_placebo"]["rank"]["reportable"], False)
        self.assertEqual(result["spatial_placebo"]["rank"]["reason"], "treated_pre_rmspe_exceeds_threshold")
        self.assertTrue(result["weekly"])
        self.assertIn("post_average_gap", result["fit"])

    def test_budget_exhaustion_is_not_optimal(self):
        donors = [[10.0, 10.0], [0.0, 10.0]]
        observed = [20.0, 10.0]
        exhausted = fit_weights(donors, observed, 0.0, max_iterations=1)
        self.assertEqual(exhausted["termination_reason"], "budget_exhausted")
        self.assertFalse(exhausted["optimal"])
        self.assertFalse(exhausted["converged"])
        solved = fit_weights(donors, observed, 0.0, max_iterations=20)
        self.assertEqual(solved["termination_reason"], "kkt_converged")
        self.assertTrue(solved["optimal"])

    def test_balanced_panel_sample(self):
        result = _analyze("data/balanced-panel.json")
        self.assertEqual(result["protocol"]["treated"], "treated")
        self.assertEqual(len(result["weekly"]), len(result["protocol"]["weeks"]))
        self.assertTrue(all(set(entry) >= {"week", "observed", "synthetic", "gap", "period"} for entry in result["weekly"]))
        self.assertEqual(result["fit"]["solver"]["optimal"], True)
        self.assertLessEqual(result["fit"]["kkt"]["max_residual"], 1e-8)
        json.dumps(result, allow_nan=False)

    def test_cli_analyze_json(self):
        completed = subprocess.run(
            [sys.executable, "-m", "analysis_data.cli", "analyze", "data/hand-example.json"],
            check=True, capture_output=True, text=True,
        )
        parsed = json.loads(completed.stdout)
        self.assertEqual(parsed["fit"]["lambda"], 0.01)
        self.assertEqual(parsed["spatial_placebo"]["included_denominator"], 3)


if __name__ == "__main__":
    unittest.main()
