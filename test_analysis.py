import math

import pytest

from analysis import chi_square_test, load_data, stratified_prevalence, wilson_ci


def test_overall_prevalence_matches_thesis():
    rows = [r for r in load_data() if r["group"] == "sex"]
    total_n = sum(int(r["examined"]) for r in rows)
    total_k = sum(int(r["positive"]) for r in rows)
    assert (total_n, total_k) == (100, 5)
    assert math.isclose(total_k / total_n, 0.05)


def test_sex_stratification():
    rows = [r for r in load_data() if r["group"] == "sex"]
    res = {s["category"]: s for s in stratified_prevalence(rows)}
    assert res["Male"]["positive"] == 2 and res["Male"]["examined"] == 55
    assert res["Female"]["positive"] == 3 and res["Female"]["examined"] == 45


def test_sex_distribution_matches_thesis_percentages():
    # thesis Table 4.1: positives are 40% male, 60% female
    rows = [r for r in load_data() if r["group"] == "sex"]
    res = {s["category"]: s for s in stratified_prevalence(rows)}
    assert res["Male"]["positive"] / 5 == pytest.approx(0.40)
    assert res["Female"]["positive"] / 5 == pytest.approx(0.60)


def test_wilson_ci_known_values():
    # Wilson CI for 5/100 is roughly 2.2% to 11.2%
    lo, hi = wilson_ci(5, 100)
    assert lo == pytest.approx(0.022, abs=0.003)
    assert hi == pytest.approx(0.112, abs=0.004)


def test_wilson_ci_zero_and_all():
    assert wilson_ci(0, 10) == (0.0, 0.0) or wilson_ci(0, 10)[1] > 0
    lo, hi = wilson_ci(10, 10)
    assert hi == pytest.approx(1.0)
    assert lo < 1.0


def test_wilson_ci_rejects_bad_input():
    with pytest.raises(ValueError):
        wilson_ci(6, 5)
    with pytest.raises(ValueError):
        wilson_ci(1, 0)


def test_pregnancy_uses_females_only_denominator():
    rows = [r for r in load_data() if r["group"] == "pregnancy_status"]
    assert sum(int(r["examined"]) for r in rows) == 45  # 45 female patients


def test_chi_square_runs_per_group():
    data = load_data()
    for group in ("sex", "age_group", "marital_status", "pregnancy_status"):
        sub = [r for r in data if r["group"] == group]
        res = chi_square_test(sub)
        assert 0 <= res["p_value"] <= 1
