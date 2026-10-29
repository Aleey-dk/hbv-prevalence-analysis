"""
Prevalence of Hepatitis B Virus (HBsAg) among patients attending
Dutse General Hospital, Jigawa State, Nigeria (N = 100).

Companion analysis for the BSc Microbiology thesis by Aliyu Haruna Ubale
(Federal University Dutse, December 2024).

Functions
---------
- wilson_ci: 95% Wilson score confidence interval for a proportion
- stratified_prevalence: overall and per-group prevalence with CIs
- chi_square_test: chi-square test of independence for a 2xK table
"""
import csv
import math
from pathlib import Path

import numpy as np
from scipy import stats

DATA_PATH = Path(__file__).parent / "data" / "patients_summary.csv"
Z_95 = 1.959963984540054  # two-sided 95% normal quantile


def wilson_ci(successes: int, total: int, z: float = Z_95) -> tuple[float, float]:
    """Wilson score interval for a binomial proportion."""
    if total <= 0:
        raise ValueError("total must be positive")
    if not 0 <= successes <= total:
        raise ValueError("successes must be between 0 and total")
    p = successes / total
    denom = 1 + z**2 / total
    centre = (p + z**2 / (2 * total)) / denom
    half = z * math.sqrt(p * (1 - p) / total + z**2 / (4 * total**2)) / denom
    return (max(0.0, centre - half), min(1.0, centre + half))


def stratified_prevalence(rows: list[dict]) -> list[dict]:
    """Prevalence and 95% Wilson CI for each group/category row."""
    out = []
    for r in rows:
        n, k = int(r["examined"]), int(r["positive"])
        lo, hi = wilson_ci(k, n)
        out.append({
            "group": r["group"],
            "category": r["category"],
            "examined": n,
            "positive": k,
            "prevalence": k / n,
            "ci95_low": lo,
            "ci95_high": hi,
        })
    return out


def chi_square_test(rows: list[dict]) -> dict:
    """Chi-square test of independence on a 2xK positives/negatives table."""
    positives = np.array([int(r["positive"]) for r in rows])
    examined = np.array([int(r["examined"]) for r in rows])
    negatives = examined - positives
    table = np.vstack([positives, negatives])
    if (table.sum(axis=0) == 0).any():
        raise ValueError("each category must have at least one patient")
    chi2, p, dof, _ = stats.chi2_contingency(table)
    return {"chi2": chi2, "p_value": p, "dof": dof}


def load_data(path: Path = DATA_PATH) -> list[dict]:
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def main() -> None:
    rows = load_data()
    total_n = sum(int(r["examined"]) for r in rows if r["group"] == "sex")
    total_k = sum(int(r["positive"]) for r in rows if r["group"] == "sex")
    lo, hi = wilson_ci(total_k, total_n)
    print(f"Overall prevalence: {total_k}/{total_n} = {total_k/total_n:.1%} "
          f"(95% CI {lo:.1%}-{hi:.1%})\n")
    for group in ("sex", "age_group", "marital_status", "pregnancy_status"):
        sub = [r for r in rows if r["group"] == group]
        print(f"== {group} ==")
        for s in stratified_prevalence(sub):
            print(f"  {s['category']}: {s['positive']}/{s['examined']} = "
                  f"{s['prevalence']:.1%} (95% CI {s['ci95_low']:.1%}-{s['ci95_high']:.1%})")
        res = chi_square_test(sub)
        print(f"  chi-square = {res['chi2']:.3f}, dof = {res['dof']}, p = {res['p_value']:.4f}\n")


if __name__ == "__main__":
    main()
