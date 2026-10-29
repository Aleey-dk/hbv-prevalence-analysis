# HBV Seroprevalence Analysis — Dutse General Hospital (N = 100)

Companion computational analysis for the undergraduate thesis
**"Prevalence of Hepatitis B Virus Among Patients Attending Dutse General
Hospital, Jigawa State"** by Aliyu Haruna Ubale (FSC/MCB/19/1020),
Department of Microbiology and Biotechnology, Federal University Dutse,
December 2024. Supervisor: Dr. Mas'ud Abdullahi.

## What this repo does

- Computes overall HBsAg seroprevalence with 95% Wilson score confidence intervals
- Computes prevalence stratified by sex, age group, marital status and pregnancy status
- Runs chi-square tests of independence for each stratification
- Unit tests verify every headline figure against the thesis tables

## Data

`data/patients_summary.csv` reconstructs the grouped counts from the thesis
results tables (Tables 4.1, 4.3, 4.4 and 4.5):

| Stratification | Groups (examined / positive) |
|---|---|
| Sex | Male 55/2, Female 45/3 |
| Age | 17-25: 26/1, 25-30: 39/2, 30-45: 35/2 |
| Marital status | Married 30/1, Single 60/4, Widowed 10/0 |
| Pregnancy (females only) | Pregnant 15/2, Not pregnant 30/1 |

Overall prevalence: **5/100 = 5.0%** (95% Wilson CI ≈ 2.2%–11.2%).

## Running

```
pip install numpy scipy pytest
python analysis.py
pytest
```

Requires Python 3.10+.
