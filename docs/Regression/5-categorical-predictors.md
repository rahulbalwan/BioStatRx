# Categorical Predictors (Dummy Coding)

Many biostat predictors are categorical:
- sex (female/male)
- treatment group (placebo/drug)
- smoking status (never/former/current)
- genotype category (AA/AG/GG)

---

## 1. Dummy coding (one-hot with reference group)

For a 3-level variable (never/former/current):
- choose a reference (e.g., never)
- create indicators:
  - former = 1 if former else 0
  - current = 1 if current else 0

Model:
$$
Y = \beta_0 + \beta_1 \cdot I(\text{former}) + \beta_2 \cdot I(\text{current}) + \cdots
$$

Interpretation:
- $\beta_1$ = difference (former − never), adjusted
- $\beta_2$ = difference (current − never), adjusted

---

## 2. Interactive simulation: treatment groups

!!! interactive "Python"
    ```python
    import numpy as np
    import pandas as pd
    import statsmodels.formula.api as smf

    np.random.seed(5)
    n = 300

    # Categorical treatment variable
    trt = np.random.choice(["placebo", "low_dose", "high_dose"], size=n, p=[0.4, 0.3, 0.3])
    age = np.random.uniform(30, 75, n)

    # True effects on SBP change (negative is improvement)
    # placebo: 0
    # low_dose: -4
    # high_dose: -8
    effect = np.where(trt == "placebo", 0, np.where(trt == "low_dose", -4, -8))

    y = 2 + 0.05*(age-50) + effect + np.random.normal(0, 4, n)

    df = pd.DataFrame({"sbp_change": y, "treatment": trt, "age": age})

    # By default, statsmodels chooses an alphabetical reference.
    # Force placebo as reference:
    model = smf.ols("sbp_change ~ C(treatment, Treatment(reference='placebo')) + age", data=df).fit()
    print(model.summary())
    ```

Interpret:
- treatment[T.low_dose] = low_dose − placebo (adjusted)
- treatment[T.high_dose] = high_dose − placebo (adjusted)

---

## 3. Logistic regression with categories (OR interpretation)

Same idea, but exponentiate coefficients to get OR vs reference group.

!!! interactive "Python"
    ```python
    import statsmodels.formula.api as smf
    import numpy as np

    np.random.seed(6)
    n = 1000

    trt = np.random.choice(["placebo", "drug"], size=n, p=[0.5, 0.5])
    age = np.random.uniform(20, 80, n)

    # True log-odds: drug reduces risk
    beta0 = -3.0
    beta_age = 0.04
    beta_drug = -0.6   # OR ~ 0.55

    lp = beta0 + beta_age*age + beta_drug*(trt=="drug")
    p = 1/(1 + np.exp(-lp))
    event = np.random.binomial(1, p, n)

    df = pd.DataFrame({"event": event, "treatment": trt, "age": age})

    m = smf.logit("event ~ C(treatment) + age", data=df).fit(disp=False)
    print(m.summary())

    # ORs
    params = m.params
    print("\nOdds Ratios:")
    print(np.exp(params))
    ```

---

## 4. Practical tips

- Always state the **reference group** in reporting
- Make reference clinically meaningful (e.g., placebo, never-smoker)
- For ordered categories (mild/moderate/severe), consider trend tests or treat as ordinal

---

## Exercises

<details>
<summary>Click to try</summary>

1. Change reference group to "high_dose" and re-interpret coefficients.  
2. Add a 4-level smoking status and fit with adjustment for age and sex.  
3. Fit the same model but treat categories as numeric codes and explain why it’s risky.

</details>
