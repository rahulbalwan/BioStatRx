# Reporting Regression Results (Biostat / Clinical Style)

A regression model is only useful if you report it clearly:
- effect sizes
- uncertainty (CI)
- clinical meaning (absolute risk, predicted values)
- limitations

---

## 1. What to report (minimum)

### Linear regression
- $\hat\beta$ with 95% CI (and units!)
- interpret per clinically meaningful increment (e.g., per 10 years, per 5 BMI)
- model diagnostics summary (residuals, influential points)

### Logistic regression
- OR with 95% CI
- baseline risk or predicted probabilities for typical patient profiles
- AUC + calibration (if prediction)
- events per variable / overfitting consideration

---

## 2. Interactive: Create a clean results table (Linear)

!!! interactive "Python"
    ```python
    import numpy as np
    import pandas as pd
    import statsmodels.api as sm

    np.random.seed(77)
    n = 350
    age = np.random.uniform(30, 80, n)
    bmi = np.random.normal(27, 4, n)
    smoker = np.random.binomial(1, 0.3, n)

    sbp = 95 + 0.8*age + 1.1*bmi + 5.0*smoker + np.random.normal(0, 12, n)

    df = pd.DataFrame({"sbp": sbp, "age": age, "bmi": bmi, "smoker": smoker})

    X = sm.add_constant(df[["age", "bmi", "smoker"]])
    fit = sm.OLS(df["sbp"], X).fit()

    # Build coefficient table with CI
    ci = fit.conf_int()
    out = pd.DataFrame({
        "beta": fit.params,
        "CI_low": ci[0],
        "CI_high": ci[1],
        "p_value": fit.pvalues
    })

    # Optional: interpret age per 10 years
    out.loc["age", ["beta", "CI_low", "CI_high"]] *= 10
    out.rename(index={"age": "age (per 10 years)"}, inplace=True)

    print(out)
    ```

---

## 3. Interactive: Clean OR table (Logistic)

!!! interactive "Python"
    ```python
    import numpy as np
    import pandas as pd
    import statsmodels.api as sm

    np.random.seed(78)
    n = 1200
    age = np.random.uniform(20, 85, n)
    smoker = np.random.binomial(1, 0.25, n)
    bmi = np.random.normal(27, 4.5, n)

    lp = -7 + 0.06*age + 0.8*smoker + 0.07*bmi
    p = 1/(1 + np.exp(-lp))
    disease = np.random.binomial(1, p, n)

    df = pd.DataFrame({"disease": disease, "age": age, "bmi": bmi, "smoker": smoker})

    X = sm.add_constant(df[["age", "bmi", "smoker"]])
    fit = sm.Logit(df["disease"], X).fit(disp=False)

    ci = fit.conf_int()
    tab = pd.DataFrame({
        "beta": fit.params,
        "OR": np.exp(fit.params),
        "CI_low": np.exp(ci[0]),
        "CI_high": np.exp(ci[1]),
        "p_value": fit.pvalues
    })

    # Interpret age per 10 years
    tab.loc["age", ["beta"]] *= 10
    tab.loc["age", ["OR", "CI_low", "CI_high"]] = np.exp(tab.loc["age", "beta"]), np.exp(ci.loc["age", 0]*10), np.exp(ci.loc["age", 1]*10)
    tab.rename(index={"age": "age (per 10 years)"}, inplace=True)

    print(tab)
    print("\nOutcome prevalence:", round(df["disease"].mean(), 3))
    ```

---

## 4. Writing the interpretation (templates)

### Linear regression interpretation template
> After adjusting for **[covariates]**, a **[unit]** increase in **X** was associated with a **β** change in **Y** (95% CI: [L, U]).

### Logistic regression interpretation template
> After adjusting for **[covariates]**, **X** was associated with **OR** times the odds of **Y=1** (95% CI: [L, U]).

Add clinical translation:
- predicted probability differences for typical patients
- absolute risk is often more intuitive than OR

---

## 5. Checklist (before publishing)

- Units correct and interpretable
- Reference groups clearly stated
- Linearity / nonlinearity checked for continuous predictors
- Influence diagnostics performed
- Sensitivity analysis documented (optional but strong)
- For prediction: AUC + calibration
- Avoid “significant/non-significant” framing as the only conclusion

---

## Exercises

<details>
<summary>Click to try</summary>

1. Rewrite interpretation for smoker coefficient in both linear and logistic settings.  
2. Report results using clinically meaningful increments (e.g., BMI per 5 units).  
3. Create a predicted probability for a 60-year-old smoker with BMI=30 and compare to non-smoker.

</details>
