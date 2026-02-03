# Collinearity, Confounding, and Practical Modeling Workflow

In biostat regression, the goal is often **causal-ish interpretation** (even if imperfect).
So we must think about:
- **confounding** (bias)
- **collinearity** (unstable estimates)
- **model building strategy**

---

## 1. Confounding (bias in estimated association)

A confounder:
- is associated with exposure $X$
- independently affects outcome $Y$
- is not on the causal pathway

Example:
- Exposure: smoking
- Outcome: lung function
- Confounder: age

---

## 2. Collinearity (big SE, unstable coefficients)

When predictors are highly correlated:
- coefficients may flip sign
- SE inflate
- p-values become misleading

---

## 3. Interactive: build collinearity and watch SE explode

!!! interactive "Python"
    ```python
    import numpy as np
    import pandas as pd
    import statsmodels.api as sm

    np.random.seed(30)
    n = 400

    age = np.random.uniform(20, 80, n)

    # Two highly correlated predictors:
    bmi = 25 + 0.08*(age-50) + np.random.normal(0, 2, n)
    waist = 0.9*bmi + np.random.normal(0, 0.5, n)  # very correlated with BMI

    # Outcome depends on BMI, not waist (truth)
    y = 100 + 0.7*age + 1.2*bmi + np.random.normal(0, 10, n)

    df = pd.DataFrame({"y": y, "age": age, "bmi": bmi, "waist": waist})

    m1 = sm.OLS(df["y"], sm.add_constant(df[["age", "bmi"]])).fit()
    m2 = sm.OLS(df["y"], sm.add_constant(df[["age", "bmi", "waist"]])).fit()

    print("Model without waist:\n", m1.summary(), "\n")
    print("Model with waist (collinearity):\n", m2.summary())
    ```

Watch:
- BMI coefficient and SE change when waist enters
- waist may look “significant” or BMI may lose significance due to instability

---

## 4. VIF (Variance Inflation Factor)

Rule of thumb:
- VIF > 5 (or 10) indicates problematic collinearity

!!! interactive "Python"
    ```python
    from statsmodels.stats.outliers_influence import variance_inflation_factor

    X = sm.add_constant(df[["age", "bmi", "waist"]]).values
    vifs = [variance_inflation_factor(X, i) for i in range(X.shape[1])]
    print(["const", "age", "bmi", "waist"])
    print(vifs)
    ```

---

## 5. Practical workflow

1. Plot data (relationships + missingness patterns)
2. Decide estimand: prediction vs association vs causal interpretation
3. Pre-specify covariates using subject-matter knowledge
4. Fit baseline model
5. Check diagnostics and stability (VIF, influence)
6. Sensitivity analyses (remove influential points, alternative coding, nonlinear terms)
7. Report effect sizes with CI, not just p-values

---

## Exercises

<details>
<summary>Click to try</summary>

1. Reduce correlation between waist and BMI by changing `waist = 0.5*bmi + ...`.  
2. Add a strong confounder and show bias when it’s omitted.  
3. Use regularization (next file) to stabilize correlated predictors.

</details>
