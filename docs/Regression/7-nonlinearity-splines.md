# Nonlinearity & Splines (Biostat-friendly)

Many biological relationships are **not linear**:
- dose-response curves
- risk vs age
- biomarkers with thresholds

Linear regression can mislead if the true pattern is curved.

---

## 1. Detecting nonlinearity

Signs:
- residuals show curves vs fitted
- scatterplot suggests curvature
- strong model misspecification

---

## 2. Polynomial terms (simple but can be unstable)

Example:
$$
Y = \beta_0 + \beta_1 X + \beta_2 X^2 + \varepsilon
$$

---

## 3. Splines (preferred in biostat)

Splines allow smooth, flexible curves without forcing a global polynomial.

A common option: **natural cubic splines**.

---

## 4. Interactive: linear vs polynomial vs spline fit

!!! interactive "Python"
    ```python
    import numpy as np
    import matplotlib.pyplot as plt
    import pandas as pd
    import statsmodels.formula.api as smf
    from patsy import dmatrix

    np.random.seed(21)
    n = 250
    age = np.random.uniform(20, 80, n)

    # True nonlinear relationship (risk-like curve)
    y = 120 + 0.4*age + 0.03*(age-50)**2 + np.random.normal(0, 8, n)

    df = pd.DataFrame({"y": y, "age": age})

    # Linear model
    m_lin = smf.ols("y ~ age", data=df).fit()

    # Quadratic model
    m_quad = smf.ols("y ~ age + I(age**2)", data=df).fit()

    # Spline model (natural cubic spline with df=4)
    spline_basis = dmatrix("cr(age, df=4)", data=df, return_type="dataframe")
    df_s = pd.concat([df, spline_basis], axis=1)

    cols = spline_basis.columns.tolist()
    formula = "y ~ " + " + ".join(cols)
    m_spline = smf.ols(formula, data=df_s).fit()

    # Plot
    xg = np.linspace(age.min(), age.max(), 300)
    pred_lin = m_lin.predict(pd.DataFrame({"age": xg}))
    pred_quad = m_quad.predict(pd.DataFrame({"age": xg}))

    xg_df = pd.DataFrame({"age": xg})
    xg_spline = dmatrix("cr(age, df=4)", data=xg_df, return_type="dataframe")
    xg_s = pd.concat([xg_df, xg_spline], axis=1)
    pred_spline = m_spline.predict(xg_s)

    plt.scatter(age, y)
    plt.plot(xg, pred_lin, label="Linear")
    plt.plot(xg, pred_quad, label="Quadratic")
    plt.plot(xg, pred_spline, label="Spline (df=4)")
    plt.title("Comparing Linear vs Polynomial vs Spline")
    plt.xlabel("Age")
    plt.ylabel("Outcome")
    plt.legend()
    plt.show()

    print("R^2 linear:", round(m_lin.rsquared, 3))
    print("R^2 quad:  ", round(m_quad.rsquared, 3))
    print("R^2 spline:", round(m_spline.rsquared, 3))
    ```

Try:
- Change the true curve strength: the `0.03*(age-50)**2` term
- Change spline flexibility: `df=3` vs `df=6`

---

## 5. Logistic regression + splines (common in risk modeling)

Risk vs age is rarely linear in log-odds. Splines are standard in clinical prediction modeling.

(We can add a full logistic spline demo later if you want.)

---

## Exercises

<details>
<summary>Click to try</summary>

1. Fit linear model only and inspect residual patterns.  
2. Increase spline
