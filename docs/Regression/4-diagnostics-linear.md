# Diagnostics & Assumptions (Linear Regression)

Linear regression is powerful, but **only interpretable** when model assumptions are reasonably met.

Biostat examples:
- SBP ~ age + BMI
- LDL ~ statin dose + adherence

---

## 1. Key Assumptions (What to check)

For inference (p-values, CI) to be reliable, we usually want:

1. **Linearity**: mean of $Y$ changes linearly with predictors  
2. **Independent errors**: observations independent (common in cohort, not in repeated measures)  
3. **Constant variance (homoscedasticity)**: residual spread is roughly constant  
4. **Normality of residuals**: mainly for small samples / inference  
5. **No extreme influence**: no single point drives the whole model  

---

## 2. Interactive: Create data that violates assumptions

!!! interactive "Python"
    ```python
    import numpy as np
    import matplotlib.pyplot as plt
    import statsmodels.api as sm

    np.random.seed(1)

    n = 200
    x = np.random.uniform(0, 10, n)

    # ---- choose ONE scenario by uncommenting ----
    # (A) Well-behaved linear model
    y = 5 + 2*x + np.random.normal(0, 2, n)

    # (B) Nonlinearity
    # y = 5 + 2*x + 0.6*(x-5)**2 + np.random.normal(0, 2, n)

    # (C) Heteroscedasticity (variance increases with x)
    # y = 5 + 2*x + np.random.normal(0, 0.3*x, n)

    # (D) Outlier / influential point
    # y = 5 + 2*x + np.random.normal(0, 2, n)
    # y[0] += 40
    # x[0] += 8
    # --------------------------------------------

    X = sm.add_constant(x)
    fit = sm.OLS(y, X).fit()

    # Scatter + fitted line
    plt.scatter(x, y)
    xg = np.linspace(x.min(), x.max(), 200)
    plt.plot(xg, fit.params[0] + fit.params[1]*xg)
    plt.title("Data + fitted linear regression")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.show()

    print(fit.summary())
    ```

---

## 3. Residual plots (the #1 diagnostic)

Residuals vs fitted:
- random cloud around 0 ✅
- curve pattern → nonlinearity ❌
- funnel shape → heteroscedasticity ❌

!!! interactive "Python"
    ```python
    yhat = fit.predict(X)
    resid = y - yhat

    plt.scatter(yhat, resid)
    plt.axhline(0)
    plt.title("Residuals vs Fitted")
    plt.xlabel("Fitted values")
    plt.ylabel("Residuals")
    plt.show()
    ```

---

## 4. Normal Q-Q plot (residual normality)

Normality matters most for:
- small n
- hypothesis tests, CI

!!! interactive "Python"
    ```python
    import scipy.stats as st

    st.probplot(resid, dist="norm", plot=plt)
    plt.title("Q-Q Plot of Residuals")
    plt.show()
    ```

---

## 5. Influence: leverage & Cook’s distance

In biostat datasets, a few unusual patients can dominate.
We check **Cook’s distance**.

!!! interactive "Python"
    ```python
    infl = fit.get_influence()
    cooks = infl.cooks_distance[0]

    plt.stem(cooks, use_line_collection=True)
    plt.title("Cook's Distance by observation")
    plt.xlabel("Observation index")
    plt.ylabel("Cook's D")
    plt.show()

    # Show top 5 influential points
    top = np.argsort(cooks)[-5:][::-1]
    print("Top influential indices:", top)
    print("Top Cook's D:", cooks[top])
    ```

---

## 6. What to do when assumptions fail (biostat moves)

- **Nonlinearity** → add transforms or splines (see next files)
- **Heteroscedasticity** → robust SE (HC3), transform outcome, or model variance
- **Outliers/influence** → verify data, robust regression, sensitivity analysis
- **Non-independence** (repeated measures) → mixed models / GEE (later topic)

!!! interactive "Python"
    ```python
    # Robust standard errors (HC3) for heteroscedasticity
    robust = fit.get_robustcov_results(cov_type="HC3")
    print(robust.summary())
    ```

---

## Exercises

<details>
<summary>Click to try</summary>

1. Turn on scenario (B). Explain what residual plot shows.  
2. Turn on scenario (C). Compare standard vs robust SE.  
3. Turn on scenario (D). Remove the influential point and compare coefficients.  
4. Use a log-transform of y (if y>0) and see if residuals improve.

</details>
