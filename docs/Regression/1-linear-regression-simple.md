# Simple Linear Regression (Biostat)

Simple linear regression models a **continuous outcome** using **one predictor**.

Examples in biostat:
- Systolic BP explained by age
- Cholesterol explained by BMI
- Lung function explained by smoking pack-years

---

## 1. The Model

Outcome $Y$ (continuous) and predictor $X$:

$$
Y_i = \beta_0 + \beta_1 X_i + \varepsilon_i
$$

- $\beta_0$: expected $Y$ when $X=0$ (may be non-meaningful depending on scale)
- $\beta_1$: expected change in $Y$ for a **1-unit increase** in $X$
- $\varepsilon_i$: error term

---

## 2. Biostat Interpretation

If $Y$ = systolic BP (mmHg) and $X$ = age (years), and you estimate $\hat\beta_1 = 0.8$:

> Each additional year of age is associated with **+0.8 mmHg** higher systolic BP (on average).

---

## 3. Simulate Your Own Data (Interactive: change effect size + noise)

!!! interactive "Python"
    ```python
    import numpy as np
    import matplotlib.pyplot as plt

    np.random.seed(3)

    # ---- try editing these ----
    n = 120
    beta0 = 110
    beta1 = 0.7     # effect per year
    sigma = 10      # noise (sd)
    # --------------------------

    age = np.random.uniform(20, 80, size=n)
    eps = np.random.normal(0, sigma, size=n)
    sbp = beta0 + beta1 * age + eps

    # Fit with numpy
    X = np.column_stack([np.ones(n), age])
    beta_hat = np.linalg.lstsq(X, sbp, rcond=None)[0]
    b0, b1 = beta_hat

    # Plot
    plt.scatter(age, sbp)
    xgrid = np.linspace(age.min(), age.max(), 200)
    plt.plot(xgrid, b0 + b1 * xgrid)
    plt.title("Simple Linear Regression: SBP vs Age")
    plt.xlabel("Age (years)")
    plt.ylabel("SBP (mmHg)")
    plt.show()

    print("Estimated intercept (b0):", round(b0, 3))
    print("Estimated slope (b1):", round(b1, 3))
    ```

Try:
- Increase `sigma` → more scatter, weaker signal
- Change `beta1` → stronger/weaker association
- Reduce `n` → more uncertainty

---

## 4. Estimation + Residuals

Residual:
$$
e_i = y_i - \hat{y}_i
$$

If assumptions are reasonable, residuals should look:
- centered around 0
- no trend with $X$
- roughly constant spread

!!! interactive "Python"
    ```python
    yhat = b0 + b1 * age
    resid = sbp - yhat

    plt.scatter(yhat, resid)
    plt.axhline(0)
    plt.title("Residuals vs Fitted")
    plt.xlabel("Fitted values")
    plt.ylabel("Residuals")
    plt.show()

    plt.hist(resid, bins=20, edgecolor="black")
    plt.title("Histogram of Residuals")
    plt.xlabel("Residual")
    plt.ylabel("Count")
    plt.show()
    ```

---

## 5. Inference: Confidence Interval for Slope (concept)

A 95% CI for $\beta_1$ is:
$$
\hat\beta_1 \pm t_{0.975,\,n-2}\cdot SE(\hat\beta_1)
$$

In practice you use a stats library (below).

!!! interactive "Python"
    ```python
    import statsmodels.api as sm

    model = sm.OLS(sbp, sm.add_constant(age)).fit()
    print(model.summary())
    ```

---

## 6. Exercises (Biostat flavored)

<details>
<summary>Click to try</summary>

1. Change `beta1` to 0 and confirm the fitted slope is near 0.  
2. Increase `sigma` to 25 and observe what happens to significance / CI width.  
3. Create a binary exposure (e.g., smoker/non-smoker) and fit a line using $X \in \{0,1\}$. Interpret $\beta_1$.  
4. Add a non-linear pattern: try `sbp = beta0 + beta1*age + 0.02*(age-50)**2 + eps` and see residual patterns.  

</details>
