# Logistic Regression (Binary Outcomes)

Logistic regression is used when the outcome is **binary**:
- Disease yes/no
- Treatment response yes/no
- Hospitalization yes/no
- Adverse event yes/no

---

## 1. Why not linear regression?

A linear model can predict probabilities < 0 or > 1. Logistic regression ensures predictions stay in **[0,1]** and models **odds**.

---

## 2. Model: probability, odds, log-odds

Let $Y_i \in \{0,1\}$ and $p_i = P(Y_i=1 \mid X_i)$.

**Odds:**
$$
\text{odds}_i = \frac{p_i}{1-p_i}
$$

**Logit (log-odds):**
$$
\log\left(\frac{p_i}{1-p_i}\right) = \beta_0 + \beta_1 X_{1i} + \cdots + \beta_p X_{pi}
$$

Convert log-odds back to probability:
$$
p_i = \frac{1}{1 + e^{-(\beta_0 + \beta_1 X_{1i} + \cdots + \beta_p X_{pi})}}
$$

---

## 3. Interpretation: Odds Ratio (OR)

For a 1-unit increase in $X_j$:
$$
\text{OR}_j = e^{\beta_j}
$$

- OR > 1: higher odds of outcome
- OR < 1: lower odds
- OR = 1: no association

**Biostat example:** OR = 1.30 for smoking means **30% higher odds** of disease for smokers vs non-smokers, holding others constant.

---

## 4. Interactive Simulation: create disease data with known ORs

!!! interactive "Python"
    ```python
    import numpy as np
    import pandas as pd
    import statsmodels.api as sm

    np.random.seed(42)
    n = 1500

    # Predictors
    age = np.random.uniform(20, 80, n)
    smoker = np.random.binomial(1, 0.25, n)        # 25% smokers
    bmi = np.random.normal(27, 4.5, n)

    # ---- TRUE LOGISTIC MODEL (edit these!) ----
    beta0 = -9.0
    beta_age = 0.07          # per year
    beta_smoker = 0.9        # log-odds; OR = exp(0.9) ~ 2.46
    beta_bmi = 0.08          # per 1 BMI
    # ------------------------------------------

    linpred = beta0 + beta_age*age + beta_smoker*smoker + beta_bmi*bmi
    p = 1/(1 + np.exp(-linpred))
    disease = np.random.binomial(1, p, n)

    df = pd.DataFrame({"disease": disease, "age": age, "smoker": smoker, "bmi": bmi})

    X = sm.add_constant(df[["age", "smoker", "bmi"]])
    model = sm.Logit(df["disease"], X).fit(disp=False)

    # ORs with CI
    params = model.params
    conf = model.conf_int()
    or_table = pd.DataFrame({
        "beta": params,
        "OR": np.exp(params),
        "CI_low": np.exp(conf[0]),
        "CI_high": np.exp(conf[1]),
    })

    print(or_table)
    print("\nDisease prevalence:", round(df['disease'].mean(), 4))
    ```

Try:
- Increase `beta_smoker` and see OR increase and prevalence rise
- Make `beta_bmi` negative (protective) and watch OR < 1
- Make `beta0` less negative to increase baseline risk

---

## 5. Predicted probabilities (clinical meaning)

Odds ratios can be hard to interpret clinically. Probabilities are often easier.

!!! interactive "Python"
    ```python
    import numpy as np

    # Pull fitted coefficients
    b = model.params

    def predict_prob(age, smoker, bmi):
        lp = b["const"] + b["age"]*age + b["smoker"]*smoker + b["bmi"]*bmi
        return 1/(1 + np.exp(-lp))

    # Example patients
    p1 = predict_prob(age=40, smoker=0, bmi=24)
    p2 = predict_prob(age=40, smoker=1, bmi=24)
    p3 = predict_prob(age=65, smoker=1, bmi=31)

    print("Age 40, non-smoker, BMI 24:", round(p1, 4))
    print("Age 40, smoker, BMI 24:    ", round(p2, 4))
    print("Age 65, smoker, BMI 31:    ", round(p3, 4))
    ```

---

## 6. Model performance: ROC-AUC (discrimination)

AUC answers: “How well can the model separate cases vs non-cases?”

!!! interactive "Python"
    ```python
    import numpy as np
    import matplotlib.pyplot as plt
    from sklearn.metrics import roc_curve, roc_auc_score

    pred = model.predict(X)

    fpr, tpr, thr = roc_curve(df["disease"], pred)
    auc = roc_auc_score(df["disease"], pred)

    plt.plot(fpr, tpr)
    plt.plot([0,1],[0,1], linestyle="--")
    plt.title(f"ROC Curve (AUC = {auc:.3f})")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.show()
    ```

---

## 7. Calibration (are predicted risks too high/low?)

Discrimination (AUC) is not calibration. Calibration checks whether predicted probabilities match observed rates.

!!! interactive "Python"
    ```python
    import pandas as pd
    import numpy as np

    pred = model.predict(X)
    df_cal = df.copy()
    df_cal["pred"] = pred
    df_cal["bin"] = pd.qcut(df_cal["pred"], 10, duplicates="drop")

    cal = df_cal.groupby("bin").agg(
        mean_pred=("pred", "mean"),
        obs_rate=("disease", "mean"),
        n=("disease", "size")
    ).reset_index(drop=True)

    print(cal)
    ```

If mean_pred is much higher than obs_rate in many bins, your model is overpredicting.

---

## 8. Categorical predictors + reference groups

If you have multi-level categories (e.g., race/ethnicity groups), logistic regression uses **dummy variables** with a **reference level**.

Interpretation of an OR for a category:
> odds relative to the reference group, adjusted for other covariates.

---

## 9. Common pitfalls in biostat logistic regression

- **Odds ratio ≠ risk ratio** (especially when outcome is common)
- **Separation** (perfect prediction → unstable estimates)
- **Too many predictors for too few events** (overfitting)
- **Nonlinearity** for continuous predictors (consider splines)
- **Interactions** (effect modification) matter clinically

---

## Exercises (very applied)

<details>
<summary>Click to try</summary>

1. Make the disease outcome common by increasing `beta0`. Compare OR interpretation vs probability interpretation.  
2. Add an interaction `smoker * age` in the true model and fit a model with/without interaction.  
3. Create a confounder: make smokers older, refit models with/without age adjustment, compare smoker OR.  
4. Change the sample size `n` and observe CI width and AUC stability.

</details>
