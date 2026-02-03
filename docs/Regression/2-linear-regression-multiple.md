# Multiple Linear Regression (Adjusted Effects)

Multiple linear regression extends linear regression to **multiple predictors**.

Example:
- Outcome: SBP (mmHg)
- Predictors: age, BMI, sex, smoking
- Goal: estimate association of BMI with SBP **adjusting** for age, etc.

---

## 1. Model

$$
Y_i = \beta_0 + \beta_1 X_{1i} + \beta_2 X_{2i} + \cdots + \beta_p X_{pi} + \varepsilon_i
$$

Interpretation:
- $\beta_j$ is the expected change in $Y$ for a 1-unit increase in $X_j$, **holding the other predictors constant**.

---

## 2. Confounding (Why “adjustment” matters)

Suppose:
- older people have higher SBP
- older people also tend to have higher BMI

Then BMI–SBP association can look stronger than it truly is if you ignore age.

---

## 3. Interactive Simulation: confounding + adjustment

!!! interactive "Python"
    ```python
    import numpy as np
    import pandas as pd
    import statsmodels.api as sm

    np.random.seed(7)

    n = 400

    # Age is a confounder: affects both BMI and SBP
    age = np.random.uniform(20, 80, n)

    # BMI increases with age (confounding structure)
    bmi = 22 + 0.06*(age - 50) + np.random.normal(0, 2.5, n)

    # TRUE model: SBP depends on age strongly and BMI moderately
    sbp = 105 + 0.75*age + 0.9*bmi + np.random.normal(0, 10, n)

    df = pd.DataFrame({"sbp": sbp, "age": age, "bmi": bmi})

    # Naive model (unadjusted): SBP ~ BMI
    m_naive = sm.OLS(df["sbp"], sm.add_constant(df["bmi"])).fit()

    # Adjusted model: SBP ~ BMI + Age
    m_adj = sm.OLS(df["sbp"], sm.add_constant(df[["bmi", "age"]])).fit()

    print("Unadjusted BMI slope:", round(m_naive.params["bmi"], 3))
    print("Adjusted BMI slope:  ", round(m_adj.params["bmi"], 3))
    print("\nAdjusted model summary:\n")
    print(m_adj.summary())
    ```

Try:
- Change the BMI–age relationship (`0.06`) to make confounding stronger/weaker
- Change the true BMI effect (`0.9`) and see if adjustment recovers it

---

## 4. Categorical Predictors (Sex)

Binary coding:
- sex = 0 (female), 1 (male)

Then $\beta_{\text{sex}}$ is the **mean difference** in $Y$ between males vs females (adjusted).

!!! interactive "Python"
    ```python
    import numpy as np
    import pandas as pd
    import statsmodels.api as sm

    np.random.seed(10)
    n = 300
    age = np.random.uniform(20, 80, n)
    sex = np.random.binomial(1, 0.5, n)  # 0/1
    bmi = 23 + 0.05*(age-50) + 1.2*sex + np.random.normal(0, 2.2, n)

    sbp = 100 + 0.8*age + 0.7*bmi + 4.5*sex + np.random.normal(0, 10, n)

    df = pd.DataFrame({"sbp": sbp, "age": age, "bmi": bmi, "sex": sex})

    model = sm.OLS(df["sbp"], sm.add_constant(df[["age", "bmi", "sex"]])).fit()
    print(model.summary())
    ```

Interpret:
- sex coefficient: adjusted mean difference in SBP for sex=1 vs sex=0

---

## 5. Interaction (Effect Modification)

Sometimes the effect of BMI differs by sex:

$$
SBP = \beta_0 + \beta_1 BMI + \beta_2 Sex + \beta_3 (BMI \times Sex) + \cdots
$$

- If $\beta_3 \neq 0$, the BMI slope depends on sex.

!!! interactive "Python"
    ```python
    import numpy as np
    import pandas as pd
    import statsmodels.formula.api as smf

    np.random.seed(11)
    n = 500
    sex = np.random.binomial(1, 0.5, n)
    bmi = np.random.normal(27, 4, n)

    # True effect: BMI slope is higher in sex=1
    sbp = 110 + 0.8*bmi + 3.0*sex + 0.6*(bmi*sex) + np.random.normal(0, 10, n)

    df = pd.DataFrame({"sbp": sbp, "bmi": bmi, "sex": sex})

    m_no_int = smf.ols("sbp ~ bmi + sex", data=df).fit()
    m_int = smf.ols("sbp ~ bmi * sex", data=df).fit()

    print("No interaction model:\n", m_no_int.params, "\n")
    print("With interaction model:\n", m_int.params)
    ```

---

## Exercises

<details>
<summary>Click to try</summary>

1. Create a strong confounder and show how unadjusted slope is biased.  
2. Add collinearity: make BMI highly correlated with another predictor; watch standard errors inflate.  
3. Fit an interaction model and interpret slopes separately for sex=0 and sex=1.  

</details>
