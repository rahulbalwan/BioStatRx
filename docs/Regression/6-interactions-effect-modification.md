# Interactions (Effect Modification) in Biostatistics

In biostat, we often ask:
- Does treatment work differently by sex?
- Does exposure effect differ by age group?
- Do comorbidities modify risk?

This is **effect modification** → modeled with interactions.

---

## 1. Concept

Model with interaction:
$$
Y = \beta_0 + \beta_1 X + \beta_2 Z + \beta_3 (X \times Z) + \varepsilon
$$

Interpretation:
- When $Z=0$, slope of $X$ is $\beta_1$
- When $Z=1$, slope of $X$ is $\beta_1 + \beta_3$

---

## 2. Interactive: treatment effect differs by sex (linear regression)

!!! interactive "Python"
    ```python
    import numpy as np
    import pandas as pd
    import statsmodels.formula.api as smf

    np.random.seed(12)
    n = 500

    sex = np.random.binomial(1, 0.5, n)  # 0=female, 1=male
    trt = np.random.binomial(1, 0.5, n)  # 0=control, 1=drug

    # True mean outcome (e.g., SBP reduction; more negative = better)
    # Drug works better in males:
    base = -2
    trt_effect_female = -3
    trt_effect_male = -7

    effect = np.where(sex == 0, trt_effect_female, trt_effect_male)
    y = base + trt*effect + np.random.normal(0, 4, n)

    df = pd.DataFrame({"y": y, "trt": trt, "sex": sex})

    m0 = smf.ols("y ~ trt + sex", data=df).fit()
    m1 = smf.ols("y ~ trt * sex", data=df).fit()  # includes interaction

    print("No interaction:\n", m0.params, "\n")
    print("With interaction:\n", m1.params)
    print("\nFull summary:\n", m1.summary())
    ```

Interpret:
- `trt` = treatment effect when sex=0
- `trt:sex` = difference in treatment effect between sex=1 vs sex=0

---

## 3. Visualizing interaction

!!! interactive "Python"
    ```python
    import numpy as np

    # Predicted means by group:
    b = m1.params
    # y = b0 + b1*trt + b2*sex + b3*(trt*sex)
    def pred(trt, sex):
        return b["Intercept"] + b["trt"]*trt + b["sex"]*sex + b["trt:sex"]*trt*sex

    print("Female control:", round(pred(0,0), 3))
    print("Female drug:   ", round(pred(1,0), 3))
    print("Male control:  ", round(pred(0,1), 3))
    print("Male drug:     ", round(pred(1,1), 3))
    ```

---

## 4. Logistic regression interactions (OR depends on Z)

Logistic model:
$$
\log\left(\frac{p}{1-p}\right)=\beta_0+\beta_1 X+\beta_2 Z+\beta_3 XZ
$$

The OR for $X$ is:
- when $Z=0$: $\exp(\beta_1)$
- when $Z=1$: $\exp(\beta_1 + \beta_3)$

---

## Exercises

<details>
<summary>Click to try</summary>

1. Reduce interaction strength by changing trt_effect_male closer to trt_effect_female.  
2. Add age as an additional confounder and see if interaction remains.  
3. Fit a logistic interaction: make event probabilities differ by treatment-by-sex.

</details>
