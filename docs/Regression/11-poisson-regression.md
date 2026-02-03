# Poisson Regression (Counts & Rates) 

Poisson regression is used for **count outcomes**:
- number of infections
- number of ER visits
- number of adverse events
- number of hospitalizations

It’s especially common when modeling **rates** (counts per person-time), e.g. events per patient-year.

---

## 1. When Poisson regression is appropriate

Use Poisson regression when:
- $Y$ is a **count**: $0,1,2,\dots$
- events happen independently in a time/space interval (approx.)
- mean–variance relationship is roughly Poisson:
  $$
  \mathbb{E}[Y] \approx \mathrm{Var}(Y)
  $$

If variance is much larger than the mean → consider **Negative Binomial** (overdispersion).

---

## 2. Model: log link (and why)

Poisson regression models the **log of the expected count**:

$$
Y_i \sim \mathrm{Poisson}(\mu_i),\quad
\log(\mu_i) = \beta_0 + \beta_1 X_{1i} + \cdots + \beta_p X_{pi}
$$

Because of the log link:
- $\mu_i$ is always positive
- coefficients become **multiplicative** on the original count scale

---

## 3. Interpretation: Rate Ratio / Incidence Rate Ratio (IRR)

Exponentiating coefficients gives an **incidence rate ratio** (IRR):

$$
\mathrm{IRR}_j = e^{\beta_j}
$$

Interpretation:
- IRR = 1.20 → **20% higher** event rate per 1-unit increase in predictor
- IRR = 0.80 → **20% lower** event rate

---

## 4. Modeling rates with an offset (person-time)

In cohort studies you often observe different follow-up times.

Let:
- $Y_i$ = number of events for person $i$
- $t_i$ = person-time (e.g., years of follow-up)

Model **rate** by including an **offset**:

$$
\log(\mu_i) = \beta_0 + \beta_1 X_i + \cdots + \log(t_i)
$$

where $\log(t_i)$ is included as an offset (coefficient fixed at 1).

This makes the model effectively:
- predict *events per time*
- compare rates across groups

---

## 5. Interactive simulation: events with unequal follow-up (biostat realistic)

!!! interactive "Python"
    ```python
    import numpy as np
    import pandas as pd
    import statsmodels.api as sm

    np.random.seed(123)
    n = 1500

    # Predictors
    age = np.random.uniform(20, 80, n)
    smoker = np.random.binomial(1, 0.25, n)

    # Follow-up time (person-years): varies across individuals
    t = np.random.uniform(0.5, 5.0, n)

    # ---- TRUE RATE MODEL (edit these) ----
    beta0 = -3.2
    beta_age = 0.015     # per year
    beta_smoker = 0.55   # IRR ~ exp(0.55) ≈ 1.73
    # -------------------------------------

    # log(rate) = beta0 + beta_age*age + beta_smoker*smoker
    log_rate = beta0 + beta_age*age + beta_smoker*smoker
    rate = np.exp(log_rate)  # events per person-year

    # expected count = rate * time
    mu = rate * t

    # generate counts
    y = np.random.poisson(mu, n)

    df = pd.DataFrame({"events": y, "age": age, "smoker": smoker, "py": t})

    # Fit Poisson with offset log(person-time)
    X = sm.add_constant(df[["age", "smoker"]])
    pois = sm.GLM(df["events"], X, family=sm.families.Poisson(), offset=np.log(df["py"])).fit()

    print(pois.summary())

    # IRRs with CI
    params = pois.params
    ci = pois.conf_int()
    irr = pd.DataFrame({
        "beta": params,
        "IRR": np.exp(params),
        "CI_low": np.exp(ci[0]),
        "CI_high": np.exp(ci[1]),
    })
    print("\nIRR table:\n", irr)

    print("\nMean events:", round(df["events"].mean(), 3))
    print("Mean person-years:", round(df["py"].mean(), 3))
    ```

Try:
- Increase `beta_smoker` to see smoking IRR grow
- Set `beta_age=0` to remove age effect
- Make follow-up time shorter (e.g., `t ~ Uniform(0.1, 1.0)`) and see counts drop

---

## 6. Overdispersion check (common in real biostat data)

Poisson assumes variance ≈ mean. Real data often has **overdispersion**.

A quick check:
$$
\frac{\text{Deviance}}{\text{df}} \approx 1 \text{ (Poisson ok)}, \quad \gg 1 \text{ (overdispersion)}
$$

!!! interactive "Python"
    ```python
    dev = pois.deviance
    df_resid = pois.df_resid
    print("Deviance/df:", round(dev/df_resid, 3))
    ```

If Deviance/df is large (e.g., > 1.5–2), consider:
- robust SE
- Negative Binomial regression
- zero-inflated models (if many zeros)

---

## 7. Robust standard errors (easy fix for mild overdispersion)

!!! interactive "Python"
    ```python
    pois_rob = pois.get_robustcov_results(cov_type="HC3")
    print(pois_rob.summary())

    # Robust IRR table
    params = pois_rob.params
    ci = pois_rob.conf_int()
    irr_rob = pd.DataFrame({
        "beta": params,
        "IRR": np.exp(params),
        "CI_low": np.exp(ci[0]),
        "CI_high": np.exp(ci[1]),
    })
    print("\nRobust IRR table:\n", irr_rob)
    ```

---

## 8. Practical biostat examples

- **Incidence rate of infections** by treatment group (offset = person-time)
- **Hospital visits** by exposure status adjusting for age/sex
- **Adverse event counts** in clinical trials (often person-time differs due to dropout)

---

## 9. Exercises

<details>
<summary>Click to try</summary>

1. Remove the offset and compare coefficients (should be biased when follow-up varies).  
2. Create overdispersion by mixing two subpopulations (different baseline rates) and re-check Deviance/df.  
3. Add a categorical predictor (treatment groups) and interpret IRR vs reference.  
4. Simulate many zeros by making baseline rate very small; discuss if zero-inflated models might be needed.

</details>

---

## Summary

- Poisson regression models counts using a log link.
- Exponentiated coefficients are **IRRs** (rate ratios).
- Use an **offset** to model event rates per person-time.
- Always check for **overdispersion**; consider robust SE or Negative Binomial if needed.
