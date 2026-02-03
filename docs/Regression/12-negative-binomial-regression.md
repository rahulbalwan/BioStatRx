# Negative Binomial Regression (Overdispersed Counts) 

Negative Binomial (NB) regression is used for **count outcomes** when Poisson regression fails due to **overdispersion**:

- Poisson assumption: $\mathrm{Var}(Y) \approx \mathbb{E}[Y]$
- Overdispersion: $\mathrm{Var}(Y) \gg \mathbb{E}[Y]$ (very common in health data)

Examples:
- number of ER visits (many zeros + a few heavy users)
- infections per patient-year
- asthma exacerbations per year
- adverse events per person-time

---

## 1. Why overdispersion happens (biostat reality)

Common reasons:
- unmeasured heterogeneity (frailty): some patients are “high risk”
- clustering (hospital/site/provider effects)
- event dependence (prior event increases future risk)
- excess zeros (sometimes)

Poisson under overdispersion typically gives:
- too-small standard errors
- overly optimistic p-values / narrow CI

---

## 2. Model form (same link, extra dispersion)

NB regression still uses the log link:

$$
\log(\mu_i) = \beta_0 + \beta_1 X_{1i} + \cdots + \log(t_i)
$$

but variance is larger than mean. A common parameterization:

$$
\mathrm{Var}(Y_i) = \mu_i + \alpha \mu_i^2
$$

- $\alpha > 0$ captures extra-Poisson variability
- if $\alpha \to 0$, NB behaves like Poisson

---

## 3. Interpretation: IRR stays the same idea

Exponentiated coefficients are still **incidence rate ratios (IRR)**:

$$
\mathrm{IRR}_j = e^{\beta_j}
$$

Interpretation identical to Poisson:
- IRR = 1.30 → 30% higher event rate per 1-unit increase in $X_j$

---

## 4. Interactive simulation: generate overdispersed counts

We create overdispersion by adding person-level unobserved heterogeneity (random multiplier).

!!! interactive "Python"
    ```python
    import numpy as np
    import pandas as pd
    import statsmodels.api as sm

    np.random.seed(100)
    n = 2000

    age = np.random.uniform(20, 80, n)
    smoker = np.random.binomial(1, 0.25, n)
    py = np.random.uniform(0.5, 5.0, n)  # person-years

    # True log-rate model
    beta0 = -3.4
    beta_age = 0.015
    beta_smoker = 0.55  # IRR ~ 1.73

    log_rate = beta0 + beta_age*age + beta_smoker*smoker
    base_rate = np.exp(log_rate)  # events per person-year

    # Unobserved heterogeneity (frailty-like): multiplicative random effect
    # Larger variance here -> more overdispersion
    hetero_sd = 0.9
    u = np.exp(np.random.normal(0, hetero_sd, n))

    mu = base_rate * u * py
    y = np.random.poisson(mu)

    df = pd.DataFrame({"events": y, "age": age, "smoker": smoker, "py": py})

    # Fit Poisson with offset
    X = sm.add_constant(df[["age", "smoker"]])
    pois = sm.GLM(df["events"], X, family=sm.families.Poisson(), offset=np.log(df["py"])).fit()

    print("Poisson Deviance/df:", round(pois.deviance/pois.df_resid, 3))
    print(pois.summary())
    ```

Try:
- Increase `hetero_sd` to 1.2 → more overdispersion
- Decrease to 0.2 → Poisson becomes reasonable again

---

## 5. Fit Negative Binomial model (and compare)

!!! interactive "Python"
    ```python
    # Negative Binomial GLM with offset
    nb = sm.GLM(df["events"], X, family=sm.families.NegativeBinomial(), offset=np.log(df["py"])).fit()

    print("NB Deviance/df:", round(nb.deviance/nb.df_resid, 3))
    print(nb.summary())

    # IRR tables
    import numpy as np
    import pandas as pd

    def irr_table(fit):
        ci = fit.conf_int()
        tab = pd.DataFrame({
            "beta": fit.params,
            "IRR": np.exp(fit.params),
            "CI_low": np.exp(ci[0]),
            "CI_high": np.exp(ci[1]),
        })
        return tab

    print("\nPoisson IRR:\n", irr_table(pois))
    print("\nNegative Binomial IRR:\n", irr_table(nb))
    ```

What you’ll usually see:
- Similar point estimates (sometimes)
- **NB has larger SE / wider CI** when overdispersion is present

---

## 6. When to use what (quick decision guide)

- Poisson + offset is fine when **Deviance/df ~ 1**
- If **Deviance/df >> 1**, consider:
  - Negative Binomial (best general default)
  - Poisson with robust SE (mild overdispersion)
  - Zero-inflated models (many structural zeros)

---

## Exercises

<details>
<summary>Click to try</summary>

1. Reduce `hetero_sd` until Poisson Deviance/df is near 1. Compare Poisson vs NB outputs.  
2. Simulate a 3-level treatment group and estimate IRRs vs placebo with NB.  
3. Make follow-up times very unequal and verify offset is essential.  
4. Add a second unmeasured subgroup (e.g., half the population has 3× baseline rate) and see overdispersion increase.

</details>

---

## Summary

- NB regression is designed for **overdispersed counts**.
- IRR interpretation is the same as Poisson.
- In biostat data, NB is often more realistic and yields more reliable uncertainty.
