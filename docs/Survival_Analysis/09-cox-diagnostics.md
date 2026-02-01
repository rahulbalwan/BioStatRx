# 09. Cox Model Diagnostics and Assumption Checking

A Cox model is **not complete** when you fit it.

In real biostatistical work, the correct workflow is:

**Fit Cox** → **Check assumptions** → **Fix violations** → **Interpret HRs**

If you skip diagnostics, your hazard ratios (HRs) may be:
- biased
- uninterpretable
- misleading (especially under non-proportional hazards)

This chapter is a **complete diagnostic toolkit** for Cox regression in biostatistics, with **both Python and R**.

---

## 1. What assumptions does Cox PH rely on?

The Cox proportional hazards model is:

\[
h(t|X)=h_0(t)\exp(\beta^TX)
\]

Main assumptions you must check:

### A1) Proportional hazards (PH) most important
Hazard ratios are constant over time:
\[
\frac{h(t|X_a)}{h(t|X_b)}=\exp(\beta^T(X_a-X_b)) \quad \text{does not depend on } t
\]

### A2) Correct functional form for continuous predictors
Cox assumes a linear relationship in the **log-hazard**:
\[
\log h(t|X) = \log h_0(t) + \beta_1 X_1 + \cdots
\]
So a continuous predictor effect is assumed **linear on log-hazard scale** unless you model nonlinearity.

### A3) Independent observations (or correctly handled clustering)
Patients in same hospital/site/family may be correlated.

### A4) No extreme outliers or overly influential observations
One or two subjects shouldn’t dominate your HR estimates.

---

## 2. A practical diagnostic checklist (memorize)

After fitting Cox:

1) **KM curves / visual check**
2) **PH test + Schoenfeld residual plots**
3) **Log(-log) survival plots (group PH check)**
4) **Functional form check for continuous predictors**
5) **Influential observations (dfbeta, deviance residuals)**
6) **Overall fit (concordance / C-index)**
7) **Robust SE / clustering check**
8) **Sensitivity analyses** (alternative models if needed)

---

## 3. Setup: create a working dataset (Python + R)

We will simulate data with:
- age (continuous)
- treatment (binary)
- optional PH violation example later

### 3.1 Python: simulate data

!!! interactive "Python"
    ```python
    import numpy as np
    import pandas as pd

    np.random.seed(2026)

    n = 500
    age = np.random.normal(60, 10, n)
    trt = np.random.binomial(1, 0.5, n)

    beta_age = 0.03
    beta_trt = -0.50
    base = 0.04

    hazard = base * np.exp(beta_age*age + beta_trt*trt)
    T = np.random.exponential(1/hazard)
    C = np.random.uniform(2, 25, n)

    time = np.minimum(T, C)
    event = (T <= C).astype(int)

    df = pd.DataFrame({"time": time, "event": event, "age": age, "treatment": trt})
    df.head()
    ```

### 3.2 R: simulate the same idea

!!! interactive "R"
    ```r
    set.seed(2026)

    n <- 500
    age <- rnorm(n, mean=60, sd=10)
    trt <- rbinom(n, 1, 0.5)

    beta_age <- 0.03
    beta_trt <- -0.50
    base <- 0.04

    hazard <- base * exp(beta_age*age + beta_trt*trt)
    T <- rexp(n, rate = hazard)
    C <- runif(n, min=2, max=25)

    time <- pmin(T, C)
    event <- as.integer(T <= C)

    df <- data.frame(time=time, event=event, age=age, treatment=trt)
    head(df)
    ```

---

## 4. Fit Cox model (baseline step)

### 4.1 Python: fit Cox (lifelines)

!!! interactive "Python"
    ```python
    from lifelines import CoxPHFitter

    cph = CoxPHFitter()
    cph.fit(df, duration_col="time", event_col="event")
    cph.print_summary()
    ```

### 4.2 R: fit Cox (survival)

!!! interactive "R"
    ```r
    library(survival)

    fit <- coxph(Surv(time, event) ~ age + treatment, data=df)
    summary(fit)
    ```

---

## 5. PH assumption checks (the critical part)

### 5.1 Visual pre-check: KM curves (group separation and crossing)

If curves strongly cross → PH likely violated.

#### Python (KM curves)

!!! interactive "Python"
    ```python
    import matplotlib.pyplot as plt
    from lifelines import KaplanMeierFitter

    km0 = KaplanMeierFitter()
    km1 = KaplanMeierFitter()

    ax = plt.subplot(111)
    km0.fit(df.loc[df.treatment==0, "time"], df.loc[df.treatment==0, "event"], label="Control").plot(ax=ax)
    km1.fit(df.loc[df.treatment==1, "time"], df.loc[df.treatment==1, "event"], label="Treatment").plot(ax=ax)

    plt.title("KM Curves by Treatment (quick PH visual)")
    plt.xlabel("Time")
    plt.ylabel("S(t)")
    plt.show()
    ```

#### R (KM curves)

!!! interactive "R"
    ```r
    fit_km <- survfit(Surv(time, event) ~ factor(treatment), data=df)

    plot(fit_km, col=c("black","red"), lty=1,
         xlab="Time", ylab="S(t)", main="KM Curves by Treatment")
    legend("topright", legend=c("Control","Treatment"), col=c("black","red"), lty=1)
    ```

- If curves are roughly separated without major crossing, PH is plausible.  
- If curves cross substantially, investigate PH violation formally.

---

### 5.2 Schoenfeld residual test (formal PH test)

#### R: `cox.zph()` (gold standard)

!!! interactive "R"
    ```r
    z <- cox.zph(fit)
    z
    plot(z)
    ```

Interpretation:
- Each covariate gets a test.
- Global test also reported.
- **Small p-value (<0.05)** → evidence PH violation.

The plot shows residual trend over time:
- Flat trend around 0 → good
- Clear slope / curve → violation

#### Python: lifelines `check_assumptions()`

!!! interactive "Python"
    ```python
    cph.check_assumptions(df, show_plots=True)
    ```

Interpretation:
- lifelines prints warnings when PH is violated
- shows plots (Schoenfeld-type diagnostics)

---

### 5.3 Log(-log) survival plots (classic PH visual check for groups)

If PH holds, log(-log(S(t))) curves for groups should be **roughly parallel**.

#### R: log(-log) plot

!!! interactive "R"
    ```r
    library(survival)

    fit_km <- survfit(Surv(time, event) ~ factor(treatment), data=df)

    # Extract survival curves
    s <- summary(fit_km)

    # Compute log(-log(S))
    lls <- log(-log(s$surv))

    plot(s$time, lls, type="n",
         xlab="Time", ylab="log(-log(S(t)))",
         main="Log(-log) Survival Plot (PH visual check)")

    # Two groups are stored sequentially; split by strata
    strata_lengths <- s$strata
    # A simpler approach: use survminer if available:
    ```

If you have `survminer`, this is easier:

!!! interactive "R"
    ```r
    library(survminer)

    ggsurvplot(fit_km, fun="cloglog", conf.int=FALSE,
               ggtheme=theme_minimal(),
               title="Cloglog plot: log(-log(S(t)))")
    ```

Interpretation:
- Parallel lines → PH plausible
- Non-parallel / crossing → PH questionable

Python doesn’t have a one-liner in lifelines for cloglog plots, but you can approximate by extracting survival estimates and plotting `np.log(-np.log(S))`.

---

## 6. What to do if PH is violated (practical fixes)

PH violation is common. You have multiple solutions depending on the situation.

### Option 1: Add **time interaction** (time-varying effect)
If treatment effect changes over time:

\[
h(t)=h_0(t)\exp(\beta_1X + \beta_2 X\cdot g(t))
\]

A simple choice: \(g(t)=\log(t)\).

#### R example: time interaction with `tt()`

!!! interactive "R"
    ```r
    fit_tv <- coxph(
      Surv(time, event) ~ age + treatment + tt(treatment),
      data=df,
      tt = function(x, t, ...) x * log(t + 1)
    )
    summary(fit_tv)
    ```

Interpretation:
- Now treatment effect depends on time.
- You don’t report “one HR”; you interpret HR as a function of time.

#### Python approach
lifelines can do time-varying covariates using **start–stop format** (covered in the time-dependent chapter). For PH violation, you can create an interaction variable such as `treatment * log(time)` and refit, but that is less principled unless properly formatted.

---

### Option 2: Stratified Cox
If a variable violates PH but you do NOT need its HR, stratify by it.

\[
h(t|X)=h_{0,stratum}(t)\exp(\beta X)
\]

This allows different baseline hazards per stratum.

#### R: stratification

!!! interactive "R"
    ```r
    # Example: stratify by treatment if you don't need HR for it
    fit_strat <- coxph(Surv(time, event) ~ age + strata(treatment), data=df)
    summary(fit_strat)
    ```

You lose HR estimate for the stratified variable, but PH becomes less problematic.

---

### Option 3: Use alternative summaries (RMST) or parametric models
If PH is badly violated:
- Consider **restricted mean survival time (RMST)** (time-based effect)
- Consider parametric models that allow more flexible hazards

(We can add a dedicated RMST chapter later if you want.)

---

## 7. Checking functional form (linearity) for continuous predictors

Cox assumes:
\[
\log h(t|X) \text{ is linear in } X
\]

If age effect is nonlinear (common!), a linear age term is wrong.

### 7.1 Martingale residual plots (classic approach)

#### R: martingale residual vs age

!!! interactive "R"
    ```r
    mart <- residuals(fit, type="martingale")

    plot(df$age, mart,
         xlab="Age", ylab="Martingale residual",
         main="Martingale Residuals vs Age")
    lines(lowess(df$age, mart), col="red", lwd=2)
    ```

Interpretation:
- random scatter around 0 → linear okay
- curved LOWESS trend → nonlinearity

#### Python (lifelines): martingale residuals

!!! interactive "Python"
    ```python
    mart = cph.compute_residuals(df, kind="martingale")

    import matplotlib.pyplot as plt
    plt.scatter(df["age"], mart, alpha=0.6)
    plt.xlabel("Age")
    plt.ylabel("Martingale residual")
    plt.title("Martingale Residuals vs Age")
    plt.show()
    ```

*(If your lifelines version returns a DataFrame, you may need `mart.values`.)*

### 7.2 Fixing nonlinearity: splines or categories

#### R: natural spline with `splines`

!!! interactive "R"
    ```r
    library(splines)

    fit_spline <- coxph(Surv(time, event) ~ ns(age, df=4) + treatment, data=df)
    summary(fit_spline)
    ```

#### Python: transform age using piecewise terms or spline libraries
In pure lifelines, easiest is to:
- create polynomial terms (`age`, `age^2`)
- or bin age categories
For full spline modeling, many analysts use `patsy` or `statsmodels` to build spline basis and then feed columns to lifelines.

---

## 8. Outliers and influential observations

Even if PH holds, a few subjects can dominate the model.

### 8.1 Deviance residuals (outlier detection)

#### R

!!! interactive "R"
    ```r
    dev <- residuals(fit, type="deviance")

    hist(dev, breaks=30, main="Deviance residuals", xlab="Deviance residual")
    ```

Large absolute values indicate potential outliers.

#### Python

!!! interactive "Python"
    ```python
    dev = cph.compute_residuals(df, kind="deviance")

    import matplotlib.pyplot as plt
    plt.hist(dev, bins=30)
    plt.title("Deviance residuals")
    plt.xlabel("Deviance residual")
    plt.show()
    ```

### 8.2 Influence: dfbeta (how much each subject changes coefficients)

#### R: dfbeta

!!! interactive "R"
    ```r
    dfb <- residuals(fit, type="dfbeta")

    # dfbeta is a matrix with one column per covariate
    matplot(dfb, type="p", pch=20, cex=0.6,
            main="DFBETA by covariate", ylab="DFBETA")
    abline(h=0, col="grey")
    ```

Subjects with extreme DFBETA values may have high influence.

#### Python: influence is less standardized in lifelines
lifelines doesn’t offer a perfect `dfbeta` equivalent in all versions. A practical approach:
- fit model
- refit leaving-one-out for suspected points (expensive but possible)
- or inspect deviance residual extremes and investigate those IDs

---

## 9. Overall model fit and predictive ability (C-index)

Cox does not have an R² like linear regression.
Instead, common measure is:

### Concordance index (C-index)
- 0.5 = random prediction
- 1.0 = perfect discrimination

#### Python

!!! interactive "Python"
    ```python
    cph.concordance_index_
    ```

#### R
`summary(fit)` reports concordance.
Or extract:

!!! interactive "R"
    ```r
    summary(fit)$concordance
    ```

Interpretation:
- Higher is better, but values around 0.60–0.75 are common in medical data.

---

## 10. Robust standard errors and clustering

If you have clustering (hospital/site/family), independence fails.

### 10.1 Robust SE (sandwich) — quick fix

#### R
Use `cluster()` to get robust SE for clustered IDs.

!!! interactive "R"
    ```r
    # Suppose we have a cluster variable (e.g., hospital)
    set.seed(1)
    df$hospital <- sample(1:10, nrow(df), replace=TRUE)

    fit_cluster <- coxph(Surv(time, event) ~ age + treatment + cluster(hospital), data=df)
    summary(fit_cluster)
    ```

#### Python
lifelines supports `robust=True`:

!!! interactive "Python"
    ```python
    cph_robust = CoxPHFitter()
    cph_robust.fit(df, duration_col="time", event_col="event", robust=True)
    cph_robust.print_summary()
    ```

Robust SE corrects inference (SE, CI), but does not model heterogeneity explicitly (frailty models do).

---

## 11. Practical “PH violation” simulation (see it happen!)

To understand PH diagnostics, it helps to simulate data where treatment effect changes over time.

### 11.1 R: simulate crossing hazards and test PH

!!! interactive "R"
    ```r
    set.seed(99)
    n <- 400
    trt <- rbinom(n, 1, 0.5)

    # Create time-varying effect: treatment helps early, harms late
    # We'll simulate piecewise hazards:
    # early hazard: treatment reduces
    # late hazard: treatment increases
    t_change <- 8

    base1 <- 0.06
    base2 <- 0.03

    # generate event time using simple rejection / piecewise approx
    T <- rep(NA, n)
    for (i in 1:n) {
      # simulate early time
      h1 <- base1 * ifelse(trt[i]==1, 0.6, 1.0)
      t1 <- rexp(1, rate=h1)
      if (t1 < t_change) {
        T[i] <- t1
      } else {
        # survive past change point
        h2 <- base2 * ifelse(trt[i]==1, 1.6, 1.0)
        t2 <- rexp(1, rate=h2)
        T[i] <- t_change + t2
      }
    }

    C <- runif(n, 2, 20)
    time <- pmin(T, C)
    event <- as.integer(T <= C)

    d <- data.frame(time=time, event=event, trt=trt)

    fit_np <- coxph(Surv(time, event) ~ trt, data=d)
    summary(fit_np)

    z <- cox.zph(fit_np)
    z
    plot(z)
    ```

You should often see PH violation for `trt`.

### 11.2 Python: simulate similar idea and run check_assumptions

!!! interactive "Python"
    ```python
    import numpy as np
    import pandas as pd
    from lifelines import CoxPHFitter

    np.random.seed(99)
    n = 400
    trt = np.random.binomial(1, 0.5, n)

    t_change = 8.0
    base1 = 0.06
    base2 = 0.03

    T = np.zeros(n)

    for i in range(n):
        h1 = base1 * (0.6 if trt[i]==1 else 1.0)
        t1 = np.random.exponential(1/h1)
        if t1 < t_change:
            T[i] = t1
        else:
            h2 = base2 * (1.6 if trt[i]==1 else 1.0)
            t2 = np.random.exponential(1/h2)
            T[i] = t_change + t2

    C = np.random.uniform(2, 20, n)
    time = np.minimum(T, C)
    event = (T <= C).astype(int)

    d = pd.DataFrame({"time": time, "event": event, "treatment": trt})

    cph_bad = CoxPHFitter()
    cph_bad.fit(d, duration_col="time", event_col="event")
    cph_bad.print_summary()

    cph_bad.check_assumptions(d, show_plots=True)
    ```

You should observe warnings/plots indicating PH violation.

---

## 12. Reporting diagnostics in a paper 

Examples of correct writing:

### If PH holds
> “Proportional hazards assumption was assessed using Schoenfeld residuals and was not violated.”

### If PH violated and fixed
> “The proportional hazards assumption was violated for treatment (Schoenfeld p < 0.05). A time-varying treatment effect was modeled using an interaction with log(time).”

### If stratified
> “Due to PH violation for sex, models were stratified by sex.”

---

## 13. Common mistakes (very common in student work)

### Mistake 1: Interpret HR without checking PH
Fix: Always run `cox.zph()` (R) or `check_assumptions()` (Python).

### Mistake 2: Assume linear age effect automatically
Fix: martingale residual check, use splines if needed.

### Mistake 3: Ignore clustering
Fix: robust SE or frailty.

### Mistake 4: Trust extreme late follow-up
Fix: check number at risk and consider truncating analysis window.

---

## 14. Key takeaway summary

Cox diagnostics ensure:

- HRs are interpretable (PH)
- continuous covariates are modeled correctly (linearity)
- inference is valid (robust SE/clustering)
- results are not dominated by a few subjects (influence)
- model has acceptable predictive performance (C-index)

If you remember only one thing:

# Fit Cox → always check PH → then interpret.

---

## 15. Exercises (high-value practice)

<details>
<summary>Click to try</summary>

1. Fit Cox model in R and run `cox.zph()`. Which variables violate PH (if any)?  
2. Make a dataset where treatment effect changes over time and confirm PH violation.  
3. Use martingale residual plot to detect nonlinearity in age. Then refit with spline.  
4. Create clustering (hospital variable) and compare standard SE vs robust SE.  
5. Identify influential points via deviance residuals/dfbeta and refit after removing a few. Compare HR changes.

</details>
