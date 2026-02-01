# 13. Frailty Models (Unobserved Heterogeneity + Clustering) 

In many biostatistics datasets, subjects are **not independent**.

Examples:
- patients clustered within hospitals/centers
- family members share genetics and environment
- repeated patients within the same physician practice
- multi-country clinical trials with site effects

If we ignore clustering, standard Cox regression often gives:

- standard errors too small  
- p-values too optimistic  
- confidence intervals too narrow  

Sometimes we also have **unobserved heterogeneity**:
- two patients with same covariates still have different risk because of unmeasured factors.

Frailty models handle both.

This chapter covers:

- what frailty is  
- shared frailty for clustering  
- individual frailty for unobserved heterogeneity  
- interpretation  
- how to fit in **R** (best support)  
- best practices for **Python** (robust SE + stratification alternatives)  
- what to report in papers  

---

## 1. What is “frailty” in survival analysis?

Frailty is a random effect that multiplies the hazard:

\[
h_i(t \mid X_i, u_i) = u_i \, h_0(t)\exp(\beta^T X_i)
\]

Where:
- \(u_i > 0\) is the frailty (random effect)
- usually \(E[u_i]=1\)

Interpretation:
- \(u_i > 1\): higher risk than average (frailer)
- \(u_i < 1\): lower risk than average

Frailty is like “random risk multiplier” not explained by measured covariates.

---

## 2. Why frailty matters (clinical intuition)

Even after adjusting for:
- age
- sex
- treatment
- stage

some hospitals may have higher mortality because:
- different surgical skill
- ICU quality
- patient mix not fully captured
- protocols

Frailty models capture this unmeasured cluster risk.

---

## 3. Two main frailty settings

### 3.1 Shared frailty 
All subjects in a cluster share the same frailty term.

Example: hospital clustering

\[
h_{ij}(t)=u_j h_0(t)\exp(\beta^T X_{ij})
\]

- \(j\) indexes clusters (hospital)
- \(i\) indexes patients in cluster

Interpretation:
- hospitals have random baseline risk differences

### 3.2 Individual frailty
Each subject has their own frailty (unobserved heterogeneity).

Used more in:
- recurrent events
- unmeasured patient-level risk

---

## 4. Common frailty distributions

Frailty must be positive, so common choices:

### 4.1 Gamma frailty 
- mathematically convenient
- leads to closed-form marginal likelihood in some cases

### 4.2 Log-normal frailty
- random effect on log-hazard scale
- more flexible but computationally heavier

In practice:
- gamma shared frailty is a standard default.

---

## 5. What frailty changes compared to standard Cox

### 5.1 Coefficients
Fixed effects \(\beta\) may change slightly, but often similar.

### 5.2 Standard errors
Usually increase because clustering reduces effective sample size.

### 5.3 Interpretation
Hazard ratios remain conditional on frailty (cluster effect).

---

## 6. Frailty vs robust standard errors (important distinction)

### Robust SE (sandwich)
- adjusts SE for clustering
- does NOT model cluster heterogeneity
- does NOT estimate cluster variance

### Frailty model
- explicitly models clustering as random effect
- estimates variance of frailty
- can improve predictions and interpretation

Rule of thumb:
- if you just need correct inference: robust SE
- if you want to model heterogeneity: frailty

Many papers:
- use robust SE for simplicity
- use frailty when clustering is central

---

# PART A — PYTHON (practical approach)

Python currently has limited frailty support compared to R.
The most common practical options in Python:

- robust standard errors (cluster correction)  
- stratified Cox (different baseline hazards per cluster)  

We show both.

---

## 7A. Simulate clustered survival data (Python)

We create:
- 10 hospitals
- each hospital has its own frailty multiplier
- treatment reduces hazard
- age increases hazard

!!! interactive "Python"
    ```python
    import numpy as np
    import pandas as pd

    np.random.seed(2027)

    n = 600
    hospitals = np.random.choice(range(1, 11), size=n)  # 10 hospitals

    # Gamma frailty-like multipliers per hospital
    # mean 1, some variance
    hospital_effect = {h: np.random.gamma(shape=2.0, scale=0.5) for h in range(1, 11)}

    age = np.random.normal(60, 10, n)
    trt = np.random.binomial(1, 0.5, n)

    beta_age = 0.03
    beta_trt = -0.40
    base = 0.04

    u = np.array([hospital_effect[h] for h in hospitals])

    hazard = u * base * np.exp(beta_age*age + beta_trt*trt)

    T = np.random.exponential(1/hazard)
    C = np.random.uniform(2, 25, n)

    time = np.minimum(T, C)
    event = (T <= C).astype(int)

    df = pd.DataFrame({
        "time": time,
        "event": event,
        "age": age,
        "treatment": trt,
        "hospital": hospitals
    })

    df.head()
    ```

---

## 8A. Fit Cox with and without robust SE (Python)

### 8A.1 Standard Cox (ignores clustering)

!!! interactive "Python"
    ```python
    from lifelines import CoxPHFitter

    cph = CoxPHFitter()
    cph.fit(df[["time","event","age","treatment"]], duration_col="time", event_col="event")
    cph.print_summary()
    ```

### 8A.2 Cox with robust SE

!!! interactive "Python"
    ```python
    cph_robust = CoxPHFitter()
    cph_robust.fit(df[["time","event","age","treatment"]],
                   duration_col="time", event_col="event", robust=True)
    cph_robust.print_summary()
    ```

What you’ll often see:
- coefficients similar
- SE larger with robust=True
- p-values less optimistic

Note:
`robust=True` is a general sandwich correction and does not require specifying cluster ID explicitly in lifelines.

---

## 9A. Stratified Cox by hospital (Python)

Stratification allows each hospital its own baseline hazard:

\[
h(t|X, \text{hospital}) = h_{0,hospital}(t)\exp(\beta^TX)
\]

This adjusts for hospital differences **without estimating an HR for hospital**.

!!! interactive "Python"
    ```python
    cph_strat = CoxPHFitter()
    cph_strat.fit(df, duration_col="time", event_col="event", strata=["hospital"])
    cph_strat.print_summary()
    ```

Interpretation:
- age/treatment HRs are estimated by comparing subjects within the same hospital strata.

Use stratification if:
- you don’t need hospital effect estimates
- you just need to control for hospital baseline hazard differences

---

# PART B — R (true frailty models)

R has excellent frailty support via the `survival` package.

We will show:

- shared frailty (gamma)  
- robust SE alternative  
- interpretation and reporting  

---

## 10B. Simulate clustered survival data (R)

!!! interactive "R"
    ```r
    set.seed(2027)

    n <- 600
    hospital <- sample(1:10, n, replace=TRUE)

    # Gamma frailty multipliers per hospital
    hospital_u <- rgamma(10, shape=2, scale=0.5)  # mean 1
    u <- hospital_u[hospital]

    age <- rnorm(n, 60, 10)
    trt <- rbinom(n, 1, 0.5)

    beta_age <- 0.03
    beta_trt <- -0.40
    base <- 0.04

    hazard <- u * base * exp(beta_age*age + beta_trt*trt)

    T <- rexp(n, rate=hazard)
    C <- runif(n, 2, 25)

    time <- pmin(T, C)
    event <- as.integer(T <= C)

    df <- data.frame(time=time, event=event, age=age, treatment=trt, hospital=factor(hospital))
    head(df)
    ```

---

## 11B. Standard Cox (ignores clustering)

!!! interactive "R"
    ```r
    library(survival)

    fit <- coxph(Surv(time, event) ~ age + treatment, data=df)
    summary(fit)
    ```

---

## 12B. Cox with robust SE for hospital clustering (R)

!!! interactive "R"
    ```r
    fit_robust <- coxph(Surv(time, event) ~ age + treatment + cluster(hospital), data=df)
    summary(fit_robust)
    ```

This corrects SE for clustering but does not estimate frailty variance.

---

## 13B. Shared frailty Cox model (gamma frailty)

This is the real frailty model:

!!! interactive "R"
    ```r
    fit_frailty <- coxph(Surv(time, event) ~ age + treatment + frailty(hospital), data=df)
    summary(fit_frailty)
    ```

Output includes:
- HRs for covariates
- an estimate of frailty variance (often shown as “theta” or variance parameter)
- likelihood ratio test for frailty term

---

## 14B. Interpreting frailty variance

If frailty variance ≈ 0:
- little clustering/unobserved heterogeneity
- frailty may not be needed

If frailty variance is large:
- strong between-hospital heterogeneity
- ignoring clustering likely distorted inference

Often reported as:
- variance estimate
- p-value from LRT comparing with/without frailty

---

## 15. How to interpret covariate HRs in frailty model

The HR is conditional on the frailty term.

Interpretation:
> For two patients in the same hospital (same frailty), treatment HR is X.

If you want marginal population-level effects, interpretation gets more complex.

Most clinical papers interpret conditional HRs.

---

## 16. When to prefer frailty vs stratification

### Stratification:
- controls for cluster baseline differences
- no estimate of cluster variance
- does not assume random effects distribution
- can handle many strata if enough events per stratum

### Frailty:
- models random cluster effect
- estimates cluster heterogeneity
- can improve prediction
- assumes frailty distribution (gamma/lognormal)

Use frailty when:
- clustering is an important scientific feature
- you want to quantify heterogeneity between clusters

---

## 17. Reporting examples (biostat style)

### Robust SE report
> “Cox models used robust sandwich standard errors clustered by hospital.”

### Frailty model report
> “A shared frailty Cox model with hospital-level gamma frailty was fitted to account for within-hospital correlation.”

Then report:
- HRs with CI
- frailty variance estimate
- p-value for frailty term (if relevant)

---

## 18. Common mistakes

### Mistake 1: ignore clustering entirely
Leads to under-estimated SE.

### Mistake 2: treat hospital as a normal covariate with many levels
This wastes degrees of freedom and can be unstable.

Better options:
- strata(hospital)
- cluster(hospital) robust SE
- frailty(hospital)

### Mistake 3: interpret frailty as measured covariate
Frailty is unobserved heterogeneity; not directly measured.

---

## 19. Key takeaways

- Frailty models add a multiplicative random effect to hazard.
- Shared frailty handles clustering (hospital, family, center).
- Robust SE adjusts inference but does not model heterogeneity.
- Stratification controls baseline differences without estimating effects.
- R provides full frailty modeling (`frailty()` in coxph).
- Python best practice: robust SE and/or stratified Cox.

---

## 20. Exercises

<details>
<summary>Click to try</summary>

1. Simulate data with strong hospital frailty and fit: standard Cox vs robust SE vs frailty. Compare SE and p-values.  
2. Increase frailty variance by changing gamma distribution parameters; observe stronger clustering effects.  
3. Fit stratified Cox by hospital and compare treatment HR to frailty model.  
4. Interpret frailty variance in plain language.  
5. Write a short “Methods” paragraph describing a frailty model analysis for a multi-center trial.

</details>
