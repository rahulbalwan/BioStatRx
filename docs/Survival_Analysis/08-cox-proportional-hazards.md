# Cox Proportional Hazards Model (Cox PH) 

Kaplan–Meier and log-rank help you **describe** and **compare** survival curves.

But in real biostatistics, we almost always need to answer:

- How does **treatment** affect survival **after adjusting** for age/sex/stage?
- What is the effect size (hazard ratio) with confidence intervals?
- How do multiple predictors jointly influence risk?

That is exactly what the **Cox proportional hazards model** does.

# Cox PH is the workhorse model of survival analysis.

---

## 1. What Cox regression models (the hazard)

Cox models the **hazard function**:

\[
h(t \mid X) = h_0(t)\exp(\beta_1X_1+\cdots+\beta_pX_p)
\]

Where:

- \(h(t \mid X)\): hazard at time \(t\) for covariate vector \(X\)
- \(h_0(t)\): baseline hazard (unknown function of time)
- \(\exp(\beta^TX)\): multiplicative effect of covariates

### 1.1 Key idea: multiplicative risk
Covariates scale risk *multiplicatively*:

- If \(\exp(\beta)=1.5\), hazard is 1.5× higher (50% higher risk rate).
- If \(\exp(\beta)=0.7\), hazard is 0.7× (30% lower risk rate).

---

## 2. Why Cox is called “semi-parametric”

Cox does **not** assume a parametric form for survival times (like exponential or Weibull).

- It leaves \(h_0(t)\) unspecified.
- It estimates \(\beta\) without needing to specify \(h_0(t)\).

So it’s:
- **parametric in covariates** (linear predictor \(\beta^TX\))
- **nonparametric in baseline hazard** (\(h_0(t)\))

This is why Cox is flexible and widely used.

---

## 3. The hazard ratio (HR) — the main output

For a 1-unit increase in predictor \(X_j\):

\[
HR = \exp(\beta_j)
\]

### 3.1 Clinical interpretation

- \(HR = 1\): no difference in hazard
- \(HR > 1\): higher hazard → worse survival (typically)
- \(HR < 1\): lower hazard → better survival (typically)

Examples:
- **HR = 1.20**: 20% higher instantaneous risk
- **HR = 0.65**: 35% lower instantaneous risk

### 3.2 Continuous predictors (age)
If HR for age is 1.04:

> Each extra year of age increases hazard by ~4% (assuming linear effect).

For a 10-year increase:
\[
HR_{10} = 1.04^{10}\approx 1.48
\]
So ~48% higher hazard.

---

## 4. The proportional hazards (PH) assumption (critical)

Cox assumes hazard ratios are constant over time:

\[
\frac{h(t\mid X_a)}{h(t\mid X_b)} = \exp(\beta(X_a-X_b))
\]

No \(t\) in the ratio ⇒ constant over time.

### 4.1 What PH means in plain language

If treatment HR = 0.70, then:

> Treatment reduces hazard by 30% at all times (early, middle, late).

### 4.2 What PH does NOT allow
If treatment helps early but harms later (curves cross):

- PH is violated
- Cox HR becomes difficult to interpret

We handle this later using:
- time interactions
- stratified Cox
- time-varying effects

---

## 5. Where Cox comes from: partial likelihood (intuition)

At each event time \(t_j\), Cox compares:

- the subject who had the event at \(t_j\)
- to everyone in the risk set \(R(t_j)\) (those who could have had the event at that time)

The probability that subject \(i\) fails at \(t_j\) given one event occurs then:

\[
\frac{\exp(\beta^TX_i)}{\sum_{k\in R(t_j)}\exp(\beta^TX_k)}
\]

Multiply over event times:

\[
L(\beta)=\prod_j \frac{\exp(\beta^TX_{(j)})}{\sum_{k\in R(t_j)}\exp(\beta^TX_k)}
\]

This is the **partial likelihood**.

Key insight:

We estimate \(\beta\) without ever specifying \(h_0(t)\).

---

## 6. Data structure for Cox

Minimum variables:

- `time`: observed follow-up time (event or censor)
- `event`: 1 if event, 0 if censored
- predictors: age, sex, treatment, stage, etc.

Example:

| time | event | age | trt |
|-----:|------:|----:|----:|
|  4.2 | 1 | 63 | 1 |
|  7.9 | 0 | 55 | 0 |
|  2.1 | 1 | 71 | 1 |

---

## 7. A complete simulation with known truth (Python + R)

We simulate a study where:
- age increases hazard
- treatment reduces hazard
- censoring occurs

Then we fit Cox and see if we recover the true effects.

---

# 7A. Python (lifelines) — Cox PH in practice

## 7A.1 Simulate survival data

!!! interactive "Python"
    ```python
    import numpy as np
    import pandas as pd

    np.random.seed(2026)

    n = 400

    # Covariates
    age = np.random.normal(60, 10, n)          # continuous predictor
    trt = np.random.binomial(1, 0.5, n)        # 0=control, 1=treatment

    # True coefficients (we will try to recover these)
    beta_age = 0.03        # HR per year ~ exp(0.03)=1.030
    beta_trt = -0.50       # HR treatment vs control ~ exp(-0.50)=0.607

    # Baseline hazard scale (controls overall event rate)
    base = 0.04

    # Individual hazard rates (proportional hazards)
    hazard = base * np.exp(beta_age*age + beta_trt*trt)

    # True event times from exponential with individual hazards
    T = np.random.exponential(1/hazard)

    # Random censoring times
    C = np.random.uniform(2, 25, n)

    # Observed data
    time = np.minimum(T, C)
    event = (T <= C).astype(int)

    df = pd.DataFrame({"time": time, "event": event, "age": age, "treatment": trt})
    df.head()
    ```

---

## 7A.2 Fit Cox PH model

!!! interactive "Python"
    ```python
    from lifelines import CoxPHFitter

    cph = CoxPHFitter()
    cph.fit(df, duration_col="time", event_col="event")

    cph.print_summary()
    ```

### What to look for
- `coef` ≈ true beta
- `exp(coef)` ≈ true hazard ratio
- p-values often significant if effect is strong and sample size decent

---

## 7A.3 Interpret output (common columns)

Typical Cox output includes:

- `coef` = \(\beta\)
- `exp(coef)` = hazard ratio
- `se(coef)` = standard error
- `p` = p-value testing \(\beta=0\)
- CI for HR

Interpretation example:
- exp(coef) = 0.61 for treatment → ~39% hazard reduction

---

## 7A.4 Predict survival curves for profiles

Cox can produce predicted survival curves using estimated baseline + covariates.

!!! interactive "Python"
    ```python
    import matplotlib.pyplot as plt
    import pandas as pd

    profile_control = pd.DataFrame({"age":[60], "treatment":[0]})
    profile_treated = pd.DataFrame({"age":[60], "treatment":[1]})

    s0 = cph.predict_survival_function(profile_control)
    s1 = cph.predict_survival_function(profile_treated)

    plt.plot(s0, label="Control (age=60)")
    plt.plot(s1, label="Treatment (age=60)")
    plt.legend()
    plt.title("Predicted Survival Curves from Cox PH")
    plt.xlabel("Time")
    plt.ylabel("S(t)")
    plt.show()
    ```

---

## 7A.5 Multivariable Cox example (add more covariates)

!!! interactive "Python"
    ```python
    np.random.seed(7)

    df2 = df.copy()
    df2["sex"] = np.random.binomial(1, 0.5, len(df2))  # 0=female,1=male
    df2["bmi"] = np.random.normal(27, 4, len(df2))

    cph2 = CoxPHFitter()
    cph2.fit(df2, duration_col="time", event_col="event")

    cph2.print_summary()
    ```

Interpretation:
- Each coefficient is adjusted for the others.

---

# 7B. R (survival) — Cox PH in practice

## 7B.1 Simulate survival data in R 

!!! interactive "R"
    ```r
    set.seed(2026)

    n <- 400

    age <- rnorm(n, mean=60, sd=10)
    trt <- rbinom(n, 1, 0.5)

    beta_age <- 0.03
    beta_trt <- -0.50
    base <- 0.04

    hazard <- base * exp(beta_age*age + beta_trt*trt)

    # event times
    T <- rexp(n, rate = hazard)

    # censoring times
    C <- runif(n, min=2, max=25)

    time <- pmin(T, C)
    event <- as.integer(T <= C)

    df <- data.frame(time=time, event=event, age=age, treatment=trt)
    head(df)
    ```

---

## 7B.2 Fit Cox model in R

!!! interactive "R"
    ```r
    library(survival)

    fit <- coxph(Surv(time, event) ~ age + treatment, data=df)
    summary(fit)
    ```

### What to look for in `summary(fit)`
- `coef` = \(\beta\)
- `exp(coef)` = hazard ratio
- `se(coef)` = standard error
- Wald test p-values
- 95% CI for HR

---

## 7B.3 Predicted survival curves in R (for profiles)

In R, prediction is typically done using `survfit()` on a Cox model:

!!! interactive "R"
    ```r
    newdata <- data.frame(age=c(60,60), treatment=c(0,1))

    sf <- survfit(fit, newdata=newdata)

    plot(sf, col=c("black","red"), lty=1,
         xlab="Time", ylab="S(t)",
         main="Predicted Survival Curves from Cox PH")
    legend("topright", legend=c("Control age=60","Treatment age=60"),
           col=c("black","red"), lty=1)
    ```

---

## 7B.4 Multivariable Cox example in R

!!! interactive "R"
    ```r
    set.seed(7)
    df$sex <- rbinom(nrow(df), 1, 0.5)
    df$bmi <- rnorm(nrow(df), 27, 4)

    fit2 <- coxph(Surv(time, event) ~ age + treatment + sex + bmi, data=df)
    summary(fit2)
    ```

---

## 8. Model interpretation: what Cox answers (and what it doesn’t)

### Cox DOES answer:
- association between covariates and hazard  
- adjusted hazard ratios  
- hypothesis tests for covariates  
- predicted survival curves for profiles (with baseline estimation)

### Cox DOES NOT automatically answer:
- causal effect unless design/assumptions justify causality  
- what happens when hazards are non-proportional  
- competing risks probabilities (need CIF)  
- repeated events without special modeling

---

## 9. Checking proportional hazards (preview)

You must check PH assumption before reporting results.

In later chapters we will do this deeply using:

- Schoenfeld residuals
- PH tests
- plots

Quick preview:

### Python (lifelines)
!!! interactive "Python"
    ```python
    cph.check_assumptions(df, show_plots=False)
    ```

### R (survival)
!!! interactive "R"
    ```r
    cox.zph(fit)
    ```

Interpretation:
- Small p-value suggests PH violation.

---

## 10. How to report Cox results 

A standard reporting sentence:

> “In a multivariable Cox proportional hazards model adjusting for age and sex, treatment was associated with improved survival (HR 0.61, 95% CI 0.48–0.77, p < 0.001).”

Key elements:
 - model type (Cox PH)
 - adjusted covariates
 - HR, CI, p-value
 - direction and clinical interpretation

---

## 11. Common mistakes 

### Mistake 1: interpreting HR as a probability
 HR is a ratio of hazard rates, not survival probability.

 Fix: interpret as instantaneous risk rate ratio.

### Mistake 2: ignoring PH assumption
 If PH violated, HR may change over time and “one HR” can mislead.

 Fix: check PH; use time-varying effects if needed.

### Mistake 3: using Cox for competing risks outcome
 If competing risks present, standard Cox does not directly estimate cumulative incidence.

 Fix: use competing risks methods (CIF, Fine–Gray).

### Mistake 4: nonlinear continuous covariates
 Age effect might not be linear in log-hazard.

 Fix: consider splines / transformations.

---

## 12. Key takeaways

 - Cox PH models hazard: \(h(t|X)=h_0(t)\exp(\beta^TX)\).
 - Main output: hazard ratio \(\exp(\beta)\).
 - Proportional hazards assumption: HR constant over time.
 - Cox estimates \(\beta\) using partial likelihood and risk sets.
 - Always check PH before interpreting results.
 - Cox is the standard regression tool for survival analysis.

---

## 13. Exercises

<details>
<summary>Click to try</summary>

 1. Simulate data where treatment has no effect (set beta_trt = 0). Fit Cox. What HR do you estimate?  
 2. Interpret HR = 1.25 in plain language.  
 3. Compute HR for a 10-year age increase when HR per year is 1.04.  
 4. Fit Cox in both Python and R and compare coefficients/HR.  
 5. Simulate a scenario where treatment effect changes over time (early benefit, late harm) and check PH diagnostics (preview with cox.zph / check_assumptions).

</details>
