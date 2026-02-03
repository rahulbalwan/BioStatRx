# Parametric Survival Models (Exponential, Weibull, Log-normal, AFT) 

So far you’ve learned:

- Kaplan–Meier (nonparametric)
- Cox regression (semi-parametric)

Both are powerful, but sometimes we need:

- smooth survival curves  
- explicit survival formulas  
- extrapolation beyond observed follow-up  
- simulation and prediction  
- health-economic modeling  

For these, we use:

# Parametric Survival Models

Parametric survival models assume survival times follow a known distribution:

\[
T \sim \text{some distribution}
\]

Examples:
- Exponential
- Weibull
- Log-normal
- Log-logistic

This chapter gives maximum clarity + **interactive code in Python and R**.

---

## 1. Big picture: why parametric models?

### 1.1 Limitations of Cox
Cox does not specify \(h_0(t)\), so:

- baseline hazard unknown
- extrapolation beyond observed time is unreliable
- simulations require extra steps

### 1.2 Advantages of parametric models
Parametric models give:

- closed-form \(S(t)\) and \(h(t)\)  
- smoother and more stable estimates  
- extrapolation (common in HTA, cost-effectiveness)  
- efficient estimation when distribution fits well  

### 1.3 Trade-off
You gain structure but lose flexibility:

- If distribution is wrong → bias
- Must check fit carefully

---

## 2. Two modeling frameworks: PH vs AFT

Many parametric models can be expressed as:

### 2.1 Proportional Hazards (PH) form
Similar interpretation to Cox:

\[
h(t|X)=h_0(t)\exp(\beta^TX)
\]

Hazard ratio is constant over time.

Common PH parametric models:
- Exponential PH
- Weibull PH

### 2.2 Accelerated Failure Time (AFT) form
Models time directly:

\[
\log(T)=\beta^TX+\epsilon
\]

Effect is interpreted using **time ratios**:

\[
TR=\exp(\beta)
\]

Meaning:
- TR > 1 → longer survival time
- TR < 1 → shorter survival time

Common AFT models:
- Weibull AFT
- Log-normal AFT
- Log-logistic AFT

---

## 3. Exponential model (the simplest)

### 3.1 Assumption
Hazard is constant:

\[
h(t)=\lambda
\]

### 3.2 Survival function
\[
S(t)=\exp(-\lambda t)
\]

### 3.3 Interpretation
- risk does not change over time
- often unrealistic in humans

Clinical contexts where exponential might be plausible:
- device failure with constant rate
- short follow-up where hazard seems flat

---

## 4. Weibull model (most important parametric model)

Weibull is the survival distribution you will see most often in biostatistics.

### 4.1 Weibull survival
\[
S(t)=\exp\left(-\left(\frac{t}{\lambda}\right)^k\right)
\]

- \(k\) = shape
- \(\lambda\) = scale

### 4.2 Weibull hazard
\[
h(t)=\frac{k}{\lambda}\left(\frac{t}{\lambda}\right)^{k-1}
\]

### 4.3 Why Weibull is special
Weibull can model different hazard shapes:

- \(k=1\) → constant hazard (exponential)
- \(k>1\) → increasing hazard
- \(k<1\) → decreasing hazard

So Weibull is flexible and extremely useful.

---

## 5. Log-normal model

### 5.1 Assumption
\[
\log(T)\sim N(\mu,\sigma^2)
\]

### 5.2 Hazard shape
Log-normal hazard is non-monotonic:

- rises early
- peaks
- declines

Useful when:
- early risk is high (post-surgery)
- later risk declines

---

## 6. Log-logistic model

Similar to log-normal but with heavier tails.

- non-monotonic hazard
- sometimes better for very long survivors

---

## 7. Choosing among parametric models (practical approach)

### 7.1 Fit multiple models and compare
Common criteria:
- AIC
- BIC
- likelihood
- visual fit to KM

### 7.2 General rule of thumb
- start with Weibull
- compare with log-normal/log-logistic
- consider exponential only if hazard looks flat

---

# PART A — PYTHON (lifelines)

Python parametric survival modeling is commonly done using `lifelines`.

We will:
1) simulate survival data  
2) fit exponential / Weibull / log-normal  
3) compare AIC  
4) plot model fits against KM  
5) fit AFT regression model  

---

## 8A. Simulate survival data (Python)

We simulate Weibull survival with covariate effect:

- age increases hazard
- treatment improves survival

!!! interactive "Python"
    ```python
    import numpy as np
    import pandas as pd

    np.random.seed(44)

    n = 400
    age = np.random.normal(60, 10, n)
    trt = np.random.binomial(1, 0.5, n)

    # True AFT-type generation (log-time shift)
    # We'll create baseline Weibull event times then accelerate by covariates
    shape_k = 1.7
    scale_lam = 12.0

    # Generate baseline Weibull times
    U = np.random.uniform(size=n)
    T0 = scale_lam * (-np.log(U))**(1/shape_k)

    # Covariate effects on log-time (AFT): treatment increases time, age decreases
    beta_trt = 0.35   # TR = exp(0.35)=1.42
    beta_age = -0.01  # per year TR = exp(-0.01)=0.99

    T = T0 * np.exp(beta_trt*trt + beta_age*(age-60))

    # Censoring
    C = np.random.uniform(3, 20, n)

    time = np.minimum(T, C)
    event = (T <= C).astype(int)

    df = pd.DataFrame({"time": time, "event": event, "age": age, "treatment": trt})
    df.head()
    ```

---

## 9A. Fit nonparametric KM (reference curve)

!!! interactive "Python"
    ```python
    import matplotlib.pyplot as plt
    from lifelines import KaplanMeierFitter

    km = KaplanMeierFitter().fit(df["time"], df["event"], label="KM")

    ax = km.plot(ci_show=False)
    plt.title("KM Curve (Reference)")
    plt.xlabel("Time")
    plt.ylabel("S(t)")
    plt.show()
    ```

---

## 10A. Fit parametric models (Python)

### 10A.1 Exponential, Weibull, Log-normal

!!! interactive "Python"
    ```python
    from lifelines import ExponentialFitter, WeibullFitter, LogNormalFitter

    expf = ExponentialFitter().fit(df["time"], df["event"])
    wf   = WeibullFitter().fit(df["time"], df["event"])
    lnf  = LogNormalFitter().fit(df["time"], df["event"])

    print("AIC comparison (lower is better):")
    print("Exponential AIC:", expf.AIC_)
    print("Weibull AIC    :", wf.AIC_)
    print("LogNormal AIC  :", lnf.AIC_)
    ```

---

## 10A.2 Plot fitted curves vs KM

!!! interactive "Python"
    ```python
    import matplotlib.pyplot as plt

    ax = km.plot(ci_show=False)

    expf.plot_survival_function(ax=ax, ci_show=False)
    wf.plot_survival_function(ax=ax, ci_show=False)
    lnf.plot_survival_function(ax=ax, ci_show=False)

    plt.title("KM vs Parametric Fits")
    plt.xlabel("Time")
    plt.ylabel("S(t)")
    plt.show()
    ```

Interpretation:
- Best-fitting model tends to track KM curve more closely.

---

## 11A. Regression: Weibull AFT model (Python)

Lifelines supports AFT regression models.

!!! interactive "Python"
    ```python
    from lifelines import WeibullAFTFitter

    aft = WeibullAFTFitter()
    aft.fit(df, duration_col="time", event_col="event")
    aft.print_summary()
    ```

### 11A.1 Interpretation: Time Ratios
In AFT models, exponentiated coefficients are time ratios:

- TR > 1 → longer survival time
- TR < 1 → shorter survival time

Example:
- TR = 1.40 for treatment → treated survive 40% longer (on average, depending on model)

---

## 12A. Predicted survival curves for profiles (Python)

!!! interactive "Python"
    ```python
    import pandas as pd
    import matplotlib.pyplot as plt

    p_control = pd.DataFrame({"age":[60], "treatment":[0]})
    p_treated = pd.DataFrame({"age":[60], "treatment":[1]})

    s0 = aft.predict_survival_function(p_control)
    s1 = aft.predict_survival_function(p_treated)

    plt.plot(s0, label="Control (age=60)")
    plt.plot(s1, label="Treatment (age=60)")
    plt.title("Predicted Survival (Weibull AFT)")
    plt.xlabel("Time")
    plt.ylabel("S(t)")
    plt.legend()
    plt.show()
    ```

---

# PART B — R (survival + flexsurv)

In R, parametric survival models are commonly fit using:

- `survreg()` (AFT models)
- `flexsurv` package (highly recommended for multiple distributions + AIC comparison)

We provide both.

---

## 13B. Simulate data in R (same structure)

!!! interactive "R"
    ```r
    set.seed(44)

    n <- 400
    age <- rnorm(n, 60, 10)
    trt <- rbinom(n, 1, 0.5)

    shape_k <- 1.7
    scale_lam <- 12

    U <- runif(n)
    T0 <- scale_lam * (-log(U))^(1/shape_k)

    beta_trt <- 0.35
    beta_age <- -0.01

    T <- T0 * exp(beta_trt*trt + beta_age*(age-60))

    C <- runif(n, 3, 20)

    time <- pmin(T, C)
    event <- as.integer(T <= C)

    df <- data.frame(time=time, event=event, age=age, treatment=trt)
    head(df)
    ```

---

## 14B. KM curve in R (reference)

!!! interactive "R"
    ```r
    library(survival)

    km <- survfit(Surv(time, event) ~ 1, data=df)
    plot(km, xlab="Time", ylab="S(t)", main="KM Curve (Reference)")
    ```

---

## 15B. AFT regression models using survreg()

`survreg()` fits AFT models.

Common distributions:
- Weibull
- lognormal
- loglogistic

### 15B.1 Weibull AFT

!!! interactive "R"
    ```r
    fit_weib <- survreg(Surv(time, event) ~ age + treatment, data=df, dist="weibull")
    summary(fit_weib)
    ```

### Interpretation in survreg
`survreg()` uses a different parameterization:

- coefficients relate to log(time) directly
- sign interpretation can be reversed compared to hazard models

A simple practical interpretation:
- Positive coefficient → increases survival time
- Negative coefficient → decreases survival time

To compute time ratio:
\[
TR=\exp(\beta)
\]

!!! interactive "R"
    ```r
    exp(coef(fit_weib))
    ```

---

### 15B.2 Log-normal AFT

!!! interactive "R"
    ```r
    fit_logn <- survreg(Surv(time, event) ~ age + treatment, data=df, dist="lognormal")
    summary(fit_logn)
    exp(coef(fit_logn))
    ```

---

### 15B.3 Log-logistic AFT

!!! interactive "R"
    ```r
    fit_loglog <- survreg(Surv(time, event) ~ age + treatment, data=df, dist="loglogistic")
    summary(fit_loglog)
    exp(coef(fit_loglog))
    ```

---

## 16B. Model comparison using AIC (R)

!!! interactive "R"
    ```r
    AIC(fit_weib, fit_logn, fit_loglog)
    ```

Lower AIC = better tradeoff between fit and complexity.

---

## 17B. Using flexsurv for many distributions

`flexsurvreg()` gives access to:
- exponential
- Weibull
- lognormal
- loglogistic
- Gompertz, etc.

!!! interactive "R"
    ```r
    # install.packages("flexsurv") # if needed
    library(flexsurv)

    f_exp <- flexsurvreg(Surv(time, event) ~ age + treatment, data=df, dist="exp")
    f_weib <- flexsurvreg(Surv(time, event) ~ age + treatment, data=df, dist="weibull")
    f_logn <- flexsurvreg(Surv(time, event) ~ age + treatment, data=df, dist="lognormal")
    f_llog <- flexsurvreg(Surv(time, event) ~ age + treatment, data=df, dist="llogis")

    AIC(f_exp, f_weib, f_logn, f_llog)
    ```

---

## 18B. Plot fitted parametric curves vs KM (R)

With `flexsurv`, easiest is to overlay predicted survival.

!!! interactive "R"
    ```r
    plot(km, xlab="Time", ylab="S(t)", main="KM vs Parametric Fits", col="black")

    # Overlay parametric fit for "average" covariate profile:
    nd <- data.frame(age=60, treatment=0)

    lines(f_exp, newdata=nd, type="survival", col="blue")
    lines(f_weib, newdata=nd, type="survival", col="red")
    lines(f_logn, newdata=nd, type="survival", col="green")
    lines(f_llog, newdata=nd, type="survival", col="purple")

    legend("topright",
           legend=c("KM","Exp","Weibull","LogNormal","LogLogistic"),
           col=c("black","blue","red","green","purple"), lty=1)
    ```

---

## 19. How to interpret outputs (clinically)

### 19.1 PH (hazard ratio) interpretation
If using PH form (e.g., Weibull PH):

- HR < 1 → lower hazard
- HR > 1 → higher hazard

### 19.2 AFT interpretation (time ratio)
If using AFT form:

- TR = 1.40 → 40% longer survival time
- TR = 0.80 → 20% shorter survival time

AFT is often easier to interpret for clinicians.

---

## 20. When to use parametric models 

Use parametric models when you need:

- extrapolation beyond follow-up  
- smooth survival curve + clear formula  
- simulation / risk prediction  
- cost-effectiveness modeling  

Avoid if:
- hazard shape unknown and complex
- strong non-PH patterns and poor fit

---

## 21. Reporting in biostat papers

Example parametric report:

> “A Weibull AFT model adjusting for age showed that treatment prolonged survival (time ratio 1.42, 95% CI …). Model fit was assessed by AIC and visual comparison to Kaplan–Meier curves.”

Or PH form:

> “A Weibull proportional hazards model estimated treatment HR 0.68 (95% CI …).”

Always specify:
- distribution used
- model form (PH or AFT)
- interpretation (HR or TR)
- fit assessment method (AIC, plots)

---

## 22. Common mistakes

### Mistake 1: assuming exponential without checking
Exponential assumes constant hazard—often wrong.

### Mistake 2: choosing distribution only by AIC
Also check:
- KM overlay fit
- clinical plausibility of hazard shape

### Mistake 3: interpreting AFT coefficient as hazard ratio
AFT outputs time ratios, not HR.

### Mistake 4: extrapolating too far
Even parametric extrapolation can be unrealistic.
Always justify extrapolation window.

---

## 23. Key takeaways

- Parametric survival models assume a distribution for \(T\).
- Exponential = constant hazard.
- Weibull = flexible hazard (increasing/decreasing).
- Log-normal/log-logistic = non-monotonic hazards.
- Model selection via AIC + visual KM overlay.
- Regression can be PH or AFT; interpret HR vs TR correctly.

---

## 24. Exercises

<details>
<summary>Click to try</summary>

1. Fit exponential and Weibull models and compare AIC. Which fits better?  
2. Overlay fitted curves with KM and discuss which curve matches best.  
3. Fit Weibull AFT regression and interpret treatment as a time ratio.  
4. In R, fit log-normal and log-logistic and compare AIC.  
5. Simulate data with decreasing hazard (\(k<1\)) and see which parametric model fits best.

</details>
