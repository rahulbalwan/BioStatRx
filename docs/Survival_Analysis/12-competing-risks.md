# Competing Risks (CIF, Cause-Specific Hazards, Fine–Gray) 

In many medical studies, subjects can fail from **different causes**.

Example (oncology):
- Event of interest: relapse
- Competing event: death before relapse

If a patient dies, they can no longer relapse.

So death **competes** with relapse.

This changes the meaning of “risk” and breaks naive Kaplan–Meier.

This chapter gives maximum clarity on:

- why KM is wrong under competing risks  
- cumulative incidence function (CIF)  
- cause-specific hazard models  
- Fine–Gray subdistribution hazard model  
- interpretation and reporting  
- Python + R code with plots  

---

## 1. What is a competing risk?

A competing risk is an event that:

1) prevents the event of interest from happening, and  
2) is not just censoring.

### 1.1 Example (classic)
Event of interest: relapse  
Competing event: death

If death happens, relapse can never occur afterwards.

So death is not “censoring” — it is a different event that removes the possibility of relapse.

---

## 2. Why Kaplan–Meier is wrong for competing risks

KM assumes that censoring is non-informative:

> censored subjects are like those still at risk

But in competing risks:

- a competing event is not “like censoring”
- it changes the probability structure

If you treat competing events as censored:

- KM **overestimates** probability of event of interest.

Because KM assumes those who died could still relapse later, which is impossible.

---

## 3. The correct probability: cumulative incidence function (CIF)

### 3.1 Definition
For event type \(k\):

\[
F_k(t) = P(T \le t,\ \text{event}=k)
\]

Interpretation:

> Probability of experiencing event type \(k\) by time \(t\), accounting for competing events.

This is the correct function to report when competing risks exist.

### 3.2 Relationship
If there are K event types:

\[
S(t) + \sum_{k=1}^K F_k(t) = 1
\]

Where:
- \(S(t)\) = probability of no event yet
- \(F_k(t)\) = probability of event k by time t

---

## 4. Two ways to model competing risks

### 4.1 Cause-specific hazards (CSH)
Model hazard of event k while treating other causes as censored.

This answers:

> “Among those currently event-free, what is the instantaneous hazard of event k?”

Good for etiologic questions and mechanisms.

### 4.2 Fine–Gray subdistribution hazard (SDH)
Models the hazard of the CIF directly.

This answers:

> “How do covariates affect the cumulative incidence probability of event k over time?”

Good for prognosis and absolute risk prediction.

---

## 5. Cause-specific hazard model (CSH)

For event type k:

\[
h_k(t|X)=h_{0k}(t)\exp(\beta_k^TX)
\]

Interpretation of HR:

> HR compares instantaneous hazard of cause k among those still event-free.

Important:
- This HR is NOT directly a HR on CIF probability.
- HR on cause-specific hazard and Fine–Gray can differ.

---

## 6. Fine–Gray model (subdistribution hazard)

Fine–Gray models subdistribution hazard:

\[
\tilde h_k(t|X)=\tilde h_{0k}(t)\exp(\gamma^TX)
\]

Interpretation:

> HR describes relative rate of accumulating event-k incidence over time, accounting for competing events.

This is often interpreted as:

- covariate increases/decreases cumulative incidence probability.

---

## 7. When to use which model?

### Use cause-specific Cox when:
- scientific question is about biological mechanism
- “does treatment reduce relapse hazard among those alive?”
- you want hazard-based interpretation

### Use Fine–Gray when:
- scientific question is about absolute risk prediction
- “what is probability of relapse by 5 years accounting for death?”
- clinical decision-making and prognosis

Many papers report BOTH.

---

# PART A — PYTHON

Python has limited built-in competing risks support compared to R.
But we can still do CIF estimation and cause-specific Cox in Python.

For Fine–Gray regression, R is much stronger (cmprsk, riskRegression).

We will still provide best-possible Python workflow.

---

## 8A. Simulate competing risks data (Python)

We simulate two event types:

- 1 = relapse (event of interest)
- 2 = death (competing event)

!!! interactive "Python"
    ```python
    import numpy as np
    import pandas as pd

    np.random.seed(1234)

    n = 500
    age = np.random.normal(60, 10, n)
    trt = np.random.binomial(1, 0.5, n)

    # Cause-specific hazards depend on covariates
    # Relapse hazard: treatment reduces, age increases
    base_relapse = 0.04
    base_death   = 0.03

    beta_age_relapse = 0.02
    beta_trt_relapse = -0.50

    beta_age_death = 0.04
    beta_trt_death = 0.00  # assume treatment does not affect death

    h_relapse = base_relapse * np.exp(beta_age_relapse*age + beta_trt_relapse*trt)
    h_death   = base_death   * np.exp(beta_age_death*age   + beta_trt_death*trt)

    # Simulate event times for each cause (exponential here for simplicity)
    T1 = np.random.exponential(1/h_relapse)  # relapse time
    T2 = np.random.exponential(1/h_death)    # death time

    # Observed event = whichever happens first
    T = np.minimum(T1, T2)
    event_type = np.where(T1 < T2, 1, 2)

    # Add censoring
    C = np.random.uniform(2, 25, n)
    time = np.minimum(T, C)
    status = np.where(T <= C, event_type, 0)  # 0=censored, 1=relapse, 2=death

    df = pd.DataFrame({"time": time, "status": status, "age": age, "treatment": trt})
    df.head()
    ```

Status meanings:
- 0 = censored
- 1 = relapse
- 2 = death (competing)

---

## 9A. Estimate CIF in Python (nonparametric)

We’ll estimate CIF using a simple Aalen–Johansen style approach.

This code computes CIF by:
- estimating overall survival
- estimating cause-specific increments

!!! interactive "Python"
    ```python
    import numpy as np
    import pandas as pd

    def cif_aalen_johansen(time, status, cause=1):
        # time: observed times
        # status: 0=censor, 1..K causes
        df = pd.DataFrame({"time": time, "status": status}).sort_values("time")

        unique_times = np.sort(df.loc[df.status > 0, "time"].unique())

        S = 1.0
        cif = 0.0
        out = []

        for t in unique_times:
            at_risk = (df["time"] >= t).sum()
            d_all   = ((df["time"] == t) & (df["status"] > 0)).sum()
            d_cause = ((df["time"] == t) & (df["status"] == cause)).sum()

            # increment in CIF for this cause
            if at_risk > 0:
                cif += S * (d_cause / at_risk)
                # update survival after all events at time t
                S *= (1 - d_all / at_risk)

            out.append((t, cif))

        return pd.DataFrame(out, columns=["time", f"CIF_cause{cause}"])

    cif1 = cif_aalen_johansen(df["time"], df["status"], cause=1)
    cif2 = cif_aalen_johansen(df["time"], df["status"], cause=2)

    cif1.head(), cif2.head()
    ```

---

## 10A. Plot CIF curves in Python

!!! interactive "Python"
    ```python
    import matplotlib.pyplot as plt

    plt.plot(cif1["time"], cif1["CIF_cause1"], label="Relapse CIF")
    plt.plot(cif2["time"], cif2["CIF_cause2"], label="Death CIF")

    plt.title("Cumulative Incidence Functions (Competing Risks)")
    plt.xlabel("Time")
    plt.ylabel("CIF")
    plt.legend()
    plt.show()
    ```

Interpretation:
- CIF(t) is probability of event by time t.
- These curves typically add up with survival probability.

---

## 11A. Cause-specific Cox models in Python (lifelines)

### 11A.1 Relapse model (treat death as censored)

Create event indicator:
- event=1 if relapse
- event=0 if death or censoring

!!! interactive "Python"
    ```python
    from lifelines import CoxPHFitter

    df_relapse = df.copy()
    df_relapse["event_relapse"] = (df_relapse["status"] == 1).astype(int)

    cph_relapse = CoxPHFitter()
    cph_relapse.fit(df_relapse[["time","event_relapse","age","treatment"]],
                    duration_col="time", event_col="event_relapse")

    cph_relapse.print_summary()
    ```

Interpretation:
- HR applies to relapse hazard among those event-free.

### 11A.2 Death model (treat relapse as censored)

!!! interactive "Python"
    ```python
    df_death = df.copy()
    df_death["event_death"] = (df_death["status"] == 2).astype(int)

    cph_death = CoxPHFitter()
    cph_death.fit(df_death[["time","event_death","age","treatment"]],
                  duration_col="time", event_col="event_death")

    cph_death.print_summary()
    ```

---

## 12A. Why Fine–Gray is harder in Python
Fine–Gray regression is not fully supported in standard lifelines workflows.
You can:
- estimate CIF nonparametrically 
- fit cause-specific Cox 
But for Fine–Gray regression, R is recommended.

---

# PART B — R (best tool for competing risks)

In R, competing risks is commonly handled with:

- `cmprsk` (Fine–Gray)  
- `survival` (cause-specific Cox)  
- `riskRegression` (CIF prediction / model comparison)  

---

## 13B. Simulate competing risks data in R

!!! interactive "R"
    ```r
    set.seed(1234)

    n <- 500
    age <- rnorm(n, 60, 10)
    trt <- rbinom(n, 1, 0.5)

    base_relapse <- 0.04
    base_death   <- 0.03

    beta_age_relapse <- 0.02
    beta_trt_relapse <- -0.50

    beta_age_death <- 0.04
    beta_trt_death <- 0.00

    h_relapse <- base_relapse * exp(beta_age_relapse*age + beta_trt_relapse*trt)
    h_death   <- base_death   * exp(beta_age_death*age   + beta_trt_death*trt)

    T1 <- rexp(n, rate=h_relapse)
    T2 <- rexp(n, rate=h_death)

    T <- pmin(T1, T2)
    status <- ifelse(T1 < T2, 1, 2)

    C <- runif(n, 2, 25)
    time <- pmin(T, C)
    status_obs <- ifelse(T <= C, status, 0)

    df <- data.frame(time=time, status=status_obs, age=age, treatment=trt)
    head(df)
    ```

---

## 14B. CIF estimation in R (cmprsk)

!!! interactive "R"
    ```r
    # install.packages("cmprsk") # if needed
    library(cmprsk)

    # status: 0=censored, 1=relapse, 2=death
    cif <- cuminc(ftime=df$time, fstatus=df$status)

    cif
    plot(cif, lty=1, col=c("red","blue"),
         xlab="Time", ylab="CIF",
         main="Cumulative Incidence Functions (Competing Risks)")
    legend("topleft", legend=c("Relapse","Death"), col=c("red","blue"), lty=1)
    ```

---

## 15B. CIF by treatment group 

!!! interactive "R"
    ```r
    cif_g <- cuminc(ftime=df$time, fstatus=df$status, group=df$treatment)
    cif_g
    plot(cif_g, lty=1, xlab="Time", ylab="CIF",
         main="CIF by Treatment Group")
    ```

Interpretation:
- This shows relapse probability over time accounting for death,
separately by treatment group.

---

## 16B. Cause-specific Cox models in R

### 16B.1 Relapse cause-specific Cox
Treat death as censored:

!!! interactive "R"
    ```r
    library(survival)

    event_relapse <- as.integer(df$status == 1)
    fit_relapse <- coxph(Surv(df$time, event_relapse) ~ age + treatment, data=df)
    summary(fit_relapse)
    ```

### 16B.2 Death cause-specific Cox
Treat relapse as censored:

!!! interactive "R"
    ```r
    event_death <- as.integer(df$status == 2)
    fit_death <- coxph(Surv(df$time, event_death) ~ age + treatment, data=df)
    summary(fit_death)
    ```

---

## 17B. Fine–Gray model in R (subdistribution hazard)

Fine–Gray models the CIF directly.

Use `crr()` from cmprsk:

!!! interactive "R"
    ```r
    # Fine–Gray for relapse (cause=1)
    fg <- crr(ftime=df$time, fstatus=df$status, cov1=df[,c("age","treatment")], failcode=1)
    fg
    ```

Interpretation:
- exponentiated coefficients are subdistribution HRs
- relate to CIF probability over time

---

## 18. Interpretation: cause-specific HR vs Fine–Gray HR

### Cause-specific HR answers:
> "Among those still event-free, how does X affect instantaneous hazard of relapse?"

### Fine–Gray HR answers:
> "How does X affect cumulative incidence probability of relapse over time, accounting for death?"

They can differ because:
- Fine–Gray keeps competing events “in the risk set” in a special way.

---

## 19. Practical reporting examples

### CIF reporting
> “The 5-year cumulative incidence of relapse was 0.24 (accounting for competing risk of death).”

### Cause-specific Cox reporting
> “In cause-specific Cox regression, treatment reduced relapse hazard (HR 0.62, 95% CI …).”

### Fine–Gray reporting
> “In Fine–Gray regression, treatment reduced subdistribution hazard of relapse (sHR 0.70, 95% CI …).”

Always specify:
- which model
- which interpretation
- competing event definition

---

## 20. Common mistakes 

### Mistake 1: using KM to estimate relapse probability with death competing
KM overestimates relapse probability.

Fix: use CIF.

### Mistake 2: interpreting cause-specific HR as directly affecting CIF probability
Cause-specific HR affects hazard, not CIF directly.

Fix: use Fine–Gray if question is CIF probability.

### Mistake 3: calling competing risks “censoring”
Death is not censoring if relapse is endpoint.

Fix: treat as competing event.

---

## 21. Key takeaways

- Competing risks occur when another event prevents event of interest.
- KM is wrong for event probability under competing risks.
- CIF is the correct event probability estimate.
- Cause-specific Cox models hazard among those still event-free.
- Fine–Gray models the cumulative incidence (probability) directly.
- R is currently the best ecosystem for full competing risks regression.

---

## 22. Exercises

<details>
<summary>Click to try</summary>

1. Simulate competing risks data and compare KM relapse estimate (treat death as censoring) vs CIF relapse estimate. Which is higher and why?  
2. Fit cause-specific Cox for relapse and interpret HR clinically.  
3. Fit Fine–Gray for relapse and interpret subdistribution HR.  
4. Create a scenario where treatment reduces relapse but increases death. What happens to CIF of relapse?  
5. Plot CIF curves by treatment and report 5-year relapse incidence.

</details>
