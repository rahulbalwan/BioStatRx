# 10. Time-Dependent Covariates (Extended Cox Model) 

In many real clinical and epidemiological studies, covariates **change over time**.

Examples:
- treatment starts later or changes dose
- biomarkers change (CD4 count, blood pressure, viral load)
- smoking status changes
- a patient becomes infected or vaccinated during follow-up
- hospital discharge changes risk of in-hospital outcomes

Standard Cox assumes covariates are fixed:

\[
X_i \text{ is constant for the entire follow-up}
\]

When this is false, standard Cox can become biased—sometimes severely.

This chapter covers:

- what time-dependent covariates are  
- how to structure data correctly (start–stop / counting process)  
- how to fit models in **Python and R**  
- classic bias: **immortal time bias**  
- how to interpret results correctly  

---

## 1. What is a time-dependent covariate?

A time-dependent covariate is a covariate that depends on time:

\[
X_i(t)
\]

Cox becomes:

\[
h(t|X(t))=h_0(t)\exp(\beta^T X(t))
\]

The model is still Cox PH—same idea—but covariates can update.

---

## 2. Two concepts students confuse (important)

### 2.1 Time-dependent covariate
The covariate value changes over time.

Example:
- treatment: 0 before drug starts, 1 after drug starts

### 2.2 Time-varying effect (non-PH)
The effect of a covariate changes over time (PH violation).

Example:
- treatment HR = 0.6 early but 1.1 late

These are different:
- Time-dependent covariate is about **X(t)**
- Time-varying effect is about **β(t)**

---

## 3. Why time-dependent covariates matter (biostat motivation)

### 3.1 Treatment initiation happens after baseline
Example: observational study of statin use.

If a patient starts statin at month 6, and we code:

- `statin = 1` from time 0

Then we falsely give the drug credit for survival during months 0–6 even though they weren’t on it.

This creates a major bias.

---

## 4. Immortal time bias (the classic deadly mistake)

Immortal time bias occurs when:

> to be classified as “treated”, the subject must survive long enough to receive treatment.

So early time becomes “immortal” (event cannot occur before treatment starts if you require survival to receive it).

If you label them as treated from baseline, you inflate treatment benefit.

### 4.1 Example
- surgery at day 60
- patient must survive 60 days to have surgery
If you label surgery group at day 0:

- they appear to have better survival, partly because they had to survive long enough.

Correct fix: treat surgery status as time-dependent.

---

## 5. Correct data structure: start–stop (counting process)

Instead of one row per subject:

| id | time | event | treatment |

We use multiple rows:

| id | start | stop | event | treatment |

Each row represents a time interval where covariates are constant.

---

## 6. Worked example: delayed treatment

A patient:
- enters at time 0
- starts treatment at time 5
- experiences event at time 10

Correct representation:

| id | start | stop | event | trt |
|---:|------:|-----:|------:|----:|
| 1 | 0 | 5 | 0 | 0 |
| 1 | 5 | 10 | 1 | 1 |

Now the model knows:
- untreated from 0 to 5
- treated from 5 to 10

---

## 7. Practical steps to build time-dependent dataset

### Step 1
Identify covariates that change (treatment start, biomarker measurement times)

### Step 2
Split each subject’s follow-up into intervals

### Step 3
Create start–stop rows with updated covariates

### Step 4
Fit extended Cox using start–stop format

---

# PART A — PYTHON (lifelines)

Python survival modeling with time-dependent covariates is typically done using:

- `CoxTimeVaryingFitter` from lifelines

---

## 8A. Simulate a delayed treatment dataset (Python)

We simulate:
- each subject has a treatment start time
- before that they are untreated
- after that they are treated
- treatment reduces hazard

!!! interactive "Python"
    ```python
    import numpy as np
    import pandas as pd

    np.random.seed(777)

    n = 200

    # Treatment start time (some start late, some early)
    t_start = np.random.uniform(1, 8, size=n)

    # baseline covariate
    age = np.random.normal(60, 10, n)

    # true effects
    beta_age = 0.03
    beta_trt = -0.60

    base = 0.05

    rows = []

    for i in range(n):
        # hazard before treatment
        h0 = base * np.exp(beta_age*age[i] + beta_trt*0)

        # hazard after treatment
        h1 = base * np.exp(beta_age*age[i] + beta_trt*1)

        # simulate event time in two-stage piecewise exponential
        t1 = np.random.exponential(1/h0)

        if t1 < t_start[i]:
            # event occurs before treatment begins
            event_time = t1
            # one interval only: untreated
            rows.append([i, 0, event_time, 1, 0, age[i]])
        else:
            # survived untreated period
            # then simulate time after treatment
            t2 = np.random.exponential(1/h1)
            event_time = t_start[i] + t2

            # censoring
            censor = np.random.uniform(3, 15)

            if event_time <= censor:
                # two intervals: untreated then treated; event in second
                rows.append([i, 0, t_start[i], 0, 0, age[i]])
                rows.append([i, t_start[i], event_time, 1, 1, age[i]])
            else:
                # censored in second interval
                rows.append([i, 0, t_start[i], 0, 0, age[i]])
                rows.append([i, t_start[i], censor, 0, 1, age[i]])

    tv = pd.DataFrame(rows, columns=["id","start","stop","event","treatment","age"])
    tv.head(10)
    ```

This is start–stop format.

---

## 9A. Fit time-dependent Cox in Python

!!! interactive "Python"
    ```python
    from lifelines import CoxTimeVaryingFitter

    ctv = CoxTimeVaryingFitter()
    ctv.fit(tv, id_col="id", start_col="start", stop_col="stop", event_col="event")

    ctv.print_summary()
    ```

Interpretation:
- `exp(coef)` for treatment is the hazard ratio comparing treated vs untreated at the same time.

Expected:
- HR < 1 (protective)

---

## 10A. What happens if we do it WRONG (baseline treatment only)

Wrong approach:
- label treated if they ever got treatment
- then fit standard Cox

This creates immortal time bias.

!!! interactive "Python"
    ```python
    # Collapse to one row per person
    baseline = tv.groupby("id").agg(
        time=("stop","max"),
        event=("event","max"),
        age=("age","first"),
        ever_treated=("treatment","max")
    ).reset_index()

    from lifelines import CoxPHFitter
    cph_wrong = CoxPHFitter()
    cph_wrong.fit(baseline[["time","event","age","ever_treated"]], duration_col="time", event_col="event")

    cph_wrong.print_summary()
    ```

Compare HR from:
- time-dependent model (correct)
- baseline ever-treated Cox (wrong)

You will usually see the wrong model make treatment look *more protective* than it truly is.

---

# PART B — R (survival)

In R, time-dependent covariates are handled using:

- `coxph()` with `Surv(start, stop, event)` form  
(i.e., counting-process notation)

---

## 11B. Simulate delayed treatment dataset in R

!!! interactive "R"
    ```r
    set.seed(777)

    n <- 200

    t_start <- runif(n, 1, 8)
    age <- rnorm(n, 60, 10)

    beta_age <- 0.03
    beta_trt <- -0.60
    base <- 0.05

    rows <- list()
    idx <- 1

    for (i in 1:n) {

      h0 <- base * exp(beta_age*age[i] + beta_trt*0)
      h1 <- base * exp(beta_age*age[i] + beta_trt*1)

      t1 <- rexp(1, rate = h0)

      if (t1 < t_start[i]) {
        event_time <- t1
        rows[[idx]] <- data.frame(id=i, start=0, stop=event_time, event=1,
                                  treatment=0, age=age[i])
        idx <- idx + 1
      } else {

        t2 <- rexp(1, rate = h1)
        event_time <- t_start[i] + t2
        censor <- runif(1, 3, 15)

        if (event_time <= censor) {
          rows[[idx]] <- data.frame(id=i, start=0, stop=t_start[i], event=0,
                                    treatment=0, age=age[i]); idx <- idx + 1
          rows[[idx]] <- data.frame(id=i, start=t_start[i], stop=event_time, event=1,
                                    treatment=1, age=age[i]); idx <- idx + 1
        } else {
          rows[[idx]] <- data.frame(id=i, start=0, stop=t_start[i], event=0,
                                    treatment=0, age=age[i]); idx <- idx + 1
          rows[[idx]] <- data.frame(id=i, start=t_start[i], stop=censor, event=0,
                                    treatment=1, age=age[i]); idx <- idx + 1
        }
      }
    }

    tv <- do.call(rbind, rows)
    head(tv, 10)
    ```

---

## 12B. Fit time-dependent Cox in R

Use:

\[
Surv(start, stop, event)
\]

!!! interactive "R"
    ```r
    library(survival)

    fit_td <- coxph(Surv(start, stop, event) ~ age + treatment, data=tv)
    summary(fit_td)
    ```

Interpretation:
- `exp(coef)` for treatment is HR comparing treated vs untreated intervals.

---

## 13B. WRONG baseline ever-treated model in R

!!! interactive "R"
    ```r
    baseline <- aggregate(cbind(stop, event, treatment, age) ~ id, data=tv,
                          FUN=function(x) x[1])

    # baseline time is max stop per id
    time_max <- aggregate(stop ~ id, data=tv, FUN=max)

    baseline$time <- time_max$stop
    baseline$event <- aggregate(event ~ id, data=tv, FUN=max)$event
    baseline$ever_treated <- aggregate(treatment ~ id, data=tv, FUN=max)$treatment

    fit_wrong <- coxph(Surv(time, event) ~ age + ever_treated, data=baseline)
    summary(fit_wrong)
    ```

Compare:
- `fit_td` (correct)
- `fit_wrong` (biased)

---

## 14. Interpretation rules (must be clear)

### 14.1 Treatment HR in time-dependent Cox
HR compares:
- hazard when **currently treated**
vs
- hazard when **currently untreated**

It does NOT compare “ever treated” vs “never treated”.

### 14.2 Time-dependent biomarker (e.g., CD4 count)
If CD4 changes over time, you model:

\[
h(t)=h_0(t)\exp(\beta \cdot CD4(t))
\]

Interpretation:
- at time t, higher CD4 is associated with lower/higher hazard.

---

## 15. Common real biostat scenarios

### 15.1 Treatment switching
Patients may:
- stop drug
- start new drug
- change dose

Time-dependent Cox handles this by updating treatment status across intervals.

### 15.2 Longitudinal biomarkers
If measurements occur at repeated visits:
- split follow-up at visit times
- carry-forward biomarker value until next measurement

---

## 16. Common mistakes and best practices

### Mistake 1: coding “ever treated” as baseline
→ immortal time bias

Fix: use start–stop data and time-dependent treatment variable.

### Mistake 2: using event time itself to define covariate (reverse causality)
Example:
- covariate measured after event
This can create bias.

Fix: ensure covariate values are known BEFORE each interval.

### Mistake 3: too many splits without reason
Very fine splitting can become heavy computationally.

Fix: split only when covariates change or at visit times.

---

## 17. Key takeaways

- Time-dependent covariates are covariates \(X(t)\) that change over time.
- Use start–stop (counting-process) format.
- Extended Cox handles delayed treatment and longitudinal predictors correctly.
- Avoid immortal time bias by not using baseline “ever treated.”
- Interpretation is about **current covariate value at time t**.

---

## 18. Exercises (high-value practice)

<details>
<summary>Click to try</summary>

1. Simulate delayed treatment and compare HR from correct time-dependent Cox vs wrong baseline ever-treated Cox.  
2. Add a biomarker measured every 3 months; create start–stop intervals and fit model.  
3. Create treatment switching (0→1→0) and model it.  
4. In R, fit time-dependent Cox with `Surv(start, stop, event)`. Confirm matches Python results.  
5. Explain immortal time bias in one paragraph with a real clinical example.

</details>
