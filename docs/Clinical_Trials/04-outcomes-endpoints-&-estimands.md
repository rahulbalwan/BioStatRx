# 04 — Outcomes, Endpoints, and Estimands (What Exactly Are We Estimating?)

Clinical trials are not just about collecting data — they are about answering a very specific question *clearly and unambiguously*.

A common reason trials become confusing (or results become hard to interpret) is that the question is not precisely defined.

This chapter teaches you how to define:

- **Outcome** (what is measured)
- **Endpoint** (what is analyzed)
- **Estimand** (what treatment effect you want to estimate)
- How to handle **intercurrent events** (treatment switching, rescue meds, death, dropout)
- How endpoint choice affects analysis, interpretation, and even trial success

This is a practical biostatistics chapter: you will learn how to translate a clinical question into an analysis-ready estimand.

---

## 1. Outcome vs endpoint vs estimand (key definitions)

### 1.1 Outcome
An outcome is any measurement collected during the trial.

Examples:
- systolic blood pressure
- HbA1c
- tumor size
- time-to-death
- hospital admissions

Outcomes may be collected at multiple times (baseline, week 4, week 12, etc.).

### 1.2 Endpoint
The endpoint is the trial-defined quantity derived from outcomes that will be analyzed.

Examples:
- change in systolic BP from baseline to week 12
- proportion achieving remission by week 24
- time from randomization to all-cause mortality

A trial can collect many outcomes but must specify:
- one **primary endpoint**
- possibly multiple **secondary endpoints**

### 1.3 Estimand (what effect are we trying to estimate?)
An estimand defines the treatment effect in a way that is robust to complications during follow-up.

An estimand specifies:
1) Population  
2) Treatment conditions (interventions)  
3) Endpoint variable (how outcome is defined)  
4) How intercurrent events are handled  
5) Summary measure (difference in means, HR, RR, etc.)

In plain language:
> The estimand is the precise question the trial is answering.

---

## 2. Primary vs secondary endpoints

### 2.1 Primary endpoint
- drives trial success/failure
- determines sample size
- must be pre-specified

Desirable properties:
- clinically meaningful
- measurable with low missingness
- sensitive to treatment effect
- not easily biased

### 2.2 Secondary endpoints
- supportive evidence
- mechanistic outcomes
- safety-related outcomes
- exploratory biomarkers

Statistical warning:
- multiple endpoints increase false positives unless multiplicity is addressed

---

## 3. Common endpoint types in biostatistics

### 3.1 Continuous endpoints
Examples:
- blood pressure (mmHg)
- cholesterol levels
- depression score
- lung function (FEV1)

Typical analysis:
- t-test / ANCOVA / linear regression
- mixed models for repeated measures (MMRM)

### 3.2 Binary endpoints
Examples:
- remission yes/no
- event occurred yes/no
- adverse event yes/no

Typical analysis:
- risk difference, risk ratio, odds ratio
- logistic regression or log-binomial / Poisson with robust SE

### 3.3 Count endpoints
Examples:
- number of hospital visits
- infections per person-year

Typical analysis:
- Poisson or negative binomial regression (with offset if rate)

### 3.4 Time-to-event endpoints
Examples:
- time to death
- time to relapse
- time to disease progression

Typical analysis:
- Kaplan–Meier / log-rank
- Cox regression
- competing risks methods when appropriate

---

## 4. Composite endpoints

Composite endpoints combine multiple events into one endpoint.

Example:
- major adverse cardiovascular events (MACE):
  death OR MI OR stroke

### 4.1 Why use composites?
- increases event rate (more power)
- captures multiple clinically relevant outcomes

### 4.2 Pitfalls
- components may differ in importance
- treatment may affect components differently
- interpretation becomes ambiguous

Best practice:
- report effect on the composite AND on each component.

---

## 5. Intercurrent events (the main reason estimands matter)

Intercurrent events happen after randomization and affect:
- interpretation of outcomes
- existence of outcome
- whether outcome reflects the intended treatment effect

Examples:
- treatment discontinuation
- rescue medication use
- switching to another therapy
- death before measurement
- withdrawal / loss to follow-up

Traditional trials sometimes handled these inconsistently.
The estimand framework forces clarity.

---

## 6. Estimand strategies (ICH E9(R1) in plain language)

We describe the most common strategies:

### 6.1 Treatment policy strategy
Ignore intercurrent events in the definition — analyze outcomes “as observed” regardless of switching/discontinuation.

Interpretation:
> Effect of assigning the treatment, regardless of adherence.

This often corresponds to ITT.

### 6.2 Hypothetical strategy
Ask:
> What would the outcome have been if the intercurrent event had not occurred?

Example:
- What would HbA1c be if no one used rescue medication?

This requires modeling assumptions or imputation.

### 6.3 Composite strategy
Make the intercurrent event part of the endpoint.

Example:
- define “failure” as relapse OR rescue medication use

Interpretation:
> Treatment effect on a broader clinical failure concept.

### 6.4 While-on-treatment strategy
Only consider outcomes while participants are on treatment.

Interpretation:
> Effect of treatment during actual exposure.

Danger:
- can introduce bias if discontinuation relates to prognosis

### 6.5 Principal stratum strategy (advanced)
Estimate effect within a latent subgroup defined by intercurrent event behavior.

Example:
- effect among those who would adhere under either treatment

Powerful but complex and assumptions are strong.

---

## 7. Example scenario: Diabetes trial (HbA1c at 24 weeks)

Trial:
- New drug vs standard drug
- Primary outcome: HbA1c at week 24
Complication:
- some participants use rescue medication if glucose high

Different estimands:

1) Treatment policy:
- compare HbA1c at week 24 regardless of rescue use

2) Hypothetical:
- estimate HbA1c if rescue medication had not been used

3) Composite:
- define failure as HbA1c above threshold OR rescue use

Each answers a different clinical question.

---

# Part A — Practical simulations (Python)

We simulate a trial with rescue medication to show how estimand choice changes effect estimates.

---

## 8A. Python simulation: rescue medication creates complexity

!!! interactive "Python"
    ```python
    import numpy as np
    import pandas as pd

    np.random.seed(42)

    n = 600
    trt = np.random.binomial(1, 0.5, n)  # 1 = new drug, 0 = control

    # baseline HbA1c
    hba1c0 = np.random.normal(8.5, 0.8, n)

    # true treatment effect on HbA1c reduction
    true_delta = -0.6
    noise = np.random.normal(0, 0.6, n)

    # Week 24 HbA1c if fully adherent and no rescue (latent)
    hba1c24_latent = hba1c0 + true_delta*trt + noise

    # Rescue medication used if latent HbA1c remains too high
    rescue = (hba1c24_latent > 8.2).astype(int)

    # Rescue medication improves HbA1c by extra reduction (but only among those who take it)
    rescue_effect = -0.7
    hba1c24_observed = hba1c24_latent + rescue_effect*rescue

    df = pd.DataFrame({
        "trt": trt,
        "hba1c0": hba1c0,
        "hba1c24_latent": hba1c24_latent,
        "rescue": rescue,
        "hba1c24_observed": hba1c24_observed
    })

    df.head()
    ```

---

## 9A. Python: Treatment policy estimand (use observed outcomes)

!!! interactive "Python"
    ```python
    # Treatment policy: compare observed HbA1c at 24 weeks
    mean_trt = df.loc[df.trt==1, "hba1c24_observed"].mean()
    mean_ctl = df.loc[df.trt==0, "hba1c24_observed"].mean()
    mean_trt - mean_ctl
    ```

---

## 10A. Python: Hypothetical estimand (no rescue) using latent outcome

In real life, latent is unobserved; here we use it to show the concept.

!!! interactive "Python"
    ```python
    # Hypothetical: what if no rescue were used?
    mean_trt = df.loc[df.trt==1, "hba1c24_latent"].mean()
    mean_ctl = df.loc[df.trt==0, "hba1c24_latent"].mean()
    mean_trt - mean_ctl
    ```

---

## 11A. Python: Composite estimand (failure = rescue OR high HbA1c)

Define failure as:
- rescue used OR observed HbA1c > 8.2

!!! interactive "Python"
    ```python
    fail = ((df.rescue==1) | (df.hba1c24_observed > 8.2)).astype(int)
    rate_trt = fail[df.trt==1].mean()
    rate_ctl = fail[df.trt==0].mean()
    rate_trt, rate_ctl, rate_trt - rate_ctl
    ```

Interpretation:
- this answers a different question: risk of clinical “failure”.

---

# Part B — Practical simulations (R)

---

## 12B. R simulation: rescue medication

!!! interactive "R"
    ```r
    set.seed(42)

    n <- 600
    trt <- rbinom(n, 1, 0.5)

    hba1c0 <- rnorm(n, 8.5, 0.8)
    true_delta <- -0.6
    noise <- rnorm(n, 0, 0.6)

    hba1c24_latent <- hba1c0 + true_delta*trt + noise

    rescue <- as.integer(hba1c24_latent > 8.2)

    rescue_effect <- -0.7
    hba1c24_observed <- hba1c24_latent + rescue_effect*rescue

    df <- data.frame(
      trt=trt,
      hba1c0=hba1c0,
      hba1c24_latent=hba1c24_latent,
      rescue=rescue,
      hba1c24_observed=hba1c24_observed
    )

    head(df)
    ```

---

## 13B. R: Treatment policy estimand (observed)

!!! interactive "R"
    ```r
    mean(df$hba1c24_observed[df$trt==1]) - mean(df$hba1c24_observed[df$trt==0])
    ```

---

## 14B. R: Hypothetical estimand (no rescue) using latent outcome

!!! interactive "R"
    ```r
    mean(df$hba1c24_latent[df$trt==1]) - mean(df$hba1c24_latent[df$trt==0])
    ```

---

## 15B. R: Composite estimand (failure)

!!! interactive "R"
    ```r
    fail <- as.integer(df$rescue==1 | df$hba1c24_observed > 8.2)

    rate_trt <- mean(fail[df$trt==1])
    rate_ctl <- mean(fail[df$trt==0])

    rate_trt
    rate_ctl
    rate_trt - rate_ctl
    ```

---

## 16. How to write endpoints/estimands in a protocol (template)

Example template:

**Primary endpoint:**  
HbA1c measured at 24 weeks post-randomization.

**Intercurrent event handling:**  
Use of rescue medication will be handled using a treatment policy strategy (primary analysis). Sensitivity analyses using a hypothetical strategy will be conducted via multiple imputation assuming no rescue medication use.

**Population:**  
All randomized participants (ITT).

**Summary measure:**  
Difference in mean HbA1c between groups with 95% confidence interval.

---

## 17. Exercises

<details>
<summary>Click to try</summary>

1. Change the rescue threshold from 8.2 to 8.6 and see how the treatment policy estimate changes.  
2. Increase the rescue effect magnitude (e.g., -1.2) and compare estimands.  
3. Redefine the composite endpoint using a stricter cutoff (HbA1c > 7.5).  
4. Explain (in words) what clinical question each estimand answers.  
5. In a time-to-event setting, propose intercurrent events and how you would handle them.

</details>

---

## 18. Summary

- Outcomes are measurements; endpoints are analysis targets; estimands define the treatment effect question.
- Intercurrent events (switching, rescue, dropout) can change interpretation.
- Different estimand strategies answer different clinical questions.
- A well-written trial makes endpoints and estimands explicit before data is seen.

---
