# 01 — Clinical Trials Foundations & Phases 

Clinical trials are the gold standard for answering causal questions in medicine and public health.  
Unlike observational studies, trials are designed to **control bias** (especially confounding) using **randomization**, structured follow-up, and pre-specified analysis plans.

This chapter builds the foundation you’ll need for the rest of the Clinical Trials module.

---

## 1. What is a clinical trial?

A **clinical trial** is a prospective study in which participants are assigned to one or more interventions to evaluate effects on health outcomes.

Core features:
- **Prospective** follow-up (time flows forward)
- **Intervention** is assigned by the study (not chosen by the patient/doctor)
- Outcomes are measured in a structured way

In biostatistics, trials are mainly about:
- estimating treatment effects
- quantifying uncertainty (confidence intervals)
- preventing biased inference

---

## 2. Why trials are special (vs observational studies)

### 2.1 The confounding problem
In observational studies, treatment is not assigned randomly. People who receive treatment can differ systematically.

Example:
- sicker patients might be more likely to receive a new drug
- patients with more resources might access better care

These differences cause **confounding**, making treatment comparisons biased.

### 2.2 Randomization breaks confounding (in expectation)
Randomization ensures that, on average, treatment groups are comparable at baseline:
- measured covariates: age, sex, BMI, severity
- unmeasured covariates: unknown risk factors

So differences in outcome can be attributed to treatment (under assumptions).

---

## 3. Key trial terminology

### 3.1 Intervention / treatment
What is being tested:
- drug vs placebo
- new therapy vs standard care
- dose A vs dose B

### 3.2 Control group
The comparison group could be:
- placebo
- standard-of-care
- active comparator

### 3.3 Outcome and endpoint
- **Outcome**: what is measured (e.g., blood pressure)
- **Endpoint**: the specific trial-defined quantity used for analysis (e.g., change in SBP at 12 weeks)

### 3.4 Follow-up time
Most trials have:
- start time (randomization/baseline)
- fixed follow-up schedule
- final assessment window

### 3.5 Censoring (time-to-event trials)
If the endpoint is time-to-event:
- participant may not experience the event during the study
- we handle it using survival analysis methods

---

## 4. Types of clinical trials (high-level)

### 4.1 Superiority trials
Goal: show new treatment is better than control.

\[
H_0: \Delta \le 0 \quad \text{vs} \quad H_1: \Delta > 0
\]

### 4.2 Non-inferiority trials
Goal: show new treatment is not worse than control beyond a margin \(\delta\).

\[
H_0: \Delta \le -\delta \quad \text{vs} \quad H_1: \Delta > -\delta
\]

Common when:
- new treatment is cheaper, safer, or easier to deliver

### 4.3 Equivalence trials
Goal: show treatments are similar within a margin.

\[
-\delta < \Delta < \delta
\]

---

## 5. Trial phases (I–IV)

Clinical trials are often discussed in phases (drug development context).

### Phase I — Safety & dose (first-in-human)
- small sample (20–80)
- focus: safety, tolerability, dose-finding
- may include healthy volunteers or patients (oncology often patients)

### Phase II — Preliminary efficacy
- moderate sample (100–300)
- focus: efficacy signals + continued safety
- refine dose, endpoints, feasibility

### Phase III — Confirmatory trials
- large sample (hundreds to thousands)
- focus: definitive evidence of efficacy + safety
- designed for regulatory approval decisions

### Phase IV — Post-marketing surveillance
- after approval
- rare adverse events, long-term outcomes, real-world effectiveness

---

## 6. Efficacy vs Effectiveness 

### Efficacy
- “Does it work under ideal conditions?”
- often in controlled settings, strict inclusion/exclusion criteria

### Effectiveness
- “Does it work in real-world practice?”
- pragmatic trials, broader populations, messy conditions

Biostat tip:
- pragmatic designs often need careful handling of adherence and missingness.

---

## 7. Trial protocol and statistical analysis plan (SAP)

Every trial should have:
- **protocol** (clinical + operational plan)
- **SAP** (exact statistical methods before seeing data)

This reduces:
- p-hacking
- outcome switching
- selective reporting

Minimum SAP items:
- primary endpoint definition
- analysis population (ITT vs PP)
- model and covariates
- missing data strategy
- multiplicity strategy (if relevant)

---

## 8. A simple trial example (continuous endpoint)

Imagine a parallel-group trial:
- treatment vs placebo
- outcome: systolic blood pressure at 12 weeks
- compare mean change from baseline

### Key estimand (plain language)
> “Average difference in BP at 12 weeks between treatment and placebo, among all randomized patients.”

This is typically an **ITT-type** estimand.

---

## 9. Interactive simulation (Python): why randomization matters

We will simulate a scenario where:
- treatment is associated with severity in observational allocation
- randomization fixes it

!!! interactive "Python"
    ```python
    import numpy as np
    import pandas as pd

    np.random.seed(42)

    n = 600

    # Baseline severity (higher = worse)
    severity = np.random.normal(0, 1, n)

    # True treatment effect (reduces outcome)
    true_effect = -0.6

    # Outcome depends on severity + treatment + noise
    noise = np.random.normal(0, 1, n)

    # -------- Observational allocation (confounded) --------
    # People with higher severity more likely to receive treatment
    p_treat_obs = 1 / (1 + np.exp(-1.2 * severity))
    treat_obs = np.random.binomial(1, p_treat_obs)

    outcome_obs = 1.5 * severity + true_effect * treat_obs + noise

    df_obs = pd.DataFrame({"severity": severity, "treat": treat_obs, "outcome": outcome_obs})

    # naive difference in means
    naive_diff_obs = df_obs[df_obs.treat==1].outcome.mean() - df_obs[df_obs.treat==0].outcome.mean()
    naive_diff_obs
    ```

Interpretation:
- even though treatment truly helps, confounding can distort the estimate.

Now compare to randomized allocation:

!!! interactive "Python"
    ```python
    # -------- Randomized allocation --------
    treat_rct = np.random.binomial(1, 0.5, n)
    outcome_rct = 1.5 * severity + true_effect * treat_rct + noise

    df_rct = pd.DataFrame({"severity": severity, "treat": treat_rct, "outcome": outcome_rct})

    naive_diff_rct = df_rct[df_rct.treat==1].outcome.mean() - df_rct[df_rct.treat==0].outcome.mean()
    naive_diff_rct
    ```

You should see the randomized estimate is much closer to the true effect.

---

## 10. Interactive simulation (R): why randomization matters

!!! interactive "R"
    ```r
    set.seed(42)

    n <- 600
    severity <- rnorm(n, 0, 1)

    true_effect <- -0.6
    noise <- rnorm(n, 0, 1)

    # Observational allocation (confounded)
    p_treat_obs <- 1 / (1 + exp(-1.2 * severity))
    treat_obs <- rbinom(n, 1, p_treat_obs)

    outcome_obs <- 1.5 * severity + true_effect * treat_obs + noise
    df_obs <- data.frame(severity=severity, treat=treat_obs, outcome=outcome_obs)

    naive_diff_obs <- mean(df_obs$outcome[df_obs$treat==1]) - mean(df_obs$outcome[df_obs$treat==0])
    naive_diff_obs
    ```

Randomized allocation:

!!! interactive "R"
    ```r
    treat_rct <- rbinom(n, 1, 0.5)
    outcome_rct <- 1.5 * severity + true_effect * treat_rct + noise
    df_rct <- data.frame(severity=severity, treat=treat_rct, outcome=outcome_rct)

    naive_diff_rct <- mean(df_rct$outcome[df_rct$treat==1]) - mean(df_rct$outcome[df_rct$treat==0])
    naive_diff_rct
    ```

---

## 11. Exercises

<details>
<summary>Click to try</summary>

1. Modify the confounding strength (change 1.2 to 2.0 in the logistic allocation) and see how the naive observational estimate changes.  
2. Change the true treatment effect from -0.6 to -0.2 and repeat.  
3. Increase sample size to 3000 and check how the randomized estimate stabilizes.  
4. Plot severity distributions by group for observational vs randomized settings.  
5. Think: Why does randomization work “in expectation” but not necessarily perfectly in small samples?

</details>

---

## 12. Summary

- Clinical trials are designed for **causal inference**.
- Randomization reduces confounding and improves validity.
- Trial phases (I–IV) reflect safety → efficacy → confirmation → surveillance.
- Trials must define endpoints, populations, and analysis strategy **before** seeing results.
- Simulations show how observational allocation can bias treatment comparisons.

---


