# 08 — Reporting & Interpretation 

A clinical trial result is only useful if it is reported clearly and reproducibly.

Good reporting is not just “writing well” — it is about:
- transparency in design and analysis
- correct interpretation of effect sizes and uncertainty
- communicating results so they can be trusted and reused

In this chapter, we cover:
- CONSORT reporting essentials
- effect sizes + confidence intervals (what to emphasize)
- common plots (Kaplan–Meier, forest plot, effect plots)
- subgroup analyses (and how to avoid misleading conclusions)
- reproducible tables/figures in R and Python

---

## 1. CONSORT: the standard reporting framework

CONSORT (Consolidated Standards of Reporting Trials) provides guidance for reporting randomized trials.

Key components include:
- trial design and allocation ratio
- eligibility criteria and settings
- interventions and follow-up
- outcomes and sample size calculation
- randomization method + concealment + blinding
- participant flow (flow diagram)
- baseline characteristics table
- results for primary and secondary endpoints
- harms (safety)
- interpretation and generalizability
- protocol registration and funding/conflicts

Even if you don’t memorize the checklist, your trial reporting should naturally answer:
- Who was included?
- What was compared?
- How was bias controlled?
- What was found, with uncertainty?
- How robust are the results?

---

## 2. Participant flow 

A typical flow includes:
- assessed for eligibility
- excluded (with reasons)
- randomized
- allocated to each group
- follow-up completed and lost
- analyzed (ITT, PP, safety)

Even if the diagram is produced in a word processor, you should compute these counts from data.

---

## 3. Baseline characteristics table (“Table 1”)

Baseline table is descriptive.

Recommended content:
- N randomized by arm
- demographics: age, sex, race/ethnicity (context-dependent)
- disease characteristics: baseline severity, stage, risk factors
- key prognostic factors and stratification variables

Important:
- do not use p-values to “test baseline balance”
- report summaries and (optionally) standardized differences

---

## 4. Primary endpoint reporting: focus on effect size + CI

The most important reporting elements are:
- estimate of treatment effect
- 95% confidence interval
- p-value (secondary, not primary)

Examples:
- continuous: difference in means (adjusted or unadjusted) with CI
- binary: risk difference / risk ratio / odds ratio with CI
- time-to-event: hazard ratio with CI + KM curves

Interpretation should answer:
- Is the effect clinically meaningful?
- How precise is the estimate?
- Are conclusions robust to assumptions?

---

## 5. Interpreting confidence intervals correctly

A 95% CI is not:
- “there is a 95% probability the true effect is inside” (frequentist CI interpretation is different)

A practical interpretation:
> If we repeated this trial many times, 95% of such CIs would contain the true effect.

For readers, CI gives:
- direction and magnitude
- plausible range of effects
- whether results include “no effect” (0 for differences, 1 for ratios)

---

## 6. Subgroup analyses: benefit and danger

Subgroups are tempting but risky.

### 6.1 Why subgroups can mislead
- many subgroups → many false positives
- subgroups reduce sample size → unstable estimates
- “significant in subgroup A but not in subgroup B” does not prove a difference

### 6.2 Correct approach: interaction tests
To test whether treatment differs by subgroup:

Model:
\[
Y = \beta_0 + \beta_1 T + \beta_2 S + \beta_3(T \times S) + \epsilon
\]

- \(\beta_3\) is the interaction term
- interpret interaction CI and p-value, not separate subgroup p-values

### 6.3 Recommended reporting
- pre-specify subgroups
- show forest plot with estimates and CI
- interpret cautiously, emphasize exploratory unless strongly powered

---

## 7. Key plots for clinical trial reporting

### 7.1 Continuous endpoint plots
- mean change over time by arm (with CI bands)
- distribution plots (box/violin) if appropriate

### 7.2 Binary endpoint plots
- bar plots of proportions with CI
- risk difference plot

### 7.3 Time-to-event plots
- Kaplan–Meier curves with number at risk
- cumulative incidence curves for competing risks

### 7.4 Forest plots
- subgroup effect estimates with CI
- helpful visual summary, but not proof of heterogeneity without interaction testing

---

## 8. Reproducibility: minimal standards

A reproducible trial analysis should:
- fix package versions when possible
- keep raw data separate from derived datasets
- use scripted pipelines for tables/figures
- include a clear “analysis dataset” construction step
- log all outputs used in reporting

Even for an educational site, teaching reproducibility builds professional habits.

---

# Part A — Practical reporting workflow in Python

We generate a toy parallel-group trial dataset and produce:
- Table 1 baseline summary
- primary endpoint estimate + CI
- subgroup forest plot (illustration)
- CONSORT-like counts summary

---

## 9A. Python: Create toy trial dataset

!!! interactive "Python"
    ```python
    import numpy as np
    import pandas as pd

    np.random.seed(100)

    n = 400
    trt = np.random.binomial(1, 0.5, n)

    age = np.random.normal(55, 12, n)
    sex = np.random.binomial(1, 0.48, n)  # 1=male
    baseline = np.random.normal(150, 18, n)  # baseline SBP

    # primary endpoint: SBP at week 12 (lower is better)
    true_delta = -6.0
    y = baseline + true_delta*trt + np.random.normal(0, 12, n)

    # dropout indicator (for reporting)
    dropout = np.random.binomial(1, 0.08, n)  # 8% missing
    y_obs = y.copy()
    y_obs[dropout==1] = np.nan

    df = pd.DataFrame({
        "trt": trt,
        "age": age,
        "sex": sex,
        "baseline_sbp": baseline,
        "sbp12": y_obs,
        "dropout": dropout
    })

    df.head()
    ```

---

## 10A. Python: Baseline table (“Table 1”)

We compute descriptive summaries by group.

!!! interactive "Python"
    ```python
    def summarize_cont(x):
        return pd.Series({"mean": x.mean(), "sd": x.std()})

    def summarize_bin(x):
        return pd.Series({"n": x.sum(), "pct": 100*x.mean()})

    table1_age = df.groupby("trt")["age"].apply(summarize_cont).unstack()
    table1_base = df.groupby("trt")["baseline_sbp"].apply(summarize_cont).unstack()
    table1_sex = df.groupby("trt")["sex"].apply(summarize_bin).unstack()

    table1_age, table1_base, table1_sex
    ```

You can format these for display in the site as markdown tables if desired.

---

## 11A. Python: Primary endpoint effect estimate with CI (ANCOVA)

ANCOVA often improves precision:
\[
\text{SBP}_{12} = \beta_0 + \beta_1 T + \beta_2 \text{baseline} + \epsilon
\]

Here \(\beta_1\) estimates adjusted mean difference.

!!! interactive "Python"
    ```python
    import statsmodels.formula.api as smf

    dfa = df.dropna(subset=["sbp12"]).copy()
    fit = smf.ols("sbp12 ~ trt + baseline_sbp", data=dfa).fit()
    fit.summary().tables[1]
    ```

Extract estimate + 95% CI:

!!! interactive "Python"
    ```python
    est = fit.params["trt"]
    ci = fit.conf_int().loc["trt"].tolist()
    est, ci
    ```

Interpretation:
- if estimate is negative, treatment lowers SBP relative to control

---

## 12A. Python: CONSORT-style counts 

!!! interactive "Python"
    ```python
    consort = pd.DataFrame({
        "randomized": df.groupby("trt").size(),
        "missing_primary": df.groupby("trt")["sbp12"].apply(lambda s: s.isna().sum()),
        "analyzed_primary": df.groupby("trt")["sbp12"].apply(lambda s: s.notna().sum())
    })
    consort
    ```

---

## 13A. Python: Forest plot for subgroup effects (illustration)

Subgroup: sex (male vs female)

We estimate treatment effect within subgroup (descriptive) and show forest plot.
The correct formal test is interaction (shown after).

!!! interactive "Python"
    ```python
    import matplotlib.pyplot as plt
    import numpy as np

    effects = []
    labels = []
    cis = []

    for s, name in [(0, "Female"), (1, "Male")]:
        sub = dfa[dfa.sex==s]
        f = smf.ols("sbp12 ~ trt + baseline_sbp", data=sub).fit()
        est = f.params["trt"]
        lo, hi = f.conf_int().loc["trt"]
        labels.append(name)
        effects.append(est)
        cis.append((lo, hi))

    effects = np.array(effects)
    lo = np.array([c[0] for c in cis])
    hi = np.array([c[1] for c in cis])

    y = np.arange(len(labels))

    plt.figure()
    plt.hlines(y, lo, hi)
    plt.plot(effects, y, marker="o", linestyle="None")
    plt.yticks(y, labels)
    plt.axvline(0, linestyle="--")
    plt.xlabel("Treatment effect (adjusted mean difference)")
    plt.title("Subgroup effects (illustration)")
    plt.show()
    ```

---

## 14A. Python: Interaction test (correct subgroup inference)

!!! interactive "Python"
    ```python
    fit_int = smf.ols("sbp12 ~ trt * sex + baseline_sbp", data=dfa).fit()
    fit_int.summary().tables[1]
    ```

Interpretation:
- the `trt:sex` term is the interaction
- if it is near 0 with wide CI, evidence for heterogeneity is weak

---

# Part B — Practical reporting workflow in R

---

## 15B. R: Create toy trial dataset

!!! interactive "R"
    ```r
    set.seed(100)

    n <- 400
    trt <- rbinom(n, 1, 0.5)

    age <- rnorm(n, 55, 12)
    sex <- rbinom(n, 1, 0.48)
    baseline <- rnorm(n, 150, 18)

    true_delta <- -6
    y <- baseline + true_delta*trt + rnorm(n, 0, 12)

    dropout <- rbinom(n, 1, 0.08)
    y_obs <- y
    y_obs[dropout==1] <- NA

    df <- data.frame(
      trt=trt,
      age=age,
      sex=sex,
      baseline_sbp=baseline,
      sbp12=y_obs,
      dropout=dropout
    )

    head(df)
    ```

---

## 16B. R: Baseline table (“Table 1”)

!!! interactive "R"
    ```r
    # continuous summaries
    tapply(df$age, df$trt, function(x) c(mean=mean(x), sd=sd(x)))
    tapply(df$baseline_sbp, df$trt, function(x) c(mean=mean(x), sd=sd(x)))

    # binary summary: sex=1
    tapply(df$sex, df$trt, function(x) c(n=sum(x), pct=100*mean(x)))
    ```

---

## 17B. R: Primary endpoint analysis with ANCOVA

!!! interactive "R"
    ```r
    dfa <- df[!is.na(df$sbp12), ]

    fit <- lm(sbp12 ~ trt + baseline_sbp, data=dfa)
    summary(fit)

    # treatment effect estimate and CI
    est <- coef(fit)["trt"]
    ci <- confint(fit)["trt", ]
    est
    ci
    ```

---

## 18B. R: CONSORT-style counts

!!! interactive "R"
    ```r
    randomized <- table(df$trt)
    missing_primary <- tapply(is.na(df$sbp12), df$trt, sum)
    analyzed_primary <- tapply(!is.na(df$sbp12), df$trt, sum)

    data.frame(
      trt = names(randomized),
      randomized = as.numeric(randomized),
      missing_primary = as.numeric(missing_primary),
      analyzed_primary = as.numeric(analyzed_primary)
    )
    ```

---

## 19B. R: Forest plot for subgroup effects (illustration)

!!! interactive "R"
    ```r
    # subgroup by sex
    effects <- c()
    lo <- c()
    hi <- c()
    labels <- c("Female","Male")

    for (s in c(0,1)) {
      sub <- dfa[dfa$sex==s, ]
      f <- lm(sbp12 ~ trt + baseline_sbp, data=sub)
      est <- coef(f)["trt"]
      ci <- confint(f)["trt", ]
      effects <- c(effects, est)
      lo <- c(lo, ci[1])
      hi <- c(hi, ci[2])
    }

    # simple base R forest-like plot
    plot(effects, 1:2, xlim=range(c(lo,hi)), yaxt="n",
         xlab="Treatment effect (adjusted mean difference)", ylab="",
         pch=19)
    axis(2, at=1:2, labels=labels)
    segments(lo, 1:2, hi, 1:2)
    abline(v=0, lty=2)
    title("Subgroup effects (illustration)")
    ```

---

## 20B. R: Interaction test (correct subgroup inference)

!!! interactive "R"
    ```r
    fit_int <- lm(sbp12 ~ trt*sex + baseline_sbp, data=dfa)
    summary(fit_int)
    confint(fit_int)["trt:sex", ]
    ```

---

## 21. How to write results 

A high-quality results paragraph usually includes:
- point estimate
- CI
- clear directionality
- clinical meaning

Example (continuous endpoint):
> At week 12, the adjusted mean SBP was lower in the treatment group compared with control, with an adjusted mean difference of -5.8 mmHg (95% CI -8.1 to -3.5).  

Do not over-focus on p-values; include them if required, but interpret the CI primarily.

---

## 22. Common interpretation mistakes 

- “No significant difference” does not mean “no difference”
  - check CI width and clinically relevant effects
- Subgroup results should not be claimed based on within-subgroup p-values
  - interaction test matters
- Avoid interpreting secondary endpoints as confirmatory unless multiplicity controlled
- Avoid causal language outside what the design supports (e.g., for post-hoc analyses)

---

## 23. Reproducible reporting workflow 

A simple, robust structure:

- `docs/` for narrative content and code examples
- `analysis/` for scripts/notebooks that generate figures/tables
- `data/` for toy or public datasets (avoid sensitive real trial data)
- `outputs/` for generated tables/figures used in docs

Even if your site is educational, this mirrors real-world trial reporting pipelines.

---

## 24. Exercises

<details>
<summary>Click to try</summary>

1. Modify the toy dataset so the treatment effect is -3 instead of -6 and compare CI width and conclusions.  
2. Increase dropout from 8% to 25% and recompute the primary analysis (complete-case). Discuss how missingness might bias conclusions.  
3. Create an additional subgroup (age < 60 vs ≥ 60) and make a forest plot. Then test interaction.  
4. For a binary endpoint, compute risk difference and 95% CI.  
5. Write a short CONSORT-style summary of participant flow using the computed counts.

</details>

---

## 25. Summary

- CONSORT provides a structured approach to transparent trial reporting.
- Report effect sizes and confidence intervals prominently.
- Baseline tables are descriptive; baseline p-values are discouraged.
- Subgroup claims require interaction testing and cautious interpretation.
- Reproducible scripts for tables/figures improve trust and reduce errors.

---
