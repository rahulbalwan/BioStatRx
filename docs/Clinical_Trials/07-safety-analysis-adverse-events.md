# Safety Analysis and Adverse Events (AE/SAE) + Practical Reporting

Efficacy is only half the story. In clinical trials, safety evaluation is equally important.

Safety analysis answers questions like:
- Is the new treatment associated with more adverse events?
- Are there specific event types that occur more often?
- Is the safety profile acceptable relative to benefits?
- Do adverse events depend on exposure time?

This chapter gives a practical, biostatistics-focused framework for analyzing and reporting safety data:
- AE vs SAE vs AESI definitions
- safety analysis populations
- incidence proportions vs exposure-adjusted incidence rates (EAIR)
- summary tables by event category (SOC/PT style)
- lab abnormalities and shift tables
- reproducible reporting in R and Python with toy data

---

## 1. Key safety terminology

Different trials have slightly different definitions, but the concepts are consistent.

### 1.1 Adverse event (AE)
Any untoward medical occurrence after treatment starts, whether or not related to treatment.

Important:
- “AE” does not require causality.
- If it happens after exposure, it can be recorded as an AE.

### 1.2 Serious adverse event (SAE)
Typically includes events such as:
- death
- life-threatening event
- hospitalization or prolonged hospitalization
- persistent disability or incapacity
- congenital anomaly
- other medically important event

These are high-priority events for monitoring and reporting.

### 1.3 Adverse event of special interest (AESI)
Events pre-specified as important due to:
- known drug class risks
- biological plausibility
- earlier trial signals

Example:
- thrombosis with certain therapies
- liver injury markers for hepatotoxic drugs

### 1.4 Treatment-emergent adverse event (TEAE)
A TEAE is an AE that begins (or worsens) after treatment initiation within a defined risk window.

Safety reporting often focuses on TEAEs because they align with exposure.

---

## 2. Safety analysis populations

Safety analysis is usually based on an “as-treated” or “safety set” population.

### 2.1 Safety set
All participants who received at least one dose of study treatment, analyzed according to treatment actually received.

This differs from ITT because:
- safety is linked to exposure, not assignment.

### 2.2 Why not always ITT for safety?
If someone is randomized but never takes treatment, counting them in safety summaries can dilute or distort safety risk.

Safety analyses often use:
- actual exposure time
- risk windows
- on-treatment vs intended treatment definitions

---

## 3. What to summarize in safety analysis?

Safety summaries usually include:

1) Overall AE burden
- number (%) with at least one AE
- number (%) with at least one SAE
- number (%) with discontinuation due to AE
- number (%) with death

2) Most common AEs
- by system organ class (SOC)
- by preferred term (PT)
- often sorted by frequency

3) Exposure-adjusted rates
- events per person-year (often per 100 person-years)
- important if follow-up differs between arms

4) Severity and relatedness
- mild/moderate/severe
- related vs unrelated to treatment (reported, but interpret cautiously)

5) Labs and vitals
- shifts from baseline category to worst post-baseline category
- threshold-based abnormality rates (e.g., ALT > 3× ULN)

---

## 4. Incidence proportion vs exposure-adjusted incidence rate (EAIR)

### 4.1 Incidence proportion (risk)

$$
\text{Risk} = \frac{\#\text{participants with at least one event}}{\#\text{participants at risk}}
$$

This answers:
> What proportion of participants experienced at least one event?

Use when:
- follow-up time is similar
- you want a participant-level interpretation

### 4.2 Exposure-adjusted incidence rate (EAIR)

$$
\text{EAIR} = \frac{\#\text{events}}{\text{total person-time}}
$$

Often expressed per 100 person-years:

$$
\text{EAIR}_{100} = 100 \times \frac{\#\text{events}}{\text{total person-years}}
$$

This answers:
> How frequently do events occur, accounting for time under observation?

Use when:
- follow-up differs between groups
- events can repeat (multiple infections, multiple hospitalizations)
- exposure time is clinically important

### 4.3 Which numerator to use for EAIR?
Two common choices:
- number of participants with at least one event (first-event rate)
- number of events (recurrent-event rate)

You must clearly specify which.

---

## 5. Comparing safety between groups (methods)

Safety comparisons are often descriptive, but inferential methods exist.

### 5.1 Risk difference / risk ratio
For “any AE yes/no”:
- compare proportions
- compute confidence intervals for risk difference or risk ratio

### 5.2 Poisson / negative binomial models for rates
For count outcomes (multiple events), with an offset:

$$
\log\bigl(E[Y]\bigr) = \beta_0 + \beta_1 \,\mathrm{Trt} + \log(\text{person-time})
$$

Then:
- \(\exp(\beta_1)\) is the rate ratio (treatment vs control)

If overdispersion exists, use negative binomial.

### 5.3 Time-to-event for first AE
You can use survival methods to analyze time to first AE.

Important nuance:
- time to first AE can be influenced by follow-up duration and competing risks (e.g., death).

---

# Part A — Practice in Python

We generate toy AE data with:
- participant-level exposure time
- SOC and PT categories
- recurrent event counts

---

## 6A. Python: Simulate a trial safety dataset

!!! interactive "Python"
    ```python
    import numpy as np
    import pandas as pd

    np.random.seed(2026)

    n = 300
    trt = np.random.binomial(1, 0.5, n)  # 1=treatment, 0=control

    # exposure time in days (treatment group slightly shorter, e.g., discontinuations)
    exposure_days = np.where(
        trt == 1,
        np.random.gamma(shape=5, scale=18, size=n),  # mean ~90
        np.random.gamma(shape=5, scale=20, size=n)   # mean ~100
    )

    exposure_py = exposure_days / 365.25

    # AE event generation: Poisson rate depends on treatment
    base_rate = 2.0  # events per person-year in control
    rate = np.where(trt == 1, base_rate * 1.25, base_rate)  # higher AE rate on treatment
    ae_count = np.random.poisson(rate * exposure_py)

    # create SOC/PT categories for each AE event (event-level dataset)
    socs = ["Gastrointestinal disorders", "Nervous system disorders", "Infections", "Skin disorders"]
    pts_by_soc = {
        "Gastrointestinal disorders": ["Nausea", "Diarrhoea", "Abdominal pain"],
        "Nervous system disorders": ["Headache", "Dizziness"],
        "Infections": ["Upper respiratory infection", "Urinary tract infection"],
        "Skin disorders": ["Rash", "Pruritus"]
    }

    rows = []
    for pid in range(1, n + 1):
        k = ae_count[pid - 1]
        if k == 0:
            continue
        for _ in range(k):
            soc = np.random.choice(socs, p=[0.28, 0.26, 0.26, 0.20])
            pt = np.random.choice(pts_by_soc[soc])
            severity = np.random.choice(["Mild", "Moderate", "Severe"], p=[0.65, 0.30, 0.05])
            sae = np.random.binomial(1, 0.04)  # rare SAEs
            rows.append([pid, trt[pid - 1], soc, pt, severity, sae])

    ae = pd.DataFrame(rows, columns=["participant_id", "trt", "SOC", "PT", "severity", "SAE"])

    # participant-level frame
    subj = pd.DataFrame({
        "participant_id": range(1, n + 1),
        "trt": trt,
        "exposure_days": exposure_days,
        "exposure_py": exposure_py,
        "ae_count": ae_count
    })

    subj.head(), ae.head()
    ```

---

## 7A. Python: Participant-level safety summary (at least one AE, at least one SAE)

!!! interactive "Python"
    ```python
    # indicator: at least one AE
    subj["any_AE"] = (subj["ae_count"] > 0).astype(int)

    # indicator: any SAE (participant-level)
    if len(ae) > 0:
        any_sae = ae.groupby("participant_id")["SAE"].max()
        subj = subj.merge(any_sae.rename("any_SAE"), left_on="participant_id", right_index=True, how="left")
        subj["any_SAE"] = subj["any_SAE"].fillna(0).astype(int)
    else:
        subj["any_SAE"] = 0

    # proportions by arm
    summary = subj.groupby("trt")[["any_AE", "any_SAE"]].mean()
    summary
    ```

Interpretation:
- values are proportions (risk of any AE/SAE)

---

## 8A. Python: Exposure-adjusted incidence rate (EAIR)

We compute:
- total events / total person-years

!!! interactive "Python"
    ```python
    eair = subj.groupby("trt").apply(lambda d: d["ae_count"].sum() / d["exposure_py"].sum())
    eair
    ```

Convert to events per 100 person-years:

!!! interactive "Python"
    ```python
    (100 * eair).rename("events_per_100_person_years")
    ```

---

## 9A. Python: SOC and PT frequency tables

Often safety tables report:
- number (%) of participants with at least one event in a category

We compute participant-level incidence by SOC and PT.

!!! interactive "Python"
    ```python
    # participant-level SOC incidence: unique participants within each SOC
    if len(ae) > 0:
        soc_inc = (ae.groupby(["trt", "SOC", "participant_id"])
                     .size()
                     .reset_index(name="n")
                     .groupby(["trt", "SOC"])
                     .size()
                     .reset_index(name="n_participants"))

        denom = subj.groupby("trt").size().rename("N").reset_index()

        soc_inc = soc_inc.merge(denom, on="trt")
        soc_inc["pct"] = 100 * soc_inc["n_participants"] / soc_inc["N"]

        soc_inc.sort_values(["trt", "n_participants"], ascending=[True, False]).head(10)
    ```

Similarly for PT:

!!! interactive "Python"
    ```python
    if len(ae) > 0:
        pt_inc = (ae.groupby(["trt", "PT", "participant_id"])
                    .size()
                    .reset_index(name="n")
                    .groupby(["trt", "PT"])
                    .size()
                    .reset_index(name="n_participants"))

        denom = subj.groupby("trt").size().rename("N").reset_index()

        pt_inc = pt_inc.merge(denom, on="trt")
        pt_inc["pct"] = 100 * pt_inc["n_participants"] / pt_inc["N"]

        pt_inc.sort_values(["trt", "n_participants"], ascending=[True, False]).head(12)
    ```

---

## 10A. Python: Modeling AE counts with Poisson regression (rate ratio)

We fit:

$$
\log\bigl(E[Y]\bigr) = \beta_0 + \beta_1 \,\mathrm{trt} + \log(\text{person-years})
$$

so \(\exp(\beta_1)\) is the AE rate ratio.

!!! interactive "Python"
    ```python
    import statsmodels.api as sm
    import numpy as np

    X = sm.add_constant(subj["trt"])
    model = sm.GLM(
        subj["ae_count"],
        X,
        family=sm.families.Poisson(),
        offset=np.log(subj["exposure_py"])
    )
    fit = model.fit()
    fit.summary().tables[1]
    ```

Rate ratio estimate:

!!! interactive "Python"
    ```python
    rr = float(np.exp(fit.params["trt"]))
    rr
    ```

Overdispersion check (quick heuristic):
- compare deviance to df_resid

!!! interactive "Python"
    ```python
    fit.deviance / fit.df_resid
    ```

If this is much larger than ~1, consider negative binomial.

---

# Part B — Practice in R

---

## 11B. R: Simulate safety dataset

!!! interactive "R"
    ```r
    set.seed(2026)

    n <- 300
    trt <- rbinom(n, 1, 0.5)

    exposure_days <- ifelse(
      trt == 1,
      rgamma(n, shape=5, scale=18),
      rgamma(n, shape=5, scale=20)
    )
    exposure_py <- exposure_days / 365.25

    base_rate <- 2.0
    rate <- ifelse(trt == 1, base_rate * 1.25, base_rate)
    ae_count <- rpois(n, rate * exposure_py)

    subj <- data.frame(
      participant_id = 1:n,
      trt = trt,
      exposure_days = exposure_days,
      exposure_py = exposure_py,
      ae_count = ae_count
    )

    head(subj)
    ```

Create event-level AE data with SOC/PT:

!!! interactive "R"
    ```r
    socs <- c("Gastrointestinal disorders", "Nervous system disorders", "Infections", "Skin disorders")

    pts_by_soc <- list(
      "Gastrointestinal disorders" = c("Nausea", "Diarrhoea", "Abdominal pain"),
      "Nervous system disorders" = c("Headache", "Dizziness"),
      "Infections" = c("Upper respiratory infection", "Urinary tract infection"),
      "Skin disorders" = c("Rash", "Pruritus")
    )

    ae <- data.frame()

    for (pid in 1:n) {
      k <- ae_count[pid]
      if (k == 0) next

      for (j in 1:k) {
        soc <- sample(socs, 1, prob=c(0.28, 0.26, 0.26, 0.20))
        pt <- sample(pts_by_soc[[soc]], 1)
        severity <- sample(c("Mild", "Moderate", "Severe"), 1, prob=c(0.65, 0.30, 0.05))
        sae <- rbinom(1, 1, 0.04)

        ae <- rbind(ae, data.frame(
          participant_id = pid,
          trt = trt[pid],
          SOC = soc,
          PT = pt,
          severity = severity,
          SAE = sae
        ))
      }
    }

    head(ae)
    ```

---

## 12B. R: Participant-level safety summary

!!! interactive "R"
    ```r
    subj$any_AE <- as.integer(subj$ae_count > 0)

    if (nrow(ae) > 0) {
      any_sae <- aggregate(SAE ~ participant_id, data=ae, FUN=max)
      subj <- merge(subj, any_sae, by="participant_id", all.x=TRUE)
      subj$SAE[is.na(subj$SAE)] <- 0
      names(subj)[names(subj) == "SAE"] <- "any_SAE"
    } else {
      subj$any_SAE <- 0
    }

    aggregate(cbind(any_AE, any_SAE) ~ trt, data=subj, FUN=mean)
    ```

---

## 13B. R: Exposure-adjusted incidence rate (EAIR)

!!! interactive "R"
    ```r
    eair <- aggregate(cbind(ae_count, exposure_py) ~ trt, data=subj, FUN=sum)
    eair$rate <- eair$ae_count / eair$exposure_py
    eair$rate_per_100py <- 100 * eair$rate
    eair
    ```

---

## 14B. R: SOC and PT incidence (participants with at least one event)

!!! interactive "R"
    ```r
    denom <- table(subj$trt)

    if (nrow(ae) > 0) {
      soc_inc <- unique(ae[, c("trt", "SOC", "participant_id")])
      soc_tab <- as.data.frame(table(soc_inc$trt, soc_inc$SOC))
      names(soc_tab) <- c("trt", "SOC", "n_participants")
      soc_tab$N <- as.numeric(denom[as.character(soc_tab$trt)])
      soc_tab$pct <- 100 * soc_tab$n_participants / soc_tab$N

      soc_tab[order(soc_tab$trt, -soc_tab$n_participants), ][1:10, ]
    }
    ```

PT incidence:

!!! interactive "R"
    ```r
    if (nrow(ae) > 0) {
      pt_inc <- unique(ae[, c("trt", "PT", "participant_id")])
      pt_tab <- as.data.frame(table(pt_inc$trt, pt_inc$PT))
      names(pt_tab) <- c("trt", "PT", "n_participants")
      pt_tab$N <- as.numeric(denom[as.character(pt_tab$trt)])
      pt_tab$pct <- 100 * pt_tab$n_participants / pt_tab$N

      pt_tab[order(pt_tab$trt, -pt_tab$n_participants), ][1:12, ]
    }
    ```

---

## 15B. R: Poisson regression for AE rate ratio (with offset)

!!! interactive "R"
    ```r
    fit <- glm(ae_count ~ trt + offset(log(exposure_py)), family=poisson(), data=subj)
    summary(fit)

    rr <- exp(coef(fit)["trt"])
    rr
    ```

Overdispersion check:

!!! interactive "R"
    ```r
    deviance(fit) / df.residual(fit)
    ```

If much larger than 1, consider negative binomial:

!!! interactive "R"
    ```r
    # install.packages("MASS") if needed
    library(MASS)

    fit_nb <- glm.nb(ae_count ~ trt + offset(log(exposure_py)), data=subj)
    summary(fit_nb)

    rr_nb <- exp(coef(fit_nb)["trt"])
    rr_nb
    ```

---

## 16. Lab safety and shift tables

Many trials monitor labs such as ALT/AST, creatinine, and neutrophils.

A common reporting tool is a shift table:
- baseline category (e.g., normal / high)
- worst post-baseline category

Example categories:
- Normal
- Grade 1
- Grade 2+

based on clinical thresholds.

Shift tables are descriptive but very informative for safety evaluation.

---

## 17. Practical reporting guidance

Safety reporting typically includes:
- overall AE/SAE/discontinuation/death summary
- most common AEs (SOC/PT)
- exposure-adjusted rates if follow-up differs
- severity breakdown
- relationship to treatment (reported but interpret cautiously)
- lab abnormality summaries

Important reporting principle:
- always provide denominators (N per arm)
- clearly define the risk window (what counts as treatment-emergent)

---

## 18. Exercises

<details>
<summary>Click to try</summary>

1. Change `base_rate` from 2.0 to 1.0 and treatment multiplier from 1.25 to 1.10. Recompute EAIR and Poisson RR.  
2. Make exposure time identical in both arms and compare incidence proportions vs EAIR. Do conclusions change?  
3. Simulate a scenario with strong overdispersion (e.g., add subject-level frailty) and compare Poisson vs negative binomial.  
4. Create a table of the top 5 PTs by incidence for each arm.  
5. Define an AESI category (e.g., infections) and compute incidence proportion and rate ratio for that subset.

</details>

---

## 19. Summary

- Safety analysis focuses on AEs/SAEs/AESIs, often in the safety set (exposed participants).
- Incidence proportions describe how many participants experience events.
- Exposure-adjusted rates account for different follow-up and recurrent events.
- Poisson/negative binomial models allow rate ratio estimation with person-time offsets.
- Clear, denominator-aware reporting is essential for interpretable safety results.
