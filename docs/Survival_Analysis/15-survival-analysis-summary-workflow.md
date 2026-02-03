# Survival Analysis Reporting, Tables, and Publication-Quality Plots

You now know the core survival methods. The next step is to present results **like a real biostatistics paper**.

This chapter focuses on:

 - what to report (minimum standards)  
 - how to write survival results in Methods + Results  
 - publication-quality KM plots with number-at-risk tables  
 - Cox regression tables with HR, CI, p-values  
 - common reporting mistakes to avoid  
 - Python + R code for clean outputs  

---

## 1. What survival analysis results should include

### 1.1 Kaplan–Meier (KM)
A KM analysis should report:

 - definition of time origin and endpoint
 - number of subjects and number of events
 - median follow-up time (often)
 - median survival with 95% CI (or “not reached”)
 - survival probabilities at clinically relevant times (1-year, 5-year)
 - KM plot with censor ticks and number-at-risk table
 - group comparison p-value (log-rank)

### 1.2 Cox regression
A Cox analysis should report:

 - covariates included (and why)
 - HR, 95% CI, p-values
 - PH assumption checks (Schoenfeld)
 - model performance (optional but common: concordance/C-index)
 - how missing data was handled
 - robust SE or clustering strategy if applicable

---

## 2. Recommended workflow 

A typical clinical paper flow:

 1) Descriptive summary: baseline table (not covered deeply here)
 2) KM curves by group + log-rank p-value
 3) Cox model (unadjusted then adjusted)
 4) PH checks
 5) Sensitivity analyses 
 6) Report clinically meaningful estimates (5-year survival, median)

---

## 3. A “Methods” template paragraph 

### 3.1 Methods paragraph (KM + log-rank + Cox)

> Time-to-event was defined as time from [time zero] to [event definition]. Patients without an event were censored at last follow-up. Kaplan–Meier methods were used to estimate survival functions and median survival times with 95% confidence intervals. Survival distributions were compared using the log-rank test. Cox proportional hazards regression was used to estimate hazard ratios (HRs) and 95% confidence intervals, adjusting for [covariates]. Proportional hazards assumptions were assessed using Schoenfeld residuals. All analyses were conducted in [R/Python], and statistical tests were two-sided with α=0.05.

Replace:
 - [time zero]
 - [event definition]
 - [covariates]
 - [R/Python]

---

## 4. A “Results” template paragraph 

> Among N patients, E events occurred during follow-up. Median follow-up was M months. The Kaplan–Meier estimated median survival was X months (95% CI: L–U). Five-year survival was S (95% CI: L–U). Survival differed between groups (log-rank p=...). In univariable Cox regression, treatment was associated with [higher/lower] hazard (HR ..., 95% CI ..., p=...). After adjustment for age, sex, and stage, the association remained [similar/attenuated] (adjusted HR ..., 95% CI ..., p=...). No violation of proportional hazards was observed based on Schoenfeld residuals.

---

## 5. Publication-quality Kaplan–Meier plots

What makes a plot “publication quality”?

 - clear labels + units  
 - censor tick marks  
 - confidence bands (optional)  
 - number-at-risk table  
 - clean style and legible fonts  
 - not over-cluttered  
 - shows group legend clearly  

---

# PART A — R 

R ecosystem is strongest for survival plotting via `survminer`.

---

## 6A. R: KM plot + risk table (survminer)

### 6A.1 Example data simulation

!!! interactive "R"
    ```r
    set.seed(55)

    n <- 300
    group <- rbinom(n, 1, 0.5)
    age <- rnorm(n, 60, 10)

    # Different hazards for groups
    T0 <- rexp(sum(group==0), rate=1/10)
    T1 <- rexp(sum(group==1), rate=1/14)
    T <- c(T0, T1)

    C <- runif(n, 2, 18)

    time <- pmin(T, C)
    event <- as.integer(T <= C)

    df <- data.frame(time=time, event=event, group=factor(group, labels=c("Control","Treatment")), age=age)

    head(df)
    ```

### 6A.2 KM fit + log-rank

!!! interactive "R"
    ```r
    library(survival)

    fit_km <- survfit(Surv(time, event) ~ group, data=df)
    fit_km
    ```

### 6A.3 Publication KM plot (risk table)

!!! interactive "R"
    ```r
    # install.packages("survminer") # if needed
    library(survminer)

    ggsurvplot(
      fit_km, data=df,
      conf.int=TRUE,
      risk.table=TRUE,
      censor=TRUE,
      pval=TRUE,
      xlab="Time (months)",
      ylab="Survival probability",
      title="Kaplan–Meier Survival by Group",
      legend.title="Group",
      ggtheme=theme_minimal()
    )
    ```

This produces:
 - KM curve with CI
 - censor ticks
 - log-rank p-value
 - number-at-risk table

---

## 7A. R: Extracting survival at fixed times (table-ready)

!!! interactive "R"
    ```r
    times <- c(6, 12, 24)

    s <- summary(fit_km, times=times)

    out <- data.frame(
      group = rep(names(s$strata), each=length(times)),
      time = s$time,
      surv = s$surv,
      lower = s$lower,
      upper = s$upper
    )

    out
    ```

---

## 8A. R: Cox regression table (HR, CI, p-value)

!!! interactive "R"
    ```r
    fit_cox <- coxph(Surv(time, event) ~ group + age, data=df)
    summary(fit_cox)
    ```

To create a nice table:

!!! interactive "R"
    ```r
    hr <- exp(coef(fit_cox))
    ci <- exp(confint(fit_cox))
    p <- summary(fit_cox)$coefficients[,5]

    table_out <- data.frame(
      Variable = names(hr),
      HR = hr,
      CI_lower = ci[,1],
      CI_upper = ci[,2],
      p_value = p
    )

    table_out
    ```

---

## 9A. R: PH check reporting (cox.zph)

!!! interactive "R"
    ```r
    z <- cox.zph(fit_cox)
    z
    plot(z)
    ```

Methods writing:
> “Proportional hazards assumptions were assessed with Schoenfeld residuals.”

---

# PART B — PYTHON (clean plotting + tables)

Python can do excellent survival analysis with lifelines,
but risk tables require manual work or custom plotting.

We give a clean practical pipeline.

---

## 10B. Python: KM plot with CI and censor ticks

### 10B.1 Simulate data 

!!! interactive "Python"
    ```python
    import numpy as np
    import pandas as pd

    np.random.seed(55)

    n = 300
    group = np.random.binomial(1, 0.5, n)
    age = np.random.normal(60, 10, n)

    T0 = np.random.exponential(10, (group==0).sum())
    T1 = np.random.exponential(14, (group==1).sum())
    T = np.concatenate([T0, T1])

    C = np.random.uniform(2, 18, n)

    time = np.minimum(T, C)
    event = (T <= C).astype(int)

    df = pd.DataFrame({"time": time, "event": event, "group": group, "age": age})
    df.head()
    ```

### 10B.2 KM plot

!!! interactive "Python"
    ```python
    import matplotlib.pyplot as plt
    from lifelines import KaplanMeierFitter

    km0 = KaplanMeierFitter()
    km1 = KaplanMeierFitter()

    ax = plt.subplot(111)

    km0.fit(df.loc[df.group==0, "time"], df.loc[df.group==0, "event"], label="Control").plot(ax=ax, ci_show=True)
    km1.fit(df.loc[df.group==1, "time"], df.loc[df.group==1, "event"], label="Treatment").plot(ax=ax, ci_show=True)

    plt.title("Kaplan–Meier Survival by Group")
    plt.xlabel("Time (months)")
    plt.ylabel("Survival probability")
    plt.show()
    ```

---

## 11B. Python: number-at-risk table 

We compute number at risk at chosen times.

!!! interactive "Python"
    ```python
    import numpy as np
    import pandas as pd

    times = [0, 6, 12, 18]

    def n_at_risk(df_sub, times):
        return [int((df_sub["time"] >= t).sum()) for t in times]

    risk = pd.DataFrame({
        "time": times,
        "Control": n_at_risk(df[df.group==0], times),
        "Treatment": n_at_risk(df[df.group==1], times)
    })

    risk
    ```

You can print this table under the plot in MkDocs, or convert to markdown.

---

## 12B. Python: log-rank p-value

!!! interactive "Python"
    ```python
    from lifelines.statistics import logrank_test

    res = logrank_test(
        df.loc[df.group==0, "time"],
        df.loc[df.group==1, "time"],
        event_observed_A=df.loc[df.group==0, "event"],
        event_observed_B=df.loc[df.group==1, "event"]
    )

    res.p_value, res.test_statistic
    ```

---

## 13B. Python: Cox regression table

!!! interactive "Python"
    ```python
    from lifelines import CoxPHFitter

    df2 = df.copy()
    df2["group"] = df2["group"].astype(int)

    cph = CoxPHFitter()
    cph.fit(df2[["time","event","group","age"]], duration_col="time", event_col="event")

    cph.summary
    ```

To produce a clean HR table:

!!! interactive "Python"
    ```python
    out = cph.summary.copy()
    out["HR"] = np.exp(out["coef"])
    out["CI_lower"] = np.exp(out["coef lower 95%"])
    out["CI_upper"] = np.exp(out["coef upper 95%"])

    out[["HR","CI_lower","CI_upper","p"]]
    ```

---

## 14. Reporting checklist 

### KM
 - Define time zero and event precisely
 - Report N, events, censoring
 - Show KM curve with censor ticks
 - Provide median survival + CI (or not reached)
 - Provide survival at relevant times + CI
 - Provide number at risk table
 - Log-rank p-value

### Cox
 - Report covariates and coding
 - HR + 95% CI + p-value
 - Check PH (Schoenfeld)
 - Handle clustering if needed
 - Mention missing data handling
 - Provide model performance 

---

## 15. Common reporting mistakes 

### Mistake 1: “Survival differed” without p-value or effect size
 Fix: report log-rank p and Cox HR.

### Mistake 2: HR reported without CI
 Fix: always include 95% CI.

### Mistake 3: “Median survival” with no time unit
 Fix: specify months/years.

### Mistake 4: not stating censoring rules
 Fix: define censoring explicitly.

### Mistake 5: not checking PH
 Fix: report Schoenfeld-based check.

---

## 16. Key takeaways

 - Good survival reporting requires both estimates and uncertainty.
 - Use KM + risk table + log-rank for descriptive comparison.
 - Use Cox for adjusted HR with CI and p-values.
 - Always mention PH checking.
 - Use R survminer for best publication plots; Python lifelines is strong but risk tables are manual.

---

## 17. Exercises

<details>
<summary>Click to try</summary>

 1. Make a KM plot with risk table (R) and replicate the same in Python using manual risk table.  
 2. Write a Methods paragraph for a study with time-to-relapse outcome.  
 3. Write a Results paragraph reporting median survival and adjusted HR.  
 4. Fit Cox, check PH, and report how you would write it in a paper.  
 5. Create a table of 1-, 3-, 5-year survival with 95% CI for each group.

</details>
