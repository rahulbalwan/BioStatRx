# Kaplan–Meier Confidence Intervals and Median Survival 

Kaplan–Meier gives an estimate of survival:

\[
\hat S(t)
\]

But in biostatistics, **we never report an estimate without uncertainty**.

So in real papers you will see results like:

- “5-year survival = 0.72 (95% CI: 0.64–0.80)”
- “Median survival = 18.3 months (95% CI: 15.1–22.7)”
- “Median not reached”

This chapter teaches you how to compute and interpret with **both Python and R**.:

 - standard errors  
 - Greenwood variance  
 - confidence intervals (with proper transformation)  
 - median survival and its CI  
 - reporting standards in clinical papers  



---

## 1. Why we need confidence intervals

Suppose your KM estimate at 5 years is:

\[
\hat S(5)=0.72
\]

This is a sample estimate. If you repeated the study, you would not get exactly 0.72 again.

Confidence intervals quantify uncertainty:

\[
\hat S(5)=0.72,\quad 95\%\,CI=(0.64,0.80)
\]

Clinical interpretation:

> Based on the observed data, the true 5-year survival is plausibly between 64% and 80%.

---

## 2. Standard error and Greenwood’s formula

### 2.1 KM estimator reminder

At event times \(t_1, t_2, \dots\):

\[
\hat S(t)=\prod_{t_j \le t}\left(1-\frac{d_j}{n_j}\right)
\]

where:
- \(n_j\) = number at risk just before \(t_j\)
- \(d_j\) = number of events at \(t_j\)

---

### 2.2 Greenwood variance 

The classic large-sample variance estimator for KM is:

\[
\mathrm{Var}(\hat S(t)) = \hat S(t)^2 \sum_{t_j \le t}\frac{d_j}{n_j(n_j-d_j)}
\]

Standard error:

\[
SE(\hat S(t))=\sqrt{\mathrm{Var}(\hat S(t))}
\]

---

### 2.3 Intuition

Greenwood shows why uncertainty increases late in follow-up:

- As time increases, fewer people remain at risk → \(n_j\) becomes small  
- Small \(n_j\) makes \(\frac{1}{n_j(n_j-d_j)}\) large  
- So variance accumulates faster

Therefore:
 - early CI is narrow  
 - late CI is wide (often very wide)

---

## 3. Why naive CI can be wrong

A naive CI might be:

\[
\hat S(t) \pm 1.96\,SE(\hat S(t))
\]

Problem:
- It can go below 0 or above 1  
- Survival probabilities must stay in \([0,1]\)

So we typically use a transformation that keeps CI within [0,1].

---

## 4. Log–log confidence intervals 

A widely used approach is **log–log transformed CI** (often default in software).

### 4.1 Idea
Transform KM to a scale where normal approximation works better, compute CI, then transform back.

A common form:

$$
CI_{\text{lower}}(t)=\left[\hat S(t)\right]^{\exp\!\left(z\,SE^{*}(t)\right)}
$$

$$
CI_{\text{upper}}(t)=\left[\hat S(t)\right]^{\exp\!\left(-z\,SE^{*}(t)\right)}
$$


where \(z=1.96\) for 95% CI and \(SE^*(t)\) is a standard error on the transformed scale.

You do NOT need to memorize the exact formula—**know why we transform**:
- to keep CI within [0,1]
- to improve coverage when survival is near 0 or 1

Software handles it.

---

## 5. Median survival time

### 5.1 Definition

Median survival is the time when survival drops to 0.5:

\[
t_{0.5}=\inf\{t: \hat S(t)\le 0.5\}
\]

Interpretation:

> The time by which 50% of subjects have experienced the event.

This is a standard clinical summary because it’s robust and meaningful.

---

### 5.2 Why median is preferred over mean survival
Mean survival is hard because:
- censoring may prevent observing long tail  
- mean may not exist / may be unstable  
- with heavy censoring, mean is biased without strong assumptions

Median is robust and is the standard in oncology/clinical trials.

---

### 5.3 “Median not reached”
If \(\hat S(t)\) never drops below 0.5 during follow-up, then median survival is not estimable.

In papers:
> “Median survival not reached.”

This often happens when:
- follow-up is short  
- event rate is low  
- treatment is highly effective  

---

## 6. Visual intuition: CI and median

```
S(t)
1.0 |---------
    |   \    (CI narrow early)
0.5 |----\-----------------  ← median where curve crosses 0.5
    |      \__ (CI widens late)
0.0 +------------------------- time
```

---

## 7. Compute KM + CI + median in Python (lifelines)

### 7.1 Simulate survival data (same dataset used in later sections)

!!! interactive "Python"
    ```python
    import numpy as np
    import pandas as pd

    np.random.seed(123)

    n = 250

    # True event times
    T = np.random.exponential(scale=10, size=n)

    # Censoring times
    C = np.random.uniform(2, 18, size=n)

    time = np.minimum(T, C)
    event = (T <= C).astype(int)

    df = pd.DataFrame({"time": time, "event": event})
    df.head()
    ```

---

### 7.2 Fit Kaplan–Meier and plot CI

!!! interactive "Python"
    ```python
    import matplotlib.pyplot as plt
    from lifelines import KaplanMeierFitter

    km = KaplanMeierFitter()
    km.fit(df["time"], df["event"], label="KM")

    ax = km.plot_survival_function(ci_show=True)
    plt.title("Kaplan–Meier Survival with 95% CI")
    plt.xlabel("Time")
    plt.ylabel("S(t)")
    plt.show()
    ```

---

### 7.3 Extract survival at specific times (1, 5, 10 units)

!!! interactive "Python"
    ```python
    times_of_interest = [1, 5, 10]

    # lifelines provides survival function values at times
    surv = km.survival_function_at_times(times_of_interest)

    # confidence intervals at those times
    ci = km.confidence_interval_at_times(times_of_interest)

    print("Survival estimates:")
    print(surv)

    print("\n95% CI:")
    print(ci)
    ```

---

### 7.4 Median survival time

!!! interactive "Python"
    ```python
    km.median_survival_time_
    ```

If this returns `inf` or `None`, that often indicates median was not reached.

---

### 7.5 A “publication-style” summary table in Python

!!! interactive "Python"
    ```python
    summary = pd.DataFrame({
        "time": times_of_interest,
        "S(t)": km.survival_function_at_times(times_of_interest).values
    })

    ci = km.confidence_interval_at_times(times_of_interest)
    summary["CI_lower"] = ci.iloc[:, 0].values
    summary["CI_upper"] = ci.iloc[:, 1].values

    summary
    ```

---

## 8. Compute KM + CI + median in R (survival + survminer)

### 8.1 Simulate data in R

!!! interactive "R"
    ```r
    set.seed(123)

    n <- 250
    T <- rexp(n, rate = 1/10)      # mean 10
    C <- runif(n, min = 2, max = 18)

    time <- pmin(T, C)
    event <- as.integer(T <= C)

    df <- data.frame(time=time, event=event)
    head(df)
    ```

---

### 8.2 Fit Kaplan–Meier and plot with CI

!!! interactive "R"
    ```r
    library(survival)

    fit <- survfit(Surv(time, event) ~ 1, data=df)

    plot(fit, conf.int=TRUE,
         xlab="Time", ylab="S(t)",
         main="Kaplan–Meier Survival with 95% CI")
    ```

---

### 8.3 Extract survival estimates at specific times

In base `survival`, use `summary(fit, times=...)`.

!!! interactive "R"
    ```r
    times_of_interest <- c(1, 5, 10)

    s <- summary(fit, times = times_of_interest)

    out <- data.frame(
      time = s$time,
      surv = s$surv,
      lower = s$lower,
      upper = s$upper
    )

    out
    ```

---

### 8.4 Median survival and its CI in R

`survfit` stores median in the printed output. A standard way is:

!!! interactive "R"
    ```r
    # Median survival time:
    fit

    # A direct extraction (often works):
    fit$table["median"]
    ```

For median CI, many users rely on `survminer::surv_median()`:

!!! interactive "R"
    ```r
    # install.packages("survminer")  # if needed
    library(survminer)

    surv_median(fit)
    ```

If median not reached, it will show NA or similar outputs.

---

## 9. Group-specific median survival and CI (Python + R)

Many clinical studies compare groups.

### 9.1 Python: treatment vs control KM + median

!!! interactive "Python"
    ```python
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from lifelines import KaplanMeierFitter

    np.random.seed(7)

    n = 300
    group = np.random.binomial(1, 0.5, n)  # 0=control,1=treatment

    # Treatment has longer mean survival
    T0 = np.random.exponential(9, (group==0).sum())
    T1 = np.random.exponential(13, (group==1).sum())
    T = np.concatenate([T0, T1])

    C = np.random.uniform(2, 16, n)
    time = np.minimum(T, C)
    event = (T <= C).astype(int)

    df2 = pd.DataFrame({"time": time, "event": event, "group": group})

    km0 = KaplanMeierFitter()
    km1 = KaplanMeierFitter()

    ax = plt.subplot(111)
    km0.fit(df2.loc[df2.group==0, "time"], df2.loc[df2.group==0, "event"], label="Control").plot(ax=ax, ci_show=True)
    km1.fit(df2.loc[df2.group==1, "time"], df2.loc[df2.group==1, "event"], label="Treatment").plot(ax=ax, ci_show=True)

    plt.title("KM Curves with 95% CI by Group")
    plt.xlabel("Time")
    plt.ylabel("S(t)")
    plt.show()

    print("Median (Control):", km0.median_survival_time_)
    print("Median (Treatment):", km1.median_survival_time_)
    ```

---

### 9.2 R: treatment vs control KM + median

!!! interactive "R"
    ```r
    set.seed(7)

    n <- 300
    group <- rbinom(n, 1, 0.5)

    T0 <- rexp(sum(group==0), rate = 1/9)
    T1 <- rexp(sum(group==1), rate = 1/13)
    T <- c(T0, T1)

    C <- runif(n, 2, 16)
    time <- pmin(T, C)
    event <- as.integer(T <= C)

    df2 <- data.frame(time=time, event=event, group=factor(group, labels=c("Control","Treatment")))

    library(survival)

    fit_g <- survfit(Surv(time, event) ~ group, data=df2)

    plot(fit_g, col=c("black","red"), lty=1, conf.int=TRUE,
         xlab="Time", ylab="S(t)",
         main="KM Curves with 95% CI by Group")
    legend("topright", legend=levels(df2$group), col=c("black","red"), lty=1)

    # Median by group
    fit_g$table[,"median"]
    ```

---

## 10. How to report results 

### 10.1 Survival probability at fixed time
Example:

> “Five-year survival was 0.72 (95% CI: 0.64–0.80).”

### 10.2 Median survival
Example:

> “Median survival was 18.3 months (95% CI: 15.1–22.7).”

### 10.3 Median not reached
Example:

> “Median survival was not reached during follow-up.”

Always specify:
 - time unit (months/years)
 - method (Kaplan–Meier)
 - CI level (usually 95%)

---

## 11. Common mistakes and how to avoid them

### Mistake 1: CI outside [0,1]
 Fix: use log–log CI (software default).

### Mistake 2: interpreting CI overlap as hypothesis test
 CI overlap is not a formal test for group differences.
 Use log-rank or Cox for testing.

### Mistake 3: reporting mean survival without justification
Mean requires strong assumptions; median is standard.

### Mistake 4: reporting late survival estimates with tiny risk set
 Late tail can be unstable; always check number at risk.

---

## 12. Key takeaways

- Greenwood’s formula gives variance/SE for KM.
- CI should be constructed on a transformed scale (log–log) to stay within [0,1].
- Median survival is time where \(\hat S(t)\le 0.5\).
- If \(\hat S(t)\) never drops below 0.5, median is “not reached.”
- Always report survival estimates with CI in medical research.

---

## 13. Exercises

<details>
<summary>Click to try</summary>

1. Simulate a dataset with very heavy censoring and plot KM with CI. How does CI width change?  
2. Compute \(\hat S(5)\) and its CI in Python and in R. Confirm results are similar.  
3. Create data where median is not reached and explain how you would report it.  
4. Why does Greenwood variance increase late in follow-up?  
5. In group comparison, why can’t CI overlap be used as a test?

</details>
