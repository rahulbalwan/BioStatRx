# Log-Rank Test (Comparing Survival Curves) 

Kaplan–Meier curves give a **visual comparison** between groups.

But in biostatistics, we must answer:

> Are these survival curves *statistically different*, or is the difference due to chance?

The standard hypothesis test for comparing **entire survival curves** is:

# Log-rank test

This chapter explains:
- intuition (observed vs expected)
- full formulas
- what the p-value means clinically
- assumptions and limitations
- implementation in Python and R

---

## 1. Problem setup

Suppose we have two groups:

- Group 0: Control
- Group 1: Treatment

We observe:
- follow-up time
- event indicator
- group membership

We want to test:

> Do the groups have the same survival experience over time?

---

## 2. Hypotheses

### Null hypothesis
\[
H_0: S_0(t)=S_1(t)\ \text{for all } t
\]

Meaning:
- no difference in survival curves

### Alternative hypothesis
\[
H_A: S_0(t)\ne S_1(t)
\]

Meaning:
- survival differs at some point in time

---

## 3. The core intuition: observed vs expected

At each **event time**, consider the risk set:

- how many people are at risk in each group?
- how many events occur overall?
- how many events should occur in each group if \(H_0\) were true?

If the observed events differ strongly from expected events across time,
we reject \(H_0\).

---

## 4. Notation at an event time \(t_j\)

At event time \(t_j\):

- \(n_j\) = total number at risk just before \(t_j\)
- \(n_{0j}\), \(n_{1j}\) = at risk in group 0 and group 1
- \(d_j\) = total events at \(t_j\)
- \(d_{0j}\), \(d_{1j}\) = observed events in each group

---

## 5. Expected events under the null

If the groups are identical (null hypothesis true),
events should be split according to risk set proportions.

Expected events in group 1:

\[
E_{1j}=d_j\frac{n_{1j}}{n_j}
\]

Expected events in group 0:

\[
E_{0j}=d_j\frac{n_{0j}}{n_j}
\]

---

## 6. Log-rank test statistic

Define the difference:

\[
O_1 - E_1 = \sum_j (d_{1j}-E_{1j})
\]

Variance of \(O_1-E_1\) is approximately:

\[
Var(O_1-E_1)=\sum_j
\frac{n_{0j}n_{1j}d_j(n_j-d_j)}{n_j^2(n_j-1)}
\]

The test statistic:

\[
Z=\frac{O_1-E_1}{\sqrt{Var(O_1-E_1)}}
\]

Often reported as chi-square:

\[
\chi^2 = Z^2 \sim \chi^2_1
\]

So we compute:
- \(p\)-value from \(\chi^2_1\)

---

## 7. How to interpret the p-value clinically

### If p < 0.05
We have evidence that survival differs between groups.

Clinical writing example:
> “Survival curves differed significantly (log-rank p = 0.01).”

### If p ≥ 0.05
No evidence of difference.

Example:
> “No significant difference was observed (log-rank p = 0.42).”

Important:
- A non-significant p-value does NOT prove curves are equal.
- It may mean sample size is small or censoring is heavy.

---

## 8. What log-rank is sensitive to

Log-rank is most powerful when:

hazard ratio is roughly constant over time  
(i.e., proportional hazards holds approximately)

If curves cross strongly, log-rank can lose power.

---

## 9. Limitations and special cases

### 9.1 Crossing survival curves (non-proportional hazards)
If treatment helps early but harms later, curves may cross.

Log-rank (which weights all times equally) can become misleading.

In such cases consider:
- weighted log-rank tests
- Cox models with time interaction
- restricted mean survival time (RMST)

### 9.2 Heavy censoring late in follow-up
Risk sets become tiny; late information is unstable.

Always check “number at risk” tables.

### 9.3 Ties
Multiple events at the same time:
- handled in software (Efron/Breslow methods)

---

## 10. A worked mini-example (conceptual)

Suppose at one event time:

- risk set: 60 treated, 40 control (total 100)
- total events at this time: 10

Expected treated events:
\[
E=10\cdot\frac{60}{100}=6
\]

If observed treated events are 2, then:
\[
O-E = 2-6 = -4
\]

This indicates treatment may reduce event risk at that time.

Log-rank repeats this at all event times and sums.

---

## 11. Implementation in Python (lifelines)

We will:
1) simulate two-group survival data  
2) plot KM curves  
3) run log-rank test  

### 11.1 Python simulation + KM curves

!!! interactive "Python"
    ```python
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from lifelines import KaplanMeierFitter
    from lifelines.statistics import logrank_test

    np.random.seed(202)

    n = 400
    group = np.random.binomial(1, 0.5, n)  # 0=control, 1=treatment

    # Control has shorter mean survival
    T0 = np.random.exponential(8, (group==0).sum())
    # Treatment has longer mean survival
    T1 = np.random.exponential(12, (group==1).sum())

    T = np.concatenate([T0, T1])

    # Censoring
    C = np.random.uniform(2, 15, n)

    time = np.minimum(T, C)
    event = (T <= C).astype(int)

    df = pd.DataFrame({"time": time, "event": event, "group": group})

    km0 = KaplanMeierFitter()
    km1 = KaplanMeierFitter()

    ax = plt.subplot(111)
    km0.fit(df.loc[df.group==0, "time"], df.loc[df.group==0, "event"], label="Control").plot(ax=ax)
    km1.fit(df.loc[df.group==1, "time"], df.loc[df.group==1, "event"], label="Treatment").plot(ax=ax)

    plt.title("KM Curves (Control vs Treatment)")
    plt.xlabel("Time")
    plt.ylabel("S(t)")
    plt.show()
    ```

---

### 11.2 Python: log-rank test

!!! interactive "Python"
    ```python
    result = logrank_test(
        df.loc[df.group==0, "time"],
        df.loc[df.group==1, "time"],
        event_observed_A=df.loc[df.group==0, "event"],
        event_observed_B=df.loc[df.group==1, "event"]
    )

    print(result.summary)
    ```

Output includes:
- test statistic
- p-value

---

### 11.3 Sensitivity experiments (recommended)
Try:
- Make group means equal → p-value should be large
- Make treatment much better → p-value should be tiny
- Increase censoring → test loses power

---

## 12. Implementation in R (survival package)

We will:
1) simulate data  
2) fit KM curves  
3) run log-rank test via `survdiff()`  

### 12.1 R simulation + KM curves

!!! interactive "R"
    ```r
    set.seed(202)

    n <- 400
    group <- rbinom(n, 1, 0.5)

    T0 <- rexp(sum(group==0), rate = 1/8)
    T1 <- rexp(sum(group==1), rate = 1/12)
    T <- c(T0, T1)

    C <- runif(n, 2, 15)

    time <- pmin(T, C)
    event <- as.integer(T <= C)

    df <- data.frame(time=time, event=event, group=factor(group, labels=c("Control","Treatment")))

    library(survival)

    fit <- survfit(Surv(time, event) ~ group, data=df)

    plot(fit, col=c("black","red"), lty=1,
         xlab="Time", ylab="S(t)",
         main="KM Curves (Control vs Treatment)")
    legend("topright", legend=levels(df$group), col=c("black","red"), lty=1)
    ```

---

### 12.2 R: log-rank test

In R, log-rank is:

`survdiff(Surv(time, event) ~ group)`

!!! interactive "R"
    ```r
    lr <- survdiff(Surv(time, event) ~ group, data=df)
    lr
    ```

To compute p-value:

!!! interactive "R"
    ```r
    pval <- 1 - pchisq(lr$chisq, df=1)
    pval
    ```

---

## 13. Log-rank vs Cox: what’s the difference?

### Log-rank test
- compares groups only
- gives p-value only
- does NOT adjust for covariates

### Cox regression
- gives hazard ratio (effect size)
- adjusts for covariates (age, sex, biomarkers)
- provides CI, p-values, model building

Workflow in medical papers:
1. KM plot
2. log-rank p-value
3. Cox hazard ratio (adjusted)

---

## 14. Key takeaways

- Log-rank tests whether survival curves differ across the entire follow-up.
- At each event time, compare observed vs expected events using risk sets.
- Test statistic is approximately chi-square with 1 df for two groups.
- Works best when proportional hazards roughly holds.
- Use Cox regression for adjusted hazard ratios and effect size.

---

## 15. Exercises

<details>
<summary>Click to try</summary>

1. Explain why log-rank uses risk sets rather than full sample size.  
2. Simulate data where two groups have equal hazards. What p-value do you get?  
3. Simulate crossing survival curves (early benefit, late harm). What happens to log-rank?  
4. Why is log-rank not a substitute for Cox regression?  
5. Fit KM in R and Python and confirm log-rank p-values agree.

</details>
