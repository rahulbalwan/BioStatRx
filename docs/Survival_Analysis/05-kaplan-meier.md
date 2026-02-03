# Kaplan–Meier Estimator (KM)

Kaplan–Meier (KM) is the **standard nonparametric method** to estimate the survival function:

\[
S(t)=P(T>t)
\]

KM is the first survival method you should learn because it is:

- intuitive
- clinically interpretable
- handles censoring correctly
- the foundation for log-rank tests and Cox models

In most medical papers, KM curves are the first survival result shown.

---

## 1. What Kaplan–Meier estimates

KM estimates survival probability over time:

\[
\hat S(t) \approx P(\text{survive beyond } t)
\]

Typical reported outputs:
- 1-year survival probability
- 5-year survival probability
- median survival time
- survival curves by treatment group

---

## 2. Key idea: conditional survival at each event time

KM does NOT compute “survival at time \(t\) directly.”

Instead, it breaks survival into steps at each **event time**:

At each event time \(t_j\):
- \(n_j\) = number at risk just before \(t_j\)
- \(d_j\) = number of events at \(t_j\)

Conditional probability of surviving that step:

\[
P(\text{survive past } t_j \mid \text{alive just before } t_j)
= 1 - \frac{d_j}{n_j}
\]

Then multiply across event times:

\[
\hat S(t) =
\prod_{t_j \le t}\left(1-\frac{d_j}{n_j}\right)
\]

This is called the **product-limit estimator**.

---

## 3. How censoring affects KM

Censoring:
- does **not** create a drop in \(\hat S(t)\)
- but it reduces the risk set size \(n_j\) at later times

So censoring changes future steps indirectly.

---

## 4. KM curve shape: why it’s a step function

KM curve:
- starts at 1.0
- drops only at event times
- stays flat between events
- shows censoring as tick marks

Sketch:

```
S(t)
1.0 |---------
    |        |
    |        |____
    |             |___
    |
    +-------------------- time
         event times
```

---

## 5. Fully worked manual example (must understand)

Dataset:

| Subject | time | event |
|--------:|-----:|------:|
| 1 | 2 | 1 |
| 2 | 4 | 1 |
| 3 | 5 | 0 |
| 4 | 7 | 1 |
| 5 | 8 | 0 |

Event times are: 2, 4, 7

### Step at time 2
Risk set just before 2 includes all 5:

- \(n_1=5\)
- \(d_1=1\)

\[
\hat S(2)=1\cdot\left(1-\frac{1}{5}\right)=0.8
\]

### Step at time 4
After event at 2, 4 subjects remain at risk:

- \(n_2=4\)
- \(d_2=1\)

\[
\hat S(4)=0.8\cdot\left(1-\frac{1}{4}\right)=0.8\cdot0.75=0.6
\]

### Time 5 is censoring
No survival drop at 5.  
But risk set shrinks for later times.

### Step at time 7
At time 7, two subjects are at risk (the one censored at 5 is removed):

- \(n_3=2\)
- \(d_3=1\)

\[
\hat S(7)=0.6\cdot\left(1-\frac{1}{2}\right)=0.3
\]

KM estimates:

| time | \(\hat S(t)\) |
|-----:|--------------:|
| 0 | 1.0 |
| 2 | 0.8 |
| 4 | 0.6 |
| 7 | 0.3 |

---

## 6. Implementing KM by hand (Python + R)

### 6.1 Python: manual KM calculation

!!! interactive "Python"
    ```python
    import pandas as pd

    df = pd.DataFrame({
        "time":[2,4,5,7,8],
        "event":[1,1,0,1,0]
    }).sort_values("time")

    survival = 1.0
    print("time | at_risk | events | S(t)")
    print("------------------------------")

    for t in df["time"].unique():
        at_risk = (df["time"] >= t).sum()
        events  = ((df["time"] == t) & (df["event"] == 1)).sum()

        if events > 0:
            survival *= (1 - events/at_risk)

        print(f"{t:>4} | {at_risk:>7} | {events:>6} | {survival:.3f}")
    ```

This shows:
- risk set size
- event count
- updated KM survival

---

### 6.2 R: manual KM calculation

!!! interactive "R"
    ```r
    df <- data.frame(
      time = c(2,4,5,7,8),
      event = c(1,1,0,1,0)
    )

    df <- df[order(df$time), ]

    survival <- 1

    cat("time | at_risk | events | S(t)\n")
    cat("------------------------------\n")

    for (t in unique(df$time)) {
      at_risk <- sum(df$time >= t)
      events  <- sum(df$time == t & df$event == 1)

      if (events > 0) {
        survival <- survival * (1 - events/at_risk)
      }

      cat(sprintf("%4d | %7d | %6d | %.3f\n", t, at_risk, events, survival))
    }
    ```

---

## 7. KM using standard survival libraries (Python + R)

### 7.1 Python: lifelines KaplanMeierFitter

!!! interactive "Python"
    ```python
    import numpy as np
    import matplotlib.pyplot as plt
    from lifelines import KaplanMeierFitter

    np.random.seed(1)

    n = 200
    T = np.random.exponential(10, n)
    C = np.random.uniform(2, 15, n)

    time = np.minimum(T, C)
    event = (T <= C).astype(int)

    km = KaplanMeierFitter()
    km.fit(time, event, label="KM estimate")

    km.plot()
    plt.title("Kaplan–Meier Survival Curve")
    plt.xlabel("Time")
    plt.ylabel("S(t)")
    plt.show()
    ```

---

### 7.2 R: survival package (survfit)

!!! interactive "R"
    ```r
    set.seed(1)

    n <- 200
    T <- rexp(n, rate = 1/10)
    C <- runif(n, min=2, max=15)

    time <- pmin(T, C)
    event <- as.integer(T <= C)

    library(survival)

    fit <- survfit(Surv(time, event) ~ 1)

    plot(fit, xlab="Time", ylab="S(t)", main="Kaplan–Meier Survival Curve")
    ```

---

## 8. KM curves for groups (treatment vs control)

### 8.1 Python: group KM curves

!!! interactive "Python"
    ```python
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from lifelines import KaplanMeierFitter

    np.random.seed(2)

    n = 300
    group = np.random.binomial(1, 0.5, n)  # 0=control,1=treatment

    # Treatment has longer survival
    T0 = np.random.exponential(9, (group==0).sum())
    T1 = np.random.exponential(13, (group==1).sum())
    T = np.concatenate([T0, T1])

    C = np.random.uniform(2, 16, n)
    time = np.minimum(T, C)
    event = (T <= C).astype(int)

    km0 = KaplanMeierFitter()
    km1 = KaplanMeierFitter()

    ax = plt.subplot(111)

    km0.fit(time[group==0], event[group==0], label="Control").plot(ax=ax)
    km1.fit(time[group==1], event[group==1], label="Treatment").plot(ax=ax)

    plt.title("KM Curves by Group")
    plt.xlabel("Time")
    plt.ylabel("S(t)")
    plt.show()
    ```

---

### 8.2 R: group KM curves

!!! interactive "R"
    ```r
    set.seed(2)

    n <- 300
    group <- rbinom(n, 1, 0.5)

    T0 <- rexp(sum(group==0), rate = 1/9)
    T1 <- rexp(sum(group==1), rate = 1/13)
    T <- c(T0, T1)

    C <- runif(n, 2, 16)
    time <- pmin(T, C)
    event <- as.integer(T <= C)

    library(survival)

    fit <- survfit(Surv(time, event) ~ group)

    plot(fit, col=c("black","red"), lty=1,
         xlab="Time", ylab="S(t)",
         main="KM Curves by Group")
    legend("topright", legend=c("Control","Treatment"),
           col=c("black","red"), lty=1)
    ```

---

## 9. What KM can and cannot do

### KM CAN:
estimate survival probabilities over time  
handle right censoring  
visualize group differences  
estimate median survival

### KM CANNOT:
adjust for covariates (use Cox)  
handle competing risks correctly (use CIF)  
model time-varying covariates directly  
provide causal inference alone

---

## 10. Key takeaways

- KM estimates \(\hat S(t)\) as a product of conditional survival probabilities.
- KM curve drops only at event times.
- Censoring reduces future risk sets but does not create drops.
- KM is foundational to survival analysis and is widely used in biostat papers.

---

## 11. Exercises

<details>
<summary>Click to try</summary>

1. Manually compute KM for a dataset of 6 subjects with 2 censored values.  
2. Explain in plain language why KM multiplies conditional probabilities.  
3. Simulate heavy censoring and see how KM curve becomes less stable late.  
4. Plot KM curves for two groups and interpret differences clinically.  
5. What is the difference between “survival probability” and “hazard”?

</details>
