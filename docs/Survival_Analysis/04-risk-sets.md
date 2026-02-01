# 04. Risk Sets and Event Times (The Engine Behind KM, Log-rank, and Cox)

If survival analysis had one “engine” that powers everything, it would be:

# Risk Sets

Kaplan–Meier, log-rank tests, and Cox regression all work by repeatedly asking:

> “Among the people still at risk right now, who experienced the event?”

This chapter makes that idea crystal clear and shows how risk sets are computed in both Python and R.

---

## 1. What is a risk set?

### 1.1 Definition
At an event time \(t\), the **risk set** \(R(t)\) is the set of individuals who:

- have not had the event before time \(t\), AND
- are still under observation just before \(t\) (not censored earlier)

So:

\[
R(t) = \{ i : Y_i \ge t \}
\]

(where \(Y_i\) is observed follow-up time)

### 1.2 Plain language
At time \(t\):

> “Who could possibly experience the event at time \(t\)?”

Those people are the risk set.

---

## 2. Events vs censoring: how they affect risk sets

### 2.1 Event at time \(t\)
- Counted as event
- Removed from risk set after that time

### 2.2 Censoring at time \(t\)
- NOT counted as event
- Removed from risk set after that time

Key point:

Censoring does not cause a drop in survival probability,  
but it does shrink the risk set, which affects future estimates.

---

## 3. A fully worked example

Consider 5 patients:

| Patient | Time | Status |
|--------:|-----:|-------:|
| A | 2 | Event |
| B | 4 | Event |
| C | 5 | Censored |
| D | 7 | Event |
| E | 8 | Censored |

Let’s compute the risk sets at each unique time.

---

### 3.1 Timeline view

Legend:
- **X** = event
- **|** = censoring

```
A: 0 -- X
B: 0 ---- X
C: 0 ----- |
D: 0 ------- X
E: 0 -------- |
```

---

### 3.2 Risk sets step-by-step

#### Just before time 2
Everyone is at risk:
\[
R(2) = \{A,B,C,D,E\},\quad n_1=5
\]
Event occurs (A), so after time 2:
\[
\{B,C,D,E\}
\]

#### Just before time 4
\[
R(4)=\{B,C,D,E\},\quad n_2=4
\]
Event occurs (B), so after time 4:
\[
\{C,D,E\}
\]

#### Just before time 5
\[
R(5)=\{C,D,E\},\quad n_3=3
\]
C is censored at 5, so after time 5:
\[
\{D,E\}
\]

#### Just before time 7
\[
R(7)=\{D,E\},\quad n_4=2
\]
Event occurs (D), so after time 7:
\[
\{E\}
\]

#### Just before time 8
\[
R(8)=\{E\},\quad n_5=1
\]
E is censored, then risk set becomes empty.

---

## 4. Why risk sets matter (big picture)

### 4.1 Kaplan–Meier uses risk sets

Kaplan–Meier estimator:

\[
\hat S(t)=\prod_{t_j\le t}\left(1-\frac{d_j}{n_j}\right)
\]

Where:
- \(n_j\) = size of risk set at event time \(t_j\)
- \(d_j\) = number of events at that time

So the denominator in KM is always:

the risk set at that time, NOT the original sample size.

---

### 4.2 Log-rank test uses risk sets

At each event time, log-rank compares:

- observed events vs expected events
- expected is computed based on risk set proportions

---

### 4.3 Cox regression uses risk sets

The Cox partial likelihood at each event time uses:

\[
\sum_{k \in R(t_j)}\exp(\beta X_k)
\]

It literally sums over the risk set.

So Cox cannot exist without risk sets.

---

## 5. Ties at event times

Sometimes multiple events occur at exactly the same recorded time:

- 3 patients die on day 10

This is called **ties**.

Different tie-handling methods:
- Breslow
- Efron (common default)
- Exact

Software handles this automatically, but it’s good to know why output may differ slightly between packages.

---

## 6. Interactive calculation of risk sets (Python + R)

We’ll compute:

- number at risk at each unique time
- number of events
- number censored

using code.

---

### 6.1 Python: risk set table

!!! interactive "Python"
    ```python
    import pandas as pd

    df = pd.DataFrame({
        "id": ["A","B","C","D","E"],
        "time": [2,4,5,7,8],
        "event": [1,1,0,1,0]
    }).sort_values("time")

    unique_times = df["time"].unique()

    print("time | at_risk | events | censored")
    print("-----------------------------------")

    for t in unique_times:
        at_risk = (df["time"] >= t).sum()
        events = ((df["time"] == t) & (df["event"] == 1)).sum()
        cens   = ((df["time"] == t) & (df["event"] == 0)).sum()
        print(f"{t:>4} | {at_risk:>7} | {events:>6} | {cens:>8}")
    ```

**Interpretation**
- `at_risk` is the size of the risk set just before time `t`
- `events` is number of events at time `t`
- `censored` is number censored at time `t`

---

### 6.2 R: risk set table

!!! interactive "R"
    ```r
    df <- data.frame(
      id = c("A","B","C","D","E"),
      time = c(2,4,5,7,8),
      event = c(1,1,0,1,0)
    )

    df <- df[order(df$time), ]

    unique_times <- unique(df$time)

    cat("time | at_risk | events | censored\n")
    cat("-----------------------------------\n")

    for (t in unique_times) {
      at_risk <- sum(df$time >= t)
      events  <- sum(df$time == t & df$event == 1)
      cens    <- sum(df$time == t & df$event == 0)
      cat(sprintf("%4d | %7d | %6d | %8d\n", t, at_risk, events, cens))
    }
    ```

---

## 7. Visualizing risk set size over time (Python + R)

Risk sets shrink over time; this helps build intuition.

### 7.1 Python: step plot of risk set size

!!! interactive "Python"
    ```python
    import numpy as np
    import matplotlib.pyplot as plt

    times = df["time"].values
    risk_sizes = [(df["time"] >= t).sum() for t in times]

    plt.step(times, risk_sizes, where="post")
    plt.title("Risk Set Size Over Time")
    plt.xlabel("time")
    plt.ylabel("number at risk")
    plt.show()
    ```

---

### 7.2 R: step plot of risk set size

!!! interactive "R"
    ```r
    times <- df$time
    risk_sizes <- sapply(times, function(t) sum(df$time >= t))

    plot(times, risk_sizes, type="s", main="Risk Set Size Over Time",
         xlab="time", ylab="number at risk")
    ```

---

## 8. Clinical intuition

At each event time, survival analysis compares:

> “The person who failed” vs “everyone who could have failed at that time.”

That comparison group is the risk set.

This is why:
- censoring affects denominators
- late estimates are less precise (risk sets small)
- Kaplan–Meier and Cox are fundamentally conditional methods

---

## 9. Common mistakes

### Mistake 1: using original N in denominator
Fix: always use risk set size.

### Mistake 2: leaving censored people in the risk set
Fix: remove them after censoring time.

### Mistake 3: mis-ordering events and censoring at same time
Fix: software uses a convention; usually events are counted before censoring at the same time.

---

## 10. Key takeaways

- Risk set \(R(t)\) = those still event-free and observed at time \(t\).
- KM, log-rank, and Cox all depend on risk sets.
- Censoring does not count as event, but removes people from future risk sets.
- As risk sets shrink, estimates become less stable.

---

## 11. Exercises

<details>
<summary>Click to try</summary>

1. For the sample data, write out \(R(4)\) explicitly.  
2. Explain why censoring affects future KM steps even though it doesn’t cause survival drops.  
3. Create a dataset where two events happen at the same time. How do risk sets work then?  
4. Simulate 20 subjects with random censoring and compute risk sets at each event time.  
5. Explain in plain language why Cox regression is built on risk sets.

</details>
