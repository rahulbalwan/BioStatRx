# 02. Time, Event, and Censoring (Survival Data Foundations)

This chapter builds the **core foundations** of survival analysis:

1. **Time** — what exactly are we measuring?
2. **Event** — what outcome ends follow-up?
3. **Censoring** — why follow-up ends without the event

If these three are not defined clearly, **everything else (KM, Cox, log-rank) becomes wrong or misleading**.

---

## 1. What is “time” in survival analysis?

### 1.1 Definition
In survival analysis, **time** is:

> the duration from a clearly defined starting point (“time zero”) until an event or censoring.

We always measure time relative to a **time origin**.

### 1.2 Common choices of time zero in biostatistics

| Study type | Time zero examples |
|-----------|---------------------|
| Clinical trial | randomization date, treatment start |
| Oncology cohort | diagnosis date, surgery date |
| Infection study | enrollment date, discharge date |
| Hospital study | admission date, discharge date |
| Device study | implant date, installation date |

### 1.3 Important rule: define time zero consistently
If time zero differs across patients, estimates can become biased.

Good: time zero is treatment initiation for everyone  
Bad: some start at diagnosis, others at surgery, others at drug start

---

## 2. What is an event?

### 2.1 Definition
An **event** is the outcome that ends survival follow-up for the endpoint of interest.

It must be **precisely defined**.

### 2.2 Examples of events
- Death (all-cause or disease-specific)
- Relapse / recurrence
- Hospital readmission
- ICU admission
- Infection
- Recovery (yes, “positive outcomes” can be events too)
- Device failure

### 2.3 Event coding
We usually encode:

\[
\delta =
\begin{cases}
1 & \text{event occurred} \\
0 & \text{censored}
\end{cases}
\]

So your dataset has:
- `time`
- `event`

---

## 3. What is censoring?

### 3.1 Core meaning
Censoring happens when:

> we stop observing a subject before the event is observed.

So we know:
- they survived up to last contact  
but we do not know what happens after.

### 3.2 Right censoring (most common)
Right censoring means the true event time is later than what we observe:

\[
T > C
\]

Examples:
- study ends
- patient withdraws
- patient lost to follow-up
- patient still alive at last contact

### 3.3 Censoring is NOT missingness
Censored observations still contribute valuable information:

If censored at time 4 years, we know:

> event did not occur for at least 4 years.

Survival methods are built to use this information correctly.

---

## 4. Types of censoring (biostat perspective)

### 4.1 Right censoring
Most typical in clinical studies.

```
Start ----------- |  (event not observed)
```

### 4.2 Left censoring
Event occurred before observation begins.

Example:
- infection occurred before first clinic visit, but we only detect later.

Less common in standard time-to-event clinical trials.

### 4.3 Interval censoring
Event is known to occur in an interval, but exact time unknown.

Example:
- patient tested negative at month 3
- positive at month 6
→ infection occurred between 3 and 6 months

Interval censoring requires specialized methods.

---

## 5. Non-informative censoring assumption (very important)

Most KM and Cox analyses assume:

> Censoring is independent of the event process, given model covariates.

If patients with severe disease are more likely to drop out, censoring can be informative and bias results.

Practical hints to reduce bias:
- improve follow-up
- record reasons for dropout
- sensitivity analyses

---

## 6. The survival dataset: what we actually observe

Let:
- \(T\) = true event time
- \(C\) = censoring time

We observe:

\[
Y = \min(T, C)
\]
\[
\delta = I(T \le C)
\]

This is the core survival data structure.

---

## 7. Risk set intuition (preview)

At any time \(t\), the **risk set** includes people who:
- have not had event yet
- are still under observation

Censored subjects leave the risk set at censor time.

This concept powers:
- Kaplan–Meier
- log-rank test
- Cox regression

(We cover risk sets fully in a later chapter, but keep this intuition.)

---

## 8. Visual intuition: timelines

Legend:
- **X** = event
- **|** = censored

```
Patient 1: 0 -------- X
Patient 2: 0 -------------- |
Patient 3: 0 ----- X
Patient 4: 0 ---------- |
```

At time = 6:
- patient 1 already had event
- patient 3 already had event
- patient 2 still at risk (if censor > 6)
- patient 4 maybe censored earlier

---

## 9. Interactive simulation: building survival data (Python + R)

This section shows exactly how survival data are formed using:
- event times
- censoring times

### 9.1 Python: simulate event + censoring

!!! interactive "Python"
    ```python
    import numpy as np
    import pandas as pd

    np.random.seed(42)

    n = 50

    # True event times
    T = np.random.exponential(scale=8, size=n)

    # Censoring times
    C = np.random.uniform(2, 12, size=n)

    # Observed
    time = np.minimum(T, C)
    event = (T <= C).astype(int)

    df = pd.DataFrame({
        "true_event_time_T": T,
        "censor_time_C": C,
        "time": time,
        "event": event
    }).sort_values("time")

    df.head(12)
    ```

**Try:**
- Increase censoring: `C = np.random.uniform(2, 6, size=n)`
- Decrease censoring: `C = np.random.uniform(2, 25, size=n)`

Then check censoring proportion:

!!! interactive "Python"
    ```python
    censor_prop = 1 - df["event"].mean()
    censor_prop
    ```

---

### 9.2 R: simulate event + censoring

!!! interactive "R"
    ```r
    set.seed(42)

    n <- 50

    # True event times
    T <- rexp(n, rate = 1/8)

    # Censoring times
    C <- runif(n, min = 2, max = 12)

    time <- pmin(T, C)
    event <- as.integer(T <= C)

    df <- data.frame(
      true_event_time_T = T,
      censor_time_C = C,
      time = time,
      event = event
    )

    df <- df[order(df$time), ]
    head(df, 12)
    ```

Censoring proportion:

!!! interactive "R"
    ```r
    1 - mean(df$event)
    ```

---

## 10. A concrete clinical example

Suppose:
- Study duration = 5 years
- Patient followed 3 years and then lost
- No relapse observed

We record:
- time = 3
- event = 0

We do NOT record:
- event = “no relapse”

Because we **don’t know** after 3 years.

That is the key conceptual difference.

---

## 11. Common mistakes (and fixes)

### Mistake 1: treat censored as “event-free”
Fix: censoring means unknown beyond last time.

### Mistake 2: use unequal time zero
Fix: define time origin consistently.

### Mistake 3: mix different event definitions
Fix: define event precisely and consistently.

### Mistake 4: ignore censoring assumption
Fix: assess dropout reasons, consider sensitivity.

---

## 12. Key takeaways

- Survival data consists of **time** and **event indicator**.
- Censoring means **partial follow-up**, not missingness.
- Right censoring is most common in clinical research.
- We observe \(Y=\min(T,C)\) and \(\delta=I(T\le C)\).
- Correct definitions of time zero and event are critical.

---

## 13. Exercises

<details>
<summary>Click to try</summary>

1. Give 3 different valid choices of time zero in clinical studies and explain which is best in a trial.  
2. Define right censoring in one sentence.  
3. In the simulation, increase censoring and describe how analysis becomes harder.  
4. Give a real example of interval censoring in medicine.  
5. Explain non-informative censoring assumption in plain language.

</details>
