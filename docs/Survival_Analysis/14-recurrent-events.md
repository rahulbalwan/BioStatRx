# Recurrent Events (Multiple Events per Subject) 

Many biostatistics outcomes are not “time to the first event.”

Instead, subjects can experience the event **multiple times**.

Examples:
- repeated hospitalizations
- recurrent infections (UTI, malaria episodes)
- repeated asthma exacerbations
- recurrent falls in geriatrics
- repeated device malfunctions

These are **recurrent events**.

If you use standard survival analysis on only the first event, you lose information.

This chapter covers:

- what recurrent event data looks like  
- why standard Cox is not enough  
- major recurrent-event approaches  
- Andersen–Gill (AG) model  
- PWP models (total-time and gap-time)  
- robust SE and frailty  
- Python + R implementation  

---

## 1. What is a recurrent event?

A recurrent event means:

> A subject can experience the same event multiple times during follow-up.

We observe multiple event times:

\[
T_{i1}, T_{i2}, \dots
\]

And the process may stop at:
- study end
- dropout
- death (terminal event)

---

## 2. Why “time to first event” can be insufficient

If you analyze only first hospitalization:

- you ignore later hospitalizations (important clinically)
- you underestimate disease burden
- treatment might reduce recurrence even if first event is similar

So recurrent event analysis helps quantify total burden.

---

## 3. Recurrent events vs competing risks vs terminal events

Important distinction:

- Recurrent event: can happen multiple times
- Competing risk: different event prevents event of interest
- Terminal event: death stops recurrence process

In many studies you have both:
- repeated hospitalizations
- followed by death

That is a more advanced joint modeling topic, but we mention it.

---

## 4. Data structure for recurrent events (start–stop)

Recurrent event data typically uses **counting process** format:

| id | start | stop | event | covariates |
|---:|------:|-----:|------:|-----------|
| 1 | 0 | 5 | 1 | ... |
| 1 | 5 | 9 | 1 | ... |
| 1 | 9 | 14 | 0 | ... |

- each row represents a time interval at risk
- event=1 indicates event occurred at end of interval
- after an event, subject re-enters risk set (with new interval)

---

## 5. Key modeling challenge: within-subject dependence

Events for the same subject are correlated:

- some patients are “frequent flyers”
- unmeasured frailty may increase recurrence risk

So we must adjust for correlation:
- robust SE
- frailty
- event-order models

---

## 6. Main recurrent-event model families (biostat overview)

### 6.1 Andersen–Gill (AG) model (most common starter)
- treats recurrent events as a counting process
- subject is at risk again immediately after each event
- assumes same baseline hazard for all events
- uses robust SE (cluster by subject)

Great for:
- overall recurrence rate comparison

Limitations:
- does not explicitly model event order
- assumes same risk mechanism for first and later events

### 6.2 Prentice–Williams–Peterson (PWP) models
These model event order:

- PWP Total-time: time measured from baseline
- PWP Gap-time: time resets after each event

These are good when:
- event 2 depends on having event 1
- baseline hazard differs by event number

### 6.3 Frailty models for recurrence
Random effect per subject:
- frequent-event subjects have higher frailty
- accounts for unobserved heterogeneity

---

## 7. Andersen–Gill model in detail

### 7.1 Model form
\[
h_i(t)=h_0(t)\exp(\beta^TX_i(t))
\]

But each subject can have multiple events.

Key concept: counting process \(N_i(t)\):
- counts how many events by time \(t\)

AG treats each event as contributing to Cox partial likelihood with risk sets.

### 7.2 Dependence handling
Standard errors are usually computed using:

- robust sandwich estimator clustered by subject (id)

In R: `cluster(id)`  
In Python: robust=True can help, but full recurrent support is limited

---

## 8. PWP models (event-order models)

### 8.1 Total-time
Time is measured from baseline.
Event strata created by event number.

### 8.2 Gap-time
Time resets after each event.
Good when you care about time between recurrences.

---

# PART A — R (best tool for recurrent events)

R’s `survival` package is the standard for recurrent event modeling.

We will show:
- how to build recurrent-event dataset
- fit AG model
- fit PWP total-time and gap-time

---

## 9A. Simulate recurrent event data (R)

We simulate:
- 300 subjects
- treatment reduces recurrence rate
- each subject can have up to 3 events
- random censoring

!!! interactive "R"
    ```r
    set.seed(2028)

    n <- 300
    trt <- rbinom(n, 1, 0.5)
    age <- rnorm(n, 60, 10)

    # baseline recurrence rate
    base <- 0.15

    beta_trt <- -0.50
    beta_age <- 0.02

    max_events <- 3

    rows <- list()
    idx <- 1

    for (i in 1:n) {

      # subject-specific frailty-like multiplier (optional)
      u <- rgamma(1, shape=2, scale=0.5)  # mean 1

      t <- 0
      event_count <- 0

      censor <- runif(1, 4, 12)

      while (event_count < max_events) {

        # hazard for next event (exponential gap time for simplicity)
        h <- u * base * exp(beta_trt*trt[i] + beta_age*(age[i]-60))

        gap <- rexp(1, rate=h)
        t_next <- t + gap

        if (t_next <= censor) {
          # event occurs
          rows[[idx]] <- data.frame(
            id=i,
            start=t,
            stop=t_next,
            event=1,
            event_number=event_count+1,
            trt=trt[i],
            age=age[i]
          )
          idx <- idx + 1

          t <- t_next
          event_count <- event_count + 1
        } else {
          # censored before next event
          rows[[idx]] <- data.frame(
            id=i,
            start=t,
            stop=censor,
            event=0,
            event_number=event_count+1,
            trt=trt[i],
            age=age[i]
          )
          idx <- idx + 1
          break
        }
      }
    }

    d <- do.call(rbind, rows)
    head(d, 12)
    ```

---

## 10A. Andersen–Gill model in R

Andersen–Gill is just Cox with start–stop and robust SE clustered by id:

!!! interactive "R"
    ```r
    library(survival)

    fit_ag <- coxph(Surv(start, stop, event) ~ trt + age + cluster(id), data=d)
    summary(fit_ag)
    ```

Interpretation:
- HR for `trt` compares recurrence hazard at any moment between groups.
- Robust SE accounts for within-subject correlation.

---

## 11A. PWP Total-Time model (event strata)

PWP total-time uses strata by event number:
- each event number has its own baseline hazard

!!! interactive "R"
    ```r
    fit_pwp_tt <- coxph(Surv(start, stop, event) ~ trt + age + strata(event_number) + cluster(id), data=d)
    summary(fit_pwp_tt)
    ```

Interpretation:
- HRs compare covariate effects across event strata.
- baseline hazard differs by event number.

---

## 12A. PWP Gap-Time model

For gap-time, we reset time within each event interval:

Define:
- gap_start = 0 for each interval
- gap_stop = stop - start

!!! interactive "R"
    ```r
    d$gap_start <- 0
    d$gap_stop <- d$stop - d$start

    fit_pwp_gt <- coxph(Surv(gap_start, gap_stop, event) ~ trt + age + strata(event_number) + cluster(id), data=d)
    summary(fit_pwp_gt)
    ```

Interpretation:
- models time between events
- baseline hazard differs by event number

---

## 13A. Plotting cumulative mean function
Another descriptive way to visualize recurrence is cumulative events over time.

A simple approach:
- estimate mean cumulative function by group

(Advanced: use `survfit` on counting process, or specialized packages.)

---

# PART B — PYTHON (practical options)

Python currently lacks full native tooling for recurrent events compared to R.
But you can still do:

- create start–stop data  
- use `CoxTimeVaryingFitter` as an AG-style model  
- use robust SE cautiously  
- or do marginal approaches externally  

We show an AG-style approach using lifelines time-varying fitter.

---

## 14B. Simulate recurrent events and build start–stop data (Python)

!!! interactive "Python"
    ```python
    import numpy as np
    import pandas as pd

    np.random.seed(2028)

    n = 250
    trt = np.random.binomial(1, 0.5, n)
    age = np.random.normal(60, 10, n)

    base = 0.15
    beta_trt = -0.50
    beta_age = 0.02
    max_events = 3

    rows = []

    for i in range(n):
        u = np.random.gamma(shape=2.0, scale=0.5)  # mean 1
        t = 0.0
        ev = 0
        censor = np.random.uniform(4, 12)

        while ev < max_events:
            h = u * base * np.exp(beta_trt*trt[i] + beta_age*(age[i]-60))
            gap = np.random.exponential(1/h)
            t_next = t + gap

            if t_next <= censor:
                rows.append([i, t, t_next, 1, ev+1, trt[i], age[i]])
                t = t_next
                ev += 1
            else:
                rows.append([i, t, censor, 0, ev+1, trt[i], age[i]])
                break

    d = pd.DataFrame(rows, columns=["id","start","stop","event","event_number","treatment","age"])
    d.head(12)
    ```

---

## 15B. Andersen–Gill style fit in Python (lifelines)

We fit start–stop using `CoxTimeVaryingFitter`.

!!! interactive "Python"
    ```python
    from lifelines import CoxTimeVaryingFitter

    ctv = CoxTimeVaryingFitter()
    ctv.fit(d, id_col="id", start_col="start", stop_col="stop", event_col="event")
    ctv.print_summary()
    ```

Important note:
- lifelines standard errors for recurrent events may not match R’s robust cluster(id) sandwich perfectly.
- For publication-grade recurrent-event analysis, R is typically preferred.

---

## 16. Which recurrent-event method should you choose?

### 16.1 If you want a simple “overall recurrence rate” effect:
Andersen–Gill with robust SE

### 16.2 If you believe event order matters:
PWP total-time or gap-time

### 16.3 If you believe subjects have different baseline recurrence tendencies:
frailty model (random effect per subject)

Often, analysts do:
- AG as primary
- PWP as sensitivity analysis
- frailty as additional robustness check

---

## 17. Reporting recurrent-event results (biostat style)

Example:

> “Recurrent hospitalizations were analyzed using an Andersen–Gill Cox model with robust sandwich standard errors clustered by subject. Treatment reduced the hospitalization hazard (HR 0.72, 95% CI …).”

For PWP:

> “Event-order was accounted for using PWP models stratified by event number.”

Always report:
- model type (AG / PWP)
- use of robust SE / frailty
- what time scale (total-time vs gap-time)
- interpretation of HR

---

## 18. Common mistakes

### Mistake 1: treat repeated events as independent rows without clustering correction
This underestimates SE.

### Mistake 2: ignore death as terminal event
Death ends recurrence. Treating it as censoring might be informative.

### Mistake 3: interpret AG HR as first-event effect
AG HR is about **overall recurrence hazard**, not first-event survival.

### Mistake 4: too few events per subject but still using complex models
If almost everyone has 0–1 events, use standard Cox.

---

## 19. Key takeaways

- Recurrent events require start–stop counting-process data.
- Within-subject correlation must be handled (robust SE, frailty, event-strata).
- Andersen–Gill is the common baseline recurrent-event model.
- PWP models incorporate event order (total-time / gap-time).
- R has best support; Python can approximate via time-varying Cox.

---

## 20. Exercises

<details>
<summary>Click to try</summary>

1. Simulate recurrent events with strong treatment effect and fit AG model. Interpret HR.  
2. Fit PWP total-time and compare with AG. Does effect change?  
3. Fit PWP gap-time and interpret what changes.  
4. Increase within-subject heterogeneity (frailty variance) and see how robust SE changes inference.  
5. Write a “Methods” paragraph describing recurrent-event analysis of hospitalizations.

</details>
