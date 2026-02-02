# 05 — Sample Size & Power (Continuous, Binary, Time-to-Event)

A clinical trial that is too small may fail to detect a clinically meaningful effect.  
A trial that is too large can waste resources or expose unnecessary participants to risk.

Sample size planning is therefore both a scientific and ethical necessity.

This chapter focuses on practical trial power calculations for common endpoint types:

- Continuous outcomes (difference in means)
- Binary outcomes (difference in proportions / risk)
- Time-to-event outcomes (log-rank / Cox; events drive power)
- Adjustments for dropout and non-adherence
- Simulation-based power (often the most intuitive)

We provide implementations in **R** and **Python**.

---

## 1. Core concepts: Type I error, power, effect size

### 1.1 Hypotheses
Most superiority trials test:

\[
H_0: \Delta = 0 \quad \text{vs} \quad H_1: \Delta \ne 0
\]

where \(\Delta\) is the treatment effect (difference in means, difference in proportions, log hazard ratio, etc.).

### 1.2 Type I error (\(\alpha\))
Probability of rejecting \(H_0\) when \(H_0\) is true.

Typical values:
- 0.05 (two-sided)
- 0.025 (one-sided)

### 1.3 Power (1 - \(\beta\))
Probability of detecting an effect if it truly exists.

Common targets:
- 80% power
- 90% power

### 1.4 Effect size: clinically meaningful difference
You do not power a trial for “any difference,” you power it for a difference that matters clinically.

Examples:
- 5 mmHg reduction in systolic BP
- 10% absolute increase in remission
- hazard ratio 0.75 for mortality

---

## 2. General sample size ingredients (always needed)

To plan sample size you must specify:
1) primary endpoint type (continuous/binary/time-to-event)
2) effect size
3) variability (SD, baseline rate, hazard)
4) desired power
5) significance level
6) allocation ratio
7) expected dropout / loss-to-follow-up
8) design features (clustering, repeated measures, interim analyses)

---

## 3. Continuous outcomes (difference in means)

### 3.1 Typical setting
Endpoint: change from baseline or post-treatment measurement.

Common approach:
- two-sample comparison (t-test / ANCOVA-based approximation)

Assume equal SD \(\sigma\), 1:1 allocation, two-sided \(\alpha\).

The standardized effect size is:
\[
d = \frac{\mu_1 - \mu_0}{\sigma}
\]

### 3.2 Key inputs
- \(\Delta = \mu_1 - \mu_0\) (difference in means)
- \(\sigma\) (SD)
- \(\alpha\), power
- allocation ratio

Practical note:
Using ANCOVA (adjusting for baseline) often reduces variance and increases power.  
A common approximation uses an “effective SD” based on baseline correlation.

---

## 4. Binary outcomes (difference in proportions)

### 4.1 Typical endpoints
- remission yes/no
- response yes/no
- event by fixed time yes/no

Inputs:
- control risk \(p_0\)
- treatment risk \(p_1\)
- effect measure:
  - absolute risk difference \(p_1 - p_0\)
  - risk ratio \(p_1/p_0\)
  - odds ratio (less intuitive clinically)

Power formulas typically use:
- difference in proportions

---

## 5. Time-to-event outcomes (events drive power)

For survival endpoints, power depends strongly on the **number of events** rather than the total sample size.

### 5.1 Hazard ratio
Common effect measure:
\[
HR = \frac{h_1(t)}{h_0(t)}
\]

If \(HR < 1\), treatment reduces risk.

### 5.2 Event-based approximation (log-rank)
A classic result (Freedman-type approximation):
\[
D \approx \frac{(z_{1-\alpha/2} + z_{1-\beta})^2}{(\log HR)^2 \, p(1-p)}
\]
where:
- \(D\) = required number of events
- \(p\) = allocation proportion (0.5 for equal groups)

Then total sample size depends on:
- accrual time
- follow-up time
- baseline event rate
- dropout

---

## 6. Adjusting for dropout / loss-to-follow-up

If you expect \(r\) proportion dropout, inflate:

\[
n_{\text{adjusted}} = \frac{n}{1-r}
\]

Example:
- needed \(n=200\)
- expect 15% dropout
- adjusted \(n = 200 / 0.85 = 235.3 \rightarrow 236\)

Dropout may also reduce observed events in survival trials.

---

# Part A — Practical calculations in Python

Python’s `statsmodels` provides power calculators for many common cases.  
For survival, we often use event-based calculations or simulations.

---

## 7A. Python: Continuous endpoint sample size (two-sample t-test)

Example:
- target difference \(\Delta = 5\)
- SD = 12
- alpha = 0.05
- power = 0.80

!!! interactive "Python"
    ```python
    from statsmodels.stats.power import TTestIndPower

    delta = 5
    sd = 12
    effect_size = delta / sd

    analysis = TTestIndPower()
    n_per_group = analysis.solve_power(effect_size=effect_size, alpha=0.05, power=0.80, ratio=1.0, alternative='two-sided')
    n_per_group
    ```

Total n is about \(2 \times n_{\text{per group}}\).

---

## 8A. Python: Binary endpoint sample size (two proportions)

Example:
- control response = 0.40
- treatment response = 0.55
- alpha = 0.05, power = 0.80

!!! interactive "Python"
    ```python
    from statsmodels.stats.power import NormalIndPower
    from statsmodels.stats.proportion import proportion_effectsize

    p0 = 0.40
    p1 = 0.55
    es = proportion_effectsize(p1, p0)  # Cohen's h

    analysis = NormalIndPower()
    n_per_group = analysis.solve_power(effect_size=es, alpha=0.05, power=0.80, ratio=1.0, alternative='two-sided')
    n_per_group
    ```

---

## 9A. Python: Dropout inflation helper

!!! interactive "Python"
    ```python
    import math

    def inflate_dropout(n, dropout_rate):
        return math.ceil(n / (1 - dropout_rate))

    inflate_dropout(200, 0.15)
    ```

---

## 10A. Python: Event-based survival calculation (simple approximation)

Example:
- HR = 0.75
- alpha = 0.05
- power = 0.80
- equal allocation

!!! interactive "Python"
    ```python
    import numpy as np
    from scipy.stats import norm

    alpha = 0.05
    power = 0.80
    HR = 0.75
    p = 0.5

    z_alpha = norm.ppf(1 - alpha/2)
    z_beta = norm.ppf(power)

    D = ((z_alpha + z_beta)**2) / ((np.log(HR))**2 * p * (1-p))
    D
    ```

This returns required number of events \(D\).  
To convert events into sample size, you need an event probability during follow-up (or simulate).

---

## 11A. Python: Simulation-based power (continuous endpoint)

Simulation is often the clearest way to understand power.

Scenario:
- n per group = 60
- SD = 12
- true difference = 5
- test: two-sample t-test
- power estimated by repeated simulation

!!! interactive "Python"
    ```python
    import numpy as np
    from scipy.stats import ttest_ind

    np.random.seed(1)

    def sim_power_continuous(n_per_group=60, delta=5, sd=12, reps=3000, alpha=0.05):
        pvals = []
        for _ in range(reps):
            y0 = np.random.normal(0, sd, n_per_group)
            y1 = np.random.normal(delta, sd, n_per_group)
            pvals.append(ttest_ind(y1, y0, equal_var=True).pvalue)
        return np.mean(np.array(pvals) < alpha)

    sim_power_continuous(n_per_group=60, delta=5, sd=12, reps=2000)
    ```

Now you can vary n to see how power changes.

---

# Part B — Practical calculations in R

R is widely used for trial design because it has well-developed power and survival tools.

---

## 12B. R: Continuous endpoint sample size

!!! interactive "R"
    ```r
    # Example: delta = 5, sd = 12, power = 0.80, alpha = 0.05
    delta <- 5
    sd <- 12

    power.t.test(delta = delta, sd = sd, power = 0.80, sig.level = 0.05, type = "two.sample", alternative = "two.sided")
    ```

---

## 13B. R: Binary endpoint sample size

Example:
- p0 = 0.40
- p1 = 0.55

!!! interactive "R"
    ```r
    p0 <- 0.40
    p1 <- 0.55

    power.prop.test(p1 = p0, p2 = p1, power = 0.80, sig.level = 0.05, alternative = "two.sided")
    ```

---

## 14B. R: Dropout inflation

!!! interactive "R"
    ```r
    inflate_dropout <- function(n, dropout_rate) ceiling(n / (1 - dropout_rate))

    inflate_dropout(200, 0.15)
    ```

---

## 15B. R: Event-based survival calculation (simple approximation)

!!! interactive "R"
    ```r
    alpha <- 0.05
    power <- 0.80
    HR <- 0.75
    p <- 0.5

    z_alpha <- qnorm(1 - alpha/2)
    z_beta <- qnorm(power)

    D <- ((z_alpha + z_beta)^2) / ((log(HR))^2 * p * (1-p))
    D
    ```

---

## 16B. R: Simulation-based power (continuous endpoint)

!!! interactive "R"
    ```r
    set.seed(1)

    sim_power_continuous <- function(n_per_group=60, delta=5, sd=12, reps=2000, alpha=0.05) {
      pvals <- numeric(reps)
      for (i in 1:reps) {
        y0 <- rnorm(n_per_group, 0, sd)
        y1 <- rnorm(n_per_group, delta, sd)
        pvals[i] <- t.test(y1, y0, var.equal=TRUE)$p.value
      }
      mean(pvals < alpha)
    }

    sim_power_continuous(n_per_group=60, delta=5, sd=12)
    ```

---

## 17. Practical guidance: choosing inputs responsibly

### 17.1 Effect size
Effect size should be:
- clinically meaningful
- realistic (based on prior studies or pilot data)
- justified in the protocol

### 17.2 Variability / event rate
Estimate from:
- pilot studies
- published literature
- registry data
- internal historical controls (careful)

### 17.3 Sensitivity analyses
Always compute sample size under multiple plausible assumptions:
- higher SD
- lower effect
- higher dropout
- lower event rate

Present as a table (best practice).

---

## 18. Exercises

<details>
<summary>Click to try</summary>

1. Continuous endpoint: compute n per group for delta=3 and SD=12 (alpha=0.05, power=0.80).  
2. Binary endpoint: compute n per group for p0=0.30, p1=0.40.  
3. Survival: compute required events for HR=0.80 and HR=0.70. Compare.  
4. Use simulation to estimate power for n per group = 40, 60, 80 for delta=5, SD=12.  
5. Apply dropout inflation: if you need 160 total participants and expect 20% dropout, how many should you recruit?

</details>

---

## 19. Summary

- Sample size requires specifying effect size, variability/rates, alpha, power, and allocation ratio.
- Continuous and binary endpoints have standard formula-based power methods.
- Time-to-event power is driven by the number of events; sample size depends on event probability during follow-up.
- Dropout inflates required recruitment.
- Simulation-based power provides an intuitive and flexible approach.

---
