# 03 — Trial Designs (Parallel, Cross-over, Cluster, Factorial, Adaptive)

A trial design is the blueprint that determines:
- how participants are assigned to interventions
- what comparisons are made
- how outcomes are measured
- what sources of bias and variability must be controlled

In biostatistics, the design choice affects:
- what causal question you can answer
- sample size and power
- the correct analysis method
- how generalizable the results are

This chapter covers the major clinical trial designs used in biomedical research, with clear guidance on when to use each design and how to analyze them in practice.

---

## 1. Parallel-group randomized controlled trials

### 1.1 What it is
Participants are randomized to one of two or more groups and followed forward.

Example:
- Treatment vs Control
- Dose A vs Dose B vs Placebo

Each participant receives only one assigned intervention.

### 1.2 Strengths
- simple interpretation
- minimal assumptions relative to other designs
- widely accepted by regulators
- works for many endpoint types: continuous, binary, time-to-event

### 1.3 Weaknesses
- between-subject variability can be large
- may require larger sample size than within-subject designs

### 1.4 Typical analysis
- continuous outcome: two-sample t-test / linear regression
- binary outcome: chi-square / logistic regression
- time-to-event: log-rank / Cox regression

---

## 2. Cross-over trials

Cross-over designs are common when:
- the condition is stable
- the treatment effect is reversible
- outcomes return to baseline after washout

### 2.1 Basic 2×2 cross-over design
Each participant receives both treatments in different periods:

- Sequence 1: A → B
- Sequence 2: B → A

Often includes a washout period between treatments.

### 2.2 Why cross-over is efficient
Because each participant serves as their own control, variability is reduced.

A simple way to think:
- parallel design compares person-to-person differences
- cross-over compares within-person differences

This often increases power for the same \(n\).

### 2.3 Key assumptions and problems
- **Carryover effect:** treatment effect persists into the next period
- **Period effect:** outcome changes across time regardless of treatment
- **Washout adequacy:** must remove lingering effects
- **Dropout:** can bias if not random

Cross-over designs are not suitable for:
- curative treatments
- outcomes that permanently change after first treatment
- progressive diseases where baseline changes over time

### 2.4 Typical analysis (continuous outcome)
A standard model includes treatment and period effects, and accounts for repeated measures within subject:

\[
Y_{ij} = \beta_0 + \beta_T T_{ij} + \beta_P P_{ij} + b_i + \varepsilon_{ij},
\quad b_i \sim \mathcal{N}(0,\sigma_b^2),\ \varepsilon_{ij}\sim \mathcal{N}(0,\sigma^2)
\]

- \(i\) indexes subject, \(j\) indexes period
- \(T_{ij}\) indicates treatment received in that period
- \(P_{ij}\) indicates period (1 or 2)
- \(b_i\) is a random subject effect

(Sequence effects may be explored, but the above captures the core structure.)

---

## 3. Cluster randomized trials (CRT)

In a cluster randomized trial, the unit of randomization is a group (cluster), not an individual.

Examples of clusters:
- hospitals
- schools
- villages
- clinics

### 3.1 Why use cluster randomization?
- avoids contamination (e.g., education intervention would spill over)
- intervention naturally delivered at cluster level (policy, training)
- logistical reasons

### 3.2 The key statistical issue: intracluster correlation (ICC)
Within a cluster, outcomes tend to be correlated.

A common ICC definition (continuous outcomes):

\[
\rho = \frac{\sigma^2_{\text{between}}}{\sigma^2_{\text{between}} + \sigma^2_{\text{within}}}
\]

Even small ICC can inflate required sample size.

### 3.3 Design effect (effective sample size reduction)
For equal cluster size \(m\):

\[
\text{Design Effect} = 1 + (m-1)\rho
\]

Effective sample size is approximately:

\[
n_{\text{eff}} \approx \frac{n}{\text{Design Effect}}
\]

### 3.4 Typical analysis
You must account for clustering:
- mixed effects models (random intercept for cluster)
- GEEs with cluster-robust SE
- cluster-level summary analyses (simple but less efficient)

---

## 4. Factorial trials (testing two interventions at once)

A common design is a 2×2 factorial trial:

- Factor A: Drug vs No drug
- Factor B: Exercise vs No exercise

Groups:
1. \(A_0B_0\) (control)
2. \(A_1B_0\) (drug only)
3. \(A_0B_1\) (exercise only)
4. \(A_1B_1\) (both)

### 4.1 When factorial works best
Factorial designs are most efficient when the two treatments:
- act through different mechanisms
- do not interact strongly

### 4.2 Main effects vs interaction
The key question is whether there is interaction:

\[
Y = \beta_0 + \beta_A A + \beta_B B + \beta_{AB}(A\times B) + \varepsilon
\]

- \(\beta_A\) = main effect of A (averaged over B)
- \(\beta_B\) = main effect of B
- \(\beta_{AB}\) = interaction

If interaction is large, interpretation becomes more complex and power for main effects can drop.

---

## 5. Adaptive trial designs

Adaptive designs allow modifications based on interim data while controlling Type I error.

Examples:
- group sequential designs (early stopping)
- sample size re-estimation
- adaptive randomization (more weight to better arm)
- multi-arm multi-stage (MAMS) designs

Adaptive designs require:
- careful pre-specification
- an independent data monitoring committee
- appropriate statistical control

We keep this section conceptual; later pages can go deeper.

---

# Part A — Simulations and analysis in Python

---

## 6A. Python: Parallel vs cross-over efficiency simulation

We simulate:
- same treatment effect
- compare standard error of estimated effect under two designs

!!! interactive "Python"
    ```python
    import numpy as np

    np.random.seed(202)

    def simulate_parallel(n=60, effect=-1.0, sd=2.0):
        trt = np.random.binomial(1, 0.5, n)
        y = effect*trt + np.random.normal(0, sd, n)
        est = y[trt==1].mean() - y[trt==0].mean()
        return est

    def simulate_crossover(n=30, effect=-1.0, sd_within=1.0):
        # Each subject has a baseline level; both periods measured.
        subj = np.random.normal(0, 1.0, n)

        # Outcomes under A (treatment) and B (control)
        y_A = subj + effect + np.random.normal(0, sd_within, n)
        y_B = subj + 0.0    + np.random.normal(0, sd_within, n)

        # Paired difference estimates treatment effect
        est = (y_A - y_B).mean()
        return est

    reps = 2000
    par_est = np.array([simulate_parallel() for _ in range(reps)])
    cro_est = np.array([simulate_crossover() for _ in range(reps)])

    par_est.std(), cro_est.std()
    ```

Interpretation:
- Smaller SD of estimates means greater efficiency.

---

## 7A. Python: Cluster design effect calculation

!!! interactive "Python"
    ```python
    def design_effect(m, icc):
        return 1 + (m - 1) * icc

    m = 25      # cluster size
    icc = 0.02  # intracluster correlation
    de = design_effect(m, icc)
    de
    ```

If your nominal sample size is \(n\) individuals:
- effective \(n\) is about \(n/\text{DE}\)

---

## 8A. Python: Factorial trial simulation with interaction

!!! interactive "Python"
    ```python
    import numpy as np
    import pandas as pd
    import statsmodels.formula.api as smf

    np.random.seed(777)

    n = 800
    A = np.random.binomial(1, 0.5, n)
    B = np.random.binomial(1, 0.5, n)

    # true model
    beta0 = 0
    betaA = -0.6
    betaB = -0.4
    betaAB = -0.8  # interaction (synergy)

    y = beta0 + betaA*A + betaB*B + betaAB*(A*B) + np.random.normal(0, 1, n)

    df = pd.DataFrame({"y": y, "A": A, "B": B})

    fit = smf.ols("y ~ A + B + A:B", data=df).fit()
    fit.summary().tables[1]
    ```

---

# Part B — Simulations and analysis in R

---

## 9B. R: Parallel vs cross-over efficiency simulation

!!! interactive "R"
    ```r
    set.seed(202)

    simulate_parallel <- function(n=60, effect=-1, sd=2) {
      trt <- rbinom(n, 1, 0.5)
      y <- effect*trt + rnorm(n, 0, sd)
      mean(y[trt==1]) - mean(y[trt==0])
    }

    simulate_crossover <- function(n=30, effect=-1, sd_within=1) {
      subj <- rnorm(n, 0, 1)
      yA <- subj + effect + rnorm(n, 0, sd_within)
      yB <- subj + rnorm(n, 0, sd_within)
      mean(yA - yB)
    }

    reps <- 2000
    par_est <- replicate(reps, simulate_parallel())
    cro_est <- replicate(reps, simulate_crossover())

    sd(par_est)
    sd(cro_est)
    ```

---

## 10B. R: Cluster design effect calculation

!!! interactive "R"
    ```r
    design_effect <- function(m, icc) {
      1 + (m-1)*icc
    }

    m <- 25
    icc <- 0.02
    design_effect(m, icc)
    ```

---

## 11B. R: Factorial trial simulation + regression

!!! interactive "R"
    ```r
    set.seed(777)

    n <- 800
    A <- rbinom(n, 1, 0.5)
    B <- rbinom(n, 1, 0.5)

    beta0 <- 0
    betaA <- -0.6
    betaB <- -0.4
    betaAB <- -0.8

    y <- beta0 + betaA*A + betaB*B + betaAB*(A*B) + rnorm(n, 0, 1)

    df <- data.frame(y=y, A=A, B=B)

    fit <- lm(y ~ A + B + A:B, data=df)
    summary(fit)
    ```

---

## 12. Design choice guide

### Parallel-group
Choose when:
- outcomes irreversible or long-lasting
- disease progresses over time
- intervention has long-term effects

### Cross-over
Choose when:
- condition stable
- treatment effect reversible
- outcome quickly measurable  
Avoid when:
- carryover likely
- disease progressive
- long washout needed

### Cluster randomized
Choose when:
- contamination is likely
- intervention is delivered at group level  
Remember:
- account for ICC in sample size and analysis

### Factorial
Choose when:
- testing two interventions simultaneously
- interaction expected to be small/moderate  
Plan:
- pre-specify whether interaction will be tested formally

### Adaptive
Choose when:
- early stopping or efficiency is needed  
Requires:
- strong planning, governance, and pre-specified rules

---

## 13. Exercises

<details>
<summary>Click to try</summary>

1. Modify the cross-over simulation so that carryover exists (e.g., add +0.4 to period 2 outcomes in one sequence) and see how it biases estimates if ignored.  
2. Compute design effect for \(m=10, 25, 50\) across ICC \(=0.01, 0.05, 0.10\). Make a small table.  
3. In factorial simulation, set \(\beta_{AB}=0\) and compare how the interaction term behaves across repeated simulations.  
4. Write a short paragraph: which design would you use for a vaccine education program delivered at school level, and why?  
5. For a chronic stable condition with fast reversible effect, compare sample size needs of parallel vs cross-over conceptually.

</details>

---

## 14. Summary

- Parallel-group trials are the default and simplest to interpret.
- Cross-over trials gain efficiency but require strong assumptions (no carryover, stable condition).
- Cluster randomized trials need ICC-aware design and analysis.
- Factorial designs test multiple interventions efficiently but must consider interaction.
- Adaptive designs improve efficiency but require strict pre-planning and oversight.
