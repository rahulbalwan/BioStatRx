# 02 — Randomization & Allocation Concealment

Randomization is the design feature that makes a clinical trial fundamentally different from an observational study.  
It protects the treatment comparison from confounding and forms the basis for valid inference.

However, randomization only works as intended when it is paired with **allocation concealment** (to prevent selection bias) and, when feasible, **blinding** (to reduce performance and assessment bias).

This chapter explains the major randomization approaches used in clinical trials, the logic behind each, and how to implement them in both **Python** and **R**.

---

## 1. Why randomize?

### 1.1 Confounding in observational treatment comparisons
In observational studies, treatment is often related to prognosis. For example, sicker patients might be more likely to receive an aggressive therapy. If we compare treated vs untreated outcomes, the comparison mixes:
- the effect of treatment
- the baseline differences in patients

This is the classic confounding problem.

### 1.2 What randomization accomplishes
Randomization makes treatment assignment independent of baseline covariates *in expectation*.

This implies:
- groups are comparable on both measured and unmeasured prognostic factors on average
- differences in outcomes can be attributed to treatment (under the trial’s conditions)
- standard statistical inference (confidence intervals, p-values) is justified by the design

### 1.3 “In expectation” vs “in any single trial”
Randomization guarantees balance **on average** across repeated trials. In any single realized trial, chance imbalances can occur, especially when sample size is small. This is normal and not evidence of trial failure.

---

## 2. Allocation concealment vs blinding

These two concepts are often confused. They address different biases.

### 2.1 Allocation concealment
Allocation concealment means:
- the person enrolling participants cannot predict or influence the next assignment

It happens at the time of enrollment, before treatment is assigned.

If allocation is not concealed, investigators may unconsciously alter enrollment behavior:
- enrolling high-risk patients when they expect the better arm
- delaying enrollment when they expect the worse arm

This breaks comparability of groups.

Common concealment methods:
- centralized web/telephone randomization service
- pharmacy-controlled dispensing
- sequentially numbered, opaque, sealed envelopes (SNOSE) when done correctly

### 2.2 Blinding
Blinding means:
- participants, clinicians, and/or outcome assessors do not know treatment assignment

It happens after assignment.

Blinding reduces bias from:
- differential co-interventions or care
- placebo and behavioral effects
- biased outcome assessment (especially for subjective outcomes)

You can have:
- concealed allocation without blinding
- blinding without proper concealment (rare but possible)
- both (ideal in many drug trials)

---

## 3. Randomization methods

### 3.1 Simple randomization
Each participant is assigned independently using a fixed probability.

For 1:1 allocation:
\[
P(\text{Treatment}) = 0.5,\quad P(\text{Control}) = 0.5
\]

Pros:
- simplest
- strong unpredictability (good for concealment)

Cons:
- imbalance can occur in small trials (e.g., 14 vs 6 in a trial of n=20)
- imbalance within sites (multicenter trials) is possible

Best used in:
- large trials where imbalance is unlikely to be problematic
- settings where balance at interim time points is not crucial

---

### 3.2 Block randomization
Block randomization forces balance within small blocks.

Example: block size 4, 1:1 allocation  
Each block contains 2 treatment and 2 control assignments, in random order.

Pros:
- ensures near-perfect balance throughout enrollment
- useful in small/medium sample trials
- helpful when enrollment may stop early or when interim analyses occur

Cons:
- can become predictable if block size is fixed and known

Mitigation:
- use random block sizes (e.g., 4 and 6)
- keep block structure hidden (allocation concealment)

---

### 3.3 Unequal allocation ratios (e.g., 2:1)
Sometimes trials allocate more patients to treatment:
- more safety data on new therapy
- recruitment incentives (higher chance to get active drug)
- cost or ethical constraints

If ratio is Treatment:Control = 2:1, then within a block size 6:
- 4 treatment
- 2 control

Trade-off:
- power is maximized near 1:1 for fixed total n
- unequal allocation increases required sample size for same power

---

### 3.4 Stratified randomization
Stratified randomization aims to balance groups within levels of important variables (strata), typically using block randomization inside each stratum.

Common stratification factors:
- center/site in multicenter trials
- disease severity category
- sex
- prior treatment history

Pros:
- improves balance on strong prognostic variables
- useful when sample is moderate and imbalance would harm credibility

Cons:
- too many strata becomes unmanageable
- small strata might still be imbalanced due to sparse data

Practical guidance:
- stratify only on the most prognostic variables
- keep number of strata modest

---

### 3.5 Minimization (covariate-adaptive randomization)
Minimization assigns treatment based on current imbalance across multiple covariates. Many implementations include a random component.

Pros:
- can balance several covariates even with small samples
- useful when many baseline factors are critical

Cons:
- more complex to implement and audit
- requires careful documentation and concealment
- analysis often still uses standard ITT methods, but design must be described clearly

Minimization is often run via specialized randomization systems rather than simple scripts.

---

## 4. Baseline balance: what to do and what not to do

### 4.1 What to do
- Provide a baseline characteristics table 
- Look for clinically meaningful imbalances
- Pre-specify adjustment covariates

### 4.2 What not to do: baseline p-values
Testing baseline differences with p-values is widely discouraged because:
- under randomization, any imbalance is due to chance
- baseline p-values do not diagnose trial quality
- they encourage incorrect thinking about randomization

Instead:
- report descriptive summaries
- optionally report standardized differences 

---

# Part A — Implementation in Python

The goal here is to generate reproducible randomization lists that you can actually use in a study workflow.

We will show:
- simple randomization
- block randomization
- random block sizes
- stratified block randomization
- 2:1 allocation example
- quick checks for balance

---

## 5A. Python: Simple randomization (1:1)

!!! interactive "Python"
    ```python
    import numpy as np
    import pandas as pd

    np.random.seed(123)

    n = 40
    assignments = np.random.choice(["Treatment", "Control"], size=n, p=[0.5, 0.5])

    rand_list = pd.DataFrame({
        "participant_id": range(1, n+1),
        "assignment": assignments
    })

    rand_list["assignment"].value_counts()
    ```

---

## 6A. Python: Block randomization (fixed block size, 1:1)

!!! interactive "Python"
    ```python
    import numpy as np
    import pandas as pd

    def block_randomization(n, block_size=4, seed=42):
        if block_size % 2 != 0:
            raise ValueError("For 1:1 allocation, block_size must be even.")
        rng = np.random.default_rng(seed)

        assignments = []
        while len(assignments) < n:
            block = ["Treatment"]*(block_size//2) + ["Control"]*(block_size//2)
            rng.shuffle(block)
            assignments.extend(block)

        assignments = assignments[:n]
        return pd.DataFrame({"participant_id": range(1, n+1), "assignment": assignments})

    rand_block = block_randomization(n=40, block_size=4, seed=7)
    rand_block.head(12)
    ```

---

## 7A. Python: Random block sizes 

!!! interactive "Python"
    ```python
    import numpy as np
    import pandas as pd

    def random_block_randomization(n, block_sizes=(4,6), seed=123):
        rng = np.random.default_rng(seed)
        assignments = []

        while len(assignments) < n:
            b = int(rng.choice(block_sizes))
            if b % 2 != 0:
                continue

            block = ["Treatment"]*(b//2) + ["Control"]*(b//2)
            rng.shuffle(block)
            assignments.extend(block)

        assignments = assignments[:n]
        return pd.DataFrame({"participant_id": range(1, n+1), "assignment": assignments})

    rand_rb = random_block_randomization(n=40, block_sizes=(4,6), seed=77)
    rand_rb["assignment"].value_counts(), rand_rb.head(15)
    ```

---

## 8A. Python: Block randomization with 2:1 allocation

Treatment:Control = 2:1  
A convenient block size is 6 (4 treatment, 2 control).

!!! interactive "Python"
    ```python
    import numpy as np
    import pandas as pd

    def block_randomization_ratio(n, ratio=(2,1), block_size=6, seed=123):
        if block_size != sum(ratio) * (block_size // sum(ratio)):
            # simple sanity check; user should choose compatible sizes
            pass

        rng = np.random.default_rng(seed)

        t, c = ratio
        # scale ratio to block_size
        k = block_size // (t + c)
        block = ["Treatment"]*(t*k) + ["Control"]*(c*k)

        assignments = []
        while len(assignments) < n:
            b = block.copy()
            rng.shuffle(b)
            assignments.extend(b)

        assignments = assignments[:n]
        return pd.DataFrame({"participant_id": range(1, n+1), "assignment": assignments})

    rand_21 = block_randomization_ratio(n=60, ratio=(2,1), block_size=6, seed=9)
    rand_21["assignment"].value_counts()
    ```

---

## 9A. Python: Stratified block randomization

This is common in multicenter trials:
- stratify by site
- stratify by sex
- within each stratum, block randomize with random block sizes

!!! interactive "Python"
    ```python
    import numpy as np
    import pandas as pd

    np.random.seed(10)

    n = 80
    df = pd.DataFrame({
        "participant_id": range(1, n+1),
        "site": np.random.choice(["SiteA", "SiteB", "SiteC"], size=n),
        "sex": np.random.choice(["F", "M"], size=n)
    })

    def stratified_randomization(df, strata_cols, block_sizes=(4,6), seed=123):
        rng = np.random.default_rng(seed)
        out = df.copy()
        out["assignment"] = None

        for _, idx in out.groupby(strata_cols).groups.items():
            m = len(idx)
            # different seed per stratum
            sub_seed = int(rng.integers(1, 10**9))
            sub = random_block_randomization(m, block_sizes=block_sizes, seed=sub_seed)
            out.loc[list(idx), "assignment"] = sub["assignment"].values

        return out.sort_values("participant_id")

    df_rand = stratified_randomization(df, ["site", "sex"], block_sizes=(4,6), seed=5)
    df_rand.head(12)
    ```

Quick balance check:

!!! interactive "Python"
    ```python
    # Check balance by strata
    balance = df_rand.groupby(["site", "sex", "assignment"]).size().unstack(fill_value=0)
    balance
    ```

---

# Part B — Implementation in R

R is frequently used by trial statisticians for randomization lists and audit-friendly scripts.

---

## 10B. R: Simple randomization (1:1)

!!! interactive "R"
    ```r
    set.seed(123)

    n <- 40
    assignment <- sample(c("Treatment", "Control"), n, replace=TRUE)

    rand_list <- data.frame(
      participant_id = 1:n,
      assignment = assignment
    )

    table(rand_list$assignment)
    ```

---

## 11B. R: Block randomization (fixed block size)

!!! interactive "R"
    ```r
    set.seed(7)

    block_randomization <- function(n, block_size=4) {
      if (block_size %% 2 != 0) stop("For 1:1 allocation, block_size must be even.")

      assignments <- c()

      while (length(assignments) < n) {
        block <- rep(c("Treatment", "Control"), each=block_size/2)
        block <- sample(block, length(block))
        assignments <- c(assignments, block)
      }

      assignments <- assignments[1:n]
      data.frame(participant_id=1:n, assignment=assignments)
    }

    rand_block <- block_randomization(40, block_size=4)
    head(rand_block, 12)
    ```

---

## 12B. R: Random block sizes

!!! interactive "R"
    ```r
    set.seed(77)

    random_block_randomization <- function(n, block_sizes=c(4,6)) {
      assignments <- c()

      while (length(assignments) < n) {
        b <- sample(block_sizes, 1)
        if (b %% 2 != 0) next

        block <- rep(c("Treatment", "Control"), each=b/2)
        block <- sample(block, length(block))
        assignments <- c(assignments, block)
      }

      assignments <- assignments[1:n]
      data.frame(participant_id=1:n, assignment=assignments)
    }

    rand_rb <- random_block_randomization(40, block_sizes=c(4,6))
    table(rand_rb$assignment)
    ```

---

## 13B. R: 2:1 allocation with blocks

For ratio 2:1, a block size of 6 gives:
- 4 Treatment
- 2 Control

!!! interactive "R"
    ```r
    set.seed(9)

    block_randomization_ratio <- function(n, ratio=c(2,1), block_size=6) {
      t <- ratio[1]
      c <- ratio[2]
      k <- block_size / (t + c)
      if (k != floor(k)) stop("Choose a block_size divisible by sum(ratio).")

      block <- c(rep("Treatment", t*k), rep("Control", c*k))

      assignments <- c()
      while (length(assignments) < n) {
        b <- sample(block, length(block))
        assignments <- c(assignments, b)
      }

      assignments <- assignments[1:n]
      data.frame(participant_id=1:n, assignment=assignments)
    }

    rand_21 <- block_randomization_ratio(60, ratio=c(2,1), block_size=6)
    table(rand_21$assignment)
    ```

---

## 14B. R: Stratified block randomization (by site and sex)

!!! interactive "R"
    ```r
    set.seed(10)

    n <- 80
    df <- data.frame(
      participant_id = 1:n,
      site = sample(c("SiteA","SiteB","SiteC"), n, replace=TRUE),
      sex = sample(c("F","M"), n, replace=TRUE)
    )

    stratified_randomization <- function(df, strata_cols=c("site","sex"), block_sizes=c(4,6)) {
      out <- df
      out$assignment <- NA

      strata_key <- interaction(out[, strata_cols], drop=TRUE)

      for (s in levels(strata_key)) {
        idx <- which(strata_key == s)
        m <- length(idx)

        a <- random_block_randomization(m, block_sizes=block_sizes)$assignment
        out$assignment[idx] <- a
      }

      out[order(out$participant_id), ]
    }

    df_rand <- stratified_randomization(df, strata_cols=c("site","sex"), block_sizes=c(4,6))
    head(df_rand, 12)
    ```

Balance check:

!!! interactive "R"
    ```r
    with(df_rand, table(site, sex, assignment))
    ```

---

## 15. Practical notes for real trial implementation

### 15.1 Don’t treat a simple script as the full randomization system
In real regulated trials:
- randomization must be auditable
- list generation must be documented
- access to the list must be restricted
- concealment must be guaranteed operationally

A script is typically used to generate a list that is then implemented via:
- IRT/EDC system
- pharmacy
- central randomization service

### 15.2 Document decisions clearly
Include in the protocol/SAP:
- allocation ratio
- method (simple/block/stratified)
- block size policy (random sizes recommended)
- stratification factors
- who generates the list and who has access

---

## 16. Exercises

<details>
<summary>Click to try the exercises</summary>

1. Use Python or R to generate a randomization list for n=120 with 1:1 allocation using random block sizes (4 and 6).  
2. Repeat the generation 100 times and record how often Treatment count differs from Control by more than 10 participants (simple randomization vs block randomization).  
3. Create stratified randomization by site only and compare within-site imbalance to non-stratified randomization.  
4. Write a short paragraph describing your randomization and concealment approach as if for a protocol.  
5. Explain why fixed block size can create predictability and how random block sizes help.

</details>

---

## 17. Summary

Randomization protects treatment comparisons from confounding, but proper implementation matters:
- simple randomization is unpredictable but can be imbalanced in small trials
- block randomization improves balance but can be predictable if poorly implemented
- stratified randomization improves balance on important prognostic factors
- allocation concealment prevents selection bias and is essential even if blinding is not possible
- baseline p-values are discouraged; use descriptive baseline tables instead

---
