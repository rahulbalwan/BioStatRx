# Analysis Populations (ITT, Per-Protocol, As-Treated) + Missing Data

A trial’s design starts with randomization, but the credibility of results often depends on how we define:

- *Who is included in the analysis?*  
- *How do we handle missing outcomes and protocol deviations?*

This chapter covers the most important analysis populations:
- **Intention-to-treat (ITT)**
- **Per-protocol (PP)**
- **As-treated (AT)**

and provides practical strategies for handling missing data:
- complete-case (when acceptable and when not)
- multiple imputation (MI)
- model-based approaches (likelihood / mixed models concept)

We focus on biostatistical thinking:
- what estimand each approach answers
- what biases can enter
- what to report in a paper/SAP

---

## 1. Why analysis populations matter

Trials rarely run perfectly:
- participants discontinue treatment
- participants cross over or switch therapies
- outcomes are missing due to dropout or missed visits
- protocol deviations occur

How you define the analysis population determines:
- whether randomization is preserved
- what causal question is being answered
- whether results are biased

---

## 2. Intention-to-treat (ITT)

### 2.1 Definition
In the ITT principle, you analyze participants **as randomized**, regardless of:
- adherence
- switching
- dropout (in principle; missingness still needs handling)

Population:
- usually “all randomized participants” (sometimes called “full analysis set”)

### 2.2 Why ITT is the default primary analysis
- preserves the benefits of randomization
- estimates the effect of “treatment policy”:
  - assigning the treatment in practice
  - including real-world non-adherence

Interpretation:
> What happens if we assign patients to this treatment strategy in clinical practice?

### 2.3 ITT does not mean “ignore missing data”
A common misconception:
- ITT is not “analyze only those with observed outcomes”
- missingness must be addressed appropriately

---

## 3. Per-protocol (PP)

### 3.1 Definition
Per-protocol analysis includes only participants who:
- adhered sufficiently to the protocol
- had no major protocol deviations
- met key eligibility criteria
- received adequate exposure

Different trials define PP differently, so PP must be clearly specified.

### 3.2 What PP estimates
PP aims to estimate the effect of treatment **under perfect adherence**.

Interpretation:
> What would the treatment effect be if everyone adhered as intended?

### 3.3 The main danger: selection bias
Once you exclude participants post-randomization based on adherence or deviations, groups may no longer be comparable.

Adherence can relate to prognosis:
- sicker participants may discontinue more
- adverse events may cause discontinuation

So PP can be biased unless handled carefully.

PP is usually a sensitivity analysis, not primary.

---

## 4. As-treated (AT)

### 4.1 Definition
Participants are analyzed according to what they actually received, not what they were randomized to.

This is closer to observational analysis:
- randomization is partially lost if switching occurs

Interpretation:
> Effect of receiving treatment, rather than being assigned treatment.

### 4.2 When AT is used
Often used as:
- secondary analysis
- safety analysis
- exploratory analysis

Main risk:
- confounding due to post-randomization switching

---

## 5. A practical trial workflow: ITT primary + sensitivity analyses

A common robust plan is:

Primary:
- ITT analysis, aligned with the primary estimand

Sensitivity:
- PP analysis (carefully defined)
- hypothetical estimand approaches (if relevant)
- alternative missing data assumptions

This allows readers and regulators to evaluate robustness.

---

## 6. Missing data

Missing outcomes occur due to:
- dropout
- missed visits
- withdrawal
- administrative censoring
- death (sometimes a competing event rather than “missing”)

Before choosing methods, understand *why* data are missing.

---

## 7. Missing data mechanisms
### 7.1 MCAR: Missing Completely At Random
Missingness does not depend on observed or unobserved data.

Example:
- lab machine failure on random days

If MCAR holds, complete-case analysis can be unbiased (but still loses power).

### 7.2 MAR: Missing At Random
Missingness may depend on observed data, but not on unobserved outcomes after conditioning.

Example:
- older participants miss visits more, but age is recorded

Many standard methods (MI, likelihood) assume MAR.

### 7.3 MNAR: Missing Not At Random
Missingness depends on unobserved outcomes even after conditioning on observed data.

Example:
- participants with worsening symptoms are more likely to drop out, and worsening is not fully captured by observed variables

MNAR requires sensitivity analyses; it cannot be “fixed” by standard MI without assumptions.

---

## 8. Common missing-data approaches in trials

### 8.1 Complete-case analysis
Analyze only participants with observed outcomes.

Pros:
- simple

Cons:
- biased if missingness is related to outcomes (common in trials)
- reduces power

Acceptable only under strong conditions (rare).

### 8.2 Single imputation 
Examples:
- last observation carried forward (LOCF)
- baseline carried forward
- mean imputation

Problems:
- can bias treatment effect
- underestimates uncertainty
- LOCF makes strong assumptions about disease trajectory

Usually not recommended as primary for modern trials.

### 8.3 Multiple imputation (MI)
MI fills in missing values multiple times to create several completed datasets.  
Each dataset is analyzed, then results are combined (Rubin’s rules).

Key advantages:
- accounts for uncertainty in imputation
- flexible for many outcome types
- widely accepted if assumptions are justified

### 8.4 Model-based likelihood methods (e.g., mixed models)
For repeated continuous outcomes, mixed models (MMRM) can provide valid inference under MAR without explicit imputation.

These are common in longitudinal clinical trials.

---

# Part A — Practical demonstration in Python

Python has strong imputation tools; full Rubin’s-rule MI for regression is less standard than in R,
but we can still demonstrate correct workflow and make the logic clear.

---

## 9A. Python: Simulate a trial dataset with dropout

We simulate:
- baseline severity
- treatment effect
- dropout that depends on worsening (MAR/MNAR-like depending on variables included)

!!! interactive "Python"
    ```python
    import numpy as np
    import pandas as pd

    np.random.seed(123)

    n = 500
    trt = np.random.binomial(1, 0.5, n)

    severity0 = np.random.normal(0, 1, n)
    # true treatment effect on outcome
    true_delta = -0.5
    noise = np.random.normal(0, 1, n)

    # continuous outcome at follow-up
    y = 1.0*severity0 + true_delta*trt + noise

    # dropout probability depends on outcome (worse outcomes more likely to drop)
    # we don't observe y if dropped -> MNAR-ish in truth
    p_drop = 1 / (1 + np.exp(-0.8*(y - 0.5)))
    drop = np.random.binomial(1, p_drop)

    y_obs = y.copy()
    y_obs[drop==1] = np.nan

    df = pd.DataFrame({"trt": trt, "severity0": severity0, "y_true": y, "y": y_obs, "drop": drop})
    df.head()
    ```

Check missingness rate by group:

!!! interactive "Python"
    ```python
    df.groupby("trt")["y"].apply(lambda s: s.isna().mean())
    ```

---

## 10A. Python: ITT naive complete-case estimate (often biased)

!!! interactive "Python"
    ```python
    cc = df.dropna(subset=["y"])
    est_cc = cc.loc[cc.trt==1, "y"].mean() - cc.loc[cc.trt==0, "y"].mean()
    est_cc
    ```

Compare to the true ITT difference using the latent true outcome:

!!! interactive "Python"
    ```python
    true_itt = df.loc[df.trt==1, "y_true"].mean() - df.loc[df.trt==0, "y_true"].mean()
    true_itt
    ```

If missingness is related to outcome, complete-case can distort the effect.

---

## 11A. Python: Simple imputation (mean imputation) — shown only to demonstrate why it’s weak

!!! interactive "Python"
    ```python
    df_meanimp = df.copy()
    df_meanimp["y_imp"] = df_meanimp["y"].fillna(df_meanimp["y"].mean())

    est_meanimp = df_meanimp.loc[df_meanimp.trt==1, "y_imp"].mean() - df_meanimp.loc[df_meanimp.trt==0, "y_imp"].mean()
    est_meanimp
    ```

Mean imputation can bias effect and understate uncertainty.

---

## 12A. Python: Imputation using observed covariates (MI-style workflow)

We do a principled imputation using baseline severity and treatment as predictors.

This is not full Rubin’s-rule MI, but it shows the essential idea:
- use observed data to predict missing outcomes

!!! interactive "Python"
    ```python
    import statsmodels.api as sm

    df_imp = df.copy()

    # Fit a regression model on observed cases
    obs = df_imp.dropna(subset=["y"])
    X = sm.add_constant(obs[["trt", "severity0"]])
    model = sm.OLS(obs["y"], X).fit()

    # Predict missing y and add random noise based on residual SD (one stochastic imputation)
    miss = df_imp["y"].isna()
    Xmiss = sm.add_constant(df_imp.loc[miss, ["trt", "severity0"]])

    resid_sd = obs["y"].sub(model.predict(X)).std()
    y_pred = model.predict(Xmiss)
    df_imp.loc[miss, "y"] = y_pred + np.random.normal(0, resid_sd, miss.sum())

    # ITT estimate after imputation
    est_imp = df_imp.loc[df_imp.trt==1, "y"].mean() - df_imp.loc[df_imp.trt==0, "y"].mean()
    est_imp
    ```

In real MI:
- repeat this M times
- analyze each dataset
- pool estimates + SE via Rubin’s rules
R makes this easier and is standard for clinical trials.

---

# Part B — Practical multiple imputation in R (recommended for trials)

We show MI using the `mice` package.

---

## 13B. R: Simulate trial with dropout

!!! interactive "R"
    ```r
    set.seed(123)

    n <- 500
    trt <- rbinom(n, 1, 0.5)
    severity0 <- rnorm(n, 0, 1)

    true_delta <- -0.5
    noise <- rnorm(n, 0, 1)

    y_true <- 1.0*severity0 + true_delta*trt + noise

    p_drop <- 1 / (1 + exp(-0.8*(y_true - 0.5)))
    drop <- rbinom(n, 1, p_drop)

    y <- y_true
    y[drop==1] <- NA

    df <- data.frame(trt=trt, severity0=severity0, y=y, drop=drop, y_true=y_true)
    head(df)
    ```

Missingness rate:

!!! interactive "R"
    ```r
    tapply(is.na(df$y), df$trt, mean)
    ```

---

## 14B. R: Complete-case ITT estimate (naive)

!!! interactive "R"
    ```r
    cc <- df[!is.na(df$y), ]
    est_cc <- mean(cc$y[cc$trt==1]) - mean(cc$y[cc$trt==0])
    est_cc
    ```

True ITT (latent):

!!! interactive "R"
    ```r
    true_itt <- mean(df$y_true[df$trt==1]) - mean(df$y_true[df$trt==0])
    true_itt
    ```

---

## 15B. R: Multiple imputation using mice

If you don't have it installed:
`install.packages("mice")`

!!! interactive "R"
    ```r
    library(mice)

    # Use predictors: trt and severity0
    imp <- mice(df[, c("trt","severity0","y")], m=20, seed=202, printFlag=FALSE)

    # Fit the analysis model within each imputed dataset
    fit <- with(imp, lm(y ~ trt + severity0))

    # Pool results (Rubin's rules)
    pooled <- pool(fit)
    summary(pooled)
    ```

Interpretation:
- coefficient of `trt` is the adjusted treatment effect estimate
- MI accounts for uncertainty due to missingness under MAR assumptions (given included predictors)

---

## 16. Sensitivity analyses 

Because missingness may be MNAR, trials often include sensitivity analyses such as:
- different imputation models
- pattern-mixture models (e.g., shift imputed values down/up)
- tipping-point analysis

A simple “delta adjustment” idea:
- after imputation, subtract a fixed amount from imputed outcomes in one group to represent worse-than-assumed outcomes

This is not a universal method, but it illustrates sensitivity thinking.

---

## 17. Reporting: what should appear in a trial report/SAP

For analysis populations:
- define ITT (full analysis set)
- define PP criteria explicitly
- specify handling of protocol deviations

For missing data:
- describe extent and reasons for missingness
- justify assumptions (MAR vs plausible MNAR)
- specify primary method (MI or likelihood)
- include sensitivity analyses

---

## 18. Exercises

<details>
<summary>Click to try</summary>

1. Modify dropout so it depends on baseline severity only (more MAR-like). Compare CC vs MI estimates.  
2. In the R MI example, increase number of imputations from 20 to 50 and see if results stabilize.  
3. Create a PP dataset by excluding participants with severe baseline severity and compare ITT vs PP (note the bias risk).  
4. Simulate treatment switching (e.g., some control patients switch if outcome is high) and discuss whether ITT and AT differ.  
5. Write a short paragraph: which analysis population is primary and why?

</details>

---

## 19. Summary

- ITT preserves randomization and is typically primary.
- PP and as-treated analyses can be biased because they break randomization.
- Missing data is not solved by simply dropping incomplete cases.
- Multiple imputation and likelihood-based methods provide principled handling under MAR assumptions.
- Sensitivity analyses are essential when MNAR is plausible.

---
