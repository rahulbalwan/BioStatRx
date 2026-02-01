# Regression

Regression is the core toolkit for **quantifying relationships** between variables and making **predictions**—and it sits at the heart of biostatistics, epidemiology, and clinical research.

In biomedical settings, regression helps answer questions like:

- How does **blood pressure** change with **age** and **BMI**?
- Does a **treatment** change the probability of an outcome after adjusting for covariates?
- Which risk factors are associated with **disease incidence**, **counts**, or **time-to-event** outcomes?
- How do we interpret results correctly and avoid common traps like **confounding** and **collinearity**?

This module is designed for **biostatistical applications** and emphasizes:
- clean interpretation (effect sizes + uncertainty)
- practical modeling workflow
- diagnostics and assumptions
- interactive coding in **R and Python**

---

## What you will learn

By the end of this module, you will be able to:

- Choose the correct regression model for different outcome types
- Create and interpret
- Diagnose and fix model problems
- Incorporate
- Write results like a biostat paper (tables, plots, wording)

---

## Who this module is for

This module is ideal for:
- students learning regression for biostatistics / epidemiology
- researchers working with clinical or public health datasets
- anyone preparing analyses for papers, theses, or reports

You should be comfortable with:
- basic probability/statistics concepts (mean, variance, distributions)
- reading a dataset and basic plotting
- interpreting a confidence interval and p-value

---

## Recommended tools (R + Python)

You can follow everything in either language.

### R packages
- `stats` (lm, glm)
- `car` (diagnostics, VIF)
- `broom` (tidy summaries)
- `ggplot2` (plots)
- `splines` (natural splines)
- `MASS` (negative binomial: `glm.nb`)
- `survival` (Cox regression)

### Python packages
- `statsmodels` (OLS/GLM, inference)
- `scikit-learn` (regularization + prediction workflow)
- `pandas`, `numpy` (data)
- `matplotlib` (plots)

---

## Module roadmap

Follow the pages in order for a smooth learning path. Each section includes:
- clear concepts
- biostatistical examples
- code in **R and Python**
- interpretation tips and common pitfalls

1. **Overview**  
   What regression is, outcome types, choosing the right model, and a biostat workflow.

2. **Simple Linear Regression**  
   Modeling a continuous outcome with one predictor; slope interpretation and prediction.

3. **Multiple Linear Regression**  
   Adjusting for confounders; interpreting adjusted effects; partial relationships.

4. **Logistic Regression**  
   Binary outcomes; odds ratios; predicted probabilities; classification vs inference.

5. **Diagnostics (Linear)**  
   Residual plots, normality, heteroscedasticity, leverage/influence, and what to do when assumptions fail.

6. **Categorical Predictors**  
   Dummy coding, reference groups, interpretation of group effects, ANOVA-style comparisons.

7. **Interactions**  
   Effect modification; how to model and interpret interaction terms; plotting interactions.

8. **Nonlinearity & Splines**  
   Flexible modeling of continuous predictors; splines in R/Python; avoiding overfitting.

9. **Collinearity & Confounding**  
   VIF, identifiability, unstable coefficients, confounding logic, DAG intuition (practical perspective).

10. **Regularization**  
   Ridge/Lasso/Elastic Net; cross-validation; when and why regularization matters.

11. **Reporting Results**  
   Tables, figures, clinical interpretation, and writing regression results for papers.

12. **Poisson Regression**  
   Count outcomes and rates; offsets; interpreting incidence rate ratios (IRR).

13. **Negative Binomial Regression**  
   Overdispersion; why Poisson fails; NB model interpretation; model checks.

14. **Cox Regression**  
   Time-to-event regression; hazard ratios; connection to survival analysis.

---
