
# Clinical Trials

Clinical trials are the gold standard for evaluating medical interventions because they are designed to support **causal inference**. In biostatistics, clinical trials sit at the intersection of study design, probability, modeling, ethics, and reproducible reporting.

This module is written for learners targeting **biostatistics / medical statistics** work (academia, industry, and public health). Each chapter mixes:
- clear conceptual explanations
- practical trial scenarios from biomedical settings
- hands-on **Python and R** code
- analysis and reporting workflows that mirror real trial practice

---

## What you will learn

By the end of this section, you should be able to:

- Explain why randomization supports causal inference and what can go wrong operationally
- Choose an appropriate design (parallel, cross-over, cluster, factorial, adaptive) for a biomedical question
- Define **outcomes**, **endpoints**, and **estimands** clearly (and handle intercurrent events)
- Perform sample size and power calculations for common endpoint types
- Understand ITT vs per-protocol vs as-treated analyses and handle missing data responsibly
- Summarize and model adverse events using both risk and exposure-adjusted rates
- Report trial results transparently with effect sizes, confidence intervals, and CONSORT-style structure

---

## Chapters

1. **Foundations & Phases**  
   What clinical trials are, why they matter, and how phases I–IV differ.

2. **Randomization & Allocation Concealment**  
   Simple/block/stratified randomization, concealment vs blinding, and practical list generation.

3. **Trial Designs**  
   Parallel-group, cross-over, cluster randomized trials, factorial designs, and adaptive design overview.

4. **Outcomes, Endpoints, and Estimands**  
   How to translate a clinical question into an estimand; handling rescue medication, switching, discontinuation, and other intercurrent events.

5. **Sample Size & Power**  
   Continuous, binary, and time-to-event power; dropout inflation; simulation-based power intuition.

6. **Analysis Populations & Missing Data**  
   ITT, PP, and as-treated; missing data mechanisms (MCAR/MAR/MNAR); multiple imputation in R and practical strategies in Python.

7. **Safety Analysis & Adverse Events**  
   AE/SAE/AESI concepts; participant incidence vs event rates; exposure-adjusted incidence; Poisson/negative binomial modeling for AE rates.

8. **Reporting & Interpretation**  
   CONSORT flow and tables, effect sizes and confidence intervals, forest plots, subgroup interpretation, and reproducibility best practices.

---

## How to use this module

A simple learning path:
1. Read Chapter 01 → 04 to understand trial logic and definitions.
2. Use Chapter 05 to learn how design choices connect to sample size and feasibility.
3. Use Chapter 06 and 07 to build practical analysis skills.
4. Use Chapter 08 to learn how to present results professionally.

---

## Suggested prerequisites

- Basic probability and random variables
- Linear regression fundamentals
- Comfort with reading code in Python or R

If needed, review:
- *Probability & Statistics Review* (Foundations section)
- *Regression* module (especially linear and logistic regression)

---

## Next steps

If you want to extend Clinical Trials further, strong next topics are:
- interim monitoring / group sequential designs
- multiplicity and hierarchical testing
- Bayesian trial designs
- adaptive randomization and response-adaptive designs
- pragmatic trials and real-world evidence links
- estimands in survival and competing risks settings
