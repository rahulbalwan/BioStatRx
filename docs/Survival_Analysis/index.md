# Survival Analysis

Survival analysis is the set of statistical tools used to study **time-to-event outcomes**—for example, time to death, relapse, hospital readmission, device failure, or recovery. What makes survival data special is that:

- **Not everyone experiences the event** during follow-up (they are *censored*).
- Risk changes over time and we care about the entire trajectory, not just a final outcome.
- We often want both **descriptive curves** (Kaplan–Meier) and **adjusted effect estimates** (Cox regression).
- Real clinical data frequently includes complexities like **time-dependent treatments**, **competing risks**, and **recurrent events**.

This section is designed for **biostatistics and biomedical research**, with hands-on examples in **R and Python** and a clear focus on practical interpretation.

---

## What you will learn

By the end of this module, you will be able to:

- Explain **censoring**, **risk sets**, and why standard regression fails for time-to-event outcomes  
- Construct and interpret **Kaplan–Meier curves** with confidence intervals and median survival  
- Compare groups using the **log-rank test**  
- Fit and interpret **Cox proportional hazards models** (simple + multivariable)  
- Check assumptions using **Schoenfeld residuals** and other diagnostics  
- Handle **time-dependent covariates** correctly (avoid immortal time bias)  
- Fit and interpret **parametric survival models** (Weibull, log-normal, AFT)  
- Analyze **competing risks** using CIF and Fine–Gray models  
- Account for clustering using **frailty models** and robust standard errors  
- Model **recurrent events** (Andersen–Gill, PWP)  
- Produce **publication-quality plots and tables**

---

## Who this module is for

This module is ideal for:

- Students learning biostatistics / epidemiology
- Researchers analyzing clinical or public health time-to-event data
- Anyone reading medical papers and wanting to understand survival results properly

You do **not** need to be an expert in probability, but you should be comfortable with:
- basic regression ideas
- interpreting coefficients and confidence intervals
- working with simple datasets in Python or R

---

## Recommended tools (R + Python)

You can follow the entire module using either language.

### R packages
- `survival` (core survival modeling)
- `survminer` (best KM plots + risk tables)
- `cmprsk` (competing risks + Fine–Gray)
- `flexsurv` (parametric models)
- `riskRegression` (advanced competing risk prediction)

### Python packages
- `lifelines` (KM, Cox, time-varying Cox, parametric models)
- `pandas`, `numpy`, `matplotlib` (data + visualization)

---

## Module roadmap

Below is the recommended learning path. Each page includes:

- strong conceptual explanation
- worked examples
- interactive code in both R and Python
- interpretation in biostatistical language

### Foundations
1. **01 — Why Survival Analysis?** (what makes time-to-event special)  
2. **02 — Censoring and Follow-up** (right censoring, informative censoring)  
3. **03 — Survival and Hazard Functions** (S(t), h(t), H(t), intuition)  
4. **04 — Risk Sets and Event Times** (the engine behind KM + Cox)

### Kaplan–Meier + Hypothesis Testing
5. **05 — Kaplan–Meier Estimator** (step-by-step, manual + software)  
6. **06 — KM Confidence Intervals + Median Survival**  
7. **07 — Log-rank Test** (observed vs expected events)

### Regression modeling (core)
8. **08 — Cox Proportional Hazards Model** (HR interpretation, partial likelihood intuition)  
9. **09 — Cox Diagnostics** (PH checks, Schoenfeld residuals, functional form, influence)  
10. **10 — Time-dependent Covariates** (start–stop format, avoid immortal time bias)

### Beyond Cox
11. **11 — Parametric Survival Models** (Weibull, log-normal, AFT; AIC + overlays)  
12. **12 — Competing Risks** (CIF, cause-specific Cox, Fine–Gray)  
13. **13 — Frailty Models** (clustering + unobserved heterogeneity)  
14. **14 — Recurrent Events** (Andersen–Gill, PWP total-time/gap-time)

### Putting it all together
15. **15 — Reporting + Publication-Quality Plots/Tables** (KM with risk table, Cox tables, Methods writing)

---

## A note on interpretation

In biomedical research, survival analysis is not just about “statistics”—it is about answering meaningful questions:

- How quickly do events occur?
- Does a treatment reduce risk?
- How do patient factors change prognosis?
- What is the absolute probability of an outcome by a clinically meaningful time?
- What biases (immortal time, competing risks) could distort conclusions?

This module emphasizes **interpretability**, **valid assumptions**, and **realistic workflows**.

---

