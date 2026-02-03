# Cox Proportional Hazards Regression (Time-to-Event)

Cox regression is used for **survival / time-to-event outcomes**, where:
 - not everyone experiences the event (censoring)
 - time matters (not just event yes/no)

Examples:
 - time to death
 - time to relapse
 - time to hospitalization
 - time to disease onset

---

## 1. Key ideas

### Time, event, censoring
For each subject $i$:
 - $T_i$ = observed time
 - $\delta_i$ = event indicator (1=event occurred, 0=censored)

Censoring (typical):
 - lost to follow-up
 - study ends before event

---

## 2. Hazard and Hazard Ratio (HR)

**Hazard** is the instantaneous event rate at time $t$ given survival up to $t$.

Cox model:
$$
h(t \mid X) = h_0(t)\exp(\beta_1 X_1 + \cdots + \beta_p X_p)
$$

- $h_0(t)$: baseline hazard (unspecified)
- $\exp(\beta_j)$: **hazard ratio (HR)** for a 1-unit increase in $X_j$

Interpretation:
- HR = 1.30 → 30% higher instantaneous risk
- HR = 0.70 → 30% lower instantaneous risk (protective)

---

## 3. Proportional Hazards (PH) assumption

Cox assumes hazard ratios do **not change over time**:
$$
\frac{h(t \mid X=1)}{h(t \mid X=0)} = \text{constant in } t
$$

If PH is violated:
- consider time-varying effects
- stratified Cox
- alternative models (AFT, flexible parametric, etc.)

---

## 4. Interactive simulation: survival data with known HRs

We’ll simulate with an exponential baseline hazard (for simplicity), then fit a Cox model.

!!! interactive "Python"
    ```python
    import numpy as np
    import pandas as pd

    np.random.seed(999)
    n = 2000

    age = np.random.uniform(30, 85, n)
    trt = np.random.binomial(1, 0.5, n)  # 1=drug, 0=control

    # ---- TRUE COX-LIKE LOG HAZARD MODEL (edit these) ----
    beta_age = 0.035         # HR per 1 year = exp(0.035) ~ 1.036
    beta_trt = -0.45         # HR for drug vs control = exp(-0.45) ~ 0.64
    base_hazard = 0.03       # baseline hazard rate (per unit time)
    # ----------------------------------------------------

    linpred = beta_age*age + beta_trt*trt
    hazard = base_hazard * np.exp(linpred)

    # Exponential survival time: T ~ Exp(rate=hazard)
    T = np.random.exponential(scale=1/hazard)

    # Random administrative censoring
    C = np.random.uniform(0.0, 40.0, n)
    time = np.minimum(T, C)
    event = (T <= C).astype(int)

    df = pd.DataFrame({"time": time, "event": event, "age": age, "trt": trt})

    print(df.head())
    print("\nEvent fraction:", round(df["event"].mean(), 3))
    ```

---

## 5. Fit Cox model (Python)

We’ll use `lifelines` (common in survival analysis).  
If it’s not installed in your environment, install it in your venv with `pip install lifelines`.

!!! interactive "Python"
    ```python
    from lifelines import CoxPHFitter
    import numpy as np

    cph = CoxPHFitter()
    cph.fit(df, duration_col="time", event_col="event", formula="age + trt")
    cph.print_summary()

    # HR table
    hr = np.exp(cph.params_)
    print("\nEstimated HRs:\n", hr)
    ```

Interpretation:
- `trt` HR should be around ~0.64 (protective)
- `age` HR > 1 (higher risk with older age)

---

## 6. Plot survival curves (example profiles)

Cox can give survival curves for specified covariates.

!!! interactive "Python"
    ```python
    import matplotlib.pyplot as plt
    import pandas as pd

    # Two profiles: same age, different treatment
    prof_control = pd.DataFrame({"age": [60], "trt": [0]})
    prof_drug = pd.DataFrame({"age": [60], "trt": [1]})

    sf_control = cph.predict_survival_function(prof_control)
    sf_drug = cph.predict_survival_function(prof_drug)

    plt.plot(sf_control.index, sf_control.values, label="Control (age=60)")
    plt.plot(sf_drug.index, sf_drug.values, label="Drug (age=60)")
    plt.title("Predicted Survival Curves")
    plt.xlabel("Time")
    plt.ylabel("Survival probability")
    plt.legend()
    plt.show()
    ```

---

## 7. Checking PH assumption (important!)

lifelines provides checks using Schoenfeld residuals.

!!! interactive "Python"
    ```python
    cph.check_assumptions(df, p_value_threshold=0.05, show_plots=False)
    ```

If PH is violated:
- include time interaction terms
- stratify by the violating covariate
- use flexible survival models

---

## 8. Cox vs Logistic (biostat intuition)

- Logistic: models event occurrence by a fixed time (ignores exact timing)
- Cox: uses timing + handles censoring → more efficient and appropriate for follow-up studies

---

## Exercises

<details>
<summary>Click to try</summary>

1. Increase censoring (make C smaller) and see precision drop.  
2. Change `beta_trt` to 0 and confirm HR ~ 1.  
3. Add an interaction `age*trt` and interpret effect modification.  
4. Create PH violation by making treatment effect fade over time (advanced: simulate time-varying hazards).

</details>

---

## Summary

- Cox regression models time-to-event with censoring.
- Exponentiated coefficients are **hazard ratios (HR)**.
- Proportional hazards is a key assumption — always check it.
- Cox is standard in clinical trials and cohort survival studies.
