# Survival Function, Hazard Function, and Their Relationship

This chapter introduces the two most important functions in survival analysis:

1. **Survival function** \(S(t)\): “How many survive beyond time \(t\)?”
2. **Hazard function** \(h(t)\): “How risky is it to fail right now, given survival so far?”

Everything you will do later—Kaplan–Meier, log-rank, Cox regression, parametric survival—can be understood through these two concepts.

---

## 1. The random variable \(T\)

Let \(T\) be the **true event time** for a subject.

- \(T\) is random because different subjects experience the event at different times.
- In data, we often observe \(Y = \min(T,C)\) plus an event indicator \(\delta\).

In this chapter, focus on the theoretical object \(T\).

---

## 2. Survival function \(S(t)\)

### 2.1 Definition

\[
S(t) = P(T > t)
\]

Interpretation:

> Probability that the subject survives longer than time \(t\).

### 2.2 Examples

- If \(S(5)=0.80\), then ~80% survive beyond time 5.
- If time unit is **years**, then 80% survive beyond 5 years.

### 2.3 Properties of \(S(t)\)

1. **Starts at 1**:
\[
S(0)=1
\]
(Everyone is alive at time zero)

2. **Non-increasing**:
\[
S(t_2)\le S(t_1)\quad \text{if } t_2>t_1
\]

3. **Bounded**:
\[
0 \le S(t)\le 1
\]

4. As \(t\to\infty\), \(S(t)\to 0\) (often, but not always).

---

## 3. Cumulative distribution function \(F(t)\)

Another way to describe survival time is:

\[
F(t)=P(T\le t)
\]

This is the probability the event happens by time \(t\).

Relationship:

\[
F(t)=1-S(t)
\]

So survival and cumulative incidence are two sides of the same coin (when there is a single event type).

---

## 4. Probability density and hazard: intuition

In continuous time, the event does not occur at a single exact time with positive probability, so we often use:

- density \(f(t)\)
- hazard \(h(t)\)

---

## 5. Hazard function \(h(t)\)

### 5.1 Definition (formal)

\[
h(t)=\lim_{\Delta t \to 0}\frac{P(t\le T<t+\Delta t\mid T\ge t)}{\Delta t}
\]

Interpretation:

> The instantaneous event rate at time \(t\), given survival up to time \(t\).

### 5.2 Plain-language interpretation

If hazard is high at time \(t\):

- people who are alive at \(t\) are at high risk of failing immediately.

If hazard is low:

- risk at that moment is small.

### 5.3 Hazard is NOT a probability
Hazard is a **rate**, and it can exceed 1 (depending on units).

This is a common student confusion.

---

## 6. Relationship between \(S(t)\), \(f(t)\), and \(h(t)\)

### 6.1 Density

\[
f(t) = \frac{d}{dt}F(t) = -\frac{d}{dt}S(t)
\]

### 6.2 Hazard definition in terms of density and survival

\[
h(t)=\frac{f(t)}{S(t)}
\]

This formula is extremely useful.

### 6.3 The cumulative hazard \(H(t)\)

\[
H(t)=\int_0^t h(u)\,du
\]

### 6.4 The key survival-hazard relationship

\[
S(t)=\exp(-H(t))=\exp\left(-\int_0^t h(u)\,du\right)
\]

Interpretation:

- hazard accumulates over time
- higher cumulative hazard means lower survival

This is foundational for Cox regression and parametric survival models.

---

## 7. Common hazard shapes

### 7.1 Constant hazard (Exponential model)
Examples: random failures (often unrealistic for humans)

```
hazard
  |
  |---------
  |
  +-------------- time
```

### 7.2 Increasing hazard
Examples:
- aging mortality
- chronic disease progression

```
hazard
  |
  |    /
  |   /
  |  /
  +-------------- time
```

### 7.3 Decreasing hazard
Examples:
- post-surgery early risk then recovery

```
hazard
  |
  |\ 
  | \
  |  \
  +-------------- time
```

### 7.4 Non-monotonic (rises then falls)
Examples:
- complications peak then decrease
- infection risk after surgery

```
hazard
  |
  |   /\
  |  /  \
  | /    \
  +-------------- time
```

---

## 8. Connecting hazard ratios to hazards

Cox regression models hazard as:

\[
h(t|X)=h_0(t)\exp(\beta X)
\]

If you compare two subjects with covariates \(X_1\) and \(X_2\):

\[
\frac{h(t|X_1)}{h(t|X_2)}=\exp(\beta (X_1-X_2))
\]

This ratio is constant in \(t\) → the **proportional hazards** assumption.

---

## 9. Interactive exploration: survival from hazard (Python + R)

We’ll explore how different hazard shapes generate different survival curves.

### 9.1 Python: constant hazard → exponential survival

!!! interactive "Python"
    ```python
    import numpy as np
    import matplotlib.pyplot as plt

    t = np.linspace(0, 20, 400)

    lam = 0.15  # constant hazard
    H = lam * t
    S = np.exp(-H)

    plt.plot(t, S)
    plt.title("Survival curve from constant hazard (Exponential)")
    plt.xlabel("time")
    plt.ylabel("S(t)")
    plt.show()
    ```

Try changing:
- `lam = 0.05` (slower drop)
- `lam = 0.30` (faster drop)

---

### 9.2 R: constant hazard → exponential survival

!!! interactive "R"
    ```r
    t <- seq(0, 20, length.out = 400)

    lam <- 0.15
    H <- lam * t
    S <- exp(-H)

    plot(t, S, type="l", main="Survival from constant hazard (Exponential)",
         xlab="time", ylab="S(t)")
    ```

Try:
- `lam <- 0.05`
- `lam <- 0.30`

---

## 10. Interactive exploration: Weibull hazard shapes (Python + R)

Weibull model is extremely useful because it can represent increasing or decreasing hazards.

### Weibull survival:

\[
S(t)=\exp\left(-\left(\frac{t}{\lambda}\right)^k\right)
\]

- \(k=1\) → exponential (constant hazard)
- \(k>1\) → increasing hazard
- \(k<1\) → decreasing hazard

### 10.1 Python: Weibull survival curves for different \(k\)

!!! interactive "Python"
    ```python
    import numpy as np
    import matplotlib.pyplot as plt

    t = np.linspace(0.01, 20, 500)
    lam = 8  # scale

    for k in [0.6, 1.0, 2.0, 3.0]:
        S = np.exp(- (t/lam)**k )
        plt.plot(t, S, label=f"k={k}")

    plt.title("Weibull Survival Curves (different hazard shapes)")
    plt.xlabel("time")
    plt.ylabel("S(t)")
    plt.legend()
    plt.show()
    ```

Observe:
- when \(k<1\), survival drops quickly early then slows
- when \(k>1\), survival drops slowly early then faster later

---

### 10.2 R: Weibull survival curves for different \(k\)

!!! interactive "R"
    ```r
    t <- seq(0.01, 20, length.out = 500)
    lam <- 8

    ks <- c(0.6, 1.0, 2.0, 3.0)

    plot(t, exp(-(t/lam)^ks[1]), type="l",
         main="Weibull Survival Curves (different k)",
         xlab="time", ylab="S(t)", ylim=c(0,1))

    for (k in ks[-1]) {
      lines(t, exp(-(t/lam)^k))
    }

    legend("topright", legend=paste("k=", ks), lty=1)
    ```

---

## 11. Visualizing hazard directly (Python + R)

For Weibull:

\[
h(t)=\frac{k}{\lambda}\left(\frac{t}{\lambda}\right)^{k-1}
\]

### 11.1 Python: Weibull hazard curves

!!! interactive "Python"
    ```python
    import numpy as np
    import matplotlib.pyplot as plt

    t = np.linspace(0.01, 20, 500)
    lam = 8

    for k in [0.6, 1.0, 2.0, 3.0]:
        h = (k/lam) * (t/lam)**(k-1)
        plt.plot(t, h, label=f"k={k}")

    plt.title("Weibull Hazard Curves")
    plt.xlabel("time")
    plt.ylabel("h(t)")
    plt.legend()
    plt.show()
    ```

Interpret:
- \(k<1\): high early hazard, decreasing
- \(k=1\): constant
- \(k>1\): increasing

---

### 11.2 R: Weibull hazard curves

!!! interactive "R"
    ```r
    t <- seq(0.01, 20, length.out = 500)
    lam <- 8
    ks <- c(0.6, 1.0, 2.0, 3.0)

    h_weibull <- function(t, k, lam) (k/lam) * (t/lam)^(k-1)

    plot(t, h_weibull(t, ks[1], lam), type="l",
         main="Weibull Hazard Curves",
         xlab="time", ylab="h(t)")

    for (k in ks[-1]) {
      lines(t, h_weibull(t, k, lam))
    }

    legend("topleft", legend=paste("k=", ks), lty=1)
    ```

---

## 12. Practical interpretation for biostatistics

When you see a survival curve:
- it tells you **probability of surviving beyond time t**
- it is easy to interpret clinically

When you see a hazard ratio:
- it tells you **relative instantaneous risk**
- it is the main output of Cox regression

Remember:
- Survival is about probability
- Hazard is about rate

---

## 13. Key takeaways

- \(S(t)=P(T>t)\) is the survival probability beyond time \(t\).
- \(h(t)\) is instantaneous risk rate given survival up to time \(t\).
- They are linked by:
  \[
  S(t)=\exp\left(-\int_0^t h(u)\,du\right)
  \]
- Different hazard shapes lead to different survival curve shapes.
- Cox regression is built on modeling hazard and hazard ratios.

---

## 14. Exercises

<details>
<summary>Click to try</summary>

1. If \(S(3)=0.9\), interpret it in plain language.  
2. Explain why hazard is not a probability.  
3. Which Weibull shape parameter \(k\) produces increasing hazard?  
4. For Weibull with \(k<1\), describe the clinical scenario this could represent.  
5. Simulate survival curves for different hazards and describe their differences.

</details>
