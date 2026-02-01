# Expectation & Variance

Expectation and variance are **fundamental properties of random variables** that summarize their behavior mathematically.

- **Expectation (Mean):** the “center” or long-term average of a random variable  
- **Variance:** how much values spread around the mean  
- **Standard deviation:** square root of variance (same units as the variable)

---

## 1. Expectation (Mean)

For a **discrete random variable** $X$ with PMF $P(X=x)$:

$$
\mathbb{E}[X] = \sum_x x \, P(X=x)
$$

For a **continuous random variable** with PDF $f(x)$:

$$
\mathbb{E}[X] = \int_{-\infty}^{\infty} x\, f(x)\, dx
$$

**Interpretation:** expectation is the average outcome if the experiment is repeated many times.

---

## 2. Variance

Variance measures the **spread** of a random variable around its mean:

$$
\mathrm{Var}(X) = \mathbb{E}\big[(X - \mathbb{E}[X])^2\big]
$$

- Standard deviation: $\sigma_X = \sqrt{\mathrm{Var}(X)}$  
- Small variance → values close to the mean  
- Large variance → values more spread out  

---

## Example 1: Discrete Random Variable

Number of heads in **3 coin tosses**:

| $X$ (Heads) | 0   | 1   | 2   | 3   |
|------------:|:---:|:---:|:---:|:---:|
| $P(X=x)$    | 1/8 | 3/8 | 3/8 | 1/8 |

$$
\mathbb{E}[X] =
0\cdot\frac{1}{8} + 1\cdot\frac{3}{8} + 2\cdot\frac{3}{8} + 3\cdot\frac{1}{8} = 1.5
$$

$$
\mathrm{Var}(X) = \sum_x (x - 1.5)^2 P(X=x) = 0.75
$$

---

## Example 2: Continuous Random Variable

Simulate heights:

- Mean: $\mu = 170$ cm  
- Standard deviation: $\sigma = 10$ cm  

We can estimate mean and variance from simulated data.

---

## Interactive Simulation: Discrete RV (Python)

!!! interactive "Python"
    ```python
    import numpy as np

    np.random.seed(42)

    # Simulate 3 coin tosses, 1000 trials
    tosses = np.random.choice([0, 1], size=(1000, 3))
    heads_count = np.sum(tosses, axis=1)

    # Mean and variance
    mean_heads = np.mean(heads_count)
    var_heads = np.var(heads_count)

    mean_heads, var_heads
    ```

**Run online:** [Try Python Simulation](https://replit.com/~)

---

## Interactive Simulation: Continuous RV (Python)

!!! interactive "Python"
    ```python
    import numpy as np

    np.random.seed(42)

    # Simulate 1000 heights
    heights = np.random.normal(loc=170, scale=10, size=1000)

    mean_height = np.mean(heights)
    var_height = np.var(heights)

    mean_height, var_height
    ```

**Run online:** [Try Python Simulation](https://replit.com/~)

---

## Interactive Simulation: Discrete RV (R)

!!! interactive "R"
    ```r
    set.seed(42)

    tosses <- matrix(sample(0:1, 3000, replace=TRUE), ncol=3)
    heads_count <- rowSums(tosses)

    mean_heads <- mean(heads_count)
    var_heads <- var(heads_count)

    mean_heads
    var_heads
    ```

**Run online:** [Try R Simulation](https://rdrr.io/snippets/)

---

## Interactive Simulation: Continuous RV (R)

!!! interactive "R"
    ```r
    set.seed(42)

    heights <- rnorm(1000, mean=170, sd=10)

    mean_height <- mean(heights)
    var_height <- var(heights)

    mean_height
    var_height
    ```

**Run online:** [Try R Simulation](https://rdrr.io/snippets/)

---

## Exercises

<details>
<summary>Click to try!</summary>

1. Toss 5 coins, repeat 50 trials. Compute mean and variance of heads.  
2. Roll a die 20 times, repeat 50 trials. Compute mean and variance.  
3. Simulate 500 weights (mean = 70, SD = 12). Plot histogram; compute mean & variance.  
4. Compare sample variance with theoretical variance for a discrete RV.  
5. Use Python or R to verify expectation and variance from simulated data.

</details>

---

## Summary

- Expectation gives the average value of a random variable.
- Variance quantifies how spread out values are around the mean.
- Simulations make these ideas concrete and visual.
