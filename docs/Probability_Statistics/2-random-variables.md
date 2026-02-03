# Random Variables

A **random variable** is a variable that takes numerical values based on the outcome of a **random experiment**. It allows us to quantify randomness mathematically.

---

## Types of Random Variables

1. **Discrete random variables**
   - Take **countable** values (integers or counts)
   - Examples: number of heads in coin tosses; number of patients who recover

2. **Continuous random variables**
   - Take any **real value** within an interval
   - Examples: height, weight, blood pressure, time until failure

---

## Probability Functions

- **Discrete RV:** Probability Mass Function (PMF)

$$
P(X=x) = \text{Probability that } X \text{ equals } x
$$

- **Continuous RV:** Probability Density Function (PDF)

$$
f(x) \ge 0, \quad \int_{-\infty}^{\infty} f(x)\,dx = 1
$$

- **Cumulative Distribution Function (CDF):**

$$
F(x) = P(X \le x)
$$

---

## Example 1: Discrete Random Variable

Number of heads in **3 coin tosses**:

- Random Variable: $X =$ number of heads
- Possible values: $X = 0,1,2,3$
- PMF:

$$
P(X=0)=\frac{1}{8},\quad
P(X=1)=\frac{3}{8},\quad
P(X=2)=\frac{3}{8},\quad
P(X=3)=\frac{1}{8}
$$

---

## Example 2: Continuous Random Variable

- Random Variable: $Y =$ height of a person (cm)
- Continuous values: $140 \le Y \le 200$
- A common model is the **normal distribution**:

$$
f(y) = \frac{1}{\sigma \sqrt{2\pi}}
\exp\!\Big(-\frac{(y-\mu)^2}{2\sigma^2}\Big)
$$

---

## Interactive Simulation in Python: Discrete RV

!!! interactive "Python"
    ```python
    import numpy as np

    np.random.seed(42)

    # Simulate 3 coin tosses, 10 trials (0=T, 1=H)
    tosses = np.random.choice([0, 1], size=(10, 3))
    heads_count = np.sum(tosses, axis=1)
    print("Number of heads per trial:", heads_count)
    ```

**Run online:** [Try Python Simulation](https://replit.com/~)

---

## Interactive Simulation in Python: Continuous RV

!!! interactive "Python"
    ```python
    import numpy as np
    import matplotlib.pyplot as plt

    np.random.seed(42)

    # Simulate 1000 heights (Normal distribution)
    heights = np.random.normal(loc=170, scale=10, size=1000)

    plt.hist(heights, bins=30, edgecolor="black")
    plt.title("Simulated Heights of Individuals")
    plt.xlabel("Height (cm)")
    plt.ylabel("Frequency")
    plt.show()
    ```

**Run online:** [Try Python Simulation](https://replit.com/~)

---

## Interactive Simulation in R: Discrete RV

!!! interactive "R"
    ```r
    set.seed(42)

    # Simulate 3 coin tosses, 10 trials
    tosses <- matrix(sample(0:1, 30, replace=TRUE), ncol=3)
    heads_count <- rowSums(tosses)
    print(heads_count)
    ```

**Run online:** [Try R Simulation](https://rdrr.io/snippets/)

---

## Interactive Simulation in R: Continuous RV

!!! interactive "R"
    ```r
    set.seed(42)

    # Simulate 1000 heights
    heights <- rnorm(1000, mean=170, sd=10)
    hist(heights, breaks=30, main="Simulated Heights",
         xlab="Height (cm)", ylab="Frequency")
    ```

**Run online:** [Try R Simulation](https://rdrr.io/snippets/)

---

## Properties of Random Variables

| Property | Definition |
|----------|------------|
| Mean (Expectation) | Average value of the RV |
| Variance | Measure of spread of the RV |
| Standard Deviation | Square root of variance |
| PMF / PDF | Probability distribution of the RV |
| CDF | Cumulative probability function |

---

## Exercises

<details>
<summary>Click to try!</summary>

1. Simulate 5 coin tosses 20 times. Count number of heads for each trial.  
2. Roll a die 10 times. Record the sum of outcomes.  
3. Simulate 500 random heights using Normal(mean = 165, sd = 12). Plot histogram.  
4. Compare discrete and continuous RVs using Python or R.  
5. Calculate and plot the CDF of a discrete random variable.

</details>

---

## Summary

- Random variables **quantify randomness** numerically.
- Discrete vs continuous RVs use different probability functions (PMF vs PDF).
- Simulations help visualize distributions and understand their behavior.
