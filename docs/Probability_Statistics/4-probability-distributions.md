# Probability Distributions

Probability distributions describe **how probabilities are assigned** to the possible outcomes of a random variable. They are fundamental in biostatistics, modeling uncertainty, and making data-driven decisions.

---

## 1. Discrete Probability Distributions

Discrete distributions describe **random variables with countable outcomes**.

### 1.1 Bernoulli Distribution

- Two outcomes: success (1) or failure (0)
- PMF:

$$
P(X=1)=p,\quad P(X=0)=1-p
$$

- Mean: $\mathbb{E}[X]=p$
- Variance: $\mathrm{Var}(X)=p(1-p)$

**Example:** coin toss with $p=0.5$

---

### 1.2 Binomial Distribution

- Number of successes in $n$ independent Bernoulli trials
- PMF:

$$
P(X=k)=\binom{n}{k}p^k(1-p)^{n-k},\quad k=0,1,\dots,n
$$

- Mean: $\mathbb{E}[X]=np$
- Variance: $\mathrm{Var}(X)=np(1-p)$

---

### 1.3 Poisson Distribution

- Number of events in a fixed interval with rate $\lambda$
- PMF:

$$
P(X=k)=\frac{\lambda^k e^{-\lambda}}{k!},\quad k=0,1,2,\dots
$$

- Mean and variance: $\mathbb{E}[X]=\mathrm{Var}(X)=\lambda$

---

## 2. Continuous Probability Distributions

Continuous distributions describe **random variables that can take any real value** in an interval.

### 2.1 Uniform Distribution

- PDF on $[a,b]$:

$$
f(x)=\frac{1}{b-a},\quad a\le x\le b
$$

- Mean: $\mathbb{E}[X]=\frac{a+b}{2}$
- Variance: $\mathrm{Var}(X)=\frac{(b-a)^2}{12}$

---

### 2.2 Normal (Gaussian) Distribution

- PDF:

$$
f(x)=\frac{1}{\sigma\sqrt{2\pi}}\exp\!\Big(-\frac{(x-\mu)^2}{2\sigma^2}\Big)
$$

- Mean and variance: $\mathbb{E}[X]=\mu,\quad \mathrm{Var}(X)=\sigma^2$

---

### 2.3 Exponential Distribution

- PDF:

$$
f(x)=\lambda e^{-\lambda x},\quad x\ge 0
$$

- Mean: $1/\lambda$
- Variance: $1/\lambda^2$

---

## 3. Interactive Simulations in Python

!!! interactive "Python"
    ```python
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.stats import binom, norm, expon

    # ----- Binomial -----
    n, p = 10, 0.5
    x_binom = np.arange(0, n + 1)
    plt.bar(x_binom, binom.pmf(x_binom, n, p), edgecolor="black")
    plt.title("Binomial PMF (n=10, p=0.5)")
    plt.xlabel("x")
    plt.ylabel("P(X=x)")
    plt.show()

    # ----- Normal -----
    x_norm = np.linspace(140, 200, 200)
    plt.plot(x_norm, norm.pdf(x_norm, loc=170, scale=10))
    plt.title("Normal PDF (mu=170, sd=10)")
    plt.xlabel("x")
    plt.ylabel("f(x)")
    plt.show()

    # ----- Exponential -----
    x_exp = np.linspace(0, 10, 200)
    lam = 0.5
    plt.plot(x_exp, expon.pdf(x_exp, scale=1/lam))
    plt.title("Exponential PDF (lambda=0.5)")
    plt.xlabel("x")
    plt.ylabel("f(x)")
    plt.show()
    ```

**Run online:** [Try Python Simulation](https://replit.com/~)

---

## 4. Interactive Simulations in R

!!! interactive "R"
    ```r
    # ----- Binomial -----
    x_binom <- 0:10
    probs <- dbinom(x_binom, size=10, prob=0.5)
    barplot(probs, names.arg=x_binom, main="Binomial PMF (n=10, p=0.5)")

    # ----- Normal -----
    x_norm <- seq(140, 200, length=200)
    plot(x_norm, dnorm(x_norm, mean=170, sd=10), type="l",
         main="Normal PDF (mu=170, sd=10)", xlab="x", ylab="f(x)")

    # ----- Exponential -----
    x_exp <- seq(0, 10, length=200)
    plot(x_exp, dexp(x_exp, rate=0.5), type="l",
         main="Exponential PDF (lambda=0.5)", xlab="x", ylab="f(x)")
    ```

**Run online:** [Try R Simulation](https://rdrr.io/snippets/)

---

## 5. Exercises

<details>
<summary>Click to try!</summary>

1. Simulate 10 coin tosses 1000 times. Plot the histogram of number of heads (Binomial).  
2. Generate 500 random heights from Normal(170, 10). Compute mean and variance.  
3. Simulate Poisson arrivals with $\lambda=3$ for 1000 intervals. Plot histogram.  
4. Generate 1000 Uniform(0,1) values. Plot histogram and compute mean & variance.  
5. Compare Uniform vs Normal distributions visually in Python or R.

</details>

---

## 6. Summary

- Discrete distributions: Bernoulli, Binomial, Poisson  
- Continuous distributions: Uniform, Normal, Exponential  
- Distributions are core building blocks in biostatistics and modeling uncertainty.
