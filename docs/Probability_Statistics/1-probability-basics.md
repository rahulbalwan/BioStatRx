# Probability Basics

Probability measures the **likelihood of an event** occurring. It is fundamental for statistical analysis and forms the foundation for modeling uncertainty in real-world phenomena.

---

## Basic Definitions

- **Experiment:** A process that generates outcomes (e.g., toss a coin, roll a die).
- **Sample Space $S$:** The set of all possible outcomes.
- **Event $E$:** A subset of the sample space representing outcomes of interest.
- **Probability:** A measure that satisfies $0 \le P(E) \le 1$, often computed as:

$$
P(E) = \frac{\text{Number of favorable outcomes}}{\text{Total outcomes}}
$$

---

## Key Properties

1. **Non-negativity:** $P(E) \ge 0$ for any event $E$
2. **Normalization:** $P(S) = 1$
3. **Additivity:** If $A$ and $B$ are mutually exclusive, then
   $$
   P(A \cup B) = P(A) + P(B)
   $$

---

## Example 1: Coin Toss

Consider tossing a fair coin once:

- Sample Space: $S = \{H, T\}$
- Event: Getting a Head $E = \{H\}$
- Probability: $P(E) = \frac{1}{2} = 0.5$

---

## Example 2: Rolling a Die

- Sample Space: $S = \{1,2,3,4,5,6\}$
- Event: Rolling an even number $E = \{2,4,6\}$
- Probability: $P(E) = \frac{3}{6} = 0.5$

---

## Interactive Simulation in Python

!!! interactive "Python"
    ```python
    import random

    # Simulate 10 coin tosses
    coin = ["H", "T"]
    tosses = [random.choice(coin) for _ in range(10)]
    print("Coin tosses:", tosses)

    # Simulate rolling a die 10 times
    die_rolls = [random.randint(1, 6) for _ in range(10)]
    print("Die rolls:", die_rolls)
    ```

**Run this code online:** [Try Python Simulation](https://replit.com/~)

---

## Interactive Simulation in R

!!! interactive "R"
    ```r
    # Simulate 10 coin tosses
    coin <- c("H", "T")
    tosses <- sample(coin, 10, replace = TRUE)
    print(tosses)

    # Simulate rolling a die 10 times
    die_rolls <- sample(1:6, 10, replace = TRUE)
    print(die_rolls)
    ```

**Run this code online:** [Try R Simulation](https://rdrr.io/snippets/)

---

## Exercises

<details>
<summary>Click to try the exercises!</summary>

1. Toss a fair coin 20 times. Count the number of heads and tails.  
2. Roll a six-sided die 15 times. What is the probability of rolling a number greater than 4?  
3. Toss two coins simultaneously. List the sample space and calculate the probability of getting exactly one head.

</details>

---

## Summary

- Probability quantifies uncertainty in random experiments.
- Understanding sample spaces and events is crucial for statistical analysis.
- Python and R simulations help visualize probability concepts interactively.
