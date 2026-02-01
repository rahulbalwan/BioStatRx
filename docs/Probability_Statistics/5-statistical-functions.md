# Common Statistical Functions and Summaries

Descriptive statistics summarize and describe the **main features of a dataset**. They are essential before inferential or predictive analysis.

---

## 1. Measures of Central Tendency

Central tendency tells us **where the “center” of the data lies**.

### 1.1 Mean (Average)

$$
\bar{x} = \frac{1}{n}\sum_{i=1}^{n} x_i
$$

- Sensitive to outliers.

---

### 1.2 Median

The middle value when the data is sorted:

- Resistant to outliers
- Splits data into two equal halves

---

### 1.3 Mode

The most frequent value in the dataset:

- Can be multi-modal
- Useful for categorical data

---

## 2. Measures of Spread

### 2.1 Variance and Standard Deviation

$$
\mathrm{Var}(X) = \frac{1}{n}\sum_{i=1}^{n}(x_i-\bar{x})^2
$$

$$
\sigma = \sqrt{\mathrm{Var}(X)}
$$

---

### 2.2 Range

$$
\mathrm{Range} = \max(x) - \min(x)
$$

---

### 2.3 Interquartile Range (IQR)

$$
\mathrm{IQR} = Q_3 - Q_1
$$

- Resistant to extreme values  
- Useful for detecting outliers

---

## 3. Covariance and Correlation

### 3.1 Covariance

$$
\mathrm{Cov}(X,Y) = \frac{1}{n}\sum_{i=1}^{n}(X_i-\bar{X})(Y_i-\bar{Y})
$$

- Positive → variables increase together  
- Negative → one increases as the other decreases  

---

### 3.2 Correlation

$$
\rho_{X,Y} = \frac{\mathrm{Cov}(X,Y)}{\sigma_X\sigma_Y}
$$

- Ranges from -1 to 1  
- Dimensionless

---

## 4. Interactive Simulation: Python

!!! interactive "Python"
    ```python
    import numpy as np
    import pandas as pd

    np.random.seed(42)
    data = pd.DataFrame({
        "height": np.random.normal(170, 10, 100),
        "weight": np.random.normal(70, 12, 100),
        "age": np.random.randint(18, 60, 100)
    })

    # Central Tendency
    mean_height = data["height"].mean()
    median_weight = data["weight"].median()
    mode_age = data["age"].mode()[0]

    # Spread
    var_height = data["height"].var()
    sd_weight = data["weight"].std()
    range_age = data["age"].max() - data["age"].min()
    iqr_height = data["height"].quantile(0.75) - data["height"].quantile(0.25)

    # Correlation
    corr_height_weight = data["height"].corr(data["weight"])

    mean_height, median_weight, mode_age, var_height, sd_weight, range_age, iqr_height, corr_height_weight
    ```

**Run online:** [Try Python Simulation](https://replit.com/~)

---

## 5. Interactive Simulation: R

!!! interactive "R"
    ```r
    set.seed(42)
    height <- rnorm(100, mean=170, sd=10)
    weight <- rnorm(100, mean=70, sd=12)
    age <- sample(18:59, 100, replace=TRUE)
    data <- data.frame(height, weight, age)

    # Central Tendency
    mean_height <- mean(data$height)
    median_weight <- median(data$weight)
    mode_age <- as.numeric(names(sort(table(data$age), decreasing=TRUE)[1]))

    # Spread
    var_height <- var(data$height)
    sd_weight <- sd(data$weight)
    range_age <- max(data$age) - min(data$age)
    iqr_height <- IQR(data$height)

    # Correlation
    corr_height_weight <- cor(data$height, data$weight)

    list(mean_height=mean_height, median_weight=median_weight, mode_age=mode_age,
         var_height=var_height, sd_weight=sd_weight, range_age=range_age,
         iqr_height=iqr_height, corr_height_weight=corr_height_weight)
    ```

---

## 6. Visualization (Python)

!!! interactive "Python"
    ```python
    import matplotlib.pyplot as plt

    # Histogram
    plt.hist(data["height"], bins=15, edgecolor="black")
    plt.title("Histogram of Heights")
    plt.xlabel("Height")
    plt.ylabel("Frequency")
    plt.show()

    # Boxplot
    plt.boxplot(data["weight"], vert=False)
    plt.title("Boxplot of Weights")
    plt.xlabel("Weight")
    plt.show()

    # Scatter plot
    plt.scatter(data["height"], data["weight"])
    plt.title("Height vs Weight")
    plt.xlabel("Height")
    plt.ylabel("Weight")
    plt.show()
    ```

---

## 7. Exercises

<details>
<summary>Click to try!</summary>

1. Generate 100 random ages between 20 and 60. Compute mean, median, mode, SD, range, and IQR.  
2. Simulate 100 weights (Normal distribution). Plot histogram and boxplot.  
3. Compute correlation between height and weight for 50 simulated individuals.  
4. Compare standard deviation vs IQR for a dataset with outliers.

</details>

---

## 8. Summary

- Descriptive statistics provide quick insights into data.  
- Central tendency summarizes location, while spread measures variability.  
- Correlation and covariance describe relationships between variables.  
- Visualization complements numeric summaries.
