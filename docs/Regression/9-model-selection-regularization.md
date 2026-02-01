# Model Selection & Regularization (Ridge / LASSO) — Biostat perspective

When you have many predictors:
- clinical + lab + genetic variables
- high-dimensional EHR features

You risk **overfitting** and unstable estimates.

Regularization helps:
- better prediction performance
- stability under collinearity

---

## 1. The idea

For linear regression:
- Ridge shrinks coefficients toward 0 (keeps all predictors)
- LASSO can set some coefficients exactly to 0 (feature selection)

For logistic regression:
- same ideas using penalized likelihood

---

## 2. Interactive: compare OLS vs Ridge vs LASSO (prediction)

!!! interactive "Python"
    ```python
    import numpy as np
    from sklearn.model_selection import train_test_split
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import Pipeline
    from sklearn.linear_model import LinearRegression, Ridge, Lasso
    from sklearn.metrics import mean_squared_error

    np.random.seed(101)

    n = 600
    p = 25

    X = np.random.normal(0, 1, size=(n, p))

    # True model uses only a few predictors
    true_beta = np.zeros(p)
    true_beta[[0, 3, 7, 10]] = [3.0, -2.0, 1.5, 2.2]
    y = X @ true_beta + np.random.normal(0, 3, size=n)

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

    models = {
        "OLS": Pipeline([("sc", StandardScaler()), ("m", LinearRegression())]),
        "Ridge(alpha=5)": Pipeline([("sc", StandardScaler()), ("m", Ridge(alpha=5))]),
        "Lasso(alpha=0.1)": Pipeline([("sc", StandardScaler()), ("m", Lasso(alpha=0.1))]),
    }

    for name, model in models.items():
        model.fit(X_train, y_train)
        pred = model.predict(X_test)
        rmse = mean_squared_error(y_test, pred, squared=False)
        print(name, "RMSE:", round(rmse, 3))

    # Show which coefficients LASSO kept
    lasso = models["Lasso(alpha=0.1)"].named_steps["m"]
    coef = lasso.coef_
    kept = np.where(np.abs(coef) > 1e-6)[0]
    print("\nLASSO kept predictors:", kept)
    ```

Try:
- Increase noise SD (currently `3`) → selection gets harder
- Increase p to 100 and see regularization become more helpful
- Change LASSO alpha to 0.05 or 0.2

---

## 3. Regularized Logistic Regression (binary outcomes)

!!! interactive "Python"
    ```python
    import numpy as np
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import train_test_split
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import Pipeline

    np.random.seed(202)

    n = 1200
    p = 30
    X = np.random.normal(0, 1, size=(n, p))

    true_beta = np.zeros(p)
    true_beta[[1, 4, 9]] = [1.0, -1.2, 0.8]

    lp = X @ true_beta - 0.3
    prob = 1/(1 + np.exp(-lp))
    y = np.random.binomial(1, prob, size=n)

    Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=0.3, random_state=1)

    # L2 penalty (ridge-like)
    clf = Pipeline([
        ("sc", StandardScaler()),
        ("m", LogisticRegression(penalty="l2", C=1.0, max_iter=2000))
    ])

    clf.fit(Xtr, ytr)
    pred = clf.predict_proba(Xte)[:, 1]
    print("AUC:", round(roc_auc_score(yte, pred), 3))
    ```

Biostat note:
- Penalization is common in **risk prediction** models
- For causal interpretation, regularization can complicate inference (needs careful reporting)

---

## Exercises

<details>
<summary>Click to try</summary>

1. Increase p to 200 and compare OLS vs ridge/lasso RMSE.  
2. Make predictors correlated and see ridge stabilize performance.  
3. In logistic simulation, change baseline event rate by shifting intercept (-0.3).  

</details>
