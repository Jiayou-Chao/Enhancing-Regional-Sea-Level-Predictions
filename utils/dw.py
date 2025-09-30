import numpy as np
from scipy.optimize import minimize

def logistic(t, alpha, beta, t0):
    return 1 / (1 + alpha * np.exp(-beta * (t - t0)))

def objective(params, lambda_):
    alpha, beta = params
    t0 = 2020

    w_2021 = logistic(2021, alpha, beta, t0)
    w_2100 = logistic(2100, alpha, beta, t0)

    error = (w_2021 - 0.0)**2 + (w_2100 - 1.0)**2

    regularization = lambda_ * (alpha**2 + beta**2)

    return error + regularization

lambda_ = 0.001

initial_guess = [1.0, 0.1]

result = minimize(objective, initial_guess, args=(lambda_,), method='Nelder-Mead')

alpha_opt, beta_opt = result.x

print(f"Optimized alpha: {alpha_opt}")
print(f"Optimized beta: {beta_opt}")

t0 = 2020
w_2021 = logistic(2021, alpha_opt, beta_opt, t0)
w_2100 = logistic(2100, alpha_opt, beta_opt, t0)

print(f"w(2021) = {w_2021:.6f}")
print(f"w(2100) = {w_2100:.6f}")