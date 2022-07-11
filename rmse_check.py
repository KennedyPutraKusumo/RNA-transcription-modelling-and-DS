import numpy as np

def RMSE(data, predictions, std_error=1):
    rss = residual_sum_of_squares(data, predictions, std_error)
    rmse = rss / len(data)
    rmse = np.sqrt(rmse)
    return rmse


def residual_sum_of_squares(data, predictions, std_error=1):
    rsq = ((predictions - data) / std_error) ** 2
    rss = np.sum(rsq)
    return rss


if __name__ == '__main__':
    n = 10
    data = np.random.normal(0, 1, size=n)
    predictions = np.zeros_like(data)
    rmse = RMSE(data, predictions, std_error=1)
    print(f"RMSE with n={n} is {rmse}")

    n = 100
    data = np.random.normal(0, 1, size=n)
    predictions = np.zeros_like(data)
    rmse = RMSE(data, predictions, std_error=1)
    print(f"RMSE with n={n} is {rmse}")

    n = 100000
    data = np.random.normal(0, 1, size=n)
    predictions = np.zeros_like(data)
    rmse = RMSE(data, predictions, std_error=1)
    print(f"RMSE with n={n} is {rmse}")

