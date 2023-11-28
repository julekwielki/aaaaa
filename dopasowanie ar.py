import numpy as np

def func(x, c, a, n):
    return c * (1 - np.exp(-a * np.power(x, n)))


DR = [3, 6, 12, 24, 60]
HR = [1.02, 0.66, 1.61, 1.57, 4.86]
