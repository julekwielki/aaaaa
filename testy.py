import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad


def integrand(x):
    return x**2 * np.exp(-x)


dose = []
dose_rate = []
dose_rate2 = []
dose_rate_val = 1
dose_rate_val2 = 2
time = []
divide = 10
maxx = 300
T1 = 1
T2 = 10

T21 = 12
T22 = 16
for t in range(0, maxx):
    time.append(t/divide)
    if T1*divide <= t < T2*divide:
        dose_rate.append(dose_rate_val)
    else:
        dose_rate.append(0)

    if T21*divide <= t < T22*divide:
        dose_rate2.append(dose_rate_val2)
    else:
        dose_rate2.append(0)


pa = []
pa2 = [0]
pa3 = []
paz = []

"""for t in time:
    p = dose_rate[t]**2 * np.exp(-dose_rate[t]) * t**2 * np.exp(-t)

    pz = dose_rate[t]**2 * np.exp(-dose_rate[t]) * 2 * (1 - (t**2/2 + t + 1) * np.exp(-t))
    pa.append(p)

    paz.append(pz)"""

"""
for t in time:
    p = t ** 2 * np.exp(-t)
    pa.append(p)
    pz = 2 * (1 - (t**2/2 + t + 1) * np.exp(-t))
    # pa2.append(pa2[-1]+p/divide)
    # a = quad(integrand, 0, t)
    # pa3.append(a[0])
"""

for t in time:
    pz = 0
    if t > T1:

        if dose_rate[time.index(t)] > 0:
            #pz = dose_rate_val**2 * np.exp(-dose_rate_val) * 2 * (1 - (t ** 2 / 2 + t + 1) * np.exp(-t))
            pz = dose_rate_val**2 * np.exp(-dose_rate_val) * (2- (T1**2 - (2*t + 2) *T1+t**2 + 2*t +2)* np.exp(-t+T1))
        else:
            pz = dose_rate_val**2 * np.exp(-dose_rate_val) *((T2**2 - (2*t + 2)*T2 + t**2 + 2*t + 2) * np.exp(-t+T2) - (T1**2 - (2*t + 2)*T1 + t**2 + 2*t + 2) * np.exp(-t+T1))


    if t > T21:
        if dose_rate2[time.index(t)] > 0:
            pz += dose_rate_val2**2 * np.exp(-dose_rate_val2) * (2- (T21**2 - (2*t + 2) *T21+t**2 + 2*t +2)* np.exp(-t+T21))
        else:
            pz += dose_rate_val2**2 * np.exp(-dose_rate_val2) *((T22**2 - (2*t + 2)*T22 + t**2 + 2*t + 2) * np.exp(-t+T22) - (T21**2 - (2*t + 2)*T21 + t**2 + 2*t + 2) * np.exp(-t+T21))

    # pa.append(p)
    paz.append(pz)

# pa2.pop(0)
# print(max(pa))
plt.plot(time, paz, label="Pa")
# plt.plot(time, pa2, label="pa2")
# plt.plot(time, pa3, label="integral")
# plt.plot(time, paz, label="paz")

plt.plot(time, [dose_rate[i] + dose_rate2[i] for i in range(len(dose_rate))], label="dose")
plt.legend()
plt.show()
