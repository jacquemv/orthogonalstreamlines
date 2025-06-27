import numpy as np
import matplotlib.pyplot as plt

error_msg = """
11: -8.18575 2.06207
31: -7.14853 2.58594
26: -7.37485 3.21746
16: -7.59353 3.74442
15: -8.15148 3.92484
24: -8.62031 4.80814
25: -8.86906 5.49058
4:  -9.74214 4.04289
14: -8.95985 3.62529
15: -8.15148 3.92484
10: -8.53055 2.65283
"""

x = np.array([float(s) for s in error_msg.replace(':', ' ').split()]).reshape((-1, 3))[:, 1:]

plt.plot(x[:, 0], x[:, 1])
for i in range(x.shape[0]):
    plt.text(x[i, 0], x[i, 1], str(i))
plt.show()