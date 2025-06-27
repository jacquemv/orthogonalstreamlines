import sys
sys.path.append('../data')

import numpy as np
import matplotlib.pyplot as plt
from orthogonalstreamlines.intersection import create_cable_network, tri_normals, unpack_cables
from testcases import read_case, plot_case



data = read_case('../data/nonconnected1.case')
#plot_case(*data)
normal = tri_normals(*data[:2])

out = create_cable_network(normal, *data[-4:], remove_isolated_regions=True, cut_loose_ends=True)
X = out.vertices
cables = unpack_cables(out.cables, out.cables_len)
print(X.shape[0])
plt.plot(X[:, 0], X[:, 1], '.', c='C2')
for i, x in enumerate(X):
    plt.text(x[0], x[1], str(i))
for c in cables:
    plt.plot(X[c, 0], X[c, 1], c='C2', lw=2)

plt.show()
