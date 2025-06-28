import numpy as np
import matplotlib.pyplot as plt
from orthogonalstreamlines import create_orthogonal_streamlines_mesh

ver = np.array([[0, 0, 0], [1, 0, 0], [0.5, 1, 0]])
tri = np.array([[0, 1, 2]])
orient = np.array([[np.cos(0.2), np.sin(0.2), 0]])
mesh = create_orthogonal_streamlines_mesh(ver, tri, orient, 0.05)

plt.triplot(ver[:, 0], ver[:, 1], triangles=tri, c='0.75')
for cable in mesh.cables:
    plt.plot(mesh.vertices[cable, 0], mesh.vertices[cable, 1], c='C0')
plt.show()