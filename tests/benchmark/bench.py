import sys
sys.path.append('../data')

import numpy as np
from testcases import read_off
from orthogonalstreamlines import create_orthogonal_streamlines_mesh


ver, tri = read_off('../data/finemesh2d.off')
orient = np.fromfile('../data/finemesh2d.orient', sep=' ').reshape((-1, 3))

dx = 0.05
mesh = create_orthogonal_streamlines_mesh(ver, tri, orient, dx, 
                                          verbose=True)
