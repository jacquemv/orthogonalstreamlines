from itertools import chain
import numpy as np
import matplotlib.pyplot as plt


def read_off(file):
    with open(file, 'rt') as file:
        file.readline()
        nv, nt, _ = [int(x) for x in file.readline().split()]
        ver = np.fromfile(file, sep=' ', dtype=np.float32, count=nv*3)
        tri = np.fromfile(file, sep=' ', dtype=np.int32, count=nt*4)
    return ver.reshape((-1, 3)), tri.reshape((-1, 4))[:, 1:]

def write_case(fname, ver, tri, orient, lines1, faces1, lines2, faces2, fmt='%.5f'):
    with open(fname, 'wt') as file:
        file.write(f'{ver.shape[0]} {tri.shape[0]} {len(lines1)} {len(lines2)}\n')
        ver.tofile(file, sep=' ', format=fmt)
        file.write('\n')
        tri.tofile(file, sep=' ')
        file.write('\n')
        orient.tofile(file, sep=' ', format=fmt)
        file.write('\n')
        for line, face in chain(zip(lines1, faces1), zip(lines2, faces2)):
            line.tofile(file, sep=' ', format=fmt)
            file.write('\n')
            face.tofile(file, sep=' ')
            file.write('\n')
            
def read_case(fname):
    with open(fname, 'rt') as file:
        s = file.read().splitlines()
        nv, nt, nc1, nc2 = [int(k) for k in s[0].split()]
        assert len(s) == 4 + 2*(nc1 + nc2)
        ver = np.fromstring(s[1], dtype=np.float64, sep=' ').reshape((-1, 3))
        tri = np.fromstring(s[2], dtype=np.int32, sep=' ').reshape((-1, 3))
        orient = np.fromstring(s[3], dtype=np.float64, sep=' ').reshape((-1, 3))
        lines, faces = [], []
        for i in range(nc1+nc2):
            lines.append(np.fromstring(s[2*i+4], dtype=np.float64, sep=' ').reshape((-1, 3)))
            faces.append(np.fromstring(s[2*i+5], dtype=np.int32, sep=' '))
        lines1 = lines[:nc1]
        lines2 = lines[nc1:]
        faces1 = faces[:nc1]
        faces2 = faces[nc1:]
    return ver, tri, orient, lines1, faces1, lines2, faces2


def plot_case(ver, tri, orient, lines1, faces1, lines2, faces2):
    plt.triplot(ver[:, 0], ver[:, 1], tri, c='0.8')
    for x in lines1 + lines2:
        plt.plot(x[:, 0], x[:, 1], '.-', c='C0', alpha=0.3)
    plt.axis('equal')