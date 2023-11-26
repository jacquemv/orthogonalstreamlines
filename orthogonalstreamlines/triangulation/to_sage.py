import sys
import numpy as np

n = int(sys.argv[1])

t = np.arange(n) * 2*np.pi/n
x = np.round(np.column_stack((1000*np.cos(t), 1000*np.sin(t)))).astype(np.int32)

print("PointConfiguration.set_engine('internal')")
s = repr(x)[6:-14].replace(" ", "").replace("\n", "")
print("p = PointConfiguration("+s+")")
print("list(p.triangulations())")

# run the code on: https://sagecell.sagemath.org/

