import numpy as np
from scipy import linalg

e = 1e-8
L = 10.0
N = 1000
h = (2.0*L) / N
h_square = h * h

b = np.empty(N)
a = np.empty(N-1)

for i in range(0, N, 1):
    b[i] = ((2.0 / h_square) + 4 - (6.0 / np.cosh(-L + (i * h)) ** 2))

a.fill(-(1.0 / h_square))

e, v = linalg.eigh_tridiagonal(b, a, tol=e, lapack_driver='stebz')


for i in range(0, 4):
    print(e[i])

for i in range(0, 4):
    print(v[i])
