import numpy as np

e = 1e-8

def qr_method(A):
    A_copy = A
    while True:
        x = A_copy.item((0, 0))
        Q, R = np.linalg.qr(A_copy)
        A_copy = R * Q
        if abs(x - A_copy.item(0, 0)) < e:
            break
    return A_copy.diagonal()


A = np.matrix([[1.0, 2.0, 3.0], [2.0, 4.0, 5.0], [3.0, 5.0, -1.0]])
print(qr_method(A))
