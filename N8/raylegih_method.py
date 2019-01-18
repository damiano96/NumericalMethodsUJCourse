import numpy as np

e = 1e-8


def rayleigh(A):
    vectors_b = [np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]), np.array([0.0, 0.0, 1.0])]
    eigenvalues = []
    for b in vectors_b:
        l_new = np.dot(b.T, np.dot(A, b))
        while True:
            l_old = l_new
            w = np.linalg.solve(A - l_old * np.eye(3), b)
            b = w / np.linalg.norm(w)
            l_new = np.dot(b.T, np.dot(A, b))
            err = np.linalg.norm(w - l_new * b) / np.linalg.norm(w)
            if abs(1.0 - err) < e:
                break
        eigenvalues.append(l_new)
    return eigenvalues


A = np.array([[1.0, 2.0, 3.0], [2.0, 4.0, 5.0], [3.0, 5.0, -1.0]])
print(rayleigh(A))
