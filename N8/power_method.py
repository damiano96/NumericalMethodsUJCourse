import numpy as np

e = 1e-8

def power_method(A):
    eigenvalues = []
    A_copy = np.array(A)

    for _ in range(3):
        v = np.array([1, 1, 1])
        Av = A_copy.dot(v)
        ev = v.dot(Av)
        while True:
            Av = A_copy.dot(v)
            v_new = Av / np.linalg.norm(Av)
            Av = A_copy.dot(v_new)
            ev_new = v_new.dot(Av)
            if np.abs(ev - ev_new) < e:
                break
            v = v_new
            ev = ev_new
        largest_eigenvalue = ev
        eigenvalues.append(largest_eigenvalue)
        eigenvector = np.array([v_new])
        A_copy -= largest_eigenvalue*(eigenvector * eigenvector.T)

    return eigenvalues


A = np.array([[1.0, 2.0, 3.0], [2.0, 4.0, 5.0], [3.0, 5.0, -1.0]])
print(power_method(A))
