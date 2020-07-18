import numpy as np
import matplotlib.pyplot as plt
from GaussianElimination import gaussian_elimination


def form_system(xs, ys, degree):

    N = len(xs)
    first_equ = [[N] + [sum(xs**i) for i in range(1, degree+1)] + [sum(ys)]]
    other_equs = []

    for i in range(1, degree+1):
        equ = [sum(xs**p) for p in range(i, i + degree + 1)]
        equ.append(sum(ys * xs**i))
        other_equs.append(equ)

    system = np.array(first_equ + other_equs)

    return system


def polynomial(thetas, x):
    return sum([a * x**i for i, a in enumerate(thetas)])


def squared_errors(xs, ys, thetas):
    return sum((polynomial(thetas, x) - y)**2 for x, y in zip(xs, ys))


def present_solution(degree, xs, ys):

    system = form_system(xs, ys, degree)
    thetas = gaussian_elimination(system)

    xs_plotting = [0.2 * x for x in range(5 * (int(xs[0]) - 1), 5 * (int(xs[-1]) + 1))]
    ys_plotting = [polynomial(thetas, x) for x in xs_plotting]

    plt.plot(xs_plotting, ys_plotting, label=f'The power = {degree}, the error is {squared_errors(xs, ys, thetas)}')
    plt.scatter(xs, ys, color="red", alpha=0.3)

    plt.legend()
    plt.show()


xs = np.array([0.2, 0.7, 0.8, 1.6, 7, 10], dtype='float64')
ys = np.array([10, 4, 1, 0.6, 0.746, 0.676], dtype='float64')

present_solution(2, xs, ys)
