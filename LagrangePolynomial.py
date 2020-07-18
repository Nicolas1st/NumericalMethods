from matplotlib import pyplot as plt
import random
import math


def w_prime_exp(points, for_point):
    product = 1
    for index, point in enumerate(points):
        if index != for_point:
            product *= points[for_point] - point
    return product


def w_exp(points, for_point, x):
    """divided by (x - xi)"""
    product = 1
    for index, point in enumerate(points):
        if index != for_point:
            product *= x - point
    return product


def lagrange_pol(func, points, x):
    total = 0
    for index, point in enumerate(points):
        total += func(point) / w_prime_exp(points, index) * w_exp(points, index, x)
    return total


def present_solution(true_func, points, evaluate_err_at=None):

    xs = [0.1*x for x in range(int(points[0]-1)*10, int(points[-1]+1)*10)]
    my_ys = [lagrange_pol(my_func, points, x) for x in xs]
    true_ys = [true_func(x) for x in xs]

    plt.plot(xs, my_ys, label='Interpolation polynomial')
    plt.plot(xs, true_ys, label='The function itself')
    plt.scatter(points, list(map(true_func, points)))
    plt.legend()

    if evaluate_err_at is not None:
        error = abs(newton_pol(my_func, points, evaluate_err_at) - my_func(evaluate_err_at))
        plt.title(f'The error at X* is equal to {round(error, 2)}, X* = {evaluate_err_at}.')

    plt.show()


def get_random_ascending_array(how_many, start, end):
    array = [random.uniform(start, end) for i in range(how_many)]
    return sorted(array)


def my_func(x):
    return x**5 - x**7 * math.cos(x)


present_solution(my_func, get_random_ascending_array(30, 1, 20))
