from matplotlib import pyplot as plt
import math
import random


def get_random_ascending_array(how_many, start, end):
    array = [random.uniform(start, end) for i in range(how_many)]
    return sorted(array)


def differences(func, points):

    # these are differences from the Newtons' polynomial formula
    first_level_diffs = [[func(point)] for point in points]
    second_level_diffs = []
    for i in range(len(first_level_diffs) - 1):
        numerator = first_level_diffs[i + 1][0] - first_level_diffs[i][0]
        denominator = points[i + 1] - points[i]
        new_diff = numerator / denominator
        second_level_diffs.append([new_diff, points[i + 1], points[i]])
    divided_diffs = [first_level_diffs, second_level_diffs]

    while len(divided_diffs) != len(points):
        next_level_diffs = []

        for i in range(len(divided_diffs[-1]) - 1):
            prev_level = divided_diffs[-1]
            numerator = prev_level[i + 1][0] - prev_level[i][0]
            denominator = prev_level[i + 1][1] - prev_level[i][2]
            new_diff = numerator / denominator
            next_level_diffs.append([new_diff, prev_level[i + 1][1], prev_level[i][2]])

        divided_diffs.append(next_level_diffs)
    divided_diffs = list(map(lambda x: x[0][0], divided_diffs))

    return divided_diffs


def products(points, x):
    products = [1]
    product = 1
    for i in range(len(points)-1):
        product *= x - points[i]
        products.append(product)
    return products


def newton_pol(func, points, x):
    parts = list(zip(differences(func, points), products(points, x)))
    res = 0
    for part in parts:
        res += part[0] * part[1]
    return res


def present_solution(true_func, points, evaluate_err_at=None):

    xs = [0.1*x for x in range(int(points[0]-1)*10, int(points[-1]+1)*10)]
    my_ys = [newton_pol(my_func, points, x) for x in xs]
    true_ys = [true_func(x) for x in xs]

    plt.plot(xs, my_ys, label='Interpolation polynomial')
    plt.plot(xs, true_ys, label='The function itself')
    plt.scatter(points, list(map(true_func, points)))
    plt.legend()

    if evaluate_err_at is not None:
        error = abs(newton_pol(my_func, points, evaluate_err_at) - my_func(evaluate_err_at))
        plt.title(f'The error at X* is equal to {round(error, 2)}, X* = {evaluate_err_at}.')

    plt.show()


# here you can experiment as you wish)
def my_func(x):
    return x**2 + x**3 - x + x**2 * math.sin(x)


present_solution(my_func, get_random_ascending_array(30, 1, 25))
