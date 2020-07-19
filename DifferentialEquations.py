""""
This code aims at solving second order differential equations
using four methods.
The substitution z = y' was used to convert the second order diff equ
into a system of two first order ones
"""


from matplotlib import pyplot as plt


interval = [3, 4]
h = 0.1
cond1 = 2
cond2 = 2

def true_func(x):
    return (x - 2)**3 + 1 / (x - 2)


def f(z):
    return z


def g(x, y, z):
    return ((x - 2) * z + 3 * y) / ((x - 2)**2)


def euler_method(f, g, interval, step, cond1, cond2):

    x_start = interval[0]
    y_start = cond1
    z_start = cond2
    y_true = true_func(x_start)

    # saves x, y, z, delta_y, delta_z, y_true, error
    vals = [[x_start, y_start, z_start, None, None, y_true, 0]]

    # number of points to compute values at
    points_overall = int((interval[1] - interval[0]) // step) + 1

    for point_number in range(points_overall):

        old_x = vals[-1][0]
        old_y = vals[-1][1]
        old_z = vals[-1][2]

        new_x = old_x + step
        new_y = old_y + step * f(old_z)
        new_z = old_z + step * g(old_x, old_y, old_z)

        delta_y = new_y - old_y
        delta_z = new_z - old_z
        true_y = true_func(new_x)
        error = abs(new_y - true_y)

        vals.append([new_x, new_y, new_z, delta_y, delta_z, true_y, error])

    return vals


def euler_cauchy_method(f, g, interval, step, cond1, cond2):

    x_start = interval[0]
    y_start = cond1
    z_start = cond2
    y_true = true_func(x_start)

    # saves x, y, z, y_bar, z_bar, delta_y, delta_z, y_true, error
    vals = [[x_start, y_start, z_start, None, None, None, None, y_true, 0]]

    # number of points to compute value at
    points_overall = int((interval[1] - interval[0]) // step) + 1

    for point_number in range(points_overall):

        old_x = vals[-1][0]
        old_y = vals[-1][1]
        old_z = vals[-1][2]

        new_x = old_x + step

        new_y_bar = old_y + step * f(old_z)
        new_z_bar = old_z + step * g(old_x, old_y, old_z)

        new_y = old_y + step/2 * (f(old_z) + f(new_z_bar))
        new_z = old_z + step/2 * (g(old_x, old_y, old_z) + g(new_x, new_y_bar, new_z_bar))

        delta_y = new_y - old_y
        delta_z = new_z - old_z
        true_y = true_func(new_x)
        error = abs(new_y - true_y)

        vals.append([new_x, new_y, new_z, new_y_bar, new_z_bar, delta_y, delta_z, true_y, error])

    return vals


def runge_kutt_method(f, g, interval, step, cond1, cond2, for_adams_method=False):

    x_start = interval[0]
    y_start = cond1
    z_start = cond2
    y_true = true_func(x_start)

    # saves x, y, z, delta_y, delta_z, y_true, error
    vals = [[x_start, y_start, z_start, None, None, y_true, 0]]

    # number of points to compute values at
    points_overall = int((interval[1] - interval[0]) // step) + 1
    if for_adams_method:
        points_overall = 4

    for point_number in range(points_overall):

        old_x = vals[-1][0]
        old_y = vals[-1][1]
        old_z = vals[-1][2]

        K1 = step * f(old_z)
        L1 = step * g(old_x, old_y, old_z)

        K2 = step * f(old_z + 1/2 * L1)
        L2 = step * g(old_x + step/2, old_y + K1/2, old_z + L1/2)

        K3 = step * f(old_z + L2/2)
        L3 = step * g(old_x + step/2, old_y + K2/2, old_z + L2/2)

        K4 = step * f(old_z + L3)
        L4 = step * g(old_x + step, old_y + K3, old_z + L3)

        delta_y = 1/6 * (K1 + 2*K2 + 2*K3 + K4)
        delta_z = 1/6 * (L1 + 2*L2 + 2*L3 + L4)

        new_x = old_x + step
        new_y = old_y + delta_y
        new_z = old_z + delta_z

        true_y = true_func(new_x)
        error = abs(new_y - true_y)

        vals.append([new_x, new_y, new_z, delta_y, delta_z, true_y, error])

    return vals


def adams_method(f, g, interval, step, cond1, cond2):
    # taking vals from the runge-kutt
    initial_points = runge_kutt_method(f, g, interval, step, cond1, cond2, for_adams_method=True)
    # z = y' from the first equation
    vals = []
    for point in initial_points:
        x = point[0]
        y = point[1]
        y_prime = point[2]
        y_double_prime = g(x, y, y_prime)
        true_y = point[-2]
        error = point[-1]
        vals.append([x, y, y_prime, y_double_prime, true_y, error])

    # number of points to compute value at
    points_overall = int((interval[1] - interval[0]) // step) - 3
    for point in range(points_overall):
        old_x = vals[-1][0]
        old_y = vals[-1][1]

        first_d1 = vals[-1][2]
        first_d2 = vals[-2][2]
        first_d3 = vals[-3][2]
        first_d4 = vals[-4][2]

        second_d1 = vals[-1][3]
        second_d2 = vals[-2][3]
        second_d3 = vals[-3][3]
        second_d4 = vals[-4][3]

        new_x = old_x + step
        new_y = old_y + step/24 * (55*first_d1 - 59*first_d2 + 37 * first_d3 - 9*first_d4)
        new_y_prime = first_d1 + step/24 * (55*second_d1 - 59*second_d2 + 37 * second_d3 - 9*second_d4)
        new_y_double_prime = g(new_x, new_y, new_y_prime)
        true_y = true_func(new_x)
        error = abs(true_y - new_y)

        vals.append([new_x, new_y, new_y_prime, new_y_double_prime, true_y, error])

    return vals


def get_points(points):
    xs = []
    ys = []
    for iteration in points:
        xs.append(iteration[0])
        ys.append(iteration[1])
    return xs, ys


def get_true_points(points):
    xs = []
    true_ys = []
    for point in points:
        xs.append(point[0])
        true_ys.append(point[-2])
    return xs, true_ys


def runge_romberg_error(method, p, k=2):
    h1 = method(f, g, interval, h, cond1, cond2)
    h2 = method(f, g, interval, h/2, cond1, cond2)

    errors = []
    for i in range(len(h1)):
        true_y = h1[i][-2]
        error = true_y - h1[i][1] + (h1[i][1] - h2[i*2][1]) / (k**p - 1)
        errors.append(error)
    return errors


points1 = euler_method(f, g, interval, h, cond1, cond2)
points2 = euler_cauchy_method(f, g, interval, h, cond1, cond2)
points3 = runge_kutt_method(f, g, interval, h, cond1, cond2)
points4 = adams_method(f, g, interval, h, cond1, cond2)

print('Euler\'s method')
runge_romberg_errors = runge_romberg_error(euler_method, 1)
for i, point in enumerate(points1):
    print('i', i)
    print('x', point[0])
    print('y', point[1])
    print('z', point[2])
    print('delta_y', point[3])
    print('delta_z', point[4])
    print('y_true', point[5])
    print('error', point[6])
    print('Runge-Rombergs error', runge_romberg_errors[i])
    print()

print('\n\n')

print('Euler-Cauchy\'s method')
runge_romberg_errors = runge_romberg_error(euler_cauchy_method, 2)
for i, point in enumerate(points2):
    print('i', i)
    print('x', point[0])
    print('y', point[1])
    print('z', point[2])
    print('y_bar', point[3])
    print('z_bar', point[4])
    print('delta_y', point[5])
    print('delta_z', point[6])
    print('y_true', point[7])
    print('error', point[8])
    print('Runge-Rombers error', runge_romberg_errors[i])
    print()

print('\n\n')

print('Runge-Kutt\' method')
runge_romberg_errors = runge_romberg_error(runge_kutt_method, 4)
for i, point in enumerate(points3):
    print('i', i)
    print('x', point[0])
    print('y', point[1])
    print('z', point[2])
    print('delta_y', point[3])
    print('delta_z', point[4])
    print('y_true', point[5])
    print('error', point[6])
    print('Runge-Rombergs error', runge_romberg_errors[i])
    print()

print('\n\n')

print('Adams method')
runge_romberg_errors = runge_romberg_error(adams_method, 4)
for i, point in enumerate(points4):
    print('i', i)
    print('x', point[0])
    print('y', point[1])
    print('y', point[2])
    print('y', point[3])
    print('y_true', point[4])
    print('error', point[5])
    print('Runge-Rombergs error', runge_romberg_errors[i])
    print()

plt.plot(*get_true_points(points1), label='The Truth')
plt.plot(*get_points(points1), label='Euler')
plt.plot(*get_points(points2), label='Euler-Cauchy')
plt.plot(*get_points(points3), label='Runge_Kutt')
plt.plot(*get_points(points4), label='Adams\'')
plt.legend()
plt.show()
