def rectangle_integral(func, start, end, step, side):
    
    ys = []

    if side != 'middle':
        if side == 'right':
            start += step
        elif side == 'left':
            end -= step
        x = start
        while x <= end:
            y = func(x)
            ys.append(y)
            x += step

    else:
        x1 = start
        x2 = start + step
        while x2 <= end:
            y = func((x1 + x2) / 2)
            ys.append(y)
            x1 += step
            x2 += step

    return step * sum(ys)


def trapezoid_integral(func, start, end, step, side=None):

    ys = []
    x = start

    while x <= end:
        y = func(x)
        ys.append(y)
        x += step

    ys[0] /= 2
    ys[-1] /= 2

    return step * sum(ys)


def simpson_integral(func, start, end, step, side=None):

    ys = []
    x = start

    while x <= end:
        y = func(x)
        ys.append(y)
        x += step

    multipliers = [2, 4]
    for i in range(1, len(ys)-1):
        ys[i] *= multipliers[i % 2]

    return step/3 * sum(ys)


def compute_error(integral1, integral2, step1, step2, p, true_value=0):
    k = step1 / step2
    error = true_value - integral1
    error -= (integral1 - integral2) / (k**p - 1)
    return error
