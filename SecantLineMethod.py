def will_converge(equation, a, b):
    if equation(a) * equation(b) < 0:
        return True
    return False


def secant_lines_method(equation, x1, x2, precision, iterations):
    if not will_converge(equation, x1, x2):
        print("Failed to converge")
        return None

    iteration = 0
    for i in range(iterations):

        a = equation(x1)
        b = equation(x2)
        next_x = x2 - b / ((b - a) / (x2 - x1))
        x1, x2 = x2, next_x

        err = abs(x2 - x1)
        print(f'Итерация{iteration}, x: {next_x}, err: {err}')

        if abs(x2 - x1) < precision:
            return next_x

    print("Converges too slow")
    return None
