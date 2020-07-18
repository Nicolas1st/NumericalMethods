def gaussian_elimination(A, b=None):

    dimensions = A.shape

    if b is not None:
        A = np.hstack((A, b))

    # nullify everything that is below the main diagonal
    for j in range(0, dimensions[0] - 1):
        for i in range(j + 1, A.shape[0]):
            coef = A[i, j] / A[j, j]
            A[i, :] -= coef * A[j, :]

    # nullify everything that is above the main diagonal
    for j in range(dimensions[0]-1, 0, -1):
        for i in range(j-1, -1, -1):
            coef = A[i, j] / A[j, j]
            A[i, :] -= coef * A[j, :]

    answers = []
    for i in range(dimensions[0]):
        answers.append(A[i, -1] / A[i, i])

    return answers
