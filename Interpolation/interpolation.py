from data_structures import *

# LAGRANGE

def get_lagrange_basis_polynomial_j(x_j, xs):
    numerator = Polynomial([1])
    denominator = 1
    for x in xs:
        numerator = numerator * Polynomial([1, -x])
        denominator *= (x_j - x)
    return numerator * (1 / denominator)


def get_lagrange_basis_polynomials(points):
    polynomials = []
    for point in points:
        x = point.x
        rest = points[:]
        rest.remove(point)
        xs = [point.x for point in rest]
        polynomials.append(get_lagrange_basis_polynomial_j(x, xs))
    return polynomials


def lagrange_interpolation_function(points):
    result = Polynomial([0])
    lagrange_basis_polynomials = get_lagrange_basis_polynomials(points)
    for point, lagrange_basis_polynomial in zip(points, lagrange_basis_polynomials):
        result = result + (lagrange_basis_polynomial * point.y)
    return result




# SPLINE

def create_linear_equations_spline(points):
    # 2*(len-1) + 2*(len-2) + 2
    num_equations = 4*(len(points)-1)
    mat = SqMatrix(num_equations)
    y = Vector([0]*num_equations)
    for i in range(len(points[:-1])):
        h = points[i + 1].x - points[i].x

        # S_i(x_i) = y_i
        mat.rows[i*4][i*4 + 0] = 1
        y[i*4] = points[i].y

        # S_i(x_(i+1)) = y_(i+1)
        # a_i, b_i, c_i, d_i
        coefficients = [1, h, h**2, h**3]
        for j, c in enumerate(coefficients):
            mat.rows[i*4+1][i*4 + j] = c
        y[i*4+1] = points[i+1].y

        if i != len(points[:-1]) - 1:
            # S_i'(x_i) = S_(i+1)'(x_(i+1))
            coefficients = [0, 1, 2*h, 3*h**2]
            for j, c in enumerate(coefficients):
                mat.rows[i*4+2][i*4 + j] = c
            mat.rows[i*4+2][i*4 + 4 + 1] = -1
            y[i * 4 + 2] = 0

            # S_i''(x_i) = S_(i+1)''(x_(i+1))
            coefficients = [0, 0, 2 , 6 * h]
            for j, c in enumerate(coefficients):
                mat.rows[i * 4 + 3][i * 4 + j] = c
            mat.rows[i * 4 + 3][i * 4 + 4 + 2] = -2
            y[i * 4 + 3] = 0
    # S_0''(x0) = 0
    mat.rows[num_equations-2][2] = 1
    y[num_equations-2] = 0

    # S_(n-1)''(x_n) = 0
    h = points[-1].x - points[-2].x
    mat.rows[num_equations-1][num_equations-2] = 2
    mat.rows[num_equations-1][num_equations-1] = 6*h
    y[num_equations - 1] = 0
    return mat, y


def spline_interpolation_functions(points):
    A, b = create_linear_equations_spline(points)
    L, U, P = lu_fact_with_pivoting(A)
    x = solve_with_LU(L, U, P*b)

    interpolating_functions = {}
    for i in range(len(points[:-1])):
        rang = (points[i].x, points[i+1].x)
        pol = Polynomial(x[i*4:i*4+4][::-1])
        interpolating_functions[rang] = pol

    return interpolating_functions


def eval_for_spline_functions(functions, x):
    ranges = list(functions)
    func = None
    if x <= ranges[0][0]:
        func = functions[ranges[0]]
        return func.eval(x-ranges[0][0])
    elif x >= ranges[-1][1]:
        func = functions[ranges[-1]]
        return func.eval(x-ranges[-1][0])

    for rang in ranges:
        if x <= rang[1]:
            func = functions[rang]
            return func.eval(x-rang[0])
    return None
