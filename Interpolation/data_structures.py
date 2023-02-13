import copy

class Polynomial:
    def __init__(self, coefficients):
        self.coefficients = coefficients
        self.degree = len(coefficients)-1

    def __add__(self, other):
        if isinstance(other, Polynomial):
            result: []
            to_add: []
            if self.degree >= other.degree:
                result = self.coefficients[::-1]
                to_add = other.coefficients[::-1]
            else:
                result = other.coefficients[::-1]
                to_add = self.coefficients[::-1]

            for idx, c in enumerate(to_add):
                result[idx] += c

            result.reverse()
            return Polynomial(result)
        else:
            return NotImplemented

    def __mul__(self, other):
        if isinstance(other, Polynomial):
            product_coeff = [0] * (self.degree+other.degree+1)

            for i in range(self.degree+1):
                for j in range(other.degree+1):
                    product_coeff[i+j] += self.coefficients[i] * other.coefficients[j]
            return Polynomial(product_coeff)
        elif isinstance(other, float) or isinstance(other, int):
            product_coeff = [other*c for c in self.coefficients]
            return Polynomial(product_coeff)
        else:
            return NotImplemented

    def __rmul__(self, other):
        return self.__mul__(other)

    def __repr__(self):
        return repr(self.coefficients)

    def eval(self, x):
        evaluation = 0
        for degree, coeff in enumerate(self.coefficients[::-1]):
            evaluation += x**degree * coeff
        return evaluation

class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __repr__(self):
        return '({x}, {y})'.format(x=self.x, y=self.y)


class Vector():
    def __init__(self, l: list):
        self.l = l.copy()

        # how many numbers to print before replacing them with dots
        self.print_befor_dots = 5
        # how wide will the numbers be
        self.num_len = 6

    def __len__(self):
        return len(self.l)

    def __getitem__(self, item):
        return self.l[item]

    def __setitem__(self, key, value):
        self.l[key] = value

    def __neg__(self):
        return Vector([-elem for elem in self.l])

    def __add__(self, other):
        if isinstance(other, Vector) and len(self) == len(other):
            return Vector([e1 + e2 for e1, e2 in zip(self.l, other.l)])
        else:
            raise ValueError("Addition operand is not a proper vector")

    def __sub__(self, other):
        if isinstance(other, Vector) and len(self) == len(other):
            return self + (-other)
        else:
            raise ValueError("Substraction operand is not a proper vector")

    # dot product
    def __mul__(self, other):
        if isinstance(other, Vector) and len(self) == len(other):
            vect = other
            sum = 0
            for i in range(len(self.l)):
                sum += self[i] * vect[i]
            return sum
        else:
            raise ValueError("Dot product operand is not a proper vector")

    def __str__(self):
        l_len = len(self.l)
        if l_len < self.print_befor_dots * 2 + 1:
            string = '| '
            for i in range(l_len):
                round_to = self.num_len - 2 if float(self.l[i]) > 0 else self.num_len - 3
                fl_str = str(round(float(self.l[i]), round_to))
                zeroes = '0' * (self.num_len - len(fl_str))
                fl_str = fl_str + zeroes
                string += fl_str
                if i < l_len - 1:
                    string += ', '
            string += ' |'
            return string
        string = '| '
        for i in range(self.print_befor_dots):
            round_to = self.num_len - 2 if float(self.l[i]) > 0 else self.num_len - 3
            fl_str = str(round(float(self.l[i]), round_to))
            zeroes = '0' * (self.num_len - len(fl_str))
            fl_str = fl_str + zeroes
            string += fl_str
            string += ', '
        string += '... '
        for i in range(l_len - self.print_befor_dots, l_len):
            round_to = self.num_len - 2 if float(self.l[i]) > 0 else self.num_len - 3
            fl_str = str(round(float(self.l[i]), round_to))
            zeroes = '0' * (self.num_len - len(fl_str))
            fl_str = fl_str + zeroes
            string += fl_str
            if i < l_len - 1:
                string += ', '
        string += ' |'
        return string


class SqMatrix():
    def __init__(self, n: int):
        self.rows = [Vector([0 for _ in range(n)]) for _ in range(n)]
        self.n = n

        # how many rows to print before replacing them with dots
        self.print_befor_dots = 5

    def __getitem__(self, item: int):
        return self.rows[item]

    def get_column(self, colnum: int):
        column = []
        for row in self.rows:
            column.append(row[colnum])
        return Vector(column)

    def __mul__(self, other):
        if isinstance(other, Vector) and len(other) == self.n:
            vect = other
            result = []
            for row in self.rows:
                curr_result = 0
                for idx, elem_row in enumerate(row):
                    curr_result += vect[idx] * elem_row
                result.append(curr_result)
            return Vector(result)

        elif isinstance(other, SqMatrix) and other.n == self.n:
            n = self.n
            result = SqMatrix(n)
            for rowidx in range(n):
                for columnidx in range(n):
                    row = self.rows[rowidx]
                    column = other.get_column(columnidx)
                    # dot product
                    result[rowidx][columnidx] = row * column
            return result

        else:
            raise ValueError("Multiplication operand is not a proper vector/matrix")

    def __str__(self):
        string = ''
        l_rows = len(self.rows)
        if (l_rows < self.print_befor_dots * 2 + 1):
            for row in self.rows:
                string += row.__str__() + '\n'
            return string
        for i in range(self.print_befor_dots):
            string += self.rows[i].__str__() + '\n'
        len_row = len(self.rows[0].l)
        bef_dots_row = self.rows[0].print_befor_dots
        num_len = self.rows[0].num_len
        if len_row < 2 * bef_dots_row + 1:
            string += '| '
            for i in range(len_row):
                elem = ' ' * ((num_len + 1) // 2 - 1)
                elem += ':'
                elem += ' ' * (num_len // 2)
                # accounts for ', '
                elem += '  '
                string += elem
            # removes last two spaces that were suppose to be below ', '
            string = string[:-2]
        else:
            string += '| '
            for _ in range(bef_dots_row):
                elem = ' ' * ((num_len + 1) // 2 - 1)
                elem += ':'
                elem += ' ' * (num_len // 2)
                # accounts for ', '
                elem += '  '
                string += elem
            # accounts for '... '
            string += ' :  '
            for _ in range(bef_dots_row):
                elem = ' ' * ((num_len + 1) // 2 - 1)
                elem += ':'
                elem += ' ' * (num_len // 2)
                # accounts for ', '
                elem += '  '
                string += elem
            # removes last two spaces that were suppose to be below ', '
            string = string[:-2]
        string += ' |\n'

        for i in range(l_rows - self.print_befor_dots, l_rows):
            string += self.rows[i].__str__() + '\n'
        return string


def lu_fact_with_pivoting(A):
    n = A.n
    U = copy.deepcopy(A)

    L = SqMatrix(n)
    for i in range(n):
        L[i][i] = 1
    P = SqMatrix(n)
    for i in range(n):
        P[i][i] = 1

    for k in range(n-1):
        # find pivot
        idx = k
        max_val = U.rows[k][k]
        for i in range(k, n):
            if U.rows[i][k] > max_val:
                max_val = U.rows[i][k]
                idx = i

        # interchange rows
        if k != idx:
            # change U
            new_values_k_U = [val for val in U.rows[idx][k:]]
            old_values_k_U = [val for val in U.rows[k][k:]]
            for i in range(k, n):
                U.rows[k][i] = new_values_k_U[i-k]
            for i in range(k, n):
                U.rows[idx][i] = old_values_k_U[i-k]

            #change L
            new_values_k_L = [val for val in L.rows[idx][:k]]
            old_values_k_L = [val for val in L.rows[k][:k]]
            for i in range(0, k):
                L.rows[k][i] = new_values_k_L[i]
            for i in range(0, k):
                L.rows[idx][i] = old_values_k_L[i]

            #change P
            k_row_copy = copy.deepcopy(P.rows[k])
            P.rows[k] = copy.deepcopy(P.rows[idx])
            P.rows[idx] = k_row_copy

        # standard LU
        for j in range(k+1, n):
            L[j][k] = U[j][k]/U[k][k]
            for z in range(k, n):
                U[j][z] = U[j][z] - L[j][k]*U[k][z]

    return L, U, P

def forward_substitution(L, b):
    n = L.n
    x = Vector([0]*n)
    for i in range(n):
        x[i] = b[i]
        for j in range(i):
            x[i] = x[i] - L[i][j] * x[j]
    return x

def backward_substitution(U, b):
    n = U.n
    x = Vector([0]*n)

    # iterates for i in <n-1, 0>
    for i in range(n-1, -1, -1):
        x[i] = b[i]
        #iterates for j in <n-1, i)
        for j in range(n-1, i, -1):
            x[i] = x[i] - U[i][j]*x[j]
        x[i] = x[i]/U[i][i]
    return x

def solve_with_LU(L,U, b):
    y = forward_substitution(L, b)
    return backward_substitution(U, y)
