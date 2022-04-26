# a qbit is [a, b] where a and b are complex numbers
from math import atan, pi
from random import random

class comp(object):
    def __init__(self, a : float, b: float):
        self.a = a
        self.b = b
        self.r = pow(( pow(a, 2) + pow(b, 2) ), 0.5)
        try:
            self.theta = atan((b/a))
        except ZeroDivisionError:
            self.theta = pi/2 * (1 if self.b >= 0 else -1)
    
    def __repr__(self) -> str:
        return ("%.2f %s i%.2f" % (self.a, "+" if self.b >= 0 else "-", self.b if self.b >= 0 else (-1 * self.b)))
    
    def __add__(self, other):
        return comp(self.a + other.a, self.b + other.b)
    
    def __pow__(self, p):
        power = comp(1, 0)
        for each in range(p):
            power = power * self
        return power

    def __mul__(self, other):
        return comp((self.a * other.a) - (self.b * other.b), (self.a * other.b) + (self.b * other.a))
    
    def hack(self, other):
        return comp((self.a * other.a) - (self.b * other.b), (self.a * other.b) + (self.b * other.a))

    def __mul__(self, real : float):
        try:
            return comp(self.a * real, self.b * real)
        except TypeError:
            return self.hack(real)

class Matrix(object):
    def __init__(self, r : int, c : int):
        self.nrows = r
        self.ncols = c

        self.rows = [[comp(1, 0) for i in range(self.nrows)] for x in range(self.ncols)]
    #!todo:
    # def __mul__(self, other : Matrix):
    #     prod = Matrix(self.nrows, other.ncols)
    #     for row  in range(prod.nrows):
    #         for col in range(prod.ncols):
    #             prod.rows[row][col] = 

    def __pow__(self, vector):
        prod = []
        for each in range(self.nrows):
            dot = comp(0, 0)
            for i in range(self.ncols):
                dot += self.rows[each][i] * vector[i]
            prod.append(dot)
        return prod

    def __mul__(self, scalar):
        prod = Matrix(self.nrows, self.ncols)
        for each in range(self.nrows):
            for i in range(self.ncols):
                prod.rows[each][i]  = prod.rows[each][i] * scalar
        return prod
    
    def __repr__(self) -> str:
        string = ""
        for each in self.rows:
            string += each.__repr__() + "\n"
        return string.strip()

def mod(x, squared = False):
    if type(x) == type(comp(0, 0)):
        return pow(x.r, 2) if squared else x.r
    else:
        return (x if x >= 0 else (-1 * x)) if not squared else pow(x, 2)

def qbit(bit : int) -> list:
    return [0, 1] if bit else [1, 0]

def HGATE() -> Matrix:
    HADAMARDGATE = Matrix(2, 2)
    HADAMARDGATE = HADAMARDGATE * (1/pow(2, 0.5))
    HADAMARDGATE.rows[-1][-1] =  HADAMARDGATE.rows[-1][-1] * -1
    return HADAMARDGATE


def NOT(qbit: list) -> list:
    return qbit[::-1]

def IDEN(qbit : list) -> list:
    return qbit[::]

def HAD(qbit : list) -> list:
    # print(qbit)
    return HGATE() ** qbit

# def MEASURE(qbit : list) -> list:
#     if type(qbit[0]) == type(comp(0, 0)):
#         ar = pow(qbit[0].r, 2)
#     else:
#         ar = pow(qbit[0], 2)
#     if random() <= ar: return [1, 0]
#     return [0, 1]

def MEASURE(tensor: list) -> list:
    prev = 0
    ranges = list( map( lambda x: mod(x, True), tensor ) )
    m = random()
    res = [0 for i in range(len(tensor))]
    for each in range(len(tensor)):
        ranges[each] += prev
        prev = ranges[each]
        if m <= ranges[each]:
            res[each] = 1
            return res
    return res

def tensor(q1 : list, q2 : list) -> list:
    tensorproduct = []
    for each in q1:
        for every in q2:
            try:
                tensorproduct.append(each * every)
            except TypeError:
                tensorproduct.append(every * each)
    return tensorproduct

def CNOTGATE() -> Matrix:
    cg = Matrix(4, 4)
    cg = cg * 0
    for each in range(4):
        for i in range(4):
            if each == i: cg.rows[each][i] = comp(1, 0)
    buffer = cg.rows[-1]
    cg.rows[-1] = cg.rows[-2]
    cg.rows[-2] = buffer
    return cg

def CNOT(qcontrol : list, qtarget : list) ->list:
    tens = tensor(qcontrol, qtarget)
    result = CNOTGATE() ** tens
    return result