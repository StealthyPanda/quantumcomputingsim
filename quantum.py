# a qbit is [a, b] where a and b are complex numbers
from math import atan, pi, log, cos, sin
from random import random
from typing import Type
from copy import deepcopy
mple = True

try:
    import matplotlib.pyplot as plt
except ModuleNotFoundError:
    mple = False
    print("Module `matplotlib` is not installed, will use terminal workarounds for graphing instead.")


PI = pi
pi = pi

class comp(object):
    def __init__(self, a : float = 0, b : float = 0):
        if type(a) == comp:
            self.a, self.b = a.a, a.b
        else:
            self.a = a
            self.b = b
        self.r = pow(( pow(self.a, 2) + pow(self.b, 2) ), 0.5)
        try:
            self.theta = atan((self.b/self.a))
        except ZeroDivisionError:
            self.theta = pi/2 * (1 if self.b >= 0 else -1)

    #theta is measured in radians; use quantum.pi for π
    def polar(r : float, theta : float):
        return comp(r * cos(theta), r * sin(theta))


    def __eq__(self, o: object) -> bool:
        return ((self.a == o.a) and (self.b == o.b))
    
    def __repr__(self) -> str:
        return ("%.2f %s i%.2f" % (self.a, "+" if self.b >= 0 else "-", self.b if self.b >= 0 else (-1 * self.b)))
    
    def __add__(self, other : object):
        return comp(self.a + other.a, self.b + other.b)
    
    def __pow__(self, p : float):
        return comp.polar( pow(self.r, p) , self.theta * p)
    
    def __sub__(self, other : object):
        return (self + (other * -1))

    def getcomplex(number : object or float or int):
        if type(number) == comp: return number
        return comp(number)


    def __truediv__(self, other : object or float or int):
        assert (
            (type(other) == comp) or (type(other) == float) or (type(other) == int)
        ), f"Comp cannot be divided by a {type(other)}!"
        if type(other) == comp:
            return (self * other.inverse())
        else:
            return comp(self.a/other, self.b/other)

    def __mul__(self, other : object or float or int):
        assert (
            (type(other) == comp) or (type(other) == float) or (type(other) == int)
        ), f"Comp cannot be multiplied with a {type(other)}!"
        if type(other) == comp:
            return comp((self.a * other.a) - (self.b * other.b), (self.a * other.b) + (self.b * other.a))
        else:
            return comp(self.a * other, self.b * other)
    
    def inverse(self):
        return comp( self.a/pow(self.r, 2) , (-self.b)/pow(self.r, 2))
    
    def conjugate(self):
        conj = comp(self.a, (-1 * self.b))
        return conj

class Matrix(object):

    def __eq__(self, o: object) -> bool:
        if (self.nrows != o.nrows) or (self.ncols != o.ncols): return False
        for each in range(self.nrows):
            for i in range(self.ncols):
                if self.rows[each][i] != o.rows[each][i]: return False
        return True

    def __init__(self, r : int, c : int):
        self.nrows = r
        self.ncols = c
        self.gateid = None
        self.rows = [[(comp(1) if i == x else comp()) for i in range(self.ncols)] for x in range(self.nrows)]
        self.span = int(log(r, 2))

    # def __pow__(self, vector : list):
    #     prod = []
    #     for each in range(self.nrows):
    #         dot = comp(0, 0)
    #         for i in range(self.ncols):
    #             dot += self.rows[each][i] * vector[i]
    #         prod.append(dot)
    #     return prod

    # def __mul__(self, scalar : float):
    #     prod = Matrix(self.nrows, self.ncols)
    #     for each in range(self.nrows):
    #         for i in range(self.ncols):
    #             prod.rows[each][i]  = self.rows[each][i] * scalar
    #     return prod

    def internaldot(a : list, b : list) -> float:
        assert len(a) == len(b)
        d = 0
        for each in range(len(a)):
            d += (comp.getcomplex(a[each]) * b[each])
        return d

    def __mul__(self, other : object or list or float or int or comp):
        assert (
            (type(other) == Matrix) or (type(other) == list) or (type(other) == int) or (type(other) == comp) or (type(other) == float)
        ), f"Matrix cannot be multiplied by a {type(other)}!"
        if type(other) == Matrix:
            assert self.c == other.r, f"Dimensions of matrices don't match: {self.r}x{self.c} and {other.r}x{other.c}"
            product = Matrix(self.r, other.c)
            ot = other.transpose()
            for each in range(self.r):
                for i in range(other.c):
                    product[each][i] = Matrix.internaldot(self.rows[each], other.rows[i])
            return product
        elif type(other) == list:
            buffer = Matrix(1, len(other))
            buffer.rows[0] = other
            buffer = buffer.transpose()
            return (self * buffer)
        else:
            product = deepcopy(self)
            for each in range(product.r):
                for i in range(product.c):
                    product.rows[each][i] = product.rows[each][i] * other
            return product
    
    def __truediv__(self, other : float or int or comp):
        assert (
            (type(other) == int) or (type(other) == comp) or (type(other) == float)
        ), f"Matrix cannot be divided by a {type(other)}!"
        if type(other) == comp: return(self * other.inverse())
        return (self * (1/other))


    def __repr__(self) -> str:
        string = ""
        for each in self.rows:
            string += each.__repr__() + "\n"
        return string.strip()
    
    def transpose(self):
        trans = Matrix(self.ncols, self.nrows)
        for r in range(trans.nrows):
            for c in range(trans.ncols):
                trans.rows[r][c] = self.rows[c][r]
        return trans
    
    #this is elementwise conjugate of the matrix;
    def conjugate(self):
        conj = Matrix(self.nrows, self.ncols)
        for r in range(conj.nrows):
            for c in range(conj.ncols):
                conj.rows[r][c] = self.rows[r][c].conjugate()
        return conj

def mod(x, squared = False):
    if type(x) == type(comp(0, 0)):
        return pow(x.r, 2) if squared else x.r
    else:
        return (x if x >= 0 else (-1 * x)) if not squared else pow(x, 2)

def qbit(bit : int) -> list:
    return [0, 1] if bit else [1, 0]



def HGATE(m : int = 1) -> Matrix:
    if m == 0: return 1
    inner = HGATE(m-1)
    matrix = Matrix(int(pow(2, m)), int(pow(2, m)))
    matrix.gateid = 'h'
    if m == 1:
        matrix = matrix * (1/pow(2, 0.5))
        matrix.rows[-1][-1] =  matrix.rows[-1][-1] * -1
        matrix.gateid = 'h'
        return matrix
    for each in range(int(pow(2, m))):
        for i in range(int(pow(2, m))):
            matrix.rows[each][i] = inner.rows[each % int(pow(2, m-1))][i % int(pow(2, m-1))]
            if (pow(2, m-1) <= each) and (pow(2, m-1) <= i): matrix.rows[each][i] *= -1
    matrix = matrix * (1/pow(2, 0.5))
    matrix.gateid = 'h'
    return matrix

def NGATE(n : int = 1) -> Matrix:
    gate = Matrix(pow(2, n), pow(2, n))
    gate = gate * 0
    # buffer = Matrix(pow(2, n), pow(2, n))
    for each in range(pow(2, n)):
        for i in range(pow(2, n)):
            if (each + i) == (pow(2, n) - 1): gate.rows[each][i] = comp(1, 0)
    gate.gateid = 'n'
    return gate

def IGATE(n : int = 1) -> Matrix:
    gate = Matrix(pow(2, n), pow(2, n))
    gate = gate * 0
    for each in range(pow(2, n)):
        for i in range(pow(2, n)):
            if each == i:
                gate.rows[each][i] = comp(1, 0)
    gate.gateid = 'i'
    return gate

def SETTOGATE(value : int, n : int = 1):
    if n == 1:
        matrix = Matrix(2, 2) * 0
        matrix.rows[value][0] = comp(1, 0)
        matrix.rows[value][1] = comp(1, 0)
        return matrix
    else:
        return mtensor(SETTOGATE(value = value), SETTOGATE(value = value, n =  n - 1))

def NOT(qbit: list) -> list:
    return (NGATE(int(log(len(qbit), 2))) ** qbit)

def IDEN(qbit : list) -> list:
    return qbit

def HAD(qbit : list) -> list:
    # print(qbit)
    return HGATE(int(log(len(qbit), 2))) ** qbit

def FLIP(state: list) -> list:
    flipped = []
    for each in state:
        each = comp.getcomplex(each)
        each.b *= -1
        flipped.append(each)
    return flipped

def flipgate(gate : Matrix) -> Matrix:
    flipped = Matrix(gate.nrows, gate.ncols)
    flipped.gateid = gate.gateid
    flipped.rows = deepcopy(gate.rows)
    flipped.rows = flipped.rows[::-1]
    return flipped

def SHIFTGATE(state : list, phase : float) -> list:
    shifted = []
    rotor = comp.polar(1, phase)
    for each in state:
        shifted.append(comp.getcomplex(each) * rotor)
    return shifted

def SHIFT(phase : float) -> Type[lambda x: x]:
    return lambda x: SHIFTGATE(x, phase)


def RGATE(angle : float):
    matrix = Matrix(2, 2)
    matrix.rows[0][0] = comp(sin(angle), 0)
    matrix.rows[0][1] = comp(cos(angle), 0)
    matrix.rows[1][0] = comp(cos(angle), 0)
    matrix.rows[1][1] = comp(-sin(angle), 0)
    matrix.gateid = 'r'
    return matrix


def XGATE() -> Matrix:
    matrix = Matrix(2, 2)
    matrix.rows[0][0] = 0
    matrix.rows[1][1] = 0
    matrix.gateid = 'x'
    return matrix

def YGATE() -> Matrix:
    matrix = Matrix(2, 2) * comp(0, 1)
    matrix.rows[0][0] = 0
    matrix.rows[1][1] = 0
    matrix.rows[0][1] *= -1
    matrix.gateid = 'y'
    return matrix

def ZGATE() -> Matrix:
    matrix = IGATE()
    matrix.rows[1][1] *= -1
    matrix.gateid = 'z'
    return matrix


def PHASEGATE(phase : float) -> Matrix:
    matrixgate = Matrix(2, 2) * 0
    matrixgate.rows[0][0] = 1
    matrixgate.rows[1][1] = comp.polar(1, phase)
    matrixgate.gateid = 'phase'
    return matrixgate

def dagger(gate : Matrix) -> Matrix:
    return gate.conjugate().transpose()


def mtensor(m1 : Matrix, m2 : Matrix) -> Matrix:
    r = m1.nrows * m2.nrows
    c = m1.ncols * m2.ncols
    product = Matrix(r, c)

    for each in range(m1.nrows):
        for i in range(m1.ncols):
            for x in range(m2.nrows):
                for y in range(m2.ncols):
                    product.rows[(m2.nrows * each) + x][(m2.ncols * i) + y] = comp.getcomplex(m1.rows[each][i]) * comp.getcomplex(m2.rows[x][y])
    return product


class basischange:
    def __init__(self, transformationmatrix : Matrix) -> None:
        r, c = transformationmatrix.nrows, transformationmatrix.ncols
        assert (r == 2) and (c == 2) , f"transformationmatrix must be a 2x2 matrix; got {r}x{c}"
        self.transmat = transformationmatrix
    
    def get(self, nqubitsinsystem : int = 1):
        if nqubitsinsystem == 1: return self.transmat
        else:
            retter = deepcopy(self.transmat)
            for _ in range(nqubitsinsystem - 1):
                retter = mtensor(self.transmat, retter)
            return retter

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

def printreal(matrix : Matrix) -> None:
    thing = ""
    for each in matrix.rows:
        line = ""
        for i in each: line +=(str(i.a) + ' ')
        line += '\n'
        thing += line
    print(thing)

def printimg(matrix : Matrix) -> None:
    thing = ""
    for each in matrix.rows:
        line = ""
        for i in each: line +=(str(i.b) + ' ')
        line += '\n'
        thing += line
    print(thing)

def tensor(q1 : list, q2 : list) -> list:
    tensorproduct = []
    for each in q1:
        for every in q2:
            try:
                tensorproduct.append(each * every)
            except TypeError:
                tensorproduct.append(every * each)
    return tensorproduct

def CNOTGATEOLD(controlindex : int = 0, targetindex : int = 1) -> Matrix:
    cg = Matrix(4, 4)
    cg = cg * 0
    for each in range(4):
        for i in range(4):
            if each == i: cg.rows[each][i] = comp(1, 0)
    if controlindex == 0 and targetindex == 1:
        buffer = cg.rows[-1]
        cg.rows[-1] = cg.rows[-2]
        cg.rows[-2] = buffer
    elif controlindex == 1 and targetindex == 0:
        buffer = cg.rows[-1]
        cg.rows[-1] = cg.rows[1]
        cg.rows[1] = buffer
    cg.gateid = 'cnot'
    return cg

def CNOTGATE(nqbits : int = 2, controlindex : int = 0, targetindex : int = 1) -> Matrix:
    if controlindex == targetindex: raise BaseException("control and target cannot be the same")
    if nqbits <= 1: raise BaseException("nqbits cannot be less than 2")
    if nqbits == 2: return CNOTGATEOLD(controlindex = controlindex, targetindex = targetindex)
    t = mtensor(CNOTGATEOLD(controlindex = controlindex, targetindex = targetindex), CNOTGATE(nqbits - 1, controlindex= controlindex, targetindex=targetindex))
    t.gateid = 'cnot' + str(controlindex)
    return t

def CNOTR(controlbitindex : int = 0) -> Matrix:
    gate = Matrix(4, 4) * 0
    gate.gateid = 'cnotr'

    gate.rows[0][0] = comp(1, 0)
    gate.rows[1][1] = comp(1, 0)
    gate.rows[2][2] = comp(1, 0)
    gate.rows[3][3] = comp(1, 0)

    if controlbitindex == 0: gate.rows[2], gate.rows[3] = gate.rows[3], gate.rows[2]
    else: gate.rows[1], gate.rows[3] = gate.rows[3], gate.rows[1]

    return gate

def FLIPPEDCNOTGATE(nqbits : int = 2) -> Matrix:
    fcnot = Matrix(4, 4) * 0
    fcnot.rows[0][0] = comp(1, 0)
    fcnot.rows[1][3] = comp(1, 0)
    fcnot.rows[2][2] = comp(1, 0)
    fcnot.rows[3][1] = comp(1, 0)
    fcnot.gateid = 'fcnot'
    return fcnot


def CNOT(qcontrol : list, qtarget : list) ->list:
    tens = tensor(qcontrol, qtarget)
    result = CNOTGATE() ** tens
    return result


def plotmeasurement(measurement : list, binary = True, name : str = None) -> None:
    eigenvectors = [('Ψ' + (('0' * (int(log(len(measurement), 2)) - len(str(bin(i))[2:]))) + str(bin(i))[2:])) for i in range(len(measurement))]
    if not binary : eigenvectors = [('Ψ' + str(i)) for i in range(len(measurement))]
    if mple:
        plt.ylim(0, 100)
        plt.bar(eigenvectors, measurement)
        plt.xlabel("Eigenstates")
        plt.ylabel("Percentage of outcomes")
        if name is not None : plt.title(name)
        plt.show()
    else:
        plotinterminal(measurement=measurement, binary= binary)

def plotinterminal(measurement : list, binary = True, name : str = None) -> None:
    eigenvectors = [('Ψ' + (('0' * (int(log(len(measurement), 2)) - len(str(bin(i))[2:]))) + str(bin(i))[2:])) for i in range(len(measurement))]
    if not binary : eigenvectors = [('Ψ' + str(i)) for i in range(len(measurement))]
    graph = f"\n{'' if name is None else name}"
    if name is not None: graph += '\n'
    for each in range(len(eigenvectors)):
        graph += eigenvectors[each] + ' '
        graph += ("█" * int(measurement[each] / 2))
        if measurement[each] > 0: graph += (' ' + str(measurement[each]) + '%')
        graph += '\n'
    print(graph)
    return graph

def run(shots : int, state : list, binary : bool = True, graph : bool = False, terminal : bool = False, name : str = None) -> None:
    sl = int(log(len(state), 2))
    res = [0 for i in range(len(state))]
    for each in range(shots):
        measurement = MEASURE(state)
        res[measurement.index(1)] += 1
    ret = ""
    for each in range(len(res)):
        s = str(bin(each))[2:]
        s = ('0' * (sl - len(s))) + s
        s = int(int('0b'+ s, 2)) if not binary else s
        ret += (f"|Ψ{s}> : {float(res[each]/shots) * 100}%, ")
    ret = ret[:-2]
    print(ret)
    if graph:
        if not terminal: plotmeasurement([((each/shots) * 100) for each in res], binary = binary, name = name)
        else: plotinterminal([((each/shots) * 100) for each in res], binary = binary, name = name)

def extract(measurement : list, qbitindex : int) -> list:
    nqbits = int(log(len(measurement), 2))
    index = bin(measurement.index(1))[2:]
    index = ('0' * (nqbits - len(index))) + index
    index = int(index[qbitindex])
    return ([0, 1] if index else [1, 0])


reprs = {
    'r' : '[ R ]',
    'i' : '[ I ]',
    'n' : '[ ~ ]',
    'x' : '[ X ]',
    'y' : '[ Y ]',
    'z' : '[ Z ]',
    'h' : '[ H ]',
    'cnot' : '[ ⛒ ]',
    'fcnot' : '[ ⛒ ]',
    'phase' : '[ θ ]',
    'pin' : '[ ☉ ]'
}


def getrepr(gate : Matrix) -> str:
    try : return reprs[gate.gateid]
    except : 
        print(f"Gate represntation not found: {gate.gateid}")
        return '[ ! ]'

class qprogram(object):
    def __init__(self, nqbits : int, custom : list = None, bchange : basischange = None, name : str = None) -> None:
        self.qbits = []
        self.nqbits = nqbits
        self.cache = None
        if not custom: self.qbits = [qbit(0) for i in range(nqbits)]
        else:
            for each in custom:
                self.qbits.append(qbit(each))
        self.gates = [[] for i in range(nqbits)]
        self.bchange = bchange
        self.name = name
        self.repr = f"The program {'' if self.name is None else self.name} is yet to be compiled"

    def addgates(self, qbitindex : int, gates : list):
        self.cache = None
        self.gates[qbitindex] += gates

    def __repr__(self) -> str:
        return self.repr

    def compile(self, verbose : bool = False, showcompilationresult : bool = True):
        if showcompilationresult : print(f"\nCompiling {'program' if self.name is None else self.name}...")
        longest = 0
        for each in range(len(self.gates)):
            for i in range(len(self.gates[each])):
                if self.gates[each][i].gateid == 'cnot':
                    pin = IGATE()
                    pin.gateid = 'pin'
                    self.gates[each - 1].insert(i, pin)
                elif self.gates[each][i].gateid == 'fcnot':
                    pin = IGATE()
                    pin.gateid = 'pin'
                    self.gates[each + 1].insert(i, pin)
        for each in range(len(self.gates)):
            if longest < len(self.gates[each]): longest = len(self.gates[each])
        
        for each in range(len(self.gates)):
            self.gates[each] += [IGATE() for i in range(longest - len(self.gates[each]))]
        
        self.calcrepr()
        
        if showcompilationresult:
            print(self)
            print("Compilation complete!\n")

        if verbose and self.cache is not None: print(self.cache)

        self.cache = None
        return self

    def calcrepr(self):
        string = f"\n{'' if self.name is None else self.name}"
        if self.name is not None: string += '\n'
        for each in range(self.nqbits):
            line = f"q{str(each)} ({self.qbits[each][1]}) ⮕  ---"
            for i in range(len(self.gates[each])):
                line += getrepr(self.gates[each][i])
                line += "---"
            string += line + '\n'
        self.repr = string
    
    def getstate(self, verbose : bool = False, usecache : bool = True) -> list:
        if self.cache != None and usecache:
            print("Using cache...")
            return self.cache
        length = len(self.gates[0])

        state = self.qbits[0]
        for each in range(1, self.nqbits): state = tensor(state, self.qbits[each])

        for each in range(0, length):
            gates = []
            for i in range(self.nqbits): 
                if verbose: print("ngates", len(gates))
                if self.gates[i][each].gateid == 'cnot': 
                    if verbose: print('activated')
                    del gates[-1]
                if i == 0: gates.append(self.gates[i][each])
                elif self.gates[i - 1][each].gateid != 'fcnot': gates.append(self.gates[i][each])
            if verbose: print(len(gates))
            tens = gates[0]
            # printreal(tens)
            for each in range(1, len(gates)): 
                if verbose: printreal(tens)
                tens = mtensor(tens, gates[each])
            if verbose: printreal(tens)
            state = tens ** state
            if verbose: print(state)
        
        if self.bchange is not None: state = self.bchange.get(self.nqbits) ** state
        
        self.cache = state
        return state

    def measure(self, verbose : bool = False, usecache : bool = True) -> list:
        return MEASURE(self.getstate(verbose = verbose, usecache=usecache))
    
    def run(self, shots : int = 1600, verbose : bool = False, binary : bool  = True, graph : bool = False, usecache : bool = True, terminal : bool = False) -> None:
        state = self.getstate(verbose = verbose, usecache=usecache)
        
        run(shots, state, binary=binary, graph=graph, terminal=terminal, name = self.name)