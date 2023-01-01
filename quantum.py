from math import atan, pi, log, cos, sin
from random import random
from typing import List
from copy import deepcopy
mple = True

try:
    import matplotlib.pyplot as plt
except ModuleNotFoundError:
    mple = False
    print("Module `matplotlib` is not installed, will use terminal workarounds for graphing instead.")


PI = pi
Pi = pi


class comp(object):
    """
    Class for dealing with complex numbers.
    """
    def __init__(self, a : float = 0, b : float = 0):
        """
        Defines a complex number z = a + ib
        To define a complex number in polar form, use comp.polar().
        If `a` is of type comp, it simply defines a new identical complex number.
        """
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

    def polar(r : float, theta : float):
        """
        Returns a comp in the given polar coordinates.
        Theta is measured in radians; use quantum.pi or Pi or PI for π
        """
        return comp(r * cos(theta), r * sin(theta))


    def __eq__(self, o: object) -> bool:
        if type(o) == int or type(o) == float: return ((self.a == o) and (self.b == 0))
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
        """
        Returns the comp form of whatever is passed in as input.
        If input is comp, the SAME object is returned (as opposed to comp(z), which gives a copy).
        """
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
        """
        Returns (1/z).
        """
        return comp( self.a/pow(self.r, 2) , (-self.b)/pow(self.r, 2))
    
    def conjugate(self):
        """
        Returns z̄.
        """
        conj = comp(self.a, (-1 * self.b))
        return conj

class Matrix(object):
    """
    Matrix of 2 dimensions with entries of type comp or float or int.
    A matrix on its own is a valid quantum gate that can be added to a qprogram via qprogram.addgates().

    Matrix.shape is the order of the matrix.
    Matrix.span is number of qbits the gate is supposed to input and output.
    """

    def __eq__(self, o: object) -> bool:
        if (self.nrows != o.nrows) or (self.ncols != o.ncols): return False
        for each in range(self.nrows):
            for i in range(self.ncols):
                if self.rows[each][i] != o.rows[each][i]: return False
        return True

    def __init__(self, r : int, c : int = -1, id : str = None):
        """
        Matrix(integer) returns a integer x integer Identity matrix.
        
        Matrix(a, b) returns a matrix of order a x b.
        
        Matrix(List[List[float or int or comp]]) converts a 2D matrix from list form to Matrix type.

        `id` is the gateid of the quantum gate (preferably 1 or a few chars at most).
        """
        if type(r) == list:
            assert r != [], "Empty list cannot be converted to a matrix."
            mat = r
            self.nrows = len(mat)
            self.ncols = len(mat[0])
            self.rows = mat
        else:
            self.nrows = r
            self.ncols = c if c > 0 else r
            self.rows = [[(comp(1) if i == x else comp()) for i in range(self.ncols)] for x in range(self.nrows)]
        self.gateid = id
        self.span = int(log(self.nrows, 2))
        self.shape = (self.nrows, self.ncols)

    def internaldot(a : list, b : list) -> float:
        assert len(a) == len(b)
        d = comp()
        for each in range(len(a)):
            d += (comp.getcomplex(a[each]) * b[each])
        return d

    def __mul__(self, other : object or list or float or int or comp):
        assert (
            (type(other) == Matrix) or (type(other) == list) or (type(other) == int) or (type(other) == comp) or (type(other) == float)
        ), f"Matrix cannot be multiplied by a {type(other)}!"
        if type(other) == Matrix:
            assert self.ncols == other.nrows, f"Dimensions of matrices don't match: {self.nrows}x{self.ncols} and {other.nrows}x{other.ncols}"
            product = Matrix(self.nrows, other.ncols)
            ot = other.transpose()
            for each in range(self.nrows):
                for i in range(other.ncols):
                    product.rows[each][i] = Matrix.internaldot(self.rows[each], ot.rows[i])
            return product
        elif type(other) == list:
            buffer = Matrix(1, len(other))
            buffer.rows[0] = other
            buffer = buffer.transpose()
            return (self * buffer)
        else:
            product = deepcopy(self)
            for each in range(product.nrows):
                for i in range(product.ncols):
                    product.rows[each][i] = product.rows[each][i] * other
            return product
    
    def __truediv__(self, other : float or int or comp):
        assert (
            (type(other) == int) or (type(other) == comp) or (type(other) == float)
        ), f"Matrix cannot be divided by a {type(other)}!"
        if type(other) == comp: return(self * other.inverse())
        return (self * (1/other))
    
    def determinant(matrix) -> float or comp:
        """
        Returns the determinant of this matrix.
        """
        assert matrix.nrows == matrix.ncols, f"Matrix must a square matrix; got dimensions {matrix.nrows}x{matrix.ncols}"
        if matrix.ncols == 1: return matrix.rows[0][0]
        if matrix.ncols == 2:
            return (matrix.rows[0][0] * matrix.rows[1][1]) - (matrix.rows[0][1] * matrix.rows[1][0])
        det = comp()
        for each in range(matrix.ncols):
            subdet = matrix.rows[0][each] * ((-1)**each)
            submat = [x[:each] + x[each + 1:] for x in matrix.rows[1:]]
            submatrix = Matrix(len(submat), len(submat[0]))
            submatrix.rows = submat
            # print(f'submatrix{each}\n', submatrix, '\n\n')
            subdet *= Matrix.determinant(submatrix)
            det += subdet
        return det

    def adjoint(matrix):
        """
        Returns the adjoint of this matrix.
        """
        assert matrix.nrows == matrix.ncols, f"Matrix must a square matrix; got dimensions {matrix.nrows}x{matrix.ncols}"
        adj = deepcopy(matrix)
        if matrix.nrows == 2:
            adj.rows[0][0] = matrix.rows[1][1]
            adj.rows[0][1] = matrix.rows[1][0] * -1
            adj.rows[1][0] = matrix.rows[0][1] * -1
            adj.rows[1][1] = matrix.rows[0][0]
            return adj.transpose()
        for each in range(adj.nrows):
            for i in range(adj.ncols):
                submat = [x[:i] + x[i + 1:] for x in (matrix.rows[:each] + matrix.rows[each + 1:])]
                submatrix = Matrix(len(submat), len(submat[0]))
                submatrix.rows = submat
                print(f"submatrix{each},{i}:")
                print(submatrix)
                print()
                adj.rows[each][i] = submatrix.determinant() * ((-1) ** (each + i))
        adj = adj.transpose()
        return adj
    
    def inverse(self):
        """
        Returns the inverse of this matrix.
        """
        return (self.adjoint() / self.determinant())


    def __repr__(self) -> str:
        string = ""
        for each in self.rows:
            string += each.__repr__() + "\n"
        return string.strip()
    
    def transpose(self):
        """
        Returns the transpose of this matrix.
        """
        trans = Matrix(self.ncols, self.nrows)
        for r in range(trans.nrows):
            for c in range(trans.ncols):
                trans.rows[r][c] = self.rows[c][r]
        return trans
    
    def conjugate(self):
        """
        Returns the elementwise conjugate of the matrix.
        """
        conj = Matrix(self.nrows, self.ncols)
        for r in range(conj.nrows):
            for c in range(conj.ncols):
                conj.rows[r][c] = self.rows[r][c].conjugate()
        return conj
    
    def dagger(matrix):
        """
        Returns the conjugate of the matrix.
        """
        return matrix.conjugate().transpose()
    
    def getlist(self) -> list:
        """
        If the matrix is a column or row matrix, returns a list form of the matrix.
        """
        if self.shape[0] == 1: return self.rows[0]
        elif self.shape[1] == 1: return self.transpose().rows[0]
        else: print(f"Cannot meaningfully convert to list; at least 1 dimension must be 1; shape : {self.shape}")

def mod(x, squared = False) -> float:
    """
    Returns the mod value of x.
    If squared is True, returns the square of the mod value.
    """
    if type(x) == type(comp(0, 0)):
        return pow(x.r, 2) if squared else x.r
    else:
        return (x if x >= 0 else (-1 * x)) if not squared else pow(x, 2)

def qbit(bit : int = 0) -> List[int]:
    """
    Returns a Qubit set to value bit.
    Qubit 0 is [1, 0], and 1 is [0, 1].
    """
    return [0, 1] if bit else [1, 0]


def qtensor(q1 : List[int or float or comp], q2 : List[int or float or comp]) -> List[int or float or comp]:
    """
    Returns the tensor product Q1⛒Q2.
    """
    tensorproduct = []
    for each in q1:
        for every in q2:
            tensorproduct.append(comp(each) * comp(every))
    return tensorproduct


def mtensor(m1 : Matrix, m2 : Matrix) -> Matrix:
    """
    Returns the Matrix tensor product M1⛒M2.
    """
    r = m1.nrows * m2.nrows
    c = m1.ncols * m2.ncols
    product = Matrix(r, c)

    for each in range(m1.nrows):
        for i in range(m1.ncols):
            for x in range(m2.nrows):
                for y in range(m2.ncols):
                    product.rows[(m2.nrows * each) + x][(m2.ncols * i) + y] = comp(m1.rows[each][i]) * comp(m2.rows[x][y])
    return product

def MEASURE(tensor: List[int or float or comp]) -> List[int]:
    """
    Performs a measurement on a system of qubits.
    `tensor` is a tensor in list form that represents the state of the system.
    Returns a sparse list with result of measurement.
    """
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
    """
    Prints the matrix array with only real parts of the entries.
    """
    thing = ""
    for each in matrix.rows:
        line = ""
        for i in each: line +=(str(i.a) + ' ')
        line += '\n'
        thing += line
    print(thing)

def printimg(matrix : Matrix) -> None:
    """
    Prints the matrix array with only imaginary parts of the entries.
    """
    thing = ""
    for each in matrix.rows:
        line = ""
        for i in each: line +=(str(i.b) + ' ')
        line += '\n'
        thing += line
    print(thing)


def plotmeasurement(measurement : list, binary = True, name : str = None) -> None:
    """
    Plots the measurement result graphically.
    `name` is the title for the graph.
    `binary` controls if the states are represented in binary or decimal.
    If matplotlib has not been installed, the terminal will automatically be used for plotting.
    """
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
    """
    Plots the measurement result graphically in terminal.
    `name` is the title for the graph.
    `binary` controls if the states are represented in binary or decimal.
    """
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

def run(shots : int, state : List[int or float or comp], binary : bool = True, graph : bool = False, terminal : bool = False, name : str = None) -> None:
    """
    Runs measurements on a state of the system.
    The state is measured `shots` number of times, and percentages of each outcome is recorded.
    If `graph` is set to True, the result is also plotted graphically.
    """
    
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

def extract(measurement : List[int], qbitindex : int) -> List[int]:
    """
    Returns the state of Qubit at index `qbitindex` in `measurement`.
    NOTE: returns qubits in list form ([0, 1] or [1, 0]).
    """
    nqbits = int(log(len(measurement), 2))
    index = bin(measurement.index(1))[2:]
    index = ('0' * (nqbits - len(index))) + index
    index = int(index[qbitindex])
    return ([0, 1] if index else [1, 0])

def validgate(gate : Matrix) -> bool:
    """
    Verifies if a Matrix object is a valid quantum gate.
    This aint working for some reason dont use it
    """
    print((gate.nrows == gate.ncols))
    print((gate * gate) == Matrix(gate.nrows))
    return (
        (gate.nrows == gate.ncols) and ((gate * gate) == Matrix(gate.nrows))
    )


HGATE = Matrix([
    [comp(1),  comp(1)],
    [comp(1), comp(-1)]
], id = 'h') * pow(2, -0.5)

IGATE = Matrix(2, id = 'i')

CNOT0 = Matrix(4, id = 'c0')
CNOT0.rows[2], CNOT0.rows[3] = CNOT0.rows[3], CNOT0.rows[2]

CNOT1 = Matrix(4, id = 'c1')
CNOT1.rows[1], CNOT1.rows[3] = CNOT1.rows[3], CNOT1.rows[1]

XGATE = Matrix([
    [comp(), comp(1)],
    [comp(1), comp()]
], id = 'x')

YGATE = Matrix([
    [comp(), comp(0, -1)],
    [comp(0, 1), comp( )]
], id = 'y')

ZGATE = Matrix([
    [comp(1), comp( )],
    [comp(), comp(-1)]
], id = 'z')

TGATE = Matrix([
    [comp(1), comp()],
    [comp(), comp.polar(1, pi/4)]
], id = 't')

SWAPGATE = Matrix(4, id = 'swap')
SWAPGATE.rows[1], SWAPGATE.rows[2] = SWAPGATE.rows[2], SWAPGATE.rows[1]

MTENSORIDENTITY = Matrix([
    [comp(1)]
], id = 'm')

TOFFOLIGATE = Matrix(8, id = 'tof')
TOFFOLIGATE.rows[6], TOFFOLIGATE.rows[7] = TOFFOLIGATE.rows[7], TOFFOLIGATE.rows[6]

def RGATE(angle : float) -> Matrix:
    """
    Returns a matrix of the form:

    [ sin(angle)   cos(angle) ]
    [ cos(angle)  -sin(angle) ]

    `angle` is in radians.
    """
    matrix = Matrix(2)
    matrix.rows[0][0] = comp(sin(angle), 0)
    matrix.rows[0][1] = comp(cos(angle), 0)
    matrix.rows[1][0] = comp(cos(angle), 0)
    matrix.rows[1][1] = comp(-sin(angle), 0)
    matrix.gateid = f'r({round(angle, 2)})'
    return matrix


class basischange:
    """
    A class for dealing with basis change.
    Instance of this object can be passed into qprogram for changing the basis.
    """
    def __init__(self, transformationmatrix : Matrix) -> None:
        """
        Transformation matrix is the matrix for changing basis of the system.
        Must be a 2x2 matrix.
        """
        r, c = transformationmatrix.nrows, transformationmatrix.ncols
        assert (r == 2) and (c == 2) , f"transformationmatrix must be a 2x2 matrix; got {r}x{c}"
        self.transmat = transformationmatrix
    
    def get(self, nqubitsinsystem : int = 1):
        """
        Returns the transformation matrix for an n-qubit system.
        """
        if nqubitsinsystem == 1: return self.transmat
        else:
            retter = deepcopy(self.transmat)
            for _ in range(nqubitsinsystem - 1):
                retter = mtensor(self.transmat, retter)
            return retter

class Block:
    pass

class qprogram(object):
    """
    Class for a Quantum Program.
    This class deals in the gates and circuits model of quantum computing.

    The normal pipeline for running a qprogram is to instantiate it, add gates, compile and then run.
    
    Example:
    `
    entangler = qprogram(2, name = 'quantum entanglement')
    entangler.addgates(0, [HGATE, CNOT0])
    entangler.compile()
    entangler.run(graph = True)
    `

    NOTE: It is assumed that measurement takes place at the end of the qprogram in the circuit diagram.
    """
    def __init__(self, nqbits : int, custom : list = None, bchange : basischange = None, name : str = None) -> None:
        """
        `nqbits` -> no. of qubits in the system.
        `custom` -> a list of integers that sets the initial state of each of the qubits (by default 0 for all).
        `bchange` -> a BasisChange object for changing the basis.
        `name` -> name of the Quantum program.
        """
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
        self.repr = None
        self.programmat = None

    def addgates(self, qbitindex : int, gates : list):
        """
        Adds gates to the qprogram.

        `qbitindex` is the index of the bit to which gates are added.
        `gates` is a list containing the sequence of gates to be added.
        """
        self.cache = None
        self.programmat = None
        self.gates[qbitindex] += gates
    
    def removegate(self, qbitindex : int, gateindex : int = -1) -> None:
        """
        Removes the gate at the specified position.
        
        If gate index is not provided, by default the last gate in the sequence is removed.
        """
        del self.gates[qbitindex][gateindex]

    def cleargates(self, qbitindex : int) -> None:
        """
        Removes all gates on the circuit line corresponding to Qubit at `qbitindex`.
        """
        self.gates[qbitindex] = []

    def __repr__(self) -> str:
        if self.repr is None : return f"The program {'' if self.name is None else self.name} is yet to be compiled"
        return self.repr
    
    def compile(self, verbose : bool = False, showcompilationresult : bool = True):
        """
        Compiles the qprogram. This step is necessary before running the qprogram, and everytime the circuit-gate diagram is changed.

        If `verbose`, debug info is showed as well.
        If `showcompilationresult`, the final block diagram of the program is displayed.
        """

        if showcompilationresult : print(f"\nCompiling {'program' if self.name is None else self.name}...")

        for each in self.gates:
            for i in each:
                if type(i) == Block : i.compile(showcompilationresult = False)

        longest = 0
        for each in self.gates:
            if len(each) > longest : longest = len(each)
        if verbose : print("longest set of gates :",longest)
        for each in range(len(self.gates)):
            self.gates[each] += [IGATE for _ in range(longest - len(self.gates[each]))]

        for each in range(len(self.gates)):
            for i in range(len(self.gates[each])):
                try:
                    for x in range(self.gates[each][i].span - 1):
                        x += 1
                        if self.gates[each + x][i].gateid != 'm': self.gates[each + x].insert(i, MTENSORIDENTITY)
                except IndexError:
                    print(f"Invalid position ({each}, {i}) for gate ({self.gates[each][i].gateid}) of span {self.gates[each][i].span}\nAborting compilation...")
                    return
        
        if verbose:
            for each in self.gates:
                for i in each:
                    print(i.gateid, end=' ')
                print()

        longest = 0
        for each in self.gates:
            if len(each) > longest : longest = len(each)
        if verbose : print("longest set of gates :",longest)
        for each in range(len(self.gates)):
            self.gates[each] += [IGATE for _ in range(longest - len(self.gates[each]))]        
        
        if verbose:
            print()
            for each in self.gates:
                for i in each:
                    print(i.gateid, end=' ')
                print()
        
        archways = []
        for each in range(longest):
            allidentities = True
            for i in range(len(self.gates)):
                allidentities = allidentities and (self.gates[i][each].gateid == 'i')
            if allidentities:
                if verbose: print(f"Skipping tensor product at {each} due to all I gates...")
                continue
            prod = MTENSORIDENTITY
            for i in range(len(self.gates)):
                prod = mtensor(prod, self.gates[i][each])
            archways.append(prod)
        
        if verbose:
            print("Archway shapes: ", end = "")
            for each in archways:
                print(each.shape, end = ' ')
        
        finalmat = Matrix(pow(2, len(self.gates)), pow(2, len(self.gates)))
        for each in archways[::-1]:
            finalmat = finalmat * each
        
        if verbose:
            print(f"\nFinal mat shape : {finalmat.shape}")
        
        self.programmat = finalmat

        self.calcrepr()

        if showcompilationresult:
            print(self.repr)
            print(f"\nCompilation{(' of ' + self.name) if self.name is not None else ''} complete!\n")

    
    
    def calcrepr(self):
        """
        Makes and returns the string representation of the block diagram of the qprogram.
        """
        string = f"\n{'' if self.name is None else self.name}\n"

        reprs = [f'q{i}({self.qbits[i].index(1)}) ⮕  ---' for i in range(len(self.gates))]

        for i in range(len(self.gates[0])):
            for each in range(len(self.gates)):
                gate = self.gates[each][i]
                # print('o', gate.gateid)
                if gate.gateid in ['i', 'm']:
                    reprs[each] += '--'
                    continue
                if gate.span == 1: reprs[each] += f'--[ {gate.gateid} ]--'
                else:
                    ori = len(reprs[each])
                    reprs[each] += ('-' * (len(max(reprs, key=len)) - ori))
                    ori = len(reprs[each])
                    reprs[each] += '⌈'
                    for _ in range(gate.span) : reprs[each] += f' {gate.gateid}'
                    reprs[each] += ' ⌉'
                    for x in range(gate.span - 2):
                        reprs[each + x + 1] += ('-' * (ori - len(reprs[each + x + 1])))
                        reprs[each + x + 1] += '|'
                        for _ in range(gate.span) : reprs[each + x + 1] += f' {gate.gateid}'
                        reprs[each + x + 1] += ' |'
                    reprs[each + gate.span - 1] += ('-' * (ori - len(reprs[each + gate.span - 1])))
                    reprs[each + gate.span - 1] += '⌊'
                    for _ in range(gate.span) : reprs[each + gate.span - 1] += f' {gate.gateid}'
                    reprs[each + gate.span - 1] += ' ⌋'

        longest = 0
        for each in reprs:
            if len(each) > longest: longest = len(each)

        for each in reprs:
            string += each
            string += ('-' * (longest - len(each)))
            string += '---\n'
        self.repr = string
    

    def measure(self) -> list:
        """
        Performs and returns the result of a single measurement.
        """
        state = self.qbits[0]
        for each in range(1, self.nqbits): state = qtensor(state, self.qbits[each])
        state = (self.programmat * state).getlist()

        return MEASURE(state)
    
    def run(self, shots : int = 1600, binary : bool  = True, graph : bool = False, terminal : bool = False) -> None:
        """
        Runs and records the outcome of the qprogram.

        `shots` -> No. of times the measurement of the final state of the system is repeated.
        `binary` -> Format of the numeric states displayed in the output.
        `graph` -> If set to True, the output of the qprogram is plotted graphically. May use matplotlib if installed.
        `terminal` -> Force all the output (including graphs) to be displayed only on the terminal window.
        """
        
        state = self.qbits[0]
        for each in range(1, self.nqbits): state = qtensor(state, self.qbits[each])

        state = (self.programmat * state).getlist()

        run(shots, state, binary=binary, graph=graph, terminal=terminal, name = self.name)
    
    def getblock(self, blockid : str = None) -> Block:
        """
        Converts the program to a block and returns it.

        If a block id is not provided, program name is used in its stead.
        """
        self.compile(showcompilationresult=False)
        block = Block(self.nqbits, self.name)
        block.programmat = deepcopy(self.programmat)
        return block


class Block(qprogram):
    """
    A group of gates' arrangement.

    Blocks can be used in exactly the same way as normal gates, and can contain Blocks themselves.

    A qprogram containing Blocks recursively and automatically compiles all the child blocks, hence blocks can be changed dynamically.
    """
    def __init__(self, nqbits: int, blockid: str) -> None:
        self.gateid = blockid
        super().__init__(nqbits, name = f"{blockid}_block")
    
    def __repr__(self) -> str:
        if self.repr is None : self.compile(showcompilationresult = False)
        return self.repr

    def compile(self, verbose: bool = False, showcompilationresult: bool = True):
        """
        Compiles the block (not necessary to be done by the user).
        """
        super().compile(verbose, showcompilationresult)
        self.span = self.programmat.span
        self.shape = self.programmat.shape
        self.ncols = self.programmat.ncols
        self.nrows = self.programmat.nrows
        self.rows = self.programmat.rows
    
    def getmat(self) -> Matrix:
        """
        Returns the plain corresponding to this Block.
        """
        self.compile(showcompilationresult = False)
        blockmat = self.programmat
        blockmat.gateid = self.blockid
        return blockmat

def controlledU(u : Matrix) -> Matrix:
    """
    Returns the controlled version of gate U, assuming the top qubit (q0) is the control qubit and U is applied on the
    subsequent set of qubits.
    """
    cu = Matrix(pow(2, u.span + 1), id = f"C_{u.gateid}")
    
    start = pow(2, u.span)

    for each in range(start):
        for i in range(start):
            cu.rows[each + start][i + start] = u.rows[each][i]

    return cu