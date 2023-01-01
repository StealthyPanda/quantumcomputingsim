def legacy_HGATE(m : int = 1) -> Matrix:
    if m == 0: return 1
    inner = legacy_HGATE(m-1)
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

def legacy_NGATE(n : int = 1) -> Matrix:
    gate = Matrix(pow(2, n), pow(2, n))
    gate = gate * 0
    # buffer = Matrix(pow(2, n), pow(2, n))
    for each in range(pow(2, n)):
        for i in range(pow(2, n)):
            if (each + i) == (pow(2, n) - 1): gate.rows[each][i] = comp(1, 0)
    gate.gateid = 'n'
    return gate

def legacy_IGATE(n : int = 1) -> Matrix:
    gate = Matrix(pow(2, n), pow(2, n))
    gate = gate * 0
    for each in range(pow(2, n)):
        for i in range(pow(2, n)):
            if each == i:
                gate.rows[each][i] = comp(1, 0)
    gate.gateid = 'i'
    return gate

def legacy_SETTOGATE(value : int, n : int = 1):
    if n == 1:
        matrix = Matrix(2, 2) * 0
        matrix.rows[value][0] = comp(1, 0)
        matrix.rows[value][1] = comp(1, 0)
        return matrix
    else:
        return mtensor(legacy_SETTOGATE(value = value), legacy_SETTOGATE(value = value, n =  n - 1))

def legacy_NOT(qbit: list) -> list:
    return (legacy_NGATE(int(log(len(qbit), 2))) ** qbit)

def legacy_IDEN(qbit : list) -> list:
    return qbit

def legacy_HAD(qbit : list) -> list:
    # print(qbit)
    return legacy_HGATE(int(log(len(qbit), 2))) ** qbit

def legacy_FLIP(state: list) -> list:
    flipped = []
    for each in state:
        each = comp.getcomplex(each)
        each.b *= -1
        flipped.append(each)
    return flipped

def legacy_flipgate(gate : Matrix) -> Matrix:
    flipped = Matrix(gate.nrows, gate.ncols)
    flipped.gateid = gate.gateid
    flipped.rows = deepcopy(gate.rows)
    flipped.rows = flipped.rows[::-1]
    return flipped

def legacy_SHIFTGATE(state : list, phase : float) -> list:
    shifted = []
    rotor = comp.polar(1, phase)
    for each in state:
        shifted.append(comp.getcomplex(each) * rotor)
    return shifted

def legacy_SHIFT(phase : float) -> Type[lambda x: x]:
    return lambda x: legacy_SHIFTGATE(x, phase)


def legacy_RGATE(angle : float):
    matrix = Matrix(2, 2)
    matrix.rows[0][0] = comp(sin(angle), 0)
    matrix.rows[0][1] = comp(cos(angle), 0)
    matrix.rows[1][0] = comp(cos(angle), 0)
    matrix.rows[1][1] = comp(-sin(angle), 0)
    matrix.gateid = 'r'
    return matrix


def legacy_XGATE() -> Matrix:
    matrix = Matrix(2, 2)
    matrix.rows[0][0] = 0
    matrix.rows[1][1] = 0
    matrix.gateid = 'x'
    return matrix

def legacy_YGATE() -> Matrix:
    matrix = Matrix(2, 2) * comp(0, 1)
    matrix.rows[0][0] = 0
    matrix.rows[1][1] = 0
    matrix.rows[0][1] *= -1
    matrix.gateid = 'y'
    return matrix

def legacy_ZGATE() -> Matrix:
    matrix = legacy_IGATE()
    matrix.rows[1][1] *= -1
    matrix.gateid = 'z'
    return matrix


def legacy_PHASEGATE(phase : float) -> Matrix:
    matrixgate = Matrix(2, 2) * 0
    matrixgate.rows[0][0] = 1
    matrixgate.rows[1][1] = comp.polar(1, phase)
    matrixgate.gateid = 'phase'
    return matrixgate

def legacy_CNOTGATEOLD(controlindex : int = 0, targetindex : int = 1) -> Matrix:
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

def legacy_CNOTGATE(nqbits : int = 2, controlindex : int = 0, targetindex : int = 1) -> Matrix:
    if controlindex == targetindex: raise BaseException("control and target cannot be the same")
    if nqbits <= 1: raise BaseException("nqbits cannot be less than 2")
    if nqbits == 2: return legacy_CNOTGATEOLD(controlindex = controlindex, targetindex = targetindex)
    t = mtensor(legacy_CNOTGATEOLD(controlindex = controlindex, targetindex = targetindex), legacy_CNOTGATE(nqbits - 1, controlindex= controlindex, targetindex=targetindex))
    t.gateid = 'cnot' + str(controlindex)
    return t

def legacy_CNOTR(controlbitindex : int = 0) -> Matrix:
    gate = Matrix(4, 4) * 0
    gate.gateid = 'cnotr'

    gate.rows[0][0] = comp(1, 0)
    gate.rows[1][1] = comp(1, 0)
    gate.rows[2][2] = comp(1, 0)
    gate.rows[3][3] = comp(1, 0)

    if controlbitindex == 0: gate.rows[2], gate.rows[3] = gate.rows[3], gate.rows[2]
    else: gate.rows[1], gate.rows[3] = gate.rows[3], gate.rows[1]

    return gate

def legacy_FLIPPEDCNOTGATE(nqbits : int = 2) -> Matrix:
    fcnot = Matrix(4, 4) * 0
    fcnot.rows[0][0] = comp(1, 0)
    fcnot.rows[1][3] = comp(1, 0)
    fcnot.rows[2][2] = comp(1, 0)
    fcnot.rows[3][1] = comp(1, 0)
    fcnot.gateid = 'fcnot'
    return fcnot


def legacy_CNOT(qcontrol : list, qtarget : list) ->list:
    tens = qtensor(qcontrol, qtarget)
    result = legacy_CNOTGATE() ** tens
    return result








legacy_reprs = {
    'r' : '[ R ]',
    'i' : '-----',
    'n' : '[ ~ ]',
    'x' : '[ X ]',
    'y' : '[ Y ]',
    'z' : '[ Z ]',
    'h' : '[ H ]',
    'cnot' : '[ ⛒ ]',
    'fcnot' : '[ ⛒ ]',
    'phase' : '[ θ ]',
    'pin' : '[ ☉ ]',
    'c0' : '[ 0 ]',
    'c1' : '[ 1 ]',
    'm' :  '-----',
    's' : '[ s ]'
}







def legacy_getrepr(gate : Matrix) -> str:
    try : return legacy_reprs[gate.gateid]
    except : 
        print(f"Gate represntation not found: {gate.gateid}")
        return '[ ! ]'

def legacycompile(self, verbose : bool = False, showcompilationresult : bool = True):
    if showcompilationresult : print(f"\nCompiling {'program' if self.name is None else self.name}...")
    longest = 0
    for each in range(len(self.gates)):
        for i in range(len(self.gates[each])):
            if self.gates[each][i].gateid == 'cnot':
                pin = legacy_IGATE()
                pin.gateid = 'pin'
                self.gates[each - 1].insert(i, pin)
            elif self.gates[each][i].gateid == 'fcnot':
                pin = legacy_IGATE()
                pin.gateid = 'pin'
                self.gates[each + 1].insert(i, pin)
    for each in range(len(self.gates)):
        if longest < len(self.gates[each]): longest = len(self.gates[each])
    
    for each in range(len(self.gates)):
        self.gates[each] += [legacy_IGATE() for i in range(longest - len(self.gates[each]))]
    
    self.calcrepr()
    
    if showcompilationresult:
        print(self)
        print("Compilation complete!\n")

    if verbose and self.cache is not None: print(self.cache)

    self.cache = None
    return self

def legacycalcrepr(self):
    string = f"\n{'' if self.name is None else self.name}"
    if self.name is not None: string += '\n'
    for each in range(self.nqbits):
        line = f"q{str(each)} ({self.qbits[each][1]}) ⮕ ---"
        for i in range(len(self.gates[each])):
            line += legacy_getrepr(self.gates[each][i])
            line += "---"
        string += line + '\n'
    self.repr = string

def getstate(self, verbose : bool = False, usecache : bool = True) -> list:
    if self.cache != None and usecache:
        print("Using cache...")
        return self.cache
    length = len(self.gates[0])

    state = self.qbits[0]
    for each in range(1, self.nqbits): state = qtensor(state, self.qbits[each])

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