from quantum import *


state = tensor(HAD(qbit(0)), qbit(0))
state = tensor(state, HAD(qbit(1)))

print(state, "\n\n")
print(HAD(state))
