from quantum import *

# print(CNOTR(controlbitindex=1))
# print(CNOTR(controlbitindex=1).span)

q = qprogram(2, name = 'refactor testing')

q.addgates(0, [HGATE(), HGATE()])
q.addgates(1, [IGATE(), CNOTR()])

q.compile()

q.run(graph = True, terminal = True)