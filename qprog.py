from quantum import *
import threading

newbasis = basischange(HGATE())

q = qprogram(2, bchange = None, name = 'No Basis')
r = qprogram(2, bchange = newbasis, name = 'H Basis')

q.addgates(0, [HGATE()])
q.addgates(1, [IGATE(), CNOTGATE()])
r.addgates(0, [HGATE()])
r.addgates(1, [IGATE(), CNOTGATE()])
# q.addgates(1, [IGATE(), CNOTGATE()])
# q.addgates(2, [IGATE(), IGATE(), CNOTGATE()])

q.compile()
r.compile()

q.run(graph=True)
r.run(graph=True)