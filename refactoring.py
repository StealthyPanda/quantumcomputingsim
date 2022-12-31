from quantum import *



bell = Block(2, 'B')
bell.addgates(0, [hgate, cnot0])
bell = bell.getmat()


q = qprogram(2)

q.addgates(0, [bell, bell, bell, bell, bell, bell, bell, bell])

q.compile()

q.run(graph = True)