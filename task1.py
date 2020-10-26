from sympy import *
import numpy as np
import matplotlib.pyplot as plt
l = 0.1 # in task length of link 2 not given so i decide to make it equal link 1
q0, q1= symbols('q0 q1')
m1, m2, I1, I2, g, L2, d2, h = 2,2,1,2,9.81, 1, 0.5, 10
L1= symbols('L1')


def cross(a, b): return Matrix([a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]])

O0 = Matrix([[0],[0],[0]])
O1 = Matrix([[l*cos(q0)], [l*sin(q0)], [0]])
O2 = Matrix([[(q1 + l)*cos(q0)], [(q1 + l)*sin(q0)], [0]])

R1 = Matrix([[cos(q0), -sin(q0), 0], [sin(q0), cos(q0), 0], [0, 0,1]])
R2 = R1
Jw1 = Matrix([[0,0,1],[0,0,0]]).transpose()
Jw2 = Jw1
Jv1 = Matrix([[0,0,0],[0,0,0]]).transpose()
Jv2 = Matrix([[-(q1+d2)*sin(q0),(q1+d2)*cos(q0),0],[cos(q0),sin(q0),0]]).transpose()
#Jv1 = Matrix([[cross(Matrix([[0],[0],[1]]), (O1-O0))],[zeros(3, 1)]]).reshape(3,2)
#Jv2 = Matrix([[cross(Matrix([[0],[0],[1]]), (O2-O0))],[cross(Matrix([[0],[0],[1]]), (O2-O1))]]).reshape(3,2)


D1 = m1*Jv1.transpose()*Jv1+Jw1.transpose()*R1*I1*R1.transpose()*Jw1
D2 = m2*Jv2.transpose()*Jv2+Jw2.transpose()*R2*I2*R2.transpose()*Jw2

D = sympify(D1+D2)
print(D[0,0], ", ", D[0,1])
print(D[1,0], ", ", D[1,1])
P1 = m1*g*h
P2 = m2*g*(h+(q1+d2)*sin(q0))

P = sympify(P1+P2)

G = Matrix([[diff(P, q0)],[diff(P, q1)]])

dq1, dq2, ddq1, ddq2 = symbols('dq1 dq2 ddq1 ddq2')

q, dq, ddq = [q0, q1], Matrix([dq1, dq2]),Matrix([ddq1, ddq2])

C = zeros(2)
for k in range(2):
    for j in range(2):
        for i in range(2): C[k,j] += 0.5*(diff(D[k,j],q[i]) + diff(D[k,i],q[j])-diff(D[i,j],q[k]))*dq[i]
sympify(C)
tor = D*ddq+C*dq+G

q00, q10, dq00, ddq00, dq10, ddq10 = pi / 2, 0, 0,0,0,0
dt = 0.01
U = Matrix([0, 5.81*2]) # stable vertical position with fi = pi/2
n = 300
t = [i/10 for i in range(300)]
q0p, q1p, dq0p, dq1p, ddq0p, ddq1p = [],[],[],[],[],[]
for i in range(300):
    q0p.append(q00)
    q1p.append(q10)
    dq0p.append(dq00)
    dq1p.append(dq10)
    ddq0p.append(ddq00)
    ddq1p.append(ddq10)
    #dd = D.subs(list(zip([q0, q1], [q00, q10]))).inv()
    #cc = C.subs(list(zip([q0, q1, dq00, dq10], [q00, q10, dq00, dq10])))*dq
    #gg = G.subs(list(zip([q0, q1], [q00, q10])))
    #print(dd, cc, gg)
    ddqt = D * (U - C*dq - G)
    ddq = ddqt.subs(list(zip([q0, q1, dq1, dq2], [q00, q10, dq00, dq10])))
    ddq00 = ddq[0]
    ddq11 = ddq[1]
    dq00 += float(ddq[0] * dt)
    dq10 += float(ddq[1] * dt)
    q00 += dq00 * dt
    q10 += dq10 * dt
colors = ['r', 'g', 'b']
coords = ['x', 'y', 'z']

#plt.figure(1)
#plt.plot(t, q0p, color = 'r', linestyle = 'solid', label = 'fi')
#plt.plot(t, q1p, color = 'g', linestyle = 'solid', label = 's')
#plt.figure(2)
#plt.plot(t, dq0p, color = 'r', linestyle = 'solid', label = 'fi')
#plt.plot(t, dq1p, color = 'g', linestyle = 'solid', label = 's')
#plt.figure(3)
plt.plot(t, ddq0p, color = 'r', linestyle = 'solid', label = 'fi')
plt.plot(t, ddq1p, color = 'g', linestyle = 'solid', label = 's')
plt.legend(loc='upper right')
plt.show()
