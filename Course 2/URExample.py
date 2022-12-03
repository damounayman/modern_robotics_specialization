import modern_robotics as mr
import numpy as np
import math
from math import cos,sin
from matplotlib import pyplot as plt
from numpy import matmul
from modern_robotics import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib 
matplotlib.use('tkagg')
'''
UR5e
Kinematics	theta [rad]	a [m]		d [m]		alpha [rad]	Dynamics	Mass [kg]	Center of Mass [m]
Joint 1		0		0		0.1625		pi/2			Link 1		3.761		[0, -0.02561, 0.00193]
Joint 2		0		-0.425		0		0			Link 2		8.058		[0.2125, 0, 0.11336]
Joint 3		0		-0.3922		0		0			Link 3		2.846		[0.15, 0.0, 0.0265]
Joint 4		0		0		0.1333		pi/2			Link 4		1.37		[0, -0.0018, 0.01634]
Joint 5		0		0		0.0997		-pi/2			Link 5		1.3		[0, 0.0018,0.01634]
Joint 6		0		0		0.0996		0			Link 6		0.365		[0, 0, -0.001159]
'''
L1 = 0.425
L2 = 0.392
W1 = 0.109
W2 = 0.082
H1 = 0.089
H2 = 0.095
PI = math.pi
th1 = PI/2
th2 = PI/2
th3 = PI/2

M  = np.array(
[[-1,0,0,L1+L2],
[0,0,1,W1+W2],
[0,1,0,H1-H2],
[0,0,0,1]])
p0 =np.array([0,0,0,1])

for i in range(180):

	S1 = [0,0,1,0,0,0]
	S2 = [0,1,0,-H1,0,0]
	S3 = [0,1,0,-H1,0,L1]
	S4 = [0,1,0,-H1,0,L1+L2]
	S5 = [0,0,-1,-W1,L1+L2,0]
	S6 = [0,1,0,H2-H1,0,L1+L2]
	th1 = i*PI/180
	th2 = i*PI/180
	th3 = i*PI/180
	th4 = i*PI/180/2
	th5 = i*PI/180/2
	th6 = i*PI/180/2

	Slist = np.array([S1,S2,S3,S4,S5,S6]).T
	thetalist = np.array([th1, th2, th3,th4,th5,th6])
	M  = np.array(
	[[-1,0,0,L1+L2],
	[0,0,1,W1+W2],
	[0,1,0,H1-H2],
	[0,0,0,1]])
	T07 = FKinSpace(M,Slist,thetalist)

	Slist = np.array([S1,S2,S3,S4,S5]).T
	thetalist = np.array([th1,th2, th3,th4,th5])
	M  = np.array(
	[[-1,0,0,L1+L2],
	[0,0,1,W1],
	[0,1,0,H1-H2],
	[0,0,0,1]])
	T06 = FKinSpace(M,Slist,thetalist)

	Slist = np.array([S1,S2,S3,S4]).T
	thetalist = np.array([th1,th2, th3,th4])
	M  = np.array(
	[[-1,0,0,L1+L2],
	[0,0,1,W1],
	[0,1,0,H1],
	[0,0,0,1]])
	T05 = FKinSpace(M,Slist,thetalist)

	Slist = np.array([S1,S2,S3]).T
	thetalist = np.array([th1,th2, th3])
	M  = np.array(
	[[0,0,1,L1],
	[-1,0,0,W1],
	[0,1,0,H1],
	[0,0,0,1]])
	T04 = FKinSpace(M,Slist,thetalist)

	Slist = np.array([S1,S2]).T
	thetalist = np.array([th1,th2])
	M  = np.array(
	[[1,0,0,0],
	[0,1,0,W1],
	[0,0,1,H1],
	[0,0,0,1]])
	T03 = FKinSpace(M,Slist,thetalist)

	Slist = np.array([S1]).T
	thetalist = np.array([th1])
	M  = np.array(
	[[1,0,0,0],
	[0,1,0,0],
	[0,0,1,H1],
	[0,0,0,1]])
	T02 = FKinSpace(M,Slist,thetalist)


	Slist = np.array([]).T
	thetalist = np.array([])
	M  = np.array(
	[[1,0,0,0],
	[0,1,0,0],
	[0,0,1,0],
	[0,0,0,1]])
	T01 = FKinSpace(M,Slist,thetalist)

	p1 = matmul(T01,p0)
	p2 = matmul(T02,p0)
	p3 = matmul(T03,p0)
	p4 = matmul(T04,p0)
	p5 = matmul(T05,p0)
	p6 = matmul(T06,p0)
	p7 = matmul(T07,p0)
	print(p1)
	print(p2)
	print(p3)

	l1 = np.linalg.norm(p1[:-1]-p2[:-1])
	l2 = np.linalg.norm(p2[:-1]-p3[:-1])
	l3 = np.linalg.norm(p3[:-1]-p4[:-1])
	l4 = np.linalg.norm(p4[:-1]-p5[:-1])
	l5 = np.linalg.norm(p5[:-1]-p6[:-1])
	l6 = np.linalg.norm(p6[:-1]-p7[:-1])
	print(T02)
	print(l1,'\t',l2,'\t',l3,'\t',l4,'\t',l5,'\t',l6)
	fig = plt.gcf()
	ax = plt.axes(projection='3d')
	fig.set_size_inches(5, 5)
	plt.rcParams["figure.figsize"] = (14,14)
	plt.rcParams['lines.linewidth'] = 4
	plt.rcParams['lines.color'] = 'r'
	plt.rcParams['axes.grid'] = True 
	ax.set_xlim(-0.5,0.5)
	ax.set_ylim(-0.5,0.5)
	ax.set_zlim(-0.5,0.5)
	plt.title('Screw')

	plt.plot([p1[0],p2[0],p3[0],p4[0],p5[0],p6[0],p7[0]], [p1[1],p2[1],p3[1],p4[1],p5[1],p6[1],p7[1]],[p1[2],p2[2],p3[2],p4[2],p5[2],p6[2],p7[2]],'.-')
	plt.plot([p1[0],p2[0],p3[0],p4[0],p5[0],p6[0],p7[0]], [p1[1],p2[1],p3[1],p4[1],p5[1],p6[1],p7[1]],[p1[2],p2[2],p3[2],p4[2],p5[2],p6[2],p7[2]],'r.')

	plt.draw()
	plt.pause(0.05)
	plt.clf()