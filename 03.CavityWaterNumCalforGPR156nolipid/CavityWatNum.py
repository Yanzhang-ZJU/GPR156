"""
MIT License
Copyright (c) 2023 Kun, Xi @ Zhejiang University
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

#part 1
import numpy as np
import mdtraj as md
#from Confs import *
import os

def digits(s1):
	s2 = "%.5d" % s1
	return s2

# Function to count the nunber of water molecules in the cavity #
def cavWat(ftrj, Nframe, Natom1, Natom2, Natom3, Natom4, Natom5):
	disminref=-1.0

	if (float(ftrj[Nframe][Natom5][2]) > 4.9) and \
		(float(ftrj[Nframe][Natom5][2]) < 6.5) :
		disminref = 0.0

		if (float(ftrj[Nframe][Natom5][0]) >= 4.7) and \
			( float(ftrj[Nframe][Natom5][0]) <= 8.7):
			if (float(ftrj[Nframe][Natom5][1]) >= 6.0) and \
				( float(ftrj[Nframe][Natom5][1]) <= 7.5):
				disminref = 1.0
	return disminref

#Load traj;

dirRoot = os.getcwd()

trjName = 'equiallfit.xtc'	# Merged trajectory file in xtc format with a order of energy minimization / \
							# NVT ensemble simulation / NPT ensemble simulation # 
tmpName = 'protein.pdb'		# Template file #


f = md.load(dirRoot + '/' + trjName, top=dirRoot + '/' + tmpName)
fa = f.xyz

#Atom index (@CA): F84 / I97 / F382 / I397 	# Atoms to define the cavity size #
atomInx = ['1333', '1520', '6065', '6285']

#Water Index # All water molecules #
wat_ini = 74359
wat_fin = 223193

fw1 = open('Cavwat.dat', 'w')

time = 0.0
for i in range(len(fa)):
	watnum = 0.0
	for j in range(wat_ini-1, wat_fin, 2):
		dismin = cavWat(fa, i, int(atomInx[0])-1, int(atomInx[1])-1, int(atomInx[2])-1, int(atomInx[3])-1, j)
		if float(dismin) > 0.0:
			watnum+=1.0

	#Energy minimization#
	if i < 5:
		time += 50.0
		line = 'EM1k ' + str(time) + ' Numwat ' + str(watnum)
		print(line, file=fw1)

	#NVT ensemble simulation#
	if i >= 5 and i < 55:
		time += 10.0
		line = 'NVT1 ' + str(time) + ' Numwat ' + str(watnum)
		print(line, file=fw1)

	#NPT ensemble simulation with timestep of 0.001 fs#
	if i>=55 and i < 80:
		time += 5.0
		line = 'NPT1 ' + str(time) + ' Numwat ' + str(watnum)
		print(line, file=fw1)		
	#NPT ensemble simulation with timestep of 0.002 fs#
	if i>=80 and i < len(fa):
		time += 10.0
		line = 'NPT2 ' + str(time) + ' Numwat ' + str(watnum)
		print(line, file=fw1)

fw1.close()


