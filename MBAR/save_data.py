#!/usr/bin/env python
import numpy as np

def distances_from_dump(infile, style='beadspring', POLY_END_TYPE = 1, POLY_MID_TYPES = [2], COLLOID_TYPE = 4, TYPE_IGNORES = [3], id1 = None, id2 = None,startStep=0):
	if type(POLY_MID_TYPES) == int:
		POLY_MID_TYPES = [POLY_MID_TYPES]
	if POLY_END_TYPE in POLY_MID_TYPES:
		print 'ERROR: You specified that the end type of the polymer was in POLY_MID_TYPES.'
		raise ValueError
	if type(TYPE_IGNORES) == int:
		TYPE_IGNORES = [TYPE_IGNORES]

	# Create a class of atoms so that I can 
	class Atom:
		def __init__(self, id_, type_):
			# the _ is there because id and type are built in to python and can't be overridden
			self.Pos = []
			self.Box = []
			self.neighList = []
			self.atomType = int(type_)
			self.atomID = int(id_)

		def addCoord(self, pos, box):
			'''Add a coordinate to this atom
			The expected box is of the form [(xlo,xhi),(ylo,yhi),(zlo,zhi)]'''
			for i in range(3):
				pos[i] = float(pos[i])
				for j in range(2):
					box[i][j] = float(box[i][j])	
			self.Pos.append(pos)
			self.Box.append(box)

		def addNeighbor(self, neighbor):
			'''Specify an atom that is bonded to this one'''
			self.neighList.append(neighbor)


	#First must generate a list of atoms with their respective IDs and atom types
	##print 'Creating atom list'
	with open(infile, 'rb') as f:
		dumplist = f.read().split('ITEM: TIMESTEP')

	atoms = []
	del dumplist[0] # this is an empty string


	for i in dumplist[0].split('ITEM: ATOMS id type xs ys zs')[1].split('\n'): 
	# This is a terrible way of looping through the lines that have the initial position info I need
		#print repr(i)
		line = i.split()
		if len(line) < 1:
			continue
		id_ = line[0]
		type_ = line[1]
		atoms.append(Atom(id_, type_))
	atoms.sort(key=lambda atom: atom.atomID) #Sort atoms by ID


	#Fill atoms with position data
	##print 'Filling position values'
	skipCount = 0
	for indx,timestepData in enumerate(dumplist):
		timestep = dumplist[indx].split('\n')[1]
		if int(timestep)<int(startStep):
			skipCount+=1
			continue
		temp = timestepData.split('ITEM: ATOMS id type xs ys zs')[1].split('\n')
		temp2 = timestepData.split('ITEM: ATOMS id type xs ys zs')[0].split('ITEM: BOX BOUNDS pp pp pp')[1].split('\n')
		box = []
		for i in temp2:
			if i != '':
				box.append(i.split())
				for i in range(len(box)):
					for j in range(len(box[i])):
						box[i][j] = float(box[i][j])

		if [] in box:
			del box[0]

		for atom_data in temp:
			# print repr(atom_data)
			atom_data = atom_data.split()
			if len(atom_data) == 0:
				continue
			id_ = int(atom_data[0]) # id
			Pos = [float(atom_data[2]), float(atom_data[3]), float(atom_data[4])]

			for atom in atoms:
				if atom.atomID == id_:
					for i in range(3):
						Pos[i] = (box[i][1] - box[i][0])*Pos[i] + box[i][0]
					atom.addCoord(Pos, box)
					break
			else:
				print "ID not found ", id_
				print timestepData.split('ITEM: ATOMS id type xs ys zs')[0]
				raise ValueError

	##print 'Number of data points per atom (timesteps recorded):', len(atoms[0].Pos)
	##print 'Number of points skipped:', skipCount
	#list of atoms has been filled with ID, type and position data
	#Now to add neighbors
	##print 'Adding neighbors'
	if style == 'beadspring':
		polyEnd = False
		for i in range(len(atoms)):
			if atoms[i].atomType == POLY_END_TYPE:
				if not polyEnd:
				# This atom is on the end of a polymer, and it is the beginning of a polymer
					atoms[i].addNeighbor(atoms[i+1])
					polyEnd = True
				else:
					atoms[i].addNeighbor(atoms[i-1])
					polyEnd = False
			elif atoms[i].atomType in POLY_MID_TYPES:
				atoms[i].addNeighbor(atoms[i+1])
				# atoms[i].addNeighbor(atoms[i-1])
			elif atoms[i].atomType not in TYPE_IGNORES:
				print "WARNING: Atom of unknown type encountered."
				print "atom ID ", atoms[i].atomID

	elif style == 'colloid':
		#distance between atoms of one type
		colloids = []
		for i in range(len(atoms)):
			if atoms[i].atomType != COLLOID_TYPE and atoms[i].atomType not in TYPE_IGNORES:
				print "WARNING: Atom of unknown type encountered."
				print "Atom ID ", atoms[i].atomID
			else:
				colloids.append(atoms[i])
		for i in range(len(colloids)):
			for j in range(i+1, len(colloids)):
				atoms[i].addNeighbor(atoms[j])
		atoms = colloids[:]

	elif style == 'id':
		#distance between the atoms with specified ids
		temp = []
		for atom in atoms:
			if atom.atomID == id1:
				temp.append( atom )
				id1 = None
			elif atom.atomID == id2:
				temp.append( atom )
				id2 == None
		atoms = temp

		for i in range(len(atoms)):
			for j in range(i+1, len(atoms)):
				atoms[i].addNeighbor(atoms[j])


	##print 'Calculating distance data of neighbors'
	#generate distance data of neighbors
	dists = []

	def nint(x):
		for i,j in enumerate(x):
			x[i] = round(j)
		return x

	for i in atoms:
		for j in i.neighList:
			for k, x1 in enumerate(i.Pos):
				#Find minimum distance - because periodic boundaries are a thing the bond might be across two walls
				# dist = math.sqrt(sum(min( (i.Pos[k][l]-j.Pos[k][l])**2. , (i.Pos[k][l]+j.Box[k][l][0]-j.Pos[k][l]-j.Box[k][l][1])**2., (i.Pos[k][l]+j.Box[k][l][1]-j.Pos[k][l]-j.Box[k][l][0])**2. ) for l in range(3)))
				# Pretty sure the above formula is right, it's the distance between atom1 and the wall + the distance between atom2 and the opposite wall
				x2 = j.Pos[k]
				L = np.array([0.,0.,0.])
				for box_indx,box_boundary in enumerate(i.Box[k]):
					L[box_indx] = box_boundary[1]-box_boundary[0]
				#print L
				r_ij = np.array(x2)-np.array(x1)
				r0_ij = r_ij-L*(nint(r_ij/L))
				dist = np.linalg.norm(r0_ij)
				# if dist <= min(L)/2.:
				dists.append(dist)
	return dists

def save_data(infiles, r0, k ,id1, id2, outfile, beta=1, correl_time = 50):
	allDist2D = {}
	for i,f in enumerate(infiles):
	  dists = distances_from_dump(f, style='id', id1=id1, id2=id2)
	  #remove tails from distance data
	  freq = 20
	  c,b = np.histogram(dists,100)
	  startLen= len(dists)
	  for j,c_v in enumerate(c):
	    if c_v < freq:
	      dists = [ x for x in dists if not (b[j] <= x <= b[j+1]) ]
	  print(startLen-len(dists),'points removed')
	  if (k[i],r0[i]) in allDist2D.keys():
	    allDist2D[(k[i],r0[i])] += dists
	  else:
	    allDist2D[(k[i],r0[i])] = dists



	print '\nUmbrellas are given by (k,r0) pairs:'
	print allDist2D.keys()
	print( list( (key,len(item)) for key,item in allDist2D.iteritems() ) ) 

	#allDist2D was defined as a dictionary because it was convenient, but now we want sorted data
	#remake r0,k lists (to remove duplicate umbrellas)
	r0=[]
	k=[]
	for key in sorted(allDist2D.keys()):
	  r0.append(key[1])
	  k.append(key[0])


	temp = []
	for i,_ in enumerate(r0):
	  temp.append( list(allDist2D[(k[i],r0[i])])[0::int(correl_time)] )
	allDist2D = np.array(temp)


	allDists1D = []
	for i in temp:
	  allDists1D += i

	allDists1D = np.array(allDists1D)

	#And number of independent samples for each state
	numSamples = np.zeros(len(r0), dtype=int)
	for i,temp in enumerate(allDist2D):
	  numSamples[i] = len(temp)

	print('len(numSamples) = ',len(numSamples))

	Umat = np.zeros( (len(r0), sum(numSamples)) )
	
	print('sum numSamples: ', sum(numSamples))
	print('allDist2D: ', np.shape(allDist2D))
	print('r0: ', np.shape(r0))
	print('Umat: ', np.shape(Umat))

	for i,_ in enumerate(r0):
	  Umat[i,:] =  beta*k[i]*(allDists1D - r0[i])**2


	np.savez(outfile, allDists1D=allDists1D, Umat=Umat, numSamples=numSamples)
