#This is the automated one
import shutil, os, sys
from distutils.dir_util import copy_tree

pseudo = "psps"
#this will be where all the pseudofiles are


#znucl code:

a = os.path.join("/flash", "jtv1", "OCEAN", "bin", pseudo)
aalist = os.listdir(a)
alist = sorted(aalist)
numa = map(int, alist)
#lists everything in pseudo directory    

commonz = open("Common/znucl", "r")
searchz = commonz.readlines()
commonz.close()

#looks at znucl file in Common


print "Available elements (by atomic number):"
print alist

#may want to make the system of atomic numbers to names a dictionary for easier use.

for number in numa:
        if number == 1:
                name = "H"
        elif number == 2:
                name = "He"
        elif number == 3:
                name = "Li"
        elif number == 4:
                name = "Be"
        elif number == 5:
                name = "B"
        elif number == 6:
                name = "C"
        elif number == 7:
                name = "N"
        elif number == 8:
                name = "O"
        elif number == 9:
                name = "F"
        elif number == 10:
                name = "Ne"
        elif number == 11:
                name = "Na"
        elif number == 12:
                name = "Mg"
        elif number == 13:
                name = "Al"
        elif number == 14:
                name = "Si"
        elif number == 15:
                name = "P"
        elif number == 16:
                name = "S"
        elif number == 17:
                name = "Cl"
        elif number == 18:
                name = "Ar"
        elif number == 19:
                name = "K"
        elif number == 20:
                name = "Ca"
        elif number == 21:
                name = "Sc"
	elif number == 22:
                name = "Ti"
        elif number == 23:
                name = "V"
        elif number == 24:
                name = "Cr"
        elif number == 25:
                name = "Mn"
        elif number == 26:
                name = "Fe"
        elif number == 27:
                name = "Co"
        elif number == 28:
                name = "Ni"
        elif number == 29:
                name = "Cu"
        elif number == 30:
                name = "Zn"
        #need to enter more just in case he enters a lot in the future
        else:
                name == "not inputed yet."
        #by atomic number it lists abbreviated name

        print "Atomic number " + str(number) + " is " + name
#prints what atomic number is what symbol for the user





#the final pathmaking code:

def pathmaker1():
	startdict = {}
	lastdict = {}

	#print "pathmaker1"
	#qual = fakequalitylist() #change fakequalitylist() to quality list when ready

	for object in zdict:
		start = os.path.join(a, zdict[object])
		#print start #get rid of this later
		startdict[object] = start
	return startdict

def pathmaker2():
	nextdict = {}
	semicore = semic()	
		
	#print "in pathmaker2"
	for object in pathmaker1():
		next = os.path.join(pathmaker1()[object], semicore[object])
		#print next #get rid of this later
		nextdict[object] = next
	return nextdict	

def pathmaker3():
        lastdict = {}
	qual = qualitylist()
	
	pathlist = []
	for object in pathmaker2():
		last = os.path.join(pathmaker2()[object], qual[object])
		#print last #get rid of this later
		lastdict[object] = last		

		pathlist.append(last)
		#need to grab the letter associated with the object of zdict for the dict semicore
	print "\n"
	print pathlist
#uses the dictionaries made for znucl, semicore, and quality to make paths and then it prints the paths        





#the semicore code:

def semic():
	commons = open("Common/semicore", "r")
	searchs = commons.readlines()
	commons.close()

	y = 0
	sdict = {}
	#print "\nFrom the semicore file:"
	#ssslist = []	
	for index, line in enumerate(searchs):
		semicorestring = index, line
		#print semicorestring[1] 
		letters = semicorestring[1].split()
		#splits string of letters automatically and gets rid of line return and white spaces
		#will need to be editted to allow for one semicore input also. Atm just does them for as many znucls there are		
		for bl in letters:
			#print bl 
			y = y + 1
			sdict[y] = bl
	print "\nSemicore dictionary?"
	print sdict
	return sdict
		
		
		
		#for semicore in sdict:
		#	if semicore in slist:
				
		#		print semicore
	    	#		if semicore.lower() in "t":	
		#			ssslist.append('T')
		#		        
		#		else:
                #			ssslist.append('F')
		#		return ssslist
		#for now i'm not going to check that it's actually one of the options, im just going to pick it			
#looks at semicore file in Common and for each znucl it sees if for it's semicore options which one it should
#pick based on the semicores that correspond to each znucl in Common 





#the quality code:

def qualitylist():
       	print "\nQuality requested in Common file:"

	qdict = {}

	commonq = open("Common/pp.quality", "r")
	searchq = commonq.readlines()
	commonq.close()

	for index, line in enumerate(searchq):
		qualitystring = index, line 
		print qualitystring[1]
		quality_asked_for = int(qualitystring[1])
		#grabs the one quality listed in Common and sets it equal to quality_asked_for (which should only be one value)
		for key in pathmaker2():
			path = pathmaker2()[key]
			qqlist = os.listdir(path)
        		qlist = sorted(qqlist)
        		lengthq = len(qqlist)

			print "\nThis a znucl's options:"
			qnumber_list = map(int, qlist)
			print qnumber_list
			take_closest = lambda num,collection:min(collection,key=lambda x:abs(x-num))
			#should be callable and find the closest value if give a value and a list
			closest = take_closest(quality_asked_for, qnumber_list)
			#finds out closest value to the quality listed in Common so that can be picked
			print "\nThe closest matching option to the quality requested is:"
			print closest
			
			qdict[key] = str(closest)		
			#then for each key in zdict, add a key to qdict with the closest value found above
		#pathmaker2() returns dict of paths that can be used to list of options of quality for each znucl	
					
	print "\n"
	print qdict
	return qdict 
#returns a dictionary to be used in pathmaker3() so it can match each quality with it's znucl and path 




#the znucl code (I may make this a method later but it's not necessary):

zlist = []

for index, line in enumerate(searchz):
	something = index, line
	print "\nList of everything found in znucl file (in Common):"
        print something[1]
	for znucl in alist:
        	if znucl in something[1]:
                	#print znucl #not needed atm
			zlist.append(znucl)
print "The complete list of znucls"
print zlist
#looks through znucl file and grabs the second item in the list it returns (the string with all the data)
#then if any of the stuff in that list is one of the options listed in the psps directory, it grabs that stuff,
#and makes a list of znucls from the stuff

print "\nConverting list of znucls to numbers for sorting:"
print zlist
znumber_list= map(int, zlist)
print znumber_list
znumber_list= sorted(znumber_list)
print znumber_list

zlist = map(lambda x:str(x), znumber_list)
print zlist
#turns znucls in zlist into ints and makes it znumber_list, sorts znumber_list, then changes znucls back to strings for zlist
#This makes sure that the list given to zdict is sorted numerically so keys correspond to order of semicores

print "\nDictionary for znucls:"
x = 0
zdict = {}
for znucl in zlist:
	x = x + 1
	zdict[x] = znucl
#makes the list of znucls into a dictionary so that it's easier to keep track of what's what 

print zdict

	

#print "\nList of semicores:"
#semicorelist = []
#for item in znucllist:
#	s = item
#	sslist = os.listdir(s)
#	slist = sorted(sslist)
#	print "List of options in znucl directories:"
#	print slist
#	lengths = len(sslist)
#	#makes path and lists all of semicore options, sorts it, and sees how long it is
#
#	if lengths == 1:		
#        	sssc = 'F'
#		semicorelist.append(sssc)
#		pathmaker()
#        	#if there's only one option for a znucl, that option is picked
#	elif lengths == 2:
#		print "List of semicores from Common (for ti in this case):"
#	     	for sssc in semic():
#			print sssc
#			semicorelist.append(sssc)
#			pathmaker()
#  		#if there's two options for a znucl, code looks to see if it's T or F 
#	else:
#       		print "Something went wrong."
#       		sys.exit(0)

pathmaker3()	
	
						
