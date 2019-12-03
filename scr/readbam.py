#!usr/bin/python
import os
import re
import sys
import time
import scipy.stats
import pysam

def hashRef(ref_path,ref_name):
##build reference dictionary
	REF = open(ref_path,'rb').readlines()
	Ref, Rev, For, One, Two = {}, {}, {}, {}, {}
	for lines in REF:
		lines = lines.rstrip()
		if re.match('^>' + ref_name,lines):
			continue
		else:
			i = 1
			bases = list(lines)
			for b in bases:
				Ref[i] = b
				Rev[i] = {'A':0, 'T':0, 'G':0, 'C':0, 'N':0}
				For[i] = {'A':0, 'T':0, 'G':0, 'C':0, 'N':0}
				One[i] = {'A':0, 'T':0, 'G':0, 'C':0, 'N':0}
				Two[i] = {'A':0, 'T':0, 'G':0, 'C':0, 'N':0}
#				print i
				i = i + 1
	return Ref, Rev, For, One, Two
Ref, Rev, For, One, Two = hashRef(sys.argv[2],'MT')


totalReadNumber = 0	#total seq read number
ttMapReadNumber = 0	#total map read number
qcMapReadNumber = 0	#up qc map read number
unMapReadNumber = 0	#isnot map read number
rvMapReadNumber = 0	#reverse map read number
fwMapReadNumber = 0	#forward map read number
r1MapReadNumber = 0	#read1 map read number
r2MapReadNumber = 0	#read2 map read number

BAM = pysam.AlignmentFile(sys.argv[1],'rb')
for read in BAM.fetch('MT',1,16596):
#for read in BAM.fetch('DL',1,1000):
	totalReadNumber = totalReadNumber + 1
	seq = read.query_sequence
	key = read.query_name
	pos = read.pos
	if not read.is_unmapped and len(read.cigar) == 1:
		ttMapReadNumber = ttMapReadNumber + 1
#		print ttMapReadNumber
		if not read.is_qcfail and len(seq) == 50:
			qcMapReadNumber = qcMapReadNumber + 1
			if read.is_reverse:
				rvMapReadNumber = rvMapReadNumber + 1
				for base in seq:
					if pos in Rev.keys():
						Rev[pos][base] = Rev[pos][base] + 1
						pos = pos + 1
					else:
						Rev[pos] = {'A':0, 'T':0, 'G':0, 'C':0, 'N':0}
						Rev[pos][base] = Rev[pos][base] + 1
						pos = pos + 1
#						print 'reverse', Rev[pos], pos
				if read.is_read1:
					r1MapReadNumber = r1MapReadNumber + 1
					for base in seq:
						if pos in One.keys():
							One[pos][base] = One[pos][base] + 1
							pos = pos + 1
						else:
							One[pos] = {'A':0, 'T':0, 'G':0, 'C':0, 'N':0}
							One[pos][base] = One[pos][base] + 1
							pos = pos + 1	
#						print 'readone', One[pos], One[pos][base], pos
				if read.is_read2:
#						print read.pos, read.query_sequence
					r2MapReadNumber = r2MapReadNumber + 1
					for base in seq:
						if pos in Two.keys():
							Two[pos][base] = Two[pos][base] + 1
							pos = pos + 1
						else:
							Two[pos] = {'A':0, 'T':0, 'G':0, 'C':0, 'N':0}
							Two[pos][base] = Two[pos][base] + 1
							pos = pos + 1
#						print 'readtwo', Two[pos], Two[pos][base], pos
			if not read.is_reverse:
				fwMapReadNumber = fwMapReadNumber + 1
				for base in seq:
					if pos in For.keys():
						For[pos][base] = For[pos][base] + 1
						pos = pos + 1
					else:
						For[pos] = {'A':0, 'T':0, 'G':0, 'C':0, 'N':0}
						For[pos][base] = For[pos][base] + 1
						pos = pos + 1
#						print 'forward', For[pos], For[pos][base], pos
				if read.is_read1:
					r1MapReadNumber = r1MapReadNumber + 1
					for base in seq:
						if pos in One.keys():
							One[pos][base] = One[pos][base] + 1
							pos = pos + 1
						else:
							One[pos] = {'A':0, 'T':0, 'G':0, 'C':0, 'N':0}
							One[pos][base] = One[pos][base] + 1
							pos = pos + 1	
#						print 'readone', One[pos], One[pos][base], pos
				if read.is_read2:
#					print read.pos, read.query_sequence
					r2MapReadNumber = r2MapReadNumber + 1
					for base in seq:
						if pos in Two.keys():
							Two[pos][base] = Two[pos][base] + 1
							pos = pos + 1
						else:
							Two[pos] = {'A':0, 'T':0, 'G':0, 'C':0, 'N':0}
							Two[pos][base] = Two[pos][base] + 1
							pos = pos + 1
#						print 'readtwo', Two[pos], Two[pos][base], pos
	else:
		unMapReadNumber = unMapReadNumber + 1
def strandBiasScore(a,b,c,d):
	depth = a + b + c + d
	first = a + b + .0
	later = c + d + .0
	major = a + c + .0
	minor = b + d + .0
	if major > 0 and minor > 0:
		SBS = abs(b/first - d/later)/(minor/depth)	#strand bias score
		GSB = max( (b/first) * (c/later) / (major/depth), (a/first) * (d/later) / (major/depth))	#gatk strand bias score		
		return round(SBS,3),round(GSB,3)
	elif major == 0 and minor > 0:
		return 100, 100
	elif major >0 and minor == 0:
		return 100, 100
	elif depth == 0:
		return 0, 0
	else:
		return None, None
def topTwoMmax(list):
	firstmax = max(list)
	list.remove(firstmax)
	secondmax = max(list)
	return firstmax, secondmax
#ATGCN = ['A', 'T', 'G', 'C', 'N']
for pos, base in Ref.items():
	if pos in For.keys() and pos in Rev.keys():
		Fvalues, Rvalues = For[pos].values(), Rev[pos].values()
		Fkeys,Rkeys = For[pos].keys(),Rev[pos].keys()
		Fa, Fb = topTwoMmax(Fvalues)
		Rc, Rd = topTwoMmax(Rvalues)
		print pos, strandBiasScore(Fa, Fb, Rc, Rd)
	if pos in One.keys() and pos in Two.keys():
		Ovalues, Tvalues = One[pos].values(),Two[pos].values()
		Okeys, Tkeys = One.keys(),Two.keys()
		Oa, Ob = topTwoMmax(Ovalues)
		Tc, Td = topTwoMmax(Tvalues)
		print pos, strandBiasScore(Oa, Ob, Tc, Td)

print 'totalReadNumber: %s' %totalReadNumber
print 'ttMapReadNumber: %s' %ttMapReadNumber
print 'qcMapReadNumber: %s' %qcMapReadNumber
print 'unMapReadNumber: %s' %unMapReadNumber
print 'rvMapReadNumber: %s' %rvMapReadNumber
print 'fwMapReadNumber: %s' %fwMapReadNumber
print 'r1MapReadNumber: %s' %r1MapReadNumber
print 'r2MapReadNumber: %s' %r2MapReadNumber

