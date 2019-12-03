#!usr/bin/python
import os
import re
import sys
import time
import math
import getopt
#import matplotlib as mpl
#mpl.use('Agg')
#import matplotlib.pyplot as plt
#from matplotlib.font_manager import FontProperties

###############################################################################################################################
class HETEROPLASMY:
	'Detecting Heteroplasmy .pileup'
	wholeDR = []
	reverse = ['A','T','C','G','a','t','c','g']
	Reverse = {'A':'a','T':'t','C':'c','G':'g','N':'n'}
	def __init__(self, pos, refBase, base, qual, qmin, qval, htpr):
		self.pos = pos
		self.ref = refBase
		self.base = base
		self.qual = qual
		self.qmin = qmin	#quality format lower
		self.qval = qval	#quality filter value
		self.htpr = htpr
#	def __enter__(self):
#		return self
#	def __exit__(self, type, value, trace):
#		print "type:", type
#		print "value", value
#		print "trace", trace
	def removeIndel(self):		#remove indel characters
		if re.match(r'.+(-|\+).+',self.base):
			mode = re.split(r'(\-|\+)',self.base)
			for i in range(len(mode)):
				inum = re.compile(r'\d+')
				if re.match(inum, mode[i]):
					indel = re.match(inum, mode[i]).group()
					INDEL = int(indel) + len(indel)
					mode[i] = mode[i][INDEL:]
			self.base = ''.join(mode)
	def removeLowQualityBase(self):		#remove low quality base
		module = re.compile(r'(\^.)|\$|\+|-')
		self.base = re.sub(module, '', self.base)
		base = list(self.base)
		qual = list(self.qual)
		time = len(base)
		n = 0
		if (len(base)) == len(qual):
			for q in qual:
				qabs = ord(q) - qmin
				if qabs < qval:
					base[n] = ''
					n += 1
				else:
					n += 1
			self.base = ''.join(base)
		else:
			pass
	def seqbias(self):
		def strandBiasScore(a,b,c,d):
			depth = a + b + c + d + 0.
			first = a + b + .0
			later = c + d + .0
			major = a + c + .0
			minor = b + d + .0
			if major > 0 and minor > 0 and first > 0 and later > 0:
				SBS = abs(b/first - d/later)/(minor/depth)	#strand bias score
				GSB = max( (b/first) * (c/later) / (major/depth), (a/first) * (d/later) / (major/depth))	#gatk strand bias score
				return round(SBS,3),round(GSB,3)
			elif major == 0 or minor == 0:
				return 100, 0
			elif later == 0 or first == 0:
				return 0, 100
			elif depth == 0:
				return 0, 0
			else:
				return None, None
		def topTwoMmax(list):
			firstmax = max(list)
			list.remove(firstmax)
			secondmax = max(list)
			return firstmax, secondmax
		Reverse = {'a':0, 't':0, 'g':0, 'c':0, 'n':0}
		Forward = {'A':0, 'T':0, 'G':0, 'C':0, 'N':0}
		for s in self.base:
			if s in Reverse.keys():
				Reverse[s] += 1
			elif s in Forward.keys():
				Forward[s] += 1
			else:
				continue
		Fvalues, Rvalues = Forward.values(), Reverse.values()
		Fkeys,Rkeys = Forward.keys(),Reverse.keys()
		Fa, Fb = topTwoMmax(Fvalues)
		Rc, Rd = topTwoMmax(Rvalues)
		SEQBScore, GATKScore = strandBiasScore(Fa, Fb, Rc, Rd)
		return Fa, Fb, Rc, Rd, SEQBScore, GATKScore
#		maps = histogram(self.base)
#		return maps
	def statistics(self):
		def histogram(seq):
			base = {}
			for s in seq:
				if s in base.keys():
					base[s] += 1
				else:
					base[s]  = 1
			return base
		def coefficent(hrate,depth):
			if depth > 0:
				factor = 1.96 * math.sqrt((hrate * (1 - hrate)) / depth)
				lower = 0
				upper = 1
				if (hrate - factor) > 0:
					lower = hrate - factor
				if (hrate - factor) < 1:
					upper = hrate + factor
				return ('%0.4f' %lower),('%0.4f' %upper)
			else:
				return 0,0
		base = []
		dict = histogram(self.base)
		for b in self.reverse:
			if b in dict:
				base.append(int(dict[b]))
			else:
				base.append(0)
		forwardBase = base[0:4]
		reverseBase = base[4:8]
		Fusion = map(lambda x: x[0]+x[1],zip(forwardBase,reverseBase))
		fusion = sorted(Fusion,reverse = True)
		major = fusion[0]
		minor = fusion[1]
		majorBase = Fusion.index(major)
		minorBase = Fusion.index(minor)
		shangZhi = sum([-math.log(i)*i for i in filter(lambda x:x != 0, [float(i)/sum(fusion) for i in fusion])])
		
		if major > 0:
			hrate = (minor + 0.0)/(major + minor)
			depth = sum(Fusion)
			lower, upper = coefficent(hrate,depth)
			self.heteroplasmyStatsitic = [self.reverse[majorBase],self.reverse[minorBase],str(depth),str(major),str(minor),str(hrate),str(lower),str(upper),str(shangZhi)]
		else:
			self.heteroplasmyStatsitic = ['0',self.reverse[majorBase], self.reverse[minorBase]]
		return self.heteroplasmyStatsitic

name, input, output, htpr = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])

pile = open(output + name + '.pileup','rb').readlines()

qmin, qval, = 33, 30,
hetnum, bisnum, entropy = 0, 0, 0
for p in pile:
	p = p.rstrip()
	chr, pos, ref, dep, base, qual, = p.split()
	het = HETEROPLASMY(pos, ref, base, qual, qmin, qval, htpr)
	het.removeIndel()
	het.removeLowQualityBase()
	shv = het.statistics()
	entropy += float(shv[-1])
	if float(shv[5]) > htpr/100. and float(shv[6]) > htpr/100.:
		if het.seqbias():
			Fa, Fb, Rc, Rd, SEQBScore, GATKScore = het.seqbias()
			if SEQBScore <= 1:
				BIAS = 0
				hetnum = hetnum + 1
				print(','.join([str(i) for i in [name, pos, ','.join(shv), BIAS, Fa, Fb, Rc, Rd, SEQBScore, GATKScore]]))
			else:
				BIAS = 1
				bisnum = bisnum + 1
				print(','.join([str(i) for i in [name, pos, ','.join(shv), BIAS, Fa, Fb, Rc, Rd, SEQBScore, GATKScore]]))
print "#", name, hetnum, bisnum, entropy
