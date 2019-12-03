#!usr/bin/python
import os
import re
import sys
import time
import math
import getopt


class ANNOTATION:
	'Annotation Mutations by GenBank and OMIM database'
	REFERENCE, GENEBANK, OMIMINFO, CODONTAB, = {}, {}, {}, {}
	ANNO = []
#	def __init__(self, reffasta, omimdata, genebank):#
	def __init__(self, position, majorBase, minorBase, reffasta, omimdata, genebank):
		self.pos 	= 	int(position)
		self.major 	= 	majorBase
		self.minor 	= 	minorBase
		self.ref 	= 	reffasta
		self.omim 	= 	omimdata
		self.bank 	= 	genebank
	def refDict(self):
		self.REFERENCE = re.split('\n', open(self.ref).read())[1]
	def omimDict(self):
		for omim in open(self.omim).readlines():
			omim = omim.rstrip()
			if not re.match(r'^site',omim):
				om = re.split('\t',omim)
				if re.match(r'A|G|T|C',om[1]):
					self.OMIMINFO[om[0]] = 	{
								'base' : om[1],
								'gene' : om[2],
								'omim' : om[3],
								'disease' : om[4],
								}
				else:
					continue
	def bankDict(self):
		GBText = open(self.bank).read()
		PAT = re.compile(r'\n\s{5}[a-zA-Z]+\s{12}')
		VAL = re.split(PAT,GBText)
		KEY = PAT.findall(GBText)
		if len(VAL) > len(KEY):
			for i in range(len(KEY)):
				FEATURES = KEY[i] + VAL[i+1]
				FEATURES = re.sub(r'\s\s+','',FEATURES.strip())	#del '\n' and '\s{2+}'
				FEATURES = re.sub(r'\s|\.|-','_',FEATURES)
				[startPos, abortPos, geneName, geneSyn, geneID, chainDirec, proteinID, product, translation] = ['','','','','','+','','','']
				[	patternPOS,
					patternGeneName,
					patternGeneSyn,
					patternGeneID,
					petternChainDirec,
					petternProteinID,
					petternProduct,
					petternTranslation, 	]=[
					re.compile(r'\d+'),
					re.compile(r'gene="(\w+)"'),
					re.compile(r'gene_synonym="(\w+)"'),
					re.compile(r'db_xref="GeneID:(\d+)"'),
					re.compile(r'.+complement.+'),
					re.compile(r'protein_id="(\w+)"'),
					re.compile(r'product="(\w+)"'),
					re.compile(r'translation="(\w+)"'),
				]
				if re.match(petternChainDirec, FEATURES):#c/complement '-'/
					position = patternPOS.findall(FEATURES)
					startPos = position[0]
					abortPos = position[1]
					chainDirec = '-'
					if re.match(patternGeneSyn, FEATURES):
						geneSyn     =   patternGeneSyn.findall(FEATURES)
					else:
						geneSyn     =   patternGeneName.findall(FEATURES)
					if re.match(r'^rRNA',FEATURES):
						geneName    = patternGeneName.findall(FEATURES)
						geneSyn     = patternGeneSyn.findall(FEATURES)
						geneID      = patternGeneID.findall(FEATURES)
						proteinID   = petternProteinID.findall(FEATURES)
						product     = petternProduct.findall(FEATURES)
						self.GENEBANK[startPos] = {	'abortPos' : abortPos,
									'chainDirec' : chainDirec[0],
									'geneName' : geneName[0],
									'geneSyn' : geneSyn[0],
									'geneID' : geneID[0],
									'proteinID' : proteinID,
									'product' : product[0],}
					elif re.match(r'^tRNA',FEATURES):
						geneName    = patternGeneName.findall(FEATURES)
						geneSyn     = patternGeneSyn.findall(FEATURES)
						geneID      = patternGeneID.findall(FEATURES)
						proteinID   = petternProteinID.findall(FEATURES)
						product     = petternProduct.findall(FEATURES)
						self.GENEBANK[startPos] = {	'abortPos' : abortPos,
									'chainDirec' : chainDirec[0],
									'geneName' : geneName[0],
									'geneSyn' : geneSyn[0],
									'geneID' : geneID[0],
									'proteinID' : proteinID,
									'product' : product[0],}
					elif re.match(r'^CDS', FEATURES):
						geneName    = patternGeneName.findall(FEATURES)
						geneSyn     = patternGeneSyn.findall(FEATURES)
						geneID      = patternGeneID.findall(FEATURES)
						proteinID   = petternProteinID.findall(FEATURES)
						product     = petternProduct.findall(FEATURES)
						self.GENEBANK[startPos] = 	{
									'abortPos' : abortPos,
									'chainDirec' : chainDirec[0],
									'geneName' : geneName[0],
									'geneSyn' : geneSyn[0],
									'geneID' : geneID[0],
									'proteinID' : proteinID[0],
									'product' : product[0],}
					else:
						continue
				else:
					position = patternPOS.findall(FEATURES)
					startPos = position[0]
					abortPos = position[1]
					if re.match(patternGeneSyn, FEATURES):
						geneSyn     =   patternGeneSyn.findall(FEATURES)
					else:
						geneSyn     =   patternGeneName.findall(FEATURES)
					if re.match(r'^rRNA',FEATURES):
						geneName    = patternGeneName.findall(FEATURES)
						geneID      = patternGeneID.findall(FEATURES)
						proteinID   = petternProteinID.findall(FEATURES)
						product     = petternProduct.findall(FEATURES)
						self.GENEBANK[startPos] = {	'abortPos' : abortPos,
									'chainDirec' : chainDirec[0],
									'geneName' : geneName[0],
									'geneSyn' : geneSyn[0],
									'geneID' : geneID[0],
									'proteinID' : proteinID,
									'product' : product[0],}
					elif re.match(r'^tRNA',FEATURES):
						geneName    = patternGeneName.findall(FEATURES)
						geneID      = patternGeneID.findall(FEATURES)
						proteinID   = petternProteinID.findall(FEATURES)
						product     = petternProduct.findall(FEATURES)
						self.GENEBANK[startPos] = {	'abortPos' : abortPos,
									'chainDirec' : chainDirec[0],
									'geneName' : geneName[0],
									'geneSyn' : geneSyn[0],
									'geneID' : geneID[0],
									'proteinID' : proteinID,
									'product' : product[0],}
					elif re.match(r'^CDS', FEATURES):
						geneName    = patternGeneName.findall(FEATURES)
						geneID      = patternGeneID.findall(FEATURES)
						proteinID   = petternProteinID.findall(FEATURES)
						product     = petternProduct.findall(FEATURES)
						self.GENEBANK[startPos] = 	{
									'abortPos' : abortPos,
									'chainDirec' : chainDirec[0],
									'geneName' : geneName[0],
									'geneSyn' : geneSyn[0],
									'geneID' : geneID[0],
									'proteinID' : proteinID[0],
									'product' : product[0],}
					else:
						continue
			return self.GENEBANK
	def codonChart(self):		
		self.CODONTAB	= 	{
					'TTT':'F','TTC':'F',
					'CTT':'L','CTC':'L','CTA':'L','CTG':'L','TTA':'L','TTG':'L',
					'ATT':'I','ATC':'I','ATA':'I',
					'ATG':'M',
					'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
					'TCT':'S','TCC':'S','TCA':'S','TCG':'S','AGT':'S','AGC':'S',
					'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
					'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
					'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
					'TAT':'Y','TAC':'Y',
					'CAT':'H','CAC':'H',
					'CAA':'Q','CAG':'Q',
					'AAT':'N','AAC':'N',
					'AAA':'K','AAG':'K',
					'GAT':'D','GAC':'D',
					'GAA':'E','GAG':'E',
					'TGT':'C','TGC':'C',
					'TGG':'W',
					'CGT':'R','CGC':'R','CGA':'R','CGG':'R','AGA':'R','AGG':'R',
					'GGT':'G','GGC':'G','GGA':'G','GGG':'G',
					'TAA':'_','TAG':'_','TGA':'_',
					}
		return self.CODONTAB
	
	def annoCard(self):
		for pos1, val in self.GENEBANK.items():
			start = int(pos1)
			abort = int(val['abortPos'])
			cdseq = self.REFERENCE[start:abort]
			
			if self.pos >= 1 and self.pos <= 576 or self.pos >= 16024 and self.pos <= 16569:
				self.ANNO = [str(self.pos),"D-loop",self.major,self.minor]
				return self.ANNO
				break
			elif start <= self.pos and abort >= self.pos:
				ANNO = val.values()[1:]
				if 'translation' in val.keys():
					location = self.pos - start + 1
					div, mod = divmod(location,3)
					if mod == 0:
						abc = cdseq[(div-1)*3 : (div)*3]
						base1 = abc[0] + abc[1] + self.major
						base2 = abc[0] + abc[1] + self.minor
						aa1 = self.CODONTAB[base1]
						aa2 = self.CODONTAB[base2]
					elif mod == 1:
						abc = cdseq[div*3 : (div+1)*3]
						base1 = self.major + abc[1] + abc[2]
						base2 = self.minor + abc[1] + abc[2]
						aa1 = self.CODONTAB[base1]
						aa2 = self.CODONTAB[base2]
					elif mod == 2:
						abc = cdseq[div*3 : (div+1)*3]
						base1 = abc[0] + self.minor + abc[2]
						base2 = abc[0] + self.minor + abc[2]
						aa1 = self.CODONTAB[base1]
						aa2 = self.CODONTAB[base2]
					else:
						continue
					if aa1 == aa2:
						self.ANNO = ANNO
						self.ANNO = self.ANNO + ['syn', self.major + '>' + self.minor, base1 + '>' + base2, aa1 + '>' + aa2]
					else:
						self.ANNO = ANNO
						self.ANNO = self.ANNO + ['nsy', self.major + '>' + self.minor, base1 + '>' + base2, aa1 + '>' + aa2]
					for pos2,omim in self.OMIMINFO.items():
						if self.pos == int(pos2) and self.minor == omim['base']:
							self.ANNO.append(omim['omim'])
							self.ANNO.append(omim['disease'])
							return self.ANNO
							break
						else:
							continue
				else:
					self.ANNO.append(str(self.pos))
					self.ANNO.extend(ANNO)
					return self.ANNO
					break
			else:
				pass
position = int(sys.argv[1])
refmajor = str(sys.argv[2])
altminor = str(sys.argv[3])
curr_dir = str(sys.argv[4])
CLASS=ANNOTATION(position, refmajor, altminor, curr_dir + 'data/mitochondria.fa',curr_dir + 'data/omim_pathogenic',curr_dir + 'data/mitochondria.gb')
CLASS.refDict()
CLASS.omimDict()
CLASS.bankDict()
CLASS.codonChart()
ANNOT = CLASS.annoCard()
if ANNOT:
	print(','.join([str(i) for i in ANNOT]))
else:
	print('introgenic')
