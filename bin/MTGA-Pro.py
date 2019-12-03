#!usr/bin/python

import re
import os
import sys
import math
import time
import subprocess
import getopt


pwd = sys.path[0] + '/../'

python = pwd + 'soft/python '
SOAPnuke = pwd + 'soft/SOAPnuke ' #soapnuke software abs dirs
Picard = pwd + 'soft/picard.jar '
SOAPdnvt = pwd + 'soft/SOAPdenovoTrans '

BWA = {	
	'bwa': pwd + 'soft/bwa ',
	'aln':' aln -t 1 -L -i 5 -M 9 -R 40  ',
	'samse':' samse -r \"@RG\tID:foo\tSM:bar\" ',
#	'sampe':' sampe -r \"@RG\\tID:foo\\tSM:bar\\tLB:lib\" ',
	'sampe':' sampe -r \"@RG\tID:YANWEI\tLB:AGING\tSM:HEALTH\tPL:ILLUMINA\" ',
#	'sampe': 'sampe -r \"@RG\tID:Sample1\tSM:Sample1\tCN:SequencingCenter\tDS:Project1\tDT:2015-02-09\tPL:ILLUMINA\" ',
	'indexr': pwd + 'ref/mtgenome_rline.fa ',
	'indexrfai': pwd + 'ref/mtgenome_rline.fa.fai ',
	'indexd': pwd + 'ref/mtgenome_dloop.fa ',
	'indexdfai': pwd + 'ref/mtgenome_dloop.fa.fai ',}
SAM = {	'sam': pwd + 'soft/samtools-pro ',#last batch samtools 
	'view':' view -b -S -t ',#view paratmeter
	'sort':' sort -m 300000000 -o ',#sort parameter
	'stats':' stats ',#statistics parameter
	'flagstat':' flagstat ',#flag statistics paratmeter
	'mpileup':' mpileup -d500000 ',#genarete the high-depth mpileup file
	}

def TimeReport(Notice, Iserror):
	if Iserror:
		print('#TimeNode:' + time.strftime('%c',time.localtime())+ '\n#ErrorsNode:' + Notice)
	else:
		print('#TimeNode:' + time.strftime('%c',time.localtime())+ '\n#ReportNode:' + Notice)

class SoapFilter:#using soapnuke filter PE raw fastq file
	'SOAPnuke filter' #class's function describe
	def __init__(self, Input, Output, Sample, Parameter):
		self.Input = Input #input file abs dirs
		self.Output = Output #output file abs dirs
		self.Sample = Sample #launch sample name
#		self.parameter = ' filter -l 15 -m 15 -d -q 0.5 -n 0.1 ' #filter parameter
		self.parameter = ' filter -m 5 '
	def filter(self): #main function name
		if os.path.exists(self.Input + self.Sample + '_1.fq.gz') and os.path.exists(self.Input + self.Sample + '_2.fq.gz'): #determine the input file wether exists or not
			infq1 = self.Input + self.Sample + '_1.fq.gz ' #raw fastq1
			infq2 = self.Input + self.Sample + '_2.fq.gz ' #raw fastq2 
			outfq1 = self.Output + self.Sample + '_1.fq.gz ' #clean fastq1
			outfq2 = self.Output + self.Sample + '_2.fq.gz ' #clean fastq2
			os.system(SOAPnuke + self.parameter + ' -1 ' + infq1 + ' -2 '  + infq2 + ' -Q 2 -G -o ' + self.Output + ' -C ' + outfq1 + ' -D ' + outfq2)#system command launch the initial step: filter
			TimeReport('SOAPnuke filter the PE data is done!', 0)#print the finished time management log
		elif os.path.exists(self.Input + self.sample + '.fq.gz'): #launched SE filter step
			infq0 = self.Input + self.Sample + '.fq.gz' #raw single fastq
			outfq0 = self.Output + self.Sample + '.fq.gz '#clean single fastq
			os.system(SOAPnuke + self.parameter + ' -1 ' + infq0 + ' -Q 2 -G -o ' + self.Output + ' -C ' + outfq0) #launched the initial step to SE rawdata
			TimeReport('SOAPnuke filter the SE data is done!', 0) #print the finished time management log
		else:
			TimeReport('Pelease check Rawdata whether exists or not!', 1) #print Error log
	def stats(self):
		if os.path.exists(self.Output+'') and os.path.exists(self.Output+''):
			TimeReport('SOAPnuke filter the SE data is done!', 0)
		else:
			TimeReport('Pelease check Rawdata whether exists or not!', 1)	
class FastQC:
	'Fastqc testing' #class function  describe
	def __init__(self, Fastqc, Input, Output, Sample, parameter):#function input argv
		self.Fastqc = Fastqc #software name
		self.Input = Input #import dirs
		self.Output = Output #output dirs 
		self.Sample = Sample #sample name
		if len(parameter)==0:
			self.parameter = ' ' #initial parameter
		else:
			self.parameter = parameter # import parameter
	def fastqc(self):
		if os.path.exists(self.Output + self.Sample + '_1.fq.gz ') and os.path.exists(self.Output + self.Sample + '_2.fq.gz '): #
			infq1 = self.Input + self.Sample + '_1.fq.gz '
			infq2 = self.Input + self.Sample + '_2.fq.gz '
			outfq1 = self.Output + self.Sample + '_1.fq.gz '
			outfq2 = self.Output + self.Sample + '_2.fq.gz '
			os.system(self.Fastqc + self.parameter + outfq1)
			os.system(self.Fastqc + ' ' + outfq2)
			TimeReport('Fastqc QC the clean PE data is done!', 0)
			
		elif os.path.exists(self.Output + self.Sample + '.fq.gz '):
			outfq0 = self.Output + self.Sample + '.fq.gz '
			os.system(self.Fastqc + self.parameter + outfq0)
			TimeReport('Fastqc QC the clean SE data is done!', 0)
		else:
			TimeReport('Fastqc is Failed!Please check cleandata whether exists or not!', 1)
class BwaMapper:
	'BWA alingnment' #using BWA mapping onto the reference fasta file
	def __init__(self, Input, Output, Sample, parameter):
		self.Input = Input
		self.Output = Output
		self.Sample = Sample
		self.Parameter = parameter
	def Mapper(self):
		if os.path.exists(self.Output + self.Sample + '_1.fq.gz') and os.path.exists(self.Output + self.Sample + '_2.fq.gz'):
			outfq1 = self.Output + self.Sample + '_1.fq.gz '
			outfq2 = self.Output + self.Sample + '_2.fq.gz '
			fifo1 = self.Output + self.Sample + '_1.fifo '
			fifo2 = self.Output + self.Sample + '_2.fifo '
			prefix = self.Output + self.Sample + '.bam '
			if self.Parameter:
				os.system('mkfifo ' + fifo1 + fifo2)
				os.system(BWA['bwa'] + BWA['aln'] + BWA['indexr'] + outfq1 + '>' + fifo1 + ' & ' + 
					  BWA['bwa'] + BWA['aln'] + BWA['indexr'] + outfq2 + '>' + fifo2 + ' & ' + 
					  BWA['bwa'] + BWA['sampe'] + BWA['indexr'] + fifo1 + fifo2 + outfq1 + outfq2 + ' | ' + 
#					  SAM['sam'] + SAM['view'] + BWA['indexrfai'] + '>' + prefix)
					  SAM['sam'] + SAM['view'] + BWA['indexrfai'] + ' - | ' +
					  SAM['sam'] + SAM['sort'] + prefix) #using pipe finished the PE aligment to genarate .bam file
			else:
				os.system('mkfifo ' + fifo1 + fifo2)
				os.system(BWA['bwa'] + BWA['aln'] + BWA['indexd'] + outfq1 + '>' + fifo1 + ' & ' + 
					  BWA['bwa'] + BWA['aln'] + BWA['indexd'] + outfq2 + '>' + fifo2 + ' & ' + 
					  BWA['bwa'] + BWA['sampe'] + BWA['indexd'] + fifo1 + fifo2 + outfq1 + outfq2 + ' | ' + 
					  SAM['sam'] + SAM['view'] + BWA['indexdfai'] + ' - | ' +
					  SAM['sam'] + SAM['sort'] + prefix) #using pipe finished the PE aligment to genarate .bam file
				
			TimeReport('BWAmem the ref is done!', 0)
		else:
			TimeReport('BWAmem is Failed!Please check cleandata whether exists or not!', 1)
class bqsrBam:
	def __init__(self, Input, Output, Sample, Parameter):
		self.Input = Input
		self.Output = Output
		self.Sample = Sample
		self.Parameter = Parameter
	def bqsr(self):
		if self.Parameter == 1:
#			prefBam = self.Output + self.Sample + '.bam'
			os.system(python + pwd + 'scr/gatkbqsr.py ' + ' '.join([self.Output, self.Output, self.Sample]))
			os.system('rm ' + self.Output + self.Sample + '_gatk*')
		else:
			print('GATK BQSR step is skipped')
class pileupBam:
	'filter multi-mismatch reads'
	def __init__(self, Input, Output, Sample, Parameter):
		self.Input = Input
		self.Output = Output
		self.Sample = Sample
		self.Parameter = Parameter
	def pileup(self):
		if self.Parameter == 1 :
			prefBam = self.Output + self.Sample + '.bam'
			os.system(SAM['sam'] + SAM['stats']     + prefBam + ' > ' + self.Output + self.Sample + '_rawmapstat.txt')
			os.system(SAM['sam'] + SAM['flagstat']  + prefBam + ' > ' + self.Output + self.Sample + '_rawflgstat.txt')
			file1 = self.Output + self.Sample + '_flag1.txt'
			file2 = self.Output + self.Sample + '_flag2.txt'
			file3 = self.Output + self.Sample + '_flag3.txt'
			os.system('awk -F\' \' \'{print $1}\' ' + self.Output + self.Sample + '_rawflgstat.txt > ' + file1)
			os.system('grep ^SN ' + self.Output + self.Sample + '_rawmapstat.txt' + '| cut -f 2-|awk -F\'\t\' \'{print $2}\' > ' + file2 )
			os.system('awk  -F\'(\' \'{print $2}\' ' + self.Output + self.Sample + '_rawflgstat.txt' + '|awk -F\':\' \'{print $1}\'|grep % > ' + file3)
			os.system('cat ' + file1 + ' ' + file2 + ' ' + file3 + ' > ' + self.Output + self.Sample  + '_stats.txt' )
			os.system('rm ' + file1 + ' ' + file2 + ' ' + file3)
		
			rpileup = self.Output + self.Sample + '_r.pileup'
			os.system(SAM['sam'] + ' view -bS -h -F 12 -T ' + BWA['indexr'] + prefBam + ' | ' + SAM['sam'] + ' sort -o ' + self.Output + self.Sample + '_rr.bam')
			os.system(SAM['sam'] + ' rmdup ' + self.Output + self.Sample + '_rr.bam ' + self.Output + self.Sample + '_r.bam')
			os.system('rm ' + self.Output + self.Sample + '_rr.bam')
			os.system(SAM['sam'] + SAM['mpileup'] + self.Output + self.Sample + '_r.bam ' + ' > ' + rpileup)
			dpileup = self.Output + self.Sample + '_d.pileup'
	
			prefBam = self.Output + self.Sample + '.bam'
			unmap0 = self.Output + self.Sample + '_unmap0.bam'
			unmap1 = self.Output + self.Sample + '_unmap1.bam'
			unmap2 = self.Output + self.Sample + '_unmap2.bam'
			os.system(SAM['sam'] + ' view  -b -f 12 ' + prefBam + ' > ' + unmap0)
			os.system(SAM['sam'] + ' view  -b -f 4 -F 8 ' + prefBam + ' > ' + unmap1)
			os.system(SAM['sam'] + ' view  -b -f 8 -F 4 ' + prefBam + ' > ' + unmap2)
			os.system(SAM['sam'] + ' merge -f ' + self.Output + self.Sample + '_unmap.bam ' + ' '.join([unmap0,unmap1,unmap2]))
			os.system(SAM['sam'] + ' sort -n -O SAM ' + self.Output + self.Sample + '_unmap.bam |grep -v "^@" > ' + self.Output + self.Sample + '_unmap.sam')

			os.system(python + pwd + ' scr/bam2fastq.py ' + ' '.join([self.Output, self.Output, self.Sample + '_unmap']))
		
			BM = BwaMapper(self.Input, self.Output, self.Sample + '_unmap', 0)
			BM.Mapper()

			os.system(SAM['sam'] + ' rmdup ' + self.Output + self.Sample + '_unmap.bam -| ' + SAM['sam'] + ' sort -o ' + self.Output + self.Sample + '_d.bam')
			os.system(SAM['sam'] + SAM['mpileup'] + self.Output + self.Sample + '_d.bam' + ' > ' + dpileup)
			os.system('rm ' + self.Output + self.Sample + '_unmap*')
		else:
			prefBam = self.Output + self.Sample + '.bam'
			os.system(SAM['sam'] + SAM['stats']     + prefBam + ' > ' + self.Output + self.Sample + '_rawmapstat.txt')
			os.system(SAM['sam'] + SAM['flagstat']  + prefBam + ' > ' + self.Output + self.Sample + '_rawflgstat.txt')

			file1 = self.Output + self.Sample + '_flag1.txt'
			file2 = self.Output + self.Sample + '_flag2.txt'
			file3 = self.Output + self.Sample + '_flag3.txt'
			os.system('awk -F\' \' \'{print $1}\' ' + self.Output + self.Sample + '_rawflgstat.txt > ' + file1)
			os.system('grep ^SN ' + self.Output + self.Sample + '_rawmapstat.txt' + '| cut -f 2-|awk -F\'\t\' \'{print $2}\' > ' + file2 )
			os.system('awk  -F\'(\' \'{print $2}\' ' + self.Output + self.Sample + '_rawflgstat.txt' + '|awk -F\':\' \'{print $1}\'|grep % > ' + file3)
			os.system('cat ' + file1 + ' ' + file2 + ' ' + file3 + ' > ' + self.Output + self.Sample  + '_stats.txt' )
			os.system('rm ' + file1 + ' ' + file2 + ' ' + file3)

			rpileup = self.Output + self.Sample + '.pileup'
			os.system(SAM['sam'] + ' view -bS -h -F 12 -T ' + BWA['indexr'] + prefBam + ' | ' + SAM['sam'] + ' sort -o ' + self.Output + self.Sample + '_rr.bam')
			os.system(SAM['sam'] + ' rmdup ' + self.Output + self.Sample + '_rr.bam ' + self.Output + self.Sample + '.bam')
			os.system('rm ' + self.Output + self.Sample + '_rr.bam')
			os.system(SAM['sam'] + SAM['mpileup'] + self.Output + self.Sample + '.bam ' + ' > ' + rpileup)
def coveragePic(path, name, position, coverage):
	import matplotlib as mpl
	mpl.use('Agg')
	import matplotlib.pyplot as plt
	from matplotlib.font_manager import FontProperties

	x,y = position, coverage
	maxdepth, mindepth, avgdepth = min(y), max(y), round(sum(y)/len(y),2)

	ttfont = {
                'family' : 'serif',
                'color'  : 'black',
                'weight' : 'normal',
                'size'   : 20,
        }
	xyfont = {
                'family' : 'serif',
                'color'  : 'black',
                'weight' : 'normal',
                'size'   : 15,
                }
	txfont = font = {'family': 'serif',
                'color':  'darkred',
                'weight': 'normal',
                'size': 5,
                }

	fig = plt.gcf()
	fig.set_size_inches(18,10)
	plt.style.use('ggplot')	
	plt.plot(position,coverage,color="slategrey",markersize=1)
	plt.xlim((-10, 16570))
	plt.text(13000,  1.5 * avgdepth, 'MaxDepth: ' + str(maxdepth) + '\nMinDepth: ' + str(mindepth) + '\nAverageDepth: ' + str(avgdepth) ,style='italic', bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
	plt.suptitle('Mitochondrial Genome Coverage\n' + name , fontdict = ttfont)
	plt.ylabel('Depth',fontdict = xyfont)
	plt.xlabel('Position',fontdict = xyfont)
	fig.savefig( path + name + '.png',dpi=100)

class mergePileup:
	'merge two loop*.pileup file'
	def __init__(self, Input, Output, Sample, Parameter):
		self.Input = Input
		self.Output = Output
		self.Sample = Sample
		self.Parameter = Parameter
	def merge(self):
		rpileup = self.Output + self.Sample + '_r.pileup'
		dpileup = self.Output + self.Sample + '_d.pileup'
		R = open(rpileup,'rb').readlines()
		D = open(dpileup,'rb').readlines()
		M = open(self.Output + self.Sample + '.pileup','wb')
		d = {}
		for i in D:
			if len(i.split()) == 6:
				mt, pos, ref, depth, base, qual = i.split()
				pos, depth = int(pos), int(depth)
				if int(pos) > 546:
					pos = pos - 546
				else:
					pos = pos + 16023
				d[pos] = [depth, base, qual]
			else:
				continue
		position, coverage = [],[]
		for j in R:
			if (len(j.split())) == 6:
				mt, pos, ref, depth, base, qual = j.split()
				pos, depth = int(pos), int(depth)
				if pos in d.keys():
					depth = d[pos][0] + depth
					base = d[pos][1] + base
					qual = d[pos][2] + qual
					M.write( '\t'.join([mt, str(pos), ref, str(depth), base, qual]) + '\n')
					position.append(pos)
					coverage.append(depth)
				else:
					position.append(pos)
					coverage.append(depth)
					M.write( '\t'.join([mt, str(pos), ref, str(depth), base, qual]) + '\n')
			else:
				continue
		
		M.close()
		coveragePic(self.Output, self.Sample, position, coverage)
		
class heterCall:
	'using .pileup file to calling heteroplasmy'
	def __init__(self, Input, Output, Sample, Parameter,):
		self.Input = Input
		self.Output = Output
		self.Sample = Sample
		self.Parameter = Parameter
	def filterSBias(self):
		OSP = self.Output + self.Sample
		SHV = self.Output + self.Sample + '.shv '
		os.system(python + pwd + 'scr/readpileup.py ' + ' '.join([self.Sample, self.Input, self.Output, self.Parameter]) + ' > ' + SHV)
class mutationAnnot:
	'using OMIM and Genebank and Pathegenic database to annotate the variation'
	def __init__(self, Input, Output, Sample, Parameter):
		self.Input = Input
		self.Output = Output
		self.Sample = Sample
		self.Parameter = Parameter
	def annotate(self):
		SHV = self.Output + self.Sample + '.shv'
		ANN = self.Output + self.Sample + '.ann'
		for i in open(SHV,'rb').readlines():
			if not re.match("^#",i):
				I = i.strip().split(',')
				os.system(python + pwd + 'scr/findanno.py ' + ' '.join(I[1:4]) + ' ' + pwd + '>>' + ANN)
		os.system('paste -d\',\' ' + SHV + ' ' + ANN + '> ' + self.Output + self.Sample + self.Parameter+'.hetvar')
		os.system('rm ' + SHV + ' ' + ANN)
class TransAssembly:
	'using SOAPdenovo to assembly the Ultr-Depth MT sequencing'
	def __init__(self, Input, Output, Sample, Parameter):
		self.Input = Input
		self.Output = Output
		self.Sample = Sample
		self.Parameter = Parameter
	def mtAssembly(self):
		infq1 = self.Output + self.Sample + '_1.fq.gz '
		infq2 = self.Output + self.Sample + '_2.fq.gz '
		config = self.Output + self.Sample + '.trans.config'
		assemb = self.Output + self.Sample + '.trans'
		CFG = open(config,'wb')
		CFG.write('#mal read length\nmax_rd_len=50\n[LIB]\n#average insert size\navg_ins=200\n' + 'q1='+infq1 + '\nq2=' + infq2 + '\n')
		CFG.close()
		os.system(SOAPdnvt + ' all -s '+ config + ' -o ' + assemb + ' -K ' + self.Parameter)
class PhylotreeMer:
	'using Phy-Mer to call the mitochondria haplotype'
	def __init__(self, Input, Output, Sample, Parameter):
		self.Input = Input
		self.Output = Output
		self.Sample = Sample
		self.Parameter = Parameter
	def haplocall(self):
		outbam = self.Output + self.Sample + '.bam'
		script = pwd + 'scr/Phy-Mer.py '
		haplodatabase = pwd + 'data/Haplogroup_motifs.csv '
		phylodatabase = pwd + 'data/PhyloTree_b16_k12.txt '
		mt_haplogroup = self.Output + self.Sample + '.haplogroup'
		os.system(python + script + '--verbose --def-snps=' + haplodatabase + phylodatabase + outbam + '>' +  mt_haplogroup)


def print_help():
        print('''
MTGA-Pro - MiTochondrial Genome Analysis Project  (2016 Oct 15, compiled Thu 17:10:58 Aug 10 2017)

Introduction:
	1.Function: Variantion Calling & Haplopgroup Calling & Genome Assembly for BGISEQ500 Sequencing Data
	2.Copy Right: BGI Research Institute 2017
	3.Auther: Yanwei QI [qiyanwei@genomics.cn [FirstChoice] & qiyanweii@icloud.com]
	4.Usage: python MTGA-Pro.py 
			--input      = ../test/ 
			--output     = ../test/ 
			--sample     = CL100020393_L02_18 
			--isfilter   = 1 
			--ismapping  = 1 
			--isbqsr     = 1 
			--ispileup   = 1,1
			--isheter    = 1,1
			--isassembly = 1,45
			--ishaplo    = 1,45
''')

def print_usage():
	print('''
Argument:
	-i,--input		[string /dir/] 	input rawdata path
	-o,--output		[string /dir/] 	output result path
	-s,--sample		[string Name]  	sample name
	-f,--isfilter		[int 0/1] 	Running SOAPnuke Filter
	-m,--ismapping		[int 0/1]	Mapping on reference
	-b,--isbqsr		[int 0/1]	Running GATK to BQSR the BAM file
	-l,--ispileup		[int 0/1,int 0/1]	The first number is make .pileup file, 
								the second number is merge dloop and rCRS reference.
	-t,--isheter		[int 0/1,float 1/5]	The first number is calling heteroplasmy, 
								the second number is defind the frequency of heteroplasmy.
	-a,--isassembly		[int 0/1,int 45/47]	Running SOAPdenovo-Trans to assembly the Mitochondrial genome.
	-p,--ishaplo		[int 0/1,int 12/45]	Running Phy-Mer to calling the haplogroup of sample.
	-h,--help
	-u,--usage
''')

try:
	opts, args = getopt.getopt(sys.argv[1:],"vhu:i:o:s:f:m:c:t:a:p",["version","help","usage","input=","output=","sample=","isfilter=","ismapping=","isbqsr=","ispileup=","isheter=","isassembly=","ishaplo="])
except getopt.GetoptError as err:
	print_help()
	print_usage()
	sys.exit(1)

for arg, par in opts:
	if arg in ('-i','--input'):
		Input = par
	if arg in ('-o','--output'):
		output = par
	if arg in ('-s','--sample'):
		Sample = par
	if arg in ('-f','--isfilter'):
		Isfilter = int(par)
	if arg in ('-m','--ismapping'):
		Ismapping = int(par)
	if arg in ('-b','--isbqsr'):
		Isbqsr = int(par)
	if arg in ('-c','--ispileup'):
		Ispileup,circle = par.split(',')
	if arg in ('-t','--isheter'):
		Isheter,freq = par.split(',')
	if arg in ('-a','--isassembly'):
		Isassembly,akmer = par.split(',')
	if arg in ('-p','--ishaplo'):
		Ishaplo,hkmer = par.split(',')
	if arg in ('-u','--usage'):
		print_usage()
		sys.exit(1)
	if arg in ('-h','-help'):
		print_help()
		sys.exit(1)
def main(Input,output,Sample,Isfilter,Ismapping,Isbqsr,Ispileup,circle,Isheter,freq,Isassembly,akmer,Ishaplo,hkmer):
	Output = output + Sample + '/'
	if not os.path.exists(Output):
		os.system('mkdir ' + Output)
	if Isfilter:
		SF = SoapFilter(Input, Output, Sample, '')
		SF.filter()
	if Ismapping:
		BM = BwaMapper(Input, Output, Sample, 1 )
		BM.Mapper()
	if Isbqsr:
		BB = bqsrBam(Input, Output, Sample, 1)
		BB.bqsr()
	if int(Ispileup):
		SB = pileupBam(Input, Output, Sample, int(circle))
		SB.pileup()
	if int(circle):
		MP = mergePileup(Input, Output, Sample, int(circle))
		MP.merge()
	if int(Isheter):
		HC = heterCall(Input, Output, Sample, freq )
		HC.filterSBias()
		SA = mutationAnnot(Input, Output, Sample, '')
		SA.annotate()
	if int(Isassembly):
		TA = TransAssembly(Input, Output, Sample, akmer)
		TA.mtAssembly()
	if int(Ishaplo):
		PM = PhylotreeMer(Input, Output, Sample, hkmer)
		PM.haplocall()
	print(Input,output,Sample,Isfilter,Ismapping,Isbqsr,Ispileup,circle,Isheter,freq,Isassembly,akmer,Ishaplo,hkmer)

if __name__ == '__main__':
	main(Input,output,Sample,Isfilter,Ismapping,Isbqsr,Ispileup,circle,Isheter,freq,Isassembly,akmer,Ishaplo,hkmer)
