#!/usr/bin/python
import sys
import pysam
import array
from collections import Counter
from gzip import GzipFile
import pickle
import os
import argparse
from Bio import SeqIO

class MappingStats:
	class VectorData:
		def __init__(self,typ,data):
			self.typ=typ
			if type(data)==str:
				self.data=array.array(self.typ)
				self.data.fromstring(data)
			elif type(data)==list:
				self.data=array.array(self.typ,data)
			else:
				self.data=data

			self.softclip=0
			self.ceil=0

		def ceilMode(self):
			self.ceil=Counter(self.data).most_common()
			if self.ceil[0][0]: self.ceil=self.ceil[0][0]
			else: 
				self.ceil.pop(0)
				if not self.ceil: self.ceil=1
				else: self.ceil=self.ceil[0][0]

		def ceilAvg(self):
			if not self.data: return 0
			self.ceil=sum(self.data)/len(self.data)

		def ceilMax(self):
			if not self.data: return 0
			self.ceil=max(self.data)


		def vectorSummary(self,l,normalized=True,binary=False):
			'''Normalized will transform data to values from 0 to 1
			   binary will interpret original positions as 0 or 1, 1 meaning different than 0 (so any value)
			   the binary flag will convert the vector to 1's and 0's and then averages will be performed 
			   leading to values between 0 and 1 '''

			if binary: normalized=False
			
			res=[i for i in (float(sum(i))/len(i) for i in self.chunks(self.data,l,binary))]
			if not normalized: return res
			return [min(float(i)/self.ceil,1) for i in res]


		def toList(self):
			return [self.typ,self.data.tostring(),self.ceil]

		@staticmethod
		def chunks(l, n,binary):
			"""Yield successive n chunks from l."""
			chunksize = float(len(l))/n
			intchunk=int(chunksize)
			if not chunksize.is_integer(): intchunk+=1

			if binary:
				for i in range(n):
					pos=int(chunksize*i)
					yield [1 if j else 0 for j in l[pos:pos+intchunk]]
			else:
				for i in range(n):
					pos=int(chunksize*i)
					yield l[pos:pos+intchunk]

		@staticmethod
		def fromList(l):
			typ,data,ceil=l
			v=MappingStats.VectorData(typ,data)
			v.ceil=ceil
			return v


		def mean(self):
			return sum(self.data)/float(len(self.data))
		
		def median(self):
			lst=self.data
			lst = sorted(lst)
			if len(lst) < 1:
					return None
			if len(lst) %2 == 1:
					return lst[((len(lst)+1)/2)-1]
			else:
					return float(sum(lst[(len(lst)/2)-1:(len(lst)/2)+1]))/2.0
		
		def mode(self):
			return Counter(self.data).most_common(1)[0][0]

		def __gt__(self,val):
			assert type(val) in [int,float]
			lst=self.data
			vals=[i for i in self.data if i>val]
			return float(len(vals))/len(lst)*100


		def __lt__(self,val):
			assert type(val) in [int,float]
			lst=self.data
			vals=[i for i in self.data if i<val]
			return float(len(vals))/len(lst)*100



	def __init__(self,sampleid="NA"):
		self.data={'coverages':{},'mapq':{},'highcov':{},'called':{},'mutated':{},'hz':{}}
		self.ID=sampleid;
		self.ref=None

	def processBam(self,bam):
		if self.ID=="NA": self.ID=os.path.basename(bam)
		self.bam=pysam.Samfile(bam)
		self.readBam()

	def covStats(self,ofile):
		ofile=open(ofile,"w")
                covs=[1,20,50,100,150,200]
		ofile.write("ID\tREF\tmean\tmedian\tmode\tsoftclip\t")
                ofile.write('\t'.join([str(i) for i in covs])+'\n')

		for ch,data in self.data["coverages"].items():
                    strwr="{0}\t{1}\t{2:.1f}\t{3:.1f}\t{4}\t{5:.1f}\t"+"\t".join(["{"+str(i)+":.1f}" for i in range(6,len(covs)+6)])+'\n'
                    vals=[self.ID,ch,data.mean(),data.median(),data.mode(),self.data['softclip'][ch]/float(len(data.data))]
                    for i in covs:
                        vals.append(data>(i-1))
                    ofile.write(strwr.format(*vals))


#		for ch,data in self.data["coverages"].items():
#		    ofile.write("{0}\t{1}\t{2:.1f}\t{3:.1f}\t{4}\t{5:.1f}\t{6:.1f}\t{7:.1f}\n".format(self.ID,ch,data.mean(),data.median(),data.mode(),data>0,data>4,data>9))
		ofile.close()

	def varStats(self,ofile):
		pass
		

	def processVcf(self,vcf,ref):
		if self.ID=="NA": self.ID=os.path.basename(vcf)
		if ref: 
			self.ref=dict([(i.name,str(i.seq).upper()) for i in SeqIO.parse(ref,"fasta")])
		if vcf.split(".")[-1].lower() in ["gz",'gzip']:
			self.vcf=GzipFile(vcf)
		else:
			self.vcf=open(vcf)
		self.readVcf()

	def vSum(self,dtype,chromosome,size,norm=True,binary=False):
		'''Normalized will transform data to values from 0 to 1
		   binary will interpret original positions as 0 or 1, 1 meaning different than 0 (so any value)
		   the binary flag will convert the vector to 1's and 0's and then averages will be performed 
		   leading to values between 0 and 1, this is useful for when you are not interested in values but 
		   just PRESENCE of value'''

		assert dtype in self.data and chromosome in self.data[dtype]
		return self.data[dtype][chromosome].vectorSummary(size,norm,binary)

	def switchData(self):	
		for datatype,data in self.data.items():
			for chromosome in data.keys():
				if type(data[chromosome])!=list: data[chromosome]=data[chromosome].toList()
				else: data[chromosome]=MappingStats.VectorData.fromList(data[chromosome])

	def store(self,f):
		self.switchData()
		f=GzipFile(f,"w")
		pickle.dump([self.ID,self.data],f)
		f.close()
		self.switchData()  

	def load(self,f):
		f=GzipFile(f)
		self.ID,self.data=pickle.load(f)
		f.close()
		self.switchData()

	def chromosomes(self):
	    return set(self.data["coverages"].keys()+self.data["mapq"].keys()+self.data["highcov"].keys()+self.data["called"].keys())

	@staticmethod
	def zeros(n):
		for i in xrange(n):
			yield 0

	@staticmethod
	def mapq2array(d,length):
		buckets=[(pos,sum(l)/len(l)) for pos,l in d.items()]
		data=array.array('b',MappingStats.zeros(length))
		for pos,val in buckets:
			for i in xrange(pos,min(pos+100,length)):
				data[i]=val
	
		return data

	def extractField(self,f,string,sep):
		fpos=string.index(f)+len(f)
		endpos=string.find(sep,fpos)
		return string[fpos:endpos]


	def readVcf(self):
	    #return set(self.data["coverages"].keys()+self.data["mapq"].keys()+self.data["highcov"].keys()+self.data["called"].keys())
		for i in self.vcf:
			if i[0]=="#": continue
			ch=i[:i.index('\t')]

			if ch not in self.data["coverages"]:
			    self.data["coverages"][ch]=MappingStats.VectorData('H',[])
			    self.data["mapq"][ch]=MappingStats.VectorData('b',[])
			    self.data["highcov"][ch]=MappingStats.VectorData('H',[])
			    self.data["called"][ch]=MappingStats.VectorData('b',[])
			    self.data["mutated"][ch]=MappingStats.VectorData('b',[])
			    self.data["hz"][ch]=MappingStats.VectorData('b',[])

			bc=self.extractField("BaseCounts=",i,';')
			bc4=self.extractField("BaseCounts4=",i,';')
			bc4=sum([int(j) for j in bc4.split(",")])
			self.data["highcov"][ch].data.append(bc4)
			bc=sum([int(j) for j in bc.split(",")])
			self.data["coverages"][ch].data.append(bc)
			mq=int(self.extractField("MQ=",i,';'))
			self.data["mapq"][ch].data.append(mq)
			bcall=self.extractField("BCALL=",i,';')
			if "PASS" in i and bcall in "ACTG": 
			    self.data["called"][ch].data.append(1)
			    called=True
			else: 
			    self.data["called"][ch].data.append(0)
			    called=False
			if 'z\t' in i: self.data["hz"][ch].data.append(1)
			else: self.data["hz"][ch].data.append(0)

			
			if self.ref:
				if not called:
				    self.data["mutated"][ch].data.append(0)
				    continue
				pos=i.find("\t",len(ch)+1)
				pos=int(i[len(ch)+1:pos])
				if bcall!=self.ref[ch][pos-1]: self.data["mutated"][ch].data.append(1)
				else: self.data["mutated"][ch].data.append(0)

		for i in self.data["coverages"].values():
			i.ceilMode()
		for i in self.data["mapq"].values():
			i.ceilMax()
		for i in self.data["highcov"].values():
			i.ceilMode()
		for i in self.data["mapq"].values():
			i.ceilMax()
		for i in self.data["called"].values():
			i.ceilMax()
		for i in self.data["mutated"].values():
			i.ceilMax()


	def readBam(self):
		coverages=[array.array('H',self.zeros(i)) for i in self.bam.lengths]
		mapq={}
		
		softclip={}

		for al in self.bam:
			if al.is_unmapped: continue
			ref=al.tid
			pos=al.pos
                        minref=99999999999
                        maxref=0

                        for _,rpos in al.aligned_pairs:
                                if rpos and rpos<minref:
                                    minref=rpos
                                if rpos and rpos>maxref:
                                    maxref=rpos
                        for i in xrange(minref,maxref+1):
                            coverages[ref][i]+=1

			refname=self.bam.references[ref]
			softclip.setdefault(refname,0)
			for code,l in al.cigar:
				if code==4:
					softclip[refname]+=l
					

		
			pos=pos-(pos%100)
			if not ref in mapq: mapq[ref]={}
			if pos not in mapq[ref]: mapq[ref][pos]=[]
			mapq[ref][pos].append(al.mapq)
		
		for i,j in zip(self.bam.references,coverages):
			self.data["coverages"][i]=MappingStats.VectorData("H",j)
			self.data["coverages"][i].ceilMode()
		
		mapq=[i[1] for i in sorted(mapq.items())]
		for i in range(len(mapq)):
			mapq[i]=self.mapq2array(mapq[i],self.bam.lengths[i])
		
		for i,j in zip(self.bam.references,mapq):
			self.data["mapq"][i]=MappingStats.VectorData("b",j)
			self.data["mapq"][i].ceilMax()

		self.data['softclip']=softclip


if __name__=="__main__":
	parser = argparse.ArgumentParser(description='Mapping Statistics')
	parser.add_argument("-s", dest="SampleID",help="SampleId",required=False,default='NA')
	parser.add_argument("-b", dest="bam",help="Bam file",required=False)
	parser.add_argument("-v", dest="Vcf",help="Vcf file",required=False)
	parser.add_argument("-r", dest="reference",help="fastafile (if you want to get mutations)",required=False)
	parser.add_argument("-o", dest="heatout",help="HeatMaps output file",required=False)
	parser.add_argument("-cs", dest="covStats",help="Coverage Stats (bam or VCF)",required=False)
	parser.add_argument("-vs", dest="variantstats",help="Variant statistics (only vcf)",required=False)
	options = parser.parse_args()

	if not options.bam and not options.Vcf: 
	    print "You must specifiy at least a bam or a vcf file!"
	    sys.exit(-1)

	if options.bam and options.Vcf:
	    print "You can specify either bam or a vcf but not both!"
	    sys.exit(-1)
	
	if not options.heatout and not options.covStats and not options.variantstats:
	    print "You must specifiy some output!"
	    sys.exit(-1)

	if (options.bam and options.variantstats):
	    print "Wrong parameters! you need vcf for variantStats"
	    sys.exit(-1)

	a=MappingStats(options.SampleID)

	if options.bam: a.processBam(options.bam)
	if options.Vcf: a.processVcf(options.Vcf,options.reference)

	if options.heatout: a.store(options.heatout)
	if options.covStats: a.covStats(options.covStats)
	if options.variantstats: a.varStats(options.variantstats)
    
#	a.load("deleteme")
#
#	#print vectorSummary(data["unitig_2|quiver"]["mapq"],100)
#	import code
#	code.interact(local=locals())

