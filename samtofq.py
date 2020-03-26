#!/usr/bin/python
import pysam
import argparse
import subprocess
from fqreader import *
from string import maketrans
from threading import Thread

class SamToFq:
	def __init__(self,inbam,fout1,fout2=None):
		self.out1=FastQWriter(fout1)	
		if not fout2:
			self.out2=self.out1
		else:
			self.out2=FastQWriter(fout2)	
	
		self.inp=pysam.Samfile(inbam)

	def run(self):
		self.thread=Thread(target=self.worker,kwargs={})
		self.thread.start()
		

	def join(self):
		self.thread.join()

	def worker(self):
		transdna=maketrans("ACTGactg","TGACtgac")

		reads={}
		for i in self.inp:
			if i.flag&2048: continue
	
			r1=i
	
			if not r1.flag&16:
				self.out1.writeRecord(["","@"+r1.qname.split("/")[0]+"/1",r1.seq,"+",r1.qual])
#			else:
#				self.out1.writeRecord(["","@"+r1.qname.split("/")[0]+"/1",r1.seq.translate(transdna)[::-1],"+",r1.qual[::-1]])
			
	
		self.out1.close()
		if self.out2!=self.out1: self.out2.close()

if __name__=="__main__":
	parser = argparse.ArgumentParser(description='Fix FastqFiles')
	parser.add_argument("-i", dest="inp",help="Sam/Bam input (default stdin)",default="-")
	parser.add_argument("-o", dest="out1",help="Output file (default stdout)",default="-")
	parser.add_argument("-o2", dest="out2",help="Second output file (if not shuffled file)",default=None)
	parser.add_argument("-hd", dest="header",help="Text file where you want to dump the header information",default=None)
	args = parser.parse_args()


	q=SamToFq(args.inp,args.out1,args.out2)
	q.run()
	q.join()
