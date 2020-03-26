#!/usr/bin/python

from fqreader import *
import sys
import os
import argparse
import uuid

class FqMerge:
	class GeneratorWrapper:
		def __init__(self,fqr):
			self.gen=fqr.__iter__()
			self.nx=self.gen.next()
			self.name=self.nx[1]
		def pop(self):
			ret=self.nx
			try: 
				self.nx=self.gen.next()
				self.name=self.nx[1]
			except: self.nx=None
			return ret

	def __init__(self,nameList):
		self.fils=[FastQReader(i) for i in nameList]
		self.sortedlist=[FqMerge.GeneratorWrapper(i) for i in self.fils]
		self.sortedlist.sort(key=lambda x: x.name)

	def __iter__(self):
		while self.sortedlist:
			el=self.sortedlist[0].pop()
			if not el: 
				self.sortedlist.pop(0)
			else:
				self.sortedlist.sort(key=lambda x: x.name)
				yield el

	def close(self):
		for i in self.fils:
			i.close()
		

def fqSort(fqfil,out):
	READSXFILE=200000
	fils=[]
	nfiles=0

	if fqfil.name=="-":
		sourcename=str(uuid.uuid4())
	else:
		sourcename=fqfil.name

	rds=fqfil.__iter__()


	while True:
		reads=[]
		k=0
		for i in rds:
			reads.append(i)
			k+=1
			if k>=READSXFILE:
				break
		
		reads.sort(key=lambda x: x[1])
		nf=FastQWriter(sourcename+".sort."+str(nfiles)+".tmp")
		nfiles+=1
		fils.append(nf)
		for i in reads:
			nf.writeRecord(i)
		if len(reads)!=READSXFILE:
			break

	fqfil.close()
	for i in fils:
		i.close()
	fm=FqMerge([i.name for i in fils])
	out=FastQWriter(out)
	for i in fm:
		out.writeRecord(i)

	out.close()
	fm.close()
	for i in fils:
		os.remove(i.name)


if __name__=="__main__":
	parser = argparse.ArgumentParser(description='Sort FastqFile')
	parser.add_argument("-i", dest="fin",help="Input file1 (default stdin)",default="-")
	parser.add_argument("-o", dest="fout",help="Output file (default stdout)",default="-")
	args = parser.parse_args()

	fqSort(FastQReader(args.fin),args.fout)
