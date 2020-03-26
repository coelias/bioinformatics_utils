from gzip import GzipFile
import sys
import itertools
import StringIO
import re

class FqDeshuffler:
	def __init__(self,inp,fill=True):
		# If fill=True, reads that are single reads will be mated with a new read containing an N
		self.fqr=FastQReader(inp)

	def __iter__(self):
		reads={}
		for i in self.fqr:
			name=i[0]
			l=reads.setdefault(name,[])
			l.append(i)
			if len(l)==2:
				del reads[name]
				yield l

		for name,i in reads.items():
			if len(i)==1:
				r1=i[0]
				r2=list(i[0][:])
				r2[2]='N'
				r2[4]='A'
				yield r1,r2

class FastQReader:
	'''FastQ file parser
	   detects if it is gzip or not
	   returns for every record a tuple(redname(without read number),first line,dna,thrdline,qual)'''

	# FASTA FORMATs
	# 1: @GAII01:1:1:0:1704#0/1
	# 2: @M01746:45:000000000-ADAHU:1:1101:13881:2101 1:N:0:10
	# 3: @SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36
	# 4: @ERR532261.1.1 HWI-ST858_56:5:1101:1304:1950 length=100




	def __init__(self,ifile):
		self.name=ifile
		if ifile=="-":
			self.inp=sys.stdin
		elif type(ifile)==file:
			self.inp=ifile
		else:
			try: self.inp=GzipFile(ifile)
			except: self.inp=open(ifile,"r")

		self.inibuffer=StringIO.StringIO()
		k=100
		for i in self.inp:
			self.inibuffer.write(i)
			k-=1
			if not k: break

		self.inibuffer.seek(0)
		first=self.inibuffer.readline()
		self.CRlen=len(first)-len(first.strip())
		self.inibuffer.seek(0)

		self.detectFormat()

		#extravalues={}

		#for line in self.inibuffer:
		#	if line[0]=='@':
		#		line=line[:-self.CRlen]
		#		if line[-2:-1]=='/':
		#			readnumberlen=len(line)-len(line.split("/")[0])
		#		else:
		#			readnumberlen=len(line)-len(line.split(" ")[0])
		#		extravalues.setdefault(readnumberlen,0)
		#		extravalues[readnumberlen]+=1

		#if extravalues:
		#	extravalues=[i[::-1] for i in extravalues.items()]
		#	extravalues.sort(reverse=True)
		#	self.readnumberlength=extravalues[0][1]
		#else:
		#	self.readnumberlength=0
		#self.inibuffer.seek(0)

	def detectFormat(self):
		bf=[i.strip() for i in self.inibuffer]
		self.inibuffer.seek(0)

		header=bf[16]   #picking 3rd header to avoid 1st and 2nd
		if header[-2]=="/" and header[-1] in "12":
			self.fastqformat=1
		elif re.findall("[12]:[YN]:[^:]+:[^: ]+$",header):
			self.fastqformat=2
#		elif re.findall("[0-9]+:[12]:[0-9]+:[0-9]+ length=[0-9]+$"): 
#			self.fastqformat=3
		else:
			header=header.split(" ")[0]
			if header[-2]=="." and header[-1] in "1234": self.fastqformat=4
			else:
				raise Exception("Cannot detect fastq format")


	def seekInit(self,*args):
		self.inp.seek(0)
		self.inibuffer.seek(0)

	def seek(self,pos):
		self.inp.seek(pos)

	def close(self):
		self.inp.close()

	def __iter__(self):
		NAME=0
		READ=1
		NAME2=2
		state=None

		entry=tuple()

		for i in itertools.chain(self.inibuffer,self.inp):
			i=i[:-self.CRlen]
			if not i:
				entry=tuple()
				state=None
				continue
			if state==None:
				if i[0]=="@":
					if self.fastqformat==1: rname=i[:-2]
					elif self.fastqformat==2: rname=i.split(" ")[0]
					elif self.fastqformat==4: rname=i.split(" ")[0][:-2]
					
					entry+=(rname[1:],i,)
					state=NAME
				else: 
					entry=tuple()
					state=None
			elif state==NAME:
				state=READ
				entry+=(i,)
			elif state==READ :
				if i[0]=="+":
					state=NAME2
					entry+=(i,)
				else:
					entry=tuple()
					state=None
			elif state==NAME2:
				if len(entry[2])==len(i):
					entry+=(i,)
					yield entry
				entry=tuple()
				state=None
	

class FastQWriter:
	'''Every record you write expects the same format as FastQReader returns'''

	BUFFLEN=1024*1
	def __init__(self,ofile):
		self.name=ofile
		if ofile=="-":
			self.output=sys.stdout
		else:
			if type(ofile)==file:
				self.output=ofile
			else:
				self.output=open(ofile,"w")
		self.closed=False
		self.buffer=""

	def close(self):
		if not self.closed:
			if self.buffer: self.output.write(self.buffer)
			self.buffer=""
			self.output.close()
			self.closed=True

	def writeRecord(self,i):
		self.buffer+=i[1]+"\n"+i[2]+"\n"+i[3]+"\n"+i[4]+"\n"
		# Next line is a faster way of doing len(), oh tricky carlos...	(old line: if len(self.buffer)>=FastQWriter.BUFFLEN:)
		if self.buffer[FastQWriter.BUFFLEN:FastQWriter.BUFFLEN+1]:
			self.output.write(self.buffer[:FastQWriter.BUFFLEN])
			self.buffer=self.buffer[FastQWriter.BUFFLEN:]
	write=writeRecord

	def writeGen(self,gen):
		for i in gen:
			self.buffer+=i[1]+"\n"+i[2]+"\n"+i[3]+"\n"+i[4]+"\n"
			if self.buffer[FastQWriter.BUFFLEN:FastQWriter.BUFFLEN+1]:
				self.output.write(self.buffer[:FastQWriter.BUFFLEN])
				self.buffer=self.buffer[FastQWriter.BUFFLEN:]

	def writeRecord1(self,i):
		if i[0][-2:] not in ["/1","/2"]:
			self.buffer+=i[1]+"/2\n"+i[2]+"\n"+i[3]+"\n"+i[4]+"\n"
		else:
			self.buffer+=i[0]+"\n"+i[1]+"\n"+i[2]+"\n"+i[3]+"\n"
		if self.buffer[FastQWriter.BUFFLEN:FastQWriter.BUFFLEN+1]:
			self.output.write(self.buffer[:FastQWriter.BUFFLEN])
			self.buffer=self.buffer[FastQWriter.BUFFLEN:]

	def writeRecord2(self,i):
		if i[0][-2:] not in ["/1","/2"]:
			self.buffer+=i[1]+"/2\n"+i[2]+"\n"+i[3]+"\n"+i[4]+"\n"
		else:
			self.buffer+=i[1]+"\n"+i[2]+"\n"+i[3]+"\n"+i[4]+"\n"
		if self.buffer[FastQWriter.BUFFLEN:FastQWriter.BUFFLEN+1]:
			self.output.write(self.buffer[:FastQWriter.BUFFLEN])
			self.buffer=self.buffer[FastQWriter.BUFFLEN:]

	def writeStr(self,s):
		self.output.write(s)

if __name__=='__main__':
	import sys
	FastQReader(sys.argv[1])
