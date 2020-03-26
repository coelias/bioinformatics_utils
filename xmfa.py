#!/usr/bin/python


from Bio import SeqIO
import sys
import re
import argparse

rename=re.compile("^([^:]+):([0-9]+)-([0-9]+)( [+-])?")

class XMFA:
    class Seq:
	def __init__(self,name,start,end,seq,rc=False):
	    self.name=name
	    self.start=start
	    self.end=end
	    self.seq=str(seq)

	def fastaid(self):
	    return "{0}:{1}-{2}".format(self.name,self.start,self.end)
	
	def __str__(self):
	    return "{0} ({4}) , ({1}:{2}), [{3}...]".format(self.fastaid(),self.start,self.end,self.seq[:50],len(self.seq))

	def __len__(self):
	    return len(self.seq)
	
	def toString(self):
	    return ">{0}\n{1}".format(self.fastaid(),self.seq)


    class Alignment:
	def __init__(self,seqs):
	    assert len(set([len(i) for i in seqs]))==1
	    self.length=len(seqs[0])
	    self.seqs=dict([(i.name,i) for i in seqs])

	def __iter__(self):
	    q=sorted(self.seqs.items())
	    q=[i[1] for i in q]
	    for i in itertools.izip(q):
		yield i

	def sequences(self):
	    return self.seqs.values()

	def __len__(self):
		return self.length

	def toString(self):
	    return "\n".join([i.toString() for i in self.seqs.values()])


	def __str__(self):
	    return "Alignment of seqs {0} ({1} bp long)".format([i.fastaid() for i in self.seqs.values()],self.length)

    def __init__(self,filename):
	self.alignments=[]
	self.uniqueSequences={}

	cur=[]
	
	for i in SeqIO.parse(filename,"fasta"):
		try: sqid,start,end,rc=rename.findall(i.name)[0]
		except: raise Exception ("Wrong XFMA name: {0}".format(repr(i.name)))
		if not rc: rc="+"
		rc=rc.strip()

		if i.seq[-1]=="=":
			cur.append(XMFA.Seq(sqid,int(start),int(end),str(i.seq)[:-1],rc!='-'))
			if len(cur)>1:
				self.alignments.append(XMFA.Alignment(cur))
			else:
				self.uniqueSequences.setdefault(cur[0].name,[]).append(cur[0])
			cur=[]
		else:
			cur.append(XMFA.Seq(sqid,int(start),int(end),str(i.seq),rc!='-'))

	self.alignments.sort(key=lambda x: len(x),reverse=True)

	self.fastaid2object={}
	for i,j in self.uniqueSequences.items():
	    for k in j:
		self.fastaid2object[k.fastaid()]=k

	for i in self.alignments:
	    for s in i.sequences():
		self.fastaid2object[s.fastaid()]=i

    def __getitem__(self,key):
	return self.fastaid2object[key]

    def trim(self,lentgh):
	self.alignments=[i for i in self.alignments if len(i)>lentgh]
	for i in self.uniqueSequences.keys():
	    self.uniqueSequences[i]=[j for j in self.uniqueSequences[i] if len(j)>lentgh]

    def __str__(self):
	cad="\n".join([str(i) for i in self.alignments])
	alnlen=sum([len(i) for i in self.alignments])
	cad+="\n\nAlignment size={0}\n".format(alnlen)

	cad+="\nUnique Sequences:\n"

	for i,j in self.uniqueSequences.items():
		j.sort(key=lambda x: len(x),reverse=True)
		for seq in j:
			cad+="{0} : {1}\n".format(i,seq)
		alnlen=sum([len(k) for k in j])
		cad+="Total unique bs for {1}: {0}\n\n".format(alnlen,i)



	return cad

if __name__=="__main__":
	parser = argparse.ArgumentParser(description='XFMA extraction tool')
	parser.add_argument("-i",dest="info",help="Shows general information",required=False,default=False,action="store_true")
	parser.add_argument("-t",dest="tlength",help="Ignores everything smaller than tlength bp",required=False,default=0,type=int)
	parser.add_argument("-u",dest="uniqseqs",help="Dumps unique sequences for a given id",required=False,default=None)
	parser.add_argument("-e",dest="fastaid",help="Extracts Unique-sequence or Alignment given a seq ID",required=False,default=None,nargs="+")
	parser.add_argument('FILE')
	args = parser.parse_args()

	a=XMFA(args.FILE)
	if args.tlength:
	    a.trim(args.tlength)
	if args.info:
	    print a
	if args.uniqseqs:
	    for i in a.uniqueSequences[args.uniqseqs]:
		print i.toString()
	if args.fastaid:
	    for i in args.fastaid:
		print a[i].toString()
