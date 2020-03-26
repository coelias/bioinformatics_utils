import sys
import re
import uuid
import os
from random import randint

def fqparser(fqfile):
    with open(fqfile) as fqf:
        lst = []
        for i in fqf:
            if i[0] != '@' and not lst:
                continue
            lst.append(i.strip())
            if len(lst) == 4:
                if lst[0][0] == '@' and lst[2][0] == '+':
                    yield lst
                lst = []

reuuid=re.compile('[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}')
validrecord=re.compile('^@[0-9a-z]{8}-[0-9a-z]{4}-[0-9a-z]{4}-[0-9a-z]{4}-[0-9a-z]{12}(?:\s+(?:\S+)=(?:\S+))+')


try:
    assert os.path.isfile(sys.argv[1])
    blocksize=int(sys.argv[2])
    fname=os.path.basename(sys.argv[1][:-2])
except:
    print "Usage: {} FILE.fq blocksize".format(sys.argv[0])
    sys.exit(-1)

buff=[]
nblock=1
knt=0

for record in fqparser(sys.argv[1]):
    if not validrecord.match(record[0]):
        uuid_list=reuuid.findall(record[0])
        if not uuid_list:
            guid=uuid.uuid4()
        else:
            guid=uuid_list[0]

        knt+=1
        record[0]="@{} runid=789979b907d213a575d472468164d65bf02669e8 read={} ch={} start_time=2017-05-12T16:11:59Z\n".format(guid,knt,randint(0,511))

    buff.append(record)

    if len(buff)==blocksize:
        z=open("{}_{}.fastq".format(sys.argv[1].rsplit('.',1)[0],nblock),'w')
        nblock+=1
        z.write('\n'.join(['\n'.join(i) for i in buff])+'\n')
        z.close()
        buff=[]

if buff:
    z=open("fastq_{}_runid_789979b907d213a575c472468164d65bf02669e8.fastq".format(nblock),'w')
    z.write(''.join([''.join(i) for i in buff]))
    z.close()
