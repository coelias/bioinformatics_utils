#!/usr/bin/python
import sys
import re
from uuid import uuid4
import datetime
import random
import os

random_runid=''.join([random.choice("abcfed1234567890") for i in range(40)])

re_uuid = re.compile("[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}",re.I)
re_chandread = re.compile("_ch([0-9]+)_read([0-9]+)")
re_date = re.compile('_(201[0-9]{5})_')
HEADER_REGEX = re.compile(
        '(?:^[>@]([0-9a-z]{8}-[0-9a-z]{4}-[0-9a-z]{4}-[0-9a-z]{4}-[0-9a-z]{12})|(\S+)=(\S+))')
now=datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ")


def process(header,seq,*others):
#    if r[0][0] != '@' or r[2][0] != '+':
#        raise Exception('Invalid FastQ format')

    fields = dict([(j, k) for _, j, k in HEADER_REGEX.findall(header)[1:]])

    uuid = re_uuid.findall(header)
    if uuid:
        uuid = uuid[0]
    else:
        uuid = str(uuid4())

    chanread = re_chandread.findall(header)
    if chanread:
        ch, read = chanread[0]
    else:
        ch = read = 'NA'

    date = re_date.findall(header)
    if date:
        date = "{}T00:00:00Z".format(date[0])
    else:
        date = now

    fields.setdefault('ch', ch)
    fields.setdefault('read', read)
    fields.setdefault('start_date', date)
    fields.setdefault('runid', random_runid)

    header = "{} {}".format(uuid,' '.join(["{}={}".format(i,j) for i,j in fields.items()]))
    if not others:
        sys.stdout.write("@{}\n{}\n+\n{}\n".format(header,seq,'I'*len(seq)))
    else:
        sys.stdout.write("@{}\n{}\n+\n{}\n".format(header,seq,others[0]))

record = []

if len(sys.argv)!=2:
    print 'Usage:',os.path.basename(sys.argv[0]),'FILE.[fasta/fa/fq/fastq] --> (stdout) ONT compatible FastQ'
    sys.exit(1)

fmt=sys.argv[1].lower().rsplit('.',1)[-1]

if fmt in ['fa','fasta']:
    rlen=2
if fmt in ['fq','fastq']:
    rlen=4


for i in open(sys.argv[1]):
    record.append(i.strip())
    if len(record) == rlen:
        if rlen==2:
            process(record[0],record[1])
        if rlen==4:
            process(record[0],record[1],record[3])
        record=[]
