#!/usr/bin/python

import argparse
import subprocess
import uuid
import os
import re
import shutil
import sys

def findinlist(l,cad):
    pos=0
    for i in l:
        if re.findall(cad,i): 
            return pos
        pos+=1
    return -1

class Mummer():
    def __init__(self,r,q,rn,qn,k=31):
        self.k=k
        self.ref=r
        self.que=q
        self.ref_n=rn
        self.que_n=qn
        self.devnull=open('/dev/null','w')
        self.prefix=str(uuid.uuid4())[:8]
    
        if not self.ref_n:
            self.ref_n=re.sub("\.[^.]+$","",os.path.basename(self.ref))
        if not self.que_n:
            self.que_n=re.sub("\.[^.]+$","",os.path.basename(self.que))
    
    def do_mummer(self):
        with open(self.prefix+'.mum','w') as mumout:
            mum=subprocess.Popen(['mummer','-mum','-n','-l',str(self.k),'-b','-c',self.ref,self.que],stderr=self.devnull,stdout=mumout)
            mum.wait()
    
    def do_plot(self,outpng='out.png',outmum=None):
        mumplt=subprocess.Popen(['mummerplot','--large','-t','png','-p',self.prefix,self.prefix+'.mum'],stderr=self.devnull,stdout=self.devnull)
        mumplt.wait()
        with open(self.prefix+'.gp') as a:
            gnuplot_inp=[i for i in a if 'set mouse' not in i]

        gnuplot_inp[0]="set terminal png large size 2000,2000 font \"Helvetica,14\"\n"
        xlabel=findinlist(gnuplot_inp,'^set xlabel')
        ylabel=findinlist(gnuplot_inp,'^set ylabel')

        gnuplot_inp[xlabel]='set xlabel "{}"\n'.format(self.ref_n)
        gnuplot_inp[ylabel]='set ylabel "{}"\n'.format(self.que_n)

        gnuplt=subprocess.Popen(['gnuplot'],stderr=self.devnull,stdout=self.devnull,stdin=subprocess.PIPE)
        gnuplt.stdin.write(''.join(gnuplot_inp))
        gnuplt.stdin.close()
        gnuplt.wait()

        if outmum:
            shutil.move(self.prefix+'.mum',outmum)
        else:
            os.unlink(self.prefix+'.mum')

        shutil.move(self.prefix+'.png',outpng)

        try:
            os.unlink(self.prefix+'.gp')
            os.unlink(self.prefix+'.fplot')
            os.unlink(self.prefix+'.rplot')
        except:
            pass


        
    
    
if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Mummer wrapper')

    parser.add_argument("-r", dest="ref", help="ref fasta", required=True)
    parser.add_argument("-q", dest="query", help="query fasta", required=True)

    parser.add_argument("-rn", dest="refname", help="ref name", required=False)
    parser.add_argument("-qn", dest="queryname", help="query name", required=False)

    parser.add_argument("-mo", dest="outmum", help="Mumer output", required=False)
    parser.add_argument("-o", dest="outpng", help="PNG output [out.png]", required=False,default='out.png')

    parser.add_argument("-k", dest="k", help="Kmer size", required=False, default=31, type=int)

    options = parser.parse_args()

    m=Mummer(options.ref,options.query,options.refname,options.queryname,k=options.k)

    m.do_mummer()
    m.do_plot(outpng=options.outpng,outmum=options.outmum)


