import pysam
import sys
import argparse
import random
import re


class Subsample:
	def __init__(self, inp, nsamp, mem):
            self.mem = mem
            self.inp = inp
            nsamp = nsamp.split(',')
            self.reads = []
            self.header = None

            self.percsamples=[]
            self.nsamples=[]

            for i in nsamp:
                if re.findall('^[0-9]+%$', i):
                    self.percsamples.append(int(i[:-1]))
                elif re.findall('^[0-9]+$', i):
                    self.nsamples.append([int(i),int(i)])

            self.nreads = self.countReads()

            for i in self.percsamples:
                self.nsamples.append([int(self.nreads*(float(i)/100)),i])

        def run(self):
            for n,name in self.nsamples:
                outf=self.inp[:-3]+str(name)+".bam"
                wr=pysam.Samfile(outf,'wb',header=self.header)
                if self.mem:
                    random.shuffle(self.reads)
                    while n:
                        wr.write(self.reads[n])
                        n-=1
                else:
                    posics=range(self.nreads)
                    random.shuffle(posics)
                    posics=set(posics[:n])
                    k=0
                    with pysam.Samfile(self.inp) as f:
                        for i in f:
                            if k in posics:
                                wr.write(i)
                            k+=1
                wr.close()


                

        def countReads(self):
            a = pysam.Samfile(self.inp)

            self.header = a.header
            if self.mem:
                self.reads = list(a)
                a.close()
                return len(self.reads)
            else:
                k = 0
                for i in a:
                    k += 1
                a.close()
                return k


            
		
		
	
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Subsample Bam files')
	parser.add_argument("-n", dest="nsamp",  help="Subsample number or percentage [, separated it several needed]. Examples: -n 10000,20000,90%,35%", default=0)
	parser.add_argument("-m", dest="mem",  help="Load bam into memory (faster if multiple subsamples needed.", action='store_true' ,default=False)
	parser.add_argument(dest="input", help="SAM/BAM file")
	args = parser.parse_args()

	sub = Subsample(args.input, args.nsamp, args.mem)
	sub.run()


