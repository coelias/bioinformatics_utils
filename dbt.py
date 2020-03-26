import threading
import Queue
import subprocess
import shlex
import argparse
import uuid
import os
import sys

COMM=Queue.Queue()

def worker(cmd,outf):
    out=open(outf,'w')
    cmd=subprocess.Popen(shlex.split(cmd),stdin=subprocess.PIPE,stdout=out)
    while True:
        inp=COMM.get()
        if inp=='DIE': break
        cmd.stdin.write(inp)
        cmd.stdin.flush()
    cmd.stdin.close()
    out.flush()
    out.close()

def divide_stream(st,bs=1):
    bl=[]
    for i in st:
        bl.append(i)
        if len(bl)==bs:
            yield ''.join(bl)
            bl=[]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='FastaM tools')

    parser.add_argument(
        "-s",
        dest="inpsize",
        help="Input size (lines) [1]",
        required=False,
        default=1,
        type=int)
    parser.add_argument(
        "-t",
        dest="threads",
        help="Threads [5]",
        default=5,
        type=int)
    parser.add_argument(
        "-c",
        dest="cmd",
        help="Command to execute",
        required=True)
    parser.add_argument(
        "-o",
        dest="output",
        help="outputfile [stdout]",
        default='-')

    parser.add_argument("input", nargs='*', default=["-"])
    options = parser.parse_args()

    workers=[]
    for i in range(options.threads):
        outf=".deleteme."+str(uuid.uuid4())[:8]
        w=threading.Thread(target=worker,kwargs={"cmd":options.cmd,"outf":outf})
        w.start()
        workers.append([w,outf])

    for i in options.input:
        if i == '-':
            inp = sys.stdin
        else:
            inp = open(i)

        for bl in divide_stream(inp, options.inpsize):
            COMM.put(bl)

    for i in range(options.threads):
        COMM.put('DIE')

    for w,_ in workers:
        w.join()

    if options.output=='-':
        output_file=sys.stdout
    else:
        output_file=open(options.output,'w')

    for _,outp in workers:
        with open(outp) as q:
            for i in q:
                output_file.write(i)
        os.unlink(outp)

