import array
import sys
from Bio import SeqIO,Seq

def maskHomopolymer(arr, margin, k, pos, a, c, t, g):

    arrlen = len(arr)

    i=pos

    while i < arrlen:
        rightch = arr[i]
        if rightch == 'A': a += 1
        elif rightch == 'C': c += 1
        elif rightch == 'T': t += 1
        elif rightch == 'G': g += 1

        if a < margin and c < margin and t < margin and g < margin:
            break

        i += 1

        leftch = arr[i-k]
        if leftch == 'A': a -= 1
        elif leftch == 'C': c -= 1
        elif leftch == 'T': t -= 1
        elif leftch == 'G': g -= 1

    pos -= k

    while pos<i:
        arr[pos]='N'
        pos+=1

    return i

def detectAndMaskHomopolymer(dna, k=25, margin=.7):
    a = c = t = g = 0
    arr = array.array('c', dna.upper())

    margin = int(k*margin)

    i = 0
    winsz = 0
    arrlen = len(arr)

    while i < arrlen:
        if winsz >= k:
            if a > margin or c > margin or t > margin or g > margin:
                i = maskHomopolymer(arr, margin, k, i, a, c, t, g)
                winsz = a = c = t = g = 0
                if i >= arrlen: break

            leftch = arr[i-k]
            if leftch == 'A': a -= 1
            elif leftch == 'C': c -= 1
            elif leftch == 'T': t -= 1
            elif leftch == 'G': g -= 1

        rightch = arr[i]
        if rightch == 'A': a += 1
        elif rightch == 'C': c += 1
        elif rightch == 'T': t += 1
        elif rightch == 'G': g += 1

        i += 1
        winsz += 1

    return arr.tostring()

if len(sys.argv)!=2:
    print "usage:",sys.argv[0],"input.fasta/fastq"
    sys.exit(-1)

inp=sys.argv[1]
ext=inp.lower().split('.')[-1]

if ext in ['fa','fasta']:
    fmt='fasta'
elif ext in ['fq','fastq']:
    fmt='fastq'
else:
    print "Format not recognized"
    sys.exit(-1)

inp=SeqIO.parse(inp,fmt)

for i in inp:
    dna = detectAndMaskHomopolymer(str(i.seq))
    i.seq=Seq.Seq(dna)
    print i.format(fmt),
