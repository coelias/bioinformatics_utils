#!/usr/bin/python
from StringIO import StringIO
import re
import urllib


def _spFromAccession():
    cache = {}

    def a(acc):
        if acc in cache:
            return cache[acc]

#        try:
        b = urllib.urlopen('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={}&rettype=fasta&retmode=xml'.format(acc))

        for i in b:
            if 'TSeq_taxid' in i:
                c = re.findall('[0-9]+', i)[0]
                b.close()
                return c
#        except:
#            b = urllib.urlopen('https://www.ebi.ac.uk/ena/data/view/{}&display=xml&download=xml&filename={}.xml'.format(acc))
#
#            for i in b:
#                if 'taxId=' in i:
#                    c = re.findall('([0-9]+)">', i)[0]
#                    b.close()
#                    return c
            
    return a

acc2seqid = _spFromAccession()

if __name__ == '__main__':
    import argparse
    from Queue import Queue
    from threading import Thread
    import sys

    THREADS = 2
    parser = argparse.ArgumentParser(description='NCBI Accessions to taxid')
    parser.add_argument('-i', dest="input",help="File input (one accession per line, default stdin)",type=str,default='-')
    parser.add_argument('-o', dest="outfile",help="file output (default stdout)",type=str,default='-')
    parser.add_argument('-n', dest="ncbidbs",help="Use only if massive number of accession numbers, will download ncbidbs",action='store_true')

    args = parser.parse_args()

    finput = sys.stdin if args.input == '-' else open(args.input)
    fout = sys.stdout if args.outfile == '-' else open(args.outfile,'w')

    accessions = set([i.strip() for i in finput])

    print len(accessions), 'total accession numbers'

    results = 0

    if args.ncbidbs:
        from gzip2 import gzfile
        resources = [
          'http://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz',
          'http://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz',
          'http://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gss.accession2taxid.gz',
          'http://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_est.accession2taxid.gz'
        ]

        for resource in resources:
            sys.stderr.write('\nDownloading {}\n'.format(resource))
            resource = urllib.urlopen(resource)
            size = int(re.findall('Content-Length: ([0-9]+)', str(resource.headers))[0])
            a = gzfile(fileobj=resource, sz=size)

            for i in a:
                i = i.split()
                found=None
                if i[1] in accessions: found=i[1]
                elif i[0] in accessions: found=i[0]

                if found:
                    results += 1
                    fout.write('{}\t{}\n'.format(found, i[2]))
                    tot = '  {}/{}'.format(results, len(accessions))
                    sys.stderr.write('{}\033[{}D'.format(tot, len(tot)))
                    sys.stderr.flush()
                    if results >= len(accessions):
                        break

            if results >= len(accessions):
                break
    else:
        inacc = Queue()
        outinfo = Queue()

        def worker():
            while True:
                acc = inacc.get()
                if acc == 'DIE': break
                try: outinfo.put("{}\t{}\n".format(acc, acc2seqid(acc)))
                except: outinfo.put("{}\tERR\n".format(acc))
            outinfo.put('DONE')

        ths = []

        for i in accessions:
            inacc.put(i)

        for i in range(THREADS):
            th = Thread(target=worker)
            th.start()
            ths.append(th)
            inacc.put('DIE')

        done = 0
        while True:
            if done == THREADS: break
            d = outinfo.get()
            if d == 'DONE': done += 1
            else:
                results += 1
                fout.write(d)
                sys.stderr.write('\r{}/{}'.format(results, len(accessions)))

        for i in ths:
            i.join()
