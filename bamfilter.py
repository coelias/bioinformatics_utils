#!/usr/bin/python
import argparse
from pysam import *
import sys
import os
from math import *


class Filter:
    def __init__(self):
        pass

    def filter(self, r):
        pass


class Dedup(Filter):
    def __init__(self, dic):
        self.dedup = dic

    def filter(self,al):
        if al.flag & 2304:
            return True
        if al.qname in self.dedup:
            if self.dedup[al.qname]:
                for k, l in al.tags:
                    if k == 'AS':
                        if l == dedup[al.qname]:
                            dedup[al.qname] = None
                            return False
            return True

        return False


class Quality(Filter):
    def __init__(self, q=None):
        self.q = q

    def filter(self,al):
        for k, l in al.tags:
            if k == 'AS':
                if l>=self.q: return False
                return True

        return False

class QLen(Filter):
    def __init__(self, q=None):
        self.q = q

    def filter(self,al):
        if al.qlen<self.q: return True
        return False

class Mapped(Filter):
    def __init__(self, q=None):
        self.q = q

    def filter(self,al):
        if al.is_unmapped: return True
        return False

class Unmapped(Filter):
    def __init__(self, q=None):
        self.q = q

    def filter(self,al):
        if al.is_unmapped: return False
        return True

class TrimQual(Filter):
    def filter(self,al):
        al.qual=""
        return False

class QCScore(Filter):
    def __init__(self, q=None):
        self.q = q

    def mean_qscore(self,scores):
        """ Returns the phred score corresponding to the mean of the probabilities
        associated with the phred scores provided.
        :param scores: Iterable of phred scores.
        :returns: Phred score corresponding to the average error rate, as
        estimated from the input phred scores.
        """
        if len(scores) == 0:
            return 0.0
        sum_prob = 0.0
        for val in scores:
            sum_prob += pow(10, -0.1 * val)
        mean_prob = sum_prob / len(scores)
        return -10.0 * log10(mean_prob)

    def qstring_to_phred(self,quality):
        """ Compute standard phred scores from a quality string. """
        qscores = [ord(q) - 33 for q in quality]
        return qscores

    def filter(self,al):
        if self.mean_qscore(self.qstring_to_phred(al.qual))<self.q: return True
        return False

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Drops reads from bam files')
    parser.add_argument(
        "-r",
        dest="replace",
        help="replace original file",
        action='store_true',
        default=False)

    parser.add_argument(
        "-t",
        dest="trimqual",
        help="Trim qualities",
        action='store_true',
        default=False)

    parser.add_argument(
        "-d",
        dest="altmappings",
        help="Removes alternative mappings (only keeps max score)",
        action='store_true',
        default=False)

    parser.add_argument(
        "-a",
        dest="mapqual",
        help="Filters by minimum alignment score",
        default=0,
        type=int)

    parser.add_argument(
        "-m",
        dest="mapped",
        help="Filters unmapped reads",
        action='store_true',
        default=False)

    parser.add_argument(
        "-u",
        dest="unmapped",
        help="Filters mapped reads, leaving unmapped reads",
        action='store_true',
        default=False)

    parser.add_argument(
        "-q",
        dest="qcscore",
        help="Filters by minimum readqcscore",
        default=0,
        type=int)

    parser.add_argument(dest="input", help="SAM/BAM file")
    args = parser.parse_args()

    filters=[]

    if args.altmappings:
        dedup = {}
        with Samfile(args.input) as s:
            for al in s:
                score = None
                for i, j in al.tags:
                    if i == 'AS':
                        score = j
                        break
                if not score: continue
                dedup.setdefault(al.qname, []).append(score)

        dedup = dict([(i, max(j)) for i, j in dedup.iteritems() if len(j) > 1])

        filters.append(Dedup(dedup))

    if args.mapqual:
        filters.append(Quality(args.mapqual))

    if args.qcscore:
        filters.append(QCScore(args.qcscore))

    if args.mapped:
        filters.append(Mapped())

    if args.unmapped:
        filters.append(Unmapped())

    if args.trimqual:
        filters.append(TrimQual())

    src = Samfile(args.input)
    dst = Samfile(args.input[:-4]+".filter.bam", "wb", header=src.header)

    for i in src:
        for f in filters:
            if all([not f.filter(i) for f in filters]):
                dst.write(i)

    dst.close()

    if args.replace:
        import shutil
        shutil.move(args.input[:-4]+'.filter.bam', args.input)
