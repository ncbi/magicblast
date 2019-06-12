#!/usr/bin/env python3
#============================================================================
#
#                           PUBLIC DOMAIN NOTICE
#               National Center for Biotechnology Information
#
#  This software/database is a "United States Government Work" under the
#  terms of the United States Copyright Act.  It was written as part of
#  the author's official duties as a United States Government employee and
#  thus cannot be copyrighted.  This software/database is freely available
#  to the public for use. The National Library of Medicine and the U.S.
#  Government have not placed any restriction on its use or reproduction.
#
#  Although all reasonable efforts have been taken to ensure the accuracy
#  and reliability of the software and data, the NLM and the U.S.
#  Government do not and cannot warrant the performance or results that
#  may be obtained by using this software or data. The NLM and the U.S.
#  Government disclaim all warranties, express or implied, including
#  warranties of performance, merchantability or fitness for any particular
#  purpose.
#
#  Please cite the author in any work or product based on this material.
#
# ===========================================================================
#
# Author: Greg Boratyn boratyng@ncbi.nlm.nih.gov
#
# ---------------------------------------------------------------------------

"""Parse a SAM/BAM file and get an intron support table"""

import sam
import gff
import gtf
import txt
import gzip
import pickle
import argparse
import sys
import numpy as np
import pandas as pd
from pyfaidx import Fasta
from collections import Counter
from collections import defaultdict

def get_splice_signals(introns, fasta_filename):
    """Find splice signals for introns and return as dictionary indexed by introns"""
    sites = {}
    genome = Fasta(fasta_filename)
    for i in introns:
        if i.seqid not in genome:
            sites[i] = 'xxxx'
        else:
            # indices into pyfaidx sequences are zero-based
            sites[i] = genome[i.seqid][(i.start - 1):(i.start + 1)].seq.upper() + genome[i.seqid][(i.end - 2):i.end].seq.upper()

    return sites


def print_splice_signal_histogram(introns, signals):
    """Find and print histogram of splice signals for a given set of introns"""
    hist = Counter()
    for i in introns:
        hist[signals[i]] += 1
    for s, c in sorted([(s, hist[s]) for s in hist], key=lambda x: x[1], reverse=True):
        print('{0}\t{1}'.format(s, c))
    print('')


class EmptySpliceSignal:
    """No splice signal, class needed for defaultdict"""
    def __call__(self):
        """Return a constant representing empty splice signal"""
        return '---'


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Generate intron support table from a SAM/BAM file')
    parser.add_argument('--bam', metavar='FILE', dest='bamfile', type=str,
                        help='BAM file')
    parser.add_argument('--gff', metavar='FILE', dest='gfffile', type=str,
                        help='GFF file')
    parser.add_argument('--genome', metavar='FILE', dest='genome', type=str,
                        help='Genome FASTA file')
    parser.add_argument('--filter-by', metavar='LIST', dest='filter_by',
                        type=str, help='Filter reads')
    parser.add_argument('--filter-annot', metavar='STRING',
                        dest='filter_annot', type=str,
                        help='Filter annotated introns by transcript accession')
    parser.add_argument('--introns', metavar='FILE', dest='introns', type=str,
                        help='Output file for the intron table'
                        '(default: stdout)', default='-')
    parser.add_argument('--spec', dest='spec', action='store_true',
                        help='Show results for sensitivity and specificity'
                        ' analysis')
    parser.add_argument('--splice_histogram', dest='splice_histogram',
                        action='store_true',
                        help='Show splice signal histogram')
    parser.add_argument('--numbers', dest='numbers', action='store_true',
                        help='Show numbers of annotated and unannotated '
                        'introns')
    parser.add_argument('--sort', metavar='FILE', dest='sort', type=str,
                        help='Output sort file for ROC score computation')
    parser.add_argument('--counts', metavar='FILE', dest='counts', type=str,
                        help='Tab delimited file with numbers of read '
                        'placements for weighting read contrinution.'
                        'Format: read_id, mate_id (1, 2), number of placements')
    parser.add_argument('--max-count', metavar='FILE', dest='max_count',
                        type=int, help='Maximum count for computation of '
                        'weighted read alignment counts')
                        

    args = parser.parse_args()

#    if not args.gfffile:
#        raise InputError('Annotation file not specified, use --gff option')

    gff_introns = {}
    if args.gfffile:
        print('Reading annotatiotns', file=sys.stderr)

        if args.gfffile.endswith('pickle'):
            f = open(args.gfffile, 'rb')
            gff_introns = pickle.load(f)
            f.close()
        else:
            if args.gfffile.endswith('.gz'):
                f = gzip.GzipFile(args.gfffile, 'r')
            else:
                f = open(args.gfffile)

            if args.gfffile.endswith('.gff') or args.gfffile.endswith('.gff.gz'):
                gff_introns = gff.get_splice_sites(f,
                                                   accession=args.filter_annot)
            elif args.gfffile.endswith('.gtf') or args.gfffile.endswith('.gtf.gz'):
                gff_introns = gtf.get_splice_sites(f)
            elif args.gfffile.endswith('.sam') or args.gfffile.endswith('.sam.gz'):
                gff_introns = sam.get_introns_with_reads(args.gfffile)
            elif args.gfffile.endswith('.txt'):
                gff_introns = txt.get_introns(f)
            else:
                raise InputError('Unrecognized annotation file extension, '\
                                 'must be one of these: .gff, .gtf, .sam, '\
                                 '.pickle')
            f.close()
        print('{0} introns in the annotation'.format(len(gff_introns)))
        print('done', file=sys.stderr)
    else:
        gff_introns = {}

    placements = None
    if args.counts:
        max_count = None
        if args.max_count:
            max_count = args.max_count
        placements = sam.read_placements(args.counts, max_count)

    # parse command line filtering arguments
    filter_by = None
    if args.filter_by:
        filter_by = {}
        a = args.filter_by.split()
        for w, i in zip(a[::2], a[1::2]):
            if w == 'read_id':
                f = open(i)
                filter_by[w] = {i: 1 for i in f.read().splitlines()}
                f.close()
            else:
                filter_by[w] = int(i)
            print('   filtering reads by {0}\t{1}'.format(w, i), file=sys.stderr)


    # get introns from a SAM/BAM file
    print('Reading SAM/BAM file', file=sys.stderr)
    introns = sam.get_introns_with_reads(args.bamfile, force_single = True,
                                         filter_by = filter_by,
                                         placements = placements)
    print('done', file=sys.stderr)

    splice_signals = defaultdict(EmptySpliceSignal())
    if args.genome:
        splice_signals = get_splice_signals(introns, args.genome)

        if args.splice_histogram:
            print_splice_signal_histogram(introns, splice_signals)


    # sort introns by end position, start position, reference sequence id
    keys = sorted(introns.keys(), key=lambda x: x.end)
    keys.sort(key=lambda x: x.start)
    keys.sort(key=lambda x: x.seqid)

    # print intron coverage table
    f = sys.stdout
    if args.introns != '-':
        f = open(args.introns, 'w')
    for k in keys:
        known = 'NEW'
        if k in gff_introns:
            known = 'KNOWN'
        print('{0}\t{1}\t{2}\t{3}\t{4}\t{5}'\
              .format(k.seqid, k.start, k.end, introns[k], known,
                      splice_signals[k]), file=f)
    if args.introns != '-':
        f.close()
    

    # create a dataframe to compute sensitivity, precision, and a sort file
    if args.spec or args.sort or args.numbers:

        d = {'intron': [':'.join(str(i) for i in [k.seqid, k.start, k.end]) \
                        for k in keys],
             'coverage': [introns[k] for k in keys],
             'annotation': ['KNOWN' if k in gff_introns else 'NEW' \
                            for k in keys]}

        df = pd.DataFrame(data=d)

        num_annot = df['annotation'].apply(func=lambda x: 0 if x == 'KNOWN' \
                                           else 1)
        df['key'] = 2 * df['coverage'] + num_annot
        df.sort_values(by='key', ascending=False, inplace=True)
#        df = df.reindex(columns=['intron', 'coverage', 'annotation', 'key'])

    # print sensitivity/specificity report
    if args.spec:
        num_known = len(gff_introns)
        print('Number of known introns: {0}'.format(num_known))
        for i in [1, 2, 3, 5]:
            print('Coverage >= {0}'.format(i))
            idx = (df['coverage'] >= i)
            num_tp = df[idx & (df['annotation'] == 'KNOWN')].shape[0]
            sensitivity = float(num_tp) / float(num_known)
            specificity = float(num_tp) / float(df[idx].shape[0])
            print('Sensitivity: {0}'.format(sensitivity))
            print('Precision: {0}'.format(specificity))
            print('')
                  
    # print number of annotated and unannotated introns
    if args.numbers:
        num_known = len(gff_introns)
        print('Number of known introns: {0}'.format(num_known))
        for i in [1, 2, 3, 5]:
            print('Coverage >= {0}'.format(i))
            idx = (df['coverage'] >= i)
            num_tp = df[idx & (df['annotation'] == 'KNOWN')].shape[0]
            num_fp = df[idx & (df['annotation'] != 'KNOWN')].shape[0]
            print('Number of annotated introns: {0}'.format(num_tp))
            print('Number of unannotated introns: {0}'.format(num_fp))
            print('')

    # print the sort file
    if args.sort:
        with open(args.sort, 'w') as f:
            n = 0
            for r in df.iterrows():
                print('\t'.join(str(i) for i in [n, r[1]['coverage'], \
                              '+' if r[1]['annotation'] == 'KNOWN' else '-', \
                              r[1]['intron']]), file=f)
                n += 1



