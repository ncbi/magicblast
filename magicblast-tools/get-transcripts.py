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

"""Get transcript sequences from a genome and a GTF file"""

import gtf
import gff
from pyfaidx import Fasta
import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Get transcript sequences from a genome and a GTF file')
    parser.add_argument('--genome', metavar='FILE', dest='genome', type=str,
                        help='Reference sequence in FASTA format')
    parser.add_argument('--gff', metavar='FILE', dest='gff', type=str,
                        help='GFF or GTF file')
    parser.add_argument('--select', metavar='STRING', dest='select', type=str,
                        help='Print only sequences whose id contain provided '
                        'string')

    args = parser.parse_args()

    f = open(args.gff)
    if args.gff.endswith('.gtf'):
        transcripts = gtf.get_transcripts(f)
    elif args.gff.endswith('.gff'):
        transcripts = gff.get_mrnas(f)
    else:
        raise ValueError('Unrecognized file extension for: {}. Only GFF or GTF'
                         ' files are allowed'.format(args.gff))
    f.close()

    genome = Fasta(args.genome)

    for i in transcripts:

        strand = transcripts[i].exons[0].strand
        sequence = ''

        exons = sorted(transcripts[i].exons, key=lambda x: x.start,
                       reverse = (strand == '-'))

        seqid = transcripts[i].seqid

        for exon in exons:

            if exon.strand != strand:
                raise ValueError('Mismatched strands for transcript: {0}'.\
                                 format(i))

            if strand == '-':
                sequence += genome[seqid][(exon.start - 1):(exon.end)].\
                            reverse.complement.seq.upper()
            else:
                sequence += genome[seqid][(exon.start - 1):(exon.end)].\
                            seq.upper()

        seqid = i
        if 'Name' in transcripts[i].attributes:
            # A few mRNAs align to both X and Y chromosomes in slightly
            # different locations, so we are adding reference id to sequence
            # id to distinguish between the two alignments
            seqid = transcripts[i].attributes['Name'] + ':' + transcripts[i].seqid

        if args.select and args.select not in seqid:
            continue

        print('>{0}'.format(seqid))
        for n in range(0, len(sequence), 80):
            print('{0}'.format(sequence[n:(n + 80)]))

            
