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

"""GTF file parser"""

from collections import namedtuple
from base import Intron
from base import Exon
from base import mRNA
import sys
import re


Record = namedtuple('Record', ['seqid', 'source', 'feature', 'start',
                               'end', 'score', 'strand', 'frame',
                               'attribute'])

def parse(line):
    """Parse a single line and return Record"""
    fields = line.rstrip().split('\t')
    r = Record(
        seqid = None if fields[0] == '.' else fields[0],
        source = None if fields[1] == '.' else fields[1],
        feature = None if fields[2] == '.' else fields[2],
        start = None if fields[3] == '.' else int(fields[3]),
        end = None if fields[4] == '.' else int(fields[4]),
        score = None if fields[5] == '.' else float(fields[5]),
        strand = None if fields[6] == '.' else fields[6],
        frame = None if fields[7] == '.' else int(fields[7]),
        attribute = None if fields[8] == '.' else fields[8]
        )
    return r



def get_transcripts(stream):
    """Collect transcripts as collections of CDS/exons"""
    transcripts = {}
    
    for line__ in stream:
        if isinstance(line__, str):
            line = line__
        elif isinstance(line__, bytes):
            line = line__.decode()
        else:
            raise InputError('Unsupported stream data')

        if line.startswith('#') or not line.strip():
            continue
        
        f = parse(line)
        if f.feature == 'exon':
            m = re.search('transcript_id "([a-zA-Z0-9\.]+)";', f.attribute)
            if not m:
                raise ValueError('Gene id could not be found')
            index = m.group(1)
            exons = []
            if index not in transcripts:
                transcripts[index] = mRNA(seqid = f.seqid, start = None,
                                          end = None, strand = f.strand,
                                          exons = [], attributes = '')
                
            transcripts[index].exons.append(Exon(seqid = f.seqid,
                                                 start = f.start, end = f.end,
                                                 strand = f.strand))

    return transcripts


def get_splice_sites_from_exons(exons, use_strand):
    """Collect splice sites from list of exons and return as a dictionary"""
    sites = {}
    sorted_exons = sorted(exons, key=lambda x: x.start)
    for f, s in zip(sorted_exons, sorted_exons[1:]):
        strand = None
        if use_strand:
            strand = f.strand
        sites[Intron(seqid = f.seqid, start = f.end + 1, end = s.start - 1,
                     strand = strand)] = 1
        
    return sites

def get_splice_sites(stream, use_strand = False):
    """Get splice sites from a mRNAs and returs and a dictionary of introns"""
    sites = {}
    transcripts = get_transcripts(stream)
    for r in transcripts:
        s = get_splice_sites_from_exons(transcripts[r].exons, use_strand)
        for i in s:
            sites[i] = 1
    return sites



if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(description='Generate a list of introns from a GTF file')
    parser.add_argument('gtffile', metavar='FILE', type=str, help='GFF file')

    args = parser.parse_args()

    with open(args.gtffile) as f:
        introns = get_splice_sites(f, use_strand = True)

    for k in introns:
        print(f'{k.seqid}\t{k.start - 2}\t{k.end}\t{k.strand}')

