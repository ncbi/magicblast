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

"""GFF file parser"""

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
    try:
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
    except ValueError:
        print(line)

    return r


def get_introns(stream):
    """Collect introns from gff stream and return as a dictionary"""
    introns = {}
    for line in stream:
        if line.startswith('#'):
            continue

        r = parse(line)
        if r.feature != 'intron':
            continue

        introns[Intron(seqid = r.seqid, start = r.start, end = r.end,
                       strand = r.strand)] = 1

    return introns


def get_mrnas(stream, source = None):
    """Collect mRNA extents with exons"""
    mrnas = {}
    
    for line__ in stream:
        if isinstance(line__, str):
            line = line__
        elif isinstance(line__, bytes):
            line = line__.decode()
        else:
            raise InputError('Unsupported stream data')

        if line.startswith('#'):
            continue
        
        f = parse(line)

        if source is not None and f.source != source:
            continue
        
#        if f.feature in ['mRNA', 'transcript', 'primary_transcript', 'miRNA',
#                         'lnc_RNA', 'gene', 'snoRNA', 'antisense_RNA']:
#        if f.feature in ['mRNA', 'transcript']:
        if f.feature == 'mRNA':
#            m = re.search('ID=(rna\d\d*);', f.attribute)
            m = re.search('ID=(\w[\w:]*)', f.attribute)
            if not m:
                raise ValueError('mRNA id could not be found')
            index = m.group(1)
            exons = []
            if index in mrnas:
                if mrnas[index].start is not None:
                    raise RuntimeError('mRNA with the same id is already present')
                exons = mrna[index].exons

            attributes = {}
            for a in f.attribute.rstrip().split(';'):
                r = a.split('=')
                if (len(r) == 2):
                    attributes[r[0]] = r[1]

            mrnas[index] = mRNA(seqid = f.seqid, start = f.start, end = f.end,
                                strand = f.strand, exons = exons,
                                attributes = attributes)

        if f.feature == 'exon':
#            m = re.search('Parent=(rna\d\d*)', f.attribute)
            m = re.search('Parent=(\w[\w:]*)', f.attribute)
            if not m:
                # there seem to be exons not assigned to mRNAs
#                raise ValueError('Exon without parent: {0}'.format(f.attribute))
                print('WARNING: Exon without parent: {0}'.format(line))
                continue

            index = m.group(1)
            if index not in mrnas:
#                raise  ValueError('Parent of the exon not found: {0}'.format(line))
#                mrnas[index] = mRNA(seqid = f.seqid, start = None, end = None,
#                                    strand = None, exons = [])
                continue

            mrnas[index].exons.append(Exon(seqid = f.seqid, start = f.start,
                                            end = f.end, strand = f.strand))

        # these things appear in RNA-seq, but are not mRNA
        if f.feature == 'five_prime_UTR':
            m = re.search('ID=(id\d\d*);', f.attribute)
            if not m:
                raise ValueError("5'UTR id could not be found")
            
            index = m.group(1)
            if index not in mrnas:
                mrnas[index] = mRNA(seqid = f.seqid, start = None, end = None,
                                    strand = None, exons = [], attributes = {})

            mrnas[index].exons.append(Exon(seqid = f.seqid, start = f.start,
                                           end = f.end, strand = f.strand))

        # this is a hack
        # there are introns with no exons in the gff file, we create fake
        # exons for easier processing
        if f.feature == 'intron':
            m = re.search('ID=(id\d\d*);', f.attribute)
            if not m:
                raise ValueError('Intron id could not be found')
            
            index = m.group(1)
            if index in mrnas:
                raise RuntimeError('mRNA element already present')

            exon1 = Exon(seqid = f.seqid, start = f.start - 2,
                         end = f.start - 1, strand = f.strand)
            exon2 = Exon(seqid = f.seqid, start = f.end + 1, end = f.end + 2,
                         strand = f.strand)
            mrnas[index] = mRNA(seqid = f.seqid, start = None, end = None,
                                strand = None, exons = [exon1, exon2],
                                attributes = {})


    return mrnas

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

def get_splice_sites(stream, source = None, use_strand = False, accession = None):
    """Get splice sites from a mRNAs and returs and a dictionary of introns"""
    sites = {}
    mrnas = get_mrnas(stream, source)
    for r in mrnas:
        if accession is not None and 'Name' in mrnas[r].attributes and \
           not mrnas[r].attributes['Name'].startswith(accession):
            continue

        s = get_splice_sites_from_exons(mrnas[r].exons, use_strand)
        for i in s:
            sites[i] = 1
    return sites



if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(description='Generate a list of introns from a GFF file')
    parser.add_argument('gfffile', metavar='FILE', type=str, help='GFF file')

    args = parser.parse_args()

    with open(args.gfffile) as f:
        introns = get_splice_sites(f, use_strand = True)

    for k in introns:
        print(f'{k.seqid}\t{k.start}\t{k.end}\t{k.strand}')

            
