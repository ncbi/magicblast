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

"""Useful functions for getting information from a SAM/BAM file that work on
top of pysam"""

import pysam
from collections import Counter
from collections import namedtuple
from collections import defaultdict
from base import Intron
import re


CIGAR_MATCH = 0
CIGAR_INSERTION = 1
CIGAR_DELETION = 2
CIGAR_INTRON = 3
CIGAR_SOFT_CLIP = 4


def open_sam_or_bam(filename):
    """Open a SAM or BAM file and return file handle"""
    is_bam = ""
    if filename.endswith('.bam'):
        is_bam = 'b'
    return pysam.AlignmentFile(filename, 'r' + is_bam)


def get_standard_read_name(r, force_single = False, trim = False):
    """Standardize read names. Some programs change read names.
       If force_single is true, .1 and .2 will be added to read names so that
       they can be treated as single. If trim is true, then the last 2
       characters of read name will be trimmed. This is to remove .1 and .2
       from a SAM or BAM file."""
    read_name = r.query_name
    # standardise read accessions
    # hisat puts .R. before read number and removes .1 and .2
    read_name = re.sub(r'\.R\.', '.', read_name)
    if force_single and r.flag & 1:
        if r.flag & 64:
            read_name += '.1'
        else:
            read_name += '.2'

    if trim:
        read_name = read_name[:-2]

    return read_name



def get_introns_from_cigar(position, cigar):
    """Get introns start and stop positions from SAM alignment position and
    cigar string"""

    # pysam cigartuple op codes
    introns = []
    kMatch = 0
    kIns = 1
    kDel = 2
    kIntron = 3
    kSoftClip = 4

    kMatchOnly = 7
    kMismatch = 8

    s_offset = position

    for op, num in cigar:
        if op in [kMatch, kMatchOnly, kMismatch, kDel]:
            s_offset += num
        elif op == kIntron:
            introns.append((s_offset, s_offset + num - 1))
            s_offset += num

    return introns


def get_exons_from_cigar(position, cigar):
    """Get exon start and stop positions from SAM alignment position and CIGAR
    string"""
    exons = []
    s_offset = position
    start = s_offset
    for op, num in cigar:
        if op in [0, 2, 7, 8]:
            s_offset += num
        elif op == 3:
            exons.append((start + 1, s_offset))
            s_offset += num
            start = s_offset

    exons.append((start + 1, s_offset))

    return exons


def get_exons(line):
    """Get a list of exons from a single SAM alignment"""
    return get_exons_from_cigar(line.reference_start, line.cigartuples)


def do_filter(r, filter_by):
    """Apply filter to an alignment, return True if alignment passes"""

    if not isinstance(filter_by, dict):
        raise ValueError('filter_by argument must be a dictionary')

    for k in filter_by:
        if k == 'edit_distance':
            edist = r.get_tag('NM')
            if edist > filter_by[k]:
                return False

        elif k == 'score':
            score = get_score(r)
            if score < filter_by[k]:
                return False

        elif k == 'count':
            count = r.get_tag('NH')
            if count > filter_by[k]:
                return False

        elif k == 'edit_distance_clip':
            edist = 0
            for op, num in r.cigartuples:
                if op == 4:
                    edist += num
            edist += r.get_tag('NM')
            if edist > filter_by[k]:
                return False

        elif k == 'read_id':
            read_name = r.query_name.replace('.R.', '.')
            if read_name not in filter_by[k]:
                return False

        else:
            raise ValueError('Unrecognised filter name: {0}'.format(k))
            
                    
    return True


def get_introns_(stream, force_single, filter_by):
    """Get intron postitions from a SAM stream"""
    introns = Counter()

    for read in stream:

        # skip unaligned reads
        if read.flag & 4:
            continue

        # apply alignment filters
        if filter_by is not None:
            if not do_filter(r, filter_by):
                continue

        # read.reference_start is zero based
        for (f, t) in get_introns_from_cigar(read.reference_start + 1,
                                             read.cigartuples):

            subject = stream.getrname(read.reference_id)
            introns[Intron(seqid = subject, start = f, end = t,
                           strand = None)] += 1

    return introns


def get_introns(filename, force_single, filter_by, trim):
    """Get intron positions"""
    f = open_sam_or_bam(filename)
    introns = get_introns_(f, force_single, filter_by)
    f.close()
    return introns


def get_introns_with_reads_(stream, force_single, filter_by, trim, with_reads,
                            placements):
    """Get intron postitions from a SAM stream"""
    introns = None
    if with_reads:
        introns = defaultdict(set)
    else:
        if placements is None:
            introns = defaultdict(int)
        else:
            introns = defaultdict(float)

    for r in stream:

        # skip unalined reads
        if r.cigartuples is None:
            continue

        # apply alignment filters
        if filter_by is not None:
            if not do_filter(r, filter_by):
                continue

        # read.reference_start is zero based
        for (f, t) in get_introns_from_cigar(r.reference_start + 1, r.cigartuples):
            strand = '+'
            if r.flag & 16 != 0:
                strand = '-'

            subject = stream.getrname(r.reference_id)
            i = Intron(seqid = subject, start = f, end = t, strand = None)
            # standardise read accessions
            # hisat puts .R. before read number and removes .1 and .2
            read_name = get_standard_read_name(r, force_single = force_single,
                                               trim = trim)

            num_clipped = 0
            for op, num in r.cigartuples:
                if op == 4:
                    num_clipped += num

#            if i not in introns:
#                introns[i] = set()
            if with_reads:
                introns[i].add(read_name)
            else:
                if placements is None:
                    introns[i] += 1
                else:
                    introns[i] += 1.0 / placements[read_name]
                
    return introns


def get_introns_with_reads(filename, force_single = False, filter_by = None,
                           trim = False, with_reads = False, placements = None):
    """Get intron positions with reads"""
    f = open_sam_or_bam(filename)
    introns = get_introns_with_reads_(f, force_single = force_single,
                                      filter_by = filter_by, trim = trim,
                                      with_reads = with_reads,
                                      placements = placements)
    f.close()
    return introns


def get_score(line):
    """Compute alignment score from CIGAR and edit distance for a single line
    of SAM/BAM file. This is prefered to relaying on AS tag to be able to
    compare mappers with different scoring schemes"""
    penalty = 8
    score = 0

    # get edit distance
    edit_dist = line.get_tag('NM')

    for op, num in line.cigar:
        # score matches (and mismatches as matches)
        if op == 0:
            score += num
        # score gaps
        elif op in [1, 2]:
            score -= num * penalty
            edit_dist -= num

    # add penalty for mismatches and subtract match scores for them
    score -= (edit_dist * penalty) + edit_dist
    return score
