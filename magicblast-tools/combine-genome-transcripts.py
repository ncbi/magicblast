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

"""Compare mapping to genome with mapping to transcripts"""

import sam
import pysam
import gff
import gtf
import argparse
import bisect
import sys
import contextlib


def is_equal(query_name_1, query_name_2):
    """Determine whether two read ids represent the same read"""
    q1 = query_name_1
    q2 = query_name_2
    if q1[-1] in ['a', 'b']:
        q1 = q1[:-1]
    if q2[-1] in ['a', 'b']:
        q2 = q2[:-1]

    return q1 == q2

def get_aligns(stream):
    """Get all alignments for a single reads. The SAM or BAM file must be
    sorted by read name"""

    result = []
    for align in stream:
        if len(result) == 0 or \
               is_equal(result[0].query_name, align.query_name):
            result.append(align)
        else:
            yield result
            result = [align]

    if len(result) > 0:
        yield result


def index_aligns(aligns):
    """Index alignments by reference name and reference start position"""
    result = {}
    for i in aligns:

        k = (i.reference_name, i.reference_start)
        result[k] = i

    return result


def get_score(aligns):
    """Get alignment score, composite for paired reads"""
    # for single reads return score of the first alignment
    # we assume that only top scoring alignments are reported
    if aligns[0].flag & 1 == 0:
        return aligns[0].get_tag('AS') if aligns[0].flag & 4 == 0 else 0

    # for paired reads report sum of scores for properly paired alignments
    # (bit 2 set), and single read score for other cases
    scores = []
    forward = {}
    reverse = {}
    for i in aligns:
        # if alignment is not properly paired, save score
        if i.flag & 2 == 0:
            scores.append(i.get_tag('AS') if i.flag & 4 == 0 else 0)
            continue
        
        # for properly paired alignments, find the mate alignment and add
        # scores
        d = None
        if i.flag & 64:
            d = forward
        elif i.flag & 128:
            d = reverse
        else:
            raise ValueError(f'Unrecognised paired flags for alignment: {i}')

        k = (i.reference_name, i.reference_start)
        d[k] = i

    for i in forward.values():
        k = (i.next_reference_name, i.next_reference_start)
        if k not in reverse:
            raise ValueError(f'Missing mate for {i}')
        scores.append(i.get_tag('AS') +  reverse[k].get_tag('AS'))

    return max(scores)
            

def compare_alignments(align_1, align_2):
    """Compare alignmet scores.
    Return 1 if align_1 score is larger than align_2 score, zero if
    align_1 score == align_2 score, and -1 if align_1 score is smaller than
    align_2 score"""
    score_1 = get_score(align_1)
    score_2 = get_score(align_2)

    if score_1 == score_2:
        return 0
    elif score_1 > score_2:
        return 1
    return -1


def get_alignment_end(align):
    """Find end position of a SAM alignment"""
    if align.flag & 4:
        raise ValueError(f'Read unaligned: {align}')

    pos = align.reference_start
    for op, num in align.cigartuples:
        if op in [sam.CIGAR_MATCH, sam.CIGAR_DELETION, sam.CIGAR_INTRON]:
            pos += num

    return pos
    

def transcript2genome(positions, align, transcript):
    """Translate a position on a transcript to position on a genome"""

    # find cumulative exon lengths
    exon_lens = [e.end - e.start + 1 for e in transcript.exons]
    cum_exon_lens = [sum(exon_lens[:i + 1]) for i in range(0, len(exon_lens))]
        
    for pos in positions:

        # if transcript is annotated on the negative negative strand of the
        # genome, go from the end of the read alignment to the transcript:
        # find alignment position from sequence end and reverse CIGAR
        if transcript.strand == '-':
            align_start = sum(exon_lens) - get_alignment_end(align) + 1

        ind = bisect.bisect_left(cum_exon_lens, align_start)
        if ind > 0 and cum_exon_lens[ind] == align_start:
            ind += 1
        if ind == 0:
            yield transcript.exons[ind].start + align_start - 1
        else:
            yield transcript.exons[ind].start + align_start - 1 - \
                      cum_exon_lens[ind - 1]


def reverse_complement(sequence):
    """Reverse complement a sequence"""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    result = []
    for i in reversed(sequence):
        result.append(complement.get(i, 'N'))
    return ''.join(result)


def remap_to_genome(align, annot, mates = None):
    """Remap alignment to a transctript to alignment on a genome. Generates
    a string, single line SAM output"""
    # if a read is unaligned retrun the same SAM line
    if align.flag & 4:
        return align.to_string()

    if align.reference_name not in annot:
        raise ValueError(f'{align.reference_name} not present in annotation')

    transcript = annot[align.reference_name]

    flag = align.flag
    sequence = align.query_sequence
    if transcript.strand == '-':
        flag ^= 16
        sequence = reverse_complement(sequence)

    exon_lens = [e.end - e.start + 1 for e in transcript.exons]

    # alignment of read to transcript left-most position 
    align_start = align.reference_start + 1
    cigartuples = align.cigartuples

    # if transcript is annotated on the negative negative strand of the genome
    # go from the end of the read alignment to the transcript:
    # find alignment position from sequence end and reverse CIGAR
    if transcript.strand == '-':
        align_start = sum(exon_lens) - get_alignment_end(align) + 1
        cigartuples = reversed(cigartuples)

    # find read alignment to the genome start position
    cum_exon_lens = [sum(exon_lens[:i + 1]) for i in range(0, len(exon_lens))]
    ind = bisect.bisect_left(cum_exon_lens, align_start)
    if ind > 0 and cum_exon_lens[ind] == align_start:
        ind += 1
    if ind == 0:
        start = transcript.exons[ind].start + align_start - 1
    else:
        start = transcript.exons[ind].start + align_start - 1 - \
                cum_exon_lens[ind - 1]

    op2cigar = {sam.CIGAR_MATCH: 'M', sam.CIGAR_INSERTION: 'I',
               sam.CIGAR_DELETION: 'D', sam.CIGAR_SOFT_CLIP: 'S'}

    # generate new CIGAR
    cigar = ''
    g_pos = start
    for op, num in cigartuples:
        if op == sam.CIGAR_INTRON:
            raise RuntimeError('Bad alignment')

        if op in [sam.CIGAR_MATCH, sam.CIGAR_DELETION]:

            bases_left = num
            while bases_left > 0:

                # if rthe current alignment segment ends within current exon
                if g_pos + bases_left <= transcript.exons[ind].end:
                    cigar += f'{bases_left}{op2cigar[op]}'
                    g_pos += bases_left
                    break
                else:

                    # if the current alignment segment spans an intron
                    exon_bases = transcript.exons[ind].end - g_pos + 1
                    bases_left -= exon_bases
                    g_pos += exon_bases
                    cigar += f'{exon_bases}{op2cigar[op]}'
                    if bases_left > 0:
                        intron = transcript.exons[ind + 1].start - transcript.exons[ind].end - 1
                        cigar += f'{intron}N'
                        g_pos += intron
                        ind += 1

        elif op in [sam.CIGAR_INSERTION, sam.CIGAR_SOFT_CLIP]:
            cigar += f'{num}{op2cigar[op]}'
        else:
            raise ValueError(f'Unsupported CIGAR operation: {op}')

    # compute mate start for paired reads
    mate_name = '*'
    mate_start = 0
    if align.flag & 1 and align.next_reference_name is not None :
        if align.next_reference_name == align.reference_name:
            mate_name = '='
            mate_start = transcript2genome([align.next_reference_start], align,
                                           transcript)

        else:
            # get alignment start for a new transcript
            mate_transcript = annot[align.next_reference_name]
            mate_name = annot[align.next_reference_name].seqid

            mate_align = mates[(align.next_reference_name,
                                align.next_reference_start)]

            mate_start = transcript2genome([align.next_reference_start],
                                           mate_align,
                                           mate_transcript)


    result = f'{align.query_name}\t{flag}\t{annot[align.reference_name].seqid}'\
             f'\t{start}\t255\t{cigar}\t{mate_name}\t{mate_start}\t0' \
             f'\t{sequence}\t*'\
             f'\tNH:i:{align.get_tag("NH")}\tAS:i:{align.get_tag("AS")}'\
             f'\tNM:i:{align.get_tag("NM")}'

    return result


@contextlib.contextmanager
def wopen(filename = None, mode = 'w', template = None):
    """Open stream for writing either to stdout or a file"""
    if 'b' in mode:
        ff = pysam.AlignmentFile(filename, mode, template = template)
    else:
        if filename and filename != '-' and filename != 'stdout':
            ff = open(filename, 'w')
        else:
            ff = sys.stdout
        print(template.text, file=ff)

    try:
        yield ff
    finally:
        if ff is not sys.stdout:
            ff.close()
            

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Add annotation to Magic-BLAST mapping')
    parser.add_argument('--to-genome', metavar='FILE', dest='genome', type=str,
                        help='BAM file with mapping to a genome')
    parser.add_argument('--to-transcripts', metavar='FILE', dest='transcripts',
                        type=str, help='BAM file with mapping to a transcripts')
    parser.add_argument('--gff', metavar='FILE', dest='gfffile', type=str,
                        help='Genome annotations file')
    parser.add_argument('--out', metavar='FILE', dest='outfile', type=str,
                        help='Output SAM file', default='-')
    parser.add_argument('-b', dest='isbam', action='store_true',
                        help='Output BAM')
                        

    args = parser.parse_args()
                            
    # read annotation file
    with open(args.gfffile) as f:
        if args.gfffile.endswith('gff'):
            m = gff.get_mrnas(f)

            # index mRNAs by accession
            mrnas = {}
            for k in m.keys():
                new_key = ':'.join([m[k].attributes['Name'], m[k].seqid])
                if not new_key.startswith('NM_'):
                    continue

                if new_key in mrnas:
                    raise ValueError(f'{new_key} already present in mrnas')

                m[k].exons.sort(key = lambda x: x.start)
                mrnas[new_key] = m[k]
        else:
            m = gtf.get_transcripts(f)
            mrnas = m

    mode = 'w'
    isbam = args.isbam or args.outfile.endswith('bam')
    if isbam:
        mode = mode + 'b'

    # read and compare BAM files
    with sam.open_sam_or_bam(args.genome) as fg, sam.open_sam_or_bam(args.transcripts) as ft, wopen(args.outfile, mode = mode, template = fg) as out:
        for genome, transcript in zip(get_aligns(fg), get_aligns(ft)):
            
#            print(f'genome:\t{genome[0]}')
#            print(f'transcript\t{transcript[0]}')
#            print('')

            duplicates = {}
            if compare_alignments(genome, transcript) >= 0:
                for i in genome:
                    if isbam:
                        out.write(i)
                    else:
                        print(i.to_string(), file=out)
            else:

                if transcript[0].flag & 1:
                    indexed = index_aligns(transcript)

                for i in transcript:
                    remapped = remap_to_genome(i, mrnas, indexed)
                    fields = remapped.split()
                    key = (fields[2], fields[3])
                    # check if a similar alignment was already produced;
                    # alignments to RNA variants may result in the
                    # same alignment on the genome
                    if key not in duplicates:
                        duplicates[key] = 1
                        if isbam:
                            out.write(pysam.AlignedSegment.fromstring(
                                remapped, fg.header))
                        else:
                            print(remapped, file=out)

