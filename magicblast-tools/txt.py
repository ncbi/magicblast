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

"""Parse a text junctions file"""

from base import Intron

def get_introns(stream):
    """Get introns locations and return them as a dictionary"""
    introns = {}
    for line in stream:
        fields = line.rstrip().split()[1].split(':')
        seqid = fields[0]
        ff = fields[1].split('-')
        start = int(ff[0])
        end = int(ff[1])
        introns[Intron(seqid = seqid, start = start + 1, end = end - 1,
                       strand = None)] = 1

    return introns

