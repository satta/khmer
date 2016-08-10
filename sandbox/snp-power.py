#! /usr/bin/env python
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2016, The Regents of the University of California.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
#     * Neither the name of the Michigan State University nor the names
#       of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written
#       permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Contact: khmer-project@idyll.org
# pylint: disable=invalid-name,missing-docstring,no-member

from __future__ import print_function
import argparse
import screed
from khmer import khmer_args
from khmer import khmer_logger
import sys

#import cProfile


def get_parser():
    parser = khmer_args.build_nodegraph_args(descr='Calculate')
    parser.add_argument('-l', '--limit', type=int, default=None,
                        help='stop after this many sequence windows of length '
                        '2k-1 (for debugging or profiling)')
    parser.add_argument('input_sequence_filename', help='The name of the input'
                        ' FASTA sequence file.')
    return parser


def get_windows(sequence, ksize):
    """
    Grab all windows of size 2k-1.

    Ignores the k-1 nucleotides at the end of each sequence.
    """
    for i in range(len(sequence)):
        if i < ksize - 1 or i + ksize > len(sequence):
            continue
        minpos = i - ksize + 1
        maxpos = i + ksize
        yield sequence[minpos:maxpos]


def get_mutations(record, ksize, limit=None):
    """
    Generate all mutations of the middle nucleotide of the sequence.
    """
    i = ksize - 1
    for j, window in enumerate(get_windows(record.sequence, ksize)):
        assert len(window) == (2 * ksize) - 1
        if j > 0 and j % 100000 == 0:
            perc = float(j) / 1000000.0
            khmer_logger.log_info('  processed {:.1f} Mb of sequence {}'.format(perc, record.name.split()[0]))
            if limit is not None and j >= limit:
                break
        for nucl in 'ACGT':
            if nucl != window[i]:
                yield window[:i] + nucl + window[i+1:]


def count_mutation_collisions(args, nodegraph):
    total = 0
    hits = 0
    for i, record in enumerate(screed.open(args.input_sequence_filename)):
        for mutatedseq in get_mutations(record, args.ksize, args.limit):
            for kmer in nodegraph.get_kmers(mutatedseq):
                total += 1
                if nodegraph.get(kmer) > 0:
                    hits += 1
    print(total, hits, float(hits) / float(total))


def main(args):
    khmer_logger.log_info('allocating nodegraph')
    nodegraph = khmer_args.create_nodegraph(args)

    khmer_logger.log_info('consuming input')
    nodegraph.consume_fasta(args.input_sequence_filename)

    khmer_logger.log_info('generating mutations')
    count_mutation_collisions(args, nodegraph)
    #cProfile.runctx('func(r, ng)', {'r': args, 'func': count_mutation_collisions, 'ng': nodegraph}, {});



if __name__ == '__main__':
    main(get_parser().parse_args())
