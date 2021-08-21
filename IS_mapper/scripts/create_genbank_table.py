#!/usr/bin/env python

from argparse import (ArgumentParser, FileType)
from Bio import SeqIO
from Bio import SeqFeature
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.Blast.Applications import NcbiblastnCommandline
from operator import itemgetter
import os, sys, re, collections, operator
import numpy as np
from collections import OrderedDict

def parse_args():
    '''
    Parse the input arguments, use -h for help.
    '''

    parser = ArgumentParser(description="Create a table of features for ISMapper (improvement pathway)")
    parser.add_argument('--five_bed', type=str, required=True, help='five prime end bed file')
    parser.add_argument('--three_bed', type=str, required=True, help='three prime end bed file')
    parser.add_argument('--assembly', type=str, required=True, help='assembly file for annotation (can be genbank or fasta)')
    parser.add_argument('--type', type=str, required=True, help='the type of assembly file (is either a genbank or fasta)')
    parser.add_argument('--output', type=str, required=True, help='prefix for output file')
    return parser.parse_args()

def create_feature(hit, end):
    '''
    Creates a feature from the hit, setting colours
    based on orientation and noting whether it is a
    left or right end hit.
    '''
    # Set up coordinates
    start = int(hit[1])
    stop = int(hit[2])
    quals = {}

    location = SeqFeature.FeatureLocation(start, stop)
    # Annotate with end information and colour accordingly
    if end == 'five':
        quals['colour'] = '2'
        quals['end'] = 'left_end'
        feat_type = 'left end'
    elif end == 'three':
        quals['colour'] = '7'
        quals['end'] = 'right end'
        feat_type = 'right_end'
    # Create feature
    feature = SeqFeature.SeqFeature(location, type=feat_type, qualifiers=quals)
    
    return feature

def main():

    args = parse_args()

    header = ['contig', 'end', 'x', 'y']
    results = collections.defaultdict(dict)
    # If both bed files empty, then no hits found
    if os.stat(args.five_bed)[6] == 0 and os.stat(args.three_bed)[6] == 0:
        output = open(args.output + '_table.txt', 'w')
        output.write('\t'.join(header) + '\n')
        output.write('No hits found')
        output.close()
        # Exit ISMapper
        sys.exit()

    # Start counting hits
    hit_no = 1
    # Check first that there are left end hits
    if os.stat(args.five_bed)[6] != 0:
        with open(args.five_bed) as five_bed:
            for line in five_bed:
                info = line.strip().split('\t')
                results[info[0]]['hit_' + str(hit_no)] = ['five', info[1], info[2]]
                hit_no += 1
    # Check first that there are right end hits
    if os.stat(args.three_bed)[6] != 0:
        with open(args.three_bed) as three_bed:
            for line in three_bed:
                info = line.strip().split('\t')
                results[info[0]]['hit_' + str(hit_no)] = ['three', info[1], info[2]]
                hit_no += 1

    # Open up table for writing output into
    output = open(args.output + '_table.txt', 'w')
    output.write('\t'.join(header) + '\n')
    # Convert the assembly to a genbank file first
    if args.type == 'fasta':
        print('Creating multi entry genbank for annotation...')
        count = SeqIO.convert(args.assembly, 'fasta', args.output + '.gbk', 'genbank', generic_dna)
        print('Successfully converted %i records' % count)
        # read in file
        record_list = SeqIO.parse(args.output + '.gbk', 'genbank')
    # Otherwise read in the genbank
    elif args.type == 'genbank':
        record_list = SeqIO.parse(args.assembly, 'genbank')
    # Initalise a list where any edited contigs will go
    new_record_list = []
    # Intialise number of features added
    feature_count = 0
    # Go through each contig and see if there is a left or right end hit
    for record in record_list:
        if record.name in results:
            for hit in results[record.name]:
                # If there is, then annotate the hit
                new_feature = create_feature(results[record.name][hit], results[record.name][hit][0])
                record.features.append(new_feature)
                feature_count += 1
                # Add the hit to the output table
                output.write(record.name + '\t' + '\t'.join(results[record.name][hit]) + '\n')
            new_record_list.append(record)
        # Otherwise just pass over this contig
        else:
            new_record_list.append(record)
    # Write out the new genbank file
    SeqIO.write(new_record_list, args.output + '_annotated.gbk', 'genbank')
    print('Added ' + str(feature_count) + ' features to ' + args.output + '_annotated.gbk')
    output.close()

if __name__ == '__main__':
    main()
