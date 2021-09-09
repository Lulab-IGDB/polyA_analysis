#!/usr/bin/env python
import sys
if len(sys.argv) != 2:
        print('\tneed: infile in bed format')
        print('\tPlease contact WJQ if you have any questions')
        quit()

infile = sys.argv[1]
tag_dict = {}
infile = open(infile, 'r')
for line in infile:
        (chr, start, end, read_id, MAPQ, strand) = line.strip('\n').split('\t')
	(read_name, sample_barcode_UMI, ccs_pass) = read_id.split(':')
	tag = sample_barcode_UMI + ':' + chr + ':' + start + ':' + end
        tag_dict[tag] = read_id
	
for tag, read_id in tag_dict.items():
        print(read_id)
