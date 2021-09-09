#!/usr/bin/env python
import os
import sys
import statistics

def toDigit(ascii_list):
	digit_list = [(ord(i)-33) for i in ascii_list]
	return digit_list

if len(sys.argv) != 3:
        msg = ['python',sys.argv[0], 'XXX.polyA_trim.out.txt','XXX.clean.filter.bam.featureCounts','>sample.polyA_tail_length.txt']
        msg = ' '.join([ str(x) for x in msg])
        print(msg)
        quit()

fc_dict = {}
fcounts = sys.argv[2]
fcounts = open(fcounts,'r')
for line in fcounts:
        (name, status, strand, Ensm_id) = line.strip('\n').split('\t')
        if status == 'Assigned':
                fc_dict[name] = Ensm_id
        else:
                fc_dict[name] = 'Unknown'

infile = sys.argv[1]
infile = open(infile, 'r')
for line in infile:
        (tag, name, New, tail_seq, tail_qual, Old, old_tail_seq, old_tail_qual) = line.strip('\n').split('\t')
	if tag == 'TRUE':
	        (read_name, sample, read_pass) = name.split(':')

        	A = tail_seq.count('A')
	        T = tail_seq.count('T')
	        C = tail_seq.count('C')
	        G = tail_seq.count('G')
	        nonA = T + C + G

        	Ensm_id = ''
	        if name in fc_dict:
        	        Ensm_id = fc_dict[name]
	        else:
	                sys.stderr.write('Please confirm that the input file is corresponding')
			quit()
		if tail_qual != 'No':
			tail_qual = toDigit(tail_qual)
			qual_avg = statistics.mean(tail_qual)
		else:
			qual_avg = 93.0

        	out = [sample, name, Ensm_id, read_pass, 1, A, T, C, G, nonA, 0, tail_seq, qual_avg]
	        out = '\t'.join([str(x) for x in out])
        	print(out)

