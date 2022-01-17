#!/usr/bin/env python
import os, sys

bedfile = sys.argv[1]
mapfile = sys.argv[2]
notefile = sys.argv[3]
outprefix = sys.argv[4]

outfile = open(os.path.join(outprefix + '.polyA_onAPA.txt'), 'w')
outfile.write('read_name\tEnsm_id\ttag\tP_D\torder\tcounts\tend\tAPA\tchrom\tstrand\tpass\tA\tT\tC\tG\tnonA\ttail_seq\tqual_avg\n')

note_dict = {}
notefile = open(notefile,'r')
for line in notefile:
	(sample, read_name, Ensm_id, _pass, _1, A, T, C, G, nonA, _0, tail_seq, qual_avg) \
		= line.strip('\n').split('\t')
	note_dict[read_name] = (Ensm_id, _pass, A, T, C, G, nonA, tail_seq, qual_avg)

APA_dict = {}
mapfile = open(mapfile,'r')
for line in mapfile:
	(Ensm_id, chrom, strand, total_reads_c, peaks_number, loci, loci_reads_c) \
		= line.strip('\n').split('\t')

	if loci != 'loci':
		if loci != '':
			list_loci = map(int, loci.split(','))
			if strand == '-':
				list_loci = [x * (-1) for x in list_loci]
			APA_loci = sorted(list_loci)
		else:
			APA_loci = []
		APA_dict[Ensm_id] = APA_loci

APA_sites = {}
bedfile = open(bedfile,'r')
for line in bedfile:
	(chrom, start, end, read_name, XX, strand) = line.strip('\n').split('\t')
	if read_name in note_dict:
		(Ensm_id, _pass, A, T, C, G, nonA, tail_seq, qual_avg) = note_dict[read_name]
		if Ensm_id in APA_dict:
			APA_loci = APA_dict[Ensm_id]
			APA_c = len(APA_loci)
			if strand == '-':
				end = int(start) * (-1) -1

			tag = ''
			if APA_loci == []:
				tag = 'noAPA'
				key = (Ensm_id, 'no')
				site = 'no' 
			else:
				if int(end) > APA_loci[-1]+5:
					tag = 'Longer'
					key = (Ensm_id, 'ex')
					site = 'extra'

				for i in range(APA_c):
					if APA_loci[i]-5 <= int(end) <= APA_loci[i]+5:
						tag = 'Class_1'
						key = (Ensm_id, i+1)
						site = APA_loci[i]
					if APA_loci[i]-10 <= int(end) < APA_loci[i]-5:
						tag = 'Class_2'
						key = (Ensm_id, i+1)
						site = APA_loci[i]

			if tag == '':
				tag = 'Class_3'
				key = (Ensm_id, 'de')
				site = 'degra'

			P_D = 'NA'
			if tag == 'Class_1' and APA_c >=2:
				(Ensm_id, order) = key
				if order == 1:
					P_D = 'Proximal'
				if order == APA_c:
					P_D = 'Distal'

			end = abs(int(end))
			APA_sites.setdefault(key, []).append((read_name, Ensm_id, chrom, \
				strand, _pass, A, T, C, G, nonA, tail_seq, qual_avg, end, site, tag, P_D))

for key, notes in APA_sites.items():
	[Ensm_id, order] = key
	counts = len(notes)
	for note in notes:
		(read_name, Ensm_id, chrom, strand, _pass, A, T, C, G, nonA, tail_seq, qual_avg, \
			end, site, tag, P_D) = note
		outfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % \
			(read_name, Ensm_id, tag, P_D, order, counts, end, site, chrom, strand, \
			_pass, A, T, C, G, nonA, tail_seq, qual_avg))
outfile.close()
