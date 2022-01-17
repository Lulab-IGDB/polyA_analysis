#!/usr/bin/env python
import os, sys, pysam, numpy

def getPeaks(depth):
	N = len(depth)
	w = 5
	minDist = 20
	minSupport = 10
	peaks = []
	counts = []
	while True:
		currPeaks = numpy.zeros(N)
		for i in xrange(N):
			for p in peaks:
				if abs(i-p) < minDist-1: 
					break
			else:
				if numpy.sum(depth[max(0,i-w):min(N,i+w+1)]) >= minSupport:
					currPeaks[i] = depth[i]*2 + numpy.median(depth[max(0,i-w):min(N,i+w+1)])
		if numpy.max(currPeaks) == 0:
			break
		cp = numpy.argmax(currPeaks)
		if cp not in peaks:
			peaks.append(cp)
			counts.append(int(numpy.sum(depth[max(0,cp-w):min(N,cp+w+1)])))

	return peaks, counts

notefile = sys.argv[1]
bedfile = sys.argv[2]
outprefix = sys.argv[3]

outfile = open(os.path.join(outprefix + '.APAsites.csv'), 'w')
outfile.write('ensembl_id\tchrom\tstrand\ttotal_reads\tpeaks\tloci\tloci_reads\n')

note_dict = {}
notefile = open(notefile,'r')
for line in notefile:
	(sample, read_name, Ensm_id, read_pass, _1, A, T, C, G, nonA, _0, tail_seq, qual_avg) = \
		line.strip('\n').split('\t')
	note_dict[read_name] = Ensm_id

polyAMap = {}
bedfile = sys.argv[2]
bedfile = open(bedfile,'r')
for line in bedfile:
	(chrom, start, end, read_name, XX, strand) = line.strip('\n').split('\t')
	if read_name in note_dict:
		Ensm_id = note_dict[read_name]
		if strand == '-':
			end = int(start) + 1
		polyAMap.setdefault((Ensm_id, chrom, strand), []).append(int(end))

for key, loci in polyAMap.items():
	(Ensm_id, chrom, strand) = key
	depth = numpy.zeros(max(loci) - min(loci) + 1, int)
	for locus in loci:
			depth[locus - min(loci)] += 1
	peaks, counts = getPeaks(depth)

	outfile.write('%s\t%s\t%s\t%d\t%d\t%s\t%s\n' % \
		(Ensm_id, chrom, strand, len(loci), len(peaks), \
		','.join([str(x + min(loci)) for x in peaks]), \
		','.join([str(y) for y in counts])))

outfile.close()
