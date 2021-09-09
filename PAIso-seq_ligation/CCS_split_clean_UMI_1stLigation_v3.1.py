#!/usr/bin/env python
import sys
import pysam
import regex
from itertools import islice

def revComp(s):
        revCDict = {'G':'C','C':'G','T':'A','A':'T','N':'N','-':'-'}
        return ''.join([revCDict[x] for x in s])[::-1]

def removeAdapterEnd(adapter, seq, qual):
        matches = regex.finditer(adapter, seq, overlapped = False)
        for it in matches:
                end = it.end()
        try:
                seq = seq[end:]
                qual = qual[end:]	
        except:
                seq = seq
                qual = qual
        return seq, qual

def removeAdapterStart(adapter, seq, qual):
        seq = revComp(seq)
        qual = qual[::-1]
        seq, qual = removeAdapterEnd(adapter, seq, qual)
        seq = revComp(seq)
        qual = qual[::-1]
        return seq, qual

def splitCCS(ccs_name, ccs_pass, seq, qual):
        split = []
        loci = [len(seq)]
        tso3 = '(GTACTCTGCGTTGATACCACTGCTT){e<=3}'
        matches = regex.finditer(tso3, seq, overlapped = False)
        for it in matches:
                loci.append(it.end())

        tso5 = '(AAGCAGTGGTATCAACGCAGAGTAC){e<=3}'
        matches = regex.finditer(tso5, seq, overlapped = False)
        for it in matches:
                loci.append(it.start())

        loci_sort = sorted(loci)
        for i in range(len(loci_sort)-1):
                if loci_sort[i+1]-loci_sort[i] < 10:
                        loci_sort[i] = loci_sort[i+1]
        tso_loci = list(set(loci_sort))

        for i in range(len(tso_loci)-1):
                target_seq = seq[tso_loci[i]:tso_loci[i+1]]
                target_qual = qual[tso_loci[i]:tso_loci[i+1]]
                order = str(i+1)
                read = (ccs_name, ccs_pass, order, target_seq, target_qual)
                split.append(read)
        return split

infile = sys.argv[1]
passfile = sys.argv[2]

if len(sys.argv) != 3:
        print('\tneed: infile, passfile, barcodefile')
        print('\tPlease contact WJQ if you have any questions')
        quit()

pass_dict = {}
passfile = open(passfile,'r')
for line in passfile:
        (ccs_name, ccs_pass) = line.strip('\n').split('\t')
        pass_dict[ccs_name] = ccs_pass

splits = []
ccs_count = 0
infile = pysam.FastxFile(infile)
for ccs in islice(infile, None):
        ccs_count = ccs_count + 1
        ccs_name = ccs.name
        ccs_pass = pass_dict[ccs_name]
        seq = ccs.sequence
        qual = ccs.quality
        split = splitCCS(ccs_name, ccs_pass, seq, qual)
        splits.extend(split)

        if len(split) == 0:
                sys.stderr.write(ccs_name + '\tNO_TSO\n')

for read in splits:
        end1 = end2 = end3 = end4 = end5 = ''
        (ccs_name, ccs_pass, order, seq, qual) = read
        tso3 = '(GTACTCTGCGTTGATACCACTGCTT){e<=4}'
        tso5 = '(AAGCAGTGGTATCAACGCAGAGTAC){e<=4}'
        tso5gg = '(AAGCAGTGGTATCAACGCAGAGTAC){e<=4}(ATGGG){e<=2}'

        match1 = regex.finditer(tso3, seq, overlapped = False)
        for it in match1:
                end1 = it.end()

        match2 = regex.finditer(tso5, seq, overlapped = False)
        for it in match2:
                end2 = it.end()

        if end1 != '' and end2 != '':
                seq1, qual = removeAdapterEnd(tso5gg, seq, qual)
                if len(seq1) != len(seq):
                    seq =seq1
                    seq1_1, qual = removeAdapterStart(tso5gg, seq, qual)
                    if len(seq1_1) == len(seq):
                        tag = 'GOOD'
                        strand = '+'
                        seq =seq1_1
                    else:
                        tag = 'P5_5P'
                        strand = '?'
                        seq =seq1_1
                else:
                        seq2, qual = removeAdapterStart(tso5gg, seq1, qual)
                        if len(seq2) != len(seq1):
                            seq = revComp(seq2)
                            qual = qual[::-1]
                            seq2_1, qual = removeAdapterStart(tso5gg, seq, qual)
                            if len(seq2_1) == len(seq):
                                tag = 'GOOD'
                                strand = '-'
                            else:
                                tag = 'P5_5P_rev' 
                                seq = seq2_1
                                strand = '?'

                        else:
                                seq, qual = removeAdapterEnd(tso5, seq2, qual)
                                tag = 'noATGGG'
                                strand = '?'
                seq, qual = removeAdapterStart(tso5, seq, qual)

        else:
                tag = '01_tso'
                strand = '?'

        ccs_id = ccs_name + strand + order + ':noBARnoUMI:' + ccs_pass

        if len(seq) > 50:
                if tag == 'GOOD':
                        out = [ccs_name, ccs_id, 'noBARnoUMI', tag, strand, seq, qual]
                        out = '\t'.join([str(x) for x in out])
                        print(out)
                else:
                        err = [ccs_name, ccs_id, 'noBARnoUMI', tag, strand, seq, qual]
                        err = '\t'.join([str(x) for x in err])
                        sys.stderr.write(err + '\n')

logfile = sys.argv[1] + '.log.txt'
log_file = open(logfile, 'w')
log_file.write(str(ccs_count) + '\n')
log_file.close()

