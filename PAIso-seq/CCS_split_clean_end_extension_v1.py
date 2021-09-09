#!/usr/bin/env python
import re
import sys
import pysam
import regex
import parasail
from itertools import islice

def revComp(s):
        revCDict = {'G':'C','C':'G','T':'A','A':'T','N':'N','-':'-'}
        return ''.join([revCDict[x] for x in s])[::-1]

def splitCCS(ccs_name, ccs_pass, strand, seq, qual, barcode_dict, mismatch):
        split = []
        all_bc_loci = [(0, 0)]

        for barcode_name, barcode_seq in barcode_dict.items():
                bc_loci = []
                pattern = "("+ barcode_seq + "){e<=" + str(mismatch) + "}"
                matches = regex.finditer(pattern, seq, overlapped = False)
                for it in matches:
                        bc_loci.append((barcode_name, it.end()))
                all_bc_loci.extend(bc_loci)
        sort_loci = sorted(all_bc_loci, key=lambda all_bc_loci: all_bc_loci[1])

        for i in range(len(sort_loci)-1):
                barcode_name = sort_loci[i+1][0]
                barcode_seq = barcode_dict[barcode_name]
                ccs_id = ccs_name + strand + str(i+1) + ':' + barcode_name + ':' + ccs_pass
                target_seq = seq[sort_loci[i][1]:sort_loci[i+1][1]]
                target_qual = qual[sort_loci[i][1]:sort_loci[i+1][1]]
                transcript = (barcode_name, barcode_seq, ccs_name, ccs_id, strand, target_seq, target_qual)
                split.append(transcript)
        return split

def removeAdapterStart(adapter, seq, qual):
        gap_open, gap_extend = 8, 4
        score_matrix = parasail.matrix_create("ACGT", 5, -5)
        pairwise = parasail.sg_dx_trace(adapter, seq, gap_open, gap_extend, score_matrix)
        aligned_adapter = pairwise.traceback.query
        aligned_seq = pairwise.traceback.ref
        pattern = r"^(-*)([ATCG-]+?)(-*)$"
        matches = re.finditer(pattern, aligned_adapter)
        for match in matches:
                start = len(match.group(1))
        seq = seq[:start]
        qual = qual[:start]
        return seq, qual

def removeAdapterEnd(adapter, seq, qual):
        matches = regex.finditer(adapter, seq, overlapped = False)
        for it in matches:
                end = it.end()
        try:
                seq = seq[end:]
                qual = qual[end:]
        except:
                pass
        return seq, qual

def oneAdapter(adapter, seq, qual):
        matches = regex.finditer(adapter, seq, overlapped = False)
        for it in matches:
                start = it.start()
        try:
                seq = seq[start:]
                qual = qual[start:]
        except:
                pass
        return seq, qual

infile = sys.argv[1]
passfile = sys.argv[2]
barcodefile = sys.argv[3]
mismatch = int(sys.argv[4])

if len(sys.argv) != 5:
        print('\tneed: infile, passfile, barcodefile, mismatch')
        print('\tPlease contact WJQ if you have any questions')
        quit()

barcode_dict = {}
barcodefile = pysam.FastxFile(barcodefile)
for read in islice(barcodefile, None):
        barcode_dict[read.name] = read.sequence[:22]

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

        strand = '+'
        seq = ccs.sequence
        qual = ccs.quality
        split = splitCCS(ccs_name, ccs_pass, strand, seq, qual, barcode_dict, mismatch)

        strand = '-'
        seq = revComp(ccs.sequence)
        qual = ccs.quality[::-1]
        split.extend(splitCCS(ccs_name, ccs_pass, strand, seq, qual, barcode_dict, mismatch))
        splits.extend(split)

        if len(splits) == 0:
                sys.stderr.write(ccs_name + '\tNO_Barcode\n')

for read in splits:
        end = ''
        (barcode_name, barcode_seq, ccs_name, ccs_id, strand, seq, qual) = read
        seq, qual = removeAdapterStart(barcode_seq, seq, qual)

        tso3 = '(GTACTCTGCGTTGATACCACTGCTT){e<=3}'
        seq, qual =  removeAdapterEnd(tso3, seq, qual)

        tso5 = '(AAGCAGTGGTATCAACGCAGAGTAC){e<=3}'
        seq, qual = oneAdapter(tso5, seq, qual)

        p5 = '(AAGCAGTGGTATCAACGCAGAGTAC){e<=3}(ATGGG){e<=2}'
        matches = regex.finditer(p5, seq, overlapped = False)
        for it in matches:
                end = it.end()
        if end != '':
                seq = seq[end:]
                qual = qual[end:]
                if len(seq) > 50:
                        out = [ccs_name, ccs_id, barcode_name, 'GOOD', strand, seq, qual]
                        out = '\t'.join([str(x) for x in out])
                        print(out)
                else:
                        err = [ccs_name, ccs_id, barcode_name, strand, 'No_sequence', seq, qual]
                        err = '\t'.join([str(x) for x in err])
                        sys.stderr.write(err + '\n')
        else:
                tso5 = '(AAGCAGTGGTATCAACGCAGAGTAC){e<=3}'
                seq, qual =  removeAdapterEnd(tso5, seq, qual)
                if len(seq) > 50:
                        out = [ccs_name, ccs_id, barcode_name, 'NO_P5', strand, seq, qual]
                        out = '\t'.join([str(x) for x in out])
                        print(out)
                else:
                        err = [ccs_name, ccs_id, barcode_name, strand, 'No_sequence', seq, qual]
                        err = '\t'.join([str(x) for x in err])
                        sys.stderr.write(err + '\n')

logfile = sys.argv[1] + '.log.txt'
log_file = open(logfile, 'w')
log_file.write(str(ccs_count) + '\n')
log_file.close()

