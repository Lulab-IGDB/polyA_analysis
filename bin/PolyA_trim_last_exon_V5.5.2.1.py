#! /usr/bin/env python
import sys
import pysam

"""
Poly(A) tail extraction method V5.5.2.1
### 1. Extend soft-clips from minimap2 to any number of A to the 5' direction as the candidate tail.
2. Flag reads with more than 10% T and 10% C and 10% G simutaneously "HIGH_TCG".
3. Move tail from 5' to 3' end, if there is a change in nucleotide, score 1. Flag reads with score larger than 12 as "FALSE_score_12+".
4. The remaining reads are flagged as "TRUE".
5. Reads that ends within the last exon was further marked by "Last_" appended to the tag.
"""

def revComp(s):
        revCDict = {'G':'C','C':'G','T':'A','A':'T','N':'N','-':'-'}
        return ''.join([revCDict[x] for x in s])[::-1]

def toASCII(digit_list):
        ascii_list = [(chr(i+33)) for i in digit_list]
        ascii_list = ''.join([str(x) for x in ascii_list])
        return ascii_list

if len(sys.argv) != 4:
        print('\tneed: infile, passfile, barcodefile, mismatch')
        print('\tPlease contact WJQ if you have any questions')
        quit()

gene_id_dict = {}
last_exon = open(sys.argv[1],'r')
for line in last_exon:
        (chr_num, start, end, strand, gene_id, transcript_id, gene_name, exon_number) = line.strip('\n').split('\t')
        if strand == '+':
                gene_id_dict.setdefault(gene_id,[]).append((int(start), int(end)+100))
        if strand == '-':
                gene_id_dict.setdefault(gene_id,[]).append((int(start)-100, int(end)))

fc_dict = {}
fcounts = open(sys.argv[2],'r')
for line in fcounts:
        (name, status, strand, Ensm_id) = line.strip('\n').split('\t')
        if status == 'Assigned':
                fc_dict[name] = Ensm_id
        else:
                fc_dict[name] = 'Unknown'

bamfile = pysam.AlignmentFile(sys.argv[3], 'rb')
for read in bamfile:
        out = []
        Ensm_id = tag = ''
        if read.is_unmapped:
                continue
        if read.is_secondary:
                continue
        if read.is_supplementary:
                continue

        read_name = read.query_name
        read_seq = read.query_sequence
        read_block = read.get_blocks()
        read_qual = list(read.query_qualities)
        read_qual = toASCII(read_qual)

        tail_len = 0
        tail_seq = tail_qual = ''

        if read.is_reverse == False :
                cigar = read.cigartuples[-1][0]
                tail_len = read.cigartuples[-1][1]
                read_end = int(read_block[-1][-1])

        if read.is_reverse == True :
                read_seq = revComp(read_seq)
                read_qual = read_qual[::-1]
                cigar = read.cigartuples[0][0]
                tail_len = read.cigartuples[0][1]
                read_end = int(read_block[0][0])

        if cigar == 4:
                tail_seq = read_seq[-1*tail_len:]
                tail_qual = read_qual[-1*tail_len:]

        old_tail_seq = tail_seq
        old_tail_qual = tail_qual

        # extend polyA
#        for i in range(len(read_seq) - len(tail_seq))[::-1]:
#                if read_seq[i] == 'A':
#                        tail_seq = 'A' + tail_seq
#                        tail_qual =  read_qual[i] + tail_qual
#                else:
#                        break

        # filter based on T C or G cotent
        A = tail_seq.count('A')
        T = tail_seq.count('T')
        C = tail_seq.count('C')
        G = tail_seq.count('G')
        Total = float(A + T + C + G)
        if Total > 0 and T/Total >= 0.1 and C/Total >= 0.1 and G/Total >= 0.1 :
                tag = 'HIGH_TCG'

        else:
        # calculate score
                score = 0
                for i in range(len(tail_seq)-1):
                        if tail_seq[i+1] != tail_seq[i]:
                                score += 1
                if score > 12:
                        tag = 'FALSE_score_12+'
                else:
                        tag = 'TRUE'

                if old_tail_seq == '':
                        old_tail_seq = old_tail_qual = 'No'
                if tail_seq == '':
                        tail_seq = tail_qual = 'No'

        if len(tail_seq) < 10000 and read_name in fc_dict:
                gene_id = fc_dict[read_name]
                if gene_id != 'Unknown':
                        start_end = gene_id_dict[gene_id]
                        for i in range(len(start_end)):
                                start = start_end[i][0]
                                end = start_end[i][1]
                                if read_end >= start and read_end <= end:
                                        tag = 'Last_' + tag
                                        break

#        out = [tag, read_name, "New:", tail_seq, tail_qual, "Old:", old_tail_seq, old_tail_qual, score]
        out = [tag, read_name, "New:", tail_seq, tail_qual, "Old:", old_tail_seq, old_tail_qual]
        out = '\t'.join([str(x) for x in out])
        print(out)

