import os
import optparse
from collections import Counter
import re
from Bio import SeqIO

VERSION = '02.03.2018'
NAME = 'BLAST2single_copy_v2.py'
descr = """Usage:
python BLAST2single_copy_v2.py -b <blast result> -f <original query fasta> -o <output file>
"""
parser = optparse.OptionParser(description=descr)
parser.add_option('-b', '--blast', help='blast result file')
parser.add_option('-f', '--fasta', help='original query fasta')
parser.add_option('-o', '--output', help='output file name')
args = parser.parse_args()[0]

# FUNCTIONS
# Function for creating intervals based on a list of sites with a hit (z) a and a counter of z (d)

def pc(z, d, name, length):
    nums = [ z[0] ]
    result = []
    for k in range(2,max(z)+1):
        if d[k] == d[k-1]:
            nums.append(k)
        if d[k] != d[k-1]:
            result.append(name + ' ' + str(length) + ' ' + str(min(nums)) + '-' + str(max(nums)) + ": " + str(d[k-1]))
            nums = [k]
        if k == z[-1]:
            result.append(name + ' ' + str(length) + ' ' + str(min(nums)) + '-' + str(max(nums)) + ": " + str(d[k-1]))
    return(result)


def BLAST2single(blastout, fasta, outputfile):
    file = open(blastout)
    count = 0
    single_count = 0
    fasta = fasta
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    output = open(outputfile, 'w')
    print >> output, "query\tquery_length\tcopynumber\tpercent_aligned" 
    fulldata = file.read()
    blastversion = fulldata.split('\n')[0]
    for transcript in fulldata.split(blastversion)[1:]:
        count = count + 1
        hits = [] # make a list for end cutted append
        hits_2 = [] # make a list without end cuts
        name =  re.search('# Query: ([\S]+)', transcript).group(1)
        length = len(seq_dict[name].seq)
        for i in transcript.split('\n'):
            if i.split('\t')[0] == name:
                for site in range(int(i.split('\t')[6])+10, int(i.split('\t')[7])+1-10):
                    hits.append(site)
                for site in range(int(i.split('\t')[6]), int(i.split('\t')[7])+1):
                    hits_2.append(site)

        if hits != []:
            hits.sort()
            hits_2_set = set(hits_2)
            mappings = Counter(hits)
            #pc(hits, mappings, name, length)
            mapsum = set()
            for i in pc(hits, mappings, name, length):
                 mapsum.add(i.split(": ")[1])
            if list(mapsum) == ['1', '0']:
                print >>output, name+'\t'+str(length)+'\t'+"single-copy_locus"+"\t"+str(float(len(hits_2_set))/float(length))
                #print name, length, "single-copy_locus"
                single_count = single_count + 1
            else:
                print >>output, name+'\t'+str(length)+'\t'+"multi-copy_locus"+"\t"+str(float(len(hits_2_set))/float(length))
                #print name, length, "multi-copy_locus"
        else:
            print >>output, name+'\t'+str(length)+'\t'+"no_hits_found"+"\t"+"NA"
            #print name, length, "no_hits_found"
    output.close()

    print count
    print single_count

# MAIN
BLAST2single(args.blast, args.fasta, args.output)
