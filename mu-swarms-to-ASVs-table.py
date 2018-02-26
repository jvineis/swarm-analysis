#!/usr/bin/env python

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='''This script takes the output returned from swarmv2.2.2 -w flag, which is a list of lists '\n' containing the sequence ids for each swarm, and converts it into a '\n' ASV count matrix. It will also write a reduced fasta file for the representative nodes if you want to remove swarms below a particular count threshold''')
parser.add_argument('-s' , help = "The swarm output of the -o flag: a list of ids for each swarm (therefore a list of lists)  -  commonly run like this: swarm -d 1 -f -t 10 -z pooled-samples-derep.fa -s pooled-samples-derep-stats.txt -w pooled-samples-node-representatives.fa -o pooled-samples-node-table.txt")
parser.add_argument('-o' , help = "The filename for the output table with node ids as rows and samples as columns.")
parser.add_argument('-l' , help = "A single column list of the sample names. eg.  for a dereplicated sequence in your analysis that looks like this V1-33-33-5-nrfA_415|M01781:52:000000000-ARME0:1:1101:13516:15760;size=47 from the swarm table, the sample name would be V1-33-33-5-nrfA")
parser.add_argument('-n' , help = "The swarm representative nodes which is the -w flag in swarm")
parser.add_argument('-min', help = "The threshold value to keep a swarm,  anything observed fewer than this many times will be removed")
args = parser.parse_args()

outfile = open(args.o, 'w')
swarms = open(args.s, 'rU')
samples = open(args.l, 'rU')
node_seqs = open(args.n, 'rU')

smpls = []
for line in samples:
    smpls.append(line.strip())
print("Here are your samples")
print(smpls)

rep_seq_swarms = {}
ns = [] # create a list to hold the node names contained in the fasta file
for rec in SeqIO.parse(node_seqs, "fasta"):
    ns.append(rec.id.split(";")[0])
    rep_seq_swarms[rec.id.split(";")[0]] = rec.seq # we are going to need the sequences and headers again.  So we have to make a dictionary here.


ks = [] # a list for the sequence swarm node names that contain a count greater than the min
sw_dic = {} # create a dictionary for the swarm node table.  Use the node name that matches a name in the fasta file header list created on lines 25-27.  

# make a new swarms list that contains swarms with sequence counts greater than the "min"
kn = [] # an empty list the keep the "abundant" swarms 
for node in swarms:
    x = node.strip().split() # create a vaiable for a single swarm list - each sequence in a swarm looks something like this
# V2-0-0-5-nrfA_492|M01781:52:000000000-ARME0:1:1101:4451:19578;size=2. "Sample name | sequence id ; abundance retained from deduplicating"
    sizes = 0
    for line in x: # We can cull out the nodes that total less than the minimum here.
        size = line.strip(';').split('=')[1]
        sizes += int(size)
    if sizes > int(args.min):
        kn.append(x)
        # seach the swarm list for a matiching name in the fasta id file.  There should be a match in there. Othewise there is a problemo

for line in kn:
    x = line[0:len(line)]
    for id in x:
        if id.split(";")[0] in ns:
            sw_dic[id.split(";")[0]] = line # Use the name found in the swarm as the key and the sequence ids in the swarm as the values
            ks.append(id.split(";")[0]) # Append the keeper abundant node names to create a reduced fasta file

print("you have %d total nodes and %d passed the min filter of %s occurrences"%(len(ns), len(sw_dic.keys()), args.min))

# This is to check if you suspect there might be a problem between the names in the swarm-node-table and the node-fasta file
#for i in ns:
#    if i not in sw_dic.keys():
#        print(i, "not found")

def survey(swarm_dic, sample_id):
    id_sums = []
    node_total = []
    for key in swarm_dic.keys():
        tot = 0
        for line in sw_dic[key]:
            sample = line.split("_")[0]
            size = line.strip(";").split('=')[1]
            if sample == sample_id:
                tot += int(size)
        id_sums.append(str(tot))
        node_total.append(tot)
    print(sample_id,sum(node_total))
    return id_sums

print("Here are the sequence counts per sample")
## Use the survey function to look for each sample and the abundance of the sequnece in each list of swarm nodes.

outfile.write("samples"+'\t') # Write each of the swarm node ids to the output file
for key in sw_dic.keys():
    outfile.write(key+'\t')

outfile.write('\n')
for line in smpls:
    outfile.write(line+'\t'+'\t'.join(survey(sw_dic, line))+'\n')

outseq = open("reduced-node-fasta"+'-min'+args.min+'.fa', 'w')
for rec in rep_seq_swarms.keys():
    print rec
    if rec in ks:
        outseq.write(">"+rec +'\n'+str(rep_seq_swarms[rec])+'\n')



outfile.close()

            

