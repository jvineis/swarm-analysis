# swarm-analysis
## This git contains information for processing functional genes nrfA, dsrB, and nirS using swarm and the analysis in R.

1.Lets start by merging and quality filtering the data. 

# Steps for all maker genes
I downloaded each of the marker gene files nirS, dsrB, and nrfA from a (google drive) [https://drive.google.com/drive/u/1/folders/0B77oGC7t35HOTmxfQXJIa1g0WGM] to a directory called oxygen_gradient and set up a separate directory for each sample within a directory of the marker gene.  There are 48 total samples for each of the genes.  

Do something like this to create a directory for each of the samples that will contain each R1 and R2 file/sample

    mv John_nirS nirS
    cd nirS
    ls *.gz | cut -f 1 -d "_" | sort | uniq > samples.txt
    for i in `cat samples.txt`; do mkdir $i; done
    for i in `cat samples.txt`; do mv $i*.fastq $i; done

Make a file that resembles the 00_DEMULTIPLEXING_REPORT(since these reads are already demultiplexed).  The file looks like this

    sample  r1	 r2
    P1-0-0-5	 /Users/joevineis/Dropbox/oxygen_gradient/nirS/P1-0-0-5/P1-0-0-5_S25_L001_R1_001.fastq	/Users/joevineis/Dropbox/oxygen_gradient/nirS/P1-0-0-5/P1-0-0-5_S25_L001_R2_001.fastq
    P1-1-5-2	 /Users/joevineis/Dropbox/oxygen_gradient/nirS/P1-1-5-2/P1-1-5-2_S26_L001_R1_001.fastq	/Users/joevineis/Dropbox/oxygen_gradient/nirS/P1-1-5-2/P1-1-5-2_S26_L001_R2_001.fastq
    P1-16-16-6	 /Users/joevineis/Dropbox/oxygen_gradient/nirS/P1-16-16-6/P1-16-16-6_S30_L001_R1_001.fastq	/Users/joevineis/Dropbox/oxygen_gradient/nirS/P1-16-16-6/P1-16-16-6_S30_L001_R2_001.fastq

# nirS
To filter out all reads with more than three mismatches in the merged region and trim off the primers, create a config-generator.shx bash script that looks like this. 

    #!/bin/bash
    iu-gen-configs 00-file-tomake-ini.txt --r1 ^CCTA[C,T]TGGCCGCC[A,G]CA[A,G]T --r2 ^GCCGCCGTC[A,G]TG[A,C,G]AGGAA

    for i in *.ini
    do
        iu-merge-pairs $i --enforce-Q30-check
    done

    for i in *MERGED
    do
        iu-filter-merged-reads $i -m 3
    done

You may notice that many failed to merge.  I created a summary of the reads lost during this step
   
    sample	raw	filtered	percent-retained
    P1-0-0-5	26052	15446		59.28911408
    P1-1-5-2	22733	13703		60.27800994
    P1-16-16-6	28433	17061		60.00422045
    P1-22-22-5	19105	11818		61.85815232
    P1-3-3-5	17312	9465		54.67305915
    P1-33-33-5	16126	9286		57.5840258


concatenate all of the mismatches-max-3.fa files 'they are named something like this but it might not be exactly correct'

    cat *mismatches-max-3.fa > sequences-to-decompose.fa

Pad the sequences
 
    o-pad-with-gaps sequences-to-decompose.fa -o sequences-to-decompose-padded.fa
    
remove all of the merged filtered data keeping only the sequences-to-decompose-padded.fa file.  Create a new directory in the oxygen_gradient dir called nirS-MED-analysis and run the MED script with a min substative abundance of 20

    decompose ../nirS/sequences-to-decompose-padded.fa -M 20  

# nrfA
The same filtering method was applied to nrfA with the exception of the config_generator.shx.  Here it is shown below

    #!/bin/bash
    iu-gen-configs 00-file-tomake-ini.txt --r1 ^CA[A,G]TG[C,T]CA[C,T]GT[C,G,T]GA[A,G]T --r2 ^T[A,T].GGCAT[A,G]TG[A,G]CA[A,G]TC

    for i in *.ini
    do
        iu-merge-pairs $i --marker-gene-stringent --max-num-mismatches 3
    done

# dsrB
Once again its the same and here is the config_generator.shx

    iu-gen-configs 00-file-tomake-ini.txt --r1-prefix ^CA[T,C]AC.CA[A,G]GG.TGG --r2-prefix ^CAGTT[A,G,T]CC[A,G]CAG[A,T]ACAT
    for i in *.ini
    do
        iu-merge-pairs $i --enforce-Q30-check
    done

    for i in *MERGED
    do
        iu-filter-merged-reads $i -m 3
    done

# SWARM 

1.I recommend running swarm at d = 1.  This provides the highest resolution possible and we can deal with that downstream.  The ouput of swarm required some scripting to make it usable.  First concatenate your merged fasta files, then derplicate your pooled(concatenated) sequences using vsearch, and then run SWARM
    
    cat *MERGED > pooled-samples.fa
    vsearch --derep_fulllength pooled-samples.fa --sizeout --output pooled-samples-derep.fa
    swarm -d 1 -f -t 10 -z pooled-samples-derep.fa -s pooled-samples-derep-stats.txt -w pooled-samples-node-representatives.fa -o pooled-samples-node-table.txt
    
The -f specifies fastidious OTU joining of linked low aboundant nodes

2.Run the mu-swarms-to-ASVs-table.py to convert the SWARM output into a count matrix and filter out low abundant nodes.  The script will also ouput a fasta file with the representative nodes that meet the "min" threshold criteria - called reduced-node-fasta-'min'.fa.  

    python ~/scripts/mu-swarms-to-ASVs-table.py -s pooled-samples-node-table.txt -o swarm-test-count.txt -l ../samples.txt -n pooled-samples-node-representatives.fa -min 10

3.Now you can run famsa to align and trimal to remove large gaps in the alignment and then build a tree with FastTree.

    famsa reduced-node-fasta-min10.fa reduced-node-fasta-min10-famsa.fa
    trimal -in reduced-node-fasta-min10-famsa.fa -out reduced-node-fasta-min10-famsa-trimal.fa -gappyout
    FastTree -nt reduced-node-fasta-min10-famsa-trimal.fa > reduced-node-fasta-min10-famsa-trimal.tre
    vsearch --usearch_global reduced-node-fasta-min10.fa --db ~/scripts/databas/functional_fa_dbs/nrfA_ncbi_parks.fa --blast6out NODE-HITS.txt --id 0.6

4.Build a phyloseq output from these results 

    python ~/scripts/mu-swarm-to-phyloseq-objects.py -tax_ref ~/scripts/databas/functional_fa_dbs/ncbi-parks-LINEAGE-STRINGS.txt -hits NODE-HITS.txt -med swarm-d1-fastidious-min10-matrix-count.txt -fa reduced-node-fasta-min10.fa

5.You could run this as a bash script and run for multiple min count thresholds like this

    #!/bin/bash

    for i in 10 3 0
    do
	mkdir min-count"$i"
    	python ~/Documents/swarm-analysis/mu-swarms-to-ASVs-table.py -s pooled-samples-node-table.txt -o reduced-min"$i"-matrix-count.txt -l ../samples.txt -n pooled-samples-node-representatives.fa -min "$i"
    	vsearch --usearch_global reduced-node-fasta-min"$i".fa --db ~/Documents/BOWEN/functional_gene_db/dsrB/pfam_NIR_SIR/dsrB_ncbi_parks.fa --blast6out NODES-min"$i"-HITS.txt --id 0.6
    	python ~/Documents/swarm-analysis/mu-swarm-to-phyloseq-objects.py -tax_ref ~/scripts/databas/functional_fa_dbs/ncbi-parks-LINEAGE-STRINGS.txt -hits NODES-min"$i"-HITS.txt -med reduced-min"$i"-matrix-count.txt -fa reduced-node-fasta-min"$i".fa
    	mv PHYLOSEQ-MATRIX-OBJECT.txt PHYLOSEQ-TAX-OBJECT.txt reduced-node-fasta-min"$i".fa reduced-min"$i"-matrix-count.txt min-count"$i"
    	famsa min-count"$i"/reduced-node-fasta-min"$i".fa min-count"i"/reduced-node-fasta-famsa"$i".fa
    	trimal -in min-count"$i"/reduced-node-fasta-famsa"$i".fa min-count"i"/reduced-node-fasta-famsa-trimal"$i".fa -gappyout
    	FastTree -nt min-count"$i"/reduced-node-fasta-famsa-trimal"$i".fa > min-count"i"/reduced-node-fasta-famsa-trimal"$i".tre
    done


6.Load this into phyloseq (R) to make beautiful things.  Something like this

    library("phyloseq")
    library("ape")
    mat_nirS = read.table("PHYLOSEQ-MATRIX-OBJECT.txt", header = TRUE, sep = "\t", row.names = 1)
    tax_nirS = read.table("PHYLOSEQ-TAX-OBJECT.txt", header = TRUE, sep = ";", row.names = 1)
    meta_nirS = read.table("nrfA-map.txt", header = TRUE, sep = "\t", row.names = 1)
    tree_nirS = read.tree("nirS_fasttree.tre")    

    mat_nirS = as.matrix(mat_nirS)
    tax_nirS = as.matrix(tax_nirS)

    OTU = otu_table(mat_nirS, taxa_are_rows = TRUE)
    TAX = tax_table(tax_nirS)
    META = sample_data(meta_nirS)

    nrfA_physeq = phyloseq(OTU,TAX,META,tree_nirS)
  

