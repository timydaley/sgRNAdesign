# sgRNAdesign


Design rules in order of priority:

       1. Exact match of PAM (Required)

       2. No exact match of seed sequence (positions 1-8 on 3' end)

       3. If exact match of seed region, then anything must be at least an edit distance of 3 away.


Input:

	1. fastq format file of input regions, determined separately

	2. fastq format file of genome


Output:

	1. Mapped-read format file of candidate sgRNAs 


Workflow:

	1. Find cadidate sgRNAs using PAM (NGG for S. Pyongenes, but can theoretically be anything)

	   a. Use PAM as anchor for sgRNAs

	2. Construct a hash table of seed sequences to look for exact matches to the seed sequence

	3. Iteratively hash genome (no need to load the whole genome)

	   a. Avoid unknown seqs (any base !in {ACGT}: don't test)

	   b. 

	4. Only one exact match (which should be the sgRNA in the genome): return sgRNA

	5. More than one exact match, then do exact alignment of candidates against proposed sgRNA 