#laSV

laSV (Local Assembly based Structural Variation discovery tools)  
version 1.0.3  
2016-01-28

Jiali Zhuang(jiali.zhuang@umassmed.edu)  
Weng Lab  
Programs in Bioinformatics and Integrative Biology  
University of Massachusetts Medical School

laSV is a software that employs a local de novo assembly based approach to detect genomic structural variations from whole-genome high-throughput sequencing datasets.

Please send questions, suggestions and bug reports to:  
jiali.zhuang@umassmed.edu


laSV is a free software whose distribution is governed by the terms of the GNU General Public License as published by the Free Software Foundation, version 3 or any later version. 

Test dataset can be downloaded from laSV website: http://zlab.umassmed.edu/~zhuangj/laSV/  

##Requirements and installation

laSV runs on Linux X64 systems. The following are the software required by laSV:

bwa (version 0.7.x, http://sourceforge.net/projects/bio-bwa/files/);  
samtools (http://samtools.sourceforge.net/);  
bedtools (https://github.com/arq5x/bedtools2);  
twoBitToFa (part of the UCSC utilities, http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/);  
fastx_collapser (part of the fastx_toolkit, http://hannonlab.cshl.edu/fastx_toolkit/download.html); and  
Perl package BioPerl (http://www.bioperl.org/wiki/Main_Page).  

Please make sure these tools/software are properly installed and included in the path.   

Three files need to be ready in the reference genome information directory (argument of -R option):  
1. FASTA file of the reference genome sequence: GENOME_NAME.fa (GENOME_NAME is the name of the reference genome, e.g., hg19, mm10, dm3, also the argument of option -G);  
2. 2bit file of the reference genome: GENOME_NAME.2bit  
3. RepeatMasker file of the reference genome: GENOME_NAME_rpmk.bed  
4. Chromosome length informtion of the reference genome: GENOME_NAME.chromInfo  
All these file can be found and downloaded from the UCSC genome browser (https://genome.ucsc.edu/).  


For installing laSV first unzip the file. Then in the directory “laSV_v1.x.x/” execute:  
bash ./install.sh  
make MAXK=k  
where k is the the maximal k-mer size value (please supply either 31, 63, or 95).  



##Execution


To execute laSV just run the bash script run_laSV.sh. The various options and parameters are explained below:


        -i       Input FASTQ sequence file prefix with full path. 


        -I       Input compressed FASTQ sequence file prefix with full path. 


        -D       Path to laSV directory


        -o       Path to output directory. Default is the current directory


        -R       Path to directory where the reference genome information is placed


        -G       Reference genome name (e.g. hg19, mm10, dm3, etc.)


        -s       Read Phred score threshold for building de Bruijn graph. Default is 20


        -L       The height of the hash table for storing the de Bruijn graph. Default is 24


        -W       The width of the hash table for storing the de Bruijn graph. Default is 250


        -k       Kmer size for building de Bruijn graph. Default is 53


        -f       Insert size of the library. Default is 500


        -l       An integer specifying the read length. Default is 100


        -t       Number of threads used. Default is 8


        -h       Show help message


If you use option -i to supply your input, the actual file names should end with “_1.fq” and “_2.fq”. If you use option -I, the actual file names should be “_1.fq.gz” and “_2.fq.gz”. The prefix of the input file will be used as the prefix of the output files. 

Iqbal et al. implemented a very efficient way of storing the de Bruijn graph by a hash table in their CORTEX_VAR software. We adopted their approach for storing the de Bruijn graph and the -L and -W options specify the dimensions of the hash table. Note that the -L value is an exponential of 2 so when the value increases by 1 the height of the table doubles. The default value is sufficient for a high-coverage human whole genome sequencing dataset.

Please supply an odd number as K-mer size. 



##Output files

The main output of laSV is a list of detected structural variations in the VCF format with suffix “.SVs.vcf”. 

You may run “summarize_SVs.pl” in the scripts/ directory to convert the VCF file into BED format. The resulting file will have suffix “.SVs.summary”. The first three columns  of the file are the coordinate of the SV. The 4th column is the breakpoint ID and the 5th column is the inferred SV type. The 6th column indicates the length of the SV. The 7th column represents the estimated allele frequency of the SV and the number of the 8th column indicates the number of read pairs that support the SV. 

The FASTA file with suffix “.contigs.fa” contains all the maximal unambiguous contigs produced by the local assembly. 



##Example


We attached a test dataset generated by simulating SVs in the chr20 of the human genome (hg19). Users can run the following command to test if laSV is properly installed. 

cd /Your_path_to_laSV/laSV_v1.0.2/  
wget http://bib.umassmed.edu/~zhuangj/laSV_resources/test_dataset.tar.gz  
wget http://bib.umassmed.edu/~zhuangj/laSV_resources/hg19_reference.tar.gz  
tar -zxvf test_dataset.tar.gz  
tar -zxvf hg19_reference.tar.gz  
cd test_dataset  
bash ../scripts/run_laSV.sh -I chr20_sim_100_400 -D ../ -R ../reference_genomes/ -G hg19 -s 5 -k 53 -f 400 -l 100 -t 8
