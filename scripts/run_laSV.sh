#!/bin/bash -x
# laSV (Local Assembly based Structural Variation discovery tools)
# version 1.0.1
# 2015-01-27
# 
# Jiali Zhuang(jiali.zhuang@umassmed.edu)
# Weng Lab
# Programs in Bioinformatics and Integrative Biology
# University of Massachusetts Medical School
#
# laSV is a free software whose distribution is governed by
# the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

#usage function
usage() {
echo -en "\e[1;36m"
cat <<EOF

usage: $0 -i PREFIX -D laSV_dir -o Output_dir -R Ref_dir -G Ref_genome -s Score -k kmer_size -L block_height -W block_width -f insert_size -l read_length -t Threads -h 

laSV is a software package for detecting structural variations
from whole-genome sequencing data. 

Please send questions, suggestions and bug reports to:
jiali.zhuang@umassmed.edu

Options:
        -i     Input FASTQ sequence file prefix with full path. The actual file name should be PREFIX_1.fq and PREFIX_2.fq
        -I     Input compressed FASTQ sequence file prefix with full path. The actual file name should be PREFIX_1.fq.gz and PREFIX_2.fq.gz
        -D     Path to laSV directory
        -o     Path to output directory. Default is the current directory
        -R     Path to directory where the reference genome information is placed
        -G     Reference genome name (e.g. hg19, mm10, dm3, etc.)
        -s     Read Phred score threshold for building de Bruijn graph. Default is 20
        -L     The height of the memory block allocated for storing the de Bruijn graph. Default is 24
        -W     The width of the memory block allocated for storing the de Bruijn graph. Default is 250
        -k     Kmer size for building de Bruijn graph. Default is 53
        -f     An integer specifying the length of the fragments (inserts) of the library. Default is 500
        -l     An integer specifying the read length. Default is 100
        -t     Number of threads used. Default is 8
        -h     Show help message

EOF
echo -en "\e[0m"
}

# taking options
while getopts "hi:t:f:I:o:R:D:l:k:L:W:s:p:G:" OPTION
do
        case $OPTION in
                h)
                        usage && exit 1
		;;
                i)
                        FQ=$OPTARG
		;;
	        I)
		        FZ=$OPTARG
		;;
	        f)
		        INSERT=$OPTARG
		;;
	        l)
		        READLEN=$OPTARG
		;;
	        k)
		        KMER=$OPTARG
		;;
	        L)
		        MLEN=$OPTARG
		;;
	        W)
		        MWID=$OPTARG
	        ;;
	        s)
		        SCORE=$OPTARG
		;;
	        p)
		        PRUNE=$OPTARG
		;;
                o)
                        OUTDIR=$OPTARG
                ;;
                t)
                        CPU=$OPTARG
                ;;
                D)
                        DIR=$OPTARG
                ;;
	        R)
		        REF=$OPTARG
		;;
	        G)
		        GENOME=$OPTARG
		;;
                ?)
                        usage && exit 1
                ;;
        esac
done


if [[ -z $FQ ]] && [[ -z $FZ ]]
then
        echo -e "\e[1;31mPlease set input files\e[0m"
        usage && exit 1
fi

if [[ -z $DIR ]]
then 
        echo -e "\e[1;31mPlease set the path to laSV directory\e[0m"
        usage && exit 1
fi

if [[ -z $REF ]]
then 
        echo -e "\e[1;31mPlease set the path to directory where reference information is placed[0m"
        usage && exit 1
fi

[ ! -z "${CPU##*[!0-9]*}" ] || CPU=8
[ ! -z "${INSERT##*[!0-9]*}" ] || INSERT=500
[ ! -z "${READLEN##*[!0-9]*}" ] || READLEN=100
[ ! -z "${MLEN##*[!0-9]*}" ] || MLEN=24
[ ! -z "${MWID##*[!0-9]*}" ] || MWID=250
[ ! -z "${SCORE##*[!0-9]*}" ] || SCORE=20
[ ! -z "${PRUNE##*[!0-9]*}" ] || PRUNE=2
[ ! -z "${KMER##*[!0-9]*}" ] || KMER=53
[ ! -z $OUTDIR ]  || OUTDIR=$PWD

mkdir -p "${OUTDIR}" || echo -e "\e[1;31mWarning: Cannot create directory ${OUTDIR}...\e[0m"
cd ${OUTDIR} || echo -e "\e[1;31mError: Cannot access directory ${OUTDIR}... Exiting...\e[0m" || exit 1
touch ${OUTDIR}/.writting_permission && rm -rf ${OUTDIR}/.writting_permission || echo -e "\e[1;31mError: Cannot write in directory ${OUTDIR}... Exiting...\e[0m" || exit 1

function checkExist {
        echo -ne "\e[1;32m\"${1}\" is using: \e[0m" && which "$1"
        [[ $? != 0 ]] && echo -e "\e[1;36mError: cannot find software/function ${1}! Please make sure that you have installed the pipeline correctly.\nExiting...\e[0m" && exit 1
}
echo -e "\e[1;35mTesting required softwares/scripts:\e[0m"
checkExist "echo"
checkExist "cat"
checkExist "rm"
checkExist "mkdir"
checkExist "date"
checkExist "mv"
checkExist "sort"
checkExist "touch"
checkExist "awk"
checkExist "grep"
checkExist "fastx_collapser"
checkExist "twoBitToFa"
checkExist "bwa"
checkExist "samtools"
checkExist "bedtools"
echo -e "\e[1;35mDone with testing required softwares/scripts, starting pipeline...\n\n\e[0m"

if [[ -z $FQ ]]
then
    name=`basename $FZ`
    i=$FZ"_1.fq.gz"
    j=$FZ"_2.fq.gz"
else
    name=`basename $FQ`
    i=$FQ"_1.fq"
    j=$FQ"_2.fq"
fi

echo $i > left.file
echo $j > right.file


#############################################
## Build deBruijn graph and print branches ##
#############################################
echo -e "\e[1;35mConstructing de Bruijn graph from reads and print branch sequences...\n\e[0m"
$DIR/bin/dB_graph --kmer_size $KMER --mode build --pe_list left.file,right.file --mem_height $MLEN --mem_width $MWID \
--quality_score_threshold $SCORE --remove_low_coverage_supernodes $PRUNE --dump_binary $name.ctx
mv branches.fa $name.branches.fa
echo -e "\e[1;35mBranch sequences printed\n\e[0m"

rm left.file right.file


###############################
## Map reads to the branches ##
###############################
echo -e "\e[1;35mDetecting connectivity between branches...\n\e[0m"
bwa index -a bwtsw $name.branches.fa
bwa mem -t $CPU -k $KMER -r 20 -w 0 -B 100 -O 100 -L 0 -U 0 -T $KMER -v 2 -a -S $name.branches.fa $i $j | perl $DIR/scripts/process.branching.mapping.sam.pl $name $KMER
echo -e "\e[1;35mConnectivity building completed\n\e[0m"

rm $name.branches.fa*


##############################
## Print contigs from graph ##
##############################
echo -e "\e[1;35mTravering graph and printing contigs...\n\e[0m"
$DIR/bin/dB_graph --kmer_size $KMER --mode assembly --multicolour_bin $name.ctx --read_span_file $name.read.span --mem_height $MLEN --mem_width $MWID \
--max_read_len $READLEN --insert_size $INSERT --min_contig_len $READLEN --min_sup 1
echo -e "\e[1;35mTraversing completed. Processing contigs...\n\e[0m"
perl $DIR/scripts/reformat.pl contigs.fa | fastx_collapser -o $name.contigs.fa
perl $DIR/scripts/change_contig_name.pl $name.contigs.fa
echo -e "\e[1;35mContigs printing completed\n\e[0m"

#rm contigs.fa $name.read.span $name.ctx


##################################
## Estimate SV allele frequency ##
##################################
echo -e "\e[1;35mAligning contigs to the reference genome...\n\e[0m"
bwa mem -t $CPU -T 30 -v 2 $REF/$GENOME.fa $name.contigs.fa > $name.aligned.sam
echo -e "\e[1;35mAlignment completed. Detecting putative SVs...\n\e[0m"
perl $DIR/scripts/process.alignment.sam.pl $name.aligned $REF/$GENOME\_rpmk.bed
perl $DIR/scripts/get_flanking_seqs.pl $name $REF/$GENOME.chromInfo $REF/$GENOME.2bit
perl $DIR/scripts/filter_breakpoints.pl $name $REF/$GENOME\_rpmk.bed

rm $name.aligned.sam
echo -e "\e[1;35mPutative SV sequences ready.\nMapping reads to the variant alleles...\n\e[0m"
bwa index -a is $name.breakpoints.fa
bwa aln -t $CPU -n 3 -l 100 -R 10000 $name.breakpoints.fa $i -f 1.sai
bwa aln -t $CPU -n 3 -l 100 -R 10000 $name.breakpoints.fa $j -f 2.sai
bwa sampe -P -n 1 $name.breakpoints.fa 1.sai 2.sai $i $j | samtools view -bSf 0x2 - > $name.breakpoints.bam
samtools sort -@ 8 -m 3G $name.breakpoints.bam $name.breakpoints.sorted
rm $name.breakpoints.bam 1.sai 2.sai $name.breakpoints.fa.*
samtools index $name.breakpoints.sorted.bam

echo -e "\e[1;35mMapping reads to the reference1 alleles...\n\e[0m"
bwa index -a is $name.reference1.fa
bwa aln -t $CPU -n 3 -l 100 -R 10000 $name.reference1.fa $i -f 1.sai
bwa aln -t $CPU -n 3 -l 100 -R 10000 $name.reference1.fa $j -f 2.sai
bwa sampe -P -n 1 $name.reference1.fa 1.sai 2.sai $i $j | samtools view -bSf 0x2 - > $name.reference1.bam
samtools sort -@ 8 -m 3G $name.reference1.bam $name.reference1.sorted
rm $name.reference1.bam 1.sai 2.sai $name.reference1.fa.*
samtools index $name.reference1.sorted.bam

echo -e "\e[1;35mMapping reads to the reference2 alleles...\n\e[0m"
bwa index -a is $name.reference2.fa
bwa aln -t $CPU -n 3 -l 100 -R 10000 $name.reference2.fa $i -f 1.sai
bwa aln -t $CPU -n 3 -l 100 -R 10000 $name.reference2.fa $j -f 2.sai
bwa sampe -P -n 1 $name.reference2.fa 1.sai 2.sai $i $j | samtools view -bSf 0x2 - > $name.reference2.bam
samtools sort -@ 8 -m 3G $name.reference2.bam $name.reference2.sorted
rm $name.reference2.bam 1.sai 2.sai $name.reference2.fa.*
samtools index $name.reference2.sorted.bam

echo -e "\e[1;35mEstimating allele frequency...\n\e[0m"
perl $DIR/scripts/estimate_breakpoint_frequency.pl $name $DIR/scripts
perl $DIR/scripts/consolidate.freq.pl $name.breakpoints.freq
mv $name.breakpoints.freq.consolidate $name.breakpoints.freq

echo -e "\e[1;35mFrequency estimation completed\nPrinting vcf file\e[0m"
perl $DIR/scripts/classify_SVs.pl $name $REF/$GENOME.2bit

rm $name.breakpoints.fa $name.reference1.fa $name.reference2.fa $name.*.sorted.bam*

echo -e "\e[1;35mlaSV successfully completed!\e[0m"
###################################
## The end of run_laSV.sh script ##
###################################
