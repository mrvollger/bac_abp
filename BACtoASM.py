import os
import glob
from Bio import SeqIO
import re

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)
shell.executable("/bin/bash")
# set -x (exapnds vars)
# set -v (does not expand)
shell.prefix("source %s/env_BACtoASM.cfg; " % SNAKEMAKE_DIR)
#/net/eichler/vol2/home/mvollger/projects/bac_abp/env_BACtoASM.cfg

"""
Dependencies: should be taken care of by loading the env_PSV.cfg file
"""
# loocation of my additional scripts for assembly by phasing
base2="/net/eichler/vol2/home/mvollger/projects/abp"
# options for regular blasr
blasr_options ="-bestn 1 -minMapQV 30 -minAlignLength 2000 -minPctIdentity 70 -clipping subread "
blasr_options2="-bestn 10 -maxMatch 15 -minMapQV 0 -minAlignLength 1000 "
# options for the verison of blasr in pitchfork
pitchfork="source /net/eichler/vol2/home/mvollger/projects/builds/pitchfork/setup_pitchfork.sh && "
blasr_options_new="--bestn 1  --minAlnLength 2000 --minPctAccuracy 75 "
blasr_options_new2="--bestn 10 --maxMatch 15 --minAlignLength 500 "
#CANU_DIR="/net/eichler/vol5/home/mchaisso/software/canu/Linux-amd64/bin"
utils="/net/eichler/vol2/home/mvollger/projects/utility"
#CANU_DIR="/net/eichler/vol21/projects/bac_assembly/nobackups/canu/Linux-amd64/bin"
CANU_DIR="/net/eichler/vol2/home/mvollger/projects/builds/canu/Linux-amd64/bin"
vector="/net/eichler/vol4/home/jlhudd/projects/pacbio/vector/vector.fasta"



rule all:	
	input: 
		#bam="vectorTrimmed.subreads.bam",
		#plot=dynamic("{contigID}/depth.png"),
		#miropdf=dynamic("{contigID}/contig_compare.pdf"),
		#contigs=dynamic("{contigID}/ref.fasta"),
		#contigBam=dynamic("{contigID}/reads.sequal.bam"),
		#contigSortedBam=dynamic("{contigID}/reads.sorted.bam"),
		#depth=dynamic("{contigID}/depth.txt"),
		reads="reads.fasta", 
		gfaplot="assembly/asm.unitigs.gfa.png",
		miropdf=dynamic("{contigID}/contig_compare.pdf",)
		#asm='assembly/asm.unitigs.fasta',
	message: 'Running BACtoAsm'

# if the reads are form the RS2?, middle age tech
if(os.path.exists("bax.fofn")):
	rule baxToBam:
		input:
			bax="bax.fofn",
		output:
			subreads="movie.subreads.bam",
		shell:
			"""
			{pitchfork} bax2bam -f {input.bax} --subread -o movie
			rm movie.scraps.bam*
			"""
	rule inputFofn:
		input:
			subreads="movie.subreads.bam",
		output:
			fofn="input.fofn"
		shell:
			"""
			echo $(readlink -f {input.subreads} ) > {output.fofn}
			"""

elif( os.path.exists("old.bax.fofn") ):
	# if the reads are the really old reads
	rule mapVector:
		input:
			bax="old.bax.fofn",
		output:
			trimfofn="trim.fofn",
			vectorm1="vector.m1",
		threads:8
		shell:
			"""
			mkdir -p baxFiles 
			for readFile in $(cat {input.bax}); do
				echo $readFile
				cp $readFile baxFiles/.
			done
			ls $PWD/baxFiles/*.ba*.h5 > {output.trimfofn}	
			
			blasr {output.trimfofn} {vector} -nproc {threads} -out {output.vectorm1}
			"""

	rule removeVectorOld:
		input:
			vectorm1="vector.m1",
			trimfofn="trim.fofn",
		output:
			trimDone = "baxFiles/trimDone.txt",
		shell:
			"""
			for baxFile in $(cat {input.trimfofn}); do
				echo $baxFile
				python /net/eichler/vol2/home/mvollger/projects/bac_abp/TrimAdapterFromRegionTable.py $baxFile {input.vectorm1}
			done
			touch {output.trimDone}
			"""

	rule convertToFasta:
		input:
			trimDone = "baxFiles/trimDone.txt",
			trimfofn="trim.fofn",
		output:
			reads="reads.fasta"
		shell:
			"""
			> {output.reads}
			for baxFile in $(cat {input.trimfofn}); do 
				/net/eichler/vol5/home/mchaisso/projects/blasr/cpp/pbihdfutils/bin/pls2fasta \
						$baxFile fasta.temp -trimByRegion
				cat fasta.temp >> {output.reads}
				rm fasta.temp
			done
			"""

# else input fofn must exist and it must contain paths to bam files
# sequal machine

if(os.path.exists("orgin.fasta")):
	rule mapSubreadsToFasta:
		input:
			orgFasta="orgin.fasta",
			fofn="input.fofn"
		output:
			bam1="fofn_orgin.bam"
		threads: 8
		shell:
			# get the alignmetn then filter the sam to not include useless stuff
			"""
			{pitchfork} blasr {blasr_options_new} --nproc {threads} --bam --out {output.bam1} {input.fofn} {input.orgFasta}
			"""
else:
	print("no a batched pool of bacs")
	rule mapSubreadsToFasta:
		input:
			fofn="input.fofn"
		output:
			bam1="fofn_orgin.bam"
		shell:
			# get the alignmetn then filter the sam to not include useless stuff
			"""
			samtools merge -f {output.bam1} $(cat {input.fofn})
			"""


rule findVector:
    input:    
        bam1="fofn_orgin.bam",
    output:
        novector="noVectorReads.list",
    threads: 8
    shell:
        # get the alignmetn then filter the sam to not include useless stuff
        """
        {pitchfork} blasr {blasr_options_new} --nproc {threads} --bam --out containsVector.bam \
                --unaligned {output.novector} --noPrintUnalignedSeqs {input.bam1} {vector} 
        """

rule removeVector:
    input:
        novector="noVectorReads.list",
        bam1="fofn_orgin.bam",
    output:
        bam="vectorTrimmed.subreads.bam",
    shell:
        """
        {utils}/extractBams.py {input.bam1} {input.novector} --out {output.bam}
        """

rule downsample:
    input:
        bam="vectorTrimmed.subreads.bam",
        fasta="orgin.fasta",
    output:
        downbam="downsample/vectorTrimmed.subreads.bam",
    shell:
        """
        mkdir -p downsample 
        samtools view -h -b -s 0.10 {input.bam} > {output.downbam}
        """


# get fasta file to be input into canu 
# I am keeping this but actualy going to use fastq now
rule bamToFasta:
    input:
        bam="vectorTrimmed.subreads.bam",
    output:
        reads="reads.fasta", 
    shell:
        """
        samtools fasta {input.bam} > {output.reads}
        """

# get fastq file to be input into canu 
rule bamToFastq:
    input:
        bam="vectorTrimmed.subreads.bam",
        reads="reads.fasta", 
    output:
        readsq="reads.fastq", 
    shell:
        """
        bamtools convert -format fastq -in {input.bam} -out {output.readsq}
        #samtools fastq {input.bam} > {output.readsq}
        # I am not very confident that eaither of these are keeping the qual info
        """
# run the assembly
rule runAssembly:
    input: 
        readsq='reads.fasta'
    output: 
        asm='assembly/asm.unitigs.fasta',
        asm2='assembly/asm.contigs.fasta',
        gfa="assembly/asm.unitigs.gfa",
        gfa2="assembly/asm.contigs.gfa",
    threads: 16
    shell:
        '''
        #corMhapSensitivity=high corMinCoverage=0 
        #contigFilter="50 2000 0.75 0.75 50" 
        if [ -s {input} ]; then
            module load java/8u25 && {CANU_DIR}/canu -pacbio-raw {input.readsq} \
                useGrid=false gnuplotTested=true \
                maxThreads={threads} \
                stageDirectory=$TMPDIR   \
                -d assembly -p asm \
                genomeSize=200000 \
                corOutCoverage=100000 \
                correctedErrorRate=0.040 \
                MhapSensitivity=low \
                corMaxEvidenceErate=0.15 \
                obtOvlHashBits=23 \
                utgOvlHashBits=23 \
                minOverlapLength=1000 minReadLength=1000 \
                || ( >&2 echo " no real assembly" && \
                mkdir -p assembly && \
                > {output.asm} )
        else
            >&2 echo " no real assembly"
            mkdir -p assembly
            > {output.asm}
        fi
        '''

# plot GFA graph
rule plotGFA:
	input:
		gfa="assembly/asm.unitigs.gfa",
		gfa2="assembly/asm.contigs.gfa",
	output:
		gfaplot="assembly/asm.unitigs.gfa.png",
		gfaplot2="assembly/asm.contigs.gfa.png",
	shell:
		"""
		/net/eichler/vol2/home/mvollger/projects/GFA/plotGFA.py {input.gfa}
		/net/eichler/vol2/home/mvollger/projects/GFA/plotGFA.py {input.gfa2}
		"""

'''
rule mvAssembly:
	input:
		asm="assembly/asm.contigs.fasta",
	output:
		ref="ref.fasta",
	shell:
		"""
		cp {input.asm} {output.ref}
		"""

rule asmmiropeats:
    input:
        gfaplot="assembly/asm.unitigs.gfa.png",
        asm="assembly/asm.unitigs.fasta",
    output:
        asmpdf="miropeats.pdf",
    shell:
        """
        miropeats -s 1000 {input.asm} > contig_compare
        mv threshold1000 contig_compare.ps
        ps2pdf contig_compare.ps {output.asmpdf}
        """

# mask repetative sequences, this makes the depth profile possible to read
# masking is not helping, disabled
rule maskUnitigs:
    input:
        asmpdf="miropeats.pdf",
        asm="assembly/asm.unitigs.fasta"
    output:
        masked="assembly/asm.unitigs.masked.fasta"
    shell:
        """
        echo "Maksing is not helping, disabled"
        cp {input.asm} {output.masked}
        #{base2}/maskRef.py --ref {input.asm} --out {output.masked}
        """
'''

# split the asm
#import runCmd
from subprocess import call
rule splitFasta:
    input:
        masked='assembly/asm.contigs.fasta'
    output:
        contigs=dynamic("{contigID}/ref.fasta"),
    run:
        records=list(SeqIO.parse(input.masked, "fasta"))
        length = len(records)
        print("Num of contigs: {}".format(length))
        names=[] 
        for contig in records:
            name = contig.name
            names.append(name)
            #print(name)
            #runCmd.exe("mkdir -p " + name)
            call(["mkdir", "-p", name])
            SeqIO.write(contig, name + "/ref.fasta", "fasta")
        
        #storage.store("contigs", names)

rule miropeats:
    input:
        contigs="{contigID}/ref.fasta",
    output:
        miropdf="{contigID}/contig_compare.pdf",
    shell:
        """
        pushd {wildcards.contigID}
        miropeats -s 1000 -onlyintra ref.fasta > contig_compare
        mv threshold1000 contig_compare.ps
        ps2pdf contig_compare.ps contig_compare.pdf
        popd
        """

rule mapToContig:
    input:
        reads="reads.fasta",
        contigs="{contigID}/ref.fasta",
    output:
        contigBam="{contigID}/reads.sequal.bam",
    threads: 4
    shell:
        """
        {pitchfork} blasr {blasr_options_new2} --nproc {threads} --bam \
                --out {output.contigBam} \
                {input.reads} {input.contigs} 
        """

rule sortContigBam:
    input:
        contigBam="{contigID}/reads.sequal.bam",
    output:
        contigSortedBam="{contigID}/reads.sorted.bam",
    shell:
        """
        samtools sort {input.contigBam} > {output.contigSortedBam}
        """

rule depthProfiles:
    input:
        contigSortedBam="{contigID}/reads.sorted.bam",
    output:
        depth="{contigID}/depth.txt",
    shell:
        """
        samtools depth -a {input.contigSortedBam} > {output.depth} 
        """
    
rule plotDepth:
    input:
        depth="{contigID}/depth.txt",
    output:
        plot="{contigID}/depth.png",
    shell:
        """
        {utils}/plotDepth.py {input.depth} {output.plot} 
        """
    




