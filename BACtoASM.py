import os
import glob
from Bio import SeqIO
import re

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)
shell.executable("/bin/bash")
# set -x (exapnds vars)
# set -v (does not expand)
shell.prefix("source %s/env_BACtoASM.cfg; set -x " % SNAKEMAKE_DIR)
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
#CANU_DIR="/net/eichler/vol5/home/mchaisso/software/canu/Linux-amd64/bin"
utils="/net/eichler/vol2/home/mvollger/projects/utility"
#CANU_DIR="/net/eichler/vol21/projects/bac_assembly/nobackups/canu/Linux-amd64/bin"
CANU_DIR="/net/eichler/vol2/home/mvollger/projects/builds/canu/Linux-amd64/bin"
vector="/net/eichler/vol4/home/jlhudd/projects/pacbio/vector/vector.fasta"



rule all:	
    input: 
        bam1="fofn_orgin.bam",
        plot=dynamic("{contigID}/depth.png"),
        contigs=dynamic("{contigID}/ref.fasta"),
        contigSam=dynamic("{contigID}/reads.orig.sam"),
        contigBam=dynamic("{contigID}/reads.orig.bam"),
        depth=dynamic("{contigID}/depth.txt"),
        reads="reads.fasta", 
        asm='assembly/asm.unitigs.fasta',
    message: 'Running BACtoAsm'


rule baxToBam:
    input:
        bax="bax.fofn",
    output:
        subreads="movie.subreads.bam",
    shell:
        """
        {pitchfork} bax2bam -f {input.bax} --subread -o movie 
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
        readsq='reads.fastq'
    output: 
        asm='assembly/asm.unitigs.fasta'
    threads: 16
    shell:
        '''
        if [ -s {input} ]; then
            module load java/8u25 && {CANU_DIR}/canu -pacbio-raw {input.readsq} \
                genomeSize=200000 -d assembly -p asm useGrid=false \
                stageDirectory=$TMPDIR gnuplotTested=true  \
                corMhapSensitivity=high corMinCoverage=0 corOutCoverage=500 \
                MhapSensitivity=low \
                corMaxEvidenceErate=0.15 \
                minOverlapLength=750 minReadLength=2000 \
                maxThreads={threads} cnsThreads={threads} ovlThreads={threads} mhapThreads={threads} \
                contigFilter="50 2000 0.75 0.75 50" \
                || ( >&2 echo " no real assembly" && \
                mkdir -p assembly && \
                > {output.asm} )

        else
            >&2 echo " no real assembly"
            mkdir -p assembly
            > {output.asm}
        fi
        '''

# mask repetative sequences, this makes the depth profile possible to read
rule maskUnitigs:
    input:
        asm='assembly/asm.unitigs.fasta'
    output:
        masked='assembly/asm.unitigs.masked.fasta'
    shell:
        """
        {base2}/maskRef.py --ref {input.asm} --out {output.masked}
        """


# split the asm
#import runCmd
from subprocess import call
rule splitFasta:
    input:
        masked='assembly/asm.unitigs.masked.fasta'
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



rule mapToContig:
    input:
        reads="reads.fasta",
        contigs="{contigID}/ref.fasta",
    output:
        contigSam="{contigID}/reads.orig.sam"
    threads: 4
    shell:
        """
        blasr {blasr_options2} -nproc {threads} -sam -out temp.{wildcards.contigID}.sam {input.reads} {input.contigs}
        samtools view -h -F 4 temp.{wildcards.contigID}.sam > {output.contigSam}
        rm temp.{wildcards.contigID}.sam
        """

rule ContigSamToBam:
    input:
        contigSam="{contigID}/reads.orig.sam"
    output:
        contigBam="{contigID}/reads.orig.bam",
        contigBai="{contigID}/reads.orig.bam.bai",
    threads: 4
    shell:
        """
        samtools view -b -S {input.contigSam} | samtools sort -T tmp -@ {threads} -m 4G -o {output.contigBam} 
        samtools index {output.contigBam}
        """

rule depthProfiles:
    input:
        contigBam="{contigID}/reads.orig.bam",
    output:
        depth="{contigID}/depth.txt",
    shell:
        """
        samtools depth -a {input.contigBam} > {output.depth} 
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
    




