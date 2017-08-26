import os
import glob
from Bio import SeqIO
import re

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)
shell.executable("/bin/bash")
shell.prefix("source %s/env_PSV.cfg; " % SNAKEMAKE_DIR)

"""
Dependencies: should be taken care of by loading the env_PSV.cfg file
"""

blasrDir = '~mchaisso/projects/blasr-repo/blasr'
base = '/net/eichler/vol5/home/mchaisso/projects/AssemblyByPhasing/scripts/abp'
base2="/net/eichler/vol21/projects/bac_assembly/nobackups/scripts"
blasr_options ="-bestn 1 -minMapQV 30 -minAlignLength 2000 -minPctIdentity 70 -clipping subread "
blasr_options2="-bestn 10 -maxMatch 15 -minMapQV 0 -minAlignLength 1000 "
#CANU_DIR="/net/eichler/vol5/home/mchaisso/software/canu/Linux-amd64/bin"
CANU_DIR="/net/eichler/vol21/projects/bac_assembly/nobackups/canu/Linux-amd64/bin"

rule all:	
    input: 'h5_mapped_orgin.bam',
        plot=dynamic("{contigID}/depth.png"),
        contigs=dynamic("{contigID}/ref.fasta"),
        contigSam=dynamic("{contigID}/reads.orig.sam"),
        contigBam=dynamic("{contigID}/reads.orig.bam"),
        depth=dynamic("{contigID}/depth.txt"),
        reads="reads.fasta", 
        asm='assembly/asm.unitigs.fasta',
    message: 'Running BACtoAsm'

rule mapBasToFasta:
    input:
        orgFasta="orgin.fasta",
        h5fofn="bax.fofn",
    output:
        sam1="h5_mapped_orgin.sam"
    threads: 8
    shell:
        # get the alignmetn then filter the sam to not include useless stuff
        """
        blasr {blasr_options} -nproc {threads} -sam -out temp.sam {input.h5fofn} {input.orgFasta}
        samtools view -h -F 4 temp.sam > {output.sam1}
        rm temp.sam
        """

rule samToBam:
    input:
        sam1="h5_mapped_orgin.sam",
    output:
        bam1="h5_mapped_orgin.bam",
        bai1="h5_mapped_orgin.bam.bai",
    threads:16
    shell:
        """
        samtools view -b -S {input.sam1} | samtools sort -T tmp -@ {threads} -m 4G -o {output.bam1} 
        samtools index {output.bam1}
        """

# get fasta file to be input into canu 
# I am keeping this but actualy going to use fastq now
rule bamToFasta:
    input:
        bam1="h5_mapped_orgin.bam",
    output:
       reads="reads.fasta", 
    shell:
        """
        samtools fasta {input.bam1} > {output.reads}
        """
# get fastq file to be input into canu 
rule bamToFastq:
    input:
        bam1="h5_mapped_orgin.bam",
    output:
        readsq="reads.fastq", 
    shell:
        """
        samtools fastq {input.bam1} > {output.readsq}
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
                gnuplotTested=true  corMhapSensitivity=high corMinCoverage=1 corOutCoverage=300 \
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
import runCmd
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
            runCmd.exe("mkdir -p " + name)
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
        {base2}/plotDepth.py {input.depth} {output.plot} 
        """
    




