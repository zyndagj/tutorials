Parabricks Workshop
=====================

This course introduces `NVIDIA Clara™ Parabricks® <https://www.nvidia.com/en-us/clara/genomics/>`_ for read alignment and variant calling to demonstrate the benefits of GPU acceleration.

Course objectives
-----------------

As DNA sequencing has decreased in price, experiments have often been designed to use deeper coverage or more samples. This means more data is being produced and more time is spent on the bioinformatics analysis. CPU-based tools are commonly used for both read alignment and variant calling. However, significant time can be saved when GPU-accelerated tools are used. This will be demonstrated with NVIDIA Parabricks with the following steps:

* DNA alignment with fq2bam
* Calling variants with haplotypeCaller
* Potential downstream analyses

This tutorial will show you how to run our core alignment tool, FQ2BAM, which allows you to align a FASTQ file according to GATK best practices at blazing speeds. This includes the gold-standard alignment tool BWA-MEM with inbuilt co-ordinate sorting of the output file, and optionally application of base-quality-score-recalibration and marking of duplicate reads.

Requirements
------------

* NVIDIA GPU and driver greater than version 465.32.*
* Container runtime that supports NVIDIA GPUs and Docker containers:

  * `Docker <https://docs.docker.com/get-docker/>`_
  * `Singularity <https://sylabs.io/docs/>`_
  * `Apptainer <https://apptainer.org/>`_
  * `Enroot <https://github.com/NVIDIA/enroot>`_

* The following containers:

  * `nvcr.io/nvidia/clara/clara-parabricks:4.2.0-1 <https://catalog.ngc.nvidia.com/orgs/nvidia/teams/clara/containers/clara-parabricks>`_
  * `biocontainers/bwa:v0.7.15_cv4 <https://hub.docker.com/r/biocontainers/bwa/tags>`_
  * `broadinstitute/gatk:4.3.0.0 <https://hub.docker.com/r/broadinstitute/gatk/tags>`_
  * `biocontainers/seqtk:v1.3-1-deb_cv1 <https://hub.docker.com/r/biocontainers/seqtk/tags>`_
  * `gzynda/snpmatch:5.0.1 <https://hub.docker.com/r/gzynda/snpmatch/tags>`_
  * `mgibio/samtools-cwl:1.16.1 <https://hub.docker.com/r/mgibio/samtools-cwl/tags>`_

Data Used
-------------------

This tutorial will be using paired-end reads from the TDr-7 strain of *Arabidopsis thaliana* from the `1001 Genomes project <https://1001genomes.org/index.html>`_

* `TDr-7 <https://www.ebi.ac.uk/ena/browser/view/SRR519591>`_

The *A. thaliana* reference genome, `TAIR10 <https://www.arabidopsis.org/>`_, with "Chr" removed from chromosome names to match naming convention used by 1001 Genomes project.

.. code-block:: shell

    # Download reference
    curl https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas.gz | zcat | sed -e "s/^>Chr/>/" > TAIR10_chr_all.fasta
    # Index with samtools
    samtools faidx TAIR10_chr_all.fasta

The `SNPmatch database <https://figshare.com/articles/dataset/SNP_dataset_for_A_thaliana_1001_Genomes_project/5514403?backTo=/collections/_/3722743>`_ for the 1001 Genomes Project.

All downloads can be performed by running :download:`download_data.sh <assets/download_data.sh>`

Introduction to Clara Parabricks
--------------------------------

NVIDIA Clara™ Parabricks® is a software suite for the secondary analysis (alignment, variant calling) of DNA and RNA from short and long reads on GPU-accelerated hardware. On the most performant platforms, Clara Parabricks can analyze a typical whole human genome dataset in about 25 minutes, instead of 30 hours for other methods. It was designed to be easy to run while also matching the output from commonly used software, which means analyses can be reproduced even without GPUs.

.. image:: assets/analysis_steps.png

Parabricks v4 supports a variety of GPU-accelerated secondary analyses: alignment, preprocessing, variant calling, QC, and even some variant annotation (which lies in tertiary analysis). The chart below shows each pipeline supported in Parabricks v4.2.0-1, but take a look at `the documentation <https://docs.nvidia.com/clara/parabricks/4.2.0/documentation/tooldocs/standalonetools.html>`_ for more information.

.. image:: assets/pb_pipelines.png

Getting Clara Parabricks
########################

Clara Parabricks is available as a container on NGC. 

https://catalog.ngc.nvidia.com/orgs/nvidia/teams/clara/containers/clara-parabricks

This tutorial will be using ``v4.2.0-1`` of the container, and it can be pulled to your system as follows:

.. code-block:: shell

    # Docker
    docker pull nvcr.io/nvidia/clara/clara-parabricks:4.2.0-1

    # singularity
    singularity pull docker://nvcr.io/nvidia/clara/clara-parabricks:4.2.0-1

To streamline this workshop on `ACES <https://hprc.tamu.edu/aces/>`_, an `LMOD <https://lmod.readthedocs.io/en/latest/>`_ modulefile was created with helper variables and functions for calling singularity containers. Please load it with:

.. code-block:: shell

    module use /scratch/user/u.gz28467/workshop/modules/modulefiles/
    module load parabricks_workshop/4.2.0-1

.. note::

    Since there are many ways to invoke the parabricks container, all example commands will omit container runtime commands. If you are using apptainer/singularity, make sure to include the ``--nv`` flag to `mount NVIDIA libraries and devices <https://docs.sylabs.io/guides/3.11/user-guide/gpu.html>`_.

DNA alignment with fq2bam
-------------------------

Unless you're starting with pre-aligned reads in a ``.bam`` file, the first step to many bioinformatics pipelines is alignment. Parabricks has the `fq2bam <https://docs.nvidia.com/clara/parabricks/4.2.0/documentation/tooldocs/man_fq2bam.html#man-fq2bam>`_ pipeline for DNA (based on bwa mem) and the `rna_fq2bam <https://docs.nvidia.com/clara/parabricks/4.2.0/documentation/tooldocs/man_rna_fq2bam.html#man-rna-fq2bam>`_ pipeline for RNA (based on STAR). This tutorial uses DNA inputs and will focus on fq2bam.

When fq2bam is run, reads (compressed or not) are aligned by GPU-bwa mem, alignments are sorted by coordinate, duplicates are marked, and Base Quality Score Recalibration (BQSR) is performed, and a final ``.bam`` file is output.

.. image:: https://docscontent.nvidia.com/dims4/default/b2cc884/2147483647/strip/true/crop/1230x402+0+0/resize/2460x804!/format/webp/quality/90/?url=https%3A%2F%2Fk3-prod-nvidia-docs.s3.us-west-2.amazonaws.com%2Fbrightspot%2Fsphinx%2F0000018b-6753-d717-adef-ffffd61b0000%2Fclara%2Fparabricks%2F4.2.0%2F_images%2Ffq2bam.png

Depending on how you need to process your sample, fq2bam has a `lot of options <https://docs.nvidia.com/clara/parabricks/4.2.0/documentation/tooldocs/man_fq2bam.html#fq2bam-reference>`_. Since the data used in this tutorial is paired-end, we're going to be using:

.. code-block::

  --ref REF             Path to the reference file. (default: None)
  --in-fq [IN_FQ [IN_FQ ...]]
                        Path to the pair-ended FASTQ files followed by optional read groups with quotes
  --out-bam OUT_BAM     Path of a BAM/CRAM file after Marking Duplicates. (default: None)
  --num-gpus NUM_GPUS   Number of GPUs to use for a run. (default: number detected)

Indexing the reference for fq2bam
#################################

Parabricks does require the the reference genome be indexed with bwa before it can run fq2bam. This function was not ported to GPU, so you'll need to run this with the CPU version of the code.

.. code-block:: shell

    bwa index TAIR10_chr_all.fasta

.. note::

    For the sake of time, the reference genome we'll be using was already indexed.
    This will need to be done for any other genomes you use.

Running fq2bam on one GPU
##########################

By default, Parabricks will use all GPUs detected on your system. However, that can be controlled with the ``--num-gpus`` flag. Lets start with one and see how things run.

.. literalinclude:: assets/run_fq2bam.sh
    :caption: :download:`run_fq2bam.sh <assets/run_fq2bam.sh>`

Running fq2bam on two GPUs
##########################

Now run the same pipeline on two GPUs. If you only requested a single GPU through your scheduler, you may need to start a new job.

Is there a larger effect on certain phases?

Running the CPU equivalent
##########################

For every workflow in Clara Parabricks, an equivalent CPU workflow is provided to reproduce results. `fq2bam is no exception <https://docs.nvidia.com/clara/parabricks/4.2.0/documentation/tooldocs/man_fq2bam.html#compatible-cpu-based-bwa-mem-gatk4-commands>`_, and the below example takes those commands and wraps them in a bash function.

.. literalinclude:: assets/run_fq2bam_cpu.sh
    :caption: :download:`run_fq2bam_cpu.sh <assets/run_fq2bam_cpu.sh>`

If you compare the size of the ``.bam`` files created by Parabricks and the CPU pipeline, you'll notice that the CPU pipeline are larger.

.. code-block:: shell

    $ ls -lh *bam
    -rw-rw-r-- 1 u.gz28467 u.gz28467 1.6G Mar  7 01:48 TDr-7_10M_cpu_sorted_marked.bam
    -rw-rw-r-- 1 u.gz28467 u.gz28467 1.4G Mar  7 01:35 TDr-7_10M_pb.bam

This is due to compression levels. If you compare the bam files with samtools, they contain the same number of alignments.

.. code-block:: shell

    $ samtools flagstat TDr-7_10M_pb.bam
    20022096 + 0 in total (QC-passed reads + QC-failed reads)
    20000000 + 0 primary
    0 + 0 secondary
    22096 + 0 supplementary
    3608170 + 0 duplicates
    3608170 + 0 primary duplicates
    19033358 + 0 mapped (95.06% : N/A)
    19011262 + 0 primary mapped (95.06% : N/A)

    $ samtools flagstat TDr-7_10M_cpu_sorted_marked.bam 
    20022096 + 0 in total (QC-passed reads + QC-failed reads)
    20000000 + 0 primary
    0 + 0 secondary
    22096 + 0 supplementary
    3608170 + 0 duplicates
    3608170 + 0 primary duplicates
    19033358 + 0 mapped (95.06% : N/A)
    19011262 + 0 primary mapped (95.06% : N/A)

Optional Exercises
##########################

* Have you tried using a different number of GPUs?
* What was the speedup when using a GPU?
* What happens if your input reads are compressed?
* Are the reads the same if you view them?

Calling Variants with haplotypeCaller
-------------------------------------

The Parabricks haplotypeCaller is a re-implementation of the `GATK HaplotypeCaller <https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller>`_, which is used to call germline (inherited) SNPs and indels through the local re-assembly of haplotypes and identify genotypes.

Like humans, `A. thaliana is diploid <https://www.pnas.org/doi/abs/10.1073/pnas.92.24.10831>`_ and has two copies of each chromosome, so at any given location across the genome, all aligned bases are the same (homozygous), or about half of the reads have one base and half have another (heterozygous). The chloroplast (C) is haploid, similar to the sex chromosomes in humans, and cannot be heterozygous.

.. image:: https://docscontent.nvidia.com/sphinx/0000018b-6753-d717-adef-ffffd61b0000/clara/parabricks/4.2.0/_images/parabricks-web-graphics-1259949-r2-haplotypecaller.svg

Similar to fq2bam, the haplotypeCaller pipeline in Parabricks has `many options <https://docs.nvidia.com/clara/parabricks/4.2.0/documentation/tooldocs/man_haplotypecaller.html#specifying-haplotype-caller-options>`_, but we'll only need to use these:

  --ref REF             Path to the reference file. (default: None)
  --in-bam IN_BAM       Path to the input BAM/CRAM file for variant calling. The argument may also be a local folder containing several bams; each will be processed by 1 GPU in batch mode. (default: None)
  --out-variants OUT_VARIANTS
                        Path of the vcf/g.vcf/gvcf.gz file after variant calling. The argument may also be a local folder in batch mode. (default: None)
  --num-gpus NUM_GPUS   Number of GPUs to use for a run. (default: number detected)

Indexing the reference for haplotypeCaller
##########################################

Parabricks does require the the reference genome be indexed with ``gatk CreateSequenceDictionary`` and ``samtools faidx`` before it can run haplotypeCaller. This function was not ported to GPU, so you'll need to run this with the CPU version of the code.

.. code-block:: shell

    gatk CreateSequenceDictionary -R TAIR10_chr_all.fasta; samtools faidx TAIR10_chr_all.fasta

.. note::

    For the sake of time, the reference genome we'll be using was already indexed.
    This will need to be done for any other genomes you use.

Running haplotypeCaller
#######################

HaplotypeCaller takes the ``.bam`` files created by fq2bam as input, and calls variants on the aligned reads.

.. literalinclude:: assets/run_pb_haplo.sh
    :caption: :download:`run_pb_haplo.sh <assets/run_pb_haplo.sh>`

.. note::

    Alternatively, both fq2bam and haplotypeCaller can be run with the `germline pipeline <https://docs.nvidia.com/clara/parabricks/4.2.0/documentation/tooldocs/man_germline.html#man-germline>`_.

Running the CPU equivalent
##########################

The `CPU equivalent of haplotypeCaller <https://docs.nvidia.com/clara/parabricks/4.2.0/documentation/tooldocs/man_haplotypecaller.html#compatible-gatk4-command>`_ only requires a single call to GATK, but it's much more time intensive than the Parabricks version. Once again, the commands have been wrapped in a bash function for easy usage.

.. literalinclude:: assets/run_cpu_haplo.sh
    :caption: :download:`run_cpu_haplo.sh <assets/run_cpu_haplo.sh>`

Optional Exercises
##########################

* Run the germline pipeline
* How much is runtime affected by the number of GPUs?
* Try running `DeepVariant <https://docs.nvidia.com/clara/parabricks/4.2.0/documentation/tooldocs/man_deepvariant.html#man-deepvariant>`_
* Try `re-training DeepVariant <https://docs.nvidia.com/clara/parabricks/4.2.0/tutorials/dvtraining.html>`_

Genotyping sample
-------------------------

Now that we have a ``.vcf`` file of all variants, we're ready to perform some "tertiary" analyses. One common one is sample identification based on genotype. In humans, your specific alleles, or variants, could be used to determine your ancestry. In plants, you could determine what variety of a crop you sequenced.

For our example, we're going to verify the identity of the sample we've been working with. The 1001 Genomes project actually has a web portal called `AraGeno <http://arageno.gmi.oeaw.ac.at/>`_ for identifying samples based on the called SNPs, but we're going to run `SNPmatch <https://github.com/Gregor-Mendel-Institute/SNPmatch>`_ manually.

.. literalinclude:: assets/run_snpmatch.sh
    :caption: :download:`run_snpmatch.sh <assets/run_snpmatch.sh>`

    
This will then create a ``.matches.json`` file, which will have an accession ID. Find this ID in the `1001 Genomes accessions list <https://1001genomes.org/accessions.html>`_ and see if it matches TDr-7.

Optional Exercises
##########################

* Try genotyping another sample from `HERE <https://www.ebi.ac.uk/ena/browser/view/SRP012869>`_
* Run your own data

Next Steps
----------

Try the RNA-seq pipelines:

* `rna_fq2bam <https://docs.nvidia.com/clara/parabricks/4.2.0/documentation/tooldocs/man_rna_fq2bam.html#man-rna-fq2bam>`_
* `starfusion <https://docs.nvidia.com/clara/parabricks/4.2.0/documentation/tooldocs/man_starfusion.html#man-starfusion>`_

Register for NVIDIA Deep Learning Institutes:

* `Training DeepVariant Models using Parabricks* [DLIT52115] <https://www.nvidia.com/gtc/session-catalog/?tab.catalogallsessionstab=16566177511100015Kus&search=parabricks#/session/1669934478047001d6Ot>`_
* `Variant Calling on Whole Exome Data using Parabricks <https://www.nvidia.com/en-us/on-demand/session/gtcfall22-dlit41350/?playlistId=playList-c395267f-7c85-4a96-90bb-574392cbd162>`_

Other packages for genomics analyses:

* `Rapids_singlecell <https://github.com/scverse/rapids_singlecell>`_: A GPU-accelerated tool for scRNA analysis.
* `RAPIDS <https://rapids.ai/>`_: GPU accelerated data science
* `Metacache-GPU <https://github.com/muellan/metacache/blob/master/docs/gpu_version.md>`_: memory efficient, fast & precise taxnomomic classification system for metagenomic read mapping

`Register for GTC 2024 <https://www.nvidia.com/gtc/?ncid=GTC-NVGZYNDA>`_ and look out for genomics talks!