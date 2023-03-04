Parabricks Workshop
=====================

This course introduces `NVIDIA Clara Parabricks <https://www.nvidia.com/en-us/clara/genomics/>`_ for read alignment and variant calling to demonstrate the benefits of GPU acceleration.

Course objectives
-----------------

As DNA sequencing has decreased in price, experiments have often been designed to use deeper coverage or more samples. This means more data is being produced and more time is spent on the bioinformatics analysis. CPU-based tools are commonly used for both read alignment and variant calling. However, significant time can be saved when GPU-accelerated tools are used. This will be demonstrated with NVIDIA Parabricks with the following steps:

* DNA alignment with fq2bam
* Calling variants with haplotypecaller
* Potential downstream analyses

This tutorial will show you how to run our core alignment tool, FQ2BAM, which allows you to align a FASTQ file according to GATK best practices at blazing speeds. This includes the gold-standard alignment tool BWA-MEM with inbuilt co-ordinate sorting of the output file, and optionally application of base-quality-score-recalibration and marking of duplicate reads.

Introduction to Clara Parabricks
--------------------------------

Parabricks is a software suite for performing secondary analysis of next generation sequencing (NGS) DNA and RNA data. It delivers results at blazing fast speeds and low cost. Clara Parabricks can analyze 30x WGS (whole human genome) data in about 25 minutes, instead of 30 hours for other methods. Its output matches commonly used software, making it fairly simple to verify the accuracy of the output.

.. note::
**Bioinformatics Analysis Stages**
* **Primary:** Production of reads and quality scores (basecalling)
* **Secondary:** Read filtering, Alignment/Assembly, and variant calling
* **Tertiary:** Multi-sample processing, annotation and filtering of variants, association analysis, ...


DNA alignment with fq2bam
-------------------------

Data we'll be using
###################
