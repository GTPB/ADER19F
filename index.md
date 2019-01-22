---
layout: page
title: ADER18S
tagline: Analysis of Differential Expression with RNAseq (Second course in 2018)
description: Introductory course that covers practical aspects of the analysis of differential gene expression by RNAseq
---
![](./pages/images/Entry_index_image.jpg)

## Course Description
This introductory course covers practical aspects of the **analysis of differential gene expression by RNAseq**, from planning the gathering of sequence data to the generation of tables of differentially expressed gene lists and visualization of results. For this edition of the course, we will also explore some specificities of **single-cell RNAseq data analysis**. Towards the end, we will cover some of the initial steps of secondary analysis, such as **functional enrichment** of the obtained gene lists. Participants will first start learning the concepts using small example datasets, and then will apply the learned concepts in the training room using real sized examples. At the end of the course, participants should be able to autonomously apply most of the learned methods to their own data. 


## Target Audience
Life Scientists who want to be able to use NGS data (RNAseq) to infer genes differentially expressed between different conditions. Computational researchers that wish to get acquainted with the concepts and methodologies used in RNAseq are also welcome.

---

## Course Documentation

**Note -** All the datasets used for this training course is available in the following button. You need to unzip this file and follow the instructions throughout the documentation.

[**Download ADER18S Datasets**](https://github.com/GTPB/ADER18S/archive/data.zip){: .btn} <sub><i>File Size: 482,6MB</i></sub>

<br/>

### Day 1

#### 1 - Plan your experiment using NGS technologies:
+ **1.1 -** [List possibilities and limitations of NGS sequencing technologies](pages/L01_L02.md)

      What choices do you have when sending your samples to the sequencing facility

+ **1.2 -** [Choose adequate sequencing for your biological question](pages/L01_L02.md#LO1)

      How do the sequencing choices influence the kind of questions you can answer


#### 2 - [List steps in the analysis of RNAseq differential expression experiments](pages/L01_L02.md#LO2) 
        What are the steps in RNAseq data analysis


#### 3 - Assess the general quality of the raw data from the sequencing facility
+ **3.1 -** [Interpret what are fastq files and what is their content](pages/L03.md)
		    
      What information is in fastq files, and how is it organized
    
+ **3.2 -** [Use software like FastQC to process fastq files and produce QC reports](pages/L03.md#LO3.2)
		    
      Detect low quality bases in the QC reports  
      Detect sequence bias and possible presence of adaptors and other contaminants


#### 4 - Do simple processing operations in the raw data to improve its quality
+ **4.1 -** [Use  trimmomatic to remove low quality bases from your reads](pages/L04.md)
		    
      Use trimmomatic to filter/trim low quality bases from your reads

+ **4.2 -** [Use trimmomatic to remove adaptors and other artefactual sequences from your reads](pages/L04.md#LO4.2)
		
      Remove adaptors (eg. illumina adaptors) and unwanted sequences (eg. polyA tails) from your reads
      Check results using FastQC on filtered data



#### 5 - Generate alignments of processed reads against a reference genome
+ **5.1 -** [What is a reference genome, versioning and where to obtain genomes](pages/L05.md)
		
      Are genomes constant?
      Obtain genome fasta from Ensembl
       
+ **5.2 -** [Alignment software: hisat2](pages/L05.md#LO5.2)
		  
      What are the requisites for using burrows-wheeler approaches?
      Prepare a reference genome to use with hisat2
      
+ **5.3 -** [Run an alignment: the SAM/BAM alignment format](pages/L05.md#LO5.3)
		
      Run hisat2 in an example dataset
      What is the SAM/BAM format

<br/>

### Day 2

#### 6 - Assess the general quality of the alignments and detect possible problems
+ **6.1 -** [What is a reference gene annotation, versioning and where to obtain](pages/L06.md)
		  
      What is the GFF/GTF format
      Obtain genome GTF from Ensembl

+ **6.2 -** [Visualizing alignments in IGV for single genes](pages/L06.md#LO6.2)

+ **6.3 -** [Use Qualimap to assess quality of alignments](pages/L06.md#LO6.3)
		
      Interpret general alignment statistics such as percentage of aligned reads
      Check the reports to assess RNA integrity and diversity


#### 7 - Generate tables of counts
+ **7.1 -** [The process of generating gene counts from genome aligments](pages/L07.md)
		
      What parameters we need to consider when counting

+ **7.2 -** [Use featurecounts to generate table of gene counts](pages/L07.md#LO7.2)
		  
      Interpret results from featurecounts
      
+ **7.3 -** [Using Salmon to generate counts only with the transcriptome](pages/L07.md#LO7.3)
		
      Interpret results from salmon

<br/>

### Day 3
#### 8 - Generate lists of differentially expressed genes, at least for a simple pairwise comparison
+ **8.1 -** [Execute a pairwise differential expression analysis](pages/L08.md)
		  
      Use Galaxy to produce differentially expressed genes with DESeq2

+ **8.2 -** [Interpretation and visualization of results](pages/L08.md#LO8.2)
		  
      PCA plots comparing all samples: detection of outliers, and batch effects
      Heatmaps and other plots

+ **8.3 -** [Use more complex settings than simple pairwise comparisons](pages/L08.md#LO8.3)
      
      Account for batch effects and paired data

+ **8.4 -** [Gain control over your analysis using R and Rstudio](pages/L08.md#LO8.4)
		  
      Use R in Rstudio to make a pairwise comparison using DESeq2 and edgeR
      Use edgeR to perform more complex analysis such as ANOVA-like all versus all comparisons

<br/>

### Day 4
#### 9 - Understand specificies of differential gene expression in single-cell RNAseq
+ **9.1 -** [Overview of Single Cell RNA-seq (scRNA-seq)](pages/L09.md)
		  
      Specificities of single-cell RNAseq, using the Chromium system as example
      Diferences in raw data preprocessing and counting

+ **9.2 -** [Generate a count matrix for a single-cell RNAseq dataset](pages/L09.md#LO9.2)
		
      Use Cell Ranger to preprocess a Chromium (10x Genomics) datasets
      Use Dropseq tools to obtain an UMI count matrix for a non-standard dataset

+ **9.3 -** [Identification and characterization of cell subpopulations in a UMI count matrix](pages/L09.md#LO9.3)
		
      Quality checking and filtering of the count table
      Interpreting PCA plots and dimensionality reduction
      Indentify genes that distinguish the different groups

<br/>

### Day 5 
#### 10 - Perform simple functional enrichment analysis and understand the concepts involved [(**Slides (pdf)**)](./assets/ADER18_-_Functional_Enrichment_Analysis.pdf)
+ **10.1 -** [How to extract meaning from a list of genes](pages/L10.md)

      What are functional annotations, what types exist, and where to get them

+ **10.2 -** [Understand the concept of functional enrichment analysis, and the statistics involved](pages/L10.md#LO10.2)
		
      What is enrichment analysis and how is it performed
      How to define sample and population sets
      Why do we need multiple test corrections

+ **10.3 -** [Interpret the results of functional enrichment analysis](pages/L10.md#LO10.3)
		  
      What can we get from enrichment analysis results
      Using functional enrichment analysis with real lists of genes

---

### [Learning Objectives and Course Pre-requisites](pages/objectives_prerequisites.md)

### [Instructors](pages/instructors.md)

---

The source for this course webpage is [on github](https://github.com/GTPB/ADER18S).
