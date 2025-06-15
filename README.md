# clonal phasing based on Long sequencing reads

 

# 

# title

## summary

## Dependencies

- python 3.7
- perl
- R version 4.4.3

- RAxML
- bwa  
- samtools
- Spades

- 

## Repo Contents
### data: genome sequences and meta data of sample
  genomes: genomes of each sample(fasta format)
  public genomes:
     List of public strains for phylogenetic analysis 

### 1.genome_assembly: 
     The scripts for genome assembly and gene annotation

### Phylogenetic analysis
      
### 2.Alignment_to_representative_genome: 
     The scripts for genome assembly

### 3.iSNV_SNP_calling: 
   -- Step1:   Extract nucleotide composition at each site from mpileup file
             Usage: bash step1_mpileup2ntfreq.sh  $Numb_of_Thread  $path_of_mpileup
   -- Step2:  iSNV and SNP calling for samples of each person
            Usage: snakemake -s iSNV_calling-pipeline.py -p
            input file:  config file: List file of nucleotide composition files from multiple samples(One sample per line), file name : {person}.ntfreq_file_list.txt
            Modify the corresponding file paths in the script according to the specific project

### 4. Genome_phasing_CCS_reads
   -- Step1: filter.SNV.locus




## repeat region 

- cut the ref genome to 150bp simulation sequences:

  python  cut 

- map the  simulation reads to the other bacteria genomes

- 



## Phylogenetic analysis

- prepare SNPs sequences:

​     ` python XXX`

- Phylogenetic analysis was 

`raxmlHPC-PTHREADS -s $align_trim_file -n raxmltree_result -m GTRGAMMAI -f a -x 12345 -N 1000 -p 123456 -T 24 -k`



## Genome assembly

`genome assembly : Spade XXX`



## variations annotation

`snpeff XXX`

1. To convert iSNV table into vcf format(VCFv4.1) that can be recognized  by SnpEff software: `python ./scripts/Variantions_annotion/iSNVTable_2_vcf_allsample.py  -i  ./data/   -o $allsamplesVcfPath   -r   MN908947.3`
2. Annoting the mutation by using SnpEff software ,for example: `java -jar snpEff.jar ann  MN908947   $$allsamplesVcfPath/${sample}.vcf    >$outSnpeff/{sample}.snpeffAnno.vcf`

## Citation

If you use data, results or conclusion from this work, please cite:



## Acknowledgement
