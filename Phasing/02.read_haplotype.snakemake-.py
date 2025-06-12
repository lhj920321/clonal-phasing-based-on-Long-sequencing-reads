###################################
#######phasing 
###################################
configfile: "/home/try/Documents/projects/liuhj-Other/phasing-2024/P2-t.config.yaml" 


##sample params
person=config["person"]
CCS_Samps=config["CCS_Sample"]
Reference_samp=config["REF_Sample"]

##path params
software_Path=config["software_Path"]
iSNVtable_P=config["iSNVtable_P"]
refFasta_path=config["refFasta_path"]
bam_path=config["bam_path"]
OUT_path=config["OUT_path"]


rule all:
	input:
	##read_LongReadbam_file
	#1 read bam -read haplotype
	 	 expand("{OUT_path}/readHap.{samp}_2_{Ref_samp}/2.read-haplotype.Stat.txt",\
	 		samp={CCS_Samps},Ref_samp=Reference_samp,OUT_path=OUT_path),



#1
rule read_LongReadbam_file:
	input:
		reliablePosiVcf=expand("{iSNVtable_P}/{person}.reliable.locus.CCS.table.txt",\
			iSNVtable_P=iSNVtable_P,person=person),
		BamF=expand("{bam_path}/{samp}/{samp}_2_{Ref_samp}.sort.bam",bam_path=bam_path,samp={CCS_Samps},Ref_samp=Reference_samp),
		refSampGnm=expand("{refFasta_path}/{Ref_samp}.fasta",refFasta_path=refFasta_path,Ref_samp=Reference_samp),
	output:	
		SNV_type=expand("{OUT_path}/readHap.{samp}_2_{Ref_samp}/2.read-haplotype.txt",\
			OUT_path=OUT_path,samp={CCS_Samps},Ref_samp=Reference_samp),
		SNV_Type_Stat=expand("{OUT_path}/readHap.{samp}_2_{Ref_samp}/2.read-haplotype.Stat.txt",\
			OUT_path=OUT_path,samp={CCS_Samps},Ref_samp=Reference_samp)		
	params:
		scriptPath=config["software_Path"],
		outP=expand("{OUT_path}/readHap.{samp}_2_{Ref_samp}",\
			OUT_path=OUT_path,samp={CCS_Samps},Ref_samp=Reference_samp),
		ONT_PosiIndex_step="200",
	shell:
		"python {params.scriptPath}/function_read_ONTbam.py  -b {input.BamF} -v {input.reliablePosiVcf} \
		-r {input.refSampGnm} -o {params.outP} -s {params.ONT_PosiIndex_step}"  










