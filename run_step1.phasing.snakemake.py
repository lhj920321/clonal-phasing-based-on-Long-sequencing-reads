###################################
#######phasing 
###################################


##sample params
CCS_Samps=config["CCS_Sample"]
Reference_samp=config["REF_Sample"]
##SNV params
iSNVminFreq=config["iSNVminFreq"]
iSNVmaxFreq=config["iSNVmaxFreq"]
##path params
SNV_path=config["SNV_path"]
refFasta_path=config["refFasta_path"]
bam_path=config["bam_path"]
software_Path=config["software_Path"]
OUT_path=config["OUT_path"]
##out params
minPerct_covered=config["minPerct_covered"]
minhapNumNeed_forOut=config["minhapNumNeed_forOut"]
minDepth_forHapOut=config["minDepth_forHapOut"]
print(OUT_path)



rule all:
	input:
	##read_LongReadbam_file
	#0 filt SNV posi
		expand("{OUT_path}/phasing_{samp}_2_{Ref_samp}--{iSNVminFreq}-{iSNVmaxFreq}/1.out-Filted_reliablePosi.txt",\
			samp={CCS_Samps},Ref_samp=Reference_samp,\
			SNV_path=SNV_path,OUT_path=OUT_path,software_Path=software_Path,iSNVminFreq=iSNVminFreq,iSNVmaxFreq=iSNVmaxFreq),
	#1 read bam -reads haplotype
	 	 expand("{OUT_path}/phasing_{samp}_2_{Ref_samp}--{iSNVminFreq}-{iSNVmaxFreq}/1.out-SNVtype_of_Filted_Cites.txt",\
	 		samp={CCS_Samps},Ref_samp=Reference_samp,\
	 		OUT_path=OUT_path,software_Path=software_Path,iSNVminFreq=iSNVminFreq,iSNVmaxFreq=iSNVmaxFreq),
	# #2 phasing 
	 	 expand("{OUT_path}/phasing_{samp}_2_{Ref_samp}--{iSNVminFreq}-{iSNVmaxFreq}/step2_{samp}_2_{Ref_samp}_phasing.log",\
	 		samp={CCS_Samps},Ref_samp=Reference_samp,\
	 		OUT_path=OUT_path,software_Path=software_Path,iSNVminFreq=iSNVminFreq,iSNVmaxFreq=iSNVmaxFreq),



#0,filter the SNV sites
rule Filt_Posi_inPerson:
	input:
		CCS_iSNVtable = expand("{SNV_path}/all.iSNV_with_SNP.pyResults.txt", SNV_path=SNV_path),
	output:
		out_reliableCite_CCS = expand("{OUT_path}/phasing_{samp}_2_{Ref_samp}--{iSNVminFreq}-{iSNVmaxFreq}/1.out-Filted_reliablePosi.txt",\
			OUT_path=OUT_path,software_Path=software_Path,iSNVminFreq=iSNVminFreq,iSNVmaxFreq=iSNVmaxFreq,samp={CCS_Samps},Ref_samp=Reference_samp)
	params:
		scriptPath=config["software_Path"],
		minFreq=config["iSNVminFreq"],
		maxFreq=config["iSNVmaxFreq"],
		Sample={CCS_Samps}
	shell:
		"python {params.scriptPath}/Step0__read_vcf_FiltPosi.py   \
		-C {input.CCS_iSNVtable}  -m   {params.minFreq} -M   {params.maxFreq} \
		-o {output.out_reliableCite_CCS} -s {params.Sample} "



#1
rule read_LongReadbam_file:
	input:
		#configF="/shared/liuhj/HP/process/phasing_CCS/config.txt",
		reliablePosiVcf=expand("{OUT_path}/phasing_{samp}_2_{Ref_samp}--{iSNVminFreq}-{iSNVmaxFreq}/1.out-Filted_reliablePosi.txt",\
			OUT_path=OUT_path,iSNVminFreq=iSNVminFreq,iSNVmaxFreq=iSNVmaxFreq,samp={CCS_Samps},Ref_samp=Reference_samp),
		BamF=expand("{bam_path}/{samp}_2_{Ref_samp}.sort.bam",bam_path=bam_path,samp={CCS_Samps},Ref_samp=Reference_samp),
		refSampGnm=expand("{refFasta_path}/{Ref_samp}.fasta",refFasta_path=refFasta_path,Ref_samp=Reference_samp),
	output:	
		SNV_type=expand("{OUT_path}/phasing_{samp}_2_{Ref_samp}--{iSNVminFreq}-{iSNVmaxFreq}/1.out-SNVtype_of_Filted_Cites.txt",\
			OUT_path=OUT_path,iSNVminFreq=iSNVminFreq,iSNVmaxFreq=iSNVmaxFreq,samp={CCS_Samps},Ref_samp=Reference_samp),
		SNV_Type_Stat=expand("{OUT_path}/phasing_{samp}_2_{Ref_samp}--{iSNVminFreq}-{iSNVmaxFreq}/1.out-SNVtype_of_Filted_Cites.Stat.txt",\
			OUT_path=OUT_path,iSNVminFreq=iSNVminFreq,iSNVmaxFreq=iSNVmaxFreq,samp={CCS_Samps},Ref_samp=Reference_samp)		
	params:
		scriptPath=config["software_Path"],
		outP=expand("{OUT_path}/phasing_{samp}_2_{Ref_samp}--{iSNVminFreq}-{iSNVmaxFreq}",\
			OUT_path=OUT_path,iSNVminFreq=iSNVminFreq,iSNVmaxFreq=iSNVmaxFreq,samp={CCS_Samps},Ref_samp=Reference_samp),
		ONT_PosiIndex_step="200",
	shell:
		"python {params.scriptPath}/Step1__read_ONTbam_f.py  -b {input.BamF} -v {input.reliablePosiVcf} -r {input.refSampGnm} -o {params.outP} -s {params.ONT_PosiIndex_step}"  


#2
rule phasing:
	input:
		reliablePosiVcf=expand("{OUT_path}/phasing_{samp}_2_{Ref_samp}--{iSNVminFreq}-{iSNVmaxFreq}/1.out-Filted_reliablePosi.txt",\
			OUT_path=OUT_path,iSNVminFreq=iSNVminFreq,iSNVmaxFreq=iSNVmaxFreq,samp={CCS_Samps},Ref_samp=Reference_samp),
		SNV_Type_Stat=expand("{OUT_path}/phasing_{samp}_2_{Ref_samp}--{iSNVminFreq}-{iSNVmaxFreq}/1.out-SNVtype_of_Filted_Cites.Stat.txt",\
			OUT_path=OUT_path,iSNVminFreq=iSNVminFreq,iSNVmaxFreq=iSNVmaxFreq,samp={CCS_Samps},Ref_samp=Reference_samp),
		refSampGnm=expand("{refFasta_path}/{Ref_samp}.fasta",refFasta_path=refFasta_path,Ref_samp=Reference_samp)
	output:
		phasingLog=temp(expand("{OUT_path}/phasing_{samp}_2_{Ref_samp}--{iSNVminFreq}-{iSNVmaxFreq}/step2_{samp}_2_{Ref_samp}_phasing.log",\
			OUT_path=OUT_path,iSNVminFreq=iSNVminFreq,iSNVmaxFreq=iSNVmaxFreq,samp={CCS_Samps},Ref_samp=Reference_samp)),
	params:
		scriptPath=config["software_Path"],
		outP=expand("{OUT_path}/phasing_{samp}_2_{Ref_samp}--{iSNVminFreq}-{iSNVmaxFreq}",\
			OUT_path=OUT_path,iSNVminFreq=iSNVminFreq,iSNVmaxFreq=iSNVmaxFreq,samp={CCS_Samps},Ref_samp=Reference_samp),
		startLociNum="1",
		minCovPerct=config["minPerct_covered"],
		minhapNum_forOut=config["minhapNumNeed_forOut"],
		minDepthOut=config["minDepth_forHapOut"],
		sample={CCS_Samps},
		refSamp=Reference_samp,
	shell:
		"python {params.scriptPath}/Step2_phasing.py  -f {input.reliablePosiVcf} -F {input.SNV_Type_Stat} \
		-o {params.outP}  -s {params.startLociNum} -m {params.minhapNum_forOut} -c {params.minCovPerct} \
		-l {output.phasingLog}  -S {params.sample} -D {params.minDepthOut}  -R {input.refSampGnm} -r {params.refSamp}"





