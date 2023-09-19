#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys,os,time,re

def install_package(package_name):
        package_name = package_name.replace("_", "-")  # 下载pip fake_useragent 包时  包名是:fake-useragent
        p = os.popen("pip list --format=columns")  # 获取所有包名 直接用 pip list 也可获取
        pip_list = p.read()  # 读取所有内容
        print(pip_list)
        if package_name in pip_list:
            print("{} had been installed".format(package_name))
            return True
        else:
            print("{}had not been installed! Starting to install this package.".format(package_name))
            p = os.popen("pip install {}".format(package_name))
            if "Success" in p.read():
                print("{} installed success!".format(package_name))
                return True if "Success" in p.read() else False


# 调用执行检测 如果没安装 则自动安装
install_package('optparse')


import sys,os,time,re 
from optparse import OptionParser
print("ooooooo")



def readVcf(VcfF):
    PosiFreqD,vcflD = {},{}
    for VcfFl in open(VcfF).readlines():
        if "#" not in VcfFl and VcfFl != "\n":
            Tags = VcfFl.split("\n")[0].split("\t")
            Cite = Tags[1]
            Freq = float(Tags[9].split(":")[5].split("%")[0])
            if Cite not in PosiFreqD and Freq >= minFreq * 100 and Freq < maxFreq * 100:
                PosiFreqD[Cite] = Freq
                vcflD[Cite] = VcfFl
    return PosiFreqD,vcflD               



def iSNV_SNP_tableRead(iSNV_SNP_table):
    iSNVSNPtablels = open(iSNV_SNP_table,'r').readlines()
    HeaderL = iSNVSNPtablels[0].split("\n")[0].split("\t")
    lie = 0
    lieD = {}
    for Header in HeaderL:
        lieD[Header] = lie
        lie += 1
    lsD = {}
    for iSNVSNPtablel in iSNVSNPtablels:
        if "#" not in iSNVSNPtablel and iSNVSNPtablel != "\n" :
            lL = iSNVSNPtablel.split("\n")[0].split("\t")
            Pos = lL[0]
            lsD[Pos] = lL
    return lieD,lsD



def samp_Freq_Dic(samples,lieD,lsD):
    samp_FreqDic = {}
    sampNA_PosiDic = {}
    #print(lieD)
    for sample in samples:
        Lie = lieD[sample]
        samp_FreqDic[sample] = {}
        sampNA_PosiDic[sample] = []
        for Posi in lsD:
            Freq = lsD[Posi][Lie]
            if Freq != "NA" :
                if Freq != "NO":                
                    Freq = float(lsD[Posi][Lie])
                    #if Freq >= 0.02:
                    samp_FreqDic[sample][Posi] = Freq
                else:
                    Freq = 0
                    samp_FreqDic[sample][Posi] = Freq
            else:
                sampNA_PosiDic[sample].append(Posi)
    return samp_FreqDic,sampNA_PosiDic



def Person_NA_Posi(samples,sampNA_PosiDic):
    personNAPosi = {}
    for sample in samples:
        sampleNAPosi = sampNA_PosiDic[sample]
        for posi in sampleNAPosi:
            if posi not in personNAPosi:
                personNAPosi[posi] = ''
    return personNAPosi



def PersoniSNVLoci(sampIDs,NGS_FreqDic,NGSlsD,NGSpersNAPosi,CCS_FreqDic,CCSlsD,CCSpersNAPosi,iSNVNum):
    PersiSNVLoci = {}    
    for posi in NGSlsD.keys():
        if posi in CCSlsD.keys() and posi not in NGSpersNAPosi and posi not in CCSpersNAPosi:
            iSNVcount = 0
            PosiFreqD = {}
            for sampID in sampIDs:
                NGS_Freq = NGS_FreqDic[sampID][posi]
                CCS_Freq = CCS_FreqDic[sampID][posi]
                PosiFreqD[sampID] = CCS_Freq
                if (NGS_Freq >= iSNVminFreq and NGS_Freq < iSNVmaxFreq) and (CCS_Freq >= iSNVminFreq and CCS_Freq < iSNVmaxFreq):
                    iSNVcount += 1
            if iSNVcount >= iSNVNum:
                PersiSNVLoci[posi] = PosiFreqD
    return  PersiSNVLoci



def OutPosi(Sample,sampIDs,PersiSNVLoci,out_reliableCite_CCS):
    Head = "\t".join(["#iSNVLoci"] + sampIDs)
    lsL = [Head]
    for iSNVLoci in PersiSNVLoci:
        lL = [iSNVLoci]
        for SampID in sampIDs:
            lL.append(str(PersiSNVLoci[iSNVLoci][SampID]))
            #print(lL)
        lsL.append("\t".join(lL))

    O = "\n".join(lsL)
    #print(SampID)
    for SampID in sampIDs:
        print(Sample)
        if Sample == SampID.split("_2_")[0] :
        #if Sample.split("phasing_")[1] == SampID.split("_2_")[0] :
            print(len(O))
            print("Filt reliable Posi Num :" + str(len(O)-1))
            OutFO = open(out_reliableCite_CCS,'a')
            OutFO.write(O)
            OutFO.close()
            print(out_reliableCite_CCS.split("/")[-2].split("phasing_")[1].split("_2_")[0])
            personID = out_reliableCite_CCS.split("/")[-2].split("phasing_")[1].split("_2_")[0].split("-")[0]    
            
            personoutF="/".join(out_reliableCite_CCS.split("/")[0:-2]) + "/" + personID + ".Filted_reliablePosi.txt"
            if (os.path.exists(personoutF)):
                os.remove(personoutF)
            personoutFO = open(personoutF,'a')
            personoutFO.write(O)
            personoutFO.close()



def main():
################################################################################################################################
# program
################################################################################################################################
    startTime = time.time()
    print("Filt reliable Posi")
    print("min Freq of iSNV :" + str(iSNVminFreq))
    print("max Freq of iSNV :" + str(iSNVmaxFreq))

    
    NGSlieD,NGSlsD = iSNV_SNP_tableRead(NGS_iSNVtable)
    CCSlieD,CCSlsD = iSNV_SNP_tableRead(CCS_iSNVtable)
    sampIDs = []
    for sampID in CCSlieD:
        if "#Posi" not in sampID and "snv/SNP" not in sampID and sampID not in ["P24-C_sC-j_2_P24-C_sC-t","P8-E-t_2_P8-E-x","P22-C_sL-t_2_P22-C_sL-d"]:  #
            sampIDs.append(sampID)
    NGS_FreqDic,NGSsampNA_PosiDic = samp_Freq_Dic(sampIDs,NGSlieD,NGSlsD)
    NGSpersNAPosi = Person_NA_Posi(sampIDs,NGSsampNA_PosiDic)

    CCS_FreqDic,CCSsampNA_PosiDic = samp_Freq_Dic(sampIDs,CCSlieD,CCSlsD)
    CCSpersNAPosi = Person_NA_Posi(sampIDs,CCSsampNA_PosiDic)



    PersiSNVLoci = PersoniSNVLoci(sampIDs,NGS_FreqDic,NGSlsD,NGSpersNAPosi,CCS_FreqDic,CCSlsD,CCSpersNAPosi,iSNVNum)
    print(out_reliableCite_CCS)
    OutPosi(Sample,sampIDs,PersiSNVLoci,out_reliableCite_CCS)


 
   



    endTime = time.time()
    sys.stdout.write("Total time taken: "+str(endTime-startTime)+" seconds\n")






################################################################################################################################
# Parameters
################################################################################################################################

if __name__ == "__main__":
    usage = "usage:python  %prog [option]"
    parser = OptionParser(usage=usage)

    parser.add_option("-C","--CCS_iSNVtable",
                        dest = "CCS_iSNVtable",
                        default = "",
                        metavar = "file",
                        help = "CCS iSNV table  [required]")
    parser.add_option("-m","--iSNVminFreq",
                        dest = "iSNVminFreq",
                        default = "",
                        metavar = "float",
                        help = "min Freq of iSNV [required]")
    parser.add_option("-M","--iSNVmaxFreq",
                        dest = "iSNVmaxFreq",
                        default = "",
                        metavar = "float",
                        help = "max Freq of iSNV [required]")

    parser.add_option("-o","--out_reliableCite_CCS",
                        dest = "out_reliableCite_CCS",
                        default = "",
                        metavar = "file",
                        help = "Output file. [required],name: 1.out-Filted_reliablePosi.txt")
    parser.add_option("-s","--Sample",
                        dest = "Sample",
                        default = "",
                        metavar = "str",
                        help = "Sample [required]")

    (options,args) = parser.parse_args()
    #NGS_iSNVtable       = os.path.abspath(options.NGS_iSNVtable)
    CCS_iSNVtable     = os.path.abspath(options.CCS_iSNVtable)
    iSNVminFreq = float(options.iSNVminFreq)
    iSNVmaxFreq      = float(options.iSNVmaxFreq)
    out_reliableCite_CCS  = os.path.abspath(options.out_reliableCite_CCS)
    Sample = options.Sample
    iSNVNum = 1
    NGS_iSNVtable = CCS_iSNVtable



################################################################################################################################
# Check output file
################################################################################################################################

    if (os.path.exists(out_reliableCite_CCS)):
        os.remove(out_reliableCite_CCS)


    main()








