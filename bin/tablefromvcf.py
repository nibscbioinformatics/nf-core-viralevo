#This is a script to read possibly annotated VCF files from lofreq and from ivar and extract the variants from them into a nice joint table
#script relies on filenames which must be of the form
#annotated files:
#${sampleID}_${caller}_annotated.vcf
#unannotated ivar and lofreq files:
#${sampleID}_ivar.vcf or ${sampleprefix}_lofreq.vcf

import sys
import os

varcallsdir = sys.argv[1]
infiles = os.listdir(varcallsdir)
#infiles = sys.argv[1]
fileout = open(sys.argv[2], "w")

basicpassalt = int(sys.argv[3]) #100
basicpassproportion = float(sys.argv[4]) #0.01


#ivar vcf has lines like:
#NC_045512.2     45      .       G       A       .       FALSE   IVAR_DP=245;IVAR_GFF=NA;IVAR_REFAA=NA;IVAR_ALTAA=NA;ANN=A|intergenic_region|MODIFIER|CHR_START-ORF1ab|CHR_START-GU280_gp01|intergenic_region|CHR_START-GU280_gp01|||n.45G>A||||||       GT:PVAL:AQ:DP:AF        G/A:0.482283:37,32:244,1:0.00408163
#lofreq has lines like:
#NC_045512.2     241     .       C       T       49314.0 PASS    DP=5101;AF=0.993923;SB=46;DP4=0,21,2067,3003;ANN=T|intergenic_region|MODIFIER|CHR_START-ORF1ab|CHR_START-GU280_gp01|intergenic_region|CHR_START-GU280_gp01|||n.241C>T||||||

fileout.write("Sample,Caller,Region,Position,Ref,Alt,Ref_Reads,Alt_Reads,Proportion,Basic_Pass,Gene\n")

for infile in infiles:
    if ("_lofreq_annotated.vcf" in infile) or ("_lofreq.vcf" in infile):
        filein = open(varcallsdir+"/"+infile)
        caller = "lofreq"
        samplename = infile.replace("_lofreq_annotated.vcf","").replace("_lofreq.vcf","")
        vcfout = open(samplename+"_"+caller+"_filtered.vcf", "w")
        for line in filein:
            if line[0] == "#":
                vcfout.write(line)
                continue
            collect = line.rstrip().split("\t")
            chromosome = collect[0]
            position = collect[1]
            ref = collect[3]
            alt = collect[4]
            truevar = (collect[6]=="PASS")
            refdepth = int(collect[-1].split("DP4=")[1].split(";")[0].split(",")[0]) + int(collect[-1].split("DP4=")[1].split(";")[0].split(",")[1])
            altdepth = int(collect[-1].split("DP4=")[1].split(";")[0].split(",")[2]) + int(collect[-1].split("DP4=")[1].split(";")[0].split(",")[3])
            proportion = collect[7].split(";AF=")[1].split(";")[0]
            basicpass = (int(altdepth) >= basicpassalt) and (float(proportion) >= basicpassproportion) and (truevar)
            if basicpass:
                vcfout.write(line)
            if ("_lofreq_annotated.vcf" in infile):
              gene = collect[7].split(";ANN=")[1].split("|")[3]
            else:
              gene = "NA"
            fileout.write(",".join([samplename,caller,chromosome,position,ref,alt,str(refdepth),str(altdepth),str(proportion),str(basicpass),gene])+"\n")
        filein.close()
        vcfout.close()
    if ("_ivar_annotated.vcf" in infile) or ("_ivar.vcf" in infile):
        filein = open(varcallsdir+"/"+infile)
        caller = "ivar"
        samplename = infile.replace("_ivar_annotated.vcf","").replace("_ivar.vcf","")
        vcfout = open(samplename+"_"+caller+"_filtered.vcf", "w")
        for line in filein:
            if line[0] == "#":
                vcfout.write(line)
                continue
            if len(line.rstrip()) < 1:
                continue
            collect = line.rstrip().split("\t")
            chromosome = collect[0]
            position = collect[1]
            ref = collect[3]
            alt = collect[4]
            truevar = (collect[6]=="TRUE")
            refdepth = collect[-1].split(":")[3].split(",")[0]
            altdepth = collect[-1].split(":")[3].split(",")[1]
            proportion = collect[-1].split(":")[4]
            basicpass = (int(altdepth) >= basicpassalt) and (float(proportion) >= basicpassproportion) and (truevar)
            if basicpass:
                vcfout.write(line)
            if ("_ivar_annotated.vcf" in infile):
              gene = collect[7].split(";ANN=")[1].split("|")[3]
            else:
              gene = "NA"
            fileout.write(",".join([samplename,caller,chromosome,position,ref,alt,str(refdepth),str(altdepth),str(proportion),str(basicpass),gene])+"\n")
        filein.close()
        vcfout.close()

fileout.close()
