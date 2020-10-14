
import sys

vcf_infile = open(sys.argv[1])
vcf_outfile = open(sys.argv[2], "w")
min_depth = int(sys.argv[4])
min_AF = float(sys.argv[5])
annovar_outfile = open(sys.argv[3], "w")

header = True
for line in vcf_infile :
    if header:
        if line.find("#CHROM") != -1 :
            sample = line.strip().split("\t")[-1]
            header = False
        vcf_outfile.write(line)
        continue
    lline = line.split("\t")
    ref = lline[3]
    alt = lline[4]
    if len(ref) + len(alt) == 2 :
        continue
    filter = lline[6]
    if not (filter == "PASS" or filter == "clustered_events" or filter == "germline") :
        continue
    chrom = lline[0]
    if len(ref) > 1 :
        ref = ref[1:]
        alt = "-"
        start_pos = str(int(lline[1])+1)
        end_pos = str(int(start_pos) + len(ref) - 1)
    else :
        ref = "-"
        alt = alt[1:]
        start_pos = lline[1]
        end_pos = start_pos
    Format = lline[8].split(":")
    AD_i = 0
    AF_i = 0
    DP_i = 0
    F1R2_i = 0
    F2R1_i = 0
    i = 0
    for f in Format :
        if f == "AD" :
            AD_i = i
        if f == "AF" :
            AF_i = i
        if f == "DP" :
            DP_i = i
        if f == "F1R2" :
            F1R2_i = i
        if f == "F2R1" :
            F2R1_i = i
        i += 1
    Data = lline[9].split(":")
    vRatio = Data[AF_i]
    alleleFreq = Data[AD_i]
    readDepth = Data[DP_i]
    strandInfo = "-"
    F1R2 = Data[F1R2_i]
    F2R1 = Data[F2R1_i]
    if int(readDepth) < min_depth :
        continue
    if float(vRatio) < min_AF :
        continue
    vcf_outfile.write(line)
    if alt == "-" :
        annovar_outfile.write(str(chrom[3:]) + "\t" + str(start_pos) + "\t" + str(end_pos) + "\t" + str(ref) + "\t" + str(alt) + "\t" + "comments: sample=" + sample + " variantAlleleRatio=" + str(vRatio) + " alleleFreq=" + str(alleleFreq) + " readDepth=" + str(readDepth) + " Tumor_Del=+" + F1R2.split(",")[1] + "|-" + F2R1.split(",")[1] + " Tumor_var_plusAmplicons=- Tumor_var_minusAmplicons=- Tumor_ref_plusAmplicons=- Tumor_ref_minusAmplicons=-" + "\n")
    else :
        annovar_outfile.write(str(chrom[3:]) + "\t" + str(start_pos) + "\t" + str(end_pos) + "\t" + str(ref) + "\t" + str(alt) + "\t" + "comments: sample=" + sample + " variantAlleleRatio=" + str(vRatio) + " alleleFreq=" + str(alleleFreq) + " readDepth=" + str(readDepth) + " Tumor_Ins=+" + F1R2.split(",")[1] + "|-" + F2R1.split(",")[1] + " Tumor_var_plusAmplicons=- Tumor_var_minusAmplicons=- Tumor_ref_plusAmplicons=- Tumor_ref_minusAmplicons=-" + "\n")
    
vcf_outfile.close()
vcf_infile.close()
annovar_outfile.close()
