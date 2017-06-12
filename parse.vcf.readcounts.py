from sys import argv
script, vcf, output =argv


'''
script that parses a VCF file to a file with only 012 genotype status and ref and alt read counts
ex input VCF field : T   A   0/1:41,13:41:99:0,99,1485
   output : 1:41:13
parses multliallelic sites to their own line
ex input VCF field:  T   A/AC    0/2:34,0,11:21:99:0,0,1291
   output:  T   A   0:34:0
            T   AC  1:34:11
for missing genotypes output "NA:0:0" to distinguish missing from hom. ref
'''

OUT=open(output, 'w')
v=open(vcf)

def parse_geno(vcf_field, alt_num=1):
    line=vcf_field.split(":")
    geno=line[0]
    if geno == "./.":
        format="NA:0:0"
        return format
    else:
        allele_count=geno.count(str(alt_num))
        ref_reads=line[1].split(',')[0]
        alt_reads=line[1].split(',')[alt_num]
        format="%s:%s:%s" % (allele_count, ref_reads, alt_reads)
        return format

for v_line in v:
# find and write sample file names to output file
    if "#CHROM" in v_line:
        sample_name=v_line.strip().split("\t")[9:]
        format="%s\t%s\t%s\t%s\t%s\t%s\n" % ("CHR","POS","REF","ALT","AC",'\t'.join(map(str,sample_name)))
        OUT.write(format)
	continue
    elif "#" in v_line:
	continue
    parsed=v_line.strip().split("\t")
    chr=parsed[0]
    pos=parsed[1]
    ref=parsed[3]
    alt=parsed[4]
    alt_list=alt.split(",")
    AC=parsed[7].split(";")[0].split("=")[1]
    AC_list=AC.split(",")
    samples=parsed[9:]

# If multilple alternate aleles, parse each alt allele separately
    if len(alt_list)>1:
        for i in range(len(alt_list)):
            alti=alt_list[i]
            ACi=AC_list[i]
            altnum=int(i)+1
# remove sites with GATK * notation (sites not present due to a nearby deletion in that sample)
            if alti == "<*:DEL>" or alti == "*":
                continue
# parse snps coded as pseudo indels at multiallelic sites
# ex: ATGC/TTGC to A/T
            if (len(ref )>1) & (len(ref) == len(alti)):
                if ref[1:] == alti[1:]:
                    refn=ref[0]
                    alti=alti[0]
                    parsed_fields=[parse_geno(j,alt_num=altnum) for j in samples]
                    format="%s\t%s\t%s\t%s\t%s\t%s\n" % (chr,pos,ref,alti,ACi,'\t'.join(map(str,parsed_fields)))
                    OUT.write(format)
            else:
                parsed_fields=[parse_geno(j,alt_num=altnum) for j in samples]        
                format="%s\t%s\t%s\t%s\t%s\t%s\n" % (chr,pos,ref,alti,ACi,'\t'.join(map(str,parsed_fields)))
                OUT.write(format)
    else:
        parsed_fields=[parse_geno(j) for j in samples] 
        format="%s\t%s\t%s\t%s\t%s\t%s\n" % (chr,pos,ref,alt,AC,'\t'.join(map(str,parsed_fields)))
        OUT.write(format)
