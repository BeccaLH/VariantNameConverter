##VCF to Text variant converter
##
## Written by: Rebecca Haines (contact rebecca.haines@nuh.nhs.uk)
## Date: 2015/07/20
## Version: 1.0
##
##Convert variants in vcf files to a format that Mutalyzer can handle
##Variants in vcf files are not in a suitable format for Mutalyzer. 
##The variants need to be converted from columns to the format:
##Chr#:g.#GenomicCoord#REF>ALT, eg chr17:g.41223094T>C
##
##use: python VCFtoTXTconverter.py "input vcf name" "output txt file name"

import sys
import vcf #PyVCF library for handling vcf files in python

#open the vcf using the vcf.Reader method
vcf_reader = vcf.Reader(open("%s" %sys.argv[1], 'r')) 

#Open the output file
with open("%s" %sys.argv[2], 'w') as f:

#write the results to the output file
    for record in vcf_reader:
        ALT = str(record.ALT[0]) #for some reason record.ALT is a list not a string
        variant = record.CHROM+':g.'+str(record.POS)+record.REF+'>'+ALT
        f.write(variant)
        f.write('\n')

#Confirm the job is complete
print "Job complete"