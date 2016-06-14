## Rebecca Haines rebecca.haines@nuh.nhs.uk
## 2016/06/14
## Mutalyzer Position Converter
## Version: 1.0
##
## Takes a vcf file, converts the variants to the correct format for the 
## Mutalyzer batch query web service then queries the web service to annotate 
## the variants with HGVS nomenclature.
##
## Use: python vcf_variant_converter_v1.py "VCF_FileName"

import sys
import vcf #PyVCF library for handling vcf files in python
from suds.client import Client
import base64
import datetime
import time

VCF_FILE = sys.argv[1]

# make sure the input file exists before proceeding. 
# don't want to generate audit file etc if the input filename is incorrect.
try:
    fhand = open(sys.argv[1])
except:
    print "VCF file cannot be opened. Please check that the file exists and is in vcf format then try again."
    exit()


def audittrail(audit_file, audit_info):
    '''To write relevant audit trail information to a file.
    '''
    with open(audit_file, "a") as audit:
        audit.write(audit_info)
        audit.write("\n")


def vcftotxt(vcf_file, output_name, indel_filename, audit_file):
    '''Converts the variants in  VCF to a format suitable for Mutalyzer Batch
    Queries (Chr#:g.#GenomicCoord#REF>ALT, eg chr17:g.41223094T>C).
    Output is a text file containing the variants in a list.'''

    # check input file is a vcf
    vcf_filename = str(vcf_file) #converts filename to string
    assert ".vcf" in vcf_filename, "File chosen is not a vcf. Please select a file ending '.vcf'."

    # open the vcf using PyVCF methods
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))


    # Open the output file to write variants to
    with open(output_name, "w") as variants_txtfile:
        
        # the record is the full variant line in the vcf
        for record in vcf_reader:
            for ALT in record.ALT: # some records contain more than one alt
                if len(record.REF) >1 or len(ALT) >1: # program for SNVs only
                    # call the indelvariants function to write these variants to a different file.
                    indelvariants(indel_filename, audit_file, record.CHROM, record.POS, record.REF, ALT)
                else:
                    variant = record.CHROM+":g."+str(record.POS)+record.REF+">"+str(ALT)
                    # Some vcfs contain "chr" in the CHROM field but others do not.
                    # the Mutalyzer input format requires "chr" in the name.
                    if "chr" not in variant: 
                        variant = "chr" + variant
                variants_txtfile.write(variant) # write the variant to the output text file
                variants_txtfile.write("\n")

    # sanity check for user.
    print "Variants text file has been generated."
    # write file name to audit file
    audit = "Variants text file has been generated. File name %s" %output_name
    audittrail(audit_file, audit)
    
    # convert file to Base64 for input to mutalyzer (calling the next function)
    converted_file = txtBase64(output_name) 
    return converted_file 


def indelvariants(indelFilename,audit_file,chrom,pos,ref,alt):
    '''to handle insertion/deletion variants. 
    This output file will need further work before being suitable for Mutalyzer 
    submission because the input format required for indels is more complex.
    A file containing insertion/deletion variants in the form: chr#:g.###REF>ALT
    '''
    
    # Open the output file to write variants to
    with open(indelFilename, "a") as variants_indelfile:
        
        variant = chrom+":g."+str(pos)+ref+">"+str(alt)
                # Some vcfs contain "chr" in the CHROM field but others do not.
                # the Mutalyzer input format requires "chr" in the name.
                
        if "chr" not in variant: 
            variant = "chr" + variant
        
        variants_indelfile.write(variant) # write the variant to the output text file
        variants_indelfile.write("\n")


def txtBase64(output_text):
    '''Open a text file and convert to Base64 encoding (required for Mutalyzer).
    '''
    with open(output_text, "rb") as textfile:
        Base64File = textfile.read()
        MutalyzerInputfile = Base64File.encode("base64")

    # functions need to return something I thought- this should be holding
    # the encoded file in memory ready for use in the next function
    return MutalyzerInputfile


def MutalyzerBatchSubmission(MutalyzerInputfile, MutalyzerProcess, GenomeBuild, audit_file):
    '''Submit the batch job to Mutalyzer
    '''
    URL = "https://mutalyzer.nl/services/?wsdl"
    # This is taken from the suds example from the mutalyzer website.
    # Not clear what it's for, but it seems necessary
    c = Client(URL, cache=None)
    o = c.service
    
    # write Mutalyzer version to audit file
    softwareversion = o.info().version
    MutalyzerVersion = "Mutalyzer version: " + softwareversion
    audittrail(audit_file, MutalyzerVersion) 
    
    # check Mutalyzer version is as expected. If an update has occured the 
    # program must be re-validated
    if softwareversion != "2.0.20.dev":
        warning = "Mutalyzer version has been updated. Contact the bioinformatician. \
Position Converter program has not been run."
        print warning
        audittrail(audit_file,warning)
        exit()
    
    # add HGVS version info to the audit trail file   
    HGVSVersion = "HGVS version: " + o.info().nomenclatureVersion
    audittrail(audit_file, HGVSVersion)

    # This submits the batch job, providing the required arguments
    BatchJobID = o.submitBatchJob(MutalyzerInputfile, MutalyzerProcess, GenomeBuild)
    # print type(BatchJobID)  # uncomment this line to investigate issues with 
    # the webservice. Type should be Class.

    MutalyzerBatchID = "Mutalyzer Batch job ID: %s" %BatchJobID
    
    #write batch job info to the audit file
    audittrail(audit_file, MutalyzerBatchID) 

    print "Job submitted to Mutalyzer"
    print MutalyzerBatchID

    # monitor the progress of the batch job

    r = o.monitorBatchJob(BatchJobID)

    # checking whether the batch job is complete.
    while r >0:
        r = o.monitorBatchJob(BatchJobID)
        time.sleep(10) # wait for 10 seconds before checking status of job again

    # get the result of the batch job once the job is complete (the result of 
    # the monitor the job query == 0)
    MutalyzerOutput = o.getBatchJob(BatchJobID)

    return MutalyzerOutput


def decodeBase64(MutalyzerOutput, NewVarFilename):
    '''decode a text file from Base64
    '''
    #need to decode the output result to make it human readable
    ResultDecode = base64.b64decode(MutalyzerOutput)

    #output the results file
    with open(NewVarFilename, 'w') as newf:
        newf.write(ResultDecode)


def MutalyzerPositionConverter(vcf_file):
    '''call all of the functions to run the program
    '''
    # output file names generated based on name of input VCF
    PosConvAudit_file = vcf_file + "_audit.txt"
    output_name = vcf_file + ".txt"
    indel_var_filename = vcf_file + "_indel_vars.text"
    
    # write version, start date and time and input file details to audit trail file
    software_version = "VCF variant converter v1.0"
    print_start_time = datetime.datetime.now().strftime("%y-%m-%d %H:%M")
    job_started = "Job started at: " + print_start_time
    input_vcf = "VCF file for analysis %s" %vcf_file
    audittrail(PosConvAudit_file, software_version)
    audittrail(PosConvAudit_file, job_started)
    audittrail(PosConvAudit_file, input_vcf)
    
    # convert the variants in the VCF to a text file for submission to Mutalyzer
    fileformutalyzer = vcftotxt(vcf_file, output_name, indel_var_filename, PosConvAudit_file)

    # Run the Mutalyzer query
    MutalyzerProcess = "PositionConverter"
    GenomeBuild = "GRCh37"
    MutalyzerResult = MutalyzerBatchSubmission(fileformutalyzer, MutalyzerProcess, GenomeBuild, PosConvAudit_file)

    # convert the Mutalyzer output to make it human readable
    NewVarFilename = str(vcf_file) + "_converted_variants.txt"
    decodeBase64(MutalyzerResult, NewVarFilename)

    # write time job finished to audit file
    print_stop_time = datetime.datetime.now().strftime("%y-%m-%d %H:%M")
    job_finished = "Job finished at: " + print_stop_time
    audittrail(PosConvAudit_file, job_finished)
   
    # Sanity check for the user
    print "Job complete."


# call the last function to run everything
MutalyzerPositionConverter(VCF_FILE)
