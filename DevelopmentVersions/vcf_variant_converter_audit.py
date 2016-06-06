## Rebecca Haines
## 2016/06/01
## Mutalyzer Position Converter
## Version: 1.0
##
## Takes a vcf file, converts the variants to the correct format for the Mutalyzer
## batch query web service then queries the web service to annotate the variants
## with HGVS nomenclature.
##
## Use: python vcf_variant_converter_v1.py "VCF_FileName"

import sys
import vcf #PyVCF library for handling vcf files in python
from suds.client import Client
import base64
import datetime

VCF_FILE = sys.argv[1]


def audittrail(audit_file, audit_info):
    '''To write relevant audit trail information to a file.
    '''
    with open(audit_file, "a") as audit:
        audit.write(audit_info)
        audit.write("\n")


def vcftotxt(vcf_file, output_name, audit_file):
    '''Converts the variants in  VCF to a format suitable for Mutalyzer Batch
    Queries (Chr#:g.#GenomicCoord#REF>ALT, eg chr17:g.41223094T>C).
    Output is a text file containing the variants in a list.'''

    #check input file is a vcf
    vcf_filename = str(vcf_file) #converts filename to string
    assert ".vcf" in vcf_filename, "File chosen is not a vcf. Please select a file ending '.vcf'."

    #open the vcf using PyVCF methods
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))

    #Open the output file to write variants to
    with open(output_name, "w") as variants_txtfile:

        for record in vcf_reader:
            ALT = str(record.ALT[0]) #record.ALT is a list
            variant = record.CHROM+":g."+str(record.POS)+record.REF+">"+ALT
            variants_txtfile.write(variant)
            variants_txtfile.write("\n")

    #sanity check for user.
    print "Variants text file has been generated."
    # write file name to audit file
    audit = "Variants text file has been generated. File name %s" %output_name
    audittrail(audit_file, audit)
    
    converted_file = txtBase64(output_name) # nest this here to save a step
    return converted_file # you need to return the file handle to pass between functions


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
    
    # write Mutalyzer and HGVS version to audit file
    MutalyzerVersion = "Mutalyzer version: " + o.info().version
    HGVSVersion = "HGVS version: " + o.info().nomenclatureVersion
    audittrail(audit_file, MutalyzerVersion) 
    audittrail(audit_file, HGVSVersion)

    # This submits the batch job, providing the required arguments
    BatchJobID = o.submitBatchJob(MutalyzerInputfile, MutalyzerProcess, GenomeBuild)

    MutalyzerBatchID = "Mutalyzer Batch job ID: %s" %BatchJobID

    audittrail(audit_file, MutalyzerBatchID) #write batch job info to the audit file

    print "Job submitted to Mutalyzer"
    print "Batch job identifier: " + BatchJobID

    # monitor the progress of the batch job

    r = o.monitorBatchJob(BatchJobID)

    # Definitely inefficient but it works...
    while r >0:
        r = o.monitorBatchJob(BatchJobID)

    # get the result of the batch job once the job is complete (the result of the
    # monitor the job query == 0)
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
    PosConvAudit_file = vcf_file + "audit.txt"
    output_name = vcf_file + ".txt"
    
    # write version, start date and time and input file details to audit trail file
    software_version = "VCF variant converter v1.0"
    print_date_time = datetime.datetime.now().strftime("%y-%m-%d %H:%M")
    job_started = "Job started at: " + print_date_time
    input_vcf = "VCF file for analysis %s" %vcf_file
    audittrail(PosConvAudit_file, software_version)
    audittrail(PosConvAudit_file, job_started)
    audittrail(PosConvAudit_file, input_vcf)
    
    # convert the variants in the VCF to a text file for submission to Mutalyzer
    fileformutalyzer = vcftotxt(vcf_file, output_name, PosConvAudit_file)

    # Run the Mutalyzer query
    MutalyzerProcess = "PositionConverter"
    GenomeBuild = "GRCh37"
    MutalyzerResult = MutalyzerBatchSubmission(fileformutalyzer, MutalyzerProcess, GenomeBuild, PosConvAudit_file)

    # convert the Mutalyzer output to make it human readable
    NewVarFilename = str(vcf_file) + "_converted_variants.txt"
    decodeBase64(MutalyzerResult, NewVarFilename)

    # write time job finished to audit file
    job_finished = "Job finished at: " + print_date_time
    audittrail(PosConvAudit_file, job_finished)
   
    # Sanity check for the user
    print "Job complete."


# call the last function to run everything
MutalyzerPositionConverter(VCF_FILE)
