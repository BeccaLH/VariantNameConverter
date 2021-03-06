## Rebecca Haines
## 2016/06/01
## Mutalyzer Position Converter
## Version: 1.0
##
## Takes a vcf file, converts the variants to the correct format for the Mutalyzer
## batch query web service then queries the web service to annotate the variants
## with HGVS nomenclature.

import sys
import vcf #PyVCF library for handling vcf files in python
from suds.client import Client
import base64
import datetime

VCF_FILE = sys.argv[1]


print software_version
print "Job started at: " + print_date_time


def audittrail(audit_file, audit_info):
    '''To write relevant audit trail information to a file.
    '''
    with open(audit_file, "a") as audit:
        audit.write(audit_info)


# this function was DEFINTELY working 2015/05/31
def vcftotxt(vcf_file, output_name):
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
            ALT = str(record.ALT[0]) #for some reason record.ALT is a list not a string; probably as ALT can have multiple alleles?
            variant = record.CHROM+":g."+str(record.POS)+record.REF+">"+ALT
            variants_txtfile.write(variant)
            variants_txtfile.write("\n")

    #sanity check for user.
    print "Variants text file has been generated."
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

def MutalyzerBatchSubmission(MutalyzerInputfile, MutalyzerProcess, GenomeBuild):
    '''Submit the batch job to Mutalyzer
    '''
    URL = "https://mutalyzer.nl/services/?wsdl"
    # This is taken from the suds example from the mutalyzer website.
    # Not clear what it's for, but it seems necessary
    c = Client(URL, cache=None)
    o = c.service

    # This submits the batch job, providing the required arguments
    BatchJobID = o.submitBatchJob(MutalyzerInputfile, MutalyzerProcess, GenomeBuild)

    print "Job submitted to Mutalyzer"
    print "Batch job identifier: " + BatchJobID

    # monitor the progress of the batch job

    r = o.monitorBatchJob(BatchJobID)

    # Definitely inefficient but it works...
    while r >0:
        r = o.monitorBatchJob(BatchJobID)

    # get the result of the batch job once the job is complete (the result of the
    #   monitor the job query == 0)
    MutalyzerOutput = o.getBatchJob(BatchJobID)

    # again because functions need to return something, right?
    return MutalyzerOutput

def decodeBase64(MutalyzerOutput, NewVarFilename):
    '''decode a text file from Base64
    '''
    #need to decode the output result to make it human readable
    ResultDecode = base64.b64decode(MutalyzerOutput)

    #output the results file
    with open('%s.txt' %NewVarFilename, 'w') as newf:
        newf.write(ResultDecode)


def MutalyzerPositionConverter(vcf_file):
    '''call all of the functions to run the program
    '''
    PosConvAudit_file = vcf_file + "audit.txt"
    output_name = vcf_file + ".txt"
    software_version = "VCF variant converter v1.0 \n"
    print_date_time = datetime.datetime.now().strftime("%y-%m-%d %H:%M")
    job_started = "Job started at: " + print_date_time + "\n"
    
    audittrail(PosConvAudit_file, software_version)
    audittrail(PosConvAudit_file, job_started)
    
    fileformutalyzer = vcftotxt(vcf_file, output_name)

    MutalyzerProcess = "PositionConverter"
    GenomeBuild = "GRCh37"

    MutalyzerResult = MutalyzerBatchSubmission(fileformutalyzer, MutalyzerProcess, GenomeBuild)

    NewVarFilename = str(vcf_file) + "_converted_variants"

    decodeBase64(MutalyzerResult, NewVarFilename)

    # Sanity check for the user
    print "Job complete."

# call the last function to run everything
MutalyzerPositionConverter(VCF_FILE)
