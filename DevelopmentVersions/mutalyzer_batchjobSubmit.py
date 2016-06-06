## Mutalyzer batch job submission
##
## Written by: Rebecca Haines (contact rebecca.haines@nuh.nhs.uk)
## Date: 2015/07/20
## Version: 1.0
##
## Input files are txt format files produced by VCFtoTXTconverter.py
##
## batch job submit takes two arguments (or 3 for the position converter):
##     1. input file
##     2. process (NameChecker (default), SyntaxChecker, 
##        PositionConverter, SnpConverter)
##     3. genome build (only required for PositionConverter)
## add additional arugment to name output file
##     4. output file name 
##
## use: python mutalyzer_batchjobSubmit.py "inputfilename" "PositionChecker" "GRCh37" "Mutalyzer020715"

import sys
from suds.client import Client

URL = 'https://mutalyzer.nl/services/?wsdl'

#convert the input file to base64 encoding (required for submission to mutalyzer)
with open(sys.argv[1], "rb") as f:
    data = f.read()
    inputfile = data.encode("base64")

#checking correct number of input arguments
if len(sys.argv) < 5:
    print 'Please provide a input file, process, genome build and output filename'
    sys.exit(1)
    
#This is taken from the suds example from the mutalyzer website. Not clear what it's for, but it seems necessary
c = Client(URL, cache=None)
o = c.service

#This submits the batch job, providing the required arguments
BatchJobID = o.submitBatchJob(inputfile,sys.argv[2], sys.argv[3])

print "Job submitted"
print 'Batch job identifier: ' + BatchJobID

#monitor the progress of the batch job
r = o.monitorBatchJob(BatchJobID)

while r >0:
    r = o.monitorBatchJob(BatchJobID)

#get the result of the batch job once the job is complete (the result of the 
#   monitor the job query == 0)
result = o.getBatchJob(BatchJobID)

#need to decode the output result to make it human readable
import base64
result_decode = base64.b64decode(result)

#output the results file
with open('%s.txt' %sys.argv[4], 'w') as newf:
    newf.write(result_decode)

#Confirm that the job has completed
print "Job complete"