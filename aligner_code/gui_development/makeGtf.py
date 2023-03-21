#/usr/bin
import re

import re
import operator
import sys,string,random
import sys,os
from subprocess import Popen,PIPE,STDOUT
#from string import maketrans


comma=re.compile(",")
colon=re.compile(":")
equal=re.compile("=")
dash=re.compile("-")
under=re.compile("_")
s=re.compile("\s+")
accep=re.compile("ACCEPTABLE")
isdig=re.compile("\d")
isarrow=re.compile("^>")
isname=re.compile("NM_1")
notarrow=re.compile("[a-z]")
parens=re.compile("\)")
tab=re.compile("\t")
tmpFile = sys.argv[1]


count=0
udo="udo"
with open(tmpFile, 'r') as fileIn1:
  lineA = [ line[:-1] for line in fileIn1.readlines() ]

  
#fileIn1=open(tmpFile,"r")
lineB = list(set(lineA))


#Check for any spaces - abort if spaces found
for var in lineA:
  if(" " in var):
    sys.stderr.write("A space was detected. Please remove. \n")
    sys.exit()
#Check for duplicate gene names/sequences - abort if duplicates identified
if(len(lineA) != len(lineB)):
  sys.stderr.write("Duplicate detected. Please fix. \n")
  sys.exit()



name_dict={}
k=0
for j in range(0,len(lineA)):
  if(isarrow.findall(lineA[j])):
    namer= lineA[j][1:]
    lamer=namer+"\t"+"protien_coding"+"\t"+"exon"+"\t"+"1"+"\t"+"50"+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+'gene_id '+'"'+namer+'"'+';transcript_id '+'"'+namer+'";'
    print(lamer)


