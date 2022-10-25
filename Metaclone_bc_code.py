# This program counts the number of occurrences of each specified barcode in a fastq file to quantify the relative abundances of different bacterial strains in a sample that was sequenced. 
# There is currently no correction for reading both strands or potentially reading reverse, each directionality is done separately but they are mostly on FWD.
# This program requires your .fastq.gz files to be unzipped, for that I used gunzip -dk *_paired.fastq.gz, should keep the original file and make an unzipped version in the same directory. 
# To get to the point of running this script you need to trim and unzip, get all your _paired_R*.fastq files in the same directory.
# For trimming, I used Trimmomatic trimmomatic-0.38 with the following parameters: 
# java -jar trimmomatic-0.38.jar PE -phred33 CROP:150 ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:150:30 MINLEN:100 
# Make a list of all the files to process with ls *.fastq | awk -F ".fastq" '{print $1}'
# This should give you the file names in a list without the extension. These will be used with the sys.arg[1] to call the basename of the file and use it across the script to avoid having to make one script per file. 
# To call the script do: python nameofscript filebasename
# e.g. python Barcode_Seq_Trimming_Mc35.py D0-35-1_S7_L001_R1_paired


#Import modules. 
import re
import os
import numpy as np
import pandas as pd
from StringIO import StringIO
from io import StringIO
#from itertools import islice
from Bio.Seq import Seq #Biopython module
from Bio import SeqIO
pd.options.mode.chained_assignment = None
import sys

#sys.argv[0]

csv = '/path_to_output/' + sys.argv[1] + ".csv"
fastqfile = '/path_to_trimmed_data/' + sys.argv[1] + ".fastq"


#Open the input file. 
out = open(csv, 'w')
out.seek(0,2)
out.write('Effector_label,'+'Barcode,'+'Directionality,'+'Read,'+'Matches,'+'\n')

inbarcode = open('/path_to_bc_file/', 'r')

#########################	R1

for line in inbarcode:
	barcodecount_read = 0
	#barcodecount_readrev = 0
	#barcodecount_readcom = 0
	#barcodecount_readrevcom = 0
	line = line.strip('\n')
	label = re.sub(r'(.+)\t.+', r'\1', line)
	#print(label)
	barcode = re.sub(r'.+\t(.+)', r'\1', line)
	#print(barcode)
	infastq = open(fastqfile, 'r')
	#read_number = (sum(1 for line in infastq))/4
	#print(read_number)
	for read in SeqIO.parse(infastq, "fastq"): #Loop through each read in your fastq file by parsing with SeqIO. 
		read = read.seq
		readcomplement = read.complement() #Get complement of read.
		readrevcomplement = read.reverse_complement() #Get reverse complement of read.
		readrev = readrevcomplement.complement() #Get reverse of read.
		if re.search(str(barcode), str(read)): #Search for barcode in the read.
			barcodecount_read = barcodecount_read + 1 #Count the number of reads that contain the barcode.
		#if re.search(str(barcode), str(readrev)):
		#	barcodecount_readrev = barcodecount_readrev + 1
		#if re.search(str(barcode), str(readcomplement)):
		#	barcodecount_readcom = barcodecount_readcom + 1
		#if re.search(str(barcode), str(readrevcomplement)):
		#	barcodecount_readrevcom = barcodecount_readrevcom + 1
	
	print(barcode)
	out.seek(0,2)
	out.write(label+ "," + barcode +','+'Forward' +","+sys.argv[1]+","+ str(barcodecount_read)+"\n")
	infastq.close()
