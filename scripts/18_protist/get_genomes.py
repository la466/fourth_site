#!/usr/bin/python

# Script number:				18.2
# File:							1 of 2
# Prerequisite script(s):
# Prerequisite file(s):
# Description:					Download protists genomes from embl protists


import sys
import imp
import urllib2
import HTMLParser
import time
from datetime import timedelta
from ftplib import FTP
import shutil
import urllib2
from contextlib import closing
import gzip
import glob
import os.path
import os


#############
# Variables #
#############

url = "ftp://ftp.ensemblgenomes.org/pub/protists/release-36/fasta/"



#############
# Functions #
#############


# Get the script descriptions
def script_misc():
	print(script_name)
	print(script_description)


# Get a list of the accession numbers
def read_accessions(file_path):

	accession_list = []

	accession_file = open(file_path, 'U')
	for line in accession_file:
		line = line.replace("\n", "")
		accession_list.append(line)

	return accession_list, len(accession_list)


# Download the genomes from EMBL
def download_genomes(start_time):


	download_files= []


	ensembl = FTP('ftp.ensemblgenomes.org')
	ensembl.login()

	wd = '/pub/release-36/protists/fasta/'
	ensembl.cwd(wd)
	# ensembl.retrbinary('RETR ' + gbname, open(gbname, 'wb').write)

	files = []

	files = list_folders(ensembl)

	for f in files:

		# if len(download_files):
		path = wd + f
		ensembl.cwd(path)

		dirs = list_folders(ensembl)

		if 'cds' in dirs:
			# print(path)
			path = path + '/cds'
			ensembl.cwd(path)

			files = list_folders(ensembl)
			for file in files:
				if file.endswith('.gz'):
					download_files.append(path + '/' + file)

		else:
			for g in dirs:

				path = wd + f + '/' + g

				ensembl.cwd(path)

				dirs = list_folders(ensembl)

				if 'cds' in dirs:
					# print(path)
					path = path + '/cds'
					ensembl.cwd(path)

					files = list_folders(ensembl)

					for file in files:
						if file.endswith('.gz'):
							download_files.append(path + '/' + file)


	ensembl.quit()

	if os.path.exists('protists/'):
		shutil.rmtree('protists/')
		os.mkdir('protists/')
	else:
		os.mkdir('protists/')


	for path in download_files:

		file_url = 'ftp://ftp.ensemblgenomes.org' + path

		with closing(urllib2.urlopen(file_url)) as r:

			new_path = 'protists/' + path.split("/")[-1]
			with open(new_path, 'wb') as f:
				print 'copying %s' % file_url
				shutil.copyfileobj(r, f)


	print('extracting files')
	for src_name in glob.glob(os.path.join('protists/', '*.gz')):
	    base = os.path.basename(src_name)
	    dest_name = os.path.join('protists/', base[:-3])
	    with gzip.open(src_name, 'rb') as infile:
	        with open(dest_name, 'wb') as outfile:
	            for line in infile:
	                outfile.write(line)

		os.remove('protists/' + base)





def list_folders(ensembl):

	files = []

	try:
	    files = ensembl.nlst()
	except ftplib.error_perm, resp:
	    if str(resp) == "550 No files found":
	        print "No files in this directory"
	    else:
	        raise
	return(files)

###############
# Run script #
###############


def run():

	start_time = time.time()


	download_genomes(start_time)	# Download bacteria EMBL files



if __name__ == "__main__":
	run()
