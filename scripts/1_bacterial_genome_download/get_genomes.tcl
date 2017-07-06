#!/bin/sh
# the next line restarts using tclsh \
exec tclsh "$0" "$@"

# Script number: 			1.1
# File: 					1 of 2
# Prerequisite script(s):	N/A
# Prerequisite file(s):		Bacterial Accession number file (.txt)
# Description: 				Download bacterial genomes from EMBL

proc getSeqBatch {accFileName} {



	

	# read the accession numbers:
	if [catch {open $accFileName r} accFile] {
		puts stdout "Error: cannot open $accFileName for reading"
		return -1
	}
	set accessNumbers [split [read $accFile]]
	set num [llength $accessNumbers]
	close $accFile

	# pass the accession numbers to getSeq one at a time:
set j 1
	foreach accNo $accessNumbers {
getSeq $accNo $j
puts "$j of $num"
set j [expr $j +1]
}
	
}

proc getSeq {accNo j1} {


##################################################################
# PARAMETER DEFINITIONS START:

# name of the main web page:
#set url {http://www.ebi.ac.uk/Tools/dbfetch/emblfetch}

# name of the HTML document called by the query that starts the search 
# (get this from an inspection of the web page source in your browser:
#  <form method="post" action="$queryPage"> ):
#set queryPage {/cgi-bin/emblfetch}

# name of the input form variable for the queryPage
# (get this from an inspection of the web page source in your browser:
# <input type="text" name="$queryName" value="" size=30> ):
#set queryName {id}

set timeout 10000 	;# timeout in msec for web access



	# check for valid accession number:
	if  {[regexp {[A-Z]+.[0-9]+} $accNo]} { puts stdout "extracting $accNo"} else {
puts stdout "Invalid accession number <$accNo>"} 


	# load the http procedures:
	package require http

	# encode the query properly:

	#	set query [http::formatQuery db embl style raw $queryName $accNo]

		set url http://www.ebi.ac.uk/Tools/dbfetch/emblfetch?db=embl&id=$accNo&format=embl&style=raw&Retrieve=Retrieve

	# POST the query and get the result page:
set token [http::geturl $url]
set data1 [http::data $token]

if {[regexp {CDS} $data1]} {

			set out [open $accNo.embl w] 	
			puts $out $data1
			close $out
			
				}


http::cleanup $token
	

}

if {$argc ==0} {return -code error "I need a file name"}

if {$argc >0}  {
set accFileName [lindex $argv 0]
if {[file exists $accFileName]} {
	set t0 [clock seconds]
				getSeqBatch $accFileName
				puts "getSeq: done with processing $accFileName  in [expr [clock seconds] - $t0] sec."
					} else {return -code error "File $accFileName doesn't exist\n"}

}


if {![file exists bacterial_genomes]} {file mkdir bacterial_genomes}

set flist [glob *embl]

foreach f $flist {
file rename $f bacterial_genomes
}



