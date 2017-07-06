#!/bin/sh
# the next line restarts using tclsh \
exec tclsh "$0" "$@"


#read in the data file

set in [open Mesoplasma_florum_l1.ASM830v1.cds.all.fa r]
set data [read $in]
close $in


#break file into individual genes

regsub -all {>AAT} $data "£>AAT" data
set genes [split $data £]

set nucs [list A C T G]

#count the number of 2nd codons starting with each nucleotide
foreach n $nucs {
set cnt($n) 0

set cntstart($n) 0
}


foreach gn $genes  {

if {[regexp {>} $gn]} {

regexp {>.*?\n(.*?)$} $gn all seq


set seq [string trim $seq]

set 4th [string index $seq 3]

incr cnt($4th)
}
}


#find the total number of 4th site nucleotides

puts "Proportions of each base at the fourth site:"
set total 0.0

foreach n $nucs {
set  total [expr $total + $cnt($n)]
}


#calculate the proportion of 4th sites for each nucleotide
foreach n $nucs {

set prop($n) [expr ($cnt($n)*1.0)/($total*1.0)]

puts "\tProportion 4th site $n: $prop($n)"
puts ""
}


#calculate the total number of codons in the genome starting with each nucleotide
set totalcdn 0
set totalgenes 0
foreach gn $genes  {

if {[regexp {>} $gn]} {
incr totalgenes
regexp {>.*?\n(.*?)$} $gn all seq


set seq [string trim $seq]
regsub -all {[^A-Za-z]} $seq {} seq

for {set i 0} {$i < [expr [string length $seq]-3]} {incr i +3} {

set base [string index $seq $i]

if {[regexp {[^A|T|C|G]} $base]} {puts "odd base: $base"}

incr cntstart($base)
incr totalcdn

}
}
}



#...and calculate the proportion...
puts "Proportions of each base at the start of all codons in the genome:"
foreach n $nucs {

set propall($n) [expr ($cntstart($n)*1.0)/($totalcdn*1.0)]

puts "\tProportion all 1st codon positions $n: $propall($n)"
puts ""

}


#now calculateN4 ratio - proportion N at ste 4 / proportion of all codons starting N

foreach n $nucs {

set normalprop($n) [expr ($prop($n)*1.0)/($propall($n)*1.0)]

puts "Normalized proportion 4th site $n: $normalprop($n)"
puts ""


}

#report counts
puts "Total codons in genome: $totalcdn"
puts "Number of A starting codons at codon 2: $cnt(A)"
puts "Number of codon starting A: $cntstart(A)"
puts "\a"

set out [open results4thsite.xls w]
puts $out "Nucloetide\tProp4th\tProp_all\tN4"

foreach n $nucs {
puts $out "$n\t$prop($n)\t$propall($n)\t$normalprop($n)"

}

close $out
