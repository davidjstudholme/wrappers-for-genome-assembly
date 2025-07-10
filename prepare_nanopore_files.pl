#!/usr/bin/perl

use strict ;
use warnings ;
use Bio::SearchIO;
use Bio::SeqIO;

my $dirname = ".";

opendir(DIR, $dirname) or die "can't opendir $dirname: $!";
while (defined(my $file = readdir(DIR))) {

    if ($file =~ m/(\d+)_([\w\d]+)_([\w\d]+)_R1_(\d+)\.fastq\.gz/) {

	warn "\n$2\n";

	### Look for nanopore reads data
	#my $nanopore_file1 = "/data/djs217/PROJECTS/22_project-11319/nanopore/P00143_2024_07_01_11319_run1/01_raw_reads//11319_${2}_PAW68016_rebasecalled_pass.fastq.gz";
	#my $nanopore_file2 = "/data/djs217/PROJECTS/22_project-11319/nanopore/P00143_2024_07_01_11319_run2/01_raw_reads//11319_${2}_PAW68016_rebasecalled_pass.fastq.gz";

	my $nanopore_file1 = "/data/djs217/PROJECTS/22_project-11319/nanopore/P00146_20240710_11319_run1/01_raw_reads//11319_${2}_PAW67997_rebasecalled_pass.fastq.gz";
        my $nanopore_file2 = "/data/djs217/PROJECTS/22_project-11319/nanopore/P00146_20240710_11319_run2/01_raw_reads//11319_${2}_PAW67997_rebasecalled_pass.fastq.gz";
	
	my $combined_nanopore_file = "${2}.ONT.combined.fastq.gz";
	
	### Make symlinks to the nanopore files - this will facilitate upload to SRA
	if (-s $nanopore_file1 > 1000  and -s $nanopore_file2 > 1000) {
	    system("ln -s $nanopore_file1 $2.nanopore_run-1.fastq.gz");
	    system("ln -s $nanopore_file2 $2.nanopore_run-2.fastq.gz");
	} else {
	    die "Missing files or empty files: $nanopore_file1 and $nanopore_file2";
	}
	    
	### Trim the nanopore reads
	my $trimmed_nanopore_file = "11319_${2}_PAW68016_rebasecalled_pass.runs1-and-2.filtlong.fastq.gz";
	if (-s $trimmed_nanopore_file > 1000) {
	    warn "$trimmed_nanopore_file already exists\n";
	} else {
	    ### Concatenate nanopore run files 1 and 2 into a single combined file
	    if (-s $combined_nanopore_file > 1000) {
		warn "$combined_nanopore_file already exists\n";
	    } else {
		my $combined_files_command = "zcat $nanopore_file1 $nanopore_file2 | gzip > $combined_nanopore_file";
		warn "$combined_files_command\n";
		system($combined_files_command);
	    }
	}
	
	### Perform the filtering with filtlong
	#my $filtlong_command = "filtlong --min_length 2000 --keep_percent 50 $combined_nanopore_file | gzip > $trimmed_nanopore_file";	
	#my $filtlong_command = "filtlong --min_length 5000 --keep_percent 90 $combined_nanopore_file | gzip > $trimmed_nanopore_file";
	#if (-s $trimmed_nanopore_file > 1000) {
	#    warn "$trimmed_nanopore_file already exists\n";
	#} else {
	#    warn "$trimmed_nanopore_file is missing or empty\n";
	#    warn "$filtlong_command\n";
	#    system($filtlong_command);
	#}
	
	    
	### Make symlink to nanopore reads file
	if (-e $trimmed_nanopore_file) {
	    system("ln -s $trimmed_nanopore_file .");
	} else {
	    warn "$trimmed_nanopore_file does not exist\n";
	}
	
    }	    
}

closedir(DIR);





	
