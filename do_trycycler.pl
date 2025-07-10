#!/usr/bin/perl

use strict ;
use warnings ;
use Bio::SearchIO;
use Bio::SeqIO;

my $dirname = ".";

opendir(DIR, $dirname) or die "can't opendir $dirname: $!";
while (defined(my $trimmed_reads_file1 = readdir(DIR))) {
    
    #if ($trimmed_reads_file1 =~ m/([\w\d]+?)_([\w\d]+)_([\w\d]+)_R1_(\d+)_val_1\.fq\.gz/) {
    if ($trimmed_reads_file1 =~ /(SRR\d+)_1_val_1\.fq\.gz/) {
	
	my $trimmed_reads_file2 = "${1}_2_val_2.fq.gz";
     	
	if (-s $trimmed_reads_file2) {
	    warn "Found $trimmed_reads_file1 and $trimmed_reads_file2\n";
	} else {
	    die "Failed to find $trimmed_reads_file2";
	}
		
	### Sub-sample the reads
	my $subsample_dir = "${1}.nanopore.subsamples";
	if (-d $subsample_dir) {
	    warn "$subsample_dir already exists\n";
	} else {
	    my $subsample_command = "trycycler subsample --genome_size 5m --reads *.filtlong.fastq.gz --out_dir $subsample_dir";
	    system( $subsample_command);
	}
	die "Subsampling failed\n" unless my @files = glob("$subsample_dir");

	### Step 1: Do Flye assembly on each subsample
	my $assemblies_dir = "assemblies";
        if (-d $assemblies_dir) {
            warn "$assemblies_dir/ already exists\n";
        } else {
            my $create_assemblies_dir_command = "mkdir $assemblies_dir";
            warn "$create_assemblies_dir_command\n";
            system($create_assemblies_dir_command);
        }
	
	foreach my $i ( '01', '02', '03', '04', '05', '06', '07', '08', '09', 10, 11, 12) {
	    my $sample_fastq_file = "sample_${i}.fastq";
	    my $flye_output_dir = "sample_${i}.flye";
	    my $flye_command = "flye --nano-hq $subsample_dir/$sample_fastq_file --out-dir $flye_output_dir --threads 16";
	    if (-d  $flye_output_dir and -s "$flye_output_dir/assembly.fasta") {
		warn " $flye_output_dir already exists\n";
	    } else {
		warn "$flye_command\n";
		system($flye_command);
	    }
	    my $assembly_file = "$flye_output_dir/assembly.fasta";
	    my $assembly_symlink = "${i}.assembly.fasta"; 
	    if (-l "$assemblies_dir/$assembly_symlink") {
		warn " $assembly_symlink already exists\n";
	    } else {
		if (-s $assembly_file) {
		    warn "OK, $assembly_file exists`n";
		    my $assembly_symlink_command = "cd $assemblies_dir && ln -s ../$assembly_file $assembly_symlink && cd -";
		    warn "$assembly_symlink_command\n";
		    system($assembly_symlink_command);
		}
		die "Assembly $assembly_file does not exist" unless -d  $flye_output_dir and -s "$flye_output_dir/assembly.fasta";;
	    }
	    
	}
		
	### Step 2: Clustering the contigs
	### Assemblies directory is $assemblies_dir
	### Sequence reads file is $combined_nanopore_file or *filtlong.fastq.gz

	my $trycycler_output_dir = "trycycler";
	if ( -d $trycycler_output_dir) {
	    warn "$trycycler_output_dir already exists, so will not perform clustering step\n";
	} else {
	    warn "$trycycler_output_dir does not exist";
	    my $trycycler_cluster_command = "trycycler cluster --assemblies $assemblies_dir/*.fasta --reads *filtlong.fastq.gz --out_dir $trycycler_output_dir";
	    warn "$trycycler_cluster_command\n";
	    system($trycycler_cluster_command);
	}

	### Step 3: Reconciling the contigs
	### Reconcile all the clusters in folder $trycycler_output_dir

	opendir(DH, $trycycler_output_dir) or die $!;
	my @files = readdir(DH);
	closedir(DH);
	
	foreach my $file (@files) {
	    if ($file =~ m/cluster_\d+$/ and -d "$trycycler_output_dir/$file") {
		my $reconcile_output_file = "$trycycler_output_dir/$file/2_all_seqs.fasta";
		if (-s $reconcile_output_file) {
		    warn "$reconcile_output_file already exists, so do not perform reconcile\n";
		} else {
		    warn "Will reconcile $trycycler_output_dir/$file";
		    my $command = "trycycler reconcile --max_add_seq 1600 --reads *filtlong.fastq.gz --cluster_dir $trycycler_output_dir/$file";
		    warn "$command\n";
		    system($command);
		    die "Failed Step 3 (Reconcile): $command" unless -s $reconcile_output_file;
		}
	    } else {
		warn "Ignoring $file\n";
	    }
	}
	
	### Step 4: Multiple sequence alignment
	### Will perform this step on file $reconcile_output_file

	foreach my $file (@files) {
            if ($file =~ m/cluster_\d+$/ and -d "$trycycler_output_dir/$file") {
		my $reconcile_output_file = "$trycycler_output_dir/$file/2_all_seqs.fasta";
		my $alignment_output_file = "$trycycler_output_dir/$file/3_msa.fasta";
		if (-s $reconcile_output_file) {
		    if (-s $alignment_output_file) {
			warn "$alignment_output_file already exists, so no need to perform alignment step\n";
		    } else {
			warn "perform alignment on $reconcile_output_file, which already exists\n";
			my $command = "trycycler msa --cluster_dir $trycycler_output_dir/$file";
			warn "$command\n";
			system($command);
			die "Failed Step 4 (Align): $command" unless -s $alignment_output_file;
		    }
		} else {
		    warn "$reconcile_output_file does not exist, so cannot perform alignment";
		}
            } else {
                warn "Ignoring $file\n";
            }
        }
    	
	### Step 5: Partitioning reads
	### Perform this once for the whole directory $trycycler_output_dir
	### Sequence reads file is $combined_nanopore_file or *filtlong.fastq.gz    

	if (my @files = glob("$trycycler_output_dir/cluster_*/4_reads.fastq")) {
	    warn "It looks like partitioning has been done already. File(s) exist(s): $trycycler_output_dir/cluster_*/4_reads.fastq\n";
	} else {
	    warn "$trycycler_output_dir/cluster_*/4_reads.fastq does not exist";
	    warn "performing partition step on $trycycler_output_dir/cluster_*\n";
	    my $command = "trycycler partition --reads  --cluster_dirs $trycycler_output_dir/cluster_*";
	    warn "$command\n";
	    system($command);
	    die "Failed Step 5 (Partition reads)" unless @files = glob("$trycycler_output_dir/cluster_*/4_reads.fastq");;
	}


	### Step 6: Generating a consensus

	foreach my $file (@files) {
	    
	    if ($file =~ m/cluster_\d+$/ and -d "$trycycler_output_dir/$file") {
		
		my $consensus_output_file = "$trycycler_output_dir/$file/7_final_consensus.fasta";
		if (-s $consensus_output_file) {
		    warn "$consensus_output_file already exists, so do not perform consensus step\n";
		} else {
		    my $reconcile_output_file = "$trycycler_output_dir/$file/2_all_seqs.fasta";
		    my $alignment_output_file = "$trycycler_output_dir/$file/3_msa.fasta";
		    my $partition_output_file = "$trycycler_output_dir/$file/4_reads.fastq";
		    if (-s $partition_output_file) {
			warn "perform consensus on $partition_output_file, which already exists\n";
			my $command = "trycycler consensus --cluster_dir $trycycler_output_dir/$file";
			warn "$command\n";
			system($command);
			die "Failed step 6 (Consensus)" unless -s $consensus_output_file;
                    }
		}

	    } else {
		warn "Ignoring $file\n";
	    }

	}

    	### Step 7: polishing

	### Long-read polishing with medaka
	
	foreach my $file (@files) {
	    
            if ($file =~ m/cluster_\d+$/ and -d "$trycycler_output_dir/$file") {
		
                my $medaka_output_file = "$trycycler_output_dir/$file/8_medaka.fasta";

		if (-s $medaka_output_file) {
		    warn "$medaka_output_file already exists, so do not need to do medaka polishing\n";
		} else {

		    if (-s $medaka_output_file) {
			warn "$medaka_output_file already exists, so do not perform medaka polishing step\n";
		    } else {
			my $medaka_output_file = "$trycycler_output_dir/$file/8_medaka.fasta";
			warn "perform medaka polishing on $trycycler_output_dir/$file, which already exists\n";
			my $command = "medaka_consensus -i $trycycler_output_dir/$file/4_reads.fastq -d $trycycler_output_dir/$file/7_final_consensus.fasta -o $trycycler_output_dir/$file/medaka -m r1041_e82_400bps_sup_variant_v4.3.0  -t 12; mv $trycycler_output_dir/$file/medaka/consensus.fasta $trycycler_output_dir/$file/8_medaka.fasta; rm $trycycler_output_dir/$file/*.mmi $trycycler_output_dir/$file/*.fai";
			warn "$command\n";
			system($command);
			die "Failed step 7 (Medaka polishing)" unless -s $medaka_output_file;
		    }
		}
		
            } else {
                warn "Ignoring $file\n";
            }
	    
        }
	
	### concatenate
	my $concat_command = "cat $trycycler_output_dir/cluster_*/8_medaka.fasta > $trycycler_output_dir/medaka.fasta";
        warn "$concat_command\n";
        system($concat_command);
	
	### Do BWA aligment ahead of polishing with short reads
	### Use files $trimmed_reads_file1 and $trimmed_reads_file2 against $trycycler_output_dir/medaka.fasta
	my $sam_file1 = "$trycycler_output_dir/alignments_1.sam";
	my $sam_file2 =	"$trycycler_output_dir/alignments_2.sam";
	if (-s $sam_file1 and -s $sam_file1) {
	    warn "$sam_file1 and $sam_file1 a;ready exist, so do not need to run bwa mem\n";
	} else {
	    ### BWA alignment
	    my $index_command = "bwa index $trycycler_output_dir/medaka.fasta";
	    warn "$index_command\n";
	    system($index_command);
	    my $bwa_command1 = "bwa mem -t 16 -a $trycycler_output_dir/medaka.fasta $trimmed_reads_file1 > $sam_file1";
	    warn "$bwa_command1\n";
	    system($bwa_command1);
	    my $bwa_command2 = "bwa mem -t 16 -a $trycycler_output_dir/medaka.fasta $trimmed_reads_file2 > $sam_file2";
	    warn "$bwa_command2\n";
	    system($bwa_command2);
	}

	### Short-read polishing with Polypolish                                                                                                                                                      
        ### Use files $trimmed_reads_file1 and $trimmed_reads_file2 against $trycycler_output_dir/medaka.fasta
	my $filtered_sam_file1 = "$trycycler_output_dir/filtered_1.sam";
	my $filtered_sam_file2 = "$trycycler_output_dir/filtered_2.sam";
	my $final_polished_assembly_file = "$trycycler_output_dir/medaka.polypolish.fasta";
	if (-s $filtered_sam_file1 and -s $filtered_sam_file2) {
	    warn "$filtered_sam_file1 and $filtered_sam_file2 already exist, so do not run polypolish filter";
	} else {
	    my $polypolish_filter_command = "polypolish filter --in1 $sam_file1 --in2 $sam_file2 --out1 $filtered_sam_file1 --out2 $filtered_sam_file2";
	    warn "$polypolish_filter_command\n";
	    system($polypolish_filter_command);
	}

	if (-s $final_polished_assembly_file) {
	    warn "$final_polished_assembly_file already exists, do do not perform polypolish polishing\n";
	} else {
	    my $polypolish_polish_command = "polypolish polish $trycycler_output_dir/medaka.fasta $filtered_sam_file1 $filtered_sam_file2 > $final_polished_assembly_file";
	    warn "$polypolish_polish_command\n";
	    system($polypolish_polish_command);
	    die "Failed step 7 (Polypolish polishing)" unless -s $final_polished_assembly_file;
	}


	warn "\nFinal assembly is at $final_polished_assembly_file\n\n";
    }
    
}
closedir(DIR);





	
