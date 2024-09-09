#!/usr/bin/env perl

# Description: This script converts bedGraph to fixedStep wig format with defined step size. Input file may be compressed as .gz.
# Coordinates in bedGraph input are assumed to be 0-based (http://genome.ucsc.edu/goldenPath/help/bedgraph.html).
# Coordinates in wig output are 1-based (http://genome.ucsc.edu/goldenPath/help/wiggle.html). 

# Usage: bedgraph_to_wig.pl --bedgraph input.bedgraph --wig output.wig --step step_size [--compact]
# --bedgraph : specify input file in bedGraph format.
# --wig : specify output file in fixedStep format.
# --step : specify step size. Note that span is set to be identical to step.
# --compact : if selected, steps with value equal to 0 will not be printed. This saves space but was not allowed in original wig format, thus some scripts using wig file as input may not understand it.

# Credits: This script was written by Sebastien Vigneau (sebastien.vigneau@gmail.com) in Alexander Gimelbrant lab (Dana-Farber Cancer Institute). The inspiration for this script comes from Dave Tang's own version (http://davetang.org/wiki/tiki-index.php?page=wig).


use strict;
use warnings;
use Getopt::Long;
use List::Util qw[min max];

my $usage = "Usage: $0 --bedgraph <infile.bedgraph> --wig <outfile.wig> --step <step_size> [--compact]\n";

# Parse command line arguments

my $infile; # bedGraph input file name
my $outfile; # wig output file name
my $step; # wig step size
my $compact; # if selected, steps with value equal to 0 will not be printed

GetOptions (
  "bedgraph=s" => \$infile,
  "wig=s" => \$outfile,
  "step=i" => \$step,
  "compact" => \$compact,
) or die ("Error in command line arguments!\n$usage\n");

# Open input file. If it is compressed with gunzip, uncompress it.

if ($infile =~ /\.gz$/){
  open(IN,'-|',"gunzip -c $infile") || die "Could not open $infile: $!\n";
} else {
  open(IN,'<',$infile) || die "Could not open $infile: $!\n";
}

# Open output file.

open(OUT,'>',$outfile) || die "Could not open $outfile: $!\n";


# bedGraph to wig conversion starts here.

# Print main header.

print OUT "track type=wiggle_0 name=\"$infile\" description=\"$infile\" visibility=full\n";

# Initialize variables.

my $cur_chr = 0; # chromosome being processed
my $cur_pos = 0; # position of current step
my $next_pos = $cur_pos + $step; # position of next step
my $exp_pos = 0; # expected position if no step was skipped; used with --compact option
my $cur_val = 0; # value for current step

while (<IN>) {

  chomp;

  # Skip comment lines
  next if (/^track/);
  next if (/^#/);

  # Parse relevant information in current line 
  # e.g: chr1 3000400 3000500 2
  my ($chr, $start, $end, $val) = split(/\t/);
  
  # Print header for new chromosome and initialize variables.
  if ($chr ne $cur_chr) {
    $cur_chr = $chr;
    $cur_pos = 0;
    $next_pos = $cur_pos + $step;
    $cur_val = 0;
    if (!$compact) { # If --compact option selected, header will be printed immediately before non-null value.
      print OUT "fixedStep chrom=$chr start=", $cur_pos + 1, " step=$step span=$step\n";
      # +1 was added to convert from 0-based to 1-based coordinates.
    }
  }
  
  # Print values when gap in bedGraph file is greater than step.
  while ($start >= $next_pos) {
    print_wig_line($cur_chr, \$cur_pos, \$next_pos, \$exp_pos, \$cur_val, $chr, $start, $end, $val, $step);
  }
  
  # Print values when step overlaps with bedGraph interval and bedGraph interval is longer than step.
  while ($end >= $next_pos) {
    $cur_val += $val * ($next_pos - max($cur_pos, $start));
    print_wig_line($cur_chr, \$cur_pos, \$next_pos, \$exp_pos, \$cur_val, $chr, $start, $end, $val, $step);
  }

  # Update value when end of bedGraph interval is contained within step.
  if ($end < $next_pos) {
    $cur_val += $val * ($end - max($cur_pos, $start));
  }

}

close(IN);
close(OUT);

exit(0);


# Print or skip line in wig file depending on --compact option and value; update variables for next step.
sub print_wig_line {
  my ($cur_chr, $cur_pos, $next_pos, $exp_pos, $cur_val, $chr, $start, $end, $val, $step) = @_;
  if (!$compact) { # Always print if --compact option was not selected.
    my $cur_ave_val = $$cur_val / $step;
    print OUT "$cur_ave_val\n";
  } elsif ($$cur_val != 0) { # Skips printing if --compact option selected and value is null.
    if ($$cur_pos == 0 || $$cur_pos != $$exp_pos) {
      # Adds header if first step in chromosome, or if previous step had null value and was skipped in print out.
      print OUT "fixedStep chrom=$chr start=", $$cur_pos + 1, " step=$step span=$step\n";
      # +1 was added to convert from 0-based to 1-based coordinates.
    }
    my $cur_ave_val = $$cur_val / $step;
    print OUT "$cur_ave_val\n";
    $$exp_pos = $$next_pos;
  }
  $$cur_pos = $$next_pos;
  $$next_pos = $$cur_pos + $step;
  $$cur_val = 0;
}
