#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Data::Dumper;


# globals

my $EXE = basename($0);
my $VERSION = "0.1.1";
my $AUTHOR = 'Torsten Seemann (@torstenseemann) and Pablo Manuel Gonzalez de la Rosa';
my $HOMEPAGE = "undefined";

# SAM file TSV columns
use constant {
  SAM_RID   => 0,
  SAM_FLAG  => 1,
  SAM_RNAME => 2,
  SAM_POS   => 3,
  SAM_MAPQ  => 4,
  SAM_CIGAR => 5,
  SAM_TLEN  => 8,
  SAM_SEQ   => 9,
};

#----------------------------------------------------------------------
# command line parameters

my $telomere  = "TTAGGC";
my $debug     = 0;

#----------------------------------------------------------------------
sub usage {
  my($exitcode) = @_;
  $exitcode=0 if !defined($exitcode) or $exitcode eq 'help';
  my $fh = $exitcode ? \*STDERR : \*STDOUT;
  print $fh
    "SYNOPSIS\n  Get coordinates of soft clipped telomeric sequences \n",
    "OUTPUT\n  TSV printed to STDOUT. Columns: 1) reference seq ID, 2) position of softclipped telomere,\n",
    "  3) R or L indicating if telomere occurs at right or left of alignment and 4) MAPQ\n",
    "AUTHOR\n  $AUTHOR\n",
    "USAGE\n",
    "  % minimap2 ref.fa telomeric_reads.fa | softclipped_telo_coords.pl > clippedTeloPos.tsv\n",
    "OPTIONS\n",
    "  --help           This help\n",
    "  --version        Print version and exit\n",
    "  --telomere STR   Telomeric repeat for species (default is nematodes' $telomere)\n",
    "  --debug          Print verbose debug info to stderr\n",
    "HOMEPAGE\n  $HOMEPAGE\n",
    "";
  exit($exitcode);
}

#----------------------------------------------------------------------
# getopts

GetOptions(
  "help"       => \&usage,
  "version"    => \&version,
  "telomere=s" => \$telomere,
  "debug"      => \$debug,
) or usage(1);
             
!@ARGV and -t STDIN and err("Please provide or pipe a SAM file to $EXE");

#----------------------------------------------------------------------
# main
msg("$EXE $VERSION by $AUTHOR");


my $total=0;
my $leftend=0;
my $rightend=0;
my $addedlen=0;
# at most this number of bases can lack telomeric repeat
# intended for dealing with read errors
my $telosearchspace =  3 * length($telomere); 
my $revcomptelomere=reverse $telomere;
$revcomptelomere =~ tr/ACGTacgt/TGCAtgca/;

# read SAM one line ar a time
while (my $line = <ARGV>) {
  # skip SAM header
  if ($line =~ m/^@/) {
    next;
  }
  $total++;
  my @sam = split m/\t/, $line;
  # ensure there is softclipped alignment before heavyweight parsing
  my $isSoftClipped = ($sam[SAM_CIGAR] =~ /\dS/);
  my $isHardClipped = ($sam[SAM_CIGAR] =~ /\dH/);
  my $isPrimaryAlignment = !($sam[SAM_FLAG] & 2048 or $sam[SAM_FLAG] & 256);
  if ($isSoftClipped and not $isHardClipped and $isPrimaryAlignment ) {
    my $forwardStrand = !($sam[SAM_FLAG] & 16);
    my $start = $sam[SAM_POS];
    my $readlen = length($sam[SAM_SEQ]);
    my $end = $start + $readlen  - 1;
    my $contigname = $sam[SAM_RNAME];
    my ($SL, undef, $SR) 
      = ($sam[SAM_CIGAR] =~ m/ ^ (?:(\d+)S)? (.*?) (?:(\d+)S)? $/x);
    $SL ||= 0; $SR ||= 0;

    ## Check telomere on the right side of read
    if ( $SR > $telosearchspace ) {
      my $potentialSequence = substr($sam[SAM_SEQ], $readlen - $SR, $readlen + 1);
      my $endswithtelomere = ($potentialSequence =~ /\Q$telomere\E\w{0,\Q$telosearchspace\E}$/);
      if ($endswithtelomere){
        # adjust end of alignment by indels and deletion length
        # this is obtained by parsing the CIGAR string
        # code adapted from https://github.com/holmeso/adamaperl/blob/56fee71f0c431e3b98ef5f3fc1a2a7213e755e14/lib/QCMG/SamTools/Bam/Alignment.pm#L569
        my $adjust = 0;
        my @cigar  = $sam[SAM_CIGAR] =~ /(\d+)(\w)/g;
        while (@cigar) {
          my ($len,$op) = splice(@cigar,0,2);
          $adjust += $len if $op eq 'I';
          $adjust -= $len if $op eq 'D';
        }
        $adjust += $SL + $SR;
        my $adjustedend = $start + $readlen - $adjust - 1; 
        print join("\t", $sam[SAM_RNAME], $adjustedend, "R", $sam[SAM_MAPQ]), "\n";
        $rightend++;
      }

    } 
    if ( $SL > $telosearchspace ){
      ## Check telomere on the left side of read
      my $potentialSequence = substr($sam[SAM_SEQ], 0, $SL + 1);
      my $startswithtelomere = ($potentialSequence =~ /^\w{0,\Q$telosearchspace\E}\Q$revcomptelomere/);
      if($startswithtelomere){
        print join("\t", $sam[SAM_RNAME], $start, "L", $sam[SAM_MAPQ]), "\n";
        $leftend++;
      }
    }
  }
}



# stats
msg("Total SAM records $total");
msg("Left end softclipped telomeres $leftend");
msg("Right end softclipped telomeres $rightend");


msg("Done.");


#----------------------------------------------------------------------
sub version {
  print "$EXE $VERSION\n";
  exit(0);
}

#----------------------------------------------------------------------
sub msg {
  print STDERR "[$EXE] @_\n";
}

#----------------------------------------------------------------------
sub err {
  msg("ERROR:", @_);
  exit(1);
}