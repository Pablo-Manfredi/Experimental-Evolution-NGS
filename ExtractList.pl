#!/usr/bin/perl -w
use strict; use warnings; use Cwd qw(cwd);
#Declare vriables
my ($Current_StringSCV, $Contig, $SCVpos, %PositionsList);
my ($VcfDIR) = cwd; # get the path of the current directory
my ($fileOUT) = 'SitesOfInterest.tsv';
opendir (DIR1, $VcfDIR)|| die "can't opendir $VcfDIR: $!";my (@Files_in_Dir) = readdir(DIR1);closedir(DIR1);
foreach my $Files(@Files_in_Dir)
  {
  chomp $Files;
  if ($Files =~ /^.*polyploid\.1\.m\.A\.vcf$/)
    {
    open(SCVfile, $Files) || die ("Couldn't open the $Files");
    while (<SCVfile>)
      {
      $Current_StringSCV = $_;
      chomp $Current_StringSCV;
      if ($Current_StringSCV =~ /^([^\#][^\t]+)\t([0-9]+)\t[^\t]+\t[^\t]+\t[^\t]+\t([^\t]+)\t.+$/)
         {if ($3 > 85){$PositionsList{$1}{$2} = 1;}
      }
    }
    close (SCVfile);
  }
}
###########   PRINT THE FILE OUT
open(OUTFILE, ">$fileOUT") || die ("Couldn't open the OUTFILE");
foreach $Contig (sort (keys (%PositionsList)))
    {
    foreach $SCVpos (sort {$a <=> $b}(keys (%{$PositionsList{$Contig}})))
      {
      print OUTFILE $Contig."\t".$SCVpos."\n";
      }
    }
close (OUTFILE);
exit;