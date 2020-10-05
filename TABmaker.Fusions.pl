#!/usr/bin/perl -w
#Nom du sub.       : FusionReads.pl
#Author            : Dr. Pablo Manfredi
#Goal              : identify fusion reads in genome assemblies
#Inpout            : an assembly .sam file
#Outpout           : a .tsv file containing the cordinates of patentially recombined sites
#Outpout           : a .tsv file containing the cordinates of the reads detected as fusion reads
#Outpout           : a .fastq file containing all the reads identified as fusion reads for de novo assembly purpose
use warnings;use strict;use Cwd qw(cwd);

##########      START OF ADAPT TO YOUR NEEDS
my ($FusionsPerPositionCutoff)= 4;
##########      ENDO OF ADAPT TO YOUR NEEDS

#Declaring variables needed
my ($Sample,%SampleList,$file_name_OUT, $VCFfile,$Current,$phrase,@word,$countStartReads,$countStopReads,$countBULLshit,%Reference,$InsertLength,$pattern,%Comparison,$NaCounter,$CoverageConfirmation);
my ($currentPosit,%FusPoints, %SCVdata,$SCVpos,$replicon,$LinkTo,$LinkToreplicon,$QuantifLink);
my ($lastPA,$cov,$posChromosome,$Contig,$PaRefs);
my ($VcfDIR) = cwd; # get the path of the current directory
opendir (DIR1, $VcfDIR)|| die "can't opendir $VcfDIR: $!";my (@Files_in_Dir) = readdir(DIR1);closedir(DIR1);

###    1 
###########################     TABLE containing the cordinates of patentially recombined sites
foreach my $Files(@Files_in_Dir)
  {
  chomp $Files;
  if ($Files =~ /^(.*)\.sam\_FUSION\.POINTS\.tsv$/)
    {
    $Sample = $1; $SampleList{$Sample}=1;
    open(FILEIN1, $Files) || die ("Couldn't open the $Files");
    Loop1:while (<FILEIN1>)
      {
      $Current = $_;
      chomp $Current;
      if ($Current =~ /^([^\t]+)\t([0-9]+)\t([0-9]+)$/)
        {
        $replicon = $1;
        $SCVpos = $2;
        if ($3 > $FusionsPerPositionCutoff)
          {
          $FusPoints{$replicon}{$SCVpos}{$Sample}=$3;
          }
        }
      }
    close (FILEIN1);
    }
  }

###########################     Print OUT the file
$file_name_OUT = 'PositionFusions.tsv';
open(OUTFILE, ">$file_name_OUT") || die ("Couldn't open the OUTFILE");
print OUTFILE "Contig\tPosition";
foreach $Sample(sort (keys (%SampleList)))
   {
   print OUTFILE "\t".$Sample;
   }
print OUTFILE "\n";
foreach $replicon (sort {$a cmp $b} (keys (%FusPoints)))
  {
  foreach my $currentPosit (sort {$a <=> $b} (keys (%{$FusPoints{$replicon}})))
      {
      print OUTFILE $replicon."\t".$currentPosit;
      foreach $Sample(sort (keys (%SampleList)))
        {
        if (!($FusPoints{$replicon}{$currentPosit}{$Sample})){$FusPoints{$replicon}{$currentPosit}{$Sample}=0;}
        print OUTFILE "\t".$FusPoints{$replicon}{$currentPosit}{$Sample};
        }
      print OUTFILE "\n";
      }
  }
close (OUTFILE);


###    2
###########################     TABLE containing the cordinates of the reads detected as fusion reads
foreach my $Files(@Files_in_Dir)
  {
  chomp $Files;
  if ($Files =~ /^(.*)\.sam\_FUSION\.LINKING\.tsv$/)
    {
    $Sample = $1; $SampleList{$Sample}=1;
    open(FILEIN2, $Files) || die ("Couldn't open the $Files");
    Loop2:while (<FILEIN2>)
      {
      $Current = $_;
      chomp $Current;
      if ($Current =~ /^([^\t]+)\t([0-9]+)\t([^\t]+)\t([0-9]+)\t([0-9]+)$/)
        {
        $replicon = $1;
        $SCVpos = $2;
        $LinkToreplicon = $3;
        $LinkTo = $4;
        $QuantifLink=$5;
        if ($QuantifLink > $FusionsPerPositionCutoff)
          {
          $SCVdata{$replicon}{$SCVpos}{$LinkToreplicon}{$LinkTo}{$Sample} = $QuantifLink;
          }
        }
      }
    close (FILEIN2);
    }
  }
###########################     ASSOCIATING COORDINATES linked by a fusion read
foreach $replicon(keys %SCVdata)
                   {
                   foreach $SCVpos(sort (keys %{$SCVdata{$replicon}}))
                           {
                           foreach $LinkToreplicon(sort (keys %{$SCVdata{$replicon}{$SCVpos}}))
                               {
                               foreach my $LinkTo(sort (keys %{$SCVdata{$replicon}{$SCVpos}{$LinkToreplicon}}))
                                       {
                                       foreach $Sample(sort (keys (%SampleList)))
                                           {
                                           if (($SCVdata{$replicon}{$SCVpos}{$LinkToreplicon}{$LinkTo}{$Sample})&&($SCVdata{$replicon}{$SCVpos}{$LinkToreplicon}{$LinkTo}{$Sample} > 0)){}
                                           else{$SCVdata{$replicon}{$SCVpos}{$LinkToreplicon}{$LinkTo}{$Sample} = '0';}
                                           }
                                       }
                               }
                           }
                   }
###########################     Print OUT the file
$file_name_OUT= 'PositionLinks.tsv';
open(OUTFILE, ">$file_name_OUT") || die ("Couldn't open the OUTFILE");
print OUTFILE "LinkFROM\tLinkTO";
foreach $Sample(sort (keys (%SampleList)))
   {
   print OUTFILE "\t".$Sample;
   }
print OUTFILE "\n";
foreach my $replicon (sort {$a cmp $b} (keys (%SCVdata)))
  {
  foreach $currentPosit (sort {$a <=> $b} (keys (%{$SCVdata{$replicon}})))
      {
      foreach my $currentLink (sort {$a cmp $b} (keys (%{$SCVdata{$replicon}{$currentPosit}})))
              {
              foreach my $currentLinkPosi (sort {$a <=> $b} (keys (%{$SCVdata{$replicon}{$currentPosit}{$currentLink}})))
                      {
                      print OUTFILE $replicon."\t".$currentPosit."\t".$currentLink."\t".$currentLinkPosi;
                      foreach $Sample(sort (keys (%SampleList)))
                             {
                             print OUTFILE "\t".$SCVdata{$replicon}{$currentPosit}{$currentLink}{$currentLinkPosi}{$Sample};
                             }
                      $Contig = $replicon;
                      print OUTFILE "\n";
                      }
              }
      }
  }
close (OUTFILE);
exit;