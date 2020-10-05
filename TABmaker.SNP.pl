#!/usr/bin/perl -w
use warnings; use Cwd qw(cwd); use strict;
my ($Sample, $Current, $Contig, $SNPpos, %SNPdata, %SampleList, %PositionsList, $compareVar, $AltBase, $readDepth, $IsINDEL, %Altcov, $ref, $AltBaseQ, $HighestPercent);
my ($VcfDIR) = cwd; # get the path of the current directory
opendir (DIR1, $VcfDIR)|| die "can't opendir $VcfDIR: $!";my (@Files_in_Dir) = readdir(DIR1);closedir(DIR1);
foreach my $Files(@Files_in_Dir)
  {
  chomp $Files;
  if ($Files =~ /^(.*)\.AllBase\.vcf$/)
    {
    $Sample = $1; $SampleList{$Sample}=1;
    open(FILEIN, $Files) || die ("Couldn't open the $Files");
    Loop1:while (<FILEIN>)
      {
      $Current = $_;
      chomp $Current;
      undef(%Altcov);
      $IsINDEL = 'N';
      if ($Current =~ /^\#.*$/){next Loop1;}
      if ($Current =~ /^.+INDEL\;IDV\=([^\;]+)\;IMF\=[^\;]+\;DP.*$/)
        {if ($1 > 4){$IsINDEL = 'Y';}else{next Loop1;}}
      if ($Current =~ /^([^\#][^\t]+)\t([0-9]+)\t[^\t]+\t([^\t]+)\t([^\t]+)\t([^\t]+)\t[^\t]+\t.*DP\=([0-9]+)\;.*DP4\=([0-9]+)\,([0-9]+)\,([0-9]+)\,([0-9]+)\;.*$/){
        $Contig = $1; $SNPpos = $2;$ref = $3;$AltBase = $4;$AltBaseQ = $5;$readDepth = $6;
        $Altcov{'AltFract'} =($9+$10)/$readDepth;
        $Altcov{'AltFract'} = sprintf("%.1f", $Altcov{'AltFract'});
        $PositionsList{$Contig}{$SNPpos}{$IsINDEL}{'PrintDif'} = 0;
        $PositionsList{$Contig}{$SNPpos}{$IsINDEL}{'Print10Percent'} = 0;
        if ($AltBase eq '.'){$AltBase = $ref;}
        if ($SNPdata{$Contig}{$SNPpos}{$IsINDEL}{$Sample}{'Quant'})
           {if ($SNPdata{$Contig}{$SNPpos}{$IsINDEL}{$Sample}{'Quant'} > $Altcov{'AltFract'}){next Loop1;}}
        $SNPdata{$Contig}{$SNPpos}{$IsINDEL}{$Sample}{'Quant'}=$Altcov{'AltFract'};
        $SNPdata{$Contig}{$SNPpos}{$IsINDEL}{$Sample}{'Qual'}=$AltBaseQ;
        $SNPdata{$Contig}{$SNPpos}{$IsINDEL}{$Sample}{'Event'} = $ref.'/'.$AltBase;
        }
       else {print $Current."\n";}
       }
      close (FILEIN);
      } 
    }
#################     Exclude Positions with constant variants
foreach $Contig (sort (keys (%PositionsList)))
    {
    foreach $SNPpos (sort {$a <=> $b}(keys (%{$PositionsList{$Contig}})))
      {
      foreach $IsINDEL (sort (keys (%{$PositionsList{$Contig}{$SNPpos}})))
        {
        $compareVar = 'ND';$HighestPercent = 0;
        $PositionsList{$Contig}{$SNPpos}{$IsINDEL}{'PrintDif'} = 0;
        $PositionsList{$Contig}{$SNPpos}{$IsINDEL}{'Print10Percent'} = 0;
        foreach $Sample (sort (keys (%SampleList)))
          {
          if ($SNPdata{$Contig}{$SNPpos}{$IsINDEL}{$Sample}){}
          else{$SNPdata{$Contig}{$SNPpos}{$IsINDEL}{$Sample}{'Event'}= 'NO';$SNPdata{$Contig}{$SNPpos}{$Sample}{$IsINDEL}{'Quant'}= 0;}
          if ($compareVar eq 'ND'){$compareVar = $SNPdata{$Contig}{$SNPpos}{$IsINDEL}{$Sample}{'Event'};}
          else{if ((split /\,/, $compareVar)[0] ne (split /\,/, $SNPdata{$Contig}{$SNPpos}{$IsINDEL}{$Sample}{'Event'})[0]){$PositionsList{$Contig}{$SNPpos}{$IsINDEL}{'PrintDif'} = 1;}}
          if ($SNPdata{$Contig}{$SNPpos}{$IsINDEL}{$Sample}{'Quant'}){}else{$SNPdata{$Contig}{$SNPpos}{$IsINDEL}{$Sample}{'Quant'} = 0;}
          if ($HighestPercent < $SNPdata{$Contig}{$SNPpos}{$IsINDEL}{$Sample}{'Quant'}){$HighestPercent = $SNPdata{$Contig}{$SNPpos}{$IsINDEL}{$Sample}{'Quant'};}
          if (!($SNPdata{$Contig}{$SNPpos}{$IsINDEL}{$Sample}{'Qual'})){$SNPdata{$Contig}{$SNPpos}{$IsINDEL}{$Sample}{'Qual'}=0;}
          }
        if ($HighestPercent > 0.1){$PositionsList{$Contig}{$SNPpos}{$IsINDEL}{'Print10Percent'} = 1;}
        } 
      } 
    }
#################     Print the output files
my ($file_name_OUT) = 'Evolution.tsv';
open(OUTFILE, ">$file_name_OUT") || die ("Couldn't open the OUTFILE");
print OUTFILE "Contig\tPosition";
foreach $Sample (sort (keys (%SampleList))){print OUTFILE "\t".'Ref/'.$Sample;}
foreach $Sample (sort (keys (%SampleList))){print OUTFILE "\t".$Sample.'.Ratio';}
foreach $Sample (sort (keys (%SampleList))){print OUTFILE "\t".$Sample.'.Qual';}
print OUTFILE "\n";
foreach $Contig (sort (keys (%PositionsList)))
    {
    foreach $SNPpos (sort {$a <=> $b}(keys (%{$PositionsList{$Contig}})))
      {
      foreach $IsINDEL (sort (keys (%{$PositionsList{$Contig}{$SNPpos}})))
        {
        if ($PositionsList{$Contig}{$SNPpos}{$IsINDEL}{'PrintDif'} == 1 && $PositionsList{$Contig}{$SNPpos}{$IsINDEL}{'Print10Percent'} == 1)
           {
           print OUTFILE $Contig."\t".$SNPpos;
           foreach $Sample (sort (keys (%SampleList)))
             {print OUTFILE "\t".$SNPdata{$Contig}{$SNPpos}{$IsINDEL}{$Sample}{'Event'};}
           foreach $Sample (sort (keys (%SampleList)))
             {print OUTFILE "\t".$SNPdata{$Contig}{$SNPpos}{$IsINDEL}{$Sample}{'Quant'};}
           foreach $Sample (sort (keys (%SampleList)))
             {print OUTFILE "\t".$SNPdata{$Contig}{$SNPpos}{$IsINDEL}{$Sample}{'Qual'};}
           print OUTFILE "\n";
           } 
        }
      }
    }
close(OUTFILE);
exit;