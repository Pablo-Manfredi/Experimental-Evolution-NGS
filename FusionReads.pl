#!/usr/bin/perl -w
#Nom du sub.       : FusionReads.pl
#Author            : Dr. Pablo Manfredi
#Goal              : identify fusion reads in genome assemblies
#Inpout            : an assembly .sam file
#Outpout           : a .tsv file containing the cordinates of patentially recombined sites
#Outpout           : a .tsv file containing the cordinates of the reads detected as fusion reads
#Outpout           : a .fastq file containing all the reads identified as fusion reads for de novo assembly purpose
use strict; use warnings;use Cwd qw(cwd);

##########      START OF ADAPT TO YOUR NEEDS
my ($PartialMatchCutoff)= 20;
##########      ENDO OF ADAPT TO YOUR NEEDS

#Declaring variables needed
my ($Sample,$out_file,$phrase,@word,%Readsword,$Match,$SoftClip,$readRef,$Pair,$replicon,$PositionToReccord,%FullMatch,%ReadsMemory,%CordinatesFromSecondary,$possibleReplicon, $readLength);

# The script will look for all files with a .sam extention in its directory
my ($VcfDIR) = cwd; # get the path of the current directory
opendir (DIR1, $VcfDIR)|| die "can't opendir $VcfDIR: $!";my (@Files_in_Dir) = readdir(DIR1);closedir(DIR1);
foreach my $File_in(@Files_in_Dir)
  {
  chomp $File_in;
  if ($File_in =~ /^(.*)\.sam$/)
    {
    $Sample = $1;
    undef (%Readsword);undef (%FullMatch);undef (%ReadsMemory);undef (%CordinatesFromSecondary);
##########      START PARSING OF THE SAM FILE
=example_data
BS-DSFCONTROL06:65:C6E8JANXX:8:1101:8066:2427   99      gi|110227054|gb|AE004091.2|     408824  255     55M71S  =       408824  -197    CGAACAGATCGGGTACAGGCATGAATGCACTCTCAGGCTAGGCGGCTTCGGGCTCAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGTTCTAATCTCGTATGCCGTCTTCTGCTTGAAAAA  CCCCBEGGGGGGGGGGGGGGGGGG>FGGGGGGGGGGGGGGGGGGEG<EFEBGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFG=0FGGGGGGGGGGGGEG/DGEGGGGGGG  AS:i:110        XN:i:0
XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:55 YS:i:110        YT:Z:CP
BS-DSFCONTROL06:65:C6E8JANXX:8:1101:8066:2427   147     gi|110227054|gb|AE004091.2|     408824  255     71S55M  =       408824  197     TTCTTTTTTTTTTAATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTC
GAACAGATCGGGTACAGGCATGAATGCACTCTCAGGCTAGGCGGCTTCGGGCTC  ...GGGGGGGGGGGGGGGG>GD>CEFFGFGGGGGCGGGGGFGDGGEGGGGFFFEGBEE?GGGFGGFGFGGGGGFDGGGGGGGGGGFGGGFFFB<1EGGGFGGFBGFGEGGBGGGGGGGGGGCCCBB  AS:i:110        XN:i:0
NB501515:29:H7TGCAFXY:1:11101:18364:1057        83      AE004091.2_PAO1 5133034 255     150M    =       5132759 -425    GCGCTCGCCGGTCGCTTCAAGGCGATCACCGACAAGCACCCGAGGCTGGAGAAGACCAACAACCGCCAGCACCTGGGAACCCTCGGCACCGGCAACCACTTCATCGAGGTCTGCCTGGACGAGGCCGACCGCGTCTGGTTCATGCTGCAC        /<A<AAA<EAAEEEEEEAAE/AAEEEEEEEAEEAEA<EEEE<AEEEEEEEEEE/EEEEEEEEEEAEA/AAEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEAEAEEEE/EEEEEEEEEEEAA6AA  AS:i:300        XN:i:0  XM:i:0
        XO:i:0  XG:i:0  NM:i:0  MD:Z:150        YS:i:298        YT:Z:CP
NB501515:29:H7TGCAFXY:2:11204:5601:20105        99      AE004091.2_PAO1 3991681 11      57S93M  =       3991681 208     ACCAGGGCCCGGACCTCGGCGGCGCTCGCCAGGCCCAGGGCTTTGCGGGCCAGGCCGCGCGCAGCTCGTCCACTCCCAGGCCGACCAGCAGCGGCAACGCCAGCGGGTCGGCCGCCAGTTCGCCGCAGACGCCGCCCCACTTGCCCTCGG    AAAAAEEEEAEEAAEEEEEEEEEE//E/EEAAEEEE//EEE/EEE</AEAEAAEEE/EAE/E/////A//</A/<EE//EE<//AE//E//</AE/<A<AE//A/<E//AE<E</E///A/<<E/E///6/<6///6</<///AAA////   AS:i:176        XS:i:130        XN:i:0  XM:i:2  XO:i:0  XG:i:0  NM:i:2  MD:Z:3T73A15    YS:i:213 YT:Z:CP
=cut
    open(DATA_FILE,$File_in) || die ("Couldn't open the DATA_FILE");
    while (<DATA_FILE>)
                {
                $phrase = $_;
                chomp $phrase;
                unless ($phrase=~/^\@.+/)
                  {
                  @word = split(/\t/, $phrase);
                  $replicon = $word[2];
                  $readRef =$word[0];
                  # Interpreting the sam format tags in terms of read pairs for a pair end sequencing
                  # This is required as the pair number modifies the rest of the tag
                  if ($word[1] < 128)
                     {
                     $Pair = 1;
                     }
                  elsif (($word[1] > 127)&&($word[1] < 256))
                     {
                     $Pair = 2;
                     }
                  elsif (($word[1] > 255)&&($word[1] <384))
                     {
                     $Pair = 1;
                     }
                  elsif ($word[1] > 383)
                     {
                     $Pair = 2;
                     }
                  else {print 'Big Problem'."\n";}
                  $readLength = scalar($word[1]);
                  $readRef = $readRef.'.'.$Pair;
                  # Interpreting the CIGAR coded alignement
                  if ($word[5] eq $readLength.'M')
                     {
                     $Match = $1;
                     $FullMatch{$readRef}=1;
                     }
                  elsif ($word[5] =~ /^([0-9]+)M.*$/)
                     {
                     $Match = $1;
                     if ($word[5] =~ /.*[A-Z]([0-9]+)S$/)
                        {
                        $SoftClip = $1;
                        if ($Match > $PartialMatchCutoff && $SoftClip > $PartialMatchCutoff)
                            {
                            $PositionToReccord = $Match+$word[3];#1-based offset into the forward reference strand where leftmost character of the alignment occurs
                            $Readsword{$readRef}{'MS'}{$replicon}{$PositionToReccord}++;
                            if ($word[1] > 255){$ReadsMemory{$word[0]}{$Pair}{$replicon}{'Secondary'}='MS';}
                            elsif ($word[1] < 256){$ReadsMemory{$word[0]}{$Pair}{$replicon}{'Primary'}='MS';$ReadsMemory{$word[0]}{$Pair}{$replicon}{'PrimaryPos'}=$PositionToReccord;}
                            $ReadsMemory{$word[0]}{$Pair}{$replicon}{'Seq'}=$word[9];
                            $ReadsMemory{$word[0]}{$Pair}{$replicon}{'Qual'}=$word[10];
                            foreach $possibleReplicon (keys %{$ReadsMemory{$word[0]}{$Pair}})
                                    {
                                    if (
                                       ($ReadsMemory{$word[0]}{$Pair}{$replicon}{'Secondary'})
                                       &&
                                       ($ReadsMemory{$word[0]}{$Pair}{$possibleReplicon}{'Primary'})
                                       &&
                                       (($ReadsMemory{$word[0]}{$Pair}{$possibleReplicon}{'Primary'} ne $ReadsMemory{$word[0]}{$Pair}{$replicon}{'Secondary'})||($possibleReplicon ne $replicon))
                                       )
                                          {
                                          $ReadsMemory{$word[0]}{$Pair}{$replicon}{'Print'} =1;
                                          $CordinatesFromSecondary{$possibleReplicon}{$ReadsMemory{$word[0]}{$Pair}{$possibleReplicon}{'PrimaryPos'}}{$replicon}{$PositionToReccord}++;
                                          }
                                    }
                            }
                        }
                     }
                  elsif ($word[5] =~ /^([0-9]+)S.*/)
                     {
                     $SoftClip = $1;
                     if ($word[5] =~ /.*[A-Z]([0-9]+)M$/)
                        {
                        $Match = $1;
                        if ($Match > $PartialMatchCutoff && $SoftClip > $PartialMatchCutoff)
                            {
                            $PositionToReccord = $word[3];
                            $Readsword{$readRef}{'SM'}{$replicon}{$PositionToReccord}++;
                            if ($word[1] > 255){$ReadsMemory{$word[0]}{$Pair}{$replicon}{'Secondary'}='SM';}
                            elsif ($word[1] < 256){$ReadsMemory{$word[0]}{$Pair}{$replicon}{'Primary'}='SM';$ReadsMemory{$word[0]}{$Pair}{$replicon}{'PrimaryPos'}=$PositionToReccord;}
                            $ReadsMemory{$word[0]}{$Pair}{$replicon}{'Seq'}=$word[9];
                            $ReadsMemory{$word[0]}{$Pair}{$replicon}{'Qual'}=$word[10];
                            
                            foreach $possibleReplicon (keys %{$ReadsMemory{$word[0]}{$Pair}}) 
                                    {
                                    if (
                                    ($ReadsMemory{$word[0]}{$Pair}{$replicon}{'Secondary'})
                                     &&
                                     ($ReadsMemory{$word[0]}{$Pair}{$possibleReplicon}{'Primary'})
                                     &&
                                     (($ReadsMemory{$word[0]}{$Pair}{$possibleReplicon}{'Primary'} ne $ReadsMemory{$word[0]}{$Pair}{$replicon}{'Secondary'})||($possibleReplicon ne $replicon))
                                        )
                                          {
                                          $ReadsMemory{$word[0]}{$Pair}{$replicon}{'Print'} =1;
                                          $CordinatesFromSecondary{$possibleReplicon}{$ReadsMemory{$word[0]}{$Pair}{$possibleReplicon}{'PrimaryPos'}}{$replicon}{$PositionToReccord}++;
                                          }
                                    }
                            }
                        }
                     }
                  }
                }
    close (DATA_FILE);
##########      STOP PARSING OF THE SAM FILE


##########      START PRINTING FUSION.POINTS.fastq
    my ($out_file_Seq) = $File_in.'_FUSION.POINTS.fastq';
    open(FILEOUT,">$out_file_Seq") || die ("Couldn't open the FILEOUT");
    my ($CurrentReadName);
    foreach $CurrentReadName (keys %ReadsMemory)
        {
        foreach $replicon (keys %{$ReadsMemory{$CurrentReadName}{1}})
          {
          if($ReadsMemory{$CurrentReadName}{1}{$replicon}{'Print'})
            {
            if ($ReadsMemory{$CurrentReadName}{1}{$replicon}{'Print'} ==1)
              {
              print FILEOUT '@'.$CurrentReadName."\n".$ReadsMemory{$CurrentReadName}{1}{$replicon}{'Seq'}."\n".'+'."\n".$ReadsMemory{$CurrentReadName}{1}{$replicon}{'Qual'}."\n";
              }
            }
          }
        foreach $replicon (keys %{$ReadsMemory{$CurrentReadName}{2}})
          {
          if($ReadsMemory{$CurrentReadName}{2}{$replicon}{'Print'})
            {
            if ($ReadsMemory{$CurrentReadName}{2}{$replicon}{'Print'} ==1)
              {
              print FILEOUT '@'.$CurrentReadName."\n".$ReadsMemory{$CurrentReadName}{2}{$replicon}{'Seq'}."\n".'+'."\n".$ReadsMemory{$CurrentReadName}{2}{$replicon}{'Qual'}."\n";
              }
            }
          }
        }
    close (FILEOUT);
##########      STOP PRINTING FUSION.POINTS.fastq

##########      START PRINTING FFUSION.LINKING.tsv
    my ($out_file_Seq_Cordinates) = $File_in.'_FUSION.LINKING.tsv';
    open(FILEOUT2,">$out_file_Seq_Cordinates") || die ("Couldn't open the FILEOUT2");
    my ($CurrentSeconPOs);
# Here $CurrentReadName is the first position
    foreach $possibleReplicon (keys %CordinatesFromSecondary)
        {
        #print $replicon."\n";
        foreach $CurrentReadName (keys %{$CordinatesFromSecondary{$possibleReplicon}})
          {
          #print $replicon."\t".$CurrentReadName."\n";
          foreach $replicon (keys %{$CordinatesFromSecondary{$possibleReplicon}{$CurrentReadName}})
            {
            foreach $CurrentSeconPOs (keys %{$CordinatesFromSecondary{$possibleReplicon}{$CurrentReadName}{$replicon}})
                    {
                    if ($CordinatesFromSecondary{$possibleReplicon}{$CurrentReadName}{$replicon}{$CurrentSeconPOs})
                        {
                        print FILEOUT2 $possibleReplicon."\t".$CurrentReadName."\t".$replicon."\t".$CurrentSeconPOs."\t".$CordinatesFromSecondary{$possibleReplicon}{$CurrentReadName}{$replicon}{$CurrentSeconPOs}."\n";
                        }
                    }
            }
          }
        }
    close (FILEOUT2);
##########      STOP PRINTING FUSION.LINKING.tsv

##########      START Counting fusion reads per position
    my (%Counter);
    my ($CurrentRead);
    my ($Supatut);
    foreach $CurrentRead (keys %Readsword)
        {
        foreach $replicon (keys %{$Readsword{$CurrentRead}{'SM'}})
          {
          if ((keys %{$Readsword{$CurrentRead}{'SM'}{$replicon}})&&(keys %{$Readsword{$CurrentRead}{'MS'}{$replicon}})&&(! $FullMatch{$readRef}))
             {
             foreach $Supatut (keys %{$Readsword{$CurrentRead}{'SM'}{$replicon}})
                     {
                     $Counter{$replicon}{$Supatut}++;
                     }
             foreach $Supatut (keys %{$Readsword{$CurrentRead}{'MS'}{$replicon}})
                     {
                     $Counter{$replicon}{$Supatut}++;
                     }
             }
          }
        }
##########      STOP Counting fusion reads per position

##########      START PRINTING FUSION.POINTS.tsv
    $out_file = $File_in.'_FUSION.POINTS.tsv';
    open(FILEOUT3,">$out_file") || die ("Couldn't open the FILEOUT3");
    foreach $replicon (keys %Counter)
        {
        foreach $Supatut (keys %{$Counter{$replicon}})
          {
          if ($Counter{$replicon}{$Supatut} > 1)
            {
            print FILEOUT3 $replicon."\t".$Supatut."\t".$Counter{$replicon}{$Supatut}."\n";
            }
          }
        }
    close (FILEOUT3);
##########      STOP PRINTING FUSION.POINTS.tsv

    }
}
exit;