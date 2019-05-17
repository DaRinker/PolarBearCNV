#!/usr/bin/perl -w
use strict;

##### USAGE
##### ./overlap.pl gene_coordinate_file cnv_file > output_file
#
# 	"gene_coordinate_file" is tab delimited with each gene/feature on a separate line
#	"gene_coordinate_file" is 5 columns consisting of: LOCUS	START	STOP	STRAND	CHROMOSOME
#	
#	Example "gene_coordinate_file":
#	LOCUS	START	STOP	STRAND	CHROMOSOME
#	ID=Uma_R000001Source=ENSP00000155926Function="TRIB2"	384179	406202	-	scaffold79
#	ID=Uma_R000002Source=ENSP00000379406Function="LPIN1"	1140529	1192082	-	scaffold79
#
#	"cnv_file" is tab delimited output file of control-freec ending in "CNVs"	
#	"cnv_file" consists of 5 columns: scaffold/chromosome	start	stop	copy_number	gain or loss
#
#	Example "cnv_file" (no header line necessary):
#	scaffold79      60000   77999   4       gain
#	scaffold79      492000  509999  5       gain
#


#initial variants
my @gene_raw;
my @cnv_raw;
my @final_1;
my @final_2;


#open gene file                                      
@gene_raw = openfile_array($ARGV[0]); 

#open cnv_table file
@cnv_raw = openfile_array($ARGV[1]);
 
my $size_cnv = scalar(@cnv_raw);
my $size_gene = scalar(@gene_raw);

#search through cnv file for each gene

for ( my $i=0; $i<$size_gene; $i++){   
  my @gene_raw_line = split(/\s/,$gene_raw[$i]);  
  my $count_hit=0;
  my $count_mishit=0;
  my $line;
   foreach my $cnv_raw_line(@cnv_raw){
   my @cnv_raw_line=split(/\s/,$cnv_raw_line);
   if ($gene_raw_line[4] eq $cnv_raw_line[0] && $gene_raw_line[2]>= $cnv_raw_line[1] && $gene_raw_line[1]<= $cnv_raw_line[2]){
     $count_hit= $count_hit + 1;
     if ($count_hit==1){
      $line = $gene_raw_line[0]." ".$gene_raw_line[4]." ".$gene_raw_line[1]." ".$gene_raw_line[2]." ".$cnv_raw_line[1]." ".$cnv_raw_line[2]." ".$cnv_raw_line[3]."\n";
      $final_1[$i]=$line;
      }elsif($count_hit == 2){
      $line=~ s/\n/ $cnv_raw_line[3] flag\n/g; 
      $final_1[$i]=$line;}elsif($count_hit>2){
        $line=~ s/flag\n/$cnv_raw_line[3] flag\n/g;
        $final_1[$i]=$line;}
        } else {
     $count_mishit = $count_mishit + 1;
     if ($count_mishit == $size_cnv){
     my $line =  $gene_raw_line[0]." ".$gene_raw_line[4]." ".$gene_raw_line[1]." ".$gene_raw_line[2]." * * 2"."\n";
     push (@final_1, $line);
     }
    }
  }
}

foreach my $final_1_line(@final_1){
  my @final_1_line=split(/\s/,$final_1_line);
  my $size_line = scalar (@final_1_line);
  if ($final_1_line[$size_line-1] ne 'flag') {
    push (@final_2,$final_1_line);
  }else{
     my $sum=0;
     for (my $i=6 ; $i< $size_line-1; $i++){
        $sum = $sum + $final_1_line[$i];
     }
     my $average = $sum/($size_line - 7);
     my $line= $final_1_line[0]." ".$final_1_line[1]." ".$final_1_line[2]." ".$final_1_line[3]." ".$final_1_line[4]." ".$final_1_line[5]." "."$average flag\n";
    push (@final_2,$line);
  }
}

shift @final_2;
unshift (@final_2,'GENE_ID SCAFFOLD GENE_START GENE_STOP CNV_START CNV_STOP CN NOTE'."\n");
 
print @final_2;



########################openfile_array#################################
sub openfile_array {
  my $filename= $_[0];
  unless (open(FILE, $filename)){
  print "cannot find the file \"$filename\"!!!";
  }
  open (FILE, $filename);
  my @data=<FILE>;
  close FILE;
  return @data;
}

