=cut
	function:
	1. Perform statistic analysis on the discordant-aligned fragments from each sample:
    	(1) Proportions of fragments derived from one chromosome, or >1 chromosomes; 
    	(2) Proportion of fragments with identical 6-mer end motifs in each type
    2. Perform  identical 6-mer end motifs analysis on the discordant-aligned fragments from each group:
		(1) Identical 6-mer end motifs profile of discordant-aligned fragments



=cut


use strict;
my $list=shift;
my %K6_motif=();
my $N=0;  
my $Ends_equal_k6=0;
#Fragments from one chromosome 
my $N_S_T=0; #Count of Total
#Fragments from two  chromosomes with one junction
my $N_T_T=0; #Count of Total
my $out="statics_discordant_fragments_".$list;
open(LI,$list)or die $!;
open(OUT,">>$out")or die $!;
	print OUT "file\tFragment Count\tFragments from one chromosome\tFragments from more than one chromosomes\tIdentical 6-mer in discordantly-aligned fragments\n";
while(my $file=<LI>){
	chomp $file;
    $N=0;
    $Ends_equal_k6=0;
    $N_S_T=0; 
    $N_T_T=0;
    my %read=();

    open(IN,$file)or die $!;
	while(my $line=<IN>){
		if($line=~/^\@/){
		# print $line;		
		}else{
			my @temp=split /\t/,$line;
			if(!exists($read{$temp[0]})){
				$read{$temp[0]}=$line;
			}else{
				$read{$temp[0]}.=$line;
			}
		}
	}
	close(IN);

	my @read_id=keys(%read);
	
	for(my $i=0;$i<@read_id;$i++){
   	 	my @temp_k6=();
   	 	my $k6="";
   	 	my %chr=();
   	 	
		my @pair=split /\n/, $read{$read_id[$i]};
		for my $each(@pair){
			my @temp2=split /\t/,$each;
			my $test=$temp2[1] & 16;
			if($test == 0){
				$k6=substr($temp2[9],0,6);
				$chr{$temp2[2]}=1;
			}elsif($test == 16){
				$k6=reverse(substr($temp2[9],-6));
        		my @tempkmer=split //,$k6;
       	 		my $tempk6="";
        		for(@tempkmer){
        	    	if($_ eq "A"){
        	    		$tempk6.="T"
         	  		 }elsif($_ eq "T"){
          		  		$tempk6.="A"
            		}elsif($_ eq "G"){
            			$tempk6.="C"
            		}elsif($_ eq "C"){
            			$tempk6.="G"
            		}else{  # for N
						$tempk6.=$_;
            			#print "error...$k6....$tempk6....$temp2[9]...\n";
            		}
        		}
        		$k6=$tempk6;
        		$chr{$temp2[2]}=1;
			}else{
					print "error1...\n";
			}
			if($k6 ne ""){
				push @temp_k6,$k6;
			}
     	
    	 } 
    
        if(@temp_k6==2){
        	$N++;
        	my @chr_key=keys(%chr);
        	if(@chr_key==1){$N_S_T++;}else{$N_T_T++;}
       		my $k5_0_0=substr($temp_k6[0],0,5);
     		my $k5_0_1=substr($temp_k6[0],-5);
     		my $k5_1_0=substr($temp_k6[1],0,5);
     		my $k5_1_1=substr($temp_k6[1],-5);
     		if($temp_k6[0] eq $temp_k6[1] or $k5_0_0 eq $k5_1_1 or $k5_0_1 eq $k5_1_0){
     			$Ends_equal_k6++;
     			if(exists($K6_motif{$k6})){
        			$K6_motif{$k6}++;
    			}else{
       		 		$K6_motif{$k6}=1; 
   				}
     		}
     	}
	}
	open(OUT,">>$out")or die $!;
	print OUT "$file\t$N\t".($N_S_T/$N)."\t".($N_T_T/$N)."\t". ($Ends_equal_k6/$N)."\n";
}
close();

my $kmer="";
my @k6mer=sort{$a<=>$b}keys(%K6_motif);
foreach(@k6mer){   
   	$kmer.=$_."\t".$K6_motif{$_}."\n"; 
}
my $out="k6mer_equal_discordant_".$list;
open(OUTT,">$out")or die $!;
print OUTT $kmer;
close(OUTT);
	
