#!/usr/bin/perl -w
BEGIN{  push (@INC, '/home/fredjiang/fredProgs/PMs')  };


use strict;
use warnings;



############################################################################################################
#

#         
############################################################################################################

use GD;
use Cwd;
use File::Basename;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::SeqFeature::Lite;
use TimeWork;
use Interproscan;
use Bio::Graphics::Panel;
use DirFileHandle;
use OnlyOnePidWork;


use Bio::DB::EUtilities;


                      
package EUtilitieswork;

my $JL_ProteinGB_DB_efatchproteinGBworking_Flag="/home/fredjiang/EightT/fredjiang.2018.04.02/GB_database20181213/EfatchproteinGB_WorkingFlag.20190326.txt";
my $Longest_waiting_time='1h';
my $largst_number_onetimeFeach_limit=5000;

sub efatch_proteinGB {  #EUtilitieswork::efatch_proteinGB ($proteinID_or_Array, $outFILE);
	my ($proteinID_or_Array, $outFILE)=@_;
  
  my $warnMsgBody="\nIn package  EUtilitieswork,\tIn sub efatch_proteinGB,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $subCalInfom=DirFileHandle::print_SubCallerInform;
  
  my $number_limit=$largst_number_onetimeFeach_limit;
  
  my $nowWorkTimePidInformMark=OnlyOnePidWork::CheckAndStart_step_forOnlyOnePidWork($JL_ProteinGB_DB_efatchproteinGBworking_Flag, $Longest_waiting_time);	
  
  if (   ( defined ( $nowWorkTimePidInformMark )  ) && ( $nowWorkTimePidInformMark=~/\S+/ )   ){
  	
  	
  	if (   (  defined ( $outFILE )  ) && ( $outFILE=~m/\S+/)   ){			}
	  else{
	  	DieWork::Just_dieWork( $die_MsgHead."\n\$outFILE=$outFILE is not right!!!: $!".$subCalInfom );
	  }
	  
	  my @proteinID_Array;
	  if (   (  defined ( $proteinID_or_Array )  ) &&  ( $proteinID_or_Array=~m/\S+/)   ){
	  	if (  ref ( $proteinID_or_Array ) eq 'ARRAY'  ){
	  		@proteinID_Array=@{ $proteinID_or_Array };
	  	}
	  	else {
	  		@proteinID_Array=( $proteinID_or_Array );
	  	}
	  }
	  my $new_proID_array=[ @proteinID_Array ];
	  
	  my $proteinID_array_SIZE=@proteinID_Array;
	  
	  if ($proteinID_array_SIZE > $number_limit){
	  	my $Segmented_array=ArrayHashChange::Segment_bigArray_into_smallArrayes($new_proID_array, $number_limit); 
	  	if (   (  defined ( $Segmented_array )  ) &&  (  ref ( $Segmented_array ) eq 'ARRAY'  )   ){
	  		
	  		open (IN, ">$outFILE") or DieWork::Just_dieWork( $die_MsgHead."\n\cannot create \$outFILE=$outFILE !!!: $!".$subCalInfom );
	  		for (  my $i=0; $i<@{ $Segmented_array }; $i++  ){
	  			if (   (  defined ( $Segmented_array->[$i] )  ) &&  (  ref ( $Segmented_array->[$i] ) eq 'ARRAY'  )   ){
	  				my $partFileOut=$outFILE.".part.".$i.".txt";
	  				EUtilitieswork::efatch_proteinGB_inOnlyOnePidwork ($Segmented_array->[$i], $partFileOut); 
	  				my $temp_String=InFileHandle::readAllfileIntoAstring($partFileOut); 
	  				print IN $temp_String;
	  			}
	  			else{
	  		    DieWork::Just_dieWork( $die_MsgHead."\n\$Segmented_array->[$i]=$Segmented_array->[$i] should be a ARRAY ref!!!: $!".$subCalInfom );
	  	    }
	  		}
	  		close (IN);
	  	}
	  	else{
	  		DieWork::Just_dieWork( $die_MsgHead."\n\$Segmented_array=$Segmented_array should be a ARRAY ref!!!: $!".$subCalInfom );
	  	}
	  }
	  else{
	  	EUtilitieswork::efatch_proteinGB_inOnlyOnePidwork ($proteinID_or_Array, $outFILE);  	
	  }
	  
  	
  	
  	
  	
  	OnlyOnePidWork::CheckAnd_Done_step_forOnlyOnePidWork( $JL_ProteinGB_DB_efatchproteinGBworking_Flag, $nowWorkTimePidInformMark );  	
  
  }
  
}


sub efatch_proteinGB_inOnlyOnePidwork {  #EUtilitieswork::efatch_proteinGB_inOnlyOnePidwork ($proteinID_or_Array, $outFILE);
	my ($proteinID_or_Array, $outFILE)=@_;
	
	my $warnMsgBody="\nIn package  EUtilitieswork,\tIn sub efatch_proteinGB_inOnlyOnePidwork,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $subCalInfom=DirFileHandle::print_SubCallerInform;
	
	
	
	my $reTryMaxTime=10;
	my $sleepTime=5;
	
	if (   (  defined ( $outFILE )  ) && ( $outFILE=~m/\S+/)   ){			}
	else{
		DieWork::Just_dieWork( $die_MsgHead."\n\$outFILE=$outFILE is not right!!!: $!".$subCalInfom );
	}
	
	my @proteinID_Array;
	if (   (  defined ( $proteinID_or_Array )  ) &&  ( $proteinID_or_Array=~m/\S+/)   ){
		if (  ref ( $proteinID_or_Array ) eq 'ARRAY'  ){
			@proteinID_Array=@{ $proteinID_or_Array };
		}
		else {
			@proteinID_Array=( $proteinID_or_Array );
		}
	}
	       
	my $factory = Bio::DB::EUtilities->new(-eutil   => 'efetch',
                                         -db      => 'protein',
                                         -id      => \@proteinID_Array,
                                         -email   => 'fredjiang240@126.com',
                                         #-rettype => 'acc'
                                         -rettype => 'gbwithparts'
                                         
                                         );
                                         
  if ( -e ( $outFILE ) ){ system ("rm -f $outFILE"); }
                                           
  my $Work_done_flag=0;
  FORMK: for ( my $i=0; $i<$reTryMaxTime; $i++ ) {                                  
    my $pidKey=fork();    
	  if (  !defined (  $pidKey )  ) {                  DieWork::Just_dieWork( $die_MsgHead."\n Error in fork: $!".$subCalInfom );       }
        
    if ($pidKey == 0) {    warn "$warnMsgBody\n  Child   fork: My pid = $$\t\t \$outFILE=$outFILE\n \$i=$i\n\n";  print "$warnMsgBody\n Child   fork: My pid = $$\t\t \$outFILE=$outFILE\n\n";
                                        
      $factory->get_Response(-file => $outFILE);
      exit 0;
    } 
    
    waitpid($pidKey, 0);
    
    if ( -e ( $outFILE ) ){
    	$Work_done_flag=1;
    	last FORMK;    
    }
    else{
    	my $ciShu=$i+1;
    	my $sleep_Seconds=$sleepTime*$ciShu;
    	my $warnMsg="\n\n\$outFILE=$outFILE\n Try to efatch the $ciShu time!!\nSleep for $sleep_Seconds seconds!!!\n\n\n";
    	warn  $warnMsg;
    	print $warnMsg;
    	sleep( $sleep_Seconds );
    }
  } 
  if ( $Work_done_flag==0 ){
  	DieWork::Just_dieWork( $die_MsgHead."\n\$outFILE=$outFILE maybe cannot  connect to eutils.ncbi.nlm.nih.gov !!!: $!".$subCalInfom );
  }   	
  
  
}

sub efatch_proteinGB_Insub{  #EUtilitieswork::efatch_proteinGB_Insub ($proteinID_or_Array, $outFILE);
	my ($proteinID_or_Array, $outFILE)=@_;
	
	my $warnMsgBody="\nIn package  EUtilitieswork,\tIn sub efatch_proteinGB_Insub,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  DirFileHandle::print_SubCallerInform;
	
	my @proteinID_Array;
	if (   (  defined ( $proteinID_or_Array )  ) &&  ( $proteinID_or_Array=~m/\S+/)   ){
		if (  ref ( $proteinID_or_Array ) eq 'ARRAY'  ){
			@proteinID_Array=@{ $proteinID_or_Array };
		}
		else {
			@proteinID_Array=( $proteinID_or_Array );
		}
	}
	       
	my $factory = Bio::DB::EUtilities->new(-eutil   => 'efetch',
                                         -db      => 'protein',
                                         -id      => \@proteinID_Array,
                                         -email   => 'fredjiang240@126.com',
                                         #-rettype => 'acc'
                                         -rettype => 'gbwithparts'
                                         
                                         );
  $factory->get_Response(-file => $outFILE);
}


sub test0{
  my @ids = qw(1621261 89318838 68536103 20807972 730439);
  
  my $factory = Bio::DB::EUtilities->new(-eutil   => 'efetch',
                                         -db      => 'protein',
                                         -id      => \@ids,
                                         -email   => 'fredjiang240@123.com',
                                         #-rettype => 'acc'
                                         -rettype => 'gbwithparts'
                                         
                                         );
  
  #my @accs = split(m{\n},$factory->get_Response->content);
  
  #print join(',',@accs), "\n";
  
  my $file = 'myseqs.gb';
  
  # dump <HTTP::Response> content to a file (not retained in memory)
  $factory->get_Response(-file => $file);
  
  my $seqin = Bio::SeqIO->new(-file => $file,
                              -format => 'genbank');
  
  while (my $seq = $seqin->next_seq) {
     my $seqence= $seq->();
     print "\$seqence=$seqence\n";
     # do whatever....
  }
}




sub test1{
my $id = 27479347;

  # No sequence, just CONTIG
  my $factory = Bio::DB::EUtilities->new(-eutil => 'efetch',
                                         -db      => 'nucleotide',
                                         -id      => $id,
                                         -email   => 'mymail@foo.bar',
                                         -rettype => 'gb');
  
  $factory->get_Response(-file => 'contigfile.gb');
  
  # Get sequence and all features
  $factory->set_parameters(-rettype => 'gbwithparts');
  
  # file with sequence
  $factory->get_Response(-file => 'full_contig.gb');
}  
                                         
sub test2{
	
	my @ids = qw(CAB02640 EAS10332 YP_250808 NP_623143 P41007);

  my $factory = Bio::DB::EUtilities->new(-eutil => 'efetch',
                                       -db      => 'protein',
                                       -id      => \@ids,
                                       -email   => 'mymail@foo.bar',
                                       -rettype => 'gi');

  my @gis = split(m{\n},$factory->get_Response->content);

  print join(',',@gis), "\n";

}            


1;