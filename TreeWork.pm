#!/usr/bin/perl -w
BEGIN{  push (@INC, '/home/fredjiang/fredProgs/PMs')  };

use strict;
use warnings;
use Bio::TreeIO;
use DieWork;
use InFileHandle;

                      
package TreeWork;


##################################################################################################
#
#   
#
#
##################################################################################################


#TreeWork::Get_nodeNB_to_orgSeqNAME_HASH_nexus($in_nex_file);
sub Get_nodeNB_to_orgSeqNAME_HASH_nexus{
	my ($in_nex_file)=@_;
	
	my $warnMsgBody="\nIn package  TreeWork,\tIn sub Get_nodeNB_to_orgSeqNAME_HASH_nexus,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	DieWork::Check_FileDirExist_or_DIE( $in_nex_file, "\$in_nex_file", $die_MsgHead, $caller_inform  );
	my $outString=InFileHandle::readAllfileIntoAstring($in_nex_file);
	
	
	if ( $outString=~m/
    	                         begin\s+trees;\n
    	                         \s+translate\n
                               (                                  
                                 (?:\s+\'\d+\'\s+\'\S+\',\n)*   #'751'	'lcl|JLDATANB_5230',
                                 (?:\s+\'\d+\'\s+\'\S+\'\n)
                               #   
                               )
                               \s+;  
                         
                            /x
	    )
	{ #DieWork::Print_and_warn( "\$1=$1 \$2=$2\n" );  
		my $mainString=$1;                                                     
    my @node_to_name_line_array=($mainString=~m/                       
                                                               
                                                  (\s+\'\d+\'\s+\'\S+\')   #'751'	'lcl|JLDATANB_5230',
                                                              
                                                                      
  	                                  	       /xg 
  	                             );
  	foreach my $each_node_to_name_line (@node_to_name_line_array) {   
  	  if ( $each_node_to_name_line=~ m/^\s+\'(\d+)\'\s+\'(\S+)\'/ ){
  	  	DieWork::Print_and_warn( "\$1=$1 \$2=$2\n" );
  	  }
  	  else{
  	  	DieWork::Just_dieWork( $die_MsgHead."\n \$each_node_to_name_line=$each_node_to_name_line should be fit the regular expression !!  $!\n\n\n".$caller_inform );	
  	  }
  	}
	}
	else{
		DieWork::Just_dieWork( $die_MsgHead."\n \$outString=\$outString should be fit the regular expression !!  $!\n\n\n".$caller_inform );	
	}
	
} 
 


1;