
#!/usr/bin/perl -w

BEGIN{  push (@INC, '/home/fredjiang/fredProgs/PMs')  };

use strict;
use warnings;

use Storable ;
use Storable qw(nstore store_fd nstore_fd freeze thaw dclone);

use PrintSubArrayHash;
use FastaFileHandle;
use SeqSegmentsTools;
use BlastHandle;
use DirFileHandle;
use DieWork;
use InFileHandle;
use GenBankHandle;
use EUtilitieswork;
package  DiamondWork_Find_CHAR;



#statistic the specific char presentation number in the diamond out hash file.

sub Statistics_of_Diamond_chunkRst_fromFILE{ #my $outHash_suit_forStat=DiamondWork_Find_CHAR::Statistics_of_Diamond_chunkRst_fromFILE( $dmd_chunk_hash_file );
	my ($dmd_chunk_hash_file, $dmd_db_fasta_idx_File, $cutOff_hitNb, $cutOff_Evalue, $cutOff_LcScore )=@_;
	
	my $warnMsgBody="\nIn package  DiamondWork_Find_CHAR,\tIn sub Statistics_of_Diamond_chunkRst_fromFILE,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	if   (   (  defined ( $dmd_chunk_hash_file )  ) && ( $dmd_chunk_hash_file=~m/\S+/ ) && (  ( -e $dmd_chunk_hash_file )  )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$dmd_chunk_hash_file=$dmd_chunk_hash_file should be a defined not empty file  !!  $!\n\n\n".$caller_inform ); 	}
	warn "\n\n\$dmd_chunk_hash_file=$dmd_chunk_hash_file\n";
	my $dmd_chunk_hash_Hash=Storable::retrieve( $dmd_chunk_hash_file );
	my $outHash_suit_forStat=DiamondWork_Find_CHAR::Build_Hash_for_Statistics_of_Diamond_chunkRst( $dmd_chunk_hash_Hash, $dmd_db_fasta_idx_File  );
	my $Sorted_hit_HASH     =DiamondWork_Find_CHAR::Sort_hitGroup_for_each_CodonPositon          ( $outHash_suit_forStat                         );
	my $out_statistic_hash  =DiamondWork_Find_CHAR::Doing_Statistics_of_Diamond_chunkRst         ( $Sorted_hit_HASH, $cutOff_hitNb, $cutOff_Evalue, $cutOff_LcScore  );
	return [ $Sorted_hit_HASH, $out_statistic_hash ];
	  
	
}

sub Build_Hash_for_Statistics_of_Diamond_chunkRst{  #my $outHash_suit_forStat=DiamondWork_Find_CHAR::Build_Hash_for_Statistics_of_Diamond_chunkRst( $dmd_chunk_hash );
	my ( $dmd_chunk_hash, $dmd_db_fasta_idx_File )=@_;
	
	my $warnMsgBody="\nIn package  DiamondWork_Find_CHAR,\tIn sub Build_Hash_for_Statistics_of_Diamond_chunkRst,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	if  (   (  defined ( $dmd_chunk_hash )  ) && (  ref ( $dmd_chunk_hash ) eq 'HASH'  )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$dmd_chunk_hash=$dmd_chunk_hash should be a HASH  !!  $!\n\n\n".$caller_inform ); 	}
	
	my $reGetThePep_Char_by_index_orNot=0;
	if   (   (  defined ( $dmd_db_fasta_idx_File )  ) && ( $dmd_db_fasta_idx_File=~m/\S+/ )   ){
	  if ( -e $dmd_db_fasta_idx_File ){
	  	$reGetThePep_Char_by_index_orNot=1;
	  }
	  else{		DieWork::Just_dieWork( $die_MsgHead."\n \$dmd_db_fasta_idx_File=$dmd_db_fasta_idx_File should be a defined not empty file  !!  $!\n\n\n".$caller_inform ); 	}
	}
	
	
	
	my $stp_rst_Key_hash=ArrayHashChange::Change_2d_hashINTO_3d_hash ($dmd_chunk_hash, '0_2_a_DmdMrStpHASH');  #'0_2_a_DmdMrStpHASH' Õâ¸ö³£Á¿×Ö·û´® ÔÚDiamond_multipleCore_working.plÖÐ¶¨Òå
  if  (   (  defined ( $stp_rst_Key_hash )  ) && (  ref ( $stp_rst_Key_hash ) eq 'HASH'  )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$stp_rst_Key_hash=$stp_rst_Key_hash should be a HASH  !!  $!\n\n\n".$caller_inform ); 	}
	
	
	my $outHash_suit_forStat;          #the out put hash , holding all the information to show the stop codon readthrough events
	my $tempHash_to_del_repeat_result; #an temp hash ,holding same postion 3-DNA-codon from different hits, the best hit will fit into this temp hash
	foreach my $eachStpRstFile (    sort {$a cmp $b} (   keys (  %{ $stp_rst_Key_hash }  )   )    ){  DieWork::Print_and_warn( "\n20191009-0-0-0 \$eachStpRstFile=\$eachStpRstFile\n" );
		
		if   (   (  defined ( $eachStpRstFile )  ) && ( $eachStpRstFile=~m/\S+/ ) && (  ( -e $eachStpRstFile )  )   ){}																																																																								#DieWork::Print_and_warn( "\n20191009-0-0-0-0 \$eachStpRstFile=$eachStpRstFile\n" ); 
	  else{		 DieWork::Just_dieWork( $die_MsgHead."\n \$eachStpRstFile=$eachStpRstFile should be a defined not empty file  !!  $!\n\n\n".$caller_inform ); 	}
	  my $eachStpRst_HASH=Storable::retrieve( $eachStpRstFile );                           																																																																								 #DieWork::Print_and_warn( "\n20191009-0-0-0-1 \$eachStpRst_HASH=$eachStpRst_HASH\n" );
	  if  (   (  defined ( $eachStpRst_HASH )  ) && (  ref ( $eachStpRst_HASH ) eq 'HASH'  )   ){}
	  else{		DieWork::Just_dieWork( $die_MsgHead."\n \$eachStpRst_HASH=$eachStpRst_HASH should be a HASH  !!  $!\n\n\n".$caller_inform ); 	}
		
		if  (   (  defined ( $eachStpRst_HASH->{'0_0_1_Query_array'} )  ) && (  ref ( $eachStpRst_HASH->{'0_0_1_Query_array'} ) eq 'ARRAY'  )   ){  																																																																								#DieWork::Print_and_warn( "\n20191009-0-0-0-2 \$eachStpRst_HASH->{'0_0_1_Query_array'}=$eachStpRst_HASH->{'0_0_1_Query_array'}\n" );
	    for (  my $qry_idx=0; $qry_idx < @{ $eachStpRst_HASH->{'0_0_1_Query_array'} }; $qry_idx++  ) {     
	    	my $ovlCk_name=$eachStpRst_HASH->{'0_0_1_Query_array'}->[$qry_idx]->{'2_0_0_0_query_name'};      DieWork::Print_and_warn( "\n20191009-0-0-1 \$ovlCk_name=\$eachStpRst_HASH->{'0_0_1_Query_array'}->[$qry_idx]->{'2_0_0_0_query_name'}=$ovlCk_name\n" );
	    	
	    	if  (   (  defined ( $eachStpRst_HASH->{'0_0_1_Query_array'}->[$qry_idx]->{'2_0_Z_hitArray'} )  ) && (  ref ( $eachStpRst_HASH->{'0_0_1_Query_array'}->[$qry_idx]->{'2_0_Z_hitArray'} ) eq 'ARRAY'  )   ){
	    	  for (  my $hit_idx=0; $hit_idx < @{ $eachStpRst_HASH->{'0_0_1_Query_array'}->[$qry_idx]->{'2_0_Z_hitArray'} }; $hit_idx++  ) {
	    	  my $hit___name=$eachStpRst_HASH->{'0_0_1_Query_array'}->[$qry_idx]->{'2_0_Z_hitArray'}->[$hit_idx]->{'3_0_1__name'};  DieWork::Print_and_warn( "\n20191009-0-0-2 \$hit___name=\$eachStpRst_HASH->{'0_0_1_Query_array'}->[$qry_idx]->{'2_0_Z_hitArray'}->[$hit_idx]->{'3_0_1__name'}=$hit___name\n" );
	    		    	 
	    	  
	    	    if  (      (  defined ( $eachStpRst_HASH->{'0_0_1_Query_array'}->[$qry_idx]->{'2_0_Z_hitArray'}->[$hit_idx]->{'4_0_0_hspArray'} )             ) 
	    	            && (  ref     ( $eachStpRst_HASH->{'0_0_1_Query_array'}->[$qry_idx]->{'2_0_Z_hitArray'}->[$hit_idx]->{'4_0_0_hspArray'} ) eq 'ARRAY'  )      ) 
	    	    {
	    	      for (  my $hsp_idx=0; $hsp_idx < @{ $eachStpRst_HASH->{'0_0_1_Query_array'}->[$qry_idx]->{'2_0_Z_hitArray'}->[$hit_idx]->{'4_0_0_hspArray'} }; $hsp_idx++  ) {  DieWork::Print_and_warn( "\n20191009-0-0-3 \$hsp_idx=$hsp_idx\n" );
	    	
	    	      	
	    	      	
	    	      	if  (      (  defined ( $eachStpRst_HASH->{'0_0_1_Query_array'}->[$qry_idx]->{'2_0_Z_hitArray'}->[$hit_idx]->{'4_0_0_hspArray'}->[$hsp_idx]->{'5_6_0_Hsp_AA_match'} )            ) 
	    	                && (  ref     ( $eachStpRst_HASH->{'0_0_1_Query_array'}->[$qry_idx]->{'2_0_Z_hitArray'}->[$hit_idx]->{'4_0_0_hspArray'}->[$hsp_idx]->{'5_6_0_Hsp_AA_match'} ) eq 'HASH'  )      ) 
	    	        {
	    	      	   foreach my $eachCharKey (    sort {$a cmp $b} (   keys (  %{ $eachStpRst_HASH->{'0_0_1_Query_array'}->[$qry_idx]->{'2_0_Z_hitArray'}->[$hit_idx]->{'4_0_0_hspArray'}->[$hsp_idx]->{'5_6_0_Hsp_AA_match'} }  )   )    ){  DieWork::Print_and_warn( "\n20191009-0-0-3 \$eachCharKey=$eachCharKey\n" );
	    	      	   	 
	    	      	   	 if  (      (  defined ( $eachStpRst_HASH->{'0_0_1_Query_array'}->[$qry_idx]->{'2_0_Z_hitArray'}->[$hit_idx]->{'4_0_0_hspArray'}->[$hsp_idx]->{'5_6_0_Hsp_AA_match'}->{$eachCharKey} )            ) 
	    	                     && (  ref     ( $eachStpRst_HASH->{'0_0_1_Query_array'}->[$qry_idx]->{'2_0_Z_hitArray'}->[$hit_idx]->{'4_0_0_hspArray'}->[$hsp_idx]->{'5_6_0_Hsp_AA_match'}->{$eachCharKey} ) eq 'ARRAY'  )      ) 
	    	             {
	    	      	   	   for (  my $chr_idx=0; $chr_idx < @{ $eachStpRst_HASH->{'0_0_1_Query_array'}->[$qry_idx]->{'2_0_Z_hitArray'}->[$hit_idx]->{'4_0_0_hspArray'}->[$hsp_idx]->{'5_6_0_Hsp_AA_match'}->{$eachCharKey} }; $chr_idx++  ) {  DieWork::Print_and_warn( "\n20191009-0-0-4 \$chr_idx=$chr_idx\n" );
	    	      	   	     
	    	      	   	     if  (      (  defined ( $eachStpRst_HASH->{'0_0_1_Query_array'}->[$qry_idx]->{'2_0_Z_hitArray'}->[$hit_idx]->{'4_0_0_hspArray'}->[$hsp_idx]->{'5_6_0_Hsp_AA_match'}->{$eachCharKey}->[$chr_idx] )            ) 
	    	                         && (  ref     ( $eachStpRst_HASH->{'0_0_1_Query_array'}->[$qry_idx]->{'2_0_Z_hitArray'}->[$hit_idx]->{'4_0_0_hspArray'}->[$hsp_idx]->{'5_6_0_Hsp_AA_match'}->{$eachCharKey}->[$chr_idx] ) eq 'HASH'  )      ) 
	    	                 {
	    	                 	 
	    	                 	 
	    	                 	 #begin to build the statistic hash
	    	                 	 
	    	                 	 my $tempShortChar_HASH=$eachStpRst_HASH->{'0_0_1_Query_array'}->[$qry_idx]->{'2_0_Z_hitArray'}->[$hit_idx]->{'4_0_0_hspArray'}->[$hsp_idx]->{'5_6_0_Hsp_AA_match'}->{$eachCharKey}->[$chr_idx];
	    	                 	 
	    	                 	 my $genome_file=$dmd_chunk_hash->{$ovlCk_name}->{'0_0_0_orgFile_path'}; #the genomic file path, to reprensent the organism
	    	                 	 my $genome___id=$dmd_chunk_hash->{$ovlCk_name}->{'0_0_2_orgParimayID'}; #DNA file ID
	    	                 	 my $AA_DNA_type=$tempShortChar_HASH->{'6_0_3_Qur_DNAchar'};             #STOP codon type: TAA TAG TGA
	    	                 	 
	    	                 	 my $MatchAA_Pos=$tempShortChar_HASH->{'6_1_1_hit_pos_stt'};             #Hit postion start
	    	                 	 my $MatchPeP_AA=$tempShortChar_HASH->{'6_1_3_hit_PEPchar'};             #Hit protein amino acid residue
	    	                 	 ## if we need reget the PEP char , then do the work below
	    	                 	 if ( $reGetThePep_Char_by_index_orNot==1){
	    	                 	 	 $MatchPeP_AA= FastaFileHandle::feach_seqString_from_idx_File($hit___name, $dmd_db_fasta_idx_File, $MatchAA_Pos, $MatchAA_Pos);
	    	                 	 }
	    	                 	 
	    	                 	 
	    	                 	 my $realSttPost=$dmd_chunk_hash->{$ovlCk_name}->{'0_3_4_Sgm______stt'}-1+$tempShortChar_HASH->{'6_0_1_Qur_pos_stt'};
	    	                 	 my $realEndPost=$dmd_chunk_hash->{$ovlCk_name}->{'0_3_4_Sgm______stt'}-1+$tempShortChar_HASH->{'6_0_2_Qur_pos_end'};
	    	                 	 my $Pos_Stn_KEY=$tempShortChar_HASH->{'6_0_0_Qur_stn_ZoF'}."_".$realSttPost."_".$realEndPost;
	    	                 	                                                                         #postion 
	    	                 	                                                                         
                           my $Real_querySttPos=$dmd_chunk_hash->{$ovlCk_name}->{'0_3_4_Sgm______stt'}-1+$eachStpRst_HASH->{'0_0_1_Query_array'}->[$qry_idx]->{'2_0_Z_hitArray'}->[$hit_idx]->{'4_0_0_hspArray'}->[$hsp_idx]->{'5_2_0_Hsp_2StEd_1Qr_1St'};
                           my $Real_queryEndPos=$dmd_chunk_hash->{$ovlCk_name}->{'0_3_4_Sgm______stt'}-1+$eachStpRst_HASH->{'0_0_1_Query_array'}->[$qry_idx]->{'2_0_Z_hitArray'}->[$hit_idx]->{'4_0_0_hspArray'}->[$hsp_idx]->{'5_2_1_Hsp_2StEd_1Qr_2Ed'};
                           
                           	    	                 	                                                                         
	    	                 	 
	    	                 	 my $Hit_nameKEY=$hit___name;                                            #hit id
	    	                 	 
	    	                 	 DieWork::Print_and_warn( "\n20191009-0-0-4 \$Pos_Stn_KEY=\$Pos_Stn_KEY\n" );
	    	                 	 
	    	                 	 #detail information hash
	    	                 	 my $final_point_hash;
	    	                 	 
	    	                 	 $final_point_hash->{'3_0_5_hitDescription'   }=$eachStpRst_HASH->{'0_0_1_Query_array'}->[$qry_idx]->{'2_0_Z_hitArray'}->[$hit_idx]->{'3_0_5_hitDescription'   };
	    	                 	 $final_point_hash->{'3_0_6_hitEvalue'        }=$eachStpRst_HASH->{'0_0_1_Query_array'}->[$qry_idx]->{'2_0_Z_hitArray'}->[$hit_idx]->{'3_0_6_hitEvalue'        };
	    	                 	 $final_point_hash->{'3_0_7_hitScore'         }=$eachStpRst_HASH->{'0_0_1_Query_array'}->[$qry_idx]->{'2_0_Z_hitArray'}->[$hit_idx]->{'3_0_7_hitScore'         };
	    	                 	 $final_point_hash->{'3_2_0_hitTotalIdentical'}=$eachStpRst_HASH->{'0_0_1_Query_array'}->[$qry_idx]->{'2_0_Z_hitArray'}->[$hit_idx]->{'3_2_0_hitTotalIdentical'};
	    	                 	 $final_point_hash->{'3_2_1_hitTotalConserved'}=$eachStpRst_HASH->{'0_0_1_Query_array'}->[$qry_idx]->{'2_0_Z_hitArray'}->[$hit_idx]->{'3_2_1_hitTotalConserved'};
	    	                 	 
	    	                 	 $final_point_hash->{'5_2_0_0_Hsp_qrREALstt'  }=$Real_querySttPos;
	    	                 	 $final_point_hash->{'5_2_1_0_Hsp_qrREALend'  }=$Real_queryEndPos;
	    	                 	 
	    	                 	 
	    	                 	 $final_point_hash->{'6_1_3_hit_PEPchar'      }=$MatchPeP_AA;
	    	                 	 
	    	                 	 $final_point_hash->{'6_2_0_locAlnScore'      }=$tempShortChar_HASH->{'6_2_0_locAlnScore'};
	    	                 	 $final_point_hash->{'6_2_1_ChrShow_Str'      }=$tempShortChar_HASH->{'6_2_1_ChrShow_Str'};
	    	                 	 $final_point_hash->{'6_2_2_Qur_sub_Str'      }=$tempShortChar_HASH->{'6_2_2_Qur_sub_Str'};
	    	                 	 $final_point_hash->{'6_2_3_Mch_sub_Str'      }=$tempShortChar_HASH->{'6_2_3_Mch_sub_Str'};
	    	                 	 $final_point_hash->{'6_2_4_Shw_Sco_Str'      }=$tempShortChar_HASH->{'6_2_4_Shw_Sco_Str'};
	    	                 	 $final_point_hash->{'6_2_5_Hit_sub_Str'      }=$tempShortChar_HASH->{'6_2_5_Hit_sub_Str'};
	    	                 	 
	    	                 	 $final_point_hash->{'7_1_0_BLThshP'          }=$eachStpRstFile;
	    	                 	 $final_point_hash->{'7_1_1_qry_idx'          }=$qry_idx;
	    	                 	 $final_point_hash->{'7_1_2_hit_idx'          }=$hit_idx;
	    	                 	 $final_point_hash->{'7_1_3_hsp_idx'          }=$hsp_idx;
	    	                 	                                                                     
	    	                 	 
	    	                 	 if (       (  defined ( $tempHash_to_del_repeat_result->{$genome_file}->{$AA_DNA_type}->{$genome___id}->{$Pos_Stn_KEY}->{$Hit_nameKEY}->{'3_0_6_hitEvalue'} )   )
	    	                 	        &&  (  $tempHash_to_del_repeat_result->{$genome_file}->{$AA_DNA_type}->{$genome___id}->{$Pos_Stn_KEY}->{$Hit_nameKEY}->{'3_0_6_hitEvalue'} =~ m/\S+/     )
	    	                 	        &&  (  $tempHash_to_del_repeat_result->{$genome_file}->{$AA_DNA_type}->{$genome___id}->{$Pos_Stn_KEY}->{$Hit_nameKEY}->{'3_0_6_hitEvalue'} < $final_point_hash->{'3_0_6_hitEvalue'}     )
	    	                 	    )
	    	                 	 {
	    	                 	 	 #do nothing # when the postion was recorded, and the recorded hit evalue was smaller , then the current one would not be recorded 
	    	                 	 }
	    	                 	 else{ # when the codintion above was not satisfied, then do the upgrade work
	    	                 	 	 
	    	                 	 	 #delete the old and worse recored hit
	    	                 	 	 if (       (  defined ( $tempHash_to_del_repeat_result->{$genome_file}->{$AA_DNA_type}->{$genome___id}->{$Pos_Stn_KEY}->{$Hit_nameKEY}->{'6_1_3_hit_PEPchar'} )   )
	    	                 	          &&  (  $tempHash_to_del_repeat_result->{$genome_file}->{$AA_DNA_type}->{$genome___id}->{$Pos_Stn_KEY}->{$Hit_nameKEY}->{'6_1_3_hit_PEPchar'} =~ m/\S+/     )
	    	                 	    )
	    	                 	   {
	    	                 	 	   delete (
	    	                 	 	     $outHash_suit_forStat -> {$genome_file}     # the genomic file path, to reprensent the organism # # '/home/fredjiang/work/Algae/20150522NewDATA/AllSpecies/Nannochloropsis_oculata20150522/Est/Nannochloropsis_oculata20150522'
	    	                 	                             -> {$AA_DNA_type}     #STOP codon type: TAA TAG TGA                       # # 'TGA'  »ò 'TAA'  
	    	                 	                             -> {$genome___id}     #DNA sequence ID                                    # # 'scaffold123'  or 'lcl|1ctg0'                         
	    	   	    	                   	                 -> {$Pos_Stn_KEY}     #Postion string                                     # # '+_23334_23890' or  '-_1268_1266'   
	    	                 	                             -> {$Hit_nameKEY}     #hit id                                             # # 'gi|953442847|ref|XP_014564976.1|'
	    	                 	 	   
	    	                 	 	   );
	    	                 	 	 }
	    	                 	 	 
	    	                 	 	 # recorded the current better one
	    	                 	 	 $outHash_suit_forStat -> {$genome_file}     # the genomic file path, to reprensent the organism # # '/home/fredjiang/work/Algae/20150522NewDATA/AllSpecies/Nannochloropsis_oculata20150522/Est/Nannochloropsis_oculata20150522'
	    	                 	                         -> {$AA_DNA_type}     #STOP codon type: TAA TAG TGA                       # # 'TGA'  »ò 'TAA'                                                                                                            
	    	                 	                         -> {$genome___id}     #DNA sequence ID                                    # # 'scaffold123'  or 'lcl|1ctg0'                                                                                              
	    	                 	                         -> {$Pos_Stn_KEY}     #Postion string                                     # # '+_23334_23890' or  '-_1268_1266'                                                                                          
	    	                 	                         -> {$Hit_nameKEY}     #hit id                                             # # 'gi|953442847|ref|XP_014564976.1|'                                                                                         
	    	                 	   =$final_point_hash;                                     
	    	                 	   
	    	                 	   #recorded the best hit of a particalur postion
	    	                 	   $tempHash_to_del_repeat_result 
	    	                 	                         -> {$genome_file}     # the genomic file path, to reprensent the organism # # '/home/fredjiang/work/Algae/20150522NewDATA/AllSpecies/Nannochloropsis_oculata20150522/Est/Nannochloropsis_oculata20150522'
	    	                 	                         -> {$AA_DNA_type}     #STOP codon type: TAA TAG TGA                       # # 'TGA'  »ò 'TAA'                                                                                                            
	    	                 	                         -> {$genome___id}     #DNA sequence ID                                    # # 'scaffold123'  or 'lcl|1ctg0'                                                                                              
	    	                 	                         -> {$Pos_Stn_KEY}     #Postion string                                     # # '+_23334_23890' or  '-_1268_1266'                                                                                          
	    	                 	                         -> {$Hit_nameKEY}     #hit id                                             # # 'gi|953442847|ref|XP_014564976.1|'                                                                                         
	    	                 	   =$final_point_hash;                                     
	    	                 	 }
	    	                 	                      	    	                 	                       	    	                 	                       	    	                 	                                              
	    	                 } 
 	    	      	   	   }
  	    	      	   }	 
	    	      	   }
	    	        }
	    	      }
	    	    
	    	    }
	    	
	    	  
	    	  
	    	  }
	    	}
	    	
	    	
	    }  
	  }
	  
		
	}
	
	return $outHash_suit_forStat;
}

sub Sort_hitGroup_for_each_CodonPositon{  # my $new_Hash= DiamondWork_Find_CHAR::Sort_hitGroup_for_each_CodonPositon( $inHash_suit_forStat );
	my ($inHash_suit_forStat)=@_;
	my ($warnMsgBody, $warnMsgHead, $die_MsgHead, $caller_inform)=@{ DieWork::BuildWarnDieHeadInform( 'DiamondWork_Find_CHAR', 'Sort_hitGroup_for_each_CodonPositon' ) };
  
  DieWork::Check_Hash_or_DIE( $inHash_suit_forStat,     "\$inHash_suit_forStat",    $die_MsgHead, $caller_inform  );  
  
  
  my $temp_hash_for_information_add; #build a temp hash to add information from blast prased hash
  
  my $new_Hash;
  foreach my $fnaPATH (    sort {$a cmp $b}(   keys(  %{ $inHash_suit_forStat }  )   )    ){
  	DieWork::Check_Hash_or_DIE( $inHash_suit_forStat->{$fnaPATH},     "\$inHash_suit_forStat->{\$fnaPATH}=\$inHash_suit_forStat->{$fnaPATH}",    $die_MsgHead, $caller_inform  );  
  	foreach my $CodonTYPE (    sort {$a cmp $b}(   keys(  %{ $inHash_suit_forStat->{$fnaPATH} }  )   )    ){
  		DieWork::Check_Hash_or_DIE( $inHash_suit_forStat->{$fnaPATH}->{$CodonTYPE},     "\$inHash_suit_forStat->{\$fnaPATH}->{\$CodonTYPE}=\$inHash_suit_forStat->{$fnaPATH}->{$CodonTYPE}",    $die_MsgHead, $caller_inform  );  
  		foreach my $GenomicID (    sort {$a cmp $b}(   keys(  %{ $inHash_suit_forStat->{$fnaPATH}->{$CodonTYPE} }  )   )    ){
  			DieWork::Check_Hash_or_DIE( $inHash_suit_forStat->{$fnaPATH}->{$CodonTYPE}->{$GenomicID},     "\$inHash_suit_forStat->{\$fnaPATH}->{\$CodonTYPE}->{\$GenomicID}=\$inHash_suit_forStat->{$fnaPATH}->{$CodonTYPE}->{$GenomicID}",    $die_MsgHead, $caller_inform  );  
  		  foreach my $PositonKEY (    sort {$a cmp $b}(   keys(  %{ $inHash_suit_forStat->{$fnaPATH}->{$CodonTYPE}->{$GenomicID} }  )   )    ){
  		    DieWork::Check_Hash_or_DIE( $inHash_suit_forStat->{$fnaPATH}->{$CodonTYPE}->{$GenomicID}->{$PositonKEY},     "\$inHash_suit_forStat->{\$fnaPATH}->{\$CodonTYPE}->{\$GenomicID}->{\$PositonKEY}=\$inHash_suit_forStat->{$fnaPATH}->{$CodonTYPE}->{$GenomicID}->{$PositonKEY}",    $die_MsgHead, $caller_inform  );  
  		  	my $Pos_mutiHIT_HASH=$inHash_suit_forStat->{$fnaPATH}->{$CodonTYPE}->{$GenomicID}->{$PositonKEY};
  		  	
  		  	my $inGroupSortNUMBER=0;
  		  	my $realPEPresidue;
  		  	foreach my $HitIDkey (    
  		  	                          sort {  (      (  $Pos_mutiHIT_HASH->{$b}->{'6_2_0_locAlnScore'}         <=> $Pos_mutiHIT_HASH->{$a}->{'6_2_0_locAlnScore'}  ) 
  		  	                                     ||  (  $Pos_mutiHIT_HASH->{$a}->{'3_0_6_hitEvalue'}           <=> $Pos_mutiHIT_HASH->{$b}->{'3_0_6_hitEvalue'}    ) 
  		  	                                     ||  (  $Pos_mutiHIT_HASH->{$b}->{'3_2_1_hitTotalConserved'}   <=> $Pos_mutiHIT_HASH->{$a}->{'3_2_1_hitTotalConserved'}    ) 
  		  	                                     ||  (  $Pos_mutiHIT_HASH->{$b}->{'3_2_0_hitTotalIdentical'}   <=> $Pos_mutiHIT_HASH->{$a}->{'3_2_0_hitTotalIdentical'}    ) 
  		  	                                  )
  		  	                               } (   keys(  %{ $inHash_suit_forStat->{$fnaPATH}->{$CodonTYPE}->{$GenomicID}->{$PositonKEY} }  )   )    ){
  		      DieWork::Check_Hash_or_DIE( $inHash_suit_forStat->{$fnaPATH}->{$CodonTYPE}->{$GenomicID}->{$PositonKEY}->{$HitIDkey},     "\$inHash_suit_forStat->{\$fnaPATH}->{\$CodonTYPE}->{\$GenomicID}->{\$PositonKEY}->{\$HitIDkey}=\$inHash_suit_forStat->{$fnaPATH}->{$CodonTYPE}->{$GenomicID}->{$PositonKEY}->{$HitIDkey}",    $die_MsgHead, $caller_inform  );  
  		  	  if ( $inGroupSortNUMBER == 0 ){
  		      	$realPEPresidue=$inHash_suit_forStat->{$fnaPATH}->{$CodonTYPE}->{$GenomicID}->{$PositonKEY}->{$HitIDkey}->{'6_1_3_hit_PEPchar'} 
  		      }
  		      $new_Hash->{$fnaPATH}->{$CodonTYPE}->{$realPEPresidue}->{$GenomicID}->{$PositonKEY}->{$HitIDkey}=Storable::dclone( $inHash_suit_forStat->{$fnaPATH}->{$CodonTYPE}->{$GenomicID}->{$PositonKEY}->{$HitIDkey} );
  		      
  		      $new_Hash->{$fnaPATH}->{$CodonTYPE}->{$realPEPresidue}->{$GenomicID}->{$PositonKEY}->{$HitIDkey}->{'8_0_0_sort_order'}=$inGroupSortNUMBER;
  		      
  		      #add path, query index, hit index
  		      if ( $inGroupSortNUMBER==0){
  		      	my $tempRef=$new_Hash->{$fnaPATH}->{$CodonTYPE}->{$realPEPresidue}->{$GenomicID}->{$PositonKEY}->{$HitIDkey};
  		      	$temp_hash_for_information_add->{ $tempRef->{'7_1_0_BLThshP'} }->{ $tempRef->{'7_1_1_qry_idx'} }->{ $tempRef->{'7_1_2_hit_idx'} }->{ $tempRef->{'7_1_3_hsp_idx'} }
  		      	=[ $fnaPATH, $CodonTYPE, $realPEPresidue, $GenomicID, $PositonKEY, $HitIDkey ];
  		      }
  		      
  		      $inGroupSortNUMBER++;
  		    }
  		    
  		    
  		  }
  		}
  	}
  }
  
  
  #add  more information only for the best hit
  DieWork::Check_Hash_or_DIE( $temp_hash_for_information_add,     "\$temp_hash_for_information_add",    $die_MsgHead, $caller_inform  ); 
  foreach my $BLThshP (    sort {$a cmp $b}(   keys(  %{ $temp_hash_for_information_add }  )   )    ){
  	
  	if   (   (  defined ( $BLThshP )  ) && ( $BLThshP=~m/\S+/ ) && (  ( -e $BLThshP )  )   ){}																																																																								#DieWork::Print_and_warn( "\n20191009-0-0-0-0 \$eachStpRstFile=$eachStpRstFile\n" ); 
	  else{		 DieWork::Just_dieWork( $die_MsgHead."\n \$BLThshP=$BLThshP should be a defined not empty file  !!  $!\n\n\n".$caller_inform ); 	}
	  my $eachStpRst_HASH=Storable::retrieve( $BLThshP );                           																																																																								 #DieWork::Print_and_warn( "\n20191009-0-0-0-1 \$eachStpRst_HASH=$eachStpRst_HASH\n" );
	  if  (   (  defined ( $eachStpRst_HASH )  ) && (  ref ( $eachStpRst_HASH ) eq 'HASH'  )   ){}
	  else{		DieWork::Just_dieWork( $die_MsgHead."\n \$eachStpRst_HASH=$eachStpRst_HASH should be a HASH  !!  $!\n\n\n".$caller_inform ); 	}
  	
  	DieWork::Check_Hash_or_DIE( $temp_hash_for_information_add->{$BLThshP},     "\$temp_hash_for_information_add->{\$BLThshP}=\$temp_hash_for_information_add->{$BLThshP}",    $die_MsgHead, $caller_inform  );  
  	foreach my $qry_idx (    sort {$a <=> $b}(   keys(  %{ $temp_hash_for_information_add->{$BLThshP} }  )   )    ){
  		DieWork::Check_Hash_or_DIE( $temp_hash_for_information_add->{$BLThshP}->{$qry_idx},     "\$temp_hash_for_information_add->{\$BLThshP}->{\$qry_idx}=\$temp_hash_for_information_add->{$BLThshP}->{$qry_idx}",    $die_MsgHead, $caller_inform  );  
  		foreach my $hit_idx (    sort {$a <=> $b}(   keys(  %{ $temp_hash_for_information_add->{$BLThshP}->{$qry_idx} }  )   )    ){
  	    DieWork::Check_Hash_or_DIE( $temp_hash_for_information_add->{$BLThshP}->{$qry_idx}->{$hit_idx},     "\$temp_hash_for_information_add->{\$BLThshP}->{\$qry_idx}->{\$hit_idx}=\$temp_hash_for_information_add->{$BLThshP}->{$qry_idx}->{$hit_idx}",    $die_MsgHead, $caller_inform  );  
  		  foreach my $hsp_idx (    sort {$a <=> $b}(   keys(  %{ $temp_hash_for_information_add->{$BLThshP}->{$qry_idx}->{$hit_idx} }  )   )    ){
  	      DieWork::Check_Array_or_DIE( $temp_hash_for_information_add->{$BLThshP}->{$qry_idx}->{$hit_idx}->{$hsp_idx},     "\$temp_hash_for_information_add->{\$BLThshP}->{\$qry_idx}->{\$hit_idx}->{$hsp_idx}=\$temp_hash_for_information_add->{$BLThshP}->{$qry_idx}->{$hit_idx}->{$hsp_idx}",    $die_MsgHead, $caller_inform  );  
  		 
  		    my ( $fnaPATH_2, $CodonTYPE_2, $realPEPresidue_2, $GenomicID_2, $PositonKEY_2, $HitIDkey_2  )= @{ $temp_hash_for_information_add->{$BLThshP}->{$qry_idx}->{$hit_idx}->{$hsp_idx} } ;
  		    my $tempHASHhere=$eachStpRst_HASH->{'0_0_1_Query_array'}->[$qry_idx]->{'2_0_Z_hitArray'}->[$hit_idx]->{'4_0_0_hspArray'}->[$hsp_idx];
  		    
  		    $new_Hash->{$fnaPATH_2}->{$CodonTYPE_2}->{$realPEPresidue_2}->{$GenomicID_2}->{$PositonKEY_2}->{$HitIDkey_2}->{'5_4_0_Hsp_3WhSting_0AA'}=$tempHASHhere->{'5_4_0_Hsp_3WhSting_0AA'};
  		    $new_Hash->{$fnaPATH_2}->{$CodonTYPE_2}->{$realPEPresidue_2}->{$GenomicID_2}->{$PositonKEY_2}->{$HitIDkey_2}->{'5_4_1_Hsp_3WhSting_1Qr'}=$tempHASHhere->{'5_4_1_Hsp_3WhSting_1Qr'};
  		    $new_Hash->{$fnaPATH_2}->{$CodonTYPE_2}->{$realPEPresidue_2}->{$GenomicID_2}->{$PositonKEY_2}->{$HitIDkey_2}->{'5_4_2_Hsp_3WhSting_2Ho'}=$tempHASHhere->{'5_4_2_Hsp_3WhSting_2Ho'};
  		    $new_Hash->{$fnaPATH_2}->{$CodonTYPE_2}->{$realPEPresidue_2}->{$GenomicID_2}->{$PositonKEY_2}->{$HitIDkey_2}->{'5_4_3_Hsp_3WhSting_3Ht'}=$tempHASHhere->{'5_4_3_Hsp_3WhSting_3Ht'};
  		    
  		    
  		    #expand orf
  		    my $fnaPath_idx=$fnaPATH_2."Index.txt"; my $zhengfu=$PositonKEY_2; $zhengfu=~s/^(\S).*/$1/;  
  		    my $smOrfStt=$new_Hash->{$fnaPATH_2}->{$CodonTYPE_2}->{$realPEPresidue_2}->{$GenomicID_2}->{$PositonKEY_2}->{$HitIDkey_2}->{'5_2_0_0_Hsp_qrREALstt'};
  		    my $smOrfEnd=$new_Hash->{$fnaPATH_2}->{$CodonTYPE_2}->{$realPEPresidue_2}->{$GenomicID_2}->{$PositonKEY_2}->{$HitIDkey_2}->{'5_2_1_0_Hsp_qrREALend'};
  		    my $expandWithOutZF_hash=FastaFileHandle::ExpandORFforBacteriaGENOME  ($fnaPath_idx, $GenomicID_2, $smOrfStt, $smOrfEnd,  $zhengfu);
  		    $new_Hash->{$fnaPATH_2}->{$CodonTYPE_2}->{$realPEPresidue_2}->{$GenomicID_2}->{$PositonKEY_2}->{$HitIDkey_2}->{'9_0_z_expand_orf_hash'}=$expandWithOutZF_hash;
  		    #$new_Hash->{$fnaPATH_2}->{$CodonTYPE_2}->{$realPEPresidue_2}->{$GenomicID_2}->{$PositonKEY_2}->{$HitIDkey_2}->{'9_0_7_1_new_0_stt'}=$expandWithOutZF_hash->{'7_1_new_0_stt'};
  		    #$new_Hash->{$fnaPATH_2}->{$CodonTYPE_2}->{$realPEPresidue_2}->{$GenomicID_2}->{$PositonKEY_2}->{$HitIDkey_2}->{'9_0_7_1_new_1_end'}=$expandWithOutZF_hash->{'7_1_new_1_end'};
  		    $new_Hash->{$fnaPATH_2}->{$CodonTYPE_2}->{$realPEPresidue_2}->{$GenomicID_2}->{$PositonKEY_2}->{$HitIDkey_2}->{'9_0_8_0_ups_0_stt'}=$expandWithOutZF_hash->{'8_0_ups_0_stt'};
  		    $new_Hash->{$fnaPATH_2}->{$CodonTYPE_2}->{$realPEPresidue_2}->{$GenomicID_2}->{$PositonKEY_2}->{$HitIDkey_2}->{'9_0_8_0_ups_1_end'}=$expandWithOutZF_hash->{'8_0_ups_1_end'};
  		    $new_Hash->{$fnaPATH_2}->{$CodonTYPE_2}->{$realPEPresidue_2}->{$GenomicID_2}->{$PositonKEY_2}->{$HitIDkey_2}->{'9_0_8_0_ups_2_DNA'}=$expandWithOutZF_hash->{'8_0_ups_2_DNA'};
  		    $new_Hash->{$fnaPATH_2}->{$CodonTYPE_2}->{$realPEPresidue_2}->{$GenomicID_2}->{$PositonKEY_2}->{$HitIDkey_2}->{'9_0_8_0_ups_3_PEP'}=$expandWithOutZF_hash->{'8_0_ups_3_PEP'};
  		    $new_Hash->{$fnaPATH_2}->{$CodonTYPE_2}->{$realPEPresidue_2}->{$GenomicID_2}->{$PositonKEY_2}->{$HitIDkey_2}->{'9_0_8_1_dws_0_stt'}=$expandWithOutZF_hash->{'8_1_dws_0_stt'};
  		    $new_Hash->{$fnaPATH_2}->{$CodonTYPE_2}->{$realPEPresidue_2}->{$GenomicID_2}->{$PositonKEY_2}->{$HitIDkey_2}->{'9_0_8_1_dws_1_end'}=$expandWithOutZF_hash->{'8_1_dws_1_end'};
  		    $new_Hash->{$fnaPATH_2}->{$CodonTYPE_2}->{$realPEPresidue_2}->{$GenomicID_2}->{$PositonKEY_2}->{$HitIDkey_2}->{'9_0_8_1_dws_2_DNA'}=$expandWithOutZF_hash->{'8_1_dws_2_DNA'};
  		    $new_Hash->{$fnaPATH_2}->{$CodonTYPE_2}->{$realPEPresidue_2}->{$GenomicID_2}->{$PositonKEY_2}->{$HitIDkey_2}->{'9_0_8_1_dws_3_PEP'}=$expandWithOutZF_hash->{'8_1_dws_3_PEP'};
  		    
  		    
  		    
  		    
  		  }
  	  }
    }
  }
  
  return $new_Hash;
}

sub Doing_Statistics_of_Diamond_chunkRst{  #my $out_statistic_hash=DiamondWork_Find_CHAR::Doing_Statistics_of_Diamond_chunkRst( $In_Hash_suit_forStat¡¾,  $cutOff_hitNb, $cutOff_Evalue ¡¿ );
	my ( $In_Hash_suit_forStat, $cutOff_hitNb, $cutOff_Evalue, $cutOff_LcScore )=@_;
	
	my $warnMsgBody="\nIn package  DiamondWork_Find_CHAR,\tIn sub Doing_Statistics_of_Diamond_chunkRst,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	if  (   (  defined ( $In_Hash_suit_forStat )  ) && (  ref ( $In_Hash_suit_forStat ) eq 'HASH'  )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$In_Hash_suit_forStat=$In_Hash_suit_forStat should be a HASH  !!  $!\n\n\n".$caller_inform ); 	}
  
  my $cutOff_hitNb_or_not=0;  
  if  (   (  defined ( $cutOff_hitNb )  ) && (  $cutOff_hitNb =~ m/\S+/  )  ){
    if  (   (  $cutOff_hitNb =~ m/\d+/  ) && ( $cutOff_hitNb > 0 )  ){     	$cutOff_hitNb_or_not=1;    }
    else{		DieWork::Just_dieWork( $die_MsgHead."\n \$cutOff_hitNb=$cutOff_hitNb should be a positive integer > 0  !!  $!\n\n\n".$caller_inform ); 	}    
  }
  
	my $cutOff_Evalue_or_not=0;
  if  (   (  defined ( $cutOff_Evalue )  ) && (  $cutOff_Evalue =~ m/\S+/  )  ){
    if  ( $cutOff_Evalue > 0  ){     	$cutOff_Evalue_or_not=1;    }
    else{		DieWork::Just_dieWork( $die_MsgHead."\n \$cutOff_Evalue=$cutOff_Evalue should be a positive number > 0  !!  $!\n\n\n".$caller_inform ); 	}    
  }
	
	my $cutOff_LcScore_or_not=0;
  if  (   (  defined ( $cutOff_LcScore )  ) && (  $cutOff_LcScore =~ m/\S+/  )  ){
    #if  ( $cutOff_LcScore > 0  ){     	
    	$cutOff_LcScore_or_not=1;    
    #}
    #else{		DieWork::Just_dieWork( $die_MsgHead."\n \$cutOff_LcScore=$cutOff_LcScore should be a  number > 0  !!  $!\n\n\n".$caller_inform ); 	}    
  }
  
  
  
  my $out_statistic_hash;
  my $out_statistic_hash_for_cutoff;
  
  my $inHash=$In_Hash_suit_forStat;
  
  if (  ( ref ($inHash) ) eq 'HASH' ){
    foreach my $keyLev_0 (    sort { $a cmp $b } (   keys (  %{ $inHash } )   )    ){   print "    \$keyLev_0=$keyLev_0\n";
    	my $valLev_0=$inHash->{$keyLev_0};                                                print "    \$valLev_0=$valLev_0\n";
    	my $refLev_0=ref ($valLev_0);                                                     print "    \$refLev_0=$refLev_0\n\n";  
    	
      if ( $refLev_0 eq 'HASH' ){
        foreach my $keyLev_1 (    sort { $a cmp $b } (   keys (  %{ $valLev_0 } )   )    ){   print "      \$keyLev_1=$keyLev_1\n";
        	my $valLev_1=$valLev_0->{$keyLev_1};                                                print "      \$valLev_1=$valLev_1\n";                
        	my $refLev_1=ref ($valLev_1);                                                       print "      \$refLev_1=$refLev_1\n\n";
        	
      	  if ( $refLev_1 eq 'HASH' ){
            foreach my $keyLev_2 (    sort { $a cmp $b } (   keys (  %{ $valLev_1 } )   )    ){   print "        \$keyLev_2=$keyLev_2\n";
            	my $valLev_2=$valLev_1->{$keyLev_2};                                                print "        \$valLev_2=$valLev_2\n";
            	my $refLev_2=ref ($valLev_2);                                                       print "        \$refLev_2=$refLev_2\n\n";
              
               
              
            	if ( $refLev_2 eq 'HASH' ){
                foreach my $keyLev_3 (    sort { $a cmp $b } (   keys (  %{ $valLev_2 } )   )    ){   print "          \$keyLev_3=$keyLev_3\n";
                	my $valLev_3=$valLev_2->{$keyLev_3};                                                print "          \$valLev_3=$valLev_3\n";
                	my $refLev_3=ref ($valLev_3);                                                       print "          \$refLev_3=$refLev_3\n\n";
       
            	    if ( $refLev_3 eq 'HASH' ){
                    foreach my $keyLev_4 (    sort { $a cmp $b } (   keys (  %{ $valLev_3 } )   )    ){   print "            \$keyLev_4=$keyLev_4\n";
                    	
                    	
                    	
                    	my $satisfy_hit_number=0;
                    	
                    	my $valLev_4=$valLev_3->{$keyLev_4};                                                print "            \$valLev_4=$valLev_4\n";
                    	my $refLev_4=ref ($valLev_4);                                                       print "            \$refLev_4=$refLev_4\n\n";
       
            	        if ( $refLev_4 eq 'HASH' ){
                        foreach my $keyLev_5 (    sort { $a cmp $b } (   keys (  %{ $valLev_4 } )   )    ){   print "               \$keyLev_5=$keyLev_5\n";
                        	my $valLev_5=$valLev_4->{$keyLev_5};                                                print "               \$valLev_5=$valLev_5\n";
                        	my $refLev_5=ref ($valLev_5);                                                       print "               \$refLev_5=$refLev_5\n\n";
                                                                              
                          if ( $refLev_5 eq 'HASH' ){
                            my $Evalue=$valLev_5->{'3_0_6_hitEvalue'};               DieWork::Print_and_warn( "\$Evalue=\$valLev_5->{'3_0_6_hitEvalue'}=$Evalue" );
                            my $LcScor=$valLev_5->{'6_2_0_locAlnScore'};             DieWork::Print_and_warn( "\$LcScor=\$valLev_5->{'6_2_0_locAlnScore'}=$LcScor" );
                            
                            DieWork::Print_and_warn( "\$cutOff_Evalue_or_not=$cutOff_Evalue_or_not" );  DieWork::Print_and_warn( "\$cutOff_Evalue=$cutOff_Evalue" );
                            DieWork::Print_and_warn( "\$cutOff_LcScore_or_not=$cutOff_LcScore_or_not" );DieWork::Print_and_warn( "\$cutOff_LcScore=$cutOff_LcScore" );
                            DieWork::Print_and_warn( "START\$satisfy_hit_number=$satisfy_hit_number" );
                            
                            
                            
                            if    (  ( $cutOff_Evalue_or_not == 1 ) && ( $cutOff_LcScore_or_not == 0 ) && ($Evalue <= $cutOff_Evalue )  ){ $satisfy_hit_number++; }
                            elsif (  ( $cutOff_Evalue_or_not == 0 ) && ( $cutOff_LcScore_or_not == 1 ) && ($LcScor >= $cutOff_LcScore)  ){ $satisfy_hit_number++; }
                            elsif (  ( $cutOff_Evalue_or_not == 1 ) && ( $cutOff_LcScore_or_not == 1 ) && ($LcScor >= $cutOff_LcScore ) && ($Evalue <= $cutOff_Evalue )  ){ $satisfy_hit_number++; }
                            elsif (  ( $cutOff_Evalue_or_not == 0 ) && ( $cutOff_LcScore_or_not == 0 )  ) { $satisfy_hit_number++; }
                            
                            DieWork::Print_and_warn( "END  \$satisfy_hit_number=$satisfy_hit_number" );
                            
                            #my $cutOff_hitNb_or_not=0;  #cutOff_Evalue_or_not
                          
                            #foreach my $keyLev_6 (    sort { $a cmp $b } (   keys (  %{ $valLev_5 } )   )    ){   print "                \$keyLev_6=$keyLev_6\n";
                        	  #  my $valLev_6=$valLev_5->{$keyLev_6};                                                print "                \$valLev_6=$valLev_6\n";
                        	    
                                                 	    
                        	                               
                            #}
                          } 
                           			      
                        }
                      }  
                      
                      
                      my $when_to_Count=0; #set a number=0, if the satisfied hit number >= cutoff hit number, then the $when_to_Count=1.
                      if ( $cutOff_hitNb_or_not == 1 ){
                        if ( $satisfy_hit_number >= $cutOff_hitNb ){                        	$when_to_Count=1;                        }
                        else                                       {                        	$when_to_Count=0;                        }
                      }
                      else{
                      	if ( $satisfy_hit_number > 0              ){                        	$when_to_Count=1;                        }
                        else                                       {                        	$when_to_Count=0;                        }                      	
                      }
                      DieWork::Print_and_warn( " \$cutOff_hitNb_or_not=$cutOff_hitNb_or_not" );
                      DieWork::Print_and_warn( " \$cutOff_hitNb=$cutOff_hitNb" );
                      
                      DieWork::Print_and_warn( " \$when_to_Count=$when_to_Count" );
                      
                      
                      #when   $when_to_Count == 1, means the Evalue and Local score cutoff was satisfied, and staisfied hits number >= the hit number cutoff, then count the  $out_statistic_hash->{'organisms_genomeFIle'}->{'TGA/TAA/TAG'->{'X'}
                    	if (  $when_to_Count == 1  ){
                    		if   (      (  defined ( $out_statistic_hash->{$keyLev_0}->{$keyLev_1}->{$keyLev_2} )  ) 
                    	           && (  $out_statistic_hash->{$keyLev_0}->{$keyLev_1}->{$keyLev_2}=~m/\d+/      ) 
                    	           && (  $out_statistic_hash->{$keyLev_0}->{$keyLev_1}->{$keyLev_2} > 0          )   
                    	       )
                    	  { 
                    	  	$out_statistic_hash->{$keyLev_0}->{$keyLev_1}->{$keyLev_2}++;            DieWork::Print_and_warn( "++ : \$out_statistic_hash->{$keyLev_0}->{$keyLev_1}->{$keyLev_2}=$out_statistic_hash->{$keyLev_0}->{$keyLev_1}->{$keyLev_2}" );
                    	  }
                    	  else{
                    	  	$out_statistic_hash->{$keyLev_0}->{$keyLev_1}->{$keyLev_2}=1;            DieWork::Print_and_warn( "=1 : \$out_statistic_hash->{$keyLev_0}->{$keyLev_1}->{$keyLev_2}=$out_statistic_hash->{$keyLev_0}->{$keyLev_1}->{$keyLev_2}" );
                    	  }		
                    	}
                    	         
                    	
                      			          	    
            	      }
            	    }
            	  }
            	}
          	}
          }
       	}
      }
    }
  }
	    	      	
  
  return $out_statistic_hash;
}




1;

##########################################################################################################################################
# 
