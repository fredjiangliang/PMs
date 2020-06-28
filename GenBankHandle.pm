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
use Bio::Index::GenBank;
use Bio::Index::Abstract;
use Storable;
use Storable qw(nstore store_fd nstore_fd freeze thaw dclone);
use ArrayHashChange;

use DirFileHandle;


                      
package GenBankHandle;

my $JL_ProteinGB_DB                                    ="/home/fredjiang/EightT/fredjiang.2018.04.02/GB_database20181213/Pep20181213";
my $JL_ProteinGB_DB_INDEXFILE                          ="/home/fredjiang/EightT/fredjiang.2018.04.02/GB_database20181213/Pep20181213.idx.txt";
my $JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASHFile ="/home/fredjiang/EightT/fredjiang.2018.04.02/GB_database20181213/Pep20181213.Bls.acc_to_Lcs.Hsh";

my $JL_ProteinGB_DB_upgrade_database_and_acc2lcs_Flag  ="/home/fredjiang/EightT/fredjiang.2018.04.02/GB_database20181213/EfatchproteinGB_WorkingFlag.20190326.txt";
my $Longest_waiting_time                               ='20f';
                                         
                                         
#This mathod is used to download genbank information for proteins. 
#The proteins accession IDs are holded by a hash ref. The Acc id are the keys of this hash.
#the gi number ID are the values of this hash. ( get by  $inString=~m/^(gi\|\d+)\|.*/ ) If the "$hit->name() eg:  gi|676388881|ref|XP_009037719.1| " did not fit the regular expression, this gi number id could be wrong!!!!!
#this hash is usually obtained from blast result prase hash.

#then all the accession ID of those proteins will checked in our local GenBank database, if it was downloaded before, then no web woking will trying to do
#If some proteins was not found in our local database, then use the Eutility efach method to do the downloading work from ncbi
#after donwloading a temp file holding the newly downloaded Genbank information will be wrighten in the file we set
#then the local database will upgrade with this file containing all the new proteins
#at last all the GenBank information of those protein will get from our new local database!!
sub Upgrade_db_and_geting_GB_inform_from_ProtHASH {  #GenBankHandle::Upgrade_db_and_geting_GB_inform_from_ProtHASH($wholeProtHASH, $TempQrGBFile, $EachQrGBFile);
	
	my ($wholeProtHASH, $TempQrGBFile, $EachQrGBFile)=@_;  #In this HASH, the key is the '3_0_3_accessionNB' => 'XP_009037719', the value is the '3_0_2_giNub' => 'gi|676388881',
	                         # In sometimes, using the accession id cannot fatch the GenBank information correctly, then the Gi ID will do it rightly! The gi number is get by  $inString=~m/^(gi\|\d+)\|.*/
	                         
	my $warnMsgBody="\nIn package  GenBankHandle,\tIn sub Upgrade_db_and_geting_GB_inform_from_ProtHASH,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $subCallereIfm=DirFileHandle::print_SubCallerInform;
	
	if  (   (  defined ( $wholeProtHASH )  ) && (  ref ( $wholeProtHASH ) eq 'HASH'  )   ){}
  else{		DieWork::Just_dieWork( $die_MsgHead."\n \$wholeProtHASH=$wholeProtHASH should be a HASH  !!  $!\n\n\n".$subCallereIfm ); 	}
	
	if   (   (  defined ( $TempQrGBFile )  ) && ( $TempQrGBFile=~m/\S+/ )   )   {}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$TempQrGBFile=$TempQrGBFile should be a defined not empty string which could be a path to hold genbank information (just new ones)  !!  $!\n\n\n".$subCallereIfm ); 	}
	
	if   (   (  defined ( $EachQrGBFile )  ) && ( $EachQrGBFile=~m/\S+/ )   )   {}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$EachQrGBFile=$EachQrGBFile should be a defined not empty string which could be a path to hold genbank information (for all)  !!  $!\n\n\n".$subCallereIfm ); 	}
	
	
	my $wholeProtArray=ArrayHashChange::Change_Hash_to_Array ($wholeProtHASH) ;#把hash变成数组，key是数组中的元素，
	
	
	#1st, check the $JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASHFile hash, to find those Lucus ID need to be changed
	my $correctHASH;
  my $TParH=GenBankHandle::CorectTheLocusName_for_BlastHitAccName ($wholeProtArray);
  if (   (  defined ( $TParH )  ) && (  ref ( $TParH ) eq 'ARRAY'  )   ) {  ( $wholeProtArray, $correctHASH )=@{ $TParH };  }   # $wholeProtArray changed! and a $correctHASH built!! 
  else { DieWork::Just_dieWork( $die_MsgHead."\n  \$TParH=$TParH something wrong with  GenBankHandle::CorectTheLocusName_for_BlastHitAccName  $!\n\n\n".$subCallereIfm ); }
  
  #2nd, using the $correctHASH, to change the " key '3_0_3_accessionNB' => 'XP_009037719' ID "  to right ID.
  my $Acc_to_Gid_HASH=$wholeProtHASH;  #$Acc_to_Gid_HASH->{  '3_0_3_accessionNB' => 'XP_009037719'   }=  '3_0_2_giNub' => 'gi|676388881' ;
  if (   (  defined ( $correctHASH )  ) && (  ref ( $correctHASH ) eq 'HASH'  )   ) {
  	foreach my $orgAcc (    sort { $a cmp $b } (   keys (  %{ $correctHASH } )   )    ){
  	  my $newAcc=$correctHASH->{$orgAcc};
  	  $Acc_to_Gid_HASH->{$newAcc}=$Acc_to_Gid_HASH->{$orgAcc};
  	  delete ( $Acc_to_Gid_HASH->{$orgAcc} );
  	}  	
  }
  
  #3rd, building the reverse hash $Gid_to_Acc_HASH, and fill the array $EfcProteinArray for intelnet search genbank information
  my $Gid_to_Acc_HASH;                #$Gid_to_Acc_HASH->{  '3_0_2_giNub' => 'gi|676388881'   }=   '3_0_3_accessionNB' => 'XP_009037719';
  my $EfcProteinArray;                #those protein didnot found in local GB databases will put into this array
  if (   (  defined ( $wholeProtArray )  ) && (  ref ( $wholeProtArray ) eq 'ARRAY'  )   ) {
    foreach my $ProtAcc (  @{ $wholeProtArray }  ){ 
    	my $ProtGid=$Acc_to_Gid_HASH->{$ProtAcc};
      $Gid_to_Acc_HASH->{ $ProtGid }=$ProtAcc;
      my $found_in_local_db=0; 
      $found_in_local_db=GenBankHandle::Check_present_in_JL_ProteinGB_DB ($ProtAcc);    #my $msg2= "20181214-0-2\$found_in_local_db=$found_in_local_db\n"; print $msg2; warn $msg2;
      if ( $found_in_local_db == 0 ){    
      	push @{ $EfcProteinArray}, $ProtGid;   
      }  
    }
  }
  
  
  	
  #4th, acturely doing the efach process
  if (   (  defined ( $EfcProteinArray )  ) && (  ref ( $EfcProteinArray ) eq 'ARRAY'  )   ) {
  	my $pt_to_fatcch_NB=@{ $EfcProteinArray };
    my $time_3=`date`; chomp $time_3;
    my $msg4= "\n Now $time_3 , do the EUtilitieswork :: efatch_proteinGB( \$TempQrGBFile=$TempQrGBFile )\n\$pt_to_fatcch_NB=$pt_to_fatcch_NB protein to fatch!!\n\n";  
    DieWork::Print_and_warn( "\n$msg4\n\n"); 
    
    EUtilitieswork::efatch_proteinGB( $EfcProteinArray, $TempQrGBFile );  
  
    #5th, Upgrade the JL PorteinGB_database, input all new genbank files from $TempQrGBFile
    GenBankHandle::Build_upgrade_JL_ProteinGB_DB ($TempQrGBFile);
   
    #6th, Upgrade the Acc_to_realGB_lucus data hashfile     			    
    GenBankHandle::Build_upgrade_JL_DimondAcc_to_realGB_lucus_HASHFile ($TempQrGBFile, $EfcProteinArray, $Gid_to_Acc_HASH);
  
    #7th, re correct the locus name 
    $TParH=GenBankHandle::CorectTheLocusName_for_BlastHitAccName ($wholeProtArray);
    if (   (  defined ( $TParH )  ) && (  ref ( $TParH ) eq 'ARRAY'  )   ) {  ( $wholeProtArray, $correctHASH )=@{ $TParH };  }
  }
  
  
  
  #8th, build  the Genbank information holding file with all protines
  my $outHASH=GenBankHandle::Building_GB_inform_HASH_for_many_locusIDs ( $JL_ProteinGB_DB_INDEXFILE, $wholeProtArray );
  
  
  #my $outHASH=GenBankHandle::Extract_GB_inform_base_from_TOTALidxOBJ($JL_ProteinGB_DB_INDEXFILE, $inLocusARRAY);
	#my $GbkEfachWkOutARRAY=GenBankHandle::BuildTreesHash_from_geneBankFiles($GenBankFDir_0_2_a, $opt_c);
      
	return $outHASH;
}




sub CorectTheLocusName_for_BlastHitAccName{  #  GenBankHandle::CorectTheLocusName_for_BlastHitAccName ($org_nameArray);
	my ($org_nameArray)=@_;
	
	my $warnMsgBody="\nIn package  GenBankHandle,\tIn sub CorectTheLocusName_for_BlastHitAccName,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $subCallereIfm=DirFileHandle::print_SubCallerInform;
	
	my $new_nameArray;  my $correctHASH; my $changed_or_not=0;
	if  (   (  defined ( $org_nameArray )  ) && (  ref ( $org_nameArray ) eq 'ARRAY'  )   ) {
		if (   (  defined ( $JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASHFile )  ) && ( $JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASHFile=~m/\S+/ ) && (  -e ( $JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASHFile )  )   ) {
	    my $JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASH=Storable::retrieve ($JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASHFile);
	    if  (   (  defined ( $JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASH )  ) && (  ref ( $JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASH ) eq 'HASH'  )   ) {
		    for (  my $i=0; $i<@{ $org_nameArray }; $i++  ){
			    my $eachNm=$org_nameArray->[$i]; 
			    if (   ( $eachNm=~m/\S+/ ) && (  defined ( $JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASH->{$eachNm} )  ) && ( $JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASH->{$eachNm}=~m/\S+/ )   ){
			    	$new_nameArray->[$i]     =$JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASH->{$eachNm};
			    	$correctHASH->{ $eachNm }=$JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASH->{$eachNm};
			    	$changed_or_not=1;
			    }
			    else {
			    	$new_nameArray->[$i]=$eachNm; 
			    }
			  }
			}			
		}		
	}
	if ( $changed_or_not==0 ){
		return [$org_nameArray, $correctHASH ] ;
	}
	else {
		return [$new_nameArray, $correctHASH ];
	}
	
}


sub Check_present_in_JL_ProteinGB_DB{  #  GenBankHandle::Check_present_in_JL_ProteinGB_DB ($LocusNAME);
	my ($LocusNAME)=@_;
	
	my $warnMsgBody="\nIn package  GenBankHandle,\tIn sub Check_present_in_JL_ProteinGB_DB,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $subCallereIfm=DirFileHandle::print_SubCallerInform;                           #warn  $warnMsgBody ; 
	
	my $outPut_True_or_False=0;
	$outPut_True_or_False=GenBankHandle::CHeck_present_of_a_locus ($JL_ProteinGB_DB_INDEXFILE, $LocusNAME);
	
	#if ($outPut_True_or_False==0){
	#  if (   (  defined ( $JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASHFile )  ) && ( $JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASHFile=~m/\S+/ ) && (  -e ( $JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASHFile )  )   ) {
	#    my $JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASH=Storable::retrieve ($JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASHFile);
	#    if  (   (  defined ( $JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASH )  ) && (  ref ( $JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASH ) eq 'HASH'  )   ) {
	#    	if (   (  defined ( $JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASHFile->{$LocusNAME} )  ) && ( $JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASHFile->{$LocusNAME}=~m/\S+/ )   ){
	#    		my $realLcsName=$JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASHFile->{$LocusNAME};
	#    		$outPut_True_or_False=GenBankHandle::CHeck_present_of_a_locus ($JL_ProteinGB_DB_INDEXFILE, $realLcsName);
	#    	}
	#    }
	#  }
	#}	
	return $outPut_True_or_False;
}

sub Upgrade_JL_proteinGB_and_Acc_to_Lcs_HASHs_onlyOnePid{  #  GenBankHandle::Upgrade_JL_proteinGB_and_Acc_to_Lcs_HASHs_onlyOnePid ($new_GB_file, $proteinGiArray, $gi_to_old_lucs_Hash);
	# 根据 新下载的$new_GB_file,以及 蛋白名的 array $proteinGiArray， 和 gi到lucus的映射hash，对 GenBank格式的ProteinDB数据库和Acc2Lcs的数据hash，进行更新，此更新在这台电脑上，只能进行 排它性地 单独运行，不能并行运算
	
	my ($new_GB_file, $proteinGiArray, $gi_to_old_lucs_Hash)=@_;
	
	my $warnMsgBody="\nIn package  GenBankHandle,\tIn sub Upgrade_JL_proteinGB_and_Acc_to_Lcs_HASHs_onlyOnePid,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $subCalInfom=DirFileHandle::print_SubCallerInform;
  
  
  my $nowWorkTimePidInformMark=OnlyOnePidWork::CheckAndStart_step_forOnlyOnePidWork($JL_ProteinGB_DB_upgrade_database_and_acc2lcs_Flag, $Longest_waiting_time);	
  
  if (   ( defined ( $nowWorkTimePidInformMark )  ) && ( $nowWorkTimePidInformMark=~/\S+/ )   ){
  	
  	GenBankHandle::Build_upgrade_JL_ProteinGB_DB ($new_GB_file);
  	
  	if (   (  defined ( $gi_to_old_lucs_Hash )  ) && (  ref ( $gi_to_old_lucs_Hash ) eq 'HASH'  )   ){
  	  GenBankHandle::Build_upgrade_JL_DimondAcc_to_realGB_lucus_HASHFile ($new_GB_file, $proteinGiArray, $gi_to_old_lucs_Hash);
  	}
  	
  	
  	OnlyOnePidWork::CheckAnd_Done_step_forOnlyOnePidWork( $JL_ProteinGB_DB_upgrade_database_and_acc2lcs_Flag, $nowWorkTimePidInformMark );  	
  
  }
  
}



sub Build_upgrade_JL_ProteinGB_DB{ #  GenBankHandle::Build_upgrade_JL_ProteinGB_DB ($new_GB_file);
	my ($new_GB_file)=@_;
	
	my $warnMsgBody="\nIn package  GenBankHandle,\tIn sub Build_upgrade_JL_ProteinGB_DB,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $subCallereIfm=DirFileHandle::print_SubCallerInform;                    #warn  $warnMsgBody  ;
  
  my $upGrate_time=GenBankHandle::Build_and_Update_GeneBankFile_from_new_GenBankFIle ($JL_ProteinGB_DB, $new_GB_file);
  return $upGrate_time;
}

#   Build_GB_file_from_IDarray_GBidxFILE
sub Build_GB_file_from_IDarray_JL_GBidxFILE{   #GenBankHandle::Build_GB_file_from_IDarray_JL_GBidxFILE ( $ID_ARRAY, $outFIle );
	my ($ID_ARRAY, $outFIle)=@_;
	
	my $warnMsgBody="\nIn package  GenBankHandle,\tIn sub Build_GB_file_from_IDarray_JL_GBidxFILE,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $subCallereIfm=DirFileHandle::print_SubCallerInform;                  #warn  $warnMsgBody ; 
	
	GenBankHandle::Build_GB_file_from_IDarray_GBidxFILE ( $ID_ARRAY, $JL_ProteinGB_DB_INDEXFILE, $outFIle );
  
}


sub Build_upgrade_JL_DimondAcc_to_realGB_lucus_HASHFile{ #  GenBankHandle::Build_upgrade_JL_DimondAcc_to_realGB_lucus_HASHFile ($new_GB_file, $proteinGiArray, $gi_to_old_lucs_Hash);
	my ($new_GB_file, $proteinGiArray, $gi_to_old_lucs_Hash)=@_;
	
	my $warnMsgBody="\nIn package  GenBankHandle,\tIn sub Build_upgrade_JL_DimondAcc_to_realGB_lucus_HASHFile,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $subCallereIfm=DirFileHandle::print_SubCallerInform;                 #warn  $warnMsgBody ; 
  
  my $GB_file_HASH=GenBankHandle::ReadGeneBankFILE( $new_GB_file, $proteinGiArray );
  if  (   (  defined ( $GB_file_HASH )  ) && (  ref ( $GB_file_HASH ) eq 'HASH'  )   ){
  	my $ne_found=0;
  	my $ChangeHash;
  	foreach my $GB_lcs (    sort { $a cmp $b } (   keys (  %{ $GB_file_HASH } )   )    ){ 
  		my $giNb=$GB_file_HASH->{$GB_lcs}->{'002_org_inID'};
  		my $org_lcs_name=$gi_to_old_lucs_Hash->{ $giNb };
  		if ( $GB_lcs ne $org_lcs_name){
  			$ne_found=1; 
  			$ChangeHash->{$org_lcs_name}=$GB_lcs;
  		}  		
  	}  
  	if ($ne_found==1 ){
  		if (   (  defined ( $JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASHFile )  ) && ( $JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASHFile=~m/\S+/ ) && (  -e ( $JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASHFile )  )   ) {
	                                                           
	      my $JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASH=Storable::retrieve ($JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASHFile);
	      if  (   (  defined ( $JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASH )  ) && (  ref ( $JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASH ) eq 'HASH'  )   ) {
	      	foreach my $eachLcs (   keys (  %{ $ChangeHash }  )   ){
            $JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASH->{$eachLcs}=$ChangeHash->{$eachLcs};
	      	}
	      	DirFileHandle::PrintDumper ( $JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASHFile,  $JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASH   ) if (   (  defined ( $JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASH  )  ) && (  ref ( $JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASH  ) eq 'HASH'  )    );  
	      }
	      
	    }   
	    else{
	    	my $JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASH ;
	    	foreach my $eachLcs (   keys (  %{ $ChangeHash }  )   ){
          $JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASH->{$eachLcs}=$ChangeHash->{$eachLcs};
	      }
	      DirFileHandle::PrintDumper ( $JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASHFile,  $JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASH   ) if (   (  defined ( $JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASH  )  ) && (  ref ( $JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASH  ) eq 'HASH'  )    );  
	      
	    }
		  
  	}
  }   
  
  
}


########################################################################################


sub Build_and_Update_GeneBankFile_from_new_GenBankFIle{ #  GenBankHandle::Build_and_Update_GeneBankFile_from_new_GenBankFIle ($org_genBankFIle, $new_genBankFile);
	my ($org_genBankFIle, $new_genBankFile)=@_;
	my $warnMsgBody="\nIn package  GenBankHandle,\tIn sub Build_and_Update_GeneBankFile_from_new_GenBankFIle,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $subCallereIfm=DirFileHandle::print_SubCallerInform;
  
  my $fisrtWright=0;
  if ( -e ( $org_genBankFIle ) ){
  	
  }
  else{
  	$fisrtWright=1;
  	system ("cp -f $new_genBankFile $org_genBankFIle") 
  }
  
  my $org_genBank_idx=$org_genBankFIle.".idx.txt";
  
  my $inx_obj;
  #if ( -e ( $org_genBank_idx ) ){
  #	
  #	my $absIdx=Bio::Index::Abstract->new(-filename    => $org_genBank_idx,
  #	   	   	   	   	   	   	   	   	   	 -verbose     => 1);
  #	warn $absIdx->_check_file_sizes ."\n\n"; sleep (2);
  #                    
  #	 
  #	#$inx_obj=Bio::Index::GenBank->new( -filename   => $org_genBank_idx  );
  #}
  #else{
    $inx_obj=GenBankHandle::BuildIdxFile_for_GenBankFIle ($org_genBankFIle, $org_genBank_idx);
    #print "\n\nFrom \$org_genBankFIle=$org_genBankFIle build a index file:\$org_genBank_idx=$org_genBank_idx!!!!\n\n";
  #}
  my $org_count_records  =$inx_obj->count_records();
  
  
  my $in_FileOBJ  = Bio::SeqIO->new( -file   => $new_genBankFile ,
                                     -format => 'genbank'   );
                                    
  my $outFileOBJ  = Bio::SeqIO->new( -file   => ">>$org_genBankFIle" ,
                                     -format => 'genbank'    ,                                  
                                     #-flush  => 0               
                                   ); # go as fast as we can!
  my $upgrate_or_not=0;  my $upGrate_time=0;                               
  while( my $seq = $in_FileOBJ->next_seq ) { 
  	
  	my $sequnece_ID= $seq->id(); 
  	my $outPut=$inx_obj->get_Seq_by_id($sequnece_ID);
  	#print ( ref ( $outPut ) ), "xxxx\n";
  	if (   (  defined ( $outPut ) ) && (  $outPut=~/\S+/ )   ){
  		#print "The \$sequnece_ID=$sequnece_ID found in \$new_genBankFile=$new_genBankFile, is found in  \$org_genBankFIle=$org_genBankFIle \$org_genBank_idx=$org_genBank_idx\n";
  	}
  	
  	else {
  		$upGrate_time++;
  		if ($upGrate_time==1){
  			my $old_orgDatabaseFile=$org_genBankFIle.".old.txt";
  			system ("cp -f $org_genBankFIle $old_orgDatabaseFile")
  		}
  		$outFileOBJ->write_seq($seq);  #print "\n   NEW!!! The \$sequnece_ID=$sequnece_ID Cannot found in  \$org_genBankFIle=$org_genBankFIle \$org_genBank_idx=$org_genBank_idx\n\n";
  		$upgrate_or_not=1; 
  	}
  	#
  }
  
  if (  $upgrate_or_not==1 ){
  	my $new_idx=GenBankHandle::BuildIdxFile_for_GenBankFIle ($org_genBankFIle, $org_genBank_idx);
  	my $new_count_records  =$new_idx->count_records();
  	#print "\n\n\ $upGrate_time seq was newly added into $org_genBankFIle=$org_genBankFIle \$org_genBank_idx=$org_genBank_idx Upgrated!!!!\n There are $new_count_records now!!!\n\n";
  }
  else {
  	#if ($fisrtWright==1){ print "\n\n\ $org_count_records seq was newly added into $org_genBankFIle=$org_genBankFIle \$org_genBank_idx=$org_genBank_idx Upgrated!!!!\n There are $org_count_records now!!!\n\n"; }
  	#else { print "\n\nIn $new_genBankFile=$new_genBankFile , Nothing new, no need to upgrate database!!!\nThere are $org_count_records  still!!!\n"; }
  }
  
  if ($fisrtWright==1){
  	return $org_count_records;
  }
  return $upGrate_time;
  
}

sub BuildGenbankDatabase_from_moreGenBankFIle{ #  GenBankHandle::BuildGenbankDatabase_from_moreGenBankFIle ($in_FIle, $outFIle);
	my ($in_FIle, $outFIle)=@_;
	my $warnMsgBody="\nIn package  GenBankHandle,\tIn sub BuildGenbankDatabase_from_moreGenBankFIle,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $subCallereIfm=DirFileHandle::print_SubCallerInform;
  my $in_FileOBJ  = Bio::SeqIO->new( -file   => $in_FIle ,
                                     -format => 'genbank'   );
                                    
  my $outFileOBJ  = Bio::SeqIO->new( -file   => ">>$outFIle" ,
                                     -format => 'genbank'    ,                                  
                                     -flush  => 0               ); # go as fast as we can!
  while( my $seq = $in_FileOBJ->next_seq) { 
  	$outFileOBJ->write_seq($seq);
  }
  
  
}


sub Build_GB_file_from_IDarray_GBidxFILE{ #  GenBankHandle::Build_GB_file_from_IDarray_GBidxFILE ( $ID_ARRAY, $GBidxFIle, $outFIle );
	my ( $ID_ARRAY, $GBidxFIle, $outFIle )=@_;  #, $in_GBFIle
	my $warnMsgBody="\nIn package  GenBankHandle,\tIn sub Build_GB_file_from_IDarray_GBidxFILE,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $subCallereIfm=DirFileHandle::print_SubCallerInform;
  
  #my $in_FileOBJ  = Bio::SeqIO->new( -file   => $in_FIle ,
   #                                  -format => 'genbank'   );
  my $GBidxFIleO  = Bio::Index::GenBank->new( -filename   => $GBidxFIle );
                                    
  my $outFileOBJ  = Bio::SeqIO->new( -file   => ">$outFIle" ,
                                     -format => 'genbank'    ,                                  
                                     #-flush  => 0                   # go as fast as we can!
                                   );
                                   
  
  if (   (  defined ( $ID_ARRAY )  ) && (  ( ref ($ID_ARRAY) ) eq 'ARRAY'  )   ){
  	foreach my $eachProID (  @{ $ID_ARRAY }  ){
  		my $seq = $GBidxFIleO->fetch($eachProID); # Returns Bio::Seq object
                                     warn "20181214-1 \$eachProID=$eachProID \$seq=$seq\n"; print "20181214-1 \$eachProID=$eachProID \$seq=$seq\n";
  		$outFileOBJ->write_seq($seq);  warn "20181214-2 \$eachProID=$eachProID \$seq=$seq\n"; print "20181214-2 \$eachProID=$eachProID \$seq=$seq\n";
    }
  }
   
  
}



sub BuildIdxFile_for_GenBankFIle{  #  GenBankHandle::BuildIdxFile_for_GenBankFIle ($in_GenBank_FIle, $out_Index_FIle);
	my ($in_GenBank_FIle, $out_Index_FIle)=@_;
	my $warnMsgBody="\nIn package  GenBankHandle,\tIn sub BuildIdxFile_for_GenBankFIle,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $subCallereIfm=DirFileHandle::print_SubCallerInform;
  
  if ( -e ( $out_Index_FIle ) ){
    system ("mv -f $out_Index_FIle $out_Index_FIle.old");
  }
  my $inx = Bio::Index::GenBank->new( -filename   => $out_Index_FIle, 
                                      -write_flag => 'WRITE'    );
  $inx->make_index($in_GenBank_FIle); 
  return $inx;
}





sub CHeck_present_of_a_locus{  #  GenBankHandle::CHeck_present_of_a_locus ($in_GenBank_index_FIle, $LocusNAME);
	my ($in_GenBank_index_FIle, $LocusNAME)=@_;
	my $warnMsgBody="\nIn package  GenBankHandle,\tIn sub CHeck_present_of_a_locus,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $subCallereIfm=DirFileHandle::print_SubCallerInform;
  
  my $inx = Bio::Index::GenBank->new( -filename   => $in_GenBank_index_FIle  );
  my $outPut=$inx->get_Seq_by_id($LocusNAME);
  
  my $found_no_not=0;  #warn "20181214-3 \$LocusNAME=$LocusNAME \$outPut=$outPut\n"; print "20181214-3 \$LocusNAME=$LocusNAME \$outPut=$outPut\n";
  if (   (  defined ( $outPut ) ) && (  $outPut=~/\S+/ )   ){
  	$found_no_not=1;
  	#print "The \$sequnece_ID=$sequnece_ID found in \$new_genBankFile=$new_genBankFile, is found in  \$org_genBankFIle=$org_genBankFIle \$org_genBank_idx=$org_genBank_idx\n";
  }
  #warn "20181214-4 \$found_no_not=$found_no_not\n"; print "20181214-4 \$found_no_not=$found_no_not\n";
  
  return $found_no_not;
}

# Using the LocusNames array to search and build Genbank Hash
# a global value $JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASHFile was used
sub Building_GB_inform_HASH_for_many_locusIDs{    #my $outHASH=GenBankHandle::Building_GB_inform_HASH_for_many_locusIDs($in_GenBank_index_FIle, $inLocusARRAY);
	
	
	my ($in_GenBank_index_FIle, $inLocusARRAY)=@_;
	my $warnMsgBody="\nIn package  GenBankHandle,\tIn sub Building_GB_inform_HASH_for_many_locusIDs,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $subCallereIfm=DirFileHandle::print_SubCallerInform;
  
  if  (   (  defined ( $inLocusARRAY )  ) && (  ref ( $inLocusARRAY ) eq 'ARRAY'  )   ){}
  else{		DieWork::Just_dieWork( $die_MsgHead."\n \$inLocusARRAY=$inLocusARRAY should be a ARRAY  !!  $!\n\n\n".$subCallereIfm ); 	}
	
	if   (   (  defined ( $in_GenBank_index_FIle )  ) && ( $in_GenBank_index_FIle=~m/\S+/ ) && (  ( -e $in_GenBank_index_FIle )  )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$in_GenBank_index_FIle=$in_GenBank_index_FIle should be a defined not empty file  !!  $!\n\n\n".$subCallereIfm ); 	}
	
	DieWork::Check_FileDirExist_or_DIE( $JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASHFile, "\$JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASHFile", $die_MsgHead, $subCallereIfm  );
	my $JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASH=Storable::retrieve( $JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASHFile );
	my $JL_ProteinGB_DB_Dimond_realGB_lucus_to_Acc_HASH=ArrayHashChange::BuildReverseHash($JL_ProteinGB_DB_DimondAcc_to_realGB_lucus_HASH);
	
	my $inx = Bio::Index::GenBank->new( -filename   => $in_GenBank_index_FIle  );
	
	my $OutHASH;
	foreach my $eachLcus_ID (  @{ $inLocusARRAY }  ){
	  my $subHASH=GenBankHandle::Extract_GB_inform_base_from_TOTALidxOBJ($inx, $eachLcus_ID);
	  
	  #DieWork::Check_DfdNoEmptString_or_DIE( $JL_ProteinGB_DB_Dimond_realGB_lucus_to_Acc_HASH->{$eachLcus_ID}, "\$JL_ProteinGB_DB_Dimond_realGB_lucus_to_Acc_HASH->{\$eachLcus_ID}=\$JL_ProteinGB_DB_Dimond_realGB_lucus_to_Acc_HASH->{$eachLcus_ID}", $die_MsgHead, $subCallereIfm  );
	  my $changed_back_Lucus_name=$eachLcus_ID;
	  if (  DieWork::Check_DfdNoEmptString_or_NOT( $JL_ProteinGB_DB_Dimond_realGB_lucus_to_Acc_HASH->{$eachLcus_ID} )  ){
	  	$changed_back_Lucus_name=$JL_ProteinGB_DB_Dimond_realGB_lucus_to_Acc_HASH->{$eachLcus_ID};
	  }
	  $OutHASH->{$changed_back_Lucus_name}=Storable::dclone( $subHASH );
	  #$OutHASH->{$eachLcus_ID}=Storable::dclone( $subHASH );
	}
	
	return $OutHASH;
  
}

# Using the LocusName to search and build Genbank Hash
sub Extract_GB_inform_base_from_TOTALidxOBJ{      #my $outHASH=GenBankHandle::Extract_GB_inform_base_from_TOTALidxOBJ($in_GenBank_index_FIle, $LocusNAME);
	
	my ($in_GenBank_index_OBJ, $LocusNAME)=@_;
	my $warnMsgBody="\nIn package  GenBankHandle,\tIn sub Extract_GB_inform_base_from_TOTALidxOBJ,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $subCallereIfm=DirFileHandle::print_SubCallerInform;
  
  #my $inx = Bio::Index::GenBank->new( -filename   => $in_GenBank_index_FIle  );
  my $seqObj=$in_GenBank_index_OBJ->get_Seq_by_id($LocusNAME);
  
  my $outHASH;
  
  my $sequnece_ID  = $outHASH->{'000_bioPerlSeqObj_ID'}=$seqObj->id(); print "\$sequnece_ID=$sequnece_ID\n";
  	
  	                   #$outHASH->{'000_gbFileNm'}=$inFIle;          
  	                   #$outHASH->{'001_inFileNB'}=$step_nb+1;   
  	                   #$outHASH->{'002_org_inID'}=$inProtIDarray->[$step_nb] if (   (  defined ( $inProtIDarray )  ) && (  ref ( $inProtIDarray )  eq 'ARRAY'  ) && (  defined ( $inProtIDarray->[$step_nb] )  ) && (  $inProtIDarray->[$step_nb]=~m/\S+/ )   );          
  	                            	                   
  my $sequnece     = $outHASH->{'0_0_sequence'}=$seqObj->seq(); print "\$sequnece=$sequnece\n";
  if ($sequnece=~m/U/){
  	                 $outHASH->{'0_1_UfoundIN'}=1;               
  }
  
  my $sequnece_Desc= $outHASH->{'1_0_descript'}=$seqObj->desc(); print "\$sequnece_Desc=$sequnece_Desc\n";
  
  ##DirFileHandle::PrintAndWarnDumper($seqObj);
  #for my $feat_object ($seqObj->get_SeqFeatures) {
   #print "primary tag: ", $feat_object->primary_tag, "\n\n";
   #
   #for my $tag ($feat_object->get_all_tags) {
   #    print "  tag: ", $tag, "\n";
   #    for my $value ($feat_object->get_tag_values($tag)) {
   #        print "    value: ", $value, "\n";
   #    }
   #}
   #print "\n\n";
   
   # legible and long
   #my $species_object = $seqObj->species;  
   #my $species_string = $species_object->node_name;   print "11 \$species_string=$species_string\n";

   # Perlish
   my $species_string     = $outHASH->{'2_0__orgaism'}= $seqObj->species->node_name;    print "\$species_string=$species_string\n";
   my $species_txid       = $outHASH->{'3_0_txnomiID'}= $seqObj->species->ncbi_taxid;   print "\$species_txid=$species_txid\n";
   my $species_common_name= $outHASH->{'4_0_specComN'}= $seqObj->species->common_name;  #warn "\$species_common_name=$species_common_name\n";


   # either way, $species_string is "Homo sapiens"

   # get all taxa from the ORGANISM section in an array
   my @classification = $seqObj->species->classification;  print "\@classification="; 
   @classification = reverse ( @classification );
   my $foundSPEhash;
   my @newARRAY;
   FOREACHMK: foreach my $eachPHY ( @classification ){
   	if (   (  defined ( $foundSPEhash->{$eachPHY} )  ) && ( $foundSPEhash->{$eachPHY}=~m/\S+/ )   ){
   		push @newARRAY, $species_string;
   		last FOREACHMK;
   	}
   	else {
   		push @newARRAY, $eachPHY;
   	}
   	$foundSPEhash->{$eachPHY}=1;
   	print "$eachPHY     ";
   }
   print "\n";
   $outHASH->{'5_0_phyClsif'}=[  @newARRAY  ];
   
   return $outHASH;
}

sub BuildFastaFile{
	my ($InHASH, $SpeCountHASH)=@_;
	my $warnMsgBody="\nIn package  GenBankHandle,\tIn sub BuildFastaFile,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $subCallereIfm=DirFileHandle::print_SubCallerInform;
  
  my $representCUToff=0.3;
  my $allHitNb=(   keys (  %{ $InHASH }  )   ) ;
  my $cutOffNb=$allHitNb*$representCUToff;
  
  my $outInfm_HASH;
  
  my $outString; my $stepNB=1;
  if (   (  defined ( $InHASH )  ) && (  ref ( $InHASH ) eq 'HASH' )   ){
  	foreach my $eachPRPid (    sort {  ( $InHASH->{$a}->{'000_gbFileNm'} cmp $InHASH->{$b}->{'000_gbFileNm'} ) || ( $InHASH->{$a}->{'001_inFileNB'} cmp $InHASH->{$b}->{'001_inFileNB'} )   } (   keys (  %{ $InHASH }  )   )    ){
  		
  		#my $theShownPhy='';
  		#if (   (  defined ( $SpeCountHASH )  ) && (  ref ( $SpeCountHASH ) eq 'HASH' )   ){
  		#	
  		#	if  (   (  defined ( $InHASH->{$eachPRPid}->{'5_0_phyClsif'} )  ) && (  ref ( $InHASH->{$eachPRPid}->{'5_0_phyClsif'} ) eq 'ARRAY' )   ){
  		#    foreach my $eachPHY (  @{ $InHASH->{$eachPRPid}->{'5_0_phyClsif'} }  ){
      #      if (   (  defined ( $SpeCountHASH->{$eachPHY} )  ) && ( $SpeCountHASH->{$eachPHY} >0 )   ){
      #      	if ( $SpeCountHASH->{$eachPHY} >= $cutOffNb){
      #      		$theShownPhy=$eachPHY;
      #      	}
      #      }
  		#    }
  		#  }
  		#}
  		                                                              
      my $UfindMk=''; 
      #my $Prot_ID=$theShownPhy."_".$InHASH->{$eachPRPid}->{'2_0__orgaism'}."_".$InHASH->{$eachPRPid}->{'1_0_descript'};           
      my $Prot_ID=$eachPRPid;
      if (   (  defined ( $InHASH->{$eachPRPid}->{'0_1_UfoundIN'} )  ) && (  $InHASH->{$eachPRPid}->{'0_1_UfoundIN'} == 1 )   ){
  			#$thisString="u".$thisString;  
  			$UfindMk='u';
  		}
      my $addifProteinID="B".$UfindMk."_".$Prot_ID;
      my $Prot_Sq       =$InHASH->{$eachPRPid}->{'0_0_sequence'};
  		my $thisString=">".$addifProteinID."\n".$Prot_Sq."\n\n"; #warn "\$thisString=$thisString\n";
  		
  		
  		$outString.=$thisString;
  		
  		
  		###############
  		$outInfm_HASH->{'0_0_0_numberKey_hash'}->{$stepNB}->{'0_1_0_Protein___Type'}='3_the_3rd_round_result';
    	$outInfm_HASH->{'0_0_0_numberKey_hash'}->{$stepNB}->{'0_1_1_only_proteinID'}=$Prot_ID;
    	$outInfm_HASH->{'0_0_0_numberKey_hash'}->{$stepNB}->{'0_1_2_adifmProteinID'}=$addifProteinID;
    	$outInfm_HASH->{'0_0_0_numberKey_hash'}->{$stepNB}->{'0_1_3_found_U_in_Seq'}=$UfindMk;
    	$outInfm_HASH->{'0_0_0_numberKey_hash'}->{$stepNB}->{'0_1_4_lucus_or_accse'}=$eachPRPid;
    	$outInfm_HASH->{'0_0_0_numberKey_hash'}->{$stepNB}->{'0_1_a_only__Sequence'}=$Prot_Sq;
  		##################
  		
  		$stepNB++;
  	}
  }
  return [ $outString, $outInfm_HASH ];
}

sub ReadGeneBankFILE {  #GenBankHandle::ReadGeneBankFILE
	my ($inFIle, $inProtIDarray, $top_howMany)=@_;
	
	my $warnMsgBody="\nIn package  GenBankHandle,\tIn sub ReadGeneBankFILE,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $subCallereIfm=DirFileHandle::print_SubCallerInform;
  
  my $inFileOBJ  = Bio::SeqIO->new( -file   => $inFIle ,
                                    -format => 'genbank');
  my $outHASH; my $step_nb=0;
  WHILEMK: while ( my $seqObj = $inFileOBJ->next_seq() ) {
  	
  	print "\n\n\n\n";
  	my $sequnece_ID= $seqObj->id(); print "\$sequnece_ID=$sequnece_ID\n";
  	
  	                   $outHASH->{$sequnece_ID}->{'000_gbFileNm'}=$inFIle;          
  	                   $outHASH->{$sequnece_ID}->{'001_inFileNB'}=$step_nb+1;   
  	                   $outHASH->{$sequnece_ID}->{'002_org_inID'}=$inProtIDarray->[$step_nb] if (   (  defined ( $inProtIDarray )  ) && (  ref ( $inProtIDarray )  eq 'ARRAY'  ) && (  defined ( $inProtIDarray->[$step_nb] )  ) && (  $inProtIDarray->[$step_nb]=~m/\S+/ )   );          
  	                            	                   
  	my $sequnece     = $outHASH->{$sequnece_ID}->{'0_0_sequence'}=$seqObj->seq(); print "\$sequnece=$sequnece\n";
  	if ($sequnece=~m/U/){
  		                 $outHASH->{$sequnece_ID}->{'0_1_UfoundIN'}=1; print "\$sequnece=$sequnece\n";
  	}
  	
  	my $sequnece_Desc= $outHASH->{$sequnece_ID}->{'1_0_descript'}=$seqObj->desc(); print "\$sequnece_Desc=$sequnece_Desc\n";
  	
  	##DirFileHandle::PrintAndWarnDumper($seqObj);
  	#for my $feat_object ($seqObj->get_SeqFeatures) {
    #print "primary tag: ", $feat_object->primary_tag, "\n\n";
    #
    #for my $tag ($feat_object->get_all_tags) {
    #    print "  tag: ", $tag, "\n";
    #    for my $value ($feat_object->get_tag_values($tag)) {
    #        print "    value: ", $value, "\n";
    #    }
    #}
    #print "\n\n";
    
    # legible and long
    #my $species_object = $seqObj->species;  
    #my $species_string = $species_object->node_name;   print "11 \$species_string=$species_string\n";

    # Perlish
    my $species_string     = $outHASH->{$sequnece_ID}->{'2_0__orgaism'}= $seqObj->species->node_name;  print "\$species_string=$species_string\n";
    my $species_txid       = $outHASH->{$sequnece_ID}->{'3_0_txnomiID'}= $seqObj->species->ncbi_taxid;  print "\$species_txid=$species_txid\n";
    my $species_common_name= $outHASH->{$sequnece_ID}->{'4_0_specComN'}= $seqObj->species->common_name;  #warn "\$species_common_name=$species_common_name\n";


    # either way, $species_string is "Homo sapiens"

    # get all taxa from the ORGANISM section in an array
    my @classification = $seqObj->species->classification;  print "\@classification="; 
    @classification = reverse ( @classification );
    my $foundSPEhash;
    my @newARRAY;
    FOREACHMK: foreach my $eachPHY ( @classification ){
    	if (   (  defined ( $foundSPEhash->{$eachPHY} )  ) && ( $foundSPEhash->{$eachPHY}=~m/\S+/ )   ){
    		push @newARRAY, $species_string;
    		last FOREACHMK;
    	}
    	else {
    		push @newARRAY, $eachPHY;
    	}
    	$foundSPEhash->{$eachPHY}=1;
    	print "$eachPHY     ";
    }
    print "\n";
    $outHASH->{$sequnece_ID}->{'5_0_phyClsif'}=[  @newARRAY  ];
    # "sapiens Homo Hominidae Catarrhini Primates Eutheria Mammalia
    #  Euteleostomi Vertebrata Craniata Chordata Metazoa Eukaryota"
    
    $step_nb++;
    
    if  (   (  defined ( $top_howMany )  ) && ( $top_howMany=~m/\d+/ ) && ( $top_howMany > 0 )   ){
      if ( $step_nb >= $top_howMany ){
      	last WHILEMK;
      }
    }
    
  }
  	
  return $outHASH;
  
  
}





sub BuildTreesHash_from_geneBankHASH{ # my ($new_outTreeWord, $Tree_outHASH )=@{ GenBankHandle::BuildTreesHash_from_geneBankHASH ($AllGeneBkHASH) };
	my ($AllGeneBkHASH, $showing_taxonomy_Ifm_HASH, $PtNbID_toChar_HASH, $ptWd)=@_;
	
	my $warnMsgBody="\nIn package  GenBankHandle,\tIn sub BuildTreesHash_from_geneBankHASH,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $callerMSG=DirFileHandle::print_SubCallerInform;
	
	if  (   (  defined ( $AllGeneBkHASH )  ) && (  ref ( $AllGeneBkHASH ) eq 'HASH' )   ){
		
	 
	  my $speciseHASH=GenBankHandle::BuildAll_StfcName_HASH            ($AllGeneBkHASH);  #DirFileHandle::PrintDumper("0_0_2_speciseHASH.hash",$speciseHASH) if (   (  defined ( $speciseHASH )  ) && (  ref ( $speciseHASH ) eq 'HASH' )   );
	  my $AllPHY_HASH=GenBankHandle::BuildAllTxnomy_phy_HASH           ($AllGeneBkHASH);  #DirFileHandle::PrintDumper("0_0_3_AllPHY_HASH.hash",$AllPHY_HASH) if (   (  defined ( $AllPHY_HASH )  ) && (  ref ( $AllPHY_HASH ) eq 'HASH' )   );
	  my $DadSon_HASH=GenBankHandle::BuildAll_Dad_son_relationShip_HASH($AllGeneBkHASH);  #DirFileHandle::PrintDumper("0_0_4_DadSon_HASH.hash",$DadSon_HASH) if (   (  defined ( $DadSon_HASH )  ) && (  ref ( $DadSon_HASH ) eq 'HASH' )   );
	  my $SonDad_HASH=GenBankHandle::BuildAll_son_Dad_relationShip_HASH($AllGeneBkHASH);  #DirFileHandle::PrintDumper("0_0_5_SonDad_HASH.hash",$SonDad_HASH) if (   (  defined ( $SonDad_HASH )  ) && (  ref ( $SonDad_HASH ) eq 'HASH' )   );
	  
	  
	  
    
	  #my $SpecCountHASH; $SpecCountHASH=GenBankHandle::CountAllPointFormRoot_to_leaf_for_many_SPECIES($SonDad_HASH, $AllGeneBkHASH, $SpecCountHASH);  #DirFileHandle::PrintDumper("0_0_6_SpecCountHASH.hash",$SpecCountHASH) if (   (  defined ( $SpecCountHASH )  ) && (  ref ( $SpecCountHASH ) eq 'HASH' )   );
	  
	 
	  
	  my $Tree_outHASH; my $new_outTreeWord;  
	  my $new_outSmallerTreeTXT; my $new_outSmallerTreeTXT_with_leave;
	  my $new_outSmallerTreeHSH; my $new_outSmallerTreeHSH_with_leave;
	  if (   (  defined ( $AllPHY_HASH )  ) && (  ref ( $AllPHY_HASH ) eq 'HASH' )  && (  defined ( $DadSon_HASH )  ) && (  ref ( $DadSon_HASH ) eq 'HASH' )   ){
	  	
	  	
	    $Tree_outHASH=ClustalwRun::MakeTreeFrom_DadSonHashs_onlyMultiLevesTree($AllPHY_HASH, $DadSon_HASH);
	    
	    #my $realDad_son_hash=$SonDad_HASH; 
	    #$realDad_son_hash=GenBankHandle::Build_son_dad_hash_from_TREEHASH  ($Tree_outHASH, $realDad_son_hash);  DirFileHandle::PrintDumper("realDad_son_hash.hash",$realDad_son_hash) if (   (  defined ( $realDad_son_hash )  ) && (  ref ( $realDad_son_hash ) eq 'HASH' )   );
	    #my $realson_Dad_hash; $realson_Dad_hash=GenBankHandle::Build_dad_son_hash_from_TREEHASH  ($Tree_outHASH, $realson_Dad_hash);  DirFileHandle::PrintDumper("realson_Dad_hash.hash",$realson_Dad_hash) if (   (  defined ( $realson_Dad_hash )  ) && (  ref ( $realson_Dad_hash ) eq 'HASH' )   );
	    #my $real_SpecCountHASH; $real_SpecCountHASH=GenBankHandle::CountAllPointFormRoot_to_leaf_for_many_SPECIES($realDad_son_hash, $AllGeneBkHASH, $real_SpecCountHASH);  #DirFileHandle::PrintDumper("0_0_7_real_SpecCountHASH.hash",$real_SpecCountHASH) if (   (  defined ( $real_SpecCountHASH )  ) && (  ref ( $real_SpecCountHASH ) eq 'HASH' )   );
	    
	    my $real_SpecCountHASH; $real_SpecCountHASH=GenBankHandle::CountAllPointFormRoot_to_leaf_for_many_SPECIES      ($SonDad_HASH, $AllGeneBkHASH, $real_SpecCountHASH);
	    my $SpecCountHASH_wthU; $SpecCountHASH_wthU=GenBankHandle::CountAllPointFormRoot_to_leaf_for_many_SPECIES_for_U($SonDad_HASH, $AllGeneBkHASH, $SpecCountHASH_wthU);
	    
	    #my $outTreeWord; $outTreeWord=ClustalwRun::PrintTreePoints ($outTreeWord, $Tree_outHASH);
	    $new_outTreeWord        =GenBankHandle::new_PrintTreePoints_more_Inform              ($new_outTreeWord, $Tree_outHASH,                              $real_SpecCountHASH, $SpecCountHASH_wthU);
	    #DirFileHandle::PrintDumper ("1.$ptWd.real_SpecCountHASH.hash", $real_SpecCountHASH) if (   ref ($real_SpecCountHASH)    eq 'HASH' );
      #DirFileHandle::PrintDumper ("2.$ptWd.SpecCountHASH_wthU.hash", $SpecCountHASH_wthU) if (   ref ($SpecCountHASH_wthU)    eq 'HASH' );

	    #$new_outSmallerTreeTXT             =GenBankHandle::new_PrintTreePoints_more_Inform_less_showing                  ($new_outSmallerTreeTXT,            $Tree_outHASH,  $showing_taxonomy_Ifm_HASH,               $real_SpecCountHASH, $SpecCountHASH_wthU);
	    #$new_outSmallerTreeTXT_with_leave  =GenBankHandle::new_PrintTreePoints_more_Inform_less_showing_with_leave_point ($new_outSmallerTreeTXT_with_leave, $Tree_outHASH,  $showing_taxonomy_Ifm_HASH, $speciseHASH, $real_SpecCountHASH, $SpecCountHASH_wthU);
	    
	    #$new_outSmallerTreeHSH             =GenBankHandle::new_GetTreePoints_Hash_more_Inform_less_showing                  ($new_outSmallerTreeHSH,            $Tree_outHASH,  $showing_taxonomy_Ifm_HASH,               $real_SpecCountHASH, $SpecCountHASH_wthU);
	    #$new_outSmallerTreeHSH_with_leave  =GenBankHandle::new_GetTreePoints_Hash_more_Inform_less_showing_with_leave_point ($new_outSmallerTreeHSH_with_leave, $Tree_outHASH,  $showing_taxonomy_Ifm_HASH, $speciseHASH, $real_SpecCountHASH, $SpecCountHASH_wthU);
	    
	    my $MutCharSepcCountHASH=           GenBankHandle::CountAllPointFormRoot_to_leaf_for_many_SPECIES_for_many_Char    ($SonDad_HASH, $AllGeneBkHASH, $PtNbID_toChar_HASH); 
	    
	    $new_outSmallerTreeTXT             =GenBankHandle::Group_PrintTreePoints_more_Inform_less_showing                  ($new_outSmallerTreeTXT,            $Tree_outHASH,  $showing_taxonomy_Ifm_HASH,               $MutCharSepcCountHASH);
	    $new_outSmallerTreeTXT_with_leave  =GenBankHandle::Group_PrintTreePoints_more_Inform_less_showing_with_leave_point ($new_outSmallerTreeTXT_with_leave, $Tree_outHASH,  $showing_taxonomy_Ifm_HASH, $speciseHASH, $MutCharSepcCountHASH);
	    
	    $new_outSmallerTreeHSH             =GenBankHandle::Group_GetTreePoints_Hash_more_Inform_less_showing                  ($new_outSmallerTreeHSH,            $Tree_outHASH,  $showing_taxonomy_Ifm_HASH,               $MutCharSepcCountHASH);
	    $new_outSmallerTreeHSH_with_leave  =GenBankHandle::Group_GetTreePoints_Hash_more_Inform_less_showing_with_leave_point ($new_outSmallerTreeHSH_with_leave, $Tree_outHASH,  $showing_taxonomy_Ifm_HASH, $speciseHASH, $MutCharSepcCountHASH);
	    
	    
	    
      #print "\n\n\$new_outTreeWord=\n\n$new_outTreeWord\n\n\n"; 
      
      
      
    }
    
    
    
    if (   (  defined ( $Tree_outHASH )  ) && (  ref ( $Tree_outHASH ) eq 'HASH' )   ){
      #return $Tree_outHASH;
      
      return [$new_outTreeWord, $new_outSmallerTreeTXT, $new_outSmallerTreeTXT_with_leave, $new_outSmallerTreeHSH, $new_outSmallerTreeHSH_with_leave, $Tree_outHASH];
    }
    else{
    	my $warnMsg=$warnMsgHead.$callerMSG."\$Tree_outHASH=$Tree_outHASH is not right!! : $!\n\n\n";
    }
		
		
	}
	
	
	
	
  
  
}


sub BuildTreesHash_from_geneBankFiles{ #GenBankHandle::BuildTreesHash_from_geneBankFiles
	my ($inDIR, $top_howMany)=@_;
	
	my $warnMsgBody="\nIn package  GenBankHandle,\tIn sub BuildTreesHash_from_geneBankFiles,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $callerMSG=DirFileHandle::print_SubCallerInform;
	
	
	
	my ($AllGeneBkHASH, $AllRepeatRptH, $all_Pt_homologesHASH) =@{ GenBankHandle::Find_repeat_and_doing_statistic ($inDIR, $top_howMany) };
	#DirFileHandle::PrintDumper("0_0_0_AllGeneBkHASH.hash",$AllGeneBkHASH);
	#DirFileHandle::PrintDumper("0_0_1_AllRepeatRptH.hash",$AllRepeatRptH);
	
	
	#my $fileARRAY=DirFileHandle::getDirArray($inDIR);
	my $speciseHASH=GenBankHandle::BuildAllSpecise_HASH              ($AllGeneBkHASH);  #DirFileHandle::PrintDumper("0_0_2_speciseHASH.hash",$speciseHASH) if (   (  defined ( $speciseHASH )  ) && (  ref ( $speciseHASH ) eq 'HASH' )   );
	my $AllPHY_HASH=GenBankHandle::BuildAllTxnomy_phy_HASH           ($AllGeneBkHASH);  #DirFileHandle::PrintDumper("0_0_3_AllPHY_HASH.hash",$AllPHY_HASH) if (   (  defined ( $AllPHY_HASH )  ) && (  ref ( $AllPHY_HASH ) eq 'HASH' )   );
	my $DadSon_HASH=GenBankHandle::BuildAll_Dad_son_relationShip_HASH($AllGeneBkHASH);  #DirFileHandle::PrintDumper("0_0_4_DadSon_HASH.hash",$DadSon_HASH) if (   (  defined ( $DadSon_HASH )  ) && (  ref ( $DadSon_HASH ) eq 'HASH' )   );
	my $SonDad_HASH=GenBankHandle::BuildAll_son_Dad_relationShip_HASH($AllGeneBkHASH);  #DirFileHandle::PrintDumper("0_0_5_SonDad_HASH.hash",$SonDad_HASH) if (   (  defined ( $SonDad_HASH )  ) && (  ref ( $SonDad_HASH ) eq 'HASH' )   );
	
	
	

	my $SpecCountHASH; $SpecCountHASH=GenBankHandle::CountAllPointFormRoot_to_leaf_for_many_SPECIES($SonDad_HASH, $AllGeneBkHASH, $SpecCountHASH);  #DirFileHandle::PrintDumper("0_0_6_SpecCountHASH.hash",$SpecCountHASH) if (   (  defined ( $SpecCountHASH )  ) && (  ref ( $SpecCountHASH ) eq 'HASH' )   );
	
	my $outFASTAstring;
	my $outPepIfm_HASH;
	
	my $Tree_outHASH; my $new_outTreeWord;
	if (   (  defined ( $AllPHY_HASH )  ) && (  ref ( $AllPHY_HASH ) eq 'HASH' )  && (  defined ( $DadSon_HASH )  ) && (  ref ( $DadSon_HASH ) eq 'HASH' )   ){
		
		
	  $Tree_outHASH=ClustalwRun::MakeTreeFrom_DadSonHashs_onlyMultiLevesTree($AllPHY_HASH, $DadSon_HASH);
	  
	  my $realDad_son_hash=$SonDad_HASH; 
	  #$realDad_son_hash=GenBankHandle::Build_son_dad_hash_from_TREEHASH  ($Tree_outHASH, $realDad_son_hash);  DirFileHandle::PrintDumper("realDad_son_hash.hash",$realDad_son_hash) if (   (  defined ( $realDad_son_hash )  ) && (  ref ( $realDad_son_hash ) eq 'HASH' )   );
	  #my $realson_Dad_hash; $realson_Dad_hash=GenBankHandle::Build_dad_son_hash_from_TREEHASH  ($Tree_outHASH, $realson_Dad_hash);  DirFileHandle::PrintDumper("realson_Dad_hash.hash",$realson_Dad_hash) if (   (  defined ( $realson_Dad_hash )  ) && (  ref ( $realson_Dad_hash ) eq 'HASH' )   );
	  my $real_SpecCountHASH; $real_SpecCountHASH=GenBankHandle::CountAllPointFormRoot_to_leaf_for_many_SPECIES($realDad_son_hash, $AllGeneBkHASH, $real_SpecCountHASH);  #DirFileHandle::PrintDumper("0_0_7_real_SpecCountHASH.hash",$real_SpecCountHASH) if (   (  defined ( $real_SpecCountHASH )  ) && (  ref ( $real_SpecCountHASH ) eq 'HASH' )   );
	  
	  #my $outTreeWord; $outTreeWord=ClustalwRun::PrintTreePoints ($outTreeWord, $Tree_outHASH);
	  $new_outTreeWord=GenBankHandle::new_PrintTreePoints ($new_outTreeWord, $Tree_outHASH, $real_SpecCountHASH);
    print "\n\n\$new_outTreeWord=\n\n$new_outTreeWord\n\n\n"; 
    
    
    my $tempARy=GenBankHandle::BuildFastaFile($AllGeneBkHASH, $real_SpecCountHASH); 
    ($outFASTAstring, $outPepIfm_HASH)=@{ $tempARy } if (   (  defined ( $tempARy )  ) && (  ref ( $tempARy ) eq 'ARRAY' )   );
  }
  
  
  
  if (   (  defined ( $Tree_outHASH )  ) && (  ref ( $Tree_outHASH ) eq 'HASH' )   ){
    #return $Tree_outHASH;
    
    return [$outFASTAstring, $outPepIfm_HASH, $new_outTreeWord, $Tree_outHASH, $AllGeneBkHASH, $AllRepeatRptH, $all_Pt_homologesHASH];
  }
  else{
  	my $warnMsg=$warnMsgHead.$callerMSG."\$Tree_outHASH=$Tree_outHASH is not right!! : $!\n\n\n";
  }
  
  
}


sub new_PrintTreePoints{  #GenBankHandle::new_PrintTreePoints  #打印出树形结构
  my ($OutLine, $InHash, $SpecCountHASH, $levNb)=@_;
  
  my $warnMsgBody="\nIn package  GenBankHandle,\tIn sub new_PrintTreePoints,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  if (defined ($levNb)){$levNb++;}else {$levNb=0;print "\n";}
  
  if (  ( ref ( $InHash ) ) eq 'HASH' ){
    foreach my $ecKey (    sort {$a cmp $b} (   keys (  %{ $InHash }  )   )    ){ 
      if ($ecKey eq '-!End-+-Point!-'){}
      else { 
      	for (my $i=0; $i<$levNb; $i++){ print "    ";$OutLine.="    "; } 
#      	print     "+--- $ecKey:$SpecCountHASH->{$ecKey}\n";  
      	$OutLine.="+--- $ecKey:$SpecCountHASH->{$ecKey}\n";
      }   
      $OutLine=&new_PrintTreePoints( $OutLine, $InHash->{$ecKey}, $SpecCountHASH, $levNb );
      
    }
  }
  else {
    
  }
  return $OutLine;
}

sub new_PrintTreePoints_more_Inform{  #  my $OutLine=GenBankHandle::new_PrintTreePoints_more_Inform ($OutLine, $InHash, $SpecCountHASH_1, $SpecCountHASH_2, $levNb);  #打印出树形结构
  my ($OutLine, $InHash, $SpecCountHASH_1, $SpecCountHASH_2, $levNb)=@_;
  
  my $warnMsgBody="\nIn package  GenBankHandle,\tIn sub new_PrintTreePoints_more_Inform,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  if (defined ($levNb)){$levNb++;}else {$levNb=0;print "\n";}
  
  if (  ( ref ( $InHash ) ) eq 'HASH' ){
    foreach my $ecKey (    sort {$a cmp $b} (   keys (  %{ $InHash }  )   )    ){ 
      if ($ecKey eq '-!End-+-Point!-'){}
      else { 
      	if (  1   ){
      		
      		for (my $i=0; $i<$levNb; $i++){ 
      	  	$OutLine.="    ";  print "    ";
          } 
          
          my $pt_1=0; $pt_1=$SpecCountHASH_1->{$ecKey} if (   (  defined ( $SpecCountHASH_1 )  ) && (  ref ( $SpecCountHASH_1 ) eq 'HASH'  ) && (  defined ( $SpecCountHASH_1->{$ecKey} )  ) && ( $SpecCountHASH_1->{$ecKey}=~m/\S+/  )   ); $pt_1=sprintf ("%4d",$pt_1); $pt_1="All".$pt_1;
          my $pt_2=0; $pt_2=$SpecCountHASH_2->{$ecKey} if (   (  defined ( $SpecCountHASH_2 )  ) && (  ref ( $SpecCountHASH_2 ) eq 'HASH'  ) && (  defined ( $SpecCountHASH_2->{$ecKey} )  ) && ( $SpecCountHASH_2->{$ecKey}=~m/\S+/  )   ); $pt_2=sprintf ("%4d",$pt_2); $pt_2="U  ".$pt_2;
          
      	  $OutLine.="+--- $ecKey: $pt_1,   $pt_2\n";
      	  	
      	  	
      	  }
      	
      }   
      $OutLine=&new_PrintTreePoints_more_Inform( $OutLine, $InHash->{$ecKey}, $SpecCountHASH_1, $SpecCountHASH_2, $levNb );
      
    }
  }
  else {
    
  }
  return $OutLine;
}

#PtNbID_toChar_HASH
sub Group_PrintTreePoints_more_Inform_less_showing {  #  my $OutLine=GenBankHandle::Group_PrintTreePoints_more_Inform_less_showing ($OutLine, $InHash, $showing_taxonomy_Ifm_HASH, $In_MutCharSepcCountHASH, $levNb);  #打印出树形结构
  my ($OutLine, $InHash, $showing_taxonomy_Ifm_HASH, $In_MutCharSepcCountHASH, $levNb)=@_;
  
  my $warnMsgBody="\nIn package  GenBankHandle,\tIn sub Group_PrintTreePoints_more_Inform_less_showing,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  if (defined ($levNb)){$levNb++;}else {$levNb=0;print "\n";}
  
  if (  ( ref ( $InHash ) ) eq 'HASH' ){
    foreach my $ecKey (    sort {$a cmp $b} (   keys (  %{ $InHash }  )   )    ){   print "\n20181224-0-0  \$levNb=$levNb \$ecKey=$ecKey\n";
      if ($ecKey eq '-!End-+-Point!-'){  print "\n20181224-0-1  \$levNb=$levNb \$ecKey=$ecKey\n";  }
      else {    print "\n20181224-0-2  \$levNb=$levNb \$ecKey=$ecKey\n"; 
      	if (   (  defined ( $showing_taxonomy_Ifm_HASH )  ) && (  ref ( $showing_taxonomy_Ifm_HASH ) eq 'HASH' )&& (  ref ( $showing_taxonomy_Ifm_HASH->{ $ecKey } ) eq 'HASH' )   ){
      		print "\n20181224-0-3  \$levNb=$levNb \$ecKey=$ecKey \$showing_taxonomy_Ifm_HASH->{ $ecKey }=$showing_taxonomy_Ifm_HASH->{ $ecKey }\n"; 
      		
      		
      		
      		if (   (  defined ( $In_MutCharSepcCountHASH )  ) && (  ref ( $In_MutCharSepcCountHASH ) eq 'HASH'  )   ){
      			
      			my $need_to_work=0;
      			foreach my $eachChar (    sort {$a cmp $b} (   keys (  %{ $In_MutCharSepcCountHASH }  )   )    ){
      				if (      (  defined ( $In_MutCharSepcCountHASH->{$eachChar} )  ) && (  ref ( $In_MutCharSepcCountHASH->{$eachChar} ) eq 'HASH'  ) 
      				       && ( $In_MutCharSepcCountHASH->{$eachChar}->{$ecKey}=~m/\d+/  ) && ( $In_MutCharSepcCountHASH->{$eachChar}->{$ecKey} > 0  )   
      				   )
      				{
      				  $need_to_work=1;   	
      				}
      			}
      			
      			if ( $need_to_work==1 ){
      				for (my $i=0; $i<$levNb; $i++){ 
      	    	  $OutLine.="    ";  #print "\n20181224-0-4 \$levNb=$levNb \$ecKey=$ecKey  \$ecKey=$ecKey \$showing_taxonomy_Ifm_HASH->{ $ecKey }=$showing_taxonomy_Ifm_HASH->{ $ecKey }  \$i=$i \$OutLine=\n$OutLine \n\n\n"; 
              } 
              $OutLine.="+--- $ecKey:\t";
              foreach my $eachChar (    sort {$a cmp $b} (   keys (  %{ $In_MutCharSepcCountHASH }  )   )    ){
      			  	my $printText='     ';
      			  	if (      (  defined ( $In_MutCharSepcCountHASH->{$eachChar} )  ) && (  ref ( $In_MutCharSepcCountHASH->{$eachChar} ) eq 'HASH'  ) 
      			  	       && ( $In_MutCharSepcCountHASH->{$eachChar}->{$ecKey}=~m/\d+/  ) && ( $In_MutCharSepcCountHASH->{$eachChar}->{$ecKey} > 0  )   
      			  	   )
      			  	{ 
      			  		my $CharNB=$In_MutCharSepcCountHASH->{$eachChar}->{$ecKey}; $CharNB=sprintf ("%4d",$CharNB); 
      			  	  $printText=$eachChar.$CharNB;
      			  	}
      			  	$OutLine.=$printText."\t";
      			  }  
      			  $OutLine.="\n";	
      			}
      			
      		}   
      	  	
      	}
      	
      }   
      $OutLine=&Group_PrintTreePoints_more_Inform_less_showing( $OutLine, $InHash->{$ecKey}, $showing_taxonomy_Ifm_HASH, $In_MutCharSepcCountHASH, $levNb );
      
    }
  }
  else {
    
  }
  return $OutLine;
}


sub Group_PrintTreePoints_more_Inform_less_showing_with_leave_point {  #  my $OutLine=GenBankHandle::Group_PrintTreePoints_more_Inform_less_showing_with_leave_point ($OutLine, $InHash, $showing_taxonomy_Ifm_HASH, $speciseHASH, $In_MutCharSepcCountHASH, $levNb);  #打印出树形结构
  my ($OutLine, $InHash, $showing_taxonomy_Ifm_HASH, $speciseHASH, $In_MutCharSepcCountHASH, $levNb)=@_;
  
  my $warnMsgBody="\nIn package  GenBankHandle,\tIn sub Group_PrintTreePoints_more_Inform_less_showing_with_leave_point,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  if (defined ($levNb)){$levNb++;}else {$levNb=0;print "\n";}
  
  if (  ( ref ( $InHash ) ) eq 'HASH' ){
    foreach my $ecKey (    sort {$a cmp $b} (   keys (  %{ $InHash }  )   )    ){   print "\n20181224-0-0  \$levNb=$levNb \$ecKey=$ecKey\n";
      if ($ecKey eq '-!End-+-Point!-'){  print "\n20181224-0-1  \$levNb=$levNb \$ecKey=$ecKey\n";  }
      else {    print "\n201812251443-0-2  \$levNb=$levNb \$ecKey=$ecKey\n"; 
      	
      	if (   (  defined ( $In_MutCharSepcCountHASH )  ) && (  ref ( $In_MutCharSepcCountHASH ) eq 'HASH'  )   ){
      	  my $need_to_work=0;
      	  foreach my $eachChar (    sort {$a cmp $b} (   keys (  %{ $In_MutCharSepcCountHASH }  )   )    ){
      	  	if (      (  defined ( $In_MutCharSepcCountHASH->{$eachChar} )  ) && (  ref ( $In_MutCharSepcCountHASH->{$eachChar} ) eq 'HASH'  ) 
      		       && ( $In_MutCharSepcCountHASH->{$eachChar}->{$ecKey}=~m/\d+/  ) && ( $In_MutCharSepcCountHASH->{$eachChar}->{$ecKey} > 0  )   
      	  	   )
        		{
      	  	  $need_to_work=1;   	
      	  	}
        	}
      			
      	  if (        (      (  defined ( $showing_taxonomy_Ifm_HASH )  ) && (  ref ( $showing_taxonomy_Ifm_HASH ) eq 'HASH' ) 
      	                  && (  ref ( $showing_taxonomy_Ifm_HASH->{ $ecKey } ) eq 'HASH' )   
      	              )
      	          ||  
      	              (      (  defined ( $speciseHASH )  ) && (  ref ( $speciseHASH ) eq 'HASH' ) && ( $speciseHASH->{ $ecKey }=~m/\S+/ )   
      	              )
      	     )
      	  {
      		  print "\n201812251443-0-3  \$levNb=$levNb \$ecKey=$ecKey \$speciseHASH->{ $ecKey }=$speciseHASH->{ $ecKey }\n"; 
      		
      			if ( $need_to_work==1 ){
      				for (my $i=0; $i<$levNb; $i++){ 
      	    	  $OutLine.="    ";  #print "\n20181224-0-4 \$levNb=$levNb \$ecKey=$ecKey  \$ecKey=$ecKey \$showing_taxonomy_Ifm_HASH->{ $ecKey }=$showing_taxonomy_Ifm_HASH->{ $ecKey }  \$i=$i \$OutLine=\n$OutLine \n\n\n"; 
              } 
              $OutLine.="+--- $ecKey:\t";
              foreach my $eachChar (    sort {$a cmp $b} (   keys (  %{ $In_MutCharSepcCountHASH }  )   )    ){
      			  	my $printText='     ';
      			  	if (      (  defined ( $In_MutCharSepcCountHASH->{$eachChar} )  ) && (  ref ( $In_MutCharSepcCountHASH->{$eachChar} ) eq 'HASH'  ) 
      			  	       && ( $In_MutCharSepcCountHASH->{$eachChar}->{$ecKey}=~m/\d+/  ) && ( $In_MutCharSepcCountHASH->{$eachChar}->{$ecKey} > 0  )   
      			  	   )
      			  	{ 
      			  		my $CharNB=$In_MutCharSepcCountHASH->{$eachChar}->{$ecKey}; $CharNB=sprintf ("%4d",$CharNB); 
      			  	  $printText=$eachChar.$CharNB;
      			  	}
      			  	$OutLine.=$printText."\t";
      			  }  
      			  $OutLine.="\n";	
      			}
      			
      		} 
      	  ##################
      		     	  	
      	}
      	
      	
      	  
      	
      	
      }   
      $OutLine=&Group_PrintTreePoints_more_Inform_less_showing_with_leave_point( $OutLine, $InHash->{$ecKey}, $showing_taxonomy_Ifm_HASH, $speciseHASH, $In_MutCharSepcCountHASH, $levNb );
      
    }
  }
  else {
    
  }
  return $OutLine;
}


sub Group_GetTreePoints_Hash_more_Inform_less_showing {  #  my $OutHash=GenBankHandle::Group_GetTreePoints_Hash_more_Inform_less_showing ($OutHash, $InHash, $showing_taxonomy_Ifm_HASH, $In_MutCharSepcCountHASH, $levNb);  #打印出树形结构
  my ($OutHash, $InHash, $showing_taxonomy_Ifm_HASH, $In_MutCharSepcCountHASH, $levNb)=@_;
  
  my $warnMsgBody="\nIn package  GenBankHandle,\tIn sub Group_GetTreePoints_Hash_more_Inform_less_showing,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  if (defined ($levNb)){$levNb++;}else {$levNb=0;print "\n";}
  
  if (  ( ref ( $InHash ) ) eq 'HASH' ){
    foreach my $ecKey (    sort {$a cmp $b} (   keys (  %{ $InHash }  )   )    ){   #print "\n20181224-0-0  \$levNb=$levNb \$ecKey=$ecKey\n";
      if ($ecKey eq '-!End-+-Point!-'){  print "\n20181224-0-1  \$levNb=$levNb \$ecKey=$ecKey\n";  }
      else {    #print "\n20181224-0-2  \$levNb=$levNb \$ecKey=$ecKey\n"; 
      	if (   (  defined ( $showing_taxonomy_Ifm_HASH )  ) && (  ref ( $showing_taxonomy_Ifm_HASH ) eq 'HASH' )&& (  ref ( $showing_taxonomy_Ifm_HASH->{ $ecKey } ) eq 'HASH' )   ){
      		#print "\n20181224-0-3  \$levNb=$levNb \$ecKey=$ecKey \$showing_taxonomy_Ifm_HASH->{ $ecKey }=$showing_taxonomy_Ifm_HASH->{ $ecKey }\n"; 
      		
      		
      		######################
      		if (   (  defined ( $In_MutCharSepcCountHASH )  ) && (  ref ( $In_MutCharSepcCountHASH ) eq 'HASH'  )   ){      			
      			
      			foreach my $eachChar (    sort {$a cmp $b} (   keys (  %{ $In_MutCharSepcCountHASH }  )   )    ){
      				if (      (  defined ( $In_MutCharSepcCountHASH->{$eachChar} )  ) && (  ref ( $In_MutCharSepcCountHASH->{$eachChar} ) eq 'HASH'  ) 
      				       && ( $In_MutCharSepcCountHASH->{$eachChar}->{$ecKey}=~m/\d+/  ) && ( $In_MutCharSepcCountHASH->{$eachChar}->{$ecKey} > 0  )   
      				   )
      				{
      				  my $CharNB=$In_MutCharSepcCountHASH->{$eachChar}->{$ecKey}; 
      				  
      				  #$OutHash->{$ecKey}->{'0_level'}=$levNb;
      	  	    $OutHash->{$ecKey}->{$eachChar}=$CharNB;    	        	
      				  
      				  
      				}
      			}
      			
      			
      		} 
      		######################
      		      		   		
      	}      	
      }   
      $OutHash=&Group_GetTreePoints_Hash_more_Inform_less_showing( $OutHash, $InHash->{$ecKey}, $showing_taxonomy_Ifm_HASH, $In_MutCharSepcCountHASH, $levNb );
      
    }
  }
  else {
    
  }
  return $OutHash;
}

sub Group_GetTreePoints_Hash_more_Inform_less_showing_with_leave_point {  #  my $OutHash=GenBankHandle::Group_GetTreePoints_Hash_more_Inform_less_showing_with_leave_point ($OutHash, $InHash, $showing_taxonomy_Ifm_HASH, $speciseHASH, $In_MutCharSepcCountHASH, $levNb);  #打印出树形结构
  my ($OutHash, $InHash, $showing_taxonomy_Ifm_HASH, $speciseHASH, $In_MutCharSepcCountHASH, $levNb)=@_;
  
  my $warnMsgBody="\nIn package  GenBankHandle,\tIn sub Group_GetTreePoints_Hash_more_Inform_less_showing_with_leave_point,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  if (defined ($levNb)){$levNb++;}else {$levNb=0;print "\n";}
  
  if (  ( ref ( $InHash ) ) eq 'HASH' ){
    foreach my $ecKey (    sort {$a cmp $b} (   keys (  %{ $InHash }  )   )    ){   #print "\n20181224-0-0  \$levNb=$levNb \$ecKey=$ecKey\n";
      if ($ecKey eq '-!End-+-Point!-'){  print "\n20181224-0-1  \$levNb=$levNb \$ecKey=$ecKey\n";  }
      else {    #print "\n20181224-0-2  \$levNb=$levNb \$ecKey=$ecKey\n"; 
      	
      	if (        (      (  defined ( $showing_taxonomy_Ifm_HASH )  ) && (  ref ( $showing_taxonomy_Ifm_HASH ) eq 'HASH' ) 
      	                && (  ref ( $showing_taxonomy_Ifm_HASH->{ $ecKey } ) eq 'HASH' )   
      	            )
      	        ||  
      	            (      (  defined ( $speciseHASH )  ) && (  ref ( $speciseHASH ) eq 'HASH' ) && ( $speciseHASH->{ $ecKey }=~m/\S+/ ) 
      	            )
      	   )
      	{
      		#print "\n20181224-0-3  \$levNb=$levNb \$ecKey=$ecKey \$showing_taxonomy_Ifm_HASH->{ $ecKey }=$showing_taxonomy_Ifm_HASH->{ $ecKey }\n"; 
      		
      		######################
      		if (   (  defined ( $In_MutCharSepcCountHASH )  ) && (  ref ( $In_MutCharSepcCountHASH ) eq 'HASH'  )   ){      			
      			
      			foreach my $eachChar (    sort {$a cmp $b} (   keys (  %{ $In_MutCharSepcCountHASH }  )   )    ){
      				if (      (  defined ( $In_MutCharSepcCountHASH->{$eachChar} )  ) && (  ref ( $In_MutCharSepcCountHASH->{$eachChar} ) eq 'HASH'  ) 
      				       && ( $In_MutCharSepcCountHASH->{$eachChar}->{$ecKey}=~m/\d+/  ) && ( $In_MutCharSepcCountHASH->{$eachChar}->{$ecKey} > 0  )   
      				   )
      				{
      				  my $CharNB=$In_MutCharSepcCountHASH->{$eachChar}->{$ecKey}; 
      				  
      				  #$OutHash->{$ecKey}->{'0_level'}=$levNb;
      	  	    $OutHash->{$ecKey}->{$eachChar}=$CharNB;    	        	
      				  
      				  
      				}
      			}
      			
      			
      		} 
      		######################
      		     	  	
      	}
      	
      }   
      $OutHash=&Group_GetTreePoints_Hash_more_Inform_less_showing_with_leave_point( $OutHash, $InHash->{$ecKey}, $showing_taxonomy_Ifm_HASH, $speciseHASH, $In_MutCharSepcCountHASH, $levNb );
      
    }
  }
  else {
    
  }
  return $OutHash;
}




sub new_PrintTreePoints_more_Inform_less_showing {  #  my $OutLine=GenBankHandle::new_PrintTreePoints_more_Inform_less_showing ($OutLine, $InHash, $showing_taxonomy_Ifm_HASH, $SpecCountHASH_1, $SpecCountHASH_2, $levNb);  #打印出树形结构
  my ($OutLine, $InHash, $showing_taxonomy_Ifm_HASH, $SpecCountHASH_1, $SpecCountHASH_2, $levNb)=@_;
  
  my $warnMsgBody="\nIn package  GenBankHandle,\tIn sub new_PrintTreePoints_more_Inform_less_showing,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  if (defined ($levNb)){$levNb++;}else {$levNb=0;print "\n";}
  
  if (  ( ref ( $InHash ) ) eq 'HASH' ){
    foreach my $ecKey (    sort {$a cmp $b} (   keys (  %{ $InHash }  )   )    ){   print "\n20181224-0-0  \$levNb=$levNb \$ecKey=$ecKey\n";
      if ($ecKey eq '-!End-+-Point!-'){  print "\n20181224-0-1  \$levNb=$levNb \$ecKey=$ecKey\n";  }
      else {    print "\n20181224-0-2  \$levNb=$levNb \$ecKey=$ecKey\n"; 
      	if (   (  defined ( $showing_taxonomy_Ifm_HASH )  ) && (  ref ( $showing_taxonomy_Ifm_HASH ) eq 'HASH' )&& (  ref ( $showing_taxonomy_Ifm_HASH->{ $ecKey } ) eq 'HASH' )   ){
      		print "\n20181224-0-3  \$levNb=$levNb \$ecKey=$ecKey \$showing_taxonomy_Ifm_HASH->{ $ecKey }=$showing_taxonomy_Ifm_HASH->{ $ecKey }\n"; 
      		
      		if (   (  defined ( $SpecCountHASH_2 )  ) && (  ref ( $SpecCountHASH_2 ) eq 'HASH'  ) && (  defined ( $SpecCountHASH_2->{$ecKey} )  ) && ( $SpecCountHASH_2->{$ecKey}=~m/\d+/  ) && ( $SpecCountHASH_2->{$ecKey}>0  )   ){
      			
      			for (my $i=0; $i<$levNb; $i++){ 
      	    	$OutLine.="    ";  print "\n20181224-0-4 \$levNb=$levNb \$ecKey=$ecKey  \$ecKey=$ecKey \$showing_taxonomy_Ifm_HASH->{ $ecKey }=$showing_taxonomy_Ifm_HASH->{ $ecKey }  \$i=$i \$OutLine=\n$OutLine \n\n\n"; 
            } 
            
            my $pt_1=0; $pt_1=$SpecCountHASH_1->{$ecKey} if (   (  defined ( $SpecCountHASH_1 )  ) && (  ref ( $SpecCountHASH_1 ) eq 'HASH'  ) && (  defined ( $SpecCountHASH_1->{$ecKey} )  ) && ( $SpecCountHASH_1->{$ecKey}=~m/\S+/  )   ); $pt_1=sprintf ("%4d",$pt_1); $pt_1="All\t".$pt_1;
            my $pt_2=0; $pt_2=$SpecCountHASH_2->{$ecKey} if (   (  defined ( $SpecCountHASH_2 )  ) && (  ref ( $SpecCountHASH_2 ) eq 'HASH'  ) && (  defined ( $SpecCountHASH_2->{$ecKey} )  ) && ( $SpecCountHASH_2->{$ecKey}=~m/\S+/  )   ); $pt_2=sprintf ("%4d",$pt_2); $pt_2="U\t"  .$pt_2;
            
      	    $OutLine.="+--- $ecKey:\t$pt_1\t$pt_2\n"; print "\n20181224-0-5 \$levNb=$levNb \$ecKey=$ecKey \$showing_taxonomy_Ifm_HASH->{ $ecKey }=$showing_taxonomy_Ifm_HASH->{ $ecKey }  \$OutLine=\n$OutLine \n\n\n\n\n"; 
      	  	
      			
      			
      		}
      		
      		
      	  	
      	}
      	
      }   
      $OutLine=&new_PrintTreePoints_more_Inform_less_showing( $OutLine, $InHash->{$ecKey}, $showing_taxonomy_Ifm_HASH, $SpecCountHASH_1, $SpecCountHASH_2, $levNb );
      
    }
  }
  else {
    
  }
  return $OutLine;
}


sub new_PrintTreePoints_more_Inform_less_showing_with_leave_point {  #  my $OutLine=GenBankHandle::new_PrintTreePoints_more_Inform_less_showing_with_leave_point ($OutLine, $InHash, $showing_taxonomy_Ifm_HASH, $speciseHASH, $SpecCountHASH_1, $SpecCountHASH_2, $levNb);  #打印出树形结构
  my ($OutLine, $InHash, $showing_taxonomy_Ifm_HASH, $speciseHASH, $SpecCountHASH_1, $SpecCountHASH_2, $levNb)=@_;
  
  my $warnMsgBody="\nIn package  GenBankHandle,\tIn sub new_PrintTreePoints_more_Inform_less_showing_with_leave_point,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  if (defined ($levNb)){$levNb++;}else {$levNb=0;print "\n";}
  
  if (  ( ref ( $InHash ) ) eq 'HASH' ){
    foreach my $ecKey (    sort {$a cmp $b} (   keys (  %{ $InHash }  )   )    ){   print "\n20181224-0-0  \$levNb=$levNb \$ecKey=$ecKey\n";
      if ($ecKey eq '-!End-+-Point!-'){  print "\n20181224-0-1  \$levNb=$levNb \$ecKey=$ecKey\n";  }
      else {    print "\n201812251443-0-2  \$levNb=$levNb \$ecKey=$ecKey\n"; 
      	
      	if (        (      (  defined ( $showing_taxonomy_Ifm_HASH )  ) && (  ref ( $showing_taxonomy_Ifm_HASH ) eq 'HASH' ) 
      	                && (  ref ( $showing_taxonomy_Ifm_HASH->{ $ecKey } ) eq 'HASH' )   
      	            )
      	        ||  
      	            (      (  defined ( $speciseHASH )  ) && (  ref ( $speciseHASH ) eq 'HASH' ) && ( $speciseHASH->{ $ecKey }=~m/\S+/ )   
      	            )
      	   )
      	{
      		print "\n201812251443-0-3  \$levNb=$levNb \$ecKey=$ecKey \$speciseHASH->{ $ecKey }=$speciseHASH->{ $ecKey }\n"; 
      		
      		if (   (  defined ( $SpecCountHASH_2 )  ) && (  ref ( $SpecCountHASH_2 ) eq 'HASH'  ) && (  defined ( $SpecCountHASH_2->{$ecKey} )  ) && ( $SpecCountHASH_2->{$ecKey}=~m/\d+/  ) && ( $SpecCountHASH_2->{$ecKey}>0  )   ){
      			
      			for (my $i=0; $i<$levNb; $i++){ 
      	    	$OutLine.="    ";  #print "\n20181224-0-4 \$levNb=$levNb \$ecKey=$ecKey  \$ecKey=$ecKey \$showing_taxonomy_Ifm_HASH->{ $ecKey }=$showing_taxonomy_Ifm_HASH->{ $ecKey }  \$i=$i \$OutLine=\n$OutLine \n\n\n"; 
            } 
            
            my $pt_1=0; $pt_1=$SpecCountHASH_1->{$ecKey} if (   (  defined ( $SpecCountHASH_1 )  ) && (  ref ( $SpecCountHASH_1 ) eq 'HASH'  ) && (  defined ( $SpecCountHASH_1->{$ecKey} )  ) && ( $SpecCountHASH_1->{$ecKey}=~m/\S+/  )   ); $pt_1=sprintf ("%4d",$pt_1); $pt_1="All\t".$pt_1;
            my $pt_2=0; $pt_2=$SpecCountHASH_2->{$ecKey} if (   (  defined ( $SpecCountHASH_2 )  ) && (  ref ( $SpecCountHASH_2 ) eq 'HASH'  ) && (  defined ( $SpecCountHASH_2->{$ecKey} )  ) && ( $SpecCountHASH_2->{$ecKey}=~m/\S+/  )   ); $pt_2=sprintf ("%4d",$pt_2); $pt_2="U\t"  .$pt_2;
            
      	    $OutLine.="+--- $ecKey:\t$pt_1\t$pt_2\n"; print "\n201812251443-0-5 \$levNb=$levNb \$ecKey=$ecKey \$speciseHASH->{ $ecKey }=$speciseHASH->{ $ecKey }  \$OutLine=\n$OutLine \n\n\n\n\n"; 
      	  	     			
      		}
      		     	  	
      	}
      	
      }   
      $OutLine=&new_PrintTreePoints_more_Inform_less_showing_with_leave_point( $OutLine, $InHash->{$ecKey}, $showing_taxonomy_Ifm_HASH, $speciseHASH, $SpecCountHASH_1, $SpecCountHASH_2, $levNb );
      
    }
  }
  else {
    
  }
  return $OutLine;
}


sub new_GetTreePoints_Hash_more_Inform_less_showing {  #  my $OutHash=GenBankHandle::new_GetTreePoints_Hash_more_Inform_less_showing ($OutHash, $InHash, $showing_taxonomy_Ifm_HASH, $SpecCountHASH_1, $SpecCountHASH_2, $levNb);  #打印出树形结构
  my ($OutHash, $InHash, $showing_taxonomy_Ifm_HASH, $SpecCountHASH_1, $SpecCountHASH_2, $levNb)=@_;
  
  my $warnMsgBody="\nIn package  GenBankHandle,\tIn sub new_GetTreePoints_Hash_more_Inform_less_showing,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  if (defined ($levNb)){$levNb++;}else {$levNb=0;print "\n";}
  
  if (  ( ref ( $InHash ) ) eq 'HASH' ){
    foreach my $ecKey (    sort {$a cmp $b} (   keys (  %{ $InHash }  )   )    ){   #print "\n20181224-0-0  \$levNb=$levNb \$ecKey=$ecKey\n";
      if ($ecKey eq '-!End-+-Point!-'){  print "\n20181224-0-1  \$levNb=$levNb \$ecKey=$ecKey\n";  }
      else {    #print "\n20181224-0-2  \$levNb=$levNb \$ecKey=$ecKey\n"; 
      	if (   (  defined ( $showing_taxonomy_Ifm_HASH )  ) && (  ref ( $showing_taxonomy_Ifm_HASH ) eq 'HASH' )&& (  ref ( $showing_taxonomy_Ifm_HASH->{ $ecKey } ) eq 'HASH' )   ){
      		#print "\n20181224-0-3  \$levNb=$levNb \$ecKey=$ecKey \$showing_taxonomy_Ifm_HASH->{ $ecKey }=$showing_taxonomy_Ifm_HASH->{ $ecKey }\n"; 
      		
      		if (   (  defined ( $SpecCountHASH_2 )  ) && (  ref ( $SpecCountHASH_2 ) eq 'HASH'  ) && (  defined ( $SpecCountHASH_2->{$ecKey} )  ) && ( $SpecCountHASH_2->{$ecKey}=~m/\d+/  ) && ( $SpecCountHASH_2->{$ecKey}>0  )   ){
      			
            my $pt_1=0; $pt_1=$SpecCountHASH_1->{$ecKey} if (   (  defined ( $SpecCountHASH_1 )  ) && (  ref ( $SpecCountHASH_1 ) eq 'HASH'  ) && (  defined ( $SpecCountHASH_1->{$ecKey} )  ) && ( $SpecCountHASH_1->{$ecKey}=~m/\S+/  )   );
            my $pt_2=0; $pt_2=$SpecCountHASH_2->{$ecKey} if (   (  defined ( $SpecCountHASH_2 )  ) && (  ref ( $SpecCountHASH_2 ) eq 'HASH'  ) && (  defined ( $SpecCountHASH_2->{$ecKey} )  ) && ( $SpecCountHASH_2->{$ecKey}=~m/\S+/  )   );
            
      	  	$OutHash->{$ecKey}->{'0_level'}=$levNb;
      	  	$OutHash->{$ecKey}->{'1_allNB'}=$pt_1;
      	  	$OutHash->{$ecKey}->{'2_U__NB'}=$pt_2;   	  	
      			
      		}      		
      	}      	
      }   
      $OutHash=&new_GetTreePoints_Hash_more_Inform_less_showing( $OutHash, $InHash->{$ecKey}, $showing_taxonomy_Ifm_HASH, $SpecCountHASH_1, $SpecCountHASH_2, $levNb );
      
    }
  }
  else {
    
  }
  return $OutHash;
}

sub new_GetTreePoints_Hash_more_Inform_less_showing_with_leave_point {  #  my $OutHash=GenBankHandle::new_GetTreePoints_Hash_more_Inform_less_showing_with_leave_point ($OutHash, $InHash, $showing_taxonomy_Ifm_HASH, $speciseHASH, $SpecCountHASH_1, $SpecCountHASH_2, $levNb);  #打印出树形结构
  my ($OutHash, $InHash, $showing_taxonomy_Ifm_HASH, $speciseHASH, $SpecCountHASH_1, $SpecCountHASH_2, $levNb)=@_;
  
  my $warnMsgBody="\nIn package  GenBankHandle,\tIn sub new_GetTreePoints_Hash_more_Inform_less_showing_with_leave_point,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  if (defined ($levNb)){$levNb++;}else {$levNb=0;print "\n";}
  
  if (  ( ref ( $InHash ) ) eq 'HASH' ){
    foreach my $ecKey (    sort {$a cmp $b} (   keys (  %{ $InHash }  )   )    ){   #print "\n20181224-0-0  \$levNb=$levNb \$ecKey=$ecKey\n";
      if ($ecKey eq '-!End-+-Point!-'){  print "\n20181224-0-1  \$levNb=$levNb \$ecKey=$ecKey\n";  }
      else {    #print "\n20181224-0-2  \$levNb=$levNb \$ecKey=$ecKey\n"; 
      	
      	if (        (      (  defined ( $showing_taxonomy_Ifm_HASH )  ) && (  ref ( $showing_taxonomy_Ifm_HASH ) eq 'HASH' ) 
      	                && (  ref ( $showing_taxonomy_Ifm_HASH->{ $ecKey } ) eq 'HASH' )   
      	            )
      	        ||  
      	            (      (  defined ( $speciseHASH )  ) && (  ref ( $speciseHASH ) eq 'HASH' ) && ( $speciseHASH->{ $ecKey }=~m/\S+/ ) 
      	            )
      	   )
      	{
      		#print "\n20181224-0-3  \$levNb=$levNb \$ecKey=$ecKey \$showing_taxonomy_Ifm_HASH->{ $ecKey }=$showing_taxonomy_Ifm_HASH->{ $ecKey }\n"; 
      		
      		if (   (  defined ( $SpecCountHASH_2 )  ) && (  ref ( $SpecCountHASH_2 ) eq 'HASH'  ) && (  defined ( $SpecCountHASH_2->{$ecKey} )  ) && ( $SpecCountHASH_2->{$ecKey}=~m/\d+/  ) && ( $SpecCountHASH_2->{$ecKey}>0  )   ){
      			
      			my $pt_1=0; $pt_1=$SpecCountHASH_1->{$ecKey} if (   (  defined ( $SpecCountHASH_1 )  ) && (  ref ( $SpecCountHASH_1 ) eq 'HASH'  ) && (  defined ( $SpecCountHASH_1->{$ecKey} )  ) && ( $SpecCountHASH_1->{$ecKey}=~m/\S+/  )   );
            my $pt_2=0; $pt_2=$SpecCountHASH_2->{$ecKey} if (   (  defined ( $SpecCountHASH_2 )  ) && (  ref ( $SpecCountHASH_2 ) eq 'HASH'  ) && (  defined ( $SpecCountHASH_2->{$ecKey} )  ) && ( $SpecCountHASH_2->{$ecKey}=~m/\S+/  )   );
                     
            $OutHash->{$ecKey}->{'0_level'}=$levNb;
      	  	$OutHash->{$ecKey}->{'1_allNB'}=$pt_1;
      	  	$OutHash->{$ecKey}->{'2_U__NB'}=$pt_2;    			
      		}
      		     	  	
      	}
      	
      }   
      $OutHash=&new_GetTreePoints_Hash_more_Inform_less_showing_with_leave_point( $OutHash, $InHash->{$ecKey}, $showing_taxonomy_Ifm_HASH, $speciseHASH, $SpecCountHASH_1, $SpecCountHASH_2, $levNb );
      
    }
  }
  else {
    
  }
  return $OutHash;
}

sub Build_dad_son_hash_from_TREEHASH{  #    GenBankHandle::Build_dad_son_hash_from_TREEHASH  
  my ($treeHASH, $DadSonHash, $levNb)=@_;
  
  my $warnMsgBody="\nIn package  GenBankHandle,\tIn sub Build_dad_son_hash_from_TREEHASH,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  #my $callerMSG=DirFileHandle::print_SubCallerInform;
	
  if  (   (   defined ( $levNb )   ) && ( $levNb=~m/\d+/ ) && ( $levNb>=0 )   ){  	
  }
  else {
  	$levNb=0;
  }
  
  if (   (   defined ( $treeHASH )   ) && (  ref ( $treeHASH ) eq 'HASH'  )    ){
    foreach my $mid_1_Point (    sort { $a cmp $b } (   keys (  %{ $treeHASH }  )   )    ){
    	if (   (   defined ( $treeHASH->{$mid_1_Point} )   ) && (  ref ( $treeHASH->{$mid_1_Point} ) eq 'HASH'  )    ){
    		foreach my $mid_2_Point (    sort { $a cmp $b } (   keys (  %{ $treeHASH->{$mid_1_Point} }  )   )    ){
    		  if (   (   defined ( $treeHASH->{$mid_1_Point}->{$mid_2_Point} )   ) && (  ref ( $treeHASH->{$mid_1_Point}->{$mid_2_Point} ) eq 'HASH'  )    ){
    		    $DadSonHash->{$mid_1_Point}->{$mid_2_Point}=1; #my $MSG= "201811271329-1 \$levNb=$levNb \$DadSonHash->{$mid_1_Point}->{$mid_2_Point}=$DadSonHash->{$mid_1_Point}->{$mid_2_Point}\n"; warn $MSG; print $MSG;
    		    $DadSonHash=GenBankHandle::Build_dad_son_hash_from_TREEHASH($treeHASH->{$mid_1_Point}, $DadSonHash, $levNb);
    		    $levNb++;
    		  }
    		}
    		
    	}
    }
  }
  else{
  	
  }
  return $DadSonHash;
 
}


sub Build_son_dad_hash_from_TREEHASH{  #    GenBankHandle::Build_son_dad_hash_from_TREEHASH
  my ($treeHASH, $SonDadHash, $levNb)=@_;
  
  my $warnMsgBody="\nIn package  GenBankHandle,\tIn sub Build_son_dad_hash_from_TREEHASH,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  #my $callerMSG=DirFileHandle::print_SubCallerInform;
	
  if  (   (   defined ( $levNb )   ) && ( $levNb=~m/\d+/ ) && ( $levNb>=0 )   ){  	
  }
  else {
  	$levNb=0;
  }
  
  if (   (   defined ( $treeHASH )   ) && (  ref ( $treeHASH ) eq 'HASH'  )    ){
    foreach my $mid_1_Point (    sort { $a cmp $b } (   keys (  %{ $treeHASH }  )   )    ){
    	if (   (   defined ( $treeHASH->{$mid_1_Point} )   ) && (  ref ( $treeHASH->{$mid_1_Point} ) eq 'HASH'  )    ){
    		foreach my $mid_2_Point (    sort { $a cmp $b } (   keys (  %{ $treeHASH->{$mid_1_Point} }  )   )    ){
    		  if (   (   defined ( $treeHASH->{$mid_1_Point}->{$mid_2_Point} )   ) && (  ref ( $treeHASH->{$mid_1_Point}->{$mid_2_Point} ) eq 'HASH'  )    ){
    		    $SonDadHash->{$mid_2_Point}=$mid_1_Point;   #my $MSG= "201811271329-2 \$levNb=$levNb \$SonDadHash->{$mid_2_Point}=$mid_1_Point=$SonDadHash->{$mid_2_Point}\n"; warn $MSG; print $MSG;
    		    $SonDadHash=GenBankHandle::Build_son_dad_hash_from_TREEHASH($treeHASH->{$mid_1_Point}, $SonDadHash, $levNb);
    		    $levNb++;
    		  }
    		}
    		
    	}
    }
  }
  else{
  	
  }
  return $SonDadHash;
 
}

# my $MutCharSepcCountHASH=           GenBankHandle::CountAllPointFormRoot_to_leaf_for_many_SPECIES_for_many_Char    ($SonDad_HASH, $AllGeneBkHASH, $PtNbID_toChar_HASH); 
sub CountAllPointFormRoot_to_leaf_for_many_SPECIES_for_many_Char{
	my ($SonDadHASH, $in_AllGeneBkHASH, $In_PtNbID_toChar_HASH)=@_;
  
  my $warnMsgBody="\nIn package  GenBankHandle,\tIn sub CountAllPointFormRoot_to_leaf_for_many_SPECIES_for_many_Char,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $callerMSG=DirFileHandle::print_SubCallerInform;
  
  my $outPutHASH;
  if (   (  defined ( $In_PtNbID_toChar_HASH )  ) && (   ref ( $In_PtNbID_toChar_HASH ) eq 'HASH' )   ){
  	my $char_type_hash;
  	foreach my $eachPtNbID (    sort { $a cmp $b } (   keys (  %{ $In_PtNbID_toChar_HASH }  )   )    ){
  		$char_type_hash->{ $In_PtNbID_toChar_HASH->{$eachPtNbID} }=1;
  	}
  	
  	if (   (  defined ( $in_AllGeneBkHASH )  ) && (   ref ( $in_AllGeneBkHASH ) eq 'HASH' ) && (  defined ( $char_type_hash )  ) && (   ref ( $char_type_hash ) eq 'HASH' )   ){
  	  foreach my $eachChar (    sort { $a cmp $b } (   keys (  %{ $char_type_hash }  )   )    ){  
  		
  		
  			
  			my  $CountHASH;
      	foreach my $eachProID (    sort { $a cmp $b } (   keys (  %{ $in_AllGeneBkHASH }  )   )    ){  
      		
      		
      		
      		if  (      (  defined ( $in_AllGeneBkHASH->{$eachProID}                   )  )
      		        && (  defined ( $in_AllGeneBkHASH->{$eachProID}->{'2_0__orgaism'} )  )
      		        && ( $in_AllGeneBkHASH->{$eachProID}->{'2_0__orgaism'}=~m/\S+/       )
      		        && (  defined ( $in_AllGeneBkHASH->{$eachProID}->{'a_0_InfmlIdx'} )  )
      		        && ( $in_AllGeneBkHASH->{$eachProID}->{'a_0_InfmlIdx'}=~m/\d+/       )   
      		    )
      		{
      		  my $eachSpece=$in_AllGeneBkHASH->{$eachProID}->{'2_0__orgaism'};
      		  my $eachSequn=$in_AllGeneBkHASH->{$eachProID}->{'0_0_sequence'};   
      		  my $eachPtINb=$in_AllGeneBkHASH->{$eachProID}->{'a_0_InfmlIdx'};   
      		  my $theCharHr;  #&PtKLook( 1, "\$theCharHr=$theCharHr" );
      		  if (   (  defined ( $In_PtNbID_toChar_HASH->{$eachPtINb} )  ) && ( $In_PtNbID_toChar_HASH->{$eachPtINb}=~m/\S+/ )   ){
      		  	$theCharHr=$In_PtNbID_toChar_HASH->{$eachPtINb}; &PtKLook( 2, "\$theCharHr=$theCharHr=\$In_PtNbID_toChar_HASH->{$eachPtINb}" );
      		  } 
      		  else {
      		  	my $In_PtNbID_toChar_HASH_dump=DirFileHandle::ReturnDumperInform ( $In_PtNbID_toChar_HASH->{$eachPtINb} );
      		  	DieWork::Just_dieWork( $die_MsgHead."\n\$In_PtNbID_toChar_HASH->{$eachPtINb}=$In_PtNbID_toChar_HASH->{$eachPtINb} is not right!!! \n\$In_PtNbID_toChar_HASH_dump=\n$In_PtNbID_toChar_HASH_dump\n\n$! \n\n\n" );
      		  }
      		 
      		  &PtKLook( 3, "\$eachChar=$eachChar" );
      		  if (   (  defined ( $eachChar )  ) && ( $eachChar=~m/\S+/ ) && (  defined ( $theCharHr )  ) && ( $theCharHr=~m/\S+/ ) &&  ($eachChar eq  $theCharHr)   ){
      		    if (   (   defined ( $CountHASH->{$eachSpece} )   ) && ( $CountHASH->{$eachSpece}=~m/\S+/ )   ){  
      		      $CountHASH->{$eachSpece}++; &PtKLook( 4, "\$CountHASH->{$eachSpece}=$CountHASH->{$eachSpece}" );
      		    }
      		    else {
      		    	$CountHASH->{$eachSpece}=1; &PtKLook( 5, "\$CountHASH->{$eachSpece}=$CountHASH->{$eachSpece}" );
      		    }    
      		    my $EachCountHASH; $EachCountHASH=GenBankHandle::CountAllPointFromRoot_to_leaf($SonDadHASH, $eachSpece, $EachCountHASH);
      	      $CountHASH=GenBankHandle::Join_2_hash__sameKeyOF1_willBE_added_by_HASH2_2_LEVEL($CountHASH, $EachCountHASH);
      		  }
      	    
      		
      		}
      		
      	}
      	
      	$outPutHASH->{$eachChar}=Storable::dclone( $CountHASH ) if (   (  defined ( $CountHASH )  ) && (   ref ( $CountHASH ) eq 'HASH' )   );
      	DirFileHandle::PrintAndWarnDumper ($outPutHASH->{$eachChar}, "\n20181230-2-0 \$outPutHASH->{$eachChar}= \n");
      }
  		
  		
  		
  	}
  }
  
  
    
  return $outPutHASH;

}

sub PtKLook{	my ($a, $b)=@_;	print "\n20181230-2-".$a." ".$b."\n";}

sub CountAllPointFormRoot_to_leaf_for_many_SPECIES_for_U{
	my ($SonDadHASH, $in_AllGeneBkHASH, $CountHASH)=@_;
  
  my $warnMsgBody="\nIn package  GenBankHandle,\tIn sub CountAllPointFormRoot_to_leaf_for_many_SPECIES_for_U,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $callerMSG=DirFileHandle::print_SubCallerInform;
  
  if (   (  defined ( $in_AllGeneBkHASH )  ) && (   ref ( $in_AllGeneBkHASH ) eq 'HASH' )   ){
  	foreach my $eachProID (    sort { $a cmp $b } (   keys (  %{ $in_AllGeneBkHASH }  )   )    ){  
  		
  		if  (      (  defined ( $in_AllGeneBkHASH->{$eachProID}                   )  )
  		        && (  defined ( $in_AllGeneBkHASH->{$eachProID}->{'2_0__orgaism'} )  )
  		        && ( $in_AllGeneBkHASH->{$eachProID}->{'2_0__orgaism'}=~m/\S+/       )   
  		    )
  		{
  		  my $eachSpece=$in_AllGeneBkHASH->{$eachProID}->{'2_0__orgaism'};
  		  my $eachSequn=$in_AllGeneBkHASH->{$eachProID}->{'0_0_sequence'};   
  		  my $U_pos_array=FastaFileHandle::Find_U_pos_inAminoAcid( $eachSequn ); 
    		if (   (  defined ( $U_pos_array )  ) && (  ref( $U_pos_array ) eq 'ARRAY'  )   ){
  		    if (   (   defined ( $CountHASH->{$eachSpece} )   ) && ( $CountHASH->{$eachSpece}=~m/\S+/ )   ){
  		      $CountHASH->{$eachSpece}++;
  		    }
  		    else {
  		    	$CountHASH->{$eachSpece}=1;
  		    }    
  		    my $EachCountHASH; $EachCountHASH=GenBankHandle::CountAllPointFromRoot_to_leaf($SonDadHASH, $eachSpece, $EachCountHASH);
  	      $CountHASH=GenBankHandle::Join_2_hash__sameKeyOF1_willBE_added_by_HASH2_2_LEVEL($CountHASH, $EachCountHASH);
  		  }
  		  
  	    
  		
  		}
  		
  	}
  }
  return $CountHASH;

}

sub CountAllPointFormRoot_to_leaf_for_many_SPECIES{
	my ($SonDadHASH, $in_AllGeneBkHASH, $CountHASH)=@_;
  
  my $warnMsgBody="\nIn package  GenBankHandle,\tIn sub CountAllPointFormRoot_to_leaf_for_many_SPECIES,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $callerMSG=DirFileHandle::print_SubCallerInform;
  
  if (   (  defined ( $in_AllGeneBkHASH )  ) && (   ref ( $in_AllGeneBkHASH ) eq 'HASH' )   ){
  	foreach my $eachProID (    sort { $a cmp $b } (   keys (  %{ $in_AllGeneBkHASH }  )   )    ){  
  		
  		if  (      (  defined ( $in_AllGeneBkHASH->{$eachProID}                   )  )
  		        && (  defined ( $in_AllGeneBkHASH->{$eachProID}->{'2_0__orgaism'} )  )
  		        && ( $in_AllGeneBkHASH->{$eachProID}->{'2_0__orgaism'}=~m/\S+/       )   
  		    )
  		{
  		  my $eachSpece=$in_AllGeneBkHASH->{$eachProID}->{'2_0__orgaism'};
  		  if (   (   defined ( $CountHASH->{$eachSpece} )   ) && ( $CountHASH->{$eachSpece}=~m/\S+/ )   ){
  		    $CountHASH->{$eachSpece}++;
  		  }
  		  else {
  		  	$CountHASH->{$eachSpece}=1;
  		  }
  		  
  	    my $EachCountHASH; $EachCountHASH=GenBankHandle::CountAllPointFromRoot_to_leaf($SonDadHASH, $eachSpece, $EachCountHASH);
  	    $CountHASH=GenBankHandle::Join_2_hash__sameKeyOF1_willBE_added_by_HASH2_2_LEVEL($CountHASH, $EachCountHASH);
  		
  		}
  		
  	}
  }
  return $CountHASH;

}

sub CountAllPointFromRoot_to_leaf{  #    GenBankHandle::CountAllPointFromRoot_to_leaf
  my ($SonDadHASH, $InSpece, $CountHASH, $levNb)=@_;
  
  my $warnMsgBody="\nIn package  GenBankHandle,\tIn sub CountAllPointFromRoot_to_leaf,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $callerMSG=DirFileHandle::print_SubCallerInform;
	
  if  (   (   defined ( $levNb )   ) && ( $levNb=~m/\d+/ ) && ( $levNb>=0 )   ){  	
  }
  else {
  	$levNb=0;
  }
  if (   (   defined ( $SonDadHASH->{$InSpece} )   ) && ( $SonDadHASH->{$InSpece}=~m/\S+/ )   ){
    my $newSpec=$SonDadHASH->{$InSpece};
    if (   (   defined ( $CountHASH->{$newSpec} )   ) && ( $CountHASH->{$newSpec}=~m/\d+/ )   ){
    	$CountHASH->{$newSpec}++;    print "20181127-1 \$levNb=$levNb \$CountHASH->{$newSpec}=$CountHASH->{$newSpec}\n";
    }
    else {
    	$CountHASH->{$newSpec}=1;    print "20181127-2 \$levNb=$levNb \$CountHASH->{$newSpec}=$CountHASH->{$newSpec}\n";
    }
    
    #if ( $levNb==0) { 
    #	if (   (   defined ( $CountHASH->{$InSpece} )   ) && ( $CountHASH->{$InSpece}=~m/\d+/ )   ){
    #	  $CountHASH->{$InSpece}++;  print "20181127-01 \$CountHASH->{$InSpece}=$CountHASH->{$InSpece}\n";
    #  }
    #  else {
    #  	$CountHASH->{$InSpece}=1;  print "20181127-02 \$CountHASH->{$InSpece}=$CountHASH->{$InSpece}\n";
    #  }
    #}
    $levNb++;
    $CountHASH=CountAllPointFromRoot_to_leaf($SonDadHASH, $newSpec, $CountHASH, $levNb);   
    
  }
  else {    
  }
  return $CountHASH;
}



sub Find_repeat_and_doing_statistic{  #GenBankHandle::Find_repeat_and_doing_statistic
	my ($inDIR, $in_top_howMany)=@_;
	
	my $warnMsgBody="\nIn package  GenBankHandle,\tIn sub Find_repeat_and_doing_statistic,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $callerMSG=DirFileHandle::print_SubCallerInform;
	
	
	my $fileARRAY=DirFileHandle::getDirArray($inDIR);
	my $bigGBKhash;
	my $repeatHASH;
	my $PtFindHash;
	foreach my $ecFile (  @{ $fileARRAY }  ){ 
		my $realFilePath=$inDIR."/".$ecFile; warn "\$realFilePath=$realFilePath\n\n\n";
		my $GeBkHASH=GenBankHandle::ReadGeneBankFILE($realFilePath, '', $in_top_howMany); #DirFileHandle::PrintDumper($ecFile."hash",$GeBkHASH);
		if ( ref ($GeBkHASH) eq 'HASH' ){
			foreach my $eachPRPid (    sort { $a cmp $b } (   keys (  %{ $GeBkHASH }  )   )    ){
				$bigGBKhash->{$eachPRPid}=$GeBkHASH->{$eachPRPid};
				$repeatHASH->{$eachPRPid}->{'0_0_2_repeatFileHASH'}->{$ecFile}=1;
				$PtFindHash->{$ecFile}->{$eachPRPid}=1;
			}
		}
	}
	my $realRptHASH;
	if ( ref ($repeatHASH) eq 'HASH' ){
		foreach my $eachPRPid (    sort { $a cmp $b } (   keys (  %{ $repeatHASH }  )   )    ){
			if (   (  defined ( $repeatHASH->{$eachPRPid}->{'0_0_2_repeatFileHASH'} )  ) && (  ref ( $repeatHASH->{$eachPRPid}->{'0_0_2_repeatFileHASH'} ) eq 'HASH' )   ){
				my $repeatNB=(   keys (  %{ $repeatHASH->{$eachPRPid}->{'0_0_2_repeatFileHASH'} }  )   );
				if ( $repeatNB>1){
					$realRptHASH->{$eachPRPid}->{'0_0_1_repeat__number'}=$repeatNB;
					$realRptHASH->{$eachPRPid}->{'0_0_2_repeatFileHASH'}=$repeatHASH->{$eachPRPid}->{'0_0_2_repeatFileHASH'};
				}
			}
		}
	}
	return [$bigGBKhash, $realRptHASH, $PtFindHash]
}


sub BuildAll_StfcName_HASH{  #   my $outHASH=GenBankHandle::BuildAll_StfcName_HASH($InHASH);
	my ($InHASH)=@_;
	
	my $warnMsgBody="\nIn package  GenBankHandle,\tIn sub BuildAll_StfcName_HASH,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  DirFileHandle::print_SubCallerInform;
  
  my $outHASH;
  if (   (  defined ( $InHASH )  ) && (  ref ( $InHASH ) eq 'HASH' )   ){
  	foreach my $eachPRPid (    sort { $a cmp $b } (   keys (  %{ $InHASH }  )   )    ){  print "\n201812251638-1 \$eachPRPid=$eachPRPid\n";
  		if  (   (  defined ( $InHASH->{$eachPRPid}->{'2_0__orgaism'} )  ) && ( $InHASH->{$eachPRPid}->{'2_0__orgaism'}=~m/\S+/ )   ){ print "\n201812251638-2 \$InHASH->{$eachPRPid}->{'2_0__orgaism'}=$InHASH->{$eachPRPid}->{'2_0__orgaism'}\n";
  		  if (    (   defined (  $outHASH->{ $InHASH->{$eachPRPid}->{'2_0__orgaism'} }  )   )&& ( $outHASH->{ $InHASH->{$eachPRPid}->{'2_0__orgaism'} }=~m/\d+/ )    ){
  		    $outHASH->{ $InHASH->{$eachPRPid}->{'2_0__orgaism'} }++;   print "\n201812251638-3 \$outHASH->{ \$InHASH->{$eachPRPid}->{'2_0__orgaism'} }=\$outHASH->{ $InHASH->{$eachPRPid}->{'2_0__orgaism'} }=$outHASH->{ $InHASH->{$eachPRPid}->{'2_0__orgaism'} }\n";
  		  }
  		  else{
  		  	$outHASH->{ $InHASH->{$eachPRPid}->{'2_0__orgaism'} }=1;   print "\n201812251638-4 \$outHASH->{ \$InHASH->{$eachPRPid}->{'2_0__orgaism'} }=\$outHASH->{ $InHASH->{$eachPRPid}->{'2_0__orgaism'} }=$outHASH->{ $InHASH->{$eachPRPid}->{'2_0__orgaism'} }\n";
  		  }  		  
  		}
  	}
  }
  return $outHASH;
}


sub BuildAllSpecise_HASH{  #GenBankHandle::BuildAllSpecise_HASH
	my ($InHASH)=@_;
	
	my $warnMsgBody="\nIn package  GenBankHandle,\tIn sub BuildAllSpecise_HASH,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  DirFileHandle::print_SubCallerInform;
  
  my $outHASH;
  if (   (  defined ( $InHASH )  ) && (  ref ( $InHASH ) eq 'HASH' )   ){
  	foreach my $eachPRPid (    sort { $a cmp $b } (   keys (  %{ $InHASH }  )   )    ){
  		if  (   (  defined ( $InHASH->{$eachPRPid}->{'4_0_specComN'} )  ) && ( $InHASH->{$eachPRPid}->{'4_0_specComN'}=~m/\S+/ )   ){
  		  if (    (   defined (  $outHASH->{ $InHASH->{$eachPRPid}->{'4_0_specComN'} }  )   )&& ( $outHASH->{ $InHASH->{$eachPRPid}->{'4_0_specComN'} }=~m/\d+/ )    ){
  		    $outHASH->{ $InHASH->{$eachPRPid}->{'2_0__orgaism'} }++;
  		  }
  		  else{
  		  	$outHASH->{ $InHASH->{$eachPRPid}->{'2_0__orgaism'} }=1;
  		  }  		  
  		}
  	}
  }
  return $outHASH;
}

sub BuildAllTxnomy_phy_HASH{  #GenBankHandle::BuildAllTxnomy_phy_HASH
	my ($InHASH)=@_;
	
	my $warnMsgBody="\nIn package  GenBankHandle,\tIn sub BuildAllTxnomy_phy_HASH,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  DirFileHandle::print_SubCallerInform;
  
  my $outHASH;
  if (   (  defined ( $InHASH )  ) && (  ref ( $InHASH ) eq 'HASH' )   ){
  	foreach my $eachPRPid (    sort { $a cmp $b } (   keys (  %{ $InHASH }  )   )    ){
  		if  (   (  defined ( $InHASH->{$eachPRPid}->{'5_0_phyClsif'} )  ) && (  ref ( $InHASH->{$eachPRPid}->{'5_0_phyClsif'} ) eq 'ARRAY' )   ){
  		  foreach my $eachPHY (  @{ $InHASH->{$eachPRPid}->{'5_0_phyClsif'} }  ){
  		  	$outHASH->{$eachPHY}=1;
  		  }
  		}
  	}
  }
  return $outHASH;
}

sub BuildAll_Dad_son_relationShip_HASH{  #GenBankHandle::BuildAll_Dad_son_relationShip_HASH
	my ($InHASH)=@_;
	
	my $warnMsgBody="\nIn package  GenBankHandle,\tIn sub BuildAll_Dad_son_relationShip_HASH,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  DirFileHandle::print_SubCallerInform;
  
  my $outHASH;
  if (   (  defined ( $InHASH )  ) && (  ref ( $InHASH ) eq 'HASH' )   ){
  	foreach my $eachPRPid (    sort { $a cmp $b } (   keys (  %{ $InHASH }  )   )    ){
  		if  (   (  defined ( $InHASH->{$eachPRPid}->{'5_0_phyClsif'} )  ) && (  ref ( $InHASH->{$eachPRPid}->{'5_0_phyClsif'} ) eq 'ARRAY' )   ){
  		  my $lastPHY; my $stepNB=0;
  		  foreach my $eachPHY (  @{ $InHASH->{$eachPRPid}->{'5_0_phyClsif'} }  ){
  		  	if ( $stepNB >0 ){
  		  		$outHASH->{$lastPHY}->{$eachPHY}=1;
  		  	}
  		  	$lastPHY=$eachPHY;
  		  	$stepNB++;
  		  }
  		}
  	}
  }
  return $outHASH;
}

sub BuildAll_son_Dad_relationShip_HASH{  #GenBankHandle::BuildAll_son_Dad_relationShip_HASH
	my ($InHASH)=@_;
	
	my $warnMsgBody="\nIn package  GenBankHandle,\tIn sub BuildAll_son_Dad_relationShip_HASH,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  DirFileHandle::print_SubCallerInform;
  
  my $outHASH;
  if (   (  defined ( $InHASH )  ) && (  ref ( $InHASH ) eq 'HASH' )   ){
  	foreach my $eachPRPid (    sort { $a cmp $b } (   keys (  %{ $InHASH }  )   )    ){
  		if  (   (  defined ( $InHASH->{$eachPRPid}->{'5_0_phyClsif'} )  ) && (  ref ( $InHASH->{$eachPRPid}->{'5_0_phyClsif'} ) eq 'ARRAY' )   ){
  		  my $lastPHY; my $stepNB=0;
  		  foreach my $eachPHY (  @{ $InHASH->{$eachPRPid}->{'5_0_phyClsif'} }  ){
  		  	if ( $stepNB >0 ){
  		  		$outHASH->{$eachPHY}=$lastPHY;
  		  	}
  		  	$lastPHY=$eachPHY;
  		  	$stepNB++;
  		  }
  		}
  	}
  }
  return $outHASH;
}

sub Join_2_hash__sameKeyOF1_willBE_overwrite_by_HASH2_2_LEVEL{ #GenBankHandle::Join_2_hash__sameKeyOF1_willBE_overwrite_by_HASH2_2_LEVEL
	my ($InHASH_1, $InHASH_2)=@_;
	
	my $warnMsgBody="\nIn package  GenBankHandle,\tIn sub Join_2_hash__sameKeyOF1_willBE_overwrite_by_HASH2,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  DirFileHandle::print_SubCallerInform;
  
  
  if (   (  defined ( $InHASH_2 )  ) && (  ref ( $InHASH_2 ) eq 'HASH' )   ){
  	foreach my $eachKey_lev1 (    sort { $a cmp $b } (   keys (  %{ $InHASH_2 }  )   )    ){
  		if  (   (  defined ( $InHASH_2->{$eachKey_lev1} )  ) && (  ref ( $InHASH_2->{$eachKey_lev1} ) eq 'HASH' )   ){  		  
  		  foreach my $eachKey_lev2 (    sort { $a cmp $b } (   keys (  %{ $InHASH_2->{$eachKey_lev1} }  )   )    ){
  		  	$InHASH_1->{$eachKey_lev1}->{$eachKey_lev2}=$InHASH_2->{$eachKey_lev1}->{$eachKey_lev2};
  		  }
  		}
  		else {
  			$InHASH_1->{$eachKey_lev1}=$InHASH_2->{$eachKey_lev1};
  		}
  	}
  }
  return $InHASH_1;
}

sub Join_2_hash__sameKeyOF1_willBE_added_by_HASH2_2_LEVEL{ #GenBankHandle::Join_2_hash__sameKeyOF1_willBE_added_by_HASH2_2_LEVEL
	my ($InHASH_1, $InHASH_2)=@_;
	
	my $warnMsgBody="\nIn package  GenBankHandle,\tIn sub Join_2_hash__sameKeyOF1_willBE_added_by_HASH2_2_LEVEL,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  DirFileHandle::print_SubCallerInform;
  
  
  if (   (  defined ( $InHASH_2 )  ) && (  ref ( $InHASH_2 ) eq 'HASH' )   ){
  	foreach my $eachKey_lev1 (    sort { $a cmp $b } (   keys (  %{ $InHASH_2 }  )   )    ){
  		if  (   (  defined ( $InHASH_2->{$eachKey_lev1} )  ) && (  ref ( $InHASH_2->{$eachKey_lev1} ) eq 'HASH' )   ){  		  
  		  foreach my $eachKey_lev2 (    sort { $a cmp $b } (   keys (  %{ $InHASH_2->{$eachKey_lev1} }  )   )    ){
  		  	$InHASH_1->{$eachKey_lev1}->{$eachKey_lev2}+=$InHASH_2->{$eachKey_lev1}->{$eachKey_lev2};
  		  }
  		}
  		else {
  			$InHASH_1->{$eachKey_lev1}+=$InHASH_2->{$eachKey_lev1};
  		}
  	}
  }
  return $InHASH_1;
}




1;