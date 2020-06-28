#!/usr/bin/perl -w
BEGIN{  push (@INC, '/home/fredjiang/fredProgs/PMs')  };

use strict;
use warnings;
use FastaFileHandle;
use DirFileHandle;
use DieWork;
use DiamondWork;
use BlastPraser;
use BlastHandle;
use Storable ;
use Storable qw(nstore store_fd nstore_fd freeze thaw dclone);
use ClustalwRun;
use ExcelHandle;
use File::Basename;
use Interproscan;
use InFileHandle;
use Interproscan_MultiplAnalysis;
use CDhitPraser;
use File::Spec;        
                      
package SelenoProteinData;


##################################################################################################
#
#

#working steps:
#1st    #  my $outHash_holding_blstInform= SelenoProteinData::Blast_and_CDhit_to_getInformation_forSelDB($in_old_DataBaseFastaFILE, $blastAnaysis_OF_oldSELENODB_DIR);
        #  DirFileHandle::PrintDumper_plReadFile($outHash_holding_blstInform_FILE, $outHash_holding_blstInform)
        
#2nd    #  hand make working of family identification
#open and copy all information of 1m_BlsChDIR/0_6_0_SelDB__Table.txt into 1m_BlsChDIR/0_4_0_CDhitDIR/4m_CdHit_WKingDIR010.ips.xlsx
#do the hand checking and family identification working in the 4m_CdHit_WKingDIR010.ips.xlsx, usually it was copied into the windows working computer
#then it was recopied into the HPC, and name-changed into 2m2m_CdHit_WKingDIR010.ips.handWorking.xlsx


#copy the txt infomation of   2m2m_CdHit_WKingDIR010.ips.handWorking.xlsx       and pasted into 3m1m_handMade.divFamily.excel.txt
#then do the recheck working below

#3rd    re divided of selenoprotein families
#my $inExcelTXTfile="old_xls_informaltil_add_group.excel.txt";
#SelenoProteinData::With_manual_HANDLE_excel_build_family_HASH( $hand_made_familyDivided_excel_1, 60, $hand_made_familyDivided_DIR_1, 46, 1, 'SF'  );

#4th    re divided of selenoprotein families
#my $inExcelTXTfile="old_xls_needtoRediv.txt";
#SelenoProteinData::With_manual_HANDLE_excel_build_family_HASH( $hand_made_familyDivided_excel_2, 60, $hand_made_familyDivided_DIR_2, 47, 1, 'BG'  ); #BG means bigGroup
   
#
#
##################################################################################################

#The first original Selenoprotein database file, all the id was ordered with JL_number
my $orginal_JL_SelproDATABASE="/home/fredjiang/BlastDB/SelDB_20121109/SelDB_20121107/SelDB_1351570402/org_selenoproteinDB";

#excel file including different information of those  first original Selenoprotein data, all the id was ordered with JL_number
my $fake_excel_selProInform_1="/home/fredjiang/BlastDB/SelDB_20121109/SelDB_20121107/SelDB_1351570402/selenoproteinDB_inform.xls"; #text file

#excel file including some selenoproteins from ZhangYan
my $fake_excel_selProInform_2="/home/fredjiang/BlastDB/SelDB_20121109/SelDB_20121107/SelDB_1352294452/SelDB_fromZY.txt_inform.xls"; #text file

#the Hash holding informations from the 2 excel files above
my $big_raw_inform_Hash_file='/home/fredjiang/OneTSSD/fredjiang20190321/SELENODB_20191217/0m_old_xls_inform.table.hsh';



#   #the last new selenoprotein DB, in which something wrong with U were delete, the cis-elements were also deleted and added the selenoproteins predicted from my work----metazone and algae selenoproteins, and added the archare Loki selenorpteins.
#    my $Last_new__20170322__SelDB="/home/fredjiang/BlastDB/SelDB_2017032220170322__SelDB";

# the selenoprotein database file to extrac the fasta sequences and do the addtional anaysis work.
my $in_old_DataBaseFastaFILE      ='/home/fredjiang/BlastDB/SelDB_20170119/20170119__SelDB';
# the directory to hold the anlysis result
my $blastAnaysis_OF_oldSELENODB_DIR='/home/fredjiang/OneTSSD/fredjiang20190321/SELENODB_20191217/1m_BlsChDIR';
# the hash dumpfile to hold the anlysis result
my $outHash_holding_blstInform_FILE='/home/fredjiang/OneTSSD/fredjiang20190321/SELENODB_20191217/2m_BlsCdHit.hsh';


#The #
my $hand_made_familyDivided_excel_1=   '/home/fredjiang/OneTSSD/fredjiang20190321/SELENODB_20191217/3m1m_handMade.divFamily.excel.txt';
my $hand_made_familyDivided_DIR_1=     '/home/fredjiang/OneTSSD/fredjiang20190321/SELENODB_20191217/3m1m_DIR';
my $hand_made_familyDivided_excel_2=   '/home/fredjiang/OneTSSD/fredjiang20190321/SELENODB_20191217/3m2m_handMade.divSupFml.excel.txt';
my $hand_made_familyDivided_DIR_2=     '/home/fredjiang/OneTSSD/fredjiang20190321/SELENODB_20191217/3m2m_DIR';



#          my $SelenoProteinHASH;
#          my $SELproIDkey__HASH;
#          
#          
#          
#          my $EachSELpro___HASH;
#          
#          $SELproIDkey__HASH->{'JLDB_xxxxxxxx'    }=$EachSELpro___HASH;
#          
#          
#          
#          
#          
#          my $OutDbId_HASH_KEY='';   #Íâ²¿idÐÎ³ÉÒ»¸öhash,keyÊÇid£¬Öµ¿ÉÒÔÊÇ1£¬Ò²¿ÉÒÔÊÇÕâ¸öidµÄ¸½¼ÓÐÅÏ¢
#          
#          
#          my $species_name_KEY='';   #ÎïÖÖÃû
#          my $taxonomy_ID__KEY='';   #ÎïÖÖµÄtaxnomy ID
#          my $StrainIDorNamKEY='';   #¾úÖêÐÅÏ¢
#          
#          my $AbsNoChngePepKEY='';   #Ã¿¸öÎøµ°°×µÄ È·¶¨µÄ Î¨Ò»µÄ °±»ùËáÐòÁÐ£¨¿ÉÒÔºÍºóÃæµÄ»ùÒò×é»òEstÖÐ»ñµÃµÄPEPÐòÁÐÏàÍ¬»òÓÐ²îÒì£©
#          
#          my $dataResousRefKEY='';   #¸ÃÎøµ°°×µÄÎÄÏ×À´Ô´
#          
#          my $GenomicID____KEY='';   #»ùÒò×é fastaÎÄ¼þµÄid
#          my $GM_CdsPosListKEY='';   #´Ó»ùÒò×éÖÐ »ñµÃµÄ ¸ÃÎøµ°°×»ùÒòµÄ CDSµÄ¸÷¸ö ×é³ÉµÄÎ»ÖÃ£¬ÀàÐÍµÈÐÅÏ¢
#          my $GM_SegPosHsh_KEY='';   #´Ó»ùÒò×éÖÐ »ñµÃµÄ ¸ÃÎøµ°°×»ùÒòµÄ ¸÷¸ö ×é³ÉµÄÎ»ÖÃ£¬ÀàÐÍµÈÐÅÏ¢hash
#          my $GM_CdsDNAseq_KEY='';   #´Ó»ùÒò×éÖÐ »ñµÃµÄ ¸ÃÎøµ°°×»ùÒòµÄ CDS µÄDNAÐòÁÐ
#          my $GM_CdsPEPseq_KEY='';   #´Ó»ùÒò×éÖÐ »ñµÃµÄ ¸ÃÎøµ°°×»ùÒòµÄ CDS µÄPEPÐòÁÐ
#          my $Gm_SECIS_Ary_KEY='';   #´Ó»ùÒò×éÖÐ »ñµÃµÄ ¸ÃÎøµ°°×»ùÒòµÄ SECIS µÄÊý×éµÄÐÅÏ¢£¨¿ÉÄÜÓÐ¶à¸öSECIS£©
#          
#          my $ESTcontigID__KEY='';   #EstÆ´½ÓµÄContig fastaÎÄ¼þµÄid
#          my $EstCdsPosListKEY='';   #´ÓEstÆ´½ÓContigÖÐ »ñµÃµÄ ¸ÃÎøµ°°×»ùÒòµÄ CDSµÄ¸÷¸ö ×é³ÉµÄÎ»ÖÃ£¬ÀàÐÍµÈÐÅÏ¢
#          my $EstSegPosHsh_KEY='';   #´ÓEstÆ´½ÓContigÖÐ »ñµÃµÄ ¸ÃÎøµ°°×»ùÒòµÄ ¸÷¸ö ×é³ÉµÄÎ»ÖÃ£¬ÀàÐÍµÈÐÅÏ¢hash
#          my $EstCdsDNAseq_KEY='';   #´ÓEstÆ´½ÓContigÖÐ »ñµÃµÄ ¸ÃÎøµ°°×»ùÒòµÄ CDS µÄDNAÐòÁÐ
#          my $EstCdsPEPseq_KEY='';   #´ÓEstÆ´½ÓContigÖÐ »ñµÃµÄ ¸ÃÎøµ°°×»ùÒòµÄ CDS µÄPEPÐòÁÐ
#          my $EstSECIS_Ary_KEY='';   #´ÓEstÆ´½ÓContigÖÐ »ñµÃµÄ ¸ÃÎøµ°°×»ùÒòµÄ SECIS µÄÊý×éµÄÐÅÏ¢£¨¿ÉÄÜÓÐ¶à¸öSECIS£©
#          
#          my $Est_ID__HASH_KEY='';   # EstÆ´½ÓContig,Õâ¸öhash¼ÇÂ¼ÁËÆäÊÇÓÉÄÄÐ©ESTÆ´½Ó¶øÀ´£¬ÆäIDÎªkey£¬ÖµÎª¸üÏêÏ¸ÐÅÏ¢
#          
#          
#          #$EachSELpro___HASH->{'0_0_0_'};
#          
#          
#          $SelenoProteinHASH->{'0_0_0_ID_key_HASH'}=$SELproIDkey__HASH;



sub Step_1 {    # SelenoProteinData::Step_1();
	
	my $outHash_holding_blstInform= SelenoProteinData::Blast_and_CDhit_to_getInformation_forSelDB($in_old_DataBaseFastaFILE, $blastAnaysis_OF_oldSELENODB_DIR);
  DirFileHandle::PrintDumper_plReadFile($outHash_holding_blstInform_FILE, $outHash_holding_blstInform)

}        

sub Step_2 {    # SelenoProteinData::Step_2();
	
	SelenoProteinData::With_manual_HANDLE_excel_build_family_HASH( $hand_made_familyDivided_excel_1, 60, $hand_made_familyDivided_DIR_1, 46, 1, 'SF'  );

}   

sub Step_3 {    # SelenoProteinData::Step_3();
	
	SelenoProteinData::With_manual_HANDLE_excel_build_family_HASH( $hand_made_familyDivided_excel_2, 60, $hand_made_familyDivided_DIR_2, 47, 1, 'BG'  ); #BG means bigGroup

}   







sub Blast_and_CDhit_to_getInformation_forSelDB{    # my $outHash= SelenoProteinData::Blast_and_CDhit_to_getInformation_forSelDB($oldDataBasePath, $workingDIRpath);
	my ($oldDataBasePath, $workingDIRpath)=@_;
	
	my ($warnMsgBody, $warnMsgHead, $die_MsgHead, $caller_inform)=@{ DieWork::BuildWarnDieHeadInform( 'SelenoProteinData', 'Blast_and_CDhit_to_getInformation_forSelDB' ) };
  
  DieWork::Check_FileDirExist_or_DIE   ( $oldDataBasePath,     "\$oldDataBasePath",    $die_MsgHead, $caller_inform  );
  DieWork::Check_DfdNoEmptString_or_DIE( $workingDIRpath,      "\$workingDIRpath",     $die_MsgHead, $caller_inform  );  
	
	my $mkDir_wkingDirPath="mkdir -p $workingDIRpath";  	DieWork::Print_and_warn( "\n$warnMsgBody$mkDir_wkingDirPath\n\n");
	system ( "$mkDir_wkingDirPath");
	
	
	my $oldDataFastaHash=FastaFileHandle::BuildHashFromFastaFile_Name_as_Key_checkRepeat_WithSeq ( $oldDataBasePath ) ; 
	
	
	my $InputSelProFile_noSuffix=$workingDIRpath."/0_0_0_in_fasta";
	
	my $oldDataFastaHash_HshFile=$InputSelProFile_noSuffix."\.hsh";  DirFileHandle::PrintDumper_plReadFile($oldDataFastaHash_HshFile, $oldDataFastaHash) if ( ref($oldDataFastaHash) eq 'HASH' ) ;
	
	my $InputSelProFile_FastaTxt=$InputSelProFile_noSuffix."\.txt";
	my $InputSelProFile_IndexTxt=$InputSelProFile_noSuffix."\.idx";
	
	
	my $inputSelProFastaString;   
	
  my ($inHash)=$oldDataFastaHash;
  
  my $seqTotalNumber=0;
  if (  ( ref ($inHash) ) eq 'HASH' ){
    foreach my $keyLev_0 (    sort { $a cmp $b } (   keys (  %{ $inHash } )   )    ){   print "    \$keyLev_0=$keyLev_0\n";
    	my $valLev_0=$inHash->{$keyLev_0};                                                print "    \$valLev_0=$valLev_0\n";
    	my $refLev_0=ref ($valLev_0);                                                     print "    \$refLev_0=$refLev_0\n\n";  
    	
    	#if (   (  defined ( $segNameHASH )  ) && (   ref ( $segNameHASH ) eq 'HASH' )   ){
      if (   ( $refLev_0 eq 'HASH' ) && (  defined ( $valLev_0->{ '0_0_6_____sequence' } )  ) && (   $valLev_0->{ '0_0_6_____sequence' } =~ m/\S+/  )   ) {
      	
      	$inputSelProFastaString.=">$keyLev_0\n$valLev_0->{ '0_0_6_____sequence' }\n\n";
      	
      }
      $seqTotalNumber++;
    }
  }
        
	FastaFileHandle::BuildFastaFile_withFastaString ($InputSelProFile_FastaTxt, $inputSelProFastaString  );  #
	FastaFileHandle::BuildIdxFile_for_fastaFile     ($InputSelProFile_FastaTxt, $InputSelProFile_IndexTxt);  #
	
	
	
	
	my $outPutDiamondBlastp__Out=$workingDIRpath."/0_1_0_diamondRst.txt";
	
	my $blastRsltPrase_Hash_FILE=$workingDIRpath."/0_1_1_blstPraseOUT.hsh"; 
	
	my $blast_result_Prase__HASH;
	my $diamond_work_is_already_done =0;
	if (  DieWork::Check_FileDirExist_or_NOT( $blastRsltPrase_Hash_FILE )  ){
		$diamond_work_is_already_done = DiamondWork::Check_diamond_done_EasyVersion( $oldDataFastaHash, $blastRsltPrase_Hash_FILE );
	}	
	if ( $diamond_work_is_already_done == 0 ){
		DiamondWork::RunDiamond_vs_NRdb( $InputSelProFile_FastaTxt, $outPutDiamondBlastp__Out, 'blastp', 5, 0, 99 );
		$blast_result_Prase__HASH=BlastPraser::BioPerlBlastPraser20191021_justPrase ($outPutDiamondBlastp__Out, 'blastp', $InputSelProFile_IndexTxt,  $seqTotalNumber ) ;
	}
	else{
		$blast_result_Prase__HASH=Storable::retrieve( $blastRsltPrase_Hash_FILE );
	}
	
	
	 
	
	DirFileHandle::PrintDumper_plReadFile($blastRsltPrase_Hash_FILE, $blast_result_Prase__HASH) if ( ref($blast_result_Prase__HASH) eq 'HASH' ) ; 
  my $blstQryBestHit_Hash_FILE=$workingDIRpath."/0_1_2_blsQryBstHit.hsh";
  
  #my $blast_result_Prase__HASH=Storable::retrieve( $blastRsltPrase_Hash_FILE );
  
  my $All_Qry_Hit_Hash        =BlastPraser::New_ChangeBlastResult_into_Qry_Hit_HASH( $blast_result_Prase__HASH );  
  my $Total_Best_blastQueryHit_hash=BlastPraser::GetTheBestHit_for_eachQuery_basedON( $All_Qry_Hit_Hash, '3_0_6_hitEvalue' );
  DirFileHandle::PrintDumper_plReadFile($blstQryBestHit_Hash_FILE, $Total_Best_blastQueryHit_hash)      if ( ref ($Total_Best_blastQueryHit_hash) eq 'HASH' ) ;
  
  my $blstQryBestHit_tempGKB_FILE=$workingDIRpath."/0_1_3_blsQryBstHit.tempGbk.txt";
  my $blstQryBestHit_realGKB_FILE=$workingDIRpath."/0_1_4_blsQryBstHit.realGbk.txt";
  my $Total_Bshit_genank_InfomHASH_FILE=$workingDIRpath."/0_1_5_blsQryBstHit.gbk.ifm.hsh";
  
  my $All_Acc_to_Gi_HASH   =BlastPraser::Get_Acc_to_Gi_HASH_from_Qry_bestHit_HASH( $Total_Best_blastQueryHit_hash );
  my $Total_Bshit_genank_InfomHASH=Storable::retrieve( $Total_Bshit_genank_InfomHASH_FILE );  # delete
#my $Total_Bshit_genank_InfomHASH=GenBankHandle::Upgrade_db_and_geting_GB_inform_from_ProtHASH($All_Acc_to_Gi_HASH, $blstQryBestHit_tempGKB_FILE, $blstQryBestHit_realGKB_FILE);
#DirFileHandle::PrintDumper_plReadFile($Total_Bshit_genank_InfomHASH_FILE, $Total_Bshit_genank_InfomHASH)      if ( ref ($Total_Bshit_genank_InfomHASH) eq 'HASH' ) ;
  
  
  
  
  my $blastRslt_Cov1_Idt1_HASH=BlastHandle::Extract_subQueryHit_HASH( $blast_result_Prase__HASH, 1, 1 );
  my $blstRst_Cov1Idt1_HshFILE=$workingDIRpath."/0_2_0_0_blt_Cov1Idt1.hsh";
  DirFileHandle::PrintDumper_plReadFile($blstRst_Cov1Idt1_HshFILE, $blastRslt_Cov1_Idt1_HASH)      if ( ref ($blastRslt_Cov1_Idt1_HASH) eq 'HASH' ) ;
  
  
  my $blastRslt_Cov1_Idt1_QryHit_HASH=BlastPraser::New_ChangeBlastResult_into_Qry_Hit_HASH( $blastRslt_Cov1_Idt1_HASH );
  my $bstQrHt_Cov1Idt1_HshFILE=$workingDIRpath."/0_2_0_1_bstQrHitCId1.hsh";
  my $Cov1Id1_Best_blstQeryHit_hash=BlastPraser::GetTheBestHit_for_eachQuery_basedON( $blastRslt_Cov1_Idt1_QryHit_HASH, '3_0_6_hitEvalue' );
  DirFileHandle::PrintDumper_plReadFile($bstQrHt_Cov1Idt1_HshFILE, $Cov1Id1_Best_blstQeryHit_hash)      if ( ref ($Cov1Id1_Best_blstQeryHit_hash) eq 'HASH' ) ;
  
#  my $blastRslt_Cov1_Idt99HASH=BlastHandle::Extract_subQueryHit_HASH( $blast_result_Prase__HASH, 1, 0.99 );
#  my $blstRst_Cov1Idt99HshFILE=$workingDIRpath."/0_2_1_bltCov1Idt99.hsh";
#  DirFileHandle::PrintDumper_plReadFile($blstRst_Cov1Idt99HshFILE, $blastRslt_Cov1_Idt99HASH)      if ( ref ($blastRslt_Cov1_Idt99HASH) eq 'HASH' ) ;
  
  my $Temp_all_hit_gebank_File=$workingDIRpath."/0_3_0_temp_HitGBak.txt";
  my $Real_all_hit_gebank_File=$workingDIRpath."/0_3_1_real_HitGBak.txt";
  
                                                                                                  #0_2_0_0_blt_Cov1Idt1.hsh   0_3_0_temp_HitGBak.txt     0_3_1_real_HitGBak.txt
#my $Cov1ID1_hit_gebank_InfomHASH=SelenoProteinData::GetGenBank_Inform_from_BLASTpraseFile_HASH ($blastRslt_Cov1_Idt1_HASH, $Temp_all_hit_gebank_File, $Real_all_hit_gebank_File);
  my $AllHitGebankInfomHASHfil=$workingDIRpath."/0_3_2_real_HitGBak.hsh";
  my $Cov1ID1_hit_gebank_InfomHASH=Storable::retrieve( $AllHitGebankInfomHASHfil );  # delete
  DirFileHandle::PrintDumper_plReadFile($AllHitGebankInfomHASHfil, $Cov1ID1_hit_gebank_InfomHASH)      if ( ref ($Cov1ID1_hit_gebank_InfomHASH) eq 'HASH' ) ;
  
  
  
  
  
  my $CDhit_working_DIR=$workingDIRpath."/0_4_0_CDhitDIR";  
  my ($Cd_hit_SelPID_to_ClusterID_HASH_FILE, $Cd_hit_ClusterID_to_Aln_HASH_FILE )=@{ CDhitPraser::RunCdhit_and_Prase  ($InputSelProFile_FastaTxt, $InputSelProFile_IndexTxt, $oldDataFastaHash_HshFile, $CDhit_working_DIR)  };
  my $cluster_sepKeyHASH    =Storable::retrieve( $Cd_hit_SelPID_to_ClusterID_HASH_FILE );  
  my $cluster_cluID_AlnHASH =Storable::retrieve( $Cd_hit_ClusterID_to_Aln_HASH_FILE    );     
  
  
  my $big_raw_inform_Hash=SelenoProteinData::Build_the_raw_informHASH_from_oldEXCELfile( $fake_excel_selProInform_1, $fake_excel_selProInform_2, $big_raw_inform_Hash_file );
  my $SelProteinHash=SelenoProteinData::Build_SelenoProteinHASH( $oldDataFastaHash, $big_raw_inform_Hash, $Total_Bshit_genank_InfomHASH, $Cov1ID1_hit_gebank_InfomHASH, $Total_Best_blastQueryHit_hash, $Cov1Id1_Best_blstQeryHit_hash, $cluster_sepKeyHASH, $cluster_cluID_AlnHASH);
  
  
  my $SelProteinHash_FILE=$workingDIRpath."/0_5_0_out____SelDB.hsh";
  DirFileHandle::PrintDumper_plReadFile($SelProteinHash_FILE, $SelProteinHash)      if ( ref ($SelProteinHash) eq 'HASH' ) ;
  
  
  my $Table_from_OutHash=MatrixCsvChange::Change2dHashIntoCsvString_withSortFunc($SelProteinHash,0,'','','','');


  my $Table_from_OutHash_FILE=$workingDIRpath."/0_6_0_SelDB__Table.txt";
  FastaFileHandle::BuildFastaFile_withFastaString ($Table_from_OutHash_FILE, $Table_from_OutHash);
  
  #system ("tar -cf test.txt.clstr_dir.tar $CDhit_out_WorkingDIR");
  my $relativeWkDIRpath= File::Spec->abs2rel( $workingDIRpath );
  ExcelHandle::TarGZ_cmd ( $relativeWkDIRpath ); 
  
  
	return $SelProteinHash;
	
	
}


sub GetGenBank_Inform_from_BLASTpraseFile{ #SelenoProteinData::GetGenBank_Inform_from_BLASTpraseFile ($inBlastPrsHASH, $TempQrGBFile, $EachQrGBFile);
	my ($inBlastPrsFile, $TempQrGBFile, $EachQrGBFile)=@_;
	
	my $warnMsgBody="\nIn package  SelenoProteinData,\tIn sub GetGenBank_Inform_from_BLASTpraseFile,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	if   (   (  defined ( $inBlastPrsFile )  ) && ( $inBlastPrsFile=~m/\S+/ ) && (  ( -e $inBlastPrsFile )  )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$inBlastPrsFile=$inBlastPrsFile should be a defined not empty file  !!  $!\n\n\n".$caller_inform ); 	}
	my $inBlastPrsHASH=Storable::retrieve( $inBlastPrsFile );
	my $Genk_HASH=SelenoProteinData::GetGenBank_Inform_from_BLASTpraseFile_HASH($inBlastPrsHASH, $TempQrGBFile, $EachQrGBFile);
	return $Genk_HASH;
}


sub GetGenBank_Inform_from_BLASTpraseFile_HASH{  #SelenoProteinData::GetGenBank_Inform_from_BLASTpraseFile_HASH ($inBlastPrsHASH, $TempQrGBFile, $EachQrGBFile);
	my ($inBlastPrsHASH, $TempQrGBFile, $EachQrGBFile)=@_;
	
	my $warnMsgBody="\nIn package  SelenoProteinData,\tIn sub GetGenBank_Inform_from_BLASTpraseFile_HASH,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	if  (   (  defined ( $inBlastPrsHASH->{'0_0_1_Query_array'} )  ) && (  ref ( $inBlastPrsHASH->{'0_0_1_Query_array'} ) eq 'ARRAY'  )   ){}
  else{		DieWork::Just_dieWork( $die_MsgHead."\n \$inBlastPrsHASH->{'0_0_1_Query_array'}=$inBlastPrsHASH->{'0_0_1_Query_array'} should be a ARRAY  !!  $!\n\n\n".$caller_inform ); 	}
	
	my $AllHit_HASH=BlastPraser::New_GetAllHit_HASH_fromBlastRstArrayHash( $inBlastPrsHASH );
	my $Genk_HASH=GenBankHandle::Upgrade_db_and_geting_GB_inform_from_ProtHASH($AllHit_HASH, $TempQrGBFile, $EachQrGBFile);
	return $Genk_HASH;
}


# build a hash holding information for the 2 old text excel files
#input 3  value   $fake_excel_selProInform_1 $fake_excel_selProInform_2   $big_raw_inform_Hash_file
#                     my $big_raw_inform_Hash=SelenoProteinData::Build_the_raw_informHASH_from_oldEXCELfile();
sub Build_the_raw_informHASH_from_oldEXCELfile{
	my ($fake_excel_selProInform_1, $fake_excel_selProInform_2, $big_raw_inform_Hash_file)=@_;
	
	my ($warnMsgBody, $warnMsgHead, $die_MsgHead, $caller_inform)=@{ DieWork::BuildWarnDieHeadInform( 'SelenoProteinData', 'Build_the_raw_informHASH_from_oldEXCELfile' ) };
  
  DieWork::Check_FileDirExist_or_DIE   ( $fake_excel_selProInform_1,     "\$fake_excel_selProInform_1",    $die_MsgHead, $caller_inform  );  
	DieWork::Check_FileDirExist_or_DIE   ( $fake_excel_selProInform_2,     "\$fake_excel_selProInform_2",    $die_MsgHead, $caller_inform  );  
	DieWork::Check_DfdNoEmptString_or_DIE( $big_raw_inform_Hash_file,      "\$big_raw_inform_Hash_file",     $die_MsgHead, $caller_inform  );  
	
	
	my $fake_excel_HASH_selProInform_1=InFileHandle::FormTableToHash  ($fake_excel_selProInform_1, 0);
	my $fake_excel_HASH_selProInform_2=InFileHandle::FormTableToHash  ($fake_excel_selProInform_2, 0);
	
	my $outReal=ArrayHashChange::CheckHashKeyPreat ($fake_excel_HASH_selProInform_1, $fake_excel_HASH_selProInform_2);
	
	my $big_raw_inform_Hash=ArrayHashChange::PushSmallHash_into_bigHash($fake_excel_HASH_selProInform_1, $fake_excel_HASH_selProInform_2);
	
	DirFileHandle::PrintDumper_plReadFile($big_raw_inform_Hash_file, $big_raw_inform_Hash) if ( ref($big_raw_inform_Hash) eq 'HASH' ) ; 

	
	return  $big_raw_inform_Hash;
}



#use 5 hashs to build selenoprotein hash
#my $SelProteinHash=SelenoProteinData::Build_SelenoProteinHASH( $oldDataFastaHash, $big_raw_inform_Hash, $Total_Bshit_genank_InfomHASH,  $Cov1ID1_hit_gebank_InfomHASH, $Total_Best_blastQueryHit_hash, $Cov1Id1_Best_blstQeryHit_hash, $cluster_sepKeyHASH, $cluster_cluID_AlnHASH);
sub Build_SelenoProteinHASH{
	my ( $oldDataFastaHash, $big_raw_inform_Hash, $Total_Bshit_genank_InfomHASH,  $Cov1ID1_hit_gebank_InfomHASH, $Total_Best_blastQueryHit_hash, $Cov1Id1_Best_blstQeryHit_hash, $cluster_sepKeyHASH, $cluster_cluID_AlnHASH)=@_;
	
	my $warnMsgBody="\nIn package  SelenoProteinData,\tIn sub Build_SelenoProteinHASH,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	DieWork::Check_Hash_or_DIE( $oldDataFastaHash,              "\$oldDataFastaHash",              $die_MsgHead, $caller_inform  );
	DieWork::Check_Hash_or_DIE( $big_raw_inform_Hash,           "\$big_raw_inform_Hash",           $die_MsgHead, $caller_inform  );
	DieWork::Check_Hash_or_DIE( $Total_Bshit_genank_InfomHASH,  "\$Total_Bshit_genank_InfomHASH",  $die_MsgHead, $caller_inform  );
	DieWork::Check_Hash_or_DIE( $Cov1ID1_hit_gebank_InfomHASH,  "\$Cov1ID1_hit_gebank_InfomHASH",  $die_MsgHead, $caller_inform  );
	DieWork::Check_Hash_or_DIE( $Total_Best_blastQueryHit_hash, "\$Total_Best_blastQueryHit_hash", $die_MsgHead, $caller_inform  );
	DieWork::Check_Hash_or_DIE( $Cov1Id1_Best_blstQeryHit_hash, "\$Cov1Id1_Best_blstQeryHit_hash", $die_MsgHead, $caller_inform  );
	DieWork::Check_Hash_or_DIE( $cluster_sepKeyHASH,            "\$cluster_sepKeyHASH",            $die_MsgHead, $caller_inform  );
	DieWork::Check_Hash_or_DIE( $cluster_cluID_AlnHASH,         "\$cluster_cluID_AlnHASH",         $die_MsgHead, $caller_inform  );
	
	my $SelProteinHash;
	
	foreach my $SelDB_ProID (    sort { $a cmp $b } (   keys (  %{ $oldDataFastaHash } )   )    ){ 
	  DieWork::Check_DfdNoEmptString_or_DIE( $oldDataFastaHash->{$SelDB_ProID}->{'0_0_4_orgShort__ID'}, "\$oldDataFastaHash->{\$SelDB_ProID}->{'0_0_4_orgShort__ID'}=\$oldDataFastaHash->{$SelDB_ProID}->{'0_0_4_orgShort__ID'}", $die_MsgHead, $caller_inform  );
	  my $shortJLdb_ID=$oldDataFastaHash->{$SelDB_ProID}->{'0_0_4_orgShort__ID'};
	  $SelProteinHash->{ $shortJLdb_ID }=Storable::dclone( $oldDataFastaHash->{$SelDB_ProID} );
	  
	  DieWork::Check_DfdNoEmptString_or_DIE( $cluster_sepKeyHASH->{$SelDB_ProID}, "\$cluster_sepKeyHASH->{\$SelDB_ProID}=\$cluster_sepKeyHASH->{$SelDB_ProID}", $die_MsgHead, $caller_inform  );
	  my $clusterNB=$cluster_sepKeyHASH->{$SelDB_ProID};
	  $SelProteinHash->{ $shortJLdb_ID }->{'000000__cluster_NB'            }=$clusterNB;
	  
	  DieWork::Check_DfdNoEmptString_or_DIE( $cluster_cluID_AlnHASH->{$clusterNB}, "\$cluster_cluID_AlnHASH->{\$clusterNB}=\$cluster_cluID_AlnHASH->{$clusterNB}", $die_MsgHead, $caller_inform  );
	  foreach my $cluKeyWord (    sort { $a cmp $b } (   keys (  %{ $cluster_cluID_AlnHASH->{$clusterNB} } )   )    ){  DieWork::Print_and_warn( "\n201910291639-0-0-0 \$SelProteinHash->{ \$shortJLdb_ID }->{'1_0_0_'.$cluKeyWord}=\$cluster_cluID_AlnHASH->{$clusterNB}->{$cluKeyWord}=\$SelProteinHash->{ $shortJLdb_ID }->{'1_0_0_'.$cluKeyWord}=$cluster_cluID_AlnHASH->{$clusterNB}->{$cluKeyWord} \n\n" );
	    $SelProteinHash->{ $shortJLdb_ID }->{'1_0_0_'.$cluKeyWord}=$cluster_cluID_AlnHASH->{$clusterNB}->{$cluKeyWord}; 
	  }
	  
	  
	  if ( DieWork::Check_Hash_or_NOT( $big_raw_inform_Hash->{$shortJLdb_ID} ) ){ 
	  	foreach my $OldExcelKey (    sort { $a cmp $b } (   keys (  %{ $big_raw_inform_Hash->{$shortJLdb_ID} } )   )    ){  DieWork::Print_and_warn( "\n20191029-0-1-0 No \$big_raw_inform_Hash->{$shortJLdb_ID}->{$OldExcelKey}=$big_raw_inform_Hash->{$shortJLdb_ID}->{$OldExcelKey} \n\n" );
	  	  $SelProteinHash->{ $shortJLdb_ID }->{'2_0_0_'.$OldExcelKey}=$big_raw_inform_Hash->{$shortJLdb_ID}->{$OldExcelKey}; 
	  	}	  	
	  } 
	  else{
	  	DieWork::Print_and_warn( "\n20191029-0-0-0 No \$SelDB_ProID=$SelDB_ProID found in \$big_raw_inform_Hash=$big_raw_inform_Hash\n\n" );
	  }
	  
	  if   ( DieWork::Check_Hash_or_NOT( $Cov1Id1_Best_blstQeryHit_hash->{$SelDB_ProID} ) ){
	  	$SelProteinHash->{ $shortJLdb_ID }->{'3_0_0_3_0_1__name'            }=$Cov1Id1_Best_blstQeryHit_hash->{$SelDB_ProID}->{'3_0_1__name'            };
	  	$SelProteinHash->{ $shortJLdb_ID }->{'3_0_0_3_0_2_giNub'            }=$Cov1Id1_Best_blstQeryHit_hash->{$SelDB_ProID}->{'3_0_2_giNub'            };
	  	my $AccID_lucsID=
	  	$SelProteinHash->{ $shortJLdb_ID }->{'3_0_0_3_0_3_accessionNB'      }=$Cov1Id1_Best_blstQeryHit_hash->{$SelDB_ProID}->{'3_0_3_accessionNB'      };
	  	$SelProteinHash->{ $shortJLdb_ID }->{'3_0_0_3_0_4_hitSeqLength'     }=$Cov1Id1_Best_blstQeryHit_hash->{$SelDB_ProID}->{'3_0_4_hitSeqLength'     };
	  	$SelProteinHash->{ $shortJLdb_ID }->{'3_0_0_3_0_5_hitDescription'   }=$Cov1Id1_Best_blstQeryHit_hash->{$SelDB_ProID}->{'3_0_5_hitDescription'   };
	  	$SelProteinHash->{ $shortJLdb_ID }->{'3_0_0_3_0_6_hitEvalue'        }=$Cov1Id1_Best_blstQeryHit_hash->{$SelDB_ProID}->{'3_0_6_hitEvalue'        };
	  	$SelProteinHash->{ $shortJLdb_ID }->{'3_0_0_3_2_0_hitTotalIdentical'}=$Cov1Id1_Best_blstQeryHit_hash->{$SelDB_ProID}->{'3_2_0_hitTotalIdentical'};
	  	$SelProteinHash->{ $shortJLdb_ID }->{'3_0_0_3_2_1_hitTotalConserved'}=$Cov1Id1_Best_blstQeryHit_hash->{$SelDB_ProID}->{'3_2_1_hitTotalConserved'};
	  	$SelProteinHash->{ $shortJLdb_ID }->{'3_0_0_3_2_5_total_qry_covrAte'}=$Cov1Id1_Best_blstQeryHit_hash->{$SelDB_ProID}->{'3_2_5_total_qry_covrAte'};
	  	$SelProteinHash->{ $shortJLdb_ID }->{'3_0_0_3_2_7_hitCoverRate'     }=$Cov1Id1_Best_blstQeryHit_hash->{$SelDB_ProID}->{'3_2_7_hitCoverRate'     };
	  	
	  	#DieWork::Check_Hash_or_DIE( $Cov1ID1_hit_gebank_InfomHASH->{$AccID_lucsID}, "\$Cov1ID1_hit_gebank_InfomHASH->{\$AccID_lucsID}=\$Cov1ID1_hit_gebank_InfomHASH->{$AccID_lucsID}", $die_MsgHead, $caller_inform );
	  	  
	  	if   ( DieWork::Check_Hash_or_NOT( $Cov1ID1_hit_gebank_InfomHASH->{$AccID_lucsID} ) ){
	  	  foreach my $GenbankKey (    sort { $a cmp $b } (   keys (  %{ $Cov1ID1_hit_gebank_InfomHASH->{$AccID_lucsID} } )   )    ){
	  	  	if   ( DieWork::Check_Array_or_NOT( $Cov1ID1_hit_gebank_InfomHASH->{$AccID_lucsID}->{$GenbankKey} ) ){
	  	  		my $outString=ArrayHashChange::ChangArrayIntoSimpleString ( $Cov1ID1_hit_gebank_InfomHASH->{$AccID_lucsID}->{$GenbankKey} );
	  	  		$SelProteinHash->{ $shortJLdb_ID }->{'4_0_0_'.'5_0_phyClsif'     }=$outString;
	  	  		$SelProteinHash->{ $shortJLdb_ID }->{'4_0_0_'.'6_0_phyClsifArray'}=Storable::dclone( $Cov1ID1_hit_gebank_InfomHASH->{$AccID_lucsID}->{$GenbankKey} );
	  	  		
	  	  	}
	  	  	else{
	  	  		$SelProteinHash->{ $shortJLdb_ID }->{'4_0_0_'.$GenbankKey}=$Cov1ID1_hit_gebank_InfomHASH->{$AccID_lucsID}->{$GenbankKey};
	  	  	}
	  	  }
	  	}
	  	else{
	  		DieWork::Print_and_warn( "\n20191029-0-0-1-0 No \$AccID_lucsID=$AccID_lucsID found in \$Cov1ID1_hit_gebank_InfomHASH=$Cov1ID1_hit_gebank_InfomHASH\n\n" );
	  	}
	  	
	  }
	  elsif( DieWork::Check_Hash_or_NOT( $Total_Best_blastQueryHit_hash->{$SelDB_ProID} ) ){
	  	$SelProteinHash->{ $shortJLdb_ID }->{'3_0_0_3_0_1__name'            }=$Total_Best_blastQueryHit_hash->{$SelDB_ProID}->{'3_0_1__name'            };
	  	$SelProteinHash->{ $shortJLdb_ID }->{'3_0_0_3_0_2_giNub'            }=$Total_Best_blastQueryHit_hash->{$SelDB_ProID}->{'3_0_2_giNub'            };
	  	my $AccID_lucsID=
	  	$SelProteinHash->{ $shortJLdb_ID }->{'3_0_0_3_0_3_accessionNB'      }=$Total_Best_blastQueryHit_hash->{$SelDB_ProID}->{'3_0_3_accessionNB'      };
	  	$SelProteinHash->{ $shortJLdb_ID }->{'3_0_0_3_0_4_hitSeqLength'     }=$Total_Best_blastQueryHit_hash->{$SelDB_ProID}->{'3_0_4_hitSeqLength'     };
	  	$SelProteinHash->{ $shortJLdb_ID }->{'3_0_0_3_0_5_hitDescription'   }=$Total_Best_blastQueryHit_hash->{$SelDB_ProID}->{'3_0_5_hitDescription'   };
	  	$SelProteinHash->{ $shortJLdb_ID }->{'3_0_0_3_0_6_hitEvalue'        }=$Total_Best_blastQueryHit_hash->{$SelDB_ProID}->{'3_0_6_hitEvalue'        };
	  	$SelProteinHash->{ $shortJLdb_ID }->{'3_0_0_3_2_0_hitTotalIdentical'}=$Total_Best_blastQueryHit_hash->{$SelDB_ProID}->{'3_2_0_hitTotalIdentical'};
	  	$SelProteinHash->{ $shortJLdb_ID }->{'3_0_0_3_2_1_hitTotalConserved'}=$Total_Best_blastQueryHit_hash->{$SelDB_ProID}->{'3_2_1_hitTotalConserved'};
	  	$SelProteinHash->{ $shortJLdb_ID }->{'3_0_0_3_2_5_total_qry_covrAte'}=$Total_Best_blastQueryHit_hash->{$SelDB_ProID}->{'3_2_5_total_qry_covrAte'};
	  	$SelProteinHash->{ $shortJLdb_ID }->{'3_0_0_3_2_7_hitCoverRate'     }=$Total_Best_blastQueryHit_hash->{$SelDB_ProID}->{'3_2_7_hitCoverRate'     };
	  	
	  	
	  	if   ( DieWork::Check_Hash_or_NOT( $Total_Bshit_genank_InfomHASH->{$AccID_lucsID} ) ){
	  	  foreach my $GenbankKey (    sort { $a cmp $b } (   keys (  %{ $Total_Bshit_genank_InfomHASH->{$AccID_lucsID} } )   )    ){
	  	  	if   ( DieWork::Check_Array_or_NOT( $Total_Bshit_genank_InfomHASH->{$AccID_lucsID}->{$GenbankKey} ) ){
	  	  		my $outString=ArrayHashChange::ChangArrayIntoSimpleString ( $Total_Bshit_genank_InfomHASH->{$AccID_lucsID}->{$GenbankKey} );
	  	  		$SelProteinHash->{ $shortJLdb_ID }->{'4_0_0_'.'5_0_phyClsif'     }=$outString;
	  	  		$SelProteinHash->{ $shortJLdb_ID }->{'4_0_0_'.'6_0_phyClsifArray'}=Storable::dclone( $Total_Bshit_genank_InfomHASH->{$AccID_lucsID}->{$GenbankKey} );
	  	  		
	  	  	}
	  	  	else{
	  	  		$SelProteinHash->{ $shortJLdb_ID }->{'4_0_0_'.$GenbankKey}=$Total_Bshit_genank_InfomHASH->{$AccID_lucsID}->{$GenbankKey};
	  	  	}
	  	  }
	  	}
	  	else{
	  		DieWork::Print_and_warn( "\n20191029-0-0-1-1 No \$AccID_lucsID=$AccID_lucsID found in \$Total_Bshit_genank_InfomHASH=$Total_Bshit_genank_InfomHASH\n\n" );
	  	}
	  	
	  	
	  	
	  }
	  else{
	  	#DieWork::Just_dieWork( $die_MsgHead."\n \$SelDB_ProID=$SelDB_ProID cannot be found in \$Total_Best_blastQueryHit_hash=$Total_Best_blastQueryHit_hash  !!  $!\n\n\n".$caller_inform ); 
	    DieWork::Print_and_warn( "\n20191029-0-0-1 No \$SelDB_ProID=$SelDB_ProID found in \$Total_Best_blastQueryHit_hash=$Total_Best_blastQueryHit_hash\n\n" );
	  }
	  
	  
	}
	return $SelProteinHash;
}













#SelenoProteinData::With_manual_HANDLE_excel_build_family_HASH( $inExcelTXTfile, $mutiple_core, $outRstDir, $key_col_1d, $key_col_2d, $GB_id_HEAD  );
sub With_manual_HANDLE_excel_build_family_HASH{  
	my ( $inExcelTXTfile, 
	     $mutiple_core,
	     $outRstDir,
	     $key_col_1d,
	     $key_col_2d,
	     $GB_id_HEAD   #short name for Selenoprotein Family
	     
	   )=@_;  #my $inExcelTXTfile="old_xls_informaltil_add_group.excel.txt";
	
	my ($warnMsgBody, $warnMsgHead, $die_MsgHead, $caller_inform)=@{ DieWork::BuildWarnDieHeadInform( 'SelenoProteinData', 'With_manual_HANDLE_excel_build_family_HASH' ) };
  
  DieWork::Check_FileDirExist_or_DIE( $inExcelTXTfile,     "\$inExcelTXTfile",    $die_MsgHead, $caller_inform  );  
  if ( DieWork::Check_DfdNoEmptNUMBER_or_NOT( $mutiple_core )  ){ } else { $mutiple_core=10;   } 	  
	if ( DieWork::Check_DfdNoEmptString_or_NOT( $outRstDir    )  ){ } else { $outRstDir   =$inExcelTXTfile.'_outDIR';   } 
	my $mkoutRstDirCMD=" mkdir -p $outRstDir"; if (  -d ($outRstDir)  ) {} else {   DieWork::Print_and_warn( $mkoutRstDirCMD."\n" );  system ( "$mkoutRstDirCMD" ); }
	
	
	if ( DieWork::Check_DfdNoEmptNUMBER_or_NOT( $key_col_1d )  ){ } else { $key_col_1d=46;   } 	  
	if ( DieWork::Check_DfdNoEmptNUMBER_or_NOT( $key_col_2d )  ){ } else { $key_col_2d=1;    } 	 
  my $raw_groupedHash=InFileHandle::FormTableToHash_to_2D_hash  ($inExcelTXTfile, $key_col_1d, $key_col_2d );
  DieWork::Check_Hash_or_DIE( $raw_groupedHash,            "\$raw_groupedHash",    $die_MsgHead, $caller_inform  );  
  my $raw_groupedHash_readFILE=$outRstDir."/0m_in_excel_changed_HASH.txt";   DirFileHandle::PrintDumper_plReadFile($raw_groupedHash_readFILE, $raw_groupedHash) if ( ref($raw_groupedHash) eq 'HASH' ) ; 

  my $Wei_shu=8;
  
  if ( DieWork::Check_DfdNoEmptString_or_NOT( $GB_id_HEAD    )  ){ } else { $GB_id_HEAD   ='SF';   } 
 

  
  my $KEY_nmuberID_VAL_wordID___HASH;
  my $KEY_wordID___VAL_nmuberID_HASH;
  
  my $newDB_GrpID_SedID_SEQ_HASH;
  my $Group_NB=1;
  foreach my $SelFml_NAME (    sort { $a cmp $b } (   keys (  %{ $raw_groupedHash } )   )    ){                                       print "    \$SelFml_NAME=$SelFml_NAME\n";
    if ( DieWork::Check_DfdNoEmptString_or_NOT( $SelFml_NAME    )  ){
    
      my $showGpNB=  MatrixCsvChange::SprintfKeyHead_inPut_weishu( $Wei_shu, $Group_NB );  $showGpNB=$GB_id_HEAD.'_'.$showGpNB;         print "    \$showGpNB=$showGpNB\n"; 
      
      $KEY_nmuberID_VAL_wordID___HASH->{$showGpNB}=$SelFml_NAME;
      $KEY_wordID___VAL_nmuberID_HASH->{$SelFml_NAME}=$showGpNB;    
      
      DieWork::Check_Hash_or_DIE( $raw_groupedHash->{$SelFml_NAME},                      "\$raw_groupedHash->{$SelFml_NAME}",         $die_MsgHead, $caller_inform  );   
      foreach my $SelProJLDB_ID (    sort { $a cmp $b } (   keys (  %{ $raw_groupedHash->{$SelFml_NAME} }  )   )    ){  
        $newDB_GrpID_SedID_SEQ_HASH->{$showGpNB}->{$SelProJLDB_ID}=$raw_groupedHash->{$SelFml_NAME}->{$SelProJLDB_ID}->{'0_0_6_____sequence'};
      }    
      $Group_NB++; 
    }
  }
  
  my $KEY_nmuberID_VAL_wordID___HASH_file=$outRstDir."/".'1m_1m_SFid_SFname_HASH';    DirFileHandle::PrintDumper_plReadFile($KEY_nmuberID_VAL_wordID___HASH_file, $KEY_nmuberID_VAL_wordID___HASH) if ( ref($KEY_nmuberID_VAL_wordID___HASH) eq 'HASH' ) ;  
  my $KEY_wordID___VAL_nmuberID_HASH_file=$outRstDir."/".'1m_2m_SFname_SFid_HASH';    DirFileHandle::PrintDumper_plReadFile($KEY_wordID___VAL_nmuberID_HASH_file, $KEY_nmuberID_VAL_wordID___HASH) if ( ref($KEY_nmuberID_VAL_wordID___HASH) eq 'HASH' ) ;
	my $newDB_GrpID_SedID_SEQ_HASH_file    =$outRstDir."/".'2m_SFid_SqID_Seq_HASH';     DirFileHandle::PrintDumper_plReadFile($newDB_GrpID_SedID_SEQ_HASH_file,     $newDB_GrpID_SedID_SEQ_HASH    ) if ( ref($newDB_GrpID_SedID_SEQ_HASH    ) eq 'HASH' ) ;
  
  
  my $reDIVstepDIR=$outRstDir."/"."0m_reDivStepDIR"; my $mkdirCMD=" mkdir -p $reDIVstepDIR"; if (  -d ($reDIVstepDIR)  ) {} else {   DieWork::Print_and_warn( $mkdirCMD."\n" );  system ( "$mkdirCMD" ); }

  my $PutIntoDIRHASH=DirFileHandle::BuildTreeLevel_DIR_forAGroup_ofFILES($newDB_GrpID_SedID_SEQ_HASH, $reDIVstepDIR);

  my $PutIntoDIRHASH_file=$outRstDir."/"."3m_1m_SFid_to_FastaPATH_HASH";           DirFileHandle::PrintDumper_plReadFile ($PutIntoDIRHASH_file, $PutIntoDIRHASH);
  
  my $fastaFileHASH;
  foreach my $eachSFid (   keys  ( %{ $newDB_GrpID_SedID_SEQ_HASH }  )   ){
  	DieWork::Check_FileDirExist_or_DIE( $PutIntoDIRHASH->{ $eachSFid }, "\$newDB_GrpID_SedID_SEQ_HASH->{ \$eachSFid }=$newDB_GrpID_SedID_SEQ_HASH->{ $eachSFid }", "", ""  );
  	my $fastaString;
  	foreach my $eachSepID (    sort { $a cmp $b } (   keys  ( %{ $newDB_GrpID_SedID_SEQ_HASH->{ $eachSFid } }  )   )    ){
  	  $fastaString.=">$eachSepID\n$newDB_GrpID_SedID_SEQ_HASH->{ $eachSFid }->{ $eachSepID }\n\n";  #DieWork::Print_and_warn( "\$newDB_GrpID_SedID_SEQ_HASH->{ \$eachSFid }->{ $\eachSepID }=\$newDB_GrpID_SedID_SEQ_HASH->{ $eachSFid }->{ $eachSepID }=$newDB_GrpID_SedID_SEQ_HASH->{ $eachSFid }->{ $eachSepID }\n" ); 
  	}
  	my $fastaFILEpath=$PutIntoDIRHASH->{ $eachSFid }."/0_0_fasta.txt";
  	
  	FastaFileHandle::BuildFastaFile_withFastaString ($fastaFILEpath, $fastaString);
  	$fastaFileHASH->{$fastaFILEpath}=$eachSFid;
  }
  my $fastaFileHASH_file=$outRstDir."/"."3m_2m_fastaPath_to_SFid_HASH";  DirFileHandle::PrintDumper_plReadFile ($fastaFileHASH_file, $fastaFileHASH);
  
  
  my $indexFileHASH_file  =$outRstDir."/"."4m_1m_fastaPath_to_idxFilePath_HASH"; 
  my $musFileHASH_file    =$outRstDir."/"."4m_2m_fastaPath_to_muscle_HASH";    
  my $itprscFileHASH_file =$outRstDir."/"."4m_3m_fastaPath_to_interProScan_HASH";
  my $domPNGFileHASH_file =$outRstDir."/"."4m_4m_fastaPath_to_PNG_HASH";
  my $nexTreeFileHASH_file=$outRstDir."/"."4m_5m_fastaPath_to_Mrbayes_HASH";
  
  my $Mutilple_buildIndex_CMD  ="perl /home/fredjiang/fredProgs/PMs/MutThrdPLS/Mutilple_buildIndex.pl   -h $fastaFileHASH_file                          -o $indexFileHASH_file   -c $mutiple_core";                     
  my $Mutilple_muscle_CMD      ="perl /home/fredjiang/fredProgs/PMs/MutThrdPLS/Mutilple_muscle.pl       -h $fastaFileHASH_file                          -o $musFileHASH_file     -c $mutiple_core";                     
  my $Mutilple_InterproScan_CMD="perl /home/fredjiang/fredProgs/PMs/MutThrdPLS/Mutilple_InterproScan.pl -h $fastaFileHASH_file                          -o $itprscFileHASH_file  -c $mutiple_core";                    
  my $Mutilple_PrintPNG_CMD    ="perl /home/fredjiang/fredProgs/PMs/MutThrdPLS/Mutilple_PrintPNG.pl     -h $musFileHASH_file    -p $itprscFileHASH_file -o $domPNGFileHASH_file  -c $mutiple_core";  
  my $Mutilple_MrBays_CMD      ="perl /home/fredjiang/fredProgs/PMs/MutThrdPLS/Mutilple_MrBays.pl       -h $musFileHASH_file                            -o $nexTreeFileHASH_file -c $mutiple_core";                    
  
  DieWork::Print_and_warn( $Mutilple_buildIndex_CMD    ); system ( $Mutilple_buildIndex_CMD   ); #  system ( "perl /home/fredjiang/fredProgs/PMs/MutThrdPLS/Mutilple_buildIndex.pl   -h $fastaFileHASH_file -o $indexFileHASH_file -c $mutiple_core");                       
  DieWork::Print_and_warn( $Mutilple_muscle_CMD        ); system ( $Mutilple_muscle_CMD       ); #  system ( "perl /home/fredjiang/fredProgs/PMs/MutThrdPLS/Mutilple_muscle.pl       -h $fastaFileHASH_file -o $musFileHASH_file -c $mutiple_core");                             
  DieWork::Print_and_warn( $Mutilple_InterproScan_CMD  ); system ( $Mutilple_InterproScan_CMD ); #  system ( "perl /home/fredjiang/fredProgs/PMs/MutThrdPLS/Mutilple_InterproScan.pl -h $fastaFileHASH_file -o $itprscFileHASH_file -c 50");                     
  DieWork::Print_and_warn( $Mutilple_PrintPNG_CMD      ); system ( $Mutilple_PrintPNG_CMD     ); #  system ( "perl /home/fredjiang/fredProgs/PMs/MutThrdPLS/Mutilple_PrintPNG.pl     -h $musFileHASH_file   -p  $itprscFileHASH_file  -o $domPNGFileHASH_file -c 50");     
  DieWork::Print_and_warn( $Mutilple_MrBays_CMD        ); system ( $Mutilple_MrBays_CMD       ); #  system ( "perl /home/fredjiang/fredProgs/PMs/MutThrdPLS/Mutilple_MrBays.pl       -h $musFileHASH_file   -o $nexTreeFileHASH_file  -c 30");                             
                                                                                                                                                                                                                                                                
  
  my $information_hash_putIntoExcel;                                                                           
  $information_hash_putIntoExcel->{$indexFileHASH_file  }->{'1m1m_indexF'}=1;                                  
  $information_hash_putIntoExcel->{$musFileHASH_file    }->{'2m3m_delnum'}=1;                                  
  $information_hash_putIntoExcel->{$itprscFileHASH_file }->{'3m2m_iPsRed'}=1;                                  
  $information_hash_putIntoExcel->{$domPNGFileHASH_file }->{'4m1m_domPng'}=1;                                  
  $information_hash_putIntoExcel->{$nexTreeFileHASH_file}->{'5m1m_MrBays'}=1;                                  
                                                                                                               
  my $new_big_Path_hash=ArrayHashChange::GetBigHash_fromLotsOfSubHASH_FILE( $information_hash_putIntoExcel ); 
  my $new_big_Path_hash_file=$outRstDir."/"."5m_1m_fastaPath_to_TotalPathHASH.txt";                    DirFileHandle::PrintDumper_plReadFile ($new_big_Path_hash_file, $new_big_Path_hash);                        
                                                                                                  
  my $ToHyperLink_HASH=SelenoProteinData::BuildPathExcel_for_SelID_SelFamily($fastaFileHASH_file, $newDB_GrpID_SedID_SEQ_HASH_file, $new_big_Path_hash_file);
  my $ToHyperLink_HASH_file=$outRstDir."/"."5m_2m_TotalPathHASH_needTOchangeHYPER.txt";                DirFileHandle::PrintDumper_plReadFile ($ToHyperLink_HASH_file, $ToHyperLink_HASH);                                                        
  
  my $hylinkHASH;                                                                                              
  $hylinkHASH->{'1m1m_fmlFst'}='1_fmlFst';                                                                     
  $hylinkHASH->{'2m3m_delnum'}='2_delnum';                                                                     
  $hylinkHASH->{'3m2m_iPsRed'}='3_iPsRed';                                                                     
  $hylinkHASH->{'4m1m_domPng'}='4_domPng';                                                                     
  $hylinkHASH->{'5m1m_MrBays'}='5_MrBays';                                                                     
                                                                                                               
  my $hyperLinkedHASH=ExcelHandle::ChangeALLpathTOhyperLINK_in_specifc_colKEYS($ToHyperLink_HASH, $hylinkHASH);
  my $hyperLinkedHASH_file=$outRstDir."/"."6m_1m_TotalPathHASH_HYPERlink.txt";            
  DirFileHandle::PrintDumper_plReadFile ($hyperLinkedHASH_file, $hyperLinkedHASH);                             
                                                                                                               
  my $outExcel=$outRstDir."/"."6m_2m_TotalPathHASH_HYPERlinkExcel.txt";                   
  my $outString= MatrixCsvChange::PrintExcelTableFileFromHASH( $outExcel, $hyperLinkedHASH);                   
  
  my $relativeOutRstDIR= File::Spec->abs2rel( $outRstDir );
                                                                                                                 
  ExcelHandle::TarGZ_cmd ( $relativeOutRstDIR );                                                                      	    
	
}



# my $outHASH=SelenoProteinData::BuildPathExcel_for_SelID_SelFamily($inFasta_SelFaml_HASH_FILE, $inSelFaml_SelID_HASH_FILE, $inFastaKEY_other_path_HASH_FILE);
sub BuildPathExcel_for_SelID_SelFamily{  # input 3 hashDumpFile, 1st fastaPath to selFamily , 2nd SelFamily to Selprotein ID , 3rd fastaPath to other result file path(msf, png, nex .et al). Output the excel input hash, showing all the path to be changed into hyperlink
	my ($inFasta_SelFaml_HASH_FILE, $inSelFaml_SelID_HASH_FILE, $inFastaKEY_other_path_HASH_FILE)=@_;
	   
	my ($warnMsgBody, $warnMsgHead, $die_MsgHead, $caller_inform)=@{ DieWork::BuildWarnDieHeadInform( 'SelenoProteinData', 'BuildPathExcel_for_SelID_SelFamily' ) };

  DieWork::Check_FileDirExist_or_DIE( $inFasta_SelFaml_HASH_FILE,       "\$inFasta_SelFaml_HASH_FILE",        $die_MsgHead, $caller_inform  );
	DieWork::Check_FileDirExist_or_DIE( $inSelFaml_SelID_HASH_FILE,       "\$inSelFaml_SelID_HASH_FILE",        $die_MsgHead, $caller_inform  );
	DieWork::Check_FileDirExist_or_DIE( $inFastaKEY_other_path_HASH_FILE, "\$inFastaKEY_other_path_HASH_FILE",  $die_MsgHead, $caller_inform  );
	
	my $inFasta_SelFaml_HASH       =Storable::retrieve ( $inFasta_SelFaml_HASH_FILE );
	my $inSelFaml_SelID_HASH       =Storable::retrieve ( $inSelFaml_SelID_HASH_FILE );
	my $inFastaKEY_other_path_HASH =Storable::retrieve ( $inFastaKEY_other_path_HASH_FILE );
	
	DieWork::Check_Hash_or_DIE( $inFasta_SelFaml_HASH,       "\$inFasta_SelFaml_HASH",        $die_MsgHead, $caller_inform  );
	DieWork::Check_Hash_or_DIE( $inSelFaml_SelID_HASH,       "\$inSelFaml_SelID_HASH",        $die_MsgHead, $caller_inform  );
	DieWork::Check_Hash_or_DIE( $inFastaKEY_other_path_HASH, "\$inFastaKEY_other_path_HASH",  $die_MsgHead, $caller_inform  );
	
	my $inSelFaml_Fasta_HASH=ArrayHashChange::BuildReverseHash($inFasta_SelFaml_HASH);
	my $Selid_to_SelFML_HASH=ArrayHashChange::Reverse_level1KEY_Level2Key_HASH($inSelFaml_SelID_HASH);
	
	my $outHASH;
	
	foreach my $SelProID (    sort { $a cmp $b } (  keys ( %{ $Selid_to_SelFML_HASH } )  )    ){ 
		
		DieWork::Check_Hash_or_DIE( $Selid_to_SelFML_HASH->{$SelProID}, "\$Selid_to_SelFML_HASH->{\$SelProID}=\$Selid_to_SelFML_HASH->{$SelProID}", $die_MsgHead, $caller_inform  );
		
		foreach my $SfamilyName (    sort { $a cmp $b } (  keys ( %{ $Selid_to_SelFML_HASH->{$SelProID} } )  )    ){ 
		  
		  $outHASH->{$SelProID}->{'0m1m_family'}=$SfamilyName;
		  
		  DieWork::Check_FileDirExist_or_DIE( $inSelFaml_Fasta_HASH->{$SfamilyName}, "\$inSelFaml_Fasta_HASH->{\$SfamilyName}=\$inSelFaml_Fasta_HASH->{$SfamilyName}", $die_MsgHead, $caller_inform  );
	
		  my $fastaFilePath=$inSelFaml_Fasta_HASH->{$SfamilyName};
		  
		  $outHASH->{$SelProID}->{'1m1m_fmlFst'}=$fastaFilePath;
		  
		  DieWork::Check_Hash_or_DIE( $inFastaKEY_other_path_HASH->{$fastaFilePath}, "\$inFastaKEY_other_path_HASH->{\$fastaFilePath}=\$inFastaKEY_other_path_HASH->{$fastaFilePath}", $die_MsgHead, $caller_inform  );
		  
		  foreach my $rst_type_key (    sort { $a cmp $b } (  keys ( %{ $inFastaKEY_other_path_HASH->{$fastaFilePath} } )  )    ){ 
		  	$outHASH->{$SelProID}->{$rst_type_key}=$inFastaKEY_other_path_HASH->{$fastaFilePath}->{$rst_type_key};
		  }
		
		}
		
		
	}
	
	return $outHASH;
}


1;