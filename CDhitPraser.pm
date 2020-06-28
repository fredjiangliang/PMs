#!/usr/bin/perl -w
BEGIN{  push (@INC, '/home/fredjiang/fredProgs/PMs')  };

use strict;
use warnings;
use DieWork;
use DirFileHandle;
use InFileHandle;
                      
package CDhitPraser;
                                            
                                             
my $cd_hit_path="/home/fredjiang/tools/cd-hit-v4.8.1-2019-0228/cd-hit" ;
my $coservationCutoff=0.4;
my $wordLength       =2  ;

                                            
# /home/fredjiang/tools/cd-hit-v4.8.1-2019-0228/cd-hit -i fastafile.txt -o clusterResult.txt -c 0.4 -n 2                # the cluster result will be clusterResult.txt.clstr                            


sub RunCdhit_and_Prase{  # CDhitPraser::RunCdhit_and_Prase  ($inFastaPath, $InputFastaIndexFl, $fasta_Hash_file, $CdHitOutDir);
	
	my ($inFastaPath, $InputFastaIndexFl, $fasta_Hash_file, $CdHitOutDir)=@_;
	
	my ($warnMsgBody, $warnMsgHead, $die_MsgHead, $caller_inform)=@{ DieWork::BuildWarnDieHeadInform( 'CDhitPraser', 'RunCdhit_and_Prase' ) };
	
	DieWork::Check_FileDirExist_or_DIE   ( $inFastaPath,       "\$inFastaPath",       $die_MsgHead, $caller_inform  );  
	DieWork::Check_FileDirExist_or_DIE   ( $InputFastaIndexFl, "\$InputFastaIndexFl", $die_MsgHead, $caller_inform  );  
	DieWork::Check_FileDirExist_or_DIE   ( $fasta_Hash_file,   "\$fasta_Hash_file",   $die_MsgHead, $caller_inform  );  
	DieWork::Check_DfdNoEmptString_or_DIE( $CdHitOutDir,       "\$CdHitOutDir",       $die_MsgHead, $caller_inform  );  
	my $mkoutRstDirCMD=" mkdir -p $CdHitOutDir"; if (  -d ($CdHitOutDir)  ) {} else {   DieWork::Print_and_warn( $mkoutRstDirCMD."\n" );  system ( "$mkoutRstDirCMD" ); }
	
	my $CdHitOutFile                     =$CdHitOutDir."/1m_CdHit_c${coservationCutoff}_n${wordLength}.txt";
	
	CDhitPraser::Run_CdHit( $inFastaPath, $CdHitOutFile, $coservationCutoff, $wordLength );
	
	my $CdHitOutFile_cltrFile=$CdHitOutFile.".clstr";
	DieWork::Check_FileDirExist_or_DIE   ( $CdHitOutFile_cltrFile,       "\$CdHitOutFile_cltrFile",       $die_MsgHead, $caller_inform  );  
	
	
  my $Cdhit_prase_Hsh                  =CDhitPraser::Prase_CDhit_clstr ($CdHitOutFile_cltrFile);
  my $Cdhit_prase_Hsh_FILE             =$CdHitOutDir."/2m_CH_prsOut.hsh";
  DirFileHandle::PrintDumper_plReadFile ($Cdhit_prase_Hsh_FILE, $Cdhit_prase_Hsh);
	
	my $CD_hit_out_ProteinID_keyHASH     =CDhitPraser::Build_ProteinID_KEY_HASH_fromCDhit_outHASH($Cdhit_prase_Hsh);                                          
  my $CD_hit_out_ProteinID_keyHASH_file=$CdHitOutDir."/3m_ProID_to_CluID.hsh";
  DirFileHandle::PrintDumper_plReadFile ($CD_hit_out_ProteinID_keyHASH_file, $CD_hit_out_ProteinID_keyHASH);
  
  #my $InputFastaIndexFl='/home/fredjiang/OneTSSD/fredjiang20190321/TestSELENODBwork/1_test/0_0_0_in_fasta.idx';
  
  my $CDhit_out_WorkingDIR             =$CdHitOutDir."/4m_CdHit_WKingDIR";  
  
  
  #my $fasta_Hash_file="/home/fredjiang/OneTSSD/fredjiang20190321/TestSELENODBwork/1_test/0_0_0_in_fasta.hsh";
  my $fasta_Hash=Storable::retrieve($fasta_Hash_file);
  DieWork::Check_Hash_or_DIE            ( $fasta_Hash,                    "\$fasta_Hash",                     $die_MsgHead, $caller_inform  );
	
  my $outAlignmentHASH;
  $outAlignmentHASH     = CDhitPraser::Build_mutiple_alignment_based_on_CDhit( $Cdhit_prase_Hsh_FILE, $InputFastaIndexFl, $fasta_Hash, $CDhit_out_WorkingDIR );
  
  
  
  my $outAlignmentHASH_file            =$CdHitOutDir.'/5m_CdHit_almPath.hsh';
  DirFileHandle::PrintDumper_plReadFile($outAlignmentHASH_file, $outAlignmentHASH) ; 
  
  my $rltv_CDhit_out_WorkingDIR=File::Spec->abs2rel( $CDhit_out_WorkingDIR );
  #ExcelHandle::TarGZ_cmd ( $rltv_CDhit_out_WorkingDIR); 
  
  
  
  
  my $Table_from_OutHash=MatrixCsvChange::Change2dHashIntoCsvString_withSortFunc($outAlignmentHASH,0,'','','','');
  
  
  
  #my $outTable                         =$CdHitOutDir.'/6m_CDhit_IPS.xlsx';
  
  #FastaFileHandle::BuildFastaFile_withFastaString ($outTable, $Table_from_OutHash);
	
	return [ $CD_hit_out_ProteinID_keyHASH_file, $outAlignmentHASH_file ];
}


sub Run_CdHit{  # CDhitPraser::Run_CdHit( $inFastaPath, $outFileName, $coservationCutoff, $wordLength );
	my ($inFastaPath, $outFileName, $coservationCutoff, $wordLength)=@_;
	
	my ($warnMsgBody, $warnMsgHead, $die_MsgHead, $caller_inform)=@{ DieWork::BuildWarnDieHeadInform( 'CDhitPraser', 'Run_CdHit' ) };
  
  DieWork::Check_FileDirExist_or_DIE   ( $inFastaPath,       "\$inFastaPath",       $die_MsgHead, $caller_inform  );  
	DieWork::Check_DfdNoEmptString_or_DIE( $outFileName,       "\$outFileName",       $die_MsgHead, $caller_inform  );  
	DieWork::Check_DfdNoEmptString_or_DIE( $coservationCutoff, "\$coservationCutoff", $die_MsgHead, $caller_inform  );  
	DieWork::Check_DfdNoEmptString_or_DIE( $wordLength,        "\$wordLength",        $die_MsgHead, $caller_inform  );  
	
	my $CdHit_CMD= "$cd_hit_path -i $inFastaPath -o $outFileName -c $coservationCutoff -n $wordLength ";
	DieWork::Print_and_warn( $CdHit_CMD );
	system ($CdHit_CMD);
	
}

sub Prase_CDhit_clstr{  #CDhitPraser::Prase_CDhit_clstr ($inFile);
	my ($inFile)=@_;
	
  my $warnMsgBody="\nIn package  CDhitPraser,\tIn sub Prase_CDhit_clstr,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	if   (   (  defined ( $inFile )  ) && ( $inFile=~m/\S+/ ) && (  ( -e $inFile )  )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$inFile=$inFile should be a defined not empty file  !!  $!\n\n\n".$caller_inform ); 	}
  
  my $outArray=InFileHandle::readAllfileInto_a_Array($inFile, "\n>Cluster");
  #DirFileHandle::PrintDumper($inFile."\.hsh", $outArray) ;
  
  my $outHASH;
  my $total_cluster=0;
  
  if  (   (  defined ( $outArray )  ) && (  ref ( $outArray ) eq 'ARRAY'  )   ){
    for ( my $i=0; $i<@{ $outArray }; $i++ ){    #DieWork::Print_and_warn( "\n20191022-0-0-0 \$outArray->[$i]=$outArray->[$i]\n\n" ); 
    	
    	if ($outArray->[$i]=~m/
    	                         \s+(\d+)\n
                               (?:
                                  
                                  #(?:
#                                  (\d+\s+\d+aa,\s+>\S+.*\n\+#\.\.\.\s)
#                                  )
                                  ( #?:
                                  \d+\s+\d+aa,\s+>\S+.*\n\S+
                                    (?:\n
                                      \S+
                                      .*
                                    )?
                                    \.\.\.\s
                                  )
                                  (?:(?:\s+at\s+\d+%\n)?|(?:\*\n)?)
                                  
                                  
                               )+  
                         
                            /x
         )
      { my $cluster_NB=$1;        DieWork::Print_and_warn( "\n20200604-0-0-0 \$2=$2\n\n" );                                               
        my @eachSingle=( $outArray->[$i]=~m/                         
                                                (                   
                                                  (?:\d+\s+\d+aa,\s+>\S+.*\n\S+(?:\n\S+.*)?\.\.\.\s)
                                                  (?:(?:\s+at\s+\d+%\n)?|(?:\*\n)?)
            
                                                )                    
                                                                      
  	                                  	    /xg 
  	                     );
  	    my $clustal_in_idx=0;                             	                           
        foreach my $eachSigLeElmt (@eachSingle) {          
         	if ($eachSigLeElmt=~m/
                                  (\d+)\s+(\d+)aa,\s+>(\S+.*\n\S+)(?:\n\S+.*)?\.\.\.\s
                                  (?:\s+at\s+(\d+)%\n)?
                                  
                         
                            /x
             )
         	{
         		#DieWork::Print_and_warn( "\n20191022-0-0-0 \$1=$1 \$2=$2 \$3=$3 \$4=$4\n\n" );  
         		#DieWork::Print_and_warn( "\n20191022-0-0-0 \$5=$5\n\n" )  if ( defined ($5) );
         		
         		$outHASH->{'0_1_0_clusterHASH'}->{$cluster_NB}->{'1_1_0_inClusterArray'}->[$clustal_in_idx]->{'0_1_0_outCluster_NB'}=$cluster_NB;
         		
         		$outHASH->{'0_1_0_clusterHASH'}->{$cluster_NB}->{'1_1_0_inClusterArray'}->[$clustal_in_idx]->{'0_1_1_inCluster_pepIdx'}=$1;
         		$outHASH->{'0_1_0_clusterHASH'}->{$cluster_NB}->{'1_1_0_inClusterArray'}->[$clustal_in_idx]->{'0_1_3_pep_length'}=$2;
         		$outHASH->{'0_1_0_clusterHASH'}->{$cluster_NB}->{'1_1_0_inClusterArray'}->[$clustal_in_idx]->{'0_1_2_pep_____ID'}=$3;
         		#$outHASH->{'0_1_0_clusterHASH'}->{$cluster_NB}->{'1_1_0_inClusterArray'}->[$clustal_in_idx]->{'0_1_4_pep_ShortSeq'}=$4;
         		if ( defined ($5) ){
         			my $Coserv_pct=$5."%";
         			my $Coserv_nub=$5/100;
         			$outHASH->{'0_1_0_clusterHASH'}->{$cluster_NB}->{'1_1_0_inClusterArray'}->[$clustal_in_idx]->{'0_1_5_pep_Coserv_nub'}=$Coserv_nub;
         		  $outHASH->{'0_1_0_clusterHASH'}->{$cluster_NB}->{'1_1_0_inClusterArray'}->[$clustal_in_idx]->{'0_1_6_pep_Coserv_pct'}=$Coserv_pct;
         		}
         		
         	}
         	else{
         		DieWork::Just_dieWork( $die_MsgHead."\n \$eachSigLeElmt=$eachSigLeElmt didnot fit for the regular expression  !!  $!\n\n\n".$caller_inform ); 
         	}
          $clustal_in_idx++;	 
          $outHASH->{'0_1_0_clusterHASH'}->{$cluster_NB}->{'1_0_0_total_pep_in_cluster'}=$clustal_in_idx; 
      	}
      }
      else{
      	DieWork::Just_dieWork( $die_MsgHead."\n \$outArray->[$i]=$outArray->[$i] didnot fit for the regular expression  !!  $!\n\n\n".$caller_inform ); 	
      }
    	$total_cluster++;
    	
    }
  }
  
  $outHASH->{'0_0_0_total_cluster_nb'}=$total_cluster;
  return $outHASH;
}


#build the ProteinID_KEY to cluster_ID hash from the CD-hit prased out hash
#my $ProteinID_KEY_HASH_fromCDhit_outHASH=CDhitPraser::Build_ProteinID_KEY_HASH_fromCDhit_outHASH( $CdHitOutHash );
sub Build_ProteinID_KEY_HASH_fromCDhit_outHASH{
	my ( $CdHitOutHash )=@_;
	
	my ($warnMsgBody, $warnMsgHead, $die_MsgHead, $caller_inform)=@{ DieWork::BuildWarnDieHeadInform( 'CDhitPraser', 'Build_ProteinID_KEY_HASH_fromCDhit_outHASH' ) };
	
	DieWork::Check_Hash_or_DIE( $CdHitOutHash->{'0_1_0_clusterHASH'}, "\$CdHitOutHash->{'0_1_0_clusterHASH'}", $die_MsgHead, $caller_inform  );
	my $ProteinID_KEY_HASH_fromCDhit_outHASH;
	foreach my $cluster_nb (    sort { $a <=> $b } (   keys (  %{ $CdHitOutHash->{'0_1_0_clusterHASH'} } )   )    ){  #DieWork::Print_and_warn( "\n201910291554-0-0-0  \$cluster_nb=$cluster_nb \n\n" );
	  
	  DieWork::Check_Array_or_DIE( $CdHitOutHash->{'0_1_0_clusterHASH'}->{$cluster_nb}->{'1_1_0_inClusterArray'}, "\$CdHitOutHash->{'0_1_0_clusterHASH'}->{$cluster_nb}->{'1_1_0_inClusterArray'}", $die_MsgHead, $caller_inform  );
	  
	  my $inClusterArray=$CdHitOutHash->{'0_1_0_clusterHASH'}->{$cluster_nb}->{'1_1_0_inClusterArray'};
	  
	  
	  for (  my $i=0; $i < @{ $inClusterArray } ; $i++  ){  #DieWork::Print_and_warn( "\n201910291554-0-0-1  \$i=$i \n\n" );
	  	
	  	$ProteinID_KEY_HASH_fromCDhit_outHASH->{ $inClusterArray->[$i]->{'0_1_2_pep_____ID'} }=$inClusterArray->[$i]->{'0_1_0_outCluster_NB'}; 
	  	#DieWork::Print_and_warn( "\n201910291554-0-0-2  \$ProteinID_KEY_HASH_fromCDhit_outHASH->{ \$inClusterArray->[$i]->{'0_1_2_pep_____ID'} }=\$ProteinID_KEY_HASH_fromCDhit_outHASH->{ $inClusterArray->[$i]->{'0_1_2_pep_____ID'} }=\$inClusterArray->[$i]->{'0_1_0_outCluster_NB'}=$inClusterArray->[$i]->{'0_1_0_outCluster_NB'} \n\n" );
	  }
	  
	}
	#DirFileHandle::PrintAndWarnDumper ($ProteinID_KEY_HASH_fromCDhit_outHASH, "\$ProteinID_KEY_HASH_fromCDhit_outHASH=\n");
	return $ProteinID_KEY_HASH_fromCDhit_outHASH;
	
}


sub Build_mutiple_alignment_based_on_CDhit{  #      CDhitPraser::Build_mutiple_alignment_based_on_CDhit( $CdHitOutHashFile, $InputSelProFile_IndexTxt, $input_fasta_Hash, $CDhit_out_WorkingDIR );
	my ( $CdHitOutHashFile, $InputSelProFile_IndexTxt, $input_fasta_Hash, $CDhit_out_WorkingDIR )=@_;	
	
	my ($warnMsgBody, $warnMsgHead, $die_MsgHead, $caller_inform)=@{ DieWork::BuildWarnDieHeadInform( 'CDhitPraser', 'Build_mutiple_alignment_based_on_CDhit' ) };
		
	DieWork::Check_FileDirExist_or_DIE    ( $CdHitOutHashFile,                    "\$CdHitOutHashFile",                    $die_MsgHead, $caller_inform  );
  DieWork::Check_FileDirExist_or_DIE    ( $InputSelProFile_IndexTxt,            "\$InputSelProFile_IndexTxt",             $die_MsgHead, $caller_inform  );
  DieWork::Check_Hash_or_DIE            ( $input_fasta_Hash,                    "\$input_fasta_Hash",                     $die_MsgHead, $caller_inform  );
	DieWork::Check_DfdNoEmptString_or_DIE ( $CDhit_out_WorkingDIR,                "\$CDhit_out_WorkingDIR",                 $die_MsgHead, $caller_inform  );
  
  
  my $CdHitOutHash=Storable::retrieve($CdHitOutHashFile);	
	
	my $out_cluster_HASH=CDhitPraser::Build_mutiple_alignment_based_on_CDhit_HASH( $CdHitOutHash, $InputSelProFile_IndexTxt, $input_fasta_Hash, $CDhit_out_WorkingDIR );
	
	return $out_cluster_HASH;
}


sub Build_mutiple_alignment_based_on_CDhit_HASH{   #   CDhitPraser::Build_mutiple_alignment_based_on_CDhit_HASH( $CdHitOutHashFile, $InputSelProFile_IndexTxt, $input_fasta_Hash, $CDhit_out_WorkingDIR );
  my ( $CdHitOutHash, $InputSelProFile_IndexTxt, $input_fasta_Hash, $CDhit_out_WorkingDIR )=@_;
	
	my ($warnMsgBody, $warnMsgHead, $die_MsgHead, $caller_inform)=@{ DieWork::BuildWarnDieHeadInform( 'CDhitPraser', 'Build_mutiple_alignment_based_on_CDhit_HASH' ) };
	
	DieWork::Check_Hash_or_DIE            ( $CdHitOutHash->{'0_1_0_clusterHASH'}, "\$CdHitOutHash->{'0_1_0_clusterHASH'}",  $die_MsgHead, $caller_inform  );
	DieWork::Check_FileDirExist_or_DIE    ( $InputSelProFile_IndexTxt,            "\$InputSelProFile_IndexTxt",             $die_MsgHead, $caller_inform  );
  DieWork::Check_Hash_or_DIE            ( $input_fasta_Hash,                    "\$input_fasta_Hash",                     $die_MsgHead, $caller_inform  );
	DieWork::Check_DfdNoEmptString_or_DIE ( $CDhit_out_WorkingDIR,                "\$CDhit_out_WorkingDIR",                 $die_MsgHead, $caller_inform  );
  
  
	
	my $clustalw_or_not=1;   #   $clustalw_or_not=1 do running the clustaw, 0 or nothing don't run
	
	my $mkDir_wkingDirPath="mkdir -p $CDhit_out_WorkingDIR";  
	DieWork::Print_and_warn( "\n$warnMsgBody$mkDir_wkingDirPath\n\n");
	system ( "$mkDir_wkingDirPath");
	
	my $Basename_CDhit_out_WorkingDIR=File::Basename::basename $CDhit_out_WorkingDIR;
	my $dirname__CDhit_out_WorkingDIR=File::Basename::dirname  $CDhit_out_WorkingDIR;
	
	my $total_IPS_hash;
	
	my $out_cluster_HASH;
	foreach my $cluster_nb (    sort { $a <=> $b } (   keys (  %{ $CdHitOutHash->{'0_1_0_clusterHASH'} } )   )    ){  
	  
	  my $each_cluster_Path=$CDhit_out_WorkingDIR."/".$cluster_nb;
	  my $HyperLink_cluPath=$Basename_CDhit_out_WorkingDIR."/".$cluster_nb;
	   
	  my $mkDir_wkingDir_each_cluster_Path="mkdir -p $each_cluster_Path";  
	  DieWork::Print_and_warn( "\n$warnMsgBody$mkDir_wkingDir_each_cluster_Path\n\n");
	  system ( "$mkDir_wkingDir_each_cluster_Path");
	  
	  my $total_pep_in_cluster=$CdHitOutHash->{'0_1_0_clusterHASH'}->{$cluster_nb}->{'1_0_0_total_pep_in_cluster'};
	  $out_cluster_HASH->{$cluster_nb}->{'000_totalPep_incluster'}=$total_pep_in_cluster;
	  
	  if  (   (  defined ( $CdHitOutHash->{'0_1_0_clusterHASH'}->{$cluster_nb}->{'1_1_0_inClusterArray'} )  ) && (  ref ( $CdHitOutHash->{'0_1_0_clusterHASH'}->{$cluster_nb}->{'1_1_0_inClusterArray'} ) eq 'ARRAY'  )   ){}
    else{		DieWork::Just_dieWork( $die_MsgHead."\n \$CdHitOutHash->{'0_1_0_clusterHASH'}->{$cluster_nb}->{'1_1_0_inClusterArray'}=$CdHitOutHash->{'0_1_0_clusterHASH'}->{$cluster_nb}->{'1_1_0_inClusterArray'} should be a ARRAY  !!  $!\n\n\n".$caller_inform ); 	}
	  
	  my $inClusterArray=$CdHitOutHash->{'0_1_0_clusterHASH'}->{$cluster_nb}->{'1_1_0_inClusterArray'};
	  
	  my $each_Cluster_Fasta_string='';
	  for (  my $i=0; $i < @{ $inClusterArray } ; $i++  ){
	  	my $pep_____ID=$inClusterArray->[$i]->{'0_1_2_pep_____ID'};
	  	my $seq_string=FastaFileHandle::feach_seqString_from_idx_File( $pep_____ID, $InputSelProFile_IndexTxt );
	  	$each_Cluster_Fasta_string.=">$pep_____ID\n$seq_string\n\n";
	  }
	  
	  
	  my $each_cluster_fasta_txt=$each_cluster_Path."/0_0_0_cluster_fasta.txt";    #$out_cluster_HASH->{$cluster_nb}->{'001_fastaFile'}=$each_cluster_fasta_txt;
	  my $each_cluster_fasta_idx=$each_cluster_Path."/0_0_1_cluster_fasta.idx";    #$out_cluster_HASH->{$cluster_nb}->{'002_fastaIndx'}=$each_cluster_fasta_idx;
	  FastaFileHandle::BuildFastaFile_withFastaString ($each_cluster_fasta_txt, $each_Cluster_Fasta_string);  #
    FastaFileHandle::BuildIdxFile_for_fastaFile     ($each_cluster_fasta_txt, $each_cluster_fasta_idx   );  #
    
    my $cluster_TEMP_fasta_fas=$each_cluster_Path."/0_0_2_clustTEMPfast.fas";    
    my $cluster_TEMP_fasta_msf=$each_cluster_Path."/0_0_3_clustTEMPfast.msf";   
    
    my $each_cluster_fasta_msf=$each_cluster_Path."/0_0_4_cluster_fasta.msf";        
    my $Hplk_cluster_fasta_msf=$HyperLink_cluPath."/0_0_4_cluster_fasta.msf";
    
    ClustalwRun::RunMuscle_and_obtain_varous_Format_alnFILE  ($each_cluster_fasta_txt, $cluster_TEMP_fasta_fas, $cluster_TEMP_fasta_msf, $each_cluster_fasta_msf);
    
    #my $inPathHead="/test";
    #my $msf_hyperLinkFileAddrs=ExcelHandle::HylikFileAdd($Hplk_cluster_fasta_msf, $inPathHead, 'msf');  
    
    #$Hplk_cluster_fasta_msf=ExcelHandle::Lnx2Dos    ( $Hplk_cluster_fasta_msf );
    #$Hplk_cluster_fasta_msf=ExcelHandle::makeHyperLk( $Hplk_cluster_fasta_msf,  'Msf_file' ); 
    $Hplk_cluster_fasta_msf=ExcelHandle::HylikFileAddress_easyVersion($Hplk_cluster_fasta_msf, 'Msf_file'); 
    
    $out_cluster_HASH->{$cluster_nb}->{'003_msfFile'}=$Hplk_cluster_fasta_msf;
    
    Interproscan::Doing_all_the_Ipsc_Search_MutipleThread_fromFASTAfile($each_cluster_fasta_txt, 10); #This step should be used for the total fasta file!
    
    my $path_cluster_IPSC_HASH=$each_cluster_Path."/0_0_4_inter_ProScan.hsh";
    my $Hplk_cluster_IPSC_HASH=$HyperLink_cluPath."/0_0_4_inter_ProScan.hsh";  #$Hplk_cluster_fasta_msf=ExcelHandle::Lnx2Dos( $Hplk_cluster_IPSC_HASH ); $Hplk_cluster_IPSC_HASH=ExcelHandle::makeHyperLk( $Hplk_cluster_IPSC_HASH,  'Msf_file' );
    $Hplk_cluster_IPSC_HASH=ExcelHandle::HylikFileAddress_easyVersion($Hplk_cluster_IPSC_HASH, 'IpSCAN.hash', 1);  
	  $out_cluster_HASH->{$cluster_nb}->{'004_IPShash'}=$Hplk_cluster_IPSC_HASH; 
	  
	  my $CInterProScan_out_HASH=Interproscan::Eatract_interProscanHash_from_JL_ipsDatabase_inputFastaFile ($each_cluster_fasta_txt);
	  DirFileHandle::PrintDumper_plReadFile($path_cluster_IPSC_HASH, $CInterProScan_out_HASH) if ( ref($CInterProScan_out_HASH) eq 'HASH' ) ; 
  
    $total_IPS_hash=ArrayHashChange::PushSmallHash_into_bigHash($total_IPS_hash, $CInterProScan_out_HASH);
	  

	  my $path_clusterCluIPs_PNG=$each_cluster_Path."/0_0_5_clust_msf_IPS.png";
    my $Hplk_clusterCluIPs_PNG=$HyperLink_cluPath."/0_0_5_clust_msf_IPS.png";  #$Hplk_cluster_fasta_msf=ExcelHandle::Lnx2Dos( $Hplk_cluster_IPSC_HASH ); $Hplk_cluster_IPSC_HASH=ExcelHandle::makeHyperLk( $Hplk_cluster_IPSC_HASH,  'Msf_file' );
    $Hplk_clusterCluIPs_PNG=ExcelHandle::HylikFileAddress_easyVersion($Hplk_clusterCluIPs_PNG, 'Msf_PNG'); 
	  $out_cluster_HASH->{$cluster_nb}->{'006_IPSCluPNG'}=$Hplk_clusterCluIPs_PNG; 
	  
	  ClustalwRun::PrintPngForClustalw  ($each_cluster_fasta_msf, $CInterProScan_out_HASH, $path_clusterCluIPs_PNG);  #, $msf_new_showNAMEhash);
	  
	  my $dirname__CDhit_out_WorkingDIR=File::Basename::dirname  $CDhit_out_WorkingDIR;

	  
	  
	}
	
	my $excel_file_path_basename  =$Basename_CDhit_out_WorkingDIR."010.ips.xlsx";
	my $excel_file_path           =$dirname__CDhit_out_WorkingDIR."/".$excel_file_path_basename;
	my $middle_hperLkCLUSTEDrstDIR=$Basename_CDhit_out_WorkingDIR."011.cluDIR";
	my $middle_clustalRSTdir      =$dirname__CDhit_out_WorkingDIR."/".$middle_hperLkCLUSTEDrstDIR;
	
	Interproscan_MultiplAnalysis::Build_excel_for_mutipl_IPShash(   
	                                                                $total_IPS_hash,                #   Total hash holding all interproscan prased out information                   
	                                                                $input_fasta_Hash,              #   fasta hash holding fasta hash
	                                                                $excel_file_path,               #   output excel file path
	                                                                '',                             #   sheet name 1 , no need to set if you can use the default name sh1_no_order
	                                                                '',                             #   sheet name 2 , no need to set if you can use the default name sh2_in_order 
	                                                                $clustalw_or_not,               #   $clustalw_or_not=1 do running the clustaw, 0 or nothing don't run
	                                                                $middle_clustalRSTdir,          #   clustalw directory holding clustalw results
	                                                                $middle_hperLkCLUSTEDrstDIR     #   hypling clustalw directory path 
	                                                            );
	                                                                
	                                                                
	
	return $out_cluster_HASH;
}


1;