#!/usr/bin/perl -w
BEGIN{  push (@INC, '/home/fredjiang/fredProgs/PMs')  };

use strict;
use warnings;
use HTML::TreeBuilder;
use XML::Simple;
use Data::Dumper;
$Data::Dumper::Purity=1;
use Storable ;
use Storable qw(nstore store_fd nstore_fd freeze thaw dclone);
use Excel::Writer::XLSX;
use Excel::Writer::XLSX::Utility;
use Bio::Seq;
use Bio::SeqIO;
use File::Basename;

use DirFileHandle;
use TimeWork;
use ExcelHandle;
use ForeachHash;
use ClustalwRun;


#sub DumpXml{          #输入: 1，“xml文件”，                           2，“pathway信息对应的需要生成的 hash Dumperfile 的名字”; 
                       #返回: 1，“Pathway 的hash的引用”                                                                                                                              ##  已被更改的返回 #############1，“xml的hash的引用”，   
                       #写入: 1，是“hash dumperfile”，2是 “hash dumperfile的 可读文件”
                       
                       
                       
#sub InterProScanRun{  #$PepSeqFastaFile, $pepInterProScan_Result, $hperLk_pepIPS_rstHtmlDir
                       #输入: 1，包含需要运行InterProScan的Fasta文件， 2，InterProScan的结果文件名，以此生成大批文件和文件夹， 3，包含输出Html的文件夹名， 4，是否真的要进行InterProScan，“1 ”是yes，其它是No
                       #返回: 无 
                       #输出: 一系列输出文件
                       
                       
#sub praseInterProScanHtml{  #解析 interproscan输出的html文件 
                             #输入: 1，interproscan输出的html文件名，
                             #输出：1，解析后的Hash
                             #写入：1，是解析后的Hash的“hash dumperfile”，2,是 解析后的Hash的“hash dumperfile的 可读文件” 






#sub SeqIn_PrasedHashOut{  #($inSeqFile, $outHashFile, $runOrNot, $outFilesDir, $delOrnot)                         
                           #输入: 1，用于interproscan分析的 fasta文件， 2，输出的Hash文件的名字，用来装本程序的输出hash  3，是否需要interproscan的运行？“1”，运行   
                           #      4，输出文件的路径，如没有填写则按输入文件构建输出文件夹   5，是否要生成临时文件并在最后删除临时文件？“1”，则生成并删除
                           #输出：1，解析后的Hash
                           #写入：1，是解析后的Hash的“hash dumperfile”，2,是 解析后的Hash的“hash dumperfile的 可读文件” 
                               
                      
package Interproscan;


my $JL_INterProScan_database         ="/home/fredjiang/EightT/fredjiang.2018.04.02/INterProScan_database20181217";
my $JL_INterProScan_tempWorkDIR=$JL_INterProScan_database."/0_0_TempWorkDIR";
my $JL_INterProScan_RealDoneDIR=$JL_INterProScan_database."/0_1_realDonkDIR";
my $JL_INterProScan_DoneWorkHsh=$JL_INterProScan_database."/0_2_realDone.Hsh";
my $JL_INterProScan_upGradeMark=$JL_INterProScan_database."/0_3_upGradeING__";
my $JL_SubDIrSize=1000;
my $JL_interProscanShowName='QUERY';
my $JL_MutipleThread_doing_Search_pl="/home/fredjiang/fredProgs/PLs/InterProscan_tempWorkDoing.pl";   


sub Doing_all_the_Ipsc_Search_MutipleThread_fromFASTAfile{  #        Interproscan::Doing_all_the_Ipsc_Search_MutipleThread_fromFASTAfile($inFastaFile, $coreToUse);
	
	my ($inFastaFile, $coreToUse)=@_;
	
	my $warnMsgBody="\nIn package  Interproscan,\tIn sub Doing_all_the_Ipsc_Search_MutipleThread_fromFASTAfile,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody"; #warn "$warnMsgHead\$inFastaFile=$inFastaFile\n";
  #my $subCallereIfm=DirFileHandle::print_SubCallerInform; #warn $subCallereIfm;
  
  my $workingHASH;
  my $InFastaSeHASH=ClustalwRun::NewReadInSeq_into_A_HASH($inFastaFile, 'fasta');
	if (   (  defined ( $InFastaSeHASH )  ) && (  ref ( $InFastaSeHASH ) eq 'HASH'  ) && (  defined ( $InFastaSeHASH->{'2_seqKeyHash'} )  ) && (  ref ( $InFastaSeHASH->{'2_seqKeyHash'} ) eq 'HASH'  )   ) {
		foreach my $eachSeq (    sort { $a cmp $b } (   keys (  %{ $InFastaSeHASH->{'2_seqKeyHash'} }  )   )    ){ 
		  $workingHASH->{$eachSeq}=$InFastaSeHASH->{'2_seqKeyHash'}->{$eachSeq};
		}
	}
  
  if (   (  defined ( $workingHASH )  ) && (  ref ( $workingHASH ) eq 'HASH'  )   ){
  	Interproscan::Doing_all_the_Ipsc_Search_MutipleThread($workingHASH, $coreToUse);
  }
  
}


sub Doing_all_the_Ipsc_Search_MutipleThread{  #   Interproscan::Doing_all_the_Ipsc_Search_MutipleThread($inHash, $coreToUse);
	my ($inHash, $coreToUse)=@_;
	
	my $warnMsgBody="\nIn package  Interproscan,\tIn sub Doing_all_the_Ipsc_Search_MutipleThread,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	#my $subCallereIfm=DirFileHandle::print_SubCallerInform;
	
	if  (   (  defined ( $inHash )  ) && (  ref ( $inHash ) eq 'HASH'  )   ){
		my $inHashTempFile_basename=TimeWork::GetTimeDirOrFileName;
		$inHashTempFile_basename="IN_searchHASH_".$inHashTempFile_basename;
		my $inHashTempFile_wholpath=$JL_INterProScan_tempWorkDIR."/".$inHashTempFile_basename;
		DirFileHandle::PrintDumper($inHashTempFile_wholpath, $inHash) if (   (  defined ( $inHash )  ) && (  ref ( $inHash ) eq 'HASH'  )   );
		
		if (  -e ( $inHashTempFile_wholpath )  ){
			system ( "perl $JL_MutipleThread_doing_Search_pl  -i $inHashTempFile_wholpath -c $coreToUse");
			
		}
		
		if (  -e ( $inHashTempFile_wholpath )  ){ system ("rm -f $inHashTempFile_wholpath" )	;	}
		my $readFile=&GetReadTxt($inHashTempFile_wholpath);
		if (  -e ( $readFile                )  ){ system ("rm -f $readFile"                )	;	}
		
	}
	
	
	
}

sub GetReadTxt{
	my ($inF)=@_;
	my $ouF=$inF.".read.txt";
	return $ouF;
}



sub Filter_out_searchedSeq{   #      Interproscan::Filter_out_searchedSeq($whole_inSeq_HASH);
	my ($whole_inSeq_HASH)=@_;
	
	my $warnMsgBody="\nIn package  Interproscan,\tIn sub Filter_out_searchedSeq,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	#my $subCallereIfm=DirFileHandle::print_SubCallerInform;
	
	my $new_HASH;  my $filger_or_not=0;
	if (   (  defined ( $JL_INterProScan_DoneWorkHsh )  ) && ( $JL_INterProScan_DoneWorkHsh=~m/\S+/ ) && (  -e ( $JL_INterProScan_DoneWorkHsh )  )   ) {
	  my $JL_INterProScan_DoneWorkHsh_HASH=Storable::retrieve  ( $JL_INterProScan_DoneWorkHsh );
	  if  (   (  defined ( $JL_INterProScan_DoneWorkHsh_HASH )  ) && (  ref ( $JL_INterProScan_DoneWorkHsh_HASH ) eq 'HASH'  ) && (  defined ( $JL_INterProScan_DoneWorkHsh_HASH->{'0_1_0_SeqKeyHASH'} )  ) && (  ref ( $JL_INterProScan_DoneWorkHsh_HASH->{'0_1_0_SeqKeyHASH'} ) eq 'HASH'  )   ) {
	    
	    if (   (  defined ( $whole_inSeq_HASH )  ) && (  ref ( $whole_inSeq_HASH ) eq 'HASH'  )   ){
	    	foreach my $eachSeq (    sort { $a cmp $b } (   keys (  %{ $whole_inSeq_HASH }  )   )    ){
	        my $newString=Interproscan::MakeGoodString_forInterProscan ($eachSeq);
	        if (   (  defined ( $JL_INterProScan_DoneWorkHsh_HASH->{'0_1_0_SeqKeyHASH'}->{$newString} )  ) && ( $JL_INterProScan_DoneWorkHsh_HASH->{'0_1_0_SeqKeyHASH'}->{$newString}=~m/\S+/ )   ){
	        	$filger_or_not=1;
	        }
	        else {
	        	$new_HASH->{$eachSeq}=$whole_inSeq_HASH->{$eachSeq};
	        }
	      }
	    }
	  }
	}	
	if ( $filger_or_not==0 ){
		$new_HASH=$whole_inSeq_HASH;
	}
	return $new_HASH;
}


sub Check_inSeq_done_orNot_in_JL_interProscan_database{   #      Interproscan::Check_inSeq_done_orNot_in_JL_interProscan_database($inString);
	my ($inString)=@_;
	
	my $warnMsgBody="\nIn package  Interproscan,\tIn sub Check_inSeq_done_orNot_in_JL_interProscan_database,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	#my $subCallereIfm=DirFileHandle::print_SubCallerInform;
  
  my $found=0;
  if (   (  defined ( $JL_INterProScan_DoneWorkHsh )  ) && ( $JL_INterProScan_DoneWorkHsh=~m/\S+/ ) && (  -e ( $JL_INterProScan_DoneWorkHsh )  )   ) {
	  my $JL_INterProScan_DoneWorkHsh_HASH=Storable::retrieve  ( $JL_INterProScan_DoneWorkHsh );
	  if  (   (  defined ( $JL_INterProScan_DoneWorkHsh_HASH )  ) && (  ref ( $JL_INterProScan_DoneWorkHsh_HASH ) eq 'HASH'  ) && (  defined ( $JL_INterProScan_DoneWorkHsh_HASH->{'0_1_0_SeqKeyHASH'} )  ) && (  ref ( $JL_INterProScan_DoneWorkHsh_HASH->{'0_1_0_SeqKeyHASH'} ) eq 'HASH'  )   ) {
	    $inString=~s/\s+//g;  $inString=Interproscan::MakeGoodString_forInterProscan ($inString);
	    if (   (  defined ( $JL_INterProScan_DoneWorkHsh_HASH->{'0_1_0_SeqKeyHASH'}->{$inString} )  ) && ( $JL_INterProScan_DoneWorkHsh_HASH->{'0_1_0_SeqKeyHASH'}->{$inString}=~m/\S+/ )   ){
	    	$found=1;
	    }
	  }
	}	
	return $found;
}

sub DOing_InterProscanWork_in_temp_DIR_fromFASTAfile{  #     Interproscan::DOing_InterProscanWork_in_temp_DIR_fromFASTAfile($inFastaFile);
	my ($inFastaFile)=@_;
	
	
	my $warnMsgBody="\nIn package  Interproscan,\tIn sub DOing_InterProscanWork_in_temp_DIR_fromFASTAfile,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody"; #warn "$warnMsgHead\$inFastaFile=$inFastaFile\n";
  #my $subCallereIfm=DirFileHandle::print_SubCallerInform; #warn $subCallereIfm;
  
  my $workingHASH;
  my $InFastaSeHASH=ClustalwRun::NewReadInSeq_into_A_HASH($inFastaFile, 'fasta');
	if (   (  defined ( $InFastaSeHASH )  ) && (  ref ( $InFastaSeHASH ) eq 'HASH'  ) && (  defined ( $InFastaSeHASH->{'2_seqKeyHash'} )  ) && (  ref ( $InFastaSeHASH->{'2_seqKeyHash'} ) eq 'HASH'  )   ) {
		foreach my $eachSeq (    sort { $a cmp $b } (   keys (  %{ $InFastaSeHASH->{'2_seqKeyHash'} }  )   )    ){ 
		  $workingHASH->{$eachSeq}=$InFastaSeHASH->{'2_seqKeyHash'}->{$eachSeq};
		}
	}
  
  if (   (  defined ( $workingHASH )  ) && (  ref ( $workingHASH ) eq 'HASH'  )   ){
  	Interproscan::DOing_InterProscanWork_in_temp_DIR_fromHASH($workingHASH);
  }
  
}

sub DOing_InterProscanWork_in_temp_DIR_fromHASH{  #     Interproscan::DOing_InterProscanWork_in_temp_DIR_fromHASH($workingHASH);
	my ($workingHASH)=@_;
	
	
	my $warnMsgBody="\nIn package  Interproscan,\tIn sub DOing_InterProscanWork_in_temp_DIR_fromHASH,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody"; #warn "$warnMsgHead\$inFastaFile=$inFastaFile\n";
  #my $subCallereIfm=DirFileHandle::print_SubCallerInform; #warn $subCallereIfm;
  
  my $JL_INterProScan_DoneWorkHsh_HASH;
  if (   (  defined ( $JL_INterProScan_DoneWorkHsh )  ) && ( $JL_INterProScan_DoneWorkHsh=~m/\S+/ ) && (  -e ( $JL_INterProScan_DoneWorkHsh )  )   ) {
	  $JL_INterProScan_DoneWorkHsh_HASH=Storable::retrieve  ( $JL_INterProScan_DoneWorkHsh );
  }
	 
  
	
	#my $InFastaSeHASH=ClustalwRun::NewReadInSeq_into_A_HASH($inFastaFile, 'fasta');
	my $InFastaSeHASH=Storable::dclone ( $workingHASH );
	if (   (  defined ( $InFastaSeHASH )  ) && (  ref ( $InFastaSeHASH ) eq 'HASH'  )   ) {
		my $TempReportHash;
		my $inSeqNb=1;
	  foreach my $eachSeq (    sort { $a cmp $b } (   keys (  %{ $InFastaSeHASH }  )   )    ){ 
	  	my $realEachSeq= Interproscan::MakeGoodString_forInterProscan ($eachSeq);
	  	
	  	my $CheckPrsent=0;
	  	if  (      (  defined ( $JL_INterProScan_DoneWorkHsh_HASH )  ) && (  ref ( $JL_INterProScan_DoneWorkHsh_HASH ) eq 'HASH'  ) 
	  	        && (  defined ( $JL_INterProScan_DoneWorkHsh_HASH->{'0_1_0_SeqKeyHASH'} )  ) && (  ref ( $JL_INterProScan_DoneWorkHsh_HASH->{'0_1_0_SeqKeyHASH'} ) eq 'HASH'  )   
	  	        && (  defined ( $JL_INterProScan_DoneWorkHsh_HASH->{'0_1_0_SeqKeyHASH'}->{$realEachSeq} )  ) && (  $JL_INterProScan_DoneWorkHsh_HASH->{'0_1_0_SeqKeyHASH'}->{$realEachSeq}=~m/\S+/  )
	  	    )
	    {
	      $CheckPrsent=1;
	    }	  
	 	
	  	
	 	
	  	if ( $CheckPrsent ==0 ){
	  		
	  		my $BigOutFilesDir=TimeWork::GetTimeDirOrFileName;
	  		 
	  		$BigOutFilesDir         =$TempReportHash->{$realEachSeq}->{'0_0_0_BigDIR_baseName'}=$BigOutFilesDir."_".$inSeqNb;
	  		my $wholepth            =$TempReportHash->{$realEachSeq}->{'0_0_1_BigDIR_wholPath'}=$JL_INterProScan_tempWorkDIR."/".$BigOutFilesDir;   system ("mkdir -p $wholepth");
	  		
	  		my $inScanFastaFile     =$TempReportHash->{$realEachSeq}->{'0_0_2_inFast_baseName'}='0_0_inPut_Fasta.txt';
	  		my $inScanFastaFile_Path=$TempReportHash->{$realEachSeq}->{'0_0_3_inFast_wholPath'}=$wholepth."/".$inScanFastaFile;  
	  		my $fastaName           =$TempReportHash->{$realEachSeq}->{'0_0_4_inFast_ShowName'}=$JL_interProscanShowName;
	  		
	  		Interproscan::BuildOneSeqFastaFile_forInterProscan ($realEachSeq, $inScanFastaFile_Path, $fastaName);
	  		
	  		##############################################################
	  		
	  		my $outFilesDir         =$TempReportHash->{$realEachSeq}->{'0_0_5_OtDIR_baseName'}='0_1_output_DIR';
	  		my $outFilesDir_path    =$TempReportHash->{$realEachSeq}->{'0_0_6_OtDIR_wholPath'}=$wholepth."/".$outFilesDir;                       system ("mkdir -p $outFilesDir_path");
	  		
	  		my $outRstFile          =$TempReportHash->{$realEachSeq}->{'0_0_7_OtRst_wholPath'}=$outFilesDir_path."/1_0_Rst";
        my $outRstHtml_baseName =$TempReportHash->{$realEachSeq}->{'0_0_8_0_OtHtm_BaseNm'}="1_1_Htm";         
        my $outRstHtml          =$TempReportHash->{$realEachSeq}->{'0_0_8_OtHtmDIRwhPath'}=$outFilesDir_path."/".$outRstHtml_baseName;       system ("mkdir -p $outRstHtml");
        my $outRstXml           =$TempReportHash->{$realEachSeq}->{'0_0_9_OtXml_wholPath'}=$outFilesDir_path."/1_0_Rst.xml";
  
    	  			  			  			  			  			  			  			  			  		    print "20181218-2011-1 \$inScanFastaFile_Path=$inScanFastaFile_Path\n";
        &InterProScanRun( $inScanFastaFile_Path, $outRstFile, $outRstHtml);     print "20181218-2011-2 \$inScanFastaFile_Path=$inScanFastaFile_Path\n";
        
        my $outRstXmlHASH=&DumpXml($outRstXml);
        my $allHtmlDirFiles=DirFileHandle::getDirArray($outRstHtml);
        
        my $allHtmlHash;
        foreach my $eachHtmlFile (  @{ $allHtmlDirFiles }  ){    #print "\$eachHtmlFile=$eachHtmlFile\n";
          if ($eachHtmlFile=~m/.html$/){                          print "20181218-2011 \$eachHtmlFile=$eachHtmlFile\n";
          	my $outHash=&praseInterProScanHtml( "$outRstHtml/$eachHtmlFile" );
          	foreach my $k0 (    sort { $a cmp $b } (   keys (  %{ $outHash } )   )    ){
          	  $allHtmlHash->{$k0}=$outHash->{$k0};
          	  if (  defined( $outRstXmlHASH->{$k0} )  ){
          	  	$allHtmlHash->{$k0}->{'Pathway'}=$outRstXmlHASH->{$k0};
          	  }
          	}
          }
        }
        
        my $outHashFile         =$TempReportHash->{$realEachSeq}->{'0_0_a_OtHsh_baseName'}='0_2_outHash.hsh';                                                             
        my $outHashFile_path    =$TempReportHash->{$realEachSeq}->{'0_0_b_OtHsh_wholPath'}=$wholepth."/".$outHashFile;        
  
        DirFileHandle::PrintDumper($outHashFile_path    , $allHtmlHash) if (   (  defined ( $allHtmlHash )  ) && (  ref ( $allHtmlHash ) eq 'HASH'  )   );
	  		
	  		my $SiglDonHSH;
	  		$SiglDonHSH->{$realEachSeq}=Storable::dclone ( $TempReportHash->{$realEachSeq} );
	  		my $inDir_infomrHashFile=$TempReportHash->{$realEachSeq}->{'0_0_d_SiglDonHSHPath'}=$wholepth."/".'9_0_SiglDon.HSH';     
	  	  DirFileHandle::PrintDumper($inDir_infomrHashFile, $allHtmlHash) if (   (  defined ( $inDir_infomrHashFile )  ) && (  ref ( $inDir_infomrHashFile ) eq 'HASH'  )   ); 
	  		##############################################################
	  		                         $TempReportHash->{$realEachSeq}->{'0_0_e_DoneAndRecPath'}=$wholepth."/".'9_1_DoAdRec.txt';
	  		
	  	}
	  	$inSeqNb++;
	  }
	  
	  if (   (  defined ( $TempReportHash )  ) && (  ref ( $TempReportHash ) eq 'HASH'  )    ){
	  	
	  	
	  	
	  	
	  	my $timeDIRWord=TimeWork::GetTimeDirOrFileName;
	  	
	  	my $tempWorkHashFile=$JL_INterProScan_tempWorkDIR."/TEMP_work_Mark_".$timeDIRWord.".HsH";
	  	
	  	
	  	
	  	DirFileHandle::PrintDumper($tempWorkHashFile,$TempReportHash) if  (   (  defined ( $TempReportHash )  ) && (  ref ( $TempReportHash ) eq 'HASH'  )   );
	  	
	  	foreach my $realEachSeq (    sort { $a cmp $b } (   keys (  %{ $TempReportHash }  )   )    ){ 
	  		if (   (  defined ( $TempReportHash->{$realEachSeq}->{'0_0_e_DoneAndRecPath'} )  ) && (  $TempReportHash->{$realEachSeq}->{'0_0_e_DoneAndRecPath'}=~m/\S+/  )    ){
	  		  my $DoneAndRecPath=$TempReportHash->{$realEachSeq}->{'0_0_e_DoneAndRecPath'};
	  		  open ( DONERCD, ">$DoneAndRecPath") or &Just_dieWork(  "cannot create \$DoneAndRecPath=$DoneAndRecPath : $!\n\n\n"  );
	  		  print DONERCD "This work is done and recorded into the $tempWorkHashFile !!\n\n";
	  		  close (DONERCD);
	  		} 
	  		
	  	}
	  	
	  	
	  	
	  	
	  }
	  
	}
	
}


sub Just_dieWork{
	my ( $dieWord )=@_;
	
	print $dieWord;
	die   $dieWord;
}


sub dieWork{
	my ( $dieWord )=@_;
	Interproscan::Del_JL_INterProScan_upGradeMark();
	print $dieWord;
	die   $dieWord;
}

sub Del_JL_INterProScan_upGradeMark{   #Interproscan::Del_JL_INterProScan_upGradeMark();
	if (  -e( $JL_INterProScan_upGradeMark )  ){
	  system ( "rm -f $JL_INterProScan_upGradeMark" );
	}	
}


#Interproscan::Find_Finished_ophan_tempDIR_toBuildTEMPWORKHASH_and_doingtheCleanWork();
sub Find_Finished_ophan_tempDIR_toBuildTEMPWORKHASH_and_doingtheCleanWork{
	my $warnMsgBody="\nIn package  Interproscan,\tIn sub Find_Finished_ophan_tempDIR_toBuildTEMPWORKHASH_and_doingtheCleanWork,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $subCallereIfm=DirFileHandle::print_SubCallerInform;
	
	my $allTemmDIRfile=DirFileHandle::getDirArray($JL_INterProScan_tempWorkDIR);
  
  my $isThere_anything_new=0;
  
  my $TEMP_work_Mark_HASH;
  if (   (  defined ( $allTemmDIRfile )  ) && (  ref ( $allTemmDIRfile ) eq 'ARRAY'  )    ){
    foreach my $each_TEMP_DIR (  @{ $allTemmDIRfile }  ){ 
	    if ( $each_TEMP_DIR=~/^D\d+T\d+P\d+_\d+$/ ){
	    	my $each_TEMP_DIR_path=$JL_INterProScan_tempWorkDIR."/".$each_TEMP_DIR;
	      my $SiglDonHSHFile=$each_TEMP_DIR_path."/9_0_SiglDon.HSH";
	      
	      
	      if (  -e ( $SiglDonHSHFile )  ){
	      	my $SiglDonHSH = Storable::retrieve ( $SiglDonHSHFile );
	      	if (   (  defined ( $SiglDonHSH )  ) && (  ref ( $SiglDonHSH ) eq 'HASH'  )    ){
	      		my $DoneAndRecFile=$each_TEMP_DIR_path."/9_1_DoAdRec.txt";   #0_0_e_DoneAndRecPath'}=$wholepth."/".'9_1_DoAdRec.txt
	      		if (  -e ( $SiglDonHSHFile )  ){
	      			
	      		}
	      		else{
	      			foreach my $eachSeq (    sort { $a cmp $b } (   keys (  %{ $SiglDonHSH } )   )    ){
	      		    $TEMP_work_Mark_HASH->{$eachSeq}=Storable::dclone ( $SiglDonHSH->{$eachSeq} );
                $TEMP_work_Mark_HASH->{$eachSeq}->{'0_0_e_DoneAndRecPath'}=$each_TEMP_DIR_path."/".'9_1_DoAdRec.txt';
	      		    $isThere_anything_new=1;
	      		  }	      		
	      		}    		
	      		
	      	}	      	
	      	else{	      	  
	      		if (  -d ( $each_TEMP_DIR_path )  ){ system ( "rm -rf $each_TEMP_DIR_path"); }	        
	      	}
	      	
	      }
	      
	    }
	  }
	}
	if ( $isThere_anything_new==1 ){
		
		 	my $timeDIRWord=TimeWork::GetTimeDirOrFileName;
	  	
	  	my $tempWorkHashFile=$JL_INterProScan_tempWorkDIR."/TEMP_work_Mark_".$timeDIRWord.".HsH";
	  	
	  	
	  	
	  	DirFileHandle::PrintDumper($tempWorkHashFile,$TEMP_work_Mark_HASH) if  (   (  defined ( $TEMP_work_Mark_HASH )  ) && (  ref ( $TEMP_work_Mark_HASH ) eq 'HASH'  )   );
	  	
	  	foreach my $realEachSeq (    sort { $a cmp $b } (   keys (  %{ $TEMP_work_Mark_HASH }  )   )    ){ 
	  		if (   (  defined ( $TEMP_work_Mark_HASH->{$realEachSeq}->{'0_0_e_DoneAndRecPath'} )  ) && (  $TEMP_work_Mark_HASH->{$realEachSeq}->{'0_0_e_DoneAndRecPath'}=~m/\S+/  )    ){
	  		  my $DoneAndRecPath=$TEMP_work_Mark_HASH->{$realEachSeq}->{'0_0_e_DoneAndRecPath'};
	  		  open ( DONERCD, ">$DoneAndRecPath") or &Just_dieWork(  "cannot create \$DoneAndRecPath=$DoneAndRecPath : $!\n\n\n"  );
	  		  print DONERCD "This work is done and recorded into the $tempWorkHashFile !!\n\n";
	  		  close (DONERCD);
	  		} 
	  		
	  	}
	  	
		
		
	}
	    
	
}


#Interproscan::MoveEverthing_from_temp_to_DoneDIR_and_upgrade_theHASH();
sub MoveEverthing_from_temp_to_DoneDIR_and_upgrade_theHASH{  #   Interproscan::MoveEverthing_from_temp_to_DoneDIR_and_upgrade_theHASH
	
	my $warnMsgBody="\nIn package  Interproscan,\tIn sub MoveEverthing_from_temp_to_DoneDIR_and_upgrade_theHASH,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $subCallereIfm=DirFileHandle::print_SubCallerInform;
	
	if (  -e( $JL_INterProScan_upGradeMark )  ){
		my $dieMsg=" $die_MsgHead\n\nAnother upgrating pid maybe running now, please check the ps aux !!\n\n\n\n";
		print $dieMsg;  die $dieMsg;
	}
	
	
	
	open ( UPGRDMK , ">$JL_INterProScan_upGradeMark") or &dieWork ( "$die_MsgHead cannot create \$JL_INterProScan_upGradeMark=$JL_INterProScan_upGradeMark : $! \n\n\n" ); 
	my $timeNow=TimeWork::GetHumanRead_NOW_Time;
	print UPGRDMK "$timeNow\n\n\n\n$subCallereIfm\n\n\n\n"; close (UPGRDMK);
	
	
	my $SubDIrSize=$JL_SubDIrSize;
	
	
	my $JL_INterProScan_DoneWorkHsh_HASH;
	if (  -e ( $JL_INterProScan_DoneWorkHsh )  ){
		$JL_INterProScan_DoneWorkHsh_HASH=Storable::retrieve ($JL_INterProScan_DoneWorkHsh);
	}
	
	
	
	
	my $allTemmDIRfile=DirFileHandle::getDirArray($JL_INterProScan_tempWorkDIR);
  
  my $isThere_anything_new=0;
  
  my $Big_Seq_key_HASH;
  if (   (  defined ( $allTemmDIRfile )  ) && (  ref ( $allTemmDIRfile ) eq 'ARRAY'  )    ){
    foreach my $eachTEMP_Hashfile (  @{ $allTemmDIRfile }  ){ 
	    if ( $eachTEMP_Hashfile=~/^TEMP_work_Mark_\S+\.HsH$/ ){
	    	my $TEMP_FILE_wholePath=$JL_INterProScan_tempWorkDIR."/".$eachTEMP_Hashfile;
	    	if (  -e ( $TEMP_FILE_wholePath )  ){
	    	  my $TEMP_File_HASH=Storable::retrieve ($TEMP_FILE_wholePath);
	    	  if (   (  defined ( $TEMP_File_HASH )  ) && (  ref ( $TEMP_File_HASH ) eq 'HASH'  )    ){
	    	    foreach my $eachSeq (    sort { $a cmp $b } (   keys (  %{ $TEMP_File_HASH } )   )    ){
	    	    	if (      (  defined ( $JL_INterProScan_DoneWorkHsh_HASH )  ) && (  ref ( $JL_INterProScan_DoneWorkHsh_HASH ) eq 'HASH'  ) 
	    	    	       && (  defined ( $JL_INterProScan_DoneWorkHsh_HASH->{  '0_1_0_SeqKeyHASH'} )  ) && (  ref ( $JL_INterProScan_DoneWorkHsh_HASH->{'0_1_0_SeqKeyHASH'} ) eq 'HASH'  ) 
	    	    	       && (  defined ( $JL_INterProScan_DoneWorkHsh_HASH->{'0_1_0_SeqKeyHASH'}->{$eachSeq} )  )    
	    	    	   )
	    	    	{
	    	    	  #do nothing	
	    	    	}
	    	    	else{
	    	    	  
	    	    	  $isThere_anything_new=1;
	    	    	  
	    	    	  my $New_SubDir;
	    	    	  my $New_DrNumb;
	    	    	  my $New_InsubN; 	
	    	    	  my $New_totNmb;
	    	    		if (      (  defined ( $JL_INterProScan_DoneWorkHsh_HASH )  ) && (  ref ( $JL_INterProScan_DoneWorkHsh_HASH ) eq 'HASH'  ) 
	    	    		       && (  defined ( $JL_INterProScan_DoneWorkHsh_HASH->{'0_0_0_LastSubDir'} )  ) && ( $JL_INterProScan_DoneWorkHsh_HASH->{'0_0_0_LastSubDir'}=~m/\d+/ )
	    	    		       && (  defined ( $JL_INterProScan_DoneWorkHsh_HASH->{'0_0_1_LsSbDrNumb'} )  ) && ( $JL_INterProScan_DoneWorkHsh_HASH->{'0_0_1_LsSbDrNumb'}=~m/\d+/ )
	    	    		       && (  defined ( $JL_INterProScan_DoneWorkHsh_HASH->{'0_0_2_LastInsubN'} )  ) && ( $JL_INterProScan_DoneWorkHsh_HASH->{'0_0_2_LastInsubN'}=~m/\d+/ )
	    	    		       && (  defined ( $JL_INterProScan_DoneWorkHsh_HASH->{'0_0_3_total_Numb'} )  ) && ( $JL_INterProScan_DoneWorkHsh_HASH->{'0_0_3_total_Numb'}=~m/\d+/ )
	    	    		   )
	    	    	  {
	    	    		  	    	    		  
	    	    		  $New_SubDir=$JL_INterProScan_DoneWorkHsh_HASH->{'0_0_0_LastSubDir'};
	    	    		  $New_DrNumb=$JL_INterProScan_DoneWorkHsh_HASH->{'0_0_1_LsSbDrNumb'};
	    	    		  $New_InsubN=$JL_INterProScan_DoneWorkHsh_HASH->{'0_0_2_LastInsubN'}+1;
	    	    		  $New_totNmb=$JL_INterProScan_DoneWorkHsh_HASH->{'0_0_3_total_Numb'}+1;
	    	    		  
	    	    		  if ( $New_InsubN > $SubDIrSize ){
	    	    		  	$New_InsubN=1;
	    	    		  	$New_DrNumb=$New_DrNumb+1;
	    	    		  	$New_SubDir=sprintf ("%05d",$New_DrNumb);
	    	    		  }
	    	    		  
	    	    		  $JL_INterProScan_DoneWorkHsh_HASH->{'0_0_0_LastSubDir'}=$New_SubDir;
	    	    		  $JL_INterProScan_DoneWorkHsh_HASH->{'0_0_1_LsSbDrNumb'}=$New_DrNumb;
	    	    		  $JL_INterProScan_DoneWorkHsh_HASH->{'0_0_2_LastInsubN'}=$New_InsubN; 
	    	    		  $JL_INterProScan_DoneWorkHsh_HASH->{'0_0_3_total_Numb'}=$New_totNmb;  
	    	    		}
	    	    		else{
	    	    			$New_SubDir=1;     
	    	    		  $New_DrNumb=1;    $New_SubDir=sprintf ("%05d",$New_DrNumb);
	    	    		  $New_InsubN=1;
	    	    			
	    	    			$JL_INterProScan_DoneWorkHsh_HASH->{'0_0_0_LastSubDir'}=$New_SubDir;
	    	    		  $JL_INterProScan_DoneWorkHsh_HASH->{'0_0_1_LsSbDrNumb'}=$New_DrNumb;
	    	    		  $JL_INterProScan_DoneWorkHsh_HASH->{'0_0_2_LastInsubN'}=$New_InsubN;
	    	    		  $JL_INterProScan_DoneWorkHsh_HASH->{'0_0_3_total_Numb'}=1;
	    	    		}
	    	    		
	    	    		
	    	    		
	    	    		my $new_tmpHASH=Storable::dclone( $TEMP_File_HASH->{$eachSeq} );
	    	    		my $new_whole_parent_DIR =$new_tmpHASH->{'0_1_1_new_BigDIR_wholPath'}=$JL_INterProScan_RealDoneDIR."/".$New_SubDir; system ("mkdir -p $new_whole_parent_DIR");
	    	    		my $new_whole_DIR        =$new_tmpHASH->{'0_1_1_new_BigDIR_wholPath'}=$JL_INterProScan_RealDoneDIR."/".$New_SubDir."/".$TEMP_File_HASH->{$eachSeq}->{'0_0_0_BigDIR_baseName'}; 
	    	    		 
	    	    		system ( "mv -f   $TEMP_File_HASH->{$eachSeq}->{'0_0_1_BigDIR_wholPath'}  $new_whole_parent_DIR/" );
	    	    	
#	    	    		system ("mkdir -p $new_whole_DIR");   system ( "cp -rf   $TEMP_File_HASH->{$eachSeq}->{'0_0_1_BigDIR_wholPath'}  $new_whole_parent_DIR/" );
	    	    		
	    	    		                          $new_tmpHASH->{'0_1_3_new_inFast_wholPath'}=$new_whole_DIR."/".$TEMP_File_HASH->{$eachSeq}->{'0_0_2_inFast_baseName'};
	    	    		                          $new_tmpHASH->{'0_1_6_new_OtDIR_wholPath' }=$new_whole_DIR."/".$TEMP_File_HASH->{$eachSeq}->{'0_0_5_OtDIR_baseName' };  
	    	    		                          $new_tmpHASH->{'0_1_8_new_OtHtmDIRwhPath' }=$new_tmpHASH->{'0_1_6_new_OtDIR_wholPath' }."/".$TEMP_File_HASH->{$eachSeq}->{'0_0_8_0_OtHtm_BaseNm' };
	    	    		                          $new_tmpHASH->{'0_1_b_new_OtHsh_wholPath' }=$new_whole_DIR."/".$TEMP_File_HASH->{$eachSeq}->{'0_0_a_OtHsh_baseName' }; 
                
                warn "\$new_whole_DIR=$new_whole_DIR\n";
                my $newTempInformHashFile=$new_tmpHASH->{'1_0_0_new_tpHsh_wholPath' }=$new_whole_DIR."/".'1_0_tmpIfom.hsh';
                
                DirFileHandle::PrintDumper($newTempInformHashFile, $new_tmpHASH) if  (   (  defined ( $new_tmpHASH )  ) && (  ref ( $new_tmpHASH ) eq 'HASH'  )   );
                
                $JL_INterProScan_DoneWorkHsh_HASH->{'0_1_0_SeqKeyHASH'}->{$eachSeq}->{'0_0_0_OtHsh_wholPath'}=$new_tmpHASH->{'0_1_b_new_OtHsh_wholPath' };
                
                if (  -e ( $new_tmpHASH->{'0_1_b_new_OtHsh_wholPath' } )  ){
                	my $new_InterProScan_HASH=Storable::retrieve( $new_tmpHASH->{'0_1_b_new_OtHsh_wholPath' } );
                	if (      (  defined ( $new_InterProScan_HASH )  ) && (  ref ( $new_InterProScan_HASH ) eq 'HASH'  ) 
                	       && (  defined ( $new_InterProScan_HASH->{$JL_interProscanShowName} )  ) && (  ref ( $new_InterProScan_HASH->{$JL_interProscanShowName} ) eq 'HASH'  )    
                	       && (  defined ( $new_InterProScan_HASH->{$JL_interProscanShowName}->{'file'} )  ) && ( $new_InterProScan_HASH->{$JL_interProscanShowName}->{'file'}=~m/\S+/  )
                	   )
                	{
                	  my $filePath=$new_InterProScan_HASH->{$JL_interProscanShowName}->{'file'};
                	  my $oldFiPth=$new_tmpHASH->{'0_0_8_OtHtmDIRwhPath'}."/".$JL_interProscanShowName.".html";
                	  if ( $filePath eq $oldFiPth ){
                	  	my $new_filePath=$new_tmpHASH->{'0_1_8_new_OtHtmDIRwhPath' }."/".$JL_interProscanShowName.".html";
                	  	$new_InterProScan_HASH->{$JL_interProscanShowName}->{'file'}=$new_filePath;
                	  }
                	  else{
                	  	my $dieMsg4="$die_MsgHead The file \$new_tmpHASH->{'0_1_b_new_OtHsh_wholPath' }=$new_tmpHASH->{'0_1_b_new_OtHsh_wholPath' } \n\n\$filePath=$filePath \nwhould eq \n\$oldFiPth=$oldFiPth\n !!\n\n\n\n";
		                  print $dieMsg4; Interproscan::Del_JL_INterProScan_upGradeMark(); die $dieMsg4;
                	  }
                	}
                	else {
                	  my $dieMsg3="$die_MsgHead The file \$new_tmpHASH->{'0_1_b_new_OtHsh_wholPath' }=$new_tmpHASH->{'0_1_b_new_OtHsh_wholPath' } \$new_InterProScan_HASH->{$JL_interProscanShowName}->{'file'}=$new_InterProScan_HASH->{$JL_interProscanShowName}->{'file'} is not right !!\n\n\n\n";
		                print $dieMsg3; Interproscan::Del_JL_INterProScan_upGradeMark(); die $dieMsg3;
                  }
                  system  ("rm -f $new_tmpHASH->{'0_1_b_new_OtHsh_wholPath'}");
                  my $rdFl=$new_tmpHASH->{'0_1_b_new_OtHsh_wholPath'}.".read.txt"; if (  -e ( $rdFl )  ){ system  ("rm -f $rdFl"); }
                
                  DirFileHandle::PrintDumper($new_tmpHASH->{'0_1_b_new_OtHsh_wholPath' }, $new_InterProScan_HASH) if (   (  defined ( $new_InterProScan_HASH )  ) && (  ref ( $new_InterProScan_HASH ) eq 'HASH'  )   );
                
                  
                }
                else {
                	my $dieMsg2="$die_MsgHead The file \$new_tmpHASH->{'0_1_b_new_OtHsh_wholPath' }=$new_tmpHASH->{'0_1_b_new_OtHsh_wholPath' } is not right !!\n\n\n\n";
		              print $dieMsg2; Interproscan::Del_JL_INterProScan_upGradeMark(); die $dieMsg2;
                }
                
                
	    	    		
	    	    		##############################################
	    	    		
	    	    		
	    	    		
	    	    		
	    	    	}
	    	    	
	    	    	my $eachDIRpath=$TEMP_File_HASH->{$eachSeq}->{'0_0_2_inFast_baseName'};
	    	    	if (   ( -d ( $eachDIRpath ) ) && ( $eachDIRpath=~m/\d{8}\w\d{6}/ )   ){
	    	    		system  ("rm -rf $eachDIRpath");
	    	    	}
	    	    	
	    	    }
	    	  }
	    	  
	    	  system ( "rm -rf $TEMP_FILE_wholePath");
	    	  my $rdFl2=$TEMP_FILE_wholePath.".read.txt"; if (  -e ( $rdFl2 )  ){ system  ("rm -f $rdFl2"); }
	    	  
	    	}
	    	
	    	
	    }
	  } 
	
  }
	
	if ( $isThere_anything_new==1 ){
	  DirFileHandle::PrintDumper($JL_INterProScan_DoneWorkHsh,$JL_INterProScan_DoneWorkHsh_HASH) if (   (  defined ( $JL_INterProScan_DoneWorkHsh_HASH )  ) && (  ref ( $JL_INterProScan_DoneWorkHsh_HASH ) eq 'HASH'  )   );	
	  
	  
	  
	  
	}
	
	
	Interproscan::Del_JL_INterProScan_upGradeMark();
	
}


sub Eatract_interProscanHash_from_JL_ipsDatabase_inputFastaFile{  #  my $finalOutHash=Interproscan::Eatract_interProscanHash_from_JL_ipsDatabase_inputFastaFile ($inFastaFile);
	my ( $inFastaFile )=@_;
	
	my $warnMsgBody="\n\n\n   In package  Interproscan,\tIn sub Eatract_interProscanHash_from_JL_ipsDatabase_inputFastaFile,\n\n";
  my $warnMsgHead="\n\n\n$warnMsgBody";
  my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $InFastaSeHASH=ClustalwRun::NewReadInSeq_into_A_HASH($inFastaFile, 'fasta');
  
  my $theRealSeqHash;
  if (   (  defined ( $InFastaSeHASH )  ) && (  ref ( $InFastaSeHASH ) eq 'HASH'  ) && (  defined ( $InFastaSeHASH->{'1_id_seqHash'} )  ) && (  ref ( $InFastaSeHASH->{'1_id_seqHash'} ) eq 'HASH'  )   ) {
		foreach my $eachFastaName (    sort { $a cmp $b } (   keys (  %{ $InFastaSeHASH->{'1_id_seqHash'} }  )   )    ){
	  	#$outHash->{'1_id_seqHash'}->{ $seqObj->primary_id }->{'3_seque_Order'}=$seqNb; 
	  	$theRealSeqHash->{$eachFastaName}=$InFastaSeHASH->{'1_id_seqHash'}->{$eachFastaName}->{'4____sequence'} if (   (  defined ( $InFastaSeHASH->{'1_id_seqHash'}->{$eachFastaName}->{'4____sequence'}  )  ) && ( $InFastaSeHASH->{'1_id_seqHash'}->{$eachFastaName}->{'4____sequence'}=~m/\S+/ )   ) ; 
	  }
	}
  
  my $finalOutHash;
  if (   (  defined ( $theRealSeqHash )  ) && (  ref ( $theRealSeqHash ) eq 'HASH'  )   ){
  	$finalOutHash=Interproscan::Extract_interProscanHASH_from_JL_ipsDatabase ($theRealSeqHash);
  }
  
  return $finalOutHash;
	
}

sub Extract_interProscanHASH_from_JL_ipsDatabase{   #  Interproscan::Extract_interProscanHASH_from_JL_ipsDatabase ($seqHash);
	
	my ( $seqHash )=@_;
	
	my $warnMsgBody="\n\n\n   In package  Interproscan,\tIn sub Extract_interProscanHASH_from_JL_ipsDatabase,\n\n";
  my $warnMsgHead="\n\n\n$warnMsgBody";
  my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	
	
	
	my $JL_INterProScan_DoneWorkHsh_HASH;
  if (   (  defined ( $JL_INterProScan_DoneWorkHsh )  ) && ( $JL_INterProScan_DoneWorkHsh=~m/\S+/ ) && (  -e ( $JL_INterProScan_DoneWorkHsh )  )   ) {
	  $JL_INterProScan_DoneWorkHsh_HASH=Storable::retrieve  ( $JL_INterProScan_DoneWorkHsh );
	}
	
	my $newOut_ipsHASH;
	
	if (   (  defined ( $seqHash )  ) && (  ref ( $seqHash ) eq 'HASH'  )    ){
		foreach my $fastaName (    sort { $a cmp $b } (   keys (  %{ $seqHash } )   )    ){
		  my $seqString=$seqHash->{$fastaName};
		  my $relString=Interproscan::MakeGoodString_forInterProscan ($seqString);
		  
		  warn "20181218-1 ".$JL_INterProScan_DoneWorkHsh_HASH."\n";
		  warn "20181218-1 ".$JL_INterProScan_DoneWorkHsh_HASH->{'0_1_0_SeqKeyHASH'}."\n";
		  warn "20181218-2 ".$JL_INterProScan_DoneWorkHsh_HASH->{'0_1_0_SeqKeyHASH'}->{$relString}."\n";
		  warn "20181218-3 ".$JL_INterProScan_DoneWorkHsh_HASH->{'0_1_0_SeqKeyHASH'}->{$relString}->{'0_0_0_OtHsh_wholPath'}."\n";
		  
      if  (      (  defined ( $JL_INterProScan_DoneWorkHsh_HASH )  ) && (  ref ( $JL_INterProScan_DoneWorkHsh_HASH ) eq 'HASH'  ) 
              && (  defined ( $JL_INterProScan_DoneWorkHsh_HASH->{'0_1_0_SeqKeyHASH'} )  ) && (  ref ( $JL_INterProScan_DoneWorkHsh_HASH->{'0_1_0_SeqKeyHASH'} ) eq 'HASH'  )   
              && (  defined ( $JL_INterProScan_DoneWorkHsh_HASH->{'0_1_0_SeqKeyHASH'}->{$relString} )  ) && (  ref ( $JL_INterProScan_DoneWorkHsh_HASH->{'0_1_0_SeqKeyHASH'}->{$relString} ) eq 'HASH'  )  
              && (  defined ( $JL_INterProScan_DoneWorkHsh_HASH->{'0_1_0_SeqKeyHASH'}->{$relString}->{'0_0_0_OtHsh_wholPath'} )       ) 
              && (            $JL_INterProScan_DoneWorkHsh_HASH->{'0_1_0_SeqKeyHASH'}->{$relString}->{'0_0_0_OtHsh_wholPath'}=~m/\S+/ ) 
              && (  -e      ( $JL_INterProScan_DoneWorkHsh_HASH->{'0_1_0_SeqKeyHASH'}->{$relString}->{'0_0_0_OtHsh_wholPath'} )       )   
          )
      { 
      	my $single_ipsHASH=Storable::retrieve  ( $JL_INterProScan_DoneWorkHsh_HASH->{'0_1_0_SeqKeyHASH'}->{$relString}->{'0_0_0_OtHsh_wholPath'} );
      	
      	$newOut_ipsHASH->{$fastaName}= Storable::dclone ( $single_ipsHASH->{$JL_interProscanShowName } );
       
      }
      else {
      	my $dieMsg="$die_MsgHead in the hash file cannot found the follow string key:\n\$relString=$relString\n\n !!\n\n\n\n";
        print $dieMsg; die $dieMsg;
      }

		  
		  
		}
	}
	return $newOut_ipsHASH;
	
}

sub MakeGoodString_forInterProscan{  #  Interproscan::MakeGoodString_forInterProscan ($inString);
	my ($inString)=@_;
  
  my $warnMsgBody="\n\n\n   In package  Interproscan,\tIn sub MakeGoodString_forInterProscan,\n\n";
  my $warnMsgHead="\n\n\n$warnMsgBody";
    
  $inString=~s/\s+//g; $inString=~s/\*//g; $inString= uc ( $inString );
  return $inString;
}

sub BuildFastastring_forInterProscan{           #Interproscan::BuildFastastring_forInterProscan ($fastaName, $inString);
  my ($fastaName, $inString)=@_;
  my $warnMsgBody="\n\n\n   In package  Interproscan,\tIn sub BuildFastastring_forInterProscan,\n\n";
  my $warnMsgHead="\n\n\n$warnMsgBody";
  my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  if (   (  defined ( $fastaName )  ) && ( $fastaName=~/\S+/ )  ){} 
  else {
  	$fastaName='FASTA'; 
    warn "$die_MsgHead\nNo FASTA name set, use FASTA as name\n\n"; print  "$die_MsgHead\nNo FASTA name set, use FASTA as name\n\n";
  }
  #open (IN,">$FilePath")  or die "DIE!!!!!\n$warnMsgHead cannot create \$FilePath=$FilePath:$!\n\n";
  $fastaName=~s/^\s+//;  $fastaName=~s/\s+$//;  $inString=~s/\s+//g; $inString=~s/\*//g;
  #print IN ">$fastaName\n$inString\n\n";
  #close (IN);
  my $outPUT=">$fastaName\n$inString\n\n";
  return $outPUT;
}



sub GetSegLimit_size{   #   my $outLimit=Interproscan::GetSegLimit_size ($in_Seq_HASH, $core_number);    
	my ($in_Seq_HASH, $core_number)=@_;
	
	my $warnMsgBody="\n\n\n   In package  Interproscan,\tIn sub GetSegLimit_size,\n\n";
  my $warnMsgHead="\n\n\n$warnMsgBody";
  my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $subCallereIfm=DirFileHandle::print_SubCallerInform;
	
	my $biggest_chuck_seqLimit=10;
	
	my $total_length=0; my $total_number=0;
	if  (   (  defined ( $in_Seq_HASH )  ) && (  ref ( $in_Seq_HASH ) eq 'HASH'  )   ){
		foreach my $eachSeq(   keys (  %{ $in_Seq_HASH }  )   ){
			$total_length+=$in_Seq_HASH->{$eachSeq} if (   (  defined ( $in_Seq_HASH->{$eachSeq} )  ) && ( $in_Seq_HASH->{$eachSeq}=~m/\d+/ )   );
			$total_number++;
		}
	}
	if ( $total_number <=0 ){
		&Just_dieWork( "$die_MsgHead\n$subCallereIfm\n\nThe\$in_Seq_HASH=$in_Seq_HASH\nThe \$total_number=$total_number\n Something is wrong!!!\n" );
	}
	my $avrg_10_seq_length=$biggest_chuck_seqLimit*$total_length/$total_number;
	
	my $outLimit;
	if ( ($total_length > 0) && ($core_number > 0) ){
		$outLimit=$total_length/$core_number;
		if ( $outLimit >= $avrg_10_seq_length){
			$outLimit=$avrg_10_seq_length;
		}
	}
	else{
		&Just_dieWork( "$die_MsgHead\n$subCallereIfm\n\nThe\$in_Seq_HASH=$in_Seq_HASH\nThe \$core_number=$core_number\n Something is wrong!!!\n" );
	}
	return $outLimit;
	
}






sub dis_seq_hash_IntoGrpuped_HASH{   #   my $outPutSegKeyGroups=Interproscan::dis_seq_hash_IntoGrpuped_HASH ($org_inputSegSizeHash, $sizeLimit);                #将很多序列，按照序列长度大小的要求，放入一系列小于该大小要求的数组中，该数组可用于下一步的文件生成
  my ($org_inputSegSizeHash, $sizeLimit)=@_;
  
  my $warnMsgBody="\n\n\n   In package  Interproscan,\tIn sub dis_seq_hash_IntoGrpuped_HASH,\n\n";
  my $warnMsgHead="\n\n\n$warnMsgBody";
  my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $outPutSegKeyGroups;
  my $oPsKgIdx=0;
  
  my $inputSegSizeHash=Storable::dclone ( $org_inputSegSizeHash );
  foreach my $BigNumberKey (    sort { $inputSegSizeHash->{$b} <=> $inputSegSizeHash->{$a} } (   keys (  %{ $inputSegSizeHash }  )   )    ){                                #print "\n\n\cl\cl\n\nOutSide of if\n";
    if (   (  defined ( $inputSegSizeHash->{$BigNumberKey} )  ) && ( $inputSegSizeHash->{$BigNumberKey}=~m/\d+/ ) && ( $inputSegSizeHash->{$BigNumberKey} > 0 )  ){    
      
      my $bigNumberSize=$inputSegSizeHash->{$BigNumberKey};               #print "\$oPsKgIdx=$oPsKgIdx\t\$bigNumberSize=$bigNumberSize=\$inputSegSizeHash->{\$BigNumberKey}=\$inputSegSizeHash->{$BigNumberKey}\n";
      
      push @{ $outPutSegKeyGroups->[$oPsKgIdx] },$BigNumberKey;           #print "\n\n\nprint Dumper (\$outPutSegKeyGroups)\n****************Big-Out*********************\n\cl";    #print Dumper ($outPutSegKeyGroups);    #print "\n\n\n";        #print "\n\n\nprint Dumper (\$inputSegSizeHash)\n*****************Big__In***BeforeDelete*****\n\cl";    #print Dumper ($inputSegSizeHash);    #print "\n";
      delete $inputSegSizeHash->{$BigNumberKey};
      
      if ($bigNumberSize > $sizeLimit){           
      	#die "In sub disIDsIntoGrpups, \$bigNumberSize=$bigNumberSize > \$sizeLimit=$sizeLimit\n";
      }
      
      else{
        #print "\nprint Dumper (\$inputSegSizeHash)\n*****************Big__In****After-Delete****\n\cl";    #print Dumper ($inputSegSizeHash);    #print "\n\n\n";
        my $addedNumber=$bigNumberSize;                                #print "\$addedNumber=\$bigNumberSize=$addedNumber\n";#这个数字用来将 长度加和 以计算是否超过限度
        SFOREACHMARK: foreach my $SmallNBKey (    sort { $inputSegSizeHash->{$a} <=> $inputSegSizeHash->{$b} } (   keys (  %{ $inputSegSizeHash }  )   )    ){
          $addedNumber+=$inputSegSizeHash->{$SmallNBKey};  #print "\$addedNumber=$addedNumber \t\$addedNumber+=\$inputSegSizeHash->{\$SmallNBKey}=\$inputSegSizeHash->{$SmallNBKey}=$inputSegSizeHash->{$SmallNBKey}\n";
          if ($addedNumber > $sizeLimit){                  #print ">>>>>>>>>>\$addedNumber=$addedNumber > \$sizeLimit=$sizeLimit\n\cl\n";     
            last SFOREACHMARK;
          }
          else{                                            #print "<=<=<=<=<=<=\$addedNumber=$addedNumber <= \$sizeLimit=$sizeLimit\n";           
          	push @{ $outPutSegKeyGroups->[$oPsKgIdx] },$SmallNBKey;
          	#print "\n\n\nprint Dumper (\$outPutSegKeyGroups)\n*****************Sml-Out********************\n\cl";        #print Dumper ($outPutSegKeyGroups);        #print "\n\n\n";
          	delete $inputSegSizeHash->{$SmallNBKey};
          	#print "\n\n\nprint Dumper (\$inputSegSizeHash)\n******************Sml__In*******************\n\cl";        #print Dumper ($inputSegSizeHash);        #print "\n\n\n";
          }
          
        }
      
      }
      $oPsKgIdx++;
      
    }
    
  }
  
  return $outPutSegKeyGroups;
  
}









###########################################

sub BuildOneSeqFastaFile_forInterProscan{   #Interproscan::BuildOneSeqFastaFile_forInterProscan ($inString, $FilePath, $fastaName);
	my ($inString, $FilePath, $fastaName)=@_;
  my $warnMsgBody="\n\n\n   In package  Interproscan,\tIn sub BuildOneSeqFastaFile_forInterProscan,\n\n";
  my $warnMsgHead="\n\n\n$warnMsgBody";
  my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  if (   (  defined ( $fastaName )  ) && ( $fastaName=~/\S+/ )  ){} 
  else {
  	$fastaName='FASTA'; 
    warn "$die_MsgHead\nNo FASTA name set, use FASTA as name\n\n"; print  "$die_MsgHead\nNo FASTA name set, use FASTA as name\n\n";
  }
  
  open (IN,">$FilePath")  or die "$die_MsgHead cannot create \$FilePath=$FilePath:$!\n\n";
  
  $fastaName=~s/^\s+//;  $fastaName=~s/\s+$//;  $inString= Interproscan::MakeGoodString_forInterProscan ($inString);
  my $outPUT=">$fastaName\n$inString\n\n";
  print IN $outPUT;
  close (IN);
  
  return $outPUT;
}


sub BuildDomainOrderHash{
  my $domainTypeInformOrderHash;                              
  $domainTypeInformOrderHash->{'Protein family membership'}=0;
  $domainTypeInformOrderHash->{'Domains and repeats'      }=1;
  $domainTypeInformOrderHash->{'Molecular Function'       }=2;
  $domainTypeInformOrderHash->{'Biological Process'       }=3;
  $domainTypeInformOrderHash->{'Cellular Component'       }=4;
  $domainTypeInformOrderHash->{'Pathway'                  }=5;
  return $domainTypeInformOrderHash; 
}

sub BuildDomainShortNameHash{
  my $domainTypeShortNMhash; 
  $domainTypeShortNMhash->{'Protein family membership'}='Fml';
  $domainTypeShortNMhash->{'Domains and repeats'      }='Dom';
  $domainTypeShortNMhash->{'Molecular Function'       }='GoM';
  $domainTypeShortNMhash->{'Biological Process'       }='GoB';
  $domainTypeShortNMhash->{'Cellular Component'       }='GoC';
  $domainTypeShortNMhash->{'Pathway'                  }='Paw';
  return $domainTypeShortNMhash; 
}


sub InterProScanRun{
	

  my ($PepSeqFastaFile, $pepInterProScan_Result, $hperLk_pepIPS_rstHtmlDir)=@_;   warn  "20181218-2011-3 Interproscan::InterProScanRun($PepSeqFastaFile, $pepInterProScan_Result, $hperLk_pepIPS_rstHtmlDir)\n"; 	
                                                                                  print "20181218-2011-4 Interproscan::InterProScanRun($PepSeqFastaFile, $pepInterProScan_Result, $hperLk_pepIPS_rstHtmlDir)\n";
  #/home/fredjiang/
  my $ipScan_command="/home/fredjiang/tools/interProScan/interproscan-5-44.0/interproscan.sh  -i $PepSeqFastaFile -f xml -f tsv -f svg -f html -f gff3  -iprlookup -pa -goterms -b $pepInterProScan_Result";
  my $mkHTMLdir_command="mkdir -p $hperLk_pepIPS_rstHtmlDir";
  my $unzipHtml_in2_HTMLdir_Command="tar xzfv ${pepInterProScan_Result}.html.tar.gz -C $hperLk_pepIPS_rstHtmlDir";  print "20181218-2011-5 $ipScan_command\n";                 warn "\n$ipScan_command\n\n";
    
  system ("$ipScan_command");                                                                                       print "20181218-2011-6  $mkHTMLdir_command\n";              warn "\n$mkHTMLdir_command\n\n";
    
  system ("$mkHTMLdir_command");                                                                                    
                                                                                                                    print "20181218-2011-7 $unzipHtml_in2_HTMLdir_Command\n";  warn "\n$unzipHtml_in2_HTMLdir_Command\n\n";
  system ("$unzipHtml_in2_HTMLdir_Command");
  
}

sub SeqIn_PrasedHashOut{
  my ($inSeqFile, $outHashFile, $runOrNot, $outFilesDir, $delOrnot)=@_;
  my $outResDIR;
 
  if (  ( defined ($outFilesDir) )  &&  ($outFilesDir=~m/\S+/)  ){$outResDIR=$outFilesDir;}
  else {
  	if (  ( defined ($delOrnot) )  &&  ($delOrnot == 1)  ){   	$outFilesDir=TimeWork::GetTimeDirOrFileName; }
  	else                                                  {   	$outFilesDir="${inSeqFile}.RstDir"; }
  }
  
  my $mkdirCMD="mkdir -p $outFilesDir";  
  system ("$mkdirCMD");
  
  my $outRstFile="$outFilesDir/Rst";
  my $outRstHtml="$outFilesDir/Htm";
  my $outRstXml="$outFilesDir/Rst.xml";
  
  if (  ( defined ($runOrNot) )  &&  ($runOrNot ==1)  ){
  	if (   ( -e ($inSeqFile) ) && ( $outFilesDir eq "${inSeqFile}.RstDir" ) && ( -d ($outFilesDir) )   ){
  		system ("rm -rf $outFilesDir");
  		system ("$mkdirCMD");
  		
  	}
    &InterProScanRun( $inSeqFile, $outRstFile, $outRstHtml);
    if (  ( defined ($delOrnot) )  &&  ($delOrnot == 1)  ){   	
    	if ($outFilesDir =~ m/\d{8}\w\d{6}/){
    	  system ("rm -rf $outFilesDir") ;
    	  
      }
    }
  }
  
  my $outRstXmlHASH=&DumpXml($outRstXml);
  
  my $allHtmlDirFiles=DirFileHandle::getDirArray($outRstHtml);
  
  my $allHtmlHash;
  foreach my $eachHtmlFile (  @{ $allHtmlDirFiles }  ){ print "\$eachHtmlFile=$eachHtmlFile\n";
    if ($eachHtmlFile=~m/.html$/){ print "\$eachHtmlFile=$eachHtmlFile\n";
    	my $outHash=&praseInterProScanHtml( "$outRstHtml/$eachHtmlFile" );
    	foreach my $k0 (    sort { $a cmp $b } (   keys (  %{ $outHash } )   )    ){
    	  $allHtmlHash->{$k0}=$outHash->{$k0};
    	  if (  defined( $outRstXmlHASH->{$k0} )  ){
    	  	$allHtmlHash->{$k0}->{'Pathway'}=$outRstXmlHASH->{$k0};
    	  }
    	}
    }
  }
  
  if (  ( defined ($outHashFile) )  &&  ($outHashFile=~m/\S+/)  ){$outHashFile=$outHashFile;}
  else {$outHashFile=$inSeqFile."HASH";}
  
  DirFileHandle::PrintDumper($outHashFile,$allHtmlHash);
  return $allHtmlHash;
  
}


sub DumpXml{
  my ($xmlFile, $pathWayHashFile)=@_;
  my $xmlHashFile="${xmlFile}.Hash";
  
  
  if (  ( defined ($pathWayHashFile) )  &&  ($pathWayHashFile=~m/\w+/)  ){$pathWayHashFile=$pathWayHashFile;} 
  else {$pathWayHashFile="${xmlFile}.PathWay.Hash";}
  
  my $config = XML::Simple::XMLin($xmlFile);
  DirFileHandle::PrintDumper($xmlHashFile,$config);
  my $pathWayHash;  my $pathWayFound=0;
  foreach my $pepDna (   keys (  %{ $config }  )   ){    print "\$pepDna=$pepDna\n";
  	if (  ( ref ($config->{$pepDna}) ) eq 'ARRAY' ){
  		foreach my $eachPtArray (  @{ $config->{$pepDna} }  ){  	    #print "  \$eachPtArray->{'xref'}->{'id'}=$eachPtArray->{'xref'}->{'id'}\n";
  			
  			if (  ( ref ($eachPtArray->{'matches'}) ) eq 'HASH' ){
  			  foreach my $keyLev_0 (    sort { $a cmp $b } (   keys (  %{ $eachPtArray->{'matches'} } )   )    ){   print "    \$keyLev_0=$keyLev_0\n";
  			  	my $valLev_0=$eachPtArray->{'matches'}->{$keyLev_0};                                                print "    \$valLev_0=$valLev_0\n";
  			  	my $refLev_0=ref ($valLev_0);                                                                       print "    \$refLev_0=$refLev_0\n\n";  
  			  	
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
  			                  	my $valLev_4=$valLev_3->{$keyLev_4};                                                print "            \$valLev_4=$valLev_4\n";
  			                  	my $refLev_4=ref ($valLev_4);                                                       print "            \$refLev_4=$refLev_4\n\n";
  			                  	
  			                  	
  			                  	if ($keyLev_3 eq 'pathway-xref'){
  			              	      my $acId=$valLev_2->{'ac'};   
  			              	      print "              \$eachPtArray->{'xref'}->{'id'}=$eachPtArray->{'xref'}->{'id'}\t";
  			              	      print "\$acId=$acId\t";
  			              	      #print "\$valLev_4->{'id'}=$valLev_4->{'id'}\t";
  			              	      #print "\$valLev_4->{'db'}=$valLev_4->{'db'}\t";
  			              	      #print "\$keyLev_4=$keyLev_4\n\n\n";
  			              	      $pathWayFound=1;
  			              	      
  			              	      if ($keyLev_4 eq 'id'){
  			              	        $pathWayHash->{ $eachPtArray->{'xref'}->{'id'} }->{ $valLev_4 }->{'Definition'    }=$valLev_3->{'db'};
  			              	        $pathWayHash->{ $eachPtArray->{'xref'}->{'id'} }->{ $valLev_4 }->{'DomainDataBase'}=$valLev_3->{'name'};
  			              	      }
  			              	      $pathWayFound=2;
  			              	    }
  			     
  			          	        if ( $refLev_4 eq 'HASH' ){
  			                      foreach my $keyLev_5 (    sort { $a cmp $b } (   keys (  %{ $valLev_4 } )   )    ){   print "              \$keyLev_5=$keyLev_5\n";
  			                      	my $valLev_5=$valLev_4->{$keyLev_5};                                                print "              \$valLev_5=$valLev_5\n";
  			                      	my $refLev_5=ref ($valLev_5);                                                       print "              \$refLev_5=$refLev_5\n\n";
  			      
  			                        
  			                        if ( ($keyLev_3 eq 'pathway-xref') && ($keyLev_5 eq 'id') ){
  			                          $pathWayFound=3;
  			                          
  			                          $pathWayHash->{ $eachPtArray->{'xref'}->{'id'} }->{ $valLev_5 }->{'Definition'    }=$keyLev_4;
  			              	          $pathWayHash->{ $eachPtArray->{'xref'}->{'id'} }->{ $valLev_5 }->{'DomainDataBase'}=$valLev_4->{'db'};
  			                          
  			                        }
  			                        
  			                              
  			      
  			      
  			                        if ( $refLev_5 eq 'HASH' ){
  			                          foreach my $keyLev_6 (    sort { $a cmp $b } (   keys (  %{ $valLev_5 } )   )    ){   print "              \$keyLev_6=$keyLev_6\n";
  			                      	    my $valLev_6=$valLev_5->{$keyLev_6};                                                print "              \$valLev_6=$valLev_6\n";
  			                      	    my $refLev_6=ref ($valLev_6);                                                       print "              \$refLev_6=$refLev_6\n\n\n\n";
  			                          
  			                            
  			                          
  			                          
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
  			}
  	  }
  	}
  }
  
  if (  ($pathWayFound==2) || ($pathWayFound==3)  ){
    DirFileHandle::PrintDumper($pathWayHashFile,$pathWayHash);
  }
  elsif ($pathWayFound==1){
    die "In sub \$pathWayFound==1,\n \$xmlFile=$xmlFile, \$pathWayHashFile=$pathWayHashFile\nWrong Wrong Wrong!\n\n\n";
  }
  
  
  
  
  return $pathWayHash;
}









##sub8.3          praseInterProScanHtml{  #解析 interproscan输出的html文件
sub praseInterProScanHtml{  #解析 interproscan输出的html文件
  my ( $file )=@_;
  
  #my $Protein_name_regular_Expression="\\d+\\w\\d+\\w\\.\\S";
  
  
  my $finalBigHash;
  
  
  #$finalBigHash->{$proTeinNm}->{'file'}=$file;  
  #$finalBigHash->{$proTeinNm}->{'length'}=$proLength; ##################^^^^^^^^^^____!____^^^^^^^^^^###################
  #$finalBigHash->{$proTeinNm}->{$ClassType}->{$signID}=[$signNm,$signUrl,$subArray];##################^^^^^^^^^^____!____^^^^^^^^^^###################
  #                                                                       $subArray->[$j]->{$subSigID}=[$subSigUrl,$subSigNm,$subSubArray];
  #                                                                                                                          $subSubArray->[$realK]=[$segDataBase, $segID, $segNm, $seqUrl, $segCav];
  #$finalBigHash->{$proTeinNm}->{$noPir}->{$noPirID}=[$noPirNm,$noPirUrl,$subNoPirArray,$noPirDataBase];##################^^^^^^^^^^____!____^^^^^^^^^^###################
  #                                                                      $subNoPirArray->[$realJ]=[$subNoPirDataBase,$subNoPirID,$subNoPirNm,$subNoPirUrl,$subNoPirCav];

  
  my $familyClassHash;
  my $domainClassHash;  #mkdir "noFml";            print "sdsfsdfsdfds\n";

  
  my $domainRelation;
  my $noInformation;

  my $clsFyHash;

  my $definiionHash;
  
  
  if ($file=~m/\.html$/){ 	  print "\$file =$file \n";     print "\n";
    my $tree=HTML::TreeBuilder->new;      $tree->parse_file($file);
    DirFileHandle::PrintDumper("$file.Tree",$tree);  
      
    my $treeHd=                    $tree->{'_content'}[1]{'_content'}[0]{'_content'}[0]{'_content'}[1]{'_content'}[1];    
    
    my $treeHd_sbtp_1=             $tree->{'_content'}[1]{'_content'}[0]{'_content'}[0]{'_content'}[0]{'_content'}[1];   #DirFileHandle::PrintDumper("treeHd_sbtp_1.HS",$treeHd_sbtp_1); 
    
    my $proTeinNm;                                                            #print ref ($treeHd_sbtp_1), "\n\n" ;     #print ref ($treeHd_sbtp_1->{'_content'}), "\n\n" ;     #print ref ($treeHd_sbtp_1->{'_content'}[0]), "\n\n" ;     #print ref ($treeHd_sbtp_1->{'_content'}[0]{'_content'}), "\n\n" ; 

    my $interProInformExist=0;

    if (  ( ref ($treeHd_sbtp_1->{'_content'}[0]) ) eq 'HTML::Element'  ){
      my $proTeinNm_sbtp_1= $treeHd_sbtp_1->{'_content'}[0]{'_content'}[0]  ;    
      $proTeinNm_sbtp_1=~s/\s//g;       print "\$proTeinNm_sbtp_1=$proTeinNm_sbtp_1\n";
     	$treeHd=  $treeHd_sbtp_1; 
     	$proTeinNm=  $proTeinNm_sbtp_1;           	
    }
    else{ 
    	$proTeinNm=               $treeHd->{'_content'}[0]{'_content'}[0]; $proTeinNm=~s/\s//g;  print "\$proTeinNm=$proTeinNm\n";
    	$interProInformExist=1;
    }
    
   
    
    $finalBigHash->{$proTeinNm}->{'file'}=$file;  
    #-----------------------------------------------------+++--------------------------------------------------------------------------
    ###########输入蛋白名                                         14m23641A.anophagefferens;Grx1;ACJI01001294_+7             
    print "\n";
    printf "%-20s\t%-20s\n",  "1--", $treeHd->{'_content'}[0]{'_content'}[0];                                ## 固定的
    
    #-----------------------------------------------------+++--------------------------------------------------------------------------
    ###########输入蛋白有多少个amino acids                        217 amino acids                                                  
    print "\n";
    
    my $WholeProLength=$treeHd->{'_content'}[1]{'_content'}[0]{'_content'}[1]{'_content'}[0]; my $proLength;
    if ($WholeProLength=~m/(\d+)\s+amino\s+acid/){$proLength=$1;} else {die "\$file=$file\n\$WholeProLength=$WholeProLength\n";}
    #printf "%-20s\t%-20s\n",  "2--", $treeHd->{'_content'}[1]{'_content'}[0]{'_content'}[1]{'_content'}[0];  ## 固定的
    printf "%-20s\t%-20s\n",  "2--", $WholeProLength;
    $finalBigHash->{$proTeinNm}->{'length'}=$proLength; ##################^^^^^^^^^^____!____^^^^^^^^^^###################
    #-----------------------------------------------------+++--------------------------------------------------------------------------
    
    if ($interProInformExist==1){
    
      my $idexNb=2;
      print "2222 \$idexNb=$idexNb\n";
      #-----------------------------------------------------+++--------------------------------------------------------------------------
      ###########家族 头文字                                        Protein family membership                                  
      my $prtFmlMbshp=$treeHd->{'_content'}[$idexNb]{'_content'}[0];    print "\n";    print "\n";    printf "%-20s\t%-20s\n","2.5--", $prtFmlMbshp;    print "\n";
      if ($prtFmlMbshp=~m/Protein family membership/){} else{die "\$file=$file\n\$prtFmlMbshp=$prtFmlMbshp\n";}
      #-----------------------------------------------------+++--------------------------------------------------------------------------
      $idexNb++;  #!!!!!!!!!!!!!!!!      如果不die 这里 就变成了  3  
      print "3333 \$idexNb=$idexNb\n";
      #-----------------------------------------------------+++--------------------------------------------------------------------------
      
      my $fmlClsed=0;    # 这个也没有用到
      
      #-----------------------------------------------------+++--------------------------------------------------------------------------
      for (my $i=0; $i<@{ $treeHd->{'_content'}[$idexNb]{'_content'} }; $i++){
      	if (  ( ref ($treeHd->{'_content'}[$idexNb]{'_content'}[$i]) )  eq 'HTML::Element'  ){
          ###########家族 id           
          my $fmlID=$treeHd->{'_content'}[$idexNb]{'_content'}[$i]{'_content'}[0]{'_content'}[1]{'_content'}[0] ;          printf "%-20s\t%-20s\n",  "$i)3--", $fmlID;      #$finalBigHash->{$proTeinNm}->{$prtFmlMbshp}->{$fmlID}->{} ##########################################################
          
          ###########家族 名字       
          my $fmlDefin=$treeHd->{'_content'}[$idexNb]{'_content'}[$i]{'_content'}[0]{'_content'}[0]{'_content'}[0];       printf "%-20s\t%-20s\n",  "$i)4--", $fmlDefin;      #$finalBigHash->{$proTeinNm}->{$prtFmlMbshp}->{$fmlDefin}->{} ##########################################################
          ###########家族 网址
          my $fmlUrl=$treeHd->{'_content'}[$idexNb]{'_content'}[$i]{'_content'}[0]{'_content'}[0]{'href'};                printf "%-20s\t%-20s\n",  "$i)5--", $fmlUrl;       print "\n";
          
          $fmlID=~s/^\s*\(\s*//; $fmlID=~s/\s*\)\s*$//;
          $finalBigHash->{$proTeinNm}->{$prtFmlMbshp}->{$fmlID}->{'Definition'}=$fmlDefin;
          $finalBigHash->{$proTeinNm}->{$prtFmlMbshp}->{$fmlID}->{'Url'       }=$fmlUrl;
          
          
          
          ############################################################
          #############  下面的 没有用在我们的程序中  ################
          if ($fmlID=~m/\(\s*(\S+)\s*\)/){ $fmlID=$1;
            if (defined ($familyClassHash->{$fmlID}) ){              print "cccccccccaaaaaaaaaaaaaaa\$familyClassHash->{$fmlID},$file\n\n";
              push @{$familyClassHash->{$fmlID}}, $file;
            }
            else {
              $familyClassHash->{$fmlID}=[$file];                    print "ccccccccccbbbbbbbbbbbbb\$familyClassHash->{$fmlID}=[$file]\n\n";
            }
            #print "&print_all_sub_array( \$familyClassHash->{$fmlID} )\n;"; &print_all_sub_array( $familyClassHash->{$fmlID} );
          }
          else {
            $fmlClsed=1;
          }      
          #############  上面的 没有用在我们的程序中  ################
          ############################################################     
        }
      }                                                            print "\n\n";
      #if ($fmlClsed == 1) {    	system ("cp -f $file noFml/$file");    }
      
      
      #-----------------------------------------------------+++--------------------------------------------------------------------------
      $idexNb++;  #!!!!!!!!!!!!!!!!         这里  变成了  4  
      print "4444 \$idexNb=$idexNb\n";
      #-----------------------------------------------------+++--------------------------------------------------------------------------
      
      #-----------------------------------------------------+++--------------------------------------------------------------------------
      ###########保守区域 头文字                                    Domains and repeats                                                                           Domains and repeats
      my $DomAndRpHd=$treeHd->{'_content'}[$idexNb]{'_content'}[0];    print "\n";    printf "%-20s\t%-20s\n","5.5--", $DomAndRpHd;    print "\n";
      if ($DomAndRpHd=~m/Domains and repeats/){} else{die "\$file=$file\n\$DomAndRpHd=$DomAndRpHd\n";}
      #-----------------------------------------------------+++--------------------------------------------------------------------------
      $idexNb++;  #!!!!!!!!!!!!!!!!      如果不die 这里 就变成了  5  
      print "5555 \$idexNb=$idexNb\n";
      #-----------------------------------------------------+++--------------------------------------------------------------------------
      
      #-----------------------------------------------------+++--------------------------------------------------------------------------
      for (my $i=0; $i<@{ $treeHd->{'_content'}[$idexNb]{'_content'} }; $i++){
        if (  ( ref ($treeHd->{'_content'}[$idexNb]{'_content'}[0]{'_content'}[1]{'_content'}[$i]) eq 'HTML::Element' )  ){  
          ###########保守区域 标记                                                                                                 
          my $DomainWords =$treeHd->{'_content'}[$idexNb]{'_content'}[0]{'_content'}[1]{'_content'}[$i]{'_content'}[0]{'_content'}[0];               #print ref $treeHd->{'_content'}[$idexNb]{'_content'}[0]{'_content'}[1]{'_content'}[$i], "\n";                             
          my $treeDomainHd=$treeHd->{'_content'}[$idexNb]{'_content'}[0]{'_content'}[1]{'_content'}[$i]{'_content'}[1]{'_content'}[0];        
          #printf "%-20s\t%-20s\n",  "$i)6--", $treeDomainHd;      print "\n"; 
          printf "%-20s\t%-20s\n",  "$i)6--", $DomainWords;      print "\n"; 
                        
          for (my $j=0; $j<@{ $treeDomainHd->{'_content'} }; $j++){  
          	#DirFileHandle::PrintDumper("treeDomainHd.$j.HS",$treeDomainHd->{'_content'}[$j]);        	#print "RefHere \$j=$j ", ref ($treeDomainHd->{'_content'}[$j]), "\n";
          	if (  ( ref ($treeDomainHd->{'_content'}[$j]) eq 'HTML::Element' )  ){    #print "\$j=$j\t$treeDomainHd->{'_content'}[$j]{'_content'}[1]{'_content'}[0]{'_content'}[0]{'_content'}[0]{'_content'}[0]\n\n";
              if (defined ($treeDomainHd->{'_content'}[$j]{'_content'}[1]{'_content'}[0]{'_content'}[0]{'_content'}[0]{'_content'}[0])){
                if (       $treeDomainHd->{'_content'}[$j]{'_content'}[1]{'_content'}[0]{'_content'}[0]{'_content'}[0]{'_content'}[0]=~m/\w+/){
                
                  ###########保守区域 覆盖范围 比如 10 - 19                                                                                
                  printf "%-20s\t%-20s\n",  "$i)($j)7--",  $treeDomainHd->{'_content'}[$j]{'_content'}[0]{'_content'}[1];                print "\n";
                  my $domainRegion=$treeDomainHd->{'_content'}[$j]{'_content'}[0]{'_content'}[1];
                  printf "%-20s\t%-20s\n",  "$i)($j)7--",  $treeDomainHd->{'_content'}[$j]{'_content'}[1]{'_content'}[0]{'_content'}[0]{'_content'}[0]{'_content'}[0];                print "\n";  
                  my $domainDefine=$treeDomainHd->{'_content'}[$j]{'_content'}[1]{'_content'}[0]{'_content'}[0]{'_content'}[0]{'_content'}[0];
                  printf "%-20s\t%-20s\n",  "$i)($j)7--",  $treeDomainHd->{'_content'}[$j]{'_content'}[1]{'_content'}[0]{'_content'}[0]{'_content'}[0]{'href'};                print "\n";  
                  my $domainUrl=$treeDomainHd->{'_content'}[$j]{'_content'}[1]{'_content'}[0]{'_content'}[0]{'_content'}[0]{'href'};
                  printf "%-20s\t%-20s\n",  "$i)($j)7--",  $treeDomainHd->{'_content'}[$j]{'_content'}[1]{'_content'}[0]{'_content'}[0]{'_content'}[1]{'_content'}[0];                print "\n";       
                  my $domainID=$treeDomainHd->{'_content'}[$j]{'_content'}[1]{'_content'}[0]{'_content'}[0]{'_content'}[1]{'_content'}[0];
                  $domainID=~s/^\s*\(\s*//; $domainID=~s/\s*\)\s*$//; 
                  my $start; my $end;
                  if ($domainRegion=~m/^\s*((\d+(?:,\d+)?)\s*-\s*(\d+(?:,\d+)?))\s*$/){$start=$2; $end=$3; $domainRegion=$1; $start=~s/,//g; $start=~s/\s+//g; $end=~s/,//g; $end=~s/\s+//g;} else {die"\$domainRegion=$domainRegion !!!!!!!!!!!!!!\n\n\n\n\n";}
                  
                  $finalBigHash->{$proTeinNm}->{$DomAndRpHd}->{'DomainLineNBhash'}->{$i}->{'DomainIDhash'}->{$domainID}->{'Definition'  }=$domainDefine;
                  $finalBigHash->{$proTeinNm}->{$DomAndRpHd}->{'DomainLineNBhash'}->{$i}->{'DomainIDhash'}->{$domainID}->{'Url'         }=$domainUrl;
                  $finalBigHash->{$proTeinNm}->{$DomAndRpHd}->{'DomainLineNBhash'}->{$i}->{'DomainIDhash'}->{$domainID}->{'CavRegHash'  }->{$domainRegion}->{'RegionStart'}=$start;
                  $finalBigHash->{$proTeinNm}->{$DomAndRpHd}->{'DomainLineNBhash'}->{$i}->{'DomainIDhash'}->{$domainID}->{'CavRegHash'  }->{$domainRegion}->{'RegionEnd'  }=$end;
                  
                  $finalBigHash->{$proTeinNm}->{$DomAndRpHd}->{'DomainLineNBhash'}->{$i}->{'CavReg2DomainIDHash'}->{$domainRegion}->{'RegionStart'}=$start;
                  $finalBigHash->{$proTeinNm}->{$DomAndRpHd}->{'DomainLineNBhash'}->{$i}->{'CavReg2DomainIDHash'}->{$domainRegion}->{'RegionEnd'  }=$end;
                  $finalBigHash->{$proTeinNm}->{$DomAndRpHd}->{'DomainLineNBhash'}->{$i}->{'CavReg2DomainIDHash'}->{$domainRegion}->{'DomainID'   }=$domainID;
                  
                               #$treeDomainHd->{'_content'}[$j]{'_content'}[1]{'_content'}[0]{'_content'}[0]{'_content'}[1]{'_content'}[0];
                  if (  ( ref ($treeDomainHd->{'_content'}[$j]{'_content'}[1]{'_content'}[0]{'_content'}[0]{'_content'}[2]{'_content'}[0]) eq 'HTML::Element' )  ){
                  	#warn "\$proTeinNm=$proTeinNm\n";
                  	printf "%-20s\t%-20s\n",  "$i)($j)7.2--1",  $treeDomainHd->{'_content'}[$j]{'_content'}[1]{'_content'}[0]{'_content'}[0]{'_content'}[2]{'_content'}[0];  print "\n";
                  	printf "%-20s\t%-20s\n",  "$i)($j)7.2--2",  $treeDomainHd->{'_content'}[$j]{'_content'}[1]{'_content'}[0]{'_content'}[0]{'_content'}[2]{'_content'}[0]{'_content'}[0];  print "\n";
                  	printf "%-20s\t%-20s\n",  "$i)($j)7.2--3",  $treeDomainHd->{'_content'}[$j]{'_content'}[1]{'_content'}[0]{'_content'}[0]{'_content'}[2]{'_content'}[0]{'_content'}[0]{'_content'}[0];  print "\n";
                    #printf "%-20s\t%-20s\n",  "$i)($j)7.2--4",  $treeDomainHd->{'_content'}[$j]{'_content'}[1]{'_content'}[0]{'_content'}[0]{'_content'}[2]{'_content'}[0]{'_content'}[0]{'_content'}[0]{'_content'}[0];  print "\n"; 
                    #my $subDef=$treeDomainHd->{'_content'}[$j]{'_content'}[1]{'_content'}[0]{'_content'}[0]{'_content'}[2]{'_content'}[0]{'_content'}[0]{'_content'}[0]{'_content'}[0];
                    my $subDef=$treeDomainHd->{'_content'}[$j]{'_content'}[1]{'_content'}[0]{'_content'}[0]{'_content'}[2]{'_content'}[0]{'_content'}[0]{'_content'}[0];                                                    #warn "\$subDef=$subDef\n";
                    #printf "%-20s\t%-20s\n",  "$i)($j)7.2--",  $treeDomainHd->{'_content'}[$j]{'_content'}[1]{'_content'}[0]{'_content'}[0]{'_content'}[2]{'_content'}[0]{'_content'}[0]{'_content'}[0]{'href'};  print "\n"; 
                    printf "%-20s\t%-20s\n",  "$i)($j)7.2--",  $treeDomainHd->{'_content'}[$j]{'_content'}[1]{'_content'}[0]{'_content'}[0]{'_content'}[2]{'_content'}[0]{'_content'}[0]{'href'};  print "\n"; 
                    #my $subUrl=$treeDomainHd->{'_content'}[$j]{'_content'}[1]{'_content'}[0]{'_content'}[0]{'_content'}[2]{'_content'}[0]{'_content'}[0]{'_content'}[0]{'href'};                                                           #warn "\$subUrl=$subUrl\n";
                    my $subUrl=$treeDomainHd->{'_content'}[$j]{'_content'}[1]{'_content'}[0]{'_content'}[0]{'_content'}[2]{'_content'}[0]{'_content'}[0]{'href'};      
                    
                    #printf "%-20s\t%-20s\n",  "$i)($j)7.2--",  $treeDomainHd->{'_content'}[$j]{'_content'}[1]{'_content'}[0]{'_content'}[0]{'_content'}[2]{'_content'}[0]{'_content'}[0]{'_content'}[1]{'_content'}[0];  print "\n"; 
                    #my $subID=$treeDomainHd->{'_content'}[$j]{'_content'}[1]{'_content'}[0]{'_content'}[0]{'_content'}[2]{'_content'}[0]{'_content'}[0]{'_content'}[1]{'_content'}[0];                                                     #warn "\$subID=$subID\n";
                    printf "%-20s\t%-20s\n",  "$i)($j)7.2--",  $treeDomainHd->{'_content'}[$j]{'_content'}[1]{'_content'}[0]{'_content'}[0]{'_content'}[2]{'_content'}[0]{'_content'}[1]{'_content'}[0];  print "\n"; 
                    my $subID=$treeDomainHd->{'_content'}[$j]{'_content'}[1]{'_content'}[0]{'_content'}[0]{'_content'}[2]{'_content'}[0]{'_content'}[1]{'_content'}[0];
                    
                    $subID=~s/^\s*\(\s*//; $subID=~s/\s*\)\s*$//; 
                    $finalBigHash->{$proTeinNm}->{$DomAndRpHd}->{'DomainLineNBhash'}->{$i}->{'DomainIDhash'}->{$domainID}->{'subDomainsHash'}->{$subID}->{'Definition'} =$subDef;
                    $finalBigHash->{$proTeinNm}->{$DomAndRpHd}->{'DomainLineNBhash'}->{$i}->{'DomainIDhash'}->{$domainID}->{'subDomainsHash'}->{$subID}->{'Url'       } =$subUrl;
                    
                  
                  
                  }
                  
                  
                  
                  
                  #my $DomainRoot=$treeDomainHd->{'_content'}[$j]{'_content'}[1]{'_content'}[0]{'_content'}[0];
                  #&printSubDomainInf ($DomainRoot, 0, $file, $proTeinNm); #$inRoot, $rootLv, $fileNm, $proNm, $upperNm
                 }
              }
            }
	        } 
	      }
	      print "\n\n";
	    }
	    print "\n\n\n";
	    
	    #-----------------------------------------------------+++--------------------------------------------------------------------------
      $idexNb++;  #!!!!!!!!!!!!!!!!         这里  变成了  6  
      print "6666 \$idexNb=$idexNb\n";
      #-----------------------------------------------------+++--------------------------------------------------------------------------
      
      #-----------------------------------------------------+++--------------------------------------------------------------------------
	    #                                                              Detailed signature matches                                
	    my $DetailSignatureMatches=$treeHd->{'_content'}[$idexNb]{'_content'}[0];	  print "\n";	  printf "%-20s\t%-20s\n","11.0--",$DetailSignatureMatches;    print "\n";
      #-----------------------------------------------------+++--------------------------------------------------------------------------
	    if ($DetailSignatureMatches=~m/Detailed signature matches/){ 
        $idexNb++;  #!!!!!!!!!!!!!!!!      如果符合它 这里 就变成了  7  
        print "7777 \$idexNb=$idexNb\n";
      }
      #-----------------------------------------------------+++--------------------------------------------------------------------------
      
      #-----------------------------------------------------+++--------------------------------------------------------------------------
      #my $treeHd=                    $tree->{'_content'}[1]{'_content'}[0]{'_content'}[0]{'_content'}[1]{'_content'}[1];    
      #my $treeHd_sbtp_1=             $tree->{'_content'}[1]{'_content'}[0]{'_content'}[0]{'_content'}[0]{'_content'}[1];
	    #$sigRoot=$tree->{'_content'}[1]{'_content'}[0]{'_content'}[0]{'_content'}[0]{'_content'}[1]{'_content'}[7]{'_content'}[0];
	    my $isThereIPRifom=0;
	    my $sigRoot=$treeHd->{'_content'}[$idexNb]{'_content'}[0];            #	  print "\cl\cl"; &print_all_sub_array ($sigRoot);  print "\cl\cl";
	    if (  ( ref ($sigRoot->{'_content'}) eq 'ARRAY' )  ){
	      for (my $i=0; $i<@{ $sigRoot->{'_content'} }; $i++ ){
	        #                                                                                                                        
	        #标明是 family 还是 domain 还是别的                   
	        my $ClassType=$sigRoot->{'_content'}[$i]{'_content'}[0]{'_content'}[0]{'_content'}[0]{'title'};
	        $ClassType=~s/^\s+//; $ClassType=~s/\s+$//;
	        printf "%-20s\t%-20s\n", "$i)title:", $sigRoot->{'_content'}[$i]{'_content'}[0]{'_content'}[0]{'_content'}[0]{'title'};
          print "\n";
          
	        #这个是family或这domain的id，显示在左边，还有小图标 ，以及该id对应的网址
	        my $signID=$sigRoot->{'_content'}[$i]{'_content'}[0]{'_content'}[0]{'_content'}[2]{'_content'}[0];   $signID=~s/^\s+//;  $signID=~s/\s+$//;
	        my $signUrl=$sigRoot->{'_content'}[$i]{'_content'}[0]{'_content'}[0]{'_content'}[2]{'href'};        $signUrl=~s/^\s+//; $signUrl=~s/\s+$//;                     
	        printf "%-20s\t%-20s\n", "$i)13-2-", $sigRoot->{'_content'}[$i]{'_content'}[0]{'_content'}[0]{'_content'}[2]{'_content'}[0];
	        printf "%-20s\t%-20s\n", "$i)14-2-", $sigRoot->{'_content'}[$i]{'_content'}[0]{'_content'}[0]{'_content'}[2]{'href'};
	        print "\n";
	        
	        #这个是family或这domain的名字，显示在上面 ，以及该名字对应的网址，和上面id的网址相同  
	        my $signNm=$sigRoot->{'_content'}[$i]{'_content'}[0]{'_content'}[1]{'_content'}[0]{'_content'}[0];   $signNm=~s/^\s+//;  $signNm=~s/\s+$//;                         
	        printf "%-20s\t%-20s\n", "$i)15-2-", $sigRoot->{'_content'}[$i]{'_content'}[0]{'_content'}[1]{'_content'}[0]{'_content'}[0];
	        printf "%-20s\t%-20s\n", "$i)16-2-", $sigRoot->{'_content'}[$i]{'_content'}[0]{'_content'}[1]{'_content'}[0]{'href'};
	        print "\n";
	        
	        #$finalBigHash->{$proTeinNm}->{$ClassType}->{$signID}=[$signNm,$signUrl]
	        
	        my $subArray;
	        for (my $j=0; $j<@{ $sigRoot->{'_content'}[$i]{'_content'}[1]{'_content'}[1]{'_content'} }; $j++){ 
	          my $sigInerRoot_1=$sigRoot->{'_content'}[$i]{'_content'}[1]{'_content'}[1];
	          #
	          my $subSigID=$sigInerRoot_1->{'_content'}[$j]{'_content'}[0]{'_content'}[0]{'_content'}[0]; $subSigID=~s/^\s+//;  $subSigID=~s/\s+$//;
	          my $subSigUrl=$sigInerRoot_1->{'_content'}[$j]{'_content'}[0]{'_content'}[0]{'href'};       $subSigUrl=~s/^\s+//;  $subSigUrl=~s/\s+$//;  
	          my $subSigNm=$sigInerRoot_1->{'_content'}[$j]{'_content'}[0]{'_content'}[0]{'title'};       $subSigNm=~s/^\s+//;   $subSigNm=~s/\s+$//;                              
	          printf "%-20s\t%-20s\n", "$i)$j)17-2-",   $sigInerRoot_1->{'_content'}[$j]{'_content'}[0]{'_content'}[0]{'_content'}[0];
	          printf "%-20s\t%-20s\n", "$i)$j)18-2-",   $sigInerRoot_1->{'_content'}[$j]{'_content'}[0]{'_content'}[0]{'href'};
	          printf "%-20s\t%-20s\n", "$i)$j)18.5-2-", $sigInerRoot_1->{'_content'}[$j]{'_content'}[0]{'_content'}[0]{'title'};
	          printf "%-20s\t%-20s\n", "$i)$j)19-2-",   $sigInerRoot_1->{'_content'}[$j]{'_content'}[0]{'_content'}[2]{'_content'}[0]  if defined ($sigInerRoot_1->{'_content'}[$j]{'_content'}[0]{'_content'}[2]{'_content'}[0]);
	          #$subArray->[$j]->{$subSigID}=[$subSigUrl,$subSigNm,]
	          print "\n";
	          $isThereIPRifom=1;
	          
	          #下面是具体的各个以图形和浮动气泡表示的各domain的信息
	          my $subSubArray; my $realK=0;
	          for (my $k=0; $k< @{                               $sigInerRoot_1->{'_content'}[$j]{'_content'}[1]{'_content'}[0]{'_content'} }; $k++){
	          	my $sigInerRoot_2=$sigInerRoot_1->{'_content'}[$j]{'_content'}[1]{'_content'}[0];
	          	if (  ( ref ($sigInerRoot_2->{'_content'}[$k]) )  eq 'HTML::Element'  ){	      	  
                if (defined (                                 $sigInerRoot_2->{'_content'}[$k]{'_content'}[1]{'_content'}[4]{'_content'}[0] )){
                  if (                                        $sigInerRoot_2->{'_content'}[$k]{'_content'}[1]{'_content'}[4]{'_content'}[0] =~m/\w+/){
                    
                    my $segCav=$sigInerRoot_2->{'_content'}[$k]{'_content'}[0]{'_content'}[1];                     $segCav=~s/^\s+//;     $segCav=~s/\s+$//;
                    my $segDataBase=$sigInerRoot_2->{'_content'}[$k]{'_content'}[1]{'_content'}[0]{'_content'}[0]; $segDataBase=~s/^\s+//;$segDataBase=~s/\s+$//;
                    my $segID=$sigInerRoot_2->{'_content'}[$k]{'_content'}[1]{'_content'}[4]{'_content'}[0];       $segID=~s/^\s+//;      $segID=~s/\s+$//;
                    my $segNm=$sigInerRoot_2->{'_content'}[$k]{'_content'}[1]{'_content'}[6]{'_content'}[0];       $segNm=~s/^\s+//;      $segNm=~s/\s+$//;
	                  my $seqUrl=$sigInerRoot_2->{'_content'}[$k]{'_content'}[1]{'_content'}[4]{'href'};             $seqUrl=~s/^\s+//;     $seqUrl=~s/\s+$//;
	                  printf "%-20s\t%-20s\n","$i)$j)$k)20-2-",  $sigInerRoot_2->{'_content'}[$k]{'_content'}[0]{'_content'}[1];
	                  printf "%-20s\t%-20s\n","$i)$j)$k)21-2-",  $sigInerRoot_2->{'_content'}[$k]{'_content'}[1]{'_content'}[0]{'_content'}[0];
	                  printf "%-20s\t%-20s\n","$i)$j)$k)22-2-",  $sigInerRoot_2->{'_content'}[$k]{'_content'}[1]{'_content'}[4]{'_content'}[0];
	                  printf "%-20s\t%-20s\n","$i)$j)$k)23-2-",  $sigInerRoot_2->{'_content'}[$k]{'_content'}[1]{'_content'}[6]{'_content'}[0];
                    printf "%-20s\t%-20s\n","$i)$j)$k)24-2-",  $sigInerRoot_2->{'_content'}[$k]{'_content'}[1]{'_content'}[4]{'href'};
	                  #$subSubArray->[$realK]=[$segDataBase, $segID, $segNm, $seqUrl, $segCav];
	                  my $start; my $end;
                    if ($segCav=~m/^\s*((\d+(?:,\d+)?)\s*-\s*(\d+(?:,\d+)?))\s*$/){$start=$2; $end=$3; $segCav=$1;  $start=~s/,//g; $start=~s/\s+//g; $end=~s/,//g; $end=~s/\s+//g;} else {die"\$segCav=$segCav !!!!!!!!!!!!!!\n\n\n\n\n";}
	                  $segNm=~s/^\s*\(\s*//; $segNm=~s/\s*\)\s*$//; 
	                  
	                  $subSubArray->[$realK]->{'DomainDataBase'}=$segDataBase;
	                  $subSubArray->[$realK]->{'DomainID'      }=$segID;
	                  $subSubArray->[$realK]->{'Definition'    }=$segNm;
	                  $subSubArray->[$realK]->{'Url'           }=$seqUrl;
	                  
	                  
	                  $subSubArray->[$realK]->{'CavRegHash'    }->{$segCav}->{'RegionStart'   }=$start;
	                  $subSubArray->[$realK]->{'CavRegHash'    }->{$segCav}->{'RegionEnd'     }=$end;
	                  
	                 
	                  
	                  $realK++;
	                  
	                  print "\n";
	                }
	              } 
	            }
	          }
	          #$subArray->[$j]->{$subSigID}=[$subSigUrl,$subSigNm,$subSubArray];
	          
	          $subArray->[$j]->{$subSigID}->{'subSubArray'   }=$subSubArray;
	          $subArray->[$j]->{$subSigID}->{'Definition'    }=$subSigNm;
	          $subArray->[$j]->{$subSigID}->{'Url'           }=$subSigUrl;
	         
	          
	          print "\n\n";
	        }
	        #$finalBigHash->{$proTeinNm}->{$ClassType}->{$signID}=[$signNm,$signUrl,$subArray];##################^^^^^^^^^^____!____^^^^^^^^^^###################
	        $finalBigHash->{$proTeinNm}->{$ClassType}->{$signID}->{'Definition'    }=$signNm;
	        $finalBigHash->{$proTeinNm}->{$ClassType}->{$signID}->{'Url'           }=$signUrl;
	        $finalBigHash->{$proTeinNm}->{$ClassType}->{$signID}->{'subArray'      }=$subArray;
	        
	        
	        $definiionHash->{$signID}=$signNm;
	        if ( defined ($clsFyHash->{$ClassType}->{$signID} ) ) { push @{ $clsFyHash->{$ClassType}->{$signID} }, $proTeinNm; }
	        else                                                  { $clsFyHash->{$ClassType}->{$signID}=[$proTeinNm];          }
	        
	        print "\n\n\n";
	      }
	    }
	    print "\n\n\n\n";
	    
      #-----------------------------------------------------+++--------------------------------------------------------------------------	  
	    if ($isThereIPRifom==1){
	      $idexNb++;  #!!!!!!!!!!!!!!!!      如果$isThereIPRifom==1 这里 就变成了  8  
	      print "8888 \$idexNb=$idexNb\n";
	    }
      #-----------------------------------------------------+++--------------------------------------------------------------------------	  
      
      
      
      
      
      my $noIPRiformYesOrNo=0;   
	                               
	                               
	    my $noIprRoot=$treeHd->{'_content'}[$idexNb];
	    #warn "1 ref ($noIprRoot)\n";   
	    
	                                 
      if (    (   (  ref ($noIprRoot )  ) eq 'HTML::Element'   ) && (   (  ref ($noIprRoot->{'_content'}[0] )  ) eq 'HTML::Element'   ) && (   (  ref ($noIprRoot->{'_content'}[0]{'_content'}[0] )  ) eq 'HTML::Element'   ) && (   (  ref ($noIprRoot->{'_content'}[0]{'_content'}[0]{'_content'}[0] )  ) eq 'HTML::Element'   ) && ( $noIprRoot->{'_content'}[0]{'_content'}[0]{'_content'}[0]{'title'} eq 'Unintegrated signatures' )    ){
      	
      	#warn "2 ref ($noIprRoot->{'_content'}[0])\n"; 
      	#warn "3 ref ($noIprRoot->{'_content'}[0]{'_content'}[0])\n";  
	      #warn "4 ref ($noIprRoot->{'_content'}[0]{'_content'}[0]{'_content'}[0])\n"; 
      	#-----------------------------------------------------+++--------------------------------------------------------------------------
	       
	      
	      print $noIprRoot->{'_content'}[0],"\n";                                      
	      print $noIprRoot->{'_content'}[0]{'_content'}[0],"\n";                       
	      print $noIprRoot->{'_content'}[0]{'_content'}[0]{'_content'}[0],"\n";        
	      #DirFileHandle::PrintDumper("ttt",$noIprRoot->{'_content'}[0]{'_content'}[0]{'_content'}) ;
	       
	      #标明是 family 还是 domain 还是别的                   
	      printf "%-20s\t%-20s\n", "  )title:", $noIprRoot->{'_content'}[0]{'_content'}[0]{'_content'}[0]{'title'};
        
        print "\n";
        #这个是family或这domain的id，显示在左边，还有小图标 ，
        my $noPir=$noIprRoot->{'_content'}[0]{'_content'}[0]{'_content'}[1]; $noPir=~s/^\s+//;     $noPir=~s/\s+$//;
        printf "%-20s\t%-20s\n", "  )25--",   $noIprRoot->{'_content'}[0]{'_content'}[0]{'_content'}[1];
        #这个是family或这domain的名字，显示在上面 ，                               
        printf "%-20s\t%-20s\n", "  )26--",   $noIprRoot->{'_content'}[0]{'_content'}[1]{'_content'}[0];
        print "\n";
        
        print ref ($noIprRoot->{'_content'}[1]{'_content'}[1]{'_content'}), "\n";
        if (  ( ref ($noIprRoot->{'_content'}[1]{'_content'}[1]{'_content'}) )  eq 'ARRAY'  ){  print "aaaaaaaaaaaa\n\n";
	        for (my $i=0; $i<@{ $noIprRoot->{'_content'}[1]{'_content'}[1]{'_content'} }; $i++ ){    print "\$i=$i\n";
	          my $noIprInerRoot_1=$noIprRoot->{'_content'}[1]{'_content'}[1] ;
          
            printf "%-20s\t%-20s\n", "$i))27--",   $noIprInerRoot_1->{'_content'}[$i]{'_content'}[0]{'_content'}[0]{'_content'}[0];
            printf "%-20s\t%-20s\n", "$i))28--",   $noIprInerRoot_1->{'_content'}[$i]{'_content'}[0]{'_content'}[0]{'href'};
            printf "%-20s\t%-20s\n", "$i))29--",   $noIprInerRoot_1->{'_content'}[$i]{'_content'}[0]{'_content'}[0]{'title'};
            printf "%-20s\t%-20s\n", "$i))30--",   $noIprInerRoot_1->{'_content'}[$i]{'_content'}[0]{'_content'}[2]{'_content'}[0]  if ( defined ($noIprInerRoot_1->{'_content'}[$i]{'_content'}[0]{'_content'}[2]{'_content'}[0]) );
            print "\n";
          
            my $subNoPirArray; my $realJ=0;
            #下面是具体的各个以图形和浮动气泡表示的各domain的信息
            for (my $j=0; $j< @{                        $noIprInerRoot_1->{'_content'}[$i]{'_content'}[1]{'_content'}[0]{'_content'} }; $j++){
            	my $noIprInerRoot_2=$noIprInerRoot_1->{'_content'}[$i]{'_content'}[1]{'_content'}[0];
            	if (  ( ref ($noIprInerRoot_2->{'_content'}[$j]) )  eq 'HTML::Element'  ){
                if (defined (                                 $noIprInerRoot_2->{'_content'}[$j]{'_content'}[1]{'_content'}[4]{'_content'}[0] )){
                  if (                                        $noIprInerRoot_2->{'_content'}[$j]{'_content'}[1]{'_content'}[4]{'_content'}[0] =~m/\w+/){
                    
                    
                    my $subNoPirCav=     $noIprInerRoot_2->{'_content'}[$j]{'_content'}[0]{'_content'}[1];                $subNoPirCav=~     s/^\s+//;$subNoPirCav=~     s/\s+$//;
                    my $subNoPirDataBase=$noIprInerRoot_2->{'_content'}[$j]{'_content'}[1]{'_content'}[0]{'_content'}[0]; $subNoPirDataBase=~s/^\s+//;$subNoPirDataBase=~s/\s+$//;
                    my $subNoPirID=      $noIprInerRoot_2->{'_content'}[$j]{'_content'}[1]{'_content'}[4]{'_content'}[0]; $subNoPirID=~      s/^\s+//;$subNoPirID=~      s/\s+$//;
                    my $subNoPirNm=      $noIprInerRoot_2->{'_content'}[$j]{'_content'}[1]{'_content'}[6]{'_content'}[0]; $subNoPirNm=~      s/^\s+//;$subNoPirNm=~      s/\s+$//;
                    my $subNoPirUrl=     $noIprInerRoot_2->{'_content'}[$j]{'_content'}[1]{'_content'}[4]{'href'};        $subNoPirUrl=~     s/^\s+//;$subNoPirUrl=~     s/\s+$//;
                    printf "%-20s\t%-20s\n","$i)$j))31--",  $noIprInerRoot_2->{'_content'}[$j]{'_content'}[0]{'_content'}[1];
                    printf "%-20s\t%-20s\n","$i)$j))32--",  $noIprInerRoot_2->{'_content'}[$j]{'_content'}[1]{'_content'}[0]{'_content'}[0];
                    printf "%-20s\t%-20s\n","$i)$j))33--",  $noIprInerRoot_2->{'_content'}[$j]{'_content'}[1]{'_content'}[4]{'_content'}[0];
                    printf "%-20s\t%-20s\n","$i)$j))34--",  $noIprInerRoot_2->{'_content'}[$j]{'_content'}[1]{'_content'}[6]{'_content'}[0];
                    printf "%-20s\t%-20s\n","$i)$j))35--",  $noIprInerRoot_2->{'_content'}[$j]{'_content'}[1]{'_content'}[4]{'href'};
                    
                    #$subNoPirArray->[$realJ]=[$subNoPirDataBase,$subNoPirID,$subNoPirNm,$subNoPirUrl,$subNoPirCav];
                    
                    
                    
                    
                    my $start; my $end;
                    if ($subNoPirCav=~m/^\s*((\d+(?:,\d+)?)\s*-\s*(\d+(?:,\d+)?))\s*$/){$start=$2; $end=$3; $subNoPirCav=$1;$start=~s/,//g; $start=~s/\s+//g; $end=~s/,//g; $end=~s/\s+//g;} else {die"\$subNoPirCav=$subNoPirCav !!!!!!!!!!!!!!\n\n\n\n\n";}
	                  $subNoPirNm=~s/^\s*\(\s*//; $subNoPirNm=~s/\s*\)\s*$//; 
	                  
	                  $subNoPirArray->[$realJ]->{'DomainDataBase'}=$subNoPirDataBase;
	                  $subNoPirArray->[$realJ]->{'DomainID'      }=$subNoPirID;
	                  $subNoPirArray->[$realJ]->{'Definition'    }=$subNoPirNm;
	                  $subNoPirArray->[$realJ]->{'Url'           }=$subNoPirUrl;
	                 
	                  $subNoPirArray->[$realJ]->{'CavRegHash'    }->{$subNoPirCav}->{'RegionStart'   }=$start;
	                  $subNoPirArray->[$realJ]->{'CavRegHash'    }->{$subNoPirCav}->{'RegionEnd'     }=$end;
                    
                    
                    
                    
                    
                    
                    
                    $realJ++;
                    $noIPRiformYesOrNo=1;
                    print "\n";
                  }
                }    
              }
	          }
	          my $noPirDataBase=$subNoPirArray->[0]->{'DomainDataBase'};
	          my $noPirID=      $subNoPirArray->[0]->{'DomainID'      };
	          my $noPirNm=      $subNoPirArray->[0]->{'Definition'    };  $noPirNm=~s/^\s*\(\s*//; $noPirNm=~s/\s*\)\s*$//; 
	          my $noPirUrl=     $subNoPirArray->[0]->{'Url'           };
	          print "\n\n";
	          #$finalBigHash->{$proTeinNm}->{$noPir}->{$noPirID}=[$noPirNm,$noPirUrl,$subNoPirArray,$noPirDataBase];##################^^^^^^^^^^____!____^^^^^^^^^^###################
	          
	          $finalBigHash->{$proTeinNm}->{$noPir}->{$noPirID}->{'DomainDataBase'}=$noPirDataBase;
	          $finalBigHash->{$proTeinNm}->{$noPir}->{$noPirID}->{'Definition'    }=$noPirNm;
	          $finalBigHash->{$proTeinNm}->{$noPir}->{$noPirID}->{'Url'           }=$noPirUrl;
	          $finalBigHash->{$proTeinNm}->{$noPir}->{$noPirID}->{'subNoPirArray' }=$subNoPirArray;
	          
	          
	          $definiionHash->{$noPirID}=$noPirNm;
	          if ( defined ($clsFyHash->{$noPir}->{$noPirID} ) ) { push @{ $clsFyHash->{$noPir}->{$noPirID} }, $proTeinNm; }
	          else                                             { $clsFyHash->{$noPir}->{$noPirID}=[$proTeinNm];          }	    
	          
	        }
	      }
	      	  
	      print "\n\n\n";
        #-----------------------------------------------------+++--------------------------------------------------------------------------
      
      	
      	
      }
      
      
      
      
      
      #-----------------------------------------------------+++--------------------------------------------------------------------------	  
	    if ($noIPRiformYesOrNo==1){
	      $idexNb++;  #!!!!!!!!!!!!!!!!      如果$noIPRiformYesOrNo==1 这里 就变成了  9  
	      print "9999 \$idexNb=$idexNb\n";
	    }
      #-----------------------------------------------------+++--------------------------------------------------------------------------	  
      
      
      
      my $GoTitle=$treeHd->{'_content'}[$idexNb]{'_content'}[0]; $GoTitle=~s/^\s+//; $GoTitle=~s/\s+$//;    printf "%-20s\t%-20s\n","GO title:--",  $GoTitle;    print "\n\n"; 
      
      #-----------------------------------------------------+++--------------------------------------------------------------------------	  
      $idexNb++;  #!!!!!!!!!!!!!!!!       这里 就变成了  10  
      print "1010 \$idexNb=$idexNb\n";
      #-----------------------------------------------------+++--------------------------------------------------------------------------	  
      
      
      my $GoClassType;
      for (my $i=0; $i<@{ $treeHd->{'_content'}[$idexNb]{'_content'}[0]{'_content'} }; $i++){ print "\$i=$i\n";
        my $GoRoot=$treeHd->{'_content'}[$idexNb]{'_content'}[0];
        my $GoClass=$GoRoot->{'_content'}[$i]{'_content'}[0];     printf "%-20s\t%-20s\n","GoClass",  $GoClass, "\n";      
        if ( (defined ($GoClass)) && ($GoClass=~m/^\w+(\s+\w+){1,2}$/) ){
          printf "%-20s\t%-20s\n","GO_sub_class",  $GoClass;
          $GoClassType=$GoClass;  $GoClassType=~s/^\s+//; $GoClassType=~s/\s+$//;
          print "\n";
          $i++;
        
          #printf "%-20s\t%-20s\n","GO  1111:--", $GoRoot->{'_content'}[$i]{'_content'}[0];   
          if (  ( ref ($GoRoot->{'_content'}[$i]{'_content'}[0]{'_content'}) )  eq 'ARRAY'  ){
            for (my $j=0; $j<@{ $GoRoot->{'_content'}[$i]{'_content'}[0]{'_content'} }; $j++){   print "\$j=$j\n";
            	my $GoIdRoot=$GoRoot->{'_content'}[$i]{'_content'}[0];
            	if ( (defined ($GoIdRoot->{'_content'}[$j]{'_content'}[0])) && ($GoIdRoot->{'_content'}[$j]{'_content'}[0]=~m/\w+/) ){
            	  my $GoID=$GoIdRoot->{'_content'}[$j]{'_content'}[0];
            	  my $GoUrl=$GoIdRoot->{'_content'}[$j]{'href'};
            	  $j++;  # 这里需要加一,因为后面的一个数据，也就是$GoIDdef,是下一个 $j 作为序号 指引的数据
            	  my $GoIDdef=$GoIdRoot->{'_content'}[($j)];
            	  
            	  
            	  printf "%-20s\t%-20s\n","GO__ID:",       $GoID;
            	  printf "%-20s\t%-20s\n","GO__describ:",  $GoIDdef;
            	  printf "%-20s\t%-20s\n","GO__Url:",      $GoUrl;
            	  print "\n";
            	  if ($GoID=~m/GO:\d+/){
            	    #$finalBigHash->{$proTeinNm}->{$GoTitle}->{$GoClassType}->{$GoID}=[$GoIDdef,$GoUrl];##################^^^^^^^^^^____!____^^^^^^^^^^###################
            	    
            	    $finalBigHash->{$proTeinNm}->{$GoTitle}->{$GoClassType}->{$GoID}->{'Definition'    }=$GoIDdef;
            	    $finalBigHash->{$proTeinNm}->{$GoTitle}->{$GoClassType}->{$GoID}->{'Url'           }=$GoUrl;
            	    
            	    $definiionHash->{$GoID}=$GoIDdef;
            	  }
            	  $j++;  #这里需要加一, 因为，下面的结构中 又有一个需要加1的结构
            	}
            }
            print "\n\n";  
          }
        }
      }
      print "\n\n\n";
      
      #$tree=$tree->delete;
    }
	  
	  
	}
	  
	print "\n\n\n\n";

  
  DirFileHandle::PrintDumper("${file}.HS",$finalBigHash) ;
  return $finalBigHash;
}




sub MakeDomainExcelFromIterProScan{
  my ( $in_fastaFile,
       $in_clustalwOrNot,
       $in_interProScanRunOrNot,                  #$in_interProscanRstHash,
       $out_Dir,
       $out_excelFile,                            #$in_clustalwFile,
       $level_1_KeyHash,
       $level_1_KeyOrderHash,
       $level_2_KeyHash,       
       $level_2_KeyOrderHash,
       $level_3_KeyHash,
       $level_3_KeyOrderHash
       
       
     )=@_;
     
   my $seqIDhash=ClustalwRun::readInSeq($in_fastaFile, 'fasta');
 	
  
  
  #1  首先，生成一个Hash，以指定的 关键字作为多层 hash的 不同层的 Key。 这里是用 species作为第一层， 蛋白的功能分类作为第二层，每一个蛋白ID，必须只能对应 一种 功能分类 和 家族。
  my $multiLevelHash;
  my $all_level1_Hash;  #这个是Level1 的所有Key的hash
  my $all_level2_Hash;  #这个是Level2 的所有Key的hash
  my $all_level3_Hash;  #这个是Level3 的所有Key的hash
  
  if (  ( ref ( $seqIDhash->{'id_seqHash'} ) ) eq 'HASH' ){
    foreach my $seqID (    sort { $seqIDhash->{'id_seqHash'}->{$a}->{'seqOrder'} <=> $seqIDhash->{'id_seqHash'}->{$b}->{'seqOrder'} } (   keys (  %{ $seqIDhash->{'id_seqHash'} }  )   )    ){  #'seqIDhash' 这个key，指向的是该SeqID在原来那个fas文件中的顺序
      my $level_1_Matrix_key='UnSpecified_1'; if(  defined ( $level_1_KeyHash->{$seqID} )  ) { $level_1_Matrix_key=$level_1_KeyHash->{$seqID}; }
      my $level_2_Matrix_key='UnSpecified_2'; if(  defined ( $level_2_KeyHash->{$seqID} )  ) { $level_2_Matrix_key=$level_2_KeyHash->{$seqID}; }
      my $level_3_Matrix_key='UnSpecified_3'; if(  defined ( $level_3_KeyHash->{$seqID} )  ) { $level_3_Matrix_key=$level_3_KeyHash->{$seqID}; }
      $multiLevelHash->{$level_1_Matrix_key}->{$level_2_Matrix_key}->{$level_3_Matrix_key}->{$seqID}=1; print "\$multiLevelHash->{$level_1_Matrix_key}->{$level_2_Matrix_key}->{$level_3_Matrix_key}->{$seqID}=$multiLevelHash->{$level_1_Matrix_key}->{$level_2_Matrix_key}->{$level_3_Matrix_key}->{$seqID}\n";
      $all_level1_Hash->{$level_1_Matrix_key}=1;
      $all_level2_Hash->{$level_2_Matrix_key}=1;
      $all_level3_Hash->{$level_3_Matrix_key}=1;
    }
  }
  
  if (    (  defined ( $level_1_KeyOrderHash )  )  || (  ( ref ( $level_1_KeyOrderHash ) ) eq 'HASH' )    ) {} else { $level_1_KeyOrderHash=ForeachHash::buildOrderHash( $all_level1_Hash,'string'); }
  if (    (  defined ( $level_2_KeyOrderHash )  )  || (  ( ref ( $level_2_KeyOrderHash ) ) eq 'HASH' )    ) {} else { $level_2_KeyOrderHash=ForeachHash::buildOrderHash( $all_level2_Hash,'string'); }
  if (    (  defined ( $level_3_KeyOrderHash )  )  || (  ( ref ( $level_3_KeyOrderHash ) ) eq 'HASH' )    ) {} else { $level_3_KeyOrderHash=ForeachHash::buildOrderHash( $all_level3_Hash,'string'); }
  






  
  
  my $level1_ItScClu_Hash;
  
  my $mkDirCmd="mkdir -p $out_Dir"; if ( -d ($out_Dir) ){} else { warn "$mkDirCmd\n"; print "$mkDirCmd\n";  system ($mkDirCmd); } 
  my $BsNmOut_Dir=File::Basename::basename $out_Dir; warn "\$BsNmOut_Dir=$BsNmOut_Dir\n";
  	
  
  my $lv1Idx=0;
  if (  ( ref ( $multiLevelHash ) ) eq 'HASH' ){
    foreach my $level_1_key1 (    sort { $level_1_KeyOrderHash->{$a} <=> $level_1_KeyOrderHash->{$b} } (   keys (  %{ $multiLevelHash }  )   )    ){ 
    	my $lv1_ecDir="$out_Dir/lv1_${lv1Idx}";
    	my $mkDirCmd_lv1="mkdir -p $lv1_ecDir"; if ( -d ($mkDirCmd_lv1) ){} else { warn "$mkDirCmd_lv1\n"; print "$mkDirCmd_lv1\n";  system ($mkDirCmd_lv1); } 
      
      #对每一个 $level_1_key1, 比如说，这是一个家族的类别，那么这个key下面的 level2 lv3的key都属于这个家族。那么，对这一个家族的蛋白，就需要 进行比对和 interproscan的操作
      
      my $lv1_seqFile="$out_Dir/lv1_${lv1Idx}/lv1_seqFile.txt";     $level1_ItScClu_Hash->{$level_1_key1}->{'lv1_seqFile'}=$lv1_seqFile;    my $lv1_seqFile_lk="$BsNmOut_Dir/lv1_${lv1Idx}/lv1_seqFile.txt";     $level1_ItScClu_Hash->{$level_1_key1}->{'lv1_seqFile_lk'}=$lv1_seqFile_lk;        
      my $lv1_cluFile="$out_Dir/lv1_${lv1Idx}/lv1_cluFile.msf";     $level1_ItScClu_Hash->{$level_1_key1}->{'lv1_cluFile'}=$lv1_cluFile;    my $lv1_cluFile_lk="$BsNmOut_Dir/lv1_${lv1Idx}/lv1_cluFile.msf";     $level1_ItScClu_Hash->{$level_1_key1}->{'lv1_cluFile_lk'}=$lv1_cluFile_lk;  
      my $lv1_ItScHsh="$out_Dir/lv1_${lv1Idx}/lv1_ItScHsh";         $level1_ItScClu_Hash->{$level_1_key1}->{'lv1_ItScHsh'}=$lv1_ItScHsh;      
      my $lv1_ISHRead="$out_Dir/lv1_${lv1Idx}/lv1_ItScHsh.read.txt";$level1_ItScClu_Hash->{$level_1_key1}->{'lv1_ISHRead'}=$lv1_ISHRead;    my $lv1_ISHRead_lk="$BsNmOut_Dir/lv1_${lv1Idx}/lv1_ItScHsh.read.txt";$level1_ItScClu_Hash->{$level_1_key1}->{'lv1_ISHRead_lk'}=$lv1_ISHRead_lk;  
      my $lv1_ItScDir="$out_Dir/lv1_${lv1Idx}/lv1_ItScDir";         $level1_ItScClu_Hash->{$level_1_key1}->{'lv1_ItScDir'}=$lv1_ItScDir;    my $lv1_ItScDir_lk="$BsNmOut_Dir/lv1_${lv1Idx}/lv1_ItScDir";         $level1_ItScClu_Hash->{$level_1_key1}->{'lv1_ItScDir_lk'}=$lv1_ItScDir_lk;  
      my $Lv1_ItClPng="$out_Dir/lv1_${lv1Idx}/Lv1_ItClPng.png";     $level1_ItScClu_Hash->{$level_1_key1}->{'Lv1_ItClPng'}=$Lv1_ItClPng;    my $Lv1_ItClPng_lk="$BsNmOut_Dir/lv1_${lv1Idx}/Lv1_ItClPng.png";     $level1_ItScClu_Hash->{$level_1_key1}->{'Lv1_ItClPng_lk'}=$Lv1_ItClPng_lk;  
      
      my $lv1_seqIOobj_OUT =Bio::SeqIO->new(-file   => ">$lv1_seqFile",    
                                            -format => 'fasta'         );
      
      my $lv2Idx=0;
      if (  ( ref ( $multiLevelHash->{$level_1_key1} ) ) eq 'HASH' ){
        foreach my $level_2_key1 (    sort { $level_2_KeyOrderHash->{$a} <=> $level_2_KeyOrderHash->{$b} } (   keys (  %{ $multiLevelHash->{$level_1_key1} }  )   )    ){ 
          
          
          
          #以后可以启用# my $lv2_ecDir="$out_Dir/lv1_${lv1Idx}/lv2_${lv2Idx}";          
    	    #以后可以启用# my $mkDirCmd_lv2="mkdir -p $lv2_ecDir"; if ( -d ($out_Dir) ){} else { warn "$mkDirCmd_lv2\n"; print "$mkDirCmd_lv2\n";  system ($mkDirCmd_lv2); } 
          
          my $lv3Idx=0;
          if (  ( ref ( $multiLevelHash->{$level_1_key1}->{$level_2_key1} ) ) eq 'HASH' ){
            foreach my $level_3_key1 (    sort { $level_3_KeyOrderHash->{$a} <=> $level_3_KeyOrderHash->{$b} } (   keys (  %{ $multiLevelHash->{$level_1_key1}->{$level_2_key1} }  )   )    ){ 
     
              if (  ( ref ( $multiLevelHash->{$level_1_key1}->{$level_2_key1}->{$level_3_key1} ) ) eq 'HASH' ){
                foreach my $SequnceID (    sort { $seqIDhash->{'id_seqHash'}->{$a}->{'seqOrder'} <=> $seqIDhash->{'id_seqHash'}->{$b}->{'seqOrder'} } (   keys (  %{ $multiLevelHash->{$level_1_key1}->{$level_2_key1}->{$level_3_key1} }  )   )    ){ 
                	
                	my $SeqObj = Bio::Seq->new( -display_id => $SequnceID,
                                                     -seq => $seqIDhash->{'id_seqHash'}->{$SequnceID}->{'sequence'} );
                	$lv1_seqIOobj_OUT->write_seq($SeqObj); 
                	
                }
              }
              $lv3Idx++;
              
            }
          }
          $lv2Idx++;
        }
      }
      
      if ($in_clustalwOrNot == 1){ClustalwRun::runingClustalw($lv1_seqFile, $lv1_cluFile);}      
     
      my $lv1_ItScanHash=Interproscan::SeqIn_PrasedHashOut($lv1_seqFile, $lv1_ItScHsh,$in_interProScanRunOrNot,$lv1_ItScDir );      
      $level1_ItScClu_Hash->{$level_1_key1}->{'lv1_ItScanHash'}=$lv1_ItScanHash;       
      
      ClustalwRun::PrintPngForClustalw($lv1_cluFile, $lv1_ItScanHash, $Lv1_ItClPng);                                                
      
      $level1_ItScClu_Hash->{$level_1_key1}->{'lv1_CluSeqOrdr'}=ClustalwRun::readInSeqAndUseClustalwOrder($lv1_seqFile, 'fasta', $lv1_cluFile, 'msf');
      
      
      $lv1Idx++;
      
    }
  }
  
  


   
  
   
   
  
  #3 建立3个hash， 主要是 两个类型的导航Matrix
  my $mainKey_MatrixHash;     #建立一个hash，用来装各种二维表
  my @matrix_type=('MainKey_Lv1', 'MainKey_Lv2');
  my $matrix_type_OrderHash;  
  $matrix_type_OrderHash->{'MainKey_Lv1'}=$level_1_KeyOrderHash;
  $matrix_type_OrderHash->{'MainKey_Lv2'}=$level_2_KeyOrderHash;
  
 
  my $fmtKeyWordsHash;             #这个hash,主要用来生成 excel的格式关键字 hash。其中装下不同的 格式关键字，然后在后面的 步骤中，给每个关键字 对应生成一种 格式。
  
  my $allIDHash_for_EachLevle1Type;#这个hash,用来生成 每个蛋白功能分类中的 相应(Family,Domain,GO,Pathway)每一个类型的ID的 全集的
  my $ordHash=&BuildDomainOrderHash();                         #这个hash,用来生成 各个 类别的 DOmain kegg family Go 的排序
  my $snmHash=&BuildDomainShortNameHash();                     #这个hash,用来生成 各个 类别的 DOmain kegg family Go 的缩写 
  
  
  
  
  if (  ( ref ( $multiLevelHash ) ) eq 'HASH' ){
    foreach my $level_1_key (    sort { $level_1_KeyOrderHash->{$a} <=> $level_1_KeyOrderHash->{$b} } (   keys (  %{ $multiLevelHash }  )   )    ){ 
    
      if (  ( ref ( $multiLevelHash->{$level_1_key} ) ) eq 'HASH' ){
        foreach my $level_2_key (    sort { $level_2_KeyOrderHash->{$a} <=> $level_2_KeyOrderHash->{$b} } (   keys (  %{ $multiLevelHash->{$level_1_key} }  )   )    ){ 
          
          if (  ( ref ( $multiLevelHash->{$level_1_key}->{$level_2_key} ) ) eq 'HASH' ){
            foreach my $level_3_key (    sort { $level_3_KeyOrderHash->{$a} <=> $level_3_KeyOrderHash->{$b} } (   keys (  %{ $multiLevelHash->{$level_1_key}->{$level_2_key} }  )   )    ){ 
     
              if (  ( ref ( $multiLevelHash->{$level_1_key}->{$level_2_key}->{$level_3_key} ) ) eq 'HASH' ){  #In_level1_ItScClu_Hash->{$level_1_key}->{'lv1_CluSeqOrdr'}->{'id_seqHash'}->{$a}->{'seqOrder'}
                foreach my $SequnceID (    sort { $level1_ItScClu_Hash->{$level_1_key}->{'lv1_CluSeqOrdr'}->{'id_seqHash'}->{$a}->{'seqOrder'} <=>$level1_ItScClu_Hash->{$level_1_key}->{'lv1_CluSeqOrdr'}->{'id_seqHash'}->{$b}->{'seqOrder'} } (   keys (  %{ $multiLevelHash->{$level_1_key}->{$level_2_key}->{$level_3_key} }  )   )    ){ 
                
                  print "\$SequnceID=$SequnceID\n";  #warn "\$SequnceID=$SequnceID\n\n\n"; #sleep (5);
                  
                  my $domTempHash=&GetFmlDomGoKeggHash($level1_ItScClu_Hash->{$level_1_key}->{'lv1_ItScanHash'}, $SequnceID);
                  my $domHash=$domTempHash->{'realHash'};
                  
                  
                  
                  
                  
                  if (  ( ref ( $domHash ) ) eq 'HASH' ){
                    foreach my $domType (    sort { $ordHash->{$a} <=> $ordHash->{$b} } (   keys (  %{ $domHash }  )   )    ){ 
                      
                      if (  ( ref ( $domHash->{$domType} ) ) eq 'HASH' ){
                        foreach my $domID (    sort { $a cmp $b } (   keys (  %{ $domHash->{$domType} }  )   )    ){ 
                        	
                        	
                        	#将有dom信息的这些对应的 节点，记入'DoHavePointMatrix'中；
                        	$mainKey_MatrixHash->{'DoHavePointMatrix'}->{$level_1_key}->{$level_2_key}->{$level_3_key}->{$SequnceID}->{$domType}->{$domID}=1;
                        	
                        	
                        	
                        	#记录下 一个蛋白功能类别中，相应的类型的 Domain信息中的 所有ID的全集， 这里 $level_3_key  $level3AllKey, 都是为了照顾 有第三维信息的情况下，加入的中间结构 变量
                        	$allIDHash_for_EachLevle1Type->{$level_1_key}->{$level_3_key }->{$domType}->{$domID}=1;   # warn "\$allIDHash_for_EachLevle1Type->{$level_1_key}->{$level_3_key }->{$domType}->{$domID}=$allIDHash_for_EachLevle1Type->{$level_1_key}->{$level_3_key }->{$domType}->{$domID}\n";
                        	
                        	
                        	# 装填 以  物种 和 蛋白功能 为 行标识符 的 matirx             		      	
                        	
                          
                        	my @keyArray= ($level_1_key, $level_2_key);  # 要注意，这个数组，和 前面的 数组 @matrix_type=('MainKey_Lv1', 'MainKey_Lv1')一一对应
                        	for (my $key_i=0; $key_i<@keyArray; $key_i++){
                        		
                        		if ( defined ($mainKey_MatrixHash->{'matrix'}->{ $matrix_type[$key_i] }->{$domType}->{ $keyArray[$key_i] }->{$domID}->{$level_3_key}) ) {
              		      		  $mainKey_MatrixHash->{'matrix'}->{ $matrix_type[$key_i] }->{$domType}->{ $keyArray[$key_i] }->{$domID}->{$level_3_key}++;    		     
              		      	  }
              		      	  else{
              		      	  	$mainKey_MatrixHash->{'matrix'}->{ $matrix_type[$key_i] }->{$domType}->{ $keyArray[$key_i] }->{$domID}->{$level_3_key}=1;             
              		      	  }
              		      	  
              		          $mainKey_MatrixHash->{'DoMainHash'}->{$domType}->{$domID}->{'IDhyperLink' }=ExcelHandle::makeHyperLk($domHash->{$domType}->{$domID}->{'Url'       }, $domID);
              		      	  $mainKey_MatrixHash->{'DoMainHash'}->{$domType}->{$domID}->{'DEF'         }=                         $domHash->{$domType}->{$domID}->{'Definition'};
              		      	  $mainKey_MatrixHash->{'DoMainHash'}->{$domType}->{$domID}->{'DefHyperLink'}=ExcelHandle::makeHyperLk($domHash->{$domType}->{$domID}->{'Url'       }, $domHash->{$domType}->{$domID}->{'Definition'});
              		        
                        	}
                        	
                        	                        	
                        	
                        	#装填 格式数组hash
              		        
              		        
              		        $fmtKeyWordsHash->{'level_1_key_fmt'}->{$level_1_key}=1;
              		        $fmtKeyWordsHash->{'level_2_key_fmt'}->{$level_2_key}=1;
              		        $fmtKeyWordsHash->{'level_3_key_fmt'}->{$level_3_key}=1;
              		        $fmtKeyWordsHash->{'level_4_key_fmt'}->{$domType    }=1;
              		        $fmtKeyWordsHash->{'level_5_key_fmt'}->{$domID      }=1;
                        	
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
  
  # 构建 excel文件对象  构建 format hash
  my $outXlsBook= Excel::Writer::XLSX->new($out_excelFile);
  my $emptyHash;
  my $excelCellsFmtHash=ExcelHandle::MakeFormatHash_1($emptyHash,$fmtKeyWordsHash); 
  
  
  my $stNmHr='AllLv3'; 
  
  #对所有leve3，构建 全蛋白 Domain excel的hash，以及对应的 导航matrix的 定位HyperLink 的hash
  my $AllDmImfomExcelHash;
  ($AllDmImfomExcelHash,  $mainKey_MatrixHash->{'AllLv3_matrixLinkHash'}->{'AllLv3'} )
  =@{ &make3ImportantMatrix( $level1_ItScClu_Hash                ,   #1
        	                   $all_level1_Hash                    ,   #2
	                           $level_1_KeyOrderHash               ,   #3
	                           $all_level2_Hash                    ,   #4
	                           $level_2_KeyOrderHash               ,   #5
	                           $all_level3_Hash                    ,   #6
	                           $level_3_KeyOrderHash               ,   #7
	                           $multiLevelHash                     ,   #8
	                           $allIDHash_for_EachLevle1Type       ,   #9
	                           $stNmHr                             ,   #10
	                           $mainKey_MatrixHash                 ,   #11
	                           \@matrix_type                       ,   #12
	                           $ordHash                            ,   #13
                             $snmHash                                #14
	                          )
	 };
	
	#打印出包含所有Level 3的key的 excel sheet                       
	  
  my $OutSheetObj=$outXlsBook->add_worksheet($stNmHr);
  
  #warn "1\n";
  ExcelHandle::PtExcelBasedOnFmtHash($outXlsBook, $OutSheetObj, $AllDmImfomExcelHash ,$excelCellsFmtHash );
  	
  $OutSheetObj->set_column(0 ,0   ,4  );  # 第0    列宽度为   第0   
  $OutSheetObj->set_column(1 ,1   ,40 );
  $OutSheetObj->set_column(2 ,2   ,30 );
  $OutSheetObj->set_column(3 ,3   ,10 );
  $OutSheetObj->set_column(4 ,4   ,15 );
  $OutSheetObj->set_column(5 ,5   ,4  );
  $OutSheetObj->set_column(6 ,300 ,2  );
  
  
  #然后，对每一个 Level3 的key，再制作一个单独的 key的 excel sheet，以及对应的 导航matrix 的 hash
  if (   (  keys %{ $all_level3_Hash }  ) > 1   ){
  	
    if (  ( ref ( $all_level3_Hash ) ) eq 'HASH' ){
      foreach my $level_3_key (    sort { $level_3_KeyOrderHash->{$a} <=> $level_3_KeyOrderHash->{$b} } (   keys (  %{ $all_level3_Hash }  )   )    ){    
        
        #建立一个 只有单个 $level_3_key 的Hash
        my $oneLv3KeyHash; 
        $oneLv3KeyHash->{$level_3_key}=1;
        
        my $EcStNmHr=$level_3_key; 
        
        my $EachDmImfomExcelHash;
        ($EachDmImfomExcelHash,  $mainKey_MatrixHash->{'EachLv3_matrixLinkHash'}->{$level_3_key} )
        =@{ &make3ImportantMatrix( $level1_ItScClu_Hash                ,   #1
        	                         $all_level1_Hash                    ,   #2
	                                 $level_1_KeyOrderHash               ,   #3
	                                 $all_level2_Hash                    ,   #4
	                                 $level_2_KeyOrderHash               ,   #5
	                                 $oneLv3KeyHash                      ,   #6
	                                 $level_3_KeyOrderHash               ,   #7
	                                 $multiLevelHash                     ,   #8
	                                 $allIDHash_for_EachLevle1Type       ,   #9
	                                 $EcStNmHr                           ,   #10
	                                 $mainKey_MatrixHash                 ,   #11
	                                 \@matrix_type                       ,   #12
	                                 $ordHash                            ,   #13
                                   $snmHash                                #14
	                                )
	      };
	                            
	                           
          
        my $EcOutSheetObj=$outXlsBook->add_worksheet($EcStNmHr);
        #warn "2\n";
        ExcelHandle::PtExcelBasedOnFmtHash($outXlsBook, $EcOutSheetObj, $EachDmImfomExcelHash ,$excelCellsFmtHash );
  	
        $OutSheetObj->set_column(0 ,0   ,4  );  # 第0    列宽度为   第0   
        $OutSheetObj->set_column(1 ,1   ,40 );
        $OutSheetObj->set_column(2 ,2   ,30 );
        $OutSheetObj->set_column(3 ,3   ,10 );
        $OutSheetObj->set_column(4 ,4   ,15 );
        $OutSheetObj->set_column(5 ,5   ,4  );
        $OutSheetObj->set_column(6 ,300 ,2  );
        
        
      }
    } 
  
  
  } 
  
  
  for ( my $key_j=0; $key_j < @matrix_type; $key_j++ ){
    if (  ( ref ( $mainKey_MatrixHash->{'DoMainHash'} ) ) eq 'HASH' ){  #$mainKey_MatrixHash->{'DoMainHash'}->{$domType}->{$domID}->{'DefHyperLink'}=ExcelHandle::makeHyperLk($domHash->{$domType}->{$domID}->{'Url'       }, $domHash->{$domType}->{$domID}->{'Definition'});
      foreach my $domTpKey (    sort { $ordHash->{$a} <=>  $ordHash->{$b} } (   keys (  %{ $mainKey_MatrixHash->{'DoMainHash'} }  )   )    ){
      	my $inMatrixHash;
      	my $SttLineNumbe=0;
      	( $inMatrixHash, $SttLineNumbe )=@{ 
      	  	                                  &FromHash2Excel_domainMatrix( $mainKey_MatrixHash->{'matrix'}->{ $matrix_type[$key_j] }->{$domTpKey},                  #1      #$mainKey_MatrixHash->{'matrix'}->{ $matrix_type[$key_i] }->{$domType}->{ $keyArray[$key_i] }->{$domID}->{$level_3_key}
      	  	                                                                $mainKey_MatrixHash->{'AllLv3_matrixLinkHash'}->{'AllLv3'}->{$domTpKey},                 #2
      	  	                                                                $mainKey_MatrixHash->{'DoMainHash'}->{$domTpKey},                                        #3
      	  	                                                                $matrix_type_OrderHash->{ $matrix_type[$key_j] },                                        #4
      	  	                                                                $all_level3_Hash,                                                                        #5
      	  	                                                                $SttLineNumbe,                                                                           #6
      	  	                                                                $inMatrixHash                                                                            #7
      	  	                                                               ) 
      	  	                               };
        
      	if (    (   keys (  %{ $all_level3_Hash }  )   ) > 1    ){
      	  foreach my $level_3_key (    sort { $level_3_KeyOrderHash->{$a} <=> $level_3_KeyOrderHash->{$b} } (   keys (  %{ $all_level3_Hash }  )   )    ){  
      	    
      	    #建立一个 只有单个 $level_3_key 的Hash
            my $oneLv3KeyHash; 
            $oneLv3KeyHash->{$level_3_key}=1;
      	    
      	    ( $inMatrixHash, $SttLineNumbe )=@{ 
      	  	                                     &FromHash2Excel_domainMatrix( $mainKey_MatrixHash->{'matrix'}->{ $matrix_type[$key_j] }->{$domTpKey},       #1    #$mainKey_MatrixHash->{'matrix'}->{ $matrix_type[$key_i] }->{$domType}->{ $keyArray[$key_i] }->{$domID}->{$level_3_key}
      	  	                                                                   $mainKey_MatrixHash->{'EachLv3_matrixLinkHash'}->{$level_3_key}->{$domTpKey}, #2
      	  	                                                                   $mainKey_MatrixHash->{'DoMainHash'}->{$domTpKey},                             #3
      	  	                                                                   $matrix_type_OrderHash->{ $matrix_type[$key_j] },                             #4
      	  	                                                                   $oneLv3KeyHash,                                                               #5
      	  	                                                                   $SttLineNumbe,                                                                #6
      	  	                                                                   $inMatrixHash                                                                 #7
      	  	                                                                  ) 
      	  	                                  };
      	  
      	  
      	  }      	
      	}
      	my $TNmH_Key=$snmHash->{$domTpKey};
      	my $DmMtxShtNM="${TNmH_Key}_$key_j"; #warn "\$DmMtxShtNM=$DmMtxShtNM\n"; 
    	  my $DmMtxOutSheetObj=$outXlsBook->add_worksheet($DmMtxShtNM); 
    	  #warn "3\n";
    	  ExcelHandle::PtExcelBasedOnFmtHash($outXlsBook, $DmMtxOutSheetObj, $inMatrixHash ,$excelCellsFmtHash );
    	  $DmMtxOutSheetObj->set_column(0 ,0 ,30 );  # 第0    列宽度为   第0   
        $DmMtxOutSheetObj->set_column(1 ,150 ,2);  # 第1    列宽度为   第1 
      	
      }
    }
  }
  
  
  
  
#  
  
}


sub make3ImportantMatrix{
	my (
	  $In_level1_ItScClu_Hash                ,   #1
	  $In_all_level1_Hash                    ,   #2
	  $In_level_1_KeyOrderHash               ,   #3
	  $In_all_level2_Hash                    ,   #4
	  $In_level_2_KeyOrderHash               ,   #5
	  $In_all_level3_Hash                    ,   #6
	  $In_level_3_KeyOrderHash               ,   #7
	  $In_multiLevelHash                     ,   #8
	  $In_allIDHash_for_EachLevle1Type       ,   #9
	  $In_DomInformSheetName                 ,   #10
	  $In_mainKey_MatrixHash                 ,   #11
	  $In_mainkeyTypeArray                   ,   #12
	  $In_ordHash                            ,   #13
    $In_snmHash                                #14
	  
	  
	)=@_;
	
  my $Out_AllDmImfomExcelHash;   #这个HASH用来　装用来打印　各个蛋白的　Ｄｏｍａｉｎ信息的　Ｅｘｃｅｌ的HASH
  my $out_MainKey_matrixHash;
  
  
  my $col_i=0;
  if (  ( ref ( $In_all_level1_Hash ) ) eq 'HASH' ){
    foreach my $level_1_key (    sort { $In_level_1_KeyOrderHash->{$a} <=> $In_level_1_KeyOrderHash->{$b} } (   keys (  %{ $In_all_level1_Hash }  )   )    ){
    	
      if (  ( ref ( $In_all_level2_Hash ) ) eq 'HASH' ){
        foreach my $level_2_key (    sort { $In_level_2_KeyOrderHash->{$a} <=> $In_level_2_KeyOrderHash->{$b} } (   keys (  %{ $In_all_level2_Hash }  )   )    ){ 
    
          if (  ( ref ( $In_all_level3_Hash ) ) eq 'HASH' ){
            foreach my $level_3_key (    sort { $In_level_3_KeyOrderHash->{$a} <=> $In_level_3_KeyOrderHash->{$b} } (   keys (  %{ $In_all_level3_Hash }  )   )    ){    
                
                
              if (  ( ref ( $In_multiLevelHash->{$level_1_key}->{$level_2_key}->{$level_3_key} ) ) eq 'HASH' ){  
                foreach my $SequnceID (    sort { $In_level1_ItScClu_Hash->{$level_1_key}->{'lv1_CluSeqOrdr'}->{'id_seqHash'}->{$a}->{'seqOrder'} <=> $In_level1_ItScClu_Hash->{$level_1_key}->{'lv1_CluSeqOrdr'}->{'id_seqHash'}->{$b}->{'seqOrder'} } (   keys (  %{ $In_multiLevelHash->{$level_1_key}->{$level_2_key}->{$level_3_key} }  )   )    ){ 
                
                  
                  my $j=0;                                                                                                                                                                                   #$level1_ItScClu_Hash->{$level_1_key1}->{'lv1_seqFile_lk'}=$lv1_seqFile_lk;
  				                                                                                                                                                                                                   #$level1_ItScClu_Hash->{$level_1_key1}->{'lv1_cluFile_lk'}=$lv1_cluFile_lk;
  				        push @{ $Out_AllDmImfomExcelHash->{ $level_1_key } }, [ $col_i,$j, ExcelHandle::makeHyperLk( $In_level1_ItScClu_Hash->{$level_1_key}->{'lv1_seqFile_lk'}, $SequnceID   ) ]; $j++;              #                                                                          
  				        push @{ $Out_AllDmImfomExcelHash->{ $level_2_key } }, [ $col_i,$j, ExcelHandle::makeHyperLk( $In_level1_ItScClu_Hash->{$level_1_key}->{'lv1_cluFile_lk'}, $level_1_key ) ]; $j++;              #$level1_ItScClu_Hash->{$level_1_key1}->{'lv1_ISHRead_lk'}=$lv1_ISHRead_lk;
  				        push @{ $Out_AllDmImfomExcelHash->{ $level_1_key } }, [ $col_i,$j, ExcelHandle::makeHyperLk( $In_level1_ItScClu_Hash->{$level_1_key}->{'lv1_ISHRead_lk'}, $level_2_key ) ]; $j++;              #$level1_ItScClu_Hash->{$level_1_key1}->{'lv1_ItScDir_lk'}=$lv1_ItScDir_lk; 
  				        push @{ $Out_AllDmImfomExcelHash->{ $level_3_key } }, [ $col_i,$j, ExcelHandle::makeHyperLk( $In_level1_ItScClu_Hash->{$level_1_key}->{'lv1_ItScDir_lk'}, $level_3_key ) ]; $j++;              #$level1_ItScClu_Hash->{$level_1_key1}->{'Lv1_ItClPng_lk'}=$Lv1_ItClPng_lk;
  				        push @{ $Out_AllDmImfomExcelHash->{ $level_1_key } }, [ $col_i,$j, ExcelHandle::makeHyperLk( $In_level1_ItScClu_Hash->{$level_1_key}->{'Lv1_ItClPng_lk'}, $SequnceID   ) ]; $j++;              #ExcelHandle::makeHyperLk($domHash->{$domType}->{$domID}->{'Url'       }, $domHash->{$domType}->{$domID}->{'Definition'});
              		        
  				        
  				        $j++;
  				        
  				        #print "\$Out_AllDmImfomExcelHash=$Out_AllDmImfomExcelHash\n";
  				        
  				        if (  ( ref ( $In_allIDHash_for_EachLevle1Type->{$level_1_key}->{$level_3_key} ) ) eq 'HASH' ){
  				          foreach my $eachDomType(   sort { $In_ordHash->{$a} <=> $In_ordHash->{$b} }(   keys(  %{ $In_allIDHash_for_EachLevle1Type->{$level_1_key}->{$level_3_key} }  )   )    ){     #warn "\$eachDomType=$eachDomType\n";
  				            push @{ $Out_AllDmImfomExcelHash->{ 'BlackRed' } }, [$col_i,$j, "$In_snmHash->{$eachDomType}:"  ] ; $j++;
  				            push @{ $Out_AllDmImfomExcelHash->{ 'BlackRed' } }, [$col_i,$j, ''                              ] ; $j++;
  				            
  				            if (  ( ref ( $In_allIDHash_for_EachLevle1Type->{$level_1_key}->{$level_3_key}->{$eachDomType} ) ) eq 'HASH' ){
  				              foreach my $eachDomID(   sort {$a cmp $b}(   keys(  %{ $In_allIDHash_for_EachLevle1Type->{$level_1_key}->{$level_3_key}->{$eachDomType} }  )   )    ){  				    #warn "\$eachDomID=$eachDomID\n";	
  				              	
  				              	
  				              	#$mainKey_MatrixHash->{'DoHavePointMatrix'}->{$level_1_key}->{$level_2_key}->{$level_3_key}->{$SequnceID}->{$domType}->{$domID}=1
  				              	if (  defined ( $In_mainKey_MatrixHash->{'DoHavePointMatrix'}->{$level_1_key}->{$level_2_key}->{$level_3_key}->{$SequnceID}->{$eachDomType}->{$eachDomID} )  ){
  				                #if (  defined ( $In_allIDHash_for_EachLevle1Type->{$level_1_key}->{$level_3_key}->{$eachDomType}->{$eachDomID} )  ){
  				                  push @{ $Out_AllDmImfomExcelHash->{ $eachDomID } }, [$col_i,$j, $In_mainKey_MatrixHash->{'DoMainHash'}->{$eachDomType}->{$eachDomID}->{'DefHyperLink'}  ] ; 
  				                  #warn "\$In_mainKey_MatrixHash->{'DoMainHash'}->{$eachDomType}->{$eachDomID}->{'DefHyperLink'}=$In_mainKey_MatrixHash->{'DoMainHash'}->{$eachDomType}->{$eachDomID}->{'DefHyperLink'}\n";
  				                 
  				                  my $cellAddrr=Excel::Writer::XLSX::Utility::xl_rowcol_to_cell($col_i,$j);   
  				                  my $SheetLink="#$In_DomInformSheetName!$cellAddrr";             
  				                                                                                                                                #In_mainkeyTypeArray     @matrix_type=('MainKey_Lv1', 'MainKey_Lv1');    
  				                                                                                                                                
  				                                                                                                                                
  				                                                                                                                                
  				                  my @keyArray   = ($level_1_key, $level_2_key);  # 要注意，这个数组，和 前面的 数组 @matrix_type=('MainKey_Lv1', 'MainKey_Lv1')一一对应                                                                                                              
                            my @matrix_type= @{ $In_mainkeyTypeArray };     # 要注意，这个数组，和 前面的 数组 @matrix_type=('MainKey_Lv1', 'MainKey_Lv1')一一对应
                        	  for (my $key_i=0; $key_i<@keyArray; $key_i++){               		      	    
              		      	    $out_MainKey_matrixHash->{$eachDomType}->{ $keyArray[$key_i] }->{$eachDomID}->{$level_3_key}
              		      	    =ExcelHandle::makeHyperLk( $SheetLink, $In_mainKey_MatrixHash->{'matrix'}->{ $matrix_type[$key_i] }->{$eachDomType}->{ $keyArray[$key_i] }->{$eachDomID}->{$level_3_key} );      
              		      	    print "\$out_MainKey_matrixHash->{$eachDomType}->{ $keyArray[$key_i] }->{$eachDomID}->{$level_3_key}=$out_MainKey_matrixHash->{$eachDomType}->{ $keyArray[$key_i] }->{$eachDomID}->{$level_3_key}\n\n";        		          
                        	  }
                        	
                        	 
  				                  
  				                }  
  				                $j++;
  				              }
  				            }
  				            $j++;
  				            $j++;
  				            
  				          }
  				        }
  				        
  				        $col_i++;
    		          
  			        
                  
                  
                  
                                             	
                }
                
              }
            }
            
          }          
        }
        
        
          $col_i++;                                             #行数，行号
      } 
        
      
      
    }
  }
  return [$Out_AllDmImfomExcelHash, $out_MainKey_matrixHash];
}





sub GetFmlDomGoKeggHash{
  my ($interProscanHash, $ProteinID)=@_;
  
 
  
  
  if (  ( ref ( $interProscanHash->{$ProteinID} ) ) eq 'HASH' ){
    
    my $outHash;
    
    
    if (  ( ref ( $interProscanHash->{$ProteinID}->{'Protein family membership'} ) ) eq 'HASH' ){
      foreach my $fmlID (    sort { $a cmp $b } (   keys (  %{ $interProscanHash->{$ProteinID}->{'Protein family membership'} }  )   )    ){       
        $outHash->{'realHash'}->{'Protein family membership'}->{$fmlID}->{'Definition'}=$interProscanHash->{$ProteinID}->{'Protein family membership'}->{$fmlID}->{'Definition'};
        $outHash->{'realHash'}->{'Protein family membership'}->{$fmlID}->{'Url'       }=$interProscanHash->{$ProteinID}->{'Protein family membership'}->{$fmlID}->{'Url'       };      
      }
    }
    
    if (  ( ref ( $interProscanHash->{$ProteinID}->{'Domains and repeats'}->{'DomainLineNBhash'} ) ) eq 'HASH' ){
      foreach my $number (    sort { $a cmp $b } (   keys (  %{ $interProscanHash->{$ProteinID}->{'Domains and repeats'}->{'DomainLineNBhash'} }  )   )    ){     
      	if (  ( ref ( $interProscanHash->{$ProteinID}->{'Domains and repeats'}->{'DomainLineNBhash'}->{$number}->{'DomainIDhash'} ) ) eq 'HASH' ){
          foreach my $DomID (    sort { $a cmp $b } (   keys (  %{ $interProscanHash->{$ProteinID}->{'Domains and repeats'}->{'DomainLineNBhash'}->{$number}->{'DomainIDhash'} }  )   )    ){  
            $outHash->{'realHash'}->{'Domains and repeats'}->{$DomID}->{'Definition'}=$interProscanHash->{$ProteinID}->{'Domains and repeats'}->{'DomainLineNBhash'}->{$number}->{'DomainIDhash'}->{$DomID}->{'Definition'};
            $outHash->{'realHash'}->{'Domains and repeats'}->{$DomID}->{'Url'       }=$interProscanHash->{$ProteinID}->{'Domains and repeats'}->{'DomainLineNBhash'}->{$number}->{'DomainIDhash'}->{$DomID}->{'Url'       };  
          }
        }          
      }
    }
    
    if (  ( ref ( $interProscanHash->{$ProteinID}->{'GO term prediction'} ) ) eq 'HASH' ){
      foreach my $GOtype (    sort { $a cmp $b } (   keys (  %{ $interProscanHash->{$ProteinID}->{'GO term prediction'} }  )   )    ){       
      	
      	if (  ( ref ( $interProscanHash->{$ProteinID}->{'GO term prediction'}->{$GOtype} ) ) eq 'HASH' ){
          foreach my $GOid (    sort { $a cmp $b } (   keys (  %{ $interProscanHash->{$ProteinID}->{'GO term prediction'}->{$GOtype} }  )   )    ){  
            $outHash->{'realHash'}->{$GOtype}->{$GOid}->{'Definition'}=$interProscanHash->{$ProteinID}->{'GO term prediction'}->{$GOtype}->{$GOid}->{'Definition'};
            $outHash->{'realHash'}->{$GOtype}->{$GOid}->{'Url'       }=$interProscanHash->{$ProteinID}->{'GO term prediction'}->{$GOtype}->{$GOid}->{'Url'       }; 
          }
        }     
      }
    }
    
    if (  ( ref ( $interProscanHash->{$ProteinID}->{'Pathway'} ) ) eq 'HASH' ){
      foreach my $ptwID (    sort { $a cmp $b } (   keys (  %{ $interProscanHash->{$ProteinID}->{'Pathway'} }  )   )    ){       
        $outHash->{'realHash'}->{'Pathway'}->{$ptwID}->{'Definition'}=$interProscanHash->{$ProteinID}->{'Pathway'}->{$ptwID}->{'Definition'    };
        $outHash->{'realHash'}->{'Pathway'}->{$ptwID}->{'Url'       }=$interProscanHash->{$ProteinID}->{'Pathway'}->{$ptwID}->{'DomainDataBase'};      
      }
    }
    
    return $outHash;
    
  }
  else {
    die "Die here: \$interProscanHash=$interProscanHash, \$ProteinID=$ProteinID \n\n";
  }
}


#&FromHash2Excel_domainMatrix(                                                 $mainKey_MatrixHash->{'matrix'}->{ $matrix_type[$key_j] },                        #$mainKey_MatrixHash->{'matrix'}->{ $matrix_type[$key_i] }->{$domType}->{ $keyArray[$key_i] }->{$domID}->{$level_3_key}
#      	  	                                                                   $mainKey_MatrixHash->{'EachLv3_matrixLinkHash'}->{$level_3_key}, 
#      	  	                                                                   $mainKey_MatrixHash->{'DoMainHash'}, 
#      	  	                                                                   $matrix_type_OrderHash->{ $matrix_type[$key_j] },
#      	  	                                                                   $oneLv3KeyHash
#      	  	                                                                   $SttLineNumbe,
#      	  	                                                                   $inMatrixHash                                         
#      	  	                                                                  ) 



sub FromHash2Excel_domainMatrix{
  my ($inHash,           #1
      $inLinkHash,       #2
      $inDmHASH,         #3
      $IndexHash,        #4
      $lv3_hash,         #5
      $stLineNb,         #6
      $proMatrixHash     #7
      )=@_;
      
  
  #DirFileHandle::PrintDumper('IndexHash.txt',$IndexHash)  ; 
  #DirFileHandle::PrintDumper('inHash.txt',$inHash)  ; 
  
  my $AllLevel_2_KeyHash;
  my $AllLevel_3_KeyHash;  ##################
  my $level_1_sum_hash;
  my $level_2_sum_hash;
  
  foreach my $Level_1_key (   sort { $IndexHash->{$a} <=> $IndexHash->{$b} } (   keys (  %{ $inHash }  )   )    ) {
  	foreach my $Level_2_key (   sort {$a cmp $b} (   keys (  %{ $inHash->{$Level_1_key} }  )   )    ) {
  		$AllLevel_2_KeyHash->{$Level_2_key}=1;   print "MMMMMMMMMMM000011 \$AllLevel_2_KeyHash->{\$Level_2_key}=\$AllLevel_2_KeyHash->{$Level_2_key}=$AllLevel_2_KeyHash->{$Level_2_key}\n";
  		foreach my $Level_3_key (   sort {$a cmp $b} (   keys (  %{ $lv3_hash }  )   )    ) {
  			$AllLevel_3_KeyHash->{$Level_3_key}=1;   print "MMMMMMMMMMM002222 \$AllLevel_3_KeyHash->{$Level_3_key}=\$AllLevel_3_KeyHash->{$Level_3_key}=$AllLevel_3_KeyHash->{$Level_3_key}\n";
  			if (  defined ( $level_1_sum_hash->{$Level_1_key}->{$Level_3_key}->{'allSum'} )  ) { $level_1_sum_hash->{$Level_1_key}->{$Level_3_key}->{'allSum'}+=$inHash->{$Level_1_key}->{$Level_2_key}->{$Level_3_key}; }
  			else                                                                               { $level_1_sum_hash->{$Level_1_key}->{$Level_3_key}->{'allSum'}=$inHash->{$Level_1_key}->{$Level_2_key}->{$Level_3_key}; } 
  			if (  defined ( $level_1_sum_hash->{$Level_1_key}->{$Level_3_key}->{'kindSum'} )  ){ $level_1_sum_hash->{$Level_1_key}->{$Level_3_key}->{'kindSum'}++; }
  			else                                                                               { $level_1_sum_hash->{$Level_1_key}->{$Level_3_key}->{'kindSum'}=1; } 
  			if (  defined ( $level_2_sum_hash->{$Level_2_key}->{$Level_3_key}->{'allSum'} )  ) { $level_2_sum_hash->{$Level_2_key}->{$Level_3_key}->{'allSum'}+=$inHash->{$Level_1_key}->{$Level_2_key}->{$Level_3_key}; }
  			else                                                                               { $level_2_sum_hash->{$Level_2_key}->{$Level_3_key}->{'allSum'}=$inHash->{$Level_1_key}->{$Level_2_key}->{$Level_3_key}; } 
  			if (  defined ( $level_2_sum_hash->{$Level_2_key}->{$Level_3_key}->{'kindSum'} )  ){ $level_2_sum_hash->{$Level_2_key}->{$Level_3_key}->{'kindSum'}++; }
  			else                                                                               { $level_2_sum_hash->{$Level_2_key}->{$Level_3_key}->{'kindSum'}=1; }                
  		}
  	}
  }
  
  my $i=0;
  if   ($stLineNb==0) {}
  else { $i=$stLineNb+1;}
  my $outPutHash;       
   #warn "\$stLineNb=$stLineNb\n";
   #warn "\$proMatrixHash=$proMatrixHash\n";
 
  
  my $blankColNumber=1;  # 统计数据和matrix之间的列数
  #$colStartNb是matrix开始的列数  #第一列是 行标识符， 第二列开始 是各个Lv3的统计数据， 从$colStartNb列 开始 才是matrix的数据
  my $colStartNb=(   keys (  %{ $AllLevel_3_KeyHash }  )   )*2+1+$blankColNumber;   #=5;  ###   
  
  # 1st Line  打印出 列标识符的ID和超链接
  my $j=$colStartNb;
  foreach my $Level_2_key (   sort {$a cmp $b} (   keys (  %{ $AllLevel_2_KeyHash }  )   )    ) {
  	my $fisrtJ=$j;
   	my $k=0;
   	foreach my $Level_3_key (   sort {$a cmp $b} (   keys (  %{ $AllLevel_3_KeyHash }  )   )    ) {
  		$j=$j+$k;	  
  	  $k++;
  	}                                                                               #这个是 matrix导航的 列标识符的 内容
  	print "MMMMMMMMMMM110 push \@{  \$outPutHash->{ $Level_2_key }  },[$i,$fisrtJ,\$inDmHASH->{$Level_2_key}->{'IDhyperLink'}=$inDmHASH->{$Level_2_key}->{'IDhyperLink'}, $i, $j ]\n";  
  	push @{  $outPutHash->{ $Level_2_key }  },[$i,$fisrtJ,$inDmHASH->{$Level_2_key}->{'IDhyperLink'},$i,$j ];  
  	
  	$j++;
  }$i++;
  
  # 2nd Line  打印出 列标识符的定义
  $j=$colStartNb;
  foreach my $Level_2_key (   sort {$a cmp $b} (   keys (  %{ $AllLevel_2_KeyHash }  )   )    ) {
  	my $fisrtJ=$j;
  	my $k=0;
  	foreach my $Level_3_key (   sort {$a cmp $b} (   keys (  %{ $AllLevel_3_KeyHash }  )   )    ) {
  		$j=$j+$k;	  
  	  $k++;
  	}
  	push @{  $outPutHash->{ $Level_2_key }  },[$i,$fisrtJ,$inDmHASH->{$Level_2_key}->{'DEF'},$i,$j ];  
  	print "MMMMMMMMMMM111 push \@{  \$outPutHash->{ $Level_2_key }  },[1,$fisrtJ,$j,$inDmHASH->{$Level_2_key}->{'DEF'},$i,$j ]\n";

  	$j++;
  }$i++;
  
  # 3rd Line 打印出 lv3 的标示
  $j=$colStartNb;
  foreach my $Level_2_key (   sort {$a cmp $b} (   keys (  %{ $AllLevel_2_KeyHash }  )   )    ) {
  	my $k=0;
  	foreach my $Level_3_key (   sort {$a cmp $b} (   keys (  %{ $AllLevel_3_KeyHash }  )   )    ) {
  		$j=$j+$k;
  		push @{  $outPutHash->{ $Level_3_key }  },[$i,$j,$Level_3_key ];  print "MMMMMMMMMMM220 push \@{  \$outPutHash->{ $Level_3_key }  },[$i,$j,$Level_3_key ]\n";  	  
  	  $k++;
  	}
  	$j++;
  }$i++;

  # 4th Line  打印出 每一列的 覆盖数
  $j=$colStartNb;
  foreach my $Level_2_key (   sort {$a cmp $b} (   keys (  %{ $AllLevel_2_KeyHash }  )   )    ) {
  	my $k=0;
  	foreach my $Level_3_key (   sort {$a cmp $b} (   keys (  %{ $AllLevel_3_KeyHash }  )   )    ) {
  		$j=$j+$k;
  		if (  defined ( $level_2_sum_hash->{$Level_2_key}->{$Level_3_key}->{'kindSum'} )  ){
  		  push @{  $outPutHash->{ $Level_3_key }  },[$i,$j,$level_2_sum_hash->{$Level_2_key}->{$Level_3_key}->{'kindSum'} ];  print "MMMMMMMMMMM221 push \@{  \$outPutHash->{ $Level_3_key }  },[$i,$j,$level_2_sum_hash->{$Level_2_key}->{$Level_3_key}->{'kindSum'} ]\n";
  	  }
  	  $k++;
  	}
  	$j++;
  }$i++; 
  
  # 5th Line  打印出 每一列的 总数
  $j=$colStartNb;
  foreach my $Level_2_key (   sort {$a cmp $b} (   keys (  %{ $AllLevel_2_KeyHash }  )   )    ) {
  	my $k=0;
  	foreach my $Level_3_key (   sort {$a cmp $b} (   keys (  %{ $AllLevel_3_KeyHash }  )   )    ) {
  		$j=$j+$k;
      if (  defined ( $level_2_sum_hash->{$Level_2_key}->{$Level_3_key}->{'allSum'} )  ){
  		  push @{  $outPutHash->{ $Level_3_key }  },[$i,$j,$level_2_sum_hash->{$Level_2_key}->{$Level_3_key}->{'allSum'} ];  print "MMMMMMMMMMM222 push \@{  \$outPutHash->{ $Level_3_key }  },[$i,$j,$level_2_sum_hash->{$Level_2_key}->{$Level_3_key}->{'allSum'} ]\n";
  	  }  	  
  	  $k++;
  	}
  	$j++;
  }$i++;
  
  # 6th Line 开始打matrix
  
  foreach my $Level_1_key (   sort { $IndexHash->{$a}  <=> $IndexHash->{$b} } (   keys (  %{ $inHash }  )   )    ) {
  	$j=0;
  	# 1st col  打印出 该行的 第一列 行标识符
  	push @{  $outPutHash->{ $Level_1_key }  },[$i,$j, $Level_1_key ]; $j++;
  	# 2nd col  打印出 该行的 第二列到第？列 各个lv3 对应的 行 统计数据
  	foreach my $Level_3_key (   sort {$a cmp $b} (   keys (  %{ $AllLevel_3_KeyHash }  )   )    ) {
  		if (  defined ( $level_1_sum_hash->{$Level_1_key}->{$Level_3_key}->{'kindSum'} )  ){
  		  push @{  $outPutHash->{ $Level_3_key }  },[$i,$j,$level_1_sum_hash->{$Level_1_key}->{$Level_3_key}->{'kindSum'} ]; 
  	  }
  	  $j++;
  	}
  	# ？th col 打印出 该行的 第二列到第？列 各个lv3 对应的 行 统计数据
  	foreach my $Level_3_key (   sort {$a cmp $b} (   keys (  %{ $AllLevel_3_KeyHash }  )   )    ) {
  		if (  defined ( $level_1_sum_hash->{$Level_1_key}->{$Level_3_key}->{'allSum'} )  ){
  		  push @{  $outPutHash->{ $Level_3_key }  },[$i,$j,$level_1_sum_hash->{$Level_1_key}->{$Level_3_key}->{'allSum'} ]; 
  	  }
  	  $j++;
  	}
  	
  	
  	# $colStartNb col， matrix开始列
  	$j=$colStartNb;
  	foreach my $Level_2_key (   sort {$a cmp $b} (   keys (  %{ $AllLevel_2_KeyHash }  )   )    ) {
  		my $k=0;
  		foreach my $Level_3_key (   sort {$a cmp $b} (   keys (  %{ $AllLevel_3_KeyHash }  )   )    ) {
  			$j=$j+$k;
  			if ( defined ($inHash->{$Level_1_key}->{$Level_2_key}->{$Level_3_key}) ){
  				#$SpeciseDomainMatrixHash->{$eachDomType}->{'matrixLinkHash'}->{$eachSpHere}->{$eachDomID}->{$eachCUHere}->{$ecTp} $fmlDomainMatrixHash->{$eachDomType}->{'matrixLinkHash'}->{$eachPtFml}->{$eachDomID}->{$eachCUHere}->{$ecTp}
  			  push @{  $outPutHash->{ $Level_3_key }  },[$i,$j, $inLinkHash->{$Level_1_key}->{$Level_2_key}->{$Level_3_key} ];  #  $SpeciseDomainMatrixHash->{$eachDomType}->{'matrixLinkHash'}->{$eachDomID}->{$ecTp}->{'DefHyperLink'}
  			  print "MMMMMMMMMMM33 push \@{  \$outPutHash->{ $Level_3_key }  } ,[$i,$j, \$inLinkHash->{$Level_1_key}->{$Level_2_key}->{$Level_3_key}=$inLinkHash->{$Level_1_key}->{$Level_2_key}->{$Level_3_key} ]\n";
  		  }
  		  $k++;
  		}
  		$j++;
  	}
  	$i++;
  }
  my $outPut=[$outPutHash, $i];  print "\$outPut=$outPut\n\n";
  return $outPut;			   ##################
} 







1;