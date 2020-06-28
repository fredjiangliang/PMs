
#!/usr/bin/perl -w

BEGIN{  push (@INC, '/home/fredjiang/fredProgs/PMs')  };

use strict;
use warnings;


use Bio::Graphics;
use Bio::SearchIO;
use Bio::SearchIO::Writer::HTMLResultWriter;
use Storable ;
use Storable qw(nstore store_fd nstore_fd freeze thaw dclone);

use SeqSegmentsTools;
use ArrayHashChange;
use FastaFileHandle;
use MatrixCsvChange;
use BLOSUMalignSCORE;
use String_Work;

package  BlastHandle;


######################################################################################################################################################
#    
#    sub DrawBlast{      #pmAble#   #绘制Blast 的图像
#    
#    sub BioPerlBlastPraser20170325{       #pmAble#   #解析Blast文件
#    
#    sub DrawBlast2seq_20170330{     #pmAble#   #绘制Blast2seq的图像
#    
#    sub BioPerlBlastPraser_bl2seq_20170330{ #  #pmAble#    #解析blast 这个用于解析 bl2seq 的函数
#    
#    sub overLayerHadle{    #pmAble#  #BioPerlBlastPraser_bl2seq_20170330的子函数 处理overlayyer的hash
#    
#    sub GetCuttedLength{  #pmAble#  #BioPerlBlastPraser_bl2seq_20170330的子函数 获得 cutted的 片段的长度
#    
#    sub getWholeMatchedStartEnd{   #pmAble#  #BioPerlBlastPraser_bl2seq_20170330的子函数 获得 match的片段的 最开头 和 最结尾
#    
#    
#    
#######################################################################################################################################################################


#my $outStrst=BlastHandle::FastaCMD_blastdbCMD_runSub($dbFile, $entryNm,[  $blastVersion, $strand, $sttNb, $endNb]);
sub FastaCMD_blastdbCMD_runSub{
	
	my ($dbFile, $entryNm, $blastVersion, $strand, $sttNb, $endNb)=@_;
	
	my $warnMsgBody="\nIn package  BlastHandle,\tIn sub FastaCMD_blastdbCMD_runSub,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  my $subCalInfom=DirFileHandle::print_SubCallerInform;
  
  DieWork::Check_FileDirExist_or_DIE   ( $dbFile,  "\$dbFile",  $die_MsgHead, $subCalInfom  );
  DieWork::Check_DfdNoEmptString_or_DIE( $entryNm, "\$entryNm", $die_MsgHead, $subCalInfom  );
 
  
  my $blastVersion_new_or_not=0;
  if ( DieWork::Check_DfdNoEmptString_or_NOT( $blastVersion ) ){
  	$blastVersion=lc $blastVersion;   	$blastVersion=~s/\s//g;
  	if   (  ( $blastVersion eq 'blastall' ) || ( $blastVersion eq 'fastacmd'    ) || ( $blastVersion eq 'old' )  ){
  		$blastVersion_new_or_not=0;
  	}
  	elsif(  ( $blastVersion eq 'blastplus' ) || ( $blastVersion eq 'blastdbcmd' ) || ( $blastVersion eq 'new' )  ){
  		$blastVersion_new_or_not=1;
  	}
  	else{
  		DieWork::Just_dieWork( $die_MsgHead."\n\$blastVersion=$blastVersion should be blastall, blastplus, fastacmd, blastdbcmd, old, new  $!".$subCalInfom ); 
  	}  
  }
  else{
  	$blastVersion_new_or_not=1;
  }
  
  
  
  my $cmd_exeFile;
  if    ( $blastVersion_new_or_not == 0 ){ $cmd_exeFile="fastacmd";    }
  elsif ( $blastVersion_new_or_not == 1 ){ $cmd_exeFile="blastdbcmd";  }
  
  my $database_path;
  if    ( $blastVersion_new_or_not == 0 ){ $database_path=" -d  $dbFile ";  }
  elsif ( $blastVersion_new_or_not == 1 ){ $database_path=" -db $dbFile ";  }
  
  my $real_entryNm;
  if    ( $blastVersion_new_or_not == 0 ){ $real_entryNm=    " -s $entryNm ";  }
  elsif ( $blastVersion_new_or_not == 1 ){ $real_entryNm=" -entry $entryNm ";  }
  
  my $realStrandPos_inform='';
  if ( DieWork::Check_DfdNoEmptString_or_NOT( $strand ) ){
  	if (  ( $strand eq '+' ) || ( $strand eq '-' )  ){
  		DieWork::Check_DfdNoEmptString_or_NOT( $sttNb, "\$sttNb", $die_MsgHead, $subCalInfom  );
  		DieWork::Check_DfdNoEmptString_or_NOT( $endNb, "\$endNb", $die_MsgHead, $subCalInfom  );
  		if (  ( $sttNb > 0 ) && ( $endNb > 0 ) && ( $sttNb <= $endNb)  ){
  			if    (  ( $blastVersion_new_or_not == 0 ) && ( $strand eq '+' )   ){  
  				$realStrandPos_inform=" -S 1 -L $sttNb,$endNb ";  
  			}
  			elsif (  ( $blastVersion_new_or_not == 0 ) && ( $strand eq '-' )   ){  
  				$realStrandPos_inform=" -S 2 -L $sttNb,$endNb ";  
  			}
  			elsif (  ( $blastVersion_new_or_not == 1 ) && ( $strand eq '+' )   ){ 
  				$realStrandPos_inform=" -strand plus -range $sttNb-$endNb "; 
  			}
  			elsif (  ( $blastVersion_new_or_not == 1 ) && ( $strand eq '-' )   ){ 
  				$realStrandPos_inform=" -strand minus -range $sttNb-$endNb "; 
  			}
  		}
  		else{
  			DieWork::Just_dieWork( $die_MsgHead."\n\$sttNb=$sttNb and $endNb=$endNb should be numbers >  0, and  \$sttNb=$sttNb should <= $endNb=$endNb $!".$subCalInfom ); 
  		}
  		
  	}
  	else{
  		DieWork::Just_dieWork( $die_MsgHead."\n\$strand=$strand should be + or - $!".$subCalInfom ); 
  	}
  }
  else{
  	$realStrandPos_inform='';
  }
  
  my $runingCMD=$cmd_exeFile.$database_path.$real_entryNm.$realStrandPos_inform;
  my $HmReadTimeNow=TimeWork::GetHumanRead_NOW_Time();
  DieWork::Print_and_warn( "\n $HmReadTimeNow \$runingCMD=$runingCMD\n\n" );
  
  my $runingResult=`$runingCMD`;
  
  return $runingResult;
	
}

sub FormatDB_withORout_Plus_version{  #   BlastHandle::FormatDB_withORout_Plus_version ($inFasta_to_format, $DNAorPEP, $outDB, $blastVersion);
	my ($inFasta_to_format, $DNAorPEP, $outDB, $blastVersion)=@_;
	                                                                    DieWork::Print_and_warn( "\n  \$blastVersion=$blastVersion\n\n" );
	my $warnMsgBody="\nIn package  BlastHandle,\tIn sub FormatDB_withORout_Plus_version,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";  
  my $subCalInfom=DirFileHandle::print_SubCallerInform;
  
  DieWork::Check_FileDirExist_or_DIE   ( $inFasta_to_format,  "\$inFasta_to_format",  $die_MsgHead, $subCalInfom  );
  DieWork::Check_DfdNoEmptString_or_DIE( $DNAorPEP,           "\$DNAorPEP",           $DNAorPEP,    $DNAorPEP     );
  
  my $blastVersion_new_or_not=0;
  if ( DieWork::Check_DfdNoEmptString_or_NOT( $blastVersion ) ){
  	$blastVersion=lc $blastVersion;   	$blastVersion=~s/\s//g;
  	if   (  ( $blastVersion eq 'blastall' ) || ( $blastVersion eq 'fastacmd'    ) || ( $blastVersion eq 'old' )  ){
  		$blastVersion_new_or_not=0;
  	}
  	elsif(  ( $blastVersion eq 'blastplus' ) || ( $blastVersion eq 'blastdbcmd' ) || ( $blastVersion eq 'new' )  ){
  		$blastVersion_new_or_not=1;
  	}
  	else{
  		DieWork::Just_dieWork( $die_MsgHead."\n\$blastVersion=$blastVersion should be blastall, blastplus, fastacmd, blastdbcmd, old, new  $!".$subCalInfom ); 
  	}  
  }
  else{
  	$blastVersion_new_or_not=1;
  }
  
  my $cmd_exeFile;
  if    ( $blastVersion_new_or_not == 0 ){ $cmd_exeFile="formatdb";     }
  elsif ( $blastVersion_new_or_not == 1 ){ $cmd_exeFile="makeblastdb";  }
  
  my $database_path;
  if    ( $blastVersion_new_or_not == 0 ){ $database_path=" -i  $inFasta_to_format ";  }
  elsif ( $blastVersion_new_or_not == 1 ){ $database_path=" -in $inFasta_to_format ";  }
   
  
  my $DNA_PEP_option;  $DNAorPEP=lc $DNAorPEP;   	$DNAorPEP=~s/\s//g;
  if    (  ( $DNAorPEP eq 'd') || ( $DNAorPEP eq 'dna') || ( $DNAorPEP eq 'nucl') || ( $DNAorPEP eq 'nucleotide')  ){
  	if    ( $blastVersion_new_or_not == 0 ){ $DNA_PEP_option=" -p F ";          }
    elsif ( $blastVersion_new_or_not == 1 ){ $DNA_PEP_option=" -dbtype nucl ";  }
  }
  elsif (  ( $DNAorPEP eq 'p') || ( $DNAorPEP eq 'pep') || ( $DNAorPEP eq 'prot') || ( $DNAorPEP eq 'protein')  ){
  	if    ( $blastVersion_new_or_not == 0 ){ $DNA_PEP_option=" -p T ";          }
    elsif ( $blastVersion_new_or_not == 1 ){ $DNA_PEP_option=" -dbtype prot ";  }
  }
  else{
  	DieWork::Just_dieWork( $die_MsgHead."\n\$DNAorPEP=$DNAorPEP should be d dna nucl or p pep prot $!".$subCalInfom );
  }
  
  my $parse_seqids_OBTion;
  if    ( $blastVersion_new_or_not == 0 ){ $parse_seqids_OBTion=" -o T ";    }
  elsif ( $blastVersion_new_or_not == 1 ){ $parse_seqids_OBTion=" -parse_seqids ";  }
  
  
  my $outDatabase_option='';	
  if ( DieWork::Check_DfdNoEmptString_or_NOT( $outDB ) ){	
		if    ( $blastVersion_new_or_not == 0 ){ $outDatabase_option=" -t $outDB ";    }
    elsif ( $blastVersion_new_or_not == 1 ){ $outDatabase_option=" -out $outDB ";  }		
	}
  
  my $runingCMD=$cmd_exeFile.$database_path.$DNA_PEP_option.$parse_seqids_OBTion.$outDatabase_option;
  my $HmReadTimeNow=TimeWork::GetHumanRead_NOW_Time();
  DieWork::Print_and_warn( "\n $HmReadTimeNow \$runingCMD=$runingCMD\n\n" );
  
  system ( "$runingCMD" );
  
  
}

#BlastHandle::RunBlast_Old_or_New_VERSION($blastProgram, $queryFile, $database, $outPutFile[, $blastVersion, $evalue, $genetiCode]);
sub RunBlast_Old_or_New_VERSION{
	my ($blastProgram, $queryFile, $database, $outPutFile, $blastVersion, $evalue, $genetiCode, $gCodeForTblastX)=@_;
	my $warnMsgBody="\nIn package  BlastHandle,\tIn sub RunBlast_Old_or_New_VERSION,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";	
	my $subCalInfom=DirFileHandle::print_SubCallerInform;
	
	DieWork::Check_FileDirExist_or_DIE   ( $queryFile,    "\$queryFile",    $die_MsgHead, $subCalInfom  );
  DieWork::Check_FileDirExist_or_DIE   ( $database,     "\$database",     $die_MsgHead, $subCalInfom  );
  DieWork::Check_DfdNoEmptString_or_DIE( $blastProgram, "\$blastProgram", $die_MsgHead, $subCalInfom  );
  DieWork::Check_DfdNoEmptString_or_DIE( $outPutFile,   "\$outPutFile",   $die_MsgHead, $subCalInfom  );
  
	my $blastVersion_new_or_not=0;
  if ( DieWork::Check_DfdNoEmptString_or_NOT( $blastVersion ) ){
  	$blastVersion=lc $blastVersion;   	$blastVersion=~s/\s//g;
  	if   (  ( $blastVersion eq 'blastall' ) || ( $blastVersion eq 'old' )  ){
  		$blastVersion_new_or_not=0;
  	}
  	elsif(  ( $blastVersion eq 'blastplus' ) || ( $blastVersion eq 'new' )  ){
  		$blastVersion_new_or_not=1;
  	}
  	else{
  		DieWork::Just_dieWork( $die_MsgHead."\n\$blastVersion=$blastVersion should be blastall, blastplus, fastacmd, blastdbcmd, old, new  $!".$subCalInfom ); 
  	}  
  }
  else{
  	$blastVersion_new_or_not=1;
  }
  
  $blastProgram=lc $blastProgram;   	$blastProgram=~s/\s//g;
  if (  ( $blastProgram eq 'blastn') || ( $blastProgram eq 'blastp') || ( $blastProgram eq 'tblastn') || ( $blastProgram eq 'blastx') || ( $blastProgram eq 'tblastx')  ){  	  }
  else{ DieWork::Just_dieWork( $die_MsgHead."\n\$blastProgram=$blastProgram should be blastn, blastp, tblastn, blastx, or tblastx  $!".$subCalInfom );  }
  
  my $cmd_exeFile; my $query_option; my $database_option; my $outPutFile_obtion; my $Seg_option='';
  if    ( $blastVersion_new_or_not == 0 ){ $cmd_exeFile="blastall -p $blastProgram"; $query_option=" -i $queryFile ";     $database_option=" -d $database "; $outPutFile_obtion=" -o $outPutFile ";   $Seg_option=" -F F ";    }
  elsif ( $blastVersion_new_or_not == 1 ){ $cmd_exeFile=$blastProgram;               $query_option=" -query $queryFile "; $database_option=" -db $database ";$outPutFile_obtion=" -out $outPutFile "; $Seg_option=" -seg  no ";}
  
  my $eValue_option='';
  if ( DieWork::Check_DfdNoEmptString_or_NOT( $evalue ) ){
    if    ( $blastVersion_new_or_not == 0 ){ $eValue_option=" -e $evalue ";      }
    elsif ( $blastVersion_new_or_not == 1 ){ $eValue_option=" -evalue $evalue "; }
  }
  
  my $Gcd_option='';
  if ( DieWork::Check_DfdNoEmptString_or_NOT( $genetiCode ) ){
  	DieWork::Check_DfdNoEmptNUMBER_or_DIE( $genetiCode, "\$genetiCode", $die_MsgHead, $subCalInfom  );
    if    ( $blastVersion_new_or_not == 0 ){ $Gcd_option=" -D $genetiCode ";          }
    elsif ( $blastVersion_new_or_not == 1 ){ 
      if     ( $blastProgram eq 'tblastn' ){ $Gcd_option=" -db_gencode    $genetiCode "; }  #query_gencode
      elsif  ( $blastProgram eq 'blastx'  ){ $Gcd_option=" -query_gencode $genetiCode "; } 
      elsif  ( $blastProgram eq 'tblastx' ){ DieWork::Check_DfdNoEmptString_or_DIE( $gCodeForTblastX, "\$gCodeForTblastX", $die_MsgHead, $subCalInfom  ); $Gcd_option=" -query_gencode $genetiCode -db_gencode $gCodeForTblastX "; } 	 
      else{
      	DieWork::Just_dieWork( $die_MsgHead."\n\$blastProgram=$blastProgram: \$Gcd_option=$Gcd_option is not OK $!".$subCalInfom );
      }
    }
  }
  
  #my $Seg_option='';
  #if ( DieWork::Check_DfdNoEmptString_or_NOT( $Seg_option ) ){
  #	if ( $Seg_option )
  #  if    ( $blastVersion_new_or_not == 0 ){ $Seg_option=" -F F ";      }
  #  elsif ( $blastVersion_new_or_not == 1 ){ $Seg_option=" -seg  no "; }
  #}
  
  my $runingCMD=$cmd_exeFile.$query_option.$database_option.$outPutFile_obtion.$eValue_option.$Gcd_option.$Seg_option;
  my $HmReadTimeNow=TimeWork::GetHumanRead_NOW_Time();
  DieWork::Print_and_warn( "\n $HmReadTimeNow \$runingCMD=$runingCMD\n\n" );  
  system ( "$runingCMD" );
  
}

sub RunTblastN{  #   BlastHandle::RunTblastN ($queryFile, $database, $outPutFile, $evalue, $genetiCode, $blastVersion);
	my ($queryFile, $database, $outPutFile, $evalue, $genetiCode, $blastVersion)=@_;
	
	my $warnMsgBody="\nIn package  BlastHandle,\tIn sub RunTblastN,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	
	my $blastallPath="blastall";  my $blastPlusTblastNPath="tblastn";
	
  
  if ( defined ($genetiCode) ){
  	if ($genetiCode=~m/^\d+$/){  		  	}  	else {die "\n\n\nDIE!!!!!!!!!\n$warnMsgHead\n\$genetiCode=$genetiCode\nIT should be a number!!!\n\n\n\n";}  	
  }else {$genetiCode=1;}
   
  my $blastToComp_command     ="$blastallPath -i $queryFile -d $database -p tblastn -D $genetiCode -e $evalue -o $outPutFile -a 10 -F F";
	
	my $blastToComp_Plus_command="$blastPlusTblastNPath -query $queryFile -db $database -db_gencode $genetiCode -evalue $evalue   -out $outPutFile -seg  no";
	
	#$blastToComp_command=$blastToComp_Plus_command;
	
	
	my $realCMD=$blastToComp_Plus_command;
  if ( DieWork::Check_DfdNoEmptString_or_NOT( $blastVersion ) ){
  	$blastVersion=lc $blastVersion;
  	$blastVersion=s/\s//g;
  	
  	if (  ( $blastVersion eq 'blastall' ) || ( $blastVersion eq 'old' )  ){
  	  $realCMD=$blastToComp_command;	 
  	}
  }
  
  
  
  print "$warnMsgHead\n$realCMD\n";
  warn  "$warnMsgHead\n$realCMD\n";
  system ( "$realCMD" );
	
	
	
}

sub Run_blast2seq_with_fastaString{  #   BlastHandle::Run_blast2seq_with_fastaString
	#sub Map_PepSeq_BackTo_DNAseq{ #Input 1, protein sequence string, 2, DNA sequence string, 3, a directory path to hold middle step information. TARGET: map the protein sequence back to dna sequence, find the coding region and in frame tga positions 
  my ($inString_1, $inString_2, $workingDir)=@_;
  
 
  
  my $warnMsgBody="\nIn package  BlastHandle,\tIn sub Run_blast2seq_with_fastaString,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
   
  my $File_Str_1_key_word='0_File_Str_1';
  my $File_Str_2_key_word='0_File_Str_2';

  my $inSt_len_1_key_word='1_inSt_len_1';                                           
  my $inSt_len_2_key_word='1_inSt_len_2';                                         
  
  my $BlRstS1xS2_key_word='2_BlRstS1xS2';
  my $BlRstS2xS1_key_word='3_BlRstS2xS1';
  
  #my $bstBlS1xS2_key_word='4_bstBlS1xS2';
  #my $bstBlS2xS1_key_word='5_bstBlS2xS1';
  
  my $Png_S1vsS2_key_word='4_Png_S1vsS2';
  my $Png_S2vsS1_key_word='5_Png_S2vsS1';
  
  my $blIfm_S1S2_key_word='6_blIfm_S1S2';
  my $blIfm_S2S1_key_word='7_blIfm_S2S1';
  
  my $FblHS_S1S2_key_word='8_FblHS_S1S2';
  my $FblHS_S2S1_key_word='9_FblHS_S2S1';

  
  #my $blastDIRWD="bl2seq_work";

  
  my $outInformationHash;  
    
  #s$inString_1=~s/\n//g; $inString_2=~s/\s//g;
  my $inSt_len_1=FastaFileHandle::GetSeqLength($inString_1);           $outInformationHash->{$inSt_len_1_key_word}=$inSt_len_1; 
  my $inSt_len_2=FastaFileHandle::GetSeqLength($inString_2);           $outInformationHash->{$inSt_len_2_key_word}=$inSt_len_2; 
  
  
  
  my $tempTimeDirNeedDel=0; if (  ( defined($workingDir) ) && ($workingDir=~/\S+/)  ) {} else {		$workingDir=TimeWork::GetTimeDirOrFileName();	$tempTimeDirNeedDel=1; }  #If the $analysisDIR is not set, then build a temporary directory named by time #print "\$workingDir=$workingDir\n";
  system ( "mkdir -p $workingDir ") ;
  
  my $blastWorkDIr=$workingDir;    	system ( "mkdir -p $blastWorkDIr ") ;
  
  my $File_Str_1="$blastWorkDIr/$File_Str_1_key_word.txt";             $outInformationHash->{$File_Str_1_key_word}=$File_Str_1;            
  my $File_Str_2="$blastWorkDIr/$File_Str_2_key_word.txt";	           $outInformationHash->{$File_Str_2_key_word}=$File_Str_2;
  
  #my $FastaStr_1=">$File_Str_1_key_word\n$inString_1\n\n";   
  #my $FastaStr_2=">$File_Str_2_key_word\n$inString_2\n\n";     
  open (INSTR_1, ">$File_Str_1") or die "$die_MsgHead\ncannot create \$File_Str_1=$File_Str_1 :$!\n\n\n"; print INSTR_1 $inString_1;  close (INSTR_1);
  open (INSTR_2, ">$File_Str_2") or die "$die_MsgHead\ncannot create \$File_Str_2=$File_Str_2 :$!\n\n\n"; print INSTR_2 $inString_2;  close (INSTR_2);
  
  my $BlRstS1xS2="$blastWorkDIr/$BlRstS1xS2_key_word.txt";             $outInformationHash->{$BlRstS1xS2_key_word}=$BlRstS1xS2;                    	
  my $BlRstS2xS1="$blastWorkDIr/$BlRstS2xS1_key_word.txt";             $outInformationHash->{$BlRstS2xS1_key_word}=$BlRstS2xS1;                    	
  
  #my $bstBlS1xS2="$blastWorkDIr/$bstBlS1xS2_key_word.html";            $outInformationHash->{$bstBlS1xS2_key_word}=$bstBlS1xS2;
  #my $bstBlS2xS1="$blastWorkDIr/$bstBlS2xS1_key_word.html";            $outInformationHash->{$bstBlS2xS1_key_word}=$bstBlS2xS1;
  
  #my $makeblastdbPath="formatdb";                  #blastplus vesrion:  my $makeblastdbPath="~/EightT/bin/makeblastdb";
  #my $tblastnExecPath="blastall";                  #blastplus vesrion:  my $tblastnExecPath="~/EightT/bin/tblastn";
  #$bl2seq_app_Path="bl2seq";
  my $appPathNme="bl2seq";
  my $subAppName="blastn";
  my $eValueCtof=0.000001;
  
  #my $mkDatabaseCMD_for_str_1="$makeblastdbPath -i $File_Str_1 -p F -o T";                                                                       #blastplus vesrion: #my $mkDatabaseCMD_forDNA="$makeblastdbPath -in $DNAseqFastaFile -dbtype nucl -parse_seqids";
  #my $mkDatabaseCMD_for_str_2="$makeblastdbPath -i $File_Str_2 -p F -o T";
  #my $mkDatabaseCMD_for_str_2="$makeblastdbPath -i $File_Str_2 -p T -o T";
  
  #my $db_gencode=1; 	my $U_pos_array=&Find_U_pos_inAminoAcid($inPEPstring);	if ( ref($U_pos_array) eq 'ARRAY' ){ $db_gencode=10; }
  #my $blastCMD=$tblastnExecPath." -p tblastn  -i ".$PEPseqFastaFile." -d "."$DNAseqFastaFile -o $tBLASTnRslt -e 0.0000001 -D $db_gencode -F F";    #blastplus vesrion: #my $blastCMD="$tblastnExecPath –query $PEPseqFastaFile –db $DNAseqFastaFile –out $tBLASTnRslt -evalue 0.0000001 -db_gencode $db_gencode -seg no" ;
  
  #my $blast__CMD=$appPathNme." -F F -p $subAppName  -i ".$File_Str_1." -j "."$File_Str_2 -o $tBLASTnRslt -e 0.0000001 -D $db_gencode -F F";    #blastplus vesrion: #my $blastCMD="$tblastnExecPath –query $PEPseqFastaFile –db $DNAseqFastaFile –out $tBLASTnRslt -evalue 0.0000001 -db_gencode $db_gencode -seg no" ;
  
  my $bltCMDS1S2="$appPathNme -F F -p $subAppName  -i $File_Str_1 -j $File_Str_2 -e $eValueCtof -o $BlRstS1xS2";
  my $bltCMDS2S1="$appPathNme -F F -p $subAppName  -i $File_Str_2 -j $File_Str_1 -e $eValueCtof -o $BlRstS2xS1";
  
  
  
  #warn "\n\$mkDatabaseCMD_forDNA=$mkDatabaseCMD_forDNA\n"; print "\n\$mkDatabaseCMD_forDNA=$mkDatabaseCMD_forDNA\n";        system ("$mkDatabaseCMD_forDNA");  
  #warn "\n\$mkDatabaseCMD_forPEP=$mkDatabaseCMD_forPEP\n"; print "\n\$mkDatabaseCMD_forPEP=$mkDatabaseCMD_forPEP\n";        system ("$mkDatabaseCMD_forPEP");  
  warn "\n\$bltCMDS1S2=$bltCMDS1S2\n";            print "\n\$bltCMDS1S2=$bltCMDS1S2\n";       system ("$bltCMDS1S2");
  warn "\n\$bltCMDS2S1=$bltCMDS2S1\n";            print "\n\$bltCMDS2S1=$bltCMDS2S1\n";       system ("$bltCMDS2S1");
  
  my $Png_S1vsS2="$blastWorkDIr/$Png_S1vsS2_key_word.png";             $outInformationHash->{$Png_S1vsS2_key_word}=$Png_S1vsS2;                    	
  my $Png_S2vsS1="$blastWorkDIr/$Png_S2vsS1_key_word.png";             $outInformationHash->{$Png_S2vsS1_key_word}=$Png_S2vsS1;
  
  
  my $blIfm_S1S2=&DrawBlast2seq_20180920($BlRstS1xS2, $Png_S1vsS2 );   $outInformationHash->{$blIfm_S1S2_key_word}=$blIfm_S1S2;   
  my $blIfm_S2S1=&DrawBlast2seq_20180920($BlRstS2xS1, $Png_S2vsS1 );   $outInformationHash->{$blIfm_S2S1_key_word}=$blIfm_S2S1;
  
  my $FblHS_S1S2="$blastWorkDIr/$File_Str_1_key_word.hash";  if ( ref ($blIfm_S1S2)  eq 'HASH' ) { DirFileHandle::PrintDumper($FblHS_S1S2,$blIfm_S1S2); $outInformationHash->{$File_Str_1_key_word}=$FblHS_S1S2; }
  my $FblHS_S2S1="$blastWorkDIr/$File_Str_2_key_word.hash";	 if ( ref ($blIfm_S2S1)  eq 'HASH' ) { DirFileHandle::PrintDumper($FblHS_S2S1,$blIfm_S2S1); $outInformationHash->{$File_Str_2_key_word}=$FblHS_S2S1; }
  	
  
  return $outInformationHash; 

}



sub ChangeTo100PercentNB{
	my ($inNumber)=@_;
	
	my $warnMsgBody="\nIn package  BlastHandle,\tIn sub ChangeTo100PercentNB,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	my $outNumber=100*$inNumber;  $outNumber=sprintf "%.2f",$outNumber; $outNumber="$outNumber%";
	return $outNumber;
}


sub Change100PercentNB_to_pointNB{
	my ($inNumber)=@_;
	
	my $warnMsgBody="\nIn package  BlastHandle,\tIn sub Change100PercentNB_to_pointNB,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $caller_inform=DirFileHandle::print_SubCallerInform;

  my $outNumber;
	if ($inNumber=~m/^\s*(\d+(\.\d+)?)%\s*$/){
		$outNumber=$1/100;
	}
	else {
		my $dieMsg=$die_MsgHead.$caller_inform."\$inNumber=$inNumber is not a correct 100% precent number :$! \n\n\n\n";
		die $dieMsg;
	}
	
	return $outNumber;
}
    
sub DrawBlast{      #pmAble#   #绘制Blast 的图像
	my ($file, $outPNGfile)=@_;   warn "In sub &DrawBlast (\$file=$file, \$outPNGfile=$outPNGfile)\n\n";
	open (BLASTPNGOUT,">$outPNGfile") or die "cannot create \$outPNGfile=$outPNGfile : $!\n\n";

  my $searchio = Bio::SearchIO->new(-file   => $file,
                                    -format => 'blast') or die "parse failed";


  my $result = $searchio->next_result() or die "\n\nno result\n\n";

  my $panel = Bio::Graphics::Panel->new(-length    => $result->query_length,
                                        -width     => 800,
                                        -pad_left  => 10,
                                        -pad_right => 10,
                                       );

  my $full_length = Bio::SeqFeature::Generic->new(-start =>   1, 
                                                  -end   =>   $result->query_length,
                                                  -seq_id=>   $result->query_name);
  $panel->add_track($full_length,
                    -glyph   => 'arrow',
                    -tick    => 2,
                    -fgcolor => 'black',
                    -double  => 1,
                    -label   => 1,
                   );

  my $track = $panel->add_track(-glyph       => 'segments',  #'graded_segments',
                                -label       => 1,
                                -connector   => 'dashed',
                                -bgcolor     => 'red',
                                -fgcolor     => 'red',
                                -font2color  => 'black',
                                -sort_order  => 'high_score',
                                -description => sub {
                                                       my $feature = shift;
                                                       return unless $feature->has_tag('description');
                                                       my ($description) = $feature->each_tag_value('description');
                                                       my $score = $feature->score;
                                                       #my $eValue= $feature->
                                                       "score=$score $description";
                                                    }
                                );
  
  while( my $hit = $result->next_hit ) {
    #next unless $hit->significance < 1E-20;
    my $feature = Bio::SeqFeature::Generic->new(-score   => $hit->raw_score,
                                                -seq_id  => $hit->name,
                                                -tag     => {
                                                             description => $hit->description
                                                            },
                                               );
    while( my $hsp = $hit->next_hsp ) {
      $feature->add_sub_SeqFeature($hsp,'EXPAND');
    }
  
    $track->add_feature($feature);
  }
  
  print BLASTPNGOUT $panel->png;

  close (BLASTPNGOUT);
	
	
	
}

sub BuildTGAPosArrayFrom_TGAchangeHASH{
	my ($in_TGAchangeHASH, $wholeCtgLength)=@_;
	my $warnMsgHead="\n\n\n		In package BlastHandle,		\n In sub BuildTGAPosArrayFrom_TGAchangeHASH,\n"; 
	my $outArray; my $arIdx=0;
	if ( ref($in_TGAchangeHASH) eq 'HASH' ){
	  foreach my $tga_A_pos (    sort (   keys (  %{ $in_TGAchangeHASH }  )   )    ) {
		  if (   (  ref( $in_TGAchangeHASH->{$tga_A_pos} ) eq 'HASH'  ) && (  defined ( $in_TGAchangeHASH->{$tga_A_pos}->{'2_fromCha'} )  ) && (  ( $in_TGAchangeHASH->{$tga_A_pos}->{'2_fromCha'} eq 'A' ) || ( $in_TGAchangeHASH->{$tga_A_pos}->{'2_fromCha'} eq 'a' )  )   ){
		  	my $tgaPosZoF=$in_TGAchangeHASH->{$tga_A_pos}->{'0_directi'};
		  	my $ZorF_factor=1; if ($tgaPosZoF eq '-'){ $ZorF_factor=-1; }
		  	$outArray->[$arIdx]->{'0_ExonHead'}=$in_TGAchangeHASH->{$tga_A_pos}->{'1_postion'}-$ZorF_factor*2;
		  	$outArray->[$arIdx]->{'1_ExonTail'}=$in_TGAchangeHASH->{$tga_A_pos}->{'1_postion'};
		  	$outArray->[$arIdx]->{'2_ExonZorF'}=$tgaPosZoF;
		  	$outArray->[$arIdx]->{'4_SegmType'}='inFrameTGA';

		  	if (  ( defined ($wholeCtgLength) ) && ($wholeCtgLength=~m/^\d+$/) ){
		  		my $TGAsttPos=$outArray->[$arIdx]->{'0_ExonHead'};
		  		my $TGA_frame;
		  		if ($tgaPosZoF eq '+'){ $TGA_frame=$TGA_frame%3;                     if ($TGA_frame==0){ $TGA_frame=3; } }
		  		if ($tgaPosZoF eq '-'){ $TGA_frame=($wholeCtgLength-$TGA_frame+1)%3; if ($TGA_frame==0){ $TGA_frame=3; } }
		  		$TGA_frame=$tgaPosZoF.$TGA_frame; 
		  		$outArray->[$arIdx]->{'3_ExonFram'}=$TGA_frame;
		  	}
		  	
		  	$arIdx++;
		  }
		  else {
		  	my $in_TGAchangeHASH_string=DirFileHandle::PrintAndWarnDumper ($in_TGAchangeHASH);		  	die "\n\n\nDIE!!!!!$warnMsgHead\nThe \$in_TGAchangeHASH=$in_TGAchangeHASH\n\n\n";
		  }
	  }	
	}
	return $outArray;
}

sub FindAllUForQueryAndConsAAinHit_20180626{   # my $blastParsedHash= BlastHandle::FindAllUForQueryAndConsAAinHit_20180626 ($blastfile, $outBestBlastFile, $numHits, $report_type ) ;
	my ($blastfile, $outBestBlastFile, $numHits, $report_type )=@_;
	
	
	my $warnMsgBody="\nIn package  BlastHandle,\tIn sub FindAllUForQueryAndConsAAinHit_20180626,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	my $blastParsedHash=BioPerlBlastPraser20170325($blastfile, $outBestBlastFile, $numHits, $report_type);
  $blastParsedHash=FindAllUForQueryAndConsAAinHit_fromHASH_20190322 ($blastParsedHash, $outBestBlastFile, $numHits, $report_type );
	
	return $blastParsedHash;
}



sub FindAllUForQueryAndConsAAinHit_fromHASH_20190322{  # my $blastParsedHash= BlastHandle::FindAllUForQueryAndConsAAinHit_fromHASH_20190322 ($blastParsedHash, $outBestBlastFile, $numHits, $report_type ) ;
	my ($blastParsedHash, $outBestBlastFile, $numHits, $report_type )=@_;
  
  my $warnMsgBody="\nIn package  BlastHandle,\tIn sub FindAllUForQueryAndConsAAinHit_fromHASH_20190322,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	
	my $caller_inform=DirFileHandle::print_SubCallerInform;
#} 

#sub FindAllUForQueryAndConsAAinHit_20180626{  #从blast结果中，通常是tblastn，找出每个hsp中的query的U的位置，以及对应hit的位置上是什么氨基酸，并记录其DNA，并将hit中的U对应位置的TGA翻译为U，并记录所有*号，及其密码子
#	my ($blastfile, $outBestBlastFile, $numHits, $report_type )=@_;
#	my $blastParsedHash=BioPerlBlastPraser20170325($blastfile, $outBestBlastFile, $numHits, $report_type);


	if (   (ref ($blastParsedHash) eq 'HASH') && (  ref ( $blastParsedHash->{'_ResultArray'} ) eq 'ARRAY'  )   ){
		for (my $rstIdx=0; $rstIdx<@{ $blastParsedHash->{'_ResultArray'} }; $rstIdx++){
		  
		  if (   (  ref ( $blastParsedHash->{'_ResultArray'}->[$rstIdx] ) eq 'HASH'  )  &&  (  ref ( $blastParsedHash->{'_ResultArray'}->[$rstIdx]->{'_hitArray'} ) eq 'ARRAY'  )   ){
		  	for (my $hitIdx=0; $hitIdx<@{ $blastParsedHash->{'_ResultArray'}->[$rstIdx]->{'_hitArray'} }; $hitIdx++){
		  	
		  	  if (   (  ref ($blastParsedHash->{'_ResultArray'}->[$rstIdx]->{'_hitArray'}->[$hitIdx] ) eq 'HASH'  )  &&  (  ref ( $blastParsedHash->{'_ResultArray'}->[$rstIdx]->{'_hitArray'}->[$hitIdx]->{'_hspArray'} ) eq 'ARRAY'  )   ){
		  	    for (my $hspIdx=0; $hspIdx<@{ $blastParsedHash->{'_ResultArray'}->[$rstIdx]->{'_hitArray'}->[$hitIdx]->{'_hspArray'} }; $hspIdx++){
		  	    
		  	      #$blastParsedHash->{'_ResultArray'}->[$rstIdx]->{'_hitArray'}->[$hitIdx]->{'_hspArray'}->[$hspIdx]->{'_Hsp_3WhSting_4_1QrNB2Ca'}  	   
		          my $QryWholeString=$blastParsedHash->{'_ResultArray'}->[$rstIdx]->{'_hitArray'}->[$hitIdx]->{'_hspArray'}->[$hspIdx]->{'_Hsp_3WhSting_1Qr'};     #query 包括gap的比对序列string.   the string of query including gaps
		          my $HitWholeString=$blastParsedHash->{'_ResultArray'}->[$rstIdx]->{'_hitArray'}->[$hitIdx]->{'_hspArray'}->[$hspIdx]->{'_Hsp_3WhSting_3Ht'};     #hit   包括gap的比对序列string.   the string of hit   including gaps
		           
			  	    my $QryZoF        =$blastParsedHash->{'_ResultArray'}->[$rstIdx]->{'_hitArray'}->[$hitIdx]->{'_hspArray'}->[$hspIdx]->{'_Hsp_1stfm_3ZoF_1qur'};  #query Hsp的正负. Minus or Plus of the query. query + or -.     
			      	my $HitZoF        =$blastParsedHash->{'_ResultArray'}->[$rstIdx]->{'_hitArray'}->[$hitIdx]->{'_hspArray'}->[$hspIdx]->{'_Hsp_1stfm_4ZoF_1hit'};  #hit   Hsp的正负. Minus or Plus of the hit  . hit   + or -.
			  	   
			  	    my $QryStt        =$blastParsedHash->{'_ResultArray'}->[$rstIdx]->{'_hitArray'}->[$hitIdx]->{'_hspArray'}->[$hspIdx]->{'_Hsp_2StEd_1Qr_1St'};   #query Hsp的起始位置 . Start postion of the query.     
			       	my $QryEnd        =$blastParsedHash->{'_ResultArray'}->[$rstIdx]->{'_hitArray'}->[$hitIdx]->{'_hspArray'}->[$hspIdx]->{'_Hsp_2StEd_1Qr_2Ed'};   #query Hsp的终止位置 . End   postion of the query. 
			      	
			      	my $HitStt        =$blastParsedHash->{'_ResultArray'}->[$rstIdx]->{'_hitArray'}->[$hitIdx]->{'_hspArray'}->[$hspIdx]->{'_Hsp_2StEd_2Ht_1St'};   #hit   Hsp的起始位置 . Start postion of the hit.     
			      	my $HitEnd        =$blastParsedHash->{'_ResultArray'}->[$rstIdx]->{'_hitArray'}->[$hitIdx]->{'_hspArray'}->[$hspIdx]->{'_Hsp_2StEd_2Ht_2Ed'};   #hit   Hsp的终止位置 . End   postion of the hit. 
			  	
			  	    my $dtDBnm        =$blastParsedHash->{'_ResultArray'}->[$rstIdx]->{'_database_name'};
			  	    my $hitNam        =$blastParsedHash->{'_ResultArray'}->[$rstIdx]->{'_hitArray'}->[$hitIdx]->{'_accessionNB'};
			  	    
			  	    my $QryNam        =$blastParsedHash->{'_ResultArray'}->[$rstIdx]->{'_query_name'}; 
			  	    
			  	    my $idxFnm=$dtDBnm; 			  	            $idxFnm=~s/\.txt$/\.idx/;
			  	    if   (   (  defined ( $idxFnm )  ) && ( $idxFnm=~m/\S+/ ) && (  ( -e $idxFnm )  )   ){}
			  	    else{		DieWork::Just_dieWork( $die_MsgHead."\n \$idxFnm=$idxFnm should be a not empty idx file  !!  $!\n\n\n".$caller_inform ); 	}
	
			  	    
	
			  	    
			  	    print "\n 20190406-0-0-0 ", TimeWork::GetNowTimePid_microSecond ()," \$QryNam=$QryNam \n 20190406-0-0-1 ", TimeWork::GetNowTimePid_microSecond ()," \$QryWholeString=$QryWholeString \n";
			  	    print "\n 20190406-0-0-2 ", TimeWork::GetNowTimePid_microSecond ()," \$hitNam=$hitNam \n 20190406-0-0-3 ", TimeWork::GetNowTimePid_microSecond ()," \$HitWholeString=$HitWholeString \n";
			  	    
			  	    my $target_found_mark=0;
			  	    
			  	    
			  	    my ($qry_AlnPos_to_Char_HASH, $qry_AlnPos_to_OrgPos_HASH, $qry_OrgPos_to_AlnPos_HASH);
			  	    my ($hit_AlnPos_to_Char_HASH, $hit_AlnPos_to_OrgPos_HASH, $hit_OrgPos_to_AlnPos_HASH); my $tblsatnOrgHSPDNAseq;
			  	    
			  	    my ($qry_seq_including_C, $qry_seq_including_U, $hit_seq_including_star, $hit_seq_including_w, $hit_seq_including_C,$hit_seq_including_X)=
			  	       (0,                    0,                    0,                       0,                    0,                   0,                  );
			  	       
			  	    if (0){ # 下面这种 判断 query和hit中的UC*WX出现的 模式，适合于 寻找较小数据量中，所有 query和target中的UC U* UW UX UU CC C* CW CX CU。
			  	    	if ($QryWholeString=~m/C/i  ){$qry_seq_including_C=1;   }   if ($QryWholeString=~m/U/i ){$qry_seq_including_U=1;}
			  	      if ($HitWholeString=~m/\*/i ){$hit_seq_including_star=1;}   if ($HitWholeString=~m/W/i ){$hit_seq_including_w=1;}     if ($HitWholeString=~m/C/i ){$hit_seq_including_C=1;}   if ($HitWholeString=~m/X/i ){$hit_seq_including_X=1;}
		            if ( ($qry_seq_including_C) ||  ($qry_seq_including_U) ||  ($hit_seq_including_star) ||  ($hit_seq_including_w) ||  ($hit_seq_including_C) ||  ($hit_seq_including_X) ){
		              $target_found_mark=1;
		            }
			  	    }
			  	    
			  	    if (1){ # 下面这种 判断 query和hit中的UC*WX出现的 模式，适合于 寻找   仅仅同时在hit中，w,x,* 出现， ###############在query中出现了U和C时，才开始后续分析了。
			  	    	#if ($QryWholeString=~m/C/i  ){$qry_seq_including_C=1;   }   if ($QryWholeString=~m/U/i ){$qry_seq_including_U=1;}
			  	      if ($HitWholeString=~m/\*/i ){$hit_seq_including_star=1;}   if ($HitWholeString=~m/W/i ){$hit_seq_including_w=1;} if ($HitWholeString=~m/X/i ){$hit_seq_including_X=1;}
		            #if (   (  ($qry_seq_including_C) ||  ($qry_seq_including_U)  ) && (  ( ($hit_seq_including_star) ||  ($hit_seq_including_w)  ||  ($hit_seq_including_X)  )   ){
		            if (  ($hit_seq_including_star) ||  ($hit_seq_including_w) || ($hit_seq_including_X) ){
		            	  $target_found_mark=1;
		            }
			  	    }
			  	    
		          
		          print "\n 20190406-0-0-4 ", TimeWork::GetNowTimePid_microSecond ()," \$qry_seq_including_C=$qry_seq_including_C \$qry_seq_including_U=$qry_seq_including_U\n";
		          print "\n 20190406-0-0-5 ", TimeWork::GetNowTimePid_microSecond ()," \$hit_seq_including_star=$hit_seq_including_star \$hit_seq_including_w=$hit_seq_including_w \$hit_seq_including_C=$hit_seq_including_C \$hit_seq_including_X=$hit_seq_including_X \n";
		          if ( $target_found_mark == 1 ){
		          	print "\n 20190406-0-0-6 ", TimeWork::GetNowTimePid_microSecond ()," \$qry_seq_including_C=$qry_seq_including_C \$qry_seq_including_U=$qry_seq_including_U    \$hit_seq_including_star=$hit_seq_including_star \$hit_seq_including_w=$hit_seq_including_w \$hit_seq_including_C=$hit_seq_including_C \$hit_seq_including_X=$hit_seq_including_X \n";
		          	my $Blastx_or_not=0;
		          	if (   (  defined ( $report_type )  ) && ( $report_type=~m/\S+/ )   ){
		          		$report_type=lc $report_type;
		          		if  ( $report_type eq 'blastx' ){
		          			$Blastx_or_not=1;
		          		}
		          	}
		          	if ( $Blastx_or_not==0){
		          	  $qry_AlnPos_to_Char_HASH=&build_NBtoWord_HASH_for_String($QryWholeString);    ( $qry_AlnPos_to_OrgPos_HASH, $qry_OrgPos_to_AlnPos_HASH                      ) =@{  &build_NBtoWord_HASH_for_String_at_OrgSEQUENCE($qry_AlnPos_to_Char_HASH, $QryZoF,  $QryStt, $QryEnd, 1                   )  }; #DirFileHandle::PrintAndWarnDumper ($qry_AlnPos_to_Char_HASH, "20190406-0-1-0");
			  	        $hit_AlnPos_to_Char_HASH=&build_NBtoWord_HASH_for_String($HitWholeString);    ( $hit_AlnPos_to_OrgPos_HASH, $hit_OrgPos_to_AlnPos_HASH, $tblsatnOrgHSPDNAseq) =@{  &build_NBtoWord_HASH_for_String_at_OrgSEQUENCE($hit_AlnPos_to_Char_HASH, $HitZoF,  $HitStt, $HitEnd, 3, $idxFnm, $hitNam )  };	#DirFileHandle::PrintAndWarnDumper ($hit_AlnPos_to_Char_HASH, "20190406-0-1-1");				  	      	
		          	}
		          	else {
		          		$qry_AlnPos_to_Char_HASH=&build_NBtoWord_HASH_for_String($QryWholeString);    ( $qry_AlnPos_to_OrgPos_HASH, $qry_OrgPos_to_AlnPos_HASH                      ) =@{  &build_NBtoWord_HASH_for_String_at_OrgSEQUENCE($qry_AlnPos_to_Char_HASH, $QryZoF,  $QryStt, $QryEnd, 3, $idxFnm, $QryNam )  };
			  	        $hit_AlnPos_to_Char_HASH=&build_NBtoWord_HASH_for_String($HitWholeString);    ( $hit_AlnPos_to_OrgPos_HASH, $hit_OrgPos_to_AlnPos_HASH, $tblsatnOrgHSPDNAseq) =@{  &build_NBtoWord_HASH_for_String_at_OrgSEQUENCE($hit_AlnPos_to_Char_HASH, $HitZoF,  $HitStt, $HitEnd, 1                   )  };	  	      
		          	}
		          	
		          	my $ChangTGAintoSECforHitString=$HitWholeString;
			  	      my $otherSTOPchangeString      =$HitWholeString;
			  	      
			  	      
			  	      
			  	      my $QryCysPosArray; if ($qry_seq_including_C){$QryCysPosArray=&findingIdexOfAword ( 'C', $QryWholeString );}  #找到query中所有的C的位置  find the postions of all of the C in query
			  	      my $QrySecPosArray; if ($qry_seq_including_U){$QrySecPosArray=&findingIdexOfAword ( 'U', $QryWholeString );}  #找到query中所有的U的位置  find the postions of all of the U in query
			  	      my $QryCysPosArrayString='';$QryCysPosArrayString=ArrayHashChange::ChangArrayIntoSimpleString ($QryCysPosArray) if (  ( ref ($QryCysPosArray) ) eq 'ARRAY'  );
			  	      my $QrySecPosArrayString='';$QrySecPosArrayString=ArrayHashChange::ChangArrayIntoSimpleString ($QrySecPosArray) if (  ( ref ($QrySecPosArray) ) eq 'ARRAY'  );
			  	      print "\n 20190406-0-2-0 ", TimeWork::GetNowTimePid_microSecond ()," \$QrySecPosArray=$QrySecPosArray=$QryCysPosArrayString \$QryCysPosArray=$QryCysPosArray=$QrySecPosArrayString\n";
			  	      if (    (  ( ref ($QrySecPosArray) ) eq 'ARRAY'  ) || (  ( ref ($QryCysPosArray) ) eq 'ARRAY'  )    ) {   #这个ref成立，则表明在query中找到了U.    Find U in query!!
			  	      	print "\n 20190406-0-2-1 ", TimeWork::GetNowTimePid_microSecond ()," \$QrySecPosArray=$QrySecPosArray=$QryCysPosArrayString \$QryCysPosArray=$QryCysPosArray=$QrySecPosArrayString\n";
			  	      	my @all_U_C_posArray; 
			  	      	push @all_U_C_posArray,  @{ $QrySecPosArray } if (  ( ref ($QrySecPosArray) ) eq 'ARRAY'  );
			  	      	push @all_U_C_posArray,  @{ $QryCysPosArray } if (  ( ref ($QryCysPosArray) ) eq 'ARRAY'  );
			  	      		
			  	        my $findTGA=0;
			  	      	foreach my $eachQrySecPos (  @all_U_C_posArray  ){
			  	      		if (  defined ( $hit_AlnPos_to_OrgPos_HASH->{$eachQrySecPos}->{'3DNAchar'} )  ){
			  	      			my $tp3DNAcodon=uc ( $hit_AlnPos_to_OrgPos_HASH->{$eachQrySecPos}->{'3DNAchar'} );  print "\n 20190406-0-2-2 ", TimeWork::GetNowTimePid_microSecond ()," \$eachQrySecPos=$eachQrySecPos \$tp3DNAcodon=$tp3DNAcodon \n";
			  	      			if ($tp3DNAcodon eq 'TGA'){
			  	      				$findTGA++;
			  	      			  substr($ChangTGAintoSECforHitString, ($eachQrySecPos-1),1)='U';
			  	      			  
			  	      			  $blastParsedHash->{'_ResultArray'}->[$rstIdx]->{'_hitArray'}->[$hitIdx]->{'_hspArray'}->[$hspIdx]->{'_Hsp_4Hit_1Utga_MatchToQryUC'}->{ $qry_AlnPos_to_Char_HASH->{$eachQrySecPos} }->{$eachQrySecPos}=$hit_AlnPos_to_OrgPos_HASH->{$eachQrySecPos};			  	    			  
			  	      			                                                                                   print "\n 20190406-0-2-3 ", TimeWork::GetNowTimePid_microSecond ()," \$eachQrySecPos=$eachQrySecPos \$tp3DNAcodon=$tp3DNAcodon \n";                                                                                                 
			  	      			  #下面是用来在以后的分析中，如基因结构识别等步骤中，将TGA转变成TGC的
			  	      			  my $hitOrgPos=$hit_AlnPos_to_OrgPos_HASH->{$eachQrySecPos}->{'real_1_3_Pos'}->{'3'};	
			  	      			  $blastParsedHash->{'_ResultArray'}->[$rstIdx]->{'_hitArray'}->[$hitIdx]->{'_TGAchangeHASH'}->{$hitOrgPos}->{'0_directi'}=$HitZoF;	
			  	      			  $blastParsedHash->{'_ResultArray'}->[$rstIdx]->{'_hitArray'}->[$hitIdx]->{'_TGAchangeHASH'}->{$hitOrgPos}->{'1_postion'}=$hitOrgPos;			  	    			  
			  	      			  $blastParsedHash->{'_ResultArray'}->[$rstIdx]->{'_hitArray'}->[$hitIdx]->{'_TGAchangeHASH'}->{$hitOrgPos}->{'2_fromCha'}='A';
			  	      			  $blastParsedHash->{'_ResultArray'}->[$rstIdx]->{'_hitArray'}->[$hitIdx]->{'_TGAchangeHASH'}->{$hitOrgPos}->{'3_intoCha'}='C';
			  	      			  $blastParsedHash->{'_ResultArray'}->[$rstIdx]->{'_hitArray'}->[$hitIdx]->{'_TGAchangeHASH'}->{$hitOrgPos}->{'4_sourcet'}='tBlastn, \'A\' of tga matched to '.$qry_AlnPos_to_Char_HASH->{$eachQrySecPos};
			  	      			  
			  	      			  
			  	      		  }
			  	      		}
			  	      		
			  	      	}
			  	      	
			  	      	if ($findTGA>0){           print "\n 20190406-0-2-4 ", TimeWork::GetNowTimePid_microSecond ()," \$blastParsedHash->{'_ResultArray'}->[$rstIdx]->{'_hitArray'}->[$hitIdx]->{'_hspArray'}->[$hspIdx]->{'_Hsp_3WhSting_4Uh'} \n";
			  	        	
			  	        	$blastParsedHash->{'_ResultArray'}->[$rstIdx]->{'_hitArray'}->[$hitIdx]->{'_hspArray'}->[$hspIdx]->{'_Hsp_3WhSting_4Uh'}=  $ChangTGAintoSECforHitString;
			  	        	$otherSTOPchangeString=$ChangTGAintoSECforHitString;
			  	        	
			  	        	#$blastParsedHash->{'_ResultArray'}->[$rstIdx]->{'_hitArray'}->[$hitIdx]->{'_hspArray'}->[$hspIdx]->{'_Hsp_3WhSting_4_1QrAlnPosToChar'}      =$qry_AlnPos_to_Char_HASH;
			  	        	#$blastParsedHash->{'_ResultArray'}->[$rstIdx]->{'_hitArray'}->[$hitIdx]->{'_hspArray'}->[$hspIdx]->{'_Hsp_3WhSting_4_2HtAlnPosToChar'}      =$hit_AlnPos_to_Char_HASH;			  	      	
			  	        	#$blastParsedHash->{'_ResultArray'}->[$rstIdx]->{'_hitArray'}->[$hitIdx]->{'_hspArray'}->[$hspIdx]->{'_Hsp_5CoordPos_3QrAlnPosToOrgPos'}    =$qry_AlnPos_to_OrgPos_HASH;
			  	        	#$blastParsedHash->{'_ResultArray'}->[$rstIdx]->{'_hitArray'}->[$hitIdx]->{'_hspArray'}->[$hspIdx]->{'_Hsp_5CoordPos_4QrOrgPosToAlnPos'}    =$qry_OrgPos_to_AlnPos_HASH;
			  	        	#$blastParsedHash->{'_ResultArray'}->[$rstIdx]->{'_hitArray'}->[$hitIdx]->{'_hspArray'}->[$hspIdx]->{'_Hsp_5CoordPos_5HtAlnPosToOrgPos'}    =$hit_AlnPos_to_OrgPos_HASH;
			  	        	#$blastParsedHash->{'_ResultArray'}->[$rstIdx]->{'_hitArray'}->[$hitIdx]->{'_hspArray'}->[$hspIdx]->{'_Hsp_5CoordPos_6HtOrgPosToAlnPos'}    =$hit_OrgPos_to_AlnPos_HASH;
			  	        	$blastParsedHash->{'_ResultArray'}->[$rstIdx]->{'_hitArray'}->[$hitIdx]->{'_hspArray'}->[$hspIdx]->{'_Hsp_5CoordPos_7HitOrgDNASequens'}    =$tblsatnOrgHSPDNAseq
			  	        	
			  	        	
			  	        	
			  	        }
			  	      }
			  	      
			  	      
			  	      my $HitStpPosArray; if ($hit_seq_including_star){$HitStpPosArray=&findingIdexOfAword ( '*', $HitWholeString );}  #寻找Hit  中的所有*的位置  find the postions of all of the * in hit
		            my $HitTrpPosArray; if ($hit_seq_including_w   ){$HitTrpPosArray=&findingIdexOfAword ( 'W', $HitWholeString );}  #寻找Hit  中的所有W的位置  find the postions of all of the W in hit
		            my $HitXxxPosArray; if ($hit_seq_including_X   ){$HitXxxPosArray=&findingIdexOfAword ( 'X', $HitWholeString );}  #寻找Hit  中的所有X的位置  find the postions of all of the X in hit		            
		            #下面这行是为了收集 hit中的CYs，这里先去掉
		            #my $HitCysPosArray; if ($hit_seq_including_C   ){$HitCysPosArray=&findingIdexOfAword ( 'C', $HitWholeString );}  #寻找Hit  中的所有C的位置  find the postions of all of the C in hit
			  	  	  
			  	  	  my $HitStpPosArrayString='';$HitStpPosArrayString=ArrayHashChange::ChangArrayIntoSimpleString ($HitStpPosArray) if (  ( ref ($HitStpPosArray) ) eq 'ARRAY'  );
			  	      my $HitTrpPosArrayString='';$HitTrpPosArrayString=ArrayHashChange::ChangArrayIntoSimpleString ($HitTrpPosArray) if (  ( ref ($HitTrpPosArray) ) eq 'ARRAY'  );
			  	      my $HitXxxPosArrayString='';$HitXxxPosArrayString=ArrayHashChange::ChangArrayIntoSimpleString ($HitXxxPosArray) if (  ( ref ($HitXxxPosArray) ) eq 'ARRAY'  );
			  	      #下面这行是为了收集 hit中的CYs，这里先去掉
			  	      #my $HitCysPosArrayString='';$HitCysPosArrayString=ArrayHashChange::ChangArrayIntoSimpleString ($HitCysPosArray) if (  ( ref ($HitCysPosArray) ) eq 'ARRAY'  );
			  	      print "\n 20190406-0-3-0 ", TimeWork::GetNowTimePid_microSecond ()," \$HitStpPosArray=$HitStpPosArray=$HitStpPosArrayString \$HitTrpPosArray=$HitTrpPosArray=$HitTrpPosArrayString \$HitXxxPosArray=$HitXxxPosArray=$HitXxxPosArrayString\n";
			  	  	  #print "\n 20190406-0-3-0 \$HitStpPosArray=$HitStpPosArray=$HitStpPosArrayString \$HitTrpPosArray=$HitTrpPosArray=$HitTrpPosArrayString \$HitCysPosArray=$HitCysPosArray=$HitCysPosArrayString\n";
			  	  	  
			  	  	  my @all_StopKind_Hit_posArray; 
			  		    push @all_StopKind_Hit_posArray,  @{ $HitStpPosArray } if (  ( ref ($HitStpPosArray) ) eq 'ARRAY'  );
			  	  	  push @all_StopKind_Hit_posArray,  @{ $HitTrpPosArray } if (  ( ref ($HitTrpPosArray) ) eq 'ARRAY'  );
			  	  	  push @all_StopKind_Hit_posArray,  @{ $HitXxxPosArray } if (  ( ref ($HitXxxPosArray) ) eq 'ARRAY'  );
			  	  	  #下面这行是为了收集 hit中的CYs，这里先去掉
			  	  	  #push @all_StopKind_Hit_posArray,  @{ $HitCysPosArray } if (  ( ref ($HitCysPosArray) ) eq 'ARRAY'  );
			  	  	  
			  	  	  
			  	  	  
			  	  	  my $spcialSTOPfound=0;  
			  	  	  foreach my $eachHitStpPos (  @all_StopKind_Hit_posArray  ){
			  	  	  	                                                                                                                                                                                       
			  	  	  	if (   
			  	  	  	       (  defined ( $qry_AlnPos_to_Char_HASH->{$eachHitStpPos} )   ) 
			  	  	  	       && 
			  	  	  	       (  defined ( $blastParsedHash->{'_ResultArray'}->[$rstIdx]->{'_hitArray'}->[$hitIdx]->{'_hspArray'}->[$hspIdx]->{'_Hsp_4Hit_1Utga_MatchToQryUC'} )  )
			  	  	  	       && 
			  	  	  	       (  defined ( $blastParsedHash->{'_ResultArray'}->[$rstIdx]->{'_hitArray'}->[$hitIdx]->{'_hspArray'}->[$hspIdx]->{'_Hsp_4Hit_1Utga_MatchToQryUC'}->{ $qry_AlnPos_to_Char_HASH->{$eachHitStpPos} } )  )   
			  	  	  	   ){ } #表示这个*已经被记录为和query的U或C能够match上的U了
			  	  	  	else {
			  	  	  		my $tp3DNAcodon=uc ( $hit_AlnPos_to_OrgPos_HASH->{$eachHitStpPos}->{'3DNAchar'} ); 
			  	  	  		#B J Z X
			  	  	  		if ($tp3DNAcodon eq 'TGA'){   # Z
			  	  	  			substr($otherSTOPchangeString, ($eachHitStpPos-1),1)='Z';                   print "\n 20190406-0-3-1 ", TimeWork::GetNowTimePid_microSecond ()," \$eachHitStpPos=$eachHitStpPos \$tp3DNAcodon=$tp3DNAcodon \n";
			  	  	  			$blastParsedHash->{'_ResultArray'}->[$rstIdx]->{'_hitArray'}->[$hitIdx]->{'_hspArray'}->[$hspIdx]->{'_Hsp_4Hit_2Utga_noMchToQryUC'}->{ $qry_AlnPos_to_Char_HASH->{$eachHitStpPos} }->{$eachHitStpPos}=$hit_AlnPos_to_OrgPos_HASH->{$eachHitStpPos};	
			  	  	  			
			  	  	  			#下面是用来在以后的分析中，如基因结构识别等步骤中，将TGA转变成TGC的
			  	  	  			my $hitOrgPos=$hit_AlnPos_to_OrgPos_HASH->{$eachHitStpPos}->{'real_1_3_Pos'}->{'3'};
			  	      			$blastParsedHash->{'_ResultArray'}->[$rstIdx]->{'_hitArray'}->[$hitIdx]->{'_TGAchangeHASH'}->{$hitOrgPos}->{'0_directi'}=$HitZoF;	
			  	      			$blastParsedHash->{'_ResultArray'}->[$rstIdx]->{'_hitArray'}->[$hitIdx]->{'_TGAchangeHASH'}->{$hitOrgPos}->{'1_postion'}=$hitOrgPos;			  	    			  
			  	      			$blastParsedHash->{'_ResultArray'}->[$rstIdx]->{'_hitArray'}->[$hitIdx]->{'_TGAchangeHASH'}->{$hitOrgPos}->{'2_fromCha'}='A';
			  	      			$blastParsedHash->{'_ResultArray'}->[$rstIdx]->{'_hitArray'}->[$hitIdx]->{'_TGAchangeHASH'}->{$hitOrgPos}->{'3_intoCha'}='C';
			  	      			$blastParsedHash->{'_ResultArray'}->[$rstIdx]->{'_hitArray'}->[$hitIdx]->{'_TGAchangeHASH'}->{$hitOrgPos}->{'4_sourcet'}='tBlastn, a of tga matched to '.$qry_AlnPos_to_Char_HASH->{$eachHitStpPos};
			  	  	  			
			  	  	  			
			  	  	  			
                     $spcialSTOPfound++;
			  	  	  		}
			  	  	  		elsif($tp3DNAcodon eq 'TAG'){  #O
			  	  	  			substr($otherSTOPchangeString, ($eachHitStpPos-1),1)='O';                   print "\n 20190406-0-3-2 ", TimeWork::GetNowTimePid_microSecond ()," \$eachHitStpPos=$eachHitStpPos \$tp3DNAcodon=$tp3DNAcodon \n";
			  	  	  			$blastParsedHash->{'_ResultArray'}->[$rstIdx]->{'_hitArray'}->[$hitIdx]->{'_hspArray'}->[$hspIdx]->{'_Hsp_4Hit_3Otag_1SpecialStop'}->{ $qry_AlnPos_to_Char_HASH->{$eachHitStpPos} }->{$eachHitStpPos}=$hit_AlnPos_to_OrgPos_HASH->{$eachHitStpPos};
			  	  	  			$spcialSTOPfound++;
			  	  	  		}
			  	  	  		elsif($tp3DNAcodon eq 'TAA'){  #B  
			  	  	  			substr($otherSTOPchangeString, ($eachHitStpPos-1),1)='B';                   print "\n 20190406-0-3-3 ", TimeWork::GetNowTimePid_microSecond ()," \$eachHitStpPos=$eachHitStpPos \$tp3DNAcodon=$tp3DNAcodon \n";            
			  	  	  			$blastParsedHash->{'_ResultArray'}->[$rstIdx]->{'_hitArray'}->[$hitIdx]->{'_hspArray'}->[$hspIdx]->{'_Hsp_4Hit_4Btaa_2SpecialStop'}->{ $qry_AlnPos_to_Char_HASH->{$eachHitStpPos} }->{$eachHitStpPos}=$hit_AlnPos_to_OrgPos_HASH->{$eachHitStpPos};
			  	  	  			$spcialSTOPfound++;
			  	  	  		}
			  	  	  		
			  	  	  	}
			  	  	  }
			  	  	  if ($spcialSTOPfound>0){                                                          print "\n 20190406-0-4-0 ", TimeWork::GetNowTimePid_microSecond ()," \$spcialSTOPfound=$spcialSTOPfound \n";          
			  	  	    $blastParsedHash->{'_ResultArray'}->[$rstIdx]->{'_hitArray'}->[$hitIdx]->{'_hspArray'}->[$hspIdx]->{'_Hsp_3WhSting_5Sh'}=  $otherSTOPchangeString;
			  	  	  }
			  	      
		             
		            #my $HitXxxPosArray; if ($hit_seq_including_X){$HitCysPosArray=&findingIdexOfAword ( 'X', $HitWholeString );}  #寻找Hit  中的所有X的位置  find the postions of all of the X in hit
			  	  	  #my $HitXxxPosArrayString='';$HitXxxPosArrayString=ArrayHashChange::ChangArrayIntoSimpleString ($HitXxxPosArray) if (  ( ref ($HitXxxPosArray) ) eq 'ARRAY'  );
			  	      #print "\n 20190406-0-5-0 \$HitXxxPosArray=$HitXxxPosArray=$HitXxxPosArrayString \n";
			  	  	  
			  	  	  #my @all_UknwKind_Hit_posArray;
			  	  	  #push @all_UknwKind_Hit_posArray,  @{ $HitXxxPosArray } if (  ( ref ($HitXxxPosArray) ) eq 'ARRAY'  );
			  	  	  #foreach my $eachHitUkwPos (  @all_UknwKind_Hit_posArray  ){
			  	  	  #  if (  ( $hit_AlnPos_to_OrgPos_HASH->{$eachHitUkwPos} eq 'X' ) || ( $hit_AlnPos_to_OrgPos_HASH->{$eachHitUkwPos} eq 'x' )  ){
			  	  	  #	  $blastParsedHash->{'_ResultArray'}->[$rstIdx]->{'_hitArray'}->[$hitIdx]->{'_hspArray'}->[$hspIdx]->{'_Hsp_4Hit_5Xxxx_noStdrtOthCd'}->{ $qry_AlnPos_to_Char_HASH->{$eachHitUkwPos} }->{$eachHitUkwPos}=$hit_AlnPos_to_OrgPos_HASH->{$eachHitUkwPos};			  	      				
			  	  	  #  }
			  	  	  #} 
		          	
		          	
		          	
		          }
			  	    
			  	    
			  	    
			  	    

		           
		                    
		  	    
		  	    
		  	    }
		  	  }
		  	
		  	}
		  }
		  
		  
		}		
	}
	return $blastParsedHash;
	
}

sub GetGiNumber{
	my ($inString)=@_;
	my $outString;
	#gi|377656257|pdb|3QPM|A
	if ($inString=~m/^(gi\|\d+)\|.*/) {
		$outString=$1;
	} 
}


sub BioPerlBlastPraser20170325{       #    BlastHandle::BioPerlBlastPraser20170325 ($blastfile, $outBestBlastFile, $numHits, $report_type, $INPUTdatabasePATH ) #pmAble#   #解析Blast文件
  my ($blastfile, $outBestBlastFile, $numHits, $report_type, $INPUTdatabasePATH )=@_;   #my ($blastfile, $allSelProInformHASHinput, $numHits)=@_;    
  
  
  my $warnMsgBody="\nIn package  BlastHandle,\tIn sub BioPerlBlastPraser20170325,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $subCalInfom=DirFileHandle::print_SubCallerInform;
  
  #warn  "\n\nNow in Sub &BioPerlpraseBlast!\n";      warn  "Input 1:\$blastfile=$blastfile\t" if (defined ($blastfile)); warn  "Input 2:\$outBestBlastFile=$outBestBlastFile\t" if (defined ($outBestBlastFile)); warn  "Input 3:\$numHits=$numHits\t" if (defined ($numHits));  warn  "Input 4:\$report_type=$report_type\t" if (defined ($report_type));  warn  "\n\n";
  print "\n\cl\nNow in Sub &BioPerlpraseBlast!\n:";  print "Input 1:\$blastfile=$blastfile\t" if (defined ($blastfile)); print "Input 2:\$outBestBlastFile=$outBestBlastFile\t" if (defined ($outBestBlastFile)); print "Input 3:\$numHits=$numHits\t" if (defined ($numHits));  print "Input 4:\$report_type=$report_type\t" if (defined ($report_type));  print "\n\cl\n";
  
  my $in;
  if (   ( defined ($report_type) )  && ($report_type=~m/^\s*tblastn\s+$/i)   ){ 
  	$report_type='tblastn';
    $in = Bio::SearchIO->new(-format => 'blast', -file => $blastfile, -report_type => $report_type); warn "1\n";
  }
  else {
  	$in = Bio::SearchIO->new(-format => 'blast', -file => $blastfile                              ); warn "2\n";
  }
  my $prasedOUThash;
  my $resultArrayIdx=0;
  my $maxBestinformToget=5;  #这个数字是用来 限制， $summuryHere这个变量的长度的
  
  
  my $QueryHit_excel_HASH;
  
  my $SortArray;  my $SortArrayIdx=0;
  my $totalHit=0;
  
  # extraction of information for each result recursively
  while ( my $result = $in->next_result ) {
	  # the name of the query sequence    	
	  print  $result->query_name . "\t";
	  my $queryMinStart =  99999999999999999999;
    my $queryMaxEnd   = -99999999999999999999;
	  my $queryName=$result->query_name;
	  #if (  defined ( $allSelProInformHASHinput->{$queryName} )  ){ }
   	#else { my $find=0;
   	#  foreach my $qrNmKey (   keys (  %{ $allSelProInformHASHinput }  )   ){
   	#    if ($qrNmKey=~m/^$queryName/){$queryName=$qrNmKey;   $find=1;    }
   	#  }
   	#  if ($find==0){
   	#    die "\nNO \$queryName=$queryName was found and Changed:\$blastfile=$blastfile\n\n";
   	#  }
   	#}
   	$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_query_name'}=$queryName; 
   	
   	my $databaseFind=$result->database_name;                                                                                                   #print  "\n20190322-0-1 database=".$databaseFind."\n";
   	if (   (  defined ( $databaseFind )  ) && ( $databaseFind=~m/\S+/)   ){  #&& (  -e ( $databaseFind )  )    )  
   		                                                                                                                                       #print  "\n20190322-0-2 database=".$databaseFind."\n"; 		
   	} 
   	else {
   		                                                                                                                                        #print  "\n20190322-0-3 database=".$databaseFind."\t\$INPUTdatabasePATH=$INPUTdatabasePATH\n";
   	  if (   (  defined ( $INPUTdatabasePATH )  ) && ( $INPUTdatabasePATH=~m/\S+/ )    ){ # && ( $INPUTdatabasePATH=~m/\S+/ ) && ( -e $INPUTdatabasePATH  )
   	  	$databaseFind=$INPUTdatabasePATH;                                                                                                    #print  "\n20190322-0-4 database=".$databaseFind."\n";
   	  }
   	}
   	
   	
   	$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_database_name'}=$databaseFind;                                                    #print  "\n20190322-0-5 database=".$databaseFind."\n";
    
    # the length of the query sequence #    	
    print  $result->query_length;
    my $query_length=$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_query_length'}=$result->query_length;
    
    $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_num_hits'}=$result->num_hits;  print  $result->num_hits;
    # output "no hits found" if there is no hits
    $totalHit+=$result->num_hits;  print "\n\$totalHit=$totalHit\t\$result->num_hits=$result->num_hits\n";
    
    my $summuryHere='';  my $summuryHereLength=0;
    if ( $result->num_hits == 0 ) {
		  print  "\tNo hits found\n";
    } 
    else {
		  my $count = 0;
      # process each hit recursively
		  while (my $hit = $result->next_hit) {
			  print  "\t" if ($count > 0);
                          # get the accession numbers of the hits
			  print  "\t" . $hit->accession . "\t";
                          # get the lengths of the hit sequences
        print  $hit->length . "\t";
                          # get the description of the hit sequences
			  print  $hit->description . "\t";
                          # get the E value of the hit
			  print  $hit->significance . "\t";
                          #get the bit score of the hit
			  print  $hit->bits . "\t";
        
        my $HitName=$hit->accession;
	      #if (  defined ( $allSelProInformHASHinput->{$HitName} )  ){ }
   	    #else { my $find=0;
   	    #  foreach my $qrNmKey (   keys (  %{ $allSelProInformHASHinput }  )   ){
   	    #    if ($qrNmKey=~m/^$HitName/){$HitName=$qrNmKey;   $find=1;    }
   	    #  }
   	    #  if ($find==0){
   	    #    die "\nNO \$HitName=$HitName was found and Changed:\$blastfile=$blastfile\n\n";
   	    #  }
       	#}
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'0_0_0_locus'}   =$hit->locus();        
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'0_0_1__name'}   =$hit->name();  
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'0_0_2_giNub'}   =&GetGiNumber( $hit->name() );  
           
        #$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_accessionNB'}  =$hit->accession;
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_accessionNB'}   =$HitName;
        
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hitSeqLength'}  =$hit->length;
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hitDescription'}=$hit->description;
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hitEvalue'}     =$hit->significance;
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hitScore'}      =$hit->bits;
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_algorithm'}      =$hit->algorithm;
        
        
        my $accIDinShort=$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_accessionNB'};
        $accIDinShort=~s/^(\S+)\s*.*$/$1/;  #214g2702C.paradoxa20150522;thioredoxin_reductase2;Contig53901_-3

        
        my $hspcount = 0; 
        my $totalIdentical=0;
        my $totalConserved=0;
        
        my $totConsePercent=0;  #
        
        
        my $queryCoverHash;  
        my $hitCoverHash;
        
        my $queryCoverRate=0;  
        my $qureyCovRatePct=0;
        my $total_qry_covrAte=0;
        my $total_qry_covrAte_pct=0;
        my $hitCoverRate=0;
        my $hitCovRatePct=0;
        
                          # process the top HSP for the top hit
			  while (my $hsp = $hit->next_hsp) {
			  	
			  	
			  	
          print  "\t\t\t\t\t\t\t", if ($hspcount > 0);
                          	# get the frame of the query sequence
			  	print  $hsp->query->frame . "\t";
                                  # get the start and the end of the query sequence in the alignment
			  	print  $hsp->start('query') . "\t" . $hsp->end('query'). "\t";
                                  # get the start and the end of the hit sequence in the alignment
			  	print  $hsp->start('hit') . "\t" . $hsp->end('hit') . "\t";
                                  # get the similarity value
			  	printf  "%.1f" , ($hsp->frac_conserved * 100);
			  	print  "%\t";
                                  # get the identity value
			  	printf  "%.1f" , ($hsp->frac_identical * 100);
			  	
			  	#$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_HspQurFrame'}        =$hsp->frame('query');  #oldVersion   #bioperl 得到的blast的hsp的frame是不对的
			  	#$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_HspHitFrame'}        =$hsp->frame('hit');    #oldVersion   #bioperl 得到的blast的hsp的frame是不对的
			  	
			  	
			  	
			  	#$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_Hsp_1stfm_3frm_2qur'}       =$hsp->frame('query');        my $QrAbsFm=$hsp->frame('query')+1;  #bioperl 得到的blast的hsp的frame是不对的
			  	#$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_Hsp_1stfm_4frm_2hit'}       =$hsp->frame('hit');          my $HtAbsFm=$hsp->frame('hit')+1;   #bioperl 得到的blast的hsp的frame是不对的
			  	$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_Hsp_1stfm_1stn_1qur'}       =$hsp->strand('query');       my $QrZF; if ($hsp->strand('query') >=0){$QrZF='+';}else {$QrZF='-';} 
			  	$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_Hsp_1stfm_2stn_1hit'}       =$hsp->strand('hit');         my $HtZF; if ($hsp->strand('hit')   >=0){$HtZF='+';}else {$HtZF='-';} 
			  	#$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_Hsp_1stfm_1qur'}       =$QrZF.$QrAbsFm;   #bioperl 得到的blast的hsp的frame是不对的
			  	#$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_Hsp_1stfm_2hit'}       =$HtZF.$HtAbsFm;   #bioperl 得到的blast的hsp的frame是不对的
			  	$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_Hsp_1stfm_3ZoF_1qur'}       =$QrZF;       
			  	$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_Hsp_1stfm_4ZoF_1hit'}       =$HtZF;         
			  	
			  	$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_HspQueryStart'}      =$hsp->start('query');  #oldVersion
			  	$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_HspQueryEnd'}        =$hsp->end('query');    #oldVersion
			  	$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_HspHitStart'}        =$hsp->start('hit');    #oldVersion
			  	$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_HspHitEnd'}          =$hsp->end('hit');      #oldVersion
			  				  	
			  	$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_Hsp_2StEd_1Qr_1St'}        =$hsp->start('query');       $queryCoverHash->[$hspcount]->{'_start'}=$hsp->start('query'); 
			  	$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_Hsp_2StEd_1Qr_2Ed'}        =$hsp->end('query');         $queryCoverHash->[$hspcount]->{'_end'}  =$hsp->end('query');   
			  	$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_Hsp_2StEd_2Ht_1St'}        =$hsp->start('hit');         $hitCoverHash->[$hspcount]->{'_start'}  =$hsp->start('hit');     
			  	$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_Hsp_2StEd_2Ht_2Ed'}        =$hsp->end('hit');           $hitCoverHash->[$hspcount]->{'_end'}    =$hsp->end('hit');       
			  	

			  	
			  	$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_HspConserverd'}      =$hsp->frac_conserved ;
			  	$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_HspIdentical'}       =$hsp->frac_identical ;
			  	



			  	
			  	#$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_Hsp_1Qry_string'}      =$hsp->query_string ;    #oldVersion
			  	#$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_Hsp_3Hit_string'}      =$hsp->hit_string ;      #oldVersion
			  	#$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_Hsp_2Hom_string'}      =$hsp->homology_string ; #oldVersion
			  	
			  	$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_Hsp_3WhSting_1Qr'}      =$hsp->query_string ;
			  	$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_Hsp_3WhSting_2Ho'}      =$hsp->homology_string ;
			  	$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_Hsp_3WhSting_3Ht'}      =$hsp->hit_string ;
			  	
			    
			    
			    

			  


			  	if (0){   #！这里的 if(0)是将下面的代码 无效化了     	
			  	   
			  	    	#下面分别是 考虑用  seq_inds  和  get_aln 来解决 query和 hit之间位置一一对应的功能，但都不太好用。
			  	    	if(0){#！这里的 if(0)是将下面的代码 无效化了     
			  	    		# seq_ind的问题是，得到了各种match的id，但没有query和hit之间的对应关系
			  	    	  my ($conQryIns, $conHitIns, $mmsQryIns, $mmsHitIns, $gapQryIns, $gapHitIns, $fsfQryIns, $fsfHitIns );
			  	        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_Hsp_4ids_1idt_1qry'}           =[$hsp->seq_inds('query', 'identical')];
			  	        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_Hsp_4ids_1idt_2hit'}           =[$hsp->seq_inds('hit',   'identical')];  #'gap'
			  	        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_Hsp_4ids_2con_1qry'}=$conQryIns=[$hsp->seq_inds('query', 'conserved')];
			  	        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_Hsp_4ids_2con_2hit'}=$conHitIns=[$hsp->seq_inds('hit',   'conserved')];  #'gap'
			  	        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_Hsp_4ids_3msm_1qry'}=$mmsQryIns=[$hsp->seq_inds('query', 'mismatch')];
			  	        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_Hsp_4ids_3msm_2hit'}=$mmsHitIns=[$hsp->seq_inds('hit',   'mismatch')];   #
			  	        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_Hsp_4ids_4gap_1qry'}=$gapQryIns=[$hsp->seq_inds('query', 'gap')];
			  	        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_Hsp_4ids_4gap_2hit'}=$gapHitIns=[$hsp->seq_inds('hit',   'gap')];   #conserved
			  	        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_Hsp_4ids_5fsf_1qry'}=$fsfQryIns=[$hsp->seq_inds('query', 'frameshift')];
			  	        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_Hsp_4ids_5fsf_2hit'}=$fsfHitIns=[$hsp->seq_inds('hit',   'frameshift')];   #conserved
                  
                  my ($conQryInsHash, $conHitInsHash, $mmsQryInsHash, $mmsHitInsHash, $gapQryInsHash, $gapHitInsHash, $fsfQryInsHash, $fsfHitInsHash );
                  $conQryInsHash=ArrayHashChange::ChangeArrayToHash($conQryIns); $mmsQryInsHash=ArrayHashChange::ChangeArrayToHash($mmsQryIns); $gapQryInsHash=ArrayHashChange::ChangeArrayToHash($gapQryIns); $fsfQryInsHash=ArrayHashChange::ChangeArrayToHash($fsfQryIns); 
                  $conHitInsHash=ArrayHashChange::ChangeArrayToHash($conHitIns); $mmsHitInsHash=ArrayHashChange::ChangeArrayToHash($mmsHitIns); $gapHitInsHash=ArrayHashChange::ChangeArrayToHash($gapHitIns); $fsfHitInsHash=ArrayHashChange::ChangeArrayToHash($fsfHitIns); 
                  if ( ArrayHashChange::Check4HashForReapeatKey($conQryInsHash, $mmsQryInsHash, $gapQryInsHash, $fsfQryInsHash) ){
                  	$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_Hsp_4ids_5DIE_1qry'}=ArrayHashChange::Check4HashForReapeatKey_showMSG($conQryInsHash, $mmsQryInsHash, $gapQryInsHash, $fsfQryInsHash); 
                  	
                  }
                  if ( ArrayHashChange::Check4HashForReapeatKey($conHitInsHash, $mmsHitInsHash, $gapHitInsHash, $fsfHitInsHash) ){
                  	$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_Hsp_4ids_5DIE_1hit'}=ArrayHashChange::Check4HashForReapeatKey_showMSG($conHitInsHash, $mmsHitInsHash, $gapHitInsHash, $fsfHitInsHash); 
                  }
                  
			  	    	  # get_aln，得到了query和hit之间的对应关系，但是
			  	    	  $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_Hsp_aln'}              = $hsp->get_aln;
			  	    	  #$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_Hsp_4HtStWrdArr'}      = $SomWordArray;
			  	    	  
			  	    	  #foreach my $wdDB (  @{ $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_Hsp_4HtStWrdArr'} }  ){
			  	    	  #	my $qrSeqObj=$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_Hsp_aln'}->{'_seq'}->{ $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_Hsp_aln'}->{'_order'}->{0} };
			  	    	  #	my $sbSeqObj=$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_Hsp_aln'}->{'_seq'}->{ $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_Hsp_aln'}->{'_order'}->{1} };
			  	    	  #	warn "\$wdDB=$wdDB\n";
			  	    	  #	warn "\$qrSeqObj->id=", $qrSeqObj->column_from_residue_number($qrSeqObj->{'start'}+$wdDB-1), "\n";
                  #  my $qloc = $qrSeqObj->location_from_column($wdDB);
                  #  my $sloc = $sbSeqObj->location_from_column($wdDB);
			  	    	  #  push (   (  @{ $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_Hsp_4QQQQWrdArr'} }  ),  $qloc   );
			  	    	  #  push (   (  @{ $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_Hsp_4HHHHWrdArr'} }  ),  $sloc   );
			  	    	  #}
			  	    	  #$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_Hsp_aln'}->{'_seq'}->{ $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_Hsp_aln'}->{'_order'}->{0} }->column_from_residue_number()
			  	      }
			  	      
			  	      
			  	}
			  	
			  	
			  	
			  	
			  	
		           

			  	
			  	
			  	
			  	$totalIdentical+=$hsp->frac_identical ;
			  	$totalConserved+=$hsp->frac_conserved ;
			  	
			  	my $smallQNB=$hsp->start('query'); my $bigQNB=$hsp->end('query');
			  	if ($smallQNB > $bigQNB ){ $bigQNB=$hsp->start('query'); my $smallQNB=$hsp->end('query'); }
			  	
			  	if ( $smallQNB < $queryMinStart ){  $queryMinStart = $smallQNB;   }
			  	if ( $bigQNB   > $queryMaxEnd   ){  $queryMaxEnd   = $bigQNB;     }
			  	  
		      print  "%\n";
          $hspcount++;
        }
        
        $queryCoverHash=&overLayerHadle( $queryCoverHash );  
        $hitCoverHash  =&overLayerHadle( $hitCoverHash ); 
        
        my $queryCoverLength=&GetCuttedLength($queryCoverHash);
        my $hitCoverLength=&GetCuttedLength($hitCoverHash);
        
        my ($hit_queryStart, $hit_queryEnd)=@{ &getWholeMatchedStartEnd($queryCoverHash) };
        my ($hit_HitStart,   $hit_HitEnd)  =@{ &getWholeMatchedStartEnd($hitCoverHash)   };
        
        
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'queryCoverHash'}  =$queryCoverHash;        
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'queryCoverLength'}=$queryCoverLength;
        
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'hitCoverHash'}    =$hitCoverHash;
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'hitCoverLength'}  =$hitCoverLength;
       
        
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'hit_queryStart'}      =$hit_queryStart;
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'hit_queryEnd'}        =$hit_queryEnd;
        my $hit_queryRangeLength =(  (abs ($hit_queryEnd-$hit_queryStart)  ) + 1   );
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'hit_queryRangeLength'}=$hit_queryRangeLength;
        
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'hit_HitStart'}      =$hit_HitStart;
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'hit_HitEnd'}        =$hit_HitEnd;
        my $hit_HitRangeLength   =(  (abs ($hit_HitEnd-$hit_HitStart)  ) + 1   );
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'hit_HitRangeLength'}=$hit_HitRangeLength;
        
        
        
        if ($hspcount > 0){
          $totalIdentical=$totalIdentical/$hspcount;
          $totalConserved=$totalConserved/$hspcount;
          $totConsePercent=100*$totalConserved;  $totConsePercent=sprintf "%.2f",$totConsePercent; $totConsePercent="$totConsePercent%";
          
          $queryCoverRate=$queryCoverLength/$hit_queryRangeLength;  
          $qureyCovRatePct   =100*$queryCoverRate;  $qureyCovRatePct=sprintf "%.2f",$qureyCovRatePct; $qureyCovRatePct="$qureyCovRatePct%";
          
          $total_qry_covrAte=$queryCoverLength/$query_length;
          $total_qry_covrAte_pct=100*$total_qry_covrAte;  $total_qry_covrAte_pct=sprintf "%.2f",$total_qry_covrAte_pct; $total_qry_covrAte_pct="$total_qry_covrAte_pct%";
          
        
          $hitCoverRate=$hitCoverLength/$hit_HitRangeLength;
          $hitCovRatePct     =100*$hitCoverRate;  $hitCovRatePct=sprintf "%.2f",$hitCovRatePct; $hitCovRatePct="$hitCovRatePct%";
          
          
        }
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hitTotalIdentical'}=$totalIdentical;
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hitTotalConserved'}=$totalConserved;
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'totConsePercent'   }=$totConsePercent;
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'queryCoverRate'    }=$queryCoverRate;
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'qureyCovRatePct'   }=$qureyCovRatePct;
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'total_qry_covrAte' }=$total_qry_covrAte;
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'total_qry_covrAte_pct'}=$total_qry_covrAte_pct;
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'hitCoverRate'      }=$hitCoverRate;
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'hitCovRatePct'     }=$hitCovRatePct;
        
        
        if ($summuryHereLength<$maxBestinformToget){
        	my $singleSum="${accIDinShort}($prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hitEvalue'}|$totConsePercent)";
          $summuryHere.=$singleSum;
          $summuryHereLength++;
          $QueryHit_excel_HASH->{$queryName}->{$HitName}=$singleSum;
        }
        
        $SortArray->{'EvalueSort'}->[$SortArrayIdx]->{'SotVal'}=$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hitEvalue'};
        $SortArray->{'EvalueSort'}->[$SortArrayIdx]->{'QurIdx'}=$resultArrayIdx;
        $SortArray->{'EvalueSort'}->[$SortArrayIdx]->{'HitIdx'}=$count;
        
        $SortArray->{'TotCosSort'}->[$SortArrayIdx]->{'SotVal'}=$totalConserved;
        $SortArray->{'TotCosSort'}->[$SortArrayIdx]->{'QurIdx'}=$resultArrayIdx;
        $SortArray->{'TotCosSort'}->[$SortArrayIdx]->{'HitIdx'}=$count;
        
        $SortArray->{'TotIdtSort'}->[$SortArrayIdx]->{'SotVal'}=$totalIdentical;
        $SortArray->{'TotIdtSort'}->[$SortArrayIdx]->{'QurIdx'}=$resultArrayIdx;
        $SortArray->{'TotIdtSort'}->[$SortArrayIdx]->{'HitIdx'}=$count;
        $SortArrayIdx++; 
        
			  $count++;
        
        # flow control for the number of hits needed
        if (  (defined ($numHits) ) && ($numHits=~m/\d+/) && ($numHits>0)   ){
			    last if ($count == $numHits);
			  }
		  }
    }                                                       
                                                         
    $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_QueryMinStart'}=$queryMinStart;
    $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_QueryMaxEnd'}  =$queryMaxEnd;
    $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_Summury'}      =$summuryHere;
    $resultArrayIdx++;	
  }

  my $realSortArray; 
  foreach my $EachSortType (    sort {$a cmp $b} (   keys (  %{ $SortArray }  )   )    ){
  	my $realSortIdx=0;  my $totalSummury;
  	if ($EachSortType eq 'EvalueSort'){
      foreach my $eachHash  (   sort { $a->{'SotVal'} <=> $b->{'SotVal'} }  (  @{ $SortArray->{$EachSortType} }  )   ){
      	$realSortArray->{$EachSortType}->[$realSortIdx]=$eachHash;      $realSortIdx++;      	    	
      }
    }
    else{
    	foreach my $eachHash  (   sort { $b->{'SotVal'} <=> $a->{'SotVal'} }  (  @{ $SortArray->{$EachSortType} }  )   ){
      	$realSortArray->{$EachSortType}->[$realSortIdx]=$eachHash;      $realSortIdx++;
      }
    }
  } 
  

  $prasedOUThash->{'SortHash'}=$realSortArray;
  $prasedOUThash->{'totalHit'}=$totalHit;
  
  my $beastQuery=$prasedOUThash->{'_ResultArray'}->[ $prasedOUThash->{'SortHash'}->{'EvalueSort'}->[0]->{'QurIdx'} ]->{'_query_name'};  #warn "\nbbbbbbb\$beastQuery=$beastQuery\n\n"; print "\nbbbbbbb\$beastQuery=$beastQuery\n\n";
  
##  my $writerhtml = new Bio::SearchIO::Writer::HTMLResultWriter(); 
##my $outhtml = new Bio::SearchIO(-writer => $writerhtml,
##   -file   => ">searchio.html");
##
### Loop through all the results, successively adding each one to 
### the bottom of the HTML report
##while ( $result = $searchio->next_result() ) {  
##    $outhtml->write_report($result);
##}

  #DirFileHandle::PrintDumper                    ('temp.QueryHitHash', $QueryHit_excel_HASH)     if ( ref($QueryHit_excel_HASH) eq 'HASH' ) ; 
  #MatrixCsvChange::PrintExcelTableFileFromHASH  ('temp.QueryHit.excel.txt', $QueryHit_excel_HASH, 0 ) if ( ref($QueryHit_excel_HASH) eq 'HASH' ) ; 
  
  if (   (  defined ($outBestBlastFile) ) && ($outBestBlastFile=~m/\S+/)   ){
    my $writerhtml = new Bio::SearchIO::Writer::HTMLResultWriter();      #Bio::SearchIO->new(-format => 'blast', -file => ">$outBestBlastFile");
    my $out = new Bio::SearchIO( -writer => $writerhtml,
                                 -file   => ">$outBestBlastFile");
  
    my $in_2 = Bio::SearchIO->new(-format => 'blast', -file => $blastfile);
  
    OUTBEST: while ( my $result = $in_2->next_result ) {  
      if ($beastQuery eq $result->query_name){       #print "\nbbbbbbb\$beastQuery=$beastQuery\t\t\$result->query_name=",$result->query_name,"\n\n"; warn "\nbbbbbbb\$beastQuery=$beastQuery\t\t\$result->query_name=",$result->query_name,"\n\n";
  	    $out->write_report($result); last OUTBEST;
  	  }
    }
  }
  print "\n\clEnd of Sub &BioPerlpraseBlast!\n\n";
  return $prasedOUThash;
 
  
  

}

sub GetAllHit_Array_fromBlastRstArrayHash{   # my $AllHit_Array=BlastHandle::GetAllHit_Array_fromBlastRstArrayHash( $blastRstARRAYhash );
	my ($blastRstARRAYhash)=@_;
	
	my $warnMsgBody="\nIn package  BlastHandle,\tIn sub GetAllHit_Array_fromBlastRstArrayHash,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  my $subCallereIfm=DirFileHandle::print_SubCallerInform; #warn $subCallereIfm;
  
  my $AllHit_Hash; my $AllHit_Array;
  if (   (  defined ( $blastRstARRAYhash )  ) && (  ref ( $blastRstARRAYhash ) eq 'HASH'  ) && (  defined ( $blastRstARRAYhash->{'_ResultArray'} )  ) && (  ref ( $blastRstARRAYhash->{'_ResultArray'} ) eq 'ARRAY'  )   ){
  	my $QueryIdx=0;                                                                                                                                                                                                                           #print "20190325-0-0-0-1 \$QueryIdx=$QueryIdx\n";
  	foreach my $eachQueryHASH (  @{ $blastRstARRAYhash->{'_ResultArray'} }  ){                                                                                                                                                                #print "20190325-0-0-0-2 \$QueryIdx=$QueryIdx\n";
  		if (   (  defined ( $eachQueryHASH )  ) && (  ref ( $eachQueryHASH ) eq 'HASH'  ) && (  defined ($eachQueryHASH->{'_query_name'}) ) && ($eachQueryHASH->{'_query_name'}=~m/\S+/)   ){
  		  my $queryName=$eachQueryHASH->{'_query_name'};                                                                                                                                                                                        #print "20190325-0-0-0-3 \$queryName=\$eachQueryHASH->{'_query_name'}=$eachQueryHASH->{'_query_name'}\n";
  		  if (   (  defined ( $eachQueryHASH->{'_hitArray'} )  ) && (  ref ( $eachQueryHASH->{'_hitArray'} ) eq 'ARRAY'  )   ){
  		  	my $HitIdx=0;                                                                                                                                                                                                                       #print "20190325-0-0-0-4 \$HitIdx=$HitIdx\n";
  		  	foreach my $eachHitHASH (  @{ $eachQueryHASH->{'_hitArray'} }  ){                                                                                                                                                                   #print "20190325-0-0-0-5 \$HitIdx=$HitIdx\n";
  		  	  if (   (  defined ( $eachHitHASH )  ) && (  ref ( $eachHitHASH ) eq 'HASH'  ) && (  defined ($eachHitHASH->{'_accessionNB'}) ) && ($eachHitHASH->{'_accessionNB'}=~m/\S+/)   ){
  		  	  	my $hitName=$eachHitHASH->{'_accessionNB'};                                                                                                                                                                                     print "20190325-0-0-0-6 \$hitName=$hitName=\$eachHitHASH->{'_accessionNB'}=$eachHitHASH->{'_accessionNB'}\n";
  		  	  	
  		  	  	if (   (  defined ( $AllHit_Hash->{$hitName} )  ) && ( $AllHit_Hash->{$hitName} == 1 )   ){  		  	  		
  		  	  	}
  		  	  	else{
  		  	  		push @{ $AllHit_Array }, $hitName;                       print "20190325-0-0-0-7 \$hitName=$hitName\n";
  		  	  	  $AllHit_Hash->{$hitName}=1;                              print "20190325-0-0-0-8 \$AllHit_Hash->{$hitName}=$AllHit_Hash->{$hitName}\n"; 
  		  	  	}
  		  	  	
  		  	  	
  		  	  }
  		  	  $HitIdx++;
  		  	}
  		  }
  		  
  		}
  		$QueryIdx++;
  	}
  }
	return $AllHit_Array;
}

sub ChangeBlastResult_into_Qry_Hit_HASH{ # my $blastQueryHit_hash=BlastHandle::ChangeBlastResult_into_Qry_Hit_HASH( $blastRstARRAYhash );
	my ($blastRstARRAYhash)=@_;
	
	my $warnMsgBody="\nIn package  BlastHandle,\tIn sub ChangeBlastResult_into_Qry_Hit_HASH,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  my $subCallereIfm=DirFileHandle::print_SubCallerInform; #warn $subCallereIfm;
                  
  my $blastQueryHit_hash;                                                                                                                                                                                                                     #print "20190124-0-0-0-0 \$blastQueryHit_hash=$blastQueryHit_hash\n";
  if (   (  defined ( $blastRstARRAYhash )  ) && (  ref ( $blastRstARRAYhash ) eq 'HASH'  ) && (  defined ( $blastRstARRAYhash->{'_ResultArray'} )  ) && (  ref ( $blastRstARRAYhash->{'_ResultArray'} ) eq 'ARRAY'  )   ){
  	my $QueryIdx=0;                                                                                                                                                                                                                           #print "20190124-0-0-0-1 \$QueryIdx=$QueryIdx\n";
  	foreach my $eachQueryHASH (  @{ $blastRstARRAYhash->{'_ResultArray'} }  ){                                                                                                                                                                #print "20190124-0-0-0-2 \$QueryIdx=$QueryIdx\n";
  		if (   (  defined ( $eachQueryHASH )  ) && (  ref ( $eachQueryHASH ) eq 'HASH'  ) && (  defined ($eachQueryHASH->{'_query_name'}) ) && ($eachQueryHASH->{'_query_name'}=~m/\S+/)   ){
  		  my $queryName=$eachQueryHASH->{'_query_name'};                                                                                                                                                                                        #print "20190124-0-0-0-3 \$queryName=\$eachQueryHASH->{'_query_name'}=$eachQueryHASH->{'_query_name'}\n";
  		  if (   (  defined ( $eachQueryHASH->{'_hitArray'} )  ) && (  ref ( $eachQueryHASH->{'_hitArray'} ) eq 'ARRAY'  )   ){
  		  	my $HitIdx=0;                                                                                                                                                                                                                       #print "20190124-0-0-0-4 \$HitIdx=$HitIdx\n";
  		  	foreach my $eachHitHASH (  @{ $eachQueryHASH->{'_hitArray'} }  ){                                                                                                                                                                   #print "20190124-0-0-0-5 \$HitIdx=$HitIdx\n";
  		  	  if (   (  defined ( $eachHitHASH )  ) && (  ref ( $eachHitHASH ) eq 'HASH'  ) && (  defined ($eachHitHASH->{'_accessionNB'}) ) && ($eachHitHASH->{'_accessionNB'}=~m/\S+/)   ){
  		  	  	my $hitName=$eachHitHASH->{'_accessionNB'};                                                                                                                                                                                     print "20190124-0-0-0-6 \$hitName=$hitName=\$eachHitHASH->{'_accessionNB'}=$eachHitHASH->{'_accessionNB'}\n";
  		  	  	$blastQueryHit_hash->{$queryName}->{$hitName}=Storable::dclone( $eachHitHASH );
  		  	  	
  		  	  	$blastQueryHit_hash->{$queryName}->{$hitName}->{'0_0_0_0_0_0_Query___name'}=$queryName;
  		  	  	$blastQueryHit_hash->{$queryName}->{$hitName}->{'0_0_0_0_0_1_Hit_____name'}=$hitName;
  		  	  	$blastQueryHit_hash->{$queryName}->{$hitName}->{'0_0_0_0_1_Query___length'}=$eachQueryHASH->{'_query_length'};
  		  	  	$blastQueryHit_hash->{$queryName}->{$hitName}->{'0_0_0_0_2_Query__Min_Stt'}=$eachQueryHASH->{'_QueryMinStart'};
  		  	  	$blastQueryHit_hash->{$queryName}->{$hitName}->{'0_0_0_0_3_Query__Max_End'}=$eachQueryHASH->{'_QueryMaxEnd'};
  		  	  	$blastQueryHit_hash->{$queryName}->{$hitName}->{'0_0_0_0_4_Query_Database'}=$eachQueryHASH->{'_database_name'};
  		  	  	$blastQueryHit_hash->{$queryName}->{$hitName}->{'0_0_0_0_5_Query__Summury'}=$eachQueryHASH->{'_Summury'};
  		  	  	$blastQueryHit_hash->{$queryName}->{$hitName}->{'0_0_0_0_6_Query_hits_num'}=$eachQueryHASH->{'_num_hits'};
  		  	  	
  		  	  	$blastQueryHit_hash->{$queryName}->{$hitName}->{'0_0_0_1_0_Query_arreyIdx'}=$QueryIdx;
  		  	  	$blastQueryHit_hash->{$queryName}->{$hitName}->{'0_0_0_1_1_Hit___arreyIdx'}=$HitIdx;
  		  	  	
  		  	  }
  		  	  $HitIdx++;
  		  	}
  		  }
  		  
  		}
  		$QueryIdx++;
  	}
  }
	return $blastQueryHit_hash;
}
  
sub ChangeBlastResult_into_Hit_Qry_HASH{ # my $blastQueryHit_hash=BlastHandle::ChangeBlastResult_into_Hit_Qry_HASH( $blastRstARRAYhash );
	my ($blastRstARRAYhash)=@_;
	
	my $warnMsgBody="\nIn package  BlastHandle,\tIn sub ChangeBlastResult_into_Hit_Qry_HASH,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  my $subCallereIfm=DirFileHandle::print_SubCallerInform; #warn $subCallereIfm;
  
  my $blastQueryHit_hash=BlastHandle::ChangeBlastResult_into_Qry_Hit_HASH( $blastRstARRAYhash );
  my $blastHitQuery_hash;
  if (   (  defined ( $blastQueryHit_hash )  ) && (  ref ( $blastQueryHit_hash ) eq 'HASH'  )   ){
  	$blastHitQuery_hash=ArrayHashChange::Reverse_level1KEY_Level2Key_HASH($blastQueryHit_hash);
  }
  return $blastHitQuery_hash;
}

sub GetBest_Qry_or_hit_BlastResultS {
	my ($in_blast_hit_rst_file, $Qry_or_hit, $sortByWhat, $best_number, $sortCutoff)=@_;
	
	my $warnMsgBody="\nIn package  BlastHandle,\tIn sub ChangeBlastResult_into_Hit_Qry_HASH,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  my $subCallereIfm=DirFileHandle::print_SubCallerInform; #warn $subCallereIfm;
	
	my $Qry_Hit_Type;
	if (   (  defined ( $Qry_or_hit )  ) && ( $Qry_or_hit=~m/\S+/ ) && (  (  $Qry_or_hit eq 'query'  ) || (  $Qry_or_hit eq 'hit'  )  )   ){}
	else {	  	DieWork::Just_dieWork( $die_MsgHead.$subCallereIfm."\$Qry_or_hit=$Qry_or_hit is not right!!!!!\nIt is limited to: 1 query, 2 hit\n$!\n\n" );	  }
		
	
	if (       (  defined ( $sortByWhat )  ) && ( $sortByWhat=~m/\S+/ ) 
	       &&  (  (  $sortByWhat eq '_hitEvalue'  ) || (  $sortByWhat eq '_hitTotalConserved'  ) || (  $sortByWhat eq 'hitCoverRate'  ) || (  $sortByWhat eq 'total_qry_covrAte'  ) || (  $sortByWhat eq 'queryCoverRate'  )  )   
	   ){}
	else {	  	DieWork::Just_dieWork( $die_MsgHead.$subCallereIfm."\$sortByWhat=$sortByWhat is not right!!!!!\nIt is limited to: 1 _hitEvalue, 2 _hitTotalConserved, 3 hitCoverRate, 4 total_qry_covrAte, 5 queryCoverRate\n$!\n\n" );	  }
	

	
	my $sorted_cutted_HASH;
	if (   (  defined ( $in_blast_hit_rst_file )  ) && ( $in_blast_hit_rst_file=~m/\S+/ ) && (  -e ( $in_blast_hit_rst_file )  )   ){
	  my $inBlastRstArrayHash=BlastHandle::BioPerlBlastPraser20170325( $in_blast_hit_rst_file );
	  if (   (  defined ( $inBlastRstArrayHash )  ) && (  ref ( $inBlastRstArrayHash ) eq 'HASH'  )   ){  
	    my $QH_or_HQ_hash;
	    if    (  $Qry_or_hit eq 'query'  ){
	    	$QH_or_HQ_hash=BlastHandle::ChangeBlastResult_into_Qry_Hit_HASH( $inBlastRstArrayHash );  $Qry_Hit_Type='Qry_Hit';
	    }
	    elsif (  $Qry_or_hit eq 'hit'  ){
	    	$QH_or_HQ_hash=BlastHandle::ChangeBlastResult_into_Hit_Qry_HASH( $inBlastRstArrayHash );  $Qry_Hit_Type='Hit_Qry';
	    }
	    
	                                                                                                                                                                                                                  #DirFileHandle::PrintAndWarnDumper ($QH_or_HQ_hash, "\$QH_or_HQ_hash=$QH_or_HQ_hash");
	                                                                                                                                                                                                                 #print "\n20190124-0-0-1-0\n";
	    if (   (  defined ( $QH_or_HQ_hash )  ) && (  ref ( $QH_or_HQ_hash ) eq 'HASH'  )   ){                                                                                                       #print "\n20190124-0-0-1-1 \$QH_or_HQ_hash=$QH_or_HQ_hash\n";
	    	foreach my $level_1_key (    sort { $a cmp $b } (   keys (  %{ $QH_or_HQ_hash }  )   )    ){                                                                                              #print "\n20190124-0-0-1-2 \$level_1_key=$level_1_key\n";
	  	    if (   (  defined ( $QH_or_HQ_hash->{$level_1_key} )  ) && (  ref ( $QH_or_HQ_hash->{$level_1_key} ) eq 'HASH'  )   ){	                                                                 #print "\n20190124-0-0-1-2-1 \$level_1_key=$level_1_key\n";
	  	    	
	  	    	my @Sorted_keyArray;                                                                                                                                                                   #=keys (  %{ $QH_or_HQ_hash->{$level_1_key} }   );	  	my $tp=ArrayHashChange::ChangArrayIntoSimpleString (\@Sorted_keyArray); print "\n20190124-0-0-1-3 \@Sorted_keyArray=$tp\n";    	 
	  	    	if (  $sortByWhat eq '_hitEvalue'  ){   @Sorted_keyArray= sort { $QH_or_HQ_hash->{$level_1_key}->{$a}->{$sortByWhat} <=> $QH_or_HQ_hash->{$level_1_key}->{$b}->{$sortByWhat} } (   keys (  %{ $QH_or_HQ_hash->{$level_1_key} }  )   );      }
	  	    	else                                {   @Sorted_keyArray= sort { $QH_or_HQ_hash->{$level_1_key}->{$b}->{$sortByWhat} <=> $QH_or_HQ_hash->{$level_1_key}->{$a}->{$sortByWhat} } (   keys (  %{ $QH_or_HQ_hash->{$level_1_key} }  )   );      }
	  	                                                                                                                                                                                                                    #$tp=ArrayHashChange::ChangArrayIntoSimpleString (\@Sorted_keyArray); print "\n20190124-0-0-1-4 \@Sorted_keyArray=$tp\n";  
	  	      
	  	      
	  	      my $sort_order_idx=0;  #print "\n20190124-0-0-1-5-0 \n";  
	  	      FECMK: foreach my $level_2_key ( @Sorted_keyArray ){
	  	      	if (      (  defined ( $QH_or_HQ_hash->{$level_1_key}->{$level_2_key} )                 ) && (  ref ( $QH_or_HQ_hash->{$level_1_key}->{$level_2_key}                ) eq 'HASH'  )    
	  	      	       && (  defined ( $QH_or_HQ_hash->{$level_1_key}->{$level_2_key}->{$sortByWhat} )  ) && (        $QH_or_HQ_hash->{$level_1_key}->{$level_2_key}->{$sortByWhat} =~m /\S+/    )
	  	      	   )
	  	      	{
	  	      		
	  	      		
	  	      	  
	  	      		my $filterOut_or_not=0;                                                                                                                                                                             #print "\n20190124-0-0-1-5 \$level_2_key=$level_2_key\n";  
	  	      		
	  	      		my $cutVal=$QH_or_HQ_hash->{$level_1_key}->{$level_2_key}->{$sortByWhat};                                                                                                        #print "\n20190124-0-0-1-6 \$sortCutoff=$sortCutoff \$cutVal=$cutVal=\$QH_or_HQ_hash->{$level_1_key}->{$level_2_key}->{$sortByWhat}=$QH_or_HQ_hash->{$level_1_key}->{$level_2_key}->{$sortByWhat}\n";  
	  	      	  if (   (  defined ( $sortCutoff )  ) && ( $sortCutoff=~m/\d+/ )   ){
	  	      	  	if (  $sortByWhat eq '_hitEvalue'  )  {    if ( $cutVal > $sortCutoff) { $filterOut_or_not=1; }     }
	  	      	  	else                                  {    if ( $cutVal < $sortCutoff) { $filterOut_or_not=1; }     }	
	  	      	  }                                                                                                                                                                                                              #print "\n20190124-0-0-1-7 \$sortCutoff=$sortCutoff \$cutVal=$cutVal \$filterOut_or_not=$filterOut_or_not\n";  
	  	      	  if ( $filterOut_or_not==0 ){
	  	      	  	$sorted_cutted_HASH->{$level_1_key}->{$level_2_key}=Storable::dclone( $QH_or_HQ_hash->{$level_1_key}->{$level_2_key} );
	  	      	    
	  	      	    $sorted_cutted_HASH->{$level_1_key}->{$level_2_key}->{'0_0_0_2_0___Qry_Hit_type'}=$Qry_Hit_Type;    
  		  	  	    $sorted_cutted_HASH->{$level_1_key}->{$level_2_key}->{'0_0_0_2_1___sort_by_what'}=$sortByWhat;
  		  	  	    $sorted_cutted_HASH->{$level_1_key}->{$level_2_key}->{'0_0_0_2_2___sort__cutoff'}=$sortCutoff;
  		  	  	    $sorted_cutted_HASH->{$level_1_key}->{$level_2_key}->{'0_0_0_2_3_best_numberLmt'}=$best_number if (   (  defined ( $best_number )  ) && ( $best_number=~m/\d+/ ) && ( $best_number > 0 )   );
  		  	  	    $sorted_cutted_HASH->{$level_1_key}->{$level_2_key}->{'0_0_0_2_4_sort_order_idx'}=$sort_order_idx;
  		  	  	    
  		  	  	    $sorted_cutted_HASH->{$level_1_key}->{$level_2_key}->{'0_0_0_3_0_copy_________hitEvalue'}=$QH_or_HQ_hash->{$level_1_key}->{$level_2_key}->{'_hitEvalue'        } if (   (  defined ( $QH_or_HQ_hash->{$level_1_key}->{$level_2_key}->{'_hitEvalue'        } )  ) && (        $QH_or_HQ_hash->{$level_1_key}->{$level_2_key}->{'_hitEvalue'        } =~m /\S+/ )   );
  		  	  	    $sorted_cutted_HASH->{$level_1_key}->{$level_2_key}->{'0_0_0_3_0_copy_hitTotalConserved'}=$QH_or_HQ_hash->{$level_1_key}->{$level_2_key}->{'_hitTotalConserved'} if (   (  defined ( $QH_or_HQ_hash->{$level_1_key}->{$level_2_key}->{'_hitTotalConserved'} )  ) && (        $QH_or_HQ_hash->{$level_1_key}->{$level_2_key}->{'_hitTotalConserved'} =~m /\S+/ )   );
  		  	  	    $sorted_cutted_HASH->{$level_1_key}->{$level_2_key}->{'0_0_0_3_0_copy______hitCoverRate'}=$QH_or_HQ_hash->{$level_1_key}->{$level_2_key}->{'hitCoverRate'      } if (   (  defined ( $QH_or_HQ_hash->{$level_1_key}->{$level_2_key}->{'hitCoverRate'      } )  ) && (        $QH_or_HQ_hash->{$level_1_key}->{$level_2_key}->{'hitCoverRate'      } =~m /\S+/ )   );
  		  	  	    $sorted_cutted_HASH->{$level_1_key}->{$level_2_key}->{'0_0_0_3_0_copy__otal_qry_covrAte'}=$QH_or_HQ_hash->{$level_1_key}->{$level_2_key}->{'total_qry_covrAte' } if (   (  defined ( $QH_or_HQ_hash->{$level_1_key}->{$level_2_key}->{'total_qry_covrAte' } )  ) && (        $QH_or_HQ_hash->{$level_1_key}->{$level_2_key}->{'total_qry_covrAte' } =~m /\S+/ )   );
  		  	  	    $sorted_cutted_HASH->{$level_1_key}->{$level_2_key}->{'0_0_0_3_0_copy____queryCoverRate'}=$QH_or_HQ_hash->{$level_1_key}->{$level_2_key}->{'queryCoverRate'    } if (   (  defined ( $QH_or_HQ_hash->{$level_1_key}->{$level_2_key}->{'queryCoverRate'    } )  ) && (        $QH_or_HQ_hash->{$level_1_key}->{$level_2_key}->{'queryCoverRate'    } =~m /\S+/ )   );
  		  	  	    
  		  	  	                                                                                                                                                                                                                  #print "\n20190124-0-0-1-8-0 $sort_order_idx \$QH_or_HQ_hash->{$level_1_key}->{$level_2_key}->{'_hitEvalue'        }=$QH_or_HQ_hash->{$level_1_key}->{$level_2_key}->{'_hitEvalue'        }\n";   
  		  	  	                                                                                                                                                                                                                  #print "\n20190124-0-0-1-8-1 $sort_order_idx \$QH_or_HQ_hash->{$level_1_key}->{$level_2_key}->{'_hitTotalConserved'        }=$QH_or_HQ_hash->{$level_1_key}->{$level_2_key}->{'_hitTotalConserved'        }\n";   
  		  	  	                                                                                                                                                                                                                  #DirFileHandle::PrintAndWarnDumper ($sorted_cutted_HASH->{$level_1_key}->{$level_2_key}, "\$sorted_cutted_HASH->{$level_1_key}->{$level_2_key}=$sorted_cutted_HASH->{$level_1_key}->{$level_2_key}");
  		  	  	    
	  	      	    $sort_order_idx++;
	  	      	    
	  	      	    if (   (  defined ( $best_number )  ) && ( $best_number=~m/\d+/ ) && ( $best_number > 0 ) && ( $sort_order_idx >= $best_number )   ){
	  	      	    	last FECMK;
	  	      	    }
	  	      	    
	  	      	  }
	  	      	  
	  	      	}
	  	      }
	  	    }
	  	    
	  	  }
	    }
	    
	  }
	}
	
	return $sorted_cutted_HASH;
}
  



sub findingIdexOfAword{   #在一个字符串中 找到子字符串的 出现位置及次数
  my ($inWord, $inString, $capsSensitive)=@_;
  
  if (  ( defined ($capsSensitive) ) && ( $capsSensitive == 1 )  ){ 	                                              }  #大小敏感，则按照输入的 字符 和 字符串的类型 来进行idx的检索
  else                                                            { $inWord=uc ($inWord); $inString=uc ($inString); }  #大小不敏感，通常是DNA或protein的序列，则将输入的 字符 和 字符串统一转化为大写，再来进行idx的检索
  
  my $outArray;;
  
  my $startIDX=0;
  while ($startIDX >=0){
  	$startIDX=index ($inString, $inWord, $startIDX);  #print "$inWord, $inString\t\$startIDX=$startIDX\n";
  	if ($startIDX>=0){   push (  @{ $outArray }, ( $startIDX+1 )  );   $startIDX++; }
  }
  
  return $outArray;

}                                    




sub GetSubString_from_nbHASH{  #my $outString=BlastHandle::GetSubString_from_nbHASH ($fastaCMDoutHASH, $stt, $end, $ZorF);
	my ($fastaCMDoutHASH, $stt, $end, $ZorF)=@_;
	
	my $warnMsgBody="\nIn package  BlastHandle,\tIn sub GetSubString_from_nbHASH,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $subCalInfom=DirFileHandle::print_SubCallerInform;
  
  my @nbArray;
  if    (   (  defined ( $ZorF )  ) && ( $ZorF =~m /^\S+$/ ) && (  defined ( $stt )  ) && ( $stt =~m /^\d+$/ ) && (  defined ( $end )  ) && ( $end =~m /^\d+$/ )    ){
  	if    ( $ZorF eq '+' ) {  		
  		if ( $stt > $end ) { DieWork::Just_dieWork( $die_MsgHead."\n\$ZorF=$ZorF \$stt=$stt should <= $end=\$end !!!: $!".$subCalInfom ); }
  		@nbArray=($stt .. $end);
  	}
  	elsif ( $ZorF eq '-' ) {
  		if ( $stt < $end ) { DieWork::Just_dieWork( $die_MsgHead."\n\$ZorF=$ZorF \$stt=$stt should >= $end=\$end !!!: $!".$subCalInfom ); }
  		@nbArray=($end .. $stt); @nbArray=reverse @nbArray;
  	}
  	else{
  	  DieWork::Just_dieWork( $die_MsgHead."\n\$ZorF=$ZorF should be defined and equal to + or - !!!: $!".$subCalInfom );
    }
  }
  else{
  	DieWork::Just_dieWork( $die_MsgHead."\n\$ZorF=$ZorF \$stt=$stt \$end=$end   should be defined not as NULL !!!: $!".$subCalInfom );
  }
  
  if  (   (  defined ( $fastaCMDoutHASH )  ) && (  ref ( $fastaCMDoutHASH ) eq 'HASH'  )   ){  
  }
  else{
  	DieWork::Just_dieWork( $die_MsgHead."\n\$fastaCMDoutHASH=$fastaCMDoutHASH should be defined and as a HASH ref !!!: $!".$subCalInfom );
  }
  
  my $outString;
  foreach my $eachNB ( @nbArray ){
  	$outString.=$fastaCMDoutHASH->{$ZorF}->{$eachNB};  #print "20190406-2-2-0 " , TimeWork::GetNowTimePid_microSecond (),"\$outString=$outString \$fastaCMDoutHASH->{$ZorF}->{$eachNB}=$fastaCMDoutHASH->{$ZorF}->{$eachNB}\n";
  }
	return $outString;
}
  				
sub Build_Nb_2_word_HASH_from_fastacmd{ #my $outHash=BlastHandle::Build_Nb_2_word_HASH_from_fastacmd ($ctgNM, $dbNM, $sttNb, $endNb);
	my ($ctgNM, $dbNM, $sttNb, $endNb)=@_;
	
	my $warnMsgBody="\nIn package  BlastHandle,\tIn sub Build_Nb_2_word_HASH_from_fastacmd,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $subCalInfom=DirFileHandle::print_SubCallerInform;
  
  my $fastacmdPath="fastacmd";
  my $localPostion;
  if (   (  defined ( $sttNb )  ) || (  defined ( $endNb )  )   ){
  	if (   (  defined ( $sttNb )  ) && ( $sttNb =~m /^\d+$/ ) && (  defined ( $endNb )  ) && ( $endNb =~m /^\d+$/ ) && ( 1 <= $sttNb ) && ( $sttNb <= $endNb )   ){
  		$localPostion="-L $sttNb,$endNb";
  	}
  	else {
  		DieWork::Just_dieWork( $die_MsgHead."\n\$sttNb=$sttNb or  \$endNb=$endNb is not right!!!: $!".$subCalInfom );
  	}
  }
  else{
  	$sttNb=1;
  }
  
  my $fastacmdCMD="$fastacmdPath -s \"$ctgNM\"  -d $dbNM $localPostion"; warn "\n20190926 \$fastacmdCMD=$fastacmdCMD\n\n";
	my $fastacmdSEQ=`$fastacmdCMD`;
	my $seq=FastaFileHandle::GetSequence($fastacmdSEQ); my $seqLength=length $seq; if ( $seqLength < 1 )  { DieWork::Just_dieWork( $die_MsgHead."\n\$ctgNM=$ctgNM or  \$dbNM=$dbNM or $sttNb=$sttNb or  \$endNb=$endNb is not right!!!: $!\n\$seqLength=$seqLength \$seq=$seq\n\n".$subCalInfom ); }
	my $OutHash=BlastHandle::build_NBtoWord_HASH_for_String_Start_stt_ZF ($seq, $sttNb); 
	return $OutHash;
}


sub Build_Nb_2_word_HASH_from_bioPerlInDex{ #my $outHash=BlastHandle::Build_Nb_2_word_HASH_from_bioPerlInDex ($ctgNM, $dbNM, $sttNb, $endNb);
	my ($ctgNM, $idxFile, $sttNb, $endNb)=@_;
	
	my $warnMsgBody="\nIn package  BlastHandle,\tIn sub Build_Nb_2_word_HASH_from_bioPerlInDex,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $subCalInfom=DirFileHandle::print_SubCallerInform;
  
  my $seq;  
  if (   (  defined ( $sttNb )  ) || (  defined ( $endNb )  )   ){
  	if (   (  defined ( $sttNb )  ) && ( $sttNb =~m /^\d+$/ ) && (  defined ( $endNb )  ) && ( $endNb =~m /^\d+$/ ) && ( 1 <= $sttNb ) && ( $sttNb <= $endNb )   ){  		
  	  $seq= FastaFileHandle::feach_seqString_from_idx_File($ctgNM, $idxFile, $sttNb, $endNb);
  	}
  	else {
  		DieWork::Just_dieWork( $die_MsgHead."\n\$sttNb=$sttNb or  \$endNb=$endNb is not right!!!: $!".$subCalInfom );
  	}
  }
  else{
  	$sttNb=1;
  	$seq= FastaFileHandle::feach_seqString_from_idx_File($ctgNM, $idxFile);
  }
  
  my $seqLength=length $seq; 
  if ( $seqLength < 1 )  { DieWork::Just_dieWork( $die_MsgHead."\n\$ctgNM=$ctgNM or  \$idxFile=$idxFile or $sttNb=$sttNb or  \$endNb=$endNb is not right!!!: $!\n\$seqLength=$seqLength \$seq=$seq\n\n".$subCalInfom ); }
	my $OutHash=BlastHandle::build_NBtoWord_HASH_for_String_Start_stt_ZF ($seq, $sttNb); 
	return $OutHash;
}

sub build_NBtoWord_HASH_for_String_Start_stt_ZF{  #my $outHash=BlastHandle::build_NBtoWord_HASH_for_String_Start_stt_ZF ($inString, $stt);  
	#为输入的序列，建立 位置为Key(第一个位置标为1或$stt)，该位置字符为值的 hash， ACGTGTATGTCATCGTATGACTCATGACTGACTGCTAGCTACTAGCTAGCTAGCTAGCTAGCTGATCGATCGATGCATGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTGAGCTAGCTAGCATGCTAGCTAGCTAGCATGCTAGCTAGCATCGTACGAGCT
  my ($inString, $stt)=@_;
  $inString=~s/\s//g; $inString=~s/\n//g;
  my @instr=split '',$inString;
  
  my $outHash;
  for (my $i=0; $i< @instr; $i++){
    $outHash->{'+'}->{$i+$stt}= $instr[$i];  #warn "\$outHash->{\$i+1}=\$outHash->{$i+1}=\$instr[$i]=$instr[$i]\n";
    $outHash->{'-'}->{$i+$stt}= FastaFileHandle::ReverseComplementString ($instr[$i]); 
  }
  return $outHash;
}


sub build_NBtoWord_HASH_for_String_Start_stt{  #为输入的序列，建立 位置为Key(第一个位置标为1或$stt)，该位置字符为值的 hash， ACGTGTATGTCATCGTATGACTCATGACTGACTGCTAGCTACTAGCTAGCTAGCTAGCTAGCTGATCGATCGATGCATGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTGAGCTAGCTAGCATGCTAGCTAGCTAGCATGCTAGCTAGCATCGTACGAGCT
  my ($inString, $stt)=@_;
  $inString=~s/\s//g; $inString=~s/\n//g;
  my @instr=split '',$inString;
  
  my $outHash;
  for (my $i=0; $i< @instr; $i++){
    $outHash->{$i+$stt}=$instr[$i];  #warn "\$outHash->{\$i+1}=\$outHash->{$i+1}=\$instr[$i]=$instr[$i]\n";
  }
  return $outHash;
}

sub build_NBtoWord_HASH_for_String{  #为比对的序列中的每个字符，建立 位置为Key，该位置字符为值的 hash，特别是 这种序列 YASTACTVSSCSPLPATTVSSKHGPLSSSRVLKASGRTECTSCCWSRNTAITRATYRRVL-----PDLVPGPSLHHWK*LPRSRMNITVDMSFVTIKGILPT
  my ($inString)=@_;
  $inString=~s/\s//g; $inString=~s/\n//g;
  my @instr=split '',$inString;
  
  my $outHash;
  for (my $i=0; $i< @instr; $i++){
    $outHash->{$i+1}=$instr[$i];  #warn "\$outHash->{\$i+1}=\$outHash->{$i+1}=\$instr[$i]=$instr[$i]\n";
  }
  return $outHash;
}

sub build_NBtoWord_HASH_for_String_at_OrgSEQUENCE{  #为比对的序列中的每个字符找回其在原序列中的位置
  my ($inHASH, $strand, $startNB,  $end__NB, $stepNB, $dbNM, $ctgNM)=@_;
  my $strandNB=1;  if ($strand eq '-') {$strandNB=-1;}
  
  my $smNb_of_stt_end=SeqSegmentsTools::getSmallOne($startNB, $end__NB);
  my $bgNb_of_stt_end=SeqSegmentsTools::getbigOne  ($startNB, $end__NB);
  
  my $lastNB;
  if ( $strandNB == 1 )                   {  $lastNB=$smNb_of_stt_end-1; }
  else                                    {  $lastNB=$bgNb_of_stt_end+1; }
  
  
  my $extracFastaNb_2_wd_HASH;
  if ($stepNB==3){ 
  	$extracFastaNb_2_wd_HASH=BlastHandle::Build_Nb_2_word_HASH_from_bioPerlInDex ($ctgNM, $dbNM, $smNb_of_stt_end, $bgNb_of_stt_end);
  	#$extracFastaNb_2_wd_HASH=BlastHandle::Build_Nb_2_word_HASH_from_fastacmd     ($ctgNM, $dbNM, $smNb_of_stt_end, $bgNb_of_stt_end);
  }
  
  my $outHASH;
  my $reverseOUThash;
  my $DNAseq;
  if ( ref($inHASH) eq 'HASH' ){
  	my $foecStp=0;
  	foreach my $alnPos (   sort {$a<=>$b} keys (  %{ $inHASH }  )   ) {
  		my $realPos;
  		if ($inHASH->{$alnPos} eq '-'){
  			$outHASH->{$alnPos}->{'Type'}='ITSaGAP';
  			$realPos=$lastNB;
  			
  		}
  		else {
  			$outHASH->{$alnPos}->{'Type'}='NOTaGap';
  			if ($stepNB==3){
  				my $smNB=$outHASH->{$alnPos}->{'real_1_3_Pos'}->{1}=        $lastNB+$strandNB*($stepNB-2);                         $reverseOUThash->{ $smNB }->{'alnPos'}=$alnPos;  $reverseOUThash->{ $smNB }->{'PePChar'}=$inHASH->{$alnPos};   
  				my $mdNB=$outHASH->{$alnPos}->{'real_1_3_Pos'}->{2}=        $lastNB+$strandNB*($stepNB-1);                         $reverseOUThash->{ $mdNB }->{'alnPos'}=$alnPos;  $reverseOUThash->{ $mdNB }->{'PePChar'}=$inHASH->{$alnPos};   
  				my $bgNB=$outHASH->{$alnPos}->{'real_1_3_Pos'}->{3}=$lastNB=$lastNB+$strandNB*$stepNB;                             $reverseOUThash->{ $bgNB }->{'alnPos'}=$alnPos;  $reverseOUThash->{ $bgNB }->{'PePChar'}=$inHASH->{$alnPos};   
  				
  				#for Test   $dbNM, $ctgNM
  				#my $ParmL; my $ParS; 
  				#if ($strandNB ==1){$ParmL=$smNB.",".$bgNB; $ParS=1; } 
  				#else              {$ParmL=$bgNB.",".$smNB; $ParS=2; } 
  				#my $fastacmdCMD="fastacmd -s \"$ctgNM\" -d $dbNM -L $ParmL -S $ParS";
  				#my $fastacmdSEQ=`$fastacmdCMD`; 
  				#my $seq=FastaFileHandle::GetSequence($fastacmdSEQ);
  				
  				my $threeWd_stt=$smNB; my $threeWd_end=$bgNB; my $threeWd_ZorF;
  				if ($strandNB ==1){  $threeWd_ZorF='+'; }   #$threeWd_stt=$smNB; $threeWd_end=$bgNB;
  				else              {  $threeWd_ZorF='-'; }   #$threeWd_stt=$bgNB; $threeWd_end=$smNB;
  				my $seq=BlastHandle::GetSubString_from_nbHASH ($extracFastaNb_2_wd_HASH, $threeWd_stt, $threeWd_end, $threeWd_ZorF); if ( length ($seq) != 3 ) { print "20190406-2-1-0 " , TimeWork::GetNowTimePid_microSecond ()," \$seq=$seq ,", DirFileHandle::PrintAndWarnDumper ($extracFastaNb_2_wd_HASH, "20190406-2-1-1") ,"\$extracFastaNb_2_wd_HASH=$extracFastaNb_2_wd_HASH, \$threeWd_stt=$threeWd_stt, \$threeWd_end=$threeWd_end, \$threeWd_ZorF=$threeWd_ZorF\n";}
  				
  				my $pep=FastaFileHandle::TranslateDNAtoPep($seq); 
  				$outHASH->{$alnPos}->{'3DNAchar'}=$seq;  $DNAseq.=$seq;
  				$outHASH->{$alnPos}->{'TraslaAA'}=$pep;  #print "\$seq=$seq\t\$pep=$pep\t$inHASH->{$alnPos}=\$inHASH->{$alnPos}\$outHASH->{$alnPos}->{'3DNAchar'}=$outHASH->{$alnPos}->{'3DNAchar'}\$outHASH->{$alnPos}->{'TraslaAA'}=$outHASH->{$alnPos}->{'TraslaAA'}\n";
  				if    ( $pep eq $inHASH->{$alnPos} ){}
  				elsif ( 
  				            ($pep eq 'X')
  				        &&  (
  				              ( $inHASH->{$alnPos} eq 'U' ) || ( $inHASH->{$alnPos} eq 'O' ) || ( $inHASH->{$alnPos} eq 'B' ) || ( $inHASH->{$alnPos} eq 'Z' ) || ( $inHASH->{$alnPos} eq 'J' ) || ( $inHASH->{$alnPos} eq 'X' ) || ( $inHASH->{$alnPos} eq 'C' )
  				            )
  				      ){}
  			  elsif ( 
  				            ($pep eq '*')
  				        &&  (
  				              ( $inHASH->{$alnPos} eq 'U' ) || ( $inHASH->{$alnPos} eq 'O' ) || ( $inHASH->{$alnPos} eq 'B' ) || ( $inHASH->{$alnPos} eq 'Z' ) || ( $inHASH->{$alnPos} eq 'J' ) || ( $inHASH->{$alnPos} eq 'X' ) || ( $inHASH->{$alnPos} eq 'C' )
  				            )
  				      ){}
  				else {
  					$outHASH->{$alnPos}->{'wrongFound'}="\$pep=$pep ne $inHASH->{$alnPos}=\$inHASH->{$alnPos}"; 
  					#$outHASH->{$alnPos}->{'fastacmdCMD'}=$fastacmdCMD; 
  					#$outHASH->{$alnPos}->{'fastacmdSEQ'}=$fastacmdSEQ; 
  					#die "\n\n\nDIE!!!!\nin package BlastHandle\n\n\nin Sub build_NBtoWord_HASH_for_String_at_OrgSEQUENCE\n\n  \$inHASH=$inHASH, \$strand=$strand, \$startNB=$startNB,  \$end__NB=$end__NB, \$stepNB=$stepNB, \$dbNM=$dbNM, \$ctgNM=$ctgNM \$strandNB=$strandNB \$threeWd_stt=$threeWd_stt, \$threeWd_end=$threeWd_end, $threeWd_ZorF=$threeWd_ZorF\t\$seq=$seq\t\$pep=$pep\t$inHASH->{$alnPos}=\$inHASH->{$alnPos} \n\n\n\n\n";
  					print "\n\n\nDIE!!!!\nin package BlastHandle\n\n\nin Sub build_NBtoWord_HASH_for_String_at_OrgSEQUENCE\n\n  \$inHASH=$inHASH, \$strand=$strand, \$startNB=$startNB,  \$end__NB=$end__NB, \$stepNB=$stepNB, \$dbNM=$dbNM, \$ctgNM=$ctgNM \$strandNB=$strandNB \$threeWd_stt=$threeWd_stt, \$threeWd_end=$threeWd_end, $threeWd_ZorF=$threeWd_ZorF\t\$seq=$seq\t\$pep=$pep\t$inHASH->{$alnPos}=\$inHASH->{$alnPos} \n\n\n\n\n";
  				}
  				
  				$reverseOUThash->{ $smNB }->{'DNAChar'}=substr ($seq,0,1);  
          $reverseOUThash->{ $mdNB }->{'DNAChar'}=substr ($seq,1,1); 
          $reverseOUThash->{ $bgNB }->{'DNAChar'}=substr ($seq,2,1); 

  				
  			}
  			else {
  				$outHASH->{$alnPos}->{'real_1_1_Pos'}     =$lastNB=$lastNB+$strandNB*$stepNB;
  				$reverseOUThash->{ $lastNB+$strandNB*$stepNB }->{'alnPos'}=$alnPos;
  				$reverseOUThash->{ $lastNB+$strandNB*$stepNB }->{'Char'}=$inHASH->{$alnPos};
  			}
  			
  		}
  		
  		
  		$outHASH->{$alnPos}->{'AlnChar'}=$inHASH->{$alnPos};
  		
  		
  		$foecStp++;
  	}
  }
  return [$outHASH, $reverseOUThash, $DNAseq ];
  
}



sub log_10{   #   BlastHandle::log_10 
	my ($inNumber)=@_;
	
	my $warnMsgBody="\nIn package  BlastHandle,\tIn sub log_10,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	
	my $outNumber=log($inNumber)/log(10);	
	
}

sub makeKeyWordHash{ #   BlastHandle::makeKeyWordHash 
	my ($inArray)=@_;
	
	my $warnMsgBody="\nIn package  BlastHandle,\tIn sub makeKeyWordHash,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	
	my $outHash;
	if (   (  defined ( $inArray )  ) && ( ref ($inArray) eq 'ARRAY' )   ){
		my $totalNb=@{ $inArray };
		my $weiShu=&log_10($totalNb); $weiShu=int ($weiShu+1);
		for (  my $i=0;  $i<@{ $inArray }; $i++  ){
		my $fmt_i=sprintf ("%0${weiShu}d", $i);  
		  $outHash->{ $inArray->[$i] }=$fmt_i.'_'.$inArray->[$i];
	  }
	}
	else{
		my $dieMsg="$die_MsgHead\n\n\$inArray=$inArray should be a ARRAY ref!!!!\n\n\n";
		print $dieMsg; die $dieMsg;
	}
	#print "20180920\n";
	#DirFileHandle::PrintAndWarnDumper ($outHash);
	return $outHash;
}

sub DrawBlast2seq_20180920{    #   BlastHandle::DrawBlast2seq_20180920  #pmAble#   #绘制Blast2seq的图像
	my ($file, $outPNGfile, $bl2seq_1_fastaName, $bl2seq_2_fastaName, $query_realStart, $query_realEnd)=@_;   #warn "In sub &DrawBlast2seq_20180920 (\$file=$file, \$outPNGfile=$outPNGfile)\n\n";
	
	
	
	#sub Map_PepSeq_BackTo_DNAseq{ #Input 1, protein sequence string, 2, DNA sequence string, 3, a directory path to hold middle step information. TARGET: map the protein sequence back to dna sequence, find the coding region and in frame tga positions 
  my ($inString_1, $inString_2, $workingDir)=@_;
  
 
  
  my $warnMsgBody="\nIn package  BlastHandle,\tIn sub DrawBlast2seq_20180920,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	
	my $praseHash=&BioPerlBlastPraser_bl2seq_20170330($file);  
	
	my @keyWordArray=(

	  'hit___evalue'  ,
    'tot_ConsvPct'  ,
	  
	  'qrySeqLength'  ,
	  'query__start'  ,
	  'query____end'  ,
	  'qry_cov_leng'  ,
	  'qryCutCovLen'  ,
	  'qryCtCv_CvRT'  ,
	  'qrCCv_CvRtPt'  ,
	  'qry_WlCsvPct'  ,
	              
	              
	              
	              
	  'hitSeqLength'  ,
	  'hit____start'  ,
	  'hit______end'  ,
	  'hit_cov_leng'  ,
	  'hitCutCovLen'  ,
	  'hitCtCv_CvRT'  ,
	  'htCCv_CvRtPt'  ,
	  'hit_tot_cosv'  ,
	  'hit_WlCsvPct'  ,
	  
	  
              
	  'blt_out_HASH'  ,
	  
	);
	
	my $keyWordHash=&makeKeyWordHash( [@keyWordArray] );
		
		
	
	my $outHashHere;
	#$result->query_length
	
	                       $outHashHere->{ $keyWordHash->{'blt_out_HASH'} }=$praseHash;
	my $hitTotalConserved =$outHashHere->{ $keyWordHash->{'hit_tot_cosv'} }=$praseHash->{'_ResultArray'}->[ 0 ]->{'_hitArray'}->[ 0 ]->{'_hitTotalConserved'  };       
	my $hitEvalue         =$outHashHere->{ $keyWordHash->{'hit___evalue'} }=$praseHash->{'_ResultArray'}->[ 0 ]->{'_hitArray'}->[ 0 ]->{'_hitEvalue'          };               
	my $hitSeqLength      =$outHashHere->{ $keyWordHash->{'hitSeqLength'} }=$praseHash->{'_ResultArray'}->[ 0 ]->{'_hitArray'}->[ 0 ]->{'_hitSeqLength'       };
	my $hit_queryStart    =$outHashHere->{ $keyWordHash->{'query__start'} }=$praseHash->{'_ResultArray'}->[ 0 ]->{'_hitArray'}->[ 0 ]->{'hit_queryStart'      };           
  my $hit_queryEnd      =$outHashHere->{ $keyWordHash->{'query____end'} }=$praseHash->{'_ResultArray'}->[ 0 ]->{'_hitArray'}->[ 0 ]->{'hit_queryEnd'        };             
  my $hit_queryCover    =$outHashHere->{ $keyWordHash->{'qry_cov_leng'} }=$praseHash->{'_ResultArray'}->[0  ]->{'_hitArray'}->[0  ]->{'hit_queryRangeLength'}; 
  my $queryLength       =$outHashHere->{ $keyWordHash->{'qrySeqLength'} }=$praseHash->{'_ResultArray'}->[ 0 ]->{'_query_length'};   #$result->query_length;  
	my $queryCoverLength  =$outHashHere->{ $keyWordHash->{'qryCutCovLen'} }=$praseHash->{'_ResultArray'}->[ 0 ]->{'_hitArray'}->[ 0 ]->{'queryCoverLength'    };         
  my $queryCoverRate    =$outHashHere->{ $keyWordHash->{'qryCtCv_CvRT'} }=$praseHash->{'_ResultArray'}->[0  ]->{'_hitArray'}->[0  ]->{'queryCoverRate'      };        
  my $qureyCovRatePct   =$outHashHere->{ $keyWordHash->{'qrCCv_CvRtPt'} }=$praseHash->{'_ResultArray'}->[0  ]->{'_hitArray'}->[0  ]->{'qureyCovRatePct'     }; 
  my $hit_HitStart      =$outHashHere->{ $keyWordHash->{'hit____start'} }=$praseHash->{'_ResultArray'}->[ 0 ]->{'_hitArray'}->[ 0 ]->{'hit_HitStart'        };             
  my $hit_HitEnd        =$outHashHere->{ $keyWordHash->{'hit______end'} }=$praseHash->{'_ResultArray'}->[ 0 ]->{'_hitArray'}->[ 0 ]->{'hit_HitEnd'          };               
  my $hit_HitCover      =$outHashHere->{ $keyWordHash->{'hit_cov_leng'} }=$praseHash->{'_ResultArray'}->[0  ]->{'_hitArray'}->[0  ]->{'hit_HitRangeLength'  };     
  my $hitCoverLength    =$outHashHere->{ $keyWordHash->{'hitCutCovLen'} }=$praseHash->{'_ResultArray'}->[ 0 ]->{'_hitArray'}->[ 0 ]->{'hitCoverLength'      };           
  my $hitCoverRate      =$outHashHere->{ $keyWordHash->{'hitCtCv_CvRT'} }=$praseHash->{'_ResultArray'}->[0  ]->{'_hitArray'}->[0  ]->{'hitCoverRate'        };        
  my $hitCovRatePct     =$outHashHere->{ $keyWordHash->{'htCCv_CvRtPt'} }=$praseHash->{'_ResultArray'}->[0  ]->{'_hitArray'}->[0  ]->{'hitCovRatePct'       }; 
  my $totConsePercent   =$outHashHere->{ $keyWordHash->{'tot_ConsvPct'} }=$praseHash->{'_ResultArray'}->[0  ]->{'_hitArray'}->[0  ]->{'totConsePercent'     };          
  
  my $potintNb_of_totConsePercent=&Change100PercentNB_to_pointNB($totConsePercent);
  my $qryCsvPct         =$outHashHere->{ $keyWordHash->{'qry_WlCsvPct'} }=&ChangeTo100PercentNB (  $potintNb_of_totConsePercent*($queryCoverLength/$queryLength )/100  );
  my $hitCsvPct         =$outHashHere->{ $keyWordHash->{'hit_WlCsvPct'} }=&ChangeTo100PercentNB (  $potintNb_of_totConsePercent*($hitCoverLength  /$hitSeqLength)/100  );
 ##### 

 ########## 
  
	my $bestConservSumury="         Average Identity: $totConsePercent, Evalue: $hitEvalue";
	
	
	if (($query_realStart=~m/\d+/) && ($query_realStart>1)){  }else {$query_realStart=1;}                                          
	if (($query_realEnd=~m/\d+/)   && ($query_realEnd>1)  ){  }else {$query_realEnd=$queryLength;}   #这个设定也有可能是错误的,需谨慎    
	
	$bl2seq_1_fastaName="Query: $bl2seq_1_fastaName:(${hit_queryStart}-${hit_queryEnd})";
	$bl2seq_2_fastaName="Sbjct: $bl2seq_2_fastaName:(${hit_HitStart}-${hit_HitEnd})";
	
	
	my $Caverage_1="Query Cav: $qureyCovRatePct ($queryCoverLength/$hit_queryCover/$queryLength)";
	my $Caverage_2="Sbjct Cav: $hitCovRatePct ($hitCoverLength/$hit_HitCover/$hitSeqLength)";
	warn "\$bestConservSumury=$bestConservSumury\n";
	open (BLASTPNGOUT,">$outPNGfile") or die "cannot create \$outPNGfile=$outPNGfile : $!\n\n";

  my $searchio = Bio::SearchIO->new(-file   => $file,
                                    -format => 'blast') or die "parse failed";


  my $result = $searchio->next_result() or die "\n\nno result\n\n";

  my $panel = Bio::Graphics::Panel->new(-length    => $result->query_length,
                                        -width     => 800,
                                        -pad_left  => 120,
                                        -pad_right => 80,
                                        -key_style => 'left',
                                        -start     => 1 #$query_realStart#################
                                       );

  my $full_length = Bio::SeqFeature::Generic->new(-start      =>   1, #$query_realStart, #$hit_queryStart, #1, 
                                                  -end        =>   $result->query_length, #$query_realEnd,   #$hit_queryEnd,   #$result->query_length,
                                                  -seq_id     =>   $bestConservSumury,    #$result->query_name
                                                  -fontcolor  =>   'blue'
                                                  );
  $panel->add_track($full_length,
                    -glyph   => 'arrow',
                    -tick    => 2,
                    -fgcolor => 'black',
                    -double  => 1,
                    -label   => 1,
                    -key   => $bl2seq_1_fastaName
                   );

  my $track = $panel->add_track(-glyph       => 'segments', #'graded_segments',
                                -label       => 1,
                                -connector   => 'dashed',
                                -bgcolor     => 'red',
                                -fgcolor     => 'red',
                                -font2color  => 'black',
                                -sort_order  => 'high_score',
                                -key         => $bl2seq_2_fastaName , 
                                -description => sub {
                                                       my $feature = shift;  #print_all_sub_array ($feature);
                                                       return unless $feature->has_tag('description');  
                                                       my ($description) = $feature->each_tag_value('description');
                                                       my $score = $feature->score;
                                                       
                                                       #my $eValue= $feature->
                                                       #"score=$score $description";
                                                       #"score=$score";
                                                       #$bestConservSumury;
                                                       $Caverage_2;
                                                    }
                                );
  
  while( my $hit = $result->next_hit ) {
    #next unless $hit->significance < 1E-20;
    my $feature = Bio::SeqFeature::Generic->new(-score   => $hit->raw_score,
                                                                                                                               
                                                -seq_id  => $Caverage_1, #$hit->name,
                                                -tag     => {
                                                             description => $hit->description
                                                            },
                                               );
    while( my $hsp = $hit->next_hsp ) {
    	warn "1111  \$hsp->score=",$hsp->score,"\n\$hsp->bits=",$hsp->bits,"\n\$hsp->frac_conserved=",$hsp->frac_conserved,"\n\n";
    	my $cons=$hsp->frac_conserved;
    	$hsp->score($cons);
    	warn "22222  \$hsp->score=",$hsp->score,"\n\$hsp->bits=",$hsp->bits,"\n\$hsp->frac_conserved=",$hsp->frac_conserved,"\n\n";
    	#$hsp->bits=$hsp->frac_conserved;    warn "2222  \$hsp->bits=",$hsp->bits,"\n\$hsp->frac_conserved=",$hsp->frac_conserved,"\n\n";
      $feature->add_sub_SeqFeature($hsp,'EXPAND');
    }
  
    $track->add_feature($feature);
  }
  
  print BLASTPNGOUT $panel->png;

  close (BLASTPNGOUT);
	
	
	return $outHashHere ; 
}




sub DrawBlast2seq_20170330{     #pmAble#   #绘制Blast2seq的图像
	my ($file, $outPNGfile, $bl2seq_1_fastaName, $bl2seq_2_fastaName, $query_realStart, $query_realEnd)=@_;   warn "In sub &DrawBlast2seq_20170330 (\$file=$file, \$outPNGfile=$outPNGfile)\n\n";
	my $praseHash=&BioPerlBlastPraser_bl2seq_20170330($file);  
	
	my $outHashHere;
	#$result->query_length
	
	                       $outHashHere->{'99_blt_out_HASH'     }=$praseHash;
	my $hitTotalConserved =$outHashHere->{'_hitTotalConserved'  }=$praseHash->{'_ResultArray'}->[ 0 ]->{'_hitArray'}->[ 0 ]->{'_hitTotalConserved'  };       
	my $hitEvalue         =$outHashHere->{'_hitEvalue'          }=$praseHash->{'_ResultArray'}->[ 0 ]->{'_hitArray'}->[ 0 ]->{'_hitEvalue'          };               
	my $hitSeqLength      =$outHashHere->{'_hitSeqLength'       }=$praseHash->{'_ResultArray'}->[ 0 ]->{'_hitArray'}->[ 0 ]->{'_hitSeqLength'       };            
	                                                                                                                                       
	my $hit_queryStart    =$outHashHere->{'hit_queryStart'      }=$praseHash->{'_ResultArray'}->[ 0 ]->{'_hitArray'}->[ 0 ]->{'hit_queryStart'      };           
  my $hit_queryEnd      =$outHashHere->{'hit_queryEnd'        }=$praseHash->{'_ResultArray'}->[ 0 ]->{'_hitArray'}->[ 0 ]->{'hit_queryEnd'        };             
  my $hit_queryCover    =$outHashHere->{'hit_queryRangeLength'}=$praseHash->{'_ResultArray'}->[0  ]->{'_hitArray'}->[0  ]->{'hit_queryRangeLength'};     
  my $hit_HitStart      =$outHashHere->{'hit_HitStart'        }=$praseHash->{'_ResultArray'}->[ 0 ]->{'_hitArray'}->[ 0 ]->{'hit_HitStart'        };             
  my $hit_HitEnd        =$outHashHere->{'hit_HitEnd'          }=$praseHash->{'_ResultArray'}->[ 0 ]->{'_hitArray'}->[ 0 ]->{'hit_HitEnd'          };               
  my $hit_HitCover      =$outHashHere->{'hit_HitRangeLength'  }=$praseHash->{'_ResultArray'}->[0  ]->{'_hitArray'}->[0  ]->{'hit_HitRangeLength'  };     
                                                                                                                                       
  my $queryCoverLength  =$outHashHere->{'queryCoverLength'    }=$praseHash->{'_ResultArray'}->[ 0 ]->{'_hitArray'}->[ 0 ]->{'queryCoverLength'    };         
  my $queryCoverRate    =$outHashHere->{'queryCoverRate'      }=$praseHash->{'_ResultArray'}->[0  ]->{'_hitArray'}->[0  ]->{'queryCoverRate'      };        
  my $qureyCovRatePct   =$outHashHere->{'qureyCovRatePct'     }=$praseHash->{'_ResultArray'}->[0  ]->{'_hitArray'}->[0  ]->{'qureyCovRatePct'     };        
  my $hitCoverLength    =$outHashHere->{'hitCoverLength'      }=$praseHash->{'_ResultArray'}->[ 0 ]->{'_hitArray'}->[ 0 ]->{'hitCoverLength'      };           
  my $hitCoverRate      =$outHashHere->{'hitCoverRate'        }=$praseHash->{'_ResultArray'}->[0  ]->{'_hitArray'}->[0  ]->{'hitCoverRate'        };        
  my $hitCovRatePct     =$outHashHere->{'hitCovRatePct'       }=$praseHash->{'_ResultArray'}->[0  ]->{'_hitArray'}->[0  ]->{'hitCovRatePct'       };        
                                                                                                                                     
  my $totConsePercent   =$outHashHere->{'totConsePercent'     }=$praseHash->{'_ResultArray'}->[0  ]->{'_hitArray'}->[0  ]->{'totConsePercent'     };          
  
 ##### 

 ########## 
  
	my $bestConservSumury="         Average Identity: $totConsePercent, Evalue: $hitEvalue";
	
	my $queryLength      =$outHashHere->{'_query_length'        }=$praseHash->{'_ResultArray'}->[ 0 ]->{'_query_length'};   #$result->query_length;  
	
	if (($query_realStart=~m/\d+/) && ($query_realStart>1)){  }else {$query_realStart=1;}                                          
	if (($query_realEnd=~m/\d+/)   && ($query_realEnd>1)  ){  }else {$query_realEnd=$queryLength;}   #这个设定也有可能是错误的,需谨慎    
	
	$bl2seq_1_fastaName="Query: $bl2seq_1_fastaName:(${hit_queryStart}-${hit_queryEnd})";
	$bl2seq_2_fastaName="Sbjct: $bl2seq_2_fastaName:(${hit_HitStart}-${hit_HitEnd})";
	
	
	my $Caverage_1="Query Cav: $qureyCovRatePct ($queryCoverLength/$hit_queryCover/$queryLength)";
	my $Caverage_2="Sbjct Cav: $hitCovRatePct ($hitCoverLength/$hit_HitCover/$hitSeqLength)";
	warn "\$bestConservSumury=$bestConservSumury\n";
	open (BLASTPNGOUT,">$outPNGfile") or die "cannot create \$outPNGfile=$outPNGfile : $!\n\n";

  my $searchio = Bio::SearchIO->new(-file   => $file,
                                    -format => 'blast') or die "parse failed";


  my $result = $searchio->next_result() or die "\n\nno result\n\n";

  my $panel = Bio::Graphics::Panel->new(-length    => $result->query_length,
                                        -width     => 800,
                                        -pad_left  => 120,
                                        -pad_right => 80,
                                        -key_style => 'left',
                                        -start     => 1 #$query_realStart#################
                                       );

  my $full_length = Bio::SeqFeature::Generic->new(-start      =>   1, #$query_realStart, #$hit_queryStart, #1, 
                                                  -end        =>   $result->query_length, #$query_realEnd,   #$hit_queryEnd,   #$result->query_length,
                                                  -seq_id     =>   $bestConservSumury,    #$result->query_name
                                                  -fontcolor  =>   'blue'
                                                  );
  $panel->add_track($full_length,
                    -glyph   => 'arrow',
                    -tick    => 2,
                    -fgcolor => 'black',
                    -double  => 1,
                    -label   => 1,
                    -key   => $bl2seq_1_fastaName
                   );

  my $track = $panel->add_track(-glyph       => 'segments', #'graded_segments',
                                -label       => 1,
                                -connector   => 'dashed',
                                -bgcolor     => 'red',
                                -fgcolor     => 'red',
                                -font2color  => 'black',
                                -sort_order  => 'high_score',
                                -key         => $bl2seq_2_fastaName , 
                                -description => sub {
                                                       my $feature = shift;  #print_all_sub_array ($feature);
                                                       return unless $feature->has_tag('description');  
                                                       my ($description) = $feature->each_tag_value('description');
                                                       my $score = $feature->score;
                                                       
                                                       #my $eValue= $feature->
                                                       #"score=$score $description";
                                                       #"score=$score";
                                                       #$bestConservSumury;
                                                       $Caverage_2;
                                                    }
                                );
  
  while( my $hit = $result->next_hit ) {
    #next unless $hit->significance < 1E-20;
    my $feature = Bio::SeqFeature::Generic->new(-score   => $hit->raw_score,
                                                                                                                               
                                                -seq_id  => $Caverage_1, #$hit->name,
                                                -tag     => {
                                                             description => $hit->description
                                                            },
                                               );
    while( my $hsp = $hit->next_hsp ) {
    	warn "1111  \$hsp->score=",$hsp->score,"\n\$hsp->bits=",$hsp->bits,"\n\$hsp->frac_conserved=",$hsp->frac_conserved,"\n\n";
    	my $cons=$hsp->frac_conserved;
    	$hsp->score($cons);
    	warn "22222  \$hsp->score=",$hsp->score,"\n\$hsp->bits=",$hsp->bits,"\n\$hsp->frac_conserved=",$hsp->frac_conserved,"\n\n";
    	#$hsp->bits=$hsp->frac_conserved;    warn "2222  \$hsp->bits=",$hsp->bits,"\n\$hsp->frac_conserved=",$hsp->frac_conserved,"\n\n";
      $feature->add_sub_SeqFeature($hsp,'EXPAND');
    }
  
    $track->add_feature($feature);
  }
  
  print BLASTPNGOUT $panel->png;

  close (BLASTPNGOUT);
	
	
	return $outHashHere ; 
}


sub BioPerlBlastPraser_bl2seq_20170330{ #  #pmAble#    #解析blast 这个用于解析 bl2seq 的函数
  my ($blastfile, $outBestBlastFile, $numHits)=@_;   #my ($blastfile, $allSelProInformHASHinput, $numHits)=@_;    
  warn "Now in Sub &BioPerlpraseBlast!\nInput: \$blastfile=$blastfile\t\$numHits=$numHits\n"; print "Now in Sub &BioPerlpraseBlast!\nInput: \$blastfile=$blastfile\n\n\cl\n";
  my $in = Bio::SearchIO->new(-format => 'blast', -file => $blastfile);
  
  my $prasedOUThash;
  my $resultArrayIdx=0;
  my $maxBestinformToget=5;  #这个数字是用来 限制， $summuryHere这个变量的长度的
  
  
  my $SortArray;  my $SortArrayIdx=0;
  my $totalHit=0;
  
  # extraction of information for each result recursively
  while ( my $result = $in->next_result ) {
	  # the name of the query sequence    	
	  print  $result->query_name . "\t";
	  my $queryMinStart =  99999999999999999999;
    my $queryMaxEnd   = -99999999999999999999;
	  my $queryName=$result->query_name;
	  #if (  defined ( $allSelProInformHASHinput->{$queryName} )  ){ }
   	#else { my $find=0;
   	#  foreach my $qrNmKey (   keys (  %{ $allSelProInformHASHinput }  )   ){
   	#    if ($qrNmKey=~m/^$queryName/){$queryName=$qrNmKey;   $find=1;    }
   	#  }
   	#  if ($find==0){
   	#    die "\nNO \$queryName=$queryName was found and Changed:\$blastfile=$blastfile\n\n";
   	#  }
   	#}
   	$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_query_name'}=$queryName; 
   	$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_database_name'}=$result->database_name;  print  "database=".$result->database_name."\n";
    
    # the length of the query sequence #    	
    print  $result->query_length;
    $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_query_length'}=$result->query_length;
    
    $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_num_hits'}=$result->num_hits;  print  $result->num_hits;
    # output "no hits found" if there is no hits
    $totalHit+=$result->num_hits;  print "\n\$totalHit=$totalHit\t\$result->num_hits=$result->num_hits\n";
    
    my $summuryHere='';  my $summuryHereLength=0;
    if ( $result->num_hits == 0 ) {
		  print  "\tNo hits found\n";
    } 
    else {
		  my $count = 0;
      # process each hit recursively
		  while (my $hit = $result->next_hit) {
			  print  "\t" if ($count > 0);
                          # get the accession numbers of the hits
			  print  "\t" . $hit->accession . "\t";
                          # get the lengths of the hit sequences
        print  $hit->length . "\t";
                          # get the description of the hit sequences
			  print  $hit->description . "\t";
                          # get the E value of the hit
			  print  $hit->significance . "\t";
                          #get the bit score of the hit
			  print  $hit->bits . "\t";
        
        my $HitName=$hit->accession;
	      #if (  defined ( $allSelProInformHASHinput->{$HitName} )  ){ }
   	    #else { my $find=0;
   	    #  foreach my $qrNmKey (   keys (  %{ $allSelProInformHASHinput }  )   ){
   	    #    if ($qrNmKey=~m/^$HitName/){$HitName=$qrNmKey;   $find=1;    }
   	    #  }
   	    #  if ($find==0){
   	    #    die "\nNO \$HitName=$HitName was found and Changed:\$blastfile=$blastfile\n\n";
   	    #  }
       	#}
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_accessionNB'}   =$HitName;        
        #$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_accessionNB'}  =$hit->accession;
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hitSeqLength'}  =$hit->length;
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hitDescription'}=$hit->description;
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hitEvalue'}     =$hit->significance;
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hitScore'}      =$hit->bits;
        
        my $accIDinShort=$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_accessionNB'};
        $accIDinShort=~s/^(\S+)\s*.*$/$1/;  #214g2702C.paradoxa20150522;thioredoxin_reductase2;Contig53901_-3

        
        my $hspcount = 0; 
        my $totalIdentical=0;
        my $totalConserved=0;
        my $totConsePercent=0;
        
        
        my $queryCoverHash;  
        my $hitCoverHash;
        
        my $queryCoverRate=0;  
        my $qureyCovRatePct=0;
        my $hitCoverRate=0;
        my $hitCovRatePct=0;
        
        
                          # process the top HSP for the top hit
			  while (my $hsp = $hit->next_hsp) {
			  	
			  	
			  	
          print  "\t\t\t\t\t\t\t", if ($hspcount > 0);
                          	# get the frame of the query sequence
			  	print  $hsp->query->frame . "\t";
                                  # get the start and the end of the query sequence in the alignment
			  	print  $hsp->start('query') . "\t" . $hsp->end('query'). "\t";
                                  # get the start and the end of the hit sequence in the alignment
			  	print  $hsp->start('hit') . "\t" . $hsp->end('hit') . "\t";
                                  # get the similarity value
			  	printf  "%.1f" , ($hsp->frac_conserved * 100);
			  	print  "%\t";
                                  # get the identity value
			  	printf  "%.1f" , ($hsp->frac_identical * 100);
			  	$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_HspQurFrame'}     =$hsp->query->frame;
			  	$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_HspHitFrame'}     =$hsp->hit->frame;
			  	
			  	$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_HspQueryStart'}   =$hsp->start('query');   $queryCoverHash->[$hspcount]->{'_start'}=$hsp->start('query');
			  	$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_HspQueryEnd'}     =$hsp->end('query');     $queryCoverHash->[$hspcount]->{'_end'}  =$hsp->end('query');
			  	$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_HspHitStart'}     =$hsp->start('hit');     $hitCoverHash->[$hspcount]->{'_start'}=$hsp->start('hit');
			  	$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_HspHitEnd'}       =$hsp->end('hit');       $hitCoverHash->[$hspcount]->{'_end'}  =$hsp->end('hit');
			  	$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_HspConserverd'}   =$hsp->frac_conserved ;
			  	$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hspArray'}->[$hspcount]->{'_HspIdentical'}    =$hsp->frac_identical ;
			  	
			  	
			  	
			  	
			  	
			  	
			  	$totalIdentical+=$hsp->frac_identical ;
			  	$totalConserved+=$hsp->frac_conserved ;
			  	
			  	my $smallQNB=$hsp->start('query'); my $bigQNB=$hsp->end('query');
			  	if ($smallQNB > $bigQNB ){ $bigQNB=$hsp->start('query'); my $smallQNB=$hsp->end('query'); }
			  	
			  	if ( $smallQNB < $queryMinStart ){  $queryMinStart = $smallQNB;   }
			  	if ( $bigQNB   > $queryMaxEnd   ){  $queryMaxEnd   = $bigQNB;     }
			  	  
		      print  "%\n";
          $hspcount++;
        }
        
        $queryCoverHash=&overLayerHadle( $queryCoverHash );  
        $hitCoverHash  =&overLayerHadle( $hitCoverHash ); 
        
        my $queryCoverLength=&GetCuttedLength($queryCoverHash);
        my $hitCoverLength=&GetCuttedLength($hitCoverHash);
        
        my ($hit_queryStart, $hit_queryEnd)=@{ &getWholeMatchedStartEnd($queryCoverHash) };
        my ($hit_HitStart,   $hit_HitEnd)  =@{ &getWholeMatchedStartEnd($hitCoverHash)   };
        
        
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'queryCoverHash'}  =$queryCoverHash;        
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'queryCoverLength'}=$queryCoverLength;
        
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'hitCoverHash'}    =$hitCoverHash;
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'hitCoverLength'}  =$hitCoverLength;
       
        
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'hit_queryStart'}      =$hit_queryStart;
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'hit_queryEnd'}        =$hit_queryEnd;
        my $hit_queryRangeLength =(  (abs ($hit_queryEnd-$hit_queryStart)  ) + 1   );
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'hit_queryRangeLength'}=$hit_queryRangeLength;
        
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'hit_HitStart'}      =$hit_HitStart;
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'hit_HitEnd'}        =$hit_HitEnd;
        my $hit_HitRangeLength   =(  (abs ($hit_HitEnd-$hit_HitStart)  ) + 1   );
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'hit_HitRangeLength'}=$hit_HitRangeLength;
              
                   
        
        if ($hspcount > 0){
          $totalIdentical=$totalIdentical/$hspcount;
          $totalConserved=$totalConserved/$hspcount;
          $totConsePercent=100*$totalConserved;  $totConsePercent=sprintf "%.2f",$totConsePercent; $totConsePercent="$totConsePercent%";
          
          $queryCoverRate=$queryCoverLength/$hit_queryRangeLength;  
          $qureyCovRatePct   =100*$queryCoverRate;  $qureyCovRatePct=sprintf "%.2f",$qureyCovRatePct; $qureyCovRatePct="$qureyCovRatePct%";
          
        
          $hitCoverRate=$hitCoverLength/$hit_HitRangeLength;
          $hitCovRatePct     =100*$hitCoverRate;  $hitCovRatePct=sprintf "%.2f",$hitCovRatePct; $hitCovRatePct="$hitCovRatePct%";
          
          
        }
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hitTotalIdentical'}=$totalIdentical;
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hitTotalConserved'}=$totalConserved;
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'totConsePercent'   }=$totConsePercent;
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'queryCoverRate'    }=$queryCoverRate;
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'qureyCovRatePct'   }=$qureyCovRatePct;
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'hitCoverRate'      }=$hitCoverRate;
        $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'hitCovRatePct'     }=$hitCovRatePct;
          
        
        if ($summuryHereLength<$maxBestinformToget){
          $summuryHere.="${accIDinShort}($prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hitEvalue'}|$totConsePercent)";
          $summuryHereLength++;
        }
        
        $SortArray->{'EvalueSort'}->[$SortArrayIdx]->{'SotVal'}=$prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_hitArray'}->[$count]->{'_hitEvalue'};
        $SortArray->{'EvalueSort'}->[$SortArrayIdx]->{'QurIdx'}=$resultArrayIdx;
        $SortArray->{'EvalueSort'}->[$SortArrayIdx]->{'HitIdx'}=$count;
        
        $SortArray->{'TotCosSort'}->[$SortArrayIdx]->{'SotVal'}=$totalConserved;
        $SortArray->{'TotCosSort'}->[$SortArrayIdx]->{'QurIdx'}=$resultArrayIdx;
        $SortArray->{'TotCosSort'}->[$SortArrayIdx]->{'HitIdx'}=$count;
        
        $SortArray->{'TotIdtSort'}->[$SortArrayIdx]->{'SotVal'}=$totalIdentical;
        $SortArray->{'TotIdtSort'}->[$SortArrayIdx]->{'QurIdx'}=$resultArrayIdx;
        $SortArray->{'TotIdtSort'}->[$SortArrayIdx]->{'HitIdx'}=$count;
        $SortArrayIdx++; 
        
			  $count++;
        
        # flow control for the number of hits needed
        if (  (defined ($numHits) ) && ($numHits=~m/\d+/) && ($numHits>0)   ){
			    last if ($count == $numHits);
			  }
		  }
    }                                                       
                                                         
    $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_QueryMinStart'}=$queryMinStart;
    $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_QueryMaxEnd'}  =$queryMaxEnd;
    $prasedOUThash->{'_ResultArray'}->[$resultArrayIdx]->{'_Summury'}      =$summuryHere;
    $resultArrayIdx++;	
  }

  my $realSortArray; 
  foreach my $EachSortType (    sort {$a cmp $b} (   keys (  %{ $SortArray }  )   )    ){
  	my $realSortIdx=0;  my $totalSummury;
  	if ($EachSortType eq 'EvalueSort'){
      foreach my $eachHash  (   sort { $a->{'SotVal'} <=> $b->{'SotVal'} }  (  @{ $SortArray->{$EachSortType} }  )   ){
      	$realSortArray->{$EachSortType}->[$realSortIdx]=$eachHash;      $realSortIdx++;      	    	
      }
    }
    else{
    	foreach my $eachHash  (   sort { $b->{'SotVal'} <=> $a->{'SotVal'} }  (  @{ $SortArray->{$EachSortType} }  )   ){
      	$realSortArray->{$EachSortType}->[$realSortIdx]=$eachHash;      $realSortIdx++;
      }
    }
  } 
  

  $prasedOUThash->{'SortHash'}=$realSortArray;
  $prasedOUThash->{'totalHit'}=$totalHit;
  
  my $beastQuery=$prasedOUThash->{'_ResultArray'}->[ $prasedOUThash->{'SortHash'}->{'EvalueSort'}->[0]->{'QurIdx'} ]->{'_query_name'}; warn "\nbbbbbbb\$beastQuery=$beastQuery\n\n"; print "\nbbbbbbb\$beastQuery=$beastQuery\n\n";
  
##  my $writerhtml = new Bio::SearchIO::Writer::HTMLResultWriter(); 
##my $outhtml = new Bio::SearchIO(-writer => $writerhtml,
##   -file   => ">searchio.html");
##
### Loop through all the results, successively adding each one to 
### the bottom of the HTML report
##while ( $result = $searchio->next_result() ) {  
##    $outhtml->write_report($result);
##}
  if ($outBestBlastFile=~m/\S+/){
    my $writerhtml = new Bio::SearchIO::Writer::HTMLResultWriter();      #Bio::SearchIO->new(-format => 'blast', -file => ">$outBestBlastFile");
    my $out = new Bio::SearchIO( -writer => $writerhtml,
                                 -file   => ">$outBestBlastFile");
  
    my $in_2 = Bio::SearchIO->new(-format => 'blast', -file => $blastfile);
  
    OUTBEST: while ( my $result = $in_2->next_result ) {  
      if ($beastQuery eq $result->query_name){       print "\nbbbbbbb\$beastQuery=$beastQuery\t\t\$result->query_name=",$result->query_name,"\n\n"; warn "\nbbbbbbb\$beastQuery=$beastQuery\t\t\$result->query_name=",$result->query_name,"\n\n";
  	    $out->write_report($result); last OUTBEST;
  	  }
    }
  }
  print "\n\clEnd of Sub &BioPerlpraseBlast!\n\n";
  return $prasedOUThash;

}

#$queryCoverHash=&overLayerHadle( $queryCoverHash );     #$queryCoverHash->[$hspcount]->{'_start'}=$hsp->start('query');  #$queryCoverHash->[$hspcount]->{'_end'}  =$hsp->end('query');
sub overLayerHadle{    #pmAble#  #BioPerlBlastPraser_bl2seq_20170330的子函数 处理overlayyer的hash
  my ($inArrayHash)=@_;
  my $direction=SeqSegmentsTools::DirectionDetact($inArrayHash->[0]->{'_start'}, $inArrayHash->[0]->{'_end'});
  my @sortByHeadArray;
  if ( $direction==0 ){ #从小到大
  	@sortByHeadArray=sort{ $a->{'_start'} <=> $b->{'_start'} } @{ $inArrayHash };
  }
  else{  #从大到小
  	@sortByHeadArray=sort{ $b->{'_start'} <=> $a->{'_start'} } @{ $inArrayHash };
  }
  my $JointUnitsNBhash; my $jointIdx=0;  my $unitNumber=@sortByHeadArray;  my $lastMergeOrNotMergeHash;
  for (my $i=0; $i<@sortByHeadArray; $i++){  print "\$sortByHeadArray[$i]->{'_start'}=$sortByHeadArray[$i]->{'_start'}\t\$sortByHeadArray[$i]->{'_end'}=$sortByHeadArray[$i]->{'_end'}\n";
  	if ($unitNumber ==1){
  	  $JointUnitsNBhash->[0]=$sortByHeadArray[0];  print "00 \$JointUnitsNBhash->[$jointIdx]->{'_start'}=$JointUnitsNBhash->[$jointIdx]->{'_start'}\t\$JointUnitsNBhash->[$jointIdx]->{'_end'}=$JointUnitsNBhash->[$jointIdx]->{'_end'}\n";
  	}
  	else{
  		
  		if ( $i==0){   #第一项
  		  $lastMergeOrNotMergeHash=$sortByHeadArray[0];   print "11 \$lastMergeOrNotMergeHash->{'_start'}=$lastMergeOrNotMergeHash->{'_start'}\t\$lastMergeOrNotMergeHash->{'_end'}=$lastMergeOrNotMergeHash->{'_end'}\n";
  		}
  		
      if (  ($i>0) && ($i<($unitNumber-1) )   ){ #第二项到倒数第二项
       
        if ( SeqSegmentsTools::OverlayCheck($lastMergeOrNotMergeHash->{'_start'}, $lastMergeOrNotMergeHash->{'_end'}, $sortByHeadArray[$i]->{'_start'}, $sortByHeadArray[$i]->{'_end'}) ){
          my ($newHd, $newEd)=@{ SeqSegmentsTools::JointSegments_sameDirection_Up_or_down_noEqual($lastMergeOrNotMergeHash->{'_start'}, $lastMergeOrNotMergeHash->{'_end'}, $sortByHeadArray[$i]->{'_start'}, $sortByHeadArray[$i]->{'_end'}) };
          $lastMergeOrNotMergeHash->{'_start'}=$newHd; 
          $lastMergeOrNotMergeHash->{'_end'}  =$newEd; print "22 \$lastMergeOrNotMergeHash->{'_start'}=$lastMergeOrNotMergeHash->{'_start'}\t\$lastMergeOrNotMergeHash->{'_end'}=$lastMergeOrNotMergeHash->{'_end'}\n";
        }
        else {
          $JointUnitsNBhash->[$jointIdx]=$lastMergeOrNotMergeHash; print "33 \$JointUnitsNBhash->[$jointIdx]->{'_start'}=$JointUnitsNBhash->[$jointIdx]->{'_start'}\t\$JointUnitsNBhash->[$jointIdx]->{'_end'}=$JointUnitsNBhash->[$jointIdx]->{'_end'}\n"; $jointIdx++; 
          $lastMergeOrNotMergeHash=$sortByHeadArray[$i]; print "33 \$lastMergeOrNotMergeHash->{'_start'}=$lastMergeOrNotMergeHash->{'_start'}\t\$lastMergeOrNotMergeHash->{'_end'}=$lastMergeOrNotMergeHash->{'_end'}\n";
        }
        
      }
      
      if ( $i == ($unitNumber-1) ){         #最后一项
      	
      	if ($i ==0){  #既是最后一项，又是第一项
      	  $JointUnitsNBhash->[$jointIdx]=$sortByHeadArray[0];      $jointIdx++;
      	}
      	else{      	  #是非第一项的 最后一项
      	  if ( SeqSegmentsTools::OverlayCheck($lastMergeOrNotMergeHash->{'_start'}, $lastMergeOrNotMergeHash->{'_end'}, $sortByHeadArray[$i]->{'_start'}, $sortByHeadArray[$i]->{'_end'}) ){
            my ($newHd, $newEd)=@{ SeqSegmentsTools::JointSegments_sameDirection_Up_or_down_noEqual($lastMergeOrNotMergeHash->{'_start'}, $lastMergeOrNotMergeHash->{'_end'}, $sortByHeadArray[$i]->{'_start'}, $sortByHeadArray[$i]->{'_end'}) };
            $lastMergeOrNotMergeHash->{'_start'}=$newHd; 
            $lastMergeOrNotMergeHash->{'_end'}  =$newEd;   print "44 \$lastMergeOrNotMergeHash->{'_start'}=$lastMergeOrNotMergeHash->{'_start'}\t\$lastMergeOrNotMergeHash->{'_end'}=$lastMergeOrNotMergeHash->{'_end'}\n";
            
            $JointUnitsNBhash->[$jointIdx]=$lastMergeOrNotMergeHash;print "44 \$JointUnitsNBhash->[$jointIdx]->{'_start'}=$JointUnitsNBhash->[$jointIdx]->{'_start'}\t\$JointUnitsNBhash->[$jointIdx]->{'_end'}=$JointUnitsNBhash->[$jointIdx]->{'_end'}\n"; $jointIdx++;
          }
          else {
          	
          	$JointUnitsNBhash->[$jointIdx]=$lastMergeOrNotMergeHash;print "55 \$JointUnitsNBhash->[$jointIdx]->{'_start'}=$JointUnitsNBhash->[$jointIdx]->{'_start'}\t\$JointUnitsNBhash->[$jointIdx]->{'_end'}=$JointUnitsNBhash->[$jointIdx]->{'_end'}\n"; $jointIdx++;
            $JointUnitsNBhash->[$jointIdx]=$sortByHeadArray[$i];     print "55 \$JointUnitsNBhash->[$jointIdx]->{'_start'}=$JointUnitsNBhash->[$jointIdx]->{'_start'}\t\$JointUnitsNBhash->[$jointIdx]->{'_end'}=$JointUnitsNBhash->[$jointIdx]->{'_end'}\n"; $jointIdx++;
            
          }       
        }
      	
      }
    }
  } 
  return $JointUnitsNBhash;
}   



sub GetCuttedLength{  #pmAble#  #BioPerlBlastPraser_bl2seq_20170330的子函数 获得 cutted的 片段的长度
  my ($inArrayHash)=@_;
  
  my $coverLength;   print "\n\n\n\n\n";
  for (my $i=0; $i<@{ $inArrayHash }; $i++){
  	my $segMLength=(  ( abs ($inArrayHash->[$i]->{'_end'}-$inArrayHash->[$i]->{'_start'}) ) + 1  );  #print "\$inArrayHash->[$i]->{'_end'}=$inArrayHash->[$i]->{'_end'}\t\$inArrayHash->[$i]->{'_start'}=$inArrayHash->[$i]->{'_start'}\t\t\$segMLength=$segMLength\t\t";
  	$coverLength+=$segMLength;                                                                       #print "\$coverLength=$coverLength\n";
  }                                                                                                  #print "\n\n\n\n\n";
  return $coverLength;
}

sub getWholeMatchedStartEnd{   #pmAble#  #BioPerlBlastPraser_bl2seq_20170330的子函数 获得 match的片段的 最开头 和 最结尾
  my ($inArrayHash)=@_;
  
  
  my $arrayLength=@{ $inArrayHash }; my $lastElmentNB=$arrayLength-1;
  my $lastEnd =$inArrayHash->[$lastElmentNB]->{'_end'};  #warn "\$lastEnd =\$inArrayHash->[\$lastElmentNB]->{'_end'}=$lastEnd =\$inArrayHash->[$lastElmentNB]->{'_end'}\n\n";
  my $fstStart=$inArrayHash->[0]->{'_start'};            #warn "\$fstStart=\$inArrayHash->[0]->{'_start'}=$fstStart=\$inArrayHash->[0]->{'_start'}\n\n";  #sleep (10);
  
  
  return [$fstStart, $lastEnd];
}

#my $FlnTableString = BlastHandle::BuildKeyHead_txt_from_tabularFile($inFile, $outFile);
sub BuildKeyHead_txt_from_tabularFile{  #将没有head的 blast或者 diamond的 tablur文件 加上头
	my ($inFile, $outFile)=@_;
	
	my $warnMsgBody="\nIn package  InFileHandle,\tIn sub BuildKeyHead_txt_from_tabularFile,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";  
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	my $HeadKeyArray=BlastHandle::BuildTabularHeadArray( );
	my $HeadExplmAry=BlastHandle::BuildTabularHeadExplainArray( );
	
	my $FlnTableString = InFileHandle::Add_key_head_forNOheadTable($inFile, $outFile, $HeadKeyArray, $HeadExplmAry);
	return $FlnTableString;
	
}

#my $outArrayRef=BlastHandle::BuildTabularHeadArray( );
sub BuildTabularHeadArray{  # blast 或 diamond的 tablur格式，的表头信息
	
	                  # 0           1          2          3           4            5         6         7        8          9        10           11 
	my $outArrayRef=[ 'qseqid',  'sseqid',  'pident',  'length',  'mismatch',  'gapopen',  'qstart',  'qend',  'sstart',  'send',  'evalue',  'bitscore' ];
	return $outArrayRef;
}

#my $outArrayRef=BlastHandle::BuildTabularHeadExplainArray( );
sub BuildTabularHeadExplainArray{  # blast 或 diamond的 tablur格式，的表头信息的解释信息
	
	                                               
	my $outArrayRef=[ 
	                  'query (e.g., gene) sequence id',                   #0 
	                  'subject (e.g., reference genome) sequence id',     #1
	                  'percentage of identical matches',                  #2
	                  'alignment length',                                 #3 
	                  'number of mismatches',                             #4 
	                  'number of gap openings',                           #5 
	                  'start of alignment in query',                      #6 
	                  'end of alignment in query',                        #7 
	                  'start of alignment in subject',                    #8 
	                  'end of alignment in subject',                      #9 
	                  'expect value',                                     #10
	                  'bit score'                                         #11
	                ];
	return $outArrayRef;
}




#新的blast praser，这个获得特定类型氨基酸的速度快一点，
sub BioPerlBlastPraser20190928_for_AA_char{       #    my $onlyFindingAAhash=BlastHandle::BioPerlBlastPraser20190928_for_AA_char ($blastfile, $search_AminoAcid, $report_type, $idxFnm[, $numHits] ) #pmAble#   #解析Blast文件
  my ($blastfile, $search_AminoAcid, $report_type, $idxFnm, $numHits )=@_;   
  
  
  my $warnMsgBody="\nIn package  BlastHandle,\tIn sub BioPerlBlastPraser20190928_for_AA_char,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $subCalInfom=DirFileHandle::print_SubCallerInform;
  
  #my $search_AminoAcid='*';
  
  if (   (  defined ( $idxFnm )  ) && ( $idxFnm=~m/\S+/ ) && (  -e ( $idxFnm )  )    ){ 	  }
  else{		DieWork::Just_dieWork( $die_MsgHead."\n \$idxFnm=$idxFnm should be a defined index file  !!  $!\n\n\n".$subCalInfom ); 	}
  
  #warn  "\n\nNow in Sub &BioPerlpraseBlast!\n";      warn  "Input 1:\$blastfile=$blastfile\t" if (defined ($blastfile)); warn  "Input 2:\$outBestBlastFile=$outBestBlastFile\t" if (defined ($outBestBlastFile)); warn  "Input 3:\$numHits=$numHits\t" if (defined ($numHits));  warn  "Input 4:\$report_type=$report_type\t" if (defined ($report_type));  warn  "\n\n";
  print "\n\cl\nNow in Sub &BioPerlpraseBlast!\n:";  print "Input 1:\$blastfile=$blastfile\t" if (defined ($blastfile)); 
  
  print "Input 2:\$numHits=$numHits\t" if (defined ($numHits));  print "Input 3:\$report_type=$report_type\t" if (defined ($report_type));  print "\n\cl\n";
  
  my $in;
  if (   ( defined ($report_type) )  && ($report_type=~m/^\s*tblastn\s*$/i)   ){ 
  	$report_type='tblastn';
    $in = Bio::SearchIO->new(-format => 'blast', -file => $blastfile, -report_type => $report_type); warn "1 \$report_type=$report_type\n";
  }
  elsif (   ( defined ($report_type) )  && ($report_type=~m/^\s*blastx\s*$/i)   ){ 
  	$report_type='blastx';
    $in = Bio::SearchIO->new(-format => 'blast', -file => $blastfile, -report_type => $report_type); warn "2 \$report_type=$report_type\n";
  }
  else {
  	$in = Bio::SearchIO->new(-format => 'blast', -file => $blastfile                              ); warn "3 \$report_type=$report_type\n";
  }
  #sleep(3);
  
  my $prasedOUThash;
  my $resultArrayIdx=0;    
  my $maxBestinformToget=5;  #这个数字是用来 限制， $summuryHere这个变量的长度的
  
  my $onlyFindingAAhash;   my $onlyFindingAAhash_queryIdx=0;
  
  my $QueryHit_excel_HASH;
  
  my $SortArray;  my $SortArrayIdx=0;
  my $totalHit=0;
  
  my $Check_finish_HASH;
  
  # extraction of information for each result recursively
  while ( my $result = $in->next_result ) {
	  # the name of the query sequence    	
	  print  $result->query_name . "\t";
	  my $queryMinStart =  99999999999999999999;
    my $queryMaxEnd   = -99999999999999999999;
	  my $queryName=$result->query_name;
	  
	  my $onlyFindingAA_orNot_inQury=0; 
	  
	  $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_0_0_query_name'}=$queryName; 
   	
   	$Check_finish_HASH->{$queryName}=1;  #这个hash，记录blast结果中所有的query，最终用来判断 输入的fasta文件是否被blast完整计算
   	
   	my $databaseFind=$result->database_name;                                                                                                   #print  "\n20190322-0-1 database=".$databaseFind."\n";
   	#if (   (  defined ( $databaseFind )  ) && ( $databaseFind=~m/\S+/)   ){  #&& (  -e ( $databaseFind )  )    )  
   	#	                                                                                                                                       #print  "\n20190322-0-2 database=".$databaseFind."\n"; 		
   	#} 
   	
   	
   	
   	$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_0_1_query_fastaFile'}=$databaseFind;                                                    print  "\n20190322-0-5 database=".$databaseFind."\n";
  
	  $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_0_2_query_indexFile'}=$idxFnm;      
	 
	  
   
    my $query_length=$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_1_query_length'}=$result->query_length;
    
    $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_2_0_num_hits'}=$result->num_hits;  print  $result->num_hits;
    # output "no hits found" if there is no hits
    $totalHit+=$result->num_hits;  print "\n\$totalHit=$totalHit\t\$result->num_hits=$result->num_hits\n";
    
    my $summuryHere='';  my $summuryHereLength=0;
    if ( $result->num_hits == 0 ) {
		  print  "\tNo hits found\n";
    } 
    else {
		  my $count = 0;   my $onlyFindingAAhash_hitCount=0;
      # process each hit recursively
		  while (my $hit = $result->next_hit) {
			  print  "\t" if ($count > 0);
                          # get the accession numbers of the hits
			  print  "\t" . $hit->accession . "\t";
                          # get the lengths of the hit sequences
        print  $hit->length . "\t";
                          # get the description of the hit sequences
			  print  $hit->description . "\t";
                          # get the E value of the hit
			  print  $hit->significance . "\t";
                          #get the bit score of the hit
			  print  $hit->bits . "\t";
        
        my $HitName=$hit->accession;
        
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_0_0_locus'}   =$hit->locus();        
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_0_1__name'}   =$hit->name();  
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_0_2_giNub'}   =&GetGiNumber( $hit->name() );  
           
        #$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_0_3_accessionNB'}  =$hit->accession;
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_0_3_accessionNB'}   =$HitName;
        
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_0_4_hitSeqLength'}  =$hit->length;
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_0_5_hitDescription'}=$hit->description;
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_0_6_hitEvalue'}     =$hit->significance;
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_0_7_hitScore'}      =$hit->bits;
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_0_8_algorithm'}      =$hit->algorithm;
        
        
        my $accIDinShort=$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_0_3_accessionNB'};
        $accIDinShort=~s/^(\S+)\s*.*$/$1/;  #214g2702C.paradoxa20150522;thioredoxin_reductase2;Contig53901_-3

        
        my $hspcount = 0; 
        my $totalIdentical=0;
        my $totalConserved=0;
        
        my $totConsePercent=0;  #
        
        
        my $queryCoverHash;  
        my $hitCoverHash;
        
        my $queryCoverRate=0;  
        my $qureyCovRatePct=0;
        my $total_qry_covrAte=0;
        my $total_qry_covrAte_pct=0;
        my $hitCoverRate=0;
        my $hitCovRatePct=0;
        
        my $onlyFindingAA_orNot_inHit=0;  # 用来判断在一个Hit中是否 找到要找的特定的AA的char
        
                          # process the top HSP for the top hit
			  while (my $hsp = $hit->next_hsp) {
			  	
			  	
			  	
          print  "\t\t\t\t\t\t\t", if ($hspcount > 0);
                          	# get the frame of the query sequence
			  	print  $hsp->query->frame . "\t";
                                  # get the start and the end of the query sequence in the alignment
			  	print  $hsp->start('query') . "\t" . $hsp->end('query'). "\t";
                                  # get the start and the end of the hit sequence in the alignment
			  	print  $hsp->start('hit') . "\t" . $hsp->end('hit') . "\t";
                                  # get the similarity value
			  	printf  "%.1f" , ($hsp->frac_conserved * 100);
			  	print  "%\t";
                                  # get the identity value
			  	printf  "%.1f" , ($hsp->frac_identical * 100);
			  	
			  	#$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'_HspQurFrame'}        =$hsp->frame('query');  #oldVersion   #bioperl 得到的blast的hsp的frame是不对的
			  	#$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'_HspHitFrame'}        =$hsp->frame('hit');    #oldVersion   #bioperl 得到的blast的hsp的frame是不对的
			  	
			  	
			  	
			  	#$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'_Hsp_1stfm_3frm_2qur'}       =$hsp->frame('query');        my $QrAbsFm=$hsp->frame('query')+1;  #bioperl 得到的blast的hsp的frame是不对的
			  	#$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'_Hsp_1stfm_4frm_2hit'}       =$hsp->frame('hit');          my $HtAbsFm=$hsp->frame('hit')+1;   #bioperl 得到的blast的hsp的frame是不对的
			  	$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_0_0_Hsp_1stfm_1stn_1qur'}       =$hsp->strand('query');       my $QrZF; if ($hsp->strand('query') >=0){$QrZF='+';}else {$QrZF='-';} 
			  	$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_0_1_Hsp_1stfm_2stn_1hit'}       =$hsp->strand('hit');         my $HtZF; if ($hsp->strand('hit')   >=0){$HtZF='+';}else {$HtZF='-';} 
			  	#$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'_Hsp_1stfm_1qur'}       =$QrZF.$QrAbsFm;   #bioperl 得到的blast的hsp的frame是不对的
			  	#$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'_Hsp_1stfm_2hit'}       =$HtZF.$HtAbsFm;   #bioperl 得到的blast的hsp的frame是不对的
			  	$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_0_2_Hsp_1stfm_3ZoF_1qur'}       =$QrZF;       
			  	$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_0_3_Hsp_1stfm_4ZoF_1hit'}       =$HtZF;         
			  	
			  	$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_1_0_HspQueryStart'}      =$hsp->start('query');  #oldVersion
			  	$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_1_1_HspQueryEnd'}        =$hsp->end('query');    #oldVersion
			  	$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_1_2_HspHitStart'}        =$hsp->start('hit');    #oldVersion
			  	$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_1_3_HspHitEnd'}          =$hsp->end('hit');      #oldVersion
			  	
			  	my $QryHspStt=			  	
			  	$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_2_0_Hsp_2StEd_1Qr_1St'}        =$hsp->start('query');       $queryCoverHash->[$hspcount]->{'_start'}=$hsp->start('query'); 
			  	my $QryHspEnd=
			  	$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_2_1_Hsp_2StEd_1Qr_2Ed'}        =$hsp->end('query');         $queryCoverHash->[$hspcount]->{'_end'}  =$hsp->end('query');   
			  	$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_2_2_Hsp_2StEd_2Ht_1St'}        =$hsp->start('hit');         $hitCoverHash->[$hspcount]->{'_start'}  =$hsp->start('hit');     
			  	$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_2_3_Hsp_2StEd_2Ht_2Ed'}        =$hsp->end('hit');           $hitCoverHash->[$hspcount]->{'_end'}    =$hsp->end('hit');       
			  	

			  	
			  	$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_3_0_HspConserverd'}      =$hsp->frac_conserved ;
			  	$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_3_1_HspIdentical'}       =$hsp->frac_identical ;
			  	



			  	
			  	#$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'_Hsp_1Qry_string'}      =$hsp->query_string ;    #oldVersion
			  	#$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'_Hsp_3Hit_string'}      =$hsp->hit_string ;      #oldVersion
			  	#$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'_Hsp_2Hom_string'}      =$hsp->homology_string ; #oldVersion
			  	
			  	my $Qr_aln_string=
			  	$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_4_1_Hsp_3WhSting_1Qr'}      =$hsp->query_string ;
			  	my $mc_aln_string=
			  	$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_4_2_Hsp_3WhSting_2Ho'}      =$hsp->homology_string ;
			  	my $Ht_aln_string=
			  	$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_4_3_Hsp_3WhSting_3Ht'}      =$hsp->hit_string ;
			  	
			  	my $alignStrLenth=length $Qr_aln_string;
			  	
			  	#寻找blastx的 query中特定 氨基酸 在DNA中的位置
			  	#my $search_AminoAcid='*';
			  	my $findingCharArray; 
			  	my $find_char_or_nor_inHSP=0;
			  	
			  	#计算特定氨基酸局部区域时，局部的大小，由以下参数设定 $flanking_length=5设定，表示在目标 字符 ，左右两侧各取5个字符 如 54321*12345，一共11个字符
			  	my $flanking_length=5;
			  	
			  	if ($report_type eq 'blastx'){
			  		$findingCharArray=&findingIdexOfAword ( $search_AminoAcid, $Qr_aln_string ); 
			  		if (   ( ref ($findingCharArray) eq 'ARRAY' ) && (  ( @{ $findingCharArray } ) >0  )   ) { #如果找到 对应的 字符，则进行下面的操作
			  			#DirFileHandle::PrintAndWarnDumper ($findingCharArray, "\n20190928-0-0-1\n"); 
			  			
			  			my $aln_obj  = $hsp->get_aln();                 #DirFileHandle::PrintAndWarnDumper ($aln_obj,   "\n20190928-0-0-0-0 \$aln_obj=\n");
			  			my $AlnQryObj=$aln_obj->get_seq_by_pos(1);      #DirFileHandle::PrintAndWarnDumper ($AlnQryObj, "\n20190928-0-0-0-0 \$AlnQryObj=\n");  
			  			my $AlnHitObj=$aln_obj->get_seq_by_pos(2);      #DirFileHandle::PrintAndWarnDumper ($AlnHitObj, "\n20190928-0-0-0-0 \$AlnHitObj=\n");
			        
			       
			        my $realMatchCharIdx=0;
			        for (my $charIdx=0; $charIdx<@{ $findingCharArray }; $charIdx++){ #对所有的字符进行下面操作
			          
			          my $CharPOS_inAlign=$findingCharArray->[$charIdx];
			          my $QurColOBJ=$AlnQryObj->location_from_column($CharPOS_inAlign); #根据query的字符位置，找到的location obj的类型应该是excat，start和end都是相同的才行
			          if (   (  defined ( $QurColOBJ->location_type() )  ) && (  $QurColOBJ->location_type() eq 'EXACT'  ) && (  defined ( $QurColOBJ->start() )  ) && (  defined ( $QurColOBJ->end() )  ) && ( $QurColOBJ->start() == $QurColOBJ->end() )   ){
			          	my $QurCharPostion=$QurColOBJ->start();
			          	
			          	my $HitColOBJ=$AlnHitObj->location_from_column($CharPOS_inAlign);#根据query字符位置，在hit中找对应位置，类型也必须是excat，start和end也要相同
			          	if (   (  defined ( $HitColOBJ->location_type() )  ) && (  $HitColOBJ->location_type() eq 'EXACT'  ) && (  defined ( $HitColOBJ->start() )  ) && (  defined ( $HitColOBJ->end() )  ) && ( $HitColOBJ->start() == $HitColOBJ->end() )   ){
			          	  my $hitCharPostion=$HitColOBJ->start();
			          	  
			          	  #下面进行 相应的 计算 测试 验证
			          	  my $realCharPosStt;  
			          	  my $realCharPosEnd;  
			          	  my $QryDNAchar;
			          	  
			          	  if    ( $QrZF eq '+' ){
			          	  	$realCharPosStt=3*($QurCharPostion-$QryHspStt)+$QryHspStt;
			          	  	$realCharPosEnd=$realCharPosStt+2;
			          	  	$QryDNAchar= FastaFileHandle::feach_seqString_from_idx_File($queryName, $idxFnm, $realCharPosStt, $realCharPosEnd);
			          	  }
			          	  elsif ( $QrZF eq '-'){  #query为-的时候，query的$QryHspStt实际上是 数字小的那个 位置，逻辑上的真正stt应该是$QryHspEnd 
			          	  	$realCharPosStt=3*($QurCharPostion-$QryHspEnd)+$QryHspEnd;
			          	  	$realCharPosEnd=$realCharPosStt-2;
			          	  	$QryDNAchar= FastaFileHandle::feach_seqString_from_idx_File($queryName, $idxFnm, $realCharPosEnd, $realCharPosStt);			          	  	
			          	  	$QryDNAchar= FastaFileHandle::ReverseComplementString ($QryDNAchar);
			          	  }
			          	  
			          	  $QryDNAchar= uc ($QryDNAchar);
			          	  my $QryPEPchar=FastaFileHandle::CodonTable_standard($QryDNAchar);
			          	  
			          	  print "\n\n\n\n\$QurCharPostion=$QurCharPostion \$realCharPosStt=$realCharPosStt \$realCharPosEnd=$realCharPosEnd\n";
			          	  print "\$hitCharPostion=$hitCharPostion\n";
			          	  print "\n\n\$Qr_aln_string=$Qr_aln_string\n";
			  	          print "\$Ht_aln_string=$Ht_aln_string\n\n";
			  	          my $hitPosChar=FastaFileHandle::FactchSeqJustSubStr($Ht_aln_string, '+', $CharPOS_inAlign, $CharPOS_inAlign);
			          	  
			          	  print "\n20191004-0-0-0 $CharPOS_inAlign, $flanking_length, $Qr_aln_string\n";
			          	  #获得被检测位置的 特定 氨基酸 两侧固定长度的 局部区间的 序列
			          	  my $CharPos_inAln_fromZero=$CharPOS_inAlign-1;
			          	  my ( $qur_Char_flankStr, $qur_Char_inSubStrPos)= @{ String_Work::GetSubString_from_middlePoint($CharPos_inAln_fromZero, $flanking_length, $Qr_aln_string) }; print "\n20191004-0-1-0 \$qur_Char_inSubStrPos=$qur_Char_inSubStrPos \$qur_Char_flankStr=$qur_Char_flankStr\n\n";
			          	  my ( $mch_Char_flankStr, $mch_Char_inSubStrPos)= @{ String_Work::GetSubString_from_middlePoint($CharPos_inAln_fromZero, $flanking_length, $mc_aln_string) }; print "\n20191004-0-1-1 \$mch_Char_inSubStrPos=$mch_Char_inSubStrPos \$mch_Char_flankStr=$mch_Char_flankStr\n\n";
			          	  my ( $hit_Char_flankStr, $hit_Char_inSubStrPos)= @{ String_Work::GetSubString_from_middlePoint($CharPos_inAln_fromZero, $flanking_length, $Ht_aln_string) }; print "\n20191004-0-1-2 \$hit_Char_inSubStrPos=$hit_Char_inSubStrPos \$hit_Char_flankStr=$hit_Char_flankStr\n\n";
			          	  my   $Char_highLigt_Str=FastaFileHandle::Build_simple_showChar_string($search_AminoAcid, $qur_Char_inSubStrPos);   print "\n20191004-0-1-3 \$Char_highLigt_Str=$Char_highLigt_Str\n\n";
			          	  my ( $Local_Align_Score, $ShowBlSMScoreLine, $usef_but_throw_Array)=@{ BLOSUMalignSCORE::Score_for_Alignment ($qur_Char_flankStr, $hit_Char_flankStr) };   print "\n20191004-0-1-4 \$Local_Align_Score=$Local_Align_Score\n\n";
			          	  
			  	          $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_6_0_Hsp_AA_match'}->{$search_AminoAcid}->[$realMatchCharIdx]->{'6_0_0_Qur_stn_ZoF'}=$QrZF;
			          	  $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_6_0_Hsp_AA_match'}->{$search_AminoAcid}->[$realMatchCharIdx]->{'6_0_1_Qur_pos_stt'}=$realCharPosStt;
			          	  $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_6_0_Hsp_AA_match'}->{$search_AminoAcid}->[$realMatchCharIdx]->{'6_0_2_Qur_pos_end'}=$realCharPosEnd;
			          	  $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_6_0_Hsp_AA_match'}->{$search_AminoAcid}->[$realMatchCharIdx]->{'6_0_3_Qur_DNAchar'}=$QryDNAchar;
			          	  $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_6_0_Hsp_AA_match'}->{$search_AminoAcid}->[$realMatchCharIdx]->{'6_0_4_Qur_PEPchar'}=$QryPEPchar;
			          	  
			          	  $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_6_0_Hsp_AA_match'}->{$search_AminoAcid}->[$realMatchCharIdx]->{'6_1_1_hit_pos_stt'}=$hitCharPostion;
			          	  $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_6_0_Hsp_AA_match'}->{$search_AminoAcid}->[$realMatchCharIdx]->{'6_1_3_hit_PEPchar'}=$hitPosChar;
			          	  
			          	  $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_6_0_Hsp_AA_match'}->{$search_AminoAcid}->[$realMatchCharIdx]->{'6_2_0_locAlnScore'}=$Local_Align_Score;
			          	  
			          	  $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_6_0_Hsp_AA_match'}->{$search_AminoAcid}->[$realMatchCharIdx]->{'6_2_1_ChrShow_Str'}=$Char_highLigt_Str;
			          	  $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_6_0_Hsp_AA_match'}->{$search_AminoAcid}->[$realMatchCharIdx]->{'6_2_2_Qur_sub_Str'}=$qur_Char_flankStr;
			          	  $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_6_0_Hsp_AA_match'}->{$search_AminoAcid}->[$realMatchCharIdx]->{'6_2_3_Mch_sub_Str'}=$mch_Char_flankStr;
			          	  $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_6_0_Hsp_AA_match'}->{$search_AminoAcid}->[$realMatchCharIdx]->{'6_2_4_Shw_Sco_Str'}=$ShowBlSMScoreLine;
			          	  $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_6_0_Hsp_AA_match'}->{$search_AminoAcid}->[$realMatchCharIdx]->{'6_2_5_Hit_sub_Str'}=$hit_Char_flankStr;
			          	  
			          	  $realMatchCharIdx++;
			          	  $find_char_or_nor_inHSP=1;
			          	  $onlyFindingAA_orNot_inHit=1;                  #$onlyFindingAAhash
			          	  $onlyFindingAA_orNot_inQury=1; #my $onlyFindingAAhash;   my $onlyFindingAAhash_queryIdx=0;
			          	}
			          	
			          }
			          else{ #query条件不符合，则直接die
			          	DieWork::Just_dieWork( $die_MsgHead."\n\$QurColOBJ->location_type()=$QurColOBJ->location_type() should eq EXACT \$QurColOBJ->start()=$QurColOBJ->start() should == \$QurColOBJ->start()=$QurColOBJ->start() : $!".$subCalInfom );
			          }
			         
			          
			        }
			  			
			  		}
			  	}
			  	if ( $find_char_or_nor_inHSP==1 ){
			  		$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_4_0_Hsp_3WhSting_0AA'}      =FastaFileHandle::Build_show_char_string($search_AminoAcid, $findingCharArray);  
			  	
			  	}
			  	
			  	
			  	$totalIdentical+=$hsp->frac_identical ;
			  	$totalConserved+=$hsp->frac_conserved ;
			  	
			  	my $smallQNB=$hsp->start('query'); my $bigQNB=$hsp->end('query');
			  	if ($smallQNB > $bigQNB ){ $bigQNB=$hsp->start('query'); my $smallQNB=$hsp->end('query'); }
			  	
			  	if ( $smallQNB < $queryMinStart ){  $queryMinStart = $smallQNB;   }
			  	if ( $bigQNB   > $queryMaxEnd   ){  $queryMaxEnd   = $bigQNB;     }
			  	  
		      print  "%\n";
          $hspcount++;
        }
        
        $queryCoverHash=BlastHandle::overLayerHadle( $queryCoverHash );  
        $hitCoverHash  =BlastHandle::overLayerHadle( $hitCoverHash ); 
        
        my $queryCoverLength=BlastHandle::GetCuttedLength($queryCoverHash);
        my $hitCoverLength=BlastHandle::GetCuttedLength($hitCoverHash);
        
        my ($hit_queryStart, $hit_queryEnd)=@{ BlastHandle::getWholeMatchedStartEnd($queryCoverHash) };
        my ($hit_HitStart,   $hit_HitEnd)  =@{ BlastHandle::getWholeMatchedStartEnd($hitCoverHash)   };
        
        
        #$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_1_0_queryCoverHash'}  =$queryCoverHash;        
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_1_1_queryCoverLength'}=$queryCoverLength;
        
        #$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_1_3_hitCoverHash'}    =$hitCoverHash;
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_1_4_hitCoverLength'}  =$hitCoverLength;
       
        
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_1_5_hit_queryStart'}      =$hit_queryStart;
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_1_6_hit_queryEnd'}        =$hit_queryEnd;
        my $hit_queryRangeLength =(  (abs ($hit_queryEnd-$hit_queryStart)  ) + 1   );
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_1_7_hit_queryRangeLength'}=$hit_queryRangeLength;
        
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_1_8_hit_HitStart'}      =$hit_HitStart;
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_1_9_hit_HitEnd'}        =$hit_HitEnd;
        my $hit_HitRangeLength   =(  (abs ($hit_HitEnd-$hit_HitStart)  ) + 1   );
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_1_a_hit_HitRangeLength'}=$hit_HitRangeLength;
        
        
        
        if ($hspcount > 0){
          $totalIdentical=$totalIdentical/$hspcount;
          $totalConserved=$totalConserved/$hspcount;
          $totConsePercent=100*$totalConserved;  $totConsePercent=sprintf "%.2f",$totConsePercent; $totConsePercent="$totConsePercent%";
          
          $queryCoverRate=$queryCoverLength/$hit_queryRangeLength;  
          $qureyCovRatePct   =100*$queryCoverRate;  $qureyCovRatePct=sprintf "%.2f",$qureyCovRatePct; $qureyCovRatePct="$qureyCovRatePct%";
          
          $total_qry_covrAte=$queryCoverLength/$query_length;
          $total_qry_covrAte_pct=100*$total_qry_covrAte;  $total_qry_covrAte_pct=sprintf "%.2f",$total_qry_covrAte_pct; $total_qry_covrAte_pct="$total_qry_covrAte_pct%";
          
        
          $hitCoverRate=$hitCoverLength/$hit_HitRangeLength;
          $hitCovRatePct     =100*$hitCoverRate;  $hitCovRatePct=sprintf "%.2f",$hitCovRatePct; $hitCovRatePct="$hitCovRatePct%";
          
          
        }
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_2_0_hitTotalIdentical'}=$totalIdentical;
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_2_1_hitTotalConserved'}=$totalConserved;
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_2_2_totConsePercent'   }=$totConsePercent;
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_2_3_queryCoverRate'    }=$queryCoverRate;
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_2_4_qureyCovRatePct'   }=$qureyCovRatePct;
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_2_5_total_qry_covrAte' }=$total_qry_covrAte;
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_2_6_total_qry_covrAte_pct'}=$total_qry_covrAte_pct;
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_2_7_hitCoverRate'      }=$hitCoverRate;
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_2_8_hitCovRatePct'     }=$hitCovRatePct;
        
        
        if ($summuryHereLength<$maxBestinformToget){
        	my $singleSum="${accIDinShort}($prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_0_6_hitEvalue'}|$totConsePercent)";
          $summuryHere.=$singleSum;
          $summuryHereLength++;
          $QueryHit_excel_HASH->{$queryName}->{$HitName}=$singleSum;
        }
        
        $SortArray->{'0_1_0_EvalueSort'}->[$SortArrayIdx]->{'2_1_0_SotVal'}=$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_0_6_hitEvalue'};
        $SortArray->{'0_1_0_EvalueSort'}->[$SortArrayIdx]->{'2_1_1_QurIdx'}=$resultArrayIdx;
        $SortArray->{'0_1_0_EvalueSort'}->[$SortArrayIdx]->{'2_1_2_HitIdx'}=$count;
        
        $SortArray->{'0_1_1_TotCosSort'}->[$SortArrayIdx]->{'2_1_0_SotVal'}=$totalConserved;
        $SortArray->{'0_1_1_TotCosSort'}->[$SortArrayIdx]->{'2_1_1_QurIdx'}=$resultArrayIdx;
        $SortArray->{'0_1_1_TotCosSort'}->[$SortArrayIdx]->{'2_1_2_HitIdx'}=$count;
        
        $SortArray->{'0_1_2_TotIdtSort'}->[$SortArrayIdx]->{'2_1_0_SotVal'}=$totalIdentical;
        $SortArray->{'0_1_2_TotIdtSort'}->[$SortArrayIdx]->{'2_1_1_QurIdx'}=$resultArrayIdx;
        $SortArray->{'0_1_2_TotIdtSort'}->[$SortArrayIdx]->{'2_1_2_HitIdx'}=$count;
        $SortArrayIdx++; 
        
			  
			  if ( $onlyFindingAA_orNot_inHit==1){
			    $onlyFindingAAhash->{'0_0_1_Query_array'}->[$onlyFindingAAhash_queryIdx]->{'2_0_Z_hitArray'         }->[$onlyFindingAAhash_hitCount]=Storable::dclone( $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count] );
			    $onlyFindingAAhash_hitCount++;  
			    
	        $onlyFindingAAhash->{'0_0_1_Query_array'}->[$onlyFindingAAhash_queryIdx]->{'2_0_0_0_query_name'          }=$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_0_0_query_name'};
			    $onlyFindingAAhash->{'0_0_1_Query_array'}->[$onlyFindingAAhash_queryIdx]->{'2_0_0_1_query_fastaFile'     }=$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_0_1_query_fastaFile'};
			    $onlyFindingAAhash->{'0_0_1_Query_array'}->[$onlyFindingAAhash_queryIdx]->{'2_0_0_2_query_indexFile'     }=$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_0_2_query_indexFile'};
			    $onlyFindingAAhash->{'0_0_1_Query_array'}->[$onlyFindingAAhash_queryIdx]->{'2_0_1_query_length'          }=$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_1_query_length'};
			    #$onlyFindingAAhash->{'0_0_1_Query_array'}->[$onlyFindingAAhash_queryIdx]->{'2_0_2_num_hits'              }=$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_2_num_hits'};
			    
			    $onlyFindingAAhash->{'0_0_1_Query_array'}->[$onlyFindingAAhash_queryIdx]->{'2_0_2_1_findAA_char_hitNumbs'}=$onlyFindingAAhash_hitCount;
			  }
			  
			  $count++;
			  
			  
			  
			  
        
        # flow control for the number of hits needed
        if (  (defined ($numHits) ) && ($numHits=~m/\d+/) && ($numHits>0)   ){
			    last if ($count == $numHits);
			  }
		  }
    }                                                       
                                                         
    $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_3_QueryMinStart'}=$queryMinStart;
    $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_4_QueryMaxEnd'}  =$queryMaxEnd;
    $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_5_Summury'}      =$summuryHere;
    warn "\$resultArrayIdx=$resultArrayIdx\n";
    $resultArrayIdx++;	
    
    if($onlyFindingAA_orNot_inQury==1) {
    	
    	$onlyFindingAAhash_queryIdx++;
    }
    
  }
  
  if(0){
  	my $realSortArray; 
    foreach my $EachSortType (    sort {$a cmp $b} (   keys (  %{ $SortArray }  )   )    ){
    	my $realSortIdx=0;  my $totalSummury;
    	if ($EachSortType eq '0_1_0_EvalueSort'){
        foreach my $eachHash  (   sort { $a->{'2_1_0_SotVal'} <=> $b->{'2_1_0_SotVal'} }  (  @{ $SortArray->{$EachSortType} }  )   ){
        	$realSortArray->{$EachSortType}->[$realSortIdx]=$eachHash;      $realSortIdx++;      	    	
        }
      }
      else{
      	foreach my $eachHash  (   sort { $b->{'2_1_0_SotVal'} <=> $a->{'2_1_0_SotVal'} }  (  @{ $SortArray->{$EachSortType} }  )   ){
        	$realSortArray->{$EachSortType}->[$realSortIdx]=$eachHash;      $realSortIdx++;
        }
      }
    } 
    $prasedOUThash->{'0_0_2_SortHash'}=$realSortArray;
  
  }
  
  

  
  $prasedOUThash->{'0_0_0_0_totalQry'}=$resultArrayIdx;
  $prasedOUThash->{'0_0_0_1_totalHit'}=$totalHit;
  $onlyFindingAAhash->{'0_0_0_0_the_Char_to_____find'}=$search_AminoAcid;
  $onlyFindingAAhash->{'0_0_0_1_total__________Query'}=$resultArrayIdx;
  $onlyFindingAAhash->{'0_0_0_2_totalFoundChar_Query'}=$onlyFindingAAhash_queryIdx;
  
  $onlyFindingAAhash->{'0_1_0_CheckFinishHASH'       }=$Check_finish_HASH;
  #return $prasedOUThash;
  return $onlyFindingAAhash;
##my $beastQuery=$prasedOUThash->{'0_0_1_Query_array'}->[ $prasedOUThash->{'SortHash'}->{'0_1_0_EvalueSort'}->[0]->{'2_1_1_QurIdx'} ]->{'2_0_0_0_query_name'};  #warn "\nbbbbbbb\$beastQuery=$beastQuery\n\n"; print "\nbbbbbbb\$beastQuery=$beastQuery\n\n";
  
##  my $writerhtml = new Bio::SearchIO::Writer::HTMLResultWriter(); 
##my $outhtml = new Bio::SearchIO(-writer => $writerhtml,
##   -file   => ">searchio.html");
##
### Loop through all the results, successively adding each one to 
### the bottom of the HTML report
##while ( $result = $searchio->next_result() ) {  
##    $outhtml->write_report($result);
##}

  #DirFileHandle::PrintDumper                    ('temp.QueryHitHash', $QueryHit_excel_HASH)     if ( ref($QueryHit_excel_HASH) eq 'HASH' ) ; 
  #MatrixCsvChange::PrintExcelTableFileFromHASH  ('temp.QueryHit.excel.txt', $QueryHit_excel_HASH, 0 ) if ( ref($QueryHit_excel_HASH) eq 'HASH' ) ; 
  

 
  
  

}


sub Extract_subQueryHit_HASH{  #my $Filter_out_HASH=BlastHandle::Extract_subQueryHit_HASH( $inBltPrasedHASH, $Query_Covreage_cutoff, $Hit_Identity_cutoff, $Evalue_cutoff );
	my ($inBltPrasedHASH, $Query_Covreage_cutoff, $Hit_Identity_cutoff, $Evalue_cutoff)=@_;
	
	my $warnMsgBody="\nIn package  BlastHandle,\tIn sub Extract_subQueryHit_HASH,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $subCalInfom=DirFileHandle::print_SubCallerInform;
  
  if  (   (  defined ( $inBltPrasedHASH )  ) && (  ref ( $inBltPrasedHASH ) eq 'HASH'  )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$inBltPrasedHASH=$inBltPrasedHASH should be a HASH  !!  $!\n\n\n".$subCalInfom ); 	}
	
	my $Query_Covreage_cutoff_Filter=0;
  if (   (  defined ( $Query_Covreage_cutoff )  ) && ( $Query_Covreage_cutoff=~m/\S+/ )    ){ 	  
    if (  ( $Query_Covreage_cutoff=~m/\d+(\.\d+)?/ ) && ( 0 < $Query_Covreage_cutoff ) && ( $Query_Covreage_cutoff <= 1 )  ){
    	$Query_Covreage_cutoff_Filter=1;
    }
    else{		DieWork::Just_dieWork( $die_MsgHead."\n \$Query_Covreage_cutoff=$Query_Covreage_cutoff should be a defined number > 0, <=1  !!  $!\n\n\n".$subCalInfom ); 	}
  }
  
  my $Hit_Identity_cutoff_Filter=0;
  if (   (  defined ( $Hit_Identity_cutoff )  ) && ( $Hit_Identity_cutoff=~m/\S+/ )    ){ 	  
    if (  ( $Hit_Identity_cutoff=~m/\d+(\.\d+)?/ ) && ( 0 < $Hit_Identity_cutoff ) && ( $Hit_Identity_cutoff <= 1 )  ){
    	$Hit_Identity_cutoff_Filter=1;
    }
    else{		DieWork::Just_dieWork( $die_MsgHead."\n \$Hit_Identity_cutoff=$Hit_Identity_cutoff should be a defined number > 0, <=1  !!  $!\n\n\n".$subCalInfom ); 	}
  }
  
  my $Evalue_cutoff_Filter=0;
  if (   (  defined ( $Evalue_cutoff )  ) && ( $Evalue_cutoff=~m/\S+/ )    ){ 	  
    $Evalue_cutoff_Filter=1;    
  }
  
  my $Filter_out_HASH;  
  my $filter_out_totalQry=0;
  my $filter_out_totalHit=0;
  
  foreach my $Lev_0_key (    sort {$a cmp $b} (   keys (  %{ $inBltPrasedHASH }  )   )    ){  DieWork::Print_and_warn( "\n20191022-0-0-0 \$Lev_0_key=$Lev_0_key\n\n" ); 
  	if (  defined ( $inBltPrasedHASH->{$Lev_0_key} )  ) {
  		if    (  ref ( $inBltPrasedHASH->{$Lev_0_key} ) eq 'ARRAY'  ){
  			if ( $Lev_0_key eq '0_0_1_Query_array' ){
  				
  				my $found_qry_idx=0;
  				for ( my $i=0; $i<@{ $inBltPrasedHASH->{$Lev_0_key} }; $i++ ){        DieWork::Print_and_warn( "\n20191022-0-0-1 \$inBltPrasedHASH->{$Lev_0_key}->[$i]=$inBltPrasedHASH->{$Lev_0_key}->[$i]\n\n" ); 
  	        
  	        my $temp_lev_2_key_hash; # A temp hash to hold all the not hash or array value, if the $found_qry_or_not==1, then fill all these values into final outputhash
  	        
  	        my $found_qry_or_not=0;
  	        if  (   (  defined ( $inBltPrasedHASH->{$Lev_0_key}->[$i] )  ) && (  ref ( $inBltPrasedHASH->{$Lev_0_key}->[$i] ) eq 'HASH'  )   ){
  	            
  	            
  	            foreach my $Lev_2_key (    sort {$a cmp $b} (   keys (  %{ $inBltPrasedHASH->{$Lev_0_key}->[$i] }  )   )    ){  DieWork::Print_and_warn( "\n20191022-0-0-2 \$Lev_2_key=$Lev_2_key\n\n" ); 
  	              if (  defined ( $inBltPrasedHASH->{$Lev_0_key}->[$i]->{$Lev_2_key} )  ) {
  	                if    (  ref ( $inBltPrasedHASH->{$Lev_0_key}->[$i]->{$Lev_2_key} ) eq 'ARRAY'  ){
  	                	if ( $Lev_2_key eq '2_0_Z_hitArray' ){
  	                		
  	                		my $found_hit_idx=0;
  	                		for ( my $j=0; $j<@{ $inBltPrasedHASH->{$Lev_0_key}->[$i]->{$Lev_2_key} }; $j++ ){  DieWork::Print_and_warn( "\n20191022-0-0-3 \$inBltPrasedHASH->{$Lev_0_key}->[$i]->{$Lev_2_key}->[$j]=$inBltPrasedHASH->{$Lev_0_key}->[$i]->{$Lev_2_key}->[$j]\n\n" );
  	                			
  	                			my $found_hit_or_not=0;
  	                			if  (   (  defined ( $inBltPrasedHASH->{$Lev_0_key}->[$i]->{$Lev_2_key}->[$j] )  ) && (  ref ( $inBltPrasedHASH->{$Lev_0_key}->[$i]->{$Lev_2_key}->[$j] ) eq 'HASH'  )   ){
  	                			  
  	                			  ############## Check Hit to filter!!!! ###################
  	                			  
  	                			  if (   (  defined ( $inBltPrasedHASH->{$Lev_0_key}->[$i]->{$Lev_2_key}->[$j]->{'3_0_6_hitEvalue'} )  ) && (  $inBltPrasedHASH->{$Lev_0_key}->[$i]->{$Lev_2_key}->[$j]->{'3_0_6_hitEvalue'}  )   ){ }
  	                			  else{		DieWork::Just_dieWork( $die_MsgHead."\n \$inBltPrasedHASH->{$Lev_0_key}->[$i]->{$Lev_2_key}->[$j]->{'3_0_6_hitEvalue'}=$inBltPrasedHASH->{$Lev_0_key}->[$i]->{$Lev_2_key}->[$j]->{'3_0_6_hitEvalue'} should be a defined number   !!  $!\n\n\n".$subCalInfom ); 	}
  	                			  if (   (  defined ( $inBltPrasedHASH->{$Lev_0_key}->[$i]->{$Lev_2_key}->[$j]->{'3_2_0_hitTotalIdentical'} )  ) && (  $inBltPrasedHASH->{$Lev_0_key}->[$i]->{$Lev_2_key}->[$j]->{'3_2_0_hitTotalIdentical'}  )   ){ }
  	                			  else{		DieWork::Just_dieWork( $die_MsgHead."\n \$inBltPrasedHASH->{$Lev_0_key}->[$i]->{$Lev_2_key}->[$j]->{'3_2_0_hitTotalIdentical'}=$inBltPrasedHASH->{$Lev_0_key}->[$i]->{$Lev_2_key}->[$j]->{'3_2_0_hitTotalIdentical'} should be a defined number > 0, <=1  !!  $!\n\n\n".$subCalInfom ); 	}
  	                			  if (   (  defined ( $inBltPrasedHASH->{$Lev_0_key}->[$i]->{$Lev_2_key}->[$j]->{'3_2_5_total_qry_covrAte'} )  ) && (  $inBltPrasedHASH->{$Lev_0_key}->[$i]->{$Lev_2_key}->[$j]->{'3_2_5_total_qry_covrAte'}  )   ){ }
  	                			  else{		DieWork::Just_dieWork( $die_MsgHead."\n \$inBltPrasedHASH->{$Lev_0_key}->[$i]->{$Lev_2_key}->[$j]->{'3_2_5_total_qry_covrAte'}=$inBltPrasedHASH->{$Lev_0_key}->[$i]->{$Lev_2_key}->[$j]->{'3_2_5_total_qry_covrAte'} should be a defined number > 0, <=1  !!  $!\n\n\n".$subCalInfom ); 	}
  	                			  
  	                			  my $hitEvalue        =$inBltPrasedHASH->{$Lev_0_key}->[$i]->{$Lev_2_key}->[$j]->{'3_0_6_hitEvalue'};          DieWork::Print_and_warn( "\n20191022-0-0-4 \$inBltPrasedHASH->{$Lev_0_key}->[$i]->{$Lev_2_key}->[$j]->{'3_0_6_hitEvalue'}=$inBltPrasedHASH->{$Lev_0_key}->[$i]->{$Lev_2_key}->[$j]->{'3_0_6_hitEvalue'}\n\n" ); 
  	                			  my $hitTotalIdentical=$inBltPrasedHASH->{$Lev_0_key}->[$i]->{$Lev_2_key}->[$j]->{'3_2_0_hitTotalIdentical'};
  	                			  my $total_qry_covrAte=$inBltPrasedHASH->{$Lev_0_key}->[$i]->{$Lev_2_key}->[$j]->{'3_2_5_total_qry_covrAte'};
  	                			  
  	                			  
  	                			  #real Check precess# 对于每一个要检查的项目设定一个值，最初是1，如果该条件不许要检查，则该项目为1（通过），否则如需要检查则设为0，然后再判断，通过为1，不通过为0.
  	                			                     # 这样，对每个需要检查的项目都得到一个判断是否通过的值 1通过 或 0没通过
  	                			                     # 如果所有的项目都通过，则 总体判定为通过，否则 只要有一项没通过，则 总体判定为不通过
  	                			  
  	                			  my $realEvalue_pass=1; if (  $Evalue_cutoff_Filter==1         ){ $realEvalue_pass=0; 
  	                			  	                                                               if (  $hitEvalue         <= $Evalue_cutoff         ){$realEvalue_pass=1;  } 
  	                			  	                                                             }	                                                             
  	                			  my $realHitIdt_pass=1; if (  $Hit_Identity_cutoff_Filter==1   ){ $realHitIdt_pass=0; 
  	                			  	                                                               if (  $hitTotalIdentical >= $Hit_Identity_cutoff   ){$realHitIdt_pass=1;  } 
  	                			  	                                                             }
  	                			  my $realQryCov_pass=1; if (  $Query_Covreage_cutoff_Filter==1 ){ $realQryCov_pass=0; 
  	                			  	                                                               if (  $total_qry_covrAte >= $Query_Covreage_cutoff ){$realQryCov_pass=1;  } 
  	                			  	                                                             }
  	                			  
  	                			  my $AllPass=$realQryCov_pass*$realHitIdt_pass*$realEvalue_pass;
  	                			  	                                                             
  	                			  
  	                			  
  	                			  if ($AllPass==1){
  	                			  	$found_qry_or_not=1;
  	                			  	$found_hit_or_not=1;
  	                			  	
  	                			  	DieWork::Print_and_warn( "\n20191022-0-0-5 \$inBltPrasedHASH->{$Lev_0_key}->[$i]->{$Lev_2_key}->[$j]=$inBltPrasedHASH->{$Lev_0_key}->[$i]->{$Lev_2_key}->[$j]\n\n" ); 
  	                			  	my $tempHASH=Storable::dclone ( $inBltPrasedHASH->{$Lev_0_key}->[$i]->{$Lev_2_key}->[$j] );
  	                			  	$Filter_out_HASH->{$Lev_0_key}->[$found_qry_idx]->{$Lev_2_key}->[$found_hit_idx]=$tempHASH;
  	                			  	
  	                			  	
  	                			  }
  	                			  
  	                			
  	                			}
  	                			
  	                			if ( $found_hit_or_not==1 ){
  	                				
  	                				$found_hit_idx++;  
  	                				$Filter_out_HASH->{$Lev_0_key}->[$found_qry_idx]->{'2_0_2_1_FilterOUT_num_hits'}=$found_hit_idx;
  	                				$filter_out_totalHit++;
  	                			}
  	                			
  	                		}
  	                		
  	                	}
  	                }
  	                elsif (  ref ( $inBltPrasedHASH->{$Lev_0_key}->[$i]->{$Lev_2_key} ) eq 'HASH'  ){
  	                }
  	                else{ # 当 $inBltPrasedHASH->{$Lev_0_key}->[$i]->{$Lev_2_key}  既不是 HASH的ref，也不是ARRAY的ref时
  	                	$temp_lev_2_key_hash->{$Lev_2_key}=$inBltPrasedHASH->{$Lev_0_key}->[$i]->{$Lev_2_key};
  	                }
  	              }
  	            }
  	          
  	          
  	        
  	        }
  	        if ( $found_qry_or_not==1 ){
  	        	
  	        	foreach my $temP_key (    sort {$a cmp $b} (   keys (  %{ $temp_lev_2_key_hash }  )   )    ){
  	        		$Filter_out_HASH->{$Lev_0_key}->[$found_qry_idx]->{$temP_key}=$temp_lev_2_key_hash->{$temP_key};
  	        	}     		
  	          
  	        	$found_qry_idx++;  $filter_out_totalQry++;
  	        }
  	        
  	        
  	        
  	        
          }
  			}
  		}
  		elsif (  ref ( $inBltPrasedHASH->{$Lev_0_key} ) eq 'HASH'  ){
  		}
  		else{
  			$Filter_out_HASH->{$Lev_0_key}= $inBltPrasedHASH->{$Lev_0_key} ;
  		}
  	}
  }
  
  $Filter_out_HASH->{'0_0_0_6_filter_out_totalQry'}=$filter_out_totalQry;
  $Filter_out_HASH->{'0_0_0_7_filter_out_totalHit'}=$filter_out_totalHit;
  
  return $Filter_out_HASH;
  
}

1;

##########################################################################################################################################
# 


