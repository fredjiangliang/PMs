
#!/usr/bin/perl -w

BEGIN{  push (@INC, '/home/fredjiang/fredProgs/PMs')  };

use strict;
use warnings;

use Storable ;
use Storable qw(nstore store_fd nstore_fd freeze thaw dclone);

use FastaFileHandle;
use SeqSegmentsTools;
use BlastHandle;
use DirFileHandle;
use DieWork;
use InFileHandle;
use GenBankHandle;
use EUtilitieswork;

package  DiamondWork;

my $JL_MutipleThread_doing_diamond_pl="/home/fredjiang/fredProgs/PLs/Diamond_multipleCore_working.pl";   
my $Diamond_followTBLASTN_pl         ="/home/fredjiang/fredProgs/PLs/Diamond_followTBLASTN.pl";

#DiamondWork::InFastaHash_out_MutiThrd_Work($InFastHash, $workingDir, $chunKLenthLmt, $coreToUse, $dimondType, $localDiamondDatabase);  
sub InFastaHash_out_MutiThrd_Work{ #输入 以一系列fasta文件的地址为key的hash的ref，进行多线程的diamond的分析, 注意： $fisrt_or_overwrite_chunk_mark 这个标识，在第一次运行和要进行 overwrite运行的时候， 在输入的时候需要写入 ovrt
  my ($InFastHash, $workingDir, $chunKLenthLmt, $coreToUse, $dimondType, $localDiamondDatabase )=@_;
  
  my $warnMsgBody="\nIn package  DiamondWork,\tIn sub InFastaHash_out_MutiThrd_Work,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $caller_inform=DirFileHandle::print_SubCallerInform;
  
  my $abNormFsHd='JLfa';
  
  my $subCkFasHd='JLol';
  my $segMentLth=10000;
  my $overLayLth=1000 ;
  
  
  if  (   (  defined ( $InFastHash )  ) && (  ref ( $InFastHash ) eq 'HASH'  )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$InFastHash=$InFastHash should be a HASH  !!  $!\n\n\n".$caller_inform ); 	}
	
	my ( $FastaIDKeyHASH, $ChunkPath_HASH );
	my ( $ovlCkForDmd_Hash, $chkPth_to_oDHash, $UplvID_to_oDHash );
	
	my $fisrt_or_overwrite_chunk_work_warn="\n\n\n\t\t\tAttention!!\n\n\tThe chunk make process will build differnt chunks every run!!\nPlease type \"ovrt\" if you are sure to overwrite or the first time running of chunk building process!!\nPlease type \"noru\" if you are sure to no to overwrite or no the first time running of chunk building process!!\n\n";
	warn  $fisrt_or_overwrite_chunk_work_warn;
	print $fisrt_or_overwrite_chunk_work_warn;
	sleep (1);
	my $fisrt_or_overwrite_chunk_mark=<STDIN>; chomp $fisrt_or_overwrite_chunk_mark;
	if   (   (  defined ( $fisrt_or_overwrite_chunk_mark )  ) && ( $fisrt_or_overwrite_chunk_mark=~m/\S+/ ) && (  $fisrt_or_overwrite_chunk_mark eq 'ovrt'  )   ){
	  #先直接进行文件分块，这时候，不进行单个fasta序列的内部分块
    ( $FastaIDKeyHASH, $ChunkPath_HASH )=@{ FastaChunkWork::BuildChunck_fasta_work($InFastHash, $workingDir, $chunKLenthLmt, $abNormFsHd, 'blastx') }; 
    
    #因为，要做Dimond blastx，所以要进行单个fasta序列的内部分块，目的是让所有的blastx输入DNA序列长度 相仿，防止太长的序列 取太少的结果出现
    ( $ovlCkForDmd_Hash, $chkPth_to_oDHash, $UplvID_to_oDHash )=@{ FastaChunkWork::BuildChunck_forDiamondBlastx($workingDir, $segMentLth, $overLayLth, $subCkFasHd ) }; 
    
	}
	elsif   (   (  defined ( $fisrt_or_overwrite_chunk_mark )  ) && ( $fisrt_or_overwrite_chunk_mark=~m/\S+/ ) && (  $fisrt_or_overwrite_chunk_mark eq 'noru'  )   ){
	  my $no_overwrite_chunk_building_Msg="\n\n\tAttention!!! the chunk building process was not running this time!!\n\n\n";
	  warn  $no_overwrite_chunk_building_Msg;
	  print $no_overwrite_chunk_building_Msg;
	}	
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$fisrt_or_overwrite_chunk_mark=$fisrt_or_overwrite_chunk_mark should be type in STDIN value, as ovrt or noru  !!  $!\n\n\n$fisrt_or_overwrite_chunk_work_warn\n\n\n".$caller_inform ); 	}
	
	
  
  
 
  
  
  #再对 $workingDir 中所有分段的序列 进行多线程的 diamond比对，及相关的结果的解析，将各部计算的结果 输出到一个hash文件中
  my $out_ChunkKey_Hash_file=DiamondWork::Doing_RunDiamond_vs_NRdb__MutipleThread($workingDir, $coreToUse, $dimondType, $localDiamondDatabase);
	
	#输出 结果位置hash，该hash文件 可用于其它进一步分析
	return $out_ChunkKey_Hash_file;
}



#   DiamondWork::Doing_RunDiamond_vs_NRdb__MutipleThread($workingDir, $coreToUse, $dimondType, $localDiamondDatabase);
sub Doing_RunDiamond_vs_NRdb__MutipleThread{  #调用 pl脚本 运行diamond
	my ($workingDir, $coreToUse, $dimondType, $localDiamondDatabase)=@_;
	#my ($in_ChunkKey_Hash_file, $in_FastaKey_Hash_file, $coreToUse, $out_ChunkKey_Hash_file, $out_FastaKey_Hash_file, $dimondType)=@_;
	
	
	my $warnMsgBody="\nIn package  DiamondWork,\tIn sub Doing_RunDiamond_vs_NRdb__MutipleThread,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	
	$JL_MutipleThread_doing_diamond_pl=$JL_MutipleThread_doing_diamond_pl;  #模块内部 全局变量
	
	my $notNR_database_own_database_mark='';
	if   (   (  defined ( $localDiamondDatabase )  ) && ( $localDiamondDatabase=~m/\S+/ ) && (  ( -e $localDiamondDatabase )  )   ){
	  $notNR_database_own_database_mark=" -d $localDiamondDatabase";
	}
	
	
	if (   (  defined ( $workingDir )  ) && ( $workingDir=~m/^\S+$/ )   ){				}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$workingDir=$workingDir should be right path to hold all the outHash and out chunk files !!  $!\n\n\n".$caller_inform );	}
	
	my $in_FastaKey_Hash_file=$workingDir."/0_0_2_FastaIDKeyHSH.txt";
	my $in_ChunkKey_Hash_file=$workingDir."/0_0_3_ChunkIDKeyHSH.txt";
  
	
	my $ovlCkForDmd_Hash_file=$workingDir."/0_3_1_ovlCkForDmdHs.txt";
	my $chkPth_to_oDHash_file=$workingDir."/0_3_3_CkPathtoOlDHs.txt";
	
	my $out_FastaKey_Hash_file=$workingDir."/0_5_2_Dmd_FstIDKyHS.txt";
	my $out_ChunkKey_Hash_file=$workingDir."/0_5_3_Dmd_CukIDKyHS.txt";

	if   (   (  defined ( $in_ChunkKey_Hash_file )  ) && ( $in_ChunkKey_Hash_file=~m/\S+/ ) && (  ( -e $in_ChunkKey_Hash_file )  )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$in_ChunkKey_Hash_file=$in_ChunkKey_Hash_file should be a defined not empty file  !!  $!\n\n\n".$caller_inform ); 	}
	my $in_ChunkKey_Hash=Storable::retrieve( $in_ChunkKey_Hash_file );
	if  (   (  defined ( $in_ChunkKey_Hash )  ) && (  ref ( $in_ChunkKey_Hash ) eq 'HASH'  )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$in_ChunkKey_Hash=$in_ChunkKey_Hash should be a HASH  !!  $!\n\n\n".$caller_inform ); 	}
	
	if   (   (  defined ( $in_FastaKey_Hash_file )  ) && ( $in_FastaKey_Hash_file=~m/\S+/ ) && (  ( -e $in_FastaKey_Hash_file )  )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$in_FastaKey_Hash_file=$in_FastaKey_Hash_file should be a defined not empty file  !!  $!\n\n\n".$caller_inform ); 	}
	my $in_FastaKey_Hash=Storable::retrieve( $in_FastaKey_Hash_file );
	if  (   (  defined ( $in_FastaKey_Hash )  ) && (  ref ( $in_FastaKey_Hash ) eq 'HASH'  )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$in_FastaKey_Hash=$in_FastaKey_Hash should be a HASH  !!  $!\n\n\n".$caller_inform ); 	}

  if   (   (  defined ( $coreToUse )  ) && ( $coreToUse=~m/\d+/ ) && ( $coreToUse > 0 )     ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$coreToUse=$coreToUse should be a number > 0  !!  $!\n\n\n".$caller_inform ); 	}
	
  
  if   (   (  defined ( $out_ChunkKey_Hash_file )  ) && ( $out_ChunkKey_Hash_file=~m/\S+/ )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$out_ChunkKey_Hash_file=$out_ChunkKey_Hash_file should be a defined path string  !!  $!\n\n\n".$caller_inform ); 	}
	
	if   (   (  defined ( $out_FastaKey_Hash_file )  ) && ( $out_FastaKey_Hash_file=~m/\S+/ )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$out_FastaKey_Hash_file=$out_FastaKey_Hash_file should be a defined path string  !!  $!\n\n\n".$caller_inform ); 	}
	
	my $diamondBlastx_orNot=0;
	if (   (  defined ( $dimondType )  ) && ( $dimondType=~m/\S+/)   ){
  	$dimondType=lc $dimondType;
  	if (  ( $dimondType eq 'blastp' ) || ( $dimondType eq 'blastx' )  ){
  		if ( $dimondType eq 'blastx' ){ $diamondBlastx_orNot=1;}
  	}
  	else{
  		DieWork::Just_dieWork( $die_MsgHead."\n \$dimondType=$dimondType should be blastp or blastx !!  $!\n\n\n".$caller_inform ); 
  	}  	
  } 
  else{
  	DieWork::Just_dieWork( $die_MsgHead."\n \$dimondType=$dimondType should not be empty or not defined !!  $!\n\n\n".$caller_inform ); 
  } 
  
  if ( $diamondBlastx_orNot==1){
  	if   (   (  defined ( $ovlCkForDmd_Hash_file )  ) && ( $ovlCkForDmd_Hash_file=~m/\S+/ ) && (  ( -e $ovlCkForDmd_Hash_file )  )   ){}
	  else{		DieWork::Just_dieWork( $die_MsgHead."\n \$ovlCkForDmd_Hash_file=$ovlCkForDmd_Hash_file should be a defined not empty file  !!  $!\n\n\n".$caller_inform ); 	}
	  my $ovlCkForDmd_Hash=Storable::retrieve( $ovlCkForDmd_Hash_file );
	  if  (   (  defined ( $ovlCkForDmd_Hash )  ) && (  ref ( $ovlCkForDmd_Hash ) eq 'HASH'  )   ){}
	  else{		DieWork::Just_dieWork( $die_MsgHead."\n \$ovlCkForDmd_Hash=$ovlCkForDmd_Hash should be a HASH  !!  $!\n\n\n".$caller_inform ); 	}
	  
	  if   (   (  defined ( $chkPth_to_oDHash_file )  ) && ( $chkPth_to_oDHash_file=~m/\S+/ ) && (  ( -e $chkPth_to_oDHash_file )  )   ){}
	  else{		DieWork::Just_dieWork( $die_MsgHead."\n \$chkPth_to_oDHash_file=$chkPth_to_oDHash_file should be a defined not empty file  !!  $!\n\n\n".$caller_inform ); 	}
	  my $in_FastaKey_Hash=Storable::retrieve( $in_FastaKey_Hash_file );
	  if  (   (  defined ( $in_FastaKey_Hash )  ) && (  ref ( $in_FastaKey_Hash ) eq 'HASH'  )   ){}
	  else{		DieWork::Just_dieWork( $die_MsgHead."\n \$in_FastaKey_Hash=$in_FastaKey_Hash should be a HASH  !!  $!\n\n\n".$caller_inform ); 	}
  	
  	system ( "perl $JL_MutipleThread_doing_diamond_pl  -i $chkPth_to_oDHash_file  -j $ovlCkForDmd_Hash_file -c $coreToUse -o $out_FastaKey_Hash_file -t $out_ChunkKey_Hash_file -p $dimondType $notNR_database_own_database_mark");
  	
  }
	else{
	  system ( "perl $JL_MutipleThread_doing_diamond_pl  -i $in_ChunkKey_Hash_file  -j $in_FastaKey_Hash_file -c $coreToUse -o $out_FastaKey_Hash_file -t $out_ChunkKey_Hash_file -p $dimondType $notNR_database_own_database_mark");	
	}	
  
	return $out_ChunkKey_Hash_file;	
	
}


#   DiamondWork::Doing_multiple_tBlastN_work($workingDir, $coreToUse);
sub Doing_multiple_tBlastN_work{  #多线程，对 diamond的结果，进行tblastn计算
	my ($workingDir, $coreToUse)=@_;
	my $warnMsgBody="\nIn package  DiamondWork,\tIn sub Doing_RunDiamond_vs_NRdb__MutipleThread,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	
	$Diamond_followTBLASTN_pl=$Diamond_followTBLASTN_pl;  #模块内部 全局变量
	
	DieWork::Check_FileDirExist_or_DIE          ( $workingDir, "\$workingDir", $die_MsgHead, $caller_inform  );

	system ( "perl $Diamond_followTBLASTN_pl  -i $workingDir   -c $coreToUse ");
	
}

#  my $diamond_work_is_already_done = DiamondWork::Check_diamond_done( $in_fasta_hash, $diamond_out_hash_file );
sub Check_diamond_done{   # 通过检查 diamond的 prased的hash，来检查其是否完成计算,这个检查 只适用于 包括所有的结果的hash的检查（也就是说，有些结果没有hit，也会记录下来，这是blast和diamond的基本alignment格式，如果是tabular格式则不适合）
	my ($in_fasta_hash, $diamond_out_hash_file)=@_;
	
	my $warnMsgBody="\nIn package  DiamondWork,\tIn sub Check_diamond_done,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";   
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	my $diamond_work_is_already_done=0;
	
	if   (   (  defined ( $diamond_out_hash_file )  ) && ( $diamond_out_hash_file=~m/\S+/ ) && (  ( -e $diamond_out_hash_file )  )   ){  print "20190105-0-0-0-0 \$diamond_out_hash_file=$diamond_out_hash_file \n";   }
	else{	 print "20190105-0-0-0-1 \$diamond_out_hash_file=$diamond_out_hash_file \n"; 	return  $diamond_work_is_already_done; } # DieWork::Just_dieWork( $die_MsgHead."\n \$diamond_out_hash_file=$diamond_out_hash_file should be a defined not empty file  !!  $!\n\n\n".$caller_inform ); 	}
	my $diamond_out_hash=Storable::retrieve( $diamond_out_hash_file );
	if  (   (  defined ( $diamond_out_hash )  ) && (  ref ( $diamond_out_hash ) eq 'HASH'  )   ){ print "20190105-0-0-0-2 \$diamond_out_hash=$diamond_out_hash \n"; }
	else{		print "20190105-0-0-0-3 \$diamond_out_hash=$diamond_out_hash \n"; return  $diamond_work_is_already_done; } # DieWork::Just_dieWork( $die_MsgHead."\n \$diamond_out_hash=$diamond_out_hash should be a HASH  !!  $!\n\n\n".$caller_inform ); 	}
	if  (   (  defined ( $diamond_out_hash->{'_ResultArray'} )  ) && (  ref ( $diamond_out_hash->{'_ResultArray'} ) eq 'ARRAY'  )   ){}
	else{		return  $diamond_work_is_already_done; } # DieWork::Just_dieWork( $die_MsgHead."\n \$diamond_out_hash->{'_ResultArray'}=$diamond_out_hash->{'_ResultArray'} should be a HASH  !!  $!\n\n\n".$caller_inform ); 	}
	
	
	if  (   (  defined ( $in_fasta_hash )  ) && (  ref ( $in_fasta_hash ) eq 'HASH'  )   ){ print "20190105-0-0-0-4 \$in_fasta_hash=$in_fasta_hash \n"; }
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$in_fasta_hash=$in_fasta_hash should be a HASH  !!  $!\n\n\n".$caller_inform ); 	}
	
	my $diammondOutHashQueryHash;
  foreach my $eachQry (  @{ $diamond_out_hash->{'_ResultArray'} }  ){  print "20190105-1-1-0  \$eachQry->{'_query_name'}=$eachQry->{'_query_name'}\n";    #'_query_name'
  #foreach my $eachQry (   keys (  %{ $diamond_out_hash }  )   ){  print "20190105-1-0  \$eachQry=$eachQry\n";  
    #$diammondOutHashQueryHash->{ $eachQry }=1 if (   (  defined ( $eachQry )  ) && ( $eachQry=~m/\S+/ )   );
    $diammondOutHashQueryHash->{ $eachQry->{'_query_name'} }=1 if (   (  defined ( $eachQry )  ) && (  ref ( $eachQry ) eq 'HASH' ) && (  defined ( $eachQry->{'_query_name'} )  ) && ( $eachQry->{'_query_name'}=~m/\S+/ )   );
    print "20190105-1-2-0  $diamond_out_hash_file \$diammondOutHashQueryHash->{ \$eachQry->{'_query_name'} }=\$diammondOutHashQueryHash->{ $eachQry->{'_query_name'} }=$diammondOutHashQueryHash->{ $eachQry->{'_query_name'} }\n"; 
  }
	
	my $InFasta_simple_Hash;
  foreach my $eachFastaID (   keys (  %{ $in_fasta_hash }  )   ){  print "20190105-1-1-1  \$eachFastaID=$eachFastaID\n";    #'eachFastaID'
    $InFasta_simple_Hash->{ $eachFastaID }=1 if (   (  defined ( $eachFastaID )  ) && ( $eachFastaID=~m/\S+/ )   );
    print "20190105-1-2-1  $diamond_out_hash_file \$InFasta_simple_Hash->{ \$eachFastaID }=\$InFasta_simple_Hash->{ $eachFastaID }\n"; 
  }
	
	print "20190105-2-2-0 \$diammondOutHashQueryHash=$diammondOutHashQueryHash, \$InFasta_simple_Hash=$InFasta_simple_Hash CheckHashKey_exactly_theSAME starting... \n";
	my $Same_Hash_or_not=ArrayHashChange::CheckHashKey_exactly_theSAME ( $diammondOutHashQueryHash, $InFasta_simple_Hash ); print "20190105-2-2-1 \$Same_Hash_or_not=$Same_Hash_or_not\n";
  if ( $Same_Hash_or_not == 1){
    $diamond_work_is_already_done=1;
  }
	return  $diamond_work_is_already_done;
	
}

# 这里的 $diamond_out_hash_file 是利用这个函数获得的 BlastHandle::BioPerlBlastPraser20190928_for_AA_char
#  my $diamond_work_is_already_done = DiamondWork::Check_diamond_done_EasyVersion( $in_fasta_hash, $diamond_out_hash_file );
sub Check_diamond_done_EasyVersion{   # 通过检查 diamond的 prased的hash，来检查其是否完成计算,这个检查 只适用于 包括所有的结果的hash的检查（也就是说，有些结果没有hit，也会记录下来，这是blast和diamond的基本alignment格式，如果是tabular格式则不适合）
	my ($in_fasta_hash, $diamond_out_hash_file)=@_;
	
	my $warnMsgBody="\nIn package  DiamondWork,\tIn sub Check_diamond_done_EasyVersion,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";   
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	my $diamond_work_is_already_done=0;
	
	if   (   (  defined ( $diamond_out_hash_file )  ) && ( $diamond_out_hash_file=~m/\S+/ ) && (  ( -e $diamond_out_hash_file )  )   ){  print "20190105-0-0-0-0 \$diamond_out_hash_file=$diamond_out_hash_file \n";   }
	else { DieWork::Just_dieWork( $die_MsgHead."\n \$diamond_out_hash_file=$diamond_out_hash_file should be a defined not empty file  !!  $!\n\n\n".$caller_inform ); 	}
	my $diamond_out_hash=Storable::retrieve( $diamond_out_hash_file );
	if  (   (  defined ( $diamond_out_hash )  ) && (  ref ( $diamond_out_hash ) eq 'HASH'  )   ){ print "20190105-0-0-0-2 \$diamond_out_hash=$diamond_out_hash \n"; }
	else{				
		$diamond_work_is_already_done=2;   print "20190105-0-0-0-3-0 \$diamond_out_hash=$diamond_out_hash \n"; 
		return  $diamond_work_is_already_done; 
	} # DieWork::Just_dieWork( $die_MsgHead."\n \$diamond_out_hash=$diamond_out_hash should be a HASH  !!  $!\n\n\n".$caller_inform ); 	}
	if  (   (  defined ( $diamond_out_hash->{'0_1_0_CheckFinishHASH'} )  ) && (  ref ( $diamond_out_hash->{'0_1_0_CheckFinishHASH'} ) eq 'HASH'  )   ){}
	else{	
		$diamond_work_is_already_done=2;   print "20190105-0-0-0-3-1 \$diamond_out_hash->{'0_1_0_CheckFinishHASH'}=$diamond_out_hash->{'0_1_0_CheckFinishHASH'} \n"; 	
		return  $diamond_work_is_already_done; 
	} # DieWork::Just_dieWork( $die_MsgHead."\n \$diamond_out_hash->{'_ResultArray'}=$diamond_out_hash->{'_ResultArray'} should be a HASH  !!  $!\n\n\n".$caller_inform ); 	}
	
	
	if  (   (  defined ( $in_fasta_hash )  ) && (  ref ( $in_fasta_hash ) eq 'HASH'  )   ){ print "20190105-0-0-0-4 \$in_fasta_hash=$in_fasta_hash \n"; }
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$in_fasta_hash=$in_fasta_hash should be a HASH  !!  $!\n\n\n".$caller_inform ); 	}
	
	my $diammondOutHashQueryHash=$diamond_out_hash->{'0_1_0_CheckFinishHASH'};
  
	
	my $InFasta_simple_Hash;
  foreach my $eachFastaID (   keys (  %{ $in_fasta_hash }  )   ){  print "20190105-1-1-1  \$eachFastaID=$eachFastaID\n";    #'eachFastaID'
    $InFasta_simple_Hash->{ $eachFastaID }=1 if (   (  defined ( $eachFastaID )  ) && ( $eachFastaID=~m/\S+/ )   );
    print "20190105-1-2-1  $diamond_out_hash_file \$InFasta_simple_Hash->{ \$eachFastaID }=\$InFasta_simple_Hash->{ $eachFastaID }\n"; 
  }
	
	print "20190105-2-2-0 \$diammondOutHashQueryHash=$diammondOutHashQueryHash, \$InFasta_simple_Hash=$InFasta_simple_Hash CheckHashKey_exactly_theSAME starting... \n";
	my $Same_Hash_or_not=ArrayHashChange::CheckHashKey_exactly_theSAME ( $diammondOutHashQueryHash, $InFasta_simple_Hash ); print "20190105-2-2-1 \$Same_Hash_or_not=$Same_Hash_or_not\n";
  if ( $Same_Hash_or_not == 1){
    $diamond_work_is_already_done=1;
  }
  else {
  	$diamond_work_is_already_done=2;
  }
	return  $diamond_work_is_already_done;
	
}

#  my $diamond_work_is_already_done = DiamondWork::Check_diamond_done_notRealRight( $in_fasta_hash, $diamond_out_hash_file );
sub Check_diamond_done_notRealRight{   # 通过检查 diamond的 prased的hash，来检查其是否完成计算,这个检查 适用于 tabular格式,但tablur格式实际上 可能不包含所有的结果
	my ($in_fasta_hash, $diamond_out_hash_file)=@_;
	
	my $warnMsgBody="\nIn package  DiamondWork,\tIn sub Check_diamond_done_notRealRight,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";   
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	my $diamond_work_is_already_done=0;
	
	if   (   (  defined ( $diamond_out_hash_file )  ) && ( $diamond_out_hash_file=~m/\S+/ ) && (  ( -e $diamond_out_hash_file )  )   ){  print "\n 20190105-0-0-0-0 \$diamond_out_hash_file=$diamond_out_hash_file \n";   }
	else{	 print "\n 20190105-0-0-0-1 \$diamond_out_hash_file=$diamond_out_hash_file \n"; 	return  $diamond_work_is_already_done; } # DieWork::Just_dieWork( $die_MsgHead."\n \$diamond_out_hash_file=$diamond_out_hash_file should be a defined not empty file  !!  $!\n\n\n".$caller_inform ); 	}
	my $diamond_out_hash=Storable::retrieve( $diamond_out_hash_file );
	if  (   (  defined ( $diamond_out_hash )  ) && (  ref ( $diamond_out_hash ) eq 'HASH'  )   ){ print "\n 20190105-0-0-0-2 \$diamond_out_hash=$diamond_out_hash \n"; }
	else{		print "\n 20190105-0-0-0-3 \$diamond_out_hash=$diamond_out_hash \n"; return  $diamond_work_is_already_done; } # DieWork::Just_dieWork( $die_MsgHead."\n \$diamond_out_hash=$diamond_out_hash should be a HASH  !!  $!\n\n\n".$caller_inform ); 	}
	#if  (   (  defined ( $diamond_out_hash->{'_ResultArray'} )  ) && (  ref ( $diamond_out_hash->{'_ResultArray'} ) eq 'ARRAY'  )   ){}
	#else{		return  $diamond_work_is_already_done; } # DieWork::Just_dieWork( $die_MsgHead."\n \$diamond_out_hash->{'_ResultArray'}=$diamond_out_hash->{'_ResultArray'} should be a HASH  !!  $!\n\n\n".$caller_inform ); 	}
	
	
	if  (   (  defined ( $in_fasta_hash )  ) && (  ref ( $in_fasta_hash ) eq 'HASH'  )   ){ print "\n 20190105-0-0-0-4 \$in_fasta_hash=$in_fasta_hash \n"; }
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$in_fasta_hash=$in_fasta_hash should be a HASH  !!  $!\n\n\n".$caller_inform ); 	}
	
	my $diammondOutHashQueryHash;
  #foreach my $eachQry (  @{ $diamond_out_hash->{'_ResultArray'} }  ){  print "20190105-0  \$eachQry->{'_query_name'}=$eachQry->{'_query_name'}\n";    #'_query_name'
  foreach my $eachQry (   keys (  %{ $diamond_out_hash }  )   ){  print "\n 20190105-1-0  \$eachQry=$eachQry\n";  
    $diammondOutHashQueryHash->{ $eachQry }=1 if (   (  defined ( $eachQry )  ) && ( $eachQry=~m/\S+/ )   );
    #$diammondOutHashQueryHash->{ $eachQry->{'_query_name'} }=1 if (   (  defined ( $eachQry )  ) && (  ref ( $eachQry ) eq 'HASH' ) && (  defined ( $eachQry->{'_query_name'} )  ) && ( $eachQry->{'_query_name'}=~m/\S+/ )   );
    #print "20190105-1  \$diammondOutHashQueryHash->{ \$eachQry->{'_query_name'} }=\$diammondOutHashQueryHash->{ $eachQry->{'_query_name'} }=$diammondOutHashQueryHash->{ $eachQry->{'_query_name'} }\n"; 
  }
	
	my $InFasta_simple_Hash;
  foreach my $eachFastaID (   keys (  %{ $in_fasta_hash }  )   ){  print "20190105-1-1  \$eachFastaID=$eachFastaID\n";    #'eachFastaID'
    $InFasta_simple_Hash->{ $eachFastaID }=1 if (   (  defined ( $eachFastaID )  ) && ( $eachFastaID=~m/\S+/ )   );
    #print "20190105-1  \$InFasta_simple_Hash->{ \$eachFastaID }=\$InFasta_simple_Hash->{ $eachFastaID }\n"; 
  }
	
	print "20190105-2-2-0 \$diammondOutHashQueryHash=$diammondOutHashQueryHash, \$InFasta_simple_Hash=$InFasta_simple_Hash CheckHashKey_exactly_theSAME starting... \n";
	#my $Same_Hash_or_not=ArrayHashChange::CheckHashKey_exactly_theSAME ( $diammondOutHashQueryHash, $InFasta_simple_Hash ); print "20190105-2-2-1 \$Same_Hash_or_not=$Same_Hash_or_not\n";
  my $Same_Hash_or_not=ArrayHashChange::CheckHashKey_includeRelation ( $diammondOutHashQueryHash, $InFasta_simple_Hash ); print "20190105-2-2-1 \$Same_Hash_or_not=$Same_Hash_or_not\n";
  if ( $Same_Hash_or_not == 1){
    $diamond_work_is_already_done=1;
  }
	return  $diamond_work_is_already_done;
	
}

#  DiamondWork::RunDiamond_vs_NRdb( $queryFastFil, $DiamondOUTFi, $blastProgram, $maxHitNb, $outPutFmt, $eValuSet, $identityCut );
sub RunDiamond_vs_NRdb{     # 根据 需求，进行 Diamond的运行
	
	my ( $queryFastFil, $DiamondOUTFi, $blastProgram, $maxHitNb, $outPutFmt, $eValuSet, $identityCut )=@_;
	
	print " 20190415 $queryFastFil, $DiamondOUTFi \n";
	warn  " 20190415 $queryFastFil, $DiamondOUTFi \n";
	
	my $warnMsgBody="\nIn package  DiamondWork,\tIn sub RunDiamond_vs_NRdb,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	
	my $diamondProgm='diamond';
	
	
	if (   (  defined ( $maxHitNb )  ) && ( $maxHitNb=~m/\d+/) && ( $maxHitNb > 0   )   ){}
	else {		$maxHitNb=500;  	}
	
	if (   (  defined ( $eValuSet )  ) && ( $eValuSet=~m/\d+(\.\d+)?/)   ){}
	else {		$eValuSet=0.001;  	}
	
	my $threads      ="";  #$threads      =" -p  2 ";
	my $maxTarget    ="";  $maxTarget    =" -k $maxHitNb ";
	my $eValLimit    ="";  $eValLimit    =" -e $eValuSet ";
	my $outPtFormat  ='';  $outPtFormat  =" -f $outPutFmt " if (   (  defined ( $outPutFmt )  ) && ( $outPutFmt=~m/\d+/) && ( $outPutFmt >= 0   )   );
	my $idetyCutSet  ='';  $idetyCutSet  =" --id $identityCut " if (   (  defined ( $identityCut )  ) && ( $identityCut=~m/\d+(\.\d+)?/) && ( $identityCut >= 0   )   );
	my $rangeCuligLmt='';  
	
	#if (   (  defined ( $maxHitNb )  ) && ( $maxHitNb=~m/\d+/) && ( $maxHitNb >= 0   )   ){	  $outPtFormat=" -f $outPutFmt ";	}
	
  
  my $BlastDataBas='/home/fredjiang/EightT/NR.2018.05.24/NRdimaond';
  
  if (   (  defined ( $blastProgram )  ) && ( $blastProgram=~m/\S+/)   ){
  	$blastProgram=lc $blastProgram;
  	if (  ( $blastProgram eq 'blastp' ) || ( $blastProgram eq 'blastx' )  ){
  		if ( $blastProgram eq 'blastx' ){
  			#$rangeCuligLmt=" --range-culling 100 ";
  		}
  	}
  	else{
  		DieWork::Just_dieWork( $die_MsgHead."\n \$blastProgram=$blastProgram should be blastp or blastx !!  $!\n\n\n".$caller_inform ); 
  	}  	
  } 
  else{
  	DieWork::Just_dieWork( $die_MsgHead."\n \$blastProgram=$blastProgram should not be empty or not defined !!  $!\n\n\n".$caller_inform ); 
  } 
  
  
  my $diamondCMD=$diamondProgm." ". $blastProgram." -d ".$BlastDataBas." -q ".$queryFastFil." -o ".$DiamondOUTFi."  $rangeCuligLmt $eValLimit $outPtFormat $idetyCutSet $maxTarget $threads";
  #diamond blastp -d ~/EightT/NR.2018.05.24/NRdimaond -q /home/fredjiang/EightT/fredjiang.2018.04.02/AlgaeReDo.2018.04.03/FindMoreSELENO.NOUNoC.2018.11.12/tempClustalw/PDI_e/PDI_einPsc.inseq.txt -o /home/fredjiang/EightT/fredjiang.2018.04.02/AlgaeReDo.2018.04.03/FindMoreSELENO.NOUNoC.2018.11.12/tempClustalw/PDI_e/PDI_einPsc.inseq.diamond.out.txt --threads 60 
  warn  "\n\n".$diamondCMD."\n\n";
  print "\n\n".$diamondCMD."\n\n";
  system ("$diamondCMD");
  #DirFileHandle::PrintDumper ( $opt_o, BlastHandle::BioPerlBlastPraser20170325($opt_i) );     
  
  #DiamondWork::ChangeBlastXFrame_to_withPLUS( $DiamondOUTFi, $New_DiamondOUTFi );
  
  
}

#  DiamondWork::RunDiamond_vs_ownDatabase ( $queryFastFil, $BlastDataBas, $DiamondOUTFi, $blastProgram, $maxHitNb, $outPutFmt, $eValuSet, $identityCut );
sub RunDiamond_vs_ownDatabase{     # 根据 需求，进行 Diamond的运行, the database can be set by yourslef
	
	my ( $queryFastFil, $BlastDataBas, $DiamondOUTFi, $blastProgram, $maxHitNb, $outPutFmt, $eValuSet, $identityCut )=@_;
	
	print " 20190415 $queryFastFil, $DiamondOUTFi \n";
	warn  " 20190415 $queryFastFil, $DiamondOUTFi \n";
	
	my $warnMsgBody="\nIn package  DiamondWork,\tIn sub RunDiamond_vs_NRdb,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  my $caller_inform=DirFileHandle::print_SubCallerInform;
	
  if   (   (  defined ( $queryFastFil )  ) && ( $queryFastFil=~m/\S+/ ) && (  -e ( $queryFastFil )  )   )   {}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$queryFastFil=$queryFastFil should be a defined not empty file  !!  $!\n\n\n".$caller_inform ); 	}
  
  if   (   (  defined ( $BlastDataBas )  ) && ( $BlastDataBas=~m/\S+/ ) && (  -e ( $BlastDataBas )  )   )   {}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$BlastDataBas=$BlastDataBas should be a defined not empty file  !!  $!\n\n\n".$caller_inform ); 	}
  
  if   (   (  defined ( $DiamondOUTFi )  ) && ( $DiamondOUTFi=~m/\S+/ )   )   {}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$DiamondOUTFi=$DiamondOUTFi should be a defined not empty string as output file path  !!  $!\n\n\n".$caller_inform ); 	}

  
  
	
	
	my $diamondProgm='diamond';
	
	
	if (   (  defined ( $maxHitNb )  ) && ( $maxHitNb=~m/\d+/) && ( $maxHitNb > 0   )   ){}
	else {		$maxHitNb=500;  	}
	
	if (   (  defined ( $eValuSet )  ) && ( $eValuSet=~m/\d+(\.\d+)?/)   ){}
	else {		$eValuSet=0.001;  	}
	
	my $threads      ="";  #$threads      =" -p  2 ";
	my $maxTarget    ="";  $maxTarget    =" -k $maxHitNb ";
	my $eValLimit    ="";  $eValLimit    =" -e $eValuSet ";
	my $outPtFormat  ='';  $outPtFormat  =" -f $outPutFmt " if (   (  defined ( $outPutFmt )  ) && ( $outPutFmt=~m/\d+/) && ( $outPutFmt >= 0   )   );
	my $idetyCutSet  ='';  $idetyCutSet  =" --id $identityCut " if (   (  defined ( $identityCut )  ) && ( $identityCut=~m/\d+(\.\d+)?/) && ( $identityCut >= 0   )   );
	my $rangeCuligLmt='';  
	
	#if (   (  defined ( $maxHitNb )  ) && ( $maxHitNb=~m/\d+/) && ( $maxHitNb >= 0   )   ){	  $outPtFormat=" -f $outPutFmt ";	}
	
  
  
  
  if (   (  defined ( $blastProgram )  ) && ( $blastProgram=~m/\S+/)   ){
  	$blastProgram=lc $blastProgram;
  	if (  ( $blastProgram eq 'blastp' ) || ( $blastProgram eq 'blastx' )  ){
  		if ( $blastProgram eq 'blastx' ){
  			#$rangeCuligLmt=" --range-culling 100 ";
  		}
  	}
  	else{
  		DieWork::Just_dieWork( $die_MsgHead."\n \$blastProgram=$blastProgram should be blastp or blastx !!  $!\n\n\n".$caller_inform ); 
  	}  	
  } 
  else{
  	DieWork::Just_dieWork( $die_MsgHead."\n \$blastProgram=$blastProgram should not be empty or not defined !!  $!\n\n\n".$caller_inform ); 
  } 
  
  
  my $diamondCMD=$diamondProgm." ". $blastProgram." -d ".$BlastDataBas." -q ".$queryFastFil." -o ".$DiamondOUTFi."  $rangeCuligLmt $eValLimit $outPtFormat $idetyCutSet $maxTarget $threads";
  #diamond blastp -d ~/EightT/NR.2018.05.24/NRdimaond -q /home/fredjiang/EightT/fredjiang.2018.04.02/AlgaeReDo.2018.04.03/FindMoreSELENO.NOUNoC.2018.11.12/tempClustalw/PDI_e/PDI_einPsc.inseq.txt -o /home/fredjiang/EightT/fredjiang.2018.04.02/AlgaeReDo.2018.04.03/FindMoreSELENO.NOUNoC.2018.11.12/tempClustalw/PDI_e/PDI_einPsc.inseq.diamond.out.txt --threads 60 
  warn  "\n\n".$diamondCMD."\n\n";
  print "\n\n".$diamondCMD."\n\n";
  system ("$diamondCMD");
  #DirFileHandle::PrintDumper ( $opt_o, BlastHandle::BioPerlBlastPraser20170325($opt_i) );     
  
  #DiamondWork::ChangeBlastXFrame_to_withPLUS( $DiamondOUTFi, $New_DiamondOUTFi );
  
  
}


sub ChangeBlastXFrame_to_withPLUS{  #   DiamondWork::ChangeBlastXFrame_to_withPLUS( $DiamondOUTFi , $New_DiamondOUTFi  );
	my ( $DiamondOUTFi, $New_DiamondOUTFi )=@_;  
	
	my $warnMsgBody="\nIn package  DiamondWork,\tIn sub ChangeBlastXFrame_to_withPLUS,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	if (   (  defined ( $DiamondOUTFi )  ) && ( $DiamondOUTFi=~m/\S+/ ) && (  ( -e $DiamondOUTFi )  ) && (  defined ( $New_DiamondOUTFi )  ) && ( $New_DiamondOUTFi=~m/\S+/ )   ){
		open (IN,  $DiamondOUTFi       ) or DieWork::Just_dieWork( $die_MsgHead."\n cannot open   \$DiamondOUTFi=$DiamondOUTFi : $!\n\n\n".$caller_inform ); 
		open (OUT, ">$New_DiamondOUTFi") or DieWork::Just_dieWork( $die_MsgHead."\n cannot create \$New_DiamondOUTFi=$New_DiamondOUTFi : $!\n\n\n".$caller_inform ); 
		while (<IN>){
			my $line=$_;
			$line=~s/(Frame\s+=\s+)(\d+)/$1\+$2/g;  #Frame = 2
			$line=~s/^BLASTP(\s+\d+\.)/BLASTX$1/;  
			
			print OUT $line;
		}
		close (IN);
		close (OUT);
		#my $DiamondOUTString=InFileHandle::readAllfileIntoAstring($DiamondOUTFi);
		#$DiamondOUTString=~s/(Frame\s+=\s+)(\d+)/$1\+$2/g;  #Frame = 2
		#InFileHandle::PrintStringIntoFile($DiamondOUTFi, $DiamondOUTString)
	}
	else{
		DieWork::Just_dieWork( $die_MsgHead."\n  \$DiamondOUTFi=$DiamondOUTFi and \$New_DiamondOUTFi=$New_DiamondOUTFi should be all defined not NULL string, and the \$DiamondOUTFi=$DiamondOUTFi should not be empty : $!\n\n\n".$caller_inform );
	}
	
}

#DiamondWork::Build_all_FstCk_to_Sbj_HASH($workingDir);
sub Build_all_FstCk_to_Sbj_HASH{  #从diamond的全解析结果中，找出所有的 fasta id为key，以sbj为2dkey
	my ($workingDir)=@_;
	
	my $warnMsgBody="\nIn package  DiamondWork,\tIn sub Build_all_FstCk_to_Sbj_HASH,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";  
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
  if (   (  defined ( $workingDir )  ) && ( $workingDir=~m/^\S+$/ )   ){				}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$workingDir=$workingDir should be right path to hold all the outHash and out chunk files !!  $!\n\n\n".$caller_inform );	}
	
	my $Dmw_FastaKey_Hash_file=$workingDir."/0_5_2_Dmd_FstIDKyHS.txt";
	
	if   (   (  defined ( $Dmw_FastaKey_Hash_file )  ) && ( $Dmw_FastaKey_Hash_file=~m/\S+/ )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$Dmw_FastaKey_Hash_file=$Dmw_FastaKey_Hash_file should be a defined path string  !!  $!\n\n\n".$caller_inform ); 	}
	my $Dmw_FastaKey_Hash=Storable::retrieve( $Dmw_FastaKey_Hash_file );
	if  (   (  defined ( $Dmw_FastaKey_Hash )  ) && (  ref ( $Dmw_FastaKey_Hash ) eq 'HASH'  )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$Dmw_FastaKey_Hash=$Dmw_FastaKey_Hash should be a HASH  !!  $!\n\n\n".$caller_inform ); 	}
	
	
	my $DmHashPthKey_HASH=ArrayHashChange::Change_2d_hashINTO_3d_hash ($Dmw_FastaKey_Hash, '0_2_2_DmdNROutHASH');
	
	my $orgFstCk_to_SbjHash;
	my $Sbj_to_orgFstCkHash;
	my $Sbj_to_anyCkID_Hash;
	
	foreach my $eachDmCkHashFile (    sort {$a cmp $b} (   keys (  %{ $DmHashPthKey_HASH }  )   )    ){
		
		if   (   (  defined ( $eachDmCkHashFile )  ) && ( $eachDmCkHashFile=~m/\S+/ )   ){}
	  else{		DieWork::Just_dieWork( $die_MsgHead."\n \$eachDmCkHashFile=$eachDmCkHashFile should be a defined path string  !!  $!\n\n\n".$caller_inform ); 	}
	  my $each_Dm_Ck_Hash=Storable::retrieve( $eachDmCkHashFile );
	  if  (   (  defined ( $each_Dm_Ck_Hash )  ) && (  ref ( $each_Dm_Ck_Hash ) eq 'HASH'  )   ){}
	  else{		DieWork::Just_dieWork( $die_MsgHead."\n \$each_Dm_Ck_Hash=$each_Dm_Ck_Hash should be a HASH  !!  $!\n\n\n".$caller_inform ); 	}
		
	  my $blastQueryHit_hash=BlastHandle::ChangeBlastResult_into_Qry_Hit_HASH( $each_Dm_Ck_Hash );	  
	  if  (   (  defined ( $blastQueryHit_hash )  ) && (  ref ( $blastQueryHit_hash ) eq 'HASH'  )   ){}
	  else{		DieWork::Just_dieWork( $die_MsgHead."\n \$blastQueryHit_hash=$blastQueryHit_hash should be a HASH  !!  $!\n\n\n".$caller_inform ); 	}
		
		foreach my $eachChunkID (    sort {$a cmp $b}  (   keys (  %{ $blastQueryHit_hash }  )   )    ){
			
			if (   (  defined ( $Dmw_FastaKey_Hash->{$eachChunkID}->{'0_3_1_uplvl_JLfsID'} )  ) && ( $Dmw_FastaKey_Hash->{$eachChunkID}->{'0_3_1_uplvl_JLfsID'}=~m/^\S+$/ )   ){				}
	    else{		DieWork::Just_dieWork( $die_MsgHead."\n \$Dmw_FastaKey_Hash->{$eachChunkID}->{'0_3_1_uplvl_JLfsID'}=$Dmw_FastaKey_Hash->{$eachChunkID}->{'0_3_1_uplvl_JLfsID'} should be a defined not empty string !!  $!\n\n\n".$caller_inform );	}
			my $orgFsCkID=$Dmw_FastaKey_Hash->{$eachChunkID}->{'0_3_1_uplvl_JLfsID'};
			
			if  (   (  defined ( $blastQueryHit_hash->{$eachChunkID} )  ) && (  ref ( $blastQueryHit_hash->{$eachChunkID} ) eq 'HASH'  )   ){}
	    else{		DieWork::Just_dieWork( $die_MsgHead."\n \$blastQueryHit_hash->{$eachChunkID}=$blastQueryHit_hash->{$eachChunkID} should be a HASH  !!  $!\n\n\n".$caller_inform ); 	}
			
			foreach my $eachSubjID (    sort {$a cmp $b} (    keys (  %{ $blastQueryHit_hash->{$eachChunkID} }  )   )    ){
				
				if  (   (  defined ( $blastQueryHit_hash->{$eachChunkID}->{$eachSubjID} )  ) && (  ref ( $blastQueryHit_hash->{$eachChunkID}->{$eachSubjID} ) eq 'HASH'  )   ){}
	      else{		DieWork::Just_dieWork( $die_MsgHead."\n \$blastQueryHit_hash->{$eachChunkID}->{$eachSubjID}=$blastQueryHit_hash->{$eachChunkID}->{$eachSubjID} should be a HASH  !!  $!\n\n\n".$caller_inform ); 	}
				
				$orgFstCk_to_SbjHash->{$orgFsCkID}->{$eachSubjID}=1;
				$Sbj_to_orgFstCkHash->{$eachSubjID}->{$orgFsCkID}=1;
				$Sbj_to_anyCkID_Hash->{$eachSubjID}=$eachChunkID;
				
			}
			
		}
		
  }

  if  (   (  defined ( $orgFstCk_to_SbjHash )  ) && (  ref ( $orgFstCk_to_SbjHash ) eq 'HASH'  )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$orgFstCk_to_SbjHash=$orgFstCk_to_SbjHash should be a HASH  !!  $!\n\n\n".$caller_inform ); 	}
	
	if  (   (  defined ( $Sbj_to_orgFstCkHash )  ) && (  ref ( $Sbj_to_orgFstCkHash ) eq 'HASH'  )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$Sbj_to_orgFstCkHash=$Sbj_to_orgFstCkHash should be a HASH  !!  $!\n\n\n".$caller_inform ); 	}
  
  if  (   (  defined ( $Sbj_to_anyCkID_Hash )  ) && (  ref ( $Sbj_to_anyCkID_Hash ) eq 'HASH'  )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$Sbj_to_anyCkID_Hash=$Sbj_to_anyCkID_Hash should be a HASH  !!  $!\n\n\n".$caller_inform ); 	}
  
     
  my $orgFstCk_to_SbjHash_file=$workingDir."/0_3_4_FstCKidto_sbj.HSH";   DirFileHandle::PrintDumper ( $orgFstCk_to_SbjHash_file,  $orgFstCk_to_SbjHash   );
  my $Sbj_to_orgFstCkHash_file=$workingDir."/0_3_5_sbj_toFstCKid.HSH";   DirFileHandle::PrintDumper ( $Sbj_to_orgFstCkHash_file,  $Sbj_to_orgFstCkHash   );
  my $Sbj_to_AnyCuKIDHash_file=$workingDir."/0_3_6_sbj_toAnyCkid.HSH";   DirFileHandle::PrintDumper ( $Sbj_to_orgFstCkHash_file,  $Sbj_to_anyCkID_Hash   );
  
  
  
  
  
  #my $efeachDmdSbjWkDIR=$workingDir."/0_4_0_efeachDmdSbjWkDIR";
  #my $SbjFastaString=DiamondWork::GetGeneBank_diamondTabularVer( $orgFstCk_to_SbjHash, $efeachDmdSbjWkDIR );  
  
  #my $DmdSbjprFastaStri=$workingDir."/0_4_1_DmdSbjprFasta.txt";
  
	#InFileHandle::PrintStringIntoFile($DmdSbjprFastaStri, $SbjFastaString);
	
}


#DiamondWork::Make_tblastn_qryDb_pair($workingDir);
sub Make_tblastn_qryDb_pair{  #从diamond的全解析结果中，找出所有的需要efeach genbank的subject，并做成 hash
	my ($workingDir)=@_;
	
	my $warnMsgBody="\nIn package  DiamondWork,\tIn sub Make_tblastn_qryDb_pair,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";  
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	my $chuck_hold_dir_inFileLimit=1000;
	
  if (   (  defined ( $workingDir )  ) && ( $workingDir=~m/^\S+$/ )   ){				}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$workingDir=$workingDir should be right path to hold all the outHash and out chunk files !!  $!\n\n\n".$caller_inform );	}
	
	my $in_FastaKey_Hash_file =$workingDir."/0_0_2_FastaIDKeyHSH.txt";
	DieWork::Check_FileDirExist_or_DIE          ( $in_FastaKey_Hash_file, "\$in_FastaKey_Hash_file", $die_MsgHead, $caller_inform  );
	my $in_FastaKey_Hash=Storable::retrieve( $in_FastaKey_Hash_file );	  
	DieWork::Check_Hash_or_DIE                  ( $in_FastaKey_Hash,       "\$in_FastaKey_Hash",     $die_MsgHead, $caller_inform  );	
	
	my $Dmw_FastaKey_Hash_file=$workingDir."/0_5_2_Dmd_FstIDKyHS.txt";	
	if   (   (  defined ( $Dmw_FastaKey_Hash_file )  ) && ( $Dmw_FastaKey_Hash_file=~m/\S+/ ) && (  ( -e $Dmw_FastaKey_Hash_file )  )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$Dmw_FastaKey_Hash_file=$Dmw_FastaKey_Hash_file should be a defined path string  !!  $!\n\n\n".$caller_inform ); 	}
	my $Dmw_FastaKey_Hash=Storable::retrieve( $Dmw_FastaKey_Hash_file );
	if  (   (  defined ( $Dmw_FastaKey_Hash )  ) && (  ref ( $Dmw_FastaKey_Hash ) eq 'HASH'  )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$Dmw_FastaKey_Hash=$Dmw_FastaKey_Hash should be a HASH  !!  $!\n\n\n".$caller_inform ); 	}
	
	
	
	
	my $orgFstCk_to_SbjHash_file=$workingDir."/0_3_4_FstCKidto_sbj.HSH";
	if   (   (  defined ( $orgFstCk_to_SbjHash_file )  ) && ( $orgFstCk_to_SbjHash_file=~m/\S+/ ) && (  ( -e $orgFstCk_to_SbjHash_file )  )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$orgFstCk_to_SbjHash_file=$orgFstCk_to_SbjHash_file should be a defined path string  !!  $!\n\n\n".$caller_inform ); 	}
	my $orgFstCk_to_SbjHash=Storable::retrieve( $orgFstCk_to_SbjHash_file );
	if  (   (  defined ( $orgFstCk_to_SbjHash )  ) && (  ref ( $orgFstCk_to_SbjHash ) eq 'HASH'  )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$orgFstCk_to_SbjHash=$orgFstCk_to_SbjHash should be a HASH  !!  $!\n\n\n".$caller_inform ); 	}
	
	my $Sbj_to_AnyCuKIDHash_file=$workingDir."/0_3_6_sbj_toAnyCkid.HSH";
	if   (   (  defined ( $Sbj_to_AnyCuKIDHash_file )  ) && ( $Sbj_to_AnyCuKIDHash_file=~m/\S+/ ) && (  ( -e $Sbj_to_AnyCuKIDHash_file )  )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$Sbj_to_AnyCuKIDHash_file=$Sbj_to_AnyCuKIDHash_file should be a defined path string  !!  $!\n\n\n".$caller_inform ); 	}
	my $Sbj_to_AnyCuKIDHash=Storable::retrieve( $Sbj_to_AnyCuKIDHash_file );
	if  (   (  defined ( $Sbj_to_AnyCuKIDHash )  ) && (  ref ( $Sbj_to_AnyCuKIDHash ) eq 'HASH'  )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$Sbj_to_AnyCuKIDHash=$Sbj_to_AnyCuKIDHash should be a HASH  !!  $!\n\n\n".$caller_inform ); 	}
	
	my $OrgCk_Sbjs_tblastn_wkDir=$workingDir."/0_4_0_OrgCk_Sbj_tblastn";
	my $OrgCk_Sbj_TBLN_HASH_File=$workingDir."/0_4_1_OCk_Sbj_tbltn.HSH";
	
	my @orgFstCk_id_array=(   keys (  %{ $orgFstCk_to_SbjHash }  )   );
	my $howManyOrgFstCkID=@orgFstCk_id_array;
	
	my $Dir_nb=1; my $dirTotalNB=int ( $howManyOrgFstCkID/$chuck_hold_dir_inFileLimit ); my $Dir_spt_NB= MatrixCsvChange::SprintfKeyHead( $dirTotalNB+1, $Dir_nb );
	system (" mkdir -p $OrgCk_Sbjs_tblastn_wkDir/$Dir_spt_NB");
	my $inDirNb=0;
	
	my $tblastnWrDirNB=1;
	
	my $outPutTBLASTN_infmHASH;
	
	foreach my $eachFastCk_JLID (    sort {$a cmp $b} ( @orgFstCk_id_array )   ){
		
		#if (   (  defined ( $orgFstCk_to_SbjHash->{$eachFastCk_JLID} )  ) && (  ref ( $orgFstCk_to_SbjHash->{$eachFastCk_JLID} ) eq 'HASH'  )   ){ 	}
    #else{		DieWork::Just_dieWork( $die_MsgHead."\n \$orgFstCk_to_SbjHash->{$eachFastCk_JLID}=$orgFstCk_to_SbjHash->{$eachFastCk_JLID} should be a HASH ref !!  $!\n\n\n".$caller_inform );	}
		
		DieWork::Check_Hash_or_DIE( $orgFstCk_to_SbjHash->{$eachFastCk_JLID}, "\$orgFstCk_to_SbjHash->{$eachFastCk_JLID}", $die_MsgHead, $caller_inform  );
		
		my $Sptf_tblastnWrDirNB= MatrixCsvChange::SprintfKeyHead( $howManyOrgFstCkID+1, $tblastnWrDirNB );
		
		my $tblastnWorkIngDir="$OrgCk_Sbjs_tblastn_wkDir/$Dir_spt_NB/$Sptf_tblastnWrDirNB";
		system (" mkdir -p $OrgCk_Sbjs_tblastn_wkDir/$Dir_spt_NB/$Sptf_tblastnWrDirNB");
		
		$outPutTBLASTN_infmHASH->{$eachFastCk_JLID}=Storable::dclone( $in_FastaKey_Hash->{$eachFastCk_JLID} );
		$outPutTBLASTN_infmHASH->{$eachFastCk_JLID}->{'0_4_1_tBlastnWkDIR'}=$tblastnWorkIngDir;
		$outPutTBLASTN_infmHASH->{$eachFastCk_JLID}->{'0_4_2_DNA_database'}=$tblastnWorkIngDir."/0_0_0_0_DNA_database";
		$outPutTBLASTN_infmHASH->{$eachFastCk_JLID}->{'0_4_3_DNA_db_index'}=$tblastnWorkIngDir."/0_0_1_0_DNA_db_index";
		$outPutTBLASTN_infmHASH->{$eachFastCk_JLID}->{'0_4_4_ProFasta_TXT'}=$tblastnWorkIngDir."/0_0_2_0_ProFasta.txt";
		$outPutTBLASTN_infmHASH->{$eachFastCk_JLID}->{'0_4_5_ProFasta_IDX'}=$tblastnWorkIngDir."/0_0_3_0_ProFst_index";
		
		DieWork::Check_FileDirExist_or_DIE( $outPutTBLASTN_infmHASH->{$eachFastCk_JLID}->{'0_1_4_chkFile_indx'}, "\$outPutTBLASTN_infmHASH->{$eachFastCk_JLID}->{'0_1_4_chkFile_indx'}", $die_MsgHead, $caller_inform  );	
		my $DNAckIndx=$outPutTBLASTN_infmHASH->{$eachFastCk_JLID}->{'0_1_4_chkFile_indx'};
		my $DNAstring=FastaFileHandle::feach_seqString_from_idx_File($eachFastCk_JLID, $DNAckIndx);
		
		FastaFileHandle::BuildFastaFile_with_pure_string ($outPutTBLASTN_infmHASH->{$eachFastCk_JLID}->{'0_4_2_DNA_database'}, $DNAstring, $eachFastCk_JLID);	
		FastaFileHandle::BuildIdxFile_for_fastaFile      ($outPutTBLASTN_infmHASH->{$eachFastCk_JLID}->{'0_4_2_DNA_database'}, $outPutTBLASTN_infmHASH->{$eachFastCk_JLID}->{'0_4_3_DNA_db_index'} );
		
		my $SbjString='';
		foreach my $eachSbj_ID (    sort {$a cmp $b}  (   keys (  %{ $orgFstCk_to_SbjHash->{$eachFastCk_JLID} }  )   )    ){
		  
		  #if   (   (  defined ( $orgFstCk_to_SbjHash->{$eachFastCk_JLID}->{$eachSbj_ID} )  ) && ( $orgFstCk_to_SbjHash->{$eachFastCk_JLID}->{$eachSbj_ID} == 1 )   ){}
	    #else{		DieWork::Just_dieWork( $die_MsgHead."\n \$orgFstCk_to_SbjHash->{$eachFastCk_JLID}->{$eachSbj_ID}=$orgFstCk_to_SbjHash->{$eachFastCk_JLID}->{$eachSbj_ID} should be a defined value 1  !!  $!\n\n\n".$caller_inform ); 	}
		  DieWork::Check_INTNB_equal_to_aNUMBER_or_DIE( $orgFstCk_to_SbjHash->{$eachFastCk_JLID}->{$eachSbj_ID}, 1, "\$orgFstCk_to_SbjHash->{$eachFastCk_JLID}->{$eachSbj_ID}", $die_MsgHead, $caller_inform  );
		  DieWork::Check_DfdNoEmptString_or_DIE       ( $Sbj_to_AnyCuKIDHash->{$eachSbj_ID},                        "\$Sbj_to_AnyCuKIDHash->{$eachSbj_ID}",                     $die_MsgHead, $caller_inform  );
		  my $anyCkID=$Sbj_to_AnyCuKIDHash->{$eachSbj_ID};
		  DieWork::Check_FileDirExist_or_DIE          ( $Dmw_FastaKey_Hash->{$anyCkID}->{'0_2_4_DmdOSjFstIdx'},     "\$Dmw_FastaKey_Hash->{$anyCkID}->{'0_2_4_DmdOSjFstIdx'}",  $die_MsgHead, $caller_inform  );
		  my $PROstring=FastaFileHandle::feach_seqString_from_idx_File($eachSbj_ID, $Dmw_FastaKey_Hash->{$anyCkID}->{'0_2_4_DmdOSjFstIdx'});
		  $SbjString.=">".$eachSbj_ID."\n".$PROstring."\n\n";
		
		}
		InFileHandle::PrintStringIntoFile($outPutTBLASTN_infmHASH->{$eachFastCk_JLID}->{'0_4_4_ProFasta_TXT'}, $SbjString);
		FastaFileHandle::BuildIdxFile_for_fastaFile($outPutTBLASTN_infmHASH->{$eachFastCk_JLID}->{'0_4_4_ProFasta_TXT'}, $outPutTBLASTN_infmHASH->{$eachFastCk_JLID}->{'0_4_5_ProFasta_IDX'});
		
		
		
		
		
		
		if (  $inDirNb >= ( $chuck_hold_dir_inFileLimit - 1 )  ){
			$inDirNb=0;
			$Dir_nb++;
			$Dir_spt_NB= MatrixCsvChange::SprintfKeyHead( $dirTotalNB+1, $Dir_nb );
			system (" mkdir -p $OrgCk_Sbjs_tblastn_wkDir/$Dir_spt_NB");
		}
		else{
		  $inDirNb++;	
		}
		
		$tblastnWrDirNB++;	
	}
	
	
	DieWork::Check_Hash_or_DIE( $outPutTBLASTN_infmHASH, "\$outPutTBLASTN_infmHASH", $die_MsgHead, $caller_inform  );
	
	DirFileHandle::PrintDumper ( $OrgCk_Sbj_TBLN_HASH_File,  $outPutTBLASTN_infmHASH   );
	
	return $outPutTBLASTN_infmHASH;
	
	
}



#DiamondWork::Efeach_all_diamond_subject_tabularV($workingDir);
sub Efeach_all_diamond_subject_tabularV{  #从diamond的tabular结果中，找出所有的需要efeach genbank的subject，并efeach，并更新本地gbk的database和idx。
	my ($workingDir)=@_;
	
	my $warnMsgBody="\nIn package  DiamondWork,\tIn sub Efeach_all_diamond_subject_tabularV,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";  
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
  if (   (  defined ( $workingDir )  ) && ( $workingDir=~m/^\S+$/ )   ){				}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$workingDir=$workingDir should be right path to hold all the outHash and out chunk files !!  $!\n\n\n".$caller_inform );	}
	
	my $Dmw_FastaKey_Hash_file=$workingDir."/0_5_2_Dmd_FstIDKyHS.txt";
	
	if   (   (  defined ( $Dmw_FastaKey_Hash_file )  ) && ( $Dmw_FastaKey_Hash_file=~m/\S+/ )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$Dmw_FastaKey_Hash_file=$Dmw_FastaKey_Hash_file should be a defined path string  !!  $!\n\n\n".$caller_inform ); 	}
	my $Dmw_FastaKey_Hash=Storable::retrieve( $Dmw_FastaKey_Hash_file );
	if  (   (  defined ( $Dmw_FastaKey_Hash )  ) && (  ref ( $Dmw_FastaKey_Hash ) eq 'HASH'  )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$Dmw_FastaKey_Hash=$Dmw_FastaKey_Hash should be a HASH  !!  $!\n\n\n".$caller_inform ); 	}
	
	
	my $DmHashPthKey_HASH=ArrayHashChange::Change_2d_hashINTO_3d_hash ($Dmw_FastaKey_Hash, '0_2_2_DmdNROutHASH');
	
	my $orgFstCk_to_SbjHash;
	my $Sbj_to_orgFstCkHash;
	foreach my $eachDmCkHashFile (    sort {$a cmp $b} (   keys (  %{ $DmHashPthKey_HASH }  )   )    ){
		
		if   (   (  defined ( $eachDmCkHashFile )  ) && ( $eachDmCkHashFile=~m/\S+/ )   ){}
	  else{		DieWork::Just_dieWork( $die_MsgHead."\n \$eachDmCkHashFile=$eachDmCkHashFile should be a defined path string  !!  $!\n\n\n".$caller_inform ); 	}
	  my $each_Dm_Ck_Hash=Storable::retrieve( $eachDmCkHashFile );
	  if  (   (  defined ( $each_Dm_Ck_Hash )  ) && (  ref ( $each_Dm_Ck_Hash ) eq 'HASH'  )   ){}
	  else{		DieWork::Just_dieWork( $die_MsgHead."\n \$each_Dm_Ck_Hash=$each_Dm_Ck_Hash should be a HASH  !!  $!\n\n\n".$caller_inform ); 	}
		
		foreach my $eachChunkID (    sort {$a cmp $b}  (   keys (  %{ $each_Dm_Ck_Hash }  )   )    ){
			
			if (   (  defined ( $Dmw_FastaKey_Hash->{$eachChunkID}->{'0_3_1_uplvl_JLfsID'} )  ) && ( $Dmw_FastaKey_Hash->{$eachChunkID}->{'0_3_1_uplvl_JLfsID'}=~m/^\S+$/ )   ){				}
	    else{		DieWork::Just_dieWork( $die_MsgHead."\n \$Dmw_FastaKey_Hash->{$eachChunkID}->{'0_3_1_uplvl_JLfsID'}=$Dmw_FastaKey_Hash->{$eachChunkID}->{'0_3_1_uplvl_JLfsID'} should be a defined not empty string !!  $!\n\n\n".$caller_inform );	}
			my $orgFsCkID=$Dmw_FastaKey_Hash->{$eachChunkID}->{'0_3_1_uplvl_JLfsID'};
			
			if  (   (  defined ( $each_Dm_Ck_Hash->{$eachChunkID} )  ) && (  ref ( $each_Dm_Ck_Hash->{$eachChunkID} ) eq 'HASH'  )   ){}
	    else{		DieWork::Just_dieWork( $die_MsgHead."\n \$each_Dm_Ck_Hash->{$eachChunkID}=$each_Dm_Ck_Hash->{$eachChunkID} should be a HASH  !!  $!\n\n\n".$caller_inform ); 	}
			
			foreach my $eachSubjID (    sort {$a cmp $b} (    keys (  %{ $each_Dm_Ck_Hash->{$eachChunkID} }  )   )    ){
				
				if  (   (  defined ( $each_Dm_Ck_Hash->{$eachChunkID}->{$eachSubjID} )  ) && (  ref ( $each_Dm_Ck_Hash->{$eachChunkID}->{$eachSubjID} ) eq 'ARRAY'  )   ){}
	      else{		DieWork::Just_dieWork( $die_MsgHead."\n \$each_Dm_Ck_Hash->{$eachChunkID}->{$eachSubjID}=$each_Dm_Ck_Hash->{$eachChunkID}->{$eachSubjID} should be a ARRAY  !!  $!\n\n\n".$caller_inform ); 	}
				
				$orgFstCk_to_SbjHash->{$orgFsCkID}->{$eachSubjID}=1;
				$Sbj_to_orgFstCkHash->{$eachSubjID}->{$orgFsCkID}=1;
				
			}
			
		}
		
  }

  if  (   (  defined ( $orgFstCk_to_SbjHash )  ) && (  ref ( $orgFstCk_to_SbjHash ) eq 'HASH'  )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$orgFstCk_to_SbjHash=$orgFstCk_to_SbjHash should be a HASH  !!  $!\n\n\n".$caller_inform ); 	}
	
	if  (   (  defined ( $Sbj_to_orgFstCkHash )  ) && (  ref ( $Sbj_to_orgFstCkHash ) eq 'HASH'  )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$Sbj_to_orgFstCkHash=$Sbj_to_orgFstCkHash should be a HASH  !!  $!\n\n\n".$caller_inform ); 	}
     
  my $orgFstCk_to_SbjHash_file=$workingDir."/0_3_4_FstCKidto_sbj.HSH";   DirFileHandle::PrintDumper ( $orgFstCk_to_SbjHash_file,  $orgFstCk_to_SbjHash   );
  my $Sbj_to_orgFstCkHash_file=$workingDir."/0_3_5_sbj_toFstCKid.HSH";   DirFileHandle::PrintDumper ( $Sbj_to_orgFstCkHash_file,  $Sbj_to_orgFstCkHash   );
  
  
  my $efeachDmdSbjWkDIR=$workingDir."/0_4_0_efeachDmdSbjWkDIR";
  my $SbjFastaString=DiamondWork::GetGeneBank_diamondTabularVer( $orgFstCk_to_SbjHash, $efeachDmdSbjWkDIR );  
  
  my $DmdSbjprFastaStri=$workingDir."/0_4_1_DmdSbjprFasta.txt";
  
	InFileHandle::PrintStringIntoFile($DmdSbjprFastaStri, $SbjFastaString);
	
}


# my $fastaOut=DiamondWork::GetGeneBank_for_diamondOut( $DiaMondOutHash, $workingDir, $top_howMany );
sub GetGeneBank_for_diamondOut{ #从 blastprase的hash中，获得diamond结果中所有的sbject的信息，并efeach其genbank文件，并输出fasta文件
	my ( $DiaMondOutHash, $workingDir, $top_howMany )=@_;
	
	my $warnMsgBody="\nIn package  DiamondWork,\tIn sub GetGeneBank_for_diamondOut,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	if  (   (  defined ( $workingDir )  ) && ( $workingDir=~m/\S+/ )   ){
		system ( "mkdir -p $workingDir");
	}
	
	my $whlBlsArHsh; 
  my $EfcProArHsh;
  my $Lcs2accHash;
  my $Acc2LcsHash;
  my $GBfilePthHS;
  my $TPGBflPthHS;
  
  my $LcsMaxIdtHS;
	
	my $GenBankFDir_0_2_a=      $workingDir."/".'3.1.1.GenBankFDir';  system ( "mkdir -p $GenBankFDir_0_2_a");
	my $GBfileTpDir_0_2_b=      $workingDir."/".'3.1.2.GBfileTpDir';  system ( "mkdir -p $GBfileTpDir_0_2_b");
	my $TPGBflPthHS_0_2_b_0=    $workingDir."/".'2.1.1.TPGBflPthHS.txt';
	my $LcsMaxIdtHS_0_1_a=      $workingDir."/".'4.1.1.LcsMaxIdtHS';
	
	if  (   (  defined ( $DiaMondOutHash )  ) && (  ref ( $DiaMondOutHash ) eq 'HASH'  ) && (  defined ( $DiaMondOutHash->{'_ResultArray'} )  ) && (  ref ( $DiaMondOutHash->{'_ResultArray'} ) eq 'ARRAY'  )   ){      	
        	
    #system ("mkdir -p $GenBankFDir_0_2_a");  
    #system ("mkdir -p $GBfileTpDir_0_2_b");  
    
   
    my $query_nb=1;
    foreach my $eachQrHASH (  @{ $DiaMondOutHash->{'_ResultArray'} }  ){                                                                                                                                        
    	                                                                                                                                                                                                          
    	if  (   (  defined ( $eachQrHASH )  ) && (  ref ( $eachQrHASH ) eq 'HASH'  ) && (  defined ( $eachQrHASH->{'_hitArray'} )  ) && (  ref ( $eachQrHASH->{'_hitArray'} ) eq 'ARRAY'  )   ){                  
    	  print $eachQrHASH->{'_query_name'}."\n";                                                                                                                                                                
    	  
    	  my $Qurey_name=$eachQrHASH->{'_query_name'};
    	  my $Lcs_to_acc_HASH;  my $Acc_to_lcs_HASH;
    	  my $wholeProtArray;                                                                                                                                                                               
    	  my $EfcProteinArray;       my $pt_to_fatcch_NB=0;                                                                                                                                                                                
    	  foreach my $eachHitHASH (  @{ $eachQrHASH->{'_hitArray'} }  ){                                                                                                                                          
    	  	                                                                                                                                                                                                      
    	  	if  (   (  defined ( $eachHitHASH )  ) && (  ref ( $eachHitHASH ) eq 'HASH'  )   ){                                                                                                                   
    	  		#$ProteinsHASH->{ $eachHitHASH->{'_accessionNB'} }=1;                            
    	  		
    	  		
    	  		my $lcsID=$eachHitHASH->{'_accessionNB'};
    	  		my $accID=$eachHitHASH->{'0_0_2_giNub'};
    	  		my $htIdt=$eachHitHASH->{'_hitTotalIdentical'};
    	  		if (    (  defined ( $LcsMaxIdtHS->{$lcsID} )  ) && (  defined ( $LcsMaxIdtHS->{$lcsID}->{'1_maxIdt'} )  ) && ( $LcsMaxIdtHS->{$lcsID}->{'1_maxIdt'}=~/\S+/)   ){
    	  			if ( $htIdt > $LcsMaxIdtHS->{$lcsID}->{'1_maxIdt'} ){
    	  				$LcsMaxIdtHS->{$lcsID}->{'0__Query'}=$Qurey_name;  print "\n20181227-1-1 \$LcsMaxIdtHS->{$lcsID}->{'0__Query'}=$Qurey_name=$LcsMaxIdtHS->{$lcsID}->{'0__Query'}=$Qurey_name\n ";
    	  				$LcsMaxIdtHS->{$lcsID}->{'1_maxIdt'}=$htIdt;       print "\n20181227-1-2 \$LcsMaxIdtHS->{$lcsID}->{'1_maxIdt'}=$htIdt=$LcsMaxIdtHS->{$lcsID}->{'1_maxIdt'}=$htIdt\n ";
    	  				#$LcsMaxIdtHS->{$lcsID}->{'2_QryTax'}=$addifProteinID_to_StfcNM_HASH->{$Qurey_name} if (   (  defined ( $addifProteinID_to_StfcNM_HASH->{$Qurey_name} )  ) && ( $addifProteinID_to_StfcNM_HASH->{$Qurey_name}=~m/\S+/ )   );
    	  				#                                                  print "\n20181227-1-3 \$LcsMaxIdtHS->{$lcsID}->{'2_QryTax'}=$LcsMaxIdtHS->{$lcsID}->{'2_QryTax'}\n ";
    	  			}
    	  		}
    	  		else {
    	  			$LcsMaxIdtHS->{$lcsID}->{'0__Query'}=$Qurey_name;    print "\n20181227-2-1 \$LcsMaxIdtHS->{$lcsID}->{'0__Query'}=$Qurey_name=$LcsMaxIdtHS->{$lcsID}->{'0__Query'}=$Qurey_name\n ";
    	  			$LcsMaxIdtHS->{$lcsID}->{'1_maxIdt'}=$htIdt;         print "\n20181227-2-2 \$LcsMaxIdtHS->{$lcsID}->{'1_maxIdt'}=$htIdt=$LcsMaxIdtHS->{$lcsID}->{'1_maxIdt'}=$htIdt\n ";
    	  			#$LcsMaxIdtHS->{$lcsID}->{'2_QryTax'}=$addifProteinID_to_StfcNM_HASH->{$Qurey_name} if (   (  defined ( $addifProteinID_to_StfcNM_HASH->{$Qurey_name} )  ) && ( $addifProteinID_to_StfcNM_HASH->{$Qurey_name}=~m/\S+/ )   );
    	  		  #                                                     #print "\n20181227-2-3 \$LcsMaxIdtHS->{$lcsID}->{'2_QryTax'}=$LcsMaxIdtHS->{$lcsID}->{'2_QryTax'}\n ";
    	  		}
    	  		
    	  		
    	  		
    	  		$Lcs_to_acc_HASH->{$lcsID}=$accID;  
    	  		
    	      #$Acc_to_lcs_HASH->{$accID}=$lcsID;
    	      
    	  		push @{ $wholeProtArray}, $lcsID;                                                                        #my $msg1= "20181214-0-1\$accID=$accID      \$GBfileTpDir_0_2_b=$GBfileTpDir_0_2_b \n"; print $msg1; warn $msg1;
    	  		#my $found_in_local_db=0; $found_in_local_db=GenBankHandle::Check_present_in_JL_ProteinGB_DB ($lcsID);    #my $msg2= "20181214-0-2\$found_in_local_db=$found_in_local_db\n"; print $msg2; warn $msg2;
            #if ( $found_in_local_db == 0 ){                                                                           my $msg3= "20181214-0-3\$lcsID=$lcsID      \$GBfileTpDir_0_2_b=$GBfileTpDir_0_2_b \n"; print $msg3; warn $msg3;
            #	push @{ $EfcProteinArray}, $accID;   $pt_to_fatcch_NB++;    
            #}      		  		                                                                                                                   
    	  		                                                                                                                                      
    	  	}                                                                                                                                                                                                     
    	  }
    	  
    	  my $correctHASH;
    	  my $TParH=GenBankHandle::CorectTheLocusName_for_BlastHitAccName ($wholeProtArray);
    	  if (   (  defined ( $TParH )  ) && (  ref ( $TParH ) eq 'ARRAY'  )   ) {  ( $wholeProtArray, $correctHASH )=@{ $TParH };  }  
    	  else { my $dieMsgg="\nDIE!!!!! something wrong with  GenBankHandle::CorectTheLocusName_for_BlastHitAccName \n\n"; print $dieMsgg; die $dieMsgg; }
    	  
    	  if (   (  defined ( $correctHASH )  ) && (  ref ( $correctHASH ) eq 'HASH'  )   ) {
    	  	foreach my $orgLcs (    sort { $a cmp $b } (   keys (  %{ $correctHASH } )   )    ){
    	  	  my $newLcs=$correctHASH->{$orgLcs};
    	  	  $Lcs_to_acc_HASH->{$newLcs}=$Lcs_to_acc_HASH->{$orgLcs};
    	  	  delete ( $Lcs_to_acc_HASH->{$orgLcs} );
    	  	}
    	  	
    	  }
    	  
    	  if (   (  defined ( $wholeProtArray )  ) && (  ref ( $wholeProtArray ) eq 'ARRAY'  )   ) {
    	    foreach my $ProtLcs (  @{ $wholeProtArray }  ){ 
    	    	my $accID=$Lcs_to_acc_HASH->{$ProtLcs};
    	      $Acc_to_lcs_HASH->{ $accID }=$ProtLcs;
    	      my $found_in_local_db=0; 
    	      $found_in_local_db=GenBankHandle::Check_present_in_JL_ProteinGB_DB ($ProtLcs);    #my $msg2= "20181214-0-2\$found_in_local_db=$found_in_local_db\n"; print $msg2; warn $msg2;
            if ( $found_in_local_db == 0 ){    my $msg3= "20181214-0-3\$ProtLcs=$ProtLcs \n"; print $msg3; warn $msg3;
            	#push @{ $EfcProteinArray}, $accID;   $pt_to_fatcch_NB++;
            	$EfcProArHsh->{$accID}=1;    
            }  
    	    }
    	  }
    	  
    	  if (   (  defined ( $wholeProtArray )  ) && (  ref ( $wholeProtArray ) eq 'ARRAY'  )   ){
    	  	$whlBlsArHsh->{$Qurey_name}=$wholeProtArray;      		  	
    	  }
    	  #if (   (  defined ( $EfcProteinArray )  ) && (  ref ( $EfcProteinArray ) eq 'ARRAY'  )   ){
    	  #	$EfcProArHsh->{$Qurey_name}=$EfcProteinArray;      		  	
    	  #}
    	  if (   (  defined ( $Lcs_to_acc_HASH )  ) && (  ref ( $Lcs_to_acc_HASH ) eq 'HASH'  ) && (  defined ( $Acc_to_lcs_HASH )  ) && (  ref ( $Acc_to_lcs_HASH ) eq 'HASH'  )   ){
    	  	#$Lcs2accHash->{$Qurey_name}=$Lcs_to_acc_HASH; #my $bigHash; $bigHash=ArrayHashChange::PushSmallHash_into_bigHash($bigHash, $smlHash);
    	  	#$Acc2LcsHash->{$Qurey_name}=$Acc_to_lcs_HASH;
    	  	$Lcs2accHash=ArrayHashChange::PushSmallHash_into_bigHash($Lcs2accHash, $Lcs_to_acc_HASH);
    	  	$Acc2LcsHash=ArrayHashChange::PushSmallHash_into_bigHash($Acc2LcsHash, $Acc_to_lcs_HASH);
    	  }      		  
    	  my $TempQrGBFile=  $TPGBflPthHS->{$Qurey_name}=  $GBfileTpDir_0_2_b."/".$query_nb;        		        		  
    	  my $EachQrGBFile=  $GBfilePthHS->{$Qurey_name}=  $GenBankFDir_0_2_a."/".$query_nb;
    	  
          		        		  
    	  
    	  $query_nb++;                                                                                                                                                                               
    	                                                                                              
    	}                                                                                                                                                                                                         
    } 
	
	  DirFileHandle::PrintDumper ( $workingDir."/".'1.1.1.LcsMaxIdtHS.txt',  $LcsMaxIdtHS   ) if (   (  defined ( $LcsMaxIdtHS  )  ) && (  ref ( $LcsMaxIdtHS  ) eq 'HASH'  )    ); 
        	
        	
    DirFileHandle::PrintDumper ( $workingDir."/".'1.1.2.GBfilePthHS.txt',  $GBfilePthHS   ) if (   (  defined ( $GBfilePthHS  )  ) && (  ref ( $GBfilePthHS  ) eq 'HASH'  )    ); 
    DirFileHandle::PrintDumper ( $workingDir."/".'1.1.3.TPGBflPthHS.txt',  $TPGBflPthHS   ) if (   (  defined ( $TPGBflPthHS  )  ) && (  ref ( $TPGBflPthHS  ) eq 'HASH'  )    ); 
        	                       
    DirFileHandle::PrintDumper ( $workingDir."/".'1.1.4.Lcs2accHash.txt',  $Lcs2accHash   ) if (   (  defined ( $Lcs2accHash  )  ) && (  ref ( $Lcs2accHash  ) eq 'HASH'  )    ); 
    DirFileHandle::PrintDumper ( $workingDir."/".'1.1.5.Acc2LcsHash.txt',  $Acc2LcsHash   ) if (   (  defined ( $Acc2LcsHash  )  ) && (  ref ( $Acc2LcsHash  ) eq 'HASH'  )    ); 
        	                                         
    DirFileHandle::PrintDumper ( $workingDir."/".'1.1.6.whlBlsArHsh.txt',  $whlBlsArHsh   ) if (   (  defined ( $whlBlsArHsh  )  ) && (  ref ( $whlBlsArHsh  ) eq 'HASH'  )    ); 
    DirFileHandle::PrintDumper ( $workingDir."/".'1.1.7.EfcProArHsh.txt',  $EfcProArHsh   ) if (   (  defined ( $EfcProArHsh  )  ) && (  ref ( $EfcProArHsh  ) eq 'HASH'  )    ); 
	}
	
	
      
  if  (   (  defined ( $EfcProArHsh )  ) && (  ref ( $EfcProArHsh ) eq 'HASH'  )   ){
    my $AllEfcProArray;  my $Total_Qurey_number=0;
    foreach my $EfcProIDkey (    sort { $a cmp $b } (   keys (  %{ $EfcProArHsh } )   )    ){
    	push @{ $AllEfcProArray }, $EfcProIDkey; $Total_Qurey_number++;
    }
    if (   (  defined ( $TPGBflPthHS_0_2_b_0 )  ) && ( $TPGBflPthHS_0_2_b_0=~m/\S+/ ) && (  defined ( $AllEfcProArray )  ) && (  ref ( $AllEfcProArray ) eq 'ARRAY'  )   ){
    	
    
      my $time_3=`date`; chomp $time_3;
      my $msg4= "Now $time_3 , \ndo the EUtilitieswork :: efatch_proteinGB( \$AllEfcProArra=$AllEfcProArray \$TPGBflPthHS_0_2_b_0=\n\n$TPGBflPthHS_0_2_b_0 )\n \n\n\$Total_Qurey_number=$Total_Qurey_number protein to fatch!!\n\n";  print $msg4; warn $msg4;
      #下载数据
      EUtilitieswork::efatch_proteinGB( $AllEfcProArray, $TPGBflPthHS_0_2_b_0 );  
      
      #更新数据库
      GenBankHandle::Upgrade_JL_proteinGB_and_Acc_to_Lcs_HASHs_onlyOnePid ($TPGBflPthHS_0_2_b_0, $AllEfcProArray, $Acc2LcsHash);      
      
    }
    
  
  }

	
  
    
  my $LcsMaxIdtHS_Changed=0;
  
  ###################
  #if (   (  defined ( $whlBlsArHsh_0_2_e )  ) && (  -e ( $whlBlsArHsh_0_2_e )  )   ){
    #my $whlBlsArHsh=retrieve ($whlBlsArHsh_0_2_e);
  if  (   (  defined ( $whlBlsArHsh )  ) && (  ref ( $whlBlsArHsh ) eq 'HASH'  )   ){
  	foreach my $Qurey_name (    sort { $a cmp $b } (   keys (  %{ $whlBlsArHsh } )   )    ){
  		my $wholeProtArray=$whlBlsArHsh->{$Qurey_name};
  		if  (   (  defined ( $wholeProtArray )  ) && (  ref ( $wholeProtArray ) eq 'ARRAY'  )   ){
        
        my $correctHASH;
	      my $TParH=GenBankHandle::CorectTheLocusName_for_BlastHitAccName ($wholeProtArray);
	      if (   (  defined ( $TParH )  ) && (  ref ( $TParH ) eq 'ARRAY'  )   ) {  ( $wholeProtArray, $correctHASH )=@{ $TParH };  }
	      
	      #correct the lcs name of $LcsMaxIdtHS
	      if (   (  defined ( $correctHASH )  ) && (  ref ( $correctHASH ) eq 'HASH'  )   ) {
	      	foreach my $orgLcsNm (    sort { $a cmp $b } (   keys (  %{ $correctHASH } )   )    ){  #LcsMaxIdtHS->{$lcsID}->{'1_maxIdt'} 
	      		if  (      (  defined ( $LcsMaxIdtHS )  ) && (  ref ( $LcsMaxIdtHS ) eq 'HASH'  ) && (  defined ( $LcsMaxIdtHS->{$orgLcsNm} )  ) && (  ref ( $LcsMaxIdtHS->{$orgLcsNm} ) eq 'HASH'  ) 
	      		        && (  defined ( $LcsMaxIdtHS->{$orgLcsNm}->{'1_maxIdt'} )  ) && (  $LcsMaxIdtHS->{$orgLcsNm}->{'1_maxIdt'}=~m/\S+/ )    
	      		    )
	      		{
	      			my $newLcsNm=$correctHASH->{$orgLcsNm};
	      			$LcsMaxIdtHS->{$newLcsNm}=Storable::dclone ( $LcsMaxIdtHS->{$orgLcsNm} );
	      			delete ( $LcsMaxIdtHS->{$orgLcsNm} );
	      			$LcsMaxIdtHS_Changed=1;
	      		}
	      	}
	      }
	      
	      
	      if  (   (  defined ( $GBfilePthHS )  ) && (  ref ( $GBfilePthHS ) eq 'HASH'  ) && (  defined ( $GBfilePthHS->{$Qurey_name} )  ) && ( $GBfilePthHS->{$Qurey_name}=~/\S+/ )   ){
  		    my $EachQrGBFile=  $GBfilePthHS->{$Qurey_name};
	        GenBankHandle::Build_GB_file_from_IDarray_JL_GBidxFILE ( $wholeProtArray, $EachQrGBFile );
	      }
	      
	    }
	  }
	}
	
  ###################
  
  if (  $LcsMaxIdtHS_Changed==1 ){
  	DirFileHandle::PrintDumper ( $LcsMaxIdtHS_0_1_a,  $LcsMaxIdtHS   ) if (   (  defined ( $LcsMaxIdtHS  )  ) && (  ref ( $LcsMaxIdtHS  ) eq 'HASH'  )    );         	
  }
  
  
  #print "\n201901042318-0 \$selProFamName=$selProFamName\n";
  my $GbkEfachWkOutARRAY=GenBankHandle::BuildTreesHash_from_geneBankFiles($GenBankFDir_0_2_a, $top_howMany);
  #print "\n201901042318-1 \$selProFamName=$selProFamName\n";
  
  my ($fastaOut, $outPepInfHASH,  $TreeWord, $TreeHash, $AllGBHash, $reapeatHash, $PtHhomoHASH);
  if (   (  defined ( $GbkEfachWkOutARRAY )  ) && (  ref ( $GbkEfachWkOutARRAY ) eq 'ARRAY'  )   ) {
    ($fastaOut, $outPepInfHASH,  $TreeWord, $TreeHash, $AllGBHash, $reapeatHash, $PtHhomoHASH)=@{ $GbkEfachWkOutARRAY };   
    #$FastaSeqString.="\n\n";
    #$FastaSeqString.=$fastaOut;
    #open (OUT1, ">$allFmlProte_0_3_0"  ) or die "cannot create \$allFmlProte_0_3_0=$allFmlProte_0_3_0 : $! \n\n"; print OUT1 $FastaSeqString; close (OUT1);                                                               
    open (OUT1, ">$workingDir/3.2.1.TreeShowTxt.txt"  ) or DieWork::Just_dieWork  ("cannot create $workingDir/3.2.1.TreeShowTxt.txt : $! \n\n" ); print OUT1 $TreeWord;       close (OUT1);                                                               
                                                                                  
    DirFileHandle::PrintDumper ( $workingDir."/".'3.2.2.TreeIfmhash.txt',  $TreeHash       ) if (   (  defined ( $TreeHash      )  ) && (  ref ( $TreeHash      ) eq 'HASH'  )    );                                                                                                                                              
    DirFileHandle::PrintDumper ( $workingDir."/".'3.2.3.allGebkHasH.txt',  $AllGBHash      ) if (   (  defined ( $AllGBHash     )  ) && (  ref ( $AllGBHash     ) eq 'HASH'  )    );  
    DirFileHandle::PrintDumper ( $workingDir."/".'3.2.4.reapeatRort.txt',  $reapeatHash    ) if (   (  defined ( $reapeatHash   )  ) && (  ref ( $reapeatHash   ) eq 'HASH'  )    ); 
    DirFileHandle::PrintDumper ( $workingDir."/".'3.2.5.GbPrpIfHash.txt',  $outPepInfHASH  ) if (   (  defined ( $outPepInfHASH )  ) && (  ref ( $outPepInfHASH ) eq 'HASH'  )    ); 
    DirFileHandle::PrintDumper ( $workingDir."/".'3.2.6.PtHhomoHASH.txt',  $PtHhomoHASH    ) if (   (  defined ( $PtHhomoHASH   )  ) && (  ref ( $PtHhomoHASH   ) eq 'HASH'  )    ); 
    
  }
	
  
	
	
	
	
	return $fastaOut;
	
}



# my $fastaOut=DiamondWork::GetGeneBank_diamondTabularVer( $DiaMondOutHash, $workingDir, $top_howMany );
sub GetGeneBank_diamondTabularVer{ #从 tabular格式的diamond结果中，获得所有的sbject的信息，并efeach其genbank文件，并输出fasta文件
	my ( $DiaMondOutHash, $workingDir, $top_howMany )=@_;
	
	my $warnMsgBody="\nIn package  DiamondWork,\tIn sub GetGeneBank_for_diamondOut,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	if  (   (  defined ( $workingDir )  ) && ( $workingDir=~m/\S+/ )   ){
		system ( "mkdir -p $workingDir");
	}
	
	my $whlBlsArHsh; 
  my $EfcProArHsh;
  my $Lcs2accHash;
  my $Acc2LcsHash;
  my $GBfilePthHS;
  my $TPGBflPthHS;
  
  my $LcsMaxIdtHS;
	
	my $GenBankFDir_0_2_a=      $workingDir."/".'3.1.1.GenBankFDir';  system ( "mkdir -p $GenBankFDir_0_2_a");
	my $GBfileTpDir_0_2_b=      $workingDir."/".'3.1.2.GBfileTpDir';  system ( "mkdir -p $GBfileTpDir_0_2_b");
	my $TPGBflPthHS_0_2_b_0=    $workingDir."/".'2.1.1.TPGBflPthHS.txt';
	my $LcsMaxIdtHS_0_1_a=      $workingDir."/".'4.1.1.LcsMaxIdtHS';
	
	if  (   (  defined ( $DiaMondOutHash )  ) && (  ref ( $DiaMondOutHash ) eq 'HASH'  )   ){      	
        	
    #system ("mkdir -p $GenBankFDir_0_2_a");  
    #system ("mkdir -p $GBfileTpDir_0_2_b");  
    
   
    my $query_nb=1;
    foreach my $Qurey_name (    sort (   keys (   %{ $DiaMondOutHash }  )   )    ){                                                                                                                                        
    	                                                                                                                                                                                                          
    	if  (1){                  
    	 
    	  my $Lcs_to_acc_HASH;  my $Acc_to_lcs_HASH;
    	  my $wholeProtArray;                                                                                                                                                                               
    	  my $EfcProteinArray;       my $pt_to_fatcch_NB=0;    
    	                                                                                                                                                                              
    	  foreach my $eachHitHASH (    sort (   keys (   %{ $DiaMondOutHash->{$Qurey_name} }  )   )    ){                                                                                                                                
    	  	                                                                                                                                                                                                      
    	  	if  (   1   ){                                                                                                                   
    	  		#$ProteinsHASH->{ $eachHitHASH->{'_accessionNB'} }=1;                            
    	  		
    	  		
    	  		my $lcsID=$eachHitHASH;
    	  		my $accID=BlastHandle::GetGiNumber ($lcsID);
    	  		
    	  		
    	  		
    	  		
    	  		$Lcs_to_acc_HASH->{$lcsID}=$accID;  
    	  		
    	      #$Acc_to_lcs_HASH->{$accID}=$lcsID;
    	      
    	  		push @{ $wholeProtArray}, $lcsID;                                                                        #my $msg1= "20181214-0-1\$accID=$accID      \$GBfileTpDir_0_2_b=$GBfileTpDir_0_2_b \n"; print $msg1; warn $msg1;
    	  		#my $found_in_local_db=0; $found_in_local_db=GenBankHandle::Check_present_in_JL_ProteinGB_DB ($lcsID);    #my $msg2= "20181214-0-2\$found_in_local_db=$found_in_local_db\n"; print $msg2; warn $msg2;
            #if ( $found_in_local_db == 0 ){                                                                           my $msg3= "20181214-0-3\$lcsID=$lcsID      \$GBfileTpDir_0_2_b=$GBfileTpDir_0_2_b \n"; print $msg3; warn $msg3;
            #	push @{ $EfcProteinArray}, $accID;   $pt_to_fatcch_NB++;    
            #}      		  		                                                                                                                   
    	  		                                                                                                                                      
    	  	}                                                                                                                                                                                                     
    	  }
    	  
    	  my $correctHASH;
    	  my $TParH=GenBankHandle::CorectTheLocusName_for_BlastHitAccName ($wholeProtArray);
    	  if (   (  defined ( $TParH )  ) && (  ref ( $TParH ) eq 'ARRAY'  )   ) {  ( $wholeProtArray, $correctHASH )=@{ $TParH };  }  
    	  else { my $dieMsgg="\nDIE!!!!! something wrong with  GenBankHandle::CorectTheLocusName_for_BlastHitAccName \n\n"; print $dieMsgg; die $dieMsgg; }
    	  
    	  if (   (  defined ( $correctHASH )  ) && (  ref ( $correctHASH ) eq 'HASH'  )   ) {
    	  	foreach my $orgLcs (    sort { $a cmp $b } (   keys (  %{ $correctHASH } )   )    ){
    	  	  my $newLcs=$correctHASH->{$orgLcs};
    	  	  $Lcs_to_acc_HASH->{$newLcs}=$Lcs_to_acc_HASH->{$orgLcs};
    	  	  delete ( $Lcs_to_acc_HASH->{$orgLcs} );
    	  	}
    	  	
    	  }
    	  
    	  if (   (  defined ( $wholeProtArray )  ) && (  ref ( $wholeProtArray ) eq 'ARRAY'  )   ) {
    	    foreach my $ProtLcs (  @{ $wholeProtArray }  ){ 
    	    	my $accID=$Lcs_to_acc_HASH->{$ProtLcs};
    	      $Acc_to_lcs_HASH->{ $accID }=$ProtLcs;
    	      my $found_in_local_db=0; 
    	      $found_in_local_db=GenBankHandle::Check_present_in_JL_ProteinGB_DB ($ProtLcs);    #my $msg2= "20181214-0-2\$found_in_local_db=$found_in_local_db\n"; print $msg2; warn $msg2;
            if ( $found_in_local_db == 0 ){    my $msg3= "20181214-0-3\$ProtLcs=$ProtLcs \n"; print $msg3; warn $msg3;
            	#push @{ $EfcProteinArray}, $accID;   $pt_to_fatcch_NB++;
            	$EfcProArHsh->{$accID}=1;    
            }  
    	    }
    	  }
    	  
    	  if (   (  defined ( $wholeProtArray )  ) && (  ref ( $wholeProtArray ) eq 'ARRAY'  )   ){
    	  	$whlBlsArHsh->{$Qurey_name}=$wholeProtArray;      		  	
    	  }
    	  #if (   (  defined ( $EfcProteinArray )  ) && (  ref ( $EfcProteinArray ) eq 'ARRAY'  )   ){
    	  #	$EfcProArHsh->{$Qurey_name}=$EfcProteinArray;      		  	
    	  #}
    	  if (   (  defined ( $Lcs_to_acc_HASH )  ) && (  ref ( $Lcs_to_acc_HASH ) eq 'HASH'  ) && (  defined ( $Acc_to_lcs_HASH )  ) && (  ref ( $Acc_to_lcs_HASH ) eq 'HASH'  )   ){
    	  	#$Lcs2accHash->{$Qurey_name}=$Lcs_to_acc_HASH; #my $bigHash; $bigHash=ArrayHashChange::PushSmallHash_into_bigHash($bigHash, $smlHash);
    	  	#$Acc2LcsHash->{$Qurey_name}=$Acc_to_lcs_HASH;
    	  	$Lcs2accHash=ArrayHashChange::PushSmallHash_into_bigHash($Lcs2accHash, $Lcs_to_acc_HASH);
    	  	$Acc2LcsHash=ArrayHashChange::PushSmallHash_into_bigHash($Acc2LcsHash, $Acc_to_lcs_HASH);
    	  }      		  
    	  my $TempQrGBFile=  $TPGBflPthHS->{$Qurey_name}=  $GBfileTpDir_0_2_b."/".$query_nb;        		        		  
    	  my $EachQrGBFile=  $GBfilePthHS->{$Qurey_name}=  $GenBankFDir_0_2_a."/".$query_nb;
    	  
          		        		  
    	  
    	  $query_nb++;                                                                                                                                                                               
    	                                                                                              
    	}                                                                                                                                                                                                         
    } 
	
	     	
        	
    DirFileHandle::PrintDumper ( $workingDir."/".'1.1.2.GBfilePthHS.txt',  $GBfilePthHS   ) if (   (  defined ( $GBfilePthHS  )  ) && (  ref ( $GBfilePthHS  ) eq 'HASH'  )    ); 
    DirFileHandle::PrintDumper ( $workingDir."/".'1.1.3.TPGBflPthHS.txt',  $TPGBflPthHS   ) if (   (  defined ( $TPGBflPthHS  )  ) && (  ref ( $TPGBflPthHS  ) eq 'HASH'  )    ); 
        	                       
    DirFileHandle::PrintDumper ( $workingDir."/".'1.1.4.Lcs2accHash.txt',  $Lcs2accHash   ) if (   (  defined ( $Lcs2accHash  )  ) && (  ref ( $Lcs2accHash  ) eq 'HASH'  )    ); 
    DirFileHandle::PrintDumper ( $workingDir."/".'1.1.5.Acc2LcsHash.txt',  $Acc2LcsHash   ) if (   (  defined ( $Acc2LcsHash  )  ) && (  ref ( $Acc2LcsHash  ) eq 'HASH'  )    ); 
        	                                         
    DirFileHandle::PrintDumper ( $workingDir."/".'1.1.6.whlBlsArHsh.txt',  $whlBlsArHsh   ) if (   (  defined ( $whlBlsArHsh  )  ) && (  ref ( $whlBlsArHsh  ) eq 'HASH'  )    ); 
    DirFileHandle::PrintDumper ( $workingDir."/".'1.1.7.EfcProArHsh.txt',  $EfcProArHsh   ) if (   (  defined ( $EfcProArHsh  )  ) && (  ref ( $EfcProArHsh  ) eq 'HASH'  )    ); 
	}
	
	
      
  if  (   (  defined ( $EfcProArHsh )  ) && (  ref ( $EfcProArHsh ) eq 'HASH'  )   ){
    my $AllEfcProArray;  my $Total_Qurey_number=0;
    foreach my $EfcProIDkey (    sort { $a cmp $b } (   keys (  %{ $EfcProArHsh } )   )    ){
    	push @{ $AllEfcProArray }, $EfcProIDkey; $Total_Qurey_number++;
    }
    if (   (  defined ( $TPGBflPthHS_0_2_b_0 )  ) && ( $TPGBflPthHS_0_2_b_0=~m/\S+/ ) && (  defined ( $AllEfcProArray )  ) && (  ref ( $AllEfcProArray ) eq 'ARRAY'  )   ){
    	
    
      my $time_3=`date`; chomp $time_3;
      my $msg4= "Now $time_3 , \ndo the EUtilitieswork :: efatch_proteinGB( \$AllEfcProArra=$AllEfcProArray \$TPGBflPthHS_0_2_b_0=\n\n$TPGBflPthHS_0_2_b_0 )\n \n\n\$Total_Qurey_number=$Total_Qurey_number protein to fatch!!\n\n";  print $msg4; warn $msg4;
      #下载数据
      EUtilitieswork::efatch_proteinGB( $AllEfcProArray, $TPGBflPthHS_0_2_b_0 );  
      
      #更新数据库
      GenBankHandle::Upgrade_JL_proteinGB_and_Acc_to_Lcs_HASHs_onlyOnePid ($TPGBflPthHS_0_2_b_0, $AllEfcProArray, $Acc2LcsHash);      
      
    }
    
  
  }

	
  
    
  my $LcsMaxIdtHS_Changed=0;
  
  ###################
  #if (   (  defined ( $whlBlsArHsh_0_2_e )  ) && (  -e ( $whlBlsArHsh_0_2_e )  )   ){
    #my $whlBlsArHsh=retrieve ($whlBlsArHsh_0_2_e);
  if  (   (  defined ( $whlBlsArHsh )  ) && (  ref ( $whlBlsArHsh ) eq 'HASH'  )   ){
  	foreach my $Qurey_name (    sort { $a cmp $b } (   keys (  %{ $whlBlsArHsh } )   )    ){
  		my $wholeProtArray=$whlBlsArHsh->{$Qurey_name};
  		if  (   (  defined ( $wholeProtArray )  ) && (  ref ( $wholeProtArray ) eq 'ARRAY'  )   ){
        
        my $correctHASH;
	      my $TParH=GenBankHandle::CorectTheLocusName_for_BlastHitAccName ($wholeProtArray);
	      if (   (  defined ( $TParH )  ) && (  ref ( $TParH ) eq 'ARRAY'  )   ) {  ( $wholeProtArray, $correctHASH )=@{ $TParH };  }
	      
     
	      
	      if  (   (  defined ( $GBfilePthHS )  ) && (  ref ( $GBfilePthHS ) eq 'HASH'  ) && (  defined ( $GBfilePthHS->{$Qurey_name} )  ) && ( $GBfilePthHS->{$Qurey_name}=~/\S+/ )   ){
  		    my $EachQrGBFile=  $GBfilePthHS->{$Qurey_name};
	        GenBankHandle::Build_GB_file_from_IDarray_JL_GBidxFILE ( $wholeProtArray, $EachQrGBFile );
	      }
	      
	    }
	  }
	}
	
  ###################
  

  
  
  #print "\n201901042318-0 \$selProFamName=$selProFamName\n";
  my $GbkEfachWkOutARRAY=GenBankHandle::BuildTreesHash_from_geneBankFiles($GenBankFDir_0_2_a, $top_howMany);
  #print "\n201901042318-1 \$selProFamName=$selProFamName\n";
  
  my ($fastaOut, $outPepInfHASH,  $TreeWord, $TreeHash, $AllGBHash, $reapeatHash, $PtHhomoHASH);
  if (   (  defined ( $GbkEfachWkOutARRAY )  ) && (  ref ( $GbkEfachWkOutARRAY ) eq 'ARRAY'  )   ) {
    ($fastaOut, $outPepInfHASH,  $TreeWord, $TreeHash, $AllGBHash, $reapeatHash, $PtHhomoHASH)=@{ $GbkEfachWkOutARRAY };   
    #$FastaSeqString.="\n\n";
    #$FastaSeqString.=$fastaOut;
    #open (OUT1, ">$allFmlProte_0_3_0"  ) or die "cannot create \$allFmlProte_0_3_0=$allFmlProte_0_3_0 : $! \n\n"; print OUT1 $FastaSeqString; close (OUT1);                                                               
    open (OUT1, ">$workingDir/3.2.1.TreeShowTxt.txt"  ) or DieWork::Just_dieWork  ("cannot create $workingDir/3.2.1.TreeShowTxt.txt : $! \n\n" ); print OUT1 $TreeWord;       close (OUT1);                                                               
                                                                                  
    DirFileHandle::PrintDumper ( $workingDir."/".'3.2.2.TreeIfmhash.txt',  $TreeHash       ) if (   (  defined ( $TreeHash      )  ) && (  ref ( $TreeHash      ) eq 'HASH'  )    );                                                                                                                                              
    DirFileHandle::PrintDumper ( $workingDir."/".'3.2.3.allGebkHasH.txt',  $AllGBHash      ) if (   (  defined ( $AllGBHash     )  ) && (  ref ( $AllGBHash     ) eq 'HASH'  )    );  
    DirFileHandle::PrintDumper ( $workingDir."/".'3.2.4.reapeatRort.txt',  $reapeatHash    ) if (   (  defined ( $reapeatHash   )  ) && (  ref ( $reapeatHash   ) eq 'HASH'  )    ); 
    DirFileHandle::PrintDumper ( $workingDir."/".'3.2.5.GbPrpIfHash.txt',  $outPepInfHASH  ) if (   (  defined ( $outPepInfHASH )  ) && (  ref ( $outPepInfHASH ) eq 'HASH'  )    ); 
    DirFileHandle::PrintDumper ( $workingDir."/".'3.2.6.PtHhomoHASH.txt',  $PtHhomoHASH    ) if (   (  defined ( $PtHhomoHASH   )  ) && (  ref ( $PtHhomoHASH   ) eq 'HASH'  )    ); 
    
  }
	
  
	
	
	
	
	return $fastaOut;
	
}


#  my $outHash        = DiamondWork::DimondTabularFmtPrase($inFile, $otFile, $otHash_file); 
sub DimondTabularFmtPrase{  # 根据 diamond 或者 blast 的tabular输出文件， 生成 带 字段名的新文件， 并解析生成hash文件
	my ($inFile, $otFile, $otHash_file )=@_;
	
	my $warnMsgBody="\nIn package  DiamondWork,\tIn sub DimondTabularFmtPrase,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	if (   (  defined ( $inFile )  ) && ( $inFile=~m/\S+/ ) && (  ( -e $inFile )  )   ){			}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n  \$inFile=$inFile  should be a defined File path : $!\n\n\n".$caller_inform ); 	}
	
	if (   (  defined ( $otFile )  ) && ( $otFile=~m/\S+/ )   ){			}
	else{		
		$otFile=$inFile.".AdTitl.txt";
		#DieWork::Just_dieWork( $die_MsgHead."\n  \$otFile=$otFile  should be a defined string could be File path : $!\n\n\n".$caller_inform ); 	
	}
	
	my $FlnTableString = BlastHandle::BuildKeyHead_txt_from_tabularFile($inFile, $otFile);

  my $outHash        = InFileHandle::FormTableToHash_to_2D_hash_Array($otFile, 0, 1);
  
  if  (   (  defined ( $outHash )  ) && (  ref ( $outHash ) eq 'HASH'  )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$outHash=$outHash should be a HASH  !!  $!\n\n\n".$caller_inform ); 	}
  
  if (   (  defined ( $otHash_file )  ) && ( $otHash_file=~m/\S+/ )   ){			}
	else{				$otHash_file=$inFile.".hsh";			}
	
  DirFileHandle::PrintDumper ( $otHash_file, $outHash );  
  
  return $outHash;
	
}
                 


sub GetAllProteinHASH_from_diamondMutWork{
	my ($workingDIR)=@_;
	
}


1;

##########################################################################################################################################
# 
