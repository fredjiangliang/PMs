
#!/usr/bin/perl -w

BEGIN{  push (@INC, '/home/fredjiang/fredProgs/PMs')  };

use strict;
use warnings;

use FastaFileHandle;
use Bio::Index::Fasta;
use Bio::SeqIO;

package  FastaChunkWork;


#####  FASTA文件地址 hash 实例：
#####  my $inFastaFile_HASH;
#####  $inFastaFile_HASH->{"/home/fredjiang/work/Algae/20150522NewDATA/AllSpecies/Nannochloropsis_gaditanaB3120150522/OneFiledGeno/Nannochloropsis_gaditanaB3120150522.fasta"          }=1;  

#  my ( $BigFastaIDkeyHASH, $out_ChunkID_key_HASH )=@{ FastaChunkWork::BuildChunck_fasta_work($inFastaFile_HASH, $outFile_dir ) };    #$chunk_length_limit, $abNormalFastaID_head 都可以不输入
#  my ( $BigFastaIDkeyHASH, $out_ChunkID_key_HASH )=@{ FastaChunkWork::BuildChunck_fasta_work($inFastaFile_HASH, $outFile_dir, $chunk_length_limit, $abNormalFastaID_head) };  
sub BuildChunck_fasta_work{  # 输入是 多个fasta文件的地址 ; 遍历这些地址，生成总的 fasta序列 的信息hash，每个fasta生成一个新的ID

	my ($inFastaFile_HASH, $outFile_dir, $chunk_length_limit, $abNormalFastaID_head)=@_;   #, $blastxType
	
	my $warnMsgBody="\nIn package  FastaChunkWork,\tIn sub BuildChunck_fasta_work,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	my $chuck_hold_dir_inFileLimit=1000;
	
	if (   (  defined ( $outFile_dir )  ) && ( $outFile_dir=~m/^\S+$/ )   ){			
	}
	else{
		DieWork::Just_dieWork( $die_MsgHead."\n \$outFile_dir=$outFile_dir should be right path to hold all the outHash and out chunk files !!  $!\n\n\n".$caller_inform );
	}
	
	if (   (  defined ( $abNormalFastaID_head )  ) && ( $abNormalFastaID_head=~m/^\w+$/ )   ){	
	  if (   (  length ( $abNormalFastaID_head )  ) > 4   ){
	  	$abNormalFastaID_head=~s/^(\w{4}).*$/$1/;
	  }
	}
	else{
		$abNormalFastaID_head="JLfa";
	} 
	my $absNumber_long=1000000000;      #my $showKey= MatrixCsvChange::SprintfKeyHead( $bigestNb, $prNB );
	
	
	if (   (  defined ( $chunk_length_limit )  ) && ( $chunk_length_limit=~m/^\d+$/ ) && ( $chunk_length_limit > 0 )   ){		  
	}
	else{
		$chunk_length_limit=100000000;
	} 
	
	DieWork::Print_and_warn( "20191025-0-0-1 \n\$chunk_length_limit=$chunk_length_limit\n\n" );
	
	
	system (" mkdir -p $outFile_dir");
	
  
	my $all_chunks_hold_dir=$outFile_dir."/0_0_0_all_Chunks_holder";               	#my $all_idx_hold_pthHsh=$outFile_dir."/0_0_1_all_idx_path_HASH";
	my $InFastaPathHashfile=$outFile_dir."/0_0_1_inPutFstPthHS.txt";                DirFileHandle::PrintDumper ( $InFastaPathHashfile,  $inFastaFile_HASH   ) ;
	my $out_ChunkID_key_HASH;
	my $FastaIDKey_HashFile=$outFile_dir."/0_0_2_FastaIDKeyHSH.txt";
	my $ChunkIDKey_HashFile=$outFile_dir."/0_0_3_ChunkIDKeyHSH.txt";
	
	
	if (   (  defined ( $inFastaFile_HASH )  ) && (  ref( $inFastaFile_HASH ) eq 'HASH' )   ){			}
	else{ DieWork::Just_dieWork( $die_MsgHead."\n \$inFastaFile_HASH=$inFastaFile_HASH should be a HASH ref !!  $!\n\n\n".$caller_inform );   }
	
	my $BigFastaIDkeyHASH;  my $abs_ID_number=1; my $biggest_Fasta_length=-999999999999999999999999999;
	my $reverseIDkey_HASH;
	foreach my $eachFastaFile (    sort { $a cmp $b } (   keys (  %{ $inFastaFile_HASH } )   )    ){ 
	  if (   (  defined ( $eachFastaFile )  ) && ( $eachFastaFile=~m/\S+/ ) && (  -e ( $eachFastaFile )  )   ){	
	  	
	  	my $indexFile=$eachFastaFile."Index.txt";	
	    my $inx_obj = Bio::Index::Fasta->new(-filename   => $indexFile,
                                           -write_flag => 1);
	  	
	  	$inx_obj->make_index($eachFastaFile);
	  	
	  	my $fastaIDkeyHASH=FastaFileHandle::BuildHashFromFastaFile_Name_as_Key_checkRepeat ( $eachFastaFile ) ;		
	  	if (   (  defined ( $fastaIDkeyHASH )  ) && (  ref( $fastaIDkeyHASH ) eq 'HASH' )   ){			
	  	  foreach my $eachFastaOrgID (    sort { $a cmp $b } (   keys (  %{ $fastaIDkeyHASH } )   )    ){ 
	  	  	
	  	    my $sptf_ID_nb= MatrixCsvChange::SprintfKeyHead( $absNumber_long, $abs_ID_number );
	  	    my $abs_ID=$abNormalFastaID_head."_".$sptf_ID_nb; #warn "\$abs_ID=$abs_ID \$fastaIDkeyHASH->{ $eachFastaOrgID }->{'0_0_0_orgFile_path'}=$fastaIDkeyHASH->{ $eachFastaOrgID }->{'0_0_0_orgFile_path'}\n";
	  	    
	  	    $BigFastaIDkeyHASH->{$abs_ID}->{'0_0_0_orgFile_path'}=$fastaIDkeyHASH->{ $eachFastaOrgID }->{'0_0_0_orgFile_path'}; #=$inseqFile;  
          $BigFastaIDkeyHASH->{$abs_ID}->{'0_0_1_orgFile_size'}=$fastaIDkeyHASH->{ $eachFastaOrgID }->{'0_0_1_orgFile_size'}; #=$size;  
          $BigFastaIDkeyHASH->{$abs_ID}->{'0_0_2_orgParimayID'}=$fastaIDkeyHASH->{ $eachFastaOrgID }->{'0_0_2_orgParimayID'}; #=$idHere;   
          $BigFastaIDkeyHASH->{$abs_ID}->{'0_0_3_orgMiddle_ID'}=$fastaIDkeyHASH->{ $eachFastaOrgID }->{'0_0_3_orgMiddle_ID'}; #=$middleID;   
          $BigFastaIDkeyHASH->{$abs_ID}->{'0_0_4_orgShort__ID'}=$fastaIDkeyHASH->{ $eachFastaOrgID }->{'0_0_4_orgShort__ID'}; #=$shortID;   
          $BigFastaIDkeyHASH->{$abs_ID}->{'0_0_5_seq___length'}=$fastaIDkeyHASH->{ $eachFastaOrgID }->{'0_0_5_seq___length'}; #=$Seqleth;  
          $BigFastaIDkeyHASH->{$abs_ID}->{'0_0_6_orgFstIdxFil'}=$indexFile;                                                   #=$Seqleth;  
          
          $reverseIDkey_HASH->{ $fastaIDkeyHASH->{ $eachFastaOrgID }->{'0_0_2_orgParimayID'} }=$abs_ID;
          
          if ( $biggest_Fasta_length < $BigFastaIDkeyHASH->{$abs_ID}->{'0_0_5_seq___length'}) {
          	$biggest_Fasta_length=$BigFastaIDkeyHASH->{$abs_ID}->{'0_0_5_seq___length'};
          }
	  	    
	  	    $abs_ID_number++;
	  	  }
	  	}
	  	else{
	  		DieWork::Just_dieWork( $die_MsgHead."\n \$inFastaFile_HASH=$inFastaFile_HASH should be a defined HASH ref !!  $!\n\n\n".$caller_inform );
	  	}
	  }
	  else{
	  	DieWork::Just_dieWork( $die_MsgHead."\n \$eachFastaFile=$eachFastaFile should be a defined fasta file !!  $!\n\n\n".$caller_inform );
	  }
	  
	}
	
	
	if ( $chunk_length_limit < $biggest_Fasta_length ){
		$chunk_length_limit=$biggest_Fasta_length;               DieWork::Print_and_warn( "\n20191025-0-0-1 \$chunk_length_limit=$chunk_length_limit\n\n" ); 
	}
	
	
	my $outPut_Grouped_array; 
	
	#采用了新的blastx处理方案，所以不需要使用下面这种分块方法了
	#my $blastxOrNot=0;
	#if (   (  defined ( $blastxType )  ) && ( $blastxType=~/\S+/ )   ){
	#	$blastxType=~s/\s+//g; $blastxType=lc $blastxType;
	#	if ( $blastxType eq 'blastx' ){
	#		$blastxOrNot=1;			
	#	}
	#}
	#if ( $blastxOrNot==1 ){ $outPut_Grouped_array= FastaChunkWork::disIDsIntoGrpups_987654321($BigFastaIDkeyHASH, $chunk_length_limit); }
	#else                  { $outPut_Grouped_array= FastaChunkWork::disIDsIntoGrpups          ($BigFastaIDkeyHASH, $chunk_length_limit); }
	
	$outPut_Grouped_array= FastaChunkWork::disIDsIntoGrpups          ($BigFastaIDkeyHASH, $chunk_length_limit);
	
	if (   (  defined ( $outPut_Grouped_array )  ) && (  ref( $outPut_Grouped_array ) eq 'ARRAY' )   ){			}
	else{ DieWork::Just_dieWork( $die_MsgHead."\n \$outPut_Grouped_array=$outPut_Grouped_array should be a ARRAY ref !!  $!\n\n\n".$caller_inform );   }
	
	
	
	my $GroupNb=@{ $outPut_Grouped_array }; my $Dir_nb=1; my $dirTotalNB=int ( $GroupNb/$chuck_hold_dir_inFileLimit ); my $Dir_spt_NB= MatrixCsvChange::SprintfKeyHead( $dirTotalNB+1, $Dir_nb );
	system (" mkdir -p $all_chunks_hold_dir/$Dir_spt_NB");
	my $inDirNb=0;
	for (  my $i=0; $i < @{ $outPut_Grouped_array }; $i++  ){
		my $GrpFastaFileNB     = $i+1;
		my $Sptf_GrpFastaFileNB= MatrixCsvChange::SprintfKeyHead( $GroupNb+1, $GrpFastaFileNB );
		my $Sptf_GpFstFile_NAME= $Sptf_GrpFastaFileNB.".txt";
		my $Sftf_GpFstFl_inPath= $Dir_spt_NB."/".$Sptf_GpFstFile_NAME;
		my $Sptf_GpFstFile_PATH= $all_chunks_hold_dir."/".$Sftf_GpFstFl_inPath;
		my $Sptf_GpFstIdxF_PATH= $all_chunks_hold_dir."/".$Dir_spt_NB."/".$Sptf_GrpFastaFileNB.".idx";;
		
		my $GroupFastaString   = FastaChunkWork::BuildFastaStringFor_a_group($outPut_Grouped_array->[$i], $BigFastaIDkeyHASH);  #, $inx_obj
		InFileHandle::PrintStringIntoFile($Sptf_GpFstFile_PATH, $GroupFastaString);
		
		FastaFileHandle::BuildIdxFile_for_fastaFile($Sptf_GpFstFile_PATH, $Sptf_GpFstIdxF_PATH);
		
		my @args = stat ($Sptf_GpFstFile_PATH);  my $size = $args[7];
		
		foreach my $eachAbsID (  @{ $outPut_Grouped_array->[$i] }  ){
			$BigFastaIDkeyHASH->{$eachAbsID}->{'0_1_0_chkFile_path'}=$Sptf_GpFstFile_PATH; 
	  	$BigFastaIDkeyHASH->{$eachAbsID}->{'0_1_1_chkFile_size'}=$size; 
	  	$BigFastaIDkeyHASH->{$eachAbsID}->{'0_1_2_ckFl_In__dir'}=$Dir_spt_NB;
	  	$BigFastaIDkeyHASH->{$eachAbsID}->{'0_1_3_ckFl_In_path'}=$Sftf_GpFstFl_inPath;
	  	$BigFastaIDkeyHASH->{$eachAbsID}->{'0_1_4_chkFile_indx'}=$Sptf_GpFstIdxF_PATH;
	  	$BigFastaIDkeyHASH->{$eachAbsID}->{'0_1_5_chkFile_name'}=$Sptf_GrpFastaFileNB;
	  	
	  	$out_ChunkID_key_HASH->{$Sptf_GpFstFile_PATH}->{$eachAbsID}=Storable::dclone( $BigFastaIDkeyHASH->{$eachAbsID} );
	  	
		}
		
		
		if (  $inDirNb >= ( $chuck_hold_dir_inFileLimit - 1 )  ){
			$inDirNb=0;
			$Dir_nb++;
			$Dir_spt_NB= MatrixCsvChange::SprintfKeyHead( $dirTotalNB+1, $Dir_nb );
			system (" mkdir -p $all_chunks_hold_dir/$Dir_spt_NB");
		}
		else{
		  $inDirNb++;	
		}
		
	}
	
	if (   (  defined ( $BigFastaIDkeyHASH  )  ) && (  ref ( $BigFastaIDkeyHASH  ) eq 'HASH'  )    ){
		DirFileHandle::PrintDumper ( $FastaIDKey_HashFile,  $BigFastaIDkeyHASH   );
	}
	else{ DieWork::Just_dieWork( $die_MsgHead."\n \$BigFastaIDkeyHASH=$BigFastaIDkeyHASH should be a HASH ref !!  $!\n\n\n".$caller_inform );   }
	
	if (   (  defined ( $out_ChunkID_key_HASH  )  ) && (  ref ( $out_ChunkID_key_HASH  ) eq 'HASH'  )    ){
		DirFileHandle::PrintDumper ( $ChunkIDKey_HashFile,  $out_ChunkID_key_HASH   );
	}
	else{ DieWork::Just_dieWork( $die_MsgHead."\n \$out_ChunkID_key_HASH=$out_ChunkID_key_HASH should be a HASH ref !!  $!\n\n\n".$caller_inform );   }
	
	
	
	
	
	return [ $BigFastaIDkeyHASH, $out_ChunkID_key_HASH ];
	
	
	
}



#  my ( $ovlCkForDmd_Hash, $chkPth_to_oDHash, $UplvID_to_oDHash )=@{ FastaChunkWork::BuildChunck_forDiamondBlastx($outFile_dir, $segMentLth, $overLayLth, $abNormalFastaID_head) };  
sub BuildChunck_forDiamondBlastx{  # 输入是 已经完成 fasta Chunk的文件夹 #因为，要做Dimond blastx，所以要进行单个fasta序列的内部分块，目的是让所有的blastx输入DNA序列长度 相仿，防止太长的序列 取太少的结果出现

	my ( $outFile_dir, $segMentLth, $overLayLth, $abNormalFastaID_head )=@_; 
	
	my $warnMsgBody="\nIn package  FastaChunkWork,\tIn sub BuildChunck_forDiamondBlastx,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	if (   (  defined ( $outFile_dir )  ) && ( $outFile_dir=~m/^\S+$/ )   ){				}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$outFile_dir=$outFile_dir should be right path to hold all the outHash and out chunk files !!  $!\n\n\n".$caller_inform );	}
	
	my $in_FastaKey_Hash_file=$outFile_dir."/0_0_2_FastaIDKeyHSH.txt";
	my $in_ChunkKey_Hash_file=$outFile_dir."/0_0_3_ChunkIDKeyHSH.txt";
	
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
	
	
	my $chuck_hold_dir_inFileLimit=1000;		
	
	if (   (  defined ( $segMentLth )  ) && ( $segMentLth=~m/^\d+$/ ) && ( $segMentLth > 0 )   ){		                     }
	else                                                                                        {		$segMentLth=10000;	 } 
	if (   (  defined ( $overLayLth )  ) && ( $overLayLth=~m/^\d+$/ ) && ( $overLayLth > 0 )   ){		                     }
	else                                                                                        {		$overLayLth=1000;	 } 
	
	
  
	#my $all_chunks_hold_dir=$outFile_dir."/0_3_0_all_Chunks_holder";
	my $ovlChuk_for_diam__dir=$outFile_dir."/0_3_0_ovlChuk_for_diamd";
	system ( "mkdir -p $ovlChuk_for_diam__dir" );
	
	my $ovlCkForDmd_Hash_file=$outFile_dir."/0_3_1_ovlCkForDmdHs.txt";
	my $UplvID_to_oDHash_file=$outFile_dir."/0_3_2_UplvIDtoOlDHs.txt";
	my $chkPth_to_oDHash_file=$outFile_dir."/0_3_3_CkPathtoOlDHs.txt";
	
	
	my ( $ovlCkForDmd_Hash, $UplvID_to_oDHash )=@{ FastaChunkWork::BuildOvlayerSubChunk($in_FastaKey_Hash, $segMentLth, $overLayLth, $abNormalFastaID_head) };
	
	if (   (  defined ( $ovlCkForDmd_Hash  )  ) && (  ref ( $ovlCkForDmd_Hash  ) eq 'HASH'  )    )	{
		DirFileHandle::PrintDumper ( $ovlCkForDmd_Hash_file,  $ovlCkForDmd_Hash   );
	}
	else{ DieWork::Just_dieWork( $die_MsgHead."\n \$ovlCkForDmd_Hash=$ovlCkForDmd_Hash should be a HASH ref !!  $!\n\n\n".$caller_inform );   }
	
	if (   (  defined ( $UplvID_to_oDHash  )  ) && (  ref ( $UplvID_to_oDHash  ) eq 'HASH'  )    ){
		DirFileHandle::PrintDumper ( $UplvID_to_oDHash_file,  $UplvID_to_oDHash   );
	}
	else{ DieWork::Just_dieWork( $die_MsgHead."\n \$UplvID_to_oDHash=$UplvID_to_oDHash should be a HASH ref !!  $!\n\n\n".$caller_inform );   }
	
  my $chkPth_to_oDHash;
  foreach my $eachCkPath (    sort { $a cmp $b } (   keys (  %{ $in_ChunkKey_Hash } )   )    ){  #0_1_3_ckFl_In_path
  	if  (   (  defined ( $in_ChunkKey_Hash->{$eachCkPath} )  ) && (  ref ( $in_ChunkKey_Hash->{$eachCkPath} ) eq 'HASH'  )   ){}
	  else{		DieWork::Just_dieWork( $die_MsgHead."\n \$in_ChunkKey_Hash->{$eachCkPath}=$in_ChunkKey_Hash->{$eachCkPath} should be a HASH  !!  $!\n\n\n".$caller_inform ); 	}
	  
	  my $newCkFileIN_dir;
	  my $newCkFileInPath;
	  my $absNeCkFile_Dir;
		my $absNeCkFilePath;
		my $newCkFileInName;
		my $absNeCk_IdxPath;
	  
  	my $ovlyCkIDARRAY; my $idx=0;
  	foreach my $eachAbsID (    sort { $a cmp $b } (   keys (  %{ $in_ChunkKey_Hash->{$eachCkPath} } )   )    ){ 
			
			if  (   (  defined ( $UplvID_to_oDHash->{$eachAbsID} )  ) && (  ref ( $UplvID_to_oDHash->{$eachAbsID} ) eq 'HASH'  )   ){}
	    else{		DieWork::Just_dieWork( $die_MsgHead."\n \$UplvID_to_oDHash->{$eachAbsID}=$UplvID_to_oDHash->{$eachAbsID} should be a HASH  !!  $!\n\n\n".$caller_inform ); 	}
			
			$newCkFileIN_dir=$in_FastaKey_Hash->{$eachAbsID}->{'0_1_2_ckFl_In__dir'};			
			$newCkFileInPath=$in_FastaKey_Hash->{$eachAbsID}->{'0_1_3_ckFl_In_path'};
			$newCkFileInName=$in_FastaKey_Hash->{$eachAbsID}->{'0_1_5_chkFile_name'};
			$absNeCkFile_Dir=$ovlChuk_for_diam__dir."/".$newCkFileIN_dir;
			$absNeCkFilePath=$ovlChuk_for_diam__dir."/".$newCkFileInPath;
			$absNeCk_IdxPath=$ovlChuk_for_diam__dir."/".$newCkFileIN_dir."/".$newCkFileInName.".idx";
			
			foreach my $OvLCkID (    sort { $a cmp $b } (   keys (  %{ $UplvID_to_oDHash->{$eachAbsID} } )   )    ){
			  $ovlyCkIDARRAY->[$idx]=$OvLCkID;  print  "\n 20190418-0-0-0 \$ovlyCkIDARRAY->[\$idx]=\$ovlyCkIDARRAY->[$idx]=\$OvLCkID=$OvLCkID\n";
			  
			  $ovlCkForDmd_Hash->{$OvLCkID}->{'0_2_0_olyCkFl_path'}=$absNeCkFilePath;
			  $ovlCkForDmd_Hash->{$OvLCkID}->{'0_2_1_olyCkFl_idxP'}=$absNeCk_IdxPath;
			  #$ovlCkForDmd_Hash->{$OvLCkID}->{'0_2_2_olyCkFlINpth'}=$absNeCk_IdxPath;
			  
			  $chkPth_to_oDHash->{$absNeCkFilePath}->{$OvLCkID}=Storable::dclone ( $ovlCkForDmd_Hash->{$OvLCkID} );
			   
			  $idx++; 
			}
	  	
		}
		system ( " mkdir -p $absNeCkFile_Dir " );    DieWork::Print_and_warn( "\n 20190418-0-0-0-0 \$GroupFastaString=\$GroupFastaString\n" );
		my $GroupFastaString   = FastaChunkWork::BuildFastaStringFor_OlyChuk ($ovlyCkIDARRAY, $ovlCkForDmd_Hash);  DieWork::Print_and_warn( "\n 20190418-0-0-0-1 \$GroupFastaString=\$GroupFastaString\n" ); #, $inx_obj
		InFileHandle::PrintStringIntoFile($absNeCkFilePath, $GroupFastaString);
		FastaFileHandle::BuildIdxFile_for_fastaFile($absNeCkFilePath, $absNeCk_IdxPath);  
  }
  
  
  if (   (  defined ( $ovlCkForDmd_Hash  )  ) && (  ref ( $ovlCkForDmd_Hash  ) eq 'HASH'  )    ){
		DirFileHandle::PrintDumper ( $ovlCkForDmd_Hash_file,  $ovlCkForDmd_Hash   );
	}
	else{ DieWork::Just_dieWork( $die_MsgHead."\n \$ovlCkForDmd_Hash=$ovlCkForDmd_Hash should be a HASH ref !!  $!\n\n\n".$caller_inform );   }
	
	if (   (  defined ( $chkPth_to_oDHash  )  ) && (  ref ( $chkPth_to_oDHash  ) eq 'HASH'  )    ){
		DirFileHandle::PrintDumper ( $chkPth_to_oDHash_file,  $chkPth_to_oDHash   );
	}
	else{ DieWork::Just_dieWork( $die_MsgHead."\n \$chkPth_to_oDHash=$chkPth_to_oDHash should be a HASH ref !!  $!\n\n\n".$caller_inform );   }

	return [ $ovlCkForDmd_Hash, $chkPth_to_oDHash, $UplvID_to_oDHash ];
	
	
	
}

#  my ( $newOutHash, $upLvFastaID_to_lowLvID_hash )=@{ FastaChunkWork::BuildOvlayerSubChunk($fastaIDkeyHASH, $segMentLth, $overLayLth, $abNormalFastaID_head) };
sub BuildOvlayerSubChunk{ #对已经重命名的fasta，包括其hash文件，进行再次分割，分割的方法是 根据特定的重叠长度，从头分割成n个新的fasta。
	my ($fastaIDkeyHASH, $segMentLth, $overLayLth, $abNormalFastaID_head)=@_;
	
	my $warnMsgBody="\nIn package  FastaChunkWork,\tIn sub BuildOvlayerSubChunk,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  my $caller_inform=DirFileHandle::print_SubCallerInform;
  
	
	if (   (  defined ( $abNormalFastaID_head )  ) && ( $abNormalFastaID_head=~m/^\w+$/ )   ){	
	  if (   (  length ( $abNormalFastaID_head )  ) > 4   ){
	  	$abNormalFastaID_head=~s/^(\w{4}).*$/$1/;
	  }
	}
	else{
		$abNormalFastaID_head="JLol";
	} 
	my $absNumber_long=1000000000;      #my $showKey= MatrixCsvChange::SprintfKeyHead( $bigestNb, $prNB );
	
	
	if (   (  defined ( $fastaIDkeyHASH  )  ) && (  ref ( $fastaIDkeyHASH  ) eq 'HASH'  )    ){	}
	else{ DieWork::Just_dieWork( $die_MsgHead."\n \$fastaIDkeyHASH=$fastaIDkeyHASH should be a HASH ref !!  $!\n\n\n".$caller_inform );   }
	
	
	my $newOutHash; 
	my $upLvFastaID_to_lowLvID_hash;
	
	my $new_NB=1;
	foreach my $eachFastaID (    sort { $a cmp $b } (   keys (  %{ $fastaIDkeyHASH } )   )    ){ 
		my $FastaLength= $fastaIDkeyHASH->{$eachFastaID}->{'0_0_5_seq___length'}; 
		
		my $SegHdTl_Array=SeqSegmentsTools::BuildOverLayerPosPair_for_a_length($FastaLength, $segMentLth, $overLayLth);
		
		for (  my $i=0; $i<@{ $SegHdTl_Array }; $i++  ){
			
			my $sptf_ID_nb= MatrixCsvChange::SprintfKeyHead( $absNumber_long, $new_NB );
      my $abs_ID    = $abNormalFastaID_head."_".$sptf_ID_nb;
			
			$newOutHash->{$abs_ID}=Storable::dclone( $fastaIDkeyHASH->{$eachFastaID} );
			
			$newOutHash->{$abs_ID}->{'0_3_1_uplvl_JLfsID'}=$eachFastaID; 
			$newOutHash->{$abs_ID}->{'0_3_2_Sgm_Lgth_Set'}=$segMentLth; 
			$newOutHash->{$abs_ID}->{'0_3_3_ovl_Lgth_Set'}=$overLayLth; 
			$newOutHash->{$abs_ID}->{'0_3_4_Sgm______stt'}=$SegHdTl_Array->[$i]->{'0_SgmtHead'}; 
			$newOutHash->{$abs_ID}->{'0_3_5_Sgm______end'}=$SegHdTl_Array->[$i]->{'1_SgmtTail'}; 
			$newOutHash->{$abs_ID}->{'0_3_6_Sgm______lth'}=$SegHdTl_Array->[$i]->{'1_SgmtTail'}-$SegHdTl_Array->[$i]->{'0_SgmtHead'}+1; 
			
			$upLvFastaID_to_lowLvID_hash->{$eachFastaID}->{$abs_ID}=1;
			
			$new_NB++;
		}
		
	}
	
	if (   (  defined ( $newOutHash  )  ) && (  ref ( $newOutHash  ) eq 'HASH'  )    ){		}
	else{ DieWork::Just_dieWork( $die_MsgHead."\n \$newOutHash=$newOutHash should be a HASH ref !!  $!\n\n\n".$caller_inform );   }
	
	if (   (  defined ( $upLvFastaID_to_lowLvID_hash  )  ) && (  ref ( $upLvFastaID_to_lowLvID_hash  ) eq 'HASH'  )    ){		}
	else{ DieWork::Just_dieWork( $die_MsgHead."\n \$upLvFastaID_to_lowLvID_hash=$upLvFastaID_to_lowLvID_hash should be a HASH ref !!  $!\n\n\n".$caller_inform );   }
	
	return [ $newOutHash, $upLvFastaID_to_lowLvID_hash ];
	
}

#BuildFastaStringFor_OlyChuk
#my $finalOutString=FastaChunkWork::BuildFastaStringFor_OlyChuk($inGroup, $big_org_infromHASH);
sub BuildFastaStringFor_OlyChuk{ #输入 id group的array的ref，和大的信息hash，得到group的fasta string
  my ($inGroup, $big_org_infromHASH)=@_;  #, $inx_obj
	
	my $warnMsgBody="\nIn package  FastaChunkWork,\tIn sub BuildFastaStringFor_OlyChuk,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  my $caller_inform=DirFileHandle::print_SubCallerInform;
  
  if (   (  defined ( $inGroup )  ) && (  ref ( $inGroup ) eq 'ARRAY' )   ){			}
	else{ DieWork::Just_dieWork( $die_MsgHead."\n \$inGroup=$inGroup should be a ARRAY ref !!  $!\n\n\n".$caller_inform );   }
	
	if (   (  defined ( $big_org_infromHASH )  ) && (  ref ( $big_org_infromHASH ) eq 'HASH' )   ){			}
	else{ DieWork::Just_dieWork( $die_MsgHead."\n \$big_org_infromHASH=$big_org_infromHASH should be a HASH ref !!  $!\n\n\n".$caller_inform );   }
	
	#if (   (  defined ( $reverseIDkey_HASH )  ) && (  ref ( $reverseIDkey_HASH ) eq 'HASH' )   ){			}
	#else{ DieWork::Just_dieWork( $die_MsgHead."\n \$reverseIDkey_HASH=$reverseIDkey_HASH should be a HASH ref !!  $!\n\n\n".$caller_inform );   }
	
	
	
	
	my $subGroupHASH=ArrayHashChange::GetSubHASH_from_bigHASH ($inGroup, $big_org_infromHASH);

	
	
	my $finalOutString;
	
	my $inIndex_hash=ArrayHashChange::Change_2d_hashINTO_3d_hash ($subGroupHASH, '0_0_6_orgFstIdxFil');
	if (   (  defined ( $inIndex_hash )  ) && (  ref ( $inIndex_hash ) eq 'HASH' )   ){			}
	else{ DieWork::Just_dieWork( $die_MsgHead."\n \$inIndex_hash=$inIndex_hash should be a HASH ref !!  $!\n\n\n".$caller_inform );   }
	

	
	my $idx=0;
	foreach my $indexFile (    sort { $a cmp $b } (   keys (  %{ $inIndex_hash }  )   )    ){  DieWork::Print_and_warn( "\n 20190418-0-0-1-0 \$indexFile=$indexFile\n" );
	  my $inx_obj   = Bio::Index::Fasta->new(-filename   => $indexFile);
	  
	  if (   (  defined ( $inIndex_hash->{$indexFile} )  ) && (  ref ( $inIndex_hash->{$indexFile} ) eq 'HASH' )   ){			}
	  else{ DieWork::Just_dieWork( $die_MsgHead."\n \$inIndex_hash->{$indexFile}=$inIndex_hash->{$indexFile} should be a HASH ref !!  $!\n\n\n".$caller_inform );   }
	  my $Prim_ID_HASH=ArrayHashChange::Change_2d_hashINTO_3d_hash ( $inIndex_hash->{$indexFile}, '0_0_2_orgParimayID' );
	  
	  foreach my $prim_ID (    sort { $a cmp $b } (   keys (  %{ $Prim_ID_HASH }  )   )    ){  DieWork::Print_and_warn( "\n 20190418-0-0-1-1 \$prim_ID=$prim_ID\n" );
	    my $seq_IO= $inx_obj->fetch($prim_ID);
	    
	    if (   (  defined ( $Prim_ID_HASH->{$prim_ID} )  ) && (  ref ( $Prim_ID_HASH->{$prim_ID} ) eq 'HASH' )   ){			}
	    else{ DieWork::Just_dieWork( $die_MsgHead."\n \$Prim_ID_HASH->{$prim_ID}=$Prim_ID_HASH->{$prim_ID} should be a HASH ref !!  $!\n\n\n".$caller_inform );   }
	    
	    foreach my $eachJL_ID (    sort { $a cmp $b } (   keys (  %{ $Prim_ID_HASH->{$prim_ID} }  )   )    ){  DieWork::Print_and_warn( "\n 20190418-0-0-1-2 \$eachJL_ID=$eachJL_ID\n" );
	      
	      #my $prim_ID  =$subGroupHASH->{$eachJL_ID}->{'0_0_2_orgParimayID'}; 		    #my $indexFile=$subGroupHASH->{$eachJL_ID}->{'0_0_6_orgFstIdxFil'};		    #my $inx_obj  = Bio::Index::Fasta->new(-filename   => $indexFile);
		    
		    #my $seq_IO   = $inx_obj->fetch($prim_ID);		    #my $seq_IO   =$Idx_Obj_SeqIO_HASH->{$indexFile};
		    my $Seq_hd   =$Prim_ID_HASH->{$prim_ID}->{$eachJL_ID}->{'0_3_4_Sgm______stt'};    DieWork::Print_and_warn( "\n 20190418-0-0-1-3 \$Seq_hd=$Seq_hd\n" );
		    my $Seq_ed   =$Prim_ID_HASH->{$prim_ID}->{$eachJL_ID}->{'0_3_5_Sgm______end'};    DieWork::Print_and_warn( "\n 20190418-0-0-1-4 \$Seq_ed=$Seq_ed\n" );
		    my $sequence_here=$seq_IO->subseq ( $Seq_hd, $Seq_ed ); $sequence_here=~s/\s+//g;  #DieWork::Print_and_warn(  "\n 20190418-0-0-2 \$Seq_hd=$Seq_hd \$Seq_hd=$Seq_hd \$sequence_here=$sequence_here \n" );  
		    # 下面的语句是因为 某些est中含有*,所以进行了这个操作，这个操作会 改变碱基的位置，会带来潜在的错误，正确的做法是 在cap3做Est拼接的时候，就把这个*去掉。
		    $sequence_here=~s/\*//g;
		    my $Seq_length=length ($sequence_here);                                             DieWork::Print_and_warn( "\n 20190418-0-0-1-5 \$Seq_length=$Seq_length\n" );
		    $finalOutString.=">".$eachJL_ID."\n".$sequence_here."\n\n";	
	      my $finalOutString_length=length ($finalOutString);                                 DieWork::Print_and_warn( "\n 20190418-0-0-1-6 \$finalOutString_length=$finalOutString_length\n" ); 
	    }
	    
	  }
	  
		
			
	}
	
	
	
	
	#foreach my $eachJL_ID (  @{ $inGroup }  ){  DieWork::Print_and_warn(  "\n 20190418-0-0-1-2 \$eachJL_ID=$eachJL_ID\n" );#0_0_2_orgParimayID
	#	
	#	if (   (  defined ( $subGroupHASH->{$eachJL_ID}->{'0_0_2_orgParimayID'} )  ) && (  $subGroupHASH->{$eachJL_ID}->{'0_0_2_orgParimayID'}=~m/\S+/ )   ) {
	#	  my $prim_ID  =$subGroupHASH->{$eachJL_ID}->{'0_0_2_orgParimayID'};
	#	  my $indexFile=$subGroupHASH->{$eachJL_ID}->{'0_0_6_orgFstIdxFil'};
	#	  #my $inx_obj  = Bio::Index::Fasta->new(-filename   => $indexFile);
	#	  
	#	  #my $seq_IO   = $inx_obj->fetch($prim_ID);
	#	  my $seq_IO   =$Idx_Obj_SeqIO_HASH->{$indexFile};
	#	  my $Seq_hd   =$subGroupHASH->{$eachJL_ID}->{'0_3_4_Sgm______stt'};  
	#	  my $Seq_ed   =$subGroupHASH->{$eachJL_ID}->{'0_3_5_Sgm______end'};
	#	  my $sequence_here=$seq_IO->subseq ( $Seq_hd, $Seq_ed ); $sequence_here=~s/\s+//g;  print  "\n 20190418-0-0-2 \$Seq_hd=$Seq_hd \$Seq_hd=$Seq_hd \$sequence_here=$sequence_here \n";  
	#	  # 下面的语句是因为 某些est中含有*,所以进行了这个操作，这个操作会 改变碱基的位置，会带来潜在的错误，正确的做法是 在cap3做Est拼接的时候，就把这个*去掉。
	#	  $sequence_here=~s/\*//g;
	#	  
	#	  $finalOutString.=">".$eachJL_ID."\n".$sequence_here."\n\n";	
	#	}
	#	else{
	#		DieWork::Just_dieWork( $die_MsgHead."\n \$subGroupHASH->{$eachJL_ID}->{'0_0_2_orgParimayID'}=$subGroupHASH->{$eachJL_ID}->{'0_0_2_orgParimayID'} should be a fasta name id !!  $!\n\n\n".$caller_inform );
	#	}
	#	
	#	
	#}
	
	return $finalOutString;
}


#my $finalOutString=FastaChunkWork::BuildFastaStringFor_a_group($inGroup, $big_org_infromHASH);
sub BuildFastaStringFor_a_group{ #输入 id group的array的ref，和大的信息hash，得到group的fasta string
  my ($inGroup, $big_org_infromHASH)=@_;  #, $inx_obj
	
	my $warnMsgBody="\nIn package  FastaChunkWork,\tIn sub BuildFastaStringFor_a_group,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  my $caller_inform=DirFileHandle::print_SubCallerInform;
  
  if (   (  defined ( $inGroup )  ) && (  ref ( $inGroup ) eq 'ARRAY' )   ){			}
	else{ DieWork::Just_dieWork( $die_MsgHead."\n \$inGroup=$inGroup should be a ARRAY ref !!  $!\n\n\n".$caller_inform );   }
	
	if (   (  defined ( $big_org_infromHASH )  ) && (  ref ( $big_org_infromHASH ) eq 'HASH' )   ){			}
	else{ DieWork::Just_dieWork( $die_MsgHead."\n \$big_org_infromHASH=$big_org_infromHASH should be a HASH ref !!  $!\n\n\n".$caller_inform );   }
	
	#if (   (  defined ( $reverseIDkey_HASH )  ) && (  ref ( $reverseIDkey_HASH ) eq 'HASH' )   ){			}
	#else{ DieWork::Just_dieWork( $die_MsgHead."\n \$reverseIDkey_HASH=$reverseIDkey_HASH should be a HASH ref !!  $!\n\n\n".$caller_inform );   }
	
	
	
	
	my $subGroupHASH=ArrayHashChange::GetSubHASH_from_bigHASH ($inGroup, $big_org_infromHASH);

	
	
	my $finalOutString;
	
	my $inIndex_hash=ArrayHashChange::Change_2d_hashINTO_3d_hash ($subGroupHASH, '0_0_6_orgFstIdxFil');
	if (   (  defined ( $inIndex_hash )  ) && (  ref ( $inIndex_hash ) eq 'HASH' )   ){			}
	else{ DieWork::Just_dieWork( $die_MsgHead."\n \$inIndex_hash=$inIndex_hash should be a HASH ref !!  $!\n\n\n".$caller_inform );   }
	

	
	my $idx=0;
	foreach my $indexFile (    sort { $a cmp $b } (   keys (  %{ $inIndex_hash }  )   )    ){  DieWork::Print_and_warn( "\n 20190417-0-0-1-0 \$indexFile=$indexFile\n" );
	  my $inx_obj   = Bio::Index::Fasta->new(-filename   => $indexFile);
	  
	  if (   (  defined ( $inIndex_hash->{$indexFile} )  ) && (  ref ( $inIndex_hash->{$indexFile} ) eq 'HASH' )   ){			}
	  else{ DieWork::Just_dieWork( $die_MsgHead."\n \$inIndex_hash->{$indexFile}=$inIndex_hash->{$indexFile} should be a HASH ref !!  $!\n\n\n".$caller_inform );   }
	  my $Prim_ID_HASH=ArrayHashChange::Change_2d_hashINTO_3d_hash ( $inIndex_hash->{$indexFile}, '0_0_2_orgParimayID' );
	  
	  foreach my $prim_ID (    sort { $a cmp $b } (   keys (  %{ $Prim_ID_HASH }  )   )    ){  DieWork::Print_and_warn( "\n 20190417-0-0-1-1 \$prim_ID=$prim_ID\n" );
	    my $seq_IO= $inx_obj->fetch($prim_ID);
	    
	    if (   (  defined ( $Prim_ID_HASH->{$prim_ID} )  ) && (  ref ( $Prim_ID_HASH->{$prim_ID} ) eq 'HASH' )   ){			}
	    else{ DieWork::Just_dieWork( $die_MsgHead."\n \$Prim_ID_HASH->{$prim_ID}=$Prim_ID_HASH->{$prim_ID} should be a HASH ref !!  $!\n\n\n".$caller_inform );   }
	    
	    foreach my $eachJL_ID (    sort { $a cmp $b } (   keys (  %{ $Prim_ID_HASH->{$prim_ID} }  )   )    ){  DieWork::Print_and_warn( "\n 20190417-0-0-1-2 \$eachJL_ID=$eachJL_ID\n" );
	      
	      #my $prim_ID  =$subGroupHASH->{$eachJL_ID}->{'0_0_2_orgParimayID'}; 		    #my $indexFile=$subGroupHASH->{$eachJL_ID}->{'0_0_6_orgFstIdxFil'};		    #my $inx_obj  = Bio::Index::Fasta->new(-filename   => $indexFile);
		    
		    #my $seq_IO   = $inx_obj->fetch($prim_ID);		    #my $seq_IO   =$Idx_Obj_SeqIO_HASH->{$indexFile};
		    my $sequence_here=$seq_IO->seq (  ); $sequence_here=~s/\s+//g;  #DieWork::Print_and_warn(  "\n 20190418-0-0-2 \$Seq_hd=$Seq_hd \$Seq_hd=$Seq_hd \$sequence_here=$sequence_here \n" );  
		    # 下面的语句是因为 某些est中含有*,所以进行了这个操作，这个操作会 改变碱基的位置，会带来潜在的错误，正确的做法是 在cap3做Est拼接的时候，就把这个*去掉。
		    $sequence_here=~s/\*//g;
		    my $Seq_length=length ($sequence_here);                                             DieWork::Print_and_warn( "\n 20190417-0-0-1-5 \$Seq_length=$Seq_length\n" );
		    $finalOutString.=">".$eachJL_ID."\n".$sequence_here."\n\n";	
	      my $finalOutString_length=length ($finalOutString);                                 DieWork::Print_and_warn( "\n 20190417-0-0-1-6 \$finalOutString_length=$finalOutString_length\n" ); 
	    }
	    
	  }
	  
		
			
	}
	
	##################################################################################
	
	
	
	#my $finalOutString;
	#
	##my $filePathHASH; 
	#
	#my $Idx_Obj_SeqIO_HASH; 
	#foreach my $eachJL_ID (  @{ $inGroup }  ){
	#	if (   (  defined ( $big_org_infromHASH->{$eachJL_ID}->{'0_0_2_orgParimayID'} )  ) && (  $big_org_infromHASH->{$eachJL_ID}->{'0_0_2_orgParimayID'}=~m/\S+/ )   ) {
	#	  my $prim_ID   =$big_org_infromHASH->{$eachJL_ID}->{'0_0_2_orgParimayID'};
	#	  my $indexFile = $big_org_infromHASH->{$eachJL_ID}->{'0_0_6_orgFstIdxFil'};
	#	  my $inx_obj   = Bio::Index::Fasta->new(-filename   => $indexFile);
	#	  my $seq_IO_obj= $inx_obj->fetch($prim_ID);
	#	  $Idx_Obj_SeqIO_HASH->{$eachJL_ID}=$seq_IO_obj;
	#	}
	#	else{
	#		DieWork::Just_dieWork( $die_MsgHead."\n \$big_org_infromHASH->{$eachJL_ID}->{'0_0_2_orgParimayID'}=$big_org_infromHASH->{$eachJL_ID}->{'0_0_2_orgParimayID'} should be a fasta name id !!  $!\n\n\n".$caller_inform );
	#	}
	#}
	#
	#
	#my $idx=0;
	#foreach my $eachJL_ID (  @{ $inGroup }  ){  #0_0_2_orgParimayID
	#	
	#	if (   (  defined ( $big_org_infromHASH->{$eachJL_ID}->{'0_0_2_orgParimayID'} )  ) && (  $big_org_infromHASH->{$eachJL_ID}->{'0_0_2_orgParimayID'}=~m/\S+/ )   ) {
	#	  my $prim_ID  =$big_org_infromHASH->{$eachJL_ID}->{'0_0_2_orgParimayID'};
	#	  my $indexFile=$big_org_infromHASH->{$eachJL_ID}->{'0_0_6_orgFstIdxFil'};
	#	  #my $inx_obj = Bio::Index::Fasta->new(-filename   => $indexFile);		  
	#	  #my $seq_IO = $inx_obj->fetch($prim_ID);
	#	  my $seq_IO = $Idx_Obj_SeqIO_HASH->{$eachJL_ID};
	#	  my $sequence_here=$seq_IO->seq(); $sequence_here=~s/\s+//g;
	#	  $finalOutString.=">".$eachJL_ID."\n".$sequence_here."\n\n";	
	#	}
	#	else{
	#		DieWork::Just_dieWork( $die_MsgHead."\n \$big_org_infromHASH->{$eachJL_ID}->{'0_0_2_orgParimayID'}=$big_org_infromHASH->{$eachJL_ID}->{'0_0_2_orgParimayID'} should be a fasta name id !!  $!\n\n\n".$caller_inform );
	#	}
	#	
	#	#if (   (  defined ( $big_org_infromHASH->{$eachJL_ID}->{'0_0_0_orgFile_path'} )  ) && (  $big_org_infromHASH->{$eachJL_ID}->{'0_0_0_orgFile_path'}=~m/\S+/ )   ) {
	#	#  $filePathHASH->{ $big_org_infromHASH->{$eachJL_ID}->{'0_0_0_orgFile_path'} }->{$eachJL_ID}=$idx;	
	#	#}
	#	#else{
	#	#	DieWork::Just_dieWork( $die_MsgHead."\n \$big_org_infromHASH->{$eachJL_ID}->{'0_0_0_orgFile_path'}=$big_org_infromHASH->{$eachJL_ID}->{'0_0_0_orgFile_path'} should be a right file path !!  $!\n\n\n".$caller_inform );
	#	#}
	#	#$idx++;
	#}
	
	#if (   (  defined ( $filePathHASH )  ) && (  ref ( $filePathHASH ) eq 'HASH' )   ){			}
	#else{ DieWork::Just_dieWork( $die_MsgHead."\n \$filePathHASH=$filePathHASH should be a HASH ref !!  $!\n\n\n".$caller_inform );   }
	#
	#my $id_Seq_HASH;
	#foreach my $filePath (    sort  (   keys (  %{ $filePathHASH }  )   )    ){ 
	#	if (   (  defined ( $filePath )  ) && ( $filePath=~m/^\S+$/ ) && (  -e ( $filePath )  )   ){			
	#    my $AllFastaSeqHash=FastaFileHandle::BuildFastaHash_Name_Key_seq_Val ( $filePath ) ;
	#    	   
	#    if (   (  defined ( $filePathHASH->{$filePath} )  ) && (  ref ( $filePathHASH->{$filePath} ) eq 'HASH' )   ){	
	#    	foreach my $JLID_here (    sort (   keys (  %{ $filePathHASH->{$filePath} }  )   )    ){
	#    		if (   (  defined ( $big_org_infromHASH->{$JLID_here} )  ) && (  ref ( $big_org_infromHASH->{$JLID_here} ) eq 'HASH' ) && (  defined ( $big_org_infromHASH->{$JLID_here}->{'0_0_2_orgParimayID'} )  ) && ( $big_org_infromHASH->{$JLID_here}->{'0_0_2_orgParimayID'}=~m/^\S+$/ )   ){	
	#    		  my $prim_ID=$big_org_infromHASH->{$JLID_here}->{'0_0_2_orgParimayID'};
	#    		  if (   (  defined ( $AllFastaSeqHash->{$prim_ID} )  ) && ( $AllFastaSeqHash->{$prim_ID}=~m/^\S+$/ )   ){
	#    		  	$id_Seq_HASH->{$JLID_here} = $AllFastaSeqHash->{$prim_ID};
	#    		  }
	#    		  else{
	#    		  	DieWork::Just_dieWork( $die_MsgHead."\n \$AllFastaSeqHash->{$prim_ID}=$AllFastaSeqHash->{$prim_ID} should be defined and not null !!  $!\n\n\n".$caller_inform );
	#    		  }
	#    		}
	#    		else{
	#    			DieWork::Just_dieWork( $die_MsgHead."\n \$big_org_infromHASH->{$JLID_here}=$big_org_infromHASH->{$JLID_here} should be a HASH ref \nand\n \$big_org_infromHASH->{$JLID_here}->{'0_0_2_orgParimayID'}=$big_org_infromHASH->{$JLID_here}->{'0_0_2_orgParimayID'} should not be null !!  $!\n\n\n".$caller_inform );
	#    		}	    		
	#    	}
	#    }	   
	#    	    
	#    
	#     
	#  }	  
  #	else{
	#  	DieWork::Just_dieWork( $die_MsgHead."\n \$filePath=$filePath should be right file path !!  $!\n\n\n".$caller_inform );
  #	}
  #}
	#
	#if (   (  defined ( $id_Seq_HASH )  ) && (  ref ( $id_Seq_HASH ) eq 'HASH' )   ){			}
	#else{ DieWork::Just_dieWork( $die_MsgHead."\n \$id_Seq_HASH=$id_Seq_HASH should be a HASH ref !!  $!\n\n\n".$caller_inform );   }
	#
	#my $finalOutString;
	#foreach my $JL_ID_KEY (    sort { $a cmp $b } (   keys (  %{ $id_Seq_HASH }  )   )    ){ 
	#	if (   (  defined ( $JL_ID_KEY )  ) && ( $JL_ID_KEY=~m/^\S+$/ ) && (  defined ( $id_Seq_HASH->{$JL_ID_KEY} )  ) && ( $id_Seq_HASH->{$JL_ID_KEY}=~m/^\S+$/ )   ){
	#	  $finalOutString.=">".$JL_ID_KEY."\n".$id_Seq_HASH->{$JL_ID_KEY}."\n\n";	
	#	}
	#	else{
	#		DieWork::Just_dieWork( $die_MsgHead."\n \$JL_ID_KEY=$JL_ID_KEY \$id_Seq_HASH->{$JL_ID_KEY}=$id_Seq_HASH->{$JL_ID_KEY} all should be defined no empty string !!  $!\n\n\n".$caller_inform );
	#	}
	#}
	
	return $finalOutString;
}


# my $outPutSegKeyGroups = FastaChunkWork::disIDsIntoGrpups($input_SegSize_Hash, $sizeLimit);
#            disIDsIntoGrpups{   #将很多序列，按照文件大小的要求，放入一系列小于该大小要求的数组中，该数组可用于下一步的文件生成
sub disIDsIntoGrpups{   #将很多序列，按照文件大小的要求，放入一系列小于该大小要求的数组中，该数组可用于下一步的文件生成
  my ($input_SegSize_Hash, $sizeLimit)=@_;
  
  my $warnMsgBody="\nIn package  FastaChunkWork,\tIn sub disIDsIntoGrpups,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  my $caller_inform=DirFileHandle::print_SubCallerInform;
	
  if (   (  defined ( $input_SegSize_Hash )  ) && (  ref( $input_SegSize_Hash ) eq 'HASH' )   ){			}
	else{ DieWork::Just_dieWork( $die_MsgHead."\n \$input_SegSize_Hash=$input_SegSize_Hash should be a HASH ref !!  $!\n\n\n".$caller_inform );   }
	
  my $inputSegSizeHash=Storable::dclone (  $input_SegSize_Hash  );
  
  my $outPutSegKeyGroups;
  my $oPsKgIdx=0;
  foreach my $BigNumberKey (    sort { $inputSegSizeHash->{$b}->{'0_0_5_seq___length'} <=> $inputSegSizeHash->{$a}->{'0_0_5_seq___length'} } (   keys (  %{ $inputSegSizeHash }  )   )    ){   print "\n\n\cl\cl\n\nOutSide of if\n";
    if   (  defined ( $inputSegSizeHash->{$BigNumberKey} )  ){              
      if (   (  defined ( $inputSegSizeHash->{$BigNumberKey}->{'0_0_5_seq___length'} )  ) && ( $inputSegSizeHash->{$BigNumberKey}->{'0_0_5_seq___length'} > 0    )){
      
        my $bigNumberSize=$inputSegSizeHash->{$BigNumberKey}->{'0_0_5_seq___length'}; print "\$oPsKgIdx=$oPsKgIdx\t\$bigNumberSize=$bigNumberSize=\$inputSegSizeHash->{\$BigNumberKey}->{'0_0_5_seq___length'}=\$inputSegSizeHash->{$BigNumberKey}->{'0_0_5_seq___length'}\n";
        if ($bigNumberSize > $sizeLimit){ DieWork::Just_dieWork( $die_MsgHead."\n \$bigNumberSize=$bigNumberSize should <= \$sizeLimit=$sizeLimit !!  $!\n\n\n".$caller_inform );  }
        push @{ $outPutSegKeyGroups->[$oPsKgIdx] },$BigNumberKey;   
        #print "\n\n\nprint Dumper (\$outPutSegKeyGroups)\n****************Big-Out*********************\n\cl";     DirFileHandle::PrintAndWarnDumper ($outPutSegKeyGroups);
        #print "\n\n\nprint Dumper (\$inputSegSizeHash)\n*****************Big__In***BeforeDelete*****\n\cl";    DirFileHandle::PrintAndWarnDumper ($inputSegSizeHash); 
        delete $inputSegSizeHash->{$BigNumberKey};
        #print "\nprint Dumper (\$inputSegSizeHash)\n*****************Big__In****After-Delete****\n\cl";    DirFileHandle::PrintAndWarnDumper ($inputSegSizeHash);
        my $addedNumber=$bigNumberSize; print "\$addedNumber=\$bigNumberSize=$addedNumber\n";#这个数字用来将 长度加和 以计算是否超过限度
        SFOREACHMARK: foreach my $SmallNBKey (    sort { $inputSegSizeHash->{$a}->{'0_0_5_seq___length'} <=> $inputSegSizeHash->{$b}->{'0_0_5_seq___length'} } (   keys (  %{ $inputSegSizeHash }  )   )    ){
          $addedNumber+=$inputSegSizeHash->{$SmallNBKey}->{'0_0_5_seq___length'};  print "\$addedNumber=$addedNumber \t\$addedNumber+=\$inputSegSizeHash->{\$SmallNBKey}->{'0_0_5_seq___length'}=\$inputSegSizeHash->{$SmallNBKey}->{'0_0_5_seq___length'}=$inputSegSizeHash->{$SmallNBKey}->{'0_0_5_seq___length'}\n";
          if ($addedNumber > $sizeLimit){                  print ">>>>>>>>>>\$addedNumber=$addedNumber > \$sizeLimit=$sizeLimit\n\cl\n";     
            last SFOREACHMARK;
          }
          else{                                            print "<=<=<=<=<=<=\$addedNumber=$addedNumber <= \$sizeLimit=$sizeLimit\n";           
          	push @{ $outPutSegKeyGroups->[$oPsKgIdx] },$SmallNBKey;
          	#print "\n\n\nprint Dumper (\$outPutSegKeyGroups)\n*****************Sml-Out********************\n\cl";         DirFileHandle::PrintAndWarnDumper ($outPutSegKeyGroups);
          	delete $inputSegSizeHash->{$SmallNBKey};
          	#print "\n\n\nprint Dumper (\$inputSegSizeHash)\n******************Sml__In*******************\n\cl";        DirFileHandle::PrintAndWarnDumper ($inputSegSizeHash);
          }
          
        }
        $oPsKgIdx++;
      
      
      }
      else{ DieWork::Just_dieWork( $die_MsgHead."\n \$inputSegSizeHash->{$BigNumberKey}->{'0_0_5_seq___length'}=$inputSegSizeHash->{$BigNumberKey}->{'0_0_5_seq___length'} should > 0  !!  $!\n\n\n".$caller_inform );   }
    }
    else{ 
      #Just been deleted
    }
    
  }
  
  return $outPutSegKeyGroups;
  
}


# my $outPutSegKeyGroups = FastaChunkWork::disIDsIntoGrpups_987654321($input_SegSize_Hash, $sizeLimit);
#            disIDsIntoGrpups_987654321{   #将很多序列，按照文件大小的要求，放入一系列小于该大小要求的数组中，该数组可用于下一步的文件生成,这里简单的按照大小顺序排列，最大的最早放入，最小的最后放入，这样更加有利于blastx的最佳序列数量的定义
sub disIDsIntoGrpups_987654321{   #将很多序列，按照文件大小的要求，放入一系列小于该大小要求的数组中，该数组可用于下一步的文件生成,这里简单的按照大小顺序排列，最大的最早放入，最小的最后放入，这样更加有利于blastx的最佳序列数量的定义
  my ($input_SegSize_Hash, $sizeLimit)=@_;
  
  my $warnMsgBody="\nIn package  FastaChunkWork,\tIn sub disIDsIntoGrpups_987654321,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  my $caller_inform=DirFileHandle::print_SubCallerInform;
	
  if (   (  defined ( $input_SegSize_Hash )  ) && (  ref( $input_SegSize_Hash ) eq 'HASH' )   ){			}
	else{ DieWork::Just_dieWork( $die_MsgHead."\n \$input_SegSize_Hash=$input_SegSize_Hash should be a HASH ref !!  $!\n\n\n".$caller_inform );   }
	
  my $inputSegSizeHash=Storable::dclone (  $input_SegSize_Hash  );
  
  my $outPutSegKeyGroups;
  my $oPsKgIdx=0;
  foreach my $BigNumberKey (    sort { $inputSegSizeHash->{$b}->{'0_0_5_seq___length'} <=> $inputSegSizeHash->{$a}->{'0_0_5_seq___length'} } (   keys (  %{ $inputSegSizeHash }  )   )    ){   print "\n\n\cl\cl\n\nOutSide of if\n";
    if   (  defined ( $inputSegSizeHash->{$BigNumberKey} )  ){              
      if (   (  defined ( $inputSegSizeHash->{$BigNumberKey}->{'0_0_5_seq___length'} )  ) && ( $inputSegSizeHash->{$BigNumberKey}->{'0_0_5_seq___length'} > 0    )){
      
        my $bigNumberSize=$inputSegSizeHash->{$BigNumberKey}->{'0_0_5_seq___length'}; print "\$oPsKgIdx=$oPsKgIdx\t\$bigNumberSize=$bigNumberSize=\$inputSegSizeHash->{\$BigNumberKey}->{'0_0_5_seq___length'}=\$inputSegSizeHash->{$BigNumberKey}->{'0_0_5_seq___length'}\n";
        if ($bigNumberSize > $sizeLimit){ DieWork::Just_dieWork( $die_MsgHead."\n \$bigNumberSize=$bigNumberSize should <= \$sizeLimit=$sizeLimit !!  $!\n\n\n".$caller_inform );  }
        push @{ $outPutSegKeyGroups->[$oPsKgIdx] },$BigNumberKey;   
        #print "\n\n\nprint Dumper (\$outPutSegKeyGroups)\n****************Big-Out*********************\n\cl";     DirFileHandle::PrintAndWarnDumper ($outPutSegKeyGroups);
        #print "\n\n\nprint Dumper (\$inputSegSizeHash)\n*****************Big__In***BeforeDelete*****\n\cl";    DirFileHandle::PrintAndWarnDumper ($inputSegSizeHash); 
        delete $inputSegSizeHash->{$BigNumberKey};
        #print "\nprint Dumper (\$inputSegSizeHash)\n*****************Big__In****After-Delete****\n\cl";    DirFileHandle::PrintAndWarnDumper ($inputSegSizeHash);
        my $addedNumber=$bigNumberSize; print "\$addedNumber=\$bigNumberSize=$addedNumber\n";#这个数字用来将 长度加和 以计算是否超过限度
        SFOREACHMARK: foreach my $SmallNBKey (    sort { $inputSegSizeHash->{$b}->{'0_0_5_seq___length'} <=> $inputSegSizeHash->{$a}->{'0_0_5_seq___length'} } (   keys (  %{ $inputSegSizeHash }  )   )    ){
          $addedNumber+=$inputSegSizeHash->{$SmallNBKey}->{'0_0_5_seq___length'};  print "\$addedNumber=$addedNumber \t\$addedNumber+=\$inputSegSizeHash->{\$SmallNBKey}->{'0_0_5_seq___length'}=\$inputSegSizeHash->{$SmallNBKey}->{'0_0_5_seq___length'}=$inputSegSizeHash->{$SmallNBKey}->{'0_0_5_seq___length'}\n";
          if ($addedNumber > $sizeLimit){                  print ">>>>>>>>>>\$addedNumber=$addedNumber > \$sizeLimit=$sizeLimit\n\cl\n";     
            last SFOREACHMARK;
          }
          else{                                            print "<=<=<=<=<=<=\$addedNumber=$addedNumber <= \$sizeLimit=$sizeLimit\n";           
          	push @{ $outPutSegKeyGroups->[$oPsKgIdx] },$SmallNBKey;
          	#print "\n\n\nprint Dumper (\$outPutSegKeyGroups)\n*****************Sml-Out********************\n\cl";         DirFileHandle::PrintAndWarnDumper ($outPutSegKeyGroups);
          	delete $inputSegSizeHash->{$SmallNBKey};
          	#print "\n\n\nprint Dumper (\$inputSegSizeHash)\n******************Sml__In*******************\n\cl";        DirFileHandle::PrintAndWarnDumper ($inputSegSizeHash);
          }
          
        }
        $oPsKgIdx++;
      
      
      }
      else{ DieWork::Just_dieWork( $die_MsgHead."\n \$inputSegSizeHash->{$BigNumberKey}->{'0_0_5_seq___length'}=$inputSegSizeHash->{$BigNumberKey}->{'0_0_5_seq___length'} should > 0  !!  $!\n\n\n".$caller_inform );   }
    }
    else{ 
      #Just been deleted
    }
    
  }
  
  return $outPutSegKeyGroups;
  
}


1;

##########################################################################################################################################
# 
