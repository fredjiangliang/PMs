#!/usr/bin/perl -w
BEGIN{  push (@INC, '/home/fredjiang/fredProgs/PMs')  };


use strict;
use warnings;

use Bio::Seq;
use Bio::Tools::CodonTable;
use Bio::Index::Fasta;
use BlastHandle;
use DirFileHandle;
use GeneWiseHandle;
use ProSplignHandle;
############################################################################################################

#         
############################################################################################################




package FastaFileHandle;


my ($ctg_stt_key,     $ctg_end_key,    $ctg_len_key,      $old_stt_key,     $old_end_key,     $old_DNA_key,     $old_PEP_key,   $new_stt_key,     $new_end_key,     $new_DNA_key,     $new_PEP_key,     $newDNA1lkey,             $newPEP1lkey            )
  =('7_0_ctg_0_stt',  '7_0_ctg_1_end', '7_0_ctg_2_len',   '7_1_old_0_stt',  '7_1_old_1_end', '7_1_old_2_DNA',  '7_1_old_3_PEP', '7_1_new_0_stt',  '7_1_new_1_end',  '7_1_new_2_DNA',  '7_1_new_3_PEP',  '7_1_new_4_oneLine_DNA',  '7_1_new_5_oneLine_PEP' );
#   ctgStt             ctgEnd          ctgLeth            old_stt            old_end                                             new_stt           new_end           newString 
  
my ($ups_stt_key,    $ups_end_key,     $ups_DNA_key,      $ups_PEP_key,     $dws_stt_key,    $dws_end_key,     $dws_DNA_key,      $dws_PEP_key   )
  =('8_0_ups_0_stt', '8_0_ups_1_end',  '8_0_ups_2_DNA',  '8_0_ups_3_PEP',   '8_1_dws_0_stt', '8_1_dws_1_end',  '8_1_dws_2_DNA',  '8_1_dws_3_PEP' );
#                                                                     
 
 my ($utr5_st_key,    $utr5_ed_key,     $utr5_DNA_key,       $utr3_st_key,   $utr3_ed_key,     $utr3_DNA_key,     $utr3_len_key,    )
  =('9_0_utr5_0_st', '9_0_utr5_1_ed',  '9_0_utr5_2_DNA',    '9_1_utr3_0_st', '9_1_utr3_1_ed',  '9_1_utr3_2_DNA', '9_1_utr3_3_len'   );
#                                                            utrStt           utrEnd            utrSequence        putative3UTRLength      
  
my ($upsATG_fKey,         $upsATGfMakK,          $upsSTP_fKey,         $upsSTPfMakK,          $dwsSTP_fKey,         $dwsSTPfMakK         )
  =('A_0_upsATG_0_found', 'A_0_upsATG_1_fdMak',  'A_1_upsSTP_0_found', 'A_1_upsSTP_1_fdMak',  'A_2_dwsSTP_0_found', 'A_2_dwsSTP_1_fdMak',);
#    upS__ATG_found         upS__ATG_found_mark   upS_stop_found         upS_stop_found_mark   dwS_stop_found        dwS_stop_found_mark



#  FastaFileHandle::BuildIdxFile_for_fastaFile($inFastaFile, $index_file);
sub BuildIdxFile_for_fastaFile{ #为fasta文件建立 idx文件
	my ($inFastaFile, $index_file)=@_;
	
	my $warnMsgBody="\nIn package  FastaFileHandle,\tIn sub BuildIdxFile_for_fastaFile,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  my $caller_inform=DirFileHandle::print_SubCallerInform;
  
  if   (   (  defined ( $inFastaFile )  ) && ( $inFastaFile=~m/\S+/ ) && (  ( -e $inFastaFile )  )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$inFastaFile=$inFastaFile should be a defined not empty file  !!  $!\n\n\n".$caller_inform ); 	}
  
  if   (   (  defined ( $inFastaFile )  ) && ( $inFastaFile=~m/\S+/ )   ){}
  else{
  	$index_file=$inFastaFile."Index.txt";	
  } 
	
  my $inx_obj = Bio::Index::Fasta->new(-filename   => $index_file,
                                       -write_flag => 1            );
	  	
	$inx_obj->make_index($inFastaFile);	
	
}

#my $seq_IO   =  FastaFileHandle::feach_seqOBj_from_idx_File($prim_ID, $index_file);
sub feach_seqOBj_from_idx_File{  #从index文件中 获取Seq对象
	my ($prim_ID, $index_file)=@_;
	
	my $warnMsgBody="\nIn package  FastaFileHandle,\tIn sub feach_seqOBj_from_idx_File,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	if   (   (  defined ( $index_file )  ) && ( $index_file=~m/\S+/ ) && (  ( -e $index_file )  )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$index_file=$index_file should be a defined not empty file  !!  $!\n\n\n".$caller_inform ); 	}
	
	 if   (   (  defined ( $prim_ID )  ) && ( $prim_ID=~m/\S+/ )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$prim_ID=$prim_ID should be a defined not empty string  !!  $!\n\n\n".$caller_inform ); 	}
	
	
	my $inx_obj  = Bio::Index::Fasta->new(-filename   => $index_file);
		  
  my $seq_IO   = $inx_obj->fetch($prim_ID);
  
  return $seq_IO;
		  
}


#my $seq_string   =  FastaFileHandle::feach_seqString_from_idx_File($prim_ID, $index_file 【, $stt, $end】);
sub feach_seqString_from_idx_File{  #从index文件中 获取 纯序列
	my ($prim_ID, $index_file, $stt, $end)=@_;
	
	my $warnMsgBody="\nIn package  FastaFileHandle,\tIn sub feach_seqString_from_idx_File,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	if   (   (  defined ( $index_file )  ) && ( $index_file=~m/\S+/ ) && (  ( -e $index_file )  )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$index_file=$index_file should be a defined not empty file  !!  $!\n\n\n".$caller_inform ); 	}
	
	if   (   (  defined ( $prim_ID )  ) && ( $prim_ID=~m/\S+/ )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$prim_ID=$prim_ID should be a defined not empty string  !!  $!\n\n\n".$caller_inform ); 	}
	
	
	
	#warn "\n20190926-0-0-1 \$prim_ID=$prim_ID \$index_file=$index_file\n";	 print "\n20190926-0-0-1 \$prim_ID=$prim_ID \$index_file=$index_file\n";	 
  my $seq_IO    = FastaFileHandle::feach_seqOBj_from_idx_File($prim_ID, $index_file);
  my $seq_string;
  
  if   (    (  defined ( $stt )  )  ||  (  defined ( $end )  )    ){
  	if   (    (   (  defined ( $stt )  ) && ( $stt=~m/\d+/ ) && ( $stt > 0 )   ) && (   (  defined ( $end )  ) && ( $end=~m/\d+/ ) && ( $end > 0 )   )    ){
  		if ( $stt <= $end ){ #表示  两个数字 $stt 和 $end ， 都被全部正确赋值的时候, 且$stt <= $end 成立时，则进行 子序列的抽提
  			$seq_string=$seq_IO->subseq ( $stt, $end );
  		}
  		else{ #表示 两个数字 $stt 和 $end ， 都被全部正确赋值的时候，$stt <= $end 必须成立
  			DieWork::Just_dieWork( $die_MsgHead."\n \$stt=$stt should <= \$end=$end   !!  $!\n\n\n".$caller_inform );
  		}
  	}
  	else{ #表示 两个数字 $stt 和 $end ， 当至少有一个被定义的时候，如没有被全部正确赋值，则 die
  		DieWork::Just_dieWork( $die_MsgHead."\n \$stt=$stt \$end=$end should be all not defined or all be right defined   !!  $!\n\n\n".$caller_inform );
  	}
  }
  else{ #表示 两个数字 $stt 和 $end ， 都没有被定义的时候，则直接取该序列的全长
  	$seq_string=$seq_IO->seq(); 
  }
  
  $seq_string=~s/\s+//g;
  
  return $seq_string;
		  
}






sub ExpandORFforcDNA_with_ZorF{

  my ($inCdnaContig, $smallOrfStt, $smallOrfEnd, $stran_ZoF, $IN_bestProSplgnAlnHash, $printInform)=@_; 
  
  my $WarnHeadMsg="In package FastaFileHandle,		\n In sub ExpandORFforcDNA_with_ZorF,\n"; 
  my $dieHeadMsg="\n\n\nDIE!!!!!!\n$WarnHeadMsg";	$WarnHeadMsg="\n\n\n$WarnHeadMsg";	 
  print $WarnHeadMsg;                             print "     20562056-a-0 \$printInform=$printInform $smallOrfStt, $smallOrfEnd, $stran_ZoF\n\n" if ( (defined ($printInform)));
  
  
  
  $inCdnaContig=~s/\s//g; $inCdnaContig=~s/\n//g;
  my $ctgLeth=length $inCdnaContig;
  
  my $expandWithOutZF_hash;
  my $ctgStt;
  my $ctgEnd;
  if     ($stran_ZoF eq '+'){   
    $expandWithOutZF_hash=  &ExpandORFforcDNA ($inCdnaContig, $smallOrfStt, $smallOrfEnd, 'ForSEC', $printInform) ; 
    $ctgStt=1;
    $ctgEnd        =$expandWithOutZF_hash->{$ctg_end_key            };
  
  }
  elsif  ($stran_ZoF eq '-'){   
    
    my ($tpStt, $tpEnd)=($smallOrfStt, $smallOrfEnd);
    $smallOrfStt=SeqSegmentsTools::getbigOne    ($tpStt, $tpEnd);
    $smallOrfEnd=SeqSegmentsTools::getSmallOne  ($tpStt, $tpEnd);
    
    my $revreseCdnaCtg=reverse  $inCdnaContig; $revreseCdnaCtg=~tr/acgtACGT/tgcaTGCA/;
    my $newStt=$ctgLeth-$smallOrfStt+1;                                                   print "\n20562056-a-1 \$printInform=$printInform  $smallOrfStt, $smallOrfEnd, $stran_ZoF \$newStt=\$ctgLeth-\$smallOrfStt+1=$ctgLeth-$smallOrfStt+1=$newStt\n\n" if ( (defined ($printInform)));
    my $newEnd=$ctgLeth-$smallOrfEnd+1;                                                   print "\n20562056-a-2 $printInform $smallOrfStt, $smallOrfEnd, $stran_ZoF \$newEnd=\$ctgLeth-\$smallOrfEnd+1=$ctgLeth-$smallOrfEnd+1=$newEnd\n\n"if ( (defined ($printInform)));
    $expandWithOutZF_hash= &ExpandORFforcDNA ($revreseCdnaCtg, $newStt, $newEnd, 'ForSEC', $printInform);
    my $outStt     =$expandWithOutZF_hash->{$new_stt_key            };                       print "\n20562056-a-3 $printInform  $smallOrfStt, $smallOrfEnd, $stran_ZoF \$outStt=\$expandWithOutZF_hash->{$new_stt_key}=$outStt\n\n"if ( (defined ($printInform)));
    my $outEnd     =$expandWithOutZF_hash->{$new_end_key            };                       print "\n20562056-a-4 $printInform  $smallOrfStt, $smallOrfEnd, $stran_ZoF \$outEnd=\$expandWithOutZF_hash->{$new_end_key}=$outEnd\n\n"if ( (defined ($printInform)));
    $expandWithOutZF_hash->{$new_stt_key            }   =$ctgLeth-$outStt+1; print "\n20562056-a-5 $printInform $smallOrfStt, $smallOrfEnd, $stran_ZoF \$expandWithOutZF_hash->{$new_stt_key}=\$ctgLeth-\$outStt+1=$ctgLeth-$outStt+1=$expandWithOutZF_hash->{$new_stt_key}\n\n" if ( (defined ($printInform)));
    $expandWithOutZF_hash->{$new_end_key            }   =$ctgLeth-$outEnd+1; print "\n20562056-a-6 $printInform $smallOrfStt, $smallOrfEnd, $stran_ZoF \$expandWithOutZF_hash->{$new_end_key}=\$ctgLeth-\$outEnd+1=$ctgLeth-$outEnd+1=$expandWithOutZF_hash->{$new_end_key}\n\n" if ( (defined ($printInform))) ; 
    
    $ctgStt        =$expandWithOutZF_hash->{$ctg_end_key            };
    $ctgEnd=1;
    
    
  
    
    if (  defined ( $expandWithOutZF_hash->{$utr5_st_key} )  ){
      my $utr5Stt     =$expandWithOutZF_hash->{$utr5_st_key            };
      my $utr5End     =$expandWithOutZF_hash->{$utr5_ed_key            };
      $expandWithOutZF_hash->{$utr5_st_key            }   =$ctgLeth-$utr5Stt+1; 
      $expandWithOutZF_hash->{$utr5_ed_key            }   =$ctgLeth-$utr5End+1;  		
    }
    
    if (  defined ( $expandWithOutZF_hash->{$ups_stt_key} )  ){
      my $ups_stt     =$expandWithOutZF_hash->{$ups_stt_key            };
      my $ups_end     =$expandWithOutZF_hash->{$ups_end_key            };
      $expandWithOutZF_hash->{$ups_stt_key            }   =$ctgLeth-$ups_stt+1; 
      $expandWithOutZF_hash->{$ups_end_key            }   =$ctgLeth-$ups_end+1;  	
      
    }
    
    if (  defined ( $expandWithOutZF_hash->{$utr3_st_key} )  ){
      my $utr3Stt     =$expandWithOutZF_hash->{$utr3_st_key            };
      my $utr3End     =$expandWithOutZF_hash->{$utr3_ed_key            };
      $expandWithOutZF_hash->{$utr3_st_key            }   =$ctgLeth-$utr3Stt+1; 
      $expandWithOutZF_hash->{$utr3_ed_key            }   =$ctgLeth-$utr3End+1;  		
    }
    
    if (  defined ( $expandWithOutZF_hash->{$dws_stt_key} )  ){
      my $dws_stt     =$expandWithOutZF_hash->{$dws_stt_key            };
      my $dws_end     =$expandWithOutZF_hash->{$dws_end_key            };
      $expandWithOutZF_hash->{$dws_stt_key            }   =$ctgLeth-$dws_stt+1; 
      $expandWithOutZF_hash->{$dws_end_key            }   =$ctgLeth-$dws_end+1;  		
      
    }
    
  }
  else{
    my $dieMsg="$dieHeadMsg\n\$stran_ZoF=$stran_ZoF Shound be + or -:$!\n\n\n\n\n"; print $dieMsg; warn $dieMsg;
  }
  
  $expandWithOutZF_hash->{$newDNA1lkey            }=$expandWithOutZF_hash->{$new_DNA_key            } = $expandWithOutZF_hash->{$old_DNA_key            }   =$IN_bestProSplgnAlnHash->{'6_TGA_DNA1string'}; 
  $expandWithOutZF_hash->{$newPEP1lkey            }=$expandWithOutZF_hash->{$new_PEP_key            } = $expandWithOutZF_hash->{$old_PEP_key            }   =$IN_bestProSplgnAlnHash->{'6_TGA_DNA0tr_0Sq'}; 
   
  if (  defined ( $expandWithOutZF_hash->{$ups_DNA_key} )  ){
    $expandWithOutZF_hash->{$new_DNA_key            }   =$expandWithOutZF_hash->{$ups_DNA_key}."\n".$expandWithOutZF_hash->{$new_DNA_key}; 
    $expandWithOutZF_hash->{$new_PEP_key            } 	=$expandWithOutZF_hash->{$ups_PEP_key}."\n".$expandWithOutZF_hash->{$new_PEP_key}; 
     
    $expandWithOutZF_hash->{$newDNA1lkey            }   =$expandWithOutZF_hash->{$ups_DNA_key}."  ".$expandWithOutZF_hash->{$newDNA1lkey}; 
    $expandWithOutZF_hash->{$newPEP1lkey            } 	=$expandWithOutZF_hash->{$ups_PEP_key}."  ".$expandWithOutZF_hash->{$newPEP1lkey};
  }
  if (  defined ( $expandWithOutZF_hash->{$dws_DNA_key} )  ){
    $expandWithOutZF_hash->{$new_DNA_key            }   =$expandWithOutZF_hash->{$new_DNA_key}."\n".$expandWithOutZF_hash->{$dws_DNA_key}; 
    $expandWithOutZF_hash->{$new_PEP_key            } 	=$expandWithOutZF_hash->{$new_PEP_key}."\n".$expandWithOutZF_hash->{$dws_PEP_key}; 
    
    $expandWithOutZF_hash->{$newDNA1lkey            }   =$expandWithOutZF_hash->{$newDNA1lkey}."  ".$expandWithOutZF_hash->{$dws_DNA_key}; 
    $expandWithOutZF_hash->{$newPEP1lkey            } 	=$expandWithOutZF_hash->{$newPEP1lkey}."  ".$expandWithOutZF_hash->{$dws_PEP_key}; 
  }
  
  $expandWithOutZF_hash->{$ctg_stt_key            }   =$ctgStt; 
  $expandWithOutZF_hash->{$ctg_end_key            }   =$ctgEnd;  	
  
  $expandWithOutZF_hash->{$old_stt_key           }   =$smallOrfStt; 
  $expandWithOutZF_hash->{$old_end_key           }   =$smallOrfEnd;  
  
  $expandWithOutZF_hash->{'0_SgmtHead'        }   =$expandWithOutZF_hash->{$new_stt_key            };
  $expandWithOutZF_hash->{'1_SgmtTail'        }   =$expandWithOutZF_hash->{$new_end_key            };  print "\n20562056-a-7 $printInform \$expandWithOutZF_hash->{'1_SgmtTail'}=$expandWithOutZF_hash->{'1_SgmtTail'}   =\$expandWithOutZF_hash->{$new_end_key}=$expandWithOutZF_hash->{$new_end_key}\n\n" if ( (defined ($printInform))) ; 
  $expandWithOutZF_hash->{'2_SgmtZorF'        }   =$stran_ZoF;  
  $expandWithOutZF_hash->{'3_SgmtType'        }   ='CDS_region';  	  	                
  
  my $dieMsg2="$dieHeadMsg\n20562056   $smallOrfStt, $smallOrfEnd, $stran_ZoF \$expandWithOutZF_hash->{'2_SgmtZorF'}=$expandWithOutZF_hash->{'2_SgmtZorF'}"; 
  if ( $expandWithOutZF_hash->{'2_SgmtZorF'} eq '+' ){
    if ( $expandWithOutZF_hash->{$new_stt_key} >  $expandWithOutZF_hash->{$old_stt_key} ) { $dieMsg2.=" die \$expandWithOutZF_hash->{$new_stt_key}=$expandWithOutZF_hash->{$new_stt_key} > $expandWithOutZF_hash->{$old_stt_key}=\$expandWithOutZF_hash->{$old_stt_key}\n\n\n";  print $dieMsg2; warn $dieMsg2;  }
    if ( $expandWithOutZF_hash->{$new_end_key} <  $expandWithOutZF_hash->{$old_end_key} ) { $dieMsg2.=" die \$expandWithOutZF_hash->{$new_end_key}=$expandWithOutZF_hash->{$new_end_key} < $expandWithOutZF_hash->{$old_end_key}=\$expandWithOutZF_hash->{$old_end_key}\n\n\n";  print $dieMsg2; warn $dieMsg2;  }
  }
  if ( $expandWithOutZF_hash->{'2_SgmtZorF'} eq '-' ){
    if ( $expandWithOutZF_hash->{$new_stt_key} <  $expandWithOutZF_hash->{$old_stt_key} ) { $dieMsg2.=" die \$expandWithOutZF_hash->{$new_stt_key}=$expandWithOutZF_hash->{$new_stt_key} < $expandWithOutZF_hash->{$old_stt_key}=\$expandWithOutZF_hash->{$old_stt_key}\n\n\n";  print $dieMsg2; warn $dieMsg2;  }
    if ( $expandWithOutZF_hash->{$new_end_key} >  $expandWithOutZF_hash->{$old_end_key} ) { $dieMsg2.=" die \$expandWithOutZF_hash->{$new_end_key}=$expandWithOutZF_hash->{$new_end_key} > $expandWithOutZF_hash->{$old_end_key}=\$expandWithOutZF_hash->{$old_end_key}\n\n\n";  print $dieMsg2; warn $dieMsg2;  }
  }
  
  return $expandWithOutZF_hash;

}


sub ExpandORFforBacteriaGENOME{  # my $expandWithOutZF_hash=FastaFileHandle::ExpandORFforBacteriaGENOME  ($genomeFileIdxPath, $genomeID, $smallOrfStt, $smallOrfEnd,  $smallOrfZhengFu, $printInform);
	
	my ($warnMsgBody, $warnMsgHead, $die_MsgHead, $caller_inform)=@{ DieWork::BuildWarnDieHeadInform( 'FastaFileHandle', 'ExpandORFforBacteriaGENOME' ) };
	
	my ($genomeFileIdxPath, $genomeID, $smallOrfStt, $smallOrfEnd,  $smallOrfZhengFu, $printInform)=@_;  
  DieWork::Check_FileDirExist_or_DIE      ( $genomeFileIdxPath,       "\$genomeFileIdxPath",       $die_MsgHead, $caller_inform  );  
	DieWork::Check_DfdNoEmptString_or_DIE   ( $genomeID,                "\$genomeID",                $die_MsgHead, $caller_inform  );  
	DieWork::Check_DfdNoEmptNUMBER_or_DIE   ( $smallOrfStt,             "\$smallOrfStt",             $die_MsgHead, $caller_inform  );  
	DieWork::Check_DfdNoEmptNUMBER_or_DIE   ( $smallOrfEnd,             "\$smallOrfEnd",             $die_MsgHead, $caller_inform  );  
	DieWork::Check_DfdNoEmptString_or_DIE   ( $smallOrfZhengFu,         "\$smallOrfZhengFu",         $die_MsgHead, $caller_inform  );  
	
	my $inCdnaContig=FastaFileHandle::feach_seqString_from_idx_File($genomeID, $genomeFileIdxPath );  #【, $stt, $end】
	$inCdnaContig=~s/\s//g; $inCdnaContig=~s/\n//g;
  my $ctgLeth=length $inCdnaContig;
	
	my $expandWithOutZF_hash;
	if ($smallOrfZhengFu eq '+'){
		$expandWithOutZF_hash=  &ExpandORFforcDNA ($inCdnaContig, $smallOrfStt, $smallOrfEnd, 'ForSEC', $printInform) ; 
	}
	elsif($smallOrfZhengFu eq '-'){			
	  my ($tpStt, $tpEnd)=($smallOrfStt, $smallOrfEnd);
    $smallOrfStt=SeqSegmentsTools::getbigOne    ($tpStt, $tpEnd);
    $smallOrfEnd=SeqSegmentsTools::getSmallOne  ($tpStt, $tpEnd);
    my $revreseCdnaCtg=FastaFileHandle::ReverseComplementString ($inCdnaContig);
    my $newStt=$ctgLeth-$smallOrfStt+1;                                                   #print "\n20562056-a-1 \$printInform=$printInform  $smallOrfStt, $smallOrfEnd, $stran_ZoF \$newStt=\$ctgLeth-\$smallOrfStt+1=$ctgLeth-$smallOrfStt+1=$newStt\n\n" if ( (defined ($printInform)));
    my $newEnd=$ctgLeth-$smallOrfEnd+1;   
    $expandWithOutZF_hash= &ExpandORFforcDNA ($revreseCdnaCtg, $newStt, $newEnd, 'ForSEC', $printInform);
    my $outStt     =$expandWithOutZF_hash->{$new_stt_key            };                            #print "\n20562056-a-3 $printInform  $smallOrfStt, $smallOrfEnd, $stran_ZoF \$outStt=\$expandWithOutZF_hash->{$new_stt_key}=$outStt\n\n"if ( (defined ($printInform)));
    my $outEnd     =$expandWithOutZF_hash->{$new_end_key            };                            #print "\n20562056-a-4 $printInform  $smallOrfStt, $smallOrfEnd, $stran_ZoF \$outEnd=\$expandWithOutZF_hash->{$new_end_key}=$outEnd\n\n"if ( (defined ($printInform)));
    $expandWithOutZF_hash->{$new_stt_key            }   =$ctgLeth-$outStt+1;                      #print "\n20562056-a-5 $printInform $smallOrfStt, $smallOrfEnd, $stran_ZoF \$expandWithOutZF_hash->{$new_stt_key}=\$ctgLeth-\$outStt+1=$ctgLeth-$outStt+1=$expandWithOutZF_hash->{$new_stt_key}\n\n" if ( (defined ($printInform)));
    $expandWithOutZF_hash->{$new_end_key            }   =$ctgLeth-$outEnd+1;                      #print "\n20562056-a-6 $printInform $smallOrfStt, $smallOrfEnd, $stran_ZoF \$expandWithOutZF_hash->{$new_end_key}=\$ctgLeth-\$outEnd+1=$ctgLeth-$outEnd+1=$expandWithOutZF_hash->{$new_end_key}\n\n" if ( (defined ($printInform))) ; 
    
	}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$smallOrfZhengFu=$smallOrfZhengFu should be + or -  !!  $!\n\n\n".$caller_inform ); 	}

  

	return $expandWithOutZF_hash;
	
  
}

sub ExpandORFforcDNA{  #input : 1, a + strand cDNA contig; 2, a partial orf start position; 3, a partial orf end postion. Find expanded ORF
  my $WarnHeadMsg="In package FastaFileHandle,		\n In sub ExpandORFforcDNA,\n"; my $dieHeadMsg="\n\n\nDIE!!!!!!\n$WarnHeadMsg";	$WarnHeadMsg="\n\n\n$WarnHeadMsg";	 #print $WarnHeadMsg;
  
  my ($inCdnaContig,  $smallOrfStt, $smallOrfEnd, $ForSec_orNOT, $printInform)=@_;  print "\n20562056-b00 \$printInform=$printInform \$smallOrfStt=$smallOrfStt, \$smallOrfEnd=$smallOrfEnd\n\n" if ( (defined ($printInform)));
  
  my $max5UTRlength=300;
  my $max3UTRlength=300;
  
  my $maxExpandingLength=9999;  #max expanding region
  
  $inCdnaContig=~s/\s//g; $inCdnaContig=~s/\n//g;
  my $ctgLeth=length $inCdnaContig;
  
  my ($tpStt, $tpEnd)=($smallOrfStt, $smallOrfEnd);
  $smallOrfStt=SeqSegmentsTools::getSmallOne($tpStt, $tpEnd);
  $smallOrfEnd=SeqSegmentsTools::getbigOne  ($tpStt, $tpEnd);
  
  my $inATGfound=0; my $first_3_codon=substr ($inCdnaContig, $smallOrfStt-1,   3 ); if (  $first_3_codon eq 'ATG'                               ){ $inATGfound=1; } 
  my $inSTPfound=0; my $last__3_codon=substr ($inCdnaContig, $smallOrfEnd-3,   3 ); 
  if (   (  defined ( $ForSec_orNOT )  ) && ( $ForSec_orNOT eq 'ForSEC' )   ){  	if ( ($last__3_codon eq 'TAG') || ($last__3_codon eq 'TAA')                              ){ $inSTPfound=1; }   }
  else                                                                       {  	if ( ($last__3_codon eq 'TAG') || ($last__3_codon eq 'TAA') || ($last__3_codon eq 'TGA') ){ $inSTPfound=1; }   }
  
  
  print "\n20562056-b01 $printInform \$last__3_codon=$last__3_codon, \$inSTPfound=$inSTPfound\n\n" if ( (defined ($printInform)));
  
  my ($found_upSteam_0NNN_Stop_or_not, $found_upStreamStop_1NNN_or_not, $found_upstreamATG_2NNN_or_not, $found_downstream_stop_0NNN_or_not,  $found_downstream_1NNN_stop_or_not)
    =(0,                               0,                               0,                              0,                                   0                                 );
  my ($upSteam_0NNN_Stop_Pos,          $upStreamStop_1NNN_Pos,          $upstreamATG_2NNN_Pos,          $downstream_stop_0NNN_Pos,           $downstream_1NNN_stop_Pos         ); 
  my ($upSteam_0NNN_Stop_Cod,          $upStreamStop_1NNN_Cod,          $upstreamATG_2NNN_Cod,          $downstream_stop_0NNN_Cod,           $downstream_1NNN_stop_Cod         ); 
  
  my ($found_upSteamStop_or_not, $found_upstreamATG_or_not,     $found_downStreamStop_or_not )
    =(0,                         0,                             0                            );
  my ($upStreamStopPos,          $firstATGcodonPos,             $downStreamStopPos           );
  my ($upStreamStopCodon,                                       $downStreamStopCodon         );
  
  
  if ($inATGfound==0){
  	
  	                 ###
  	#1st check ------STOP------------------------------ORIGINAL_SEGMENT, find the STOP position and codon type
  	my $upStreamSTOPArrayRef   =&findTheFisrtStop__UPstream($inCdnaContig, $smallOrfStt-1, ''                      ); 
  	if ( ref ($upStreamSTOPArrayRef   ) eq 'ARRAY'  ){ ($upStreamStopPos,         $upStreamStopCodon      )=@{ $upStreamSTOPArrayRef    };   $found_upSteamStop_or_not=1; }
   	                           ###    
    #2nd check ------STOP------NNN---------------------ORIGINAL_SEGMENT, find the NNN part after STOP position 
    if ($found_upSteamStop_or_not==1){
   	  my $upStreamStop_1NNN_StPos=1; if (  ( defined ($upStreamStopPos) ) && ($upStreamStopPos=~m/^\d+$/)  ){	$upStreamStop_1NNN_StPos=$upStreamStopPos; }
      my $upStreamStop_1NNN_ArRef=&findTheFirst_NNN_UpStream ($inCdnaContig, $smallOrfStt-1, $upStreamStop_1NNN_StPos); 
      if ( ref ($upStreamStop_1NNN_ArRef) eq 'ARRAY'  ){ ($upStreamStop_1NNN_Pos,   $upStreamStop_1NNN_Cod  )=@{ $upStreamStop_1NNN_ArRef };   $found_upStreamStop_1NNN_or_not=1; }
  	} 
   	                    	                ###
   	#3rd check -----[STOP]----[NNN]-------ATG----------ORIGINAL_SEGMENT, find the ATG position 
    my $ATGfindingStart=1; 	      
    if    (  ( defined ($upStreamStop_1NNN_Pos) ) && ($upStreamStop_1NNN_Pos=~m/^\d+$/)  ){	$ATGfindingStart=$upStreamStop_1NNN_Pos; }
    elsif (  ( defined ($upStreamStopPos      ) ) && ($upStreamStopPos      =~m/^\d+$/)  ){	$ATGfindingStart=$upStreamStopPos;       }    
    my $upStm_ATG_arrRef       =&findTheFisrtATG__UPstream($inCdnaContig, $smallOrfStt-1, $ATGfindingStart);	
    if ( ref ($upStm_ATG_arrRef       ) eq 'ARRAY' ) { ($firstATGcodonPos,         $found_upstreamATG_or_not)=@{ $upStm_ATG_arrRef }; }
    
    
   	                    	                     ###
    #4th check -----[STOP]----[NNN]-------ATG--NNN----ORIGINAL_SEGMENT, find the NNN part after ATG position
    if ($found_upstreamATG_or_not==1){ 
      my $upstreamATG_2NNN_StPos=1; 
      if    (  ( defined ($firstATGcodonPos) )    && ($firstATGcodonPos     =~m/^\d+$/)  ){	$upstreamATG_2NNN_StPos=$firstATGcodonPos; }
    	my $upstreamATG_2NNN_ArRef   =&findTheFirst_NNN_UpStream($inCdnaContig, $smallOrfStt-1, $upstreamATG_2NNN_StPos); 
    	if ( ref ($upstreamATG_2NNN_ArRef) eq 'ARRAY'  ){ ($upstreamATG_2NNN_Pos,    $upstreamATG_2NNN_Cod   )=@{ $upstreamATG_2NNN_ArRef };   $found_upstreamATG_2NNN_or_not=1;  }  	  
    }
    
  	                  ###
    #5th  check ------NNN-------------------------------ORIGINAL_SEGMENT, find the NNN position ,and no STOP ATG found
    if (  ($found_upSteamStop_or_not==0) && ($found_upstreamATG_or_not==0) && ($found_upStreamStop_1NNN_or_not==0) && ($found_upstreamATG_2NNN_or_not==0)  ){
    	my $upSteam_0NNN_Stop_StPos=$smallOrfStt-1-$maxExpandingLength; if ( $upSteam_0NNN_Stop_StPos < 1) { $upSteam_0NNN_Stop_StPos=1 ;}
    	my $upSteam_0NNN_Stop_ArRef   =&findTheFirst_NNN_UpStream($inCdnaContig, $smallOrfStt-1, $upSteam_0NNN_Stop_StPos); 
    	if ( ref ($upSteam_0NNN_Stop_ArRef) eq 'ARRAY'  ){ ($upSteam_0NNN_Stop_Pos,   $upSteam_0NNN_Stop_Cod  )=@{ $upSteam_0NNN_Stop_ArRef };  $found_upSteam_0NNN_Stop_or_not=1; }  	  
    }

  }
  
  if  ($inSTPfound==0){
  	
  	
  	  	                                           ###
  	#1st check -ORIGINAL_SEGMENT-------------------STOP----------------, find the real STOP position and codon type
    my $STP_findingEND=''; 
    my $downStreamArrayRef=&findTheFisrtStopDownstream($inCdnaContig, $smallOrfEnd+1, $STP_findingEND, $printInform); 
    if ( ref ($downStreamArrayRef) eq 'ARRAY'){ ($downStreamStopPos, $downStreamStopCodon)=@{ $downStreamArrayRef }; $found_downStreamStop_or_not=1; }    print "\n20562056-b-0 $printInform \$found_downStreamStop_or_not=$found_downStreamStop_or_not\n"  if ( (defined ($printInform))); 
    
    

  	  	                                ###
  	#2nd check -ORIGINAL_SEGMENT--------NNN--------STOP----------------, find the NNN before STOP position
    if ($found_downStreamStop_or_not==1){     	    	
      my $downstream_1NNN_stop_EDpos=''; if (  ( defined ($downStreamStopPos) ) && ($downStreamStopPos=~m/^\d+$/)  ){	$downstream_1NNN_stop_EDpos=$downStreamStopPos; }
      my $downstream_1NNN_stop_ArRef=&findTheFisrt_NNN_Downstream($inCdnaContig, $smallOrfEnd+1, $downstream_1NNN_stop_EDpos); 
      if ( ref ($downstream_1NNN_stop_ArRef) eq 'ARRAY'  ){ ($downstream_1NNN_stop_Pos,   $downstream_1NNN_stop_Cod  )=@{ $downstream_1NNN_stop_ArRef };   $found_downstream_1NNN_stop_or_not=1; }  	   	  #print "\n20562056-b-1 $printInform \$found_downstream_1NNN_stop_or_not=$found_downstream_1NNN_stop_or_not\t\$found_downStreamStop_or_not=$found_downStreamStop_or_not\n";#  if ( (defined ($printInform))); 
  	}

 	  	                                  ###
  	#3rd check -ORIGINAL_SEGMENT--------NNN---------------------------, find the NNN  position, no STOP found
 	  if (  ($found_downStreamStop_or_not==0) && ($found_downstream_1NNN_stop_or_not==0)  ){     	    	
      #my $downstream_stop_0NNN_EDPos=''; 
      my $downstream_stop_0NNN_EDPos=$smallOrfEnd+1+$maxExpandingLength; if ( $downstream_stop_0NNN_EDPos > $ctgLeth) { $downstream_stop_0NNN_EDPos=$ctgLeth ;}
      my $downstream_stop_0NNN_ArRef=&findTheFisrt_NNN_Downstream($inCdnaContig, $smallOrfEnd+1, $downstream_stop_0NNN_EDPos); 
      if ( ref ($downstream_stop_0NNN_ArRef) eq 'ARRAY'  ){ ($downstream_stop_0NNN_Pos,   $downstream_stop_0NNN_Cod  )=@{ $downstream_stop_0NNN_ArRef };   $found_downstream_stop_0NNN_or_not=1; }  	   	  #print "\n20562056-b-2 $printInform \$found_downstream_stop_0NNN_or_not=$found_downstream_stop_0NNN_or_not\t\$found_downStreamStop_or_not=$found_downStreamStop_or_not\n";#  if ( (defined ($printInform))); 
  	}

  	
  }
  
 
  
  
  my ($newStt, $newEnd);
  if ($inATGfound==1){$newStt=$smallOrfStt;}
  else {
  	if    ($found_upstreamATG_2NNN_or_not    ){$newStt=$upstreamATG_2NNN_Pos+3;  }
    elsif ($found_upstreamATG_or_not         ){$newStt=$firstATGcodonPos;                                     print "\n20562056-b-1 $printInform $smallOrfStt, $smallOrfEnd ExpandORFforcDNA  \$newStt=\$firstATGcodonPos=$newStt\n"  if ( (defined ($printInform))); }
    elsif ($found_upStreamStop_1NNN_or_not   ){$newStt=$upStreamStop_1NNN_Pos+3;  }
    elsif ($found_upSteamStop_or_not         ){$newStt=$upStreamStopPos+3;                                    print "\n20562056-b-2 $printInform $smallOrfStt, $smallOrfEnd ExpandORFforcDNA  \$newStt=\$upStreamStopPos+3=$upStreamStopPos+3=$newStt\n"  if ( (defined ($printInform))); }	
    elsif ($found_downstream_stop_0NNN_or_not){$newStt=$downstream_stop_0NNN_Pos+3;  }
    else                                      {$newStt=&FindFisrtStart_withTheSameFrame( $smallOrfStt );      print "\n20562056-b-3 $printInform $smallOrfStt, $smallOrfEnd ExpandORFforcDNA  \$newStt=&FindFisrtStart_withTheSameFrame( \$smallOrfStt=$smallOrfStt )\n"  if ( (defined ($printInform))); }	
  }
  
  if ($inSTPfound==1){$newEnd=$smallOrfEnd;}
  else {
    if       ($found_downstream_1NNN_stop_or_not ){$newEnd=$downstream_1NNN_stop_Pos - 1;  } 
    elsif    ($found_downStreamStop_or_not       ){$newEnd=$downStreamStopPos+2; print "20562056-b-4 \$newEnd=$newEnd  \$inCdnaContig=$inCdnaContig\n"  if ( (defined ($printInform))); }
    elsif    ($found_downstream_stop_0NNN_or_not ){$newEnd=$downstream_stop_0NNN_Pos - 1;  } 
    else                                          {$newEnd=&FindLastEndPos_withTheSameFrame( $smallOrfEnd,  $ctgLeth);   print "\n20562056-b-5 $printInform \$newEnd=$newEnd  \$inCdnaContig=$inCdnaContig\n"  if ( (defined ($printInform)));      }	
  }
  #$newString=substr ($inCdnaContig, $newStt-1, $newEnd-$newStt+1);
  
  my $outPUThash;
  
  my $upS_stop_found_mark; if ($found_upSteamStop_or_not   ){ $upS_stop_found_mark='upStp'; if    ($found_upStreamStop_1NNN_or_not    ){$upS_stop_found_mark='upStp_NNN';} $outPUThash->{$upsSTPfMakK}=$upS_stop_found_mark; }  
  my $upS__ATG_found_mark; if ($found_upstreamATG_or_not   ){ $upS__ATG_found_mark='upATG'; if    ($found_upstreamATG_2NNN_or_not     ){$upS__ATG_found_mark='upATG_NNN';} $outPUThash->{$upsATGfMakK}=$upS__ATG_found_mark; }  
  my $dwS_stop_found_mark; if ($found_downStreamStop_or_not){ $dwS_stop_found_mark='dwStp'; if    ($found_downstream_1NNN_stop_or_not ){$dwS_stop_found_mark='dwStp_NNN';} $outPUThash->{$dwsSTPfMakK}=$dwS_stop_found_mark; }  
  if ($found_upStreamStop_1NNN_or_not    ){$found_upSteamStop_or_not   =0;   } $outPUThash->{$upsSTP_fKey     }=$found_upSteamStop_or_not;      
  if ($found_upstreamATG_2NNN_or_not     ){$found_upstreamATG_or_not   =0;   } $outPUThash->{$upsATG_fKey     }=$found_upstreamATG_or_not;     
  if ($found_downstream_1NNN_stop_or_not ){$found_downStreamStop_or_not=0;   } $outPUThash->{$dwsSTP_fKey     }=$found_downStreamStop_or_not;  
  $outPUThash->{$new_stt_key            }=$newStt;        print "\n20562056-b-6 $printInform $smallOrfStt, $smallOrfEnd ExpandORFforcDNA\$outPUThash->{$new_stt_key            }=\$newStt=$newStt\n\n"  if ( (defined ($printInform)));
  $outPUThash->{$new_end_key            }=$newEnd;        print "\n20562056-b-7 $printInform $smallOrfStt, $smallOrfEnd ExpandORFforcDNA\$outPUThash->{$new_end_key            }=\$newEnd=$newEnd\n\n"  if ( (defined ($printInform)));
  #$outPUThash->{$new_str_key            }=$newString;       
  $outPUThash->{$ctg_end_key            }=$ctgLeth;	             	             	             	             	             	         	             
  
  
  if ( $found_upstreamATG_or_not ) {
    if ($newStt>1) {
      my $utr5Stt  	             	             	             	             	         	          =$outPUThash->{$utr5_st_key             }=1;  
      my $utr5End 	             	             	             	             	         	          =$outPUThash->{$utr5_ed_key             }=$newStt-1; 
      my $utr5_len         	              	             	             	             	         	                                         =abs ($utr5End-$utr5Stt+1);    #=$outPUThash->{$utr5_len_key            }
      if ($utr5_len>$max5UTRlength){ $utr5_len 	=$max5UTRlength;}
      my $utr5Sequence 	             	             	             	             	         	      =$outPUThash->{$utr5_DNA_key            }=substr ($inCdnaContig, $utr5Stt-1, $utr5_len);
    }
    
  }
  if ($smallOrfStt>1){
    my $ups_stt    	             	             	             	         	                        =$outPUThash->{$ups_stt_key             }=$newStt;
    my $ups_end    	             	             	             	         	                        =$outPUThash->{$ups_end_key             }=$smallOrfStt-1;  
    my $ups_length    	             	             	             	         	                                                             =abs ($ups_end-$ups_stt+1);
    my $ups_DNA_seq             	             	             	         	                        =$outPUThash->{$ups_DNA_key             }=substr ($inCdnaContig, $ups_stt-1, $ups_length);
    my $ups_PEP_seq	             	             	             	         	                        =$outPUThash->{$ups_PEP_key             }=&TranslateDNAtoPep($ups_DNA_seq);  
  }  
  
  if ( $found_downStreamStop_or_not ) {
  	if ($newEnd<$ctgLeth){
      my $utr3Stt  	             	             	             	             	         	          =$outPUThash->{$utr3_st_key             }=$newEnd + 1;  
      my $utr3End 	             	             	             	             	         	          =$outPUThash->{$utr3_ed_key             }=$ctgLeth; 
      my $utr3_len 	             	             	             	             	                  	=$outPUThash->{$utr3_len_key            }=abs ($utr3End-$utr3Stt+1); 
      if ($utr3_len>$max3UTRlength){ $utr3_len 	=$max3UTRlength;}
      my $utr3Sequence 	             	             	             	             	         	      =$outPUThash->{$utr3_DNA_key            }=substr ($inCdnaContig, $utr3Stt-1, $utr3_len);
    }
  }
  if ($smallOrfEnd<$newEnd){
    my $dws_stt    	             	             	             	         	                        =$outPUThash->{$dws_stt_key             }=$smallOrfEnd+1;
    my $dws_end    	             	             	             	         	                        =$outPUThash->{$dws_end_key             }=$newEnd;  
    my $dws_length    	             	             	             	         	                                                             =abs ($dws_end-$dws_stt+1);
    my $dws_DNA_seq             	             	             	         	                        =$outPUThash->{$dws_DNA_key             }=substr ($inCdnaContig, $dws_stt-1, $dws_length);
    my $dws_PEP_seq	             	             	             	         	                        =$outPUThash->{$dws_PEP_key             }=&TranslateDNAtoPep($dws_DNA_seq);  
  }            
  
  
  
  
  return $outPUThash;


}

sub FindFisrtStart_withTheSameFrame{
  my ($inStartNB)=@_; 
  my $outStartNB;
  if ( ($inStartNB-1)%3==0 ){$outStartNB=1;}
  if ( ($inStartNB-2)%3==0 ){$outStartNB=2;}
  if ( ($inStartNB-3)%3==0 ){$outStartNB=3;}
  return $outStartNB;
}

sub FindLastEndPos_withTheSameFrame{
  my ($inEnd__NB, $inLength_NB)=@_; 
  my $outEnd__NB;
  if ( ($inLength_NB-$inEnd__NB  )%3==0 ){$outEnd__NB=$inLength_NB;  }
  if ( ($inLength_NB-$inEnd__NB-1)%3==0 ){$outEnd__NB=$inLength_NB-1;}
  if ( ($inLength_NB-$inEnd__NB-2)%3==0 ){$outEnd__NB=$inLength_NB-2;}
  return $outEnd__NB;
}


sub findTheFirst_NNN_UpStream{
	my ($instring, $upstreamStartPos, $range_lastPos)=@_;
	my $TheFirst_NNN_UpStream;
  my $findTheFirst_NNN_UpStream_or_not=0;
  my $fisrtUpStreadm_NNN_Pos=-99999999999999999999999999999999;
  my $fisrt_upstream_NNN_pos=&findTheFisrt3CodonInFrame1__UPStream($instring, 'NNN', $upstreamStartPos,$range_lastPos); if (  ( defined ($fisrt_upstream_NNN_pos) ) && ($fisrt_upstream_NNN_pos=~m/^\d+$/) ){$findTheFirst_NNN_UpStream_or_not=1; if ($fisrt_upstream_NNN_pos>=$fisrtUpStreadm_NNN_Pos){$fisrtUpStreadm_NNN_Pos=$fisrt_upstream_NNN_pos;$TheFirst_NNN_UpStream='NNN';}}
  if ($findTheFirst_NNN_UpStream_or_not){
    return [$fisrtUpStreadm_NNN_Pos, $TheFirst_NNN_UpStream];
  }
}

sub findTheFisrt_NNN_Downstream{
  my ($instring, $downstreamStartPos, $range_lastPos, $printInform)=@_;
  my $first_downstream_NNN_Codon;
  my $found_NNN_downstream_or_not=0;
  my $fisrtDownStreadm_NNN_Pos=99999999999999999999999999999999;
  my $fisrt_downstream_NNN_pos=&findTheFisrt3CodonInFrame1DownStream($instring, 'NNN', $downstreamStartPos, $range_lastPos); if (  ( defined ($fisrt_downstream_NNN_pos) ) && ($fisrt_downstream_NNN_pos=~m/^\d+$/) ){$found_NNN_downstream_or_not=1; if ($fisrt_downstream_NNN_pos<=$fisrtDownStreadm_NNN_Pos){$fisrtDownStreadm_NNN_Pos=$fisrt_downstream_NNN_pos;$first_downstream_NNN_Codon='NNN';}}
  # print "\n   20562056-c1    $printInform  \$fisrtDownStreadm_NNN_Pos=$fisrtDownStreadm_NNN_Pos, \$first_downstream_NNN_Codon=$first_downstream_NNN_Codon\n" if ( (defined ($printInform)));
  if ($found_NNN_downstream_or_not){
  return [$fisrtDownStreadm_NNN_Pos, $first_downstream_NNN_Codon];
  }
}
   
sub findTheFisrtStopDownstream{
  my ($instring, $downstreamStartPos, $range_lastPos, $printInform)=@_;
  my $first_downstream_stopCodon;
  my $foundStop_downstream_or_not=0;
  my $fisrtDownStreadmStopPos=99999999999999999999999999999999;
  my $fisrt_downstream_TGA_pos=&findTheFisrt3CodonInFrame1DownStream($instring, 'TGA', $downstreamStartPos, $range_lastPos); if (  ( defined ($fisrt_downstream_TGA_pos) ) && ($fisrt_downstream_TGA_pos=~m/^\d+$/) ){$foundStop_downstream_or_not=1; if ($fisrt_downstream_TGA_pos<=$fisrtDownStreadmStopPos){$fisrtDownStreadmStopPos=$fisrt_downstream_TGA_pos;$first_downstream_stopCodon='TGA';}}
  my $fisrt_downstream_TAG_pos=&findTheFisrt3CodonInFrame1DownStream($instring, 'TAG', $downstreamStartPos, $range_lastPos); if (  ( defined ($fisrt_downstream_TAG_pos) ) && ($fisrt_downstream_TAG_pos=~m/^\d+$/) ){$foundStop_downstream_or_not=1; if ($fisrt_downstream_TAG_pos<=$fisrtDownStreadmStopPos){$fisrtDownStreadmStopPos=$fisrt_downstream_TAG_pos;$first_downstream_stopCodon='TAG';}}
  my $fisrt_downstream_TAA_pos=&findTheFisrt3CodonInFrame1DownStream($instring, 'TAA', $downstreamStartPos, $range_lastPos); if (  ( defined ($fisrt_downstream_TAA_pos) ) && ($fisrt_downstream_TAA_pos=~m/^\d+$/) ){$foundStop_downstream_or_not=1; if ($fisrt_downstream_TAA_pos<=$fisrtDownStreadmStopPos){$fisrtDownStreadmStopPos=$fisrt_downstream_TAA_pos;$first_downstream_stopCodon='TAA';}}
  # print "\n   20562056-c1    $printInform  \$fisrtDownStreadmStopPos=$fisrtDownStreadmStopPos, \$first_downstream_stopCodon=$first_downstream_stopCodon\n" if ( (defined ($printInform)));
  if ($foundStop_downstream_or_not){
  return [$fisrtDownStreadmStopPos, $first_downstream_stopCodon];
  }
}

sub findTheFisrtStop__UPstream{
  my ($instring, $upstreamStartPos, $range_lastPos)=@_;
  my $first_upstream_stopCodon;
  my $foundStop_upstream_or_not=0;
  my $fisrtUpStreadmStopPos=-99999999999999999999999999999999;
  my $fisrt_upstream_TGA_pos=&findTheFisrt3CodonInFrame1__UPStream($instring, 'TGA', $upstreamStartPos,$range_lastPos); if (  ( defined ($fisrt_upstream_TGA_pos) ) && ($fisrt_upstream_TGA_pos=~m/^\d+$/) ){$foundStop_upstream_or_not=1; if ($fisrt_upstream_TGA_pos>=$fisrtUpStreadmStopPos){$fisrtUpStreadmStopPos=$fisrt_upstream_TGA_pos;$first_upstream_stopCodon='TGA';}}
  my $fisrt_upstream_TAG_pos=&findTheFisrt3CodonInFrame1__UPStream($instring, 'TAG', $upstreamStartPos,$range_lastPos); if (  ( defined ($fisrt_upstream_TAG_pos) ) && ($fisrt_upstream_TAG_pos=~m/^\d+$/) ){$foundStop_upstream_or_not=1; if ($fisrt_upstream_TAG_pos>=$fisrtUpStreadmStopPos){$fisrtUpStreadmStopPos=$fisrt_upstream_TAG_pos;$first_upstream_stopCodon='TAG';}}
  my $fisrt_upstream_TAA_pos=&findTheFisrt3CodonInFrame1__UPStream($instring, 'TAA', $upstreamStartPos,$range_lastPos); if (  ( defined ($fisrt_upstream_TAA_pos) ) && ($fisrt_upstream_TAA_pos=~m/^\d+$/) ){$foundStop_upstream_or_not=1; if ($fisrt_upstream_TAA_pos>=$fisrtUpStreadmStopPos){$fisrtUpStreadmStopPos=$fisrt_upstream_TAA_pos;$first_upstream_stopCodon='TAA';}}
  if ($foundStop_upstream_or_not){
  return [$fisrtUpStreadmStopPos, $first_upstream_stopCodon];
  }
}


sub findTheFisrtATG__UPstream{
  my ($instring, $in_ATGfindingStart, $range_lastPos)=@_;

  my $foundStop_upstreamATG_or_not=0;
  
  my $fisrtUpStreadmATGPos=&findTheFisrt3CodonInFrame1__UPStream($instring, 'ATG', $in_ATGfindingStart, $range_lastPos ); #print "4--20562056 \$in_ATGfindingStart=$in_ATGfindingStart, \$range_lastPos=$range_lastPos\t \$fisrtUpStreadmATGPos=&findTheFisrt3CodonInFrame1__UPStream(\$instring, 'ATG', $in_ATGfindingStart, $range_lastPos )=$fisrtUpStreadmATGPos\t\$instring=$instring\n";	
  if (  ( defined ($fisrtUpStreadmATGPos) ) && ($fisrtUpStreadmATGPos=~m/^\d+$/)  && ($fisrtUpStreadmATGPos>=$range_lastPos)  && ($fisrtUpStreadmATGPos<=$in_ATGfindingStart) ){
    $foundStop_upstreamATG_or_not=1; 
  }
  
  if ($foundStop_upstreamATG_or_not){
    return [$fisrtUpStreadmATGPos, $foundStop_upstreamATG_or_not];
  }
}


sub findTheFisrt3CodonInFrame1DownStream{
  my ($instring, $in3codon, $downstreamStartPos, $downStreamEndPos)=@_;  
  $instring=uc $instring; $in3codon=uc $in3codon;                               #print "\n\n1111111111$instring, \n$in3codon, \$downstreamStartPos=$downstreamStartPos, \$downStreamEndPos=$downStreamEndPos\n\n";  warn "\n\n1111111111$instring, \n$in3codon, \$downstreamStartPos=$downstreamStartPos, \$downStreamEndPos=$downStreamEndPos\n\n";
  my $startIDX=$downstreamStartPos-1; 
  WHMARK: while ($startIDX >=0){
    $startIDX=index ($instring, $in3codon, $startIDX);  #print "2222222222$in3codon, $instring\n\$startIDX=$startIDX\n";
    if (  ( defined ($downStreamEndPos) ) && ($downStreamEndPos=~m/^\d+$/) ){
      if ( $startIDX>($downStreamEndPos-1) ){
        last WHMARK;  
      }
    }
    if ($startIDX>=0){   
      if ( ($startIDX-$downstreamStartPos+1)%3 == 0 ){
        return ($startIDX+1);
        last WHMARK;
      
      }
      $startIDX++; 
    }
  }
}

sub findTheFisrt3CodonInFrame1__UPStream{

  my ($instring, $in3codon, $UPstreamStartPos, $UPstreamEndPos)=@_;  # Beware!! the  $UPstreamEndPos <  $UPstreamStartPos  Here!!!!
  $instring=uc $instring; $in3codon=uc $in3codon;
  my $startIDX=$UPstreamStartPos-1;
  WHMARKTWO: while ($startIDX >=0){
    $startIDX=rindex ($instring, $in3codon, $startIDX);  #print "$in3codon, $instring\t\$startIDX=$startIDX\n";
    if (  ( defined ($UPstreamEndPos) ) && ($UPstreamEndPos=~m/^\d+$/) ){
      if ( $startIDX<($UPstreamEndPos-1) ){
        last WHMARKTWO;  
      }
    }
    if ($startIDX>=0){   
      if ( ($UPstreamStartPos-2-$startIDX-1)%3 == 0 ){
      return ($startIDX+1);
      last WHMARKTWO;
      
      }
      $startIDX--; 
    }
  }
}



sub Map_PepSeq_BackTo_DNAseq{ #Input 1, protein sequence string, 2, DNA sequence string, 3, a directory path to hold middle step information. TARGET: map the protein sequence back to dna sequence, find the coding region and in frame tga positions 
  my ($inPEPstring, $inDNAstring, $analysisDIR)=@_;
  my $dieHeadMsg="\n\n\n		In package FastaFileHandle,		\n In sub Map_PepSeq_BackTo_DNAseq,\n";
  #warn  "$dieHeadMsg\n\$inPEPstring=$inPEPstring\n\n"; #\n\$inDNAstring=$inDNAstring
  print "$dieHeadMsg\nIn sub Map_PepSeq_BackTo_DNAseq,\n\$inPEPstring=$inPEPstring\n\n"; #\n\$inDNAstring=$inDNAstring
  
  $inDNAstring=~s/\n//g; $inDNAstring=~s/\s//g;
  my $DNAseqLength=length($inDNAstring);
  
  my $outInformationString;  #holding the running middle stage information
  my $outInformationHash; 
  
  my $makeblastdbPath="formatdb";                  #blastplus vesrion:  my $makeblastdbPath="~/EightT/bin/makeblastdb";
  my $tblastnExecPath="blastall";                  #blastplus vesrion:  my $tblastnExecPath="~/EightT/bin/tblastn";
  
  my $tempTimeDirNeedDel=0; if (  ( defined($analysisDIR) ) && ($analysisDIR=~/\S+/)  ) {} else {		$analysisDIR=TimeWork::GetTimeDirOrFileName();	$tempTimeDirNeedDel=1; }  #If the $analysisDIR is not set, then build a temporary directory named by time #print "\$analysisDIR=$analysisDIR\n";
  system ( "mkdir -p $analysisDIR ") ;
  
  my $blastWorkDIr="$analysisDIR/blastWorkDIr";   	system ( "mkdir -p $blastWorkDIr ") ;
  
  my $PEPseqFastaFile="$analysisDIR/blastWorkDIr/PEPseqFastaFile.txt";  $outInformationHash->{'PEPseqFastaFile'}=$PEPseqFastaFile;            
  my $DNAseqFastaFile="$analysisDIR/blastWorkDIr/DNAseqFastaFile.txt";	$outInformationHash->{'DNAseqFastaFile'}=$DNAseqFastaFile;
  
  my $PEPseqFasta=">PEP\n$inPEPstring\n\n";    open (INPEP, ">$PEPseqFastaFile") or die "$dieHeadMsg\ncannot create \$PEPseqFastaFile=$PEPseqFastaFile :$!\n\n\n"; print INPEP $PEPseqFasta;  close (INPEP);
  my $DNAseqFasta=">DNA\n$inDNAstring\n\n";    open (INDNA, ">$DNAseqFastaFile") or die "$dieHeadMsg\ncannot create \$DNAseqFastaFile=$DNAseqFastaFile :$!\n\n\n"; print INDNA $DNAseqFasta;  close (INPEP);
  
  my $tBLASTnRslt="$analysisDIR/blastWorkDIr/tBLASTnRslt.txt";        $outInformationHash->{'tBLASTnRslt'}=$tBLASTnRslt;                    	
  my $bestBlsRslt="$analysisDIR/blastWorkDIr/bestBlsRslt.html";       $outInformationHash->{'bestBlsRslt'}=$bestBlsRslt;
  my $mkDatabaseCMD_forDNA="$makeblastdbPath -i $DNAseqFastaFile -p F -o T";                                                                       #blastplus vesrion: #my $mkDatabaseCMD_forDNA="$makeblastdbPath -in $DNAseqFastaFile -dbtype nucl -parse_seqids";
  my $mkDatabaseCMD_forPEP="$makeblastdbPath -i $PEPseqFastaFile -p T -o T";
  my $db_gencode=1; 	my $U_pos_array=&Find_U_pos_inAminoAcid($inPEPstring);	if ( ref($U_pos_array) eq 'ARRAY' ){ $db_gencode=10; }
  my $blastCMD=$tblastnExecPath." -p tblastn  -i ".$PEPseqFastaFile." -d "."$DNAseqFastaFile -o $tBLASTnRslt -e 0.0000001 -D $db_gencode -F F";    #blastplus vesrion: #my $blastCMD="$tblastnExecPath –query $PEPseqFastaFile –db $DNAseqFastaFile –out $tBLASTnRslt -evalue 0.0000001 -db_gencode $db_gencode -seg no" ;
  warn "\n\$mkDatabaseCMD_forDNA=$mkDatabaseCMD_forDNA\n"; print "\n\$mkDatabaseCMD_forDNA=$mkDatabaseCMD_forDNA\n";        system ("$mkDatabaseCMD_forDNA");  
  warn "\n\$mkDatabaseCMD_forPEP=$mkDatabaseCMD_forPEP\n"; print "\n\$mkDatabaseCMD_forPEP=$mkDatabaseCMD_forPEP\n";        system ("$mkDatabaseCMD_forPEP");  
  warn "\n\$blastCMD=$blastCMD\n";                         print "\n\$blastCMD=$blastCMD\n";                                system ("$blastCMD");
  
  my $outBlastRsltHASH=BlastHandle::FindAllUForQueryAndConsAAinHit_20180626($tBLASTnRslt, $bestBlsRslt, 1, 'tblastn');    if ( ref ($outBlastRsltHASH) eq 'HASH' ) { DirFileHandle::PrintDumper ("$analysisDIR/outBlastRsltHASH.hash", $outBlastRsltHASH); }
  
  my $proSplignDir="$analysisDIR/proSplignDir";   	system ( "mkdir -p $proSplignDir ") ;
  my $porSplignHash=ProSplignHandle::RunProSplign($PEPseqFastaFile, $DNAseqFastaFile, $proSplignDir, $db_gencode);
  my $proSplign_with_AADNAposHASH=ProSplignHandle::Build_DNA_PEP_correlate_HASH($porSplignHash, $inPEPstring, $inDNAstring);
  if ( ref ($proSplign_with_AADNAposHASH) eq 'HASH' ) { DirFileHandle::PrintDumper ("$analysisDIR/proSplign_with_AADNAposHASH.hash", $proSplign_with_AADNAposHASH); }
  my $simp_proSplign_with_AADNAposHASH=ProSplignHandle::simplify_prosplign_Out_Hash($proSplign_with_AADNAposHASH); DirFileHandle::PrintDumper ("$analysisDIR/Simp_proSplign_with_AADNAposHASH.hash", $simp_proSplign_with_AADNAposHASH) if ( ref ($simp_proSplign_with_AADNAposHASH) eq 'HASH' );
  
  
  
  if (   ( ref ($proSplign_with_AADNAposHASH) eq 'HASH' ) && (  ref ( $proSplign_with_AADNAposHASH->{'0_MatchArray'} ) eq 'ARRAY'  )   ){ 
    my $howManyProSplignMatch=@{ $proSplign_with_AADNAposHASH->{'0_MatchArray'} };
    $outInformationHash->{'howManyProSplignMatch'}=$howManyProSplignMatch;
    if ( $howManyProSplignMatch==0 ){	die"\n\n\nDIE!!!!!!!!!!!!!\n$dieHeadMsg\nCheck the \$proSplignDir=$proSplignDir, something wrong is found in prosplign step, maybe no match found!!!\n\n\n";			}
    my $outControl=1; #my $outControl=$howManyProSplignMatch;
    for (  my $i=0; $i<$outControl; $i++  ){
      my $ecPrSpMatchHash=$proSplign_with_AADNAposHASH->{'0_MatchArray'}->[$i];
      my $matchPEPinputCov=  $ecPrSpMatchHash->{'5_PepCov_1dcm'};		  
      my $mcPepIptCovPecet=$outInformationHash->{'proSplign_PepCovRate'}->{$i}=$ecPrSpMatchHash->{'5_PepCov_1pct'};
      
      if (   (  defined ( $ecPrSpMatchHash->{'3_TGA_3_UCX*_MC_PCT'} )  ) && ( $ecPrSpMatchHash->{'3_TGA_3_UCX*_MC_PCT'} >= 0.5 )   ){
        
        
        if ($matchPEPinputCov<1){
          $outInformationString.="\n$dieHeadMsg\n\$matchPEPinputCov=$matchPEPinputCov\n\$mcPepIptCovPecet=$mcPepIptCovPecet\n It should be equal or bigger than 100%! 		\n\nCheck the $analysisDIR/proSplign_with_AADNAposHASH.hash.read.txt Filen\n";   
          $outInformationHash->{'Pro_quryCovRate_smaller_than_1'}->{$i}=$mcPepIptCovPecet;
        }
        
        my $ProSpChangeSeqHash=$ecPrSpMatchHash->{'6_TGAchangeHASH'};  my $ChangeSeqHash_string=DirFileHandle::PrintAndWarnDumper ($ProSpChangeSeqHash); $outInformationString.="\n\n".$ProSpChangeSeqHash."\n\n";
        
        my $TGA_Pos_array=BlastHandle::BuildTGAPosArrayFrom_TGAchangeHASH($ProSpChangeSeqHash,$DNAseqLength);
        my $ChangedContigSequ=FastaFileHandle::CHangeSpecificSequenceInAFastaSeq($inDNAstring, $ProSpChangeSeqHash );
        
        my $geneWiseWorkDIr="$analysisDIR/$i/geneWiseWorkDIr";   	system ( "mkdir -p $geneWiseWorkDIr ") ; my $geneWiseWorkDNAFIle="$geneWiseWorkDIr/GeneWiseTGAChangedDNA.txt";    my $geneWiseOUTFIle="$geneWiseWorkDIr/GeneWiseOUT.txt";  
        my $changedContig=">DNAforGeneWISE\n$ChangedContigSequ\n\n";  open (INHFASTA, ">$geneWiseWorkDNAFIle") or die "\n$dieHeadMsg\ncannot create \$changedContig=$changedContig : $!\n\n\n" ;      print INHFASTA  $changedContig;      close (INHFASTA);
        
        my $porSplignOutPEPstring=$ecPrSpMatchHash->{'6_TGA_DNA0tr_0Sq'}; print "\$porSplignOutPEPstring=$porSplignOutPEPstring\n";
        my $porSplignOutPEPseqFastaFile="$geneWiseWorkDIr/porSplignOutPEPseqFastaFile.txt";  $outInformationHash->{'porSplignOutPEPseqFastaFile'}->{$i}=$porSplignOutPEPseqFastaFile; 
        my $porSplignOutPEPseqFasta=">ProSplignOutPEP\n$porSplignOutPEPstring\n\n";    open (PROSPLINPEP, ">$porSplignOutPEPseqFastaFile") or die "$dieHeadMsg\ncannot create \$porSplignOutPEPseqFastaFile=$porSplignOutPEPseqFastaFile :$!\n\n\n"; print PROSPLINPEP $porSplignOutPEPseqFasta;  close (PROSPLINPEP);
        
        GeneWiseHandle::GeneWiseRun($porSplignOutPEPseqFastaFile, $geneWiseWorkDNAFIle, $geneWiseOUTFIle);
        
        if (-e ($geneWiseOUTFIle) ){warn "\n\n\nIn package FastaFileHandle,\nIn sub Map_PepSeq_BackTo_DNAseq,\$geneWiseOUTFIle=$geneWiseOUTFIle Exist!!\n";}else {warn "\n$dieHeadMsg\ncannot find\$geneWiseOUTFIle=$geneWiseOUTFIle:$!\n\n"; die "\n\n\nDIE!!!\n$dieHeadMsg\n\ncannot find\$geneWiseOUTFIle=$geneWiseOUTFIle:$!";}
        my $wiseHash; 
        if (   (  defined ( $ProSpChangeSeqHash )  ) && (  ref( $ProSpChangeSeqHash ) eq 'HASH'  )   ){ $wiseHash=GeneWiseHandle::NewGenewisedbPraser($geneWiseOUTFIle, '', '','', '', '', $ProSpChangeSeqHash );  }
        else                                                                                          { $wiseHash=GeneWiseHandle::NewGenewisedbPraser($geneWiseOUTFIle, '', '','','','' );  }              	print "\n$dieHeadMsg\nEnd GeneWiseHandle::NewGenewisedbPraser of \$wiseHash=$wiseHash\n";
        
        if ( ref ( $wiseHash ) eq 'HASH' ){
          my $wiseOutHashFile="$geneWiseWorkDIr/wiseOutFile.hash";      DirFileHandle::PrintDumper($wiseOutHashFile,$wiseHash) if (ref($wiseHash) eq 'HASH'); 
          
          my $bestGeneStructureARRAY;
          if ( GeneWiseHandle::CheckTheBestWiseGeneIsANotAPseudoGene($wiseHash) ){  print "$dieHeadMsg\n\nPseudo here!\$analysisDIR=$analysisDIR\n\n";     $outInformationString.="$dieHeadMsg\n\nPseudo here!\$analysisDIR=$analysisDIR\n\n\n";     #pseudogene
            $outInformationHash->{'Is_it_a_PseudoGene'}->{$i}=1;
          }
          else{                                                                 		print "$dieHeadMsg\n\nNot Pseudo here!\$analysisDIR=$analysisDIR\n\n";
            $bestGeneStructureARRAY=GeneWiseHandle::GetBestGeneStructureFromGeneWisePraserOutHash($wiseHash);
            DirFileHandle::PrintAndWarnDumper ($bestGeneStructureARRAY);                                              		#DirFileHandle::PrintDumper ( "$analysisDIR/BestGeneStructure.hash",        $bestGeneStructureARRAY ) if ( ref($bestGeneStructureARRAY) eq 'ARRAY');		
            my ($bestWiseOutDNAstring, $bestWiseOutPEPstring)=@{ GeneWiseHandle::GetBestGeneSequence_from_WisePraserOutHash($wiseHash) };
            if ($bestWiseOutPEPstring eq $inPEPstring){
              $outInformationHash->{'wiseOutPep_equl_orgInputPep'}->{$i}=1;
            }
            else {  $outInformationString.="$dieHeadMsg\nIn \$analysisDIR=$analysisDIR\nThe Input PEP string \$bestWiseOutPEPstring is different aganist the Output \$bestWiseOutPEPstring from genewise work!\n         \$inPEPstring=$inPEPstring\n\$bestWiseOutPEPstring=$bestWiseOutPEPstring\n\n\n\n"; print "$dieHeadMsg\nIn \$analysisDIR=$analysisDIR\nThe Input PEP string \$inPEPstring is different aganist the Output $bestWiseOutPEPstring from genewise work!\n\$inPEPstring=$inPEPstring\n   \$bestWiseOutPEPstring=$bestWiseOutPEPstring\n\n\n\n";
              $outInformationHash->{'wiseOutPep_equl_orgInputPep'}->{$i}=0;		  
            }
          }
          $outInformationHash->{'bestGeneStructureARRAY'}->{$i}=$bestGeneStructureARRAY;
          $outInformationHash->{'TGA_Pos_array'}->{$i}=$TGA_Pos_array;
        }
      }
    
    }
  }
  
  #the original code for blastout genewise work
  if(0){
    if (   (  defined ( $outBlastRsltHASH->{'totalHit'} )  ) && ( $outBlastRsltHASH->{'totalHit'} >=1 )   ){} else {die "\n\n\nIn package FastaFileHandle,\nIn sub Map_PepSeq_BackTo_DNAseq,\n\n\$outBlastRsltHASH->{'totalHit'}=$outBlastRsltHASH->{'totalHit'}\nIt should be bigger than 1!\nCheck the $analysisDIR\n\n\n"; }
    my $bestHitIdx= $outBlastRsltHASH->{'SortHash'}->{'EvalueSort'}->[0]->{'HitIdx'};	       my $bestQurIdx= $outBlastRsltHASH->{'SortHash'}->{'EvalueSort'}->[0]->{'QurIdx'};
    my $quryCovRate=$outBlastRsltHASH->{'_ResultArray'}->[$bestQurIdx]->{'_hitArray'}->[$bestHitIdx]->{'queryCoverRate'};    $outInformationHash->{'best_quryCovRate'}=$quryCovRate;
    if ($quryCovRate<1){ 				print "$dieHeadMsg\n\$quryCovRate=$quryCovRate\n It should be equal or bigger than 1! 		\n\nCheck the $analysisDIR/outBlastRsltHASH.hash.read.txt Filen\n";
      $outInformationString.="\n$dieHeadMsg\n\$quryCovRate=$quryCovRate\n It should be equal or bigger than 1! 		\n\nCheck the $analysisDIR/outBlastRsltHASH.hash.read.txt Filen\n";   $outInformationHash->{'best_quryCovRate_smaller_than_1'}=$quryCovRate;
    }
    
    my $ChangeSeqHash=$outBlastRsltHASH->{'_ResultArray'}->[$bestQurIdx]->{'_hitArray'}->[$bestHitIdx]->{'_TGAchangeHASH'};  my $ChangeSeqHash_string=DirFileHandle::PrintAndWarnDumper ($ChangeSeqHash); $outInformationString.="\n\n".$ChangeSeqHash_string."\n\n";
    
    my $TGA_Pos_array=BlastHandle::BuildTGAPosArrayFrom_TGAchangeHASH($ChangeSeqHash,$DNAseqLength);
    my $ChangedContigSequ=FastaFileHandle::CHangeSpecificSequenceInAFastaSeq($inDNAstring, $ChangeSeqHash );
    my $geneWiseWorkDIr="$analysisDIR/geneWiseWorkDIr";   	system ( "mkdir -p $geneWiseWorkDIr ") ; my $geneWiseWorkDNAFIle="$geneWiseWorkDIr/GeneWiseTGAChangedDNA.txt";    my $geneWiseOUTFIle="$geneWiseWorkDIr/GeneWiseOUT.txt";  
    my $changedContig=">DNAforGeneWISE\n$ChangedContigSequ\n\n";  open (INHFASTA, ">$geneWiseWorkDNAFIle") or die "\n$dieHeadMsg\ncannot create \$changedContig=$changedContig : $!\n\n\n" ;      print INHFASTA  $changedContig;      close (INHFASTA);
    GeneWiseHandle::GeneWiseRun($PEPseqFastaFile, $geneWiseWorkDNAFIle, $geneWiseOUTFIle);
    
    if (-e ($geneWiseOUTFIle) ){warn "\n\n\nIn package FastaFileHandle,\nIn sub Map_PepSeq_BackTo_DNAseq,\$geneWiseOUTFIle=$geneWiseOUTFIle Exist!!\n";}else {warn "\n$dieHeadMsg\ncannot find\$geneWiseOUTFIle=$geneWiseOUTFIle:$!\n\n"; die "\n\n\nDIE!!!\n$dieHeadMsg\n\ncannot find\$geneWiseOUTFIle=$geneWiseOUTFIle:$!";}
    my $wiseHash; 
    if (   (  defined ( $ChangeSeqHash )  ) && (  ref( $ChangeSeqHash ) eq 'HASH'  )   ){ $wiseHash=GeneWiseHandle::NewGenewisedbPraser($geneWiseOUTFIle, '', '','', '', '', $ChangeSeqHash );  }
    else                                                                                          { $wiseHash=GeneWiseHandle::NewGenewisedbPraser($geneWiseOUTFIle, '', '','','','' );  }              	print "\n$dieHeadMsg\nEnd GeneWiseHandle::NewGenewisedbPraser of \$wiseHash=$wiseHash\n";
    my $wiseOutHashFile="$geneWiseWorkDIr/wiseOutFile.hash";      DirFileHandle::PrintDumper($wiseOutHashFile,$wiseHash) if (ref($wiseHash) eq 'HASH'); 
    
    my $bestGeneStructureARRAY;
    if ( GeneWiseHandle::CheckTheBestWiseGeneIsANotAPseudoGene($wiseHash) ){  print "$dieHeadMsg\n\nPseudo here!\$analysisDIR=$analysisDIR\n\n";     $outInformationString.="$dieHeadMsg\n\nPseudo here!\$analysisDIR=$analysisDIR\n\n\n";     #pseudogene
      $outInformationHash->{'Is_it_a_PseudoGene'}=1;
    }
    else{                                                                 		print "$dieHeadMsg\n\nNot Pseudo here!\$analysisDIR=$analysisDIR\n\n";
      $bestGeneStructureARRAY=GeneWiseHandle::GetBestGeneStructureFromGeneWisePraserOutHash($wiseHash);
      DirFileHandle::PrintAndWarnDumper ($bestGeneStructureARRAY);                                              		#DirFileHandle::PrintDumper ( "$analysisDIR/BestGeneStructure.hash",        $bestGeneStructureARRAY ) if ( ref($bestGeneStructureARRAY) eq 'ARRAY');		
      my ($bestWiseOutDNAstring, $bestWiseOutPEPstring)=@{ GeneWiseHandle::GetBestGeneSequence_from_WisePraserOutHash($wiseHash) };
      if ($bestWiseOutPEPstring eq $inPEPstring){
        $outInformationHash->{'wiseOutPep_equl_orgInputPep'}=1;
      }
      else {  $outInformationString.="$dieHeadMsg\nIn \$analysisDIR=$analysisDIR\nThe Input PEP string \$bestWiseOutPEPstring is different aganist the Output \$bestWiseOutPEPstring from genewise work!\n         \$inPEPstring=$inPEPstring\n\$bestWiseOutPEPstring=$bestWiseOutPEPstring\n\n\n\n"; print "$dieHeadMsg\nIn \$analysisDIR=$analysisDIR\nThe Input PEP string \$inPEPstring is different aganist the Output $bestWiseOutPEPstring from genewise work!\n\$inPEPstring=$inPEPstring\n   \$bestWiseOutPEPstring=$bestWiseOutPEPstring\n\n\n\n";
        $outInformationHash->{'wiseOutPep_equl_orgInputPep'}=0;		  
      }
    }
  }    
  #my $seqHeadNm="WiseOUT";    #my $DNAfastaSeq=">$seqHeadNm\n$DNAptSeq\n";     my $PEPfastaSeq=">$seqHeadNm\n$bestWiseOutPEPstring\n";        
  
  
  
  
  
  
  #if ($tempTimeDirNeedDel==1){my $rmCMD="rm -f $analysisDIR";  warn "\n\$blastCMD=$blastCMD\n";                                 system ("$blastCMD");}
  return [  $outInformationString, $outInformationHash, $proSplign_with_AADNAposHASH ];

}


sub FactchSeqBy_blastdbCMD{
my ($BlsDB, $entryNM, $strand, $sttPos, $endPos)=@_;
my $blastdbcmdPath="~/EightT/bin/blastdbcmd";

my $outPutFastaString;

my $wholesqCMDFasta;
my $wholeEntryFatchCMD="$blastdbcmdPath -db $BlsDB -entry \"$entryNM\" ";
my $segMentCMDFasta;
my $getSubSegmentOrNot=0;  
my $segmentSttEndByondLimit=0;

if (  ( defined ( $sttPos ) ) || ( defined ( $sttPos ) )  ){
if (  ($sttPos=~m/^[-\+]?\d+$/) && ($endPos=~m/^[-\+]?\d+$/) && ( defined ( $strand ) ) && ( ($strand eq '+') || ($strand eq '-') )  ) {
if ($sttPos>$endPos){die "\n\n\nIn package FastaFileHandle,\nIn sub FactchSeqBy_blastdbCMD,\n\$BlsDB=$BlsDB, \$entryNM=$entryNM,\nThe \$sttPos=$sttPos should < $endPos=\$endPos\n\n\n";}
if    ($strand eq '+') { $strand ='plus';  } 
elsif ($strand eq '-') { $strand ='minus'; } 


$getSubSegmentOrNot=1;
}
else {
die "\n\n\nIn package FastaFileHandle,\nIn sub FactchSeqBy_blastdbCMD,\n\$BlsDB=$BlsDB, \$entryNM=$entryNM,\nPlease check the \$sttPos=$sttPos, \$endPos=$endPos, \$strand=$strand\n\n\n";
}
}

warn "\n\n\nNow In package FastaFileHandle,\nIn Sub FactchSeqBy_blastdbCMD\n\n"; warn "\n\$wholeEntryFatchCMD=$wholeEntryFatchCMD\n\n\n\n"; 
$wholesqCMDFasta=`$wholeEntryFatchCMD`;

if ($getSubSegmentOrNot){		  
my $wholeLgth=FastaFileHandle::GetSeqLength ($wholesqCMDFasta);    	 
if ($endPos >  $wholeLgth) {$endPos=$wholeLgth; $segmentSttEndByondLimit=1;}
if ($sttPos <= 0         ) {$sttPos=1;          $segmentSttEndByondLimit=1;}
my $segMentCMDFastaCMD="-range $sttPos-$endPos -strand $strand";  

my $SegMtblastdbCMD="$wholeEntryFatchCMD $segMentCMDFastaCMD";
warn "\n\n\nNow In package FastaFileHandle,\nIn Sub FactchSeqBy_blastdbCMD\n\n"; warn "\n\$SegMtblastdbCMD=$SegMtblastdbCMD\n\n\n\n"; 
$outPutFastaString=$segMentCMDFasta=`$SegMtblastdbCMD`;
}
else {
$outPutFastaString=$wholesqCMDFasta;
}

return [$outPutFastaString, $segmentSttEndByondLimit, $strand, $sttPos, $endPos ]
}
                              
                              
                              #  my $outString=FastaFileHandle::FactchSeqJustSubStr($LongSeqString, $strand, $sttPos, $endPos);
sub FactchSeqJustSubStr{
  my ($LongSeqString, $strand, $sttPos, $endPos)=@_;
  my $warnMstHead="\n\n\nNow In package FastaFileHandle,\nIn Sub FactchSeqJustSubStr\n\n";
  
  if ($sttPos > $endPos)  {		my $dieMsg="\n\n\n\nDIE!!!!!!!!\n$warnMstHead\n\$sttPos=$sttPos > $endPos=\$endPos\n The Stt should smaller than the End!!!!\n\n\n"; print $dieMsg; die $dieMsg; 	}
  
  $LongSeqString=~s/\n//g; $LongSeqString=~s/\s//g; $LongSeqString= uc $LongSeqString;
  my $LongLeth=length $LongSeqString;	
  if ($endPos > $LongLeth){		my $dieMsg="\n\n\n\nDIE!!!!!!!!\n$warnMstHead\n\$LongLeth=$LongLeth < $endPos=\$endPos\n The End should smaller than the whole Length!!!!\n\n\n";	 print $dieMsg; print "$strand, $sttPos, $endPos\t\$LongSeqString=$LongSeqString\n"; die $dieMsg; 	}
  
  my $outString=substr ($LongSeqString, $sttPos-1, $endPos-$sttPos+1);
  if    ($strand eq '+'){		                                                              	}
  elsif ($strand eq '-'){		$outString=reverse $outString; $outString=~tr/ACGT/TGCA/;        }
  else                  {   my $dieMsg= "\n\n\n\nDIE!!!!!!!!\n$warnMstHead\n\$strand=$strand Should be + or - !!!!\n\n\n";  print $dieMsg; die $dieMsg;  }
  
  return $outString;

}

sub FactchSeqBy_FastaCMD{
my ($BlsDB, $entryNM, $strand, $sttPos, $endPos)=@_;
my $warnMstHead="\n\n\nNow In package FastaFileHandle,\nIn Sub FactchSeqBy_FastaCMD\n\n";
#print $warnMstHead; print "1 \$BlsDB=$BlsDB\n"; print "2 \$entryNM=$entryNM\n"; print "3 \$strand=$strand\n"; print "4 \$sttPos=$sttPos\n"; print "5 \$endPos=$endPos\n"; print "\n\n\n";
#warn  $warnMstHead; warn  "1 \$BlsDB=$BlsDB\n"; warn  "2 \$entryNM=$entryNM\n"; warn  "3 \$strand=$strand\n"; warn  "4 \$sttPos=$sttPos\n"; warn  "5 \$endPos=$endPos\n"; warn  "\n\n\n";

my $fastaCMDpath="fastacmd";

my $outPutFastaString;

my $wholesqCMDFasta;
my $wholeEntryFatchCMD="$fastaCMDpath -d $BlsDB -s \"$entryNM\" ";
my $segMentCMDFasta;
my $getSubSegmentOrNot=0;  
my $segmentSttEndByondLimit=0;
my $outStrandZoF='+';

if (  ( defined ( $sttPos ) ) || ( defined ( $endPos ) )  ){
if (  ($sttPos=~m/^[-\+]?\d+$/) && ($endPos=~m/^[-\+]?\d+$/) && ( defined ( $strand ) ) && ( ($strand eq '+') || ($strand eq '-') )  ) {
if ($sttPos>$endPos){die "\n\n\nDIE!!! In package FastaFileHandle,\nIn sub FactchSeqBy_FastaCMD,\n\$BlsDB=$BlsDB, \$entryNM=$entryNM,\nThe \$sttPos=$sttPos should < $endPos=\$endPos\n\n\n";}
$outStrandZoF=$strand;
if    ($strand eq '+') { $strand =1; } 
elsif ($strand eq '-') { $strand =2; } 


$getSubSegmentOrNot=1;
}
else {
die "\n\n\nIn package FastaFileHandle,\nIn sub FactchSeqBy_FastaCMD,\n\$BlsDB=$BlsDB, \$entryNM=$entryNM,\nPlease check the \$sttPos=$sttPos, \$endPos=$endPos, \$strand=$strand\n\n\n";
}
}

#warn "\n\n\nNow In package FastaFileHandle,\nIn Sub FactchSeqBy_FastaCMD\n\n"; warn "\n\$wholeEntryFatchCMD=$wholeEntryFatchCMD\n\n\n\n"; 
print "\n\n\nNow In package FastaFileHandle,\nIn Sub FactchSeqBy_FastaCMD\n\n"; print "\n\$wholeEntryFatchCMD=$wholeEntryFatchCMD\n\n\n\n"; 
$wholesqCMDFasta=`$wholeEntryFatchCMD`;
my $wholeLgth=FastaFileHandle::GetSeqLength ($wholesqCMDFasta); 

if ($getSubSegmentOrNot){
if ($endPos >  $wholeLgth) {$endPos=$wholeLgth; $segmentSttEndByondLimit=1;}
if ($sttPos <= 0         ) {$sttPos=1;          $segmentSttEndByondLimit=1;}
my $segMentCMDFastaCMD="-L $sttPos,$endPos -S $strand";  

my $SegMtblastdbCMD="$wholeEntryFatchCMD $segMentCMDFastaCMD";
#warn "\n\n\nNow In package FastaFileHandle,\nIn Sub FactchSeqBy_FastaCMD\n\n"; warn "\n\$SegMtblastdbCMD=$SegMtblastdbCMD\n\n\n\n"; 
$outPutFastaString=$segMentCMDFasta=`$SegMtblastdbCMD`;
}
else {
$outPutFastaString=$wholesqCMDFasta;
$sttPos=1;
$endPos=$wholeLgth;
}
#print "\$segmentSttEndByondLimit=$segmentSttEndByondLimit, \$outStrandZoF=$outStrandZoF, \$sttPos=$sttPos, \$endPos=$endPos, \$wholeLgth=$wholeLgth\n\n";
return [$outPutFastaString, $segmentSttEndByondLimit, $outStrandZoF, $sttPos, $endPos, $wholeLgth ]
}



sub Find_U_pos_inAminoAcid{
  my ($inString)=@_;
  $inString=uc($inString);
  my $outArray=BlastHandle::findingIdexOfAword('U',$inString);
  return $outArray;	
}

sub Find_AA_pos_inAminoAcid{   # my $outArray=FastaFileHandle::Find_AA_pos_inAminoAcid($inString, $AA_char);
  my ($inString, $AA_char)=@_;
  
  my $warnMsgBody="\nIn package  FastaFileHandle,\tIn sub Find_AA_pos_inAminoAcid,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
  $inString=uc($inString);
  $AA_char=uc($AA_char);
  if  (   (  defined ( $inString )  ) && ( $inString=~m/\S+/ ) && (  defined ( $AA_char )  ) && ( $AA_char=~m/\S+/ )   ){
  	if (  length ( $inString ) < 1  ){ DieWork::Just_dieWork( $die_MsgHead."\$inString=$inString, it should be longer than 1\n $!\n\n\n".$caller_inform ); }
  	if (  length ( $AA_char  ) != 1  ){ DieWork::Just_dieWork( $die_MsgHead."\$AA_char=$AA_char, it should be longer than 1\n $!\n\n\n".$caller_inform ); }  	
  }
  else{
  	DieWork::Just_dieWork( $die_MsgHead."\$inString=$inString or \$AA_char=$AA_char is not right:\n $!\n\n\n".$caller_inform );
  }
  my $outArray=BlastHandle::findingIdexOfAword($AA_char,$inString);
  return $outArray;	
}

sub BuildShow_U_PEPstring{
  my ($inString, $inArray)=@_;
  my $outSting;
  my $inArrayChgedHash;
  foreach my $ec_U_Pos ( @{$inArray} ){ $inArrayChgedHash->{$ec_U_Pos}=1; print "\$inArrayChgedHash->{$ec_U_Pos}=$inArrayChgedHash->{$ec_U_Pos}\n"; }
  $inString=~s/\n//g; $inString=~s/\s//g; $inString=uc($inString);
  my $inStrLen=length $inString;
  for (my $i=0; $i<$inStrLen; $i++){
    if (   (  defined ( $inArrayChgedHash->{$i+1} )  ) && ( $inArrayChgedHash->{$i+1} == 1)   ){
      $outSting.='U';
    }
    else{
      $outSting.=' ';
    }
  }
  $outSting.="\n";
  $outSting.="$inString\n";
  return $outSting;
}

sub Build_only_Show_U_PEPstring{
  my ($inString)=@_;
  my $outSting;
  
  #$inString=~s/\n//g; $inString=~s/\s//g; 
  $inString=uc($inString);
  my $inArray=&Find_U_pos_inAminoAcid($inString);
  
  my $inArrayChgedHash;
  foreach my $ec_U_Pos ( @{$inArray} ){ $inArrayChgedHash->{$ec_U_Pos}=1; print "\$inArrayChgedHash->{$ec_U_Pos}=$inArrayChgedHash->{$ec_U_Pos}\n"; }
  
  my $inStrLen=length $inString;
  for (my $i=0; $i<$inStrLen; $i++){
    if (   (  defined ( $inArrayChgedHash->{$i+1} )  ) && ( $inArrayChgedHash->{$i+1} == 1)   ){
      $outSting.='U';
    }
    else{
      $outSting.=' ';
    }
  }
  #$outSting.="\n";
  #$outSting.="$inString\n";
  return $outSting;
}


##  my $outSting=FastaFileHandle::Build_simple_showChar_string($inWord, $inPos);  
sub Build_simple_showChar_string{ # 输入一个Char字符，和一个位置（从0开始计数），输出一个string，这个string在位置$inPos 前，都是空，到位置$inPos 打印 $inWord字符。
	my ($inWord, $inPos)=@_;
  
  my $warnMsgBody="\nIn package  FastaFileHandle,\tIn sub Build_simple_showChar_string,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  my $caller_inform=DirFileHandle::print_SubCallerInform;
  
  if   (   (  defined ( $inWord )  ) && ( $inWord=~m/^\S$/ )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$inWord=$inWord should be a defined one char string  !!  $!\n\n\n".$caller_inform ); 	}
	
	if   (   (  defined ( $inPos )  ) && ( $inPos=~m/^\d+$/ )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$inPos=$inPos should be a defined number  !!  $!\n\n\n".$caller_inform ); 	}
	
	my $outSting;
	
	for (my $i=0; $i<$inPos; $i++){
      $outSting.=' ';    
  }
  $outSting.=$inWord; 
  
  
  return $outSting;
	
	
}

#  my $outSting=FastaFileHandle::Build_show_char_string($inWord, $inArray);  
sub Build_show_char_string{  # 输入是需要寻找的某一个字符 以及含有该字符的位置数组（该数组，可由&findingIdexOfAword获得），可生成凸显该字符及位置的序列
	my ($inWord, $inArray)=@_;
  
  my $warnMsgBody="\nIn package  FastaFileHandle,\tIn sub Build_show_char_string,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  my $caller_inform=DirFileHandle::print_SubCallerInform;

  if   (   (  defined ( $inWord )  ) && ( $inWord=~m/^\S$/ )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$inWord=$inWord should be a defined one char string  !!  $!\n\n\n".$caller_inform ); 	}
	
  
  if (   (  defined ( $inArray )  ) && (  ref ( $inArray ) eq 'ARRAY' )   ){			}
	else{ DieWork::Just_dieWork( $die_MsgHead."\n \$inArray=$inArray should be a ARRAY ref !!  $!\n\n\n".$caller_inform );   }
	
  my $outSting;
    
  my $bigest_charNb=-999999999999999999999;
  my $inArrayChgedHash;
  foreach my $ec_CHAR_Pos ( @{$inArray} ){ 
  	$inArrayChgedHash->{$ec_CHAR_Pos}=1;  
  	if ( $bigest_charNb <= $ec_CHAR_Pos ){
  		$bigest_charNb = $ec_CHAR_Pos;
  	}
    #print "\$inArrayChgedHash->{$ec_CHAR_Pos}=$inArrayChgedHash->{$ec_CHAR_Pos}\n"; 
  }
  
 
 
  for (my $i=0; $i<$bigest_charNb; $i++){
    if (   (  defined ( $inArrayChgedHash->{$i+1} )  ) && ( $inArrayChgedHash->{$i+1} == 1)   ){
      $outSting.=$inWord;
    }
    else{
      $outSting.=' ';
    }
  }
 
  return $outSting;
}

sub findingIdexOfAword{   
  my ($inWord, $inString, $capsSensitive)=@_;
  
  if (  ( defined ($capsSensitive) ) && ( $capsSensitive == 1 )  ){ 	                                              }  
  else                                                            { $inWord=uc ($inWord); $inString=uc ($inString); }  
  
  my $outArray;;
  
  my $startIDX=0;
  while ($startIDX >=0){
    $startIDX=index ($inString, $inWord, $startIDX);  #print "$inWord, $inString\t\$startIDX=$startIDX\n";
    if ($startIDX>=0){   push (  @{ $outArray }, ( $startIDX+1 )  );   $startIDX++; }
  }
  
  return $outArray;

}    


sub TranslateDNAtoPep{  #翻译三联密码子
my ($inDNAstring)=@_;  #print "\$inDNAstring=$inDNAstring\n";
$inDNAstring=~s/\s//g; $inDNAstring=~s/\n//g;

my $myCodonTable   = Bio::Tools::CodonTable->new();
my $outSeq=$myCodonTable ->translate($inDNAstring);

#print "\$outSeq=$outSeq\n\n";
return $outSeq;                                           
}



#对特定序列进行翻译，翻译的时候，特定位置的密码子和翻译成的氨基酸，按输入的 $inTGAposhash中的信息来处理
sub TranslateSeqWithTGAcodonHash{
  my ($inDNAString, $inTGAposhash)=@_;  #$inTGAposhash->{33}->{'AminoAcidResidue'}={'U'};
 
  warn     "\n\nNow in Package FastaFileHandle,\nIn Sub &TranslateSeqWithTGAcodonHash!\n";      warn  "Input 1:\$inDNAString=$inDNAString\t" if (defined ($inDNAString)); warn  "Input 2:\$inTGAposhash=$inTGAposhash\n" if (defined ($inTGAposhash));
  print "\n\cl\nNow in Package FastaFileHandle,\nIn Sub &TranslateSeqWithTGAcodonHash!\n:";     print "Input 1:\$inDNAString=$inDNAString\t" if (defined ($inDNAString)); print "Input 2:\$inTGAposhash=$inTGAposhash\n" if (defined ($inTGAposhash));
  warn  "\n\n";
  print "\n\cl\n";
  DirFileHandle::PrintAndWarnDumper($inTGAposhash);
  
  $inDNAString=~s/\s//g; $inDNAString=~s/\n//g;
  my $InDNAlength=length ($inDNAString);             print "\$InDNAlength=$InDNAlength\n";  warn "\$InDNAlength=$InDNAlength\n";
  my $outPEPstring;
  if ( ($InDNAlength % 3) == 0 ){
  
    my $startIDX=0;    
    while ($startIDX <= $InDNAlength){
      my $DNA_3_codon=substr($inDNAString, $startIDX, 3);  #print "\$DNA_3_codon=$DNA_3_codon\n";  warn "\$DNA_3_codon=$DNA_3_codon\n";
      if (  defined ( $inTGAposhash->{$startIDX+1} )  ){
        $outPEPstring.=$inTGAposhash->{$startIDX+1}->{'AminoAcidResidue'};
      }                 
      else{
        my $AminoResidu=&TranslateDNAtoPep($DNA_3_codon);
        $outPEPstring.=$AminoResidu;                       #print "\$AminoResidu=$AminoResidu\n";   warn "\$AminoResidu=$AminoResidu\n";
      }                                                    
      
      $startIDX+=3; 
    
    }
  
  }
  else{ die "\n\n\nIn package FastaFileHandle,\nIn sub TranslateSeqWithTGAcodonHash,\n\n\\$inDNAString=\n$inDNAString\n\n$InDNAlength=$InDNAlength is not a integer multiple of 3!! \n\n\n";}
  print "\$outPEPstring=$outPEPstring\n";   warn "\$outPEPstring=$outPEPstring\n";
  return $outPEPstring;
}

#   TGA position (a,b,c) direction to tgc  通常的输入是 {0}->{ '0_directi'=>'+', '1_postion'=>'236877', '2_fromCha'=> 't', '3_intoCha'=> 'c'; }

sub CHangeSpecificSequenceInAFastaSeq{  #FastaFileHandle::CHangeSpecificSequenceInAFastaSeq #将fasta格式的文件中的特定位点的序列转换为 其它序列，通常用来 将TGA 转化为TGC，也就是转化为半胱氨酸，目的是用来 进行基因预测，结构识别。  转化完后，记录下相应位点，再转化回去。
  my (
          $InFastaSeq,                              #1
          $changeHash,                              #2
          $segSeqStart,                             #3
          $segSeqEndPo,                             #4
          $orgPos_in_biger_1__in_Small_2            #5
  )=@_;
  
  #warn     "\n\nNow in Package FastaFileHandle,\nIn Sub &CHangeSpecificSequenceInAFastaSeq!\n";      warn  "Input 1:\$InFastaSeq=\$InFastaSeq\t" if (defined ($InFastaSeq)); warn  "Input 2:\$changeHash=$changeHash\t" if (defined ($changeHash));
  #print "\n\cl\nNow in Package FastaFileHandle,\nIn Sub &CHangeSpecificSequenceInAFastaSeq!\n:";     print "Input 1:\$InFastaSeq=\$InFastaSeq\t" if (defined ($InFastaSeq)); print "Input 2:\$changeHash=$changeHash\t" if (defined ($changeHash));
  #warn  "Input 3:\$segSeqStart=$segSeqStart\t" if (defined ($segSeqStart));  warn  "Input 4:\$segSeqEndPo=$segSeqEndPo\t" if (defined ($segSeqEndPo));  
  #print "Input 3:\$segSeqStart=$segSeqStart\t" if (defined ($segSeqStart));  print "Input 4:\$segSeqEndPo=$segSeqEndPo\t" if (defined ($segSeqEndPo));   
  #warn  "Input 5:\$orgPos_in_biger_1__in_Small_2=$orgPos_in_biger_1__in_Small_2\t" if (defined ($orgPos_in_biger_1__in_Small_2));
  #print "Input 5:\$orgPos_in_biger_1__in_Small_2=$orgPos_in_biger_1__in_Small_2\t" if (defined ($orgPos_in_biger_1__in_Small_2));
  #warn  "\n\n";
  #print "\n\cl\n";
  
  DirFileHandle::PrintAndWarnDumper($changeHash) if (   (  defined ( $changeHash )  ) && (  ref ( $changeHash ) eq 'HASH'  )    );
  
  $InFastaSeq=~s/\s//g;  $InFastaSeq=~s/\n//g;
  my $OnlySequenc=$InFastaSeq;
  my $sequelength=length($InFastaSeq);
  
  my $situation;
  if (  ( defined ($segSeqStart) ) && ( $segSeqStart=~m/^\d+$/ )  ){
    if (   ( defined ($orgPos_in_biger_1__in_Small_2) ) && (  ( $orgPos_in_biger_1__in_Small_2 == 1 ) || ( $orgPos_in_biger_1__in_Small_2 == 2 )  )   ){  #当定义了 替换序列的位置信息后，必须明确指出 原来的需要修改的碱基，如TGA的a在大的序列段上（标注 1），还是小的序列段上（标注 2）
      if ( $orgPos_in_biger_1__in_Small_2 == 1 ){
        $situation=1; 
      }
      else{
        $situation=2; 
      }
    }
    else{
      die "\n\n\nIn package FastaFileHandle,\nIn sub CHangeSpecificSequenceInAFastaSeq,\n\n\$segSeqStart=$segSeqStart\n\$orgPos_in_biger_1__in_Small_2=$orgPos_in_biger_1__in_Small_2 must be set as 1 or 2!!!!\n\n\n";
    }
  }
  else {
    if (  ( defined ($segSeqEndPo) ) && ( $segSeqEndPo=~m/^\d+$/ )  ) {die "\n\n\nIn package FastaFileHandle,\nIn sub CHangeSpecificSequenceInAFastaSeq,\n\n\$segSeqEndPo=$segSeqEndPo\nThe segment start postion is not set, the End position should be not set too!!!!\n\n\n";}
    else {
      $situation=3;   		
    }
  }
  
  
  if ( ref($changeHash) eq 'HASH' ){
    foreach my $numBerHere (   sort (  keys %{ $changeHash }  )   ){
      my $directi=$changeHash->{$numBerHere}->{'0_directi'};
      my $postion=$changeHash->{$numBerHere}->{'1_postion'};
      my $fromCha=$changeHash->{$numBerHere}->{'2_fromCha'};
      my $intoCha=$changeHash->{$numBerHere}->{'3_intoCha'};
      
      if ($directi eq '-'){
        $fromCha=~tr/acgtACGT/tgcaTGCA/;
        $intoCha=~tr/acgtACGT/tgcaTGCA/;
      }
      
      
      if   ($situation==1){ $postion=$postion-$segSeqStart+1; }
      elsif($situation==2){ $postion=$postion+$segSeqStart-1; }
      elsif($situation==3){ $postion=$postion;                }
      
      if (  ($postion <= 0) || ($postion > $sequelength)  ){
        warn  "\n\n\nIn package FastaFileHandle,\nIn sub CHangeSpecificSequenceInAFastaSeq,\n\n\$postion=$postion <=0, or \$sequelength=$sequelength < $postion=\$postion \n\n\n";
        print "\n\n\nIn package FastaFileHandle,\nIn sub CHangeSpecificSequenceInAFastaSeq,\n\n\$postion=$postion <=0, or \$sequelength=$sequelength < $postion=\$postion \n\n\n";
      }
      else{
        my $targetInSeq=substr ( $OnlySequenc, ($postion-1), 1 ); 
        $targetInSeq=uc $targetInSeq;
        $fromCha=uc $fromCha;
        
        if ($fromCha ne $targetInSeq) {
          die "\n\n\nIn package FastaFileHandle,\nIn sub CHangeSpecificSequenceInAFastaSeq,\n\n\$fromCha=$fromCha ne $targetInSeq=\$targetInSeq\n\n\n";
        }
        else {
          $intoCha=uc $intoCha;                             print "20181120-1 $postion $fromCha $intoCha $OnlySequenc \n";
          substr ($OnlySequenc, ($postion-1), 1)=$intoCha;  print "20181120-2 $postion $fromCha $intoCha $OnlySequenc \n";
        }  		  
      }
    
    
    }
  }
  
  return $OnlySequenc;
  
}


#  my $outHash=FastaFileHandle::BuildHashFromFastaFile_Name_as_Key_checkRepeat ( $inseqFile ) ;    
sub BuildHashFromFastaFile_Name_as_Key_checkRepeat{    #建立以fasta id为key，序列长度为 值的 hash, 并检查是否有重复的 fasta的id出现, 记录了（1）具体地址，（2）原文件大小（3）fasta原id（4）middle的id（5）短id  （6）序列长度 (7) 其它信息
  my ($inseqFile)=@_;
  
  my $warnMsgBody="\nIn package  FastaFileHandle,\tIn sub BuildHashFromFastaFile_Name_as_Key_checkRepeat,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $subCallereIfm=DirFileHandle::print_SubCallerInform;
	
	if (   (  defined ( $inseqFile )  ) && ( $inseqFile=~m/\S+/ ) && (  -e ( $inseqFile )  )   ) {
	}
	else {
		DieWork::Just_dieWork( $die_MsgHead."\n\$inseqFile=$inseqFile is not exist!!! \n\n$! \n$subCallereIfm\n\n" );
	}	
	
	my @args = stat ($inseqFile);  my $size = $args[7];
	
	
  my $seqIOobj_IN =Bio::SeqIO->new(-file   => $inseqFile,    
                                   -format => 'fasta'      );
  my $outHash;
  
  
  while (my $seqObj=$seqIOobj_IN->next_seq){   #warn "\$seqObj->primary_id=", $seqObj->primary_id, "\n\$onceCountHash->{ $seqObj->primary_id }=" , $onceCountHash->{ $seqObj->primary_id } , "\n\$seqNb_2=$seqNb_2\n";
    my $idHere   =$seqObj->primary_id; 
    my $shortID  =FastaFileHandle::GetShortGenomeID   ( $idHere ) ;
    my $middleID =FastaFileHandle::GetMilddleGenomeID ( $idHere ) ;
    #my $seqHere  =$seqObj->seq;
    my $Seqleth  =$seqObj->length;            
    
    if (   (  defined ( $outHash )  ) && (  ref( $outHash ) eq 'HASH' ) && (  defined ( $outHash->{ $idHere } )  ) && (  $outHash->{ $idHere }=~/^\d+$/ )   ){
    	DieWork::Just_dieWork( $die_MsgHead."\n in \$inseqFile=$inseqFile \t\$outHash->{ $idHere }=$outHash->{ $idHere } is repeat !! : $! \n$subCallereIfm\n\n" );
    }
    else{
      $outHash->{ $idHere }->{'0_0_0_orgFile_path'}=$inseqFile;  
      $outHash->{ $idHere }->{'0_0_1_orgFile_size'}=$size;  
      $outHash->{ $idHere }->{'0_0_2_orgParimayID'}=$idHere;   
      $outHash->{ $idHere }->{'0_0_3_orgMiddle_ID'}=$middleID;   
      $outHash->{ $idHere }->{'0_0_4_orgShort__ID'}=$shortID;   
      $outHash->{ $idHere }->{'0_0_5_seq___length'}=$Seqleth;     
      #$outHash->{ $idHere }->{'0_0_6_____sequence'}=$seqHere;     	
    }
    
    
    
  }
  return $outHash;

}

#  my $outHash=FastaFileHandle::BuildHashFromFastaFile_Name_as_Key_checkRepeat_WithSeq ( $inseqFile ) ;    
sub BuildHashFromFastaFile_Name_as_Key_checkRepeat_WithSeq{    #建立以fasta id为key，序列长度为 值的 hash, 并检查是否有重复的 fasta的id出现, 记录了（1）具体地址，（2）原文件大小（3）fasta原id（4）middle的id（5）短id  （6）序列长度 (7) 其它信息
  my ($inseqFile)=@_;
  
  my $warnMsgBody="\nIn package  FastaFileHandle,\tIn sub BuildHashFromFastaFile_Name_as_Key_checkRepeat_WithSeq,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $subCallereIfm=DirFileHandle::print_SubCallerInform;
	
	if (   (  defined ( $inseqFile )  ) && ( $inseqFile=~m/\S+/ ) && (  -e ( $inseqFile )  )   ) {
	}
	else {
		DieWork::Just_dieWork( $die_MsgHead."\n\$inseqFile=$inseqFile is not exist!!! \n\n$! \n$subCallereIfm\n\n" );
	}	
	
	my @args = stat ($inseqFile);  my $size = $args[7];
	
	
  my $seqIOobj_IN =Bio::SeqIO->new(-file   => $inseqFile,    
                                   -format => 'fasta'      );
  my $outHash;
  
  
  while (my $seqObj=$seqIOobj_IN->next_seq){   #warn "\$seqObj->primary_id=", $seqObj->primary_id, "\n\$onceCountHash->{ $seqObj->primary_id }=" , $onceCountHash->{ $seqObj->primary_id } , "\n\$seqNb_2=$seqNb_2\n";
    my $idHere   =$seqObj->primary_id; 
    my $shortID  =FastaFileHandle::GetShortGenomeID   ( $idHere ) ;
    my $middleID =FastaFileHandle::GetMilddleGenomeID ( $idHere ) ;
    my $seqHere  =$seqObj->seq;
    my $Seqleth  =$seqObj->length;            
    
    if (   (  defined ( $outHash )  ) && (  ref( $outHash ) eq 'HASH' ) && (  defined ( $outHash->{ $idHere } )  ) && (  $outHash->{ $idHere }=~/^\d+$/ )   ){
    	DieWork::Just_dieWork( $die_MsgHead."\n in \$inseqFile=$inseqFile \t\$outHash->{ $idHere }=$outHash->{ $idHere } is repeat !! : $! \n$subCallereIfm\n\n" );
    }
    else{
      $outHash->{ $idHere }->{'0_0_0_orgFile_path'}=$inseqFile;  
      $outHash->{ $idHere }->{'0_0_1_orgFile_size'}=$size;  
      $outHash->{ $idHere }->{'0_0_2_orgParimayID'}=$idHere;   
      $outHash->{ $idHere }->{'0_0_3_orgMiddle_ID'}=$middleID;   
      $outHash->{ $idHere }->{'0_0_4_orgShort__ID'}=$shortID;   
      $outHash->{ $idHere }->{'0_0_5_seq___length'}=$Seqleth;     
      $outHash->{ $idHere }->{'0_0_6_____sequence'}=$seqHere;     	
    }
    
    
    
  }
  return $outHash;

}


#  my $outHash=FastaFileHandle::BuildFastaHash_Name_Key_seq_Val ( $inseqFile ) ;    
sub BuildFastaHash_Name_Key_seq_Val{    #建立以fasta id为key，序列 为 值的 hash
  my ($inseqFile)=@_;
  
  my $warnMsgBody="\nIn package  FastaFileHandle,\tIn sub BuildHashFromFastaFile_Name_as_Key,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $subCallereIfm=DirFileHandle::print_SubCallerInform;
	
	if (   (  defined ( $inseqFile )  ) && ( $inseqFile=~m/\S+/ ) && (  -e ( $inseqFile )  )   ) {
	}
	else {
		DieWork::Just_dieWork( $die_MsgHead."\n\$inseqFile=$inseqFile is not exist!!! \n\n$! \n$subCallereIfm\n\n" );
	}	
	
  my $seqIOobj_IN =Bio::SeqIO->new(-file   => $inseqFile,    
                                   -format => 'fasta'      );
  my $outHash;
  
  
  while (my $seqObj=$seqIOobj_IN->next_seq){   #warn "\$seqObj->primary_id=", $seqObj->primary_id, "\n\$onceCountHash->{ $seqObj->primary_id }=" , $onceCountHash->{ $seqObj->primary_id } , "\n\$seqNb_2=$seqNb_2\n";
    my $idHere =$seqObj->primary_id; 
    my $seqHere=$seqObj->seq;
    #my $Seqleth=$seqObj->length;            
    $seqHere=~s/\s+//;
    
    $outHash->{ $idHere }=$seqHere;    
    
  }
  return $outHash;

}
  
#  my $outHash=FastaFileHandle::BuildHashFromFastaFile_Name_as_Key ( $inseqFile ) ;    
sub BuildHashFromFastaFile_Name_as_Key{    #建立以fasta id为key，序列长度为 值的 hash
  my ($inseqFile)=@_;
  
  my $warnMsgBody="\nIn package  FastaFileHandle,\tIn sub BuildHashFromFastaFile_Name_as_Key,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $subCallereIfm=DirFileHandle::print_SubCallerInform;
	
	if (   (  defined ( $inseqFile )  ) && ( $inseqFile=~m/\S+/ ) && (  -e ( $inseqFile )  )   ) {
	}
	else {
		DieWork::Just_dieWork( $die_MsgHead."\n\$inseqFile=$inseqFile is not exist!!! \n\n$! \n$subCallereIfm\n\n" );
	}	
	
  my $seqIOobj_IN =Bio::SeqIO->new(-file   => $inseqFile,    
                                   -format => 'fasta'      );
  my $outHash;
  
  
  while (my $seqObj=$seqIOobj_IN->next_seq){   #warn "\$seqObj->primary_id=", $seqObj->primary_id, "\n\$onceCountHash->{ $seqObj->primary_id }=" , $onceCountHash->{ $seqObj->primary_id } , "\n\$seqNb_2=$seqNb_2\n";
    my $idHere =$seqObj->primary_id; 
    my $seqHere=$seqObj->seq;
    #my $Seqleth=$seqObj->length;            
    
    $outHash->{ $idHere }=1;    
    
  }
  return $outHash;

}


sub Build_short_middle_GenomeID_HASH{ #  my ( $outHash_short_to_long, $outHash_long_to_short )=@{ FastaFileHandle::Build_short_middle_GenomeID_HASH ( $inseqFile ) };
	my ($inseqFile)=@_;
  
  my $warnMsgBody="\nIn package  FastaFileHandle,\tIn sub Build_short_middle_GenomeID_HASH,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $subCallereIfm=DirFileHandle::print_SubCallerInform;
	
	if (   (  defined ( $inseqFile )  ) && ( $inseqFile=~m/\S+/ ) && (  -e ( $inseqFile )  )   ) {
	}
	else {
		DieWork::Just_dieWork( $die_MsgHead."\n\$inseqFile=$inseqFile is not exist!!! \n\n$! \n$subCallereIfm\n\n" );
	}	
	my $segNameHASH=FastaFileHandle::BuildHashFromFastaFile_Name_as_Key ( $inseqFile ) ;
	
	my $outHash_short_to_long;
	my $outHash_long_to_short;
	if (   (  defined ( $segNameHASH )  ) && (   ref ( $segNameHASH ) eq 'HASH' )   ){
		foreach my $eachNm (   keys (  %{ $segNameHASH }  )   ){
			my $midleName=FastaFileHandle::GetMilddleGenomeID ( $eachNm ) ;
			my $shortName=FastaFileHandle::GetShortGenomeID   ( $eachNm ) ;
			
			$outHash_short_to_long->{$shortName}=$midleName;
			$outHash_long_to_short->{$midleName}=$shortName;
			
		}
	}
	
	return [$outHash_short_to_long, $outHash_long_to_short];
	
}

sub Build_short_long_GenomeID_HASH{ #  my ( $outHash_short_to_long, $outHash_long_to_short )=@{ FastaFileHandle::Build_short_long_GenomeID_HASH ( $inseqFile ) };
	my ($inseqFile)=@_;
  
  my $warnMsgBody="\nIn package  FastaFileHandle,\tIn sub Build_short_long_GenomeID_HASH,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $subCallereIfm=DirFileHandle::print_SubCallerInform;
	
	if (   (  defined ( $inseqFile )  ) && ( $inseqFile=~m/\S+/ ) && (  -e ( $inseqFile )  )   ) {
	}
	else {
		DieWork::Just_dieWork( $die_MsgHead."\n\$inseqFile=$inseqFile is not exist!!! \n\n$! \n$subCallereIfm\n\n" );
	}	
	my $segNameHASH=FastaFileHandle::BuildHashFromFastaFile_Name_as_Key ( $inseqFile ) ;
	
	my $outHash_short_to_long;
	my $outHash_long_to_short;
	if (   (  defined ( $segNameHASH )  ) && (   ref ( $segNameHASH ) eq 'HASH' )   ){
		foreach my $eachNm (   keys (  %{ $segNameHASH }  )   ){
			my $shortName=FastaFileHandle::GetShortGenomeID ( $eachNm ) ;
			
			$outHash_short_to_long->{$shortName}=$eachNm;
			$outHash_long_to_short->{$eachNm   }=$shortName;
		}
	}
	
	return [$outHash_short_to_long, $outHash_long_to_short];
	
}


sub GetShortGenomeID{  #  my $outString=FastaFileHandle::GetShortGenomeID ( $inString ) ;
	my ($inString)=@_;
	
	my $warnMsgBody="\nIn package  FastaFileHandle,\tIn sub GetShortGenomeID,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $subCallereIfm=DirFileHandle::print_SubCallerInform;
	
	my $outString;
	if (   (  defined ( $inString )  ) && ( $inString=~m/\S+/ )   ) {  #gi|158269951|gb|DS496432.1|  gi|684179921|ref|NC_024811.1| gi|507107719|emb|HG002314.1|   gi|151559145|dbj|AP006502.2|
		#if    ( $inString=~m/^gi\|\d+\|gb\|(\w).\d+\|/ ){  $outString=$1;			warn "20190109-1 \$inString=$inString\t\$outString=$outString\n";	}	
		if       ( $inString=~m/^gi\|\d+\|gb\|(\w+)\.\d+/  ){  $outString=$1;	}		#warn "20190109-0-1 \$inString=$inString\t\$outString=$outString\n";	}	
		elsif    ( $inString=~m/^gi\|\d+\|ref\|(\w+)\.\d+/ ){  $outString=$1;	}		#warn "20190109-0-2 \$inString=$inString\t\$outString=$outString\n";	}	
		elsif    ( $inString=~m/^gi\|\d+\|emb\|(\w+)\.\d+/ ){  $outString=$1;	}		#warn "20190109-0-3 \$inString=$inString\t\$outString=$outString\n";	}	
		elsif    ( $inString=~m/^gi\|\d+\|dbj\|(\w+)\.\d+/ ){  $outString=$1;	}		#warn "20190109-0-4 \$inString=$inString\t\$outString=$outString\n";	}	
		elsif    ( $inString=~m/^(\w+)\s*$/                ){  $outString=$1;	}    #warn "20190109-0-a \$inString=$inString\t\$outString=$outString\n";	}	
		elsif    ( $inString=~m/^lcl\|(\d+\w+\d+)/         ){  $outString=$1;	}    #warn "20190109-0-b \$inString=$inString\t\$outString=$outString\n";	}	  #lcl|1sig0
		elsif    ( $inString=~m/^lcl\|(\d+\w+\d+)/         ){  $outString=$1;	}    #warn "20190109-0-b \$inString=$inString\t\$outString=$outString\n";	}	  #lcl|1sig0NZ_CP012871.1
		elsif    ( $inString=~m/^lcl\|(\w+_\d+)/           ){  $outString=$1;	}    #warn "20190109-0-b \$inString=$inString\t\$outString=$outString\n";	}	  #lcl|JLDATANB_1.1
		elsif    ( $inString=~m/^(\w+)\.\d+/               ){  $outString=$1; }    #warn "20190109-0-b \$inString=$inString\t\$outString=$outString\n"; }   #NZ_CP012871.1  
		
		
		
		#elsif ( $inString=~m/\S+/ ){  $outString=$1;		  		}	
		else{		DieWork::Just_dieWork( $die_MsgHead."\n\$inString=$inString is not right!!! \n\n$! \n$subCallereIfm\n\n" ); }
		
		 print "\$inString=$inString\t\$outString=$outString\n"; 
	}
	else{
		DieWork::Just_dieWork( $die_MsgHead."\n\$inString=$inString is empty!!! \n\n$! \n$subCallereIfm\n\n" );
	}
	return $outString;
}

sub GetMilddleGenomeID{  #  my $outString=FastaFileHandle::GetMilddleGenomeID ( $inString ) ;
	my ($inString)=@_;
	
	my $warnMsgBody="\nIn package  FastaFileHandle,\tIn sub GetMilddleGenomeID,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $subCallereIfm=DirFileHandle::print_SubCallerInform;
	#warn "20190109-1-0 \$inString=$inString\n";	
	my $outString;
	if (   (  defined ( $inString )  ) && ( $inString=~m/\S+/ )   ) {  #gi|158269951|gb|DS496432.1|  gi|684179921|ref|NC_024811.1| gi|507107719|emb|HG002314.1|   gi|151559145|dbj|AP006502.2|
		#if    ( $inString=~m/^gi\|\d+\|gb\|(\w).\d+\|/ ){  $outString=$1;			warn "20190109-1 \$inString=$inString\t\$outString=$outString\n";	}	
		if       ( $inString=~m/^gi\|\d+\|(gb\|\w+\.\d+\|)/  ){  $outString=$1;		}	#warn "20190109-1-1 \$inString=$inString\t\$outString=$outString\n";	}	
		elsif    ( $inString=~m/^gi\|\d+\|(ref\|\w+\.\d+\|)/ ){  $outString=$1;		}	#warn "20190109-1-2 \$inString=$inString\t\$outString=$outString\n";	}	
		elsif    ( $inString=~m/^gi\|\d+\|(emb\|\w+\.\d+\|)/ ){  $outString=$1;		}	#warn "20190109-1-3 \$inString=$inString\t\$outString=$outString\n";	}	
		elsif    ( $inString=~m/^gi\|\d+\|(dbj\|\w+\.\d+\|)/ ){  $outString=$1;		}	#warn "20190109-1-4 \$inString=$inString\t\$outString=$outString\n";	}	
		elsif    ( $inString=~m/^(\w+)\s*$/                  ){  $outString=$1;   }  #warn "20190109-1-a \$inString=$inString\t\$outString=$outString\n";	}	
		elsif    ( $inString=~m/^(lcl\|\d+\w+\d+)/           ){  $outString=$1;	  }  #warn "20190109-0-b \$inString=$inString\t\$outString=$outString\n";	}	  #lcl|1sig0
		elsif    ( $inString=~m/^(lcl\|\w+_\d+)/             ){  $outString=$1;	}    #warn "20190109-0-b \$inString=$inString\t\$outString=$outString\n";	}	  #lcl|JLDATANB_1.1
		elsif    ( $inString=~m/^(\w+\.\d+)/                 ){  $outString=$1; }    #warn "20190109-0-b \$inString=$inString\t\$outString=$outString\n"; }   #NZ_CP012871.1
		
		#elsif ( $inString=~m/\S+/ ){  $outString=$1;		  		}	
		else{		DieWork::Just_dieWork( $die_MsgHead."\n\$inString=$inString is not right!!! \n\n$! \n$subCallereIfm\n\n" ); }
		
		#warn "20190109-1-b \$inString=$inString\t\$outString=$outString\n"; 
	}
	else{
		DieWork::Just_dieWork( $die_MsgHead."\n\$inString=$inString is empty!!! \n\n$! \n$subCallereIfm\n\n" );
	}
	return $outString;
}


sub GetHead{    #pmAble#   #获得输入文本的的第一个非空无空格字符串
  my ($inString)=@_;
  $inString=~s/^\s*//; $inString=~s/\s*$//; $inString=~s/^(\S+).*$/$1/; #$inString=substr($inString,0,10);
  return $inString;
}

#            GetSeqLength{       #获取fasta格式序列的  序列全长
sub GetSeqLength{   #FastaFileHandle::GetSeqLength       #pmAble#   #获取fasta格式序列的  序列全长
  my ($fastaSeq)=@_;
  if ($fastaSeq=~m/^\s*>\s*\S+.*\n[^>]+$/){}  else{    die "in sub GetSeqLength, the input val is not in fasta format:\n \$fastaSeq=$fastaSeq\n ";  }  #待检查
  $fastaSeq=~s/^\s+//g;
  $fastaSeq=~s/\s+$//g;
  $fastaSeq=~s/>.*//;
  $fastaSeq=~s/\s+//g;
  my $SeqLength=length $fastaSeq;
  return $SeqLength;

}




#sub10.5             GetSeqFastaHeadName{   #获取fasta格式序列的  序列名
sub GetSeqFastaHeadName{   #获取fasta格式序列的  序列名
  my ($fastaSeq)=@_;
  my $fastaHeadName;
  if ($fastaSeq=~m/^\s*>(\s*\S+.*)\n[^>]+$/){$fastaHeadName=$1;}  else{    die "in sub GetSeqFastaHeadName, the input val is not in fasta format:\n \$fastaSeq=$fastaSeq\n ";  }  #待检查 
  return $fastaHeadName;  
}

#sub10.6             GetSequence{  #获取fasta格式序列的  序列
sub GetSequence{  #获取fasta格式序列的  序列
my ($fastaSeq)=@_;
  if ($fastaSeq=~m/^\s*>\s*\S+.*\n[^>]+$/){}  else{    die "in sub GetSequence, the input val is not in fasta format:\n \$fastaSeq=$fastaSeq\n ";  }  #待检查
  $fastaSeq=~s/^\s+//g;
  $fastaSeq=~s/\s+$//g;
  $fastaSeq=~s/>.*//;
  $fastaSeq=~s/\s+//g;
  return $fastaSeq;  
} 

sub BuildHashFromFastaFile_sequence_as_Key_CheckSameToDIE{  #my $outHash=FastaFileHandle::BuildHashFromFastaFile_sequence_as_Key_CheckSameToDIE($inseqFile, $inFormat);
my ($inseqFile, $inFormat)=@_;
  my $seqIOobj_IN =Bio::SeqIO->new(-file   => $inseqFile,    
                                   -format => $inFormat       );
  my $outHash;
  
  while (my $seqObj=$seqIOobj_IN->next_seq){   #warn "\$seqObj->primary_id=", $seqObj->primary_id, "\n\$onceCountHash->{ $seqObj->primary_id }=" , $onceCountHash->{ $seqObj->primary_id } , "\n\$seqNb_2=$seqNb_2\n";
    if (   defined (  $outHash->{ $seqObj->seq }  )   ){
      my $idHere=$seqObj->primary_id; my $seqHere=$seqObj->seq; my $idTP=$outHash->{ $seqObj->seq };
      die "\n\n\nIn package FastaFileHandle, in sub BuildHashFromFastaFile_sequence_as_Key\n\$inseqFile=$inseqFile\t\$inFormat=$inFormat ,\n the sequences of these two id is the same:\n$idTP and $idHere\n The sequence:\n$seqHere\n\n";
    }
    else {
      $outHash->{ $seqObj->seq } = $seqObj->primary_id;
    }
  }
  return $outHash;
  
}


sub BuildHashFromFastaFile_sequence_as_Key{
  my ($inseqFile, $inFormat)=@_;
  my $seqIOobj_IN =Bio::SeqIO->new(-file   => $inseqFile,    
                                   -format => $inFormat       );
  my $outHash;
  
  while (my $seqObj=$seqIOobj_IN->next_seq){   #warn "\$seqObj->primary_id=", $seqObj->primary_id, "\n\$onceCountHash->{ $seqObj->primary_id }=" , $onceCountHash->{ $seqObj->primary_id } , "\n\$seqNb_2=$seqNb_2\n";
    my $idHere=$seqObj->primary_id; 
    my $seqHere=$seqObj->seq;
    
    if (   defined (  $outHash->{ $seqHere }  )   ){
      push @{  $outHash->{ $seqHere }  }, $idHere; 
      warn "\n\nJust Warning.\nIn package FastaFileHandle, in sub BuildHashFromFastaFile_sequence_as_Key\n\$inseqFile=$inseqFile\t\$inFormat=$inFormat ,\n the sequences  id: $idHere is the same as others!!!\nThe sequence:\n$seqHere\n\n";
    }
    else {
      $outHash->{ $seqHere }->[0] = $idHere;
    }
  }
  return $outHash;
  
}

sub BuildHashFromFastaFile_seqID_as_Key{  #  #  my $outHash=FastaFileHandle::BuildHashFromFastaFile_seqID_as_Key( $inseqFile ); 
	my ($inseqFile)=@_;
	
	my $warnMsgBody="\nIn package  FastaFileHandle,\tIn sub BuildHashFromFastaFile_seqID_as_Key,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  #my $caller_inform=DirFileHandle::print_SubCallerInform;
  
                                                                                              #warn "\$inseqFile=$inseqFile "; sleep (5);
   
	my $outHash;
	if  (   (  defined ( $inseqFile )  ) && (  -e ( $inseqFile )  )   ){
	
	  my $seqInObj=Bio::SeqIO->new(-file   => $inseqFile,    
                                 -format => 'fasta'     );
      
                        
    while (my $seqObj=$seqInObj->next_seq){    
      $outHash->{ $seqObj->primary_id }=$seqObj->seq();                                      #warn "\$inseqFile=$inseqFile \$seqObj->primary_id=".$seqObj->primary_id;  
    }
	
	}
	
	
  return $outHash;
	
}


#使用命令如：   samtools faidx Emiliania_huxleyiCCMP1516.fasta                                           [10:35AM]
#会提示错误：   [fai_build_core] different line length in sequence 'gi|485647780|gb|KB863137.1|'.
#                zsh: segmentation fault  samtools faidx Emiliania_huxleyiCCMP1516.fasta
#下面这个函数是用来修正这个错误的


sub fixSamltool{
  my ($inputfilename, $outputfilename)=@_;
  use Bio::SeqIO;
  my $in  = Bio::SeqIO->new(-file => "$inputfilename",
                            -format => 'Fasta');
  my $out = Bio::SeqIO->new(-file => ">$outputfilename",
                            -format => 'Fasta');
  while ( my $seq = $in->next_seq() ) {$out->write_seq($seq); }

}


sub BuildAAarray{  # my $aaArray=FastaFileHandle::BuildAAarray();
	my $aaArray;
	
	$aaArray->[0 ]='A';      	 #Alanine       
  $aaArray->[1 ]='C';      	 #Cysteine      
  $aaArray->[2 ]='D';      	 #Aspartic Acid 
  $aaArray->[3 ]='E';      	 #Glutamic Acid 
  $aaArray->[4 ]='F';      	 #Phenylalanine 
  $aaArray->[5 ]='G';      	 #Glycine       
  $aaArray->[6 ]='H';      	 #Histidine     
  $aaArray->[7 ]='I';      	 #Isoleucine    
  $aaArray->[8 ]='K';      	 #Lysine        
  $aaArray->[9 ]='L';      	 #Leucine       
  $aaArray->[10]='M';      	 #Methionine    
  $aaArray->[11]='N';      	 #Asparagine    
  $aaArray->[12]='P';      	 #Proline       
  $aaArray->[13]='Q';      	 #Glutamine     
  $aaArray->[14]='R';      	 #Arginine      
  $aaArray->[15]='S';      	 #Serine        
  $aaArray->[16]='T';      	 #Threonine     
  $aaArray->[17]='U';	       #Selenocysteine
  $aaArray->[18]='V';      	 #Valine        
  $aaArray->[19]='W';      	 #Tryptophan    
  $aaArray->[20]='Y';      	 #Tyrosine      

  return $aaArray;
}

sub GetCodonTable1{
my $codonTable_1;
$codonTable_1->{'TCA'}='S';      # Serine
$codonTable_1->{'TCC'}='S';      # Serine
$codonTable_1->{'TCG'}='S';      # Serine
$codonTable_1->{'TCT'}='S';      # Serine
$codonTable_1->{'TTC'}='F';      # Phenylalanine
$codonTable_1->{'TTT'}='F';      # Phenylalanine
$codonTable_1->{'TTA'}='L';      # Leucine
$codonTable_1->{'TTG'}='L';      # Leucine
$codonTable_1->{'TAC'}='Y';      # Tyrosine
$codonTable_1->{'TAT'}='Y';      # Tyrosine
$codonTable_1->{'TAA'}="*";      # Stop
$codonTable_1->{'TAG'}="*";      # Stop
$codonTable_1->{'TGC'}='C';      # Cysteine
$codonTable_1->{'TGT'}='C';      # Cysteine
$codonTable_1->{'TGA'}='*';      # Stop
$codonTable_1->{'TGG'}='W';      # Tryptophan
$codonTable_1->{'CTA'}='L';      # Leucine
$codonTable_1->{'CTC'}='L';      # Leucine
$codonTable_1->{'CTG'}='L';      # Leucine
$codonTable_1->{'CTT'}='L';      # Leucine
$codonTable_1->{'CCA'}='P';      # Proline
$codonTable_1->{'CCC'}='P';      # Proline
$codonTable_1->{'CCG'}='P';      # Proline
$codonTable_1->{'CCT'}='P';      # Proline
$codonTable_1->{'CAC'}='H';      # Histidine
$codonTable_1->{'CAT'}='H';      # Histidine
$codonTable_1->{'CAA'}='Q';      # Glutamine
$codonTable_1->{'CAG'}='Q';      # Glutamine
$codonTable_1->{'CGA'}='R';      # Arginine
$codonTable_1->{'CGC'}='R';      # Arginine
$codonTable_1->{'CGG'}='R';      # Arginine
$codonTable_1->{'CGT'}='R';      # Arginine
$codonTable_1->{'ATA'}='I';      # Isoleucine
$codonTable_1->{'ATC'}='I';      # Isoleucine
$codonTable_1->{'ATT'}='I';      # Isoleucine
$codonTable_1->{'ATG'}='M';      # Methionine
$codonTable_1->{'ACA'}='T';      # Threonine
$codonTable_1->{'ACC'}='T';      # Threonine
$codonTable_1->{'ACG'}='T';      # Threonine
$codonTable_1->{'ACT'}='T';      # Threonine
$codonTable_1->{'AAC'}='N';      # Asparagine
$codonTable_1->{'AAT'}='N';      # Asparagine
$codonTable_1->{'AAA'}='K';      # Lysine
$codonTable_1->{'AAG'}='K';      # Lysine
$codonTable_1->{'AGC'}='S';      # Serine
$codonTable_1->{'AGT'}='S';      # Serine
$codonTable_1->{'AGA'}='R';      # Arginine
$codonTable_1->{'AGG'}='R';      # Arginine
$codonTable_1->{'GTA'}='V';      # Valine
$codonTable_1->{'GTC'}='V';      # Valine
$codonTable_1->{'GTG'}='V';      # Valine
$codonTable_1->{'GTT'}='V';      # Valine
$codonTable_1->{'GCA'}='A';      # Alanine
$codonTable_1->{'GCC'}='A';      # Alanine
$codonTable_1->{'GCG'}='A';      # Alanine
$codonTable_1->{'GCT'}='A';      # Alanine
$codonTable_1->{'GAC'}='D';      # Aspartic Acid
$codonTable_1->{'GAT'}='D';      # Aspartic Acid
$codonTable_1->{'GAA'}='E';      # Glutamic Acid
$codonTable_1->{'GAG'}='E';      # Glutamic Acid
$codonTable_1->{'GGA'}='G';      # Glycine
$codonTable_1->{'GGC'}='G';      # Glycine
$codonTable_1->{'GGG'}='G';      # Glycine
$codonTable_1->{'GGT'}='G';      # Glycine
return $codonTable_1;
}

sub CodonTable_standard {     #   my $output=FastaFileHandle::CodonTable_standard($codon);
  my($codon) = @_;
  if (  ( defined ($codon) ) && ($codon=~m/\S+/)  ){
    $codon=uc $codon ;
    
    my $CodonTable=&GetCodonTable1;
    #$CodonTable->{'TGA'}='U';      # Selenocystein
    
    
    if ( defined ($CodonTable->{$codon}) ) {          	return $CodonTable->{$codon};                                      }
    else                                   {           	return 'X' ;            print STDERR "Bad codon \"$codon\"!!\n";   }
  }
}



sub CodonTable_TGA_U {     #   FastaFileHandle::CodonTable_TGA_U
  my($codon) = @_;
  if (  ( defined ($codon) ) && ($codon=~m/\S+/)  ){
    $codon=uc $codon ;
    
    my $CodonTable=&GetCodonTable1;
    $CodonTable->{'TGA'}='U';      # Selenocystein
    
    
    if ( defined ($CodonTable->{$codon}) ) {          	return $CodonTable->{$codon};                                      }
    else                                   {           	return 'X' ;            print STDERR "Bad codon \"$codon\"!!\n";   }
  }
}

sub CodonTable_TGA_J {     #   FastaFileHandle::CodonTable_TGA_J
  my($codon) = @_;
  if (  ( defined ($codon) ) && ($codon=~m/\S+/)  ){
    $codon=uc $codon ;
    
    my $CodonTable=&GetCodonTable1;
    $CodonTable->{'TGA'}='J';      # Selenocystein
    
    
    if ( defined ($CodonTable->{$codon}) ) {          	return $CodonTable->{$codon};                                      }
    else                                   {           	return 'X' ;            print STDERR "Bad codon \"$codon\"!!\n";   }
  }
}


sub CodonTable_TGA_U_TAG_O_TAA_B {   #   FastaFileHandle::CodonTable_TGA_U_TAG_O_TAA_B
  my($codon) = @_; $codon=uc $codon;
  if (  ( defined ($codon) ) && ($codon=~m/\S+/)  ){                  
    my $CodonTable=&GetCodonTable1;
    $CodonTable->{'TGA'}='U';      # Selenocystein
    $CodonTable->{'TAG'}='O';      # Selenocystein
    $CodonTable->{'TAA'}='B';      # Selenocystein
    
    
    if ( defined ($CodonTable->{$codon}) ) {          	return $CodonTable->{$codon};                                      }
    else                                   {           	return 'X' ;            print STDERR "Bad codon \"$codon\"!!\n";   }
  }
}

sub CodonTable_TAG_O_TAA_B {   #   FastaFileHandle::CodonTable_TAG_O_TAA_B
  my($codon) = @_; $codon=uc $codon;
  if (  ( defined ($codon) ) && ($codon=~m/\S+/)  ){                  
    my $CodonTable=&GetCodonTable1;
    #$CodonTable->{'TGA'}='U';      # Selenocystein
    $CodonTable->{'TAG'}='O';      # Selenocystein
    $CodonTable->{'TAA'}='B';      # Selenocystein
    
    
    if ( defined ($CodonTable->{$codon}) ) {          	return $CodonTable->{$codon};                                      }
    else                                   {           	return 'X' ;            print STDERR "Bad codon \"$codon\"!!\n";   }
  }
}


sub BuildFastaFile_with_pure_string{           #FastaFileHandle::BuildFastaFile_with_pure_string ($FilePath, $inString, $fastaName);
  my ($FilePath, $inString, $fastaName)=@_;
  my $warnMsgHead="\n\n\n   In package  FastaFileHandle,\tIn sub BuildFastaFile_with_pure_string,\n\n";
  if (defined ($fastaName)){} else {$fastaName='FASTA'; warn "$warnMsgHead\nNo FASTA name set, use FASTA as name\n\n"; print  "$warnMsgHead\nNo FASTA name set, use FASTA as name\n\n";}
  open (IN,">$FilePath")  or die "DIE!!!!!\n$warnMsgHead cannot create \$FilePath=$FilePath:$!\n\n";
  $fastaName=~s/^\s+//;  $fastaName=~s/\s+$//;  $inString=~s/\s+//g;
  print IN ">$fastaName\n$inString\n\n";
  close (IN);
  return 1;
}

sub BuildFastaFile_with_pure_string_forInterProscan{           #FastaFileHandle::BuildFastaFile_with_pure_string_forInterProscan
  my ($FilePath, $inString, $fastaName)=@_;
  my $warnMsgHead="\n\n\n   In package  FastaFileHandle,\tIn sub BuildFastaFile_with_pure_string_forInterProscan,\n\n";
  if (defined ($fastaName)){} else {$fastaName='FASTA'; warn "$warnMsgHead\nNo FASTA name set, use FASTA as name\n\n"; print  "$warnMsgHead\nNo FASTA name set, use FASTA as name\n\n";}
  open (IN,">$FilePath")  or die "DIE!!!!!\n$warnMsgHead cannot create \$FilePath=$FilePath:$!\n\n";
  $fastaName=~s/^\s+//;  $fastaName=~s/\s+$//;  $inString=~s/\s+//g; $inString=~s/\*//g;
  print IN ">$fastaName\n$inString\n\n";
  close (IN);
  return 1;
}

sub BuildFastaFile_withFastaString{    #FastaFileHandle::BuildFastaFile_withFastaString ($FilePath, $inString);
  
  my ($FilePath, $inString)=@_;
  my $warnMsgHead="\n\n\n   In package  FastaFileHandle,\tIn sub BuildFastaFile_withFastaString,\n\n";
  open (IN,">$FilePath")  or die "DIE!!!!!\n$warnMsgHead cannot create \$FilePath=$FilePath:$!\n\n";
  
  print IN "$inString";
  close (IN);
  return 1;
}

sub ReverseComplementString{  #my $ReverseComplement_String=FastaFileHandle::ReverseComplementString ($inString);
  my ($inString)=@_;
  
  my $warnMsgBody="\nIn package  FastaFileHandle,\tIn sub ReverseComplementString,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $caller_inform=DirFileHandle::print_SubCallerInform;
  
  $inString=~s/\s//g;
  if (   ( defined ($inString) ) && ($inString=~m/^[acgtACGTrymkswhbvdnRYMKSWHBVDN]+$/)   ){
    my $ReverseComplement_String=reverse $inString;

    #my $JianBingHash;	
    #$JianBingHash->{'R'}=['A','G'];
    #$JianBingHash->{'Y'}=['C','T'];
    #$JianBingHash->{'M'}=['A','C'];
    #$JianBingHash->{'K'}=['G','T'];
    #$JianBingHash->{'S'}=['G','C'];
    #$JianBingHash->{'W'}=['A','T'];
    #$JianBingHash->{'H'}=['A','T','C'];
    #$JianBingHash->{'B'}=['G','T','C'];
    #$JianBingHash->{'V'}=['G','A','C'];
    #$JianBingHash->{'D'}=['G','A','T'];
    #$JianBingHash->{'N'}=['A','T','C','G'];

    $ReverseComplement_String=~tr/acgtACGTrymkswhbvdnRYMKSWHBVDN/tgcaTGCAyrkmswdvbhnYRKMSWDVBHN/;

    return $ReverseComplement_String;	
  }
  else {
  	#print "WARN!!!! \$inTCGAbase=$inTCGAbase is not a ACTG type char which could be reverse and complemented!!\n\n\n";
  	return $inString;
  	DieWork::Just_dieWork( $die_MsgHead."Only acgtACGTrymkswhbvdnRYMKSWHBVDN were allowed, something wrong found in \$inString=\n$inString\n\n\n" );
  }
}

sub ReverseComplementBase{
  my ($inTCGAbase)=@_;
  if (   ( defined ($inTCGAbase) ) && ($inTCGAbase=~m/[acgtACGTrymkswhbvdnRYMKSWHBVDN]/)   ){
    my $ReverseComplementBASE=reverse $inTCGAbase;

    #my $JianBingHash;	
    #$JianBingHash->{'R'}=['A','G'];
    #$JianBingHash->{'Y'}=['C','T'];
    #$JianBingHash->{'M'}=['A','C'];
    #$JianBingHash->{'K'}=['G','T'];
    #$JianBingHash->{'S'}=['G','C'];
    #$JianBingHash->{'W'}=['A','T'];
    #$JianBingHash->{'H'}=['A','T','C'];
    #$JianBingHash->{'B'}=['G','T','C'];
    #$JianBingHash->{'V'}=['G','A','C'];
    #$JianBingHash->{'D'}=['G','A','T'];
    #$JianBingHash->{'N'}=['A','T','C','G'];

    $ReverseComplementBASE=~tr/acgtACGTrymkswhbvdnRYMKSWHBVDN/tgcaTGCAyrkmswdvbhnYRKMSWDVBHN/;

    return $ReverseComplementBASE;	
  }
  else {
  	#print "WARN!!!! \$inTCGAbase=$inTCGAbase is not a ACTG type char which could be reverse and complemented!!\n\n\n";
  	return $inTCGAbase;
  }
}

1;