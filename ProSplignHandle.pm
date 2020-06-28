
#!/usr/bin/perl -w

BEGIN{  push (@INC, '/home/fredjiang/fredProgs/PMs')  };

use strict;
use warnings;

use Storable ;
use Storable qw(nstore store_fd nstore_fd freeze thaw dclone);



use XML::Simple;

use DirFileHandle;
use FastaFileHandle;
use BlastHandle;
use SeqSegmentsTools;
use MatrixCsvChange;

package  ProSplignHandle;



sub proSplign_check_selProtein_vs_DNA{ #Input 1, protein sequence string, 2, DNA sequence string, 3, a directory path to hold middle step information. TARGET: map the protein sequence back to dna sequence, find the coding region and in frame tga positions 
	my ( $inPEPstring,                         #1
	     $inDNAstring,                         #2
	     $analysisDIR,                         #3
	     $in_gnoSeg_stt,                       #4
	     $in_gnoSeg_end,                       #5
	     $in_gnoSeg_ZoF,                       #6
	     $oriGinalDNAString,                   #7
	     $in_TgaPos,                           #8
	     $GenoOrEst_type,                      #9
	     $wkInform                             #10
	     )=@_;
	my $dieHeadMsg="\n\n\n		In package FastaFileHandle,		\n In sub proSplign_check_selProtein_vs_DNA,\n";  	#warn  "$dieHeadMsg\n\$inPEPstring=$inPEPstring\n\n"; #\n\$inDNAstring=$inDNAstring 	#print "$dieHeadMsg\nIn sub Map_PepSeq_BackTo_DNAseq,\n\$inPEPstring=$inPEPstring\n\n"; #\n\$inDNAstring=$inDNAstring
	
	$inDNAstring=~s/\n//g; $inDNAstring=~s/\s//g;
	my $DNAseqLength=length($inDNAstring);	
	
	my $tempTimeDirNeedDel=0; if (  ( defined($analysisDIR) ) && ($analysisDIR=~m/\S+/)  ) {} else {		$analysisDIR=TimeWork::GetTimeDirOrFileName();	$tempTimeDirNeedDel=1; }  #If the $analysisDIR is not set, then build a temporary directory named by time #print "\$analysisDIR=$analysisDIR\n";
  system ( "mkdir -p $analysisDIR ") ;
	
	
	my $PEPseqFastaFile="$analysisDIR/0_0_0_0_PEPseqFastaFile.txt"; FastaFileHandle::BuildFastaFile_with_pure_string($PEPseqFastaFile, $inPEPstring, 'PEP');            
	my $DNAseqFastaFile="$analysisDIR/0_0_0_1_DNAseqFastaFile.txt";	FastaFileHandle::BuildFastaFile_with_pure_string($DNAseqFastaFile, $inDNAstring, 'DNA'); 
	
	my $db_gencode=1; 	my $U_pos_array=FastaFileHandle::Find_U_pos_inAminoAcid($inPEPstring);	if ( ref($U_pos_array) eq 'ARRAY' ){ $db_gencode=10; }
	
	my $XMLreportArray=&RunProSplign($PEPseqFastaFile, $DNAseqFastaFile, $analysisDIR, $db_gencode, $GenoOrEst_type);
	
	my $porSplignHash=&praseProsplignXmlOutArray( $XMLreportArray,  #1 
	                                              $in_gnoSeg_stt,   #2 
	                                              $in_gnoSeg_end,   #3 
	                                              $in_gnoSeg_ZoF,   #4 
	                                              '',               #5 
	                                              '',               #6 
	                                              '',               #7 
	                                              $wkInform         #8 
	                                             );
		
	#my $porSplignHashName='0_0_2_0_ProSplPrsXmlHsH.txt'; my $porSplignHashFile="$analysisDIR/$porSplignHashName"; 	DirFileHandle::PrintDumper($porSplignHashFile,$porSplignHash)  if ( ref($porSplignHash) eq 'HASH' );
	
	
		
		
		
	my $proSplign_with_AADNAposHASH=&Build_DNA_PEP_correlate_HASH($porSplignHash, $inPEPstring, $oriGinalDNAString, $in_TgaPos);
	#my $DNA_PEP_coRe_HASH='0_0_2_1_DNA_PEP_coReHsH.txt';  if ( ref ($proSplign_with_AADNAposHASH) eq 'HASH' ) { DirFileHandle::PrintDumper ("$analysisDIR/$DNA_PEP_coRe_HASH", $proSplign_with_AADNAposHASH); }
	#my $simpDNAPEPco_HASH='0_0_2_2_SempDNAPEPcoHsH.txt'; my $simp_proSplign_with_AADNAposHASH=&simplify_prosplign_Out_Hash($proSplign_with_AADNAposHASH); DirFileHandle::PrintDumper ("$analysisDIR/$simpDNAPEPco_HASH", $simp_proSplign_with_AADNAposHASH) if ( ref ($simp_proSplign_with_AADNAposHASH) eq 'HASH' );
  return $proSplign_with_AADNAposHASH;
}




sub InsertSmallSegmentIntoBigSegment_noCrossBigSegment{ # insert small segment into big segment, such as "exon-------tga--------exon" ==> "exon_head----exon_middle TGA--TGA Exon_middle Exon-tail"
  my ($bigSegArray, $org_smlSegArray)=@_;  
  my $bigArryDump=DirFileHandle::ReturnDumperInform($bigSegArray); my $org_smlArryDump=DirFileHandle::ReturnDumperInform($org_smlSegArray);
  my $ptOutMsg="\n\n201808061858\$bigArryDump=$bigArryDump\n\n201808061858\n\$org_smlArryDump=$org_smlArryDump\n\n"; print $ptOutMsg;
  my $Warn____Msg="		In package ProSplignHandle,		\n In sub InsertSmallSegmentIntoBigSegment_noCrossBigSegment,\n";  #InsertSmallSegmentIntoBigSegment_inOneBigSegment
  my $WarnHeadMsg="\n\n\n$Warn____Msg";   my $die_HeadMsg="\n\n\n\nDIE!!!!!\n$Warn____Msg";
  
  my ($bigSegHeadKey, $bigSegTailKey, $bigSegTypeKey,  $bigSegZorFKey)
  =  ('0_SgmtHead',   '1_SgmtTail',   '3_SgmtType',    '2_SgmtZorF'  );
  my ($smlSegHeadKey, $smlSegTailKey, $smlSegTypeKey,  $smlSegZorFKey)
  =  ('0_SgmtHead',   '1_SgmtTail',   '3_SgmtType',    '2_SgmtZorF'  );
  my ($newSegHeadKey, $newSegTailKey, $newSegTypeKey,  $newSegZorFKey, $newSegTypeKey_bkgrd   )
  =  ('0_SgmtHead',   '1_SgmtTail',   '3_SgmtType',    '2_SgmtZorF',   '4_SgmtType_backgroud');
  
  
  #step 0 Make the divided Upstream_downstream_TGA into nomal small array
  #$VAR1 = [
  #          {
  #            '0upStream_exon_part' => {
  #                                       '0_SgmtHead' => '7129',
  #                                       '1_SgmtTail' => 7130,
  #                                       '2_SgmtZorF' => '+',
  #                                       '3_SgmtType' => 'Sec_TGA'
  #                                     },
  #            '1DwStream_exon_part' => {
  #                                       '0_SgmtHead' => '7249',
  #                                       '1_SgmtTail' => '7249',
  #                                       '2_SgmtZorF' => '+',
  #                                       '3_SgmtType' => 'Sec_TGA'
  #                                     }
  #          }
  #        ];
  my $smlSegArray; my $newIdx=0;
  if ( ref ($org_smlSegArray) eq 'ARRAY'){
    for (  my $i=0; $i<@{ $org_smlSegArray }; $i++  ){
    	if ( ref ($org_smlSegArray->[$i]) eq 'HASH'){
    		if (        (   (  defined ( $org_smlSegArray->[$i]->{'0upStream_exon_part'} )  ) && (  ref( $org_smlSegArray->[$i]->{'0upStream_exon_part'} ) eq 'HASH'  )   )
    		        ||  (   (  defined ( $org_smlSegArray->[$i]->{'0upStream_exon_part'} )  ) && (  ref( $org_smlSegArray->[$i]->{'0upStream_exon_part'} ) eq 'HASH'  )   )
    		   )
    		{
    		  if (   (  defined ( $org_smlSegArray->[$i]->{'0upStream_exon_part'} )  ) && (  ref( $org_smlSegArray->[$i]->{'0upStream_exon_part'} ) eq 'HASH'  )   ){
    		  	$smlSegArray->[$newIdx]=$org_smlSegArray->[$i]->{'0upStream_exon_part'};
    		  	$newIdx++;
    		  }
    		  if (   (  defined ( $org_smlSegArray->[$i]->{'1DwStream_exon_part'} )  ) && (  ref( $org_smlSegArray->[$i]->{'1DwStream_exon_part'} ) eq 'HASH'  )   ){
    		  	$smlSegArray->[$newIdx]=$org_smlSegArray->[$i]->{'1DwStream_exon_part'};
    		  	$newIdx++;
    		  }
    		}
    		elsif(   (  defined ( $org_smlSegArray->[$i]->{'0_SgmtHead'} )  ) && ( $org_smlSegArray->[$i]->{'0_SgmtHead'}=~m/\d+/ )   ) {
    		  $smlSegArray->[$newIdx]=$org_smlSegArray->[$i];
    		  $newIdx++;
    	  }
    	  else {
    	  	my $new_smlArryDump=DirFileHandle::ReturnDumperInform($smlSegArray);
      	  my $dieMsg="$die_HeadMsg\nthe Small Segment is not right!!\n\$org_smlArryDump=$org_smlArryDump\n\$new_smlArryDump=$new_smlArryDump\n\n";
          print $dieMsg; die $dieMsg;
    	  }
    	}
    	
    }
  }

  
  
  #step 1 check each of small segment. To make sure that, all of the small segments are included in one of big segments
  if ( ref ($smlSegArray) eq 'ARRAY'){
    for (  my $j=0; $j<@{ $smlSegArray }; $j++  ){
        	
    	
      my $smlInBig_or_not=0;
      if ( ref ($bigSegArray) eq 'ARRAY'){
        FORMK: for (  my $i=0; $i<@{ $bigSegArray }; $i++  ){
          
          if ( $bigSegArray->[$i]->{$bigSegZorFKey} eq $smlSegArray->[$j]->{$smlSegZorFKey} ){
        	  my $ZorF_one=1; if ($bigSegArray->[$i]->{$bigSegZorFKey} eq '-'){ $ZorF_one=-1; }
        	  
            if (
                 ( SeqSegmentsTools::a_between_b_c( $smlSegArray->[$j]->{$smlSegHeadKey}, $bigSegArray->[$i]->{$bigSegHeadKey}, $bigSegArray->[$i]->{$bigSegTailKey} )  )
                 &&
                 ( SeqSegmentsTools::a_between_b_c( $smlSegArray->[$j]->{$smlSegTailKey}, $bigSegArray->[$i]->{$bigSegHeadKey}, $bigSegArray->[$i]->{$bigSegTailKey} )  )
               ){
              $smlInBig_or_not=1;
              last FORMK;
              
            }
          }
        }
      }
      if ( $smlInBig_or_not==0 ){
      	my $bigArryDump=DirFileHandle::ReturnDumperInform($bigSegArray); my $smlArryDump=DirFileHandle::ReturnDumperInform($smlSegArray);
      	my $dieMsg="$die_HeadMsg\nOne of the Small Segment is not included by big Segment!!\n\$smlSegArray->[$j]->{$smlSegHeadKey}=$smlSegArray->[$j]->{$smlSegHeadKey}\n\$smlSegArray->[$j]->{$smlSegTailKey}=$smlSegArray->[$j]->{$smlSegTailKey}\n\$smlSegArray->[$j]->{$smlSegTypeKey}=$smlSegArray->[$j]->{$smlSegTypeKey}\n\$smlSegArray->[$j]->{$smlSegZorFKey}=$smlSegArray->[$j]->{$smlSegZorFKey}\n\$bigArryDump=$bigArryDump\n\$smlArryDump=$smlArryDump\n\n";
        print $dieMsg; die $dieMsg;
      }
    }
  } 
  
  #step 2 build new segment array, in which the big segments were divided by small segments
  my $newArray;
  if ( ref ($bigSegArray) eq 'ARRAY'){
    for (  my $i=0; $i<@{ $bigSegArray }; $i++  ){
    	
      my $smlInBig_or_not=0;
      if ( ref ($smlSegArray) eq 'ARRAY'){
        for (  my $j=0; $j<@{ $smlSegArray }; $j++  ){
          
          if ( $bigSegArray->[$i]->{$bigSegZorFKey} eq $smlSegArray->[$j]->{$smlSegZorFKey} ){
        	  my $ZorF_one=1; if ($bigSegArray->[$i]->{$bigSegZorFKey} eq '-'){ $ZorF_one=-1; }
        	  
            if (
                 ( SeqSegmentsTools::a_between_b_c( $smlSegArray->[$j]->{$smlSegHeadKey}, $bigSegArray->[$i]->{$bigSegHeadKey}, $bigSegArray->[$i]->{$bigSegTailKey} )  )
                 &&
                 ( SeqSegmentsTools::a_between_b_c( $smlSegArray->[$j]->{$smlSegTailKey}, $bigSegArray->[$i]->{$bigSegHeadKey}, $bigSegArray->[$i]->{$bigSegTailKey} )  )
               ){
               
              $smlInBig_or_not=1;
              
              if ($smlSegArray->[$j]->{$smlSegHeadKey}-$bigSegArray->[$i]->{$bigSegHeadKey} != 0){
                my $tempHash_1;              
                $tempHash_1->{$newSegHeadKey}=$bigSegArray->[$i]->{$bigSegHeadKey};
                $tempHash_1->{$newSegTailKey}=$smlSegArray->[$j]->{$smlSegHeadKey}-$ZorF_one;
                $tempHash_1->{$newSegTypeKey}=$bigSegArray->[$i]->{$bigSegTypeKey}; #befor_TGA?                
                $tempHash_1->{$newSegZorFKey}=$bigSegArray->[$i]->{$bigSegZorFKey};
                push @{ $newArray }, $tempHash_1;
              }
              if (1){
                my $tempHash_2;              
                $tempHash_2->{$newSegHeadKey}=$smlSegArray->[$j]->{$smlSegHeadKey};
                $tempHash_2->{$newSegTailKey}=$smlSegArray->[$j]->{$smlSegTailKey};
                $tempHash_2->{$newSegTypeKey}=$smlSegArray->[$j]->{$smlSegTypeKey};   $tempHash_2->{$newSegTypeKey_bkgrd}=$bigSegArray->[$i]->{$bigSegTypeKey};            
                $tempHash_2->{$newSegZorFKey}=$smlSegArray->[$j]->{$smlSegZorFKey};
                push @{ $newArray }, $tempHash_2;
              }
              if ($smlSegArray->[$j]->{$smlSegTailKey}-$bigSegArray->[$i]->{$bigSegTailKey} != 0){
                my $tempHash_3;              
                $tempHash_3->{$newSegHeadKey}=$smlSegArray->[$j]->{$smlSegTailKey}+$ZorF_one;
                $tempHash_3->{$newSegTailKey}=$bigSegArray->[$i]->{$bigSegTailKey};
                $tempHash_3->{$newSegTypeKey}=$bigSegArray->[$i]->{$bigSegTypeKey}; #after_TGA?
                $tempHash_3->{$newSegZorFKey}=$bigSegArray->[$i]->{$bigSegZorFKey};
                push @{ $newArray }, $tempHash_3;
              }
            }
          }
        
             
        }
      }
       
      if ($smlInBig_or_not==0){
      	my $tempHash;              
        $tempHash->{$newSegHeadKey}=$bigSegArray->[$i]->{$bigSegHeadKey};
        $tempHash->{$newSegTailKey}=$bigSegArray->[$i]->{$bigSegTailKey};
        $tempHash->{$newSegTypeKey}=$bigSegArray->[$i]->{$bigSegTypeKey};
        $tempHash->{$newSegZorFKey}=$bigSegArray->[$i]->{$bigSegZorFKey};
        push @{ $newArray }, $tempHash;
      }
    
    }
  }
  my $newArrayDump=DirFileHandle::ReturnDumperInform($newArray); #my $smlArryDump=DirFileHandle::ReturnDumperInform($newArray);
  my $ptOutMsg2="\n\n201808061858\$newArrayDump=$newArrayDump\n"; print $ptOutMsg2;
  return $newArray;
}

sub CheckAll_U_uc_Match{
	my ($inArray)=@_;
	my $warnMsgHead="\n\n\n   In package  ProSplignHandle,\tIn sub CheckAll_U_uc_Match,\n\n";
	
	my @matchRIghtAminoAcid=('U','C');
	
	#print "CheckAll_U_uc_Match DirFileHandle::PrintAndWarnDumper (\$inArray)\n";
	#DirFileHandle::PrintAndWarnDumper ($inArray);
	
	my $all_match=0;
	if ( ref ($inArray) eq 'ARRAY' ){
		FMKOUT: for (my $i=0; $i<@{ $inArray }; $i++){
			
			my $matchAminoAcid=$inArray->[$i]->{'1_mcP_aminoA'}; $matchAminoAcid=uc $matchAminoAcid;
			
			my $matchRight=0;
			FMK: foreach my $AAhere (@matchRIghtAminoAcid){
				if ($matchAminoAcid eq $AAhere){ #print "\$matchAminoAcid eq $AAhere=$matchAminoAcid eq $AAhere\n";
					$matchRight=1;
					last FMK;
				}
			}
			if ($matchRight==0){	  		$all_match=0;  last FMKOUT;			}
			else               {				$all_match=1;			              } 
			print "\$all_match=$all_match";
		}
	}
	return $all_match;
}

sub Calculate_U_Match_pct{
	
	my ($inArray, $matchAminoAcidArray)=@_;
	my $warnMsgHead="\n\n\n   In package  ProSplignHandle,\tIn sub Calculate_U_Match_pct,\n\n";
	my $matchNumber=&Calculate_U_Match_number($inArray, $matchAminoAcidArray);
	my $pct=0;  
	if ( ref ($inArray) eq 'ARRAY' ){
		my $all_U_number=@{ $inArray };  
	  $pct=$matchNumber/$all_U_number;
  }
  return $pct;
}

sub Calculate_U_Match_number{
	my ($inArray, $matchAminoAcidArray)=@_;
	my $warnMsgHead="\n\n\n   In package  ProSplignHandle,\tIn sub Calculate_U_Match_number,\n\n";
	
	my @matchRIghtAminoAcid=@{ $matchAminoAcidArray };  # ('U','C');
	
	#print "Calculate_U_Match_number DirFileHandle::PrintAndWarnDumper (\$inArray)\n";
	#DirFileHandle::PrintAndWarnDumper ($inArray);
	#print "Calculate_U_Match_number DirFileHandle::PrintAndWarnDumper (\$matchAminoAcidArray)\n";
	#DirFileHandle::PrintAndWarnDumper ($matchAminoAcidArray);
	
	my $matchNumber=0;  
	if ( ref ($inArray) eq 'ARRAY' ){

		FMKOUT: for (my $i=0; $i<@{ $inArray }; $i++){
			
			my $matchAminoAcid=$inArray->[$i]->{'1_mcP_aminoA'}; $matchAminoAcid=uc $matchAminoAcid;
			
			my $matchRight=0;
			FMK: foreach my $AAhere (@matchRIghtAminoAcid){
				if ($matchAminoAcid eq $AAhere){ print "\$matchAminoAcid eq $AAhere=$matchAminoAcid eq $AAhere\n";
					$matchRight=1;
					last FMK;
				}
			}
			if ($matchRight==0){	  		                               		}
			else               {				 $matchNumber++;		              } 
			print "\$matchNumber=$matchNumber\n\n";
		}
		
		
	}
	
	return $matchNumber;
}


sub BuildPEPcoverageAnlysis{
	my ($ExonArray)=@_; 
	my $changeArray;
	if ( ref ($ExonArray) eq 'ARRAY' ){
		for (my $i=0; $i<@{ $ExonArray }; $i++){
			#my ( $pepStt, $pepEnd )=( $ExonArray->[$i]->{'1_PEP_0stt'}, $ExonArray->[$i]->{'1_PEP_1end'} );
			$changeArray->[$i]->{'0_Stt'  }=$ExonArray->[$i]->{'1_PEP_0stt_0Pos'};  print "\$changeArray->[$i]->{'0_Stt'}=\$ExonArray->[$i]->{'1_PEP_0stt_0Pos'}=$changeArray->[$i]->{'0_Stt'}\n";
			$changeArray->[$i]->{'1_End'  }=$ExonArray->[$i]->{'1_PEP_1end_0Pos'};  print "\$changeArray->[$i]->{'1_End'}=\$ExonArray->[$i]->{'1_PEP_1end_0Pos'}=$changeArray->[$i]->{'1_End'}\n";
			$changeArray->[$i]->{'2_ZoF'  }='+';   #beacuse this mathod is used for protein sequence which have no minus(-) strand
		}
	}
	my $outArray=ProSplignHandle::overLayerHadle_with_ZoF($changeArray);
	return $outArray;
}

sub overLayerHadle_with_ZoF{     # 处理overlayyer的hash
  my ($inArrayHash)=@_;
  my $warnMsgHead="\n\nIn Package ProSplignHandle,\t In sub overLayerHadle_with_ZoF,\n\n";
  my @sortByHeadArray;
  if    ( $inArrayHash->[0]->{'2_ZoF'} eq '-' ){  #从小到大
  	@sortByHeadArray=sort{ $a->{'0_Stt'} <=> $b->{'0_Stt'} } @{ $inArrayHash };
  }
  elsif ( $inArrayHash->[0]->{'2_ZoF'} eq '+' ){  #从大到小
  	@sortByHeadArray=sort{ $b->{'0_Stt'} <=> $a->{'0_Stt'} } @{ $inArrayHash };
  }
  else{
  	my $dieMsg= "\n\n\nDIE!!!!$warnMsgHead\nThe values below Should be + or -!!\n\$inArrayHash->[0]->{'2_ZoF'}=$inArrayHash->[0]->{'2_ZoF'}\n\n";print $dieMsg; die $dieMsg;
  }
  my $JointUnitsNBhash; my $jointIdx=0;  my $unitNumber=@sortByHeadArray;  my $lastMergeOrNotMergeHash;
  for (my $i=0; $i<@sortByHeadArray; $i++){  print "\$sortByHeadArray[$i]->{'0_Stt'}=$sortByHeadArray[$i]->{'0_Stt'}\t\$sortByHeadArray[$i]->{'1_End'}=$sortByHeadArray[$i]->{'1_End'}\n";
  	if ($unitNumber ==1){
  	  $JointUnitsNBhash->[0]=$sortByHeadArray[0];  print "00 \$JointUnitsNBhash->[$jointIdx]->{'0_Stt'}=$JointUnitsNBhash->[$jointIdx]->{'0_Stt'}\t\$JointUnitsNBhash->[$jointIdx]->{'1_End'}=$JointUnitsNBhash->[$jointIdx]->{'1_End'}\n";
  	}
  	else{
  		
  		if ( $i==0){   #第一项
  		  $lastMergeOrNotMergeHash=$sortByHeadArray[0];   print "11 \$lastMergeOrNotMergeHash->{'0_Stt'}=$lastMergeOrNotMergeHash->{'0_Stt'}\t\$lastMergeOrNotMergeHash->{'1_End'}=$lastMergeOrNotMergeHash->{'1_End'}\n";
  		}
  		
      if (  ($i>0) && ($i<($unitNumber-1) )   ){ #第二项到倒数第二项
       
        if ( SeqSegmentsTools::OverlayCheck($lastMergeOrNotMergeHash->{'0_Stt'}, $lastMergeOrNotMergeHash->{'1_End'}, $sortByHeadArray[$i]->{'0_Stt'}, $sortByHeadArray[$i]->{'1_End'}) ){
          my ($newHd, $newEd)=@{ SeqSegmentsTools::JointSegments_sameDirection_Up_or_down_with_ZorF_inform($lastMergeOrNotMergeHash->{'0_Stt'}, $lastMergeOrNotMergeHash->{'1_End'}, $lastMergeOrNotMergeHash->{'2_ZoF'}, $sortByHeadArray[$i]->{'0_Stt'}, $sortByHeadArray[$i]->{'1_End'}, $sortByHeadArray[$i]->{'2_ZoF'}) };
          $lastMergeOrNotMergeHash->{'0_Stt'}=$newHd; 
          $lastMergeOrNotMergeHash->{'1_End'}  =$newEd; print "22 \$lastMergeOrNotMergeHash->{'0_Stt'}=$lastMergeOrNotMergeHash->{'0_Stt'}\t\$lastMergeOrNotMergeHash->{'1_End'}=$lastMergeOrNotMergeHash->{'1_End'}\n";
        }
        else {
          $JointUnitsNBhash->[$jointIdx]=$lastMergeOrNotMergeHash; print "33 \$JointUnitsNBhash->[$jointIdx]->{'0_Stt'}=$JointUnitsNBhash->[$jointIdx]->{'0_Stt'}\t\$JointUnitsNBhash->[$jointIdx]->{'1_End'}=$JointUnitsNBhash->[$jointIdx]->{'1_End'}\n"; $jointIdx++; 
          $lastMergeOrNotMergeHash=$sortByHeadArray[$i]; print "33 \$lastMergeOrNotMergeHash->{'0_Stt'}=$lastMergeOrNotMergeHash->{'0_Stt'}\t\$lastMergeOrNotMergeHash->{'1_End'}=$lastMergeOrNotMergeHash->{'1_End'}\n";
        }
        
      }
      
      if ( $i == ($unitNumber-1) ){         #最后一项
      	
      	if ($i ==0){  #既是最后一项，又是第一项
      	  $JointUnitsNBhash->[$jointIdx]=$sortByHeadArray[0];      $jointIdx++;
      	}
      	else{      	  #是非第一项的 最后一项
      	  if ( SeqSegmentsTools::OverlayCheck($lastMergeOrNotMergeHash->{'0_Stt'}, $lastMergeOrNotMergeHash->{'1_End'}, $sortByHeadArray[$i]->{'0_Stt'}, $sortByHeadArray[$i]->{'1_End'}) ){
            my ($newHd, $newEd)=@{ SeqSegmentsTools::JointSegments_sameDirection_Up_or_down_with_ZorF_inform($lastMergeOrNotMergeHash->{'0_Stt'}, $lastMergeOrNotMergeHash->{'1_End'}, $lastMergeOrNotMergeHash->{'2_ZoF'}, $sortByHeadArray[$i]->{'0_Stt'}, $sortByHeadArray[$i]->{'1_End'}, $sortByHeadArray[$i]->{'2_ZoF'}) };
            $lastMergeOrNotMergeHash->{'0_Stt'}=$newHd; 
            $lastMergeOrNotMergeHash->{'1_End'}  =$newEd;   print "44 \$lastMergeOrNotMergeHash->{'0_Stt'}=$lastMergeOrNotMergeHash->{'0_Stt'}\t\$lastMergeOrNotMergeHash->{'1_End'}=$lastMergeOrNotMergeHash->{'1_End'}\n";
            
            $JointUnitsNBhash->[$jointIdx]=$lastMergeOrNotMergeHash;print "44 \$JointUnitsNBhash->[$jointIdx]->{'0_Stt'}=$JointUnitsNBhash->[$jointIdx]->{'0_Stt'}\t\$JointUnitsNBhash->[$jointIdx]->{'1_End'}=$JointUnitsNBhash->[$jointIdx]->{'1_End'}\n"; $jointIdx++;
          }
          else {
          	
          	$JointUnitsNBhash->[$jointIdx]=$lastMergeOrNotMergeHash;print "55 \$JointUnitsNBhash->[$jointIdx]->{'0_Stt'}=$JointUnitsNBhash->[$jointIdx]->{'0_Stt'}\t\$JointUnitsNBhash->[$jointIdx]->{'1_End'}=$JointUnitsNBhash->[$jointIdx]->{'1_End'}\n"; $jointIdx++;
            $JointUnitsNBhash->[$jointIdx]=$sortByHeadArray[$i];     print "55 \$JointUnitsNBhash->[$jointIdx]->{'0_Stt'}=$JointUnitsNBhash->[$jointIdx]->{'0_Stt'}\t\$JointUnitsNBhash->[$jointIdx]->{'1_End'}=$JointUnitsNBhash->[$jointIdx]->{'1_End'}\n"; $jointIdx++;
            
          }       
        }
      	
      }
    }
  } 
  return $JointUnitsNBhash;
}  

sub GetCuttedLength_proSplign_version{  # 获得 cutted的 片段的长度
  my ($inArrayHash)=@_;
  
  my $coverLength;   print "\n\n\n\n\n";
  for (my $i=0; $i<@{ $inArrayHash }; $i++){
  	my $segMLength=(  ( abs ($inArrayHash->[$i]->{'1_End'}-$inArrayHash->[$i]->{'0_Stt'}) ) + 1  );  print "\$inArrayHash->[$i]->{'1_End'}=$inArrayHash->[$i]->{'1_End'}\t\$inArrayHash->[$i]->{'0_Stt'}=$inArrayHash->[$i]->{'0_Stt'}\t\t\$segMLength=$segMLength\t\t";
  	$coverLength+=$segMLength;  print "\$coverLength=$coverLength\n";
  } print "\n\n\n\n\n";
  return $coverLength;
}

sub simplify_prosplign_Out_Hash{
	my ($inHash)=@_;
	my $outHash=Storable::dclone ($inHash);
	if (  ref( $outHash->{'0_MatchArray'} ) eq 'ARRAY'  ){
		for (my $i=0; $i<@{ $outHash->{'0_MatchArray'} }; $i++){
			delete ( $outHash->{'0_MatchArray'}->[$i]->{'7_DNA_to_PEP_hash'} );
			delete ( $outHash->{'0_MatchArray'}->[$i]->{'7_PEP_to_DNA_hash'} );
			delete ( $outHash->{'0_MatchArray'}->[$i]->{'5_PepCov_1ary'} );    
			                                                                   
		}                                                                    
	}                                                                      
	return $outHash;                                                       
}                                                                        
                                                                         
                                                                         
sub praseProsplignXmlOutArray{
	
	
	
	my ( $inArray,         #1 
	     $gnoSeg_stt,      #2 
	     $gnoSeg_end,      #3 
	     $gnoSeg_ZoF,      #4 
	     $proSeg_stt,      #5 
	     $proSeg_end,      #6 
	     $proSeg_ZoF,      #7 
	     $wkInform         #8 
	    )=@_;  print "\$inArray=$inArray\n";  
  
  my $warnMsgHead="\n\n\n   In package  ProSplignHandle,\nIn sub praseProsplignXmlOutHash,\n\n";
	if (    (  defined ( $wkInform )  ) && ( $wkInform=~m/\S+/ )   ){
		$warnMsgHead.=$wkInform."\n\n";
	}
  
  
  my $gno_ChgBack=0;
  if (   (  ( defined ($gnoSeg_stt) ) && ( $gnoSeg_stt=~m/^\d+$/ )  ) || (  ( defined ($gnoSeg_end) ) && ( $gnoSeg_end=~m/^\d+$/ )  ) || (  ( defined ($gnoSeg_ZoF) ) && ( $gnoSeg_ZoF=~m/^(\+|-)$/ )  )   ){
    if (     ( defined ($gnoSeg_stt) ) && ( defined ($gnoSeg_end) ) && ( defined ($gnoSeg_ZoF) ) 
          && ( $gnoSeg_stt=~m/^\d+$/ ) && ( $gnoSeg_end=~m/^\d+$/ ) && ( $gnoSeg_ZoF=~m/^(\+|-)$/ )
       ){
      $gno_ChgBack=1;
    }
    else {
    	my $dieMsg="\n\n\nDIE!!!!!\n$warnMsgHead\nThe genome Segment Information is not right!!!!\n\$gnoSeg_stt=$gnoSeg_stt, \$gnoSeg_end=$gnoSeg_end, \$gnoSeg_ZoF=$gnoSeg_ZoF\n\n";
    	print $dieMsg; die $dieMsg;
    }
  }
	
	my $pro_ChgBack=0;
  if (   (  ( defined ($proSeg_stt) ) && ( $proSeg_stt=~m/^\d+$/ )  ) || (  ( defined ($proSeg_end) ) && ( $proSeg_end=~m/^\d+$/ )  ) || (  ( defined ($proSeg_ZoF) ) && ( $proSeg_ZoF=~m/^(\+|-)$/ )  )   ){
    if (     ( defined ($proSeg_stt) ) && ( defined ($proSeg_end) ) && ( defined ($proSeg_ZoF) ) 
          && ( $proSeg_stt=~m/^\d+$/ ) && ( $proSeg_end=~m/^\d+$/ ) && ( $proSeg_ZoF=~m/^(\+|-)$/ )
       ){
      $pro_ChgBack=1;
    }
    else {
    	my $dieMsg="\n\n\nDIE!!!!!\n$warnMsgHead\nThe protein Segment Information is not right!!!!\n\$proSeg_stt=$proSeg_stt, \$proSeg_end=$proSeg_end, \$proSeg_ZoF=$proSeg_ZoF\n\n";
    	print $dieMsg; die $dieMsg;
    }
  }
  
	
	my $outHash;
	if ( ref($inArray) eq 'ARRAY'){
		for (my $i=0; $i<@{ $inArray }; $i++){
			
			
			if (  defined ( $inArray->[$i]->{'Seq-align_ext'}->{'User-object'}->{'User-object_data'}->{'User-field'}->{'User-field_data'}->{'User-field_data_int'}  )  ){
			  my $ProSplignGeneNM=$inArray->[$i]->{'Seq-align_ext'}->{'User-object'}->{'User-object_data'}->{'User-field'}->{'User-field_data'}->{'User-field_data_int'};
			  if ($ProSplignGeneNM=~m/^\d+$/){
			  	$outHash->{'0_MatchArray'}->[$i]->{'0_0_PrSpNM'}=$ProSplignGeneNM;
			  } else {my $dieMsg_1= "DIE!!!!$warnMsgHead\n\$ProSplignGeneNM=$ProSplignGeneNM\nshound be a number!!\n\n\n"; print $dieMsg_1; die $dieMsg_1;}
			} else {my $dieMsg_2= "DIE!!!!$warnMsgHead\n\$inArray->[$i]->{'Seq-align_ext'}->{'User-object'}->{'User-object_data'}->{'User-field'}->{'User-field_data'}->{'User-field_data_int'}=$inArray->[$i]->{'Seq-align_ext'}->{'User-object'}->{'User-object_data'}->{'User-field'}->{'User-field_data'}->{'User-field_data_int'}\nshound be defined!!\n\n\n"; print $dieMsg_2; die $dieMsg_2;}
			
			if (   (  defined ( $inArray->[$i]->{'Seq-align_score'}->{'Score'} )  ) && ( ref ( $inArray->[$i]->{'Seq-align_score'}->{'Score'} ) eq 'ARRAY'  )    ){
				my $scoreArrayRef=$inArray->[$i]->{'Seq-align_score'}->{'Score'};
				for (  my $ScIdx=0; $ScIdx<@{ $scoreArrayRef }; $ScIdx++  ){
					my $EachScoreHashRef=$scoreArrayRef->[$ScIdx];
					if (  defined ( $EachScoreHashRef->{'Score_id'}->{'Object-id'}->{'Object-id_str'} )  ){
						if  (  defined ( $EachScoreHashRef->{'Score_value'}->{'Score_value_int'} )  ) {
							my $Object_id_str=$EachScoreHashRef->{'Score_id'}->{'Object-id'}->{'Object-id_str'}; my $score=$EachScoreHashRef->{'Score_value'}->{'Score_value_int'};
							if    ($Object_id_str eq 'align_length'               ){     $outHash->{'0_MatchArray'}->[$i]->{'0_1_aln_0_1_aln_lengt'}=$score; } #0_1_aln_0_1_aln_lengt
							elsif ($Object_id_str eq 'num_ident'                  ){     $outHash->{'0_MatchArray'}->[$i]->{'0_1_aln_1_ident'}=$score; }
							elsif ($Object_id_str eq 'num_positives'              ){     $outHash->{'0_MatchArray'}->[$i]->{'0_1_aln_2_postv'}=$score; }							
							elsif ($Object_id_str eq 'num_negatives'              ){     $outHash->{'0_MatchArray'}->[$i]->{'0_1_aln_3_negtv'}=$score; }
							elsif ($Object_id_str eq 'genomic_gap_length'         ){     $outHash->{'0_MatchArray'}->[$i]->{'0_1_aln_4_gn_gp'}=$score; }
							elsif ($Object_id_str eq 'product_gap_length'         ){     $outHash->{'0_MatchArray'}->[$i]->{'0_1_aln_5_pd_gp'}=$score; }
							elsif ($Object_id_str eq 'product_internal_gap_length'){     $outHash->{'0_MatchArray'}->[$i]->{'0_1_aln_6_pdIgp'}=$score; }
							else  { print "warn!!!!$warnMsgHead\n\$Object_id_str=\$Object_id_str\nThis is a new Object_id_str!!\n\n\n"; }
							
						}
						else {
							DirFileHandle::PrintAndWarnDumper ($EachScoreHashRef);
							my $dieMsgMsg= "DIE!!!!$warnMsgHead\n\$EachScoreHashRef=\$scoreArrayRef->[$ScIdx]=$scoreArrayRef->[$ScIdx]\nThe key Score_value_int shound be Defined!!\n\n\n";
							print $dieMsgMsg; die $dieMsgMsg;
						}
					}
				}
			}
			
			
			
			
			if (   (  defined ( $inArray->[$i]->{'Seq-align_segs'}->{'Seq-align_segs_spliced'}->{'Spliced-seg'} )  ) && (  ref ( $inArray->[$i]->{'Seq-align_segs'}->{'Seq-align_segs_spliced'}->{'Spliced-seg'} ) eq 'HASH'  )   ){
			  my $spliced_seqHash=$inArray->[$i]->{'Seq-align_segs'}->{'Seq-align_segs_spliced'}->{'Spliced-seg'};
			  
			  
			  
			  $outHash->{'0_MatchArray'}->[$i]->{'1_0PepNm'}=$spliced_seqHash->{'Spliced-seg_product-id'}->{'Seq-id'}->{'Seq-id_local'}->{'Object-id'}->{'Object-id_str'};
			  $outHash->{'0_MatchArray'}->[$i]->{'1_1DnaNm'}=$spliced_seqHash->{'Spliced-seg_genomic-id'}->{'Seq-id'}->{'Seq-id_local'}->{'Object-id'}->{'Object-id_str'};
			  
			  my $DNA_ZoF=$spliced_seqHash->{'Spliced-seg_genomic-strand'}->{'Na-strand'}->{'value'}; 
			  if ( $DNA_ZoF eq 'minus' ) { $DNA_ZoF='-'; 	  				  } else { $DNA_ZoF='+';  }
			  if ( $gno_ChgBack==1     ) { $DNA_ZoF=SeqSegmentsTools::Chang_ZoF_BackToOrgContig_ZoF($DNA_ZoF, $gnoSeg_ZoF); }			  
			  my $DNA_ZoF_one=1; if ($DNA_ZoF eq '-') {$DNA_ZoF_one=-1; }
			  $outHash->{'0_MatchArray'}->[$i]->{'2_DNA_ZF'}=$DNA_ZoF; 
			  #$outHash->{'0_MatchArray'}->[$i]->{'3_STP_PS'}=1 if (  defined ( $spliced_seqHash->{'Spliced-seg_modifiers'}->{'Spliced-seg-modifier'}->{'Spliced-seg-modifier_stop-codon-found'}->{'value'} )  );
			  
			  if (   (  defined ( $spliced_seqHash->{'Spliced-seg_modifiers'} )  ) && (  defined ( $spliced_seqHash->{'Spliced-seg_modifiers'}->{'Spliced-seg-modifier'} )  )   ){
			  	my $Spliced_seg_modifier_REF=$spliced_seqHash->{'Spliced-seg_modifiers'}->{'Spliced-seg-modifier'};
			    if    (  ref ( $Spliced_seg_modifier_REF ) eq 'ARRAY' )	{
			    	foreach my $ecHashRef ( @{$Spliced_seg_modifier_REF} ){
			    		if (  defined ( $ecHashRef->{'Spliced-seg-modifier_start-codon-found'} ) ) { $outHash->{'0_MatchArray'}->[$i]->{'3_0_STPe'}=1; }
			    		if (  defined ( $ecHashRef->{'Spliced-seg-modifier_stop-codon-found'} ) )  { $outHash->{'0_MatchArray'}->[$i]->{'3_1_STTe'}=1; }
			    	}
			    }
			    elsif (  ref ( $Spliced_seg_modifier_REF ) eq 'HASH' )	{
			    	if (  defined ( $Spliced_seg_modifier_REF->{'Spliced-seg-modifier_start-codon-found'} ) ) { $outHash->{'0_MatchArray'}->[$i]->{'3_0_STPe'}=1; }
			      if (  defined ( $Spliced_seg_modifier_REF->{'Spliced-seg-modifier_stop-codon-found'} ) )  { $outHash->{'0_MatchArray'}->[$i]->{'3_1_STTe'}=1; }
			    }
			  }
			  	
			  
			
			
			  if (  defined ( $spliced_seqHash->{'Spliced-seg_exons'}->{'Spliced-exon'} )  ){
			  	
			  	
			  	
			    my $TmpExonArrayHASH=$spliced_seqHash->{'Spliced-seg_exons'}->{'Spliced-exon'};
			    my $exonArray;
			    if    ( ref ($TmpExonArrayHASH) eq 'HASH' ){ 			    	$exonArray->[0]=$TmpExonArrayHASH; 			    }
			    elsif ( ref ($TmpExonArrayHASH) eq 'ARRAY'){            $exonArray=$TmpExonArrayHASH;               }
			    else  { my  $dieMsg_3= "DIE!!!!$warnMsgHead\n\$TmpExonArrayHASH=$TmpExonArrayHASH\nshound be a Hash or array!!\n\n\n";  print $dieMsg_3; die $dieMsg_3; }
			    if ( ref ($exonArray) eq 'ARRAY'){
			    	for (my $j=0; $j<@{ $exonArray }; $j++){
			    		
			    		my $chuckStart;  
			    		my $xmlExoStt=$exonArray->[$j]->{'Spliced-exon_genomic-start'}+1;  
			    		my $xmlExoEnd=$exonArray->[$j]->{'Spliced-exon_genomic-end'}+1;
			    		
			    		if ($gno_ChgBack == 1){        
                $xmlExoStt=SeqSegmentsTools::Chang_pos_BackToOrgContigPositionk($xmlExoStt, $gnoSeg_stt, $gnoSeg_end, $gnoSeg_ZoF);
      	        $xmlExoEnd=SeqSegmentsTools::Chang_pos_BackToOrgContigPositionk($xmlExoEnd, $gnoSeg_stt, $gnoSeg_end, $gnoSeg_ZoF);
      	      }
      	      my $exonSmStt= SeqSegmentsTools::getSmallOne($xmlExoStt, $xmlExoEnd); my $exonBgEnd=SeqSegmentsTools::getbigOne($xmlExoStt, $xmlExoEnd);
			    		my $realStt; my $realEnd; if ($DNA_ZoF_one>0){ ($realStt, $realEnd)=($exonSmStt,$exonBgEnd);} else { ($realStt, $realEnd)=($exonBgEnd,$exonSmStt); }
			    		$outHash->{'0_MatchArray'}->[$i]->{'4_ExoAry'}->[$j]->{'0_DNA_0stt'}=   $chuckStart=$realStt;
			    		$outHash->{'0_MatchArray'}->[$i]->{'4_ExoAry'}->[$j]->{'0_DNA_1end'}               =$realEnd;
			    		
			    		$outHash->{'0_MatchArray'}->[$i]->{'4_0_GstrcA'}->[$j]->{'0_SgmtHead'}=$realStt;
			    		$outHash->{'0_MatchArray'}->[$i]->{'4_0_GstrcA'}->[$j]->{'1_SgmtTail'}=$realEnd;
			    		$outHash->{'0_MatchArray'}->[$i]->{'4_0_GstrcA'}->[$j]->{'2_SgmtZorF'}=$outHash->{'0_MatchArray'}->[$i]->{'2_DNA_ZF'};
			    		$outHash->{'0_MatchArray'}->[$i]->{'4_0_GstrcA'}->[$j]->{'3_SgmtType'}='CDS_region';
			    		
			    		
			    		my $pep_stt_pos=$exonArray->[$j]->{'Spliced-exon_product-start'}->{'Product-pos'}->{'Product-pos_protpos'}->{'Prot-pos'}->{'Prot-pos_amin'}+1;
			    		my $pep_end_pos=$exonArray->[$j]->{'Spliced-exon_product-end'  }->{'Product-pos'}->{'Product-pos_protpos'}->{'Prot-pos'}->{'Prot-pos_amin'}+1;
			    		if ($pro_ChgBack == 1){
                $pep_stt_pos=SeqSegmentsTools::Chang_pos_BackToOrgContigPositionk($pep_stt_pos, $proSeg_stt, $proSeg_end, $proSeg_ZoF);
      	        $pep_end_pos=SeqSegmentsTools::Chang_pos_BackToOrgContigPositionk($pep_end_pos, $proSeg_stt, $proSeg_end, $proSeg_ZoF);
      	      }
			    		$outHash->{'0_MatchArray'}->[$i]->{'4_ExoAry'}->[$j]->{'1_PEP_0stt_0Pos'}=$pep_stt_pos;
			    		$outHash->{'0_MatchArray'}->[$i]->{'4_ExoAry'}->[$j]->{'1_PEP_1end_0Pos'}=$pep_end_pos;
			    		
			    		
			    		
			    		$outHash->{'0_MatchArray'}->[$i]->{'4_ExoAry'}->[$j]->{'1_PEP_3stt_1pfm'}=$exonArray->[$j]->{'Spliced-exon_product-start'}->{'Product-pos'}->{'Product-pos_protpos'}->{'Prot-pos'}->{'Prot-pos_frame'};
			    		$outHash->{'0_MatchArray'}->[$i]->{'4_ExoAry'}->[$j]->{'1_PEP_4end_1pfm'}=$exonArray->[$j]->{'Spliced-exon_product-end'  }->{'Product-pos'}->{'Product-pos_protpos'}->{'Prot-pos'}->{'Prot-pos_frame'};
			    		my $partial_or_complete; if (   (  defined ( $exonArray->[$j]->{'Spliced-exon_partial'}->{'value'} )  ) && ($exonArray->[$j]->{'Spliced-exon_partial'}->{'value'} eq 'false')   ){ $partial_or_complete='Partial'; } else {$partial_or_complete='complete';}
			    		$outHash->{'0_MatchArray'}->[$i]->{'4_ExoAry'}->[$j]->{'2_Splice_0PorC'}=$partial_or_complete;
			    		$outHash->{'0_MatchArray'}->[$i]->{'4_ExoAry'}->[$j]->{'2_Splice_1Site_0BefExon'}=$exonArray->[$j]->{'Spliced-exon_acceptor-before-exon'}->{'Splice-site'}->{'Splice-site_bases'} if (  defined ( $exonArray->[$j]->{'Spliced-exon_acceptor-before-exon'}->{'Splice-site'}->{'Splice-site_bases'} )  );
			    		$outHash->{'0_MatchArray'}->[$i]->{'4_ExoAry'}->[$j]->{'2_Splice_1Site_1AftExon'}=$exonArray->[$j]->{'Spliced-exon_donor-after-exon'    }->{'Splice-site'}->{'Splice-site_bases'} if (  defined ( $exonArray->[$j]->{'Spliced-exon_donor-after-exon'    }->{'Splice-site'}->{'Splice-site_bases'} )  );
			    		
			    		
			    		my $chunckArray; my $ckIdx=0; 
			    		my $alignChunkRef=$exonArray->[$j]->{'Spliced-exon_parts'}->{'Spliced-exon-chunk'}; 
			    		if (  ( ref ($alignChunkRef) eq 'HASH' ) || ( ref ($alignChunkRef) eq 'ARRAY' )  ){
			    		  if    ( ref ($alignChunkRef) eq 'HASH' ){
			    		  	if (  defined ( $alignChunkRef->{'Spliced-exon-chunk_diag'} )  ){} else { DirFileHandle::PrintAndWarnDumper ($alignChunkRef); my $dieMsg_3="DIE!!!!$warnMsgHead\n\$alignChunkRef->{'Spliced-exon-chunk_diag'}\nShound be Difined!!\n\n\n"; print $dieMsg_3; die $dieMsg_3; }
			    		  	my $ckLth=$chunckArray->[$ckIdx]->{'0_chunkLength'}=$alignChunkRef->{'Spliced-exon-chunk_diag'};    $chunckArray->[$ckIdx]->{'1_chunk__type'}='0_mch_1DnP';   
			    		  	          $chunckArray->[$ckIdx]->{'2_chunk___stt'}=$chuckStart;                                    $chunckArray->[$ckIdx]->{'3_chunk___end'}=$chuckStart+$DNA_ZoF_one*($ckLth-1); 
 		  	          $ckIdx++;
			    		  }
			    		  elsif ( ref ($alignChunkRef) eq 'ARRAY' ){
			    		  	my $lastChckStt=$chuckStart;  my $lastChckLth=0;
			    		  	for (my $k=0; $k<@{ $alignChunkRef }; $k++){
			    		  		my $chckHash=$alignChunkRef->[$k];    my $ckLth;          #DirFileHandle::PrintAndWarnDumper ($chckHash);    print "\n\n\$ckIdx=$ckIdx\n\n";
			    		  		if    (   (  defined ( $chckHash->{'Spliced-exon-chunk_diag'} )  )        && ( $chckHash->{'Spliced-exon-chunk_diag'}=~m/^\d+$/ )   )       { print "chunk_diag\n";        $ckLth=$chunckArray->[$ckIdx]->{'0_chunkLength'}=$chckHash->{'Spliced-exon-chunk_diag'};         print "\$chunckArray->[$ckIdx]->{'0_chunkLength'}=$chunckArray->[$ckIdx]->{'0_chunkLength'}\n"; $chunckArray->[$ckIdx]->{'1_chunk__type'}='0_mch_1DnP';  $chunckArray->[$ckIdx]->{'2_chunk___stt'}=$lastChckStt+$DNA_ZoF_one*$lastChckLth;    $chunckArray->[$ckIdx]->{'3_chunk___end'}=$lastChckStt+$DNA_ZoF_one*$lastChckLth+$DNA_ZoF_one*($ckLth-1); $lastChckStt=$chunckArray->[$ckIdx]->{'2_chunk___stt'};   $lastChckLth=$ckLth; }
			    		  		elsif (   (  defined ( $chckHash->{'Spliced-exon-chunk_match'} )  )       && ( $chckHash->{'Spliced-exon-chunk_match'}=~m/^\d+$/ )   )      { print "chunk_match\n";       $ckLth=$chunckArray->[$ckIdx]->{'0_chunkLength'}=$chckHash->{'Spliced-exon-chunk_match'};        print "\$chunckArray->[$ckIdx]->{'0_chunkLength'}=$chunckArray->[$ckIdx]->{'0_chunkLength'}\n"; $chunckArray->[$ckIdx]->{'1_chunk__type'}='0_mch_0ATG';  $chunckArray->[$ckIdx]->{'2_chunk___stt'}=$lastChckStt+$DNA_ZoF_one*$lastChckLth;    $chunckArray->[$ckIdx]->{'3_chunk___end'}=$lastChckStt+$DNA_ZoF_one*$lastChckLth+$DNA_ZoF_one*($ckLth-1); $lastChckStt=$chunckArray->[$ckIdx]->{'2_chunk___stt'};   $lastChckLth=$ckLth; }
			    		  		elsif (   (  defined ( $chckHash->{'Spliced-exon-chunk_genomic-ins'} )  ) && ( $chckHash->{'Spliced-exon-chunk_genomic-ins'}=~m/^\d+$/ )   ){ print "chunk_genomic-ins\n"; $ckLth=$chunckArray->[$ckIdx]->{'0_chunkLength'}=$chckHash->{'Spliced-exon-chunk_genomic-ins'};  print "\$chunckArray->[$ckIdx]->{'0_chunkLength'}=$chunckArray->[$ckIdx]->{'0_chunkLength'}\n"; $chunckArray->[$ckIdx]->{'1_chunk__type'}='1_ins_0Dna';  $chunckArray->[$ckIdx]->{'2_chunk___stt'}=$lastChckStt+$DNA_ZoF_one*$lastChckLth;    $chunckArray->[$ckIdx]->{'3_chunk___end'}=$lastChckStt+$DNA_ZoF_one*$lastChckLth+$DNA_ZoF_one*($ckLth-1); $lastChckStt=$chunckArray->[$ckIdx]->{'2_chunk___stt'};   $lastChckLth=$ckLth; }
			    		  		elsif (   (  defined ( $chckHash->{'Spliced-exon-chunk_product-ins'} )  ) && ( $chckHash->{'Spliced-exon-chunk_product-ins'}=~m/^\d+$/ )   ){ print "chunk_product-ins\n"; $ckLth=$chunckArray->[$ckIdx]->{'0_chunkLength'}=$chckHash->{'Spliced-exon-chunk_product-ins'};  print "\$chunckArray->[$ckIdx]->{'0_chunkLength'}=$chunckArray->[$ckIdx]->{'0_chunkLength'}\n"; $chunckArray->[$ckIdx]->{'1_chunk__type'}='1_ins_1Pep';  $chunckArray->[$ckIdx]->{'2_chunk___stt'}=$lastChckStt+$DNA_ZoF_one*$lastChckLth;    $chunckArray->[$ckIdx]->{'3_chunk___end'}=$lastChckStt+$DNA_ZoF_one*$lastChckLth;                                                                                  }
			    		  		else { DirFileHandle::PrintAndWarnDumper ($chckHash);  my $dieMsg_4= "DIE!!!!$warnMsgHead\n\$alignChunkRef->[$k]=\$chckHash=$chckHash\nshound be the array of length chunk!!\n\n\n"; print $dieMsg_4; die $dieMsg_4; }         
			    		  		
			    		  		print "\$chunckArray->[$ckIdx]->{'1_chunk__type'}=$chunckArray->[$ckIdx]->{'1_chunk__type'}\t\$chunckArray->[$ckIdx]->{'0_chunkLength'}=$chunckArray->[$ckIdx]->{'0_chunkLength'}\n\n\n\n";
			    		  		$ckIdx++;
			    		  	}
			    		  }
			    	  }
			    		else { my $dieMSG_6= "DIE!!!!$warnMsgHead\n\$exonArray->[$j]->{'Spliced-exon_parts'}->{'Spliced-exon-chunk'}=$exonArray->[$j]->{'Spliced-exon_parts'}->{'Spliced-exon-chunk'}\n\$alignChunkRef=$alignChunkRef\nshound be A Hash or Array REF!!\n\n\n"; print $dieMSG_6; die $dieMSG_6; }
			    		
			    		my $dnaExonLength; #DirFileHandle::PrintAndWarnDumper ($chunckArray);
			    		if ( ref ($chunckArray) eq 'ARRAY' ){  
			    			for (  my $ckIdx_2=0; $ckIdx_2<@{ $chunckArray }; $ckIdx_2++  ){ #print "111 \$ckIdx_2=$ckIdx_2\t\t"; print "\$chunckArray->[$ckIdx_2]->{'0_chunkLength'}=$chunckArray->[$ckIdx_2]->{'0_chunkLength'}\n";
			    				my $chkLen=$chunckArray->[$ckIdx_2]->{'0_chunkLength'}; 
			    				if (  ( $chunckArray->[$ckIdx_2]->{'1_chunk__type'} eq '0_mch_0ATG' ) ||( $chunckArray->[$ckIdx_2]->{'1_chunk__type'} eq '0_mch_1DnP' ) || ( $chunckArray->[$ckIdx_2]->{'1_chunk__type'} eq '1_ins_0Dna' )  ){$dnaExonLength+=$chunckArray->[$ckIdx_2]->{'0_chunkLength'}; } #print "\$dnaExonLength+=$chunckArray->[$ckIdx_2]->{'0_chunkLength'}=$dnaExonLength+=$chunckArray->[$ckIdx_2]->{'0_chunkLength'}\n
			    				print "222 \$ckIdx_2=$ckIdx_2\n";
			    				
			    				
			    				if (   ( $chunckArray->[$ckIdx_2]->{'1_chunk__type'} eq '0_mch_1DnP' ) && (   (  defined ( $chunckArray->[$ckIdx_2-1]->{'4_chunk_fmShft_PepIns'} )  ) &&  ( $chunckArray->[$ckIdx_2-1]->{'4_chunk_fmShft_PepIns'} ==1 )   )    ){   print "000 \$ckIdx_2=$ckIdx_2\n";
			    				  			
			    				  			$chunckArray->[$ckIdx_2]->{'4_chunk_fmShft_MatchP'}=1;
			    				  			$chunckArray->[$ckIdx_2]->{'4_chunk_fmShft_MatchP_SttFrm'}=($chunckArray->[$ckIdx_2-1]->{'4_chunk_fmShft_PepIns_EndFrm'}+1)%3;   #DirFileHandle::PrintAndWarnDumper ($chunckArray->[$ckIdx_2]);       print "000 3333 \$ckIdx_2=$ckIdx_2\n";
			    				  	    #$outHash->{'0_MatchArray'}->[$i]->{'1_3_Psdg'}=1;
			    				  	    #push @{ $outHash->{'0_MatchArray'}->[$i]->{'1_3_PsdArray'} }, $chunckArray->[$ckIdx_2];
			    				}
			    				
			    				if (   ($ckIdx_2==0) || (  $ckIdx_2==( @{ $chunckArray }-1 )  )   ) {     #print "333 \$ckIdx_2=$ckIdx_2\n";                     '1_chunk__type'
			    				  
			    				  if    (      ($ckIdx_2==0)                       && (  ( $chunckArray->[$ckIdx_2]->{'1_chunk__type'} eq '1_ins_0Dna' )  )   ){	
			    				  	$chunckArray->[$ckIdx_2]->{'5_chunkDin0hd'}=1;  
			    				  	my $strgNB=$chkLen%3; if ($strgNB==0){$strgNB=3;} $strgNB=4-$strgNB; my $chkSttFrm=$strgNB;
			    				  	$chunckArray->[$ckIdx_2]->{'5_chunkDin0hd_stt_frm'}=$strgNB;
			    				  	$outHash->{'0_MatchArray'}->[$i]->{'4_ExoAry'}->[$j]->{'0_DNA_3stt_frm'}               =$strgNB;
			    				  }
			    				  elsif (  (  $ckIdx_2==( @{ $chunckArray }-1 )  ) && (  ( $chunckArray->[$ckIdx_2]->{'1_chunk__type'} eq '1_ins_0Dna' )  )   ){	
			    				  	$chunckArray->[$ckIdx_2]->{'5_chunkDin1tl'}=1; 
			    				  	my $strgNB=$chkLen%3; if ($strgNB==0){$strgNB=3;}  my $realEndFrm=$strgNB;
			    		        $chunckArray->[$ckIdx_2]->{'5_chunkDin1tl_end_frm'}=$strgNB;
			    		        $outHash->{'0_MatchArray'}->[$i]->{'4_ExoAry'}->[$j]->{'0_DNA_4end_frm'}               =$strgNB;
			    				  }
			    				}
			    				else {  print "444 \$ckIdx_2=$ckIdx_2\n";
			    				  if ( $chunckArray->[$ckIdx_2]->{'0_chunkLength'}%3 == 0 ){ print "555 \$ckIdx_2=$ckIdx_2\n";  }           
			    				  else { print "666 \$ckIdx_2=$ckIdx_2\n";  print "\$chunckArray->[$ckIdx_2]->{'1_chunk__type'}=$chunckArray->[$ckIdx_2]->{'1_chunk__type'}\n";
			    				  	if ($chunckArray->[$ckIdx_2]->{'1_chunk__type'} eq '1_ins_0Dna'){  #print "777 \$ckIdx_2=$ckIdx_2\n";
			    				  	  $chunckArray->[$ckIdx_2]->{'4_chunk_fmShft_DnaIns'}=1;
			    				  	  $chunckArray->[$ckIdx_2]->{'0_SgmtHead'}=$chunckArray->[$ckIdx_2]->{'2_chunk___stt'};
			    				  	  $chunckArray->[$ckIdx_2]->{'1_SgmtTail'}=$chunckArray->[$ckIdx_2]->{'3_chunk___end'};
			    				  	  $chunckArray->[$ckIdx_2]->{'2_SgmtZorF'}=$DNA_ZoF;
			    				  	  $chunckArray->[$ckIdx_2]->{'3_SgmtType'}='frame_shfit_insert';
			    				  	  $outHash->{'0_MatchArray'}->[$i]->{'1_3_Psdg'}=1;
			    				  	  push @{ $outHash->{'0_MatchArray'}->[$i]->{'1_3_PsdArray'} }, $chunckArray->[$ckIdx_2];
			    				  	}
			    				  	else{  print "888 \$ckIdx_2=$ckIdx_2\n";
			    				  		if    (   ( $chunckArray->[$ckIdx_2]->{'1_chunk__type'} eq '1_ins_1Pep' ) && (   (  defined ( $chunckArray->[$ckIdx_2+1]->{'1_chunk__type'} )  ) &&  ( $chunckArray->[$ckIdx_2+1]->{'1_chunk__type'}  eq '0_mch_1DnP' )   )    ){   print "999 \$ckIdx_2=$ckIdx_2\n";
			    				  			$chunckArray->[$ckIdx_2]->{'4_chunk_fmShft_PepIns'}=1;
			    				  			$chunckArray->[$ckIdx_2]->{'0_SgmtHead'}=$chunckArray->[$ckIdx_2]->{'2_chunk___stt'};
			    				  	    $chunckArray->[$ckIdx_2]->{'1_SgmtTail'}=$chunckArray->[$ckIdx_2]->{'3_chunk___end'};
			    				  	    $chunckArray->[$ckIdx_2]->{'2_SgmtZorF'}=$DNA_ZoF;
			    				  	    $chunckArray->[$ckIdx_2]->{'3_SgmtType'}='frame_shfit_insert';
			    				  	    $outHash->{'0_MatchArray'}->[$i]->{'1_3_Psdg'}=1;
			    				  	    push @{ $outHash->{'0_MatchArray'}->[$i]->{'1_3_PsdArray'} }, $chunckArray->[$ckIdx_2];
			    				  	    my $mdFrm=$chkLen%3; if ($mdFrm==0){$mdFrm=3;}  
			    				  	    $chunckArray->[$ckIdx_2]->{'4_chunk_fmShft_PepIns_EndFrm'}=$mdFrm;                                                               #DirFileHandle::PrintAndWarnDumper ($chunckArray->[$ckIdx_2]);       print "999 2222 \$ckIdx_2=$ckIdx_2\n";
			    				  		}
			    				  		elsif (   ( $chunckArray->[$ckIdx_2]->{'1_chunk__type'} eq '0_mch_1DnP' ) && (   (  defined ( $chunckArray->[$ckIdx_2]->{'4_chunk_fmShft_MatchP'} )  ) &&  ( $chunckArray->[$ckIdx_2]->{'4_chunk_fmShft_MatchP'} ==1 )   )    ){
			    				  			     
			    				  		}
			    				  		
			    				  		else{
			    				  		  my $dieMsghere="DIE!!!!$warnMsgHead\n\$chunckArray->[$ckIdx_2]->{'0_chunkLength'}=$chunckArray->[$ckIdx_2]->{'0_chunkLength'}\nshound be a 3 times number if it is not the head chunck or tail chunck, and a match or protein gap ins!!\n\n\n";	
			    				  		  print $dieMsghere; die $dieMsghere;
			    				  		}
			    				  	  
			    				  	}
			    				  }
			    				}
			    				
			    			}
			    		} 
			    		$outHash->{'0_MatchArray'}->[$i]->{'4_ExoAry'}->[$j]->{'0_DNA_3len'}=$dnaExonLength;
			    		$outHash->{'0_MatchArray'}->[$i]->{'4_ExoAry'}->[$j]->{'3_chuckArray'}=$chunckArray;
			    	}
			    	$outHash->{'0_MatchArray'}->[$i]->{'5_PepCov_1ary'}=&BuildPEPcoverageAnlysis( $outHash->{'0_MatchArray'}->[$i]->{'4_ExoAry'} );  
			    	$outHash->{'0_MatchArray'}->[$i]->{'5_PepCov_0len'}=&GetCuttedLength_proSplign_version($outHash->{'0_MatchArray'}->[$i]->{'5_PepCov_1ary'}); 
			    	
			    }	
			    else {
			    	DirFileHandle::PrintAndWarnDumper ($inArray);
			    	DirFileHandle::PrintAndWarnDumper ($exonArray);
			  	  my $dieMsghere1="DIE!!!!$warnMsgHead\n\$exonArray=$exonArray\nshound be defined!!\n\n\n";
			  	  print $dieMsghere1; die $dieMsghere1;
			    }
			  }
			  else{
			  	DirFileHandle::PrintAndWarnDumper ($inArray);
			  	my $dieMsghere2= "DIE!!!!$warnMsgHead\n\$spliced_seqHash->{'Spliced-seg_exons'}->{'Spliced-exon'}\nshound be defined!!\n\n\n";
			  	print $dieMsghere2; die $dieMsghere2;
			  }
			  
			}
			
			
			
		}
	}
	else {
		#my $dieMsghere3=  "DIE!!!!$warnMsgHead\n\$inArray=$inArray\nshound be a array ref!!\n\n\n";
		#print $dieMsghere3; die $dieMsghere3;
	}
	
		
	return $outHash;
}



sub Build_DNA_PEP_correlate_HASH{
	my $warnMsgHead="\n\n\n   In package  ProSplignHandle,\nIn sub Build_DNA_PEP_correlate_HASH,\n\n";
	my ($inProSpliHash,   $PEPSeqString, $DNASegString, $IN_tgaPos)=@_; # $PEPfastDB, $DNAfastaDB)=@_;
	
	my $DNASeqString_length=length ($DNASegString);  
	my $PEPSeqString_length=length ($PEPSeqString);
	
	my $outPutHash;
	if (   ( ref ($inProSpliHash) eq 'HASH' ) && (  ref ( $inProSpliHash->{'0_MatchArray'} ) eq 'ARRAY'  )   ){
		for (  my $i=0; $i<@{ $inProSpliHash->{'0_MatchArray'} }; $i++  ){
			
			$inProSpliHash->{'0_MatchArray'}->[$i]->{'0_1_aln_0_0_DNA_lengt'}=$DNASeqString_length;
			$inProSpliHash->{'0_MatchArray'}->[$i]->{'0_1_aln_0_2_PEP_lengt'}=$PEPSeqString_length;
			
			
			
			
			my $DNA_to_PEP_hash;
	    my $PEP_to_DNA_hash;
			
		  my $PEPnum=$inProSpliHash->{'0_MatchArray'}->[$i]->{'1_0PepNm'};  #print "\$PEPnum=$inProSpliHash->{'0_MatchArray'}->[$i]->{'1_0PepNm'}=$PEPnum=$inProSpliHash->{'0_MatchArray'}->[$i]->{'1_0PepNm'}\n";
			my $DNAnum=$inProSpliHash->{'0_MatchArray'}->[$i]->{'1_1DnaNm'};
			my $DNAZoF=$inProSpliHash->{'0_MatchArray'}->[$i]->{'2_DNA_ZF'};  my $ZoF_one=1; if ($DNAZoF eq '-'){$ZoF_one=-1;}
			
			if (  ref ( $inProSpliHash->{'0_MatchArray'}->[$i]->{'4_ExoAry'} ) eq 'ARRAY'  ){
				my $exonArray=$inProSpliHash->{'0_MatchArray'}->[$i]->{'4_ExoAry'}; #print "\$exonArray=\$inProSpliHash->{'0_MatchArray'}->[\$i]->{'4_ExoAry'}=\$inProSpliHash->{'0_MatchArray'}->[$i]->{'4_ExoAry'}=$exonArray\n";
				for (  my $j=0; $j<@{ $exonArray }; $j++  ){                        #print "\$exonArray=$exonArray\n"; DirFileHandle::PrintAndWarnDumper ($exonArray); 
					my $exonHash=$exonArray->[$j];                                    #print "\$exonHash=\$exonArray->[$j]=$exonHash\n"; DirFileHandle::PrintAndWarnDumper ($exonHash);  print "\$exonHash=\$exonArray->[$j]=$exonHash\n"; 
					#print "xxxxxxxxxxxxxxxxxxx";
					my ( $DNAstt,                   $DNAend,                   $PEPstt,                        $PEPend,                        $sttFrm,                        $endFrm,                        )
					  =( $exonHash->{'0_DNA_0stt'}, $exonHash->{'0_DNA_1end'}, $exonHash->{'1_PEP_0stt_0Pos'}, $exonHash->{'1_PEP_1end_0Pos'}, $exonHash->{'1_PEP_3stt_1pfm'}, $exonHash->{'1_PEP_4end_1pfm'}  ); 
					#print "yyyyyyyyyyyyyyyyyyyyyy";
					my  $chkAry ;                    #DirFileHandle::PrintAndWarnDumper ( $exonHash->{'3_chuckArray'} ) ;
					$chkAry=$exonHash->{'3_chuckArray'};
					if ( ref ($chkAry) eq 'ARRAY' ){
						my $ck_Number=@{ $chkAry };  
						
				
						
						
						####################### write the PepFrm value, write pep correlated frame for each DNA postion###########################
						my $chkSttFrm;  my $chkEndFrm;  my $realEndFrm=$endFrm;
						for (my $k=0; $k<$ck_Number; $k++){                                                                                        print "\$k=$k\n";
			    		my $chkHsh=$chkAry->[$k];
			    		my ( $chkLen,                    $chkTyp,                    $chkStt,                    $chkEnd                    )
			    		  =( $chkHsh->{'0_chunkLength'}, $chkHsh->{'1_chunk__type'}, $chkHsh->{'2_chunk___stt'}, $chkHsh->{'3_chunk___end'} );   print "$chkLen,                    $chkTyp,                    $chkStt,                    $chkEnd \n";
			    		
			    		if    (   (  defined( $chkHsh->{'4_chunk_fmShft_MatchP_SttFrm'} )  ) && ($chkHsh->{'4_chunk_fmShft_MatchP_SttFrm'}=~m/\d+/)   ){
			    			$chkSttFrm=$chkHsh->{'4_chunk_fmShft_MatchP_SttFrm'}; 
			    		}
			    		elsif    (   (  defined( $chkHsh->{'5_chunkDin0hd'} )  ) && ($chkHsh->{'5_chunkDin0hd'}==1)   ){ 
			    		  if ($k!=0){  DirFileHandle::PrintAndWarnDumper($chkHsh); my $dieMsg4="DIE!!!!$warnMsgHead\n\$k=$k should be 0 !!!!!\n\n\n"; print $dieMsg4; die $dieMsg4;}			    		  
			    		  $chkSttFrm=$chkHsh->{'5_chunkDin0hd_stt_frm'};  
			    		}
			    		elsif ($k==0){$chkSttFrm=$sttFrm;}
			    		else         {$chkSttFrm=1;}                                                                     print "\$chkSttFrm=$chkSttFrm\n";
			    		
			    	  if (   (  defined( $chkHsh->{'5_chunkDin1tl'} )  ) && ($chkHsh->{'5_chunkDin1tl'}==1)   ){	
			    	    if ( $k!=($ck_Number-1) ){my $dieMsg5="DIE!!!!$warnMsgHead\n\$k=$k should be equal to \$ck_Number-1=".($ck_Number-1)." !!!!!\n\n\n";print $dieMsg5; die $dieMsg5;}			    		  
			    		  $realEndFrm=$chkHsh->{'5_chunkDin1tl_end_frm'};
			    	  }
			    	  
			    		
			    		
			    		
			    		if  ( ($chkTyp eq '0_mch_0ATG') || ($chkTyp eq '0_mch_1DnP') || ($chkTyp eq '1_ins_0Dna') ){  
			    		  for (my $stepIdx=0; $stepIdx<$chkLen; $stepIdx++){
			    		  	my $strangNB=4-$chkSttFrm; if ($strangNB==3){$strangNB=0;} print "\$strangNB=$strangNB\n";
			    		  	my $orphanDNAnb=$stepIdx+1-$strangNB;  print "\$orphanDNAnb=\$stepIdx+1-\$strangNB=$stepIdx+1-$strangNB=$orphanDNAnb\n";#if ($orphanDNAnb<0){$orphanDNAnb+=3;}
			    		  	
			    		  	my $eachFrm=$orphanDNAnb%3; if ($eachFrm==0){$eachFrm=3;}      $chkEndFrm=$eachFrm;  print "\$chkEndFrm=\$eachFrm=$chkEndFrm\n";
			    		  	my $DNApos=$chkStt+$stepIdx*$ZoF_one;                                                          			    	$DNA_to_PEP_hash->{$DNApos}->{'PepFrm'}=$eachFrm;    $DNA_to_PEP_hash->{$DNApos}->{'6Ex_numb'}=$j+1;
			    		  	if (  ($k == 0)               && ($sttFrm != 1) && (   $stepIdx    < ( 3-($sttFrm-1) )  )   ){			    		  		$DNA_to_PEP_hash->{$DNApos}->{'Splice_Exon_head_or_not'}=1;			    		  	}
			    		  	if (  ($k == ($ck_Number-1) ) && ($endFrm != 3) && (  ($stepIdx+1) > ($chkLen-$endFrm)  )   ){			    		  		$DNA_to_PEP_hash->{$DNApos}->{'Splice_Exon_tail_or_not'}=1;			    		  	}
			    		  	
			    		  	
			    		  }
			    		}
			    	  			    		
			    	}
			    	if ($chkEndFrm eq $realEndFrm){}else {my $dieMsg6= "DIE!!!!$chkEndFrm\n\$chkEndFrm=$chkEndFrm == $realEndFrm=\$realEndFrm\nThese 2 frame should be equall!!!!!\n\n\n"; print $dieMsg6; die $dieMsg6;}
			    	####################### write the PepFrm value, write pep correlated frame for each DNA postion###########################
			    	
			    	####################### Add postion correlation to hash ################################################## Add postion correlation to hash ################################################## Add postion correlation to hash ###########################
						
						my $pepPos=$PEPstt;  my $lastChunkType;  print "000 \$pepPos=$PEPstt\n";
						for (my $k=0; $k<$ck_Number; $k++){
			    		my $chkHsh=$chkAry->[$k];
			    		my ( $chkLen,                    $chkTyp,                    $chkStt,                    $chkEnd                    )
			    		  =( $chkHsh->{'0_chunkLength'}, $chkHsh->{'1_chunk__type'}, $chkHsh->{'2_chunk___stt'}, $chkHsh->{'3_chunk___end'} );
			    	  
			    	  
			    	  
			    	  if  ( ($chkTyp eq '0_mch_0ATG') || ($chkTyp eq '0_mch_1DnP') ){  	 
			    	  	my $lastFrm=3;
			    		  for (my $stepIdx=0; $stepIdx<$chkLen; $stepIdx++){
			    	      my $DNApos=$chkStt+$stepIdx*$ZoF_one;
			    	      $DNA_to_PEP_hash->{$DNApos}->{'Type'}='NOTaGap';     ###################
			    	      
			    	      my $eachFrm=$DNA_to_PEP_hash->{$DNApos}->{'PepFrm'}; print "\$eachFrm=\$DNA_to_PEP_hash->{\$DNApos}->{'PepFrm'}=$eachFrm=\$DNA_to_PEP_hash->{$DNApos}->{'PepFrm'}\n"; #$lastFrm=$eachFrm;  
			    	      
			    	      #下面第一个if是 chuck的非完整3codon开头定位类型
			    	      if (   (  defined ( $DNA_to_PEP_hash->{$DNApos}->{'Splice_Exon_head_or_not'} )  ) && ( $DNA_to_PEP_hash->{$DNApos}->{'Splice_Exon_head_or_not'} ==1  )   ){
			    	      	$DNA_to_PEP_hash->{$DNApos}->{'PepCorPos'}=$pepPos;   
			    	        
			    	        
			    	        $PEP_to_DNA_hash->{$pepPos}->{'6Ex_numb'}=$j+1;
			    	        $PEP_to_DNA_hash->{$pepPos}->{'real_1_3_Pos'}->{$eachFrm}=$DNApos; print "1111111 \$DNA_to_PEP_hash->{$DNApos}->{'PepCorPos'}=\$pepPos=$pepPos;   \$PEP_to_DNA_hash->{$pepPos}->{'real_1_3_Pos'}->{$eachFrm}=\$DNApos=$DNApos\n";
			    	        $PEP_to_DNA_hash->{$pepPos}->{'Splice_Exon_head_or_not'}=1;
			    	        $PEP_to_DNA_hash->{$pepPos}->{'Type'}='NOTaGap';
			    	      }
			    	      else{
			    	        if ( ($stepIdx==0)  && ( $k>0 ) ){  #如果不是第一个chuck，那么是中间的chunck，而且是 match的类型，则在第一个位置 将peppos加一
			    	        	if ( ($k==1) && ($lastChunkType eq '1_ins_0Dna') ){   } #如果是第二个chuck，那么如果第一个chunk是 genome ins的话，则peppos也不加一
			    	        	else{ $pepPos++;  print "111  \$pepPos=$pepPos\t"; }
			    	          
			    	        }
			    	        else {
			    	        	if ($stepIdx>0){ print "\$ZoF_one=$ZoF_one\t\$eachFrm=$eachFrm\t$lastFrm=\$lastFrm\n";  #如果 既不是非完整的3codon开头，也不是chuck的第一个位置，则该frame比上一个frame小，则peppos加一
			    	        		if ( ($eachFrm) <= ($lastFrm) ){
			    	        			
			    	        			$pepPos++;  print "222  \$pepPos=$pepPos\t";
			    	        		}
			    	        	}
			    	        }
			    	        $DNA_to_PEP_hash->{$DNApos}->{'PepCorPos'}=$pepPos;   $PEP_to_DNA_hash->{$pepPos}->{'real_1_3_Pos'}->{$eachFrm}=$DNApos;  print "2222222 \$DNA_to_PEP_hash->{$DNApos}->{'PepCorPos'}=\$pepPos=$pepPos;   \$PEP_to_DNA_hash->{$pepPos}->{'real_1_3_Pos'}->{$eachFrm}=\$DNApos=$DNApos\n";
			    	      }
			    	      
			    	      
			    	      if (   (  defined ( $DNA_to_PEP_hash->{$DNApos}->{'Splice_Exon_tail_or_not'} )  ) && ( $DNA_to_PEP_hash->{$DNApos}->{'Splice_Exon_tail_or_not'} ==1  )   ){
			    	      	$PEP_to_DNA_hash->{$pepPos}->{'Splice_Exon_tail_or_not'}=1;
			    	      }
			    	      
			    	      
			    	      $lastFrm=$eachFrm;
			    	      
			    	    }
			    	    
			    	  }
			    	  
			    	  
			    	  if  ($chkTyp eq '1_ins_0Dna'){  
			    	  	
			    	  	my $FrameShiftIns=0;
			    	  	if (   (  defined ( $chkHsh->{'4_chunk_fmShft'} )  ) && ( $chkHsh->{'4_chunk_fmShft'}==1 )   ){   $FrameShiftIns=1;		    	  					    	  	}
			    	  	
			    		  for (my $stepIdx=0; $stepIdx<$chkLen; $stepIdx++){
			    		  	my $DNApos=$chkStt+$stepIdx*$ZoF_one;
			    		  	$DNA_to_PEP_hash->{$DNApos}->{'Type'}='ITSaGAP';
			    		  	if ($FrameShiftIns){
			    		  		$DNA_to_PEP_hash->{$DNApos}->{'FramShiftIns'}=1;
			    		  	}
			    		  	#my $eachFrm=$DNA_to_PEP_hash->{$DNApos};  
			    	      $DNA_to_PEP_hash->{$DNApos}->{'PepCorPos_before_gap'}=$pepPos;    print "333333 \$DNA_to_PEP_hash->{$DNApos}->{'PepCorPos_before_gap'}=\$pepPos=$pepPos\n";
			    	    }  			    	      
			    	  }
			    	  if  ($chkTyp eq '1_ins_1Pep'){ 
			    	  	#if ($chkLen%3==0){}else{die "DIE!!!!$warnMsgHead\n\$chkLen=$chkLen should be a 3 times number!!!!!\n\n\n"; } 
			    		  for (my $pPs=0; $pPs<(int($chkLen/3)); $pPs++){
			    		  	
			    		  	if ( ($k==0) && ($pPs==0) ){} #如果k=0，则$pepPos就等于开头的位置，不加
			    		  	else { $pepPos++; }
			    		  	
			    		  	$PEP_to_DNA_hash->{$pepPos}->{'Type'}='ITSaGAP';    print "333  \$pepPos=$pepPos\t";               
			    		  	print "444444 \$PEP_to_DNA_hash->{$pepPos}->{'Type'}='ITSaGAP'\n";
			    		  	
			    		  }
			    	      
			    	  }
			    	  
			    	  
			    	  $lastChunkType=$chkTyp;
			    	  ####################### Add postion correlation to hash ################################################## Add postion correlation to hash ################################################## Add postion correlation to hash ###########################
						
			    		
			    	}
			    }
			  }
			}
			
			my $matchedPEPstring;
			if ( ref ($PEP_to_DNA_hash) eq 'HASH'){
	    	foreach my $PepPos (    sort {$a<=>$b} (   keys (  %{ $PEP_to_DNA_hash }  )   )    ){ #print "\$PepPos=$PepPos\n";
	    		my $PePChar=FastaFileHandle::FactchSeqJustSubStr($PEPSeqString, '+', $PepPos, $PepPos);
	    		$PePChar=~s/\n//g; $PePChar=~s/\s//g;  $PEP_to_DNA_hash->{$PepPos}->{'PePChar'}=$PePChar;
	    		$matchedPEPstring.=$PePChar;
	    	}	    	
	    }
			
			my $proSplginMatchDNAstring;  my $proSplginMatchPEPstring;  my $proSplginMatchDNA_frm1_string; my $proSplginMatchDNA_frm2_string; my $proSplginMatchDNA_frm3_string;
			my $lastFrm; my $dnaIdx=0; my $allTHreeFrmTest=0; my @threeCodonPosArray; my $threeCodonIdx=0; my $threeCodonSTr;
			if ( ref ($DNA_to_PEP_hash) eq 'HASH'){
				my @keyArrayOfDNA=keys (  %{ $DNA_to_PEP_hash }  ); if ($DNAZoF eq '+'){  @keyArrayOfDNA= sort {$a<=>$b} @keyArrayOfDNA; } else  {  @keyArrayOfDNA= sort {$b<=>$a} @keyArrayOfDNA; } 
				my $lastDnaPos; my $lastDNAchr; my $arraySize=@keyArrayOfDNA;   my $TransPEPNumb=0;
				my $frm1FisrtPos_found=0; my $frm1FirstPos; my $frm3Last_pos;
	    	foreach my $DnaPos ( @keyArrayOfDNA ){  print "\$DnaPos=$DnaPos\n";
	    		my $DNAChar=FastaFileHandle::FactchSeqJustSubStr($DNASegString, $DNAZoF, $DnaPos, $DnaPos);
	    		$DNAChar=~s/\n//g; $DNAChar=~s/\s//g; $DNAChar=uc $DNAChar;
	    		$DNA_to_PEP_hash->{$DnaPos}->{'DNAChar'}=$DNAChar;
	    		$DNA_to_PEP_hash->{$DnaPos}->{'DNA_ZoF'}=$DNAZoF;
	    		
	    		##############DNA seq and PEP seq build step######################
	    		my $thisFrm=$DNA_to_PEP_hash->{$DnaPos}->{'PepFrm'};  print "\$DNAZoF=$DNAZoF\n\$thisFrm=\$DNA_to_PEP_hash->{\$DnaPos}->{'PepFrm'}=\$DNA_to_PEP_hash->{$DnaPos}->{'PepFrm'}=$thisFrm\n";
	    		
	    		if ( ($thisFrm==1)  && ($frm1FisrtPos_found==0) ){$frm1FirstPos=$DnaPos;$frm1FisrtPos_found++;} #get the first frame1 pos used to expand est
	    		if (  $thisFrm==3                               ){$frm3Last_pos=$DnaPos;                      } #get each frame3 pos used to expand est
	    		if ($dnaIdx>0 ){ &FrmCheck($lastFrm,$thisFrm); }
	    		
	    		if(   (  defined ( $DNA_to_PEP_hash->{$DnaPos}->{'FramShiftIns'} )  ) && ( $DNA_to_PEP_hash->{$DnaPos}->{'FramShiftIns'}==1 )   ){
	    			$DNA_to_PEP_hash->{$DnaPos}->{'3codon'}=$DNAChar;	    		     $DNA_to_PEP_hash->{$DnaPos}->{'trasAA'}=' ';	 $proSplginMatchPEPstring.=' ';  $proSplginMatchDNA_frm1_string.=' ';     $proSplginMatchDNA_frm2_string.=' ';  $proSplginMatchDNA_frm3_string.=' ';      
	    		}
	    		elsif ( ($dnaIdx==0) && ($thisFrm==3) ){   #################  nXXX kind, n means 1 dna, X means 3 dna can be coded as amino acid or stop ################# ################# ################# ################# ################# ################# ################# ################# 
	    		  $DNA_to_PEP_hash->{$DnaPos}->{'3codon'}=$DNAChar;	    		     $DNA_to_PEP_hash->{$DnaPos}->{'trasAA'}=' ';	 $proSplginMatchPEPstring.=' ';  $proSplginMatchDNA_frm1_string.=' ';     $proSplginMatchDNA_frm2_string.=' ';  $proSplginMatchDNA_frm3_string.=$DNAChar;  $DNA_to_PEP_hash->{$DnaPos}->{'1stExonHeadNotComplete'}=1;
	    		}
	    		elsif ( ($dnaIdx==0) && ($thisFrm==2) ){   ################# nnXXX kind, n means 1 dna, X means 3 dna can be coded as amino acid or stop # The 1st n     ################# ################# ################# ################# ################# ################# ################# 
	    				    		                                                     $DNA_to_PEP_hash->{$DnaPos}->{'trasAA'}=' ';	                                 $proSplginMatchDNA_frm1_string.=' ';     $proSplginMatchDNA_frm2_string.=$DNAChar;                                        $DNA_to_PEP_hash->{$DnaPos}->{'1stExonHeadNotComplete'}=1;
	    			$lastDnaPos=$DnaPos; $lastDNAchr=$DNAChar;
	    		}
	    		elsif ( ($dnaIdx==1) && ($thisFrm==3) ){   ################# nnXXX kind, n means 1 dna, X means 3 dna can be coded as amino acid or stop # The 2nd n     ################# ################# ################# ################# ################# ################# #################  
	    			$DNA_to_PEP_hash->{$lastDnaPos}->{'3codon'}="$lastDNAchr$DNAChar";
	    			$DNA_to_PEP_hash->{$DnaPos}->{'3codon'}="$lastDNAchr$DNAChar"; $DNA_to_PEP_hash->{$DnaPos}->{'trasAA'}=' ';	 $proSplginMatchPEPstring.=' ';	    		                                                                        $proSplginMatchDNA_frm3_string=$DNAChar;  $DNA_to_PEP_hash->{$DnaPos}->{'1stExonHeadNotComplete'}=1;
	    		}
	    		elsif  (  ( $dnaIdx==($arraySize-1) ) && ($thisFrm==1)  ){   #################  XXXn  kind, n means 1 dna, X means 3 dna can be coded as amino acid or stop ################# ################# ################# ################# ################# ################# ################# #################  
	    			$DNA_to_PEP_hash->{$DnaPos}->{'3codon'}=$DNAChar;	    		     $DNA_to_PEP_hash->{$DnaPos}->{'trasAA'}=' ';	 $proSplginMatchPEPstring.=' ';  $proSplginMatchDNA_frm1_string.=$DNAChar;                                                                                $DNA_to_PEP_hash->{$DnaPos}->{'lastExonTailNotComplete'}=1;
	    	  }
	    	  elsif  (  ( $dnaIdx==($arraySize-2) ) && ($thisFrm==1)  ){   #################  XXXnn kind, n means 1 dna, X means 3 dna can be coded as amino acid or stop # The 1st n       ################# ################# ################# ################# ################# ################# #################   
	    	  		    		                                                     $DNA_to_PEP_hash->{$DnaPos}->{'trasAA'}=' ';	                                 $proSplginMatchDNA_frm1_string.=$DNAChar;                                                                                $DNA_to_PEP_hash->{$DnaPos}->{'lastExonTailNotComplete'}=1;
	    			$lastDnaPos=$DnaPos; $lastDNAchr=$DNAChar;
	    	  }
	    	  elsif  (  ( $dnaIdx==($arraySize-1) ) && ($thisFrm==2)  ){   #################  XXXnn kind, n means 1 dna, X means 3 dna can be coded as amino acid or stop # The 2nd n       ################# ################# ################# ################# ################# ################# #################    
	    	    $DNA_to_PEP_hash->{$lastDnaPos}->{'3codon'}="$lastDNAchr$DNAChar";
	    	    $DNA_to_PEP_hash->{$DnaPos}->{'3codon'}="$lastDNAchr$DNAChar"; $DNA_to_PEP_hash->{$DnaPos}->{'trasAA'}=' ';	 $proSplginMatchPEPstring.=' ';                                            $proSplginMatchDNA_frm2_string.=$DNAChar;                                      $DNA_to_PEP_hash->{$DnaPos}->{'lastExonTailNotComplete'}=1;
	    	  }
	    		else{	 #################  nnn nnn nnn nnn kind, nnn means  3 dna can be coded as amino acid or stop  ################# ################# ################# ################# ################# ################# #################  ################# ################# ################# ################# ################# ################# #################       
	    		  $threeCodonPosArray[$threeCodonIdx]=$DnaPos; $threeCodonIdx++;
	    		  $threeCodonSTr.=$DNAChar; 
	    		  
	    		  if    ($thisFrm==1){ $proSplginMatchDNA_frm1_string.=$DNAChar; }
	    		  elsif ($thisFrm==2){ $proSplginMatchDNA_frm2_string.=$DNAChar; }
	    		  elsif ($thisFrm==3){ $proSplginMatchDNA_frm3_string.=$DNAChar; }
	    		  
	    		  if ($thisFrm==3){	    		  		    		  	
	    		  	my $trasAA;
	    		  	foreach my $ecDpos (@threeCodonPosArray){
	    		  		$DNA_to_PEP_hash->{$ecDpos}->{'3codon'}=$threeCodonSTr;
	    		  		if ($threeCodonSTr eq 'TGA'){  #(  (  defined ( $DNA_to_PEP_hash->{$DnaPos}->{'PepCorPos'} )  ) && ( $DNA_to_PEP_hash->{$DnaPos}->{'PepCorPos'}=~m/\d+/ ) && (  defined ( $PEP_to_DNA_hash->{ $DNA_to_PEP_hash->{$DnaPos}->{'PepCorPos'} }->{'PePChar'} )  ) && ( $PEP_to_DNA_hash->{ $DNA_to_PEP_hash->{$DnaPos}->{'PepCorPos'} }->{'PePChar'}=~m/\S+/ )    ) ;
	    		  			if (     (  defined ( $DNA_to_PEP_hash->{$DnaPos}->{'PepCorPos'} )  ) && ( $DNA_to_PEP_hash->{$DnaPos}->{'PepCorPos'}=~m/\d+/ ) 
	    		  			      && (  defined ( $PEP_to_DNA_hash->{ $DNA_to_PEP_hash->{$DnaPos}->{'PepCorPos'} }->{'PePChar'}                  )        ) 
	    		  			      && (            $PEP_to_DNA_hash->{ $DNA_to_PEP_hash->{$DnaPos}->{'PepCorPos'} }->{'PePChar'}=~m/[ucUC]/ )    
	    		  			   )
	    		  			{  $DNA_to_PEP_hash->{$ecDpos}->{'trasAA'}=FastaFileHandle::CodonTable_TGA_U ($threeCodonSTr);   }
	    		  			else{
	    		  				 $DNA_to_PEP_hash->{$ecDpos}->{'trasAA'}=FastaFileHandle::CodonTable_TGA_J ($threeCodonSTr);
	    		  			}
	    		  			
	    		  		}
	    		  		else{
	    		  			$DNA_to_PEP_hash->{$ecDpos}->{'trasAA'}=FastaFileHandle::CodonTable_TAG_O_TAA_B ($threeCodonSTr); 
	    		  		}
	    		  		
	    		  	}
	    		  	
	    		  	$proSplginMatchPEPstring.=$DNA_to_PEP_hash->{$DnaPos}->{'trasAA'}; $TransPEPNumb++;
	    		  	
	    		  	if ($DNA_to_PEP_hash->{$DnaPos}->{'trasAA'} eq 'U'){
	    		  		if ($DNAChar ne 'A'){ die "\n\n\n\nDIE!!!!!!!!!!!!$warnMsgHead\n\$threeCodonSTr=$threeCodonSTr\n\$DNAChar=\$DNAChar\nShould Be A!!!!!!!!!!!!!\n\n\n\n\n\n"; }
	    		  		my $ck3codonArrRef=[@threeCodonPosArray]; if ($DNAZoF eq '-'){$ck3codonArrRef=[ reverse  @threeCodonPosArray ]; }
  	    	  		my $ContinuityCheckFor3Codon=&ContinuityCheckFor3Codon ( $ck3codonArrRef );
  	    	  		if ( $ContinuityCheckFor3Codon==0 ){  #如果找到U的三个密码子 是被内含子分开，则DIE！！  if the TGA codon of Sec is found be divided by intron, then DIE!!!!!!!
	    		  	    $inProSpliHash->{'0_MatchArray'}->[$i]->{'0_0_TGA_in_2_exon'}=$ck3codonArrRef; 
	    		  	   }
	    		  	  
	    		  		$inProSpliHash->{'0_MatchArray'}->[$i]->{'6_TGAchangeHASH'}->{$DnaPos}->{'0_directi'}=$DNAZoF;	
	    		  		$inProSpliHash->{'0_MatchArray'}->[$i]->{'6_TGAchangeHASH'}->{$DnaPos}->{'1_postion'}=$DnaPos;	
	    		  		$inProSpliHash->{'0_MatchArray'}->[$i]->{'6_TGAchangeHASH'}->{$DnaPos}->{'2_fromCha'}='A';
	    		  		$inProSpliHash->{'0_MatchArray'}->[$i]->{'6_TGAchangeHASH'}->{$DnaPos}->{'3_intoCha'}='C';
	    		  		$inProSpliHash->{'0_MatchArray'}->[$i]->{'6_TGAchangeHASH'}->{$DnaPos}->{'4_directi'}='prosplign, A of tga matched to '.$PEP_to_DNA_hash->{ $DNA_to_PEP_hash->{$DnaPos}->{'PepCorPos'} }->{'PePChar'}   if (  (  defined ( $DNA_to_PEP_hash->{$DnaPos}->{'PepCorPos'} )  ) && ( $DNA_to_PEP_hash->{$DnaPos}->{'PepCorPos'}=~m/\d+/ ) && (  defined ( $PEP_to_DNA_hash->{ $DNA_to_PEP_hash->{$DnaPos}->{'PepCorPos'} }->{'PePChar'} )  ) && ( $PEP_to_DNA_hash->{ $DNA_to_PEP_hash->{$DnaPos}->{'PepCorPos'} }->{'PePChar'}=~m/\S+/ )    ); #warn "\$PEP_to_DNA_hash->{ \$DNA_to_PEP_hash->{$DnaPos}->{'PepCorPos'} }->{'PePChar'}=\$PEP_to_DNA_hash->{ $DNA_to_PEP_hash->{$DnaPos}->{'PepCorPos'} }->{'PePChar'}=$PEP_to_DNA_hash->{ $DNA_to_PEP_hash->{$DnaPos}->{'PepCorPos'} }->{'PePChar'}\n";
	    		  		$inProSpliHash->{'0_MatchArray'}->[$i]->{'6_TGAchangeHASH'}->{$DnaPos}->{'5_proSpMc'}=$PEP_to_DNA_hash->{ $DNA_to_PEP_hash->{$DnaPos}->{'PepCorPos'} }->{'PePChar'} if (  (  defined ( $DNA_to_PEP_hash->{$DnaPos}->{'PepCorPos'} )  ) && ( $DNA_to_PEP_hash->{$DnaPos}->{'PepCorPos'}=~m/\d+/ ) && (  defined ( $PEP_to_DNA_hash->{ $DNA_to_PEP_hash->{$DnaPos}->{'PepCorPos'} }->{'PePChar'} )  ) && ( $PEP_to_DNA_hash->{ $DNA_to_PEP_hash->{$DnaPos}->{'PepCorPos'} }->{'PePChar'}=~m/\S+/ )    ) ;
	    		  		
	    		  		
	    		  		
	    		  		my $Sec_TGA_hash; my $Sec_TGA_pos_list_array;
	    		  		foreach my $ecDpos (@threeCodonPosArray){
	    		  		  $Sec_TGA_hash->{ $DNA_to_PEP_hash->{$ecDpos}->{'PepFrm'} }->{'DNAChar'}=$DNA_to_PEP_hash->{$ecDpos}->{'DNAChar'};
	    		  		  $Sec_TGA_hash->{ $DNA_to_PEP_hash->{$ecDpos}->{'PepFrm'} }->{'DNAPosi'}=$ecDpos;
	    		  		  $Sec_TGA_hash->{ $DNA_to_PEP_hash->{$ecDpos}->{'PepFrm'} }->{'DNAZorF'}=$DNAZoF;
	    		  	  }
	    		  	  push @{ $inProSpliHash->{'0_MatchArray'}->[$i]->{'3_TGA_0_DNA_codon_Array'} }, $Sec_TGA_hash;
	    		  		
	    		  		if ( $ContinuityCheckFor3Codon==1 ){
	    		  		  $Sec_TGA_pos_list_array->{'0_SgmtHead'}=$Sec_TGA_hash->{1}->{'DNAPosi'};
	    		  	  	$Sec_TGA_pos_list_array->{'1_SgmtTail'}=$Sec_TGA_hash->{3}->{'DNAPosi'};
	    		  	  	$Sec_TGA_pos_list_array->{'2_SgmtZorF'}=$DNAZoF;
	    		  	  	$Sec_TGA_pos_list_array->{'3_SgmtType'}='Sec_TGA';    		  			    		  		
	    		    		push @{ $inProSpliHash->{'0_MatchArray'}->[$i]->{'4_1_TGAary'} }, $Sec_TGA_pos_list_array;   		  		
	    		    	}
	    		    	else {
	    		    		my $ps_1=$threeCodonPosArray[0]; my $ps_2=$threeCodonPosArray[1]; my $ps_3=$threeCodonPosArray[2];
	    		    		$Sec_TGA_pos_list_array->{'0upStream_exon_part'}->{'2_SgmtZorF'}=$DNAZoF;  
	    		    		$Sec_TGA_pos_list_array->{'1DwStream_exon_part'}->{'2_SgmtZorF'}=$DNAZoF;  
	    		    		$Sec_TGA_pos_list_array->{'0upStream_exon_part'}->{'3_SgmtType'}='Sec_TGA';
	    		    		$Sec_TGA_pos_list_array->{'1DwStream_exon_part'}->{'3_SgmtType'}='Sec_TGA'; 
	    		    		
	    		    		$Sec_TGA_pos_list_array->{'0upStream_exon_part'}->{'0_SgmtHead'}=$ps_1;
	    		    		
	    		    		if (  ( abs ($ps_1-$ps_2) ) == 1  ){ $Sec_TGA_pos_list_array->{'0upStream_exon_part'}->{'1_SgmtTail'}=$ps_2; $Sec_TGA_pos_list_array->{'1DwStream_exon_part'}->{'0_SgmtHead'}=$ps_3;		}
	    		    		else                               { $Sec_TGA_pos_list_array->{'0upStream_exon_part'}->{'1_SgmtTail'}=$ps_1; $Sec_TGA_pos_list_array->{'1DwStream_exon_part'}->{'0_SgmtHead'}=$ps_2;		}
	    		    			    		  	  	
	    		  	  	$Sec_TGA_pos_list_array->{'1DwStream_exon_part'}->{'1_SgmtTail'}=$ps_3;
	    		  	  	
	    		  	  	push @{ $inProSpliHash->{'0_MatchArray'}->[$i]->{'4_1_TGAary'} }, $Sec_TGA_pos_list_array; 
	    		    	}
	    		    	
	    		    	
	    		    	my $SecPep_Hash; 
	    		    	my $mcType; $mcType=$DNA_to_PEP_hash->{ $threeCodonPosArray[0] }->{'Type'     }   if (   (  defined ( $DNA_to_PEP_hash->{ $threeCodonPosArray[0] }->{'Type'     } )  ) && ( $DNA_to_PEP_hash->{ $threeCodonPosArray[0] }->{'Type'     }=~m/\S+/ )   );    
	    		    	my $ncPpos; $ncPpos=$DNA_to_PEP_hash->{ $threeCodonPosArray[0] }->{'PepCorPos'}   if (   (  defined ( $DNA_to_PEP_hash->{ $threeCodonPosArray[0] }->{'PepCorPos'} )  ) && ( $DNA_to_PEP_hash->{ $threeCodonPosArray[0] }->{'PepCorPos'}=~m/\S+/ )   );
	    		  		foreach my $ecDpos (@threeCodonPosArray){
	    		  			my $ecType; $ecType=$DNA_to_PEP_hash->{ $ecDpos }->{'Type'     }  if (   (  defined ( $DNA_to_PEP_hash->{ $ecDpos }->{'Type'     } )  ) && ( $DNA_to_PEP_hash->{ $ecDpos }->{'Type'     }=~m/\S+/ )   );      
	    		  			if (   (  ( defined ($mcType) ) && ($mcType=~m/\S+/)  ) && (  ( defined ($ecType) ) && ($ecType=~m/\S+/)  ) && ($mcType ne $ecType)   ){
	    		  				my $DieMsg="DIE!!!!\n\n\n$warnMsgHead\n\$DNA_to_PEP_hash->{ $threeCodonPosArray[0] }->{'Type'}=$DNA_to_PEP_hash->{ $threeCodonPosArray[0] }->{'Type'}\n\$DNA_to_PEP_hash->{ $ecDpos }->{'Type'}=$DNA_to_PEP_hash->{ $ecDpos }->{'Type'}\nThese two type should be equal!!\n\n\n";
	    		  				print $DieMsg; die $DieMsg;
	    		  			}
	    		  			my $ecCpos; $ecCpos=$DNA_to_PEP_hash->{ $ecDpos }->{'PepCorPos'}  if (   (  defined ( $DNA_to_PEP_hash->{ $ecDpos }->{'PepCorPos'} )  ) && ( $DNA_to_PEP_hash->{ $ecDpos }->{'PepCorPos'}=~m/\S+/ )   );  
	    		  			if (   (  ( defined ($ncPpos) ) && ($ncPpos=~m/\S+/)  ) && (  ( defined ($ecCpos) ) && ($ecCpos=~m/\S+/)  ) && ($ncPpos ne $ecCpos)   ){
	    		  				my $DieMsg="DIE!!!!\n\n\n$warnMsgHead\n\$DNA_to_PEP_hash->{ $threeCodonPosArray[0] }->{'PepCorPos'}=$DNA_to_PEP_hash->{ $threeCodonPosArray[0] }->{'PepCorPos'}\n\$DNA_to_PEP_hash->{ $ecDpos }->{'PepCorPos'}=$DNA_to_PEP_hash->{ $ecDpos }->{'PepCorPos'}\nThese two postion should be equal!!\n\n\n";
	    		  				print $DieMsg; die $DieMsg;
	    		  			}
	    		  		}
	    		  		$SecPep_Hash->{'0_pep_positi'}=$TransPEPNumb;
	    		  		$SecPep_Hash->{'0_Pep_3codon'}=$DNA_to_PEP_hash->{$DnaPos}->{'3codon'};
	    		  		$SecPep_Hash->{'0_Pep_aminoA'}=$DNA_to_PEP_hash->{$DnaPos}->{'trasAA'};
	    		  		$SecPep_Hash->{'1_mcP_positi'}=$DNA_to_PEP_hash->{$DnaPos}->{'PepCorPos'};
	    		  		$SecPep_Hash->{'1_mcP_aminoA'}='-'; $SecPep_Hash->{'1_mcP_aminoA'}=$PEP_to_DNA_hash->{ $SecPep_Hash->{'1_mcP_positi'} }->{'PePChar'} if ($mcType eq 'NOTaGap');
	    		  		push @{ $inProSpliHash->{'0_MatchArray'}->[$i]->{'3_TGA_1_PEP_codon_Array'} }, $SecPep_Hash; 
	    		  		
	    		  		
	    		  	}
	    		  	
	    		  	@threeCodonPosArray=();
	    		  	$threeCodonIdx=0;
	    		  	$threeCodonSTr='';
	    		  }
	    		}
	    		
	    		$proSplginMatchDNAstring.=$DNAChar;
	    		
	    		$lastFrm=$thisFrm;
	    		$dnaIdx++;
	    		##############DNA seq and PEP seq build step######################
	    	}	
	    	$inProSpliHash->{'0_MatchArray'}->[$i]->{'0_1_aln_0_4_0_Dps_1stf1'}=$frm1FirstPos; 
	    	$inProSpliHash->{'0_MatchArray'}->[$i]->{'0_1_aln_0_5_0_Dps_Lstf3'}=$frm3Last_pos; 
	    	
	    	$inProSpliHash->{'0_MatchArray'}->[$i]->{'0_1_aln_0_4_1_3cdon_1st'}=$DNA_to_PEP_hash->{$frm1FirstPos}->{'3codon'};
	    	$inProSpliHash->{'0_MatchArray'}->[$i]->{'0_1_aln_0_5_1_3cdon_Lst'}=$DNA_to_PEP_hash->{$frm3Last_pos}->{'3codon'}; 
	    	
	    	$inProSpliHash->{'0_MatchArray'}->[$i]->{'0_1_aln_0_4_2_trsAA_1st'}=$DNA_to_PEP_hash->{$frm1FirstPos}->{'trasAA'};
	    	$inProSpliHash->{'0_MatchArray'}->[$i]->{'0_1_aln_0_5_2_trsAA_Lst'}=$DNA_to_PEP_hash->{$frm3Last_pos}->{'trasAA'};
	    	
	    	#if  ($lastFrm != 3)  	{die "\n\n\n\nDIE!!!!!!!!!!!!$warnMsgHead\n\$thisFrm=\$lastFrm=$lastFrm\nShould Be 3!!!!!!!!!!!!!\n\n\n\n\n\n";}
	    }
	    $inProSpliHash->{'0_MatchArray'}->[$i]->{'5_TGA_org_0_Ushw'}= FastaFileHandle::Build_only_Show_U_PEPstring($matchedPEPstring);
			$inProSpliHash->{'0_MatchArray'}->[$i]->{'5_TGA_org_1PEPsq'}= $matchedPEPstring;
			
			$inProSpliHash->{'0_MatchArray'}->[$i]->{'6_TGA_DNA0tr_0Sq'}= $proSplginMatchPEPstring; my $noSpaceSeq=$proSplginMatchPEPstring; $noSpaceSeq=~s/\n//g; $noSpaceSeq=~s/\s//g; my $proTransPEPLength=length $noSpaceSeq;
			$inProSpliHash->{'0_MatchArray'}->[$i]->{'6_TGA_DNA0tr_1Us'}= FastaFileHandle::Build_only_Show_U_PEPstring($proSplginMatchPEPstring);
			$inProSpliHash->{'0_MatchArray'}->[$i]->{'0_1_aln_0_3_Out_lengt'}=$proTransPEPLength;
			$inProSpliHash->{'0_MatchArray'}->[$i]->{'6_TGA_DNA1_1frmS'}= $proSplginMatchDNA_frm1_string;
			$inProSpliHash->{'0_MatchArray'}->[$i]->{'6_TGA_DNA1_2frmS'}= $proSplginMatchDNA_frm2_string;
			$inProSpliHash->{'0_MatchArray'}->[$i]->{'6_TGA_DNA1_3frmS'}= $proSplginMatchDNA_frm3_string;
			$inProSpliHash->{'0_MatchArray'}->[$i]->{'6_TGA_DNA1string'}= $proSplginMatchDNAstring;
			
			$inProSpliHash->{'0_MatchArray'}->[$i]->{'7_DNA_to_PEP_hash'}= $DNA_to_PEP_hash;               print "20181121-0-0 \n";
			$inProSpliHash->{'0_MatchArray'}->[$i]->{'7_0_represtCDSpos'}= ProSplignHandle::findTheRepresentCDSposARRAY($DNA_to_PEP_hash, 10); print "20181121-0-1 \n";
			$inProSpliHash->{'0_MatchArray'}->[$i]->{'7_PEP_to_DNA_hash'}= $PEP_to_DNA_hash;
			
			my $UUmatchAA_number=0;
			$inProSpliHash->{'0_MatchArray'}->[$i]->{'3_TGA_2_UU_MatchNum'}=$UUmatchAA_number=&Calculate_U_Match_number (  $inProSpliHash->{'0_MatchArray'}->[$i]->{'3_TGA_1_PEP_codon_Array'} ,['U']);
			$inProSpliHash->{'0_MatchArray'}->[$i]->{'3_TGA_2_UC_allMatch'}=                  &CheckAll_U_uc_Match      (  $inProSpliHash->{'0_MatchArray'}->[$i]->{'3_TGA_1_PEP_codon_Array'} );
			$inProSpliHash->{'0_MatchArray'}->[$i]->{'3_TGA_3_UCX*_MC_PCT'}=                  &Calculate_U_Match_pct    (  $inProSpliHash->{'0_MatchArray'}->[$i]->{'3_TGA_1_PEP_codon_Array'} , ['U', 'C' , 'X', '*']);
		
		  ###########################################
			if (  defined ( $inProSpliHash->{'0_MatchArray'}->[$i]->{'0_1_aln_1_ident'} )  ){ 
				my $idLth=$inProSpliHash->{'0_MatchArray'}->[$i]->{'0_1_aln_1_ident'}; my $dcm =$inProSpliHash->{'0_MatchArray'}->[$i]->{'0_1_aln_1_ident_0_dcm'}=(($idLth/3)+$UUmatchAA_number)/$PEPSeqString_length; 
				                                                                                $inProSpliHash->{'0_MatchArray'}->[$i]->{'0_1_aln_1_ident_1_pct'}=BlastHandle::ChangeTo100PercentNB ($dcm); 
			}
			if (  defined ( $inProSpliHash->{'0_MatchArray'}->[$i]->{'0_1_aln_2_postv'} )  ){ 
				my $psLth=$inProSpliHash->{'0_MatchArray'}->[$i]->{'0_1_aln_2_postv'}; my $dcm =$inProSpliHash->{'0_MatchArray'}->[$i]->{'0_1_aln_2_postv_0_dcm'}=(($psLth/3)+$UUmatchAA_number)/$PEPSeqString_length; 
				                                                                                $inProSpliHash->{'0_MatchArray'}->[$i]->{'0_1_aln_2_postv_1_pct'}=BlastHandle::ChangeTo100PercentNB ($dcm); 
			}			
			
			if (  defined ( $inProSpliHash->{'0_MatchArray'}->[$i]->{'5_PepCov_0len'} )  ){ 
				my $cvLth=$inProSpliHash->{'0_MatchArray'}->[$i]->{'5_PepCov_0len'};   my $dcm =$inProSpliHash->{'0_MatchArray'}->[$i]->{'5_PepCov_1dcm'}=$cvLth/$PEPSeqString_length; 
				                                                                                $inProSpliHash->{'0_MatchArray'}->[$i]->{'5_PepCov_1pct'}=BlastHandle::ChangeTo100PercentNB ($dcm); 
			  if (  defined ( $inProSpliHash->{'0_MatchArray'}->[$i]->{'0_1_aln_1_ident_0_dcm'} )  ){ 
			    my $Idt_dcm =$inProSpliHash->{'0_MatchArray'}->[$i]->{'0_1_aln_1_ident_2_cvg_part_dcm'}=$inProSpliHash->{'0_MatchArray'}->[$i]->{'0_1_aln_1_ident_0_dcm'}/$dcm; 
			                 $inProSpliHash->{'0_MatchArray'}->[$i]->{'0_1_aln_1_ident_2_cvg_part_pct'}=BlastHandle::ChangeTo100PercentNB ($Idt_dcm); 
			  }
			  if (  defined ( $inProSpliHash->{'0_MatchArray'}->[$i]->{'0_1_aln_2_postv_0_dcm'} )  ){ 
			    my $Idt_dcm =$inProSpliHash->{'0_MatchArray'}->[$i]->{'0_1_aln_2_postv_2_cvg_part_dcm'}=$inProSpliHash->{'0_MatchArray'}->[$i]->{'0_1_aln_2_postv_0_dcm'}/$dcm; 
			                 $inProSpliHash->{'0_MatchArray'}->[$i]->{'0_1_aln_2_postv_2_cvg_part_pct'}=BlastHandle::ChangeTo100PercentNB ($Idt_dcm); 
			  }
			}	
						
			###########################################
			
		
		}
		
		@{ $outPutHash->{'0_MatchArray'} }=sort {  ( $b->{'3_TGA_2_UU_MatchNum'} <=> $a->{'3_TGA_2_UU_MatchNum'} ) or ( $b->{'0_1_aln_1_ident_0_dcm'} <=> $a->{'0_1_aln_1_ident_0_dcm'} )  } @{ $inProSpliHash->{'0_MatchArray'} };
	  
	  my $finalOutHash;
	  foreach my $eachHashUnit (  @{ $outPutHash->{'0_MatchArray'} }  ){
	  	if (   
	  	           (  defined ( $eachHashUnit                                                  )  )  && (  ref ( $eachHashUnit                                      ) eq 'HASH'  ) 
	  	        && (  defined ( $eachHashUnit->{'7_DNA_to_PEP_hash'}                           )  )  && (  ref ( $eachHashUnit->{'7_DNA_to_PEP_hash'}               ) eq 'HASH'  )
	  	        && (  defined ( $eachHashUnit->{'7_DNA_to_PEP_hash'}->{$IN_tgaPos}             )  )  && (  ref ( $eachHashUnit->{'7_DNA_to_PEP_hash'}->{$IN_tgaPos} ) eq 'HASH'  )
	  	        && (  defined ( $eachHashUnit->{'7_DNA_to_PEP_hash'}->{$IN_tgaPos}->{'3codon'} )  )  && ( $eachHashUnit->{'7_DNA_to_PEP_hash'}->{$IN_tgaPos}->{'3codon'}=~m/\S+/ )
	  	        && (  defined ( $eachHashUnit->{'7_DNA_to_PEP_hash'}->{$IN_tgaPos}->{'PepFrm'} )  )  && ( $eachHashUnit->{'7_DNA_to_PEP_hash'}->{$IN_tgaPos}->{'3codon'}=~m/\d+/ )
	  	        && (  $eachHashUnit->{'7_DNA_to_PEP_hash'}->{$IN_tgaPos}->{'3codon'} == 'TGA'     ) 
	  	        && (  $eachHashUnit->{'7_DNA_to_PEP_hash'}->{$IN_tgaPos}->{'PepFrm'} == 1         )  
	  	){
	  		push @{ $finalOutHash->{'0_MatchArray'} }, $eachHashUnit;
	  	}
	  }  
		
		
	}

	 
	
	return $outPutHash;
	
}
sub FrmCheck{
	my ($lastFrm, $thisFrm)=@_;
	my $warnMsgHead="\n\n\n   In package  ProSplignHandle,\nIn sub FrmCheck,\n\n";;
	if    ($lastFrm=3, $thisFrm=1){}
	elsif ($lastFrm=2, $thisFrm=3){}
	elsif ($lastFrm=1, $thisFrm=2){}
	else  {my $dieMsg_7= "\n\n\n\nDIE!!!!!!!!!!!!$warnMsgHead\n\$lastFrm=$lastFrm\n\$thisFrm=$thisFrm\n\nShould be 123 in order\n\n\n\n\n\n"; print $dieMsg_7; die $dieMsg_7;}
}

sub ContinuityCheckFor3Codon{
	my ($codonArray)=@_;
	my $warnMsgHead="\n\n\n   In package  ProSplignHandle,\nIn sub ContinuityCheckFor3Codon,\n\n";;
	my $right=0;
	if ( ref($codonArray) eq 'ARRAY' ){
		if ( @{ $codonArray } < 2 ){$right=3; DirFileHandle::PrintAndWarnDumper ($codonArray); my $dieMsg_8="\n\n\n\nDIE!!!!!!!!!!!$warnMsgHead\n\$codonArray=$codonArray\n\nThis array Should be larger than 1\n\n\n\n\n\n"; print $dieMsg_8; die $dieMsg_8; }		
		my $lastPos; my $stp=0; ;
		FCMK: foreach my $posit (  @{ $codonArray }  ){
			if ($stp>0){
				if ( $posit == ($lastPos+1) ){
					$right=1;
				}
				else {$right=2; last FCMK; }
			}
			$lastPos=$posit;
			$stp++;
		}
		
	}
	if ($right != 1){ print "$warnMsgHead NOTICE!!!!! \$right=$right\n"; DirFileHandle::PrintAndWarnDumper ($codonArray); #die "\n\n\n\nDIE!!!!!!!!!!!$warnMsgHead\n\$right=$right\$codonArray=$codonArray\nCheck it, it is wrong!\n\n\n\n\n\n";
	  return 0;
	}
	else {
		return 1;
	}
}

sub findTheRepresentCDSposARRAY{ #ProSplignHandle::findTheRepresentCDSposARRAY
	my ($inHASH, $number)=@_;               
	
	my $warnMsgBody="\nIn package  ProSplignHandle,\tIn sub findTheRepresentCDSposARRAY,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $caller_inform=DirFileHandle::print_SubCallerInform;
  
  my $outARRAY;
  my @RealCDSPosARRAY;
  if (   (  defined ( $inHASH )  ) && (  ref ( $inHASH ) eq 'HASH' )   ){ print "20181121-1-1\n";
  	foreach my $DNApos (    sort { $a<=>$b} (   keys (  %{ $inHASH }  )   )   ){ print "20181121-1-2\n";  
  		if (   (  defined ( $inHASH->{$DNApos} )  ) && (  ref ( $inHASH->{$DNApos} ) eq 'HASH' ) && (  defined ( $inHASH->{$DNApos}->{'Type'} )  ) && ( $inHASH->{$DNApos}->{'Type'} eq 'NOTaGap')   ){
  			push @RealCDSPosARRAY, $DNApos;                      print "20181121-1-3\n";    
  		}
  	}
  }
  else {
  	my $dieMsg1=$die_MsgHead.$caller_inform."\$inHASH=$inHASH is not right!!! :$!\n\n\n\n";
  	print $dieMsg1;
  	die $dieMsg1;
  }
  if (   (  defined ( $number )  ) && ( $number=~m/\d+/ ) && ( $number>0 )   ){ print "20181121-1-4\n"; 
  	if ($number>=@RealCDSPosARRAY){                                      print "20181121-1-5\n"; 
  	  $outARRAY=[@RealCDSPosARRAY];	
  	}
  	else {  #print "20181121-1-6\n"; 
  		my $realCDSposARRAY_size=@RealCDSPosARRAY;
  		my $FengSan_idx_list=&BuildFengSanNumberList($realCDSposARRAY_size, $number); #print "20181121-1-7\n"; 
  		foreach my $ecFsIdx ( @{$FengSan_idx_list} ){               #print "20181121-1-8 \$RealCDSPosARRAY[$ecFsIdx]=$RealCDSPosARRAY[$ecFsIdx]\n"; 
  			push @{ $outARRAY }, $RealCDSPosARRAY[$ecFsIdx] if (   (  defined ( $RealCDSPosARRAY[$ecFsIdx] )  ) && ( $RealCDSPosARRAY[$ecFsIdx]=~/\d+/ )   );
  		}
  	}
  	
  }
  else {
  	my $dieMsg2=$die_MsgHead.$caller_inform."\$number=$number is not right!!! :$!\n\n\n\n";
  	print $dieMsg2;
  	die $dieMsg2;
  }
	return $outARRAY;
}

sub BuildFengSanNumberList{
	my ($bigNb, $smlNB)=@_;
	my $warnMsgBody="\nIn package  ProSplignHandle,\tIn sub BuildFengSanNumberList,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $caller_inform=DirFileHandle::print_SubCallerInform;
  
  my $FinalList; 
  if (   (  defined ( $bigNb )  ) && ( $bigNb=~m/\d+/ ) && ( $bigNb >0 ) && (  defined ( $smlNB )  ) && ( $smlNB=~m/\d+/ ) && ( $smlNB >0 ) && ( $bigNb >= $smlNB ) && ( $smlNB >= 2 )   ){
    my $PartNB=$smlNB-1;
    my $eaChPart=($bigNb-1)/$PartNB;
    my @tempPointNBList;
    $tempPointNBList[0]=0;
    for (my $i=0; $i<=$PartNB; $i++){
    	$tempPointNBList[$i+1]=$tempPointNBList[$i]+$eaChPart;
    }
    $FinalList->[0]=$tempPointNBList[0];
    for (my $i=1; $i<@tempPointNBList; $i++){
    	$FinalList->[$i]=$tempPointNBList[$i];
    	my $tempSml=int $tempPointNBList[$i];       my $tempBig=$tempSml+1;
    	my $detlSml=$tempPointNBList[$i]-$tempSml;  my $detlBig=$tempBig-$tempPointNBList[$i];
    	if ($tempSml == $FinalList->[$i-1]){
    		$FinalList->[$i]=$tempBig;
    	}
    	else {
    		if ($detlSml<=$detlBig){
    			$FinalList->[$i]=$tempSml;
    		}
    		else {
    			$FinalList->[$i]=$tempBig;
    		}
    	}
    	
    }
    
  }
  else{
  	my $dieMsg=$die_MsgHead.$caller_inform."\$bigNb=$bigNb or \$smlNB=$smlNB is not right!!! :$!\n\n\n\n";
  	print $dieMsg;
  	die $dieMsg;
  }
	return $FinalList;
}

sub PraseProsplignXMLoutfile{
	my $warnMsgHead="\n\n\n   In package  ProSplignHandle,\nIn sub PraseProsplignXMLoutfile,\n\n";
	my ($inXmlFile)=@_;
	my $outArray;

	my $xs1 = XML::Simple->new();

	my $tempMk=$/;
	$/='<?xml version="1.0"?>';
	open (INXML,$inXmlFile) or die "$warnMsgHead\ncannot open \$inXmlFile=$inXmlFile :$!\n";
	while (<INXML>){
	  my $inHere=$_;
	  $inHere=~s/<\?xml version="1.0"\?>$//;
	  $inHere=~s/^<\?xml version="1.0"\?>//;
	  if ($inHere=~m/\S+/){
	    $inHere='<?xml version="1.0"?>'.$inHere;
	    
	    #print "\$inHere=$inHere\n\n\n\n";
	    my $doc = $xs1->XMLin($inHere);

	    push @{ $outArray }, $doc;
	  }

	}
	close (INXML);
	$/=$tempMk;
  
	return $outArray;
	
}

sub RunProSplign{    #ProSplignHandle::RunProSplign
	my ($queryFile, $database, $workIngDir, $genetiCode, $In_GenoOrEst_type, $evalue)=@_;
	my $warnMsgHead="\n\n\n   In package  ProSplignHandle,\nIn sub RunProSplign,\n\n";
	my $proSplign_Path="prosplign";
	my $procompartPath="procompart";
	my $blastallPath="blastall";
	
	my $compFileName='0_0_1_0_bltOut_CompFile.txt';
	my $proSplXmlNme='0_0_1_1_ProSpl_Xml__out.txt';
	my $proSplAliNme='0_0_1_2_ProSpl_Aln__out.txt';
	
	
	my $noIntronFlag=''; if (  ( defined($In_GenoOrEst_type) ) && ($In_GenoOrEst_type=~m/\S+/) && ($In_GenoOrEst_type eq 'Est')  ){  $noIntronFlag=" -no_introns ";  }
	
	
	my $tempTimeDirNeedDel=0; if (  ( defined($workIngDir) ) && ($workIngDir=~m/\S+/)  ) {} else {		$workIngDir=TimeWork::GetTimeDirOrFileName();	$tempTimeDirNeedDel=1; }  #If the $analysisDIR is not set, then build a temporary directory named by time #print "\$workIngDir=$workIngDir\n";
  system ( "mkdir -p $workIngDir ") ;
  
  if ( defined ($genetiCode) ){
  	if ($genetiCode=~m/^\d+$/){  		  	}  	else {my $dieMsg_9= "\n\n\nDIE!!!!!!!!!\n$warnMsgHead\n\$genetiCode=$genetiCode\nIT should be a number!!!\n\n\n\n"; print $dieMsg_9; die $dieMsg_9;}  	
  }else {$genetiCode=1;}
  
  #########
  my $makeblastdbPath="formatdb";                                      #blastplus vesrion:  my $makeblastdbPath="~/EightT/bin/makeblastdb";   	#my $tblastnExecPath="blastall";                                      #blastplus vesrion:  my $tblastnExecPath="~/EightT/bin/tblastn";
	my $mkDatabaseCMD_forDNA="$makeblastdbPath -i $database -p F -o T";  #blastplus vesrion: #my $mkDatabaseCMD_forDNA="$makeblastdbPath -in $database -dbtype nucl -parse_seqids";
	#warn "$warnMsgHead\n\$mkDatabaseCMD_forDNA=$mkDatabaseCMD_forDNA\n"; 
	print "$warnMsgHead\n\$mkDatabaseCMD_forDNA=$mkDatabaseCMD_forDNA\n";        system ("$mkDatabaseCMD_forDNA");    
  #########
  
   
  my $evalueSet=''; if (  ( defined($evalueSet) ) && ($evalueSet=~m/\S+/)  ) {$evalueSet=" -e ".$evalue." ";}
	my $compFile="$workIngDir/$compFileName";
	my $blastToComp_command="$blastallPath -i $queryFile -d $database -p tblastn $evalueSet -D $genetiCode -m 8 | $procompartPath -t > $compFile";
	
	my $proSplignXmlFile    ="$workIngDir/$proSplXmlNme";
	my $proSplignAlimentFile="$workIngDir/$proSplAliNme";  
	my $prosplign_Command="$proSplign_Path -i $compFile -fasta $database,$queryFile -o $proSplignXmlFile -eo $proSplignAlimentFile -genetic_code $genetiCode -full -nogenbank -f x $noIntronFlag ";
	
	print "$warnMsgHead\n$blastToComp_command\n\n";  system ("$blastToComp_command");
	print "$warnMsgHead\n$prosplign_Command\n\n";    system ("$prosplign_Command"  );
	                 
	my $reportArray=&PraseProsplignXMLoutfile($proSplignXmlFile);	
	# my $proSplXmHASH='0_0_1_3_ProSplPrsXmlFil.txt'; my $outArrayFile="$workIngDir/$proSplXmHASH";	DirFileHandle::PrintDumper($outArrayFile,$reportArray)  if ( ref($reportArray) eq 'ARRAY' );
	
	return $reportArray;
	
}

sub Build_Whole_Line_alignHASH{  #ProSplignHandle::Build_Whole_Line_alignHASH
	my ($best_ProSplignOutHash, $GenoOrEst)=@_; 
	
  my $warnMsgBody="\nIn package  ProSplignHandle,\tIn sub Build_Whole_Line_alignHASH,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $dmsg1="$die_MsgHead\$GenoOrEst=$GenoOrEst should be defined , and it should be Est or OneFiledGeno!!\n\n";;
  if (  ( defined ($GenoOrEst) ) && ( ($GenoOrEst eq 'Est') || ($GenoOrEst eq 'OneFiledGeno') )  ){   } else { print $dmsg1; die $dmsg1; }
  my $pExN_0; if ($GenoOrEst eq 'Est') { $pExN_0='2_6e_pExN'; }                elsif ($GenoOrEst eq 'OneFiledGeno') { $pExN_0='0_0g_pExN'; } #Est OneFiledGeno      
  my $aAum_1; if ($GenoOrEst eq 'Est') { $aAum_1='2_5e_aAum'; }                elsif ($GenoOrEst eq 'OneFiledGeno') { $aAum_1='0_1g_aAum'; } #Est OneFiledGeno      
  my $alAA_2; if ($GenoOrEst eq 'Est') { $alAA_2='2_4e_alAA'; }                elsif ($GenoOrEst eq 'OneFiledGeno') { $alAA_2='0_2g_alAA'; } #Est OneFiledGeno      
  my $ppFm_3; if ($GenoOrEst eq 'Est') { $ppFm_3='2_3e_ppFm'; }                elsif ($GenoOrEst eq 'OneFiledGeno') { $ppFm_3='0_3g_ppFm'; } #Est OneFiledGeno      
  my $Csmk_3; if ($GenoOrEst eq 'Est') { $Csmk_3='2_3e_0CsM'; }                elsif ($GenoOrEst eq 'OneFiledGeno') { $Csmk_3='0_3g_zCSm'; } #Est OneFiledGeno      
  my $tAum_4; if ($GenoOrEst eq 'Est') { $tAum_4='2_2e_tAum'; }                elsif ($GenoOrEst eq 'OneFiledGeno') { $tAum_4='0_4g_tAum'; } #Est OneFiledGeno      
  my $trAA_5; if ($GenoOrEst eq 'Est') { $trAA_5='2_1e_trAA'; }                elsif ($GenoOrEst eq 'OneFiledGeno') { $trAA_5='0_5g_trAA'; } #Est OneFiledGeno      
  my $ACTG_6; if ($GenoOrEst eq 'Est') { $ACTG_6='2_0e_ACTG'; }                elsif ($GenoOrEst eq 'OneFiledGeno') { $ACTG_6='0_6g_ACTG'; } #Est OneFiledGeno      

  my $DnPs_a; if ($GenoOrEst eq 'Est') { $DnPs_a='Not_include_e_Dna_posi'; }   elsif ($GenoOrEst eq 'OneFiledGeno') { $DnPs_a='Not_include_g_Dna_posi'; } #Est OneFiledGeno      
  my $PePs_b; if ($GenoOrEst eq 'Est') { $PePs_b='Not_include_e_Pep_posi'; }   elsif ($GenoOrEst eq 'OneFiledGeno') { $PePs_b='Not_include_g_Pep_posi'; } #Est OneFiledGeno      
  my $DnZF_a; if ($GenoOrEst eq 'Est') { $DnZF_a='Not_include_e_Dna_ZorF'; }   elsif ($GenoOrEst eq 'OneFiledGeno') { $DnZF_a='Not_include_g_Dna_ZorF'; } #Est OneFiledGeno      

  
  #0_0g_0_UM    2_0e_ACTG
  #0_0g_alAA    2_1e_trAA
  #0_1g_ppFm    2_1e_umTA
  #0_2g_trAA    2_2e_ppFm
  #0_2g_UMtA    2_3e_alAA
  #0_3g_ACTG    2_3e_uAlA
  #0_4g_pExN
  #Not_include_0Dn_posi
  #Not_include_5Pe_posi
  
  
  
  
  
  if ( ref ($best_ProSplignOutHash) eq 'HASH' ){
  	
  	my $DNA_ZF         =$best_ProSplignOutHash->{'2_DNA_ZF'};
  	my $DNA_to_PEP_hash=$best_ProSplignOutHash->{'7_DNA_to_PEP_hash'};
  	my $PEP_to_DNA_hash=$best_ProSplignOutHash->{'7_PEP_to_DNA_hash'};
  	
  	
  	if (  ( ref ($DNA_to_PEP_hash) eq 'HASH' ) && ( ref ($PEP_to_DNA_hash) eq 'HASH' )  ){
  		my @DNAposKey=keys (  %{ $DNA_to_PEP_hash }  );
  		if    ($DNA_ZF eq '+') {  @DNAposKey=sort { $a <=> $b } @DNAposKey;  }
  		elsif ($DNA_ZF eq '-') {  @DNAposKey=sort { $b <=> $a } @DNAposKey;  }
  		else  {my $dieMsg="$die_MsgHead\n\$DNA_ZF=$DNA_ZF should be + or - !!!!!!!!!\n\n\n\n"; print $dieMsg; die $dieMsg; }
  		
  		my $whoLineHash; my $wLidx=0; 
  		foreach my $DNAposKey ( @DNAposKey ){
  			$whoLineHash->[$wLidx]->{'DnaPos'}=$DNAposKey; 
  			$whoLineHash->[$wLidx]->{'DNAZoF'}=$DNA_to_PEP_hash->{$DNAposKey}->{'DNA_ZoF'  } if (   (  defined ( $DNA_to_PEP_hash->{$DNAposKey}->{'DNA_ZoF'  } )  ) && ( $DNA_to_PEP_hash->{$DNAposKey}->{'DNA_ZoF'  }=~m/\S+/ )   );;  
  			$whoLineHash->[$wLidx]->{'PepPos'}=' '; $whoLineHash->[$wLidx]->{'PepPos'}=$DNA_to_PEP_hash->{$DNAposKey}->{'PepCorPos'} if (   (  defined ( $DNA_to_PEP_hash->{$DNAposKey}->{'PepCorPos'} )  ) && ( $DNA_to_PEP_hash->{$DNAposKey}->{'PepCorPos'}=~m/\d+/ )   );
  			$whoLineHash->[$wLidx]->{'PepFrm'}=$DNA_to_PEP_hash->{$DNAposKey}->{'PepFrm'   } if (   (  defined ( $DNA_to_PEP_hash->{$DNAposKey}->{'PepFrm'   } )  ) && ( $DNA_to_PEP_hash->{$DNAposKey}->{'PepFrm'   }=~m/\d+/ )   );
  			$whoLineHash->[$wLidx]->{'type_1'}='DNA'; 
  			$whoLineHash->[$wLidx]->{'DnaHas'}=Storable::dclone ( $DNA_to_PEP_hash->{$DNAposKey} );
  			$wLidx++;
  		}
  		
  		my $finalWhoLineHash; my $fWlIdx=0;
  		
  		my @SortedPepPosKeyArray=sort { $a <=> $b } (   keys (  %{ $PEP_to_DNA_hash }  )   );
  		for (  my $i=0; $i<@{ $whoLineHash }; $i++  ){
  			
  		  if (   (  defined ( $whoLineHash->[$i]->{'PepPos'} )  ) && ( $whoLineHash->[$i]->{'PepPos'}=~m/\d+/ )   ){
  		  	my $PepPosHere=$whoLineHash->[$i]->{'PepPos'};
  		  	
  		  	while (   (  defined ( $SortedPepPosKeyArray[0] )  ) && ( $SortedPepPosKeyArray[0] < $PepPosHere )   ){
  		  		my $PepPosKey=shift @SortedPepPosKeyArray;
  		  		if (   (  defined ( $PEP_to_DNA_hash->{$PepPosKey}->{'real_1_3_Pos'} )  ) && (  ref ( $PEP_to_DNA_hash->{$PepPosKey}->{'real_1_3_Pos'} ) eq 'HASH'  )   ){
  		  			foreach my $fmNBk (    sort { $a<=>$b } (   keys (  %{ $PEP_to_DNA_hash->{$PepPosKey}->{'real_1_3_Pos'} }  )   )    ){
  		  				if (   (  defined ( $PEP_to_DNA_hash->{$PepPosKey}->{'real_1_3_Pos'}->{$fmNBk} )  ) && ( $PEP_to_DNA_hash->{$PepPosKey}->{'real_1_3_Pos'}->{$fmNBk}=~m/\d+/ )   ){
  		  				  my $dieMsg="$die_MsgHead\n111111 \$PEP_to_DNA_hash->{$PepPosKey}->{'real_1_3_Pos'}->{$fmNBk}=$PEP_to_DNA_hash->{$PepPosKey}->{'real_1_3_Pos'}->{$fmNBk} should not be defined or as a number !!!\n\n\n";
  		  				  print  $dieMsg;
  		  				  die    $dieMsg;
  		  			  }
  		  			  else {
  		  			  	$finalWhoLineHash->[$fWlIdx]->{'PepPos'}=$PepPosKey;
  		  				  $finalWhoLineHash->[$fWlIdx]->{'PepFrm'}=$fmNBk; 
  		  				  $finalWhoLineHash->[$fWlIdx]->{'type_2'}='PEP'; 
  		  				  $finalWhoLineHash->[$fWlIdx]->{'PepHas'}=Storable::dclone ( $PEP_to_DNA_hash->{$PepPosKey} ); 
  		  		      $fWlIdx++;
  		  			  }
  		  			}
  		  		}
  		  		else {  		  			
  		  			for (my $k=1; $k<=3; $k++){ 
  		  				$finalWhoLineHash->[$fWlIdx]->{'PepPos'}=$PepPosKey;
  		  				$finalWhoLineHash->[$fWlIdx]->{'PepFrm'}=$k; 
  		  				$finalWhoLineHash->[$fWlIdx]->{'type_3'}='PEP'; 
  		  				$finalWhoLineHash->[$fWlIdx]->{'PepHas'}=Storable::dclone ( $PEP_to_DNA_hash->{$PepPosKey} ); 
  		  		    $fWlIdx++;
  		  			}
  		  		}
  		  	}
  		  	if (   (  defined ( $SortedPepPosKeyArray[0] )  ) && ( $SortedPepPosKeyArray[0] == $PepPosHere )   ){
  		  		my $PepPosKey=shift @SortedPepPosKeyArray;
  		  		
  		  		if (   (  defined ( $PEP_to_DNA_hash->{$PepPosKey}->{'real_1_3_Pos'} )  ) && (  ref ( $PEP_to_DNA_hash->{$PepPosKey}->{'real_1_3_Pos'} ) eq 'HASH'  )   ){
  		  			my $howManyKey=(   keys (  %{ $PEP_to_DNA_hash->{$PepPosKey}->{'real_1_3_Pos'} }  )   ); my $inNB=1;
  		  		  foreach my $fmNBk (    sort { $a<=>$b } (   keys (  %{ $PEP_to_DNA_hash->{$PepPosKey}->{'real_1_3_Pos'} }  )   )    ){
  		  		  	if (   (  defined ( $PEP_to_DNA_hash->{$PepPosKey}->{'real_1_3_Pos'}->{$fmNBk} )  ) && ( $PEP_to_DNA_hash->{$PepPosKey}->{'real_1_3_Pos'}->{$fmNBk}=~m/\d+/ )   ){
  		  		  	  if ( $whoLineHash->[$i]->{'DnaPos'} eq $PEP_to_DNA_hash->{$PepPosKey}->{'real_1_3_Pos'}->{$fmNBk} ){
  		  		  	  	$finalWhoLineHash->[$fWlIdx]=Storable::dclone ( $whoLineHash->[$i] ); 
  		  		  	  	$finalWhoLineHash->[$fWlIdx]->{'DnaPos'}=$whoLineHash->[$i]->{'DnaPos'};  	if ($inNB <$howManyKey) { $i++; } # The $i should be ++ in the inside 循环	  			   
  		  		  	  	$finalWhoLineHash->[$fWlIdx]->{'DNAZoF'}=$whoLineHash->[$i]->{'DNAZoF'};  	  			   
  		  		  	  	$finalWhoLineHash->[$fWlIdx]->{'PepPos'}=$PepPosKey;
  		  		  	    $finalWhoLineHash->[$fWlIdx]->{'PepFrm'}=$fmNBk; 
  		  		  	    $finalWhoLineHash->[$fWlIdx]->{'type_4'}='PEP'; 
  		  		  	    $finalWhoLineHash->[$fWlIdx]->{'PepHas'}=Storable::dclone ( $PEP_to_DNA_hash->{$PepPosKey} ); 
  		  		        $fWlIdx++;
  		  		  	  }
  		  		  	  else{
  		  		  	    my $dieMsg="$die_MsgHead\n22222222 \$PEP_to_DNA_hash->{$PepPosKey}->{'real_1_3_Pos'}->{$fmNBk}=$PEP_to_DNA_hash->{$PepPosKey}->{'real_1_3_Pos'}->{$fmNBk} \$whoLineHash->[$i]->{'DnaPos'}=$whoLineHash->[$i]->{'DnaPos'} These 2 should  be the same !!!\n\n\n";
  		  		  	    print  $dieMsg;
  		  		  	    die    $dieMsg;
  		  		  	  }  		  			  
  		  		  	
  		  		    }
  		  		    else {
  		  		  	  my $dieMsg="$die_MsgHead\n333333 \$PEP_to_DNA_hash->{$PepPosKey}->{'real_1_3_Pos'}->{$fmNBk}=$PEP_to_DNA_hash->{$PepPosKey}->{'real_1_3_Pos'}->{$fmNBk} should  be a number !!!\n\n\n";
  		  		  	  print  $dieMsg;
  		  		  	  die    $dieMsg;
  		  		    }
  		  		    $inNB++;
  		  		  }	
  		  		}	  		  		  		
  		  	}
  		  }
  		  else{
  		  	$finalWhoLineHash->[$fWlIdx]=Storable::dclone ( $whoLineHash->[$i] ); 
  		  	$fWlIdx++;
  		  }
  	  }
  	  foreach my $pepPsK ( @SortedPepPosKeyArray ){
  	  	
  	  	for (my $k=1; $k<=3; $k++){ 
  		  	if (   (  defined ( $PEP_to_DNA_hash->{$pepPsK}->{'real_1_3_Pos'}->{$k} )  ) && ( $PEP_to_DNA_hash->{$pepPsK}->{'real_1_3_Pos'}->{$k}=~m/\d+/ )   ){
  		  		my $dieMsg="$die_MsgHead\n222222 \$PEP_to_DNA_hash->{$pepPsK}->{'real_1_3_Pos'}->{$k}=$PEP_to_DNA_hash->{$pepPsK}->{'real_1_3_Pos'}->{$k} should not be defined or as a number !!!\n\n\n";
  		  		print  $dieMsg;
  		  		die    $dieMsg;
  		  	}
  		  	else {
  		  		$finalWhoLineHash->[$fWlIdx]->{'PepPos'}=$pepPsK;
  		  		$finalWhoLineHash->[$fWlIdx]->{'PepFrm'}=$k; 
  		  		$finalWhoLineHash->[$fWlIdx]->{'type_5'}='PEP'; 
  		  		$finalWhoLineHash->[$fWlIdx]->{'PepHas'}=Storable::dclone ( $PEP_to_DNA_hash->{$pepPsK} ); 
  		  	  $fWlIdx++;
  		  	}
  		  }
  	  	
  	  	
  	  }    
  	  
  	  my $ptOutWholineHash;
  	  
  	  $ptOutWholineHash->{'4_PspGnPHs'}=$finalWhoLineHash;
  	  
  	  for (  my $i=0; $i<@{ $finalWhoLineHash }; $i++  ){
  	  	
  	  	 
  	  	$ptOutWholineHash->{'2_alnShwHs'}->{$i}->{$alAA_2}=' ';  $ptOutWholineHash->{'2_alnShwHs'}->{$i}->{$alAA_2}=$finalWhoLineHash->[$i]->{'PepHas'}->{'PePChar'} if (   (  defined ( $finalWhoLineHash->[$i]->{'PepHas'}->{'PePChar'} )  ) && ( $finalWhoLineHash->[$i]->{'PepHas'}->{'PePChar'}=~m/\S+/ )  && (  defined ( $finalWhoLineHash->[$i]->{'PepFrm'} )  ) && ( $finalWhoLineHash->[$i]->{'PepFrm'}==2 )  );
  	  	$ptOutWholineHash->{'2_alnShwHs'}->{$i}->{$aAum_1}=' ';  $ptOutWholineHash->{'2_alnShwHs'}->{$i}->{$aAum_1}=$ptOutWholineHash->{'2_alnShwHs'}->{$i}->{$alAA_2} if ($ptOutWholineHash->{'2_alnShwHs'}->{$i}->{$alAA_2} eq 'U');
  	  	
  	  	$ptOutWholineHash->{'2_alnShwHs'}->{$i}->{$ppFm_3}=' ';  $ptOutWholineHash->{'2_alnShwHs'}->{$i}->{$ppFm_3}=$finalWhoLineHash->[$i]->{'PepFrm'}              if (   (  defined ( $finalWhoLineHash->[$i]->{'PepFrm'}              )  ) && ( $finalWhoLineHash->[$i]->{'PepFrm'}=~m/\d+/ )  );
  	  	
  	  	
  	  	$ptOutWholineHash->{'2_alnShwHs'}->{$i}->{$trAA_5}=' ';  $ptOutWholineHash->{'2_alnShwHs'}->{$i}->{$trAA_5}=$finalWhoLineHash->[$i]->{'DnaHas'}->{'trasAA'}  if (   (  defined ( $finalWhoLineHash->[$i]->{'DnaHas'}->{'trasAA' } )  ) && ( $finalWhoLineHash->[$i]->{'DnaHas'}->{'trasAA' }=~m/\S+/ )   && (  defined ( $finalWhoLineHash->[$i]->{'PepFrm'} )  ) && ( $finalWhoLineHash->[$i]->{'PepFrm'}==2 )  );
        $ptOutWholineHash->{'2_alnShwHs'}->{$i}->{$tAum_4}=' ';  $ptOutWholineHash->{'2_alnShwHs'}->{$i}->{$tAum_4}=$ptOutWholineHash->{'2_alnShwHs'}->{$i}->{$trAA_5} if ( ($ptOutWholineHash->{'2_alnShwHs'}->{$i}->{$trAA_5} eq 'U') || ($ptOutWholineHash->{'2_alnShwHs'}->{$i}->{$trAA_5} eq 'M') );
  	  	
  	  	$ptOutWholineHash->{'2_alnShwHs'}->{$i}->{$Csmk_3}=' ';  $ptOutWholineHash->{'2_alnShwHs'}->{$i}->{$Csmk_3}='|' if (   (  defined ( $ptOutWholineHash->{'2_alnShwHs'}->{$i}->{$trAA_5} )  ) && ( $ptOutWholineHash->{'2_alnShwHs'}->{$i}->{$trAA_5}=~m/\S+/ )  && (  defined ( $ptOutWholineHash->{'2_alnShwHs'}->{$i}->{$alAA_2} )  ) && ( $ptOutWholineHash->{'2_alnShwHs'}->{$i}->{$alAA_2}=~m/\S+/ ) && ( $ptOutWholineHash->{'2_alnShwHs'}->{$i}->{$trAA_5} eq $ptOutWholineHash->{'2_alnShwHs'}->{$i}->{$alAA_2} )    );
  	  	
  	  	$ptOutWholineHash->{'2_alnShwHs'}->{$i}->{$ACTG_6}=' ';  $ptOutWholineHash->{'2_alnShwHs'}->{$i}->{$ACTG_6}=$finalWhoLineHash->[$i]->{'DnaHas'}->{'DNAChar'} if (   (  defined ( $finalWhoLineHash->[$i]->{'DnaHas'}->{'DNAChar'} )  ) && ( $finalWhoLineHash->[$i]->{'DnaHas'}->{'DNAChar'}=~m/\S+/ )   );
  	  	$ptOutWholineHash->{'2_alnShwHs'}->{$i}->{$ACTG_6}=uc ( $ptOutWholineHash->{'2_alnShwHs'}->{$i}->{$ACTG_6} ) if (   (  defined ( $finalWhoLineHash->[$i]->{'DnaHas'}->{'PepFrm'} )  ) && ( $finalWhoLineHash->[$i]->{'DnaHas'}->{'PepFrm'}=~m/\d+/ ) && ( $finalWhoLineHash->[$i]->{'DnaHas'}->{'PepFrm'}==1 )   );
  	  	$ptOutWholineHash->{'2_alnShwHs'}->{$i}->{$ACTG_6}=lc ( $ptOutWholineHash->{'2_alnShwHs'}->{$i}->{$ACTG_6} ) if (   (  defined ( $finalWhoLineHash->[$i]->{'DnaHas'}->{'PepFrm'} )  ) && ( $finalWhoLineHash->[$i]->{'DnaHas'}->{'PepFrm'}=~m/\d+/ ) && ( $finalWhoLineHash->[$i]->{'DnaHas'}->{'PepFrm'}==2 )   );
  	  	$ptOutWholineHash->{'2_alnShwHs'}->{$i}->{$ACTG_6}=lc ( $ptOutWholineHash->{'2_alnShwHs'}->{$i}->{$ACTG_6} ) if (   (  defined ( $finalWhoLineHash->[$i]->{'DnaHas'}->{'PepFrm'} )  ) && ( $finalWhoLineHash->[$i]->{'DnaHas'}->{'PepFrm'}=~m/\d+/ ) && ( $finalWhoLineHash->[$i]->{'DnaHas'}->{'PepFrm'}==3 )   );
  	  	
  	  	$ptOutWholineHash->{'2_alnShwHs'}->{$i}->{$pExN_0}=' ';  $ptOutWholineHash->{'2_alnShwHs'}->{$i}->{$pExN_0}=$finalWhoLineHash->[$i]->{'DnaHas'}->{'6Ex_numb'} if (   (  defined ( $finalWhoLineHash->[$i]->{'DnaHas'}->{'6Ex_numb'} )  ) && ( $finalWhoLineHash->[$i]->{'DnaHas'}->{'6Ex_numb'}=~m/\d+/ )   );
  	  	
  	  	
  	  	$ptOutWholineHash->{'2_alnShwHs'}->{$i}->{$DnZF_a}=$finalWhoLineHash->[$i]->{'DNAZoF'}  if (   (  defined ( $finalWhoLineHash->[$i]->{'DNAZoF' } )  ) && ( $finalWhoLineHash->[$i]->{'DNAZoF'}=~m/\S+/ )   );
  	  	$ptOutWholineHash->{'2_alnShwHs'}->{$i}->{$DnPs_a}=$finalWhoLineHash->[$i]->{'DnaPos'}  if (   (  defined ( $finalWhoLineHash->[$i]->{'DnaPos' } )  ) && ( $finalWhoLineHash->[$i]->{'DnaPos'}=~m/\d+/ )   );
  	  	$ptOutWholineHash->{'2_alnShwHs'}->{$i}->{$PePs_b}=$finalWhoLineHash->[$i]->{'PepPos'}  if (   (  defined ( $finalWhoLineHash->[$i]->{'PepPos' } )  ) && ( $finalWhoLineHash->[$i]->{'PepPos'}=~m/\d+/ )   );
  	  	
  	  	$ptOutWholineHash->{'0___alnShw'}->{$aAum_1} .= $ptOutWholineHash->{'2_alnShwHs'}->{$i}->{$aAum_1};
  	  	$ptOutWholineHash->{'0___alnShw'}->{$alAA_2} .= $ptOutWholineHash->{'2_alnShwHs'}->{$i}->{$alAA_2};
  	  	$ptOutWholineHash->{'0___alnShw'}->{$ppFm_3} .= $ptOutWholineHash->{'2_alnShwHs'}->{$i}->{$ppFm_3};
  	  	$ptOutWholineHash->{'0___alnShw'}->{$Csmk_3} .= $ptOutWholineHash->{'2_alnShwHs'}->{$i}->{$Csmk_3};
  	  	$ptOutWholineHash->{'0___alnShw'}->{$trAA_5} .= $ptOutWholineHash->{'2_alnShwHs'}->{$i}->{$trAA_5}; 	  	                                                      
  	  	$ptOutWholineHash->{'0___alnShw'}->{$tAum_4} .= $ptOutWholineHash->{'2_alnShwHs'}->{$i}->{$tAum_4}; 	  	                                                      
  	  	$ptOutWholineHash->{'0___alnShw'}->{$ACTG_6} .= $ptOutWholineHash->{'2_alnShwHs'}->{$i}->{$ACTG_6};
  	  	$ptOutWholineHash->{'0___alnShw'}->{$pExN_0} .= $ptOutWholineHash->{'2_alnShwHs'}->{$i}->{$pExN_0}
  	  	
  	  	###################

  	  	###################
  	  	
  	  	  	  	
  	  }
  	  
  	  
  	  
  	  
  	  
  	  
  	  return $ptOutWholineHash;
  	  
  		
  	}
  	
  	
  	
  	
  }
	
}




sub Compair_the_Prosplign_Wln_Hash{  #  ProSplignHandle::Compair_the_Prosplign_Wln_Hash
	my ($inProsWlnHASH_1, $inProsWlnHASH_2, $Geno_or_EST, $file_1, $file_2)=@_;
	
	my $warnMsgBody="\nIn package  ProSplignHandle,\tIn sub Compair_the_Prosplign_Wln_Hash,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	print "201811191434-0\n";
	my $inPros_ZF_PS_Pf_HASH_1=ProSplignHandle::Get_ZorF_Posi_Pfrm_HASH_from_ProSWln_HASH($inProsWlnHASH_1, $Geno_or_EST, $file_1); print "201811191434-0-1 \$file_1=$file_1 \$inPros_ZF_PS_Pf_HASH_1=$inPros_ZF_PS_Pf_HASH_1\n";
	my $inPros_ZF_PS_Pf_HASH_2=ProSplignHandle::Get_ZorF_Posi_Pfrm_HASH_from_ProSWln_HASH($inProsWlnHASH_2, $Geno_or_EST, $file_2); print "201811191434-0-2 \$file_2=$file_2 \$inPros_ZF_PS_Pf_HASH_2=$inPros_ZF_PS_Pf_HASH_2\n";
	
	my $TheSameDNAposNUMBER=0;  my $The_same_rate_1=0;  my $The_same_rate_2=0;
	if (      (  defined ( $inPros_ZF_PS_Pf_HASH_1->{'0_1_realHASH'} )  ) && (  ref( $inPros_ZF_PS_Pf_HASH_1->{'0_1_realHASH'} ) eq 'HASH'  ) 
	       && (  defined ( $inPros_ZF_PS_Pf_HASH_2->{'0_1_realHASH'} )  ) && (  ref( $inPros_ZF_PS_Pf_HASH_2->{'0_1_realHASH'} ) eq 'HASH'  )   
	   )
	{ print "201811191434-1\n";
    foreach my $ZoF (    sort {$a cmp $b}(   keys(  %{ $inPros_ZF_PS_Pf_HASH_1->{'0_1_realHASH'} }  )   )    ){ print "201811191434-2\n";
    	if (   (  defined ( $inPros_ZF_PS_Pf_HASH_1->{'0_1_realHASH'}->{$ZoF} )  ) && (  ref( $inPros_ZF_PS_Pf_HASH_1->{'0_1_realHASH'}->{$ZoF} ) eq 'HASH'  )   ){ print "201811191434-3\n";
        foreach my $Pos (    sort {$a <=> $b}(   keys(  %{ $inPros_ZF_PS_Pf_HASH_1->{'0_1_realHASH'}->{$ZoF} }  )   )    ){ print "201811191434-4\n";
          if (   (  defined ( $inPros_ZF_PS_Pf_HASH_1->{'0_1_realHASH'}->{$ZoF}->{$Pos} )  ) && ( $inPros_ZF_PS_Pf_HASH_1->{'0_1_realHASH'}->{$ZoF}->{$Pos}=~m/\d+/ )   ){  print "201811191434-5\n";
          	
          	if (   (  defined ( $inPros_ZF_PS_Pf_HASH_2->{'0_1_realHASH'}->{$ZoF} )  ) &&(  defined ( $inPros_ZF_PS_Pf_HASH_2->{'0_1_realHASH'}->{$ZoF}->{$Pos} )  ) && ( $inPros_ZF_PS_Pf_HASH_2->{'0_1_realHASH'}->{$ZoF}->{$Pos}=~m/\d+/ )   ){
          	  print "201811191434-6\n";
          	  if ( $inPros_ZF_PS_Pf_HASH_1->{'0_1_realHASH'}->{$ZoF}->{$Pos} == $inPros_ZF_PS_Pf_HASH_2->{'0_1_realHASH'}->{$ZoF}->{$Pos} ){
          	  	$TheSameDNAposNUMBER++; print "201811191434-7 \$TheSameDNAposNUMBER=$TheSameDNAposNUMBER\n";
          	  }
          	}
          	
          }         
        }
      }   
    }
    
    if (   (  defined ( $inPros_ZF_PS_Pf_HASH_1->{'0_2_totlPsNB'} )  ) && ( $inPros_ZF_PS_Pf_HASH_1->{'0_2_totlPsNB'}=~m/\d+/ ) && (  defined ( $inPros_ZF_PS_Pf_HASH_2->{'0_2_totlPsNB'} )  ) && ( $inPros_ZF_PS_Pf_HASH_2->{'0_2_totlPsNB'}=~m/\d+/ )   ){
    	$The_same_rate_1=$TheSameDNAposNUMBER/$inPros_ZF_PS_Pf_HASH_1->{'0_2_totlPsNB'} if ($inPros_ZF_PS_Pf_HASH_1->{'0_2_totlPsNB'} > 0); print "201811191434-8 \$The_same_rate_1=$TheSameDNAposNUMBER/$inPros_ZF_PS_Pf_HASH_1->{'0_2_totlPsNB'}=$The_same_rate_1\n";
    	$The_same_rate_2=$TheSameDNAposNUMBER/$inPros_ZF_PS_Pf_HASH_2->{'0_2_totlPsNB'} if ($inPros_ZF_PS_Pf_HASH_2->{'0_2_totlPsNB'} > 0); print "201811191434-9 \$The_same_rate_2=$TheSameDNAposNUMBER/$inPros_ZF_PS_Pf_HASH_2->{'0_2_totlPsNB'}=$The_same_rate_2\n";
    }
    
  }
	else {
  	my $dieMsg=$die_MsgHead.$caller_inform."\$inPros_ZF_PS_Pf_HASH_1=$inPros_ZF_PS_Pf_HASH_1 from \$file_1=$file_1 \nor\n\$inPros_ZF_PS_Pf_HASH_2=$inPros_ZF_PS_Pf_HASH_2 from \$file_2=$file_2 is not HASH ref!!! $!\n\n\n";
  	print $dieMsg;
  	die $dieMsg;
  }
  
  return [ $The_same_rate_1, $The_same_rate_2 ];
  
}

sub Get_ZorF_Posi_Pfrm_HASH_from_ProSWln_HASH{  #ProSplignHandle::Get_ZorF_Posi_Pfrm_HASH_from_ProSWln_HASH
	my ($inProsWlnHASH, $Gno_ro_Est, $file)=@_;
	
	my $warnMsgBody="\nIn package  ProSplignHandle,\tIn sub Get_ZorF_Posi_Pfrm_HASH_from_ProSWln_HASH,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $caller_inform=DirFileHandle::print_SubCallerInform;
  
  $Gno_ro_Est=uc ($Gno_ro_Est);
  
  my $DNA_ZorF_key='Not_include_g_Dna_ZorF'; if ($Gno_ro_Est eq 'EST'){ $DNA_ZorF_key='Not_include_e_Dna_ZorF'; }
  my $DNA_posi_key='Not_include_g_Dna_posi'; if ($Gno_ro_Est eq 'EST'){ $DNA_posi_key='Not_include_e_Dna_posi'; }
  my $DNA_ppFm_key='0_3g_ppFm'             ; if ($Gno_ro_Est eq 'EST'){ $DNA_ppFm_key='2_3e_ppFm'             ; }
  
  
  my $out_ZorF_Posi_Pfrm_HASH;  print "201811191401-0 $Gno_ro_Est $file\n";
  if (   (  defined ( $inProsWlnHASH )  ) && (  ref( $inProsWlnHASH ) eq 'HASH'  ) && (  defined ( $inProsWlnHASH->{'2_alnShwHs'} )  ) && (  ref( $inProsWlnHASH->{'2_alnShwHs'} ) eq 'HASH'  )   ){ print "201811191401-1 $file\n";
    foreach my $AlnKEY (    sort {$a <=> $b}(   keys(  %{ $inProsWlnHASH->{'2_alnShwHs'} }  )   )    ){  print "201811191401-2 $Gno_ro_Est $file \$AlnKEY=$AlnKEY\n";
      my $alnHS=$inProsWlnHASH->{'2_alnShwHs'}->{$AlnKEY};
      if (      (  defined ( $alnHS->{$DNA_ZorF_key} )  ) && ( $alnHS->{$DNA_ZorF_key}=~m/\S+/ ) 
             && (  defined ( $alnHS->{$DNA_posi_key} )  ) && ( $alnHS->{$DNA_posi_key}=~m/\d+/ ) 
             && (  defined ( $alnHS->{$DNA_ppFm_key} )  ) && ( $alnHS->{$DNA_ppFm_key}=~m/\d+/ )   
         )
      { 
        $out_ZorF_Posi_Pfrm_HASH->{'0_1_realHASH'}->{ $alnHS->{$DNA_ZorF_key} }->{ $alnHS->{$DNA_posi_key} }=$alnHS->{$DNA_ppFm_key}; print "201811191401-3 $file \$out_ZorF_Posi_Pfrm_HASH->{'0_1_realHASH'}->{ $alnHS->{$DNA_ZorF_key} }->{ $alnHS->{$DNA_posi_key} }=$out_ZorF_Posi_Pfrm_HASH->{'0_1_realHASH'}->{ $alnHS->{$DNA_ZorF_key} }->{ $alnHS->{$DNA_posi_key} }=\$alnHS->{$DNA_ppFm_key}\n";
      }
    }
    
    if (      (  defined ( $out_ZorF_Posi_Pfrm_HASH )  ) && (  ref( $out_ZorF_Posi_Pfrm_HASH ) eq 'HASH'  ) 
           && (  defined ( $out_ZorF_Posi_Pfrm_HASH->{'0_1_realHASH'} )  ) && (  ref( $out_ZorF_Posi_Pfrm_HASH->{'0_1_realHASH'} ) eq 'HASH'  )   
       )
    { print "201811191401-4 $file \n";
    	$out_ZorF_Posi_Pfrm_HASH->{'0_2_totlPsNB'}=0;
    	foreach my $ZoF (    sort {$a cmp $b}(   keys(  %{ $out_ZorF_Posi_Pfrm_HASH->{'0_1_realHASH'} }  )   )    ){   print "201811191401-5 $file \n";
    	  if (   (  defined ( $out_ZorF_Posi_Pfrm_HASH->{'0_1_realHASH'}->{$ZoF} )  ) && (  ref( $out_ZorF_Posi_Pfrm_HASH->{'0_1_realHASH'}->{$ZoF} ) eq 'HASH'  )   ){   print "201811191401-6 $file \n";
          foreach my $Pos (    sort {$a <=> $b}(   keys(  %{ $out_ZorF_Posi_Pfrm_HASH->{'0_1_realHASH'}->{$ZoF} }  )   )    ){    
          	$out_ZorF_Posi_Pfrm_HASH->{'0_2_totlPsNB'}++;  print "201811191401-7 $file \$out_ZorF_Posi_Pfrm_HASH->{'0_2_totlPsNB'}=$out_ZorF_Posi_Pfrm_HASH->{'0_2_totlPsNB'}\n";
          }
        }
    	}
    }
    
  }
  else {
  	my $dieMsg=$die_MsgHead.$caller_inform."\$inProsWlnHASH=$inProsWlnHASH or \$inProsWlnHASH->{'2_alnShwHs'}=$inProsWlnHASH->{'2_alnShwHs'} is not HASH ref!!! $!\n\n\n";
  	print $dieMsg;
  	die $dieMsg;
  }
  
  return $out_ZorF_Posi_Pfrm_HASH;
  
	
}


sub build_SECIS_hash{  #ProSplignHandle::build_SECIS_hash ( $AllSECISinformationHash, $ecSelPro, $speDIRNM, $GenoOrEst, $contigNM, $strand_Z_or_F, $tgaPos, $lastEndPos )
	
	
	my ( $AllSECISinformationHash, $ecSelPro, $speDIRNM, $GenoOrEst, $contigNM, $strand_Z_or_F, $tgaPos, $lastEndPos )=@_;
	
	my $warnMsgBody="\nIn package  ProSplignHandle,\tIn sub build_SECIS_hash,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	my $SECIS_out_HASH;
	
	#my $bestSECISsgmARRAY;    	        my $all_SECISsgmARRAY;    	        my $Mbk_bestSECISsgmARRAY;    	        my $Mbk_all_SECISsgmARRAY;
  ###########SECIS information
  
  if (   defined (  $AllSECISinformationHash->{$speDIRNM}->{$GenoOrEst}->{$contigNM}->{$tgaPos}  )   ){
  	my $SECISarray=$AllSECISinformationHash->{$speDIRNM}->{$GenoOrEst}->{$contigNM}->{$tgaPos} ;  	#DirFileHandle::PrintDumper($File_00_14_SECIS_____all, $SECISarray) if ( ref($SECISarray) eq 'ARRAY' ) ;
    $SECIS_out_HASH->{'0_0_0_allSECIS_ifm_array'}=$SECISarray;
    if ( ref ($SECISarray) eq 'ARRAY' ){
    	my $all_SECISsgmARRAY=SECISwork::GetAllSECISsgmentARRAY($SECISarray);     #DirFileHandle::PrintDumper($File_00_17_SECIS_AllSgmt, $all_SECISsgmARRAY) if ( ref($all_SECISsgmARRAY) eq 'ARRAY' ) ;   
  		$SECIS_out_HASH->{'0_0_1_allSECIS_sgm_ARRAY'}=$all_SECISsgmARRAY;
      
      #my $lastEndPos; 
  		#if (  ( defined ($wholGeneEnd__Pos    ) ) && ($wholGeneEnd__Pos=~m/\d+/    )  ){$lastEndPos=$wholGeneEnd__Pos    ;}  #print "2$ecSelPro 1\$lastEndPos=$lastEndPos\n";
  		#if (  ( defined ($Exp_wholGeneEnd__Pos) ) && ($Exp_wholGeneEnd__Pos=~m/\d+/)  ){$lastEndPos=$Exp_wholGeneEnd__Pos;}  #print "2$ecSelPro 2\$lastEndPos=$lastEndPos\n";
  		
      my $sorted_SECISArray=SECISwork::SortSECISarray        ($tgaPos, $strand_Z_or_F, $SECISarray, $lastEndPos, $ecSelPro);  #DirFileHandle::PrintDumper($File_00_15_SECIS_allSort, $sorted_SECISArray) if ( ref($sorted_SECISArray) eq 'ARRAY' ) ;
  		if (   ( defined ( $sorted_SECISArray )  ) && (  ref ( $sorted_SECISArray ) eq 'ARRAY'  )   ){
  		  $SECIS_out_HASH->{'0_1_0_sorted_SECIS_ARRAY'}=$sorted_SECISArray;                                                                               #my $sorted_SECISArrayDump=DirFileHandle::ReturnDumperInform($sorted_SECISArray);                       #print "a2-20180906$ecSelPro $GenoOrEst\$sorted_SECISArray=$sorted_SECISArray\n\$sorted_SECISArrayDump=$sorted_SECISArrayDump\n";
  		}
  		my $bestSECISarray   =SECISwork::GetBestSECISarray     ($tgaPos, $strand_Z_or_F, $SECISarray, $lastEndPos, $ecSelPro);  #DirFileHandle::PrintDumper($File_00_16_SECIS_SortBst,    $bestSECISarray   ) if ( ref($bestSECISarray   ) eq 'ARRAY' ) ;
  		if (   ( defined ( $bestSECISarray )  ) && (  ref ( $bestSECISarray ) eq 'ARRAY'  )   ){
  		  $SECIS_out_HASH->{'0_1_1_best___SECIS_ARRAY'}=$bestSECISarray;    		                                                #my $bestSECISarrayDump=DirFileHandle::ReturnDumperInform($bestSECISarray);                             #print "a2-20180906$ecSelPro $GenoOrEst\$bestSECISarray=$bestSECISarray\n\$bestSECISarrayDump=$bestSECISarrayDump\n";
  		  my $bestSECISsgmARRAY=SECISwork::GetAllSECISsgmentARRAY($bestSECISarray     );                                        #DirFileHandle::PrintDumper($File_00_18_SECIS_BstSgmt, $bestSECISsgmARRAY) if ( ref($bestSECISsgmARRAY) eq 'ARRAY' ) ;
  		  if (   ( defined ( $bestSECISsgmARRAY )  ) && (  ref ( $bestSECISsgmARRAY ) eq 'ARRAY'  )   ){
  		    $SECIS_out_HASH->{'0_1_1_best_SIS_sgm_ARRAY'}=$bestSECISsgmARRAY;
  	    }
  		}
  		
  	}
  }
  return $SECIS_out_HASH;
}
      
sub Map_Est_SECIS_backto_Genome_based_on_SplignOUT{  #ProSplignHandle::Map_Est_SECIS_backto_Genome_based_on_SplignOUT ($SECIS_HASH, $SplignOUTarray)
	my ($SECIS_HASH, $SplignOUTarray)=@_;
	
	my $warnMsgBody="\nIn package  ProSplignHandle,\tIn sub Map_Est_SECIS_backto_Genome_based_on_SplignOUT,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $caller_inform=DirFileHandle::print_SubCallerInform;
  
  my $SECIS_out_HASH=$SECIS_HASH;
  if (   ( defined ( $SECIS_HASH )  ) && (  ref ( $SECIS_HASH ) eq 'HASH'  ) && ( defined ( $SplignOUTarray )  ) && (  ref ( $SplignOUTarray ) eq 'ARRAY'  )   ){
  	if (   ( defined ( $SECIS_HASH->{'0_0_1_allSECIS_sgm_ARRAY'} )  ) && (  ref ( $SECIS_HASH->{'0_0_1_allSECIS_sgm_ARRAY'} ) eq 'ARRAY'  )   ){
  	  my $all_SECISsgmARRAY=$SECIS_HASH->{'0_0_1_allSECIS_sgm_ARRAY'};
  	  my $Mbk_all_SECISsgmARRAY=SplignHandle::Map_est_sgmt_2_div_array_bk_to_geno($SplignOUTarray, $all_SECISsgmARRAY) ;   	  #DirFileHandle::PrintDumper($File_00_19_SIS_allSgmMbk, $Mbk_all_SECISsgmARRAY) if ( ref($Mbk_all_SECISsgmARRAY) eq 'ARRAY' ) ;  
  	  if (   ( defined ( $Mbk_all_SECISsgmARRAY )  ) && (  ref ( $Mbk_all_SECISsgmARRAY ) eq 'ARRAY'  )   ){
  	    $SECIS_out_HASH->{'0_2_0_mpk_all_SIS_sgm_ARRAY'}=$Mbk_all_SECISsgmARRAY;	
  	  }  	  
  	  
  	}
  	if (   ( defined ( $SECIS_HASH->{'0_1_1_best_SIS_sgm_ARRAY'} )  ) && (  ref ( $SECIS_HASH->{'0_1_1_best_SIS_sgm_ARRAY'} ) eq 'ARRAY'  )   ){
  	  my $bestSECISsgmARRAY=$SECIS_HASH->{'0_1_1_best_SIS_sgm_ARRAY'};
  	  my $Mbk_bestSECISsgmARRAY=SplignHandle::Map_est_sgmt_2_div_array_bk_to_geno($SplignOUTarray, $bestSECISsgmARRAY) ;   	  #DirFileHandle::PrintDumper($File_00_19_SIS_allSgmMbk, $Mbk_all_SECISsgmARRAY) if ( ref($Mbk_all_SECISsgmARRAY) eq 'ARRAY' ) ;  
  	  if (   ( defined ( $Mbk_bestSECISsgmARRAY )  ) && (  ref ( $Mbk_bestSECISsgmARRAY ) eq 'ARRAY'  )   ){
  	    $SECIS_out_HASH->{'0_2_1_mpk_bst_SIS_sgm_ARRAY'}=$Mbk_bestSECISsgmARRAY;	
  	  }  	  
  	  
  	}
  	
  	
  }
	return $SECIS_out_HASH;
}      
      




sub BuildAllSgmHash_formArrayDumpFILE{ #                               $HASH_00_21_all_gene_Sgmt  =ProSplignHandle::BuildAllSgmHash_formArrayDumpFILE ($HASH_00_21_all_gene_Sgmt, $nextRound_inform_Hash->{$ecSelPro}->{'analysisTypeHash'}->{'OneFiledGeno'}->{'File_00_02_PrSp_2_PsdSgm'}, "frame shift insert in Genome", '0_OnlyGeno', '01_genome_FrmShtIst_Array');   print "8,4,00b$ecSelPro\$GenoOrEst_total_number=$GenoOrEst_total_number\n";
    	        	                                                                                                                               
	my ($allSgmHash, $inSgmArrayFile, $inSgm_name, $inSgm_bigKey, $inSgm_smlKey)=@_;                                                                                                                                                                                                                                                    
	
	my $warnMsgBody="\nIn package  ProSplignHandle,\tIn sub BuildAllSgmHash_formArrayDumpFILE,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	                                                                                                                                                                                                                                                                                                                                  
	#print "8.4.000.$ecSelPro INPUT:$allSgmHash, $inSgmArrayFile, $inSgm_name, $inSgm_bigKey, $inSgm_smlKey, $ecSelPro";                                                                                                                                                                                                              
	my $inSgm_Array=Storable::retrieve($inSgmArrayFile) if (   ( defined ( $inSgmArrayFile )  ) && (  -e ( $inSgmArrayFile )  )   );  #if ( ($GeEsSpln_Yes ) && ($gnoSECIS_Yes) )   { my $tpHs; $tpHs->{'0_SgmInfom'}="SECIS found in Genome    "; $tpHs->{'1_SgmArray'}=$Best_Geno_SECISsgmArry;    $HASH_00_21_all_gene_Sgmt->{'2_genoMest'}->{'0_genome_Best_SECISArray'}=$tpHs;  }    	        	
    	        	                                                                                                                                                                                                                                                                                                                    
	if (   ( defined ( $inSgm_Array )  ) && (  ref ($inSgm_Array) eq 'ARRAY'  )   ){                                                                                                                                                                                                                                                                                              
		                                                                                                                                                                                                                                                                                                                                
		my $tpHs;                                                                                                                                                                                                                                                                                                                       
		$tpHs->{'0_SgmInfom'}=$inSgm_name;                                                                                                                                                                                                                                                                                              
		$tpHs->{'1_SgmArray'}=$inSgm_Array;                                                                                                                                                                                                                                                                                             
		$allSgmHash->{$inSgm_bigKey}->{$inSgm_smlKey}=$tpHs;                                                                                                                                                                                                                                                                            
		                                                                                                                                                                                                                                                                                                                                
	}                                                                                                                                                                                                                                                                                                                                 
	                                                                                                                                                                                                                                                                                                                                  
	return $allSgmHash;                                                                                                                                                                                                                                                                                                               
	                                                                                                                                                                                                                                                                                                                                  

}

sub BuildAllSgmHash_formSgmARRAY{ #                               $HASH_00_21_all_gene_Sgmt  =ProSplignHandle::BuildAllSgmHash_formSgmARRAY ($HASH_00_21_all_gene_Sgmt, $nextRound_inform_Hash->{$ecSelPro}->{'analysisTypeHash'}->{'OneFiledGeno'}->{'File_00_02_PrSp_2_PsdSgm'}, "frame shift insert in Genome", '0_OnlyGeno', '01_genome_FrmShtIst_Array');   print "8,4,00b$ecSelPro\$GenoOrEst_total_number=$GenoOrEst_total_number\n";
    	        	                                                                                                                               
	my ($allSgmHash, $inSgm_Array, $inSgm_name, $inSgm_bigKey, $inSgm_smlKey)=@_;                                                                                                                                                                                                                                                    
	
	my $warnMsgBody="\nIn package  ProSplignHandle,\tIn sub BuildAllSgmHash_formArrayDumpFILE,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	                                                                                                                                                                                                                                                                                                                                  
	#print "8.4.000.$ecSelPro INPUT:$allSgmHash, $inSgmArrayFile, $inSgm_name, $inSgm_bigKey, $inSgm_smlKey, $ecSelPro";                                                                                                                                                                                                              
	#my $inSgm_Array=Storable::retrieve($inSgmArrayFile) if (   ( defined ( $inSgmArrayFile )  ) && (  -e ( $inSgmArrayFile )  )   );  #if ( ($GeEsSpln_Yes ) && ($gnoSECIS_Yes) )   { my $tpHs; $tpHs->{'0_SgmInfom'}="SECIS found in Genome    "; $tpHs->{'1_SgmArray'}=$Best_Geno_SECISsgmArry;    $HASH_00_21_all_gene_Sgmt->{'2_genoMest'}->{'0_genome_Best_SECISArray'}=$tpHs;  }    	        	
    	        	                                                                                                                                                                                                                                                                                                                    
	if (   ( defined ( $inSgm_Array )  ) && (  ref ($inSgm_Array) eq 'ARRAY'  )   ){                                                                                                                                                                                                                                                                                              
		                                                                                                                                                                                                                                                                                                                                
		my $tpHs;                                                                                                                                                                                                                                                                                                                       
		$tpHs->{'0_SgmInfom'}=$inSgm_name;                                                                                                                                                                                                                                                                                              
		$tpHs->{'1_SgmArray'}=$inSgm_Array;                                                                                                                                                                                                                                                                                             
		$allSgmHash->{$inSgm_bigKey}->{$inSgm_smlKey}=$tpHs;                                                                                                                                                                                                                                                                            
		                                                                                                                                                                                                                                                                                                                                
	}                                                                                                                                                                                                                                                                                                                                 
	                                                                                                                                                                                                                                                                                                                                  
	return $allSgmHash;                                                                                                                                                                                                                                                                                                               
	                                                                                                                                                                                                                                                                                                                                  

}


sub BuildAllSgmHash_fromARRAY_forSECIS{   #$alSgmtAlSIS  =ProSplignHandle::BuildAllSgmHash_fromARRAY_forSECIS($alSgmtAlSIS, $all_SECISsgmARRAY, '0_OnlyGeno', 'Genome', '00_0', 'Gno_All__SECISArray' );
	my ($in_out_ARRAY, $in_ARRAY, $big_key_word, $SECIS_ID_head, $sm_keyWod_head, $sm_keyWod_tail)=@_;
	
	my $warnMsgBody="\nIn package  ProSplignHandle,\tIn sub BuildAllSgmHash_fromARRAY_forSECIS,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	
	if ( ref ($in_ARRAY) eq 'ARRAY' ){
		
		my $ArraySize=@{ $in_ARRAY };
		for (  my $i=0; $i<@{ $in_ARRAY }; $i++  ){ my $nb=$i+1;                   #print "8,5,03-2a $ecSelPro\$nb=\$i+1=$nb=$i+1\$in_ARRAY->[$i]=$in_ARRAY->[$i]\n"; 
			#DirFileHandle::PrintAndWarnDumper( $in_ARRAY->[$i] );                    #print "8,5,03-2b $ecSelPro\$in_ARRAY->[$i]->[0]->{'5_SgmtInfm'}="; print $in_ARRAY->[$i]->[0]->{'5_SgmtInfm'},"\n";
			my $showNumber= MatrixCsvChange::SprintfKeyHead( $ArraySize, $nb );
      my $Sgm_name= $SECIS_ID_head." SECIS $showNumber: ".$in_ARRAY->[$i]->[0]->{'5_SgmtInfm'};        #print "8,5,03-3 $ecSelPro\$Sgm_name=$Sgm_name\n"; 
			my $sml_key_Word="$sm_keyWod_head$showNumber$sm_keyWod_tail";                                #print "8,5,03-4 $ecSelPro\$sml_key_Word=$sml_key_Word\n"; 
			$in_out_ARRAY=ProSplignHandle::BuildAllSgmHash_formSgmARRAY(  
			                                                                               $in_out_ARRAY, 
			                                                                               $in_ARRAY->[$i], 
			                                                                               $Sgm_name, 
			                                                                               $big_key_word,  
			                                                                               $sml_key_Word
			                                                               			                                                               
			                                                                           );
			                                                            
			                                                               
			#print "8,5,03-3 $in_out_ARRAY=$in_out_ARRAY\n";
    }
  }
  return  $in_out_ARRAY;            
}



1;                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
                                                                                                                                           
########################################################################################################################################## 
#                                                                                                                                          