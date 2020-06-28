
#!/usr/bin/perl -w

BEGIN{  push (@INC, '/home/fredjiang/fredProgs/PMs')  };

use strict;
use warnings;

use Bio::Seq;
use Bio::SeqIO;
use Bio::Graphics;
use Bio::SearchIO;

use ExcelHandle;
use SeqSegmentsTools;
use PrintSubArrayHash;
use SeqSegmentsTools;
use FastaFileHandle;
use InFileHandle;

package  GeneWiseHandle;

################################################################################################################################################################################################################################################################################################
#
#     sub NewGenewisedbPraser{     #pmAble#  #这个是用来解析GeneWise的结果的
#     sub printHeadWordTail_Chg2Org{
#     sub PrintPngForGene{        #pmAble#     #打印基因结构
#     sub PrintPngForGeneGroups{   #pmAble#     打印成组的 gene 结构图#$proteindIDHash, $InIN_BothSecCysHash, $outFileHash, $pngFile    inKeyName    ($SecCys, $ProID, $InIN_BothSecCysWiseHash, $InIN_BothSecCysHash, $inKeyName, $pngFile, $SECISshowOrNot)
#     sub GetSmallestStart20170214{   #pmAble#  #PrintPngForGene 和 PrintPngForGeneGroups  的子函数，用来获得最小的起始点     
#     sub GetBiggestEnd20170214{   #pmAble#  #PrintPngForGene 和 PrintPngForGeneGroups 的子函数，用来获得最大的起终止点
#     sub AddTGACoreImform{    #这个是老版本， 已经被 新版本 AddTGACoreImform_20170307，所取代
#     sub AddSECISinform{   #pmAble#  #PrintPngForGene 和 PrintPngForGeneGroups 的子函数， #添加SECIS信息
#     sub AddTGACoreImform_20170307{        #pmAble#   # PrintPngForGene的子函数，添加TGAcore这个Exon类型，及相应的数据
#     sub combaineForPseudoGenes_20170307 { #pmAble#   # PrintPngForGene的子函数 #添加PseudoGene的类型 及相关数据
#     sub NEWAddFeatureTrack_into_Panel_20170307 {  #pmAble#   # PrintPngForGene #PrintPngForGene 的子函数，用来将Feature信息，加入画结构的模块
#     sub NEWAddFeatureTrack3Panel{  #这个应该是个老版本的，请使用新版本 NEWAddFeatureTrack_into_Panel_20170307
#     sub PrintPngForGeneGroups_old_2be_delete{ #$proteindIDHash, $InIN_BothSecCysHash, $outFileHash, $pngFile
#
################################################################################################################################################################################################################################################################################################

sub GeneWiseRun{
	my ($PepSeqFIle,  $contigSgm, $wiseOutFile)=@_;
	my $geneWiseDBCommandHead=" genewisedb -silent -quiet -codon /usr/share/doc/wise-doc/examples/codon.table  -m /usr/share/doc/wise-doc/examples/BLOSUM62.bla  -gene /usr/share/doc/wise-doc/examples/human.gf -prodb -trans -cdna -genes -pthread -pthr_no 1  -pretty -dnas -ace -sum -gff -gener";
  my $geneWsieWholeCommand="$geneWiseDBCommandHead $PepSeqFIle $contigSgm > $wiseOutFile";
  warn  "\n\n\nIn package GeneWiseHandle,\nIn sub GeneWiseRun,\n\n\$geneWsieWholeCommand=\n$geneWsieWholeCommand\n\n\n"; 
  print "\n\n\nIn package GeneWiseHandle,\nIn sub GeneWiseRun,\n\n\$geneWsieWholeCommand=\n$geneWsieWholeCommand\n\n\n"; 
  system ("$geneWsieWholeCommand");   
}

sub GetGffHash{    #输出wise中的gff信息
	my ($inFile)=@_;
	my $outBigHash=&NewGenewisedbPraser($inFile);
	my $outGffHash;
	if ( ref($outBigHash) eq 'HASH' ){
		if (  ref( $outBigHash->{'_genewiseResultHashAd'} ) eq 'HASH'  ){
		  foreach my $proteinNAME (    sort {$a cmp $b} (   keys ( %{ $outBigHash->{'_genewiseResultHashAd'} }  )   )    ){
		    if (  ref( $outBigHash->{'_genewiseResultHashAd'}->{$proteinNAME} ) eq 'ARRAY'  ){
		    	for (  my $i=0; $i<@{ $outBigHash->{'_genewiseResultHashAd'}->{$proteinNAME} }; $i++  ){
		    		if (  ref( $outBigHash->{'_genewiseResultHashAd'}->{$proteinNAME}->[$i] ) eq 'HASH'  ){
		    			foreach my $gff_number (    sort {$a <=> $b} (   keys ( %{ $outBigHash->{'_genewiseResultHashAd'}->{$proteinNAME}->[$i]->{'_notDividedGffHash'} }  )   )    ){
		    				if (  ref( $outBigHash->{'_genewiseResultHashAd'}->{$proteinNAME}->[$i]->{'_notDividedGffHash'}->{$gff_number} ) eq 'HASH'  ){
		    					$outGffHash->{'_genewiseResultHashAd'}->{$proteinNAME}->[$i]->{'_notDividedGffHash'}->{$gff_number}= $outBigHash->{'_genewiseResultHashAd'}->{$proteinNAME}->[$i]->{'_notDividedGffHash'}->{$gff_number};
		    				}
		    			}
		    		}
		    	}
		    } 
		  }
	  }
	}
	return $outGffHash;
}


sub OutPutGFFString{  #输出wise中的gff信息，以String的形式输出
	my ($inFile)=@_;
	my $outBigHash=&NewGenewisedbPraser($inFile);
	my $outGffString;
	if ( ref($outBigHash) eq 'HASH' ){
		if (  ref( $outBigHash->{'_genewiseResultHashAd'} ) eq 'HASH'  ){
		  foreach my $proteinNAME (    sort {$a cmp $b} (   keys ( %{ $outBigHash->{'_genewiseResultHashAd'} }  )   )    ){
		    if (  ref( $outBigHash->{'_genewiseResultHashAd'}->{$proteinNAME} ) eq 'ARRAY'  ){
		    	for (  my $i=0; $i<@{ $outBigHash->{'_genewiseResultHashAd'}->{$proteinNAME} }; $i++  ){
		    		if (  ref( $outBigHash->{'_genewiseResultHashAd'}->{$proteinNAME}->[$i] ) eq 'HASH'  ){
		    			foreach my $gff_number (    sort {$a <=> $b} (   keys ( %{ $outBigHash->{'_genewiseResultHashAd'}->{$proteinNAME}->[$i]->{'_notDividedGffHash'} }  )   )    ){
		    				
		    				if (  ref( $outBigHash->{'_genewiseResultHashAd'}->{$proteinNAME}->[$i]->{'_notDividedGffHash'}->{$gff_number} ) eq 'HASH'  ){
		    					my $lineHere;
		    					foreach my $gffColKey (    sort {$a cmp $b} (   keys ( %{ $outBigHash->{'_genewiseResultHashAd'}->{$proteinNAME}->[$i]->{'_notDividedGffHash'}->{$gff_number} }  )   )    ){
		    					  $lineHere.=$outBigHash->{'_genewiseResultHashAd'}->{$proteinNAME}->[$i]->{'_notDividedGffHash'}->{$gff_number}->{$gffColKey};
		    					  $lineHere.="\t";
		    					}
		    					$lineHere=~s/\t$//;	
		    					$outGffString.="$lineHere\n";	    					
		    				}
		    				
		    				
		    			}
		    			$outGffString.="\n";
		    		}
		    		
		    	}
		    	$outGffString.="\n";
		    } 
		  }
	  }
	}
	return $outGffString;
}

sub CheckTheBestWiseGeneIsANotAPseudoGene{  
	my ($wiseHash)=@_;
	my $theBestWiseGeneIsAPseudo_or_not=0;
	my $dieHeadMsg="\n\n\n		In package GeneWiseHandle,		\n In sub CheckTheBestWiseGeneIsANotAPseudoGene,\n";
	if (ref ($wiseHash) eq 'HASH'){
    my $bestDefine_word= '_BestScoreAllGeneHASH' ; #check the NewGenewisedbPraser method for this word!!
    my $bestPEP=                          $wiseHash->{$bestDefine_word}->{'_BestScoreQueryProteinID'};
    my $bestPEPSubNumber=                 $wiseHash->{$bestDefine_word}->{'_BestScoreQuProtNUMBER'}  ; 
    if (   
             (  defined ( $wiseHash->{'_genewiseResultHashAd'}->{$bestPEP}->[$bestPEPSubNumber]->{'_PsudoGeneFound'} )  ) 
         &&  ( $wiseHash->{'_genewiseResultHashAd'}->{$bestPEP}->[$bestPEPSubNumber]->{'_PsudoGeneFound'} )
       ) { 
      $theBestWiseGeneIsAPseudo_or_not=1;
    }
  }  
  
  return $theBestWiseGeneIsAPseudo_or_not;
  
}

sub GetBestGeneSequence_from_WisePraserOutHash{
	my ($wiseHash)=@_;
	my $dieHeadMsg="\n\n\n		In package GeneWiseHandle,		\n In sub GetBestGeneSequence_from_WisePraserOutHash,\n"; 
	print "$dieHeadMsg";
	print  "Input 1:\$wiseHash=$wiseHash\t";  print  "\n\n";
	
	if (ref ($wiseHash) eq 'HASH'){
    my $bestDefine_word= '_BestScoreAllGeneHASH' ;               #check the NewGenewisedbPraser method for this word!!
    my $bestPEP=                          $wiseHash->{$bestDefine_word}->{'_BestScoreQueryProteinID'};   print "$bestPEP=$bestPEP\n";     
    my $bestPEPSubNumber=                 $wiseHash->{$bestDefine_word}->{'_BestScoreQuProtNUMBER'}  ;   print "$bestPEPSubNumber=$bestPEPSubNumber\n";
    my $bestPepSubNumberHash=             $wiseHash->{'_genewiseResultHashAd'}->{$bestPEP}->[$bestPEPSubNumber];
    if (   (  defined ( $bestPepSubNumberHash->{'_PsudoGeneFound'} )  ) && ( $bestPepSubNumberHash->{'_PsudoGeneFound'} )   )    {     	 die "\n\n\n\nDIE!!!!\n\n$dieHeadMsg\n The best gene is  a pseudogene, please use the GetBestGeneSequence_from_WisePraserOutHash_for_pseudoGENE method!!\n\n\n";     }
    my $bestPepSubNumberBestGeneNumber=   $wiseHash->{$bestDefine_word}->{'_BestScoreGeneIdx'}       ;   print "$bestPepSubNumberBestGeneNumber=$bestPepSubNumberBestGeneNumber\n";
    
    my $bestGene_Hash                 =   $bestPepSubNumberHash->{'_EachGeneArrayAd'}->[$bestPepSubNumberBestGeneNumber];   DirFileHandle::PrintAndWarnDumper ($bestGene_Hash);
   
    my ($OutDnaString, $OutPepString)=( $bestGene_Hash->{'_EachGeneDNAseq'},  $bestGene_Hash->{'_EachGenePepSeq'} );        	print "\n\$OutDnaString=$OutDnaString\n\$OutPepString=$OutPepString\n\n";
   	if (    (  defined ( $bestGene_Hash->{'changedBackSeqWithU_DNA'} )  ) && ( $bestGene_Hash->{'changedBackSeqWithU_DNA'}=~m/\S+/ )    ){    		$OutDnaString=$bestGene_Hash->{'changedBackSeqWithU_DNA'};     print "\nWithU!!! \$OutDnaString=$OutDnaString\n"; }
   	if (    (  defined ( $bestGene_Hash->{'changedBackSeqWithU_PEP'} )  ) && ( $bestGene_Hash->{'changedBackSeqWithU_PEP'}=~m/\S+/ )    ){    		$OutPepString=$bestGene_Hash->{'changedBackSeqWithU_PEP'};     print "\nWithU!!! \$OutPepString=$OutPepString\n";	}
   	$OutDnaString=~s/\s+//g;   $OutPepString=~s/\s+//g;  
   	return [$OutDnaString, $OutPepString];
  }  
  
}

sub GetBestGeneSequence_from_WisePraserOutHash_for_pseudoGENE{
	my ($wiseHash)=@_;
	my $dieHeadMsg="\n\n\n		In package GeneWiseHandle,		\n In sub GetBestGeneSequence_from_WisePraserOutHash_for_pseudoGENE,\n"; 
	if (ref ($wiseHash) eq 'HASH'){
    my $bestDefine_word= '_BestScoreAllGeneHASH' ;                                                          #check the NewGenewisedbPraser method for this word!!
    my $bestPEP=                          $wiseHash->{$bestDefine_word}->{'_BestScoreQueryProteinID'};   
    my $bestPEPSubNumber=                 $wiseHash->{$bestDefine_word}->{'_BestScoreQuProtNUMBER'}  ;   
    my $bestPepSubNumberHash=             $wiseHash->{'_genewiseResultHashAd'}->{$bestPEP}->[$bestPEPSubNumber];
    if (   (  defined ( $bestPepSubNumberHash->{'_PsudoGeneFound'} )  ) && ( $bestPepSubNumberHash->{'_PsudoGeneFound'} )   )    {     	 }
    else {   die "$dieHeadMsg\n The best gene is  a pseudogene, please use the GetBestGeneSequence_from_WisePraserOutHash method!!\n\n\n";     }
    my $bestPepSubNumberBestGeneNumber=   $wiseHash->{$bestDefine_word}->{'_BestScoreGeneIdx'}       ;   
    
    my $bestPepSubNumberGeneArray=        $bestPepSubNumberHash->{'_EachGeneArrayAd'};
    
    my ($BigDNAstringWithOutIndel, $BigPEPstringWithOutIndel);
    if ( ref ($bestPepSubNumberGeneArray) eq 'ARRAY' ){
    	for (  my $i=0; $i<@{ $bestPepSubNumberGeneArray }; $i++  ){
    		my $bestGene_Hash                 =   $bestPepSubNumberGeneArray->{'_EachGeneArrayAd'}->[$i];  
        my ($OutDnaString, $OutPepString)=( $bestGene_Hash->{'_EachGeneDNAseq'},  $bestGene_Hash->{'_EachGenePepSeq'} );        	
   	    if (    (  defined ( $bestGene_Hash->{'changedBackSeqWithU_DNA'} )  ) && ( $bestGene_Hash->{'changedBackSeqWithU_DNA'}=~m/\S+/ )    ){    		$OutDnaString=$bestGene_Hash->{'changedBackSeqWithU_DNA'};      }
       	if (    (  defined ( $bestGene_Hash->{'changedBackSeqWithU_PEP'} )  ) && ( $bestGene_Hash->{'changedBackSeqWithU_PEP'}=~m/\S+/ )    ){    		$OutPepString=$bestGene_Hash->{'changedBackSeqWithU_PEP'};    	}
      	$OutDnaString=~s/\s+//g;        $BigDNAstringWithOutIndel.=$OutDnaString;  # if ($i>0){$BigDNAstringWithOutIndel.="!!!$OutDnaString";}
      	$OutPepString=~s/\s+//g;    	  $BigPEPstringWithOutIndel.=$OutPepString;  # if ($i>0){$BigDNAstringWithOutIndel.="!!!$OutPepString";}
    	}
    }
    return [$BigDNAstringWithOutIndel, $BigPEPstringWithOutIndel];       
  }  
  
}





sub GetBestGeneStructureFromGeneWisePraserOutHash{  #build the best gene structure 
	my ($wiseHash)=@_;
	my $dieHeadMsg="\n\n\n		In package GeneWiseHandle,		\n In sub GetBestGeneStructureFromGeneWisePraserOutHash,\n"; #print "$dieHeadMsg\n\$wiseHash=$wiseHash\n";
	if (ref ($wiseHash) eq 'HASH'){
    my $bestDefine_word= '_BestScoreAllGeneHASH' ;     # print "\$bestDefine_word=$bestDefine_word\n";             #check the NewGenewisedbPraser method for this word!!
    my $bestPEP=                          $wiseHash->{$bestDefine_word}->{'_BestScoreQueryProteinID'};   #print "\$bestPEP=$bestPEP\n";
    my $bestPEPSubNumber=                 $wiseHash->{$bestDefine_word}->{'_BestScoreQuProtNUMBER'}  ;   #print "\$bestPEPSubNumber=$bestPEPSubNumber\n";
    if (   
             (  defined ( $wiseHash->{'_genewiseResultHashAd'}->{$bestPEP}->[$bestPEPSubNumber]->{'_PsudoGeneFound'} )  ) 
         &&  ( $wiseHash->{'_genewiseResultHashAd'}->{$bestPEP}->[$bestPEPSubNumber]->{'_PsudoGeneFound'} )
       ) {  die "$dieHeadMsg\n The best gene is  a pseudogene, please use the GetBestGeneStructureFromGeneWisePraserOutHash_forPseudoGene method!!\n\n\n";  }
    else {  }
    
    my $bestPepSubNumberBestGeneNumber=   $wiseHash->{$bestDefine_word}->{'_BestScoreGeneIdx'}       ;   #print "\$bestPepSubNumberBestGeneNumber=$bestPepSubNumberBestGeneNumber\n";
    
    my $bestGene_ExonArray=     $wiseHash->{'_genewiseResultHashAd'}->{$bestPEP}->[$bestPEPSubNumber]->{'_EachGeneArrayAd'}->[$bestPepSubNumberBestGeneNumber]->{'_EachGeneExonArray'};
    my $bestGeneGffCDSifomHash= $wiseHash->{'_genewiseResultHashAd'}->{$bestPEP}->[$bestPEPSubNumber]->{'_GffImformationHash'}->{'Ordered_Match_Hash'}->{$bestPepSubNumberBestGeneNumber}->{'Ordered_CDS_Hash'};
    
    
    my $OutBestGeneARRAY;
    if ( ref($bestGene_ExonArray) eq 'ARRAY' ){
    	for (my $exonIdx=0; $exonIdx < @{ $bestGene_ExonArray }; $exonIdx++ ){
      $OutBestGeneARRAY->[$exonIdx]->{'0_SgmtHead'}=$bestGene_ExonArray->[$exonIdx]->{'_ExonHead'};       #print "\$OutBestGeneARRAY->[\$exonIdx]->{'0_SgmtHead'}=\$bestGene_ExonArray->[$exonIdx]->{'_ExonHead'}=\$OutBestGeneARRAY->[$exonIdx]->{'0_SgmtHead'}=$bestGene_ExonArray->[$exonIdx]->{'_ExonHead'}\n";
      $OutBestGeneARRAY->[$exonIdx]->{'1_SgmtTail'}=$bestGene_ExonArray->[$exonIdx]->{'_ExonTail'};       #print "\$OutBestGeneARRAY->[\$exonIdx]->{'1_SgmtTail'}=\$bestGene_ExonArray->[$exonIdx]->{'_ExonTail'}=\$OutBestGeneARRAY->[$exonIdx]->{'1_SgmtTail'}=$bestGene_ExonArray->[$exonIdx]->{'_ExonTail'}\n";
      $OutBestGeneARRAY->[$exonIdx]->{'2_SgmtZorF'}=$bestGeneGffCDSifomHash->{$exonIdx}->{'6stranZoF'};   #print "\$OutBestGeneARRAY->[\$exonIdx]->{'2_SgmtZorF'}=\$OutBestGeneARRAY->[\$exonIdx]->{'2_SgmtZorF'}=\$bestGeneGffCDSifomHash->{$exonIdx}->{'6stranZoF'}=$bestGeneGffCDSifomHash->{$exonIdx}->{'6stranZoF'}\n";
      $OutBestGeneARRAY->[$exonIdx]->{'3_SgmtFram'}=$bestGeneGffCDSifomHash->{$exonIdx}->{'7GffFrame'}+1;   #print "\$OutBestGeneARRAY->[\$exonIdx]->{'3_SgmtFram'}=\$OutBestGeneARRAY->[\$exonIdx]->{'3_SgmtFram'}=\$bestGeneGffCDSifomHash->{$exonIdx}->{'7GffFrame'}=$bestGeneGffCDSifomHash->{$exonIdx}->{'7GffFrame'}\n";
      }
      return $OutBestGeneARRAY;
    }    
  }  
}

sub ChangCDSposArrayBackToOrgContigPosition{  # If a got a cds gene structure from a segment of long contig, using the cds postion array, the segment start and end postion , segment strand, and the length of the original long contig, this method return a new array within the new postion and strand and frame information on the long contig!!!
	my ($inCDSposArray, $segmentCtgStt, $segmentCtgEnd, $segmentCtgZoF, $orgLongCtgLength)=@_;
	my $warnMstHead="\n\n\nIn package GeneWiseHandle,		\n In sub ChangCDSposArrayBackToOrgContigPosition,\n";
	my $dieHeadMsg="\n\n\nDIE!!!!!\n\n$warnMstHead"; 
	print $warnMstHead; print "1 \$inCDSposArray=$inCDSposArray\n"; print "2 \$segmentCtgStt=$segmentCtgStt\n"; print "3 \$segmentCtgEnd=$segmentCtgEnd\n"; print "4 \$segmentCtgZoF=$segmentCtgZoF\n"; print "5 \$orgLongCtgLength=$orgLongCtgLength\n"; print "\n\n\n";
	warn  $warnMstHead; warn  "1 \$inCDSposArray=$inCDSposArray\n"; warn  "2 \$segmentCtgStt=$segmentCtgStt\n"; warn  "3 \$segmentCtgEnd=$segmentCtgEnd\n"; warn  "4 \$segmentCtgZoF=$segmentCtgZoF\n"; warn  "5 \$orgLongCtgLength=$orgLongCtgLength\n"; warn  "\n\n\n";
	
	print "\n\n\n$segmentCtgStt, $segmentCtgEnd, $segmentCtgZoF, $orgLongCtgLength DirFileHandle::PrintAndWarnDumper (\$inCDSposArray)";
	DirFileHandle::PrintAndWarnDumper ($inCDSposArray);
	
	my $orgWholeLengthCtg_ZorF='+';
  my $newOutArray;
  if ( ($segmentCtgZoF eq '+') || ($segmentCtgZoF eq '-') ){                                                                   }
	else  {warn $dieHeadMsg; die "$dieHeadMsg\nThe \$segmentCtgZoF=$segmentCtgZoF\nIt should be + or -, otherwise die!!!!\n\n\n";}
	if ( ref ($inCDSposArray) eq 'ARRAY' ){
		for (  my $i=0; $i<@{ $inCDSposArray };$i++  ){  print "$segmentCtgStt, $segmentCtgEnd, $segmentCtgZoF, $orgLongCtgLength \$i=$i\nDirFileHandle::PrintAndWarnDumper (\$inCDSposArray->[$i])";  DirFileHandle::PrintAndWarnDumper ($inCDSposArray->[$i]);
			if ( ref ($inCDSposArray->[$i]) eq 'HASH' ){
			  my $oldSttPos=$inCDSposArray->[$i]->{'0_SgmtHead'}; my $newSttPos=SeqSegmentsTools::GetTheCorespondPosOFDIfferentSegement( 1, ($segmentCtgEnd-$segmentCtgStt+1), $segmentCtgZoF, $oldSttPos, $segmentCtgStt, $segmentCtgEnd, $orgWholeLengthCtg_ZorF);
			  my $oldSttEnd=$inCDSposArray->[$i]->{'1_SgmtTail'}; my $newSttEnd=SeqSegmentsTools::GetTheCorespondPosOFDIfferentSegement( 1, ($segmentCtgEnd-$segmentCtgStt+1), $segmentCtgZoF, $oldSttEnd, $segmentCtgStt, $segmentCtgEnd, $orgWholeLengthCtg_ZorF);
			  my $oldExoZoF=$inCDSposArray->[$i]->{'2_SgmtZorF'};
			  my $newExoZoF; if ($oldExoZoF eq $segmentCtgZoF){$newExoZoF='+';} else {$newExoZoF='-';} #{ ++ or --,$newExoZoF='+';  }  { +- or -+,$newExoZoF='-';  }  
			  
			  $newOutArray->[$i]->{'0_SgmtHead'}=$newSttPos;
        $newOutArray->[$i]->{'1_SgmtTail'}=$newSttEnd;
        $newOutArray->[$i]->{'2_SgmtZorF'}=$newExoZoF;
        
			  #frame work
			  if (  ( defined ($orgLongCtgLength) ) && ($orgLongCtgLength=~m/^\d+$/)  ){
			  	if (   (  defined ( $inCDSposArray->[$i]->{'3_SgmtFram'} )  ) && ( $inCDSposArray->[$i]->{'3_SgmtFram'}=~m/\d+/ )   ){
			  	  my $newExoFrm;  my $oldExoFrm=$inCDSposArray->[$i]->{'3_SgmtFram'}; my $oldAbsFrm=abs ($oldExoFrm); print "\n\n\$oldExoFrm=\$inCDSposArray->[$i]->{'3_SgmtFram'}=$oldExoFrm\n\n";  
			  	  my $realSgmStt=SeqSegmentsTools::getSmallOne($segmentCtgStt, $segmentCtgEnd); 				my $realSgmEnd=SeqSegmentsTools::getbigOne  ($segmentCtgStt, $segmentCtgEnd);
			  	  if ($oldExoZoF eq $segmentCtgZoF){					$newExoFrm=(   (  $oldAbsFrm + (  ($realSgmStt -1                 ) % 3 )  ) % 3   ); 			}  # ++ or --
			  	  else                             {					$newExoFrm=(   (  $oldAbsFrm + (  ($orgLongCtgLength-$realSgmEnd  ) % 3 )  ) % 3   ); 			}  # +- or -+
			  	  if ($newExoFrm==0){$newExoFrm=3;}    #$newExoFrm=$newExoZoF.$newExoFrm;	
			  	  $newOutArray->[$i]->{'3_SgmtFram'}=$newExoFrm;
			    }
			  }
			  
			  print "$segmentCtgStt, $segmentCtgEnd, $segmentCtgZoF, $orgLongCtgLength \$i=$i\nDirFileHandle::PrintAndWarnDumper (\$newOutArray)";  DirFileHandle::PrintAndWarnDumper ($newOutArray->[$i]);
			}
		}
		
	}
	print "\n\n\n$segmentCtgStt, $segmentCtgEnd, $segmentCtgZoF, $orgLongCtgLength DirFileHandle::PrintAndWarnDumper (\$newOutArray)";
	DirFileHandle::PrintAndWarnDumper ($newOutArray);
	return $newOutArray;
}

sub ChangCDSposArrayBackToOrgContigPosition_WithOut_FrameWork{  # If a got a cds gene structure from a segment of long contig, using the cds postion array, the segment start and end postion , segment strand, and the length of the original long contig, this method return a new array within the new postion and strand and frame information on the long contig!!!
	my ($inCDSposArray, $segmentCtgStt, $segmentCtgEnd, $segmentCtgZoF, $orgLongCtgLength)=@_;
	my $warnMstHead="\n\n\nIn package GeneWiseHandle,		\n In sub ChangCDSposArrayBackToOrgContigPosition_WithOut_FrameWork,\n";
	my $dieHeadMsg="\n\n\nDIE!!!!!\n\n$warnMstHead"; 
	print $warnMstHead; print "1 \$inCDSposArray=$inCDSposArray\n"; print "2 \$segmentCtgStt=$segmentCtgStt\n"; print "3 \$segmentCtgEnd=$segmentCtgEnd\n"; print "4 \$segmentCtgZoF=$segmentCtgZoF\n"; print "5 \$orgLongCtgLength=$orgLongCtgLength\n"; print "\n\n\n";
	warn  $warnMstHead; warn  "1 \$inCDSposArray=$inCDSposArray\n"; warn  "2 \$segmentCtgStt=$segmentCtgStt\n"; warn  "3 \$segmentCtgEnd=$segmentCtgEnd\n"; warn  "4 \$segmentCtgZoF=$segmentCtgZoF\n"; warn  "5 \$orgLongCtgLength=$orgLongCtgLength\n"; warn  "\n\n\n";
	
	print "\n\n\\nDirFileHandle::PrintAndWarnDumper (\$inCDSposArray)";
	DirFileHandle::PrintAndWarnDumper ($inCDSposArray);
	
	my $orgWholeLengthCtg_ZorF='+';
  my $newOutArray;
  if ( ($segmentCtgZoF eq '+') || ($segmentCtgZoF eq '-') ){                                                                   }
	else  {warn $dieHeadMsg; die "$dieHeadMsg\nThe \$segmentCtgZoF=$segmentCtgZoF\nIt should be + or -, otherwise die!!!!\n\n\n";}
	if ( ref ($inCDSposArray) eq 'ARRAY' ){
	 	for (  my $i=0; $i<@{ $inCDSposArray };$i++  ){  print "\$i=$i";
	    if  ( ref ($inCDSposArray->[$i]) eq 'HASH' ){
	  		my $oldSttPos=$inCDSposArray->[$i]->{'0_SgmtHead'}; my $newSttPos=SeqSegmentsTools::GetTheCorespondPosOFDIfferentSegement( 1, ($segmentCtgEnd-$segmentCtgStt+1), $segmentCtgZoF, $oldSttPos, $segmentCtgStt, $segmentCtgEnd, $orgWholeLengthCtg_ZorF);
	  		my $oldSttEnd=$inCDSposArray->[$i]->{'1_SgmtTail'}; my $newSttEnd=SeqSegmentsTools::GetTheCorespondPosOFDIfferentSegement( 1, ($segmentCtgEnd-$segmentCtgStt+1), $segmentCtgZoF, $oldSttEnd, $segmentCtgStt, $segmentCtgEnd, $orgWholeLengthCtg_ZorF);
	  		my $oldExoZoF=$inCDSposArray->[$i]->{'2_SgmtZorF'};
	  		my $newExoZoF; if ($oldExoZoF eq $segmentCtgZoF){$newExoZoF='+';} else {$newExoZoF='-';} #{ ++ or --,$newExoZoF='+';  }  { +- or -+,$newExoZoF='-';  }  
	  		
	  		$newOutArray->[$i]->{'0_SgmtHead'}=$newSttPos;
        $newOutArray->[$i]->{'1_SgmtTail'}=$newSttEnd;
        $newOutArray->[$i]->{'2_SgmtZorF'}=$newExoZoF;
        
	  		#frame work
	  		
	  	}
	  	
	  }
	}
	print "\n\n\\nDirFileHandle::PrintAndWarnDumper (\$newOutArray)";
	DirFileHandle::PrintAndWarnDumper ($newOutArray);
	return $newOutArray;
}

sub Get_GeneCDSposList_from_GeneStructureArray{
	my ($geneStructureArray)=@_;
	my $dieHeadMsg="\n\n\n		In package GeneWiseHandle,		\n In sub Get_GeneCDSposList_from_GeneStructureArray,\n"; #print "$dieHeadMsg\n\$wiseHash=$wiseHash\n";
	my $outList;
	if ( ref ($geneStructureArray) eq 'ARRAY' ){
		for (  my $i=0; $i<@{ $geneStructureArray };$i++  ){
			if ($i==0){
				$outList="$geneStructureArray->[$i]->{'0_SgmtHead'}..$geneStructureArray->[$i]->{'1_SgmtTail'}";
			}
			else {
				$outList.=",$geneStructureArray->[$i]->{'0_SgmtHead'}..$geneStructureArray->[$i]->{'1_SgmtTail'}";
			}
		}
		
	}
	else {
		DirFileHandle::PrintAndWarnDumper($geneStructureArray);
		die "$dieHeadMsg\nThe input \$geneStructureArray=$geneStructureArray should be a array ref!!";
	}
	return $outList;
}


sub Is_all_TGA_inframe_Check{
  my ($tga_pos_array, $exon_pos_array)=@_; 
  my $WarnHeadMsg="\n\n\n		In package GeneWiseHandle,		\n In sub Is_all_TGA_inframe_Check,\n";
  my $tga_pos_array_string=DirFileHandle::PrintAndWarnDumper($tga_pos_array);  my $exon_pos_array_string=DirFileHandle::PrintAndWarnDumper($exon_pos_array);
  my $Total_TGA_number=@{ $tga_pos_array };
  my $All_TGA_inframe_yes_or_not=0; 
  my $howManyTGAfoundInExon=0;
  if (  ( &CheckAllSegmentZorF($tga_pos_array) )&& (&CheckAllSegmentZorF($exon_pos_array) )  ){
  	if (  ( $tga_pos_array->[0]->{'2_SgmtZorF'} ) ne ( $exon_pos_array->[0]->{'2_SgmtZorF'} )  ){ 		warn"WarnHeadMsg, The TGA and the Exon is not in the same strand!\n\n\$tga_pos_array_string=$tga_pos_array_string\n\n\$exon_pos_array_string=$exon_pos_array_string\n\n\n";  	}
  	else{
  		for (  my $i=0; $i<@{ $tga_pos_array };$i++  ){
  			my $found_sig_UGA=0;
  		  for (  my $j=0; $j<@{ $exon_pos_array };$j++  ){
  		  	if ($tga_pos_array->[$i]->{'3_SgmtFram'} eq $exon_pos_array->[$j]->{'3_SgmtFram'}){
  		  		my ($tga_Stt, $tga_end, $Exo_Stt, $Exo_end)=( $tga_pos_array->[$i]->{'0_SgmtHead'}, $tga_pos_array->[$i]->{'1_SgmtTail'}, $exon_pos_array->[$j]->{'0_SgmtHead'}, $exon_pos_array->[$j]->{'1_SgmtTail'} );
  		  		if (  ( SeqSegmentsTools::a_between_b_c($tga_Stt, $Exo_Stt, $Exo_end) ) && ( SeqSegmentsTools::a_between_b_c($tga_Stt,  $Exo_Stt, $Exo_end) )  ){
  		  			$found_sig_UGA++;
  		  		}
  		  	}
  		  }
  		  if ($found_sig_UGA==1){
  		  	$howManyTGAfoundInExon++;
  		  }
  		  if ($found_sig_UGA>1){
  		  	die "\n\n\n\n\nDIE!!!!!\$WarnHeadMsg\nOne TGA found in more than one exons!!\nIt is impossible!!\n\n\$tga_pos_array_string=$tga_pos_array_string\n\n\$exon_pos_array_string=$exon_pos_array_string\n\n\n";
  		  }
  		}
  		
  	}
  	
  }
  else {
  	
  	die "\n\n\n\n\nDIE!!!!!$WarnHeadMsg\ndifferent + or - found in those two array \n\n\$tga_pos_array=$tga_pos_array\n$tga_pos_array_string\n\n\$exon_pos_array=$exon_pos_array\n$exon_pos_array_string\n\n";
  }
  if ($howManyTGAfoundInExon==$Total_TGA_number){
  	$All_TGA_inframe_yes_or_not=1;
  }
  return $All_TGA_inframe_yes_or_not;	
}

sub CheckAllSegmentZorF{
	my ($inSegmentArray)=@_;
	my $allTheSameDirection=0;
	my $fisrtSegZoF;
	if ( ref ($inSegmentArray) eq 'ARRAY' ){
		for (  my $i=0; $i<@{ $inSegmentArray };$i++  ){
			if ($i==0){
				$fisrtSegZoF=$inSegmentArray->[$i]->{'2_SgmtZorF'};
				$allTheSameDirection=1;
			}
			else {
				if ($inSegmentArray->[$i]->{'2_SgmtZorF'} ne $fisrtSegZoF){
					$allTheSameDirection=0;
				}
			}
		}
	}	
  return $allTheSameDirection;
}

sub Get_GeneStructureArray_from_GeneCDSposList{	
	my ($GeneCDSposList)=@_;
	my $dieHeadMsg="\n\n\n		In package GeneWiseHandle,		\n In sub Get_GeneStructureArray_from_GeneCDSposList,\n"; #print "$dieHeadMsg\n\$wiseHash=$wiseHash\n";
	my $outArray;
	my @cdsListArray=split ( ',',$GeneCDSposList ); 
	my $cdsIdx=0;
	foreach my $eachCDSlist (@cdsListArray){
		if ($eachCDSlist=~m/^(\d+)\.\.(\d+)$/){
			$outArray->[$cdsIdx]->{'0_SgmtHead'}=$1;
			$outArray->[$cdsIdx]->{'1_SgmtTail'}=$2;
			if    ($outArray->[$cdsIdx]->{'0_SgmtHead'} > $outArray->[$cdsIdx]->{'1_SgmtTail'}){
				$outArray->[$cdsIdx]->{'2_SgmtZorF'}='-'; 
				if ($cdsIdx>=1){if ($outArray->[$cdsIdx-1]->{'1_SgmtTail'}<=$outArray->[$cdsIdx]->{'0_SgmtHead'}){die "$dieHeadMsg\n\$GeneCDSposList=$GeneCDSposList\t\$outArray->[$cdsIdx]->{'2_SgmtZorF'}t\\$outArray->[$cdsIdx-1]->{'1_SgmtTail'}=$outArray->[$cdsIdx-1]->{'1_SgmtTail'}<=$outArray->[$cdsIdx]->{'0_SgmtHead'}=\$outArray->[$cdsIdx]->{'0_SgmtHead'}\nIt is not OK!!!";} }
			}
			elsif ($outArray->[$cdsIdx]->{'0_SgmtHead'} < $outArray->[$cdsIdx]->{'1_SgmtTail'}){
				$outArray->[$cdsIdx]->{'2_SgmtZorF'}='+';
				if ($cdsIdx>=1){if ($outArray->[$cdsIdx-1]->{'1_SgmtTail'}>=$outArray->[$cdsIdx]->{'0_SgmtHead'}){die "$dieHeadMsg\n\$GeneCDSposList=$GeneCDSposList\t\$outArray->[$cdsIdx]->{'2_SgmtZorF'}t\\$outArray->[$cdsIdx-1]->{'1_SgmtTail'}=$outArray->[$cdsIdx-1]->{'1_SgmtTail'}>=$outArray->[$cdsIdx]->{'0_SgmtHead'}=\$outArray->[$cdsIdx]->{'0_SgmtHead'}\nIt is not OK!!!";} }
			}
			$cdsIdx++;
		}
	}
	return $outArray;
}

sub GetBestGeneStructureFromGeneWisePraserOutHash_forPseudoGene{  #build the best gene structure 
	my ($wiseHash)=@_;
	my $dieHeadMsg="\n\n\n		In package GeneWiseHandle,		\n In sub GetBestGeneStructureFromGeneWisePraserOutHash_forPseudoGene,\n";
	if (ref ($wiseHash) eq 'HASH'){
    my $bestDefine_word= '_BestScoreAllGeneHASH' ; #check the NewGenewisedbPraser method for this word!!
    my $bestPEP=                          $wiseHash->{$bestDefine_word}->{'_BestScoreQueryProteinID'};
    my $bestPEPSubNumber=                 $wiseHash->{$bestDefine_word}->{'_BestScoreQuProtNUMBER'}  ; 
    if (   
             (  defined ( $wiseHash->{'_genewiseResultHashAd'}->{$bestPEP}->[$bestPEPSubNumber]->{'_PsudoGeneFound'} )  ) 
         &&  ( $wiseHash->{'_genewiseResultHashAd'}->{$bestPEP}->[$bestPEPSubNumber]->{'_PsudoGeneFound'} )
       ) { }
    else {    	die "$dieHeadMsg\n The best gene is not a pseudogene, please use the GetBestGeneStructureFromGeneWisePraserOutHash method!!\n\n\n";    }
    my $bestPepSubNumberBestGeneNumber=   $wiseHash->{$bestPEP}->{'_BestScoreGeneIdx'}       ;
    
   #Do the code work in the future!! 2018.07.04 
    
     
  }  
}



sub AddPsdGeneInformForWiseHash{  #add psudoinformation for wise praser out hash
	my ($wiseHash)=@_;
	if ( ref ($wiseHash) eq 'HASH' ){
		if (  ref ($wiseHash->{'_genewiseResultHashAd'} ) eq 'HASH'  ){
		  foreach my $eachPEP (    sort {$a cmp $b} (   keys(  %{ $wiseHash->{'_genewiseResultHashAd'} }  )   )    ){
		  	if (  ref ( $wiseHash->{'_genewiseResultHashAd'}->{$eachPEP} ) eq 'ARRAY'  ){
		  		for (  my $pepSubNbIdx=0; $pepSubNbIdx<@{ $wiseHash->{'_genewiseResultHashAd'}->{$eachPEP} }; $pepSubNbIdx++  ){
		  			if (  ref ( $wiseHash->{'_genewiseResultHashAd'}->{$eachPEP}->[$pepSubNbIdx]->{'_EachGeneArrayAd'} ) eq 'ARRAY'  ){
		  			  
		  			  my $geneNB=@{ $wiseHash->{'_genewiseResultHashAd'}->{$eachPEP}->[$pepSubNbIdx]->{'_EachGeneArrayAd'} };
		  			  if ($geneNB>1){
		  			  	$wiseHash->{'_genewiseResultHashAd'}->{$eachPEP}->[$pepSubNbIdx]->{'_PsudoGeneFound'}=1;
		  			  }
		  			  
		  			  for (  my $subGeneIdx=0; $subGeneIdx<@{ $wiseHash->{'_genewiseResultHashAd'}->{$eachPEP}->[$pepSubNbIdx]->{'_EachGeneArrayAd'} }; $subGeneIdx++  ){
		  			  	if ($subGeneIdx>0){
		  			  		my $match_ZoF=$wiseHash->{'_genewiseResultHashAd'}->{$eachPEP}->[$pepSubNbIdx]->{'_GffImformationHash'}->{'Ordered_Match_Hash'}->{$subGeneIdx  }->{'matchInform'}->{'6stranZoF'};
		  			  		my $preEndPos=$wiseHash->{'_genewiseResultHashAd'}->{$eachPEP}->[$pepSubNbIdx]->{'_GffImformationHash'}->{'Ordered_Match_Hash'}->{$subGeneIdx-1}->{'matchInform'}->{'4endPosti'};
		  			  		my $curSttPos=$wiseHash->{'_genewiseResultHashAd'}->{$eachPEP}->[$pepSubNbIdx]->{'_GffImformationHash'}->{'Ordered_Match_Hash'}->{$subGeneIdx  }->{'matchInform'}->{'3startPos'};
		  			  		if    ($match_ZoF eq '+'){
		  			  		  $wiseHash->{'_genewiseResultHashAd'}->{$eachPEP}->[$pepSubNbIdx]->{'_PsudoArray'}->[$subGeneIdx-1]->{'_indel_start'}= $preEndPos+1;
		  			  		  $wiseHash->{'_genewiseResultHashAd'}->{$eachPEP}->[$pepSubNbIdx]->{'_PsudoArray'}->[$subGeneIdx-1]->{'_indel___end'}= $curSttPos-1;	
		  			  		}
		  			  		elsif ($match_ZoF eq '-'){
		  			  		  $wiseHash->{'_genewiseResultHashAd'}->{$eachPEP}->[$pepSubNbIdx]->{'_PsudoArray'}->[$subGeneIdx-1]->{'_indel_start'}= $preEndPos-1;
		  			  		  $wiseHash->{'_genewiseResultHashAd'}->{$eachPEP}->[$pepSubNbIdx]->{'_PsudoArray'}->[$subGeneIdx-1]->{'_indel___end'}= $curSttPos+1;	
		  			  		}
		  			  	}
		  			  }
		  			  
		  			  
		  			}
		  		}
		  	}
		  }	
		}
	  
	}
	return $wiseHash;
}

sub NewGenewisedbPraser{     #pmAble#  #这个是用来解析GeneWise的结果的
  my (
       $geneWisedbOutFile,  #1 
       $stopPos,            #2
       $orgThreeWd,         #3
       $orgResidue,         #4
       $DNAstart,           #5                     #后面是设想，但没有这么实现 #这个参数 可以写 一类：数字，强制该片段是从该位置开始， 二类：'fastacmdStart',则会读取contig的序列名，解析开头位置（fastacmd必须为正向）， 三类：留空，则从1开始
       $DNAend,             #6                     #后面是设想，但没有这么实现 #这个参数 可以写 一类：数字，强制该片段是在该位置结束， 二类：'fastacmdend',  则会读取contig的序列名，解析开头位置（fastacmd必须为正向）， 三类：留空，则计算该序列长度为该值
       $DNAchangedHash      #7
     )=@_;    #输入时 用来接卸的 genewisedb 输出结果文件，以及Stop（实际上是TGA——U）的位置, 原始Stop位置的 三联DNA密码子， 原始 stop位置的 三联密码子翻译成的氨基酸
  
  warn     "\n\nNow in Package GeneWiseHandle,\nIn Sub &NewGenewisedbPraser!\n";      warn  "Input 1:\$geneWisedbOutFile=$geneWisedbOutFile\t" if (defined ($geneWisedbOutFile)); warn  "Input 2:\$stopPos=$stopPos\t" if (defined ($stopPos));
  print "\n\cl\nNow in Package GeneWiseHandle,\nIn Sub &NewGenewisedbPraser!\n:";     print "Input 1:\$geneWisedbOutFile=$geneWisedbOutFile\t" if (defined ($geneWisedbOutFile)); print "Input 2:\$stopPos=$stopPos\t" if (defined ($stopPos)); 
  warn  "Input 3:\$orgThreeWd=$orgThreeWd\t" if (defined ($orgThreeWd));  warn  "Input 4:\$orgResidue=$orgResidue\t" if (defined ($orgResidue));  
  print "Input 3:\$orgThreeWd=$orgThreeWd\t" if (defined ($orgThreeWd));  print "Input 4:\$orgResidue=$orgResidue\t" if (defined ($orgResidue));  
  warn  "Input 5:\$DNAstart=$DNAstart\t" if (defined ($DNAstart));  warn  "Input 6:\$DNAend=$DNAend\t" if (defined ($DNAend));  warn  "Input 7:\$DNAchangedHash=$DNAchangedHash\t" if (defined ($DNAchangedHash)); 
  print "Input 5:\$DNAstart=$DNAstart\t" if (defined ($DNAstart));  print "Input 6:\$DNAend=$DNAend\t" if (defined ($DNAend));  print "Input 7:\$DNAchangedHash=$DNAchangedHash\t" if (defined ($DNAchangedHash));
  warn  "\n\n";
  print "\n\cl\n";
  
  open (WISEFILE,"$geneWisedbOutFile") or die "cannot open \$geneWisedbOutFile=$geneWisedbOutFile : $! \n";   
  
  
  my $FinalOutWiseHASH; #建立一个总输出Hash。第一级的key是 protein的ID，
  
  #my %keyHash={
  #  '_QueryProteinFilePath'           => "用于和ＤＮＡ比对的protein的文件的路径",
  #    '_targetDNAfilePath'              => "用于和Protein比对的基因组DNA文件的路径，这个基因组DNA文件，是被切割处理过了的",
  #    '_targetDNA_sgmentStart'          => "这个切割处理过的DNA文件，其 开头 位置是原基因组的多少位置",
  #    '_targetDNA_sgmentEnd'            => "这个切割处理过的DNA文件，其 结尾 位置是原基因组的多少位置",
  #    '_genewiseResultHashAd'           => "genewiseDB的输出结果， 由全部Protein为Key形成的Hash",
  #      
  #      '_proteinAlinmentProName'         => "蛋白和DNA比对全部结果中，某个Protein结果的，query蛋白名",        
  #      '_proteinAlinmentDNAname'         => "蛋白和DNA比对全部结果中，某个Protein结果的，target DNA名",
  #      '_proteinAlinmentTotalScore'      => "蛋白和DNA比对全部结果中，某个Protein结果的，总得分，如果有pseudo出现，则可能是几个gene的得分之和",
  #      '_proteinAlinmentQueryPtHead'     => "蛋白和DNA比对全部结果中，某个Protein结果的，query蛋白 比对Match区域的 开头 的位置数",
  #      '_proteinAlinmentQueryPtTail'     => "蛋白和DNA比对全部结果中，某个Protein结果的，query蛋白 比对Match区域的 结尾 的位置数",
  #      '_proteinAlinmentDNAtrgtHead'     => "蛋白和DNA比对全部结果中，某个Protein结果的，DNA target比对Match区域的 开头 的位置数",
  #      '_proteinAlinmentDNAtrgtTail'     => "蛋白和DNA比对全部结果中，某个Protein结果的，DNA target比对Match区域的 结尾 的位置数",
  #
  #      
  #      '_proteinAlinmentHashAd'          => "蛋白和DNA比对全部结果中，某个Protein结果的，Protein Seq 和 DNA seq的 alignment的所有信息",
  #      '_proteinAlinmentDNAposKeyHASH'   => "蛋白和DNA比对全部结果中，某个Protein结果的，以DNA位置为key的，获知每个DNA位置的比对情况，如该位置是CDS还是Intron还是Pseudo，等信息",
  #    
  #    '_EachGeneArrayAd'                => "蛋白和DNA比对全部结果中，某个Protein结果的，由全部比对获得得Gene，所组成的Array（通常只有1个，如果有pseudo出现，则可能是多个）",
  #      '_EachGeneHeadPos'                => "蛋白和DNA比对全部结果中，某个Protein结果的，其中某个Gene的，开头的 位置",
  #      '_EachGeneTailPos'                => "蛋白和DNA比对全部结果中，某个Protein结果的，其中某个Gene的，结尾的 位置",
  #      '_EachGenePepSeq'                 => "蛋白和DNA比对全部结果中，某个Protein结果的，其中某个Gene的，氨基酸序列",  
  #      '_EachGeneDNAseq'                 => "蛋白和DNA比对全部结果中，某个Protein结果的，其中某个Gene的，DNA的序列",
  #      '_EachGeneSubSeqScore'            => "蛋白和DNA比对全部结果中，某个Protein结果的，其中某个Gene的，subSequence的得分，也就是这个 gene的得分",
  #      '_EachGeneExonArray'              => "蛋白和DNA比对全部结果中，某个Protein结果的，其中某个Gene的，所有Exon，所责成的Array",
  #      '_ExonHead'                       => "蛋白和DNA比对全部结果中，某个Protein结果的，其中某个Gene的，其中某个Exon的，开头的 位置",
  #      '_ExonTail'                       => "蛋白和DNA比对全部结果中，某个Protein结果的，其中某个Gene的，其中某个Exon的，结尾的 位置",
  #      '_Exonphase'                      => "蛋白和DNA比对全部结果中，某个Protein结果的，其中某个Gene的，其中某个Exon的，pahse，这个通常有开头DNA有几个没有在一个完整的三联密码子中来决定",
  #      
  #      '_ExonInGeneHead'                 => "蛋白和DNA比对全部结果中，某个Protein结果的，其中某个Gene的，其中某个Exon的，开头的 位置 结算后变成 从该gene的第一个DNA碱基开始结算，到该位置的数字",
  #      '_ExonInGeneTail'                 => "蛋白和DNA比对全部结果中，某个Protein结果的，其中某个Gene的，其中某个Exon的，结尾的 位置 结算后变成 从该gene的第一个DNA碱基开始结算，到该位置的数字",
  #                                      
  #    '_ExonHeadAAcodenPos'             => "蛋白和DNA比对全部结果中，某个Protein结果的，其中某个Gene的，其中某个Exon的，开头的 位置 计算后变成 该位置在三联密码子中的位置 0 1 2 ",
  #      '_ExonHeadAApos'                  => "蛋白和DNA比对全部结果中，某个Protein结果的，其中某个Gene的，其中某个Exon的，开头的 位置 计算后变成 氨基酸序列的 位置", 
  #      '_ExonTailAACodenPos'             => "蛋白和DNA比对全部结果中，某个Protein结果的，其中某个Gene的，其中某个Exon的，结尾的 位置 计算后变成 该位置在三联密码子中的位置 0 1 2 ",
  #      '_ExonTailAApos'                  => "蛋白和DNA比对全部结果中，某个Protein结果的，其中某个Gene的，其中某个Exon的，结尾的 位置 计算后变成 氨基酸序列的 位置",
  #                                      
  #    '_TGAinformHash'                  => "这个是用来装载 TGA-U 相关的信息，包括是否存在 tga-u，包含 TGA-U的外显子的序号，以及TGA-U的位置",
  #      '_TGAtransForUexists'             => "蛋白和DNA比对全部结果中，某个Protein结果的，其中某个Gene的，TGA-U信息，这里记录该gene是否 有 TGA-U",
  #      '_TGAincludeExonNB'               => "蛋白和DNA比对全部结果中，某个Protein结果的，其中某个Gene的，TGA-U信息，这里记录该gene 中有 TGA-U 的外显子的 序数",
  #      '_TGAPosition'                    => "蛋白和DNA比对全部结果中，某个Protein结果的，其中某个Gene的，TGA-U信息，这里记录该gene  TGA-U 的 位置数",
  #      '_UaaPostionNum'                  => "记录wise结果中，各个有U的 gene的 U在pep中的位置",                                                                        
  #      '_TGAchangedPEP'                  => "记录wise结果中，各个有U的 gene的 pep序列，这里原来被转换为C的U被重新转换回来",
  #      '_TGAchangedDNA'                  => "记录wise结果中，各个有U的 gene的 DNA序列，这里原来被转换为TGC的TGA被重新转换回来",
  #                                                                             
  #    '_pseudoSegmentAr'                => "蛋白和DNA比对全部结果中，某个Protein结果的，Pseudo位点片段所组成的 Array",
  #      '_PseudoSegmentHead'              => "蛋白和DNA比对全部结果中，某个Protein结果的， 某个 Pseudo位点片段的， 开头 位置",
  #      '_PseudoSegmentTail'              => "蛋白和DNA比对全部结果中，某个Protein结果的， 某个 Pseudo位点片段的， 结尾 位置",
  #      '_PseudoSegmentLength'            => "蛋白和DNA比对全部结果中，某个Protein结果的， 某个 Pseudo位点片段的， 长度",
  #    
  #    '_BestTotalScoreWithUhash'        => "所有 query protein中，总的比对 得分（如出现pseudo 则，可能是几个 gene的比对得分之和），总得分最高，且有TGA-U的那个gene的hash。这个应该不是正确的，但以前的程序是按这个标准来获取蛋白质序列的",
  #    '_BestScoreWithUhash'             => "所有比对结果中，                      得分最高的 Gene（注意：是gene而不是比对Result），的信息Hash。该Hash内部主要有 该 gene是哪个 Protein Query,相同query protein的序数（如无相同则该数为0）， 和 哪个Gene 序数，这3个值被保存",
  #    '_BestScoreWithOutUhash'          => "所有比对结果中，      不 含有 TGA-U的 得分最高的 Gene（注意：是gene而不是比对Result），的信息Hash。该Hash内部主要有 该 gene是哪个 Protein Query,相同query protein的序数（如无相同则该数为0）， 和 哪个Gene 序数，这3个值被保存",
  #    '_BestScoreAllGeneHASH'           => "所有比对结果中，不管 含不含有 TGA-U的 得分最高的 Gene（注意：是gene而不是比对Result），的信息Hash。该Hash内部主要有 该 gene是哪个 Protein Query,相同query protein的序数（如无相同则该数为0）， 和 哪个Gene 序数，这3个值被保存",
  #      '_BestScoreQueryProteinID'        => "得分最高的 蛋白ID，                                                       这个可以用来在 _genewiseResultHashAd 中 ，找到 _proteinAlinmentProName 的 Key",
  #      '_BestScoreQuProtNUMBER'          => "得分最高的 蛋白ID，同一个蛋白的出现第几次的序数，通常为0，                这个可以用来在 _genewiseResultHashAd 中 ，找到 _proteinAlinmentProName 的 Key,后继续找到 相同Key的蛋白的结果中，第几个结果（从0开始数）",
  #      '_BestScoreGeneIdx'               => "得分最高的 蛋白ID, 同一个蛋白的出现第几次的序数，通常为0，中的gene的序数，这个可以用来在 _genewiseResultHashAd 中 ，找到 _proteinAlinmentProName 的 Key,后继续找到 相同Key的蛋白的结果中，第几个结果（从0开始数），中的_EachGeneArrayAd 中的 \$bgIdx，其实也就是 gene的序数"
  #
  #
  #};
  
  #以下，即为本 sub 获得得最终 hash的内部信息形式。
  #$FinalOutWiseHASH->{'_QueryProteinFilePath'}=$1;
  #$FinalOutWiseHASH->{'_targetDNAfilePath'}=$DNAtargetFile=$2;
  #$FinalOutWiseHASH->{'_targetDNA_sgmentStart'}=$upstrmNb=$1;  warn "\$FinalOutWiseHASH->{'_targetDNA_sgmentStart'}=$FinalOutWiseHASH->{'_targetDNA_sgmentStart'}\n";
  #$FinalOutWiseHASH->{'_targetDNA_sgmentEnd'}=$2;              warn "\$FinalOutWiseHASH->{'_targetDNA_sgmentEnd'}  =$FinalOutWiseHASH->{'_targetDNA_sgmentEnd'}\n";
  #
  #$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_proteinAlinmentDNAname'}    =$targetDNAname;   
  #$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_proteinAlinmentDirection'}  =$forwORrevCHGnb;
  #$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_proteinAlinmentTotalScore'} =$wiseScore;
  #$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_proteinAlinmentQueryPtHead'}=$queryProtHd;
  #$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_proteinAlinmentQueryPtTail'}=$queryProtTl;
  #$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_proteinAlinmentDNAtrgtHead'}=$targetDNAhd;
  #$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_proteinAlinmentDNAtrgtTail'}=$targetDNAtl;

  #$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_pseudoSegmentAr'}->[$pseudoSegMentIdx]->{'_PseudoSegmentHead'}  =$firstPseudoPos;
  #$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_pseudoSegmentAr'}->[$pseudoSegMentIdx]->{'_PseudoSegmentTail'}  =$lastPseudoPos;
  #$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_pseudoSegmentAr'}->[$pseudoSegMentIdx]->{'_PseudoSegmentLength'}=$groupEdhash->{ '_PseudoCol' }->{ $alignINFORMhash->[$colH]->{'_GroupNB'} }->{'_PseudoNUMBER'};
  #$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_proteinAlinmentHashAd'}       =$WholeAlignHash;  #这是个很复杂的hash，主要记录了 Aliment各行各列的具体信息
  #$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_proteinAlinmentDNAposKeyHASH'}=$DNAposKeyHASH;   #这个也很复杂，记录了 比对上的所有DNA位置的 具体信息，包括类型（CDA，INTRON还是Pseudo），翻译成的氨基酸，以及内部小组序数等
  #$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneHeadPos'}=$geneStart;   
  #$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneTailPos'}=$geneEnd; 
  #$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_TGAinformHash'}->{'_TGAtransForUexists'}=0; #每个geneweie预测的基因，在最初录入 总$FinalOutWiseHASH时，都其每个 gene是否有TGA-U初始化为0.
  #$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGenePepSeq'}=$pepSeq;  
  #$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneSubSeqScore'}=$subSeqScore;
  #$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneDNAseq'}=$cDNASeq; 
  #
  #$FinalOutWiseHASH->{'_BestScoreAllGeneHASH'}->{'_BestScoreQueryProteinID'}=$QueryProteinName;    print "\$FinalOutWiseHASH->{'_BestScoreAllGeneHASH'}->{'_BestScoreQueryProteinID'}=\$QueryProteinName=$QueryProteinName\n";
  #$FinalOutWiseHASH->{'_BestScoreAllGeneHASH'}->{'_BestScoreQuProtNUMBER'}=$sameProtNameNumber;    print "\$FinalOutWiseHASH->{'_BestScoreAllGeneHASH'}->{'_BestScoreQuProtNUMBER'}=\$_BestScoreQuProtNUMBER=$_BestScoreQuProtNUMBER\n";
  #$FinalOutWiseHASH->{'_BestScoreAllGeneHASH'}->{'_BestScoreGeneIdx'}=$bgIdx;	        	           print "\$FinalOutWiseHASH->{'_BestScoreAllGeneHASH'}->{'_BestScoreGeneIdx'}=\$bgIdx=$bgIdx\n\n";
  #
  #$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneExonArray'}->[$exonIdx]->{'_ExonHead'}=$orgNbA;
  #$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneExonArray'}->[$exonIdx]->{'_ExonTail'}=$orgNbB;
  #$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneExonArray'}->[$exonIdx]->{'_Exonphase'}=$3;
  #$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneExonArray'}->[$exonIdx]->{'_ExonInGeneHead'}=$exonADDlength+1;
  #$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneExonArray'}->[$exonIdx]->{'_ExonInGeneTail'}=$exonADDlength+abs($orgNbB-$orgNbA)+1;
  #$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneExonArray'}->[$exonIdx]->{'_ExonHeadAAcodenPos'}=($FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneExonArray'}->[$exonIdx]->{'_ExonInGeneHead'})%3;
  #$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneExonArray'}->[$exonIdx]->{'_ExonHeadAApos'} = int ( ($FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneExonArray'}->[$exonIdx]->{'_ExonInGeneHead'})/3 );
  #$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneExonArray'}->[$exonIdx]->{'_ExonHeadAApos'} = int ( ($FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneExonArray'}->[$exonIdx]->{'_ExonInGeneHead'})/3 ) + 1;
  #$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneExonArray'}->[$exonIdx]->{'_ExonTailAACodenPos'}=($FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneExonArray'}->[$exonIdx]->{'_ExonInGeneTail'})%3;
  #$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneExonArray'}->[$exonIdx]->{'_ExonTailAApos'} = int ( ($FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneExonArray'}->[$exonIdx]->{'_ExonInGeneTail'})/3 );
  #$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneExonArray'}->[$exonIdx]->{'_ExonTailAApos'} = int ( ($FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneExonArray'}->[$exonIdx]->{'_ExonInGeneTail'})/3 ) + 1;
  #$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_TGAinformHash'}->{'_TGAtransForUexists'}=1;
  #$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_TGAinformHash'}->{'_TGAincludeExonNB'}  =$exonIdx;
  #$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_TGAinformHash'}->{'_TGAPosition'}       =$stopPos;
  #$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_TGAinformHash'}->{'_TGAEndPs'}          =$stopPos+2*$forwORrevCHGnb;
  #
  #$FinalOutWiseHASH->{'_BestScoreWithUhash'}->{'_BestScoreQueryProteinID'}=$QueryProteinName;                                                                                                print "\$FinalOutWiseHASH->{'_BestScoreWithUhash'}->{'_BestScoreQueryProteinID'}=\$QueryProteinName=$QueryProteinName\n";  
  #$FinalOutWiseHASH->{'_BestScoreWithUhash'}->{'_BestScoreQuProtNUMBER'}=$sameProtNameNumber;                                                                                                print "\$FinalOutWiseHASH->{'_BestScoreWithUhash'}->{'_BestScoreQuProtNUMBER'}=\$_BestScoreQuProtNUMBER=$_BestScoreQuProtNUMBER\n";
  #$FinalOutWiseHASH->{'_BestScoreWithUhash'}->{'_BestScoreGeneIdx'}=$bgIdx;	        	                                                                                                      print "\$FinalOutWiseHASH->{'_BestScoreWithUhash'}->{'_BestScoreGeneIdx'}=\$bgIdx=$bgIdx\n\n";                               
  #
  #$FinalOutWiseHASH->{'_BestTotalScoreWithUhash'}->{'_BestScoreQueryProteinID'}=$QueryProteinName;                                                                    print "\$FinalOutWiseHASH->{'_BestTotalScoreWithUhash'}->{'_BestScoreQueryProteinID'}=\$QueryProteinName=$QueryProteinName\n";                                                                                                                                                                                                             
  #$FinalOutWiseHASH->{'_BestTotalScoreWithUhash'}->{'_BestScoreQuProtNUMBER'}=$sameProtNameNumber;                                                                    print "\$FinalOutWiseHASH->{'_BestTotalScoreWithUhash'}->{'_BestScoreQuProtNUMBER'}=\$_BestScoreQuProtNUMBER=$_BestScoreQuProtNUMBER\n";
  #$FinalOutWiseHASH->{'_BestTotalScoreWithUhash'}->{'_BestScoreGeneIdx'}=$bgIdx;	        	                                                                           print "\$FinalOutWiseHASH->{'_BestTotalScoreWithUhash'}->{'_BestScoreGeneIdx'}=\$bgIdx=$bgIdx\n\n";                                                                                                                                                                                                                                          
  #
  #$FinalOutWiseHASH->{'_BestScoreWithOutUhash'}->{'_BestScoreQueryProteinID'}=$QueryProteinName;                                                                     print "\$FinalOutWiseHASH->{'_BestScoreWithOutUhash'}->{'_BestScoreQueryProteinID'}=\$QueryProteinName=$QueryProteinName\n";                                                                                                                                          
  #$FinalOutWiseHASH->{'_BestScoreWithOutUhash'}->{'_BestScoreQuProtNUMBER'}=$sameProtNameNumber;                                                                     print "\$FinalOutWiseHASH->{'_BestScoreWithOutUhash'}->{'_BestScoreQuProtNUMBER'}=\$_BestScoreQuProtNUMBER=$_BestScoreQuProtNUMBER\n";
  #$FinalOutWiseHASH->{'_BestScoreWithOutUhash'}->{'_BestScoreGeneIdx'}=$bgIdx;                                                                                       print "\$FinalOutWiseHASH->{'_BestScoreWithOutUhash'}->{'_BestScoreGeneIdx'}=\$bgIdx=$bgIdx\n\n";                                                                                                                                                                       


	#$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_TGAinformHash'}->{'_UaaPostionNum'}  =$UaaPos;       #记录wise结果中，各个有U的 gene的 U在pep中的位置
	#$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_TGAinformHash'}->{'_TGAchangedPEP'}  =$ChgBakPep;    #记录wise结果中，各个有U的 gene的 pep序列，这里原来被转换为C的U被重新转换回来
	#$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_TGAinformHash'}->{'_TGAchangedDNA'}  =$ChgBakCDNA;   #记录wise结果中，各个有U的 gene的 DNA序列，这里原来被转换为TGC的TGA被重新转换回来


  
	#my $arrForSort; $aFsIdx=0;  my $noUarrForSort; $NUaFsIdx=0;  
	my $tempMk=$/; $/=">Results for "; 
  my $tepNb=0;  
	my $DNAtargetFile;  # 这个变量，记录Genewisedb输入的 被比对的DNA文件的 绝对路径
	my $upstrmNb; #  被比对的DNA文件，是原DNA文件的片段，这个变量$upstrmNb，记录该片段的起点位置
	
	my $BestTotalScoreWithU=-99999999;           #该值用来保存 有           U的       gene的最高分（但是这个分数，是一个protein query比对得到的 一个或多个 gene的总得分， 和下面的值$BestScoreGeneWithU_Score不一样）
	my $BestScoreGeneWithU_Score=-99999999;      #该值用来保存 有           U的       gene的最高分 （这个分数，是单个gene的最高分） 这个最高分和上一个最高分，可能指向不同的 protein query
	my $BestScoreGeneWithOutU_Score=-99999999;   #该值用来保存 没有         U的       gene的最高分
	my $BestScoreGene_Score=-99999999;           #该值用来保存 不管有没有   U的 所有的gene的最高分
	
	while (<WISEFILE>){ 
		my $eachWise=$_;    warn "\$tepNb=$tepNb\n";
		
		if ($tepNb == 0){  #Dna info from:         /home/fredjiang/work/Algae/20150522NewDATA/AllSpecies/Chlamydomonas_reinhardtii20150522/CodeType_4_OutPut/a_GtbnSe_OT/Tmp_ChG/10
		  if ($eachWise=~m/Protein\s+info\s+from:\s+(\S+)\nDna\s+info\s+from:\s+(\S+)\n/){
		    $FinalOutWiseHASH->{'_QueryProteinFilePath'}=$1;
		    $FinalOutWiseHASH->{'_targetDNAfilePath'}=$DNAtargetFile=$2;
		    
		    #由于很多片段序列都是用fastacmd获得的，所以可以先 解析得到其靠头和结尾 并填入hash中
		    my $fastaFileString=InFileHandle::readAllfileIntoAstring($DNAtargetFile);
		    my $fastaFileSeqHead=FastaFileHandle::GetSeqFastaHeadName($fastaFileString);
		    if ($fastaFileSeqHead=~m/>\S+.*\s+(\d+)\.\.(\d+)/){
		    	$FinalOutWiseHASH->{'_targetDNA_sgmentStart'}=$upstrmNb=$1;  warn "\$FinalOutWiseHASH->{'_targetDNA_sgmentStart'}=$FinalOutWiseHASH->{'_targetDNA_sgmentStart'}\t\$fastaFileSeqHead=$fastaFileSeqHead\n";
		    	$FinalOutWiseHASH->{'_targetDNA_sgmentEnd'}=$2;              warn "\$FinalOutWiseHASH->{'_targetDNA_sgmentEnd'}  =$FinalOutWiseHASH->{'_targetDNA_sgmentEnd'}\t\$fastaFileSeqHead=$fastaFileSeqHead\n";		    	
		    }
		    #如果不是，则将参数中的这两个值填进去
		    else{ 
		    	warn "\$fastaFileSeqHead=$fastaFileSeqHead\tNot a Extract DNA coontig with \\d+..\\d+ pattern\nIt is OK!\n\n\n";#\$DNAtargetContain=$DNAtargetContain";   
		    	if (   (  ( defined ($DNAstart) ) && ($DNAstart=~m/^\d+$/)  )  &&  (  ( defined ($DNAstart) ) && ($DNAstart=~m/^\d+$/)  )   ) { 
		    		$upstrmNb=$FinalOutWiseHASH->{'_targetDNA_sgmentStart'}=$DNAstart; 
		    		$FinalOutWiseHASH->{'_targetDNA_sgmentEnd'}=$DNAend; 
		    	}
		      else { #如果参数中没有这两个值的化，则
		      	$upstrmNb=$FinalOutWiseHASH->{'_targetDNA_sgmentStart'}=$DNAstart=1; 
		    		$FinalOutWiseHASH->{'_targetDNA_sgmentEnd'}=$DNAend=FastaFileHandle::GetSeqLength($fastaFileString); 
		      }
		      
		      
		    }
		  }
		  else {die"\$tepNb=$tepNb\n\$eachWise=$eachWise\n\n"; 		    		  	#$upstrmNb=$FinalOutWiseHASH->{'_targetDNA_sgmentStart'}=$DNAstart; $FinalOutWiseHASH->{'_targetDNA_sgmentEnd'}=$DNAend;
		  }
		}
    elsif ($tepNb>0){
	    if ($eachWise=~m/
	                         #$1            #$2        #$3  ,$forwORrev
	                     ^\s*(\S+.*)\s+vs\s+(\S+.*)\s+\((forward|reverse)\)\s+\[\d+\]\s*\n
	                     \n
	                     (?:.*\S.*+\n){5}
	                     \n
	                    
	                     (?:\s{71}\n){5}
	                     \s{21}Alignment\s+\d+\s+Score\s+\S+\s+\(Bits\)\s+\n
	                     \n\n 
	                     #  $4 $alignMent, 具体的比对情况
	                     (  
	                       (?: 
	                            (?:\S.{18}\S) \s .{49} \s \n                     #
	                         (?:(?:\s{20}   ) \s .{49} \s \n){2}                 #
	                            (?:\S.{18}\S) \s .{49} \s \n                     #
	                         (?:(?:\s{20}   ) \s .{49} \s \n){2}                 #
	                         \n\n
	                       )+  
	                     )
	                     
	                      
	                      \/\/\n
	                      Bits\s+.*\n
	                      # $5 $wiseScore, wise得分
	                      (-?\d+(?:\.\d+)?)\s+# $6 $ProID, 蛋白的id   
	                                          (\S+)\s+(\d+)\s+(\d+)\s+\S+\s+(\d+)\s+(\d+)\s+\d+\s+\d+\n  
	                                                  #$7     $8            $9      $10                     #$7 $queryProtHd  $8 $queryProtTl   $9 $targetDNAhd  $10  $targetDNAtl  #这四个量是 query Protein 和 DNA 比对时，match上的 区域的头尾位置数。 前两个属于protein，后两个属于DNA。
	                      \/\/\n
	                      
	                      #$11 $geneInfom，大部分信息
	                      (
	                       (?:Gene\s+\d+\n
	                          Gene\s+(?:\d+)\s+(?:\d+)\s*\n
	                        (?:
	                         (?:\s+Exon\s+\d+\s+\d+\s+phase\s+\S+\n)+
	                        )
	                       )+
	                      )
	                      \/\/\n
	                      #$12 $transInfom，蛋白序列信息
	                      (
	                       (?:>\S+.*\n
	                         (?:
	                           (?:\S+\n)+
	                         )
	                       )+
	                      )
	                       \/\/\n
	                       
	                      #$13 $cDNAInform，cDNA的信息
	                      (
	                       (?:>\S+.*\n
	                         (?:
	                           (?:\S+\n)+
	                         )
	                       )+
	                      )
	                      \/\/\n
	                      
	                      #$14 $subSeqInformWithSubScore,  subsequence的信息，包括每个subsequence的得分
	                      (
	                       (?:Sequence\s+\S+\n
	                          subsequence\s+\S+\.\d+\s+\d+\s+\d+\n
	                          \n
	                          Sequence\s+\S+\.\d+\n
	                          CDS\n
	                          Start_not_found\n
	                          End_not_found\n
	                          CDS_predicted_by\s+genewise\s+\S+\n
	                         (?:
	                           (?:source_Exons\s+\d+\s+\d+\n)+
	                         )
	                          \n
	                       )+
	                      )
	                     \/\/\n              
	                     
	                      #$15 $GffImformation       
	                     (
                         (?:
	                         (?:\S+\t)+   \S+\n
	                       )+
	                     )
                       \/\/ 
	                     
	                     
	                    /mx)
	    {  
	    	my ($QueryProteinName, $targetDNAname, $forwORrev, $alignMent, $wiseScore, $ProID, $queryProtHd, $queryProtTl, $targetDNAhd, $targetDNAtl, $geneInfom, $transInfom, $cDNAInform, $subSeqInformWithSubScore, $GffImformation)  #$1 和 $6重复了。都是蛋白的 id。
	        =($1,                $2,             $3,         $4,         $5,         $6,     $7,           $8,           $9,           $10,          $11,        $12,         $13,         $14,                        $15           );
        #warn "\$QueryProteinName=$QueryProteinName,\$targetDNAname=$targetDNAname, \$forwORrev= $forwORrev\n\n\$alignMent=$alignMent\n";
        #下面这行，是以比对的方向是 forward还是 reverse来确定 后面代码中进行位置数计算时，是进行  加 还是 减的
	      my $forwORrevCHGnb=1; my $ZHENGorFU='+';
	      if ($forwORrev eq 'forward'){$forwORrevCHGnb=1; $ZHENGorFU='+';} elsif ($forwORrev eq 'reverse'){$forwORrevCHGnb=-1; $ZHENGorFU='-';} else {die"\$geneWisedbOutFile=$geneWisedbOutFile\n\n\$forwORrev=$forwORrev\n\n";}
	      
	      
	      
	      my $sameProtNameNumber=0;
	      if (  defined ( $FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName} )  ){
	      	$sameProtNameNumber=@{ $FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName} };   print "ssssssssssss \$sameProtNameNumber=$sameProtNameNumber\n";
	      }
	      $FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_proteinAlinmentDNAname'}    =$targetDNAname;   
        $FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_proteinAlinmentDirection'}  =$forwORrevCHGnb;
        $FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_proteinAlinmentTotalScore'} =$wiseScore;
        $FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_proteinAlinmentQueryPtHead'}=$queryProtHd;
        $FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_proteinAlinmentQueryPtTail'}=$queryProtTl;
        $FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_proteinAlinmentDNAtrgtHead'}=$targetDNAhd;
        $FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_proteinAlinmentDNAtrgtTail'}=$targetDNAtl;
        
	        
	     
	      my $pesudoMayBe=0;   my $AlignmentPesudoMARK=0; #这两个变量，都是用来判断，这个Protein VS DNA 的Alignment结果是否有 Pseudogene出现。 只是，前者用 得到几个gene结果（多于1个，则为有Pseudo）来判断，后者用 分析Alignment中是否有!等关键字符来判断
	      
	      my @aliMentsParts= 
	                     ($alignMent =~
	                      m/	                         
	                        ( 
	                             (?:\S.{18}\S) \s .{49} \s \n                     #
	                          (?:(?:\s{20}   ) \s .{49} \s \n){2}                 #
	                             (?:\S.{18}\S) \s .{49} \s \n                     #
	                          (?:(?:\s{20}   ) \s .{49} \s \n){2}                 #
	                          \n\n
	                        )+
	                       /gmx);

     
	      
	      my @geneInfomAr=
	                     ($geneInfom=~
	                      m/
	                        (Gene\s+\d+\n
	                           Gene\s+(?:\d+)\s+(?:\d+)\s*\n
	                          (?:
	                            (?:\s+Exon\s+\d+\s+\d+\s+phase\s+\S+\n)+
	                          )
	                        )
	                       /gmx);
	      my $ToTgeneNBs=@geneInfomAr;  
	      
	      my @trnsInfomAr=
	                     ($transInfom=~
	                      m/
	                        (>\S+.*\n
	                          (?:
	                            (?:[^>]\S*\n)+
	                          )
	                        )
	                       /gmx);
	      my $ToTtransNBs=@trnsInfomAr;
	      
	      my @cDNAInfomAr=
	                     ($cDNAInform=~
	                      m/
	                        (>\S+.*\n
	                         (?:
	                           (?:[^>]\S*\n)+
	                         )
	                        )
	                       /gmx);    
	      my $ToTcDNANBs=@cDNAInfomAr;
	                         
	      my @subSeqInformWithSubScoreAr=
	                                    ($subSeqInformWithSubScore=~
	                                     m/ 
	                                        (Sequence\s+\S+\n
	                                           subsequence\s+\S+\.\d+\s+\d+\s+\d+\n
	                                           \n
	                                           Sequence\s+\S+\.\d+\n
	                                           CDS\n
	                                           Start_not_found\n
	                                           End_not_found\n
	                                           CDS_predicted_by\s+genewise\s+\S+\n
	                                          (?:
	                                            (?:source_Exons\s+\d+\s+\d+\n)+
	                                          )
	                                           \n
	                                        )
	                                      /gmx);
	                     
	      my $ToTsubSeqNBs=@subSeqInformWithSubScoreAr;
	      
	      #################################################解析GFF信息################################################################################################################################################
	      my @GffImformationAr=(
	                          $GffImformation=~m/ 
	                                             (?:((?:\S+\t)+   \S+)\n	)
	                                           /gmx
	                        );
	      my $ToTgffLineNBs=@GffImformationAr;
	      my $gffInformHash; my $subGenNB=0; my $cdsNumbe=0; my $intronNb=0; my $subSegNB=0;
	      my $notDividedGffHash; my $notDividedGffHash_idx=0;
	      foreach my $eachGffLine(@GffImformationAr){  #print "\$eachGffLine=$eachGffLine\n";
	      	my ($DNActgID, $ProgName, $segmType, $startPos, $endPosti, $scoreNum, $stranZoF, $GffFrame, $prdctIfm)=split "\t",$eachGffLine;
	      	my $tpHASH; 
	      	$tpHASH->{'0DNActgID'}=$DNActgID; 	      	    $tpHASH->{'1ProgName'}=$ProgName;	      	    $tpHASH->{'2segmType'}=$segmType;	      	    $tpHASH->{'3startPos'}=$startPos+$upstrmNb-1;	      	    $tpHASH->{'4endPosti'}=$endPosti+$upstrmNb-1;
	      	$tpHASH->{'5scoreNum'}=$scoreNum;	      	      $tpHASH->{'6stranZoF'}=$stranZoF;	      	    $tpHASH->{'7GffFrame'}=$GffFrame;	      	    $tpHASH->{'8prdctIfm'}=$prdctIfm;
	      	
	      	if ($segmType eq 'cds'){
	      		$gffInformHash->{'Ordered_Match_Hash'}->{$subGenNB}->{'Ordered_CDS_Hash'}->{$cdsNumbe}=$tpHASH;
	      		$gffInformHash->{'Ordered_Match_Hash'}->{$subGenNB}->{'Ordered_subSegm_Hash'}->{$subSegNB}->{'hashType'}='Ordered_CDS_Hash';
	      		$gffInformHash->{'Ordered_Match_Hash'}->{$subGenNB}->{'Ordered_subSegm_Hash'}->{$subSegNB}->{'TypeOrde'}=$cdsNumbe;
	      		$cdsNumbe++; $subSegNB++;
	      	}
	      	elsif($segmType eq 'intron'){
	      		$gffInformHash->{'Ordered_Match_Hash'}->{$subGenNB}->{'Ordered_Intron_hash'}->{$intronNb}=$tpHASH;
	      		$gffInformHash->{'Ordered_Match_Hash'}->{$subGenNB}->{'Ordered_subSegm_Hash'}->{$subSegNB}->{'hashType'}='Ordered_Intron_hash';
	      		$gffInformHash->{'Ordered_Match_Hash'}->{$subGenNB}->{'Ordered_subSegm_Hash'}->{$subSegNB}->{'TypeOrde'}=$intronNb;
	      		$intronNb++; $subSegNB++;
	      	}
	      	
	      	elsif ($segmType eq 'match'){
	      		if ($notDividedGffHash_idx==0){ $subGenNB=0; }
	      		else                          { $subGenNB++;}
	      	    	      	  
	      	  $gffInformHash->{'Ordered_Match_Hash'}->{$subGenNB}->{'matchInform'}=$tpHASH;
	      	  $gffInformHash->{'Ordered_Match_Hash'}->{$subGenNB}->{'Ordered_subSegm_Hash'}->{$subSegNB}->{'hashType'}='matchInform';
	      		$gffInformHash->{'Ordered_Match_Hash'}->{$subGenNB}->{'Ordered_subSegm_Hash'}->{$subSegNB}->{'TypeOrde'}=$intronNb;
	      	   $cdsNumbe=0; $intronNb=0; $subSegNB=0;
	      	  
	      	}
	      	else{
	      		die "\n\n\nIn package GeneWiseHandle,\nIn sub NewGenewisedbPraser,\n\n\$segmType=$segmType\n\n\n";
	      	}
	      	
	      	$notDividedGffHash->{$notDividedGffHash_idx}=$tpHASH; $notDividedGffHash_idx++;
	      	
	      }
	      $FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_GffImformation'}=$GffImformation;
	      $FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_GffImformationHash'}=$gffInformHash;
	      $FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_notDividedGffHash'}=$notDividedGffHash;
	      #################################################解析GFF信息################################################################################################################################################
	      
	      print "\$ToTgeneNBs=$ToTgeneNBs\t\$ToTtransNBs=$ToTtransNBs\t\$ToTcDNANBs=$ToTcDNANBs\n";   warn "\$ToTgeneNBs=$ToTgeneNBs\t\$ToTtransNBs=$ToTtransNBs\t\$ToTcDNANBs=$ToTcDNANBs\n";   
	      
	      if ( ($ToTgeneNBs == $ToTtransNBs) && ($ToTcDNANBs == $ToTtransNBs) && ($ToTsubSeqNBs == $ToTtransNBs) && ($ToTcDNANBs >= 1) ){}  #这四个数应该相等
	      else {die "\nWrongWrongWrong!!!\n\$eachWise\n=$eachWise\n\n\n\nThe three val following should be abslute euqal and not smaller than 1\n\$ToTgeneNBs=$ToTgeneNBs\t\$ToTtransNBs=$ToTtransNBs\t\$ToTcDNANBs=$ToTcDNANBs\t\$ToTsubSeqNBs=$ToTsubSeqNBs\n\n\$geneInfom=$geneInfom\n\$transInfom=$transInfom\n\$cDNAInform=$cDNAInform\n\$subSeqInformWithSubScore=$subSeqInformWithSubScore\n\n";}
	      
	      if ($ToTcDNANBs > 1) {$pesudoMayBe=1;}  #这里是判断是否是 pseudo gene的地方
	        
	      # 下面的循环中引入了 $upstrmNb   这个变量
	      for (my $bgIdx=0; $bgIdx<$ToTgeneNBs; $bgIdx++){  #对每个 genewisedb比对得到的 gene进行 一次循环。这里的一个gene通常是被 pseudogene 分开的 gene部分，  
	        my $geneStart; my $geneEnd; my $exonInform; my $pepSeq; my $cDNASeq;  my $subSeqScore;
	        if ($geneInfomAr[$bgIdx]=~m/
	                                    Gene\s+\d+\n
	                                    Gene\s+(\d+)\s+(\d+)\s*\n 
	                                    (	 (?:\s+Exon\s+\d+\s+\d+\s+phase\s+\S+\n)+  )
	                                    /mx ){ 
	        	($geneStart, $geneEnd, $exonInform) =($1,  $2,  $3); 
	        	$geneStart=$upstrmNb+$geneStart-1; $geneEnd=$upstrmNb+$geneEnd-1;  #加上了被去掉的部分的长度。这里用了 $upstrmNb 这个外来变量	        
	        	
   		      $FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneHeadPos'}=$geneStart;   
            $FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneTailPos'}=$geneEnd; 
            
            $FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_TGAinformHash'}->{'_TGAtransForUexists'}=0; #每个geneweie预测的基因，在最初录入 总$FinalOutWiseHASH时，都其每个 gene是否有TGA-U初始化为0.
	        		
	        }  
	        else {die "\n\n\nWrong Here!!!\n\n\$geneInfomAr[$bgIdx]=\n$geneInfomAr[$bgIdx]\n";}  # 这里获得了gene的位置信息  #这里 获得了 my $geneStart; my $geneEnd; my $exonInform 这三个变量的值。 $exonInform 需要进一步解析
	        
	        if ($trnsInfomAr[$bgIdx]=~m/
	                                    >\S+.*\n
	                                    ( (?:\S+\n)+ )
	                                    /mx )
	        {
	          $pepSeq = $1;      $pepSeq=~s/\s//g;
	          $FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGenePepSeq'}=$pepSeq;  
	        }  else {die "\n\n\nWrong Here!!!\n\n\$trnsInfomAr[$bgIdx]=\n$trnsInfomAr[$bgIdx]\n";}   #这里获得了 $pepSeq 的值
	         
	        if ($cDNAInfomAr[$bgIdx]=~m/
	                                    >\S+.*\n
	                                    ( (?:\S+\n)+ )
	                                    /mx )
	        { 
	        	$cDNASeq = $1;    $cDNASeq=~s/\s//g;
	        	$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneDNAseq'}=$cDNASeq; 
	        }  else {die "\n\n\nWrong Here!!!\n\n\$cDNAInfomAr[$bgIdx]=\n$cDNAInfomAr[$bgIdx]\n";}  #这里获得了 $cDNASeq 的值
	        
	        if ($subSeqInformWithSubScoreAr[$bgIdx]=~m/CDS_predicted_by genewise\s+(\S+)\n/){	
	        	$subSeqScore=$1;  
	        	$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneSubSeqScore'}=$subSeqScore;
	        	#下面这个If语句，是用来判断 没有U的 gene中，得分最高的 gene的 关键字信息
	        	if ($BestScoreGene_Score<=$subSeqScore){                                                           print "$geneWisedbOutFile If       \$BestScoreGene_Score=$BestScoreGene_Score <= $subSeqScore=\$subSeqScore\n";
	        		$BestScoreGene_Score=$subSeqScore;                                                               print "$geneWisedbOutFile Yes Then \$BestScoreGene_Score=$BestScoreGene_Score =  $subSeqScore=\$subSeqScore\n";
	        		$FinalOutWiseHASH->{'_BestScoreAllGeneHASH'}->{'_BestScoreQueryProteinID'}=$QueryProteinName;    print "\$FinalOutWiseHASH->{'_BestScoreAllGeneHASH'}->{'_BestScoreQueryProteinID'}=\$QueryProteinName=$QueryProteinName\n";
	        	  $FinalOutWiseHASH->{'_BestScoreAllGeneHASH'}->{'_BestScoreQuProtNUMBER'}=$sameProtNameNumber;    print "\$FinalOutWiseHASH->{'_BestScoreAllGeneHASH'}->{'_BestScoreQuProtNUMBER'}=\$sameProtNameNumber=$sameProtNameNumber\n";
	        	  $FinalOutWiseHASH->{'_BestScoreAllGeneHASH'}->{'_BestScoreGeneIdx'}=$bgIdx;	        	           print "\$FinalOutWiseHASH->{'_BestScoreAllGeneHASH'}->{'_BestScoreGeneIdx'}=\$bgIdx=$bgIdx\n\n";
	        	}  
	        }       #这里获得了 $subSeqScore 的值 ；；也通过比较，最终获得了 所有得分中 最高分的单个Gene的score
	        else {die "\n\n\nWrong Here!!!\n\n\$subSeqInformWithSubScoreAr[$bgIdx]=\n$subSeqInformWithSubScoreAr[$bgIdx]\n";}
	        
	        my $direction='+'; if ($geneStart > $geneEnd) { $direction='-';}
	        #print "\cl\cl\Yes\$eachWise=n$eachWise\n\cl$1\n$2\n$3\n$4\n$5\n\$direction=$direction\n\n";                 #这里再次获得 方向 $direction 的值
	        
	        my @exonsAr=split ("\n", $exonInform);    #将 $exonInform中的内容，变成数组
	        my $tgaInclude=0;  my $cDNAbeforTGANb=0; my $cds_strucNBs=''; my $geneListAr;  # 预设左边变量，$tgaInclude 判断TGA是否被包括，$cDNAbeforTGANb 用于计算 TGA前的DNA数目，$cds_strucNBs是输出的编码区位置数列，        
	        
	        my $exonADDlength=0; #建立一个变量，用来装载，从第一个外显子的第一个碱基开始计数，直到当前这个外显子结尾，共有多少个碱基 这样一个数量。开始初始化为0
	        
	        my $exonIdx=0;    my $changeDNAposIdx=0;
	        foreach my $exonLine (@exonsAr){
	      	  if ($exonLine=~m/\s+Exon\s+(\d+)\s+(\d+)\s+phase\s+(\d+)/){
	      	  	my ($orgNbA, $orgNbB)=($1, $2);
	      	  	$orgNbA=$orgNbA+$upstrmNb-1; $orgNbB=$orgNbB+$upstrmNb-1;  #加上了被去掉的部分的长度  $upstrmNb
	      	  	$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneExonArray'}->[$exonIdx]->{'_ExonHead'}=$orgNbA;
	      	  	$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneExonArray'}->[$exonIdx]->{'_ExonTail'}=$orgNbB;
	      	  	$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneExonArray'}->[$exonIdx]->{'_Exonphase'}=$3;
	      	  	my $ExonInGeneHead=$exonADDlength+1;
	      	  	my $ExonInGeneTail=$exonADDlength+abs($orgNbB-$orgNbA)+1;
	      	  	$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneExonArray'}->[$exonIdx]->{'_ExonInGeneHead'}=$ExonInGeneHead;
	      	  	$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneExonArray'}->[$exonIdx]->{'_ExonInGeneTail'}=$ExonInGeneTail;
	      	  	
	      	  	if (  ( defined ($DNAchangedHash) )&& ( ref($DNAchangedHash) eq 'HASH' )  ){
	      	  	  foreach my $eachPOSnb (   sort {$a <=> $b} (  keys %{ $DNAchangedHash }  )   ){
	      	  	  	my $orgChangePos=$DNAchangedHash->{$eachPOSnb}->{'1_postion'};
	      	  	    if ( SeqSegmentsTools::a_between_b_c($orgChangePos,$orgNbA,$orgNbB) ){
	      	  	    	my $orgChangePosZhengFu=$DNAchangedHash->{$eachPOSnb}->{'0_directi'};
	      	  	    	if ($ZHENGorFU eq $ZHENGorFU){
	      	  	    	  my $coresBondPosInGENE=SeqSegmentsTools::GetTheCorespondPosOFDIfferentSegement($orgNbA,$orgNbB,$ZHENGorFU,$orgChangePos,$ExonInGeneHead,$ExonInGeneTail,'+');
	      	  	    	  my $yuShuDivedBy3=$coresBondPosInGENE%3; #因为需要改的都是TGA的最后一位A，转化为C，所以这个位置应该是在编码的序列中的3号frame上，可以整除3
	      	  	    	  if ($yuShuDivedBy3==0){
	      	  	    	  	
	      	  	    	  	$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_changeDNAPosHash'}->{$changeDNAposIdx}->{'0_directi'}='+';
	      	        	  	$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_changeDNAPosHash'}->{$changeDNAposIdx}->{'1_postion'}=$coresBondPosInGENE;
	                	  	$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_changeDNAPosHash'}->{$changeDNAposIdx}->{'2_fromCha'}=$DNAchangedHash->{$eachPOSnb}->{'3_intoCha'};
	          	    	  	$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_changeDNAPosHash'}->{$changeDNAposIdx}->{'3_intoCha'}=$DNAchangedHash->{$eachPOSnb}->{'2_fromCha'};
	          	    	  	
	          	    	  	$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_changeDnaPEPhash'}->{$coresBondPosInGENE-2}->{'AminoAcidResidue'}='U';
	          	    	  	
	      	  	    	  	$changeDNAposIdx++;
	      	  	    	  }
	      	  	    	}
	      	  	    }
	      	  	  }
	      	  	}
	      	  	
	      	  	
	      	  	
	      	  	
	      	  	
	      	  	
	      	  	$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneExonArray'}->[$exonIdx]->{'_ExonHeadAAcodenPos'}=($FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneExonArray'}->[$exonIdx]->{'_ExonInGeneHead'})%3;
	      	  	if ($FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneExonArray'}->[$exonIdx]->{'_ExonHeadAAcodenPos'} ==0){
	      	  		$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneExonArray'}->[$exonIdx]->{'_ExonHeadAApos'} = int ( ($FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneExonArray'}->[$exonIdx]->{'_ExonInGeneHead'})/3 );
	      	  	}
	      	  	else {
	      	  		$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneExonArray'}->[$exonIdx]->{'_ExonHeadAApos'} = int ( ($FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneExonArray'}->[$exonIdx]->{'_ExonInGeneHead'})/3 ) + 1;
	      	  	}
	      	  	
	      	  	$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneExonArray'}->[$exonIdx]->{'_ExonTailAACodenPos'}=($FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneExonArray'}->[$exonIdx]->{'_ExonInGeneTail'})%3;
	      	  	if ($FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneExonArray'}->[$exonIdx]->{'_ExonTailAACodenPos'} ==0){
	      	  		$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneExonArray'}->[$exonIdx]->{'_ExonTailAApos'} = int ( ($FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneExonArray'}->[$exonIdx]->{'_ExonInGeneTail'})/3 );
	      	  	}
	      	  	else {
	      	  		$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneExonArray'}->[$exonIdx]->{'_ExonTailAApos'} = int ( ($FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneExonArray'}->[$exonIdx]->{'_ExonInGeneTail'})/3 ) + 1;
	      	  	}      	  	
      	  	
	      	  	
	      	  	$exonADDlength+=abs($orgNbB-$orgNbA)+1;
	      	  	
	      	  	

	      	  	
	      	  	if  (    (  ( defined ($stopPos) ) && ($stopPos=~m/^\d+$/)  )&& (   ( ($orgNbA<=$stopPos) &&($stopPos<=$orgNbB) ) || ( ($orgNbB<=$stopPos) &&($stopPos<=$orgNbA) )   )    ){  #这里是判断 TGA码的位置 是否出现在该exon中。此条件成立则表明 TGA码 在 该exon中
	      	  	   $cDNAbeforTGANb+=abs ($stopPos-$orgNbA);  #计算TGA前面的 碱基数目
	      	  	   $tgaInclude=1;  #print "\$exonLine =$exonLine \n\$cDNAbeforTGANb=$cDNAbeforTGANb\n\$cDNAbeforTGANb=$cDNAbeforTGANb+=abs ($stopPos-$orgNbA)\n";
	      	  	   $FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_TGAinformHash'}->{'_TGAtransForUexists'}=1;
	      	  	   $FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_TGAinformHash'}->{'_TGAincludeExonNB'}  =$exonIdx;
	      	  	   $FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_TGAinformHash'}->{'_TGAPosition'}       =$stopPos;
	      	  	   $FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_TGAinformHash'}->{'_TGAEndPs'}          =$stopPos+2*$forwORrevCHGnb;
	      	  	   #下面这个If语句，是用来判断  有U的 gene中，得分最高的 gene的 关键字信息
	      	  	   if ($BestScoreGeneWithU_Score <= $FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneSubSeqScore'} ){       print "$geneWisedbOutFile If       \$BestScoreGeneWithU_Score=$BestScoreGeneWithU_Score <= $FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneSubSeqScore'}=\$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneSubSeqScore'}\n";                                      
	        		     $BestScoreGeneWithU_Score = $FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneSubSeqScore'};            print "$geneWisedbOutFile Yes Then \$BestScoreGeneWithU_Score=$BestScoreGeneWithU_Score  = $FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneSubSeqScore'}=\$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneSubSeqScore'}\n";                                   
	        		     $FinalOutWiseHASH->{'_BestScoreWithUhash'}->{'_BestScoreQueryProteinID'}=$QueryProteinName;                                                                                                print "\$FinalOutWiseHASH->{'_BestScoreWithUhash'}->{'_BestScoreQueryProteinID'}=\$QueryProteinName=$QueryProteinName\n";  
                   $FinalOutWiseHASH->{'_BestScoreWithUhash'}->{'_BestScoreQuProtNUMBER'}=$sameProtNameNumber;                                                                                                print "\$FinalOutWiseHASH->{'_BestScoreWithUhash'}->{'_BestScoreQuProtNUMBER'}=\$sameProtNameNumber=$sameProtNameNumber\n";
	        	       $FinalOutWiseHASH->{'_BestScoreWithUhash'}->{'_BestScoreGeneIdx'}=$bgIdx;	        	                                                                                                      print "\$FinalOutWiseHASH->{'_BestScoreWithUhash'}->{'_BestScoreGeneIdx'}=\$bgIdx=$bgIdx\n\n";                               
	        	     } 
	        	     
	        	     #下面这个If语句，是用来判断  有U的 Protein query 预测的 gene 的组中，得分最高的 gene的总分（如有pseudo则有可能是多个gene得分之和） 的关键字信息
	               if ($BestTotalScoreWithU <= $FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_proteinAlinmentTotalScore'} ){                                      print "$geneWisedbOutFile If       \$BestTotalScoreWithU=$BestTotalScoreWithU <= $FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_proteinAlinmentTotalScore'}=\$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_proteinAlinmentTotalScore'}\n";
	                 $BestTotalScoreWithU =    $FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_proteinAlinmentTotalScore'};                                        print "$geneWisedbOutFile Yes Then \$BestTotalScoreWithU=$BestTotalScoreWithU  = $FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_proteinAlinmentTotalScore'}=\$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_proteinAlinmentTotalScore'}\n";
	                 $FinalOutWiseHASH->{'_BestTotalScoreWithUhash'}->{'_BestScoreQueryProteinID'}=$QueryProteinName;                                                                    print "\$FinalOutWiseHASH->{'_BestTotalScoreWithUhash'}->{'_BestScoreQueryProteinID'}=\$QueryProteinName=$QueryProteinName\n";                                                                                                                                                                                                             
                   $FinalOutWiseHASH->{'_BestTotalScoreWithUhash'}->{'_BestScoreQuProtNUMBER'}=$sameProtNameNumber;                                                                    print "\$FinalOutWiseHASH->{'_BestTotalScoreWithUhash'}->{'_BestScoreQuProtNUMBER'}=\$sameProtNameNumber=$sameProtNameNumber\n";
	                 $FinalOutWiseHASH->{'_BestTotalScoreWithUhash'}->{'_BestScoreGeneIdx'}=$bgIdx;	        	                                                                           print "\$FinalOutWiseHASH->{'_BestTotalScoreWithUhash'}->{'_BestScoreGeneIdx'}=\$bgIdx=$bgIdx\n\n";                                                                                                                                                                                                                                          
	               } 
	      	  	}
	      	  	else {         #这里是判断 TGA码的位置 是否出现在该exon中。此条件    表明 TGA码   不在 该exon中 
	      	  		if ($tgaInclude == 0) { 
	      	  			$cDNAbeforTGANb+=abs($orgNbB-$orgNbA)+1; #print "\$cDNAbeforTGANb=$cDNAbeforTGANb+=abs($orgNbB-$orgNbA+1)\n";   #这个 应该也是没有用的语句了
	      	  			#下面这个If语句，是用来判断 没有U的 gene中，得分最高的 gene的 关键字信息
	      	  			if ($BestScoreGeneWithOutU_Score <= $FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneSubSeqScore'} ){   print "$geneWisedbOutFile If       \$BestScoreGeneWithOutU_Score=$BestScoreGeneWithOutU_Score <= $FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneSubSeqScore'}=\$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneSubSeqScore'}\n";
	        		      $BestScoreGeneWithOutU_Score =    $FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneSubSeqScore'};     print "$geneWisedbOutFile Yes Then \$BestScoreGeneWithOutU_Score=$BestScoreGeneWithOutU_Score  = $FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneSubSeqScore'}=\$FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneSubSeqScore'}\n";
	        		      $FinalOutWiseHASH->{'_BestScoreWithOutUhash'}->{'_BestScoreQueryProteinID'}=$QueryProteinName;                                                                     print "\$FinalOutWiseHASH->{'_BestScoreWithOutUhash'}->{'_BestScoreQueryProteinID'}=\$QueryProteinName=$QueryProteinName\n";                                                                                                                                          
                    $FinalOutWiseHASH->{'_BestScoreWithOutUhash'}->{'_BestScoreQuProtNUMBER'}=$sameProtNameNumber;                                                                     print "\$FinalOutWiseHASH->{'_BestScoreWithOutUhash'}->{'_BestScoreQuProtNUMBER'}=\$sameProtNameNumber=$sameProtNameNumber\n";
	        	        $FinalOutWiseHASH->{'_BestScoreWithOutUhash'}->{'_BestScoreGeneIdx'}=$bgIdx;                                                                                       print "\$FinalOutWiseHASH->{'_BestScoreWithOutUhash'}->{'_BestScoreGeneIdx'}=\$bgIdx=$bgIdx\n\n";                                                                                                                                                                       
	        	       
	        	      } 
	      	  		   
	      	  		} 
	      	  	}  #如果tga不在这个Exon中，则 直接将所有Exon长度加到TGA前面的碱基数目中
	      	  	
	      	  	
	      	  	
	      	  	
	      	    $exonIdx++;
	      	  }
	      	  else {die "Wrong: \$exonLine=$exonLine\n\$eachWise=$eachWise\n";}
	        }  #以上的循环， 完成了 my $tgaInclude=0;  my $cDNAbeforTGANb=0; my $cds_strucNBs=''; 这三个变量的赋值过程。
	        
	        
	        if (  defined ( $FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_changeDNAPosHash'} )  ){
	        	if (  ref( $FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_changeDNAPosHash'}  ) eq 'HASH' ){
	        		my $changedBackSeqWithU_DNA=FastaFileHandle::CHangeSpecificSequenceInAFastaSeq( 
	        		                                                                             $FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_EachGeneDNAseq'},
	        		                                                                             $FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_changeDNAPosHash'} 
	        		                                                                           );
              $FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'changedBackSeqWithU_DNA'}=$changedBackSeqWithU_DNA;
              	        		                                                                         
              my $changedBackSeqWithU_PEP=FastaFileHandle::TranslateSeqWithTGAcodonHash(  
                                                                                          $changedBackSeqWithU_DNA,
                                                                                          $FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_changeDnaPEPhash'}
                                                                                       );
              $FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'changedBackSeqWithU_PEP'}=$changedBackSeqWithU_PEP;                                                                                 		                                                                           
	        		                                                                           
	        	}
	        }
	         
	        
	        #$cds_strucNBs=~s/,$//;  print "\$cds_strucNBs=$cds_strucNBs\n";  #去掉$cds_strucNBs 末尾的逗号，并打印
	        
	        
	        
	        
	        #$geneListAr= &transList2Array ( $cds_strucNBs, "CDS", $geneListAr    );
	        
	        if ($tgaInclude == 1){
  	        my $vrTGApos=$stopPos; if ($direction eq '-'){ $vrTGApos=-1*$stopPos;}  #print "\$direction=$direction\n\$vrTGApos=$vrTGApos\n"; 
	          #$geneListAr= &transList2Array ( $vrTGApos, "TGA", $geneListAr    );
	        }
	        #my $ImageList=&Array2List ($geneListAr); print "\$ImageList=$ImageList\n"; 
          
          my $pepSeqNoSpace=$pepSeq; $pepSeqNoSpace=~s/\s//g; my $pepLength=length ($pepSeqNoSpace);  #print "\n( ($tgaInclude == 1) && ( ($cDNAbeforTGANb%3) == 0 ) )\n\n";
	        if ( ($tgaInclude == 1) && ( ($cDNAbeforTGANb%3) == 0 ) ) { 
	      	  #print "\$orgThreeWd=$orgThreeWd\n\$orgResidue=$orgResidue\n";
	      	  my $UaaPos=($cDNAbeforTGANb/3+1);
	      	  my $ChgBakPep =&printHeadWordTail_Chg2Org ($pepSeq, ($cDNAbeforTGANb/3+1),1,$orgResidue); #print "\$ChgBakPep=$ChgBakPep\n";
	      	  my $ChgBakCDNA=&printHeadWordTail_Chg2Org ($cDNASeq,($cDNAbeforTGANb+1),  3,$orgThreeWd); #print "\$ChgBakCDNA=$ChgBakCDNA\n";
	      	  $FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_TGAinformHash'}->{'_UaaPostionNum'}  =$UaaPos;       #记录wise结果中，各个有U的 gene的 U在pep中的位置
	      	  $FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_TGAinformHash'}->{'_TGAchangedPEP'}  =$ChgBakPep;    #记录wise结果中，各个有U的 gene的 pep序列，这里原来被转换为C的U被重新转换回来
	      	  $FinalOutWiseHASH->{'_genewiseResultHashAd'}->{$QueryProteinName}->[$sameProtNameNumber]->{'_EachGeneArrayAd'}->[$bgIdx]->{'_TGAinformHash'}->{'_TGAchangedDNA'}  =$ChgBakCDNA;   #记录wise结果中，各个有U的 gene的 DNA序列，这里原来被转换为TGC的TGA被重新转换回来
	      	  
	      	  #$arrForSort->[$aFsIdx]=[$pepLength, $wiseScore, $ProID, $ChgBakPep, $ChgBakCDNA ,$ImageList ];  $aFsIdx++;  
	        }  #                       0            1           2       3           4            5 
	        #else {
	      	  #$noUarrForSort->[$NUaFsIdx]=[$pepLength, $wiseScore, $ProID, $pepSeq, $cDNASeq ,$ImageList ];  $NUaFsIdx++;  
	        #}   
 
	        #}
	      
	        if ($pesudoMayBe == $AlignmentPesudoMARK) {} else {warn"\$eachWise=$eachWise\n\n\$pesudoMayBe=$pesudoMayBe ?==?$AlignmentPesudoMARK=\$AlignmentPesudoMARK.\n\nMaybe the alignment analysis step was canceled, so these 2 value is not equal. or Some thing is wrong!!!";}
	      } 
	    }
	    else {
	    	die "Wrong: \$eachWise=$eachWise\n\n\nNot a good pattarn\n\n";
	    }
	  }
	  $tepNb++;
	} 
	$/=$tempMk;
	
	$FinalOutWiseHASH=&AddPsdGeneInformForWiseHash($FinalOutWiseHASH);  #Add pseudoInformation
	
	return $FinalOutWiseHASH;
	                       
}


sub printHeadWordTail_Chg2Org{
  my ($seq, $stNb, $WdLength, $opgWd)=@_;           #print "\n\n$WdLength from $stNb\n"; 
  $seq=~s/\s//g;
  my $seqHead=substr ($seq, 0, ($stNb-1));          #print "\n$seqHead\n";
  my $WordsHr=substr ($seq, ($stNb-1), $WdLength);  #print "$WordsHr\n";
  my $seqTail=substr ($seq, ($stNb+$WdLength-1), ); #print "$seqTail\n\n";
  my $outVal="$seqHead$opgWd$seqTail";
  return $outVal;
}

sub PrintPngForGene{        #pmAble#     #打印基因结构
  my ($SecCys, $ProID, $InIN_BothSecCysWiseHash, $InIN_BothSecCysHash, $inKeyName, $pngFile, $SECISshowOrNot)=@_;
  open (OUT, ">$pngFile") or die "2 Cannot create \$pngFile=$pngFile : $!\n";
  my $trackHash; 
  my $panelStart=99999999999999;  my $panelEnd=-99999999999999;     
  my $wiseRstHASH=$InIN_BothSecCysWiseHash->{$SecCys}->{$ProID};
  my $BestGenesArray=$wiseRstHASH
                                              ->{'_genewiseResultHashAd'}
                                              ->{ $wiseRstHASH->{'_BestScoreWithUhash'}->{'_BestScoreQueryProteinID'} }
                                              ->[ $wiseRstHASH->{'_BestScoreWithUhash'}->{'_BestScoreQuProtNUMBER'} ]
                                              ->{'_EachGeneArrayAd'};
  
  my $bigGeneList;       
  
  if (ref($BestGenesArray) eq 'ARRAY'){
    for (my $i=0; $i< @{ $BestGenesArray }; $i++){ warn "\$i=$i\n\n";
      
      my $InputCDS_AA_sementArray=$wiseRstHASH
                                                  ->{'_genewiseResultHashAd'}
                                                  ->{ $wiseRstHASH->{'_BestScoreWithUhash'}->{'_BestScoreQueryProteinID'} }
                                                  ->[ $wiseRstHASH->{'_BestScoreWithUhash'}->{'_BestScoreQuProtNUMBER'} ]
                                                  ->{'_EachGeneArrayAd'}
                                                  ->[ $i ]
                                                  ->{'_EachGeneExonArray'};
      my $strand='+'; if ($InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_SecPosFrame'}<0){$strand='-';}
      
       warn "\naaaaaaaa1111 \$InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_CoreCODONpos'}=",$InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_CoreCODONpos'},"\n\n";      
      my $ExonPosList=&AddTGACoreImform_20170307( $InputCDS_AA_sementArray, '' , $strand); warn "\naaaaaaaa \$InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_CoreCODONpos'}=",$InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_CoreCODONpos'},"\n\n";
      #print "\cl\clAAAAAAAAAA111111 PrintSubArrayHash::print_all_sub_array (\$ExonPosList)\n\n"; PrintSubArrayHash::print_all_sub_array ($ExonPosList); print "\n\n\cl\cl\n\n";
      
      $bigGeneList->{$i}=$ExonPosList;
      print "\cl\clAAAAAAAAAA333333\ PrintSubArrayHash::print_all_sub_array (\$ExonPosList)\n\n"; PrintSubArrayHash::print_all_sub_array ($ExonPosList); print "\n\n\cl\cl\n\n";
       
           
      
    
      
    }
  }                                           
   
   my $cosPosHere=$InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_CoreCODONpos'};
   my $cosPosFram=$InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_SecPosFrame'};
   my $cosPosStrand=1; if ( ($cosPosFram<0) || ($cosPosFram eq '-') || ($cosPosFram eq '-1') || ($cosPosFram eq '-2') || ($cosPosFram eq '-3') ){$cosPosStrand=-1;}
   my ($bigExonPosList,$psdPosList)=@{ &combaineForPseudoGenes_20170307($bigGeneList, $cosPosHere, $cosPosStrand) }; 
  
  #$psdPosList=&estAddCorPos($psdPosList, $InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_CoreCODONpos'});
  if (defined($InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_bestSECISstart'})){
  	if ( $InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_bestSECISstart'}=~m/\d+/ ){  print "AAAAAAAAAA222222211111\$InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_bestSECISstart'}=$InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_bestSECISstart'}\n";
  		$psdPosList=&AddSECISinform( $psdPosList, $InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_bestSECISstart'}, $InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_bestSECISend'}, $InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_bestSECISstrand'}, $InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_BestSECISinform'} );
  	  print "\cl\clAAAAAAAAAA22222222222 PrintSubArrayHash::print_all_sub_array (\$bigExonPosList)\n\n"; PrintSubArrayHash::print_all_sub_array ($bigExonPosList); print "\n\n\cl\cl\n\n";
  	}
  }
  
  
  $panelStart=&GetSmallestStart20170214($bigExonPosList,$panelStart);  $panelStart=&GetSmallestStart20170214($psdPosList,$panelStart);
  $panelEnd  =&GetBiggestEnd20170214   ($bigExonPosList,$panelEnd);    $panelEnd  =&GetBiggestEnd20170214   ($psdPosList,$panelEnd);                                          
  
  my $SECISInform=$InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_BestSECISinform'}; $SECISInform=ExcelHandle::GetLableOfHpLk($SECISInform);
  my $CorePosCode=$InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_CoreCODONtype'};
   
  if ($SECISshowOrNot==0){$SECISInform='SECIS';}
   
  $trackHash->{ 0 }->{'_SegmentList'}   =$psdPosList;
  $trackHash->{ 0 }->{'_SegmentKeyName'}='';
  $trackHash->{ 0 }->{'_SegmentProID'}  =$ProID;  
  $trackHash->{ 0 }->{'_dataType'}      ='none';
  $trackHash->{ 0 }->{'_tkBakeCol'}     ='white';#'ghostwhite'; #'white';#'whitesmoke';#'navajowhite'; #'whitesmoke';
  $trackHash->{ 0 }->{'_CorePosCode'}   =$CorePosCode; 
  $trackHash->{ 0 }->{'_SECISinform'}   =$SECISInform; 
   
    
  $trackHash->{ 1 }->{'_SegmentList'}   =$bigExonPosList;
  $trackHash->{ 1 }->{'_SegmentKeyName'}=$inKeyName;
  $trackHash->{ 1 }->{'_SegmentProID'}  =$ProID; 
  $trackHash->{ 1 }->{'_dataType'}      ='hat'; #solid quill dushed  
  $trackHash->{ 1 }->{'_tkBakeCol'}     ='white';#'ghostwhite'; #'antiquewhite';
  $trackHash->{ 1 }->{'_CorePosCode'}   =$CorePosCode; 
  $trackHash->{ 1 }->{'_SECISinform'}   =$SECISInform; 

  my $panelLength=$panelEnd-$panelStart+1;   print "AAAAAAAAAA\$panelStart=$panelStart\t\t\$panelEnd=$panelEnd\t\t\$panelLength=$panelLength\n\n";
  $panelStart=$panelStart- int($panelLength/10); if ($panelStart<1){$panelStart=1;}
  $panelEnd=  $panelEnd  + int($panelLength/10); $panelLength=$panelEnd-$panelStart+1;
  
  my $contigNm=$InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_GenomeORestContigID'}; $contigNm=~s/^(\S+)\s*.*$/$1/;  #-key     => $contigNm,
  
  my $panel = Bio::Graphics::Panel->new(
                                             -length    => $panelLength,   #图的全长宽度能表示的碱基数量总长度
                                             -width     => 1000,                #图的实际宽度
                                             -key_style => 'left',
                                             -start     => $panelStart,
                                             -pad_left  => 100,
                                             -pad_right => 100,
                                             
                                       );  

  #然后画 标线
  my $full_length = Bio::SeqFeature::Generic->new(
                                                       -start => $panelStart, 
                                                       -end   => $panelEnd
                                                      );
  # 将标尺的 track加入到panel中。
  $panel->add_track ( $full_length,
                                         -glyph   => 'arrow',
                                         -tick    => 2,
                                         -key     => $contigNm,
                                         -fgcolor => 'black',
                                         -double  => 1
                    );
                    
  foreach my $eachPos (    sort { $a <=> $b }(   keys(  %{ $trackHash }  )   )    ){ 
    
    $panel=&NEWAddFeatureTrack_into_Panel_20170307( 
                                      $trackHash->{ $eachPos }->{'_SegmentList'},                #1
                                      $panel,                                                    #2
                                      $trackHash->{ $eachPos }->{'_SegmentKeyName'},             #3 
                                      $trackHash->{ $eachPos }->{'_dataType'},                   #4
                                      $trackHash->{ $eachPos }->{'_tkBakeCol'},                  #5
                                      $trackHash->{ $eachPos }->{'_CorePosCode'},                #6
                                      $trackHash->{ $eachPos }->{'_SECISinform'}                 #7

                                    );

  } 
  print OUT  $panel->png; 
}


sub PrintPngForGeneGroups{   #pmAble#     打印成组的 gene 结构图#$proteindIDHash, $InIN_BothSecCysHash, $outFileHash, $pngFile    inKeyName    ($SecCys, $ProID, $InIN_BothSecCysWiseHash, $InIN_BothSecCysHash, $inKeyName, $pngFile, $SECISshowOrNot)
  my ($proteindIDHash, $InIN_BothSecCysWiseHash, $InIN_BothSecCysHash, $outFileHash, $inKeyNameHash, $pngFile)=@_;
  open (OUT, ">$pngFile") or die "1 Cannot create \$pngFile=$pngFile : $!\n";
  my $trackHash; 
  my $panelStart=99999999999999;  my $panelEnd=-99999999999999; 
  my $contigNm;                  
  foreach my $SecCys (    sort { $a cmp $b }(   keys(  %{ $proteindIDHash }  )   )    ){ 
    foreach my $ProID (    sort { $a cmp $b }(   keys(  %{ $proteindIDHash->{$SecCys} }  )   )    ){ 
    	
    	if (  ( $InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_PredicFromGenoEstOrMatch'} eq 'm' ) && (  $proteindIDHash->{$SecCys}->{$ProID} eq 'gm' )  ){
    	}
    	else{
        
        $contigNm=$InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_GenomeORestContigID'};
        
        my $wiseRstHASH=$InIN_BothSecCysWiseHash->{$SecCys}->{$ProID};
        
        #my $wiseResultFile="$outFileHash->{$SecCys}/$InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_GeneWiseResult'}";  print "AAAAAAAAAA\$wiseResultFile=\$outFileHash->{$SecCys}/\$InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_GeneWiseResult'}=$outFileHash->{$SecCys}/$InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_GeneWiseResult'}=$wiseResultFile\n";
        #my $wiseRstHASH=&NewGenewisedbPraser($wiseResultFile, $InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_CoreCODONpos'}); 
        #print "\cl\clAAAAAAAAAA999999 PrintSubArrayHash::print_all_sub_array (\$wiseRstHASH)\n\n"; PrintSubArrayHash::print_all_sub_array ($wiseRstHASH); print "\n\n\cl\cl\n\n";  
        
        my $BestGenesArray=$wiseRstHASH
                                              ->{'_genewiseResultHashAd'}
                                              ->{ $wiseRstHASH->{'_BestScoreWithUhash'}->{'_BestScoreQueryProteinID'} }
                                              ->[ $wiseRstHASH->{'_BestScoreWithUhash'}->{'_BestScoreQuProtNUMBER'} ]
                                              ->{'_EachGeneArrayAd'};
  
        my $bigGeneList;
        if (ref($BestGenesArray) eq 'ARRAY'){
          for (my $i=0; $i< @{ $BestGenesArray }; $i++){ warn "\$i=$i\n\n";
          	
          	my $InputCDS_AA_sementArray=$wiseRstHASH
                                                  ->{'_genewiseResultHashAd'}
                                                  ->{ $wiseRstHASH->{'_BestScoreWithUhash'}->{'_BestScoreQueryProteinID'} }
                                                  ->[ $wiseRstHASH->{'_BestScoreWithUhash'}->{'_BestScoreQuProtNUMBER'} ]
                                                  ->{'_EachGeneArrayAd'}
                                                  ->[ $i ]
                                                  ->{'_EachGeneExonArray'};
                                                  
          	my $strand='+'; if ($InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_SecPosFrame'}<0){$strand='-';}
          
            warn "\naaaaaaaa1111 \$InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_CoreCODONpos'}=",$InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_CoreCODONpos'},"\n\n";      
            my $ExonPosList=&AddTGACoreImform_20170307( $InputCDS_AA_sementArray, '' , $strand); warn "\naaaaaaaa \$InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_CoreCODONpos'}=",$InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_CoreCODONpos'},"\n\n";
            #print "\cl\clAAAAAAAAAA111111 PrintSubArrayHash::print_all_sub_array (\$ExonPosList)\n\n"; PrintSubArrayHash::print_all_sub_array ($ExonPosList); print "\n\n\cl\cl\n\n";
          
            $bigGeneList->{$i}=$ExonPosList;
            print "\cl\clAAAAAAAAAA333333 PrintSubArrayHash::print_all_sub_array (\$ExonPosList)\n\n"; PrintSubArrayHash::print_all_sub_array ($ExonPosList); print "\n\n\cl\cl\n\n";
          }
        }     
        
        #my $tgaPos=-1;
        #if ($InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_CoreCODONtrans2AA'} eq 'U'){$tgaPos=$InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_CoreCODONpos'};}
        
        my $cosPosHere=$InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_CoreCODONpos'};
        my $cosPosFram=$InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_SecPosFrame'};
        my $cosPosStrand=1; if ( ($cosPosFram<0) || ($cosPosFram eq '-') || ($cosPosFram eq '-1') || ($cosPosFram eq '-2') || ($cosPosFram eq '-3') ){$cosPosStrand=-1;}
   
        my ($bigExonPosList,$psdPosList)=@{ &combaineForPseudoGenes_20170307($bigGeneList, $cosPosHere,$cosPosStrand ) }; 
  
  
        if (defined($InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_bestSECISstart'})){
  	      if ( $InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_bestSECISstart'}=~m/\d+/ ){  print "AAAAAAAAAA222222211111\$InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_bestSECISstart'}=$InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_bestSECISstart'}\n";
  		      $psdPosList=&AddSECISinform( $psdPosList, $InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_bestSECISstart'}, $InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_bestSECISend'}, $InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_bestSECISstrand'}, $InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_BestSECISinform'} );
  	        print "\cl\clAAAAAAAAAA22222222222 PrintSubArrayHash::print_all_sub_array (\$bigExonPosList)\n\n"; PrintSubArrayHash::print_all_sub_array ($bigExonPosList); print "\n\n\cl\cl\n\n";
      	  }
        }  
  
        $panelStart=&GetSmallestStart20170214($bigExonPosList,$panelStart);  $panelStart=&GetSmallestStart20170214($psdPosList,$panelStart);
        $panelEnd  =&GetBiggestEnd20170214   ($bigExonPosList,$panelEnd);    $panelEnd  =&GetBiggestEnd20170214   ($psdPosList,$panelEnd);                                          
   
        my $SECISInform=$InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_BestSECISinform'}; $SECISInform=ExcelHandle::GetLableOfHpLk($SECISInform);
        my $CorePosCode=$InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_CoreCODONtype'};
        my $inKeyName=$inKeyNameHash->{$SecCys}->{$ProID};
        
        $trackHash->{ $InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_CoreCODONpos'} }->{0}->{'_SegmentList'}   =$psdPosList;
        $trackHash->{ $InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_CoreCODONpos'} }->{0}->{'_SegmentKeyName'}=''; #$InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_GenomeORestContigID'};
        $trackHash->{ $InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_CoreCODONpos'} }->{0}->{'_Segmentcolor'}  ='';
        $trackHash->{ $InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_CoreCODONpos'} }->{0}->{'_SegmentProID'}  =$ProID;  
        $trackHash->{ $InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_CoreCODONpos'} }->{0}->{'_dataType'}      ='none';
        $trackHash->{ $InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_CoreCODONpos'} }->{0}->{'_tkBakeCol'}     ='ghostwhite'; #'white';#'whitesmoke';#'navajowhite'; #'whitesmoke';
        $trackHash->{ $InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_CoreCODONpos'} }->{0}->{'_CorePosCode'}      =$CorePosCode;
        $trackHash->{ $InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_CoreCODONpos'} }->{0}->{'_SECISinform'}     =$SECISInform; 
 
          
        $trackHash->{ $InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_CoreCODONpos'} }->{1}->{'_SegmentList'}   =$bigExonPosList;
        $trackHash->{ $InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_CoreCODONpos'} }->{1}->{'_SegmentKeyName'}=$inKeyName;
        $trackHash->{ $InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_CoreCODONpos'} }->{1}->{'_Segmentcolor'}  ='';
        $trackHash->{ $InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_CoreCODONpos'} }->{1}->{'_SegmentProID'}  =$ProID; 
        $trackHash->{ $InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_CoreCODONpos'} }->{1}->{'_dataType'}      ='hat'; #solid quill dushed  
        $trackHash->{ $InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_CoreCODONpos'} }->{1}->{'_tkBakeCol'}     ='ghostwhite'; #'antiquewhite';
        $trackHash->{ $InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_CoreCODONpos'} }->{1}->{'_CorePosCode'}      =$CorePosCode;
        $trackHash->{ $InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_CoreCODONpos'} }->{1}->{'_SECISinform'}     =$SECISInform; 
        	
        
        
      }
    }
  }
  my $panelLength=$panelEnd-$panelStart+1;   print "AAAAAAAAAA\$panelStart=$panelStart\t\t\$panelEnd=$panelEnd\t\t\$panelLength=$panelLength\n\n";
  $panelStart=$panelStart- int($panelLength/10); if ($panelStart<1){$panelStart=1;}
  $panelEnd=  $panelEnd  + int($panelLength/10); $panelLength=$panelEnd-$panelStart+1;
  
  $contigNm=~s/^(\S+)\s*.*$/$1/;  #-key     => $contigNm,
  
  my $panel = Bio::Graphics::Panel->new(
                                             -length    => $panelLength,   #图的全长宽度能表示的碱基数量总长度
                                             -width     => 1000,                #图的实际宽度
                                             -key_style => 'left',
                                             -start     => $panelStart,
                                             -pad_left  => 100,
                                             -pad_right => 100,
                                             
                                       );  

  #然后画 标线
  my $full_length = Bio::SeqFeature::Generic->new(
                                                       -start => $panelStart, 
                                                       -end   => $panelEnd
                                                      );
  # 将标尺的 track加入到panel中。
  $panel->add_track ( $full_length,
                                         -glyph   => 'arrow',
                                         -tick    => 2,
                                         -key     => $contigNm,
                                         -fgcolor => 'black',
                                         -double  => 1
                    );
  my $idxH=0; my $tkCol='withe';                  
  foreach my $eachPos (    sort { $a <=> $b }(   keys(  %{ $trackHash }  )   )    ){ 
    if ($idxH>0){
      $panel->add_track ( $full_length,
                                         -glyph   => 'line',                                         
                                         -bgcolor   => 'whitesmoke',
                                         -fgcolor => 'white',
                                         
                    ); 
    }
    
    if (( $idxH % 2)==0) { $tkCol='white';   } 
    else                 { $tkCol='white';   }   #ghostwhite
    foreach my $eachNb (    sort { $a <=> $b }(   keys(  %{ $trackHash->{$eachPos} }  )   )    ){
      $panel=&NEWAddFeatureTrack_into_Panel_20170307( 
                                        $trackHash->{ $eachPos }->{$eachNb}->{'_SegmentList'},                #1
                                        $panel,                                                               #2
                                        $trackHash->{ $eachPos }->{$eachNb}->{'_SegmentKeyName'},             #3 
                                        $trackHash->{ $eachPos }->{$eachNb}->{'_dataType'},                   #4
                                        $tkCol,                                                               #5
                                        $trackHash->{ $eachPos }->{$eachNb}->{'_CorePosCode'},                #6
                                        'SECIS'                                                               #7
                                      ); 
                                  
      
    }
    $idxH++;
    
  } 
  print OUT  $panel->png; 
}


sub GetSmallestStart20170214{   #pmAble#  #PrintPngForGene 和 PrintPngForGeneGroups  的子函数，用来获得最小的起始点     
	my ($inLst, $InSmSt)=@_;
	my $smallStart=99999999999999;
	if ($InSmSt<=$smallStart){$smallStart=$InSmSt;}
	if (defined ($inLst->[0])){
	  for (my $i=0; $i<@{ $inLst }; $i++){ 
      if ($inLst->[$i]->{'_headPos'}<=$smallStart){$smallStart=$inLst->[$i]->{'_headPos'};}
      if ($inLst->[$i]->{'_tailPos'}<=$smallStart){$smallStart=$inLst->[$i]->{'_tailPos'};}    
    }
  }
  return $smallStart;
}
sub GetBiggestEnd20170214{   #pmAble#  #PrintPngForGene 和 PrintPngForGeneGroups 的子函数，用来获得最大的起终止点
	my ($inLst, $inBgEd)=@_;
	my $bigEnd=-99999999999999;
	if ($inBgEd>=$bigEnd){$bigEnd=$inBgEd;}
	if (defined ($inLst->[0])){
	  for (my $i=0; $i<@{ $inLst }; $i++){
      if ($inLst->[$i]->{'_headPos'}>=$bigEnd    ){$bigEnd   =$inLst->[$i]->{'_headPos'};} 
      if ($inLst->[$i]->{'_tailPos'}>=$bigEnd    ){$bigEnd   =$inLst->[$i]->{'_tailPos'};}    
    }
  }
  return $bigEnd;
}


sub AddTGACoreImform{    #这个是老版本， 已经被 新版本 AddTGACoreImform_20170307，所取代
  my ( $CDS_AA_sementArray, $corePos, $InStrand)=@_;       
 
  my %keysUsedHere={
    
    '_ExonHead'                 => "蛋白和DNA比对全部结果中，某个Protein结果的，其中某个Gene的，其中某个Exon的，开头的 位置  ",
    '_ExonTail'                 => "蛋白和DNA比对全部结果中，某个Protein结果的，其中某个Gene的，其中某个Exon的，结尾的 位置  ",

    
    '_SegmentType'                    =>  "某个DNA区段的 类型，具体类型，就是下面这几钟类型",
    
    '_CDSinCoreExon'                  =>  "某个DNA区段的 类型，和Core 如Tga或tgc在一个外显子中",
    '_CorePos'                        =>  "某个DNA区段的 类型，core的位置，如tga的t的位置",
    '_noCoreExon'                     =>  "某个DNA区段的 类型，没有Core的位置",
  };
  
  my $outHash; my $j=0;
  for (my $i=0; $i<@{ $CDS_AA_sementArray }; $i++){ 
  	my $leftPos =$CDS_AA_sementArray->[$i]->{'_ExonHead'}; if ($InStrand eq '-')  {$leftPos  =$CDS_AA_sementArray->[$i]->{'_ExonTail'};}
  	my $rightPos=$CDS_AA_sementArray->[$i]->{'_ExonTail'}; if ($InStrand eq '-')  {$rightPos =$CDS_AA_sementArray->[$i]->{'_ExonHead'};}
  	if ( ($leftPos==$corePos) && ($rightPos==$corePos) ){
  		$outHash->[$j]->{'_headPos'}=$corePos;
      $outHash->[$j]->{'_tailPos'}=$corePos;
      $outHash->[$j]->{'_SegmentType'}='_CorePos';
      $outHash->[$j]->{'_Strand'}=$InStrand;
      $j++;
      
  	}
  	elsif ( ($leftPos==$corePos) && ($rightPos>$corePos) ){
  		$outHash->[$j]->{'_headPos'}=$corePos;
      $outHash->[$j]->{'_tailPos'}=$corePos;
      $outHash->[$j]->{'_SegmentType'}='_CorePos';
      $outHash->[$j]->{'_Strand'}=$InStrand;
      $j++;
      $outHash->[$j]->{'_headPos'}=$corePos+1;
      $outHash->[$j]->{'_tailPos'}=$rightPos;
      $outHash->[$j]->{'_SegmentType'}='_CDSinCoreExon';
      $outHash->[$j]->{'_Strand'}=$InStrand;
      $j++;
      
  	}
  	elsif ( ($leftPos<$corePos) && ($rightPos==$corePos) ){
  		$outHash->[$j]->{'_headPos'}    =$leftPos;
      $outHash->[$j]->{'_tailPos'}    =$corePos-1;
      $outHash->[$j]->{'_SegmentType'}='_CDSinCoreExon';
      $outHash->[$j]->{'_Strand'}=$InStrand;
      $j++;
      $outHash->[$j]->{'_headPos'}=$corePos;
      $outHash->[$j]->{'_tailPos'}=$corePos;
      $outHash->[$j]->{'_SegmentType'}='_CorePos';
      $outHash->[$j]->{'_Strand'}=$InStrand;
      $j++;
  	}
    elsif ( ($leftPos<$corePos) && ($corePos<$rightPos) ){
      $outHash->[$j]->{'_headPos'}    =$leftPos;
      $outHash->[$j]->{'_tailPos'}    =$corePos-1;
      $outHash->[$j]->{'_SegmentType'}='_CDSinCoreExon';
      $outHash->[$j]->{'_Strand'}=$InStrand;
      $j++;
      $outHash->[$j]->{'_headPos'}=$corePos;
      $outHash->[$j]->{'_tailPos'}=$corePos;
      $outHash->[$j]->{'_SegmentType'}='_CorePos';
      $outHash->[$j]->{'_Strand'}=$InStrand;
      $j++;
      $outHash->[$j]->{'_headPos'}=$corePos+1;
      $outHash->[$j]->{'_tailPos'}=$rightPos;
      $outHash->[$j]->{'_SegmentType'}='_CDSinCoreExon';
      $outHash->[$j]->{'_Strand'}=$InStrand;
      $j++;
    }
    else {
      $outHash->[$j]->{'_headPos'}    =$leftPos;
      $outHash->[$j]->{'_tailPos'}    =$rightPos;
      $outHash->[$j]->{'_SegmentType'}='_noCoreExon';
      $outHash->[$j]->{'_Strand'}=$InStrand;
      $j++;
    }
    
  }
  
  
  return $outHash;

}

sub AddSECISinform{   #pmAble#  #PrintPngForGene 和 PrintPngForGeneGroups 的子函数， #添加SECIS信息
  my ( $inArray, $SECISstart, $SECISend, $SECISstrand, $SECISinform)=@_;        
  if ($SECISstart<0){die "\$SECISstart=$SECISstart\n";}
  if ($SECISstart>$SECISend){my $tp=$SECISend; $SECISend=$SECISstart; $SECISstart=$tp;} # $SECISstart 被确定为较小的数字
  my $realSECIShead=$SECISstart; my $realSECIStail=$SECISend;                                                 #正义链 小 到 大
  if (  ($SECISstrand eq '-') || ($SECISstrand < 0) ){ $realSECIShead=$SECISend; $realSECIStail=$SECISstart;} #负义链 大 到 小
  
  my %keysUsedHere={
    

    
    '_SegmentType'                    =>  "某个DNA区段的 类型，具体类型，就是下面这几钟类型",
    

    '_SECIS'                          =>  "某个DNA区段的 类型，SECIS",
  };
  
  my $outHash;
   
  if ( ($SECISstart > 0) && ($SECISstart=~m/^\d+$/) ){
  	
    my $SortFu_ZhHash;
    
    if (ref($inArray) eq 'ARRAY'){
      if ( ($inArray->[0]->{'_Strand'} eq '-') || ($inArray->[0]->{'_Strand'} < 0) ){
        
        my $j=0;
        for (my $i=0; $i<@{ $inArray }; $i++){
          $SortFu_ZhHash->{$j}->{ 'HashBody' }=$inArray->[$i];
          $SortFu_ZhHash->{$j}->{ 'sortNB' }  =$inArray->[$i]->{'_headPos'};
          $j++;
        }
        my $tempHash;      
        $tempHash->{'_headPos'}    =$realSECIShead;    
        $tempHash->{'_tailPos'}    =$realSECIStail;
        $tempHash->{'_SECISinform'}=$SECISinform;
        $tempHash->{'_Strand'}     =$SECISstrand;
        $tempHash->{'_SegmentType'}='_SECIS'; 
        $SortFu_ZhHash->{$j}->{ 'HashBody' }=$tempHash; 
        $SortFu_ZhHash->{$j}->{ 'sortNB' }  =$SECISend;
        my $k=0;
        foreach my $tpk (    sort { $SortFu_ZhHash->{$b}->{'sortNB'} <=> $SortFu_ZhHash->{$a}->{'sortNB'} } (   keys (  %{ $SortFu_ZhHash }  )   )    ){
          $outHash->[$k]=$SortFu_ZhHash->{$tpk}->{'HashBody'};
          $k++;
        } 
            
      }
      else {
      	
      	my $j=0;
        for (my $i=0; $i<@{ $inArray }; $i++){
          $SortFu_ZhHash->{$j}->{ 'HashBody' }=$inArray->[$i];
          $SortFu_ZhHash->{$j}->{ 'sortNB' }  =$inArray->[$i]->{'_headPos'};
          $j++;
        }
        my $tempHash;      
        $tempHash->{'_headPos'}    =$realSECIShead;    
        $tempHash->{'_tailPos'}    =$realSECIStail;
        $tempHash->{'_SECISinform'}=$SECISinform;
        $tempHash->{'_Strand'}     =$SECISstrand;
        $tempHash->{'_SegmentType'}='_SECIS'; 
        $SortFu_ZhHash->{$j}->{ 'HashBody' }=$tempHash; 
        $SortFu_ZhHash->{$j}->{ 'sortNB' }  =$SECISstart;
        my $k=0;
        foreach my $tpk (    sort { $SortFu_ZhHash->{$a}->{'sortNB'} <=> $SortFu_ZhHash->{$b}->{'sortNB'} } (   keys (  %{ $SortFu_ZhHash }  )   )    ){
          $outHash->[$k]=$SortFu_ZhHash->{$tpk}->{'HashBody'};
          $k++;
        }      	
      }      
    }
    else{    
      $outHash->[0]->{'_headPos'}    =$realSECIShead;    
      $outHash->[0]->{'_tailPos'}    =$realSECIStail;
      $outHash->[0]->{'_SECISinform'}=$SECISinform;
      $outHash->[0]->{'_Strand'}     =$SECISstrand;
    }
  }   #elsif (){   #}
  else{
    $outHash=$inArray;  print "NNNNNNNNNNNNNNNNNNNNNNNNNNNNoooo SECIS input Here\n\n";
  }
  
  
  return $outHash;
}




sub AddTGACoreImform_20170307{        #pmAble#   # PrintPngForGene的子函数，添加TGAcore这个Exon类型，及相应的数据
  my ( $CDS_AA_sementArray, $corePos, $InStrand)=@_;       
 
  my $strandNb=1; if ( ($InStrand<0) || ($InStrand eq '-') ){$strandNb=-1;}
  my %keysUsedHere={
    
    '_ExonHead'                 => "蛋白和DNA比对全部结果中，某个Protein结果的，其中某个Gene的，其中某个Exon的，开头的 位置  ",
    '_ExonTail'                 => "蛋白和DNA比对全部结果中，某个Protein结果的，其中某个Gene的，其中某个Exon的，结尾的 位置  ",

    
    '_SegmentType'                    =>  "某个DNA区段的 类型，具体类型，就是下面这几钟类型",
    
    '_CDSinCoreExon'                  =>  "某个DNA区段的 类型，和Core 如Tga或tgc在一个外显子中",
    '_CorePos'                        =>  "某个DNA区段的 类型，core的位置，如tga的t的位置",
    '_noCoreExon'                     =>  "某个DNA区段的 类型，没有Core的位置",
  };
  
  my $step=3;
  my $pos_exonIdx_hash;
  my $headTailHash;
  my $GenoPos_CDSpos_hash; my $idx=1;
  if (ref($CDS_AA_sementArray) eq 'ARRAY'){
    for (my $i=0; $i<@{ $CDS_AA_sementArray }; $i++){
    	my $left= $CDS_AA_sementArray->[$i]->{'_ExonHead'}; warn "\$left= \$CDS_AA_sementArray->[$i]->{'_ExonHead'}=$CDS_AA_sementArray->[$i]->{'_ExonHead'}\n";
    	my $right=$CDS_AA_sementArray->[$i]->{'_ExonTail'}; warn "\$right=\$CDS_AA_sementArray->[$i]->{'_ExonTail'}=$CDS_AA_sementArray->[$i]->{'_ExonTail'}\n";
    	my $strandNBHr=1; my $length; warn "\$InStrand=$InStrand\n";
    	if    (($InStrand >0) || ($InStrand eq '+')){$strandNBHr=1; $length=$right-$left+1; if ($left<=$right){} else {print "\n\n\nprint_all_sub_array(\$CDS_AA_sementArray)\n\n";print_all_sub_array($CDS_AA_sementArray); die"\$InStrand=$InStrand \n\$CDS_AA_sementArray->[$i]->{'_ExonHead'}=$CDS_AA_sementArray->[$i]->{'_ExonHead'}\n";}}	
      elsif (($InStrand <0) || ($InStrand eq '-')){$strandNBHr=-1;$length=$left-$right+1; if ($left>=$right){} else {print "\n\n\nprint_all_sub_array(\$CDS_AA_sementArray)\n\n";print_all_sub_array($CDS_AA_sementArray); die"\$InStrand=$InStrand \n\$CDS_AA_sementArray->[$i]->{'_ExonHead'}=$CDS_AA_sementArray->[$i]->{'_ExonHead'}\n";}}	
      
      $headTailHash->{'Up2Down'}->{$i}->{'_ExonHead'}=$CDS_AA_sementArray->[$i]->{'_ExonHead'};           $headTailHash->{'Up2Down'}->{$i}->{'_ExonTail'}=$CDS_AA_sementArray->[$i]->{'_ExonTail'};
      $headTailHash->{'posTypeHash'}->{'_ExonHead'}->{ $CDS_AA_sementArray->[$i]->{'_ExonHead'} }=1;      $headTailHash->{'posTypeHash'}->{'_ExonTail'}->{ $CDS_AA_sementArray->[$i]->{'_ExonTail'} }=1;
      $headTailHash->{'Down2Up'}->{ $CDS_AA_sementArray->[$i]->{'_ExonHead'} }->{'exonIdx'}=$i;           $headTailHash->{'Down2Up'}->{ $CDS_AA_sementArray->[$i]->{'_ExonTail'} }->{'exonIdx'}=$i;
      
       
      #foreach my $eachPos ($left .. $right){
      for (my $j = 0; $j < $length; $j++) {
      	my $eachPos=$left+$j*$strandNBHr;
      	$pos_exonIdx_hash->{'Up2Down'}->{$eachPos}=$i;      #warn "\$pos_exonIdx_hash->{'Up2Down'}->{$eachPos}=\$i=$i\n";
      	$pos_exonIdx_hash->{'Down2Up'}->{$i}=$eachPos;      #warn "\$pos_exonIdx_hash->{'Down2Up'}->{$i}=\$eachPos=$eachPos\n";
      	
      	$GenoPos_CDSpos_hash->{'Up2Down'}->{$eachPos}=$idx; #warn "\$GenoPos_CDSpos_hash->{'Up2Down'}->{$eachPos}=\$idx=$idx\n";
      	$GenoPos_CDSpos_hash->{'Down2Up'}->{$idx}=$eachPos; #warn "\$GenoPos_CDSpos_hash->{'Down2Up'}->{$idx}=\$eachPos=$eachPos\n";
      	$idx++;
      }
         
    }
  }
  
  my $CorPosInCDS;   my $corePosSegmentHash; 
  if ($corePos>0){   #看是不是有输入corpos的数据
  	if ($corePos=~m/^\d+$/){ #看是不是有输入corpos的数据
  	  if (defined($GenoPos_CDSpos_hash->{'Up2Down'}->{$corePos})) { #看看是不是出现在前面的 外显子遍历数据中
  	    $CorPosInCDS=$GenoPos_CDSpos_hash->{'Up2Down'}->{$corePos}; warn "\n\n\n\$CorPosInCDS=\$GenoPos_CDSpos_hash->{'Up2Down'}->{$corePos}=$CorPosInCDS\n\n"; 
  	    for (my $i=0; $i<$step; $i++){
          my $inCdsCorSegmentPos=$CorPosInCDS+$i; my $CorSegmentGenoPos=$GenoPos_CDSpos_hash->{'Down2Up'}->{$inCdsCorSegmentPos};
          $corePosSegmentHash->{'Up2Down'}->{$CorSegmentGenoPos}=$inCdsCorSegmentPos; warn "\n\$corePosSegmentHash->{'Up2Down'}->{$CorSegmentGenoPos}=\$inCdsCorSegmentPos=$inCdsCorSegmentPos\n";
          $corePosSegmentHash->{'Down2Up'}->{$inCdsCorSegmentPos}=$CorSegmentGenoPos; warn "\n\$corePosSegmentHash->{'Down2Up'}->{$inCdsCorSegmentPos}=$\CorSegmentGenoPos=CorSegmentGenoPos\n\n";
        }
  	  }
  	}
  }
   
  my $eachPosExonTypeHash;  
  #下面是用Cds中的位置作为key来排序的，所以可以用sort，如果用genopos 则不行，因为 存在 反义互补的反向链
  foreach my $i_idx (    sort { $a <=> $b } (   keys (  %{ $headTailHash->{'Up2Down'} }  )   )    ){ print  "\$i_idx=$i_idx\n\$headTailHash->{'Up2Down'}->{$i_idx}->{'_ExonHead'}=$headTailHash->{'Up2Down'}->{$i_idx}->{'_ExonHead'}\n\$headTailHash->{'Up2Down'}->{$i_idx}->{'_ExonTail'}=$headTailHash->{'Up2Down'}->{$i_idx}->{'_ExonTail'}\n";
  	my $coreInExon=0;
  	my $CdsHead=$GenoPos_CDSpos_hash->{'Up2Down'}->{ $headTailHash->{'Up2Down'}->{$i_idx}->{'_ExonHead'} }; print "\$CdsHead=\$GenoPos_CDSpos_hash->{'Up2Down'}->{ $headTailHash->{'Up2Down'}->{$i_idx}->{'_ExonHead'} }=$GenoPos_CDSpos_hash->{'Up2Down'}->{ $headTailHash->{'Up2Down'}->{$i_idx}->{'_ExonHead'} }\n";
  	my $CdsTail=$GenoPos_CDSpos_hash->{'Up2Down'}->{ $headTailHash->{'Up2Down'}->{$i_idx}->{'_ExonTail'} }; print "\$CdsTail=\$GenoPos_CDSpos_hash->{'Up2Down'}->{ $headTailHash->{'Up2Down'}->{$i_idx}->{'_ExonTail'} }=$GenoPos_CDSpos_hash->{'Up2Down'}->{ $headTailHash->{'Up2Down'}->{$i_idx}->{'_ExonTail'} }\n";
    foreach my $eachCDSPos ( $CdsHead .. $CdsTail )  {
    	if (  defined ( $corePosSegmentHash->{'Down2Up'}->{$eachCDSPos} )  ){
    		$coreInExon=1; warn "\$corePosSegmentHash->{'Down2Up'}->{$eachCDSPos}=$corePosSegmentHash->{'Down2Up'}->{$eachCDSPos}\n";
    	}
    }
    
    my $beforeCore=1;
    foreach my $eachCDSPos ( $CdsHead .. $CdsTail )  {                   #print "\$eachCDSPos=$eachCDSPos\n";
    	my $eachGenoPos=$GenoPos_CDSpos_hash->{'Down2Up'}->{$eachCDSPos};  #print "\$eachGenoPos=$eachGenoPos\n";
    	if ($coreInExon==1){
    	  if (  defined ( $corePosSegmentHash->{'Up2Down'}->{$eachGenoPos} )  ){
    	  	$eachPosExonTypeHash->{'GenoPos'}->{$eachGenoPos}='_CorePos';    $eachPosExonTypeHash->{'CDSPos'}->{$eachCDSPos}='_CorePos'; 
    		  $beforeCore=0;
    	  }
    	  else {
    	    if ($beforeCore==1){
    	    	$eachPosExonTypeHash->{'GenoPos'}->{$eachGenoPos}='_BeforCore';   $eachPosExonTypeHash->{'CDSPos'}->{$eachCDSPos}='_BeforCore';
    	    }
      	  else {
    	    	$eachPosExonTypeHash->{'GenoPos'}->{$eachGenoPos}='_AfterCore';   $eachPosExonTypeHash->{'CDSPos'}->{$eachCDSPos}='_AfterCore';
    	    }
    	  }
      }
      else{
      	$eachPosExonTypeHash->{'GenoPos'}->{$eachGenoPos}='_noCoreExon';   $eachPosExonTypeHash->{'CDSPos'}->{$eachCDSPos}='_noCoreExon';
      	#print "\$eachPosExonTypeHash->{'GenoPos'}->{$eachGenoPos}=$eachPosExonTypeHash->{'GenoPos'}->{$eachGenoPos}\t\$eachPosExonTypeHash->{'CDSPos'}->{$eachCDSPos}=$eachPosExonTypeHash->{'CDSPos'}->{$eachCDSPos}\n";
      }
    }
    
  }
  
  #print_all_sub_array($eachPosExonTypeHash);
  my $lastPos;
  MARKLASTPOS: foreach my $CdsPos (    sort { $b <=> $a } (   keys (  %{ $eachPosExonTypeHash->{'CDSPos'} }  )   )    ){$lastPos=$CdsPos; last MARKLASTPOS;}
  
  my $outHash; my $j=0;
  foreach my $CdsPos (    sort { $a <=> $b } (   keys (  %{ $eachPosExonTypeHash->{'CDSPos'} }  )   )    ){
    my $genoPos=$GenoPos_CDSpos_hash->{'Down2Up'}->{$CdsPos};
    
    my $segmentHeadorNot=0;
    if ($CdsPos == 0)                                                                                                                                                 { $segmentHeadorNot=1; }   #所有位置的开头
    if (   (  defined ( $headTailHash->{'posTypeHash'}->{'_ExonHead'}->{ $genoPos } )  ) && (  $headTailHash->{'posTypeHash'}->{'_ExonHead'}->{ $genoPos } == 1  )   ){ $segmentHeadorNot=1; } 	 #原来单个Exon的开头
    if ($segmentHeadorNot==0){ my $lastCDSPos=$CdsPos-1; if ( $eachPosExonTypeHash->{'CDSPos'}->{$CdsPos} ne $eachPosExonTypeHash->{'CDSPos'}->{$lastCDSPos} )        { $segmentHeadorNot=1; } } #不同类型segment的开始
    if ($segmentHeadorNot==1){ 
    	$outHash->[$j]->{'_headPos'}     =$genoPos;
      $outHash->[$j]->{'_SegmentType'} =$eachPosExonTypeHash->{'CDSPos'}->{$CdsPos};
      $outHash->[$j]->{'_Strand'}     =$InStrand;
    }
    
    my $segmentTailorNot=0; my $tailInChangeType=0;
    if ($CdsPos == $lastPos)                                                                                                                                          { $segmentTailorNot=1; }   #所有位置的tail
    if (   (  defined ( $headTailHash->{'posTypeHash'}->{'_ExonTail'}->{ $genoPos } )  ) && (  $headTailHash->{'posTypeHash'}->{'_ExonTail'}->{ $genoPos } == 1  )   ){ $segmentTailorNot=1; } 	 #原来单个Exon的tail
    if ($segmentTailorNot==0){ my $NextCDSPos=$CdsPos+1; if ( $eachPosExonTypeHash->{'CDSPos'}->{$CdsPos} ne $eachPosExonTypeHash->{'CDSPos'}->{$NextCDSPos} )        { $segmentTailorNot=1; $tailInChangeType=1;} } #不同类型segment的tail
    if ($segmentTailorNot==1){ 
    	$outHash->[$j]->{'_tailPos'}     =$genoPos; if ($tailInChangeType==1){$outHash->[$j]->{'_tailPos'}     =$genoPos-1*$strandNb;}
      $outHash->[$j]->{'_SegmentType'} =$eachPosExonTypeHash->{'CDSPos'}->{$CdsPos};
      $outHash->[$j]->{'_Strand'}     =$InStrand;
      $j++;
    }
    
  }
  
  return $outHash;

}

sub combaineForPseudoGenes_20170307 { #pmAble#   # PrintPngForGene的子函数 #添加PseudoGene的类型 及相关数据
  my ($In_bigExonPosList, $corPos, $strandHere)=@_;
  
  my %keysUsedHere={
    
    '_ExonHead'                 => "蛋白和DNA比对全部结果中，某个Protein结果的，其中某个Gene的，其中某个Exon的，开头的 位置  ",
    '_ExonTail'                 => "蛋白和DNA比对全部结果中，某个Protein结果的，其中某个Gene的，其中某个Exon的，结尾的 位置  ",

    
    '_SegmentType'                    =>  "某个DNA区段的 类型，具体类型，就是下面这几钟类型",
    
    '_CDSinCoreExon'                  =>  "某个DNA区段的 类型，和Core 如Tga或tgc在一个外显子中",
    '_CorePos'                        =>  "某个DNA区段的 类型，core的位置，如tga的t的位置",
    '_noCoreExon'                     =>  "某个DNA区段的 类型，没有Core的位置",
  };
  my $outHash; my $j=0; #my $lastIdx;
  my $outHash_psd; my $k=0; #my $lastIdx;
  my $strandNb=1;
  foreach my $idx (    sort {$a <=> $b} (   keys (  %{ $In_bigExonPosList }  )   )    ){
    if ($idx>0){
    	if ( ($outHash->[$j-1]->{'_Strand'}<0) || ($outHash->[$j-1]->{'_Strand'} eq '-') ){$strandNb=-1;}
    	
    	
    	$outHash_psd->[$k]->{'_headPos'}     =$outHash->[$j-1]->{'_tailPos'}+$strandNb*1;                   #$outHash->[$j-1]->{'_tailPos'}=$outHash->[$j-1]->{'_tailPos'}-$strandNb*1;  
    	$outHash_psd->[$k]->{'_tailPos'}     =$In_bigExonPosList->{$idx}->[0]->{'_headPos'}-$strandNb*1;
    	$outHash_psd->[$k]->{'_SegmentType'} ='_PseudoPos';
      $outHash_psd->[$k]->{'_Strand'}      =$outHash->[$j-1]->{'_Strand'} ;
      $outHash_psd->[$k]->{'_length'}      =$strandNb*($outHash_psd->[$k]->{'_tailPos'}-$outHash_psd->[$k]->{'_headPos'})+1;
      $k++;
    }
    if (ref($In_bigExonPosList->{$idx}) eq 'ARRAY'){
      for ( my $i=0; $i< @{ $In_bigExonPosList->{$idx} }; $i++ ){
      	if ( SeqSegmentsTools::a_between_b_c($corPos, $In_bigExonPosList->{$idx}->[$i]->{'_headPos'}, $In_bigExonPosList->{$idx}->[$i]->{'_tailPos'}) ){
      	  $outHash_psd->[$k]->{'_headPos'}     = $corPos;
          $outHash_psd->[$k]->{'_tailPos'}     = $corPos;
          $outHash_psd->[$k]->{'_SegmentType'} ='_CorePos';
          $outHash_psd->[$k]->{'_Strand'}      =$strandHere ;
          $k++;
      	}
        $outHash->[$j]=$In_bigExonPosList->{$idx}->[$i];
        #if ( ( 0==$i ) &&  ($idx>0) ){        	my $strandNb=1; if ( ($outHash->[$j]->{'_Strand'}<0) || ($outHash->[$j]->{'_Strand'} eq '-') ){$strandNb=-1;}         $outHash->[$j]->{'_headPos'}    =    $outHash->[$j]->{'_headPos'}+$strandNb*1;      }
        $j++;
      }
      #$lastIdx=$idx;
    }
  }
      
  return [$outHash,$outHash_psd];

  
}     

#NEWAddFeatureTrack_into_Panel_20170307( 
#                                      $trackHash->{ $eachPos }->{'_SegmentList'},                #1
#                                      $panel,                                                    #2
#                                      $trackHash->{ $eachPos }->{'_SegmentKeyName'},             #3 
#                                      $trackHash->{ $eachPos }->{'_dataType'},                   #4
#                                      $trackHash->{ $eachPos }->{'_tkBakeCol'},                  #5
#                                      $trackHash->{ $eachPos }->{'_CorePosCode'},                #6
#                                      $trackHash->{ $eachPos }->{'_SECISinform'}                 #7
#
#                                    );


sub NEWAddFeatureTrack_into_Panel_20170307 {  #pmAble#   # PrintPngForGene #PrintPngForGene 的子函数，用来将Feature信息，加入画结构的模块
  my ( $AAposList,                       #1
       $tpPanel,                         #2
       $idKeyHere,                       #3
       $dataType,                        #4
       $tkBgCol,                         #5
       $CorePosCode,                     #6
       $SECISinform,                     #7
       )=@_;  
  warn "\$AAposList=$AAposList, \$tpPanel=$tpPanel, \$idKeyHere=$idKeyHere, \$CorePosCode=$CorePosCode, \$SECISinform=$SECISinform\n";  
  print   "\$AAposList=$AAposList, \$tpPanel=$tpPanel, \$idKeyHere=$idKeyHere, \$CorePosCode=$CorePosCode, \$SECISinform=$SECISinform\n";
  
  my $corCol="red"; if ($CorePosCode=~m/^\s*TGA\s*$/i){} else {$corCol="black";}
    
  #先建立各个分段的feature
  #$e1 = Bio::SeqFeature::Lite->new(-start=>1,-stop=>100,-type=>'exon');
  my $allSubFeatures;  my $CdsStrand;
  if (ref($AAposList) eq 'ARRAY'){ 
  	if ( ($AAposList->[0]->{'_Strand'} eq '-') || ($AAposList->[0]->{'_Strand'} < 0) ){$CdsStrand=-1;}
  	
  }
  my @newAAposList; if (ref($AAposList) eq 'ARRAY'){@newAAposList= @{ $AAposList };}
  #if ($CdsStrand==-1){@newAAposList= reverse @{ $AAposList };}
  if (ref($AAposList) eq 'ARRAY'){
    for (my $i=0; $i< @newAAposList; $i++){   print "\$newAAposList[$i]->{'_headPos'}=$newAAposList[$i]->{'_headPos'},\$newAAposList[$i]->{'_tailPos'}=$newAAposList[$i]->{'_tailPos'},\$newAAposList[$i]->{'_SegmentType'}=$newAAposList[$i]->{'_SegmentType'},\$newAAposList[$i]->{'_Strand'}=$newAAposList[$i]->{'_Strand'}\n";
      my $strand=1; if ( ($newAAposList[$i]->{'_Strand'} eq '-') || ($newAAposList[$i]->{'_Strand'} < 0) ){$strand=-1;}  warn "111111111111\$strand=$strand\n";  print "\$strand=$strand\n";
      #if ($newAAposList[$i]->{'_SegmentType'}=~/Exon/ ){ 
      $CdsStrand=$strand; 
      #}
      $allSubFeatures->[$i]=Bio::SeqFeature::Lite->new( -start   => $newAAposList[$i]->{'_headPos'}, 
                                                        -end     => $newAAposList[$i]->{'_tailPos'}, 
                                                        -name    => $newAAposList[$i]->{'_SegmentType'}, 
                                                        -strand  => $strand                                #$newAAposList[$i]->{'_Strand'} 
                                                        );    warn "22222222222222\$strand=$strand\n";  print "\$strand=$strand\n"; #print "\cl PrintSubArrayHash::print_all_sub_array(\$allSubFeatures->[$i]=$allSubFeatures->[$i])\n"; PrintSubArrayHash::print_all_sub_array($allSubFeatures->[$i]); print "\cl";
    }
  }
  
  
  #然后建立 总Feature
  my $AlignFeature = Bio::SeqFeature::Lite->new (-segments => $allSubFeatures,
                                                   #-key        => $idKeyHere,  
                                                 #-display_name =>$idKeyHere, 
                                                 -id       => "",#$idKeyHere, 
                                                 -strand   => $CdsStrand   #$newAAposList[0]->{'_Strand'} 
                                                 ); print "\$newAAposList[0]->{'_Strand'}=$newAAposList[0]->{'_Strand'}\n\n\cl PrintSubArrayHash::print_all_sub_array(\$AlignFeature=$AlignFeature)\n"; PrintSubArrayHash::print_all_sub_array($AlignFeature); print "\cl";
  #然后画各个比对结果
  #my $colorIdx=$tpPanel->translate_color('#FFF8DC'); warn "\$colorIdx=$colorIdx\n\n\n";sleep(1);
  
  $tpPanel->add_track(   $AlignFeature ,
                                           -key        => $idKeyHere,                
                                           -glyph       => 'transcript',
                                           -tkcolor     =>  $tkBgCol, #'whitesmoke', #$colorIdx, #'#FFF8DC',#'lightyellow',#'cyan',#'gray',
                                           #-box_subparts=> 3,
                                           -strand      => $CdsStrand,
                                           -connector   => $dataType,#'none',,#'dashed',
                                           -connector_color   => 'black',
                                           -bgcolor     =>  #$colorHere,
                                                            sub {  #'green',
                                                                  my $feature = shift;  #print "\cl 111111111 PrintSubArrayHash::print_all_sub_array(\$feature)=$feature)\n"; PrintSubArrayHash::print_all_sub_array($feature); print "\cl";
                                                                  warn "\n\nInFeature :", $feature->start,"..",$feature->end,"\n"; 
                                                                  my $returnColor= 'green'; 
                                                                  my $segType=$feature->seqname;  warn "\n\$segType=$segType\n";# \t\$colorHere=$colorHere\n";
                                                                  
                                                                  #if    ($colorHere=~m/\S+/  )          {$returnColor=$colorHere;}
                                                                  if    ($segType eq '_CDSinCoreExon')  {$returnColor='lightgreen';}
                                                                  elsif ($segType eq '_noCoreExon')     {$returnColor='lightgreen';}
                                                                  elsif ($segType eq '_AfterCore')      {$returnColor='lightgreen';}
                                                                  elsif ($segType eq '_CorePos')        {$returnColor=$corCol;}
                                                                  elsif ($segType eq '_SECIS')          {$returnColor='orange';}
                                                                  elsif ($segType eq '_BeforCore')      {$returnColor='lightgreen';}    
                                                                  elsif ($segType eq '_PseudoPos')      {$returnColor='blue';}
                                                                  
                                                                  
                                                                  warn "\$segType=$segType\t\$returnColor=$returnColor \n\n";
                                                                  $returnColor;
                                                                },  
                                           
                                           
                                           
                                           -fgcolor     =>                       #'',#'yellow',    -nothing =>
                                                           sub {  #'green',
                                                                  my $feature = shift;  #print "\cl 111111111 PrintSubArrayHash::print_all_sub_array(\$feature)=$feature)\n"; PrintSubArrayHash::print_all_sub_array($feature); print "\cl";
                                                                  warn "\n\nInFeature :", $feature->start,"..",$feature->end,"\n"; 
                                                                  my $returnColor= 'green'; 
                                                                  my $segType=$feature->seqname;  warn "\n\$segType=$segType\n";  #\t\$colorHere=$colorHere\n";
                                                                  
                                                                  #if    ($colorHere=~m/\S+/  )          {$returnColor=$colorHere;}
                                                                  if    ($segType eq '_CDSinCoreExon')  {$returnColor='lightgreen';}
                                                                  elsif ($segType eq '_noCoreExon')     {$returnColor='lightgreen';}
                                                                  elsif ($segType eq '_AfterCore')      {$returnColor='lightgreen';}
                                                                  elsif ($segType eq '_CorePos')        {$returnColor=$corCol;}
                                                                  elsif ($segType eq '_SECIS')          {$returnColor='orange';}
                                                                  elsif ($segType eq '_BeforCore')      {$returnColor='lightgreen';}    
                                                                  elsif ($segType eq '_PseudoPos')      {$returnColor='blue';}
                                                                  
                                                                  
                                                                  warn "\$segType=$segType\t\$returnColor=$returnColor \n\n";
                                                                  $returnColor;
                                                                },  
                                           -font2color  => 'red',
                                           #-key         => $keyHere,
                                           -bump        =>  3, 
                                           -hbumppad    => 3,
                                           -opacity     => 2,
                                           -sort_order  => 'shortest',
                                           -height      =>  12,
                                           -label       =>  1, 
                                           -part_labels_merge  => 0,
                                           -part_labels =>#1  #
                                                            sub {
                                                                  my ($feature,undef,$partno) = @_; 
                                                                  my $tp=$feature->{segments}->[$partno]->name; warn "\n\$feature->{segments}->[$partno]->name=$tp\n\@_=@_\n\n\n"; #print "\@_=@_\n\n\cl 2222222222 PrintSubArrayHash::print_all_sub_array(\$feature)=$feature)\n"; PrintSubArrayHash::print_all_sub_array($feature); print "\cl";
                                                                  my $legth=$feature->{segments}->[$partno]->length;
                                                                  #if (defined ($_[4])){print "\@_=@_\n\n\cl 2222222222 PrintSubArrayHash::print_all_sub_array(\$_[4]=$_[4])\n"; PrintSubArrayHash::print_all_sub_array($_[4]); print "\cl";}
                                                                  my $returnLable="Psd";
                                                                  if      ($tp eq '_CorePos'){ $returnLable=$CorePosCode; }
                                                                  elsif   ($tp eq '_SECIS'){ $returnLable=$SECISinform; }
                                                                  elsif   ($tp eq '_noCoreExon'){ $returnLable=" "; }
                                                                  elsif   ($tp eq '_BeforCore'){ $returnLable=" "; }   
                                                                  elsif   ($tp eq '_AfterCore'){ $returnLable=" "; } 
                                                                  elsif   ($tp eq '_PseudoPos'){ $returnLable=" $legth bp idel"; } 
                                                                 
                                                                  	$returnLable;
                                                                  #return $returnLable;
                                                                }
                                           #1,             #显示外显子编号
                                           #-description =>  #1
                                            #              sub {
                                            #                         my $feature = shift; print "\cl 2222222222 PrintSubArrayHash::print_all_sub_array(\$feature)=$feature)\n"; PrintSubArrayHash::print_all_sub_array($feature); print "\cl";
                                            #                         #retrun $feature->name;
                                            #                         return unless $feature->has_tag('description');
                                            #                         my ($description) = $feature->name;
                                            #                         "\$description=$description";
                                            #                    },
                                            
                     );
  return $tpPanel;                 
    


}



sub NEWAddFeatureTrack3Panel{  #这个应该是个老版本的，请使用新版本 NEWAddFeatureTrack_into_Panel_20170307
  my ($AAposList, $tpPanel, $idKeyHere, $colorHere, $keyHere)=@_;  warn "\$AAposList=$AAposList, \$tpPanel=$tpPanel, \$idKeyHere=$idKeyHere, \$colorHere=$colorHere, \$keyHere=$keyHere\n";  print   "\$AAposList=$AAposList, \$tpPanel=$tpPanel, \$idKeyHere=$idKeyHere, \$colorHere=$colorHere, \$keyHere=$keyHere\n";
    
  #先建立各个分段的feature
  #$e1 = Bio::SeqFeature::Lite->new(-start=>1,-stop=>100,-type=>'exon');
  my $allSubFeatures;
  for (my $i=0; $i< @{ $AAposList }; $i++){   print "\$AAposList->[$i]->{'_headPos'}=$AAposList->[$i]->{'_headPos'},\$AAposList->[$i]->{'_tailPos'}=$AAposList->[$i]->{'_tailPos'},\$AAposList->[$i]->{'_SegmentType'}=$AAposList->[$i]->{'_SegmentType'}\n";
    $allSubFeatures->[$i]=Bio::SeqFeature::Lite->new( -start=>$AAposList->[$i]->{'_headPos'},-stop=>$AAposList->[$i]->{'_tailPos'}, -name => $AAposList->[$i]->{'_SegmentType'} );
  }
  
  #然后建立 总Feature
  my $AlignFeature = Bio::SeqFeature::Lite->new (-segments=>$allSubFeatures, -id =>$idKeyHere, -strand       => 1);
  #然后画各个比对结果
  
  $tpPanel->add_track(   $AlignFeature ,

                                           -glyph       => 'segments', #'transcript',
                                           -connector   => 'dashed',
                                           -bgcolor     =>  #$colorHere,
                                                            sub {  #'green',
                                                                  my $feature = shift;  print "\cl PrintSubArrayHash::print_all_sub_array(\$feature)=$feature)\n"; PrintSubArrayHash::print_all_sub_array($feature); print "\cl";
                                                                  warn "InFeature :", $feature->start,"..",$feature->end,"\n"; 
                                                                  my $returnColor= 'green'; 
                                                                  my $segType=$feature->seqname;  warn "\$segType=$segType\n\$colorHere=$colorHere\n\n";
                                                                  
                                                                  if    ($colorHere=~m/\S+/  )          {$returnColor=$colorHere;}
                                                                  elsif ($segType eq '_CDSinCoreExon')  {$returnColor='lightgreen';}
                                                                  elsif ($segType eq '_noCoreExon')     {$returnColor='lightgreen';}
                                                                  #elsif ($segType eq '_BtCDS')          {$returnColor='blue';}
                                                                  elsif ($segType eq '_CorePos')        {$returnColor='red';}
                                                                  elsif ($segType eq '_SECIS')          {$returnColor='orange';}
                                                                  
                                                                  warn "\$segType=$segType, \$returnColor=$returnColor \n";
                                                                  $returnColor;
                                                                },  
                                           
                                           
                                           
                                           -fgcolor     => 'black',
                                           -font2color  => 'red',
                                           -key         => $keyHere,
                                           -bump        =>  +1,
                                           -height      =>  12,
                                           -label       =>  1,  
                                           #-part_labels =>  1,             #显示外显子编号
                                           -description =>  1
                                                          #sub {
                                                          #           my $feature = shift;
                                                          #           return unless $feature->has_tag('description');
                                                          #           my ($description) = $feature->each_tag_value('description');
                                                          #           "\$description=$description";
                                                          #      },
                              
                     );
  return $tpPanel;                 
    


}














sub PrintPngForGeneGroups_old_2be_delete{ #$proteindIDHash, $InIN_BothSecCysHash, $outFileHash, $pngFile
  my ($proteindIDHash, $InIN_BothSecCysHash, $outFileHash, $pngFile)=@_;
  open (OUT, ">$pngFile") or die "3 Cannot create \$pngFile=$pngFile : $!\n";
  my $trackHash; 
  my $panelStart=99999999999999;  my $panelEnd=-99999999999999;                   
  foreach my $SecCys (    sort { $a cmp $b }(   keys(  %{ $proteindIDHash }  )   )    ){ 
    foreach my $ProID (    sort { $a cmp $b }(   keys(  %{ $proteindIDHash->{$SecCys} }  )   )    ){ 
    	if (  ( $InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_PredicFromGenoEstOrMatch'} eq 'm' ) && (  $proteindIDHash->{$SecCys}->{$ProID} eq 'gm' )  ){
    	}
    	else{
        my $wiseResultFile="$outFileHash->{$SecCys}/$InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_GeneWiseResult'}";  print "AAAAAAAAAA\$wiseResultFile=\$outFileHash->{$SecCys}/\$InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_GeneWiseResult'}=$outFileHash->{$SecCys}/$InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_GeneWiseResult'}=$wiseResultFile\n";
        my $wiseRstHASH=&NewGenewisedbPraser($wiseResultFile, $InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_CoreCODONpos'}); 
        #print "\cl\clAAAAAAAAAA999999 PrintSubArrayHash::print_all_sub_array (\$wiseRstHASH)\n\n"; PrintSubArrayHash::print_all_sub_array ($wiseRstHASH); print "\n\n\cl\cl\n\n";  
        my $InputCDS_AA_sementArray=$wiseRstHASH
                                              ->{'_genewiseResultHashAd'}
                                              ->{ $wiseRstHASH->{'_BestScoreWithUhash'}->{'_BestScoreQueryProteinID'} }
                                              ->[ $wiseRstHASH->{'_BestScoreWithUhash'}->{'_BestScoreQuProtNUMBER'} ]
                                              ->{'_EachGeneArrayAd'}
                                              ->[ $wiseRstHASH->{'_BestScoreWithUhash'}->{'_BestScoreGeneIdx'} ]
                                              ->{'_EachGeneExonArray'};
                                              
        
        print "\cl\clAAAAAAAAAA00000 PrintSubArrayHash::print_all_sub_array (\$InputCDS_AA_sementArray)\n\n"; PrintSubArrayHash::print_all_sub_array ($InputCDS_AA_sementArray); print "\n\n\cl\cl\n\n";                                      
        my $strand='+'; if ($InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_SecPosFrame'}<0){$strand='-';}
        my $ExonPosList=&AddTGACoreImform( $InputCDS_AA_sementArray, $InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_CoreCODONpos'}, $strand);
        print "\cl\clAAAAAAAAAA111111 PrintSubArrayHash::print_all_sub_array (\$ExonPosList)\n\n"; PrintSubArrayHash::print_all_sub_array ($ExonPosList); print "\n\n\cl\cl\n\n";
        if (defined($InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_bestSECISstart'})){
        	if ( $InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_bestSECISstart'}=~m/\d+/ ){  print "AAAAAAAAAA2222222\$InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_bestSECISstart'}=$InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_bestSECISstart'}\n";
        		$ExonPosList=&AddSECISinform( $ExonPosList, $InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_bestSECISstart'}, $InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_bestSECISend'}, $InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_bestSECISstrand'}, $InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_BestSECISinform'} );
        	  print "\cl\clAAAAAAAAAA2222222 PrintSubArrayHash::print_all_sub_array (\$ExonPosList)\n\n"; PrintSubArrayHash::print_all_sub_array ($ExonPosList); print "\n\n\cl\cl\n\n";
        	}
        }
        $panelStart=&GetSmallestStart20170214($ExonPosList,$panelStart);
        $panelEnd  =&GetBiggestEnd20170214   ($ExonPosList,$panelEnd);
        
        print "\cl\clAAAAAAAAAA333333 PrintSubArrayHash::print_all_sub_array (\$ExonPosList)\n\n"; PrintSubArrayHash::print_all_sub_array ($ExonPosList); print "\n\n\cl\cl\n\n";
        
        $trackHash->{ $InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_CoreCODONpos'} }->{'_SegmentList'}   =$ExonPosList;
        $trackHash->{ $InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_CoreCODONpos'} }->{'_SegmentKeyName'}="$SecCys $ProID";
        $trackHash->{ $InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_CoreCODONpos'} }->{'_Segmentcolor'}  ='';
        $trackHash->{ $InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_CoreCODONpos'} }->{'_SegmentProID'}  =$ProID;
      }
    }
  }
  my $panelLength=$panelEnd-$panelStart+1;   print "AAAAAAAAAA\$panelStart=$panelStart\t\t\$panelEnd=$panelEnd\t\t\$panelLength=$panelLength\n\n";
  $panelStart=$panelStart- int($panelLength/10); if ($panelStart<1){$panelStart=1;}
  $panelEnd=  $panelEnd  + int($panelLength/10); $panelLength=$panelEnd-$panelStart+1;
  
  #my $contigNm=$InIN_BothSecCysHash->{$SecCys}->{$ProID}->{'_GenomeORestContigID'}; $contigNm=~s/^(\S+)\s*.*$/$1/;  #-key     => $contigNm,
  
  my $panel = Bio::Graphics::Panel->new(
                                             -length    => $panelLength,   #图的全长宽度能表示的碱基数量总长度
                                             -width     => 1000,                #图的实际宽度
                                             -key_style => 'left',
                                             -start     => $panelStart,
                                             -pad_left  => 100,
                                             -pad_right => 100,
                                             
                                       );  

  #然后画 标线
  my $full_length = Bio::SeqFeature::Generic->new(
                                                       -start => $panelStart, 
                                                       -end   => $panelEnd
                                                      );
  # 将标尺的 track加入到panel中。
  $panel->add_track ( $full_length,
                                         -glyph   => 'arrow',
                                         -tick    => 2,
                                         #-key     => $contigNm,
                                         -fgcolor => 'black',
                                         -double  => 1
                    );
                    
  foreach my $eachPos (    sort { $a <=> $b }(   keys(  %{ $trackHash }  )   )    ){ 
    
    $panel=&NEWAddFeatureTrack3Panel( 
                                      $trackHash->{ $eachPos }->{'_SegmentList'}, 
                                      $panel,
                                      $trackHash->{ $eachPos }->{'_SegmentKeyName'}, 
                                      $trackHash->{ $eachPos }->{'_Segmentcolor'}, 
                                      ''
                                    );

  } 
  print OUT  $panel->png; 
}


















1;

##########################################################################################################################################
# 



