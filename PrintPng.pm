#!/usr/bin/perl -w

BEGIN{  push (@INC, '/home/fredjiang/fredProgs/PMs')  };

use strict;
use warnings;
use ClustalwRun;
use GD;
use Bio::Graphics::Panel;


use Storable ;
use Storable qw(nstore store_fd nstore_fd freeze thaw dclone);

use DirFileHandle;

package PrintPng;


sub GetStartPos_for_sgmArray{  #  PrintPng::GetStartPos_for_sgmArray
	my ($inSgmArray)=@_;
	
	my $warnMsgBody="\nIn package  PrintPng,\tIn sub GetStartPos_for_sgmArray,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	
	my $Stt_pos;
	my $last_SgmtZorF;
	if ( ref ($inSgmArray) eq 'ARRAY' ){
	  for (my $i=0; $i<@{ $inSgmArray }; $i++){ 
	  	my $SgmtZorF; $SgmtZorF=$inSgmArray->[$i]->{'2_SgmtZorF'} if (  defined ( $inSgmArray->[$i]->{'2_SgmtZorF'} )  );	  	
	  	if (  ($i>0) && ( ($SgmtZorF eq '+') || ($SgmtZorF eq '-') ) && ($SgmtZorF ne $last_SgmtZorF)  ){
	  		my $inSgmArrayDump=DirFileHandle::ReturnDumperInform($inSgmArray);
	  		my $dieMsg="$die_MsgHead\n\$last_SgmtZorF=$last_SgmtZorF \$SgmtZorF=$SgmtZorF\tThese two should be the same!\$inSgmArrayDump=inSgmArrayDump\n\n\n";
	  	  print $dieMsg; die $dieMsg;
	  	}
	  		  	
      if ( ($SgmtZorF eq '+') || ($SgmtZorF eq '-') ){ $last_SgmtZorF=$SgmtZorF; }
    }
    if ($last_SgmtZorF eq '-')   {$Stt_pos=&GetBiggest_Number_for_sgmArray($inSgmArray); return $Stt_pos; }
    else                         {$Stt_pos=&GetSmallestNumber_for_sgmArray($inSgmArray); return $Stt_pos; }
  }
  else {
  	my $dieMsg="$die_MsgHead\nThe \$inSgmArray=$inSgmArray should be a ARRAY ref!\n\n\n"; 
  	print $dieMsg; die $dieMsg;
  }
	
	
}

sub GetEnd__Pos_for_sgmArray{   #  PrintPng::GetEnd__Pos_for_sgmArray
	
	my ($inSgmArray)=@_;
	
	my $warnMsgBody="\nIn package  PrintPng,\tIn sub GetEnd__Pos_for_sgmArray,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	
	my $End_pos;
	my $last_SgmtZorF;
	if ( ref ($inSgmArray) eq 'ARRAY' ){
	  for (my $i=0; $i<@{ $inSgmArray }; $i++){ 
	  	my $SgmtZorF; $SgmtZorF=$inSgmArray->[$i]->{'2_SgmtZorF'} if (  defined ( $inSgmArray->[$i]->{'2_SgmtZorF'} )  );	  	
	  	if (  ($i>0) && ( ($SgmtZorF eq '+') || ($SgmtZorF eq '-') ) && ($SgmtZorF ne $last_SgmtZorF)  ){
	  		my $inSgmArrayDump=DirFileHandle::ReturnDumperInform($inSgmArray);
	  		my $dieMsg="$die_MsgHead\n\$last_SgmtZorF=$last_SgmtZorF \$SgmtZorF=$SgmtZorF\tThese two should be the same!\$inSgmArrayDump=inSgmArrayDump\n\n\n";
	  	  print $dieMsg; die $dieMsg;
	  	}
	  		  	
      if ( ($SgmtZorF eq '+') || ($SgmtZorF eq '-') ){ $last_SgmtZorF=$SgmtZorF; }
    }
    if ($last_SgmtZorF eq '-')   {$End_pos=&GetSmallestNumber_for_sgmArray($inSgmArray); return $End_pos; }
    else                         {$End_pos=&GetBiggest_Number_for_sgmArray($inSgmArray); return $End_pos; }
  }
  else {
  	my $dieMsg="$die_MsgHead\nThe \$inSgmArray=$inSgmArray should be a ARRAY ref!\n\n\n"; 
  	print $dieMsg; die $dieMsg;
  }
	
}


sub GetSmallestNumber_for_sgmArray{    #  PrintPng::GetSmallestNumber_for_sgmArray
	my ($inSgmArray)=@_;
	
	my $warnMsgBody="\nIn package  PrintPng,\tIn sub GetSmallestNumber_for_sgmArray,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	
	my $SmallestNumber=99999999999999;
	
	if ( ref ($inSgmArray) eq 'ARRAY' ){
	  for (my $i=0; $i<@{ $inSgmArray }; $i++){ 
      if ($inSgmArray->[$i]->{'0_SgmtHead'}<=$SmallestNumber){$SmallestNumber=$inSgmArray->[$i]->{'0_SgmtHead'};}
      if ($inSgmArray->[$i]->{'1_SgmtTail'}<=$SmallestNumber){$SmallestNumber=$inSgmArray->[$i]->{'1_SgmtTail'};}    
    }
  }
  else {
  	my $dieMsg="$die_MsgHead\nThe \$inSgmArray=$inSgmArray should be a ARRAY ref!\n\n\n"; 
  	print $dieMsg; die $dieMsg;
  }
  return $SmallestNumber;
}

sub GetBiggest_Number_for_sgmArray{     #  PrintPng::GetBiggest_Number_for_sgmArray
	my ($inSgmArray)=@_;
	
	my $warnMsgBody="\nIn package  PrintPng,\tIn sub GetBiggest_Number_for_sgmArray,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	
	my $Biggest_Number=-99999999999999;
	
	if ( ref ($inSgmArray) eq 'ARRAY' ){
	  for (my $i=0; $i<@{ $inSgmArray }; $i++){ 
      if ($inSgmArray->[$i]->{'0_SgmtHead'}>=$Biggest_Number){$Biggest_Number=$inSgmArray->[$i]->{'0_SgmtHead'};}
      if ($inSgmArray->[$i]->{'1_SgmtTail'}>=$Biggest_Number){$Biggest_Number=$inSgmArray->[$i]->{'1_SgmtTail'};}    
    }
  }
  else {
  	my $dieMsg="$die_MsgHead\nThe \$inSgmArray=$inSgmArray should be a ARRAY ref!\n\n\n"; 
  	print $dieMsg; die $dieMsg;
  }
  return $Biggest_Number;
}

#PrintPng::PrintPngForSgmGroupArrays($allGeneStructureArrayHash->{'0_OnlyGeno'}, $Png_0_OnlyGeno_File, $nextRound_inform_Hash->{$ecSelPro}->{'analysisTypeHash'}->{'OneFiledGeno'}->{'SelProContigNM'}, $ecSelPro );
    	        			
sub PrintPngForSgmGroupArrays{  # PrintPng::PrintPngForSgmGroupArrays
	my ($InSmgArrayHASH, $PngOutFile, $GenoOrEstCtgNm, $SelProName)=@_;
	
	my $warnMsgBody="\nIn package  PrintPng,\tIn sub PrintPngForSgmGroupArrays,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	
	my $tempPngFile_1="${PngOutFile}.1.png";
  my $tempPngFile_2="${PngOutFile}.2.png";
	
	my $out_1; my $out_2; 
  open ($out_1,">$tempPngFile_1") or die "cannot create \$tempPngFile_1=$tempPngFile_1 : $!";
  open ($out_2,">$tempPngFile_2") or die "cannot create \$tempPngFile_2=$tempPngFile_2 : $!";
  
  my $InSmg_united_Array; my $idx=0;
  if  ( ref ($InSmgArrayHASH) eq 'HASH' ){
  	foreach my $eachKey (    sort {$a cmp $b} (   keys (  %{ $InSmgArrayHASH }  )   )    ){
  		$InSmg_united_Array->[$idx]=$InSmgArrayHASH->{$eachKey};
  		$idx++;
  	}
  }
  elsif ( ref ($InSmgArrayHASH) eq 'ARRAY' ){
  	$InSmg_united_Array=$InSmgArrayHASH;
  }
  else {
  	my $dieMsg="$die_MsgHead\n\$InSmgArrayHASH=$InSmgArrayHASH should be a ARRAY or HASH ref!!!\n\n";
  	print $dieMsg; die $dieMsg;
  }
  
  
  my $trackArray; 
  my $panelStart=99999999999999;  my $panelEnd=-99999999999999; 
  if  ( ref ($InSmg_united_Array) eq 'ARRAY' ){
  	for (  my $i=0; $i<@{ $InSmg_united_Array }; $i++  ){
  		if (  defined ( $InSmg_united_Array->[$i]->{'1_SgmArray'} )  ){ my $allSgm_Sm_NB=&GetSmallestNumber_for_sgmArray( $InSmg_united_Array->[$i]->{'1_SgmArray'} ); if ($allSgm_Sm_NB < $panelStart  ){ $panelStart=$allSgm_Sm_NB;} }
  		if (  defined ( $InSmg_united_Array->[$i]->{'1_SgmArray'} )  ){ my $allSgm_Bg_NB=&GetBiggest_Number_for_sgmArray( $InSmg_united_Array->[$i]->{'1_SgmArray'} ); if ($panelEnd     < $allSgm_Bg_NB){ $panelEnd  =$allSgm_Bg_NB;} } 
  		
  		$trackArray->[ $i ]->{'0_SgmInfom'}   =$InSmg_united_Array->[$i]->{'0_SgmInfom'};
      $trackArray->[ $i ]->{'1_SgmArray'}   =$InSmg_united_Array->[$i]->{'1_SgmArray'};
      #$trackArray->[ $i ]->{'2_'}  ='';
       
  		 		
  	}
  }
  
  
  my $panelLength=$panelEnd-$panelStart+1;                
  $panelStart=$panelStart- int($panelLength/10); if ($panelStart<1){$panelStart=1;}
  $panelEnd=  $panelEnd  + int($panelLength/10); $panelLength=$panelEnd-$panelStart+1;
  
  #首先画 panel_1
  my $panel_1 = Bio::Graphics::Panel->new(
                                             -length    => $panelLength,        #图的全长宽度能表示的碱基数量总长度
                                             -width     => 400,                #图的实际宽度
                                             -key_style => 'right',
                                             -start     => 1,
                                             -pad_left  => 100,
                                             -pad_right => 100,
                                             -truecolor =>1,
                                             -truetype => 'Times New Roman',
                                             #-key_font => 'Times New Roman-12:Bold' 
                                             
                                         );  
  
                                       
  #画图示                                      
  if (  ( defined ($SelProName) ) && ( $SelProName=~m/\S+/ )  ){    	
    

    
      
    $panel_1->add_track ( generic => Bio::SeqFeature::Generic->new(  -start=> ($panelLength)*40/100,
                                                                     -end  => ($panelLength)*60/100,
                                                                   
                                                                    ),
                          -glyph  => 'generic',
                          -bgcolor=> 'white',
                          -fgcolor=> 'white',
                          -key    => $SelProName,
                          -label  => 1,
                          -height => 50,
                          
                          
                    );
      
      
    
  }
  
  
  #################################################
  my $panel_2 = Bio::Graphics::Panel->new(
                                             -length    => $panelLength,   #图的全长宽度能表示的碱基数量总长度
                                             -width     => 1000,                #图的实际宽度
                                             -key_style => 'left',
                                             -start     => $panelStart,
                                             -pad_left  => 100,
                                             -pad_right => 100,
                                             -truecolor=>1,
                                             -truetype => 'Times New Roman',
                                             
                                          );  

  #然后画 标线
  my $full_length = Bio::SeqFeature::Generic->new(
                                                       -start => $panelStart, 
                                                       -end   => $panelEnd
                                                      );
  # 将标尺的 track加入到panel中。
  $panel_2->add_track ( $full_length,
                                         -glyph   => 'arrow',
                                         -tick    => 2,
                                         -key     => $GenoOrEstCtgNm,
                                         -fgcolor => 'black',
                                         -double  => 1
                      );
                      
  my $idxH=0; my $tkCol='withe';                  
  for ( my $idxH=0; $idxH<@{ $trackArray }; $idxH++ ){ 
    if ($idxH>=0){
      $panel_2->add_track ( $full_length,
                                         -glyph   => 'line',                                         
                                         -bgcolor => 'whitesmoke',
                                         -fgcolor => 'whitesmoke',
                                         -height  => 3,
                                         
                          ); 
    }
    
    if (( $idxH % 2)==0) { $tkCol='white';        } 
    else                 { $tkCol='white';        }   #ghostwhite

    $panel_2=&NEWAddFeatureTrack_into_Panel_20180904(                                                     
                                                      $panel_2,                                              #1
                                                      $trackArray->[$idxH]->{'0_SgmInfom'},                #2
                                                      $trackArray->[$idxH]->{'1_SgmArray'},                #3                                                        
                                                      $tkCol,                                              #4
                                                      
                                                    ); 
                                  
      
    
    
  } 
  ##########################################
  #$panel_1=&Add_Top_Blank_to_img ($panel_1, 50);
  
  print $out_1  $panel_1->png;   #  $tempPngFile_1
  print $out_2  $panel_2->png;   #  $tempPngFile_2
  
  close ($out_1);
  close ($out_2);

  ClustalwRun::MergeImage ($tempPngFile_1, $tempPngFile_2, $PngOutFile);
  
  &Add_Top_Blank_to_imgFile($PngOutFile, 50);
  
  if (  -e ($tempPngFile_1)  ) {system ("rm -rf $tempPngFile_1"); }
  if (  -e ($tempPngFile_2)  ) {system ("rm -rf $tempPngFile_2"); }
  

    
	
}


sub NEWAddFeatureTrack_into_Panel_20180904 {    # PrintPng::NEWAddFeatureTrack_into_Panel_20180904
  my ( 
       $tpPanel,                         
       $idKeyHere,                       
       $inSgmArray,                      
       $tkBgCol,                         

       )=@_;  
  
  my $warnMsgBody="\nIn package  PrintPng,\tIn sub NEWAddFeatureTrack_into_Panel_20180904,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
 
    
  
  my $CdsStrand;  my $dush_or_none='dashed'; my $CDS_SECIS_type=0;
  
  
  
  #先建立各个分段的feature                                                                                                                                                          
  my $allSubFeatures; 
  
  #my $disPlaName;
  if (ref($inSgmArray) eq 'ARRAY'){
    for (  my $i=0; $i< @{ $inSgmArray }; $i++  ){
      my $strand=1; 
      if (    (  defined ( $inSgmArray->[$i]->{'2_SgmtZorF'} )  ) && (   ( $inSgmArray->[$i]->{'2_SgmtZorF'} eq '-' ) || (  ( $inSgmArray->[$i]->{'2_SgmtZorF'}=~m/\d+/ ) && ( $inSgmArray->[$i]->{'2_SgmtZorF'} < 0 )  )   )    ){$strand=-1;}  
      $CdsStrand=$strand;
      
      my $SgmtType_SgmtInfm="SgmtType\t$inSgmArray->[$i]->{'3_SgmtType'}" if (  defined ( $inSgmArray->[$i]->{'3_SgmtType'} )  ); 
      $SgmtType_SgmtInfm="SgmtType\t$inSgmArray->[$i]->{'3_SgmtType'}:::SgmtInfm\t$inSgmArray->[$i]->{'5_SgmtInfm'}" if (   (  defined ( $inSgmArray->[$i]->{'3_SgmtType'} )  ) && (  defined ( $inSgmArray->[$i]->{'5_SgmtInfm'} )  )   ); 
       
      #$disPlaName.=$inSgmArray->[$i]->{'5_SgmtInfm'} if (   (  defined ( $inSgmArray->[$i]->{'3_SgmtType'} )  ) && (  defined ( $inSgmArray->[$i]->{'5_SgmtInfm'} )  ) && (  $inSgmArray->[$i]->{'3_SgmtType'} eq 'SECIS'  )   ); 
      
      $allSubFeatures->[$i]=Bio::SeqFeature::Lite->new( -start   => $inSgmArray->[$i]->{'0_SgmtHead'}, 
                                                        -end     => $inSgmArray->[$i]->{'1_SgmtTail'}, 
                                                        -name    => $SgmtType_SgmtInfm, 
                                                        -strand  => $strand                                
                                                        );    
      if (      (  defined ( $inSgmArray->[$i]->{'3_SgmtType'} )  ) 
            &&  (     ( $inSgmArray->[$i]->{'3_SgmtType'} eq '0_Geo_E_Ctg'       )
                  ||  ( $inSgmArray->[$i]->{'3_SgmtType'} eq '1_Geo_N_Ctg'       )
                  ||  ( $inSgmArray->[$i]->{'3_SgmtType'} eq '0_Geo_E_Ctg_E_Est' )
                  ||  ( $inSgmArray->[$i]->{'3_SgmtType'} eq '1_Geo_E_Est_N_Ctg' )
                  ||  ( $inSgmArray->[$i]->{'3_SgmtType'} eq '2_Geo_N_Ctg_E_Est' )
                  ||  ( $inSgmArray->[$i]->{'3_SgmtType'} eq '3_Geo_E_Ctg_N_Est' )
                  ||  ( $inSgmArray->[$i]->{'3_SgmtType'} eq '4_Geo_N_Ctg_N_Est' )
                  ||  ( $inSgmArray->[$i]->{'3_SgmtType'} eq '0_Ctg_E_Est'       )
                  ||  ( $inSgmArray->[$i]->{'3_SgmtType'} eq '1_Ctg_N_Est'       )
                  ||  ( $inSgmArray->[$i]->{'3_SgmtType'} eq '5_OTHER_TYP'       )
                  ||  ( $inSgmArray->[$i]->{'3_SgmtType'} eq 'Sec_TGA'           )
                  ||  ( $inSgmArray->[$i]->{'3_SgmtType'} eq 'frame_shfit_insert')
                ) 
        ) {$dush_or_none='none';}
      
      if (   (  defined ( $inSgmArray->[$i]->{'3_SgmtType'} )  ) && (  ( $inSgmArray->[$i]->{'3_SgmtType'} eq 'CDS_region' ) || ( $inSgmArray->[$i]->{'3_SgmtType'} eq 'SECIS' )  )   ){
      	$CDS_SECIS_type=1;
      }
          
                                                              
    }
  }
  if ($CDS_SECIS_type==1) { $dush_or_none='dashed'; }
  
  
  #然后建立 总Feature
  my $All_Feature = Bio::SeqFeature::Lite->new (-segments => $allSubFeatures,
                                                 -key        => $idKeyHere,  
                                                 #-display_name =>$disPlaName, 
                                                 -id       => "",#$idKeyHere, 
                                                 -strand   => $CdsStrand   #$newinSgmArray[0]->{'_Strand'} 
                                                 );                                                 
  #然后画各个  结果
  
  
  $tpPanel->add_track(   $All_Feature ,
                                           -key        =>  $idKeyHere,                
                                           -glyph       => 'transcript',
                                           -tkcolor     =>  $tkBgCol, #'whitesmoke', #$colorIdx, #'#FFF8DC',#'lightyellow',#'cyan',#'gray',
                                           #-box_subparts=> 3,
                                           -strand      => $CdsStrand,
                                           -connector   => $dush_or_none,#'none',,#'dashed',
                                           -connector_color   => 'black',
                                           -bgcolor     =>  #$colorHere,
                                                            sub {  #'green',
                                                                  my $feature = shift;                                          #print "\cl 111111111 PrintSubArrayHash::print_all_sub_array(\$feature)=$feature)\n"; PrintSubArrayHash::print_all_sub_array($feature); print "\cl";
                                                                                                                                # warn "\n\nInFeature :", $feature->start,"..",$feature->end,"\n"; 
                                                                  my $returnColor= 'green'; 
                                                                  my $segType=$feature->seqname;                              #warn "\n\$segType=$segType\n";# \t\$colorHere=$colorHere\n";
                                                                  if ($segType=~m/SgmtType\t(\S+)(?:\:\:\:SgmtInfm\t(.*))?$/){    
                                                                  	my $inFetuSgmtType=$1;
                                                                  	my $inFetuSgmtInfm=$2;
                                                                  	
                                                                  
                                                                    if    ($inFetuSgmtType eq 'CDS_region'        )  {$returnColor='lightgreen';}
                                                                    elsif ($inFetuSgmtType eq 'Sec_TGA'           )  {$returnColor='red'       ;}
                                                                    elsif ($inFetuSgmtType eq 'frame_shfit_insert')  {$returnColor='blue'      ;}
                                                                    elsif ($inFetuSgmtType eq 'SECIS'             )  {$returnColor='orange'    ;}
                                                                    elsif ($inFetuSgmtType eq '0_Geo_E_Ctg'       )  {$returnColor='brown'     ;}
                                                                    elsif ($inFetuSgmtType eq '1_Geo_N_Ctg'       )  {$returnColor='black'     ;}
                                                                    elsif ($inFetuSgmtType eq '0_Geo_E_Ctg_E_Est' )  {$returnColor='pink'      ;}
                                                                    elsif ($inFetuSgmtType eq '1_Geo_E_Est_N_Ctg' )  {$returnColor='gray'      ;}
                                                                    elsif ($inFetuSgmtType eq '2_Geo_N_Ctg_E_Est' )  {$returnColor='gray'      ;}
                                                                    elsif ($inFetuSgmtType eq '3_Geo_E_Ctg_N_Est' )  {$returnColor='gray'      ;}
                                                                    elsif ($inFetuSgmtType eq '4_Geo_N_Ctg_N_Est' )  {$returnColor='black'     ;}
                                                                    elsif ($inFetuSgmtType eq '0_Ctg_E_Est'       )  {$returnColor='pink'      ;}
                                                                    elsif ($inFetuSgmtType eq '1_Ctg_N_Est'       )  {$returnColor='black'     ;}
                                                                    elsif ($inFetuSgmtType eq '5_OTHER_TYP'       )  {$returnColor='white'     ;}
                                                                   
                                                                  }                                                                
                                                                  
                                                                  $returnColor;
                                                                },  
                                           
                                           
                                           
                                           -fgcolor     =>                       #'',#'yellow',    -nothing =>
                                                           sub {  #'green',
                                                                  my $feature = shift;                                          #print "\cl 111111111 PrintSubArrayHash::print_all_sub_array(\$feature)=$feature)\n"; PrintSubArrayHash::print_all_sub_array($feature); print "\cl";
                                                                                                                                # warn "\n\nInFeature :", $feature->start,"..",$feature->end,"\n"; 
                                                                  my $returnColor= 'green'; 
                                                                  my $segType=$feature->seqname;                              #warn "\n\$segType=$segType\n";# \t\$colorHere=$colorHere\n";
                                                                  if ($segType=~m/SgmtType\t(\S+)(?:\:\:\:SgmtInfm\t(.*))?$/){    
                                                                  	my $inFetuSgmtType=$1;
                                                                  	my $inFetuSgmtInfm=$2;
                                                                  	
                                                                  
                                                                    if    ($inFetuSgmtType eq 'CDS_region'        )  {$returnColor='lightgreen';}
                                                                    elsif ($inFetuSgmtType eq 'Sec_TGA'           )  {$returnColor='red'     ;}
                                                                    elsif ($inFetuSgmtType eq 'frame_shfit_insert')  {$returnColor='blue'       ;}
                                                                    elsif ($inFetuSgmtType eq 'SECIS'             )  {$returnColor='orange'    ;}
                                                                    elsif ($inFetuSgmtType eq '0_Geo_E_Ctg'       )  {$returnColor='brown'     ;}
                                                                    elsif ($inFetuSgmtType eq '1_Geo_N_Ctg'       )  {$returnColor='black'     ;}
                                                                    elsif ($inFetuSgmtType eq '0_Geo_E_Ctg_E_Est' )  {$returnColor='pink'      ;}
                                                                    elsif ($inFetuSgmtType eq '1_Geo_E_Est_N_Ctg' )  {$returnColor='gray'      ;}
                                                                    elsif ($inFetuSgmtType eq '2_Geo_N_Ctg_E_Est' )  {$returnColor='gray'      ;}
                                                                    elsif ($inFetuSgmtType eq '3_Geo_E_Ctg_N_Est' )  {$returnColor='gray'      ;}
                                                                    elsif ($inFetuSgmtType eq '4_Geo_N_Ctg_N_Est' )  {$returnColor='black'     ;}
                                                                    elsif ($inFetuSgmtType eq '0_Ctg_E_Est'       )  {$returnColor='pink'      ;}
                                                                    elsif ($inFetuSgmtType eq '1_Ctg_N_Est'       )  {$returnColor='black'     ;}
                                                                    elsif ($inFetuSgmtType eq '5_OTHER_TYP'       )  {$returnColor='white'     ;}
                                                                   
                                                                  }                                                                
                                                                  
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
                                           -part_labels =>0,  #
                                                           #sub {
                                                           #      my ($feature,undef,$partno) = @_; 
                                                           #      my $tp=$feature->{segments}->[$partno]->name; warn "\n\$feature->{segments}->[$partno]->name=$tp\n\@_=@_\n\n\n"; #print "\@_=@_\n\n\cl 2222222222 PrintSubArrayHash::print_all_sub_array(\$feature)=$feature)\n"; PrintSubArrayHash::print_all_sub_array($feature); print "\cl";
                                                           #      my $legth=$feature->{segments}->[$partno]->length;
                                                           #      #if (defined ($_[4])){print "\@_=@_\n\n\cl 2222222222 PrintSubArrayHash::print_all_sub_array(\$_[4]=$_[4])\n"; PrintSubArrayHash::print_all_sub_array($_[4]); print "\cl";}
                                                           #      my $returnLable="Psd";
                                                           #      if      ($tp eq '_CorePos'){ $returnLable=$CorePosCode; }
                                                           #      elsif   ($tp eq '_SECIS'){ $returnLable=$SECISinform; }
                                                           #      elsif   ($tp eq '_noCoreExon'){ $returnLable=" "; }
                                                           #      elsif   ($tp eq '_BeforCore'){ $returnLable=" "; }   
                                                           #      elsif   ($tp eq '_AfterCore'){ $returnLable=" "; } 
                                                           #      elsif   ($tp eq '_PseudoPos'){ $returnLable=" $legth bp idel"; } 
                                                           #     
                                                           #      	$returnLable;
                                                           #      #return $returnLable;
                                                           #    }
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

sub Add_Top_Blank_to_img{   #PrintPng::Add_Top_Blank_to_img
	my ($imgObj, $add_pix_nb)=@_;
	
	my $warnMsgBody="\nIn package  PrintPng,\tIn sub Add_Top_Blank_to_img,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";

	my ( $in_x,  $in_y  ) = $imgObj->getBounds();
	
	my $fx=0; my $fy=0;
	my $imF = GD::Image->new( $fx, $fy, 0 );                 
  
  $imF->copy( $imgObj , 0, $add_pix_nb, 0, 0, $in_x,  $in_y );

  return $imF;

	
}

sub Add_Top_Blank_to_imgFile{   #PrintPng::Add_Top_Blank_to_imgFile
	my ($img_file, $add_pix_nb)=@_;
	
	my $warnMsgBody="\nIn package  PrintPng,\tIn sub Add_Top_Blank_to_imgFile,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $imgObj= GD::Image->newFromPng( $img_file ) or die"$die_MsgHead\n\$img_file=$img_file $!\n"; 
  
  my $fx=0; my $fy=0;
	my ( $in_x,  $in_y  ) = $imgObj->getBounds(); $fx=$in_x;  $fy=$in_y+$add_pix_nb;
	
	
	my $imF = GD::Image->new( $fx, $fy, 0 );                 
  
  $imF->copy( $imgObj, 0, $add_pix_nb, 0, 0, $in_x,  $in_y );

  open ALLPNG, ">$img_file" or die "$die_MsgHead cannot whrite into \$img_file=$img_file\n";
  print ALLPNG $imF->png;
  close ALLPNG;

	
}


sub MergeImage_from_file_img1_up_img2_down{  #PrintPng::MergeImage_from_file_img1_up_img2_down #将png文件$image_1和$image_2,合并成$image_out, 应该是上下合并
  my ($image_1, $image_2, $image_out)=@_;
  
  my $warnMsgBody="\nIn package  PrintPng,\tIn sub MergeImage_from_file_img1_up_img2_down,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	
  #warn "\$image_1=$image_1\t\$image_2=$image_2\t\$image_out=$image_out\n";
  
  my $fx=0; my $fy=0; my $dist=10;  my $eOn;
  my $im1 ; my( $x1,  $y1  ); if ( -e($image_1) ){ $eOn->{1 }=1; my $dmsg_1="$die_MsgHead \$image_1=$image_1\n\n$!\n"; $im1 = GD::Image->newFromPng( $image_1 ) or &PrintAndDIe ($dmsg_1);	( $x1,  $y1  ) = $im1->getBounds() ; if ($fx<=$x1 ){$fx=$x1 ;}$fy+=($y1 +$dist); }
  my $im2 ; my( $x2,  $y2  ); if ( -e($image_2) ){ $eOn->{2 }=1; my $dmsg_2="$die_MsgHead \$image_2=$image_2\n\n$!\n"; $im2 = GD::Image->newFromPng( $image_2 ) or &PrintAndDIe ($dmsg_2);  ( $x2,  $y2  ) = $im2->getBounds() ; if ($fx<=$x2 ){$fx=$x2 ;}$fy+=($y2 +$dist); }
  #warn "\$fx=$fx; \$fy=$fy\n";
  
  my $imF = GD::Image->new( $fx, $fy, 0 );       my $startY=0;         my $dist0=1;               
  if (    (   defined(  $eOn->{ 1 }  )   ) && (  $eOn->{ 1 }==1 )    ){  $imF->copy( $im1 , 0, $startY, 0, 0, $x1 , $y1 );  $startY+=($y1 +$dist0);  }
  if (    (   defined(  $eOn->{ 2 }  )   ) && (  $eOn->{ 2 }==1 )    ){  $imF->copy( $im2 , 0, $startY, 0, 0, $x2 , $y2 );  $startY+=($y2 +$dist );  }
  #warn "\$startY=$startY\n";
                            
  open ALLPNG, ">$image_out" or die "$die_MsgHead cannot create \$image_out=$image_out\n";
  print ALLPNG $imF->png;
  close ALLPNG;
                          
}


sub MergeImage_from_file_img1_left_img2_right{   #PrintPng::MergeImage_from_file_img1_left_img2_right #将png文件$image_1和$image_2,合并成$image_out, 应该是左右合并
  my ($image_1, $image_2, $image_out)=@_;
  
  my $warnMsgBody="\nIn package  PrintPng,\tIn sub MergeImage_from_file_img1_left_img2_right,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	
  #warn "\$image_1=$image_1\t\$image_2=$image_2\t\$image_out=$image_out\n";
  #print "\n20181031 \$image_1=$image_1\t\$image_2=$image_2\t\$image_out=$image_out\n";
  
  my $fx=0; my $fy=0; my $dist=50;  my $eOn;
  my $im1 ; my( $x1,  $y1  ); if ( -e($image_1) ){  $eOn->{1 }=1; my $dmsg_1="$die_MsgHead \$image_1=$image_1\n\n$!\n"; $im1 = GD::Image->newFromPng( $image_1 ) or &PrintAndDIe ($dmsg_1);	 	( $x1,  $y1  ) = $im1->getBounds() ; if ($fy<=$y1 ){$fy=$y1 ;} $fx+=($x1 +$dist); }
  my $im2 ; my( $x2,  $y2  ); if ( -e($image_2) ){  $eOn->{2 }=1; my $dmsg_2="$die_MsgHead \$image_2=$image_2\n\n$!\n"; $im2 = GD::Image->newFromPng( $image_2 ) or &PrintAndDIe ($dmsg_2);  	( $x2,  $y2  ) = $im2->getBounds() ; if ($fy<=$y2 ){$fy=$y2 ;} $fx+=($x2 +$dist); }
  #print "\n20181031 \$fx=$fx; \$fy=$fy\n"; warn "\$fx=$fx; \$fy=$fy\n";
  
  my $imF = GD::Image->new( $fx, $fy, 0 );       my $startX=0;         my $dist0=1;               
  if (    (   defined(  $eOn->{ 1 }  )   ) && (  $eOn->{ 1 }==1 )    ){  $imF->copy( $im1 , $startX, 0, 0, 0, $x1 , $y1 );  $startX+=($x1 +$dist0);  }
  if (    (   defined(  $eOn->{ 2 }  )   ) && (  $eOn->{ 2 }==1 )    ){  $imF->copy( $im2 , $startX, 0, 0, 0, $x2 , $y2 );  $startX+=($x2 +$dist );  }
  #print "\n20181031 \$startX=$startX\n"; warn "\$startX=$startX\n";
                            
  open ALLPNG, ">$image_out" or die "$die_MsgHead cannot create \$image_out=$image_out\n";
  print ALLPNG $imF->png;
  close ALLPNG;
                          
}


sub PrintAndDIe{
	my ($printInform)=@_;
	print $printInform;
	die $printInform;
}


sub PrintPNG_files_20190111{  #  PrintPng::PrintPNG_files_20190111 ($allSgmHash, $Geno_SECIS_inform_Array, $Est__SECIS_inform_Array, $pngFile_1, $pngFile_2, $pngFile_3, $Geno_ContigName, $Est__ContigName,  $SelPro_ID);
	my ($allSgmHash, $Geno_SECIS_inform_Array, $Est__SECIS_inform_Array, $pngFile_1, $pngFile_2, $pngFile_3, $Geno_ContigName, $Est__ContigName,  $SelPro_ID)=@_;
  
  my $warnMsgBody="\nIn package  PrintPng,\tIn sub PrintPNG_files_20190111,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $callerMSG=DirFileHandle::print_SubCallerInform;

  
  
  	
  if ( ref($allSgmHash) eq 'HASH' ){
  	
  	if ( ref($allSgmHash->{'0_OnlyGeno'}) eq 'HASH' ){
  		PrintPng::PrintPngForSgmGroupArrays($allSgmHash->{'0_OnlyGeno'}, $pngFile_1, $Geno_ContigName, $SelPro_ID );
  		#print "20181031 0_OnlyGeno 2$SelPro_ID\$Geno_SECIS_inform_Array=$Geno_SECIS_inform_Array\n"; 
  		if ( ref ($Geno_SECIS_inform_Array) eq 'ARRAY' ){
  			for (my $i=0; $i<@{ $Geno_SECIS_inform_Array }; $i++){
  				my $SECISpngFile=$Geno_SECIS_inform_Array->[$i]->{'20_pngPathFileKwd'};  #print "20181031 0_OnlyGeno 2$SelPro_ID\$SECISpngFile=$SECISpngFile\n";    	        		
  				#PrintPng::MergeImage_from_file_img1_left_img2_right ($pngFile_1, $SECISpngFile, $pngFile_1);
  			}
  		}
  	}
  	if ( ref($allSgmHash->{'1_Only_Est'}) eq 'HASH' ){
  		PrintPng::PrintPngForSgmGroupArrays($allSgmHash->{'1_Only_Est'}, $pngFile_2, $Est__ContigName, $SelPro_ID );
  		#print "20181031 1_Only_Est 2$SelPro_ID\$Est__SECIS_inform_Array=$Est__SECIS_inform_Array\n"; 
  		if ( ref ($Est__SECIS_inform_Array) eq 'ARRAY' ){
  			for (my $i=0; $i<@{ $Est__SECIS_inform_Array }; $i++){
  				my $SECISpngFile=$Est__SECIS_inform_Array->[$i]->{'20_pngPathFileKwd'}; #print "20181031 1_Only_Est 1 2$SelPro_ID\$SECISpngFile=$SECISpngFile\n";
  				#PrintPng::MergeImage_from_file_img1_left_img2_right ($pngFile_2, $SECISpngFile, $pngFile_2);
  				#print "20181031 1_Only_Est 2 2$SelPro_ID\$SECISpngFile=$SECISpngFile\n";
  			}
  		}
  	}  #print "\n20181031-3 -$SelPro_ID\n"; 
  	if ( ref($allSgmHash->{'2_genoMest'}) eq 'HASH' ){
  		PrintPng::PrintPngForSgmGroupArrays($allSgmHash->{'2_genoMest'}, $pngFile_3, $Geno_ContigName, $SelPro_ID );
  		#print "\n20181031-3 2_genoMest 2$SelPro_ID\$Geno_SECIS_inform_Array=$Geno_SECIS_inform_Array\n"; 
  		if ( ref ($Geno_SECIS_inform_Array) eq 'ARRAY' ){
  			for (my $i=0; $i<@{ $Geno_SECIS_inform_Array }; $i++){
  				my $SECISpngFile=$Geno_SECIS_inform_Array->[$i]->{'20_pngPathFileKwd'}; #print "\n20181031-3 2_genoMest 2$SelPro_ID\$SECISpngFile=$SECISpngFile\n";
  				#PrintPng::MergeImage_from_file_img1_left_img2_right ($pngFile_3, $SECISpngFile, $pngFile_3);
  			}
  		}
  		
  	}
  	
  }
    	        		

}


1;