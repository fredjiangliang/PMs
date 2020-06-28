#!/usr/bin/perl -w
BEGIN{  push (@INC, '/home/fredjiang/fredProgs/PMs')  };


use strict;
use warnings;

use DirFileHandle;

############################################################################################################

#         
############################################################################################################



                      
package SeqSegmentsTools;


sub a_between_b_c{      #  my $yes_or_no=SeqSegmentsTools::a_between_b_c ($a,$b,$c); #pmAble#   #按顺序输入三个数字 a b c ，检测 a是否在b和c之间，包括等于。
  my ($a,$b,$c)=@_; my $yes=0; 
  if ($b<=$c){if ( ($b<=$a) && ($a<=$c) ){$yes=1;}}
  if ($c<=$b){if ( ($c<=$a) && ($a<=$b) ){$yes=1;}}
  return $yes;
}

sub getSmallOne{ # SeqSegmentsTools::getSmallOne
	my ($a,$b)=@_;
  if ($a<=$b){return $a;}
  else       {return $b;}
}

sub getbigOne{ # SeqSegmentsTools::getbigOne
	my ($a,$b)=@_;
  if ($a>=$b){return $a;}
  else       {return $b;}
}

sub OverlayCheck {      #检查，两个数段之间 是否有重叠的地方,重叠则 返回1，修改于2016.03.08,之前错了
  my ($seg_1_hd, $seg_1_tl, $seg_2_hd, $seg_2_tl)=@_;
  my $overlayOrNot=1;
  if   ( ($seg_1_hd<=$seg_2_hd) && ($seg_2_hd<=$seg_1_tl) ){  $overlayOrNot=1;  }
  elsif( ($seg_1_hd<=$seg_2_tl) && ($seg_2_tl<=$seg_1_tl) ){  $overlayOrNot=1;  }
  elsif( ($seg_2_hd<=$seg_1_hd) && ($seg_1_hd<=$seg_2_tl) ){  $overlayOrNot=1;  }
  elsif( ($seg_2_hd<=$seg_1_tl) && ($seg_1_tl<=$seg_2_tl) ){  $overlayOrNot=1;  }
  else                                                     {  $overlayOrNot=0;  }
  return $overlayOrNot;
}

sub a_b_including_c_d {     #  my $yes_or_no=SeqSegmentsTools::a_b_including_c_d ($a,$b,$c,$d); #检查，两个数段, (a,b) 是否 包含 (c,d);
  my ($seg_1_hd, $seg_1_tl, $seg_2_hd, $seg_2_tl)=@_;
  my $overlayOrNot=1;
  if   ( ( &a_between_b_c( $seg_2_hd, $seg_1_hd, $seg_1_tl ) ) && ( &a_between_b_c( $seg_2_tl, $seg_1_hd, $seg_1_tl ) ) ){  $overlayOrNot=1;  }  
  else                                                     {  $overlayOrNot=0;  }
  return $overlayOrNot;
}

sub Find_zoom_nb_c_from_a_b_corespondTo_g_of_e_f{    #   my $c=SeqSegmentsTools::Find_zoom_nb_c_from_a_b_corespondTo_g_of_e_f ($a, $b,    $e, $f, $g );
	my ($a, $b, $e, $f, $g)=@_;
	
	my $warnMsgBody="\nIn package  SeqSegmentsTools,\tIn sub Find_zoom_nb_c_from_a_b_corespondTo_g_of_e_f,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  my $subCallereIfm=DirFileHandle::print_SubCallerInform; #warn $subCallereIfm;
  
	if ($a>$b){ DieWork::Just_dieWork( $die_MsgHead."\$a=$a should < \$b=$b\n\n\n".$subCallereIfm ); }
	if ($e>$f){ DieWork::Just_dieWork( $die_MsgHead."\$e=$e should < \$f=$f\n\n\n".$subCallereIfm ); }
	
	my $c=($b-$a+1)*($g-$e+1)/($f-$e+1)+$a-1;
	
	return $c;
	
}

##  my $outArray=SeqSegmentsTools::BuildOverLayerPosPair_for_a_length($wholeSegLth, $segMentLth, $overLayLth);
sub BuildOverLayerPosPair_for_a_length{  #输入一个 数字1，可看作是一段序列的长度，再输入一个数字2，可看作是切片的长度，再输入一个数字3，可看作是重叠区域的长度，输出一个数组，记录了成对的 片段的开头和结尾
	
	my ($wholeSegLth, $segMentLth, $overLayLth)=@_;
	
	my $warnMsgBody="\nIn package  SeqSegmentsTools,\tIn sub BuildOverLayerPosPair_for_a_length,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  my $caller_inform=DirFileHandle::print_SubCallerInform;
  
	if (   (  defined ( $wholeSegLth  )  ) && ( $wholeSegLth=~/^\d+$/ )  && ( $wholeSegLth > 0 )    ){	}
	else{ DieWork::Just_dieWork( $die_MsgHead."\n \$wholeSegLth=$wholeSegLth should be a int number > 0 !!  $!\n\n\n".$caller_inform );   }
	
	if (   (  defined ( $segMentLth  )  ) && ( $segMentLth=~/^\d+$/ )  && ( $segMentLth > 0 )    ){	}
	else{ DieWork::Just_dieWork( $die_MsgHead."\n \$segMentLth=$segMentLth should be a int number > 0 !!  $!\n\n\n".$caller_inform );   }
	
	if (   (  defined ( $overLayLth  )  ) && ( $overLayLth=~/^\d+$/ )  && ( $overLayLth >= 0 )    ){	}
	else{ DieWork::Just_dieWork( $die_MsgHead."\n \$overLayLth=$overLayLth should be a int number >= 0 !!  $!\n\n\n".$caller_inform );   }
	
	if ( $overLayLth >= $segMentLth ){	
	  DieWork::Just_dieWork( $die_MsgHead."\n \$overLayLth=$overLayLth >= $segMentLth=\$segMentLth, it is not OK !!  $!\n\n\n".$caller_inform );
	}
	
	
	
	my $outArray;
	
	my $headNB=1;
	my $tailNb=0;
	
	my $idx_NB=0;
	
	while ( $headNB <= $wholeSegLth){
		
		$tailNb=$headNB+$segMentLth-1;
		if ($tailNb > $wholeSegLth){
			$tailNb=$wholeSegLth;
		}
		
		$outArray->[$idx_NB]->{'0_SgmtHead'}=$headNB;
		$outArray->[$idx_NB]->{'1_SgmtTail'}=$tailNb;
		
		
		$headNB+=($segMentLth-$overLayLth);
		$idx_NB++;
	}
	
	if (   (  defined ( $outArray )  ) && (  ref ( $outArray ) eq 'ARRAY' )   ){			}
	else{ DieWork::Just_dieWork( $die_MsgHead."\n \$outArray=$outArray should be a ARRAY ref !!  $!\n\n\n".$caller_inform );   }
	
	return $outArray;
	
	
}

sub Get_Overlay_rate {     #  SeqSegmentsTools::Get_Overlay_rate
  #my ($seg_1_hd, $seg_1_tl, $seg_2_hd, $seg_2_tl, $ecSelPro)=@_;       print "9-1-20180908$ecSelPro\t$seg_1_hd, $seg_1_tl, $seg_2_hd, $seg_2_tl\n";
  my ($seg_1_hd, $seg_1_tl, $seg_2_hd, $seg_2_tl)=@_;
  
  my $warnMsgBody="\nIn package  SeqSegmentsTools,\tIn sub Get_Overlay_rate,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $sm_1=&getSmallOne($seg_1_hd, $seg_1_tl); 
  my $bg_1=&getbigOne  ($seg_1_hd, $seg_1_tl); 
  my $lt_1=$bg_1-$sm_1+1;                    
  my $sm_2=&getSmallOne($seg_2_hd, $seg_2_tl);
  my $bg_2=&getbigOne  ($seg_2_hd, $seg_2_tl);
  my $lt_2=$bg_2-$sm_2+1;                    
  
  
  my $overlayOrNot=1;
  if   ( ($sm_1<=$sm_2) && ($sm_2<=$bg_1) ){  $overlayOrNot=1; } #print "9-2-20180908$ecSelPro\t\$overlayOrNot=$overlayOrNot\n"; }          
  elsif( ($sm_1<=$bg_2) && ($bg_2<=$bg_1) ){  $overlayOrNot=1; } #print "9-3-20180908$ecSelPro\t\$overlayOrNot=$overlayOrNot\n"; }          
  elsif( ($sm_2<=$sm_1) && ($sm_1<=$bg_2) ){  $overlayOrNot=1; } #print "9-4-20180908$ecSelPro\t\$overlayOrNot=$overlayOrNot\n"; }          
  elsif( ($sm_2<=$bg_1) && ($bg_1<=$bg_2) ){  $overlayOrNot=1; } #print "9-5-20180908$ecSelPro\t\$overlayOrNot=$overlayOrNot\n"; }          
  else                                     {  $overlayOrNot=0; } #print "9-6-20180908$ecSelPro\t\$overlayOrNot=$overlayOrNot\n"; }          
  my $rateArray; 
  $rateArray->[0]=0;
  $rateArray->[1]=0;
  
  if ($overlayOrNot==1){
        
    my @sortArray=($sm_1, $bg_1, $sm_2, $bg_2);
    @sortArray=sort {$a<=>$b} @sortArray;
    my $overlayerLeth=$sortArray[2]-$sortArray[1]+1;                                    #print "9-7-20180908$ecSelPro\t\$overlayerLeth=\$sortArray[2]-\$sortArray[1]+1=$sortArray[2]-$sortArray[1]+1=$overlayerLeth\n";
    $rateArray->[0]=$overlayerLeth/$lt_1;                                               #print "9-8-20180908$ecSelPro\t\$rateArray->[0]=\$overlayerLeth/\$lt_1=$rateArray->[0]\n";
    $rateArray->[1]=$overlayerLeth/$lt_2;                                               #print "9-8-20180908$ecSelPro\t\$rateArray->[1]=\$overlayerLeth/\$lt_2=$rateArray->[1]\n";
  }
  return $rateArray;
  
  
  
}


sub DirectionDetact{    #pmAble#  #对于只有两个数字的数对，通过比较两个数字的大小，来判断数对可能表示的外显子等 数段的方向， 0正向 1负向 2无方向
  my ($head, $tail)=@_;
  my $direaction=0;
  if    ($head < $tail) { $direaction=0; }
  elsif ($head > $tail) { $direaction=1; }
  elsif ($head == $tail){ $direaction=2; }
  else {die "In Sub  DirectionDetact: \$head=$head, \$tail=$tail\n\n\n";}
  return $direaction;
}

sub mergeSegmentDIGUI{   #这个函数没有使用到
  my ($hd, $tl, $ArrayHash)=@_;
  my $merged=0;
  my $newArrayHash; my $newIdx=0;
  for( my $i=0; $i<@{ $ArrayHash }; $i++ ){
  	if ( SeqSegmentsTools::OverlayCheck($hd, $tl, $ArrayHash->[$i]->{'_start'}, $ArrayHash->[$i]->{'_end'}) ){
  	  $merged=1;
  	  my ($newSt, $newEd)=@{ JointSegments_sameDirection_Up_or_down_noEqual($hd, $tl, $ArrayHash->[$i]->{'_start'}, $ArrayHash->[$i]->{'_end'}) };
  	  $newArrayHash->{$newIdx}->{'_start'}=$newSt;
  	  $newArrayHash->{$newIdx}->{'_end'}  =$newEd;
  	}
  	else {
  	  $newArrayHash->{$newIdx}->{'_start'}=$ArrayHash->[$i]->{'_start'};
  	  $newArrayHash->{$newIdx}->{'_end'}  =$ArrayHash->[$i]->{'_end'};
  	}
    
  } 
}

sub Divide_Sgments_into_merged_Groups{
	my ($in_Sgms_Hash)=@_;
	my $warnMsgBody="\nIn package  SeqSegmentsTools,\tIn sub Divide_Sgments_into_merged_Groups,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  DirFileHandle::print_SubCallerInform;
  
  my $sm_bg_Sgms_Hash;
  if (   (  defined ( $in_Sgms_Hash )  ) && (  ref( $in_Sgms_Hash ) eq 'HASH'  )   ){
    foreach my $Each_sgm_key (    sort {$a cmp $b}(   keys(  %{ $in_Sgms_Hash }  )   )    ){
      $sm_bg_Sgms_Hash->{$Each_sgm_key}=Storable::dclone( $in_Sgms_Hash->{$Each_sgm_key} );
      $sm_bg_Sgms_Hash->{$Each_sgm_key}->{'0_SgmtHead'}=SeqSegmentsTools::getSmallOne( $sm_bg_Sgms_Hash->{$Each_sgm_key}->{'0_SgmtHead'},$sm_bg_Sgms_Hash->{$Each_sgm_key}->{'1_SgmtTail'} );
      $sm_bg_Sgms_Hash->{$Each_sgm_key}->{'1_SgmtTail'}=SeqSegmentsTools::getbigOne  ( $sm_bg_Sgms_Hash->{$Each_sgm_key}->{'0_SgmtHead'},$sm_bg_Sgms_Hash->{$Each_sgm_key}->{'1_SgmtTail'} );
      
    }
  }
  my $divided_group_Array; my $divided_group_Array_idx=0;
  my $lastSgmTail=-9999999999999;
  my $step=0;
  if (   (  defined ( $sm_bg_Sgms_Hash )  ) && (  ref( $sm_bg_Sgms_Hash ) eq 'HASH'  )   ){
    foreach my $Each_sgm_key (    sort { $sm_bg_Sgms_Hash->{$a}->{'0_SgmtHead'} <=> $sm_bg_Sgms_Hash->{$b}->{'0_SgmtHead'} || $sm_bg_Sgms_Hash->{$a}->{'1_SgmtTail'} <=> $sm_bg_Sgms_Hash->{$b}->{'1_SgmtTail'} }(   keys(  %{ $sm_bg_Sgms_Hash }  )   )    ){
      if ($step==0){
      	push @{ $divided_group_Array->[$divided_group_Array_idx] },$Each_sgm_key;
      }
      elsif ($step>0){
      	if     ( $sm_bg_Sgms_Hash->{$Each_sgm_key}->{'0_SgmtHead'} <= $lastSgmTail ){
      		push @{ $divided_group_Array->[$divided_group_Array_idx] },$Each_sgm_key;
      	}
      	elsif  ( $sm_bg_Sgms_Hash->{$Each_sgm_key}->{'0_SgmtHead'} > $lastSgmTail  ){
      		$divided_group_Array_idx++;
      		push @{ $divided_group_Array->[$divided_group_Array_idx] },$Each_sgm_key;
      	}
      	
      }
      $lastSgmTail=$sm_bg_Sgms_Hash->{$Each_sgm_key}->{'1_SgmtTail'} if ( $lastSgmTail<=$sm_bg_Sgms_Hash->{$Each_sgm_key}->{'1_SgmtTail'} );
      $step++;
    }
  }
  return $divided_group_Array;
  
}

#sub              JointSegments_sameDirection_Up_or_down_noEqual{   #如果两个数段是重叠的，则将两个数段进行叠加, 并且检测数段的方向，方向满足则运行，否则die
sub JointSegments_sameDirection_Up_or_down_noEqual{     #pmAble#  #如果两个数段是重叠的，则将两个数段进行叠加, 并且检测数段的方向，方向满足则运行，否则die
  my ($seg_1_hd, $seg_1_tl, $seg_2_hd, $seg_2_tl)=@_;
  
  if (&OverlayCheck($seg_1_hd, $seg_1_tl, $seg_2_hd, $seg_2_tl)){}
  else {die "\n\nIn package SeqSegmentsTools\nThe vals below in &JointSegments_sameDirection_Up_or_down_noEqual are not overlying, they should not be Joint!!\n$seg_1_hd, $seg_1_tl\n $seg_2_hd, $seg_2_tl\n\n";}

  my $dirct_1=&DirectionDetact($seg_1_hd, $seg_1_tl);
  my $dirct_2=&DirectionDetact($seg_2_hd, $seg_2_tl);
  my ( $newHead,  $newTail )
    =( $seg_1_hd, $seg_1_tl);
    
  if (  ($dirct_1 == $dirct_2) && ( ($dirct_1 == 0) || ($dirct_1 == 1) )  ){

    if ($dirct_1 == 0){
    	if ($seg_2_hd <= $newHead){$newHead=$seg_2_hd;}
      if ($seg_2_tl >= $newTail){$newTail=$seg_2_tl;}
    }
    else{
    	if ($seg_2_hd >= $newHead){$newHead=$seg_2_hd;}
      if ($seg_2_tl <= $newTail){$newTail=$seg_2_tl;}
    }
  }
  else {
    die "\n\nIn package SeqSegmentsTools\nIn Sub JointSegments_sameDirection_Up_or_down_noEqual: The segment is not in the same direction, or there is no direction segment such as 21..21 was input!!\n\$dirct_1=$dirct_1\t\$dirct_2=$dirct_2\n\$seg_1_hd=$seg_1_hd, \$seg_1_tl=$seg_1_tl, \$seg_2_hd=$seg_2_hd, \$seg_2_tl=$seg_2_tl\n\n";
  }
  
  
  return [ $newHead,  $newTail ];
  
}

sub JointSegments_sameDirection_Up_or_down_with_ZorF_inform{     #pmAble#  #如果两个数段是重叠的，则将两个数段进行叠加, 并且检测数段的方向，方向满足则运行，否则die
  my ($seg_1_hd, $seg_1_tl, $seg_1_zf, $seg_2_hd, $seg_2_tl, $seg_2_zf)=@_;
  my $warnMsgHead="\n\nIn Package SeqSegmentsTools,\t In sub JointSegments_sameDirection_Up_or_down_with_ZorF_inform,\n\n";
  
  if (
       (  ( defined ($seg_1_hd) ) && ($seg_1_hd=~m/^\d+$/)  )  &&  (  ( defined ($seg_1_tl) ) && ($seg_1_tl=~m/^\d+$/)  )  &&  (  ( defined ($seg_2_hd) ) && ($seg_2_hd=~m/^\d+$/)  )  &&  (  ( defined ($seg_2_tl) ) && ($seg_2_tl=~m/^\d+$/)  )  
       &&
       (  ( defined ($seg_1_zf) ) && ($seg_1_zf=~m/^(\+|-)$/)  )  &&  (  ( defined ($seg_2_zf) ) && ($seg_2_zf=~m/^(\+|-)$/)  )   
     ){} else {my $dieMsg= "\n\n\nDIE!!!!$warnMsgHead\nThe input is not okay!!!!\n\nInput:{(1)$seg_1_hd, (2)$seg_1_tl, (3)$seg_1_zf, (4)$seg_2_hd, (5)$seg_2_tl, (6)$seg_2_zf}\n\n\n"; print $dieMsg; die $dieMsg;}
  
  if ($seg_1_zf ne $seg_2_zf){my $dieMsg= "\n\n\nDIE!!!!$warnMsgHead\nThe + or - of those input should be the same + or the same -!!!!\n\nInput:{(3)$seg_1_zf, (6)$seg_2_zf}\n\n\n"; print $dieMsg; die $dieMsg;}
  
  if (&OverlayCheck($seg_1_hd, $seg_1_tl, $seg_2_hd, $seg_2_tl)){}
  else {my $dieMsg= "\n\n\nDIE!!!!$warnMsgHead\nThe values below in are not overlying, they should not be Joint!!\nvalues:{(1)$seg_1_hd, (2)$seg_1_tl\n (4)$seg_2_hd, (5)$seg_2_tl}\n\n";print $dieMsg; die $dieMsg;}

  my $sm_1=&getSmallOne($seg_1_hd, $seg_1_tl); my $bg_1=&getbigOne($seg_1_hd, $seg_1_tl); $seg_1_hd=$sm_1; $seg_1_tl=$bg_1; if ($seg_1_zf eq '-'){$seg_1_hd=$bg_1; $seg_1_tl=$sm_1;}
  my $sm_2=&getSmallOne($seg_2_hd, $seg_2_tl); my $bg_2=&getbigOne($seg_2_hd, $seg_2_tl); $seg_2_hd=$sm_2; $seg_2_tl=$bg_2; if ($seg_2_zf eq '-'){$seg_2_hd=$bg_2; $seg_2_tl=$sm_2;}
  
  my ( $newHead,  $newTail )
    =( $seg_1_hd, $seg_1_tl);
    
  if ($seg_1_zf eq '+'){
  	if ($seg_2_hd <= $newHead){$newHead=$seg_2_hd;}
    if ($seg_2_tl >= $newTail){$newTail=$seg_2_tl;}
  }
  else{
  	if ($seg_2_hd >= $newHead){$newHead=$seg_2_hd;}
    if ($seg_2_tl <= $newTail){$newTail=$seg_2_tl;}
  }  
  return [ $newHead,  $newTail ];
  
}


sub GetTheCorespondPosOFDIfferentSegement{  #my $mid_unKow_pos_in_2 = SeqSegmentsTools::GetTheCorespondPosOFDIfferentSegement ( $pos_st_1, $pos_ed_1, $strand_1, $mid_Known_pos_in_1, $pos_st_2, $pos_ed_2, $strand_2);
	my (
	     $pos_st_1, $pos_ed_1, $strand_1,                 #
	     $mid_Known_pos_in_1, 
	     
	     $pos_st_2, $pos_ed_2, $strand_2
	   )=@_;
	   
	my $segment_Length_1=( abs ($pos_st_1 - $pos_ed_1) )+1;
	my $segment_Length_2=( abs ($pos_st_2 - $pos_ed_2) )+1;
	my $mid_unKow_pos_in_2;
	if ( $segment_Length_1 == $segment_Length_2 ){
	  my $sm_pos_1=&getSmallOne($pos_st_1, $pos_ed_1); my $bg_pos_1=&getbigOne($pos_st_1, $pos_ed_1);
	  my $sm_pos_2=&getSmallOne($pos_st_2, $pos_ed_2); my $bg_pos_2=&getbigOne($pos_st_2, $pos_ed_2);
	  if ($strand_1 eq '+'){ ($pos_st_1, $pos_ed_1)=($sm_pos_1, $bg_pos_1); } elsif ($strand_1 eq '-'){ ($pos_st_1, $pos_ed_1)=($bg_pos_1, $sm_pos_1); } else {die "\n\n\nDIE!!!!! In pacakege SeqSegmentsTools,\nIn sub GetTheCorespondPosOFDIfferentSegement,\n\$strand_1=$strand_1 : Which should be + or -!!!\n\n\n\n";}
	  if ($strand_2 eq '+'){ ($pos_st_2, $pos_ed_2)=($sm_pos_2, $bg_pos_2); } elsif ($strand_1 eq '-'){ ($pos_st_2, $pos_ed_2)=($bg_pos_2, $sm_pos_2); } else {die "\n\n\nDIE!!!!! In pacakege SeqSegmentsTools,\nIn sub GetTheCorespondPosOFDIfferentSegement,\n\$strand_2=$strand_2 : Which should be + or -!!!\n\n\n\n";}
	  my $detal_of_md_st_1=abs($mid_Known_pos_in_1-$pos_st_1);
	  if ($strand_2 eq '+'){$mid_unKow_pos_in_2=$pos_st_2+$detal_of_md_st_1;}
	  else                 {$mid_unKow_pos_in_2=$pos_st_2-$detal_of_md_st_1;}
	}
	else{
		die"\n\n\nDIE!!!!! In pacakege SeqSegmentsTools,\nIn sub GetTheCorespondPosOFDIfferentSegement,\nThe 2 segment is not the same about their length:\n\$pos_st_1=$pos_st_1, \$pos_ed_1=$pos_ed_1, \$strand_1=$strand_1,\n\$pos_st_2=\$pos_st_2, \$pos_ed_2=$pos_ed_2, \$strand_2=$strand_2\n\$segment_Length_1=$segment_Length_1 != $segment_Length_2=\$segment_Length_2\n\n\n";
	}
	
	return $mid_unKow_pos_in_2;
	
}



sub Chang_ZoF_BackToOrgContig_ZoF{
  my $warnMstHead="\n\n\nIn package SeqSegmentsTools,		\t In sub Chang_ZoF_BackToOrgContig_ZoF,\n";  my $dieHeadMsg="\n\n\nDIE!!!!!\n\n$warnMstHead"; 
  my ($old_ZoF, $segmentCtgZoF)=@_;
  if ( ($old_ZoF eq '+') || ($old_ZoF eq '-') ){                                                                   }
	else  {warn $dieHeadMsg; die "$dieHeadMsg\nThe \$old_ZoF=$old_ZoF\nIt should be + or -, otherwise die!!!!\n\n\n";}
	if ( ($segmentCtgZoF eq '+') || ($segmentCtgZoF eq '-') ){                                                                   }
	else  {warn $dieHeadMsg; die "$dieHeadMsg\nThe \$segmentCtgZoF=$segmentCtgZoF\nIt should be + or -, otherwise die!!!!\n\n\n";}
  my $new_ZoF; if ($old_ZoF eq $segmentCtgZoF){$new_ZoF='+';} else {$new_ZoF='-';} #{ ++ or --,$newExoZoF='+';  }  { +- or -+,$newExoZoF='-';  }  
  return $new_ZoF;
}


sub Chang_pos_BackToOrgContigPositionk{  # If you got a position from a segment of long contig, using the position, the segment start and end postion , segment strand, and the length of the original long contig, this method return a new postionn the long contig!!!
	my ($inPosition, $segmentCtgStt, $segmentCtgEnd, $segmentCtgZoF)=@_;
	my $warnMstHead="\n\n\nIn package SeqSegmentsTools,		\t In sub Chang_pos_BackToOrgContigPositionk,\n";
	my $dieHeadMsg="\n\n\nDIE!!!!!\n\n$warnMstHead"; 
	#print $warnMstHead; print "1 \$inPosition=$inPosition\n"; print "2 \$segmentCtgStt=$segmentCtgStt\n"; print "3 \$segmentCtgEnd=$segmentCtgEnd\n"; print "4 \$segmentCtgZoF=$segmentCtgZoF\n"; print "\n\n\n";
	#warn  $warnMstHead; warn  "1 \$inPosition=$inPosition\n"; warn  "2 \$segmentCtgStt=$segmentCtgStt\n"; warn  "3 \$segmentCtgEnd=$segmentCtgEnd\n"; warn  "4 \$segmentCtgZoF=$segmentCtgZoF\n";  warn  "\n\n\n";
	
  if ( ($segmentCtgZoF eq '+') || ($segmentCtgZoF eq '-') ){                                                                   }
	else  {warn $dieHeadMsg; die "$dieHeadMsg\nThe \$segmentCtgZoF=$segmentCtgZoF\nIt should be + or -, otherwise die!!!!\n\n\n";}
	
	my $orgWholeLengthCtg_ZorF='+';
	my $newPos=SeqSegmentsTools::GetTheCorespondPosOFDIfferentSegement( 1, ($segmentCtgEnd-$segmentCtgStt+1), $segmentCtgZoF, $inPosition, $segmentCtgStt, $segmentCtgEnd, $orgWholeLengthCtg_ZorF);
	  		
	
	
	return $newPos;
}

1;