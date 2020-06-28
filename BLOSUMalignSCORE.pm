#!/usr/bin/perl -w
BEGIN{  push (@INC, '/home/fredjiang/fredProgs/PMs')  };

use strict;
use warnings;
                      
package BLOSUMalignSCORE;


##################################################################################################
#
#   
#
#
##################################################################################################

                       # 0   1   2   3   4   5   6   7   8   9   10  11  12  13  14  15  16  17  18  19    20  21  22  23     24     25  
                       # A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V     U   X   *   B(N/D) Z(Q/E) J(I/L) 
my $blosum62aa=     [                                                                                                                        
                       [ 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0,   -4, -4, -4, -2,    -1,     -1      ],    #A   0
     	                 [-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -2, -3, -2, -1, -1, -3, -2, -3,   -4, -4, -4, -1,     1,     -2      ],    #R   1 
     	                 [-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,   -4, -4, -4,  4,     0,     -3      ],    #N   2
     	                 [-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,   -4, -4, -4,  4,     1,     -3      ],    #D   3
     	                 [ 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -1, -1, -2, -3, -1, -1, -2, -2, -1,    0, -4, -4, -3,    -3,     -1      ],    #C   4
     	                 [-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,   -4, -4, -4,  0,     4,     -2      ],    #Q   5
     	                 [-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,   -4, -4, -4,  1,     4,     -3      ],    #E   6
     	                 [ 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3,   -4, -4, -4,  0,    -2,     -4      ],    #G   7
     	                 [-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,   -4, -4, -4,  0,     0,     -3      ],    #H   8
     	                 [-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3,   -4, -4, -4, -3,    -3,      3      ],    #I   9
     	                 [-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1,   -4, -4, -4, -3,    -2,      3      ],    #L   10
     	                 [-1,  2,  0, -1, -1,  1,  1, -2, -1, -3, -2,  5, -1,- 3, -1,  0, -1, -3, -2, -2,   -4, -4, -4,  0,     1,     -2      ],    #K   11
     	                 [-1, -2, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1,   -4, -4, -4, -2,    -1,      2      ],    #M   12
     	                 [-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1,   -4, -4, -4, -3,    -3,      0      ],    #F   13
     	                 [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2,   -4, -4, -4, -1,    -1,     -3      ],    #P   14
     	                 [ 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,   -4, -4, -4,  1,     0,     -2      ],    #S   15
     	                 [ 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0,   -4, -4, -4,  0,    -1,     -1      ],    #T   16
     	                 [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3,   -4, -4, -4, -4,    -2,     -2      ],    #W   17
     	                 [-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1,   -4, -4, -4, -2,    -1,     -1      ],    #Y   18
     	                 [ 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4,   -4, -4, -4, -3,    -2,      2      ],    #V   19
     	                 
                       [-4, -4, -4, -4,  0, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,    9,  0,  0, -4,    -4,     -4      ],    #U   20  自己加的，依据是U和所有的氨基酸的替代得分都是-4，和C是0，和自己是9
     	                 [-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,    0,  0,  0, -4,    -4,     -4      ],    #X   21  自己加的，所有和X的替代的氨基酸得分都是-4，和U及X及*的得分是0
     	                 [-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,    0,  0,  0, -4,    -4,     -4      ],    #*   22  自己加的，和X一样
     	                 [-2, -1,  4,  4, -3,  0,  1,  0,  0, -3, -3,  0, -2, -3, -1,  1,  0, -4, -2, -3,   -4, -4, -4,  6,     1,     -3      ],    #B   23  ( N or D ) 自己加的，依据N和D的替代分，取平均值，如平局值为小数，则加0.5凑整
     	                 [-1,  1,  0,  1, -3,  4,  4, -2,  0, -3, -2,  1, -1, -3, -1,  0, -1, -2, -1, -2,   -4, -4, -4,  1,     5,     -2      ],    #Z   24  ( Q or E ) 自己加的，依据Q和E的替代分，取平均值，如平局值为小数，则加0.5凑整
     	                 [-1, -2, -3, -3, -1, -2, -3, -4, -3,  3,  3, -2,  2,  0, -3, -2, -1, -2, -1,  2,   -4, -4, -4, -3,    -2,      4      ]     #J   25  ( I or J ) 自己加的，依据Q和E的替代分，取平均值，如平局值为小数，则加0.5凑整 
    	                 
     	                                                                                                                                
     	              ];  



#my $avregaeScore=BLOSUMalignSCORE::Score_for_Alignment ($qur_line, $sbj_line);
sub Score_for_Alignment{  # 输入 alignment的两个 相同长度的 可以包括Gap的 氨基酸字符串，计算其相似性得分
	my ($qur_line, $sbj_line)=@_;
	
	my $warnMsgBody="\nIn package  BLOSUMalignSCORE,\tIn sub Score_for_Alignment,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	my $gap_open_Score=-11;
	my $gap_each_Score=-1;
	
	
	if   (   (  defined ( $qur_line )  ) && ( $qur_line=~m/\S+/ )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$qur_line=$qur_line should be a defined not empty string  !!  $!\n\n\n".$caller_inform ); 	}
	
	if   (   (  defined ( $sbj_line )  ) && ( $sbj_line=~m/\S+/ )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$sbj_line=$sbj_line should be a defined not empty string  !!  $!\n\n\n".$caller_inform ); 	}
	
	my $qur_len=length ($qur_line);
	my $sbj_len=length ($sbj_line);
	
	if ($qur_len == $sbj_len){} else {		DieWork::Just_dieWork( $die_MsgHead."\n \$qur_len=$qur_len should equal \$sbj_len=$sbj_len !!  $!\n\n\n".$caller_inform ); 	}
	
	my @Qur_strArray=split ('',$qur_line);
	my @Sbj_strArray=split ('',$sbj_line);
	
	my $AlignArray;
	for ( my $i=0; $i < @Qur_strArray; $i++ ){
		$AlignArray->[$i]=BLOSUMalignSCORE::Build_small_match_hash($Qur_strArray[$i], $Sbj_strArray[$i]);
	}
	
	my $GapGroupNumber=BLOSUMalignSCORE::Count_GapGroup_number($AlignArray);  print "\n 20191004-1642-0-0-0 \$GapGroupNumber=$GapGroupNumber\n";
	
	my $final_Score=$GapGroupNumber*$gap_open_Score;                          print "\n 20191004-1642-0-0-2 \$final_Score=$final_Score\n";
	my $match_ShowL; #制作一个 显示每个位置匹配得分的string
	for ( my $j=0; $j < @{ $AlignArray }; $j++ ){
		if   (   (  defined ( $AlignArray->[$j]->{'0_3_0_match__Type'} )  ) && ( $AlignArray->[$j]->{'0_3_0_match__Type'}=~m/\S+/ )   ){}
	  else{		DieWork::Just_dieWork( $die_MsgHead."\n \$AlignArray->[$j]->{'0_3_0_match__Type'}=$AlignArray->[$j]->{'0_3_0_match__Type'} should be a defined not empty match type String  !!  $!\n\n\n".$caller_inform ); 	}
		
		my $ScoreNB_trans2_OneWordCHAR;
		my $GapMch_type=$AlignArray->[$j]->{'0_3_0_match__Type'};
		
		
		if    ( $GapMch_type eq 'Gap') {
			$final_Score+=$gap_each_Score;                                        print "\n 20191004-1642-0-0-3 \$final_Score=$final_Score\n";   
			
			$ScoreNB_trans2_OneWordCHAR=BLOSUMalignSCORE::TransScore_to_OneChar($gap_each_Score);
			
		}
		elsif ( $GapMch_type eq 'Match') {
			my $matchScore=BLOSUMalignSCORE::Caculate_Score_for_each_Match($AlignArray->[$j]->{'0_0_0_qry_aa_char'}, $AlignArray->[$j]->{'0_1_0_sbj_aa_char'});
			$AlignArray->[$j]->{'0_4_0_match_Score'}=$matchScore;
			$final_Score+=$matchScore;                                            print "\n 20191004-1642-0-0-4 \$final_Score=$final_Score\n"; 
			
			$ScoreNB_trans2_OneWordCHAR =BLOSUMalignSCORE::TransScore_to_OneChar($matchScore);
			
		}
		
		$match_ShowL.=$ScoreNB_trans2_OneWordCHAR;			
		
	}
	my $avregaeScore=$final_Score/$qur_len;
	
	#return $match_ShowL;#这个是 用来展示两列之间的match得分的 String， 如需要可输出 
	#return $AlignArray; #这个是 包含match的所有位置的 细节信息的 ARRAY，如需要 可输出
	#return $avregaeScore; #这个是 两行match的 综合得分，就是总得分除以match长度，得到的平均分 #可作为单独输出输出
	return [ $avregaeScore, $match_ShowL, $AlignArray ];
}


#my $oneWordSHOWCHAR=BLOSUMalignSCORE::TransScore_to_OneChar($matchScore);
sub TransScore_to_OneChar{ #将 blosum等矩阵的match的Score转换为单个字符的形式，便于在match line里显示
	my ($In_matchScore)=@_;
	
	my $warnMsgBody="\nIn package  BLOSUMalignSCORE,\tIn sub TransScore_to_OneChar,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	if   (   (  defined ( $In_matchScore )  ) && ( $In_matchScore=~m/(-)?\d+/ )   ){}
  else{		DieWork::Just_dieWork( $die_MsgHead."\n \$In_matchScore=$In_matchScore should be a defined number  !!  $!\n\n\n".$caller_inform ); 	}
	
	my $ChangeHash;	
	
	$ChangeHash->{-26}='z';       $ChangeHash->{26}='Z';
	$ChangeHash->{-25}='y';       $ChangeHash->{25}='Y';
	$ChangeHash->{-24}='x';       $ChangeHash->{24}='X';
	$ChangeHash->{-23}='w';       $ChangeHash->{23}='W';
	$ChangeHash->{-22}='v';       $ChangeHash->{22}='V';
	$ChangeHash->{-21}='u';       $ChangeHash->{21}='U';
	$ChangeHash->{-20}='t';       $ChangeHash->{20}='T';
	$ChangeHash->{-19}='s';       $ChangeHash->{19}='S';
	$ChangeHash->{-18}='r';       $ChangeHash->{18}='R';
	$ChangeHash->{-17}='q';       $ChangeHash->{17}='Q';
	$ChangeHash->{-16}='p';       $ChangeHash->{16}='P';
	$ChangeHash->{-15}='o';       $ChangeHash->{15}='O';
	$ChangeHash->{-14}='n';       $ChangeHash->{14}='N';
	$ChangeHash->{-13}='m';       $ChangeHash->{13}='M';
	$ChangeHash->{-12}='l';       $ChangeHash->{12}='L';
	$ChangeHash->{-11}='k';       $ChangeHash->{11}='K';
	$ChangeHash->{-10}='j';       $ChangeHash->{10}='J';
	$ChangeHash->{-9 }='i';       $ChangeHash->{9 }='9';
	$ChangeHash->{-8 }='h';       $ChangeHash->{8 }='8';
	$ChangeHash->{-7 }='g';       $ChangeHash->{7 }='7';
	$ChangeHash->{-6 }='f';       $ChangeHash->{6 }='6';
	$ChangeHash->{-5 }='e';       $ChangeHash->{5 }='5';
	$ChangeHash->{-4 }='d';       $ChangeHash->{4 }='4';
	$ChangeHash->{-3 }='c';       $ChangeHash->{3 }='3';
	$ChangeHash->{-2 }='b';       $ChangeHash->{2 }='2';
	$ChangeHash->{-1 }='a';       $ChangeHash->{1 }='1';
	               $ChangeHash->{0}='0';
	               
  if   (   (  defined ( $ChangeHash->{$In_matchScore} )  ) && ( $ChangeHash->{$In_matchScore}=~m/^\S$/ )   ){}
  else{		DieWork::Just_dieWork( $die_MsgHead."\n \$ChangeHash->{$In_matchScore}=$ChangeHash->{$In_matchScore} should be a defined One CHAR word  !!  $!\n\n\n".$caller_inform ); 	}
	              
	return $ChangeHash->{$In_matchScore};
}


#                                            1 2    3 
#数出一个alignment中有多少连续的gap段，如 MMGGMGMMGGGGGMMMMMM（G表示gap，M表示match）有8个位置是G，但连续的Gap段 只有3处。
sub Count_GapGroup_number{  # my $GapGroup_count_number=BLOSUMalignSCORE::Count_GapGroup_number($in_Alingn_Array);
	my ($in_Alingn_Array)=@_;
	
	my $warnMsgBody="\nIn package  BLOSUMalignSCORE,\tIn sub Count_GapGroup_number,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	if  (   (  defined ( $in_Alingn_Array )  ) && (  ref ( $in_Alingn_Array ) eq 'ARRAY'  )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$in_Alingn_Array=$in_Alingn_Array should be a ARRAY  !!  $!\n\n\n".$caller_inform ); 	}
	
	my $GapGroup_count_number=0; my $last_GM_type='Type_befor_start';
	for (  my $i=0; $i < @{ $in_Alingn_Array }; $i++  ){
		
		if   (   (  defined ( $in_Alingn_Array->[$i]->{'0_3_0_match__Type'} )  ) && ( $in_Alingn_Array->[$i]->{'0_3_0_match__Type'}=~m/\S+/ )   ){}
	  else{		DieWork::Just_dieWork( $die_MsgHead."\n \$in_Alingn_Array->[$i]->{'0_3_0_match__Type'}=$in_Alingn_Array->[$i]->{'0_3_0_match__Type'} should be a defined not empty match type String  !!  $!\n\n\n".$caller_inform ); 	}
		
		my $GapMch_type=$in_Alingn_Array->[$i]->{'0_3_0_match__Type'};
		if (  ( $GapMch_type eq 'Gap') && ( $last_GM_type ne 'Gap')  ) {
			$GapGroup_count_number++;
		}
		$last_GM_type=$GapMch_type;
	}
	return $GapGroup_count_number;
}


#对每一个match的两个氨基酸或gap建立一个小的hash
sub Build_small_match_hash{    # my $outHash=BLOSUMalignSCORE::Build_small_match_hash($qry_aa, $sbj_aa);
	my ($qry_aa, $sbj_aa)=@_;
	
	my $warnMsgBody="\nIn package  BLOSUMalignSCORE,\tIn sub Build_small_match_hash,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	if   (   (  defined ( $qry_aa )  ) && ( $qry_aa=~m/^\S$/ )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$qry_aa=$qry_aa should be a defined not empty amino acid  !!  $!\n\n\n".$caller_inform ); 	}
	
	if   (   (  defined ( $sbj_aa )  ) && ( $sbj_aa=~m/^\S$/ )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$sbj_aa=$sbj_aa should be a defined not empty amino acid  !!  $!\n\n\n".$caller_inform ); 	}
	
	my $outHash;
	
	$outHash->{'0_0_0_qry_aa_char'}=$qry_aa;
  $outHash->{'0_1_0_sbj_aa_char'}=$sbj_aa;
	
	if  (   ( $qry_aa eq '-') || ( $sbj_aa eq '-')   ){
		$outHash->{'0_3_0_match__Type'}='Gap';
	}
	else{
		$outHash->{'0_3_0_match__Type'}='Match';
	}
	
	return $outHash;
}


#有一个隐藏输入，就是 全局变量$blosum62aa
sub Caculate_Score_for_each_Match{  # my $outScore= BLOSUMalignSCORE::Caculate_Score_for_each_Match($qur_aa, $sbj_aa);  #输入两个氨基酸字符（qur在前 sbj在后），由BLOSUM62计算器配对得分
	my ($qur_aa, $sbj_aa)=@_;
	
	

	my $qur_aa_nb=BLOSUMalignSCORE::aAtoNb($qur_aa);                        #print "\n 20191004-1642-0-1-0 \$qur_aa_nb=$qur_aa_nb\n";  
	my $sbj_aa_nb=BLOSUMalignSCORE::aAtoNb($sbj_aa);                        #print "\n 20191004-1642-0-1-1 \$sbj_aa_nb=$sbj_aa_nb\n";  
	              ########### 
	my $outScore =$blosum62aa->[$qur_aa_nb]->[$sbj_aa_nb];                  #print "\n 20191004-1642-0-1-2 \$outScore=\$blosum62aa->[$qur_aa_nb]->[$sbj_aa_nb]=$outScore\n";  
	              ###########
	return $outScore;
	
}

#$blosum62aa->[(&aAtoNb($QrPepAa))]->[(&aAtoNb($SbPepAa))]
sub aAtoNb{  #  BLOSUMalignSCORE::aAtoNb($aaHere);  
  my ($aaHere)=@_;
  
  my $warnMsgBody="\nIn package  BLOSUMalignSCORE,\tIn sub aAtoNb,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	
       	   
  my $aa_number_Hash ;
  $aa_number_Hash->{'A'}  = 0;	      $aa_number_Hash->{'K'}  = 11;	
  $aa_number_Hash->{'R'}  = 1;	      $aa_number_Hash->{'M'}  = 12;	
  $aa_number_Hash->{'N'}  = 2;	      $aa_number_Hash->{'F'}  = 13;	
  $aa_number_Hash->{'D'}  = 3;	      $aa_number_Hash->{'P'}  = 14;	
  $aa_number_Hash->{'C'}  = 4;	      $aa_number_Hash->{'S'}  = 15;	
  $aa_number_Hash->{'Q'}  = 5;	      $aa_number_Hash->{'T'}  = 16;	
  $aa_number_Hash->{'E'}  = 6;	      $aa_number_Hash->{'W'}  = 17;	
  $aa_number_Hash->{'G'}  = 7;	      $aa_number_Hash->{'Y'}  = 18;	
  $aa_number_Hash->{'H'}  = 8;	      $aa_number_Hash->{'V'}  = 19;	
  $aa_number_Hash->{'I'}  = 9;	      $aa_number_Hash->{'U'}  = 20;	
  $aa_number_Hash->{'L'}  = 10;	      $aa_number_Hash->{'X'}  = 21;	
  $aa_number_Hash->{'*'}  = 22;	      
  $aa_number_Hash->{'B'}  = 23;	      $aa_number_Hash->{'Z'}  = 24;	
  $aa_number_Hash->{'J'}  = 25;
  
  $aa_number_Hash->{'a'}  = 0;	      $aa_number_Hash->{'k'}  = 11;	
  $aa_number_Hash->{'r'}  = 1;	      $aa_number_Hash->{'m'}  = 12;	
  $aa_number_Hash->{'n'}  = 2;	      $aa_number_Hash->{'n'}  = 13;	
  $aa_number_Hash->{'d'}  = 3;	      $aa_number_Hash->{'p'}  = 14;	
  $aa_number_Hash->{'c'}  = 4;	      $aa_number_Hash->{'s'}  = 15;	
  $aa_number_Hash->{'q'}  = 5;	      $aa_number_Hash->{'t'}  = 16;	
  $aa_number_Hash->{'e'}  = 6;	      $aa_number_Hash->{'w'}  = 17;	
  $aa_number_Hash->{'g'}  = 7;	      $aa_number_Hash->{'y'}  = 18;	
  $aa_number_Hash->{'h'}  = 8;	      $aa_number_Hash->{'v'}  = 19;	
  $aa_number_Hash->{'i'}  = 9;	      $aa_number_Hash->{'u'}  = 20;	
  $aa_number_Hash->{'l'}  = 10;	      $aa_number_Hash->{'x'}  = 21;	
  $aa_number_Hash->{'b'}  = 23;	      $aa_number_Hash->{'z'}  = 24;	
  $aa_number_Hash->{'j'}  = 25;
  
  if (   (  defined ( $aa_number_Hash->{$aaHere} )  ) && (  $aa_number_Hash->{$aaHere} =~ m/\d+/  )   ){ 	  }
  else{		DieWork::Just_dieWork( $die_MsgHead."\n \$aa_number_Hash->{$aaHere}=$aa_number_Hash->{$aaHere} \$aaHere=$aaHere should be a right Amino acid or stop *  !!  $!\n\n\n".$caller_inform ); 	}
  
  return $aa_number_Hash->{$aaHere};
  

}



1;