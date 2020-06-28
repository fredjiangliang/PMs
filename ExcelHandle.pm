#!/usr/bin/perl -w
BEGIN{  push (@INC, '/home/fredjiang/fredProgs/PMs')  };

use strict;
use warnings;
use Excel::Writer::XLSX;
use Excel::Writer::XLSX::Utility;




package ExcelHandle;




#        example   $changeCOLkeyHASH       
#        my $changeCOLkeyHASH;
#        $changeCOLkeyHASH->{'1m1m_fmlFst'}='1_fmlFst';
#        $changeCOLkeyHASH->{'2m3m_delnum'}='2_delnum';
#        $changeCOLkeyHASH->{'3m2m_iPsRed'}='3_iPsRed';
#        $changeCOLkeyHASH->{'4m1m_domPng'}='4_domPng';
#        $changeCOLkeyHASH->{'5m1m_MrBays'}='5_MrBays';


# my $outHASH=ExcelHandle::ChangeALLpathTOhyperLINK_in_specifc_colKEYS($inPutPreExcelHASH, $changeCOLkeyHASH);
sub ChangeALLpathTOhyperLINK_in_specifc_colKEYS{  # input 2 hash, 1st the pre excel hash to change its path cells, 2nd the hash holding  the col key name which will be changed into hyper link (see the example above)
	my ($inPutPreExcelHASH, $changeCOLkeyHASH)=@_;
	
	my ($warnMsgBody, $warnMsgHead, $die_MsgHead, $caller_inform)=@{ DieWork::BuildWarnDieHeadInform( 'ExcelHandle', 'ChangeALLpathTOhyperLINK_in_specifc_colKEYS' ) };

	DieWork::Check_Hash_or_DIE( $inPutPreExcelHASH, "\$inPutPreExcelHASH", $die_MsgHead, $caller_inform  );
	DieWork::Check_Hash_or_DIE( $changeCOLkeyHASH,  "\$changeCOLkeyHASH",  $die_MsgHead, $caller_inform  );
	
	my $outHASH;
	foreach my $selIDkey (    sort { $a cmp $b } (  keys ( %{ $inPutPreExcelHASH } )  )    ){
		DieWork::Check_Hash_or_DIE( $inPutPreExcelHASH->{$selIDkey}, "\$inPutPreExcelHASH->{\$selIDkey}=\$inPutPreExcelHASH->{$selIDkey}", $die_MsgHead, $caller_inform  );
	  foreach my $colKEY (    sort { $a cmp $b } (  keys ( %{ $inPutPreExcelHASH->{$selIDkey} } )  )    ){
	  	
	  	DieWork::Check_DfdNoEmptString_or_DIE( $inPutPreExcelHASH->{$selIDkey}->{$colKEY}, "\$inPutPreExcelHASH->{\$selIDkey}->{\$colKEY}=\$inPutPreExcelHASH->{$selIDkey}->{$colKEY}", $die_MsgHead, $caller_inform  );
	  	if (  DieWork::Check_DfdNoEmptString_or_NOT( $changeCOLkeyHASH->{$colKEY} )  ) {
	  		$outHASH->{$selIDkey}->{$colKEY}=ExcelHandle::HylikFileAddress_easyVersion( $inPutPreExcelHASH->{$selIDkey}->{$colKEY}, $changeCOLkeyHASH->{$colKEY}); 
	  	}
	  	else{
	  		$outHASH->{$selIDkey}->{$colKEY}=$inPutPreExcelHASH->{$selIDkey}->{$colKEY};
	  	}
	  }
	}
	
	return $outHASH;
	
}

# ExcelHandle::TarGZ_cmd ( $inTarDIR );
sub TarGZ_cmd{   # in put a dir path, tar gz it for copy
	my ($inTarDIR, $target_path_to_hold)=@_;
	
	my ($warnMsgBody, $warnMsgHead, $die_MsgHead, $caller_inform)=@{ DieWork::BuildWarnDieHeadInform( 'ExcelHandle', 'TarGZ_cmd' ) };

	DieWork::Check_FileDirExist_or_DIE( $inTarDIR, "\$inTarDIR", $die_MsgHead, $caller_inform  );
	
	my $baseName=File::Basename::basename $inTarDIR;
	
	my $target_path="$baseName.tar.gz"; if (   (  DieWork::Check_DfdNoEmptString_or_NOT( $target_path_to_hold )  ) && ( $target_path_to_hold=~m/\.tar\.gz$/)   ){ $target_path="$target_path_to_hold/$baseName.tar.gz"; }
	my $tar_cmd="tar -czf $target_path $inTarDIR";
	DieWork::Print_and_warn( $tar_cmd."\n" );
	system "$tar_cmd";
	
}

sub HylikFileAdd{  #    my $hyperLinkFileAddress=ExcelHandle::HylikFileAdd($inFileAdd, $inPathHead, $showName, $inFileRefDumpOrNot);                                         newVersion of Hyplink maker
  my ($inFileAdd, $inPathHead, $showName, $inFileRefDumpOrNot)=@_;
  
  my $warnMsgBody="\nIn package  ExcelHandle,\tIn sub HylikFileAdd,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $subCallereIfm=DirFileHandle::print_SubCallerInform;
  
  
  if (   (  defined ( $inFileRefDumpOrNot )  ) && ( $inFileRefDumpOrNot == 1 )   ){
  	$inFileAdd=$inFileAdd.".read.txt";
  	my $pl_ref_add=$inFileAdd.".read.txt.plx";
  	system ( "cp -f $inFileAdd $pl_ref_add");
  	$inFileAdd=$pl_ref_add;
  } 
  
  #if    ($inFileAdd=~m/png$/){} 
  #elsif ($inFileAdd=~m/html$/){} 
  #else { $inFileAdd=$inFileAdd.".read.txt"; }
	$inFileAdd=ExcelHandle::Lnx2Dos( $inFileAdd ); $inFileAdd=$inPathHead."\\".$inFileAdd; $inFileAdd=ExcelHandle::makeHyperLk( $inFileAdd,  $showName ); 
	return $inFileAdd;
} 

sub HylikFileAddress_easyVersion{  #    my $hyperLinkFileAddress=ExcelHandle::HylikFileAddress_easyVersion($inFileAdd, $showName, $inFileRefDumpOrNot);                                         newVersion of Hyplink maker
  my ($inFileAdd, $showName, $inFileRefDumpOrNot)=@_;
  
  my $warnMsgBody="\nIn package  ExcelHandle,\tIn sub HylikFileAddress_easyVersion,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $subCallereIfm=DirFileHandle::print_SubCallerInform;
  
  
  if (   (  defined ( $inFileRefDumpOrNot )  ) && ( $inFileRefDumpOrNot == 1 )   ){
  	$inFileAdd=$inFileAdd.".read.txt";
  	my $pl_ref_add=$inFileAdd.".read.txt.plx";
  	system ( "cp -f $inFileAdd $pl_ref_add");
  	$inFileAdd=$pl_ref_add;
  } 
  
  #if    ($inFileAdd=~m/png$/){} 
  #elsif ($inFileAdd=~m/html$/){} 
  #else { $inFileAdd=$inFileAdd.".read.txt"; }
	$inFileAdd=ExcelHandle::Lnx2Dos( $inFileAdd ); 
	$inFileAdd=ExcelHandle::makeHyperLk( $inFileAdd,  $showName ); 
	return $inFileAdd;
} 


##sub          HplkOrNot{  #��� �õ�Ԫ�� ��ֵ �ǲ��� hyperlink
sub HplkOrNot{       #��� �õ�Ԫ�� ��ֵ �ǲ��� hyperlink
  my ($inPutVal)=@_;
  my $yesOrNot=0; #==HYPERLINK("test.txt.tempdir\Acetabularia_acetabulum20140121\D1_oldOutPut20140405\b_EtbnSe_OT\Ali_Msf\1\10__prs.txt", "-1")
  if ( $inPutVal=~m/=HYPERLINK\(\"(\S+)\",\s*\"(\S+)\"\)/ ){$yesOrNot=1;}
  return $yesOrNot;
}
##sub2.2          GetLableOfHpLk{    #���Hyperlink��Ԫ��� ��ʾֵ
sub GetLableOfHpLk{    #  my $outValHr=ExcelHandle::GetLableOfHpLk($inPutHplk); ���Hyperlink��Ԫ��� ��ʾֵ
  my ($inPutHplk)=@_;                   #=HYPERLINK("test.txt.tempdir\Alexandrium_catenella20140121\D10_oldOutPut20140414\b_EtbnSe_OT\Ali_Msf\1\23__wis.txt", "459")
  my $outValHr;                         #=HYPERLINK("a_GtbnSe_OT\Ali_Msf\1\1__wis.txt","105130")
  if ($inPutHplk=~m/=HYPERLINK\(\"(\S*)\",\s*\"(\S*)\"\)/){
    $outValHr=$2;                         #print "\$1=$1\t\$2=$2\n";    
  }
  else {
    $outValHr=$inPutHplk;
  }
  return $outValHr;
}
##sub2.3          GetHpLk{          #���Hyperlink��Ԫ��� ����������
sub GetHpLk{          #���Hyperlink��Ԫ��� ����������
  my ($inPutHplk)=@_; 
  my $outValHr;                         #=HYPERLINK("a_GtbnSe_OT\Ali_Msf\1\1__wis.txt","105130")
  if ($inPutHplk=~m/=HYPERLINK\(\"(\S*)\",\s*\"(\S*)\"\)/){
    $outValHr=$1;                         #print "\$1=$1\t\$2=$2\n";    
  }
  else {
    $outValHr=$inPutHplk;
  }
  return $outValHr;
}

sub makeHyperLk{     #my $outV=ExcelHandle::makeHyperLk($InAdr, $InVal);   #�ú�������������HYPERLINK��ʽ��Ϣ������ ��ַ �� ֵ���κ��ַ����ı����� ����
  my ($InAdr, $InVal)=@_;
  if (  defined ($InAdr)  ){}else{$InAdr='Undefined';}
  if (  defined ($InVal)  ){}else{$InVal='Undefined';}
  my $outV="=HYPERLINK(\"$InAdr\", \"$InVal\")";
  return $outV;
}

##sub2.5          Lnx2Dos{     #��linux��·����ʽ��Ϊwindows��ʽ
sub Lnx2Dos{    #ExcelHandle::Lnx2Dos  #��linux��·����ʽ��Ϊwindows��ʽ
  my ($pathFile)=@_;
  $pathFile=~tr/\//\\/;
  return $pathFile;
}
##sub2.6          Dos2Lnx{     #��windows��·����ʽ��Ϊlinux��ʽ
sub Dos2Lnx{   #ExcelHandle::Dos2Lnx   #��windows��·����ʽ��Ϊlinux��ʽ
  my ($pathFile)=@_;
  $pathFile=~tr/\\/\//;
  return $pathFile;
}


#my $outPutRowColHash=ExcelHandle::FmtHashChgIn2RowColFmtHash($inPutFmtHashAddr);
##sub8.1          FmtHashChgIn2RowColFmtHash{   #�� format��Ϊkey��hash    ת��Ϊ��    row��col��key��hash, ����ײ��ֵ�� ��  val �� format ��ɵ������ ��ַ�� �� [val ,  format]��
sub FmtHashChgIn2RowColFmtHash{   #�� format��Ϊkey��hash    ת��Ϊ��    row��col��key��hash, ����ײ��ֵ�� ��  val �� format ��ɵ������ ��ַ�� �� [val ,  format]��
  my ($inPutFmtHashAddr)=@_;
  my $outPutRowColHash;
  foreach my $fmtKeyHr (keys (%{$inPutFmtHashAddr})){                       #warn "\$fmtKeyHr=$fmtKeyHr\n";
	  foreach my $AddrOf2dArrHere (@{ $inPutFmtHashAddr->{$fmtKeyHr} }){      #warn "\$AddrOf2dArrHere=$AddrOf2dArrHere\n";
	    my ($tpRowHr, $tpColHr, $tpValHr)=@{$AddrOf2dArrHere};
	    $outPutRowColHash->{$tpRowHr}->{$tpColHr}=[$tpValHr,$fmtKeyHr];         #print "\$outPutRowColHash->{\$tpRowHr}->{\$tpColHr}=\$outPutRowColHash->{$tpRowHr}->{$tpColHr}=[\$tpValHr,\$fmtKeyHr]=[$tpValHr,$fmtKeyHr]\n";
	  }
	}
	return $outPutRowColHash;
}

#ExcelHandle::PtExcelBasedOnFmtHash($inWkBook, $inSheet, $inHash, $inforMatHash);
##sub4.2          PtExcelBasedOnFmtHash{   #��������������excel�ļ��еĺ����������� Ҫ������� �ļ��Ķ���sheet�Ķ��� ���������Ϣ�Ķ�ά����ĵ�ַ���� format��hash
sub PtExcelBasedOnFmtHash{   #��������������excel�ļ��еĺ����������� Ҫ������� �ļ��Ķ���sheet�Ķ��� ���������Ϣ�Ķ�ά����ĵ�ַ���� format��hash
  my ($inWkBook, $inSheet, $inHash, $inforMatHash)=@_;         #print "\n\n\cl\&print_all_sub_array(\$inHash)\n";&print_all_sub_array($inHash); print "\cl\n\n";
  #warn "\$inWkBook=$inWkBook\t\$inSheet=$inSheet\t\$inHash=$inHash\t\$inforMatHash=$inforMatHash\n\n";
  
  foreach my $fmt_key ( keys ( %{ $inHash } ) ){  print "\$fmt_key=$fmt_key\n";
    my $inFormat_hash=$inforMatHash->{$fmt_key}; #print "&print_all_sub_array(\$inFormat_hash)";&print_all_sub_array($inFormat_hash); print "\n\n";
    my $inFormat=$inWkBook->add_format( %{$inFormat_hash} );
    for (my $i=0; $i<@{ $inHash->{$fmt_key} }; $i++){ 
      
      my @tp=@{ $inHash->{$fmt_key}->[$i] };  #print "print \@tp\t"; map {printf "%20s\t", $_;}  @tp; print "\n";
      
      if (   (  (  defined ( $tp[3] )  ) && ($tp[3]=~m/^\d+$/) ) && (  (  defined ( $tp[3] )  ) && ($tp[4]=~m/^\d+$/)  )    ){
      	if ( ($tp[0] == $tp[3])&&($tp[1]==$tp[4]) ){
      		$inSheet->write     ($tp[0], $tp[1], $tp[2], $inFormat);
      	}else{
      		$inSheet->merge_range ($tp[0], $tp[1], $tp[3], $tp[4], $tp[2], $inFormat);
      	}
      	
      } 
      if ($tp[2]=~m/\S+/){
        $inSheet->write     ($tp[0], $tp[1], $tp[2], $inFormat);
        
      }
      
    }
  }
}




###sub8.8.1          getValHplinkFmtFrom_RowColHash{                #��row col ����hash�У�ȡ����val����lable����hpyerlink����format
sub getValHplinkFmtFrom_RowColHash{                #��row col ����hash�У�ȡ����val����lable����hpyerlink����format
  my ($hash_toExtract, $ExRow, $ExCol, $V_L_H_F)=@_;
  my ($val, $fmt)=@{ $hash_toExtract->{$ExRow}->{$ExCol} };
  my $outReturn;
  if     ($V_L_H_F eq 'V'){        $outReturn= $val;                    }
  elsif ($V_L_H_F eq 'L'){        $outReturn= &GetLableOfHpLk( $val ); }
  elsif ($V_L_H_F eq 'H'){        $outReturn= &GetHpLk( $val );        }
  elsif ($V_L_H_F eq 'F'){        $outReturn= $fmt;                    } 
  else                   {        die "In sub :\n\$V_L_H_F=$V_L_H_F,\t\$hash_toExtract=$hash_toExtract,\t\$ExRow=$ExRow,\t\$ExCol=$ExCol\n\n"; }
  return $outReturn;	
}



###sub8.8.2          PushCellIntoSheeetHash{   #��Ŀ�� hash����formatΪkey���У���Ӵ� format�ؼ���key�� cell��Ϣ���������������ֵ꣬�ȣ�
sub PushCellIntoSheeetHash{   #��Ŀ�� hash����formatΪkey���У���Ӵ� format�ؼ���key�� cell��Ϣ���������������ֵ꣬�ȣ�
  my ($hash_toFill, $row_toFill, $col_toFill, $val_toFill, $fmt_toFill)=@_;
  push @{ $hash_toFill->{$fmt_toFill} } , [$row_toFill, $col_toFill, $val_toFill];
  return $hash_toFill;
}

##sub4.1          MakeFormatHash_1 {   #������� formatHash����װ�ض��� ��ʽ�ĺ����� ���е���Ҫ�ֶ�����ĸ�ʽ��Ϣ����������װ��
sub MakeFormatHash_1 {   #������� formatHash����װ�ض��� ��ʽ�ĺ����� ���е���Ҫ�ֶ�����ĸ�ʽ��Ϣ����������װ��
  my ($insub_arrayAdd, $insub_hashAdd, $newInfmtHash)=@_;
  my $tempFmtHash;  
  my @whiteFond=(12,16,17,18,19,20,21,23,24,25,28,29,30,32,33,36,37,38,39,48,53,54,55,56,57,58,59,60,61,62,63);
  my $whiteFondHash;
  my $sizeNb=8;
  foreach my $whiteNB (@whiteFond){ $whiteFondHash->{$whiteNB}=1;  	}
  my @NB10_63colorFmt; my $idx=0;
  for (my $i=12; $i<=(63); $i++){
  	$NB10_63colorFmt[$idx]={ -bold => 1, -bg_color =>  $i, -color => 0, -size => $sizeNb };
  	if ( defined ( $whiteFondHash->{$i} ) ){
  	  if($whiteFondHash->{$i} == 1){
  	  	$NB10_63colorFmt[$idx]={ -bold => 1, -bg_color =>  $i, -color => 1, -size => $sizeNb };
  	  }
  	}
  	$idx++;
  }
  foreach my $EachNameHash (@{ $insub_arrayAdd }){
    my $BgrColNb=0;
    foreach my $EachName (   sort {$a cmp $b} (   keys (  %{ $EachNameHash }  )   )    ) {
      $tempFmtHash->{$EachName}=    $NB10_63colorFmt[$BgrColNb]; print "\$tempFmtHash->{\$EachName}=\$NB10_63colorFmt[$BgrColNb]=$NB10_63colorFmt[$BgrColNb]\n";
      $BgrColNb++;
      if ($BgrColNb >= $idx) {$BgrColNb=0;}
    }
          
  }
  
  foreach my $EachNameHash  (   sort {$a cmp $b} (   keys (  %{ $insub_hashAdd }  )   )    ){
    my $BgrColNb=0;
    foreach my $EachName (   sort {$a cmp $b} (   keys (  %{ $insub_hashAdd->{$EachNameHash} }  )   )    ) {
      $tempFmtHash->{$EachName}=    $NB10_63colorFmt[$BgrColNb]; print "\$tempFmtHash->{\$EachName}=\$NB10_63colorFmt[$BgrColNb]=$NB10_63colorFmt[$BgrColNb]\n";
      $BgrColNb++;
      if ($BgrColNb >= $idx) {$BgrColNb=0;}
    }
          
  }
  
  my @CisElementNameArray=('20170218_EFsec_SelB', '20170218_PSTK', '20170218_SBP2', '20170218_SPS2', '20170218_SecS_SelA', '20170218_YbbB', '20170218_ribosomalProteinL30', '20170218_secp43', '20170218_tRNAsec');
  my $BgrColNb=0;
  foreach my $EachName (@CisElementNameArray) {
    $tempFmtHash->{$EachName}=    $NB10_63colorFmt[$BgrColNb]; print "\$tempFmtHash->{\$EachName}=\$NB10_63colorFmt[$BgrColNb]=$NB10_63colorFmt[$BgrColNb]\n";
    $BgrColNb++;
    if ($BgrColNb >= $idx) {$BgrColNb=0;}
  }
  
  
  $BgrColNb=0;
  for ( my $i=1; $i<1000; $i++ ){
    $tempFmtHash->{$i}=    $NB10_63colorFmt[$BgrColNb];
    $BgrColNb++;
    if ($BgrColNb >= $idx) {$BgrColNb=0;}
  }
  
  
  #�������ֶ���д�ĸ�ʽ
  $tempFmtHash->{'normal'}                       ={-size => $sizeNb                                          };
  $tempFmtHash->{'bold'}                         ={-bold => 1, -size => $sizeNb                              };
  $tempFmtHash->{'BlackRed'}                     ={-bold => 1, -bg_color => 8,  -color => 2, -size => $sizeNb  };
  $tempFmtHash->{'Cys'}                          ={-bold => 1, -bg_color => 52, -color => 0, -size => $sizeNb  };   
  $tempFmtHash->{'Sec'}                          ={-bold => 1, -bg_color => 11, -color => 0, -size => $sizeNb  };   
  $tempFmtHash->{'g'}                            ={-bold => 1, -bg_color => 10, -color => 1, -size => $sizeNb  };   
  $tempFmtHash->{'m'}                            ={-bold => 1, -bg_color => 34, -color => 0, -size => $sizeNb  };   
  $tempFmtHash->{'e'}                            ={-bold => 1, -bg_color => 33, -color => 1, -size => $sizeNb  }; 
  $tempFmtHash->{'+'}                            ={-bold => 1, -bg_color => 9,  -color => 0, -size => $sizeNb  };   
  $tempFmtHash->{'-'}                            ={-bold => 1, -bg_color => 8,  -color => 1, -size => $sizeNb  };  
  $tempFmtHash->{'GEN'}                          ={-bold => 1, -bg_color => 12,  -color => 1, -size => $sizeNb  };
  $tempFmtHash->{'EST'}                          ={-bold => 1, -bg_color => 10,  -color => 1, -size => $sizeNb  }; 
  $tempFmtHash->{'_SECIS_Level_0'}=              {                                           };
  $tempFmtHash->{'_SECIS_Level_1'}=              {-bold => 1, -bg_color => 34,  -color => 0, -size => $sizeNb  };   #��SECIS��TGP��λ�ù�ϵ������0ʱ��cove score����0�ֻ���Std>1ʱ��Ϊ����ɫ
  $tempFmtHash->{'_SECIS_Level_2'}=              {-bold => 1, -bg_color => 53,  -color => 1, -size => $sizeNb  };   #��SECIS��TGP��λ�ù�ϵ������1ʱ��cove score����5�ֻ���Std>2ʱ��Ϊ�Ȼ�ɫ
  $tempFmtHash->{'_SECIS_Level_3'}=              {-bold => 1, -bg_color => 23,  -color => 0, -size => $sizeNb  };   #����SECIS��                                                    Ϊ��ɫ
  
  $tempFmtHash->{'_Level_0'}=                    {                                           };
  $tempFmtHash->{'_Level_1'}=                    {-bold => 1, -bg_color => 34,  -color => 0, -size => $sizeNb  };  
  $tempFmtHash->{'_Level_2'}=                    {-bold => 1, -bg_color => 53,  -color => 1, -size => $sizeNb  };  
  $tempFmtHash->{'_Level_3'}=                    {-bold => 1, -bg_color => 23,  -color => 0, -size => $sizeNb  };  
  
  
  foreach my $eachKey (   keys (  %{ $newInfmtHash }  )   ){
    $tempFmtHash->{ $eachKey } = $newInfmtHash->{$eachKey};
  }
  
  return $tempFmtHash;
}

sub PtColorExcel{   #pmAble#    #��ӡcolore����
  my ($bookObj, $sheetNm)=@_;
  my $fmtHash;
  $fmtHash->{'normal'}=                       {                                          };
  my $colorExcelHash;
  for (my $i=0; $i<200; $i++){
    push @{ $colorExcelHash->{$i} }, [$i,0,$i]; push @{ $colorExcelHash->{'normal'} }, [$i,1,$i];
    $fmtHash->{$i}=                       { -bold => 1, -bg_color =>  $i, -color => 1 };
  }
  my $OutSheetObj=$bookObj->add_worksheet($sheetNm);
  &PtExcelBasedOnFmtHash($bookObj, $OutSheetObj, $colorExcelHash, $fmtHash );
  
}






1;