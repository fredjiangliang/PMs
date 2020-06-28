
#!/usr/bin/perl -w

BEGIN{  push (@INC, '/home/fredjiang/fredProgs/PMs')  };

use strict;
use warnings;
#use DirFileHandle;
use DieWork;

package  InFileHandle;


sub readAllfileIntoAstring{   # my $outString=InFileHandle::readAllfileIntoAstring($infile);
  my ($infile)=@_;
  
  my $warnMsgBody="\nIn package  InFileHandle,\tIn sub readAllfileIntoAstring,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
	my $caller_inform=DirFileHandle::print_SubCallerInform;
  
  my $outString;
  open (FILE, $infile) or DieWork::Just_dieWork( $die_MsgHead."\n cannot open \$infile=$infile :  $!\n\n\n".$caller_inform ); 
  local $/=undef;        #用local是为了在其他时候用到分隔符时，不会出错  
  $outString= <FILE>;  
  close FILE;  
  return $outString;
}



sub readAllfileInto_a_Array{   # my $outArray=InFileHandle::readAllfileInto_a_Array($infile, $sgm_div_mark);
  my ($infile, $sgm_div_mark)=@_;
  
  my $warnMsgBody="\nIn package  InFileHandle,\tIn sub readAllfileInto_a_Array,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
	my $caller_inform=DirFileHandle::print_SubCallerInform;
  
  if   (   (  defined ( $sgm_div_mark )  ) && ( $sgm_div_mark=~m/\S+/ )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$sgm_div_mark=$sgm_div_mark should be a defined not empty string  !!  $!\n\n\n".$caller_inform ); 	}
	
  
  my $outArray;
  my $outArray_idx=0;
  open (FILE, $infile) or DieWork::Just_dieWork( $die_MsgHead."\n cannot open \$infile=$infile :  $!\n\n\n".$caller_inform ); 
  local $/=$sgm_div_mark;        #用local是为了在其他时候用到分隔符时，不会出错  
  #warn "\nqqqqqqq $/\n";
  while (<FILE>) {
  	my $tempString=$_;
  	$tempString=~s/^$sgm_div_mark//;
  	$tempString=~s/$sgm_div_mark$//;  #warn "\nbbbbb $tempString\n";
  	if( $tempString=~m/\S+/ ){
  		
  		$outArray->[$outArray_idx]=$tempString ;
  		$outArray_idx++;  
  	}
  }
  close FILE;  
  return $outArray;
}


sub PrintStringIntoFile{  #  InFileHandle::PrintStringIntoFile($outFile, $outString);
  my ($outFile, $outString, $ptOutWrn)=@_;
  
  my $warnMsgBody="\nIn package  outFileHandle,\tIn sub PrintStringIntoFile,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
	my $caller_inform=DirFileHandle::print_SubCallerInform;
  
  print "20190329-1718-0-0 $ptOutWrn\n" if (   (  defined ( $ptOutWrn )  ) && ( $ptOutWrn=~m/\S+/ )   );
  open (FILE, ">$outFile") or DieWork::Just_dieWork( $die_MsgHead."\n cannot create \$outFile=$outFile :  $!\n\n\n".$caller_inform ); 
  print  FILE $outString;
  close FILE;  
  return 1;
}

sub GetGenomeSequencingStatues{
  my $eukaryotesTableFile="";

}


#  my $FlnTableString = InFileHandle::Add_key_head_forNOheadTable($inFile, $outFile, $inHeadArray, $inHeadExplainArray);
sub Add_key_head_forNOheadTable{ #为没有表头 关键字属性行的 表格文本文件，添加表头（表头用Array ref的形式输入）；
	my ($inFile, $outFile, $inHeadArray, $inHeadExplainArray)=@_;
	
	my $warnMsgBody="\nIn package  InFileHandle,\tIn sub Add_key_head_forNOheadTable,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";  
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	if   (   (  defined ( $inFile )  ) && ( $inFile=~m/\S+/ ) && (  ( -e $inFile )  )   ){}
  else{		DieWork::Just_dieWork( $die_MsgHead."\n \$inFile=$inFile should be a defined not empty file  !!  $!\n\n\n".$caller_inform ); 	}
	
	if  (   (  defined ( $inHeadArray )  ) && (  ref ( $inHeadArray ) eq 'ARRAY'  )   ){}
  else{		DieWork::Just_dieWork( $die_MsgHead."\n \$inHeadArray=$inHeadArray should be a ARRAY  !!  $!\n\n\n".$caller_inform ); 	}
	  
	
	my $tableArray=InFileHandle::ReadTableFile_withOut_headLine( $inFile );  #DirFileHandle::PrintAndWarnDumper ($tableArray);
	
	my $colArraySize=@{ $tableArray->[0] };
	my $inHeadArySze=@{ $inHeadArray };
	
	if   (   (  defined ( $colArraySize )  ) && ( $colArraySize=~m/\d+/ ) && ( $colArraySize > 0 )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$colArraySize=$colArraySize should be a number > 0  !!  $!\n\n\n".$caller_inform ); 	}
	if   (   (  defined ( $inHeadArySze )  ) && ( $inHeadArySze=~m/\d+/ ) && ( $inHeadArySze >0 )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$inHeadArySze=$inHeadArySze should be a number > 0  !!  $!\n\n\n".$caller_inform ); 	}	
	if ( $colArraySize != $inHeadArySze ){ DieWork::Just_dieWork( $die_MsgHead."\n \$colArraySize=$colArraySize != $inHeadArySze=\$inHeadArySze, These 2 should be equal  !!  $!\n\n\n".$caller_inform ); }
	
	
	my $FlnTableString=InFileHandle::readAllfileIntoAstring($inFile);
	
	
	
	if  (   (  defined ( $inHeadExplainArray )  ) && (  ref ( $inHeadExplainArray ) eq 'ARRAY'  )   ){
	  my $inHeadExplainArraySze=@{ $inHeadExplainArray };
	  if   (   (  defined ( $inHeadExplainArraySze )  ) && ( $inHeadExplainArraySze=~m/\d+/ ) && ( $inHeadExplainArraySze >0 )   ){}
	  else{		DieWork::Just_dieWork( $die_MsgHead."\n \$inHeadExplainArraySze=$inHeadExplainArraySze should be a number > 0  !!  $!\n\n\n".$caller_inform ); 	}	
	  if ( $colArraySize != $inHeadExplainArraySze ){ DieWork::Just_dieWork( $die_MsgHead."\n \$colArraySize=$colArraySize != $inHeadExplainArraySze=\$inHeadExplainArraySze, These 2 should be equal  !!  $!\n\n\n".$caller_inform ); }
	  my $ExplainString=ArrayHashChange::ChangArrayIntoSimpleString_divToTable ($inHeadExplainArray);
	  $FlnTableString="#".$ExplainString."\n".$FlnTableString;
	  
	}
  else{		DieWork::Just_dieWork( $die_MsgHead."\n \$inHeadExplainArray=$inHeadExplainArray should be a defined not empty file  !!  $!\n\n\n".$caller_inform ); 	}
  
  
  my $HeadArayString=ArrayHashChange::ChangArrayIntoSimpleString_divToTable ($inHeadArray);
	$FlnTableString="#Key#".$HeadArayString."\n".$FlnTableString;
  
  InFileHandle::PrintStringIntoFile($outFile, $FlnTableString);
	
	return $FlnTableString;
}


#  my $outArray = InFileHandle::ReadTableFile_withOut_headLine( $inFile );
sub ReadTableFile_withOut_headLine{  #读取 如blast tabluar格式的 文本文件（每行有多个由制表符分开的单元格，每行单元格个数相同，每行是一个记录）
	
	my ($inFile)=@_;
	
	my $warnMsgBody="\nIn package  InFileHandle,\tIn sub ReadTableFile_withOut_headLine,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";  
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	if   (   (  defined ( $inFile )  ) && ( $inFile=~m/\S+/ ) && (  ( -e $inFile )  )   ){}
  else{		DieWork::Just_dieWork( $die_MsgHead."\n \$inFile=$inFile should be a defined not empty file  !!  $!\n\n\n".$caller_inform ); 	}
	
	open (IN, $inFile) or DieWork::Just_dieWork( $die_MsgHead."\n cannot open \$inFile=$inFile: $!\n\n\n".$caller_inform );
	
	my $outArray;
	
  my $tempMakr=$/;
  $/="\n";
  my $rowNB=0;
  my $lastColNb=0;
  while (<IN>){ #print "\$_=$_\n"; sleep (1);
  	my @OldTempAr=split '\t', $_;
  	
  	my $ColNb=0;
  	my $thisColNB=0;
  	foreach my $cuttedCell  (@OldTempAr){
  		
  		chomp $cuttedCell;  		$cuttedCell=~s/^\"//; $cuttedCell=~s/\"$//; #print "\$cuttedCell=$cuttedCell\n";
  		
  		$outArray->[$rowNB]->[$ColNb]=$cuttedCell;
  		$ColNb++;
  	}
  	if (  ( $rowNB >0 ) && ( $thisColNB != $lastColNb )  ){
  		DieWork::Just_dieWork( $die_MsgHead."\n on the \$rowNB=$rowNB, \$thisColNB=$thisColNB != $lastColNb=\$lastColNb, these 2 should be equal: $!\n\n\n".$caller_inform )
  	}
  	$lastColNb=$thisColNB;
  	
  	
    $rowNB++;
  }
  close (IN);
  $/=$tempMakr;
  
  if  (   (  defined ( $outArray )  ) && (  ref ( $outArray ) eq 'ARRAY'  )   ){}
  else{		DieWork::Just_dieWork( $die_MsgHead."\n \$outArray=$outArray should be a ARRAY  !!  $!\n\n\n".$caller_inform ); 	}
  
  
  return $outArray;
	
}

sub ChangeIntoHashForCsv{   #输入 特定格式的 NCBI的基因组Assembly的table文件，解析变成2维hash，1D的key是assembly的id，2D的key是各个关键字字段名
	#这个文件时从 https://www.ncbi.nlm.nih.gov/genome/browse#!/eukaryotes/ 下载下来 的
  my ($inFile)=@_;
  open (IN,$inFile)or die "\n\nIn package InFileHandle, in sub ChangeIntoHash\ncannot open \$inFile=$inFile : $!\n\n";
  my $outHash;
  my $keyhash;
  my $rowKeyword="#Organism Name";    #这个是 可以作为 唯一不重复的key的 关键字
  my $rowKeywoNB=0;                   #这个是 可以作为 唯一不重复的key的 关键字,所在的列
  my $tempMakr=$/;
  $/="\n";
  my $rowNB=0;
  while (<IN>){ #print "\$_=$_\n"; sleep (1);
  	my @OldTempAr=split ',', $_;
  	my @tempAr;
  	foreach my $cutDbquotationCell  (@OldTempAr){
  		my $cuttedCell=$cutDbquotationCell; 
  		chomp $cuttedCell;
  		$cuttedCell=~s/^\"//; $cuttedCell=~s/\"$//; #print "\$cuttedCell=$cuttedCell\n";
  		push @tempAr, $cuttedCell;
  	}
  	
  	my $subOrgNB=0;
  	if (   defined (  $outHash->{ $tempAr[$rowKeywoNB] }  )   ){
      $subOrgNB=keys %{  $outHash->{ $tempAr[$rowKeywoNB] }  };	
    }
  	
  	
  	my $ColNb=0;
  	foreach my $eachCell (@tempAr){
  	  if ($rowNB==0){
        $keyhash->{$ColNb}=$eachCell;
        if ($ColNb == $rowKeywoNB){
        	if ($rowKeyword eq $eachCell){
        		#$rowKeyword 和 $rowKeywoNB都没问题#
        	}
        	else {
        	  die "\n\nIn package InFileHandle, in sub ChangeIntoHash\nThe \$rowKeyword=$rowKeyword and \$rowKeywoNB=$rowKeywoNB is not OK, please check the \$inFile=$inFile\n\n\n";
        	}
        }
      }
      else{
      	$outHash->{ $tempAr[$rowKeywoNB] }->{$subOrgNB}->{ $keyhash->{$ColNb} }=$eachCell;  
      	#print "\$outHash->{ \$tempAr[\$rowKeywoNB] }->{ \$subOrgNB }->{ \$keyhash->{\$ColNb} }=\$outHash->{ $tempAr[$rowKeywoNB] }->{ $subOrgNB }->{ $keyhash->{$ColNb} }=\$eachCell=$eachCell\n";
      }
  	  $ColNb++;
  	}    
    $rowNB++;
  }
  close (IN);
  $/=$tempMakr;
  return $outHash;
}


sub ChangeIntoHashForTxt{    #  my $outHash  = InFileHandle::ChangeIntoHashForTxt ($inFile);
	#输入 特定格式的 NCBI的基因组Assembly的table文件，解析变成2维hash，1D的key是assembly的id，2D的key是各个关键字字段名
	#这个文件时从 ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/ 下载下来的
  my ($inFile)=@_;
  open (IN,$inFile)or die "\n\nIn package InFileHandle, in sub ChangeIntoHash\ncannot open \$inFile=$inFile : $!\n\n";
  my $outHash;
  my $keyhash;
  my $rowKeyword="#Organism\/Name";    #这个是 可以作为 唯一不重复的key的 关键字
  my $rowKeywoNB=0;                   #这个是 可以作为 唯一不重复的key的 关键字,所在的列
  my $tempMakr=$/;
  $/="\n";
  my $rowNB=0;
  while (<IN>){ #print "\$_=$_\n"; sleep (1);
  	my @OldTempAr=split '\t', $_;
  	my @tempAr;
  	foreach my $cutDbquotationCell  (@OldTempAr){
  		my $cuttedCell=$cutDbquotationCell; 
  		chomp $cuttedCell;
  		$cuttedCell=~s/^\"//; $cuttedCell=~s/\"$//; #print "\$cuttedCell=$cuttedCell\n";
  		push @tempAr, $cuttedCell;
  	}
  	
  	my $subOrgNB=0;
  	if (   defined (  $outHash->{ $tempAr[$rowKeywoNB] }  )   ){
      $subOrgNB=keys %{  $outHash->{ $tempAr[$rowKeywoNB] }  };	
    }
  	
  	
  	my $ColNb=0;
  	foreach my $eachCell (@tempAr){
  	  if ($rowNB==0){
        $keyhash->{$ColNb}=$eachCell;
        if ($ColNb == $rowKeywoNB){
        	if ($rowKeyword eq $eachCell){
        		#$rowKeyword 和 $rowKeywoNB都没问题#
        	}
        	else {
        	  die "\n\nIn package InFileHandle, in sub ChangeIntoHash\nThe \$rowKeyword=$rowKeyword and \$rowKeywoNB=$rowKeywoNB is not OK, please check the \$inFile=$inFile\n\n\n";
        	}
        }
      }
      else{
      	$outHash->{ $tempAr[$rowKeywoNB] }->{$subOrgNB}->{ $keyhash->{$ColNb} }=$eachCell;  
      	#print "\$outHash->{ \$tempAr[\$rowKeywoNB] }->{ \$subOrgNB }->{ \$keyhash->{\$ColNb} }=\$outHash->{ $tempAr[$rowKeywoNB] }->{ $subOrgNB }->{ $keyhash->{$ColNb} }=\$eachCell=$eachCell\n";
      }
  	  $ColNb++;
  	}    
    $rowNB++;
  }
  close (IN);
  $/=$tempMakr;
  return $outHash;
}






#  my $outHash=InFileHandle::FormTableToHash  ($inFile, $keyColNb);
sub FormTableToHash{   #输入 特定格式的 table文件，解析变成2维hash，1D的key是指定输入列的值，2D的key是各个关键字字段名
	#/home/fredjiang/work/Algae/20150522NewDATA/2018.05.03.InPutDataForAlgaeAnalysis/ProteinFmlNameIndex.txt
	
	my ($warnMsgBody, $warnMsgHead, $die_MsgHead, $caller_inform)=@{ DieWork::BuildWarnDieHeadInform( 'InFileHandle', 'FormTableToHash' ) };


  my ($inFile, $keyColNb)=@_;
  open (IN,$inFile) or DieWork::Just_dieWork ( $die_MsgHead."\n  cannot open \$inFile=$inFile :  $!\n\n\n".$caller_inform );
  my $outHash;
  my $keyhash;
  
  #my $rowKeyword="#Organism\/Name";    #这个是 可以作为 唯一不重复的key的 关键字
 
  my $HowManyValueLineWereFound=0;    #是写了值的行  则会改写为大0的数字 记录了有多少个这种行
  my $HowManyKeyWorkLineFound=0;      #是关键字行    则会改写为大0的数字 记录了有多少个这种行
  my $tempMakr=$/;
  $/="\n";
  my $rowNB=0;
  while (<IN>){ #print "\$_=$_\n"; sleep (1);
  	my $lineHere=$_;
  	if (m/\S+/){
  		
  	  my @OldTempAr=split '\t', $_;
  	  my @tempAr;
  	  foreach my $cuttedCell  (@OldTempAr){
  	  	
  	  	chomp $cuttedCell;
  	  	$cuttedCell=~s/^\s+//; $cuttedCell=~s/\s+$//; #print "\$cuttedCell=$cuttedCell\n";
  	  	push @tempAr, $cuttedCell;
  	  }  	  	
  	  
  	  
  	  my $isThisLineKeyWordLine=0;  #是关键字行    则会改写为1
  	  my $isThisLineISnotCount=0;   #是解释性语句  则会改写为1
  	  my $isThisLineValueLine=0;    #是写了值的行  则会改写为大0的数字
  	  
  	  
  	  my $ColNb=0;
  	  foreach my $eachCell (@tempAr){
  	  	
  	  	#检查 第一个cell，看这一行时 关键字行 还是解释性语句。关键字行需要取 关键字，解释性语句则 标明这一行不计数
  	  	if ($ColNb==0){
  	  	  if ($eachCell=~/^#/){   #Regular NO1 ! #开头是#标明可能是 关键字行 或 解释性语句
  	  	  	#print "\n\nIn package InFileHandle, in sub FormTableToHash\n\Regular NO1 !\t\t$eachCell=$eachCell\n\$lineHere=$lineHere\n\n";
  	  	  	#因为关键字行 或 解释性语句 应该放在最上面，所以如果输入的文件 不符合这个格式的要求，也就是说 语句运行到这里的时候，$HowManyValueLineWereFound>0 则出错
  	  	  	if ($HowManyValueLineWereFound>0){
  	  	  	  DieWork::Just_dieWork ( $die_MsgHead."\n\$inFile=$inFile, The key line and explanning line should be on the top, please check the file: \n$inFile\n\n\n".$caller_inform );  	  	  	  
  	  	  	}
  	  	  	if ($eachCell=~/^#Key#/){    #Regular NO2 !
  	  	  	  #print "\n\nIn package InFileHandle, in sub FormTableToHash\nRegular NO2 !\t\t\$eachCell=$eachCell\n\$lineHere=$lineHere\n\n";
  	  	  	  $isThisLineKeyWordLine=1;         #说明此行 是关键字行
  	  	  	  $HowManyKeyWorkLineFound++;
  	  	  	  $eachCell=~s/^#Key#\s*//;
  	  	    }
  	  	    else {
  	  	      $isThisLineISnotCount=1;          #说明此行 是解释性语句
  	  	    }
  	  	  }
  	  	  else {
  	  	  	#因为关键字行 或 解释性语句 应该放在最上面，所以如果输入的文件 不符合这个格式的要求，也就是说 语句运行到这里的时候，$HowManyKeyWorkLineFound==0 则出错
  	  	  	if ($HowManyKeyWorkLineFound==0){
  	  	  	  DieWork::Just_dieWork ( $die_MsgHead."\n\$inFile=$inFile, The key line and explanning line should be on the top, please check the file: \n$inFile\n\n\n".$caller_inform );  	 	  	  	
  	  	  	}
  	  	  	$isThisLineValueLine=1;             #说明此行 是写了值的行
  	  	  	$HowManyValueLineWereFound++;       #记录至此 有多少个 写了值的行
  	  	  	
  	  	  }
  	  	}
  	  	
  	  	#接下来，对所有的单元格值进行操作
  	  	if ($isThisLineKeyWordLine==1){         #如果是关键字行, 则将该列放入关键字hash中，这一行的操作，在数值行出现之前完成
  	  		$keyhash->{$ColNb}=$eachCell;  print "\$keyhash->{$ColNb}=\$eachCell=$keyhash->{$ColNb}\n";
  	  	}
  	  	elsif ($isThisLineValueLine==1){        #如果是数值行，  则将该数值放入 输出$outHash中，1级键是 $keyColNb决定的第几列的 值， 2级键是 该列的关键字 ， 值是该行该列的值
  	  		if (   defined (  $outHash->{ $tempAr[$keyColNb] }->{ $keyhash->{$ColNb} }  )   ){ #如果该值已经被定义过，说明 $tempAr[$keyColNb] 或者 $keyhash->{$ColNb} 发生了重复
  	  		  
  	  		  my $dieMsg_NOW="$inFile=$inFile\n\n";
  	  		  $dieMsg_NOW .= "KEY1: \$tempAr[\$keyColNb]=\$tempAr[$keyColNb]=$tempAr[$keyColNb]\n";
  	  		  $dieMsg_NOW .= "KEY2: \$keyhash->{\$ColNb}=\$keyhash->{$ColNb}=$keyhash->{$ColNb}\n\n";
  	  		  $dieMsg_NOW .= "\$outHash->{ $tempAr[$keyColNb] }->{ $keyhash->{$ColNb} }=$outHash->{ $tempAr[$keyColNb] }->{ $keyhash->{$ColNb} } has defined before!!\n";
  	  		  $dieMsg_NOW .= "Now we want it to be \$eachCell=$eachCell !!!\n";
  	  		  $dieMsg_NOW .= "Please check the file: $inFile and KEY1 or KEY2!!!\n";
  	  		  
  	  		  DieWork::Just_dieWork ( $die_MsgHead.$dieMsg_NOW.$caller_inform );
  	  		  
  	  		  
  	  	  }
  	  		$outHash->{ $tempAr[$keyColNb] }->{ $keyhash->{$ColNb} }=$eachCell;  	  		 
  	  	}
  	    $ColNb++;
  	    
  	  }
  	}    
    $rowNB++;
  }
  close (IN);
  $/=$tempMakr;
  return $outHash;
}


#  my $outHash=InFileHandle::FormTableToHash_to_2D_hash($inFile, $keyColNb_1, $keyColNb_2);
sub FormTableToHash_to_2D_hash{   #输入 特定格式的 table文件，解析变成3维hash，1D的key是第1个指定输入列的值，2D的key是第2个指定输入列的值，3D的key是各个关键字字段名
	
  my ($inFile, $keyColNb_1, $keyColNb_2)=@_;
  
  my ($warnMsgBody, $warnMsgHead, $die_MsgHead, $caller_inform)=@{ DieWork::BuildWarnDieHeadInform( 'InFileHandle', 'FormTableToHash_to_2D_hash' ) };

  
  	
	 
  open (IN,$inFile) or DieWork::Just_dieWork ( $die_MsgHead."\n  cannot open \$inFile=$inFile :  $!\n\n\n".$caller_inform );
  my $outHash;
  my $keyhash;
  
  #my $rowKeyword="#Organism\/Name";    #这个是 可以作为 唯一不重复的key的 关键字
 
  my $HowManyValueLineWereFound=0;    #是写了值的行  则会改写为大0的数字 记录了有多少个这种行
  my $HowManyKeyWorkLineFound=0;      #是关键字行    则会改写为大0的数字 记录了有多少个这种行
  my $tempMakr=$/;
  $/="\n";
  my $rowNB=0;
  while (<IN>){ #print "\$_=$_\n"; sleep (1);
  	my $lineHere=$_;
  	if (m/\S+/){
  		
  	  my @OldTempAr=split '\t', $_;
  	  my @tempAr;
  	  foreach my $cuttedCell  (@OldTempAr){
  	  	
  	  	chomp $cuttedCell;
  	  	$cuttedCell=~s/^\s+//; $cuttedCell=~s/\s+$//; #print "\$cuttedCell=$cuttedCell\n";
  	  	push @tempAr, $cuttedCell;
  	  }  	  	
  	  
  	  
  	  my $isThisLineKeyWordLine=0;  #是关键字行    则会改写为1
  	  my $isThisLineISnotCount=0;   #是解释性语句  则会改写为1
  	  my $isThisLineValueLine=0;    #是写了值的行  则会改写为大0的数字
  	  
  	  
  	  my $ColNb=0;
  	  foreach my $eachCell (@tempAr){
  	  	
  	  	#检查 第一个cell，看这一行时 关键字行 还是解释性语句。关键字行需要取 关键字，解释性语句则 标明这一行不计数
  	  	if ($ColNb==0){
  	  	  if ($eachCell=~/^#/){   #Regular NO1 ! #开头是#标明可能是 关键字行 或 解释性语句
  	  	  	#print "\n\nIn package InFileHandle, in sub FormTableToHash\n\Regular NO1 !\t\t$eachCell=$eachCell\n\$lineHere=$lineHere\n\n";
  	  	  	#因为关键字行 或 解释性语句 应该放在最上面，所以如果输入的文件 不符合这个格式的要求，也就是说 语句运行到这里的时候，$HowManyValueLineWereFound>0 则出错
  	  	  	if ($HowManyValueLineWereFound>0){
  	  	  	  die "\n\nIn package InFileHandle, in sub FormTableToHash\n\$inFile=$inFile, The key line and explanning line should be on the top, please check the file: \n$inFile\n\n\n";  	  	  	
  	  	  	}
  	  	  	if ($eachCell=~/^#Key#/){    #Regular NO2 !
  	  	  	  #print "\n\nIn package InFileHandle, in sub FormTableToHash\nRegular NO2 !\t\t\$eachCell=$eachCell\n\$lineHere=$lineHere\n\n";
  	  	  	  $isThisLineKeyWordLine=1;         #说明此行 是关键字行
  	  	  	  $HowManyKeyWorkLineFound++;
  	  	  	  $eachCell=~s/^#Key#\s*//;
  	  	    }
  	  	    else {
  	  	      $isThisLineISnotCount=1;          #说明此行 是解释性语句
  	  	    }
  	  	  }
  	  	  else {
  	  	  	#因为关键字行 或 解释性语句 应该放在最上面，所以如果输入的文件 不符合这个格式的要求，也就是说 语句运行到这里的时候，$HowManyKeyWorkLineFound==0 则出错
  	  	  	if ($HowManyKeyWorkLineFound==0){
  	  	  	  
  	  	  	  DieWork::Just_dieWork ( $die_MsgHead."\n\$inFile=$inFile, The key line and explanning line should be on the top, please check the file: \n$inFile\n\n".$caller_inform );	  	  	
  	  	  	}
  	  	  	$isThisLineValueLine=1;             #说明此行 是写了值的行
  	  	  	$HowManyValueLineWereFound++;       #记录至此 有多少个 写了值的行
  	  	  	
  	  	  }
  	  	}
  	  	
  	  	#接下来，对所有的单元格值进行操作
  	  	if ($isThisLineKeyWordLine==1){         #如果是关键字行, 则将该列放入关键字hash中，这一行的操作，在数值行出现之前完成
  	  		$keyhash->{$ColNb}=$eachCell;  #print "\$keyhash->{$ColNb}=\$eachCell=$keyhash->{$ColNb}\n";
  	  	}
  	  	elsif ($isThisLineValueLine==1){        #如果是数值行，  则将该数值放入 输出$outHash中，1级键是 $keyColNb_1决定的第几列的 值，2级键是 $keyColNb_2决定的第几列的 值， 3级键是 该列的关键字 ， 值是该行该列的值
  	  		if (   defined (  $outHash->{ $tempAr[$keyColNb_1] }->{ $tempAr[$keyColNb_2] }->{ $keyhash->{$ColNb} }  )   ){ #如果该值已经被定义过，说明 $tempAr[$keyColNb_1] 或者 $tempAr[$keyColNb_2] 或者 $keyhash->{$ColNb} 发生了重复
  	  		  my $dieMsg_NOW="$inFile=$inFile\n\n";
  	  		  $dieMsg_NOW .= "KEY1: \$tempAr[\$keyColNb_1]=\$tempAr[$keyColNb_1]=$tempAr[$keyColNb_1]\n";
  	  		  $dieMsg_NOW .= "KEY2: \$tempAr[\$keyColNb_2]=\$tempAr[$keyColNb_2]=$tempAr[$keyColNb_2]\n";
  	  		  $dieMsg_NOW .= "KEY3: \$keyhash->{\$ColNb}=\$keyhash->{$ColNb}=$keyhash->{$ColNb}\n\n";
  	  		  $dieMsg_NOW .= "\$outHash->{ $tempAr[$keyColNb_1] }->{ $tempAr[$keyColNb_2] }->{ $keyhash->{$ColNb} }=$outHash->{ $tempAr[$keyColNb_1] }->{ $tempAr[$keyColNb_2] }->{ $keyhash->{$ColNb} } has defined before!!\n";
  	  		  $dieMsg_NOW .= "Now we want it to be \$eachCell=$eachCell !!!\n";
  	  		  $dieMsg_NOW .= "Please check the file: $inFile and KEY1 KEY2 or KEY3!!!\n";
  	  		  
  	  		  DieWork::Just_dieWork ( $die_MsgHead.$dieMsg_NOW.$caller_inform );
  	  		  
  	  	  }
  	  		$outHash->{ $tempAr[$keyColNb_1] }->{ $tempAr[$keyColNb_2] }->{ $keyhash->{$ColNb} }=$eachCell;  	  		 
  	  	}
  	    $ColNb++;
  	    
  	  }
  	}    
    $rowNB++;
  }
  close (IN);
  $/=$tempMakr;
  return $outHash;
}


#  my $outHash=InFileHandle::FormTableToHash_to_2D_hash_Array($inFile, $keyColNb_1, $keyColNb_2);
sub FormTableToHash_to_2D_hash_Array{   #输入 特定格式的 table文件，解析变成3维hash，1D的key是第1个指定输入列的值，2D的key是第2个指定输入列的值，3D是这两个key下的数组
	
  my ($inFile, $keyColNb_1, $keyColNb_2)=@_;
  
  my $warnMsgBody="\nIn package  InFileHandle,\tIn sub FormTableToHash_to_2D_hash,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";  
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	 
  open (IN,$inFile) or DieWork::Just_dieWork ( $die_MsgHead."\n  cannot open \$inFile=$inFile :  $!\n\n\n".$caller_inform );
  my $outHash;
  my $keyhash;
  
  #my $rowKeyword="#Organism\/Name";    #这个是 可以作为 唯一不重复的key的 关键字
 
  my $HowManyValueLineWereFound=0;    #是写了值的行  则会改写为大0的数字 记录了有多少个这种行
  my $HowManyKeyWorkLineFound=0;      #是关键字行    则会改写为大0的数字 记录了有多少个这种行
  my $tempMakr=$/;
  $/="\n";
  my $rowNB=0;
  while (<IN>){ #print "\$_=$_\n"; sleep (1);
  	my $lineHere=$_;
  	if (m/\S+/){
  		
  	  my @OldTempAr=split '\t', $_;
  	  my @tempAr;
  	  foreach my $cuttedCell  (@OldTempAr){
  	  	
  	  	chomp $cuttedCell;
  	  	$cuttedCell=~s/^\s+//; $cuttedCell=~s/\s+$//; #print "\$cuttedCell=$cuttedCell\n";
  	  	push @tempAr, $cuttedCell;
  	  }  	  	
  	  
  	  
  	  my $isThisLineKeyWordLine=0;  #是关键字行    则会改写为1
  	  my $isThisLineISnotCount=0;   #是解释性语句  则会改写为1
  	  my $isThisLineValueLine=0;    #是写了值的行  则会改写为大0的数字
  	  
  	  my $tempHash;
  	  
  	  my $ColNb=0;
  	  foreach my $eachCell (@tempAr){
  	  	
  	  	#检查 第一个cell，看这一行时 关键字行 还是解释性语句。关键字行需要取 关键字，解释性语句则 标明这一行不计数
  	  	if ($ColNb==0){
  	  	  if ($eachCell=~/^#/){   #Regular NO1 ! #开头是#标明可能是 关键字行 或 解释性语句
  	  	  	#print "\n\nIn package InFileHandle, in sub FormTableToHash\n\Regular NO1 !\t\t$eachCell=$eachCell\n\$lineHere=$lineHere\n\n";
  	  	  	#因为关键字行 或 解释性语句 应该放在最上面，所以如果输入的文件 不符合这个格式的要求，也就是说 语句运行到这里的时候，$HowManyValueLineWereFound>0 则出错
  	  	  	if ($HowManyValueLineWereFound>0){
  	  	  	  die "\n\nIn package InFileHandle, in sub FormTableToHash\n\$inFile=$inFile, The key line and explanning line should be on the top, please check the file: \n$inFile\n\n\n";  	  	  	
  	  	  	}
  	  	  	if ($eachCell=~/^#Key#/){    #Regular NO2 !
  	  	  	  #print "\n\nIn package InFileHandle, in sub FormTableToHash\nRegular NO2 !\t\t\$eachCell=$eachCell\n\$lineHere=$lineHere\n\n";
  	  	  	  $isThisLineKeyWordLine=1;         #说明此行 是关键字行
  	  	  	  $HowManyKeyWorkLineFound++;
  	  	  	  $eachCell=~s/^#Key#\s*//;
  	  	    }
  	  	    else {
  	  	      $isThisLineISnotCount=1;          #说明此行 是解释性语句
  	  	    }
  	  	  }
  	  	  else {
  	  	  	#因为关键字行 或 解释性语句 应该放在最上面，所以如果输入的文件 不符合这个格式的要求，也就是说 语句运行到这里的时候，$HowManyKeyWorkLineFound==0 则出错
  	  	  	if ($HowManyKeyWorkLineFound==0){
  	  	  	  die "\n\nIn package InFileHandle, in sub FormTableToHash\n\$inFile=$inFile, The key line and explanning line should be on the top, please check the file: \n$inFile\n\n\n";  	  	  	
  	  	  	}
  	  	  	$isThisLineValueLine=1;             #说明此行 是写了值的行
  	  	  	$HowManyValueLineWereFound++;       #记录至此 有多少个 写了值的行
  	  	  	
  	  	  }
  	  	}
  	  	
  	  	#接下来，对所有的单元格值进行操作
  	  	if ($isThisLineKeyWordLine==1){         #如果是关键字行, 则将该列放入关键字hash中，这一行的操作，在数值行出现之前完成
  	  		$keyhash->{$ColNb}=$eachCell;  #print "\$keyhash->{$ColNb}=\$eachCell=$keyhash->{$ColNb}\n";
  	  	}
  	  	elsif ($isThisLineValueLine==1){        #如果是数值行，  则将该数值放入 输出$outHash中，1级键是 $keyColNb_1决定的第几列的 值，2级键是 $keyColNb_2决定的第几列的 值， 3级键是 该列的关键字 ， 值是该行该列的值
  	  		
  	  		$tempHash->{ $keyhash->{$ColNb} }=$eachCell;  	
  	  		
  	  		  		 
  	  	}
  	  	
  	    $ColNb++;
  	    
  	  }
  	  
  	  if (   (  defined ( $tempHash  )  ) && (  ref ( $tempHash  ) eq 'HASH'  )    ){
  	  	
  	  	#if (   defined (  $outHash->{ $tempAr[$keyColNb_1] }->{ $tempAr[$keyColNb_2] } }  )   ){ #如果该值已经被定义过，说明 $outHash->{ $tempAr[$keyColNb_1] }->{ $tempAr[$keyColNb_2] }下的数组已经有值
  	  	#	  DieWork::Just_dieWork ( $die_MsgHead."$inFile=$inFile,\n ID: \$tempAr[\$keyColNb_1]=\$tempAr[$keyColNb_1]=$tempAr[$keyColNb_1] is not unique!! \n or KeyWord:\$tempAr[\$keyColNb_2]=\$tempAr[$keyColNb_2]=$tempAr[$keyColNb_2] is not unique!!  \n or KeyWord:\$keyhash->{\$ColNb}=\$keyhash->{$ColNb}=$keyhash->{$ColNb} is not unique!!\n\nPlease check the file: $inFile".$caller_inform );
  	  		  
  	    #}
  	    push @{  $outHash->{ $tempAr[$keyColNb_1] }->{ $tempAr[$keyColNb_2] }  }, $tempHash;  	
  	    	
  	  }
  	  
  	  
  	  
  	  
  	  
  	}    
    $rowNB++;
  }
  close (IN);
  $/=$tempMakr;
  return $outHash;
}


1;



##########################################################################################################################################
# 整理所有用于后续分析的输入文件