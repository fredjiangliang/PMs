#!/usr/bin/perl -w
BEGIN{  push (@INC, '/home/fredjiang/fredProgs/PMs')  };


use strict;
use warnings;
use Bio::Seq;
use Bio::SeqIO;
use File::Spec;
use Data::Dumper;
#use DieWork;
$Data::Dumper::Purity=1;
$Data::Dumper::Sortkeys = sub {
    my $data = join '', keys %{$_[0]};
    if ($data =~ /[A-Za-z]/) {           # input is not numeric
        return [ sort { $a cmp $b } keys %{$_[0]}];
    } 
    elsif($data =~ /^\s*\d+\.?\d*\s*$/) {
        return [ sort { $a <=> $b } keys %{$_[0]} ];
    }
    else{
    	  return [ sort { $a cmp $b } keys %{$_[0]}];
    }
};


package DirFileHandle;

#sub getDirArray{      #获取参数指向的文件夹地址中的所有子文件和子文件夹的名字，组成数组，并取地址，返回
                       #输入: 1，用于获取子文件的 文件夹名
                       #输出：1，子文件名 组成的 数组的引用


#sub PrintDumper{      # 将多位 数据结构保存到一个文件，用于下个程序 读取使用
                       #输入: 1，“需要生成的 hash Dumperfile 的名字”， 2，“需要被Dump的Hash的引用”           ; 
                       #返回: 无                  ,    
                       #写入: 1，是“hash dumperfile”，2是 “hash dumperfile的 可读文件”  
                       
                       
#sub UniqSeqForInFilse{      ##去掉重复序列， 去重复的依据是  $seqObj_1->primary_id
                             #输入: 1，“需要去重复的序列的文件名”， 2，“需要去重复的序列的文件的格式” ，    3，“去重复后的序列的文件名”，4，“去重复后的序列的文件的格式” ，          ; 
                             #返回: 无                  ,    
                             #写入: 1，“去重复后的序列的文件”                        

#sub BuildSeqID_to_clsassifyHash{  #把 $seqObj->primary_id提出来，建一个HASH

#sub ChangeStarToX_inSequnce{     #把 所有序列中的*变成X

#sub UniqSeq_and_ChangeStarToX_ForInFilse{  #对输入的文件中的序列进行去重复，然后将*转化成X


#sub Mkdir_and_warn{  #DirFileHandle::BuildWarnDieHeadInform ($workingDIRpath);
#	my ($workingDIRpath)=@_;
#	my ($warnMsgBody, $warnMsgHead, $die_MsgHead, $caller_inform)=@{ DieWork::BuildWarnDieHeadInform( 'DirFileHandle', 'Mkdir_and_warn' ) };
  
#  DieWork::Check_DfdNoEmptString_or_DIE( $workingDIRpath,      "\$workingDIRpath",     $die_MsgHead, $caller_inform  );  
	
#	my $mkDir_wkingDirPath="mkdir -p $workingDIRpath";  	DieWork::Print_and_warn( "\n$warnMsgBody$mkDir_wkingDirPath\n\n");
#	system ( "$mkDir_wkingDirPath");
	
#}
  


sub AbsPath_changedInto_RelativePath{
	
	my ($warnMsgBody, $warnMsgHead, $die_MsgHead, $caller_inform)=@{ DieWork::BuildWarnDieHeadInform( 'DirFileHandle', 'AbsPath_changedInto_RelativePath' ) };   
	
	
	
}


sub AutoSortKey{
	my ($inHASH)=@_; my $allNborNot=0;
	FCMKNB: foreach my $ecKey (  keys ( %{ $inHASH } )  ){
		if ($ecKey=~m/^\s*(\+|-)?\d+(\.\d+)?\s*$/){	$allNborNot=1;		}
		else { $allNborNot=0; last FCMKNB;}
	}

  if ($allNborNot) {           
      return [   sort { $a <=> $b } (  keys ( %{ $inHASH } )  )   ];
  }  
  else{  # input is not numeric
  	  return [   sort { $a cmp $b } (  keys ( %{ $inHASH } )  )   ];
  }
}

sub getDirArray{     #my $outVal=DirFileHandle::getDirArray ($inputDir); #获取参数指向的文件夹地址中的所有子文件和子文件夹的名字，组成数组，并取地址，返回
	my($inputDir)=@_;
  opendir (INDIR,$inputDir) or die "in sub getDirArray, cannot opendir \$inputDir=$inputDir : $!\n\n";
  my @outDirFileArray = sort ( grep {!/^\.{1,2}\z/} readdir INDIR);                 #my $nb=0; map {$nb++; printf "%4d\t%30s\n", $nb, $_ ; } @outDirFileArray; print "\n";
  close (INDIR);
  my $outVal= \@outDirFileArray;
  
  return $outVal;
}


sub Warn__SubCallerInform{  #DirFileHandle::WarnSubCallerInform
	
	my $space="                      ";
	warn "\n\n\n\n";
	warn "HEADHEADHEADHEADHEADHEADHEADHEADHEAD       Show the caller information         HEADHEADHEADHEADHEADHEADHEADHEADHEAD\n";
  my $i=0; while (  my @callerArray = caller ($i++)  ){my $relI=$i-1;
  	if  ($relI>=0){
  	  warn "\nCaller  level :$relI:\n${space}Package:( $callerArray[0] )\n${space}File   :( $callerArray[1] )\n${space}Line   :( $callerArray[2] )\n${space}Sub    :( $callerArray[3] )\n\n${space}calleArray left:";
  	  for (my $j=4; $j<@callerArray; $j++) {
  	  	warn "${space}[$j]$callerArray[$j]" if (   (  defined ( $callerArray[$j] )  ) && ( $callerArray[$j]=~m/\S+/ )   );
  	  }
  	  warn "\n\n";
  	}
  }
  warn "TAILTAILTAILTAILTAILTAILTAILTAILTAIL       Show the caller information         TAILTAILTAILTAILTAILTAILTAILTAILTAIL";
  warn "\n\n\n\n";
	
}

sub print_SubCallerInform{  #DirFileHandle::print_SubCallerInform
	
	my $outInform;
	my $space="                      ";
	$outInform.="\n\n\n\n"; 
	$outInform.= "HEADHEADHEADHEADHEADHEADHEADHEADHEAD       Show the caller information         HEADHEADHEADHEADHEADHEADHEADHEADHEAD\n";  
  my $i=0; while (  my @callerArray = caller ($i++)  ){my $relI=$i-1;
  	if  ($relI>=0){
  	  $outInform.= "\nCaller  level :$relI:\n${space}Package:( $callerArray[0] )\n${space}File   :( $callerArray[1] )\n${space}Line   :( $callerArray[2] )\n${space}Sub    :( $callerArray[3] )\n\n${space}calleArray left:";
  	  for (my $j=4; $j<@callerArray; $j++) {
  	  	$outInform.= "${space}[$j]$callerArray[$j]" if (   (  defined ( $callerArray[$j] )  ) && ( $callerArray[$j]=~m/\S+/ )   );
  	  }
  	  $outInform.= "\n\n";
  	}
  }
  $outInform.= "TAILTAILTAILTAILTAILTAILTAILTAILTAIL       Show the caller information         TAILTAILTAILTAILTAILTAILTAILTAILTAIL";
  $outInform.= "\n\n\n\n";
	#print $outInform;
	return $outInform;
}

##sub     PrintDumper{    # 将多位 数据结构保存到一个文件，用于下个程序 读取使用
sub PrintDumper{    # 将多位 数据结构保存到一个文件，用于下个程序 读取使用
  my ($fileName, $inputHashName)=@_;  
  
  
  #print "\n\nIn package DirFileHandle In sub PrintDumper\n\$fileName=$fileName, \$inputHashName=$inputHashName\n\n";
  #&Warn__SubCallerInform;
  my $caller_inform=DirFileHandle::print_SubCallerInform;  #warn "\$inputHashName=$inputHashName\n\$caller_inform=$caller_inform\n\n";
  
  Storable::store $inputHashName,"$fileName"; 

  open (DUMPFILE,">${fileName}.read.txt") or die "cannot create \${fileName}.read.txt=${fileName}.read.txt : $!\n\n"; 
  print DUMPFILE Data::Dumper::Dumper ($inputHashName);  
  close (DUMPFILE);
  #warn "Created \$fileName=$fileName ,\${fileName}.read.txt=${fileName}.read.txt, 2 new hash file!\n\n"; 	
  print "Created \$fileName=$fileName ,\${fileName}.read.txt=${fileName}.read.txt, 2 new hash file!\n\n"; 
}

# DirFileHandle::PrintDumper_plReadFile ($fileName, $inputHashName)  if (    ( defined ( $inputHashName ) ) && (   (  ref ( $inputHashName ) eq 'HASH'  ) || (  ref ( $inputHashName ) eq 'ARRAY'  )   )    );
# DirFileHandle::PrintDumper_plReadFile ($fileName, $inputHashName)  if (    ( defined ( $inputHashName ) ) && (  ref ( $inputHashName ) eq 'HASH'  )    );
# DirFileHandle::PrintDumper_plReadFile ($fileName, $inputHashName)  if (    ( defined ( $inputHashName ) ) && (  ref ( $inputHashName ) eq 'ARRAY'  )   );
##sub     PrintDumper_plReadFile{    # 将多位 数据结构保存到一个文件，用于下个程序 读取使用
sub PrintDumper_plReadFile{    # 将多位 数据结构保存到一个文件，用于下个程序 读取使用
  my ($fileName, $inputHashName, $readFileName)=@_;  
  
  my ($warnMsgBody, $warnMsgHead, $die_MsgHead, $caller_inform)=@{ DieWork::BuildWarnDieHeadInform( 'DirFileHandle', 'PrintDumper_plReadFile' ) };
  
  DieWork::Check_DfdNoEmptString_or_DIE( $fileName, "\fileName", $die_MsgHead, $caller_inform  );
  
  if    (  DieWork::Check_Hash_or_NOT ( $inputHashName )  ){}
  elsif (  DieWork::Check_Array_or_NOT( $inputHashName )  ){}
  else{
  	DieWork::Just_dieWork( $die_MsgHead."\n \$inputHashName=$inputHashName should be a defined HASH or ARRAY ref  !!  $!\n\n\n".$caller_inform );
  }
  
  
  Storable::store $inputHashName,"$fileName"; 
  
  my $real_readFileName="${fileName}.read.txt.plx"; if (  DieWork::Check_DfdNoEmptString_or_NOT( $readFileName )  ){ $real_readFileName=$readFileName; }
  
  
  open (DUMPFILE,">$real_readFileName") or die "cannot create \$real_readFileName=$real_readFileName : $!\n\n"; 
  print DUMPFILE Data::Dumper::Dumper ($inputHashName);  
  close (DUMPFILE);
  warn  "Created \$fileName=$fileName ,\${fileName}.read.txt=${fileName}.read.txt, 2 new hash file!\n\n"; 	
  print "Created \$fileName=$fileName ,\${fileName}.read.txt=${fileName}.read.txt, 2 new hash file!\n\n"; 
}

sub PrintAndWarnDumper{    # 将多位 数据结构 print 和 warn。    DirFileHandle::PrintAndWarnDumper ($inputHashName, $addtionInformAtion);
  my ($inputHashName, $addtionInformAtion)=@_; 
  my $warnPrintMsgHead="\n\n\n\nIn package DirFileHandle,\nIn sub PrintAndWarnDumper,\n\n";   #print $warnPrintMsgHead;    #warn $warnPrintMsgHead;  
  
  my $outDumpString=Data::Dumper::Dumper ($inputHashName);  
  my $warnMsg="$warnPrintMsgHead\n\$inputHashName=$inputHashName\n\n$outDumpString\n\n\n";
  
  if (    (  defined ( $addtionInformAtion )  ) && ( $addtionInformAtion=~m/\S+/ )   ){
  	$warnMsg=$addtionInformAtion."\n".$warnMsg;
  }
  
  print    $warnMsg;
  warn     $warnMsg;
  return   $warnMsg;
 
  
}

sub ReturnDumperInform{  #  my $warnMsg= DirFileHandle::ReturnDumperInform ($inputHashName, $addtionInformAtion);
  my ($inputHashName, $addtionInformAtion)=@_; 
  my $outDumpString=Data::Dumper::Dumper ($inputHashName);    
  
  my $warnMsg="\n$outDumpString\n\n\n";
  if (    (  defined ( $addtionInformAtion )  ) && ( $addtionInformAtion=~m/\S+/ )   ){  
  	$warnMsg=$addtionInformAtion."\n".$warnMsg;
  }
  $warnMsg="\n\nStartLine--------ReturnDumperInformReturnDumperInformReturnDumperInformReturnDumperInform\n\n$warnMsg\n\nEnd__Line--------ReturnDumperInformReturnDumperInformReturnDumperInformReturnDumperInform\n\n";
  return   $warnMsg;
}


sub UniqSeqForInFilse{   #去掉重复序列， 去重复的依据是  $seqObj_1->primary_id
  my ($inseqFile, $inFormat, $outSeqFile, $outFormat)=@_;
  
  
  my $onceCountHash;
  
  my $seqIOobj_IN_1=Bio::SeqIO->new(-file   => $inseqFile,    
                                  -format => $inFormat   );
  
  my $seqNb_1=0;                             
  while (my $seqObj_1=$seqIOobj_IN_1->next_seq){  #warn "\$seqObj_1->primary_id=",$seqObj_1->primary_id,"\n";
  	if (   defined (  $onceCountHash->{ $seqObj_1->primary_id }  )   ){}
  	else{
  	  $onceCountHash->{ $seqObj_1->primary_id }=$seqNb_1;
  	}
  	$seqNb_1++;
  }
  
  my $seqIOobj_IN_2=Bio::SeqIO->new(-file   => $inseqFile,    
                                    -format => $inFormat       );
  my $seqIOobj_OUT =Bio::SeqIO->new(-file   => ">$outSeqFile",    
                                    -format => $outFormat      );
                                 
  my $seqNb_2=0;                             
  while (my $seqObj_2=$seqIOobj_IN_2->next_seq){   #warn "\$seqObj_2->primary_id=", $seqObj_2->primary_id, "\n\$onceCountHash->{ $seqObj_2->primary_id }=" , $onceCountHash->{ $seqObj_2->primary_id } , "\n\$seqNb_2=$seqNb_2\n";
    if ($onceCountHash->{ $seqObj_2->primary_id } ==$seqNb_2){
      $seqIOobj_OUT->write_seq($seqObj_2);   
    }
    $seqNb_2++;
  }
 

}


sub BuildSeqID_to_clsassifyHash{  #把 $seqObj->primary_id提出来，建一个HASH
  my ($inseqFile, $inFormat, $classifyKey, $inHash)=@_;
  my $seqIOobj_IN =Bio::SeqIO->new(-file   => $inseqFile,    
                                   -format => $inFormat       );
  my $outHash;
  if (  defined ($inHash)  ){ $outHash = $inHash ; }
  while (my $seqObj=$seqIOobj_IN->next_seq){   #warn "\$seqObj->primary_id=", $seqObj->primary_id, "\n\$onceCountHash->{ $seqObj->primary_id }=" , $onceCountHash->{ $seqObj->primary_id } , "\n\$seqNb_2=$seqNb_2\n";
    $outHash->{ $seqObj->primary_id } = $classifyKey;
  }
  return $outHash;
}

sub ChangeStarToX_inSequnce{     #把 所有序列中的*变成X
  my ($inseqFile, $inFormat, $outSeqFile, $outFormat)=@_;
  
  my $seqIOobj_IN =Bio::SeqIO->new(-file   => $inseqFile,    
                                   -format => $inFormat       );
  my $seqIOobj_OUT=Bio::SeqIO->new(-file   => ">$outSeqFile",    
                                   -format => $outFormat      );
                                 
  while (my $seqObj=$seqIOobj_IN->next_seq){   #warn "\$seqObj->primary_id=", $seqObj->primary_id, "\n\$onceCountHash->{ $seqObj->primary_id }=" , $onceCountHash->{ $seqObj->primary_id } , "\n\$seqNb_2=$seqNb_2\n";
    my $orgSequnce=$seqObj->seq();
    $orgSequnce=~s/\*/X/g;
    $seqObj->seq($orgSequnce);
    
    $seqIOobj_OUT->write_seq($seqObj); 
  }
  
}


sub UniqSeq_and_ChangeStarToX_ForInFilse{  #对输入的文件中的序列进行去重复，然后将*转化成X
  my ($inseqFile, $inFormat, $outSeqFile, $outFormat)=@_;
  
  
  my $onceCountHash;
  
  my $seqIOobj_IN_1=Bio::SeqIO->new(-file   => $inseqFile,    
                                  -format => $inFormat   );
  
  my $seqNb_1=0;                             
  while (my $seqObj_1=$seqIOobj_IN_1->next_seq){  #warn "\$seqObj_1->primary_id=",$seqObj_1->primary_id,"\n";
  	if (   defined (  $onceCountHash->{ $seqObj_1->primary_id }  )   ){}
  	else{
  	  $onceCountHash->{ $seqObj_1->primary_id }=$seqNb_1;
  	}
  	$seqNb_1++;
  }
  
  my $seqIOobj_IN_2=Bio::SeqIO->new(-file   => $inseqFile,    
                                    -format => $inFormat       );
  my $seqIOobj_OUT =Bio::SeqIO->new(-file   => ">$outSeqFile",    
                                    -format => $outFormat      );
                                 
  my $seqNb_2=0;                             
  while (my $seqObj_2=$seqIOobj_IN_2->next_seq){   #warn "\$seqObj_2->primary_id=", $seqObj_2->primary_id, "\n\$onceCountHash->{ $seqObj_2->primary_id }=" , $onceCountHash->{ $seqObj_2->primary_id } , "\n\$seqNb_2=$seqNb_2\n";
    if ($onceCountHash->{ $seqObj_2->primary_id } ==$seqNb_2){
    	
    	my $orgSequnce=$seqObj_2->seq();
      $orgSequnce=~s/\*/X/g;
      $seqObj_2->seq($orgSequnce);    	
      $seqIOobj_OUT->write_seq($seqObj_2);   
    }
    $seqNb_2++;
  }
 

}

#my $outHASH=DirFileHandle::BuildTreeLevel_DIR_forAGroup_ofFILES($inFile_HASH_ARRAY, $inPut_DirPath);
sub BuildTreeLevel_DIR_forAGroup_ofFILES{
	my ($inFile_HASH_ARRAY, $inPut_DirPath)=@_;
	
	my ($warnMsgBody, $warnMsgHead, $die_MsgHead, $caller_inform)=@{ DieWork::BuildWarnDieHeadInform( 'DirFileHandle', 'BuildTreeLevel_DIR_forAGroup_ofFILES' ) };
	
	my $threeHERE=3;
	
	my $working_ARRAY;
	if    (  DieWork::Check_Array_or_NOT ( $inFile_HASH_ARRAY )  ){
		$working_ARRAY=$inFile_HASH_ARRAY;
	}
	elsif (  DieWork::Check_Hash_or_NOT( $inFile_HASH_ARRAY )  ){
		$working_ARRAY=ArrayHashChange::Change_Hash_to_Array ($inFile_HASH_ARRAY) ;
	}
	else{
		DieWork::Just_dieWork( $die_MsgHead."\n $inFile_HASH_ARRAY=$inFile_HASH_ARRAY be a HASH or ARRAY ref of infiles!!  $!\n\n\n".$caller_inform );	
	}
  
  DieWork::Check_DfdNoEmptString_or_DIE( $inPut_DirPath, "\$inPut_DirPath", $die_MsgHead, $caller_inform  );
  
  
  my $howManyFiles=@{ $working_ARRAY };                           DieWork::Print_and_warn( "\$howManyFiles=$howManyFiles\n" ); 
  
  my $realWeishu=MatrixCsvChange::Get_Weishu($howManyFiles);      DieWork::Print_and_warn( "\$realWeishu=$realWeishu\n" ); 
  
  my $howManyThree=int( $realWeishu/$threeHERE );                 DieWork::Print_and_warn( "\$howManyThree=$howManyThree\n" ); 
  my $RemainderYushu= $realWeishu % $threeHERE;                   DieWork::Print_and_warn( "\$RemainderYushu=$RemainderYushu\n" );  
  if ( $RemainderYushu >0 ){ $howManyThree++};                    DieWork::Print_and_warn( "\$howManyThree=$howManyThree\n" ); 
  my $sprftWeishuNumber=$howManyThree*$threeHERE;                 DieWork::Print_and_warn( "\$sprftWeishuNumber=$sprftWeishuNumber\n" ); 
  
  my $outHASH;
  
  my $index_of_file=1;
  foreach my $eachFile (  @{ $working_ARRAY }  ){  	
  	my $sptfIdx=sprintf ("%0${sprftWeishuNumber}d", $index_of_file); 
  	
  	my $threeNbsAraay=String_Work::Chang_String_into_subStr_array($sptfIdx, $threeHERE);
  	my $joined_Path=join "/", @{ $threeNbsAraay };
  	my $the_final_path=$inPut_DirPath."/$joined_Path";            #DieWork::Print_and_warn( "\$the_final_path=$the_final_path\n" );    
  	my $mkdirCMD=" mkdir -p $the_final_path";
  	if (  -d ($the_final_path)  ) {} else {   DieWork::Print_and_warn( $mkdirCMD."\n" );  system ( "$mkdirCMD" ); }
  	
  	$outHASH->{ $eachFile }=$the_final_path;
  	
    $index_of_file++;	
  }
  
  return $outHASH;
  
}



1;