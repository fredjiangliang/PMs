
#!/usr/bin/perl -w

BEGIN{  push (@INC, '/home/fredjiang/fredProgs/PMs')  };

use strict;
use warnings;

use PrintSubArrayHash;


package  Cap3handle;


sub getPosandTypeandCharForEachEst_from_contigPos{
	my $msgHead="\n\n\nNow In package  Cap3handle,\nIn Sub getPosandTypeandCharForEachEst_from_contigPos\n\n";
	my ($ctgHash, $EstID, $ctgPos, $ctgZoF)=@_;
	my $ctgSubHash=$ctgHash->{'_Consensus'};
	if (  ( defined ($ctgZoF) ) && ( ($ctgZoF eq '+') || ($ctgZoF eq '-') )  ){			}
	else {$ctgZoF = '+';}
	my $realCtgPos=$ctgPos; if ($ctgZoF eq '-'){$realCtgPos=$ctgSubHash->{'whoLength_2len'}-$ctgPos;}
	if (   (  defined ( $ctgHash->{'_Align'}->{$EstID} )  ) && (  ref ( $ctgHash->{'_Align'}->{$EstID} ) eq 'HASH'  )   ){
		my @keyAr=(  keys ( $ctgHash->{'_Align'}->{$EstID} )  );
		my $Zf_nb=@keyAr; if ($Zf_nb==1){}else{my $dieMsg="\n\nDIE!!!!\n$msgHead\n\$Zf_nb=$Zf_nb should be 1 !!!!\n\n\n";}
		my $EstZF=$keyAr[0];
		my $EstSubHASH=$ctgHash->{'_Align'}->{$EstID}->{$EstZF};
		
		my ($outType, $outPos)=@{ $EstSubHASH->{'_Wpos2DifTypePos'}->{ $ctgSubHash->{'_DifTypePos2Wpos'}->{'_NactgType'}->{$realCtgPos} } };
		my $outWord=$EstSubHASH->{'_WlPos2Word'}->{ $ctgSubHash->{'_DifTypePos2Wpos'}->{'_NactgType'}->{$realCtgPos} };
		return [$outPos, $outWord, $outType, $EstZF];
	}
}


#sub10.2             GetCtgIdAndOrgIdsFromEstContig{  #从EstID中，把 Contig的name  和 拼接的原始Est的name都拿到
sub GetCtgIdAndOrgIdsFromEstContig{  #从EstID中，把 Contig的name  和 拼接的原始Est的name都拿到
  my ($inCtgID)=@_;
  if( $inCtgID  =~m/^(\d+(:?ctg|sig)\d+)\s+([^,]+(:?,[^,]+)*)$/  ){
  	my $ctgSigID=$1;
  	my $OrgEstIds=$2;
  	my @OrgEstIds= ($OrgEstIds=~m/([^,])/g);
  	return [ $ctgSigID, [@OrgEstIds] ];
  }
  else {
    die "In sub GetCtgIdAndOrgIdsFromEstContig, The \$inCtgID is in wrong format\n\$inCtgID=$inCtgID\n\n";
  }
}
#sub10.3             GetFileNumberFromCap3EstID{       #从Cap3的拼接contig在结果表格中的id中，找到其原来cap3生成的文件的名字
sub GetFileNumberFromCap3EstID{       #从Cap3的拼接contig在结果表格中的id中，找到其原来cap3生成的文件的名字
	my ($cap3EstID)=@_;
	if ($cap3EstID=~m/^(\d+)(:?ctg|sig)(\d+)$/){
		my $contigNB=$1;
		my $contigType=$2;
		my $ctgFileNB=$3;
		my $FileNameNumber="Ctg$ctgFileNB";
		return [$FileNameNumber,$contigNB, $contigType];
	}
	else {die"In sub GetFileNumberFromCap3EstID, The \$cap3EstID is not in correct pattern,\n\$cap3EstID=$cap3EstID\n\n";}  
}

  #下面这些key的关键字，将用在以下各个hash结构中。  由于这些关键字，会出现在 文件夹名等字符串中，需要用eq等判断语句进行判断，故不要在前面加_
#  my ($OneFiledGenoKeyWord, $EstKeyWord, $OrgEstKeyWord, $ctgTypeKeyWd, $sigTypeKeyWd)
#    =('OneFiledGeno',       'Est',       'OrgEst',       'ctg',         'sig'         ); 
 
if(0) {
#sub {
#  my $UnHandledEstID=           $IPSR_sheet_rowColHash->{$IPSR_sortedRowKey}->{$Est}->[0]  ;
#  my ($EstID, $orgEstIDarr)= @{ &GetCtgIdAndOrgIdsFromEstContig ($UnHandledEstID) };          #获得  1Est的contig或Single的ID, 以及 2拼接成该contig或single的原始est的id组成的数组的地址
#  my $EstFrm  =&GetLableOfHpLk( $IPSR_sheet_rowColHash->{$IPSR_sortedRowKey}->{$Efr}->[0] ); 
#  my $EstPos  =&GetLableOfHpLk( $IPSR_sheet_rowColHash->{$IPSR_sortedRowKey}->{$Eps}->[0] );
#  $tga_1_PosHash->{$specise}->{$EstKeyWord}->{$EstID}->{$EstFrm}->{$EstPos}=1;              ######$tga_1_PosHash##########################
#  
#  ### 下面处理Est的拼接问题，找出所有拼接成contig的原始est
#  my ($CtgFile, $WhichCtg, $CtgFileType)=@{ &GetFileNumberFromCap3EstID($EstID) };  #得到  1拼接后的est的ctgFile名，  2是这个拼接文件中的第几个contig， 3类型是contig还是single
#  my $CtgFilePath=$orgEstCap3ResultPathHash->{$specise}->{$CtgFile};                #由以上的ctgFile名，直接获取其文件绝对路径，用来进行下面的 文件prase步骤               #warn "\$CtgFileType=$CtgFileType\t\t\$ctgTypeKeyWd=$ctgTypeKeyWd\t\t\$CtgFile=$CtgFile\n\$CtgFilePath=$CtgFilePath=\$orgEstCap3ResultPathHash->{$specise}->{$CtgFile}=\$orgEstCap3ResultPathHash->{\$specise}->{\$CtgFile}\n";
#  if ($CtgFileType eq $ctgTypeKeyWd){                                               #只有 类型是 contig这种，才需要进行下面的分析
#  	my $parsedCapHash=&parseCapOut($CtgFilePath, $parseCapHashNeedKeys);            #解析cap3文件，输入为（该cap3的路径，提前设定的各种关键字），解析结果放入  $parsedCapHash, 该hash的数据结构，需详细读 下面的sub   &parseCapOut #returnedHash
#  	my $wholePos;           #定义一个变量，用来装没有空格和横杆的 contig序列中，tga的位置，frame是正的时候不需要变化，为负的时候需要一步变化。
#  	if ($EstFrm>0){
#  		$wholePos=$parsedCapHash->{$WhichCtg}->{'_Consensus'}->{'_DifTypePos2Wpos'}->{'_NactgType'}->{$EstPos};
#  	}
#  	elsif ($EstFrm<0){
#  		$wholePos=$parsedCapHash->{$WhichCtg}->{'_Consensus'}->{'_NactgTypeLengthKw'}+1-$parsedCapHash->{$WhichCtg}->{'_Consensus'}->{'_DifTypePos2Wpos'}->{'_NactgType'}->{$EstPos};
#  	}
#  	else {
#  	  die "\$EstFrm should be -1 -2 -3 or +1 +2 +3! Here, it is \$EstFrm=$EstFrm\n\n";
#  	}
#  	#下面，得到所有原始orgEst的tga对应的位置和拼接序列，用于SECIS检索
#  	foreach my $eachEstFullName (    sort {$a<=>$b} (   keys (  %{ $parsedCapHash->{$WhichCtg}->{'_Align'} }  )   )    ){
#  	  foreach my $drect (   sort  {$a<=>$b} (   keys (  %{ $parsedCapHash->{$WhichCtg}->{'_Align'}->{$eachEstFullName} }  )   )    ){
#  	    my ($wordType,$wordpos)=@{ $parsedCapHash->{$WhichCtg}->{'_Align'}->{$eachEstFullName}->{$drect}->{'_Wpos2DifTypePos'}->{$wholePos} };
#  	    #下面依据Tga的T对应的各个Est位置上的 字符类型$wordType 和 字符位置$wordpos，来获得Est相应位置在 Est序列本身序列上的 位置
#  	    my $finalTpos;  #定义 作为输出的 tga位置变量
#  	    if     ($wordType eq '_NactgType'){                     	    	$finalTpos=$wordpos;              	    }  #如果该位置是一个碱基字符，则直接获得Tga位置信息  
#  	    else{  #否则继续进行下面的操作
#  	      warn "The \$wordType=$wordType is not a nucleotide!\n\n\$CtgFilePath=$CtgFilePath,\n\$WhichCtg=$WhichCtg,\$eachEstFullName=$eachEstFullName,\$drect=$drect,\$wholePos=$wholePos\n\n\n";
#  	      print "The \$wordType=$wordType is not a nucleotide!\n\n\$CtgFilePath=$CtgFilePath,\n\$WhichCtg=$WhichCtg,\$eachEstFullName=$eachEstFullName,\$drect=$drect,\$wholePos=$wholePos\n\n\n";
#  	      
#  	      #如果相应位置不是碱基字符，而是 横杠 或 空格，则先向左边移动找到左边最近的 碱基字符，或移到位置1处为止。 或者 向右移 找到右边最近的 碱基字符，或移动到最后最长的位置处为止
#  	      my ($LeftGoodPos ,$rightGoodPos); #这两个变量，用来获得向左移动或向右移动最近的那个是碱基的位置的，Est纯碱基序列的位置
#  	      #先向左走
#  	      my $leftStep=0; my $leftProcessPos=$wholePos; my $stopSign=0;  my $leftGood=0;  my $matchNactgTypeLeftStep=0;    #my $matchNactgTypeLeftStep=0;   这个变量用来记录左移过程中，空格或横杠match的consensus序列上有几个碱基
#  	      while ($stopSign==0){
#  	      	$leftStep++;    $leftProcessPos--;   #将左走的步数$leftStep记录下来，把左走后所在的位置记录在$leftProcessPos中。
#  	      	my ($tempTp, $tempPos)=@{ $parsedCapHash->{$WhichCtg}->{'_Align'}->{$eachEstFullName}->{$drect}->{'_Wpos2DifTypePos'}->{$leftProcessPos} };  #获得该位置的字符类型$tempTp和在全序列长度上的位置数$tempPos。
#  	        if ( ($tempTp eq '_NactgType') || ($leftProcessPos == 1) ){  ## 依据步进移动后的位置的字符类型，进行相应的判断和记录。看是不是处于序列的开头，或者判断是否在中间。并最终得到 移动后的具体位置。
#  	          $stopSign=1;
#  	          if ($tempTp eq '_NactgType') { $leftGood=1;  $LeftGoodPos=$tempPos;}
#  	        }
#  	        #接下来，对$matchNactgTypeLeftStep 这一值进行计算，其方法是判断每个移动步骤中，相应的consenesus位置上是不是 '_NactgType'类型，是则 +1, 否则不变
#  	        if ($parsedCapHash->{$WhichCtg}->{'_Consensus'}->{'_Wpos2DifTypePos'}->{$leftProcessPos}->[0] eq '_NactgType'){$matchNactgTypeLeftStep++;}  
#  	      }          
#  	      #再先向右走
#  	      my $rightStep=0; my $rightProcessPos=$wholePos; my $stopSign=0;  my $rightGood=0;  my $matchNactgTypeRightStep=0; #my $matchNactgTypeRightStep=0; 这个变量用来记录右移过程中，空格或横杠match的consensus序列上有几个碱基
#  	      while ($stopSign==0){
#  	      	$rightStep++;    $rightProcessPos++;
#  	      	my ($tempTp, $tempPos)=@{ $parsedCapHash->{$WhichCtg}->{'_Align'}->{$eachEstFullName}->{$drect}->{'_Wpos2DifTypePos'}->{$rightProcessPos} };
#  	        if ( ($tempTp eq '_NactgType') || ($rightProcessPos == $parsedCapHash->{$WhichCtg}->{'_Consensus'}->{'_NactgTypeLengthKw'}) ){
#  	          $stopSign=1;
#  	          if ($tempTp eq '_NactgType') { $rightGood=1;  $rightGoodPos=$tempPos;}
#  	        }
#  	        if ($parsedCapHash->{$WhichCtg}->{'_Consensus'}->{'_Wpos2DifTypePos'}->{$rightProcessPos}->[0] eq '_NactgType'){$matchNactgTypeRightStep++;}  
#  	      }
#  	      
#  	      if ( ($leftGood ==1) || ($rightGood==1) ){
#  	      	if    ( ($leftGood ==1) && ($rightGood==1) ){   #这种情况 说明，左右都有碱基，那么这个Tga的T的位置，对应于Est的中间某个 空格或者 横杠
#  	      		if ($leftStep> $rightStep){               	      		$finalTpos=$rightGoodPos;              	      		} 
#  	      		else {                                    	      		$finalTpos=$LeftGoodPos;             	      	  }
#    	      }
#  	      	elsif ( ($leftGood ==1) && ($rightGood==0) ){   #这种情况 说明，左边有碱基，右边没有，那么这个Tga的T的位置，对应于Est的末尾的外部，那么我们给它的位置定义为 一个 大于其长度的位置
#  	      		#那么该位置的计算方法是该 Est的长度（其实就是其最后边字符的位置) 再加上左移过程中match到的consensus上的碱基字符的数量 
#  	      		$finalTpos=$parsedCapHash->{$WhichCtg}->{'_Align'}->{$eachEstFullName}->{$drect}->{'_NactgTypeLengthKw'} + $matchNactgTypeLeftStep ;
#  	        }
#  	      	elsif ( ($leftGood ==0) && ($rightGood==1) ){   #这种情况 说明，右边有碱基，左边没有，那么这个Tga的T的位置，对应于Est的开端的外部，那么我们给它的位置定义为 一个 负数       	
#  	      	  #那么该位置的计算方法是0 减去右移过程中match到的consensus上的碱基字符的数量 
#  	      		$finalTpos=                                    0                                                         - $matchNactgTypeRightStep;
#  	        }
#  	      }  
#  	      else {die "left right all wrong!!  \n\$CtgFilePath=$CtgFilePath,\n\$WhichCtg=$WhichCtg,\$eachEstFullName=$eachEstFullName,\$drect=$drect,\$wholePos=$wholePos\n\n\n";}
#  	    }
#  	    #由Est信息，指向它的原始拼接orgEst的序列及位置方向等信息
#  	    $est_orgEst_informationHash->{$specise}->{$EstKeyWord}->{$EstID}->{$EstFrm}->{$EstPos}->{$eachEstFullName}->{$drect}->{$finalTpos}=1; 
#  	    #下面由OrgEst信息，指向的序列
#  	    $est_orgEst_informationHash->{$specise}->{$OrgEstKeyWord}->{$eachEstFullName}->{$drect}->{'_NactgTypeSeqKw'}   =$parsedCapHash->{$WhichCtg}->{'_Align'}->{$eachEstFullName}->{$drect}->{'_NactgTypeSeqKw'}; 
#  	    $est_orgEst_informationHash->{$specise}->{$OrgEstKeyWord}->{$eachEstFullName}->{$drect}->{'_NactgTypeLengthKw'}=$parsedCapHash->{$WhichCtg}->{'_Align'}->{$eachEstFullName}->{$drect}->{'_NactgTypeLengthKw'}; 
#  	    
#  	    $tga_1_PosHash->{$specise}->{$OrgEstKeyWord}->{$eachEstFullName}->{$drect}->{$finalTpos}=1;              #########$tga_1_PosHash#######################
#  	  }
#  	}
#  }
#}
  }

#sub10.11              FastBuildAddLinkForEst{      # 快速建立所有 cap3的输出文件和序列名之间的关系 数据结构  
sub FastBuildAddLinkForEst{      # 快速建立所有 cap3的输出文件和序列名之间的关系 数据结构
	my ($allSpeciseDir)=@_;   # 输入通常是 "$allSpeciseDir/$eachSpecise/OrgEst/Org_${eachSpecise}.sorted.Dir";，如#Prymnesium_parvum20150522/OrgEst/Org_Prymnesium_parvum20150522.sorted.Dir
	my $addLinkHash;
      

  foreach my $eachSpecise (  @{ DirFileHandle::getDirArray($allSpeciseDir) }  ){  #Prymnesium_parvum20150522/OrgEst/Org_Prymnesium_parvum20150522.sorted.Dir
	  warn "\$eachSpecise=$eachSpecise\n";
	  if ($eachSpecise=~m/^\S+_\S+$/){    
      my $eachCapSortedDir="$allSpeciseDir/$eachSpecise/OrgEst/Org_${eachSpecise}.sorted.Dir";
      foreach my $eachNOdir( @{ DirFileHandle::getDirArray($eachCapSortedDir) } ){
      	my $absDirOfEachNOdir="$eachCapSortedDir/$eachNOdir";
        foreach my $eachCapOPfile ( @{ DirFileHandle::getDirArray($absDirOfEachNOdir) } ){
          if ($eachCapOPfile=~m/^(\S+)\.cap\.out$/){
          	my $contigName=$1;
            my $absPathOfCapOPFile="$absDirOfEachNOdir/$eachCapOPfile";
            
            $addLinkHash->{$eachSpecise}->{$contigName}=$absPathOfCapOPFile; print "\$addLinkHash->{\$eachSpecise}->{\$contigName}=\$absPathOfCapOPFile=\$addLinkHash->{$eachSpecise}->{$contigName}=$absPathOfCapOPFile\n";
          }
          else {
            #print "This file or directory:\$eachCapOPfile=$eachCapOPfile is not a file we need!\n\n";
          }
        }    
      }
    }
  }
  return $addLinkHash;
}


sub runUsearchCap3_to_clustal_Est{  #新的cap3子例程，先用usearch聚类，再拼接
  my ($orgEstInPathFile, $outEstPathFile, $InSubUsearchPathFile, $inCap3Path)=@_; #( $OrgSpieceEst, $SpieceEstDb, $usearchPathFile, ${cap3Path})  
  
  $inCap3Path="cap3";                                                                  #一般不需要改了          cap3    路径
  $InSubUsearchPathFile="~/tools/USEARCH/usearch7.0.1090_i86linux32";                         #一般不需要改了          usearch 路径
  
  my $orgEstInPath=dirname ($orgEstInPathFile); 
  my $orgEstInFile=basename ($orgEstInPathFile);
  my $sortStdOutFile="$orgEstInPath/${orgEstInFile}.sort.stdout";
  my $sortErrorFile="$orgEstInPath/${orgEstInFile}.sort.error";
  my $sortedOrgEst="$orgEstInPath/${orgEstInFile}.sorted";  
  print "\$sortedOrgEst=$sortedOrgEst=\$orgEstInPath/\${orgEstInFile}.sorted=$orgEstInPath/${orgEstInFile}.sorted\n";
  my $usearchSortCommand="$InSubUsearchPathFile  -sortbylength $orgEstInPathFile -output $sortedOrgEst -minseqlength 1 1>$sortStdOutFile 2>$sortErrorFile"; 
  print "\$usearchSortCommand=$usearchSortCommand\n";
  system ("$usearchSortCommand");  #对orgEst进行排序
  
  my $NewDirForClustering="$orgEstInPath/${orgEstInFile}.sorted.Dir";
  system ("mkdir $NewDirForClustering");
  my $usearchStdOutFile="$orgEstInPath/${orgEstInFile}.usearch.stdout";
  my $usearchErrorFile="$orgEstInPath/${orgEstInFile}.usearch.error";
  my $clusterHead="Ctg";
  my $usearchCommand="$InSubUsearchPathFile -cluster_smallmem $sortedOrgEst -id 0.9 -clusters $NewDirForClustering/$clusterHead 1>$usearchStdOutFile 2>$usearchErrorFile";
  system ("$usearchCommand");
  
  opendir (CLUSTDIR, $NewDirForClustering) or die "cannot opendir \$NewDirForClustering=$NewDirForClustering :$!\n\n";
  my @ctgsArr = sort ( grep {!/^\.{1,2}\z/} readdir CLUSTDIR);                 #my $nb=0; map {$nb++; printf "%4d\t%30s\n", $nb, $_ ; } @ctgsArr; print "\n";
  #new#已经在文件夹1.sort.stdout中 @ctgArr中为每个文件夹的名字，结果文件夹名字中为Ctg1、Ctg2等，可利用正则表达式提取出名字中的数字，可以归类到一个文件夹中）
  #my $dir_no=1;
  #my $mk1cmd="mkdir $NewDirForClustering/$dir_no";      #当前位置在原始数据文件的文件夹内，需写相对地址,./为当前目录，非/开头就是相对地址
  #print "\$mk1cmd=$mk1cmd\n";      
  #system ($mk1cmd); #or die $!;                                               #新建文件夹已数字命名
  
  open (CAPEDESTFILE,">$outEstPathFile") or die "cannot create \$outEstPathFile=$outEstPathFile : $!\n";	                                    
  foreach my $eachCtg (@ctgsArr){ 
    print "\$eachCtg=$eachCtg\n\${clusterHead}=${clusterHead}\n";
    
    if ($eachCtg=~m/^${clusterHead}(\d+)$/){       #读取文件名中的数字，能整除500的时候再另外新建文件夹
      my $ctgNb=$1; 
      my $dir_no=($ctgNb/500);
      $dir_no=(int ($dir_no))+1;
      if (-d ("$NewDirForClustering/$dir_no") ){        print "\$NewDirForClustering/\$dir_no=\$NewDirForClustering/$dir_no=$NewDirForClustering/$dir_no is existed\n";
      }
      else {
        my $mkdircommend ="mkdir $NewDirForClustering/$dir_no";
        print "\$mkdircommend=$mkdircommend\n";
        system ($mkdircommend);
      }
      
      my $movecommand="mv $NewDirForClustering/$eachCtg $NewDirForClustering/$dir_no/$eachCtg";      #当前位置在原始数据文件的文件夹内，需写相对地址,./为当前目录，非/开头就是相对地址
      print "\$movecommand=$movecommand\n";
      system ("$movecommand");
     #移动完后直接cap3
      my $cap3Comand="${inCap3Path} $NewDirForClustering/$dir_no/$eachCtg 1> $NewDirForClustering/$dir_no/${eachCtg}.cap.out 2> $NewDirForClustering/$dir_no/${eachCtg}.cap.ero";
  	  print "\$cap3Comand=$cap3Comand\n";
  	  system ("$cap3Comand");
  	  
  	  #cap3完成，剩下几步的目的？
  	  open (ACEFILE,"$NewDirForClustering/$dir_no/${eachCtg}.cap.ace") or die "cannot open \$NewDirForClustering/\$dir_no/\${eachCtg}.cap.ace=$NewDirForClustering/$dir_no/${eachCtg}.cap.ace :$!\n";
      my $capCtgNb=1; #这个数字，是为了标示每个usearch聚类后，再用cap3拼接后，形成多少个 ctg或者sig的，如果只有1个，那这个值就总是1，如果有2个或3个，则会显示为2 或者3 或者更大。
      my $tmpMK=$/;
      $/="\nCO ";
      while (<ACEFILE>) { 
      	my $aceUnit=$_;  #print "\$aceUnit=$aceUnit\n";#Contig1 782 2 1 U 
      	if ($aceUnit=~m/^Contig(\d+)\s+\S+\s+\S+\s+\S+\s+\S+\n((?:\S+\n)+)\n/){
      	  my $contigNb=$1; 
          my $contigSeq=$2;   $contigSeq=~s/\s+//g; $contigSeq=~s/\*//g; print "\$ctgNb=$ctgNb\t\$capCtgNb=$capCtgNb\n";
      	  #RD gi|66860169|gb|DR039684.1|DR039684 467 0 0
      	  my @estIDs=($aceUnit=~m/\nRD\s+(\S+)\s+/g); 
          my $estIDsLine=join (",",@estIDs );
      	  my $contigHead=">lcl|${capCtgNb}ctg$ctgNb $estIDsLine"; #print  "$contigHead\n$contigSeq";
      	  print CAPEDESTFILE "$contigHead\n$contigSeq";
      	  $capCtgNb++;
      	} else {print "In file $NewDirForClustering/$dir_no/${eachCtg}.cap.ace, this is not a good \$aceUnit=\n$aceUnit \n";}
      }
      $/=$tmpMK;
      close (ACEFILE);
  	  
  	  open (SIGFILE,"$NewDirForClustering/$dir_no/${eachCtg}.cap.singlets") or die "cannot open  \$NewDirForClustering/\$dir_no/\${eachCtg}.cap.singlets=$NewDirForClustering/$dir_no/${eachCtg}.cap.singlets :$!\n";
      
      my $tmp2MK=$/;
      $/="\n>";  
      my $sigNB=1; #这个数字，是为了标示每个usearch聚类后，再用cap3拼接后，形成多少个 ctg或者sig的，如果只有1个，那这个值就总是1，如果有2个或3个，则会显示为2 或者3 或者更大。
      while (<SIGFILE>){
      	my $sigUnit=$_; $sigUnit=~s/>$//; $sigUnit=~s/^>//; $sigUnit=~s/\n$//; $sigUnit=~s/\s*$//g; $sigUnit="$sigUnit\n";
        if ( $sigUnit=~m/\S+/){ $sigUnit=">$sigUnit";
          if ($sigUnit=~m/^>(\S+.*)\n((?:\S+\n)+)$/){
          	my $sigHead=">lcl|${sigNB}sig$ctgNb $1";
          	my $sigSeq=$2;
      	  	print CAPEDESTFILE "$sigHead\n$sigSeq";
      	  	$sigNB++;
      	  }
      	  else {die "Wrong the \$sigUnit is :\n$sigUnit\nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn\n";}
      	}
      }
      $/=$tmp2MK;
      close (SIGFILE);
      $dir_no++; 
  	 }
   }
   close (CAPEDESTFILE);
    
   
}

sub parseCapSinglets{ # Cap3handle::parseCapSinglets
	my ($inCapSingletsFile)=@_;                                 #print "\n0000-201809051456\$inCapSingletsFile=$inCapSingletsFile\n";   #, $inPutHashKeys
	my $warnMsgBody="\nIn package  Cap3handle,\tIn sub parseCapSinglets,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	
	my $returnedHash;  #用于最终输出的 hash
  my $estCtgEstNmHash;  #用于抽提出每个contig 由哪些est拼装而成，est的名字放入 这个hash
  my ($AlignKw, $ConsensusKw, $JustSequnceKw,  $Wpos2DifTypePosrelationKw, $DifTypePos2WposrelationKw, $WlPos2WordKw, $DifTypePos2WordKw, $NactgTypeKw, $SpaceTypeKw, $GangTypeKw, $wholeLengthKw,  $NactgTypeLengthKw,  $SpaceTypeLengthKw,  $GangTypeLengthKw,    $NactgTypeSeqKw,  $SpaceTypeSeqKw,  $GangTypeSeqKw   )
  =  ('_Align', '_Consensus', 'JustSequnceKw', '_Wpos2DifTypePos',         '_DifTypePos2Wpos',         '_WlPos2Word', '_DifTypePos2Word', '_NactgType', '_SpaceType', '_GangType', '_wholeLengthKw','_NactgTypeLengthKw','_SpaceTypeLengthKw','_GangTypeLengthKw','_NactgTypeSeqKw','_SpaceTypeSeqKw','_GangTypeSeqKw' );
  #    0        1            2               3                            4                          5               6                  7             8             9            10               11                   12                   13                 14                15                 16            
  #  拼接Est， 拼接Contig，   只是序列的标示     全长和不同类型位置的转换    不同类型和全长位置的转换   全长位置字符    不同类型位置字符    原序列        空格            横杆        全长长度         原序列长度          空格总长度             横杠总长度      原序列           空格序列              横杠序列    
  
  open (INCAPFILE,$inCapSingletsFile) or die "$die_MsgHead\nCannot open \$inCapSingletsFile=$inCapSingletsFile : $!\n";
  my @incapout=(<INCAPFILE>);
  my $wholeIncapout=join ('',@incapout);
  if ($wholeIncapout=~m/
                           (?:
                              >\S+.*\n
                              (?:.*\S+.*\n)+
                              \n?
                           )+  
                         
                       /x
     ){                                                         #print "1-201809051456\$inCapSingletsFile=$inCapSingletsFile\n\n";
     my @eachSinglet=( $wholeIncapout=~m/                         
                                            (                   
                                               >\S+.*\n           
                                               (?:[^>].*\S+.*\n)+     
                                               \n?                
                                            )                    
                                                                  
  	                              	    /xg 
  	                 );
  	 my $sigNb=1;                            	                           
    foreach my $eachSigLet (@eachSinglet) {                              #print "2-201809051456\$inCapSingletsFile=$inCapSingletsFile \$eachSigLet=$eachSigLet\n";
      if ($eachSigLet=~m/
    	                        >(\S+.*)\n           
                             ( (?:.*\S+.*\n)+ )     
                             \n?
    	                    /x
    	   ){
    	  my $sigLtNam=$1; 	
    	  my $sigLtSeq=$2;  
    	  $sigLtNam=~s/^\s*//; $sigLtNam=~s/\s*$//;                      #print "3-201809051456\$inCapSingletsFile=$inCapSingletsFile \$sigLtNam=$sigLtNam\$2=$2 \n";
    	  $sigLtSeq=~s/\s+//g;                                           #print "4-201809051456\$inCapSingletsFile=$inCapSingletsFile \$sigLtSeq=$sigLtSeq\n";
    	  my $sigLtLen=length ($sigLtSeq);                               #print "5-201809051456\$inCapSingletsFile=$inCapSingletsFile \$sigLtLen=$sigLtLen\n";
    	  $returnedHash->{$sigNb}->{'1_OrderHash'}->{0}=[$sigLtNam,'+'];
    	  $returnedHash->{$sigNb}->{'1_Est0lsNoD'}=$sigLtNam;
    	  $returnedHash->{$sigNb}->{'1_Est1lsAdD'}=$sigLtNam.'_+';
    	  
  

    	  
    	  $returnedHash->{$sigNb}->{$ConsensusKw}->{$JustSequnceKw    }=$returnedHash->{$sigNb}->{$AlignKw}->{$sigLtNam}->{'+'}->{$JustSequnceKw    }=$sigLtSeq;
    	  $returnedHash->{$sigNb}->{$ConsensusKw}->{$GangTypeLengthKw }=$returnedHash->{$sigNb}->{$AlignKw}->{$sigLtNam}->{'+'}->{$GangTypeLengthKw }=0;
    	  $returnedHash->{$sigNb}->{$ConsensusKw}->{$GangTypeSeqKw    }=$returnedHash->{$sigNb}->{$AlignKw}->{$sigLtNam}->{'+'}->{$GangTypeSeqKw    }='';
    	  $returnedHash->{$sigNb}->{$ConsensusKw}->{$NactgTypeLengthKw}=$returnedHash->{$sigNb}->{$AlignKw}->{$sigLtNam}->{'+'}->{$NactgTypeLengthKw}=$sigLtSeq;
    	  $returnedHash->{$sigNb}->{$ConsensusKw}->{$NactgTypeSeqKw   }=$returnedHash->{$sigNb}->{$AlignKw}->{$sigLtNam}->{'+'}->{$NactgTypeSeqKw   }=$sigLtLen;
    	  $returnedHash->{$sigNb}->{$ConsensusKw}->{$SpaceTypeLengthKw}=$returnedHash->{$sigNb}->{$AlignKw}->{$sigLtNam}->{'+'}->{$SpaceTypeLengthKw}=0;
    	  $returnedHash->{$sigNb}->{$ConsensusKw}->{$SpaceTypeSeqKw   }=$returnedHash->{$sigNb}->{$AlignKw}->{$sigLtNam}->{'+'}->{$SpaceTypeSeqKw   }='';
    	  $returnedHash->{$sigNb}->{$ConsensusKw}->{$wholeLengthKw    }=$returnedHash->{$sigNb}->{$AlignKw}->{$sigLtNam}->{'+'}->{$wholeLengthKw    }=$sigLtLen;
    	                                                  
    	  $returnedHash->{$sigNb}->{$ConsensusKw}->{'whoLength_0stt'  }=$returnedHash->{$sigNb}->{$AlignKw}->{$sigLtNam}->{'+'}->{'whoLength_0stt'  }=0;
    	  $returnedHash->{$sigNb}->{$ConsensusKw}->{'whoLength_1end'  }=$returnedHash->{$sigNb}->{$AlignKw}->{$sigLtNam}->{'+'}->{'whoLength_1end'  }=$sigLtLen-1;
    	  $returnedHash->{$sigNb}->{$ConsensusKw}->{'whoLength_2len'  }=$returnedHash->{$sigNb}->{$AlignKw}->{$sigLtNam}->{'+'}->{'whoLength_2len'  }=$sigLtLen;
    	  
    	  my @ACGT_array=split ('', $sigLtSeq);
    	  for (my $i=0; $i<@ACGT_array; $i++ ){
    	     
    	  	$returnedHash->{$sigNb}->{$ConsensusKw}->{$DifTypePos2WordKw        }->{$NactgTypeKw}->{$i}=$returnedHash->{$sigNb}->{$AlignKw}->{$sigLtNam}->{'+'}->{$DifTypePos2WordKw        }->{$NactgTypeKw}->{$i}=$ACGT_array[$i];
    	  	$returnedHash->{$sigNb}->{$ConsensusKw}->{$DifTypePos2WposrelationKw}->{$NactgTypeKw}->{$i}=$returnedHash->{$sigNb}->{$AlignKw}->{$sigLtNam}->{'+'}->{$DifTypePos2WposrelationKw}->{$NactgTypeKw}->{$i}=$i;
    	  	$returnedHash->{$sigNb}->{$ConsensusKw}->{$WlPos2WordKw                             }->{$i}=$returnedHash->{$sigNb}->{$AlignKw}->{$sigLtNam}->{'+'}->{$WlPos2WordKw                             }->{$i}=$ACGT_array[$i];
    	  	$returnedHash->{$sigNb}->{$ConsensusKw}->{$Wpos2DifTypePosrelationKw                }->{$i}=$returnedHash->{$sigNb}->{$AlignKw}->{$sigLtNam}->{'+'}->{$Wpos2DifTypePosrelationKw                }->{$i}=[$NactgTypeKw, $i];
    	  	 
    	  }
    	}
    	else {
    		 my $dieMsg1="$die_MsgHead\n\$eachSigLet=$eachSigLet=\$eachSigLet didnot fit the regular expression!!\n\n";
  	    print $dieMsg1; die $dieMsg1;
    	}
    	$sigNb++;
    }  	                          	
  }
  else {
  	my $dieMsg="$die_MsgHead\n\$wholeIncapout=$wholeIncapout=\$wholeIncapout didnot fit the regular expression!!\n\n";
  	print $dieMsg; die $dieMsg;
  }
  return $returnedHash;
     
}

#            parseCapOut{      #解析 Cap3输出文件
sub parseCapOut{ # Cap3handle::parseCapOut     #解析 Cap3输出文件
	my ($inCapOutFile)=@_;  print "\$inCapOutFile=$inCapOutFile\n";   #, $inPutHashKeys
    	
  my $warnMsgBody="\nIn package  Cap3handle,\tIn sub parseCapOut,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";	
#  Number of segment pairs = 7140; number of pairwise comparisons = 3525
#'+' means given segment; '-' means reverse complement
#
#Overlaps            Containments  No. of Constraints Supporting Overlap
#  
#





#+' means given segment; '-' means reverse complement
#
#Overlaps            Containments  No. of Constraints Supporting Overlap
#
#
#DETAILED DISPLAY OF CONTIGS

  my $returnedHash;  #用于最终输出的 hash
  my $estCtgEstNmHash;  #用于抽提出每个contig 由哪些est拼装而成，est的名字放入 这个hash
  my ($AlignKw, $ConsensusKw, $JustSequnceKw,  $Wpos2DifTypePosrelationKw, $DifTypePos2WposrelationKw, $WlPos2WordKw, $DifTypePos2WordKw, $NactgTypeKw, $SpaceTypeKw, $GangTypeKw, $wholeLengthKw,  $NactgTypeLengthKw,  $SpaceTypeLengthKw,  $GangTypeLengthKw,    $NactgTypeSeqKw,  $SpaceTypeSeqKw,  $GangTypeSeqKw   )
  =  ('_Align', '_Consensus', 'JustSequnceKw', '_Wpos2DifTypePos',         '_DifTypePos2Wpos',         '_WlPos2Word', '_DifTypePos2Word', '_NactgType', '_SpaceType', '_GangType', '_wholeLengthKw','_NactgTypeLengthKw','_SpaceTypeLengthKw','_GangTypeLengthKw','_NactgTypeSeqKw','_SpaceTypeSeqKw','_GangTypeSeqKw' );
  #    0        1            2               3                            4                          5               6                  7             8             9            10               11                   12                   13                 14                15                 16            
  #  拼接Est， 拼接Contig，   只是序列的标示     全长和不同类型位置的转换    不同类型和全长位置的转换   全长位置字符    不同类型位置字符    原序列        空格            横杆        全长长度         原序列长度          空格总长度             横杠总长度      原序列           空格序列              横杠序列    
  
  open (INCAPFILE,$inCapOutFile) or die "In Sub &praseCapOut, Cannot open \$inCapOutFile=$inCapOutFile : $!\n";
  my @incapout=(<INCAPFILE>);
  my $wholeIncapout=join ('',@incapout);
  
  my $line1_Regular='^\s*Number\s+of\s+segment\s+pairs\s+=\s+\d+;\s+number\s+of\s+pairwise\s+comparisons\s+=\s+\d+\n';
  my $line2_Regular='\'\+\'\s+means\s+given\s+segment;\s+\'-\'\s+means\s+reverse\s+complement\n\n';
  my $line4_Regular='Overlaps\s+Containments\s+No\.\s+of\s+Constraints\s+Supporting\s+Overlap\n\n';
  my $Contig_Line='\*{10,20}\s+Contig\s+\d+\s+\*{10,20}\n';
  my $estNamePlus='\S+.*[\+-]';  
  my $DetaildHead='\nDETAILED\s+DISPLAY\s+OF\s+CONTIGS\n';
  my $markerLine='\s{22,22}(?:\s{4,4}\.\s{4,4}:){6,6}\n';
  my $eachAlign='\S.{19,19}[\+-]\s\s*\S{1,60}\s*\n';
  my $onlyLINE='\s{22,22}_{60,60}\n';
  my $consensus='consensus\s{13,13}\s*\S{1,60}\s*\n';
  
  #my $wholeOtherReg='(\*{19}\s+Contig\s+\d+\s+\*{19}\n(\s*\S+.*\n)+)?\nDETAILED\s+DISPLAYs+OFs+CONTIGS\n'; #******************* Contig 1 ********************  (\*{19}\s+Contig\s+\d+\s+\*{19}\n(\s*\S+.*\n)+)?DETAILED\s+DISPLAY\s+OF\s+CONTIGS/x
  #if ($wholeIncapout=~m/^${line1_Regular}${line2_Regular}(${line4_Regular})(\*{10,20}\s+Contig\s+\d+\s+\*{10,20}\n(?:\s*\S+.*\n)+)?\n(DETAILED\s+DISPLAY\s+OF\s+CONTIGS)/x){
  #if ($wholeIncapout=~m/^${line1_Regular}${line2_Regular}(${line4_Regular})($Contig_Line^$estNamePlus\n(?:(?:\s+$estNamePlus\s+is\s+in\s+$estNamePlus\n)*(?:^$estNamePlus\n)*)*)?\n(DETAILED\s+DISPLAY\s+OF\s+CONTIGS)/x){	
  if ($wholeIncapout=~m/^( ${line1_Regular}  
                           ${line2_Regular}  
                           ${line4_Regular}
                          )
                          (  
                            (?:$Contig_Line  
                               $estNamePlus\n 
                              (?: 
                                (?:
                                  (?:\s+$estNamePlus\s+is\s+in\s+$estNamePlus\n)*
                                  (?:$estNamePlus\n)* 
                                )* 
                              )  
                            )*  
                          )?
                          ( $DetaildHead )
                          (
                            (?:
                              $Contig_Line
                              (?:
                                $markerLine
                                (?:$eachAlign)+
                                $onlyLINE
                                $consensus
                                \n
                              )+
                            )+
                          )?
                       /x
       )
    {                                                                                                             	#print "\$1=$1\n\$2=$2\n\$3=$3\n\$4=$4\n$5=$5\n\n";  
  	my $contigSeqUnit=$2;
  	my $contigAlignUnit=$4;
  	
  	
  	my $contigNBhere;
  	my $shortedEstName;
  	
  	my @contigsSeqNameArray=( $contigSeqUnit=~m/   (
                                                       $Contig_Line  
                                                       $estNamePlus\n 
                                                      (?: 
                                                        (?:
                                                          (?:\s+$estNamePlus\s+is\s+in\s+$estNamePlus\n)*
                                                          (?:$estNamePlus\n)* 
                                                        )* 
                                                      )  
                                                    )
                                                /xg 
  	                         );  	                         
  	my $i=0; #create a val to show cicle number 
  	foreach my $eachContigSeqNameUnit (@contigsSeqNameArray){   		
  		print "\$i=$i\t\$eachContigSeqNameUnit=$eachContigSeqNameUnit\n\n"; 
  		if ($eachContigSeqNameUnit=~m/\*{10,20}\s+Contig\s+(\d+)\s+\*{10,20}\n/){
      	$contigNBhere=$1;
  		  my @estNames=(  $eachContigSeqNameUnit=~m/
  		                                             \s*($estNamePlus)(?:\s+is\s+in\s+$estNamePlus)?\n
  		                                           /xg
  		               
  		               );
  		  foreach my $eachEstNm (@estNames){  #print "\$eachEstNm=$eachEstNm\n";
  		  	#gi|106779132|gb|EC097509.1|EC097509+ is in gi|106779143|gb|EC097520.1|EC097520+
  		    if ($eachEstNm=~m/^(\S+)([\+-])(?:\s+is\s+in\s+\S+[\+-])?$/){   #gi|372776673|gb|HO370755.1|HO370755- is in gi|372832635|gb|HO404838.1|HO404838+
  		    	my $estNm=$1; my $direction=$2;    if ($direction eq '-'){print "CheckThisOut \$inCapOutFile=$inCapOutFile\n\n";}   #print "000000\$estNm=$estNm\t\$direction=$direction\n"; 
  		    	$shortedEstName=$estNm;       #print "111111\$shortedEstName=$shortedEstName\n";
  		    	$shortedEstName=~s/^\s+//g; $shortedEstName=~s/\s+$//g;  #print "222222\$shortedEstName=$shortedEstName\n";
  		    	if ( (length $shortedEstName) > 20 ){$shortedEstName=substr ($shortedEstName, 0, 20); #print "333333\$shortedEstName=$shortedEstName\n";
  		    	} 
  		    	$shortedEstName=~s/\s+$//g;      #print "444444\$shortedEstName=$shortedEstName\n"; #cap3_align_name_handle_epressions
  		    	
  		    	$estCtgEstNmHash->{$contigNBhere}->{$shortedEstName}->{$direction}=$estNm;  #print "\$estCtgEstNmHash->{\$contigNBhere}->{\$shortedEstName}->{\$direction}=\$estNm=\$estCtgEstNmHash->{$contigNBhere}->{$shortedEstName}->{$direction}=$estNm\n";
  		    }
  		    else {
  		      die "The \$eachEstNm is not in correct format:\n\$eachEstNm=$eachEstNm \n\n";
  		    }
  		  }
  		  
  		}
  		else {die "The \$eachContigSeqNameUnit is not in correct format:\n\$eachContigSeqNameUnit=$eachContigSeqNameUnit \n\n";}
  		$i++;   	
  	} 
  
    my @alignArray=(   $contigAlignUnit=~m/  	                          
                                            (
                                              $Contig_Line
                                              (?:
                                                $markerLine
                                                (?:$eachAlign)+
                                                $onlyLINE
                                                $consensus
                                                \n
                                              )+
                                            )
                                          /xg 
                    
                    );
                    
    my $j=0; 
    foreach my $eachContigAlignUnit (@alignArray)    {                                print "\$j=$j\t\$eachContigAlignUnit=$eachContigAlignUnit\n\n"; 
      if ($eachContigAlignUnit=~m/\*{10,20}\s+Contig\s+(\d+)\s+\*{10,20}\n/){
      	my $ContigNumber=$1;
      	my @eachAlignPartArray=($eachContigAlignUnit=~m/      	                        
      	                                               (
                                                         $markerLine
                                                         (?:$eachAlign)+
                                                         $onlyLINE
                                                         $consensus
                                                         \n
                                                       )      	                                               
      	                                               /gx      	                        
      	                        );      	                        
      	my $alingHash;               # 建立hash 数据结构：$alingHash->{$LineHead}->{$k}=$LineSequence;   #这里的语句，将每行序列信息，放入hash中。hash的key分别为，该est在consensus中出现的顺序号，在第几个段中出现。值则是该具体段的序列
      	my $xunXuHash;               #这个hash用来装，以顺序为key，alignhead为值的 hash
        my $xunXuInAlingHash=-1;     #$xunXuHash的key
      	my $consenseHash;                        
      	                        
      	my $k=0;
      	foreach my $eachAlignPart(@eachAlignPartArray){         	  print "\$k=$k\t\$eachAlignPart=$eachAlignPart\n";
      	  if ($eachAlignPart=~m/      	                        
      	                         $markerLine
                                 ((?:$eachAlign)+)
                                 $onlyLINE
                                 ($consensus)
                                 \n      	                        
      	                       /x      	                         
      	      )
      	  {
      	    my $eachAlignPart=$1;      	    
      	    my $ConsensusPart=$2;
      	    
      	    my @eachAlignArray=($eachAlignPart=~m/
      	                                           ($eachAlign)
      	                                         /gx
      	                        );
      	    my $m=0;
      	    foreach my $eachAlingLine (@eachAlignArray) {     	    	print "Foreach 1st line print:\$m=$m\t\$eachAlingLine=$eachAlingLine\n";
      	    	if ($eachAlingLine=~m/((\S.{19,19})([\+-]))\s(\s*\S{1,60}\s*)\n/){
      	    		my $LineHead=$1;
      	    		my $LineHeadName=$2;
      	    		my $LineHeadDirect=$3;
      	    		my $LineSequence=$4;                        	    		print "\$LineHeadName=$LineHeadName\t\$LineHeadDirect=$LineHeadDirect\t\$LineSequence=$LineSequence\n";
      	    	  
      	    	  #以下二行语句，是为了获得$alingHash中每个key进入hash的时候的序号。
      	    	  if ( defined($alingHash->{$LineHead}) ){                                   }
      	    	  else {                                      
      	    	  	$xunXuInAlingHash++;           
      	    	  	$xunXuHash->{$xunXuInAlingHash}=$LineHead;      	    	  print "                \$xunXuHash->{\$xunXuInAlingHash}=\$xunXuHash->{$xunXuInAlingHash}=\$LineHead=$LineHead;\n";
      	    	  }
      	    	  
      	    	  $alingHash->{$LineHead}->{$k}=$LineSequence;   #这里的语句，将每行序列信息，放入hash中。hash的key分别为，该est在consensus中出现的顺序号，在第几个段中出现。值则是该具体段的序列
      	    	  print "\$alingHash->{\$LineHead}->{\$k}=\$alingHash->{$LineHead}->{$k}=\$LineSequence=$LineSequence\n";
      	    	  
      	    	}
      	    	else {
      	    	  die "Each Align Line is not OK!!\n\$eachAlingLine=$eachAlingLine\n";
      	    	}
      	      $m++;
      	    }
      	    
      	    if ($ConsensusPart=~m/consensus\s{13,13}(\s*\S{1,60}\s*)\n/){
      	    	my $ConsesusSequence=$1;
      	    	#print "\$k=$k\t\$ConsesusSequence=$ConsesusSequence\n\n";
      	    	$consenseHash->{$k}=$ConsesusSequence;
      	    	print "                         \$consenseHash->{\$k}=\$consenseHash->{$k}=\$ConsesusSequence=$ConsesusSequence\n";
      	    }
      	    else {
      	      die "The Consensus line is not OK!!:\n\$ConsensusPart=$ConsensusPart\n\n";
      	    }
      	  } 
      	  else {
      	    die "There is something wrong here:\n\$eachAlignPart=$eachAlignPart\n\n";
      	  }
      	  $k++;
      	}
      	
      	
      	my $finalAlignSeqHash;   #建立hash，用来装全部的 比对序列
      	my $empthSeqLine='';    	for (my $r=0; $r<60; $r++){$empthSeqLine.=' ';} #建立空行，由60个空格组成
      		
      	print "\$k=$k\n";
      	#下面循环，将所有比对序列拼接成单行，并对齐
      	for (my $n=0; $n<$k; $n++){ print "\$n=$n\t$k=$k\n";
      	  foreach my $lineHeadKey (    sort {$a cmp $b} (   keys (  %{ $alingHash }  )   )    ){ print "\$lineHeadKey=$lineHeadKey\n";
      	    if ( defined($alingHash->{$lineHeadKey}->{$n}) ){
      	      $finalAlignSeqHash->{$lineHeadKey}.=$alingHash->{$lineHeadKey}->{$n};
      	    }
      	    else {
      	      $finalAlignSeqHash->{$lineHeadKey}.=$empthSeqLine;
      	    }
      	  }
      	}
      	
      	#下面的语句将 contig序列拼成一行，并对齐
      	my $finalConsSeq='';
      	foreach my $ConsKKey (    sort {$a<=>$b} (   keys (  %{ $consenseHash }  )   )    ){
      	  $finalConsSeq.=$consenseHash->{$ConsKKey};
      	}
      	
      	#下面的语句，将各行比对 对齐序列 打印出来
      	foreach my $xunXuKey (    sort {$a<=>$b} (   keys (  %{ $xunXuHash }  )   )    ){  #这个循环中，每个循环对一个Est的数据进行处理
      	  my $tempHead=$xunXuHash->{$xunXuKey};
      	  my $tempOutSeq = $finalAlignSeqHash->{$tempHead};
      	  printf "%5s\t%20s\t%100s\n", $xunXuKey, $tempHead, $tempOutSeq; 
      	  if ($tempHead=~m/^(\S+)\s*([\+-])$/){
      	    my $hdNm=$1; my $hddrct=$2;
      	    if (defined ($estCtgEstNmHash->{$ContigNumber}->{$hdNm}->{$hddrct})){
      	    	my $fullName=$estCtgEstNmHash->{$ContigNumber}->{$hdNm}->{$hddrct};
      	    	
      	    	$returnedHash->{$ContigNumber}->{'1_OrderHash'}->{$xunXuKey}=[$fullName,$hddrct];
      	    	
      	    	##1##下面的语句用于 获得所有输出的数据结构，这里主要是用于拼接的est的序列
      	    	$returnedHash->{$ContigNumber}->{$AlignKw}->{$fullName}->{$hddrct}->{$JustSequnceKw}=$tempOutSeq;                                             #结果数据结构注释1：$returnedHash->{contig的序数}->{注明是比对或被拼接的原始est}->{原est全名}->{比对或拼接中的方向}->{注明这里是序列信息}=比对中对齐的序列
      	    	my $posIdx=0;  my $NposIdx=0; my $SposIdx=0; my $GposIdx=0; my $tOsLength=length $tempOutSeq; 
      	    	my ($alignNactgTypeSeq,  $alignSpaceTypeSeq,  $alignGangTypeSeq)=('','','');
      	    	while($posIdx < $tOsLength){                               #这个循环，是以包括空格横杠的序列的每个位置作为循环体的个体进行操作的，也就是$posIdx
      	    		my $eachWord=substr($tempOutSeq,$posIdx,1);              #这里 获得包括所有空格横杠的每个位置上的字符
      	    		
      	    		#下面获得每个位置上的相应信息
      	    		$returnedHash->{$ContigNumber}->{$AlignKw}->{$fullName}->{$hddrct}->{$WlPos2WordKw}->{$posIdx}=$eachWord;                                   #结果数据结构注释2：$returnedHash->{contig的序数}->{注明是比对或被拼接的原始est}->{原est全名}->{比对或拼接中的方向}->{注明是全长的对齐序列（包括空格横杆）的位置信息}->{全长位置}=该位置的字符
      	    		if    ($eachWord eq ' '){                                                                                                       
      	    			$returnedHash->{$ContigNumber}->{$AlignKw}->{$fullName}->{$hddrct}->{$Wpos2DifTypePosrelationKw}->{$posIdx}=[$SpaceTypeKw,$SposIdx];      #结果数据结构注释3：$returnedHash->{contig的序数}->{注明是比对或被拼接的原始est}->{原est全名}->{比对或拼接中的方向}->{注明是由全长位置直接指向 字符的类型和位置信息}->{全长位置}=【该位置的字符类型，该位置字符在所属类型字符串中的位置】
      	    			$returnedHash->{$ContigNumber}->{$AlignKw}->{$fullName}->{$hddrct}->{$DifTypePos2WposrelationKw}->{$SpaceTypeKw}->{$SposIdx}=$posIdx;     #结果数据结构注释4：$returnedHash->{contig的序数}->{注明是比对或被拼接的原始est}->{原est全名}->{比对或拼接中的方向}->{注明是由不同类型字符及位置直接指向 全长位置信息}->{该位置的字符类型}->{该位置字符在所属类型字符串中的位置}=全长位置
      	    			$returnedHash->{$ContigNumber}->{$AlignKw}->{$fullName}->{$hddrct}->{$DifTypePos2WordKw}->{$SpaceTypeKw}->{$SposIdx}=$eachWord;           #结果数据结构注释5：$returnedHash->{contig的序数}->{注明是比对或被拼接的原始est}->{原est全名}->{比对或拼接中的方向}->{注明是由不同类型字符及位置直接指向 该位置字符}->{该位置的字符类型}->{该位置字符在所属类型字符串中的位置}=该位置的字符
      	    		  $alignSpaceTypeSeq.=$eachWord;
      	    		  $SposIdx++;
      	    		}
      	    		elsif ($eachWord eq '-'){   
      	    			$returnedHash->{$ContigNumber}->{$AlignKw}->{$fullName}->{$hddrct}->{$Wpos2DifTypePosrelationKw}->{$posIdx}=[$GangTypeKw,$GposIdx];       #同结果数据结构注释3
      	    			$returnedHash->{$ContigNumber}->{$AlignKw}->{$fullName}->{$hddrct}->{$DifTypePos2WposrelationKw}->{$GangTypeKw}->{$GposIdx}=$posIdx;      #同结果数据结构注释4
      	    			$returnedHash->{$ContigNumber}->{$AlignKw}->{$fullName}->{$hddrct}->{$DifTypePos2WordKw}->{$GangTypeKw}->{$GposIdx}=$eachWord;            #同结果数据结构注释5
      	    			$alignGangTypeSeq.=$eachWord;
      	    			$GposIdx++;
      	    		}
      	    		else {
      	    			$returnedHash->{$ContigNumber}->{$AlignKw}->{$fullName}->{$hddrct}->{$Wpos2DifTypePosrelationKw}->{$posIdx}=[$NactgTypeKw,$NposIdx];      #同结果数据结构注释3
      	    			$returnedHash->{$ContigNumber}->{$AlignKw}->{$fullName}->{$hddrct}->{$DifTypePos2WposrelationKw}->{$NactgTypeKw}->{$NposIdx}=$posIdx;     #同结果数据结构注释4
      	    			$returnedHash->{$ContigNumber}->{$AlignKw}->{$fullName}->{$hddrct}->{$DifTypePos2WordKw}->{$NactgTypeKw}->{$NposIdx}=$eachWord;           #同结果数据结构注释5
      	    			$alignNactgTypeSeq.=$eachWord;
      	    			$NposIdx++;
      	    		}
      	    		$posIdx++ ;
      	    		     	    		
      	    		$returnedHash->{$ContigNumber}->{$AlignKw}->{$fullName}->{$hddrct}->{$wholeLengthKw}    =$posIdx;                                           #结果数据结构注释6：$returnedHash->{contig的序数}->{原est全名}->{比对或拼接中的方向}->{注明是 对齐后包括空格横杆的序列的字符串的长度}  =对齐后包括空格横杆的序列的字符串的长度
      	    		$returnedHash->{$ContigNumber}->{$AlignKw}->{$fullName}->{$hddrct}->{$NactgTypeLengthKw}=$NposIdx;                                          #结果数据结构注释7：$returnedHash->{contig的序数}->{原est全名}->{比对或拼接中的方向}->{注明是 对齐后包括不空格横杆的序列的字符串的长度}=对齐后包括不空格横杆的序列的字符串的长度
      	    		$returnedHash->{$ContigNumber}->{$AlignKw}->{$fullName}->{$hddrct}->{$SpaceTypeLengthKw}=$SposIdx;                                          #结果数据结构注释8：$returnedHash->{contig的序数}->{原est全名}->{比对或拼接中的方向}->{注明是 对齐后只包括空格的字符串的长度}          =对齐后只包括空格的字符串的长度
      	    		$returnedHash->{$ContigNumber}->{$AlignKw}->{$fullName}->{$hddrct}->{$GangTypeLengthKw} =$GposIdx;                                          #结果数据结构注释9：$returnedHash->{contig的序数}->{原est全名}->{比对或拼接中的方向}->{注明是 对齐后只包括横杆的字符串的长度}          =对齐后只包括横杆的字符串的长度
      	    		$returnedHash->{$ContigNumber}->{$AlignKw}->{$fullName}->{$hddrct}->{$NactgTypeSeqKw}   =$alignNactgTypeSeq;                               #结果数据结构注释10：$returnedHash->{contig的序数}->{原est全名}->{比对或拼接中的方向}->{注明是 对齐后包括不空格横杆的序列的字符串}      =对齐后包括不空格横杆的序列的字符串
      	    		$returnedHash->{$ContigNumber}->{$AlignKw}->{$fullName}->{$hddrct}->{$SpaceTypeSeqKw}   =$alignGangTypeSeq;                                #结果数据结构注释11：$returnedHash->{contig的序数}->{原est全名}->{比对或拼接中的方向}->{注明是 对齐后只包括空格的字符串}                =对齐后只包括空格的字符串
      	    		$returnedHash->{$ContigNumber}->{$AlignKw}->{$fullName}->{$hddrct}->{$GangTypeSeqKw}    =$alignSpaceTypeSeq;                               #结果数据结构注释12：$returnedHash->{contig的序数}->{原est全名}->{比对或拼接中的方向}->{注明是 对齐后只包括横杆的字符串}                =对齐后只包括横杆的字符串
      	    		      	    		
      	    	}
      	    	
      	    	
      	    	my $estLength=(  keys ( $returnedHash->{$ContigNumber}->{$AlignKw}->{$fullName}->{$hddrct}->{$DifTypePos2WposrelationKw}->{$NactgTypeKw} )  );
      	    	if (  defined ( $returnedHash->{$ContigNumber}->{$AlignKw}->{$fullName}->{$hddrct}->{$DifTypePos2WposrelationKw}->{$NactgTypeKw}->{$estLength-1} )  ){}else {die"\n\n\nDIE!!!!!$warnMsgHead\$estLength=$estLength\n";}
      	    	$returnedHash->{$ContigNumber}->{$AlignKw}->{$fullName}->{$hddrct}->{'whoLength_0stt'}=$returnedHash->{$ContigNumber}->{$AlignKw}->{$fullName}->{$hddrct}->{$DifTypePos2WposrelationKw}->{$NactgTypeKw}->{0};
      	    	$returnedHash->{$ContigNumber}->{$AlignKw}->{$fullName}->{$hddrct}->{'whoLength_1end'}=$returnedHash->{$ContigNumber}->{$AlignKw}->{$fullName}->{$hddrct}->{$DifTypePos2WposrelationKw}->{$NactgTypeKw}->{$estLength-1};
      	    	$returnedHash->{$ContigNumber}->{$AlignKw}->{$fullName}->{$hddrct}->{'whoLength_2len'}=$estLength;		
      	    	##1##和上面的成对
      	    	
      	    	
      	    }
      	    else{ print_all_sub_array ($estCtgEstNmHash); die"cap3_align_name_handle_epressions is not crrect, correct it!!\n\$estCtgEstNmHash->{\$ContigNumber}->{\$hdNm}->{\$hddrct}=\$estCtgEstNmHash->{$ContigNumber}->{$hdNm}->{$hddrct}\n"; }
      	  }
      	  else {die "The \$tempHead is not in correct format:\$tempHead=$tempHead\n\n";}
      	  
      	#$estCtgEstNmHash->{$contigNBhere}->{$shortedEstName}->{$direction}=$estNm;
      	} 
      	printf "%5s\t%20s\t%100s\n", 'C','Consensus', $finalConsSeq;
  
  #  ( $AlignKw,$ConsensusKw,$JustSequnceKw, $Wpos2DifTypePosrelationKw, $DifTypePos2WposrelationKw, $WlPos2WordKw, $DifTypePos2WordKw, $NactgTypeKw, $SpaceTypeKw, $GangTypeKw, $wholeLengthKw,  $NactgTypeLengthKw,  $SpaceTypeLengthKw,  $GangTypeLengthKw,    $NactgTypeSeqKw,  $SpaceTypeSeqKw,  $GangTypeSeqKw   )
  #  ['_Align', '_Consensus','JustSequnceKw','_Wpos2DifTypePos',         '_DifTypePos2Wpos',        '_WlPos2Word', '_DifTypePos2Word', '_NactgType', '_SpaceType', '_GangType', '_wholeLengthKw','_NactgTypeLengthKw','_SpaceTypeLengthKw','_GangTypeLengthKw','_NactgTypeSeqKw','_SpaceTypeSeqKw','_GangTypeSeqKw' ];
  #    0        1            2               3                            4                          5               6                  7             8             9            10               11                   12                   13                 14                15                 16            
  #  拼接Est， 拼接Contig，   只是序列的标示     全长和不同类型位置的转换    不同类型和全长位置的转换   全长位置字符    不同类型位置字符    原序列        空格            横杆        全长长度         原序列长度          空格总长度             横杠总长度      原序列           空格序列              横杠序列    

      	
      	##2##下面的语句用于 获得所有输出的数据结构，这里主要是用于 拼接后的contig
     	  $returnedHash->{$ContigNumber}->{$ConsensusKw}->{$JustSequnceKw}=$finalConsSeq;     	                                                                   #类似 结果数据结构注释1：$returnedHash->{contig的序数}->{注明是拼接完成后的Contig}->{注明这里是序列信息}=比对中对齐的contig序列
        my $CposIdx=0;  my $CNposIdx=0; my $CSposIdx=0; my $CGposIdx=0; my $fCsLength=length $finalConsSeq; 
      	my ($consNactgTypeSeq,  $consSpaceTypeSeq,  $consGangTypeSeq)=('','','');
      	while($CposIdx<$fCsLength){
      		my $CeachWord=substr($finalConsSeq,$CposIdx,1);              #($AlignKw,$ConsensusKw, $Wpos2DifTypePosrelationKw, $DifTypePos2WposrelationKw, $WlPos2WordKw, $DifTypePos2WordKw, $NactgTypeKw, $SpaceTypeKw, $GangTypeKw, $wholeLengthKw,  $NactgTypeLengthKw,  $SpaceTypeLengthKw,  $GangTypeLengthKw  ) 
      		
      		$returnedHash->{$ContigNumber}->{$ConsensusKw}->{$WlPos2WordKw}->{$CposIdx}=$CeachWord;                                              #类似 结果数据结构注释2：$returnedHash->{contig的序数}->{原est全名}->{注明是拼接完成后的Contig}->{注明是全长的对齐序列（包括空格横杆）的位置信息}->{全长位置}=该位置的字符
      		if    ($CeachWord eq ' '){
      			$returnedHash->{$ContigNumber}->{$ConsensusKw}->{$Wpos2DifTypePosrelationKw}->{$CposIdx}=[$SpaceTypeKw,$CSposIdx];                 #类似 结果数据结构注释3
      			$returnedHash->{$ContigNumber}->{$ConsensusKw}->{$DifTypePos2WposrelationKw}->{$SpaceTypeKw}->{$CSposIdx}=$CposIdx;                #类似 结果数据结构注释4
      			$returnedHash->{$ContigNumber}->{$ConsensusKw}->{$DifTypePos2WordKw}->{$SpaceTypeKw}->{$CSposIdx}=$CeachWord;                      #类似 结果数据结构注释5
      		  $consSpaceTypeSeq.=$CeachWord;      
      		  $CSposIdx++;
      		}
      		elsif ($CeachWord eq '-'){
      			$returnedHash->{$ContigNumber}->{$ConsensusKw}->{$Wpos2DifTypePosrelationKw}->{$CposIdx}=[$GangTypeKw,$CGposIdx];                  #类似 同结果数据结构注释3
      			$returnedHash->{$ContigNumber}->{$ConsensusKw}->{$DifTypePos2WposrelationKw}->{$GangTypeKw}->{$CGposIdx}=$CposIdx;                 #类似 同结果数据结构注释4
      			$returnedHash->{$ContigNumber}->{$ConsensusKw}->{$DifTypePos2WordKw}->{$GangTypeKw}->{$CGposIdx}=$CeachWord;                       #类似 同结果数据结构注释5
      			$consGangTypeSeq.=$CeachWord;
      			$CGposIdx++;
      		}
      		else {
      			$returnedHash->{$ContigNumber}->{$ConsensusKw}->{$Wpos2DifTypePosrelationKw}->{$CposIdx}=[$NactgTypeKw,$CNposIdx];                 #类似 同结果数据结构注释3
      			$returnedHash->{$ContigNumber}->{$ConsensusKw}->{$DifTypePos2WposrelationKw}->{$NactgTypeKw}->{$CNposIdx}=$CposIdx;                #类似 同结果数据结构注释4
      			$returnedHash->{$ContigNumber}->{$ConsensusKw}->{$DifTypePos2WordKw}->{$NactgTypeKw}->{$CNposIdx}=$CeachWord;                      #类似 同结果数据结构注释5
      			$consNactgTypeSeq.=$CeachWord;
      			$CNposIdx++;
      		}
      		$CposIdx++;		
      		
 	    		$returnedHash->{$ContigNumber}->{$ConsensusKw}->{$wholeLengthKw}    =$CposIdx;                                                       #类似结果数据结构注释6：$returnedHash->{contig的序数}->{注明是拼接完成后的Contig}->{注明是 对齐后包括空格横杆的序列的字符串的长度}  =对齐后包括空格横杆的序列的字符串的长度  
 	    		$returnedHash->{$ContigNumber}->{$ConsensusKw}->{$NactgTypeLengthKw}=$CNposIdx;                                                      #类似结果数据结构注释7：$returnedHash->{contig的序数}->{注明是拼接完成后的Contig}->{注明是 对齐后包括不空格横杆的序列的字符串的长度}=对齐后包括不空格横杆的序列的字符串的长度
 	    		$returnedHash->{$ContigNumber}->{$ConsensusKw}->{$SpaceTypeLengthKw}=$CSposIdx;                                                      #类似结果数据结构注释8：$returnedHash->{contig的序数}->{注明是拼接完成后的Contig}->{注明是 对齐后只包括空格的字符串的长度}          =对齐后只包括空格的字符串的长度          
 	    		$returnedHash->{$ContigNumber}->{$ConsensusKw}->{$GangTypeLengthKw} =$CGposIdx;                                                      #类似结果数据结构注释9：$returnedHash->{contig的序数}->{注明是拼接完成后的Contig}->{注明是 对齐后只包括横杆的字符串的长度}          =对齐后只包括横杆的字符串的长度          
          $returnedHash->{$ContigNumber}->{$ConsensusKw}->{$NactgTypeSeqKw}   =$consNactgTypeSeq;                                                 #结果数据结构注释10：$returnedHash->{contig的序数}->{注明是拼接完成后的Contig}->{注明是 对齐后包括不空格横杆的序列的字符串}      =对齐后包括不空格横杆的序列的字符串      
      	  $returnedHash->{$ContigNumber}->{$ConsensusKw}->{$SpaceTypeSeqKw}   =$consGangTypeSeq;                                                  #结果数据结构注释11：$returnedHash->{contig的序数}->{注明是拼接完成后的Contig}->{注明是 对齐后只包括空格的字符串}                =对齐后只包括空格的字符串                
      	  $returnedHash->{$ContigNumber}->{$ConsensusKw}->{$GangTypeSeqKw}    =$consSpaceTypeSeq;                                                 #结果数据结构注释12：$returnedHash->{contig的序数}->{注明是拼接完成后的Contig}->{注明是 对齐后只包括横杆的字符串}                =对齐后只包括横杆的字符串                

      		
      		
      	}
     	  ##2##和上面的成对
     	  
     	  my $ctgLength=(  keys ( $returnedHash->{$ContigNumber}->{$ConsensusKw}->{$DifTypePos2WposrelationKw}->{$NactgTypeKw} )  );
      	if (  defined ( $returnedHash->{$ContigNumber}->{$ConsensusKw}->{$DifTypePos2WposrelationKw}->{$NactgTypeKw}->{$ctgLength-1} )  ){}else {die"\n\n\nDIE!!!!!$warnMsgHead\$ctgLength=$ctgLength\n";}
      	$returnedHash->{$ContigNumber}->{$ConsensusKw}->{'whoLength_0stt'}=$returnedHash->{$ContigNumber}->{$ConsensusKw}->{$DifTypePos2WposrelationKw}->{$NactgTypeKw}->{0};
      	$returnedHash->{$ContigNumber}->{$ConsensusKw}->{'whoLength_1end'}=$returnedHash->{$ContigNumber}->{$ConsensusKw}->{$DifTypePos2WposrelationKw}->{$NactgTypeKw}->{$ctgLength-1};
      	$returnedHash->{$ContigNumber}->{$ConsensusKw}->{'whoLength_2len'}=$ctgLength;
      	    	
      }
      else {
        die "The Contig head is wrong!! :\n\$eachContigAlignUnit=$eachContigAlignUnit\n\n";
      }
      $j++;
    } 
  }
  else {
    die "\$inCapOutFile=$inCapOutFile\n\$wholeIncapout=$wholeIncapout\n";
  }
  
  
  if (  ( ref ($returnedHash) ) eq 'HASH' ){
    foreach my $keyLev_0 (    sort { $a cmp $b } (   keys (  %{ $returnedHash } )   )    ){   #print "    \$keyLev_0=$keyLev_0\n";
    	my $valLev_0=$returnedHash->{$keyLev_0}->{'1_OrderHash'};                               #print "    \$valLev_0=$valLev_0\n";
    	my $refLev_0=ref ($valLev_0);                                                           #print "    \$refLev_0=$refLev_0\n\n";  
    	
      if ( $refLev_0 eq 'HASH' ){
        foreach my $keyLev_1 (    sort { $a <=> $b } (   keys (  %{ $valLev_0 } )   )    ){   #print "      \$keyLev_1=$keyLev_1\n";
        	                                          
        	my $estName=$returnedHash->{$keyLev_0}->{'1_OrderHash'}->{$keyLev_1}->[0];
        	my $estStnd=$returnedHash->{$keyLev_0}->{'1_OrderHash'}->{$keyLev_1}->[1];
        	if (   (  defined ( $returnedHash->{$keyLev_0}->{'1_Est0lsNoD'} )  ) && ( $returnedHash->{$keyLev_0}->{'1_Est0lsNoD'}=~m/\S+/ )   ){
        		$returnedHash->{$keyLev_0}->{'1_Est0lsNoD'}.=",".$estName;
        	}
        	else {
        		$returnedHash->{$keyLev_0}->{'1_Est0lsNoD'}.=$estName;
        	}    
        	
        	if (   (  defined ( $returnedHash->{$keyLev_0}->{'1_Est1lsAdD'} )  ) && ( $returnedHash->{$keyLev_0}->{'1_Est1lsAdD'}=~m/\S+/ )   ){
        		$returnedHash->{$keyLev_0}->{'1_Est1lsAdD'}.=",".$estName."_".$estStnd;
        	}
        	else {
        		$returnedHash->{$keyLev_0}->{'1_Est1lsAdD'}.=$estName."_".$estStnd;
        	} 
        
        }
      }
      
    }
  }                                                        
          
  return $returnedHash;
  
  
}

  #  ( $AlignKw,$ConsensusKw,$JustSequnceKw, $Wpos2DifTypePosrelationKw, $DifTypePos2WposrelationKw, $WlPos2WordKw, $DifTypePos2WordKw, $NactgTypeKw, $SpaceTypeKw, $GangTypeKw, $wholeLengthKw,  $NactgTypeLengthKw,  $SpaceTypeLengthKw,  $GangTypeLengthKw,    $NactgTypeSeqKw,  $SpaceTypeSeqKw,  $GangTypeSeqKw   )
  #  ['_Align', '_Consensus','JustSequnceKw','_Wpos2DifTypePos',         '_DifTypePos2Wpos',        '_WlPos2Word', '_DifTypePos2Word', '_NactgType', '_SpaceType', '_GangType', '_wholeLengthKw','_NactgTypeLengthKw','_SpaceTypeLengthKw','_GangTypeLengthKw','_NactgTypeSeqKw','_SpaceTypeSeqKw','_GangTypeSeqKw' ];
  #->{$AlignKw}->{$fullName}->{$hddrct}->{$NactgTypeLengthKw}

sub simplify_cap3praser_Out_Hash{
	my ($inHash)=@_;
	my $outHash=Storable::dclone ($inHash);
	my $ptString;
	if (  ref( $outHash ) eq 'HASH'  ){
		foreach my $ContigNumber (    sort {$a<=>$b} (   keys (  %{ $outHash }  )   )    ){   
			 #_DifTypePos2Wpos
			 if (  ref( $outHash->{$ContigNumber}->{'_Consensus'}->{'_DifTypePos2Word'} ) eq 'HASH'  ){
         delete ( $outHash->{$ContigNumber}->{'_Consensus'}->{'_DifTypePos2Word'} );
       }
       if (  ref( $outHash->{$ContigNumber}->{'_Consensus'}->{'_DifTypePos2Wpos'} ) eq 'HASH'  ){
         delete ( $outHash->{$ContigNumber}->{'_Consensus'}->{'_DifTypePos2Wpos'} );
       }
       if (  ref( $outHash->{$ContigNumber}->{'_Consensus'}->{'_WlPos2Word'} ) eq 'HASH'  ){
         #delete ( $outHash->{$ContigNumber}->{'_Consensus'}->{'_WlPos2Word'} );
         my $wlpos2wordString=PrintSubArrayHash::change_hash_into_a_printStyle_string($outHash->{$ContigNumber}->{'_Consensus'}->{'_WlPos2Word'});
         $outHash->{$ContigNumber}->{'_Consensus'}->{'_WlPos2Word'}=$wlpos2wordString;
         $ptString.="                          Consensus + $wlpos2wordString\n";
         
       }
       if (  ref( $outHash->{$ContigNumber}->{'_Consensus'}->{'_Wpos2DifTypePos'} ) eq 'HASH'  ){
         delete ( $outHash->{$ContigNumber}->{'_Consensus'}->{'_Wpos2DifTypePos'} );
       }	                                                                          

			
		  if (   (  ref( $outHash->{$ContigNumber}->{'_Align'} ) eq 'HASH'  ) && (  ref( $outHash->{$ContigNumber}->{'1_OrderHash'} ) eq 'HASH'  )   ) {
		  	foreach my $ordKey (    sort {$a <=> $b} (   keys (  %{ $outHash->{$ContigNumber}->{'1_OrderHash'} }  )   )    ){
		  	  my ($EstID, $EstZoF)=@{ $outHash->{$ContigNumber}->{'1_OrderHash'}->{$ordKey} };
		  	  
		  	  if (  ref( $outHash->{$ContigNumber}->{'_Align'}->{$EstID}->{$EstZoF}->{'_DifTypePos2Word'} ) eq 'HASH'  ){
		        delete ( $outHash->{$ContigNumber}->{'_Align'}->{$EstID}->{$EstZoF}->{'_DifTypePos2Word'}  );
		      }	
		      if (  ref( $outHash->{$ContigNumber}->{'_Align'}->{$EstID}->{$EstZoF}->{'_DifTypePos2Wpos'} ) eq 'HASH'  ){
		        delete ( $outHash->{$ContigNumber}->{'_Align'}->{$EstID}->{$EstZoF}->{'_DifTypePos2Wpos'}  );
		      }	
		      if (  ref( $outHash->{$ContigNumber}->{'_Align'}->{$EstID}->{$EstZoF}->{'_WlPos2Word'} ) eq 'HASH'  ){
		      	my $wlpos2wordString=PrintSubArrayHash::change_hash_into_a_printStyle_string($outHash->{$ContigNumber}->{'_Align'}->{$EstID}->{$EstZoF}->{'_WlPos2Word'});
		      	$outHash->{$ContigNumber}->{'_Align'}->{$EstID}->{$EstZoF}->{'_WlPos2Word'}=$wlpos2wordString;
		      	$ptString.="$EstID $EstZoF $wlpos2wordString\n";
		        #delete ( $outHash->{$ContigNumber}->{'_Align'}->{$EstID}->{$EstZoF}->{'_WlPos2Word'}  );
		      }	
		      if (  ref( $outHash->{$ContigNumber}->{'_Align'}->{$EstID}->{$EstZoF}->{'_Wpos2DifTypePos'} ) eq 'HASH'  ){
		        delete ( $outHash->{$ContigNumber}->{'_Align'}->{$EstID}->{$EstZoF}->{'_Wpos2DifTypePos'}  );
		      }	
		  	  
		  	  
		  	}
		  		
		    #foreach my $EstID (    sort {$a cmp $b} (   keys (  %{ $outHash->{$ContigNumber}->{'_Align'} }  )   )    ){
		    #  if (  ref( $outHash->{$ContigNumber}->{'_Align'}->{$EstID} ) eq 'HASH'  ){
		    #    foreach my $EstZoF (    sort {$a cmp $b} (   keys (  %{ $outHash->{$ContigNumber}->{'_Align'}->{$EstID} }  )   )    ){
		    #      if (  ref( $outHash->{$ContigNumber}->{'_Align'}->{$EstID}->{$EstZoF}->{'_DifTypePos2Word'} ) eq 'HASH'  ){
		    #        delete ( $outHash->{$ContigNumber}->{'_Align'}->{$EstID}->{$EstZoF}->{'_DifTypePos2Word'}  );
		    #      }	
		    #      if (  ref( $outHash->{$ContigNumber}->{'_Align'}->{$EstID}->{$EstZoF}->{'_DifTypePos2Wpos'} ) eq 'HASH'  ){
		    #        delete ( $outHash->{$ContigNumber}->{'_Align'}->{$EstID}->{$EstZoF}->{'_DifTypePos2Wpos'}  );
		    #      }	
		    #      if (  ref( $outHash->{$ContigNumber}->{'_Align'}->{$EstID}->{$EstZoF}->{'_WlPos2Word'} ) eq 'HASH'  ){
		    #      	my $wlpos2wordString=PrintSubArrayHash::change_hash_into_a_printStyle_string($outHash->{$ContigNumber}->{'_Align'}->{$EstID}->{$EstZoF}->{'_WlPos2Word'});
		    #      	$outHash->{$ContigNumber}->{'_Align'}->{$EstID}->{$EstZoF}->{'_WlPos2Word'}=$wlpos2wordString;
		    #      	$ptString.="$EstID $EstZoF $wlpos2wordString\n";
		    #        #delete ( $outHash->{$ContigNumber}->{'_Align'}->{$EstID}->{$EstZoF}->{'_WlPos2Word'}  );
		    #      }	
		    #      if (  ref( $outHash->{$ContigNumber}->{'_Align'}->{$EstID}->{$EstZoF}->{'_Wpos2DifTypePos'} ) eq 'HASH'  ){
		    #        delete ( $outHash->{$ContigNumber}->{'_Align'}->{$EstID}->{$EstZoF}->{'_Wpos2DifTypePos'}  );
		    #      }	
		    #    } 
		    #  }    
			  #}
			}
			
			$outHash->{$ContigNumber}->{'0_WatchPart'}="\n$ptString";
			
		}
	}
	
	return $outHash;
}



1;

##########################################################################################################################################
# 
