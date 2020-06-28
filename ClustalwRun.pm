#!/usr/bin/perl -w
BEGIN{  push (@INC, '/home/fredjiang/fredProgs/PMs')  };


use strict;
use warnings;



############################################################################################################
#
#
#          sub clustalwForLongPathFile{  #因为clusatlw在运行的时候，如输入文件的路径过长，就会出错，所以编辑了这个函数
#          sub clustalwForLongPathFileOutOrgOrder{  #因为clusatlw在运行的时候，如输入文件的路径过长，就会出错，所以编辑了这个函数,和clustalwForLongPathFile的区别在于，输出的程序是以输入的Fasta的序列顺序来排列的
#          sub readInSeq{ #读取一个多序列的文件，将该文件变成hash，hash给出了序列数，对每个序列给出了以该序列primary_id为Key的各种信息
#          sub readInSeqAndUseClustalwOrder{  #读取一个多序列的文件，将该文件变成hash，hash给出了序列数，对每个序列给出了以该序列primary_id为Key的各种信息,这个程序增加了一个功能，利用格式是$MsfFormat的alignemnt文件$orgMsfFile中的顺序来对多序列文件进行排序
#          sub runingClustalw{  #直接运行Clustalw
#          sub runingClustalwInputOrder{  #直接运行Clustalw,但比对后 输出顺序仍然按要输入的 fasta序列的顺序输出
#          sub GetClustalwOrderHash{   #输入一个algnment文件，获得一个Hash，给出了序列ID和顺序之间的 互相对应关系的hash
#          sub ReWriteClustalw{  #对msf格式的clusatlw结果文件$orgMsfFile，进行格式转换，转换为$outFormat格式的$outMutiAlgnFile文件
#          sub RewriteClustalw20180408{     #某些输入，需要依赖 2018.04.08.genomePostionAnalysis.pl 中的中间数据结构 #输入一个 多序列比对的结果，变成有编号的新结果。同时输出mega结果，并返回一个hash，其中记录了 id和编号之间的关系
#          sub sortFastaFile{ #对fastafile进行排序，使用场景比较特异性
#          sub PrintPngForClustalw{  #利用clustalw的输出msf文件$clustalwOUTfile,和interproscan模块的结果文件$Whole_htmlParseOut,绘制图片，输出到$outPng中
#          sub MakeTreeFrom_DadSonHashs {          #利用父子关系建立树形结构           #($DS_AllDadSonHash, $DS_DomDadSonRelationHash)
#          sub MakeTreeFrom_DadSonHashs_onlyMultiLevesTree {    #利用父子关系建立树形结构                 #($DS_AllDadSonHash, $DS_DomDadSonRelationHash)
#          sub DiGuiDadSon{  利用递归法 建立树
#          sub PrintTreePoints{  #打印出树形结构
#          sub treePanelWork{  #画出树的png图
#          sub MergeImage{  #将png文件$image_1和$image_2,合并成$image_out, 应该是上下合并
#          sub AddTracks2Panel{  #将各个含有feature的Track加到Panel中
#          sub AddFeature2PanelAndTrack{ #将各个含有feature的Track加到Panel中
#          sub GetSeqPosList{  #利用clustalw的结果 及column_from_residue_number 等方法，获得经过比对对齐后，重新排布的各个片段位置形成的数列
#          sub GetSeqPosList_2 {  #比GetSeqPosList多了个Type的信息输入 #利用clustalw的结果 及column_from_residue_number 等方法，获得经过比对对齐后，重新排布的各个片段位置形成的数列
#          sub DrawClustalw{  #画出clustlaw结果的 png图 #这个函数在这个程序中没有被使用到
#          sub PrintList{  #这个函数在这个程序中没有被使用到 #获得用于打印出来的 外显子位点的 数列如： 116..156,158..181,184..193,195..202,206..263,265..313,315..330
#          sub NEWAddFeatureTrack2Panel{  #给Panel加Track
#         
############################################################################################################

use GD;
use Cwd;
use File::Basename;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::SeqFeature::Lite;
use TimeWork;
use Interproscan;
use Bio::Graphics::Panel;
use DieWork;
use FastaFileHandle;
                      
package ClustalwRun;


my $muscle_fastStart_mubmer=400;

sub ChangeSeqString_to_fit_dimand{        #  ClustalwRun::ChangeSeqString_to_fit_dimand( $inSeqString );
	my ( $inSeqString )=@_;
	
	my $warnMsgBody="\nIn package  ClustalwRun,\tIn sub ChangeSeqString_to_fit_dimand,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $caller_inform=DirFileHandle::print_SubCallerInform;
  
  my $outString;
  if  (   (  defined  ( $inSeqString )  ) && (  $inSeqString=~m/\S+/  )   ){
  	$outString=$inSeqString;
  	$outString=~s/\s+//g;
  	return $outString;
  }
  else {
  	my $dieMsg=$die_MsgHead.$caller_inform."\n\n\n\$inSeqString=$inSeqString is not right!!!! : $!\n\n\n\n";
  	print $dieMsg;
  	die   $dieMsg;
  }

}


sub muscleRUN{        #  ClustalwRun::muscleRUN;
	my ($inFastaFile, $outPutFile, $outPutFmt, $fastaCutLmt )=@_;
	
	my $warnMsgBody="\nIn package  ClustalwRun,\tIn sub muscleRUN,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $caller_inform=DirFileHandle::print_SubCallerInform;
  
  if  (   (  defined ( $fastaCutLmt )  ) && ( $fastaCutLmt=~m/\d+/ ) && ( $fastaCutLmt>0 )   ){
  	$fastaCutLmt=$fastaCutLmt;
  }
  else {
  	$fastaCutLmt=$muscle_fastStart_mubmer;
  }
  
  my $realOutFmt='';
  if (  defined  ( $outPutFmt )  ){
  	$outPutFmt=lc $outPutFmt; $outPutFmt=~s/^\s+//; $outPutFmt=~s/\s+$//;
  	if    ( $outPutFmt eq 'aln'  ){ $realOutFmt='-clw';  }
  	elsif ( $outPutFmt eq 'msf'  ){ $realOutFmt='-msf';  }
  	elsif ( $outPutFmt eq 'html' ){ $realOutFmt='-html'; }
  }
  
  my $inSeqNumber=0;
  my $tempHASH=ClustalwRun::readInSeq($inFastaFile, 'fasta');
	if (   (  defined ( $tempHASH )  ) && (  ref ( $tempHASH ) eq 'HASH'  ) && (  defined ( $tempHASH->{'howManySeq'} )  ) && (  $tempHASH->{'howManySeq'}=~m/\d+/  )   ){
		$inSeqNumber=$tempHASH->{'howManySeq'};
	}	
  
  
  if (   (  defined  ( $inFastaFile )  ) && (  -e ( $inFastaFile )  ) && (  defined  ( $outPutFile )  ) && (  $outPutFile=~m/\S+/  )   ){
  	my $command_1_slow='muscle'." -in ".$inFastaFile." ".$realOutFmt." -out ".$outPutFile;
  	my $command_2_fast='muscle'." -in ".$inFastaFile." ".$realOutFmt." -out ".$outPutFile."  -maxiters 1 -diags -sv -distance1 kbit20_3";  #muscle -in 0_3_0.all.pep -out 0_3_1.msl.test.out -maxiters 1 -diags -sv -distance1 kbit20_3
  	
  	my $command=$command_1_slow;
  	if ($inSeqNumber>= $fastaCutLmt){
  		$command=$command_2_fast;
  	}
  	print $command;
  	warn  $command;
  	system ("$command");
  }
  else {
  	my $dieMsg=$die_MsgHead.$caller_inform."\$inFastaFile=$inFastaFile or \$outPutFile=$outPutFile is not right!!! :$!\n\n\n";
  	print $dieMsg;
  	die $dieMsg;
  }
  
}






sub cobaltRun{
  my ($inFastaSeqFile, $rpsDataBasePath)=@_;
  $rpsDataBasePath="/home/fredjiang/EightT/fredjiang.2018.04.02/2018.04.02.cobalt.rpsDatabase/cdd_clique_0.75";
  my $cobaltComand="cobalt -i $inFastaSeqFile -rpsdb $rpsDataBasePath -seqalign -outfmt mfasta -parse_deflines";
#cobalt -i SELENOM.txt -rpsdb /home/fredjiang/EightT/fredjiang.2018.04.02/2018.04.02.cobalt.rpsDatabase/cdd_clique_0.75 -seqalign SelM.msf -outfmt mfasta > selm.aln
}



sub clustalwForLongPathFile{  #因为clusatlw在运行的时候，如输入文件的路径过长，就会出错，所以编辑了这个函数
  my ($inPathFile, $outPathFile, $outFormat)=@_;
  my $inFile  =File::Basename::basename $inPathFile;
  my $inPath  =File::Basename::dirname ($inPathFile);   
  my $outPath =File::Basename::dirname $outPathFile;
  my $outFile =File::Basename::basename $outPathFile;
  my $InTempFl=TimeWork::GetTimeDirOrFileName()."in.txt";  
  my $otTempFl=TimeWork::GetTimeDirOrFileName()."ot.msf";
  my $otTdndFl=TimeWork::GetTimeDirOrFileName()."in.dnd";   #这个是Clustalw自动生成的一个文件 
  my $outDNDfile=$outPathFile; $outDNDfile=~s/(msf|aln)\s*$//; $outDNDfile.="dnd";
  
  my $command1="cp -f $inPathFile $InTempFl";
  my $command2="clustalw $InTempFl -OUTPUT=$outFormat -OUTFILE=$otTempFl";
  my $command3="rm -f $InTempFl";
  my $command4="mv -f $otTempFl $outPathFile";
  my $command5="mv -f $otTdndFl $outDNDfile";
  
  
  warn "\n\n\nNow in Sub clustalwForLongPathFile\n\n";
  warn "1 We are now in :";               system ("pwd");
  warn "\n\$command1=$command1\n";        system ("$command1");  
  warn "\n\$command2=$command2\n";        system ("$command2");
  warn "\n\$command3=$command3\n";        system ("$command3");
  warn "\n\$command4=$command4\n";        system ("$command4");
  warn "\n\$command5=$command5\n";        system ("$command5");
  warn "2 We are now in :";               system ("pwd");	
  warn "\nLast moment in sub clustalwForLongPathFile\n\n\n";	 	  	
  
}


sub clustalwForLongPathFileOutOrgOrder{  #因为clusatlw在运行的时候，如输入文件的路径过长，就会出错，所以编辑了这个函数,和clustalwForLongPathFile的区别在于，输出的程序是以输入的Fasta的序列顺序来排列的
  my ($inPathFile, $outPathFile, $outFormat)=@_;
  my $inFile  =File::Basename::basename $inPathFile;
  my $inPath  =File::Basename::dirname ($inPathFile);   
  my $outPath =File::Basename::dirname $outPathFile;
  my $outFile =File::Basename::basename $outPathFile;
  my $InTempFl=TimeWork::GetTimeDirOrFileName()."in.txt";  
  my $otTempFl=TimeWork::GetTimeDirOrFileName()."ot.msf";
  my $otTdndFl=TimeWork::GetTimeDirOrFileName()."in.dnd";   #这个是Clustalw自动生成的一个文件 
  my $outDNDfile=$outPathFile; $outDNDfile=~s/(msf|aln)\s*$//; $outDNDfile.="dnd";
  
  my $command1="cp -f $inPathFile $InTempFl";
  my $command2="clustalw $InTempFl -OUTPUT=$outFormat -OUTFILE=$otTempFl -OUTORDER=INPUT";
  my $command3="rm -f $InTempFl";
  my $command4="mv -f $otTempFl $outPathFile";
  my $command5="mv -f $otTdndFl $outDNDfile";
  
  
  warn "\n\n\nNow in Sub clustalwForLongPathFile\n\n";
  warn "1 We are now in :";               system ("pwd");
  warn "\n\$command1=$command1\n";        system ("$command1");  
  warn "\n\$command2=$command2\n";        system ("$command2");
  warn "\n\$command3=$command3\n";        system ("$command3");
  warn "\n\$command4=$command4\n";        system ("$command4");
  warn "\n\$command5=$command5\n";        system ("$command5");
  warn "2 We are now in :";               system ("pwd");	
  warn "\nLast moment in sub clustalwForLongPathFile\n\n\n";	 	  	
  
}

sub GetSeq_Number_in_aln_file{  #  my $Seq_Number=ClustalwRun::GetSeq_Number_in_aln_file( $Align_File, $Align_format ); 
	my ( $Align_File, $Align_format )=@_;
	
	my $warnMsgBody="\nIn package  ClustalwRun,\tIn sub GetSeq_Number_in_aln_file,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $caller_inform=DirFileHandle::print_SubCallerInform;
  
  my $Seq_Number=0; 
  if  (   (  defined ( $Align_File )  ) && (  -e ( $Align_File ) )   ){
  	
  	
	  my $in  = Bio::AlignIO->new(-file   => $Align_File ,      #
                                -format => $Align_format   );   #将输入的 aln文件，解析成 对象
  
    my $howManyAln=0;     
    while ( my $AlnObj = $in->next_aln() ) {  
    	$Seq_Number=$AlnObj->num_sequences;  	 
    	$howManyAln++;    
    }        
                     
    if ($howManyAln>1){
      my $dieMsg=$die_MsgHead.$caller_inform."\n\$Align_File=$Align_File\t\t\$howManyAln=$howManyAln\n\n\n";
  	   DieWork::Just_dieWork( $dieMsg );
      
    }                        
    
  } 
  return  $Seq_Number; 
}

sub readInSeq{   #  ClustalwRun::readInSeq($inseqFile, $inFormat);   #读取一个多序列的文件，将该文件变成hash，hash给出了序列数，对每个序列给出了以该序列primary_id为Key的各种信息
  my ($inseqFile, $inFormat)=@_;
  
  my $warnMsgBody="\nIn package  ClustalwRun,\tIn sub readInSeq,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $caller_inform=DirFileHandle::print_SubCallerInform;
  
  
  
  
  my $seqInObj=Bio::SeqIO->new(-file   => $inseqFile,    
                               -format => $inFormat );
  my $outHash;
  
  my $seqNb=0;                             
  while (my $seqObj=$seqInObj->next_seq){
    #print "primary_id: ", $seqObj->primary_id, "\n";
    #print "display_id: ", $seqObj->display_id, "\n";
    #print "accession_number: ", $seqObj->accession_number, "\n";
    #print "desc: ", $seqObj->desc, "\n\n";
    
                $outHash->{'id_seqHash'}->{ $seqObj->primary_id }->{'display_id'      }=$seqObj->display_id;
                $outHash->{'id_seqHash'}->{ $seqObj->primary_id }->{'accession_number'}=$seqObj->accession_number;
                $outHash->{'id_seqHash'}->{ $seqObj->primary_id }->{'desc'            }=$seqObj->desc;
                $outHash->{'id_seqHash'}->{ $seqObj->primary_id }->{'seqOrder'        }=$seqNb; 
    my $seqence=$outHash->{'id_seqHash'}->{ $seqObj->primary_id }->{'sequence'        }=$seqObj->seq();   
    my $U_posAR=$outHash->{'id_seqHash'}->{ $seqObj->primary_id }->{'sequence'        }=$seqObj->seq();   
    #my $U_pos_array=FastaFileHandle::Find_U_pos_inAminoAcid($inPEPstring);	if ( ref($U_pos_array) eq 'ARRAY' ){ $db_gencode=10; }
    
    $seqNb++;
  }
  
  $outHash->{'howManySeq'}=$seqNb;
  
  return $outHash;

}

sub NewReadInSeq_into_A_HASH{   #  ClustalwRun::NewReadInSeq_into_A_HASH($inseqFile, $inFormat);   #读取一个多序列的文件，将该文件变成hash，hash给出了序列数，对每个序列给出了以该序列primary_id为Key的各种信息
  my ($inseqFile, $inFormat)=@_;
  
  my $warnMsgBody="\nIn package  ClustalwRun,\tIn sub NewReadInSeq_into_A_HASH,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $caller_inform=DirFileHandle::print_SubCallerInform;
  
  
  
  warn "\$inseqFile=$inseqFile\n\$inFormat=$inFormat\n\n";
  my $seqInObj=Bio::SeqIO->new(-file   => $inseqFile,    
                               -format => $inFormat );
  my $outHash;
  
  my $seqNb=0;                             
  while (my $seqObj=$seqInObj->next_seq){
    
    
                $outHash->{'1_id_seqHash'}->{ $seqObj->primary_id }->{'0_display__id'}=$seqObj->display_id;
                $outHash->{'1_id_seqHash'}->{ $seqObj->primary_id }->{'1_accs_number'}=$seqObj->accession_number;
                $outHash->{'1_id_seqHash'}->{ $seqObj->primary_id }->{'2_description'}=$seqObj->desc;
                $outHash->{'1_id_seqHash'}->{ $seqObj->primary_id }->{'3_seque_Order'}=$seqNb; 
    my $seqence=$seqObj->seq(); 
                $outHash->{'1_id_seqHash'}->{ $seqObj->primary_id }->{'4____sequence'}=$seqence;  
    my $U_posAR=$outHash->{'1_id_seqHash'}->{ $seqObj->primary_id }->{'5_U_Pos_Array'}=FastaFileHandle::Find_U_pos_inAminoAcid($seqence);	   
    
    $seqence=~s/\s+//g;                  
    my $seqLentgh=length ( $seqence );
                $outHash->{'2_seqKeyHash'}->{ $seqence }=$seqLentgh;	
    
    
    $seqNb++;
  }
  
                $outHash->{'0_howManySeq'}=$seqNb;
                $outHash->{'0_howMayUcSq'}=(   keys (  %{ $outHash->{'2_seqKeyHash'} }  )   );
  
  #DirFileHandle::PrintAndWarnDumper ($outHash, "\n20190104-2\n");
  return $outHash;

}


# ClustalwRun::RunMuscle_and_obtain_varous_Format_alnFILE  ($in_Fasta_File, $muscle_out_fasta_File_name, $change_X_2_U_Muscle_MSF, $Del_number_befor_names_file, $Changed_name_file, $fastaCutLmt);
sub RunMuscle_and_obtain_varous_Format_alnFILE{
	my ($in_Fasta_File, $muscle_out_fasta_File_name, $change_X_2_U_Muscle_MSF, $Del_number_befor_names_file, $Changed_name_file, $fastaCutLmt)=@_;
	
  my $warnMsgBody="\nIn package  ClustalwRun,\tIn sub RunMuscle_and_obtain_varous_Format_alnFILE,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";  
  my $caller_inform=DirFileHandle::print_SubCallerInform;
  
  
  if   (   (  defined ( $in_Fasta_File )  ) && ( $in_Fasta_File=~m/\S+/ ) && (  ( -e $in_Fasta_File )  )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$in_Fasta_File=$in_Fasta_File should be a defined not empty file  !!  $!\n\n\n".$caller_inform ); 	}
	
	if   (   (  defined ( $muscle_out_fasta_File_name )  ) && ( $muscle_out_fasta_File_name=~m/\S+/ )   ){}
	else{		$muscle_out_fasta_File_name=$in_Fasta_File."\.01_mslfas.txt"; 	}
	
	if   (   (  defined ( $change_X_2_U_Muscle_MSF )  ) && ( $change_X_2_U_Muscle_MSF=~m/\S+/ )   ){}
	else{		$change_X_2_U_Muscle_MSF=$in_Fasta_File."\.02_mslX2U.msf"; 	}
	
	if   (   (  defined ( $Del_number_befor_names_file )  ) && ( $Del_number_befor_names_file=~m/\S+/ )   ){}
	else{		$Del_number_befor_names_file=$in_Fasta_File."\.03_NoNumb.msf"; 	}
	
  
  
	ClustalwRun::muscleRUN                                       ( $in_Fasta_File, $muscle_out_fasta_File_name,      '', $fastaCutLmt                          ) ;
	ClustalwRun::ChangeU_for_alignment_file                      ( $in_Fasta_File, $muscle_out_fasta_File_name, 'fasta', $change_X_2_U_Muscle_MSF, 'msf', 'X', 'U' );
  
  
  ClustalwRun::ChangeU_for_alignment_file_del_addedNB_of_muscle( $in_Fasta_File, $muscle_out_fasta_File_name, 'fasta', $Del_number_befor_names_file, 'msf', 'X', 'U' );
  
  #ClustalwRun::ChangeU_for_alignment_file_del_addedNB_of_muscle_and_change_name ( $MsfNewNmHSH, 
  #                                                               $in_Fasta_File, $muscle_out_fasta_File_name, 'fasta', $Changed_name_file, 'msf', 'X', 'U' );
	
}





sub ChangeU_for_alignment_file{ #ClustalwRun::ChangeU_for_alignment_file
	my ( $orgFastaSeqFile, $org_Align_File, $org_format, 
	                       $out_AlignFIle,  $out_format, 
	                       $U_char_in_org_Aln, $U_char_in_new_Aln, 
	                       $cutStt, $cutEnd)=@_;
	
	my $warnMsgBody="\nIn package  ClustalwRun,\tIn sub ChangeU_for_alignment_file,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $caller_inform=DirFileHandle::print_SubCallerInform;
  
  if  (   (  defined ( $orgFastaSeqFile )  ) && (  -e ( $orgFastaSeqFile )  ) && (  defined ( $org_Align_File )  ) && (  -e ( $org_Align_File ) )   ){
  	
  	
  	my $orgFastaSeq_HASH=ClustalwRun::NewReadInSeq_into_A_HASH( $orgFastaSeqFile , $org_format);  
  	###################################  column_from_residue_number
  	
  	my $in  = Bio::AlignIO->new(-file   => $org_Align_File ,      #
                                -format => $org_format   );   #将输入的 aln文件，解析成 对象
  
  	my $idRelationHash;
  
    my $howManyAln=0;                            
    while ( my $aln = $in->next_aln() ) {
    	
    	#临时文件 保存中间fas格式的文件
    	my $tempFile=TimeWork::GetTimeDirOrFileName().".fas";
    	
    	#建立临时文件的输入输出 对象
    	my $tempFast    = Bio::SeqIO->new(-file   => ">$tempFile" ,          #
                                        -format => 'fasta');               #临时文件的打开，用于写入fasta格式的文件, For Msf
    	
    	#循环读取 msf各条信息， 写入临时文件  	
    	my $seqIdx=1;
      foreach my $seq ($aln->each_seq) {
        #获取 id 和 msf的序列
        my $idHere=$seq->id(); my $seqHere=$seq->seq();                                   # 这里得到原来的id和比对的序列
        
        
        
        
        my $outIdx=$seqIdx."-".$idHere;
        #For Mega 将其中的U转换为C
        if ($out_format eq 'mega'){
        	my $changeHASH=&build_U_changeHash( $orgFastaSeq_HASH, $idHere, $U_char_in_org_Aln, $U_char_in_new_Aln, $seq);
        	$seqHere=FastaFileHandle::CHangeSpecificSequenceInAFastaSeq($seqHere, $changeHASH); warn "\n11111    $seqHere\n";
          $seqHere=~s/\./-/g;                    warn "\n22222    $seqHere\n";
        }  
        if ($out_format eq 'msf'){
        	my $changeHASH=&build_U_changeHash( $orgFastaSeq_HASH, $idHere, $U_char_in_org_Aln, $U_char_in_new_Aln, $seq);
        	$seqHere=FastaFileHandle::CHangeSpecificSequenceInAFastaSeq($seqHere, $changeHASH);                          
        } 
        
        
        $idRelationHash->{'Org2Chg'}->{$idHere}=$seqIdx;
        $idRelationHash->{'Chg2Org'}->{$seqIdx}=$idHere;
        
         
        #建立新的 SeqObj
        my $sortFasSeqObj=Bio::Seq->new( -display_id => $outIdx,   #$seqIdx,                          #序号 
                                                -seq => $seqHere                ); 
        #print "201811201803-1 \$cutStt=$cutStt\t$cutEnd=$cutEnd\n"; 
        if (  (  defined ( $cutStt )  ) && ( $cutStt=~m/\d+/) && ( $cutStt>0 ) && (  defined ( $cutEnd )  ) && ( $cutEnd=~m/\d+/) && ( $cutEnd>0 )   ){
        	my $truncateSeq=$sortFasSeqObj->trunc($cutStt, $cutEnd);  print "201811201803-2 \$truncateSeq=$truncateSeq\n";
        	#$sortFasSeqObj->seq()=$truncateSeq; 
        	$tempFast->write_seq($truncateSeq);   
        }
        
        else {
        	#写入
          $tempFast->write_seq($sortFasSeqObj);      
        }
                                                
                                                       
            
              
        $seqIdx++;
      }
      
      #建立读取临时文件的 IO对象
      my $fastain= Bio::AlignIO->new(-file   => $tempFile ,              ##  #For Msf #将fasta的临时文件作为 输入
                                     -format => 'fasta'             );   ##
      
      
      #建立 写入 新Msf，或新格式的文件的　ＩＯ对象
      
      my $out    = Bio::AlignIO->new(-file   => ">$out_AlignFIle" ,     #   #For Msf #建立新的输出对象
                                     -format => $out_format          );   #   
      
      ##For Msf #下面进行新输出 比对文件的 写入
      while   ( my $FstAln = $fastain->next_aln() ) {
      	$out->write_aln($FstAln);
      }                             
      system ("rm -f $tempFile");    ##删除临时文件
      
      
      $howManyAln++;
    }                         
    if ($howManyAln>1){
      my $dieMsg=$die_MsgHead.$caller_inform."\n\$org_Align_File=$org_Align_File\t\t\$howManyAln=$howManyAln\n\n\n";
  	  print $dieMsg;
  	  die $dieMsg;
      
    }                        
    return  $idRelationHash;
  	#####################################
  	
  	
  	
  	
  }
  else {
  	my $dieMsg=$die_MsgHead.$caller_inform."\$orgFastaSeqFile=$orgFastaSeqFile or \$org_Align_File=$org_Align_File is not right!!! :$!\n\n\n";
  	print $dieMsg;
  	die $dieMsg;
  }
  
  
}

sub Build_mostAA_clounm_MapHASH{   #   ClustalwRun::Build_mostAA_clounm_MapHASH( $In_aln___file, $In_aln_format );
	my ($In_aln___file, $In_aln_format)=@_;
	
	my $warnMsgBody="\nIn package  ClustalwRun,\tIn sub Build_mostAA_clounm_MapHASH,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	my $anyAAfoundToBeMax=0;
	my $BigestAA_Clm_Nb=-999999;
	my $outTempARRAY;
	my $allAAarray=FastaFileHandle::BuildAAarray();
	foreach my $eachAA (  @{ $allAAarray }  ){   print "\n20190201-0-0-0 \$eachAA=$eachAA \$In_aln___file=$In_aln___file\n";
		my $tempARRAY=ClustalwRun::Build_someAA_clounm_map_HASH( $In_aln___file, $In_aln_format, $eachAA );  print "\n20190201-0-0-1 \$tempARRAY=$tempARRAY \$In_aln___file=$In_aln___file\n";
		my $eachMaxClunmNB=$tempARRAY->[4] if (   (  defined ( $tempARRAY )  ) && (  ref ( $tempARRAY ) eq 'ARRAY'  )   );
		if (   (  defined ( $eachMaxClunmNB )  ) && ( $eachMaxClunmNB=~m/\d+/ ) && ( $eachMaxClunmNB > $BigestAA_Clm_Nb )   ){  
			$outTempARRAY=$tempARRAY;                                                                                     print "\n20190201-0-0-2 \$outTempARRAY=$outTempARRAY \$In_aln___file=$In_aln___file\n";
			$BigestAA_Clm_Nb=$eachMaxClunmNB;                                                      print "\n20190201-0-0-3 \$BigestAA_Clm_Nb=$BigestAA_Clm_Nb \$In_aln___file=$In_aln___file\n";
			$anyAAfoundToBeMax=1;
		}
	}
	if ( $anyAAfoundToBeMax==0 ){ DieWork::Just_dieWork( $die_MsgHead."\$In_aln___file=$In_aln___file is not right!!".$caller_inform ); }
	return $outTempARRAY;
}

sub Build_someAA_clounm_map_HASH {  #   ClustalwRun::Build_someAA_clounm_map_HASH( $In_aln___file, $In_aln_format, $AA_char );
	my ($In_aln___file, $In_aln_format, $AA_char)=@_;
	
	my $warnMsgBody="\nIn package  ClustalwRun,\tIn sub Build_someAA_clounm_map_HASH,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	my $outHash;
	my $JustCharHoldHASH;
	my $U_clum_Char_HASH;
	
	my $olnyCtHash;
	my $BgstUClmHSH;
	my $BigestU_Clm_Unb=-999999;
	my $OUT_BigestU_Clm=-99999999; 
	
	my $real___format='msf';
	if (   (  defined ( $In_aln_format )  ) && ( $In_aln_format=~m/\S+/ )   ){
		$real___format=$In_aln_format;
	}
	
	my $inAlnIO  = Bio::AlignIO->new(
	                                   -file   => $In_aln___file ,      
                                     -format => $real___format   
                                  );  
                             
  while ( my $AlnObj = $inAlnIO->next_aln() ) {
  	
  	
  	my $Seq_NBorder=$AlnObj->num_sequences;
  	
    my $BestUClmHsh;
  	my $allClumHash;
  	
  	
  	foreach my $SeqObj ($AlnObj->each_seq) {        
      my $id_Here=$SeqObj->id() ; 
      my $seqHere=$SeqObj->seq();    #print "\n\n$seqHere\n";
      my $UposAry=FastaFileHandle::Find_AA_pos_inAminoAcid($seqHere, $AA_char);
      if  (   (  defined ( $UposAry )  ) && (  ( ref ($UposAry) ) eq 'ARRAY'  )   ){
      	foreach my $ecUpos (  @{ $UposAry }  ){ #print "\$ecUpos=$ecUpos\n";      	  #my $U_pos_inAln=$SeqObj->column_from_residue_number($ecUpos);      	  
      	  $allClumHash->{$ecUpos}=1;
      	}
      	
      }      
    }
    
    my $TotalKy      ='0_0_0_TotalNub';  
    my $EachCountHead='0_1_0_NumberOf_';
    my $Seq_IdxNmHead='0_2_0_MsfOrder_';
    
    my $seqIdx=0;
    foreach my $SeqObj ($AlnObj->each_seq) {        
      my $id_Here=$SeqObj->id() ; 
      my $seqHere=$SeqObj->seq();    
      
      my $sptfIdx=MatrixCsvChange::SprintfKeyHead( $Seq_NBorder, $seqIdx );
      my $showKey=$Seq_IdxNmHead.$sptfIdx."_".$id_Here;
      
      my $Infml_abs_ID;
      if ( $id_Here =~ m/^(\d+)_/){
      	$Infml_abs_ID=$1;
      	$Infml_abs_ID= MatrixCsvChange::Change00Number_toBack( $Infml_abs_ID );
      }
      
      
      foreach my $eachClum (   keys (  %{ $allClumHash }  )   ){
      	my $char=$SeqObj->subseq($eachClum, $eachClum); 
      	my $carCtKy=$EachCountHead.$char;
      	
      	$outHash->{ $showKey }->{$eachClum}=$char;
      	
      	$JustCharHoldHASH->{$eachClum}->{$Infml_abs_ID}=$char; print "\n20181230-0 \$JustCharHoldHASH->{$eachClum}->{$Infml_abs_ID}=\$char=$char \$id_Here=$id_Here\n"; 
      	
      	if (   (  defined ( $outHash->{ $carCtKy }->{$eachClum} )  ) && ( $outHash->{ $carCtKy }->{$eachClum}=~m/\d+/ )   ){
      		$outHash->{ $carCtKy }->{$eachClum}++;   $olnyCtHash->{ $carCtKy }->{$eachClum}++;
      		
      	}
      	else {
      		$outHash->{ $carCtKy }->{$eachClum}=1;   $olnyCtHash->{ $carCtKy }->{$eachClum}=1; 
      	}
      	
      	if (   (  defined ( $outHash->{ $TotalKy }->{$eachClum} )  ) && ( $outHash->{ $TotalKy }->{$eachClum}=~m/\d+/ )   ){
      		$outHash->{ $TotalKy }->{$eachClum}++;  $olnyCtHash->{ $TotalKy }->{$eachClum}++;
      	}
      	else {
      		$outHash->{ $TotalKy }->{$eachClum}=1;  $olnyCtHash->{ $TotalKy }->{$eachClum}=1; 
      	}      	 
      }
      
      $seqIdx++;
    }
    
    my $U_key=$EachCountHead.$AA_char;  print "\n20181228-0-1 \$U_key=$U_key\n"; 
    if  (    (  defined ( $outHash )  ) && (   (  ref ( $outHash )  ) eq 'HASH'   ) && (  defined ( $outHash->{$U_key} )  ) && (   (  ref ( $outHash->{$U_key} )  ) eq 'HASH'   )    ){
    	$BigestU_Clm_Unb=-999999;            print "\n20181228-0-2 \$BigestU_Clm_Unb=$BigestU_Clm_Unb\n"; 
    	my $BigestU_Clm;
      foreach my $ecClumU_NB (   keys (  %{ $outHash->{$U_key} }  )   ){
      	my $each_clm_Unb=$outHash->{$U_key}->{$ecClumU_NB};   print "\n20181228-0-3 \$each_clm_Unb=\$outHash->{$U_key}->{$ecClumU_NB}=$each_clm_Unb\n"; 
      	if ( $BigestU_Clm_Unb < $each_clm_Unb ){
      		$BigestU_Clm_Unb=$each_clm_Unb;
      		$BigestU_Clm=$ecClumU_NB;  print "\n20181228-0-4 \$BigestU_Clm=\$U_key=$BigestU_Clm\n"; 
      		$OUT_BigestU_Clm=$BigestU_Clm;
      	}
      }
      
      if  (    (  defined ( $olnyCtHash )  ) && (   (  ref ( $olnyCtHash )  ) eq 'HASH'   ) && ( $BigestU_Clm_Unb > 0 )    ){
      	
      	foreach my $ecKey (   keys (  %{ $olnyCtHash }  )   ){
      		if (    (  defined ( $olnyCtHash->{$ecKey}->{$BigestU_Clm} )  ) && ( $olnyCtHash->{$ecKey}->{$BigestU_Clm}=~m/\S+/ )   ) {
      			$BgstUClmHSH->{$ecKey}=$olnyCtHash->{$ecKey}->{$BigestU_Clm}; 
      		}
      	  
      	}
      	
      	 $U_clum_Char_HASH=Storable::dclone( $JustCharHoldHASH->{$BigestU_Clm} );  DirFileHandle::PrintAndWarnDumper ($U_clum_Char_HASH, "\n20181230-1 ");
      	
      }
    }     
    
      	
  }                    
  #print "20190104-0-0 \$BigestU_Clm_Unb=$BigestU_Clm_Unb\n";                
	return [ $outHash, $BgstUClmHSH, $U_clum_Char_HASH, $OUT_BigestU_Clm, $BigestU_Clm_Unb ];
	
}


sub Build_U_clounm_map_HASH {  #   ClustalwRun::Build_U_clounm_map_HASH( $In_aln___file, $In_aln_format );
	my ($In_aln___file, $In_aln_format)=@_;
	
	my $warnMsgBody="\nIn package  ClustalwRun,\tIn sub Build_U_clounm_map_HASH,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	my $outHash;
	my $JustCharHoldHASH;
	my $U_clum_Char_HASH;
	
	my $olnyCtHash;
	my $BgstUClmHSH;
	my $BigestU_Clm_Unb=-999999;
	my $OUT_BigestU_Clm=-99999999; 
	
	my $real___format='msf';
	if (   (  defined ( $In_aln_format )  ) && ( $In_aln_format=~m/\S+/ )   ){
		$real___format=$In_aln_format;
	}
	
	my $inAlnIO  = Bio::AlignIO->new(
	                                   -file   => $In_aln___file ,      
                                     -format => $real___format   
                                  );  
                             
  while ( my $AlnObj = $inAlnIO->next_aln() ) {
  	
  	
  	my $Seq_NBorder=$AlnObj->num_sequences;
  	
    my $BestUClmHsh;
  	my $allClumHash;
  	
  	
  	foreach my $SeqObj ($AlnObj->each_seq) {        
      my $id_Here=$SeqObj->id() ; 
      my $seqHere=$SeqObj->seq();    #print "\n\n$seqHere\n";
      my $UposAry=FastaFileHandle::Find_U_pos_inAminoAcid($seqHere);
      if  (   (  defined ( $UposAry )  ) && (  ( ref ($UposAry) ) eq 'ARRAY'  )   ){
      	foreach my $ecUpos (  @{ $UposAry }  ){ #print "\$ecUpos=$ecUpos\n";      	  #my $U_pos_inAln=$SeqObj->column_from_residue_number($ecUpos);      	  
      	  $allClumHash->{$ecUpos}=1;
      	}
      	
      }      
    }
    
    my $TotalKy      ='0_0_0_TotalNub';  
    my $EachCountHead='0_1_0_NumberOf_';
    my $Seq_IdxNmHead='0_2_0_MsfOrder_';
    
    my $seqIdx=0;
    foreach my $SeqObj ($AlnObj->each_seq) {        
      my $id_Here=$SeqObj->id() ; 
      my $seqHere=$SeqObj->seq();    
      
      my $sptfIdx=MatrixCsvChange::SprintfKeyHead( $Seq_NBorder, $seqIdx );
      my $showKey=$Seq_IdxNmHead.$sptfIdx."_".$id_Here;
      
      my $Infml_abs_ID;
      if ( $id_Here =~ m/^(\d+)_/){
      	$Infml_abs_ID=$1;
      	$Infml_abs_ID= MatrixCsvChange::Change00Number_toBack( $Infml_abs_ID );
      }
      
      
      foreach my $eachClum (   keys (  %{ $allClumHash }  )   ){
      	my $char=$SeqObj->subseq($eachClum, $eachClum); 
      	my $carCtKy=$EachCountHead.$char;
      	
      	$outHash->{ $showKey }->{$eachClum}=$char;
      	
      	$JustCharHoldHASH->{$eachClum}->{$Infml_abs_ID}=$char; print "\n20181230-0 \$JustCharHoldHASH->{$eachClum}->{$Infml_abs_ID}=\$char=$char \$id_Here=$id_Here\n"; 
      	
      	if (   (  defined ( $outHash->{ $carCtKy }->{$eachClum} )  ) && ( $outHash->{ $carCtKy }->{$eachClum}=~m/\d+/ )   ){
      		$outHash->{ $carCtKy }->{$eachClum}++;   $olnyCtHash->{ $carCtKy }->{$eachClum}++;
      		
      	}
      	else {
      		$outHash->{ $carCtKy }->{$eachClum}=1;   $olnyCtHash->{ $carCtKy }->{$eachClum}=1; 
      	}
      	
      	if (   (  defined ( $outHash->{ $TotalKy }->{$eachClum} )  ) && ( $outHash->{ $TotalKy }->{$eachClum}=~m/\d+/ )   ){
      		$outHash->{ $TotalKy }->{$eachClum}++;  $olnyCtHash->{ $TotalKy }->{$eachClum}++;
      	}
      	else {
      		$outHash->{ $TotalKy }->{$eachClum}=1;  $olnyCtHash->{ $TotalKy }->{$eachClum}=1; 
      	}      	 
      }
      
      $seqIdx++;
    }
    
    my $U_key=$EachCountHead."U";  print "\n20181228-0-1 \$U_key=$U_key\n"; 
    if  (    (  defined ( $outHash )  ) && (   (  ref ( $outHash )  ) eq 'HASH'   ) && (  defined ( $outHash->{$U_key} )  ) && (   (  ref ( $outHash->{$U_key} )  ) eq 'HASH'   )    ){
    	$BigestU_Clm_Unb=-999999;            print "\n20181228-0-2 \$BigestU_Clm_Unb=$BigestU_Clm_Unb\n"; 
    	my $BigestU_Clm;
      foreach my $ecClumU_NB (   keys (  %{ $outHash->{$U_key} }  )   ){
      	my $each_clm_Unb=$outHash->{$U_key}->{$ecClumU_NB};   print "\n20181228-0-3 \$each_clm_Unb=\$outHash->{$U_key}->{$ecClumU_NB}=$each_clm_Unb\n"; 
      	if ( $BigestU_Clm_Unb < $each_clm_Unb ){
      		$BigestU_Clm_Unb=$each_clm_Unb;
      		$BigestU_Clm=$ecClumU_NB;  print "\n20181228-0-4 \$BigestU_Clm=\$U_key=$BigestU_Clm\n"; 
      		$OUT_BigestU_Clm=$BigestU_Clm;
      	}
      }
      
      if  (    (  defined ( $olnyCtHash )  ) && (   (  ref ( $olnyCtHash )  ) eq 'HASH'   ) && ( $BigestU_Clm_Unb > 0 )    ){
      	
      	foreach my $ecKey (   keys (  %{ $olnyCtHash }  )   ){
      		if (    (  defined ( $olnyCtHash->{$ecKey}->{$BigestU_Clm} )  ) && ( $olnyCtHash->{$ecKey}->{$BigestU_Clm}=~m/\S+/ )   ) {
      			$BgstUClmHSH->{$ecKey}=$olnyCtHash->{$ecKey}->{$BigestU_Clm}; 
      		}
      	  
      	}
      	
      	 $U_clum_Char_HASH=Storable::dclone( $JustCharHoldHASH->{$BigestU_Clm} );  DirFileHandle::PrintAndWarnDumper ($U_clum_Char_HASH, "\n20181230-1 ");
      	
      }
    }     
    
      	
  }                    
  #print "20190104-0-0 \$BigestU_Clm_Unb=$BigestU_Clm_Unb\n";                
	return [ $outHash, $BgstUClmHSH, $U_clum_Char_HASH, $OUT_BigestU_Clm ];
	
}


sub ChangeU_for_alignment_file_del_addedNB_of_muscle{ #ClustalwRun::ChangeU_for_alignment_file_del_addedNB_of_muscle
	my ( $orgFastaSeqFile, $org_Align_File, $org_format, 
	                       $out_AlignFIle,  $out_format, 
	                       $U_char_in_org_Aln, $U_char_in_new_Aln, 
	                       $cutStt, $cutEnd)=@_;
	
	my $warnMsgBody="\nIn package  ClustalwRun,\tIn sub ChangeU_for_alignment_file_del_addedNB_of_muscle,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $caller_inform=DirFileHandle::print_SubCallerInform;
  
  if  (   (  defined ( $orgFastaSeqFile )  ) && (  -e ( $orgFastaSeqFile )  ) && (  defined ( $org_Align_File )  ) && (  -e ( $org_Align_File ) )   ){
  	
  	
  	my $orgFastaSeq_HASH=ClustalwRun::NewReadInSeq_into_A_HASH( $orgFastaSeqFile , $org_format);  
  	###################################  column_from_residue_number
  	
  	my $in  = Bio::AlignIO->new(-file   => $org_Align_File ,      #
                                -format => $org_format   );   #将输入的 aln文件，解析成 对象
  
  	my $idRelationHash;
  
    my $howManyAln=0;                            
    while ( my $aln = $in->next_aln() ) {
    	
    	#临时文件 保存中间fas格式的文件
    	my $tempFile=TimeWork::GetTimeDirOrFileName().".fas";
    	
    	#建立临时文件的输入输出 对象
    	my $tempFast    = Bio::SeqIO->new(-file   => ">$tempFile" ,          #
                                        -format => 'fasta');               #临时文件的打开，用于写入fasta格式的文件, For Msf
    	
    	#循环读取 msf各条信息， 写入临时文件  	
    	my $seqIdx=1;
      foreach my $seq ($aln->each_seq) {
        #获取 id 和 msf的序列
        my $idHere=$seq->id(); my $seqHere=$seq->seq();                                   # 这里得到原来的id和比对的序列
        warn "1,".$idHere."\n";
        
        
        
        my $outIdx=$seqIdx."-".$idHere;
        #For Mega 将其中的U转换为C
        if ($out_format eq 'mega'){
        	my $changeHASH=&build_U_changeHash( $orgFastaSeq_HASH, $idHere, $U_char_in_org_Aln, $U_char_in_new_Aln, $seq);
        	$seqHere=FastaFileHandle::CHangeSpecificSequenceInAFastaSeq($seqHere, $changeHASH); warn "\n11111    $seqHere\n";
          $seqHere=~s/\./-/g;                    warn "\n22222    $seqHere\n";
        }  
        if ($out_format eq 'msf'){
        	my $changeHASH=&build_U_changeHash( $orgFastaSeq_HASH, $idHere, $U_char_in_org_Aln, $U_char_in_new_Aln, $seq);
        	$seqHere=FastaFileHandle::CHangeSpecificSequenceInAFastaSeq($seqHere, $changeHASH);                          
        } 
        
        $idHere=~s/^\d+-//; warn "2,".$idHere."\n";
        $seq->id($idHere); warn "3,".$seq->id()."\n"; #sleep (3);
        
        $idRelationHash->{'Org2Chg'}->{$idHere}=$seqIdx;
        $idRelationHash->{'Chg2Org'}->{$seqIdx}=$idHere;
        
         
        #建立新的 SeqObj
        my $sortFasSeqObj=Bio::Seq->new( -display_id => $idHere,   #$seqIdx,                          #序号 
                                                -seq => $seqHere                ); 
        #print "201811201803-1 \$cutStt=$cutStt\t$cutEnd=$cutEnd\n"; 
        if (  (  defined ( $cutStt )  ) && ( $cutStt=~m/\d+/) && ( $cutStt>0 ) && (  defined ( $cutEnd )  ) && ( $cutEnd=~m/\d+/) && ( $cutEnd>0 )   ){
        	my $truncateSeq=$sortFasSeqObj->trunc($cutStt, $cutEnd);  print "201811201803-2 \$truncateSeq=$truncateSeq\n";
        	#$sortFasSeqObj->seq()=$truncateSeq; 
        	$tempFast->write_seq($truncateSeq);   
        }
        
        else {
        	#写入
          $tempFast->write_seq($sortFasSeqObj);      
        }
                                                
                                                       
            
              
        $seqIdx++;
      }
      
      #建立读取临时文件的 IO对象
      my $fastain= Bio::AlignIO->new(-file   => $tempFile ,              ##  #For Msf #将fasta的临时文件作为 输入
                                     -format => 'fasta'             );   ##
      
      
      #建立 写入 新Msf，或新格式的文件的　ＩＯ对象
      
      my $out    = Bio::AlignIO->new(-file   => ">$out_AlignFIle" ,     #   #For Msf #建立新的输出对象
                                     -format => $out_format          );   #   
      
      ##For Msf #下面进行新输出 比对文件的 写入
      while   ( my $FstAln = $fastain->next_aln() ) {
      	$out->write_aln($FstAln);
      }                             
      system ("rm -f $tempFile");    ##删除临时文件
      
      
      $howManyAln++;
    }                         
    if ($howManyAln>1){
      my $dieMsg=$die_MsgHead.$caller_inform."\n\$org_Align_File=$org_Align_File\t\t\$howManyAln=$howManyAln\n\n\n";
  	  print $dieMsg;
  	  die $dieMsg;
      
    }                        
    return  $idRelationHash;
  	#####################################
  	
  	
  	
  	
  }
  else {
  	my $dieMsg=$die_MsgHead.$caller_inform."\$orgFastaSeqFile=$orgFastaSeqFile or \$org_Align_File=$org_Align_File is not right!!! :$!\n\n\n";
  	print $dieMsg;
  	die $dieMsg;
  }
  
  
}



sub ChangeU_for_alignment_file_del_addedNB_of_muscle_and_change_name{ #ClustalwRun::ChangeU_for_alignment_file_del_addedNB_of_muscle_and_change_name
	my (  $MsfNewNmHSH,
	      $orgFastaSeqFile, 
	                       $org_Align_File, $org_format, 
	                       $out_AlignFIle,  $out_format, 
	                       $U_char_in_org_Aln, $U_char_in_new_Aln, 
	                       $cutStt, $cutEnd)=@_;
	
	my $warnMsgBody="\nIn package  ClustalwRun,\tIn sub ChangeU_for_alignment_file_del_addedNB_of_muscle_and_change_name,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $caller_inform=DirFileHandle::print_SubCallerInform; print "\n201812241644-1 \n";
  
  if  (   (  defined ( $orgFastaSeqFile )  ) && (  -e ( $orgFastaSeqFile )  ) && (  defined ( $org_Align_File )  ) && (  -e ( $org_Align_File ) )   ){
  	
  	
  	my $orgFastaSeq_HASH=ClustalwRun::NewReadInSeq_into_A_HASH( $orgFastaSeqFile , $org_format);  
  	###################################  column_from_residue_number
  	
  	my $in  = Bio::AlignIO->new(-file   => $org_Align_File ,      #
                                -format => $org_format   );   #将输入的 aln文件，解析成 对象
  
  	my $idRelationHash;
  
    my $howManyAln=0;                            
    while ( my $aln = $in->next_aln() ) {
    	
    	#临时文件 保存中间fas格式的文件
    	my $tempFile=TimeWork::GetTimeDirOrFileName().".fas";
    	
    	#建立临时文件的输入输出 对象
    	my $tempFast    = Bio::SeqIO->new(-file   => ">$tempFile" ,          #
                                        -format => 'fasta');               #临时文件的打开，用于写入fasta格式的文件, For Msf
    	
    	#循环读取 msf各条信息， 写入临时文件  	
    	my $seqIdx=1;
      foreach my $seq ($aln->each_seq) {
        #获取 id 和 msf的序列
        my $idHere=$seq->id(); my $seqHere=$seq->seq();                                   # 这里得到原来的id和比对的序列
        warn "1,".$idHere."\n";
        
        
        
        my $outIdx=$seqIdx."-".$idHere;
        #For Mega 将其中的U转换为C
        if ($out_format eq 'mega'){
        	my $changeHASH=&build_U_changeHash( $orgFastaSeq_HASH, $idHere, $U_char_in_org_Aln, $U_char_in_new_Aln, $seq);
        	$seqHere=FastaFileHandle::CHangeSpecificSequenceInAFastaSeq($seqHere, $changeHASH); warn "\n11111    $seqHere\n";
          $seqHere=~s/\./-/g;                    warn "\n22222    $seqHere\n";
        }  
        if ($out_format eq 'msf'){
        	my $changeHASH=&build_U_changeHash( $orgFastaSeq_HASH, $idHere, $U_char_in_org_Aln, $U_char_in_new_Aln, $seq);
        	$seqHere=FastaFileHandle::CHangeSpecificSequenceInAFastaSeq($seqHere, $changeHASH);                          
        } 
        
        $idHere=~s/^\d+-//;
        if ($idHere=~m/^(\d+)_/){ 
        	my $NbKey=$1; 
        	if (   (  defined ( $MsfNewNmHSH )  ) && (  ref ( $MsfNewNmHSH ) eq 'HASH'  ) && (  defined ( $MsfNewNmHSH->{$NbKey} )  ) && ( $MsfNewNmHSH->{$NbKey}=~m/\S+/ )   ){
        		$idHere=$MsfNewNmHSH->{$NbKey};
        	}
        	else{
        		my $dieMsg_1=$die_MsgHead."\$MsfNewNmHSH->{$NbKey}=$MsfNewNmHSH->{$NbKey} is not right!!\n\n\n"; print $dieMsg_1; die $dieMsg_1; 
        	}
        }
        else{
        	my $dieMsg_2=$die_MsgHead."\$idHere=$idHere is not right!!\n\n\n"; print $dieMsg_2; die $dieMsg_2;
        }
         print  "201812241644-2,".$idHere."\n";
        $seq->id($idHere); warn "\n201812241644-3,".$seq->id()."\n"; #sleep (3);
        
        $idRelationHash->{'Org2Chg'}->{$idHere}=$seqIdx;
        $idRelationHash->{'Chg2Org'}->{$seqIdx}=$idHere;
        
         
        #建立新的 SeqObj
        my $sortFasSeqObj=Bio::Seq->new( -display_id => $idHere,   #$seqIdx,                          #序号 
                                                -seq => $seqHere                ); 
        #print "201811201803-1 \$cutStt=$cutStt\t$cutEnd=$cutEnd\n"; 
        if (  (  defined ( $cutStt )  ) && ( $cutStt=~m/\d+/) && ( $cutStt>0 ) && (  defined ( $cutEnd )  ) && ( $cutEnd=~m/\d+/) && ( $cutEnd>0 )   ){
        	my $truncateSeq=$sortFasSeqObj->trunc($cutStt, $cutEnd);  print "201811201803-2 \$truncateSeq=$truncateSeq\n";
        	#$sortFasSeqObj->seq()=$truncateSeq; 
        	$tempFast->write_seq($truncateSeq);   
        }
        
        else {
        	#写入
          $tempFast->write_seq($sortFasSeqObj);      
        }
                                                
                                                       
            
              
        $seqIdx++;
      }
      
      #建立读取临时文件的 IO对象
      my $fastain= Bio::AlignIO->new(-file   => $tempFile ,              ##  #For Msf #将fasta的临时文件作为 输入
                                     -format => 'fasta'             );   ##
      
      
      #建立 写入 新Msf，或新格式的文件的　ＩＯ对象
      
      my $out    = Bio::AlignIO->new(-file   => ">$out_AlignFIle" ,     #   #For Msf #建立新的输出对象
                                     -format => $out_format          );   #   
      
      ##For Msf #下面进行新输出 比对文件的 写入
      while   ( my $FstAln = $fastain->next_aln() ) {
      	$out->write_aln($FstAln);
      }                             
      system ("rm -f $tempFile");    ##删除临时文件
      
      
      $howManyAln++;
    }                         
    if ($howManyAln>1){
      my $dieMsg=$die_MsgHead.$caller_inform."\n\$org_Align_File=$org_Align_File\t\t\$howManyAln=$howManyAln\n\n\n";
  	  print $dieMsg;
  	  die $dieMsg;
      
    }                        
    return  $idRelationHash;
  	#####################################
  	
  	
  	
  	
  }
  else {
  	my $dieMsg=$die_MsgHead.$caller_inform."\$orgFastaSeqFile=$orgFastaSeqFile or \$org_Align_File=$org_Align_File is not right!!! :$!\n\n\n";
  	print $dieMsg;
  	die $dieMsg;
  }
  
  
}


sub build_U_changeHash { #  ClustalwRun::build_U_changeHash
	my ($FastaSeqHash, $IN_seq_ID,  $Changed_from_Char,  $Changed_into_Char, $seq_obj)=@_;
	

  my $warnMsgBody="\nIn package  ClustalwRun,\tIn sub build_U_changeHash,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $caller_inform=DirFileHandle::print_SubCallerInform;
  
  my $change_HASH;
  if  (      (  defined ( $FastaSeqHash      )  ) && (  ref ( $FastaSeqHash ) eq 'HASH'  ) 
          && (  defined ( $IN_seq_ID         )  ) && ( $IN_seq_ID=~/\S+/        ) 
          && (  defined ( $Changed_from_Char )  ) && ( $Changed_from_Char=~/\S/ ) 
          && (  defined ( $Changed_into_Char )  ) && ( $Changed_into_Char=~/\S/ )   
      )
  {
  	$Changed_into_Char=uc $Changed_into_Char;
  	if (      (  defined ( $FastaSeqHash->{'1_id_seqHash'}                )  ) && (  ref ( $FastaSeqHash->{'1_id_seqHash'}               ) eq 'HASH' ) 
  	       && (  defined ( $FastaSeqHash->{'1_id_seqHash'}->{$IN_seq_ID}  )  ) && (  ref ( $FastaSeqHash->{'1_id_seqHash'}->{$IN_seq_ID} ) eq 'HASH' ) 
  	   )
  	{
  	 	if (   (  defined ( $FastaSeqHash->{'1_id_seqHash'}->{$IN_seq_ID}->{'5_U_Pos_Array'} )  ) && (  ref ( $FastaSeqHash->{'1_id_seqHash'}->{$IN_seq_ID}->{'5_U_Pos_Array'} ) eq 'ARRAY' )    ){
  	 		 for ( my $i=0; $i<@{ $FastaSeqHash->{'1_id_seqHash'}->{$IN_seq_ID}->{'5_U_Pos_Array'} }; $i++){
  	 		 	 my $U_pos=$FastaSeqHash->{'1_id_seqHash'}->{$IN_seq_ID}->{'5_U_Pos_Array'}->[$i];
  	 		 	 my $U_pos_inAln=$seq_obj->column_from_residue_number($U_pos);
  	 		 	 $change_HASH->{$i}->{'0_directi'}='+';
  	 		 	 $change_HASH->{$i}->{'1_postion'}=$U_pos_inAln;
  	 		 	 $change_HASH->{$i}->{'2_fromCha'}=$Changed_from_Char;
  	 		 	 $change_HASH->{$i}->{'3_intoCha'}=$Changed_into_Char;
  	 		 }
  	 	}
  	}
  	else {
  		 my $dieMsg=$die_MsgHead.$caller_inform."\$FastaSeqHash->{'1_id_seqHash'}=$FastaSeqHash->{'1_id_seqHash'} or \$FastaSeqHash->{'1_id_seqHash'}->{$IN_seq_ID}=$FastaSeqHash->{'1_id_seqHash'}->{$IN_seq_ID} is not right!!! :$!\n\n\n";
  	  print $dieMsg;
  	  die $dieMsg;
  	}
    
  }
  else {
   my $dieMsg=$die_MsgHead.$caller_inform."\$FastaSeqHash=$FastaSeqHash or \$IN_seq_ID=$IN_seq_ID or \$Changed_from_Char=$Changed_from_Char or \$Changed_into_Char=$Changed_into_Char is not right!!! :$!\n\n\n";
    print $dieMsg;
    die $dieMsg;
  }
  return $change_HASH;
}



sub readInSeqAndUseClustalwOrder{  #读取一个多序列的文件，将该文件变成hash，hash给出了序列数，对每个序列给出了以该序列primary_id为Key的各种信息,这个程序增加了一个功能，利用格式是$MsfFormat的alignemnt文件$orgMsfFile中的顺序来对多序列文件进行排序
	my ($inseqFile, $inFormat, $orgMsfFile, $MsfFormat)=@_;
	my $MsfOrderHash=&GetClustalwOrderHash($orgMsfFile, $MsfFormat);
	
	my $seqInObj=Bio::SeqIO->new(-file   => $inseqFile,    
                               -format => $inFormat );
  my $outHash;
  
  my $seqNb=0;                             
  while (my $seqObj=$seqInObj->next_seq){
    #print "primary_id: ", $seqObj->primary_id, "\n";
    #print "display_id: ", $seqObj->display_id, "\n";
    #print "accession_number: ", $seqObj->accession_number, "\n";
    #print "desc: ", $seqObj->desc, "\n\n";
    
    $outHash->{'id_seqHash'}->{ $seqObj->primary_id }->{'display_id'      }=$seqObj->display_id;
    $outHash->{'id_seqHash'}->{ $seqObj->primary_id }->{'accession_number'}=$seqObj->accession_number;
    $outHash->{'id_seqHash'}->{ $seqObj->primary_id }->{'desc'            }=$seqObj->desc;
    $outHash->{'id_seqHash'}->{ $seqObj->primary_id }->{'seqOrder'        }=$MsfOrderHash->{'Org2Chg'}->{ $seqObj->primary_id };
    $outHash->{'id_seqHash'}->{ $seqObj->primary_id }->{'sequence'        }=$seqObj->seq();  
    $seqNb++;
  }
  
  $outHash->{'howManySeq'}=$seqNb;
  
  return $outHash;
	
}

sub runingClustalw{  #直接运行Clustalw
  my ($inseqFile, $outSeqFile)=@_;
  my $cmd="clustalw $inseqFile -OUTPUT=GCG -OUTFILE=$outSeqFile";  
  warn $cmd, "\n\n"; print  $cmd, "\n\n";
  system ("$cmd");
  
}

sub runingClustalwInputOrder{  #直接运行Clustalw,但比对后 输出顺序仍然按要输入的 fasta序列的顺序输出
  my ($inseqFile, $outSeqFile)=@_;
  my $cmd="clustalw $inseqFile -OUTPUT=GCG -OUTFILE=$outSeqFile -OUTORDER=INPUT";  
  warn $cmd, "\n\n"; print  $cmd, "\n\n";
  system ("$cmd");
  
}


sub GetClustalwOrderHash{   #输入一个algnment文件，获得一个Hash，给出了序列ID和顺序之间的 互相对应关系的hash
  my ($orgMsfFile, $InFormat)=@_;
  
  my $in  = Bio::AlignIO->new(-file   => $orgMsfFile ,     #
                              -format => $InFormat     );   #将输入的 aln文件，解析成 对象
  
  my $idRelationHash;
  
  my $howManyAln=0;                            
  while ( my $aln = $in->next_aln() ) {
  	  	
  	#循环读取 各条信息 	
  	my $seqIdx=1;
    foreach my $seq ($aln->each_seq) {
      #获取 id 和 msf的序列
      my $idHere=$seq->id(); my $seqHere=$seq->seq();                                   # 这里得到原来的id和比对的序列
      
      $idRelationHash->{'Org2Chg'}->{$idHere}=$seqIdx;
      $idRelationHash->{'Chg2Org'}->{$seqIdx}=$idHere;
            
      $seqIdx++;
    }
    
    $howManyAln++;
  }                         
  if ($howManyAln>1){
    die "\n\nIn package ClustalwRun    In sub GetClustalwOrderHash\n\nOnly files consisted with one Aligment could be loaded into this sub! \n\nBut :\$orgMsfFile=$orgMsfFile\t\t\$howManyAln=$howManyAln\n\n\n";
  }                        
  return  $idRelationHash;
  
}

sub clustalBackToFastaFile {  #对$inFormat格式的clusatlw结果文件$orgMsfFile，进行格式转换，转换为Fasta文件
  my ($orgMsfFile, $inFormat, $outFastaFile )=@_;
  
  my $in  = Bio::AlignIO->new(-file   => $orgMsfFile ,     #
                              -format => $inFormat    );   #将输入的 aln文件，解析成 对象
  
  my $howManyAln=0;                            
  while ( my $aln = $in->next_aln() ) {
  	
  	#建立输出 对象
  	my $tempFast    = Bio::SeqIO->new( -file   => ">$outFastaFile" ,          #
                                       -format => 'fasta'           );               #临时文件的打开，用于写入fasta格式的文件, For Msf
  	
  	#循环读取 msf各条信息， 写入临时文件  	
  	my $seqIdx=1;
    foreach my $seq ($aln->each_seq) {
      #获取 id 和 msf的序列
      my $idHere=$seq->id(); my $seqHere=$seq->seq();                                   # 这里得到原来的id和比对的序列
      $seqHere=~s/\s//g; $seqHere=~s/\n//g; 
      $seqHere=~s/\.//g;  #print "\$seqHere=$seqHere\n";
       
      #建立新的 SeqObj
      my $sortFasSeqObj=Bio::Seq->new( -display_id => $idHere,                          #序号 
                                              -seq => $seqHere                );        
      #写入
      $tempFast->write_seq($sortFasSeqObj);          
            
      $seqIdx++;
    }
    
    $howManyAln++;
  }                         
  if ($howManyAln>1){
    die "\n\$orgMsfFile=$orgMsfFile\t\t\$howManyAln=$howManyAln\n\n\n";
  }                    
  
}


sub clustalBackToFastaHash {  #对$inFormat格式的clusatlw结果文件$orgMsfFile，进行格式转换，转换为易于生成Fasta的hash，key是id，值是序列
  my ($orgMsfFile, $inFormat)=@_;
  
  my $in  = Bio::AlignIO->new(-file   => $orgMsfFile ,     #
                              -format => $inFormat    );   #将输入的 aln文件，解析成 对象
  
  my $outHash;
  
  my $howManyAln=0;                            
  while ( my $aln = $in->next_aln() ) {
  	
  	#循环读取 msf各条信息， 写入临时文件  	
  	my $seqIdx=1;
    foreach my $seq ($aln->each_seq) {
      #获取 id 和 msf的序列
      my $idHere=$seq->id(); my $seqHere=$seq->seq();                                   # 这里得到原来的id和比对的序列
      $seqHere=~s/\s//g;
      $seqHere=~s/\.//g;
      
      if (  defined ( $outHash->{$idHere} )  ){ die "\n\n\nIn package ClustalwRun, in sub clustalBackToFastaHash, \$orgMsfFile=$orgMsfFile, \$inFormat=$inFormat\n \$idHere=$idHere was found more than once!!!\n\n\n";}
      else {      $outHash->{$idHere}=$seqHere;    }      
      $seqIdx++;
    }
    
    $howManyAln++;
  }                         
  if ($howManyAln>1){
    die "\n\$orgMsfFile=$orgMsfFile\t\t\$howManyAln=$howManyAln\n\n\n";
  }                        
  return  $outHash;
  
}

sub ReWriteClustalw{  #对msf格式的clusatlw结果文件$orgMsfFile，进行格式转换，转换为$outFormat格式的$outMutiAlgnFile文件
  my ($orgMsfFile, $outMutiAlgnFile, $outFormat)=@_;
  
  my $in  = Bio::AlignIO->new(-file   => $orgMsfFile ,     #
                              -format => 'msf'        );   #将输入的 aln文件，解析成 对象
  
  my $idRelationHash;
  
  my $howManyAln=0;                            
  while ( my $aln = $in->next_aln() ) {
  	
  	#临时文件 保存中间fas格式的文件
  	my $tempFile=TimeWork::GetTimeDirOrFileName()."fas";
  	
  	#建立临时文件的输入输出 对象
  	my $tempFast    = Bio::SeqIO->new(-file   => ">$tempFile" ,          #
                                      -format => 'fasta');               #临时文件的打开，用于写入fasta格式的文件, For Msf
  	
  	#循环读取 msf各条信息， 写入临时文件  	
  	my $seqIdx=1;
    foreach my $seq ($aln->each_seq) {
      #获取 id 和 msf的序列
      my $idHere=$seq->id(); my $seqHere=$seq->seq();                                   # 这里得到原来的id和比对的序列
      
      
      
      #my $segName=$seq->id(); my $segLength=$seq->end();        
      #my $clustalwSegMentList=&GetSeqPosList_2($seq, 1 , $segLength);
      #PrintDumper('tttt.txt',$clustalwSegMentList );
      
      
      
      
      
      #For Mega 将其中的U转换为C
      if ($outFormat eq 'mega'){
        $seqHere=~s/U/C/g;   $seqHere=~s/\./-/g;                  
      } 
      
      
      
      $idRelationHash->{'Org2Chg'}->{$idHere}=$seqIdx;
      $idRelationHash->{'Chg2Org'}->{$seqIdx}=$idHere;
      
       
      #建立新的 SeqObj
      my $sortFasSeqObj=Bio::Seq->new( -display_id => $seqIdx,                          #序号 
                                              -seq => $seqHere                );        
      #写入
      $tempFast->write_seq($sortFasSeqObj);          
            
      $seqIdx++;
    }
    
    #建立读取临时文件的 IO对象
    my $fastain= Bio::AlignIO->new(-file   => $tempFile ,              ##  #For Msf #将fasta的临时文件作为 输入
                                   -format => 'fasta'             );   ##
    
    
    #建立 写入 新Msf，或新格式的文件的　ＩＯ对象
    
    my $out    = Bio::AlignIO->new(-file   => ">$outMutiAlgnFile" ,     #   #For Msf #建立新的输出对象
                                   -format => $outFormat          );   #   
    
    ##For Msf #下面进行新输出 比对文件的 写入
    while   ( my $FstAln = $fastain->next_aln() ) {
    	$out->write_aln($FstAln);
    }                             
    system ("rm -f $tempFile");    ##删除临时文件
    
    
    $howManyAln++;
  }                         
  if ($howManyAln>1){
    die "\n\$orgMsfFile=$orgMsfFile\t\t\$howManyAln=$howManyAln\n\n\n";
  }                        
  return  $idRelationHash;
  
}
#my $Whole_htmlParseOut=Interproscan::SeqIn_PrasedHashOut($clustalwInSeqFile);


sub RewriteClustalw20180408{     #某些输入，需要依赖 2018.04.08.genomePostionAnalysis.pl 中的中间数据结构 #输入一个 多序列比对的结果，变成有编号的新结果。同时输出mega结果，并返回一个hash，其中记录了 id和编号之间的关系
  my ($idx_into_protId_Hash, 
      $clutalwInfile, 
      $clutalwOutfile,  
      $clutalwOutSortfile,
      $In_fmlCluOutGdIdxNameFile, 
      $megaOutFile, 
      $inFormat, 
      $outFormat, 
      $megaFormat , 
      $bigSmallFmlNm)=@_;
   
  my $outPutHash;             #输出的hash，其中记录了 id和编号之间的关系
  
  if ($outFormat=~m/\S+/){} else {$outFormat=$inFormat;}   #如果，第四个参数没有输入的话，则不改变文件的格式
  
  
  my $in  = Bio::AlignIO->new(-file   => $clutalwInfile ,     #
                              -format => $inFormat        );  #将输入的 aln文件，解析成 对象
  my $timeNow=time();
  my $tempFile     = "${clutalwOutfile}.${timeNow}.temp.txt";                #第1个用于中间解析的 fasta格式， 因为它的存在和使用，所以这个函数不能进行并行处理。 for Msf
  my $tempFileIdx  = "${clutalwOutfile}.${timeNow}.Idx.temp.txt";            #第2个用于中间解析的 fasta格式， 因为它的存在和使用，所以这个函数不能进行并行处理。 for Idx Msf
  my $tempFileSort = "${clutalwOutfile}.${timeNow}.sot.temp.txt";            #第3个用于中间解析的 fasta格式， 因为它的存在和使用，所以这个函数不能进行并行处理。 for sort Msf  这个是没有sort 将被sort的文件
  my $sortedTempFil= "${clutalwOutfile}.${timeNow}.sorted.temp.txt";         #第4个用于中间解析的 fasta格式， 因为它的存在和使用，所以这个函数不能进行并行处理。 for sort Msf  这个是已经sort的文件
  my $tempFileTwo  = "${clutalwOutfile}.${timeNow}.2temp2.txt";              #第5个用于中间解析的 fasta格式， 因为它的存在和使用，所以这个函数不能进行并行处理。 For Mega

  ################################################################# for Msf#############################################                           
  #下面开始进行 aln的数据流解析,
  my $idx_to_newId_hash;
  my $testIdx=0;
  while ( my $aln = $in->next_aln() ) {
    
  	my $tempFast    = Bio::SeqIO->new(-file   => ">$tempFile" ,          #
                                      -format => 'fasta');               #临时文件的打开，用于写入fasta格式的文件, For Msf
                                      
    my $tempFast_Idx= Bio::SeqIO->new(-file   => ">$tempFileIdx" ,          #
                                      -format => 'fasta');               #临时文件的打开，用于写入fasta格式的文件, For Msf
                                      

  	my $tempSortFast= Bio::SeqIO->new(-file   => ">$tempFileSort" ,         #
                                   -format => 'fasta');                  #临时文件的打开，用于写入fasta格式的文件, For sort Msf

  	my $tempFastTwo = Bio::SeqIO->new(-file   => ">$tempFileTwo" ,       #
                                      -format => 'fasta');               #临时文件的打开，用于写入fasta格式的文件, For Mega
  	
  	my $testIdx2=1;  #用于标示在这个比对文件中的序数       
  	
  	#下面开始一个 aln中的 序列流解析
  	#my $sortHash;  #将msf中的序列排序 把Seq对象抽出来，放进新的hash。
  	
  	foreach my $seq ($aln->each_seq) {
      my $idHere=$seq->id(); my $seqHere=$seq->seq();                                   ##For sort Msf #重点     这里得到原来的id和比对的序列
      my $sortFasSeqObj=Bio::Seq->new( -display_id => $idHere,                          ##For sort Msf #重点     
                                              -seq => $seqHere                );        ##For sort Msf #重点     
      $tempSortFast->write_seq($sortFasSeqObj);                                         ##For sort Msf # 写入临时 fasta文件
       
      my $msfId=  $idx_into_protId_Hash->{ $idHere }->{'ShortSP'}."-"
                 .$idx_into_protId_Hash->{ $idHere }->{'SmallFml'}."-"
                 .$idx_into_protId_Hash->{ $idHere }->{'smFmlSpUCinIdx'}."-"
                 .$idx_into_protId_Hash->{ $idHere }->{'CU'};
      $idx_to_newId_hash->{$idHere}=$msfId;

      my $seqobj = Bio::Seq->new( -display_id => $msfId,                            ##For Msf #重点     重新组成fasta格式
                                         -seq => $seqHere                );         ##For Msf #重点     在这里可以进一步 进行序列的更改
      $tempFast->write_seq($seqobj);                                                 #For Msf # 写入临时 fasta文件  
      
      
      my $idxMsfId=$idHere.$msfId;           
      my $seqobj_Idx = Bio::Seq->new( -display_id => $idxMsfId,                             ##For _Idx Msf #重点     重新组成fasta格式
                                             -seq => $seqHere                );             ##For _Idx Msf #重点     在这里可以进一步 进行序列的更改
      $tempFast_Idx->write_seq($seqobj_Idx);                                                #For _IdxMsf # 写入临时 fasta文件  
      
                                       
      my $U2CSeq=$seqHere;  $U2CSeq=~s/U/C/g;   $U2CSeq=~s/\./-/g;                 ##For Mega#重点     将其中的U转换为C
      my $FmlMark;
      if ($idx_into_protId_Hash->{ $idHere }->{'SmallFml'} eq $idx_into_protId_Hash->{ $idHere }->{'BigFml'}){}
      else{ $FmlMark=$idx_into_protId_Hash->{ $idHere }->{'SmallFml'};$FmlMark=~s/^.*(\w)$/$1/; }
      my $megaId= $idx_into_protId_Hash->{ $idHere }->{'ShortSP'}
                 .$idx_into_protId_Hash->{ $idHere }->{'CU'}
                 .$idx_into_protId_Hash->{ $idHere }->{'smFmlSpUCinIdx'}
                 .$FmlMark;                                                        ##For Mega#重点     将id进行转换                                                
      my $Megaseqobj = Bio::Seq->new( -display_id => $megaId,                      ##For Mega#重点     重新组成fasta格式
                                             -seq => $U2CSeq                );     ##For Mega#重点     在这里可以进一步 进行序列的更改
                                   
                                         
      
      
      $tempFastTwo->write_seq($Megaseqobj);   #For Mega# 写入临时 fasta文件
     
      $outPutHash->{$idHere}=$testIdx2;       #最终输出的 id和 序号的 配对  
      #print "$testIdx\t$testIdx2\t\t\$idHere=$idHere\n\$seqHere=$seqHere\n\n";   #\$seq->id()=$seq->id()\n";

      $testIdx2++;
    }
    $testIdx++;
    

    
  }

    
  if ($testIdx>0){
  	
  	
  	
  	
  	
    my $fastain= Bio::AlignIO->new(-file   => $tempFile ,              ##  #For Msf #将fasta的临时文件作为 输入
                                   -format => 'fasta'             );   ##
    my $out    = Bio::AlignIO->new(-file   => ">$clutalwOutfile" ,     #   #For Msf #建立新的输出对象
                                   -format => $outFormat          );   #   
    
    ##For Msf #下面进行新输出 比对文件的 写入
    while   ( my $FstAln = $fastain->next_aln() ) {
    	$out->write_aln($FstAln);
    }                             
    system ("rm -f $tempFile");    ##删除临时文件
    
    
    my $fastain_Idx= Bio::AlignIO->new(-file   => $tempFileIdx ,                      ##  #For Msf_Idx #将fasta的临时文件作为 输入
                                       -format => 'fasta'             );              ##
    my $out_Idx    = Bio::AlignIO->new(-file   => ">$In_fmlCluOutGdIdxNameFile" ,     #   #For Msf_Idx #建立新的输出对象
                                       -format => $outFormat          );              #   
    
    ##For Msf_Idx #下面进行新输出 比对文件的 写入
    while   ( my $FstAln_Idx = $fastain_Idx->next_aln() ) {
    	$out_Idx->write_aln($FstAln_Idx);
    }                             
    system ("rm -f $tempFileIdx");    ##删除临时文件
    
    &sortFastaFile ($tempFileSort, $sortedTempFil, $idx_to_newId_hash );
    my $fastainSort= Bio::AlignIO->new(-file   => $sortedTempFil,         ##  #For sort Msf #将fasta的临时文件作为 输入
                                       -format => 'fasta'             );   ##
    my $outSort    = Bio::AlignIO->new(-file   => ">$clutalwOutSortfile" , #   #For sort Msf #建立新的输出对象
                                       -format => $outFormat          );   #
    
    
    ##For Msf #下面进行新输出 比对文件的 写入
    while   ( my $SortFstAln = $fastainSort->next_aln() ) {
    	$outSort->write_aln($SortFstAln);
    }                             
    system ("rm -f $tempFileSort");    ##删除临时文件
    system ("rm -f $sortedTempFil");    ##删除临时文件
    
    my $fastainTwo= Bio::AlignIO->new(-file   => $tempFileTwo ,           ##  #For Mega#将fasta的临时文件作为 输入
                                      -format => 'fasta'             );   ##
    my $outTwo    = Bio::AlignIO->new(-file   => ">$megaOutFile" ,        #   #For Mega#建立新的输出对象
                                      -format => $megaFormat          );  #                                  
    
    ##For Mega#下面进行新输出 比对文件的 写入
    while   ( my $FstAlnTwo = $fastainTwo->next_aln() ) {
    	$outTwo->write_aln($FstAlnTwo);
    }                             
    system ("rm -f $tempFileTwo");    ##删除临时文件
    
  }
  ####################^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^################ for Msf#############################################
  
  
  return $outPutHash;                              
  
}

sub sortFastaFile{ #对fastafile进行排序，使用场景比较特异性
	my ($inFile, $outFile, $In_idx_to_newId_hash)=@_;
	open (IN, "$inFile") or die "cannot open \$inFile=$inFile :$!\n\n";
	open (OUT,">$outFile") or die "cannot create \$outFile=$outFile :$!\n\n";
	my $sortHash;
	my $stp=$/;
	$/=">";	
	while(<IN>){
		s/^>//;s/>$//; s/\s+$//;  s/^\s+//;
	  if (/\S+/){	
	  	
	  	    
	    if (m/^(\d+)\n/){
	    	my $tpNm=$1;
	    	my $tpSeq=$_;
	    	$tpSeq=~s/^(\d+\n)//;
	    	$sortHash->{$tpNm}=$tpSeq;
	    }
	    else{
	    	die"\$_=$_\n\n";
	    }
	  }
	}
	close (IN);
		
	$/=$stp;
	foreach my $eachFSTnm (    sort { $a <=> $b } (   keys (  %{ $sortHash }  )   )    ){
		my $newId= $In_idx_to_newId_hash->{ $eachFSTnm };
	  print OUT ">$newId\n$sortHash->{$eachFSTnm}\n";
	}
	close (OUT);
}

sub PrintPngForClustalw{   #   ClustalwRun::PrintPngForClustalw  ($clustalwOUTfile, $Whole_htmlParseOut, $outPng, $msf_new_showNAMEhash);    #利用clustalw的输出msf文件$clustalwOUTfile,和interproscan模块的结果文件$Whole_htmlParseOut,绘制图片，输出到$outPng中
  
  my ($clustalwOUTfile, $Whole_htmlParseOut, $outPng, $msf_new_showNAMEhash, $uPosClumNub )=@_;
  
  
  my $warnMsgBody="\nIn package  ClustalwRun,\tIn sub PrintPngForClustalw,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $caller_inform=DirFileHandle::print_SubCallerInform; print "\n201812241644-1 \n";
  
  
  my $shownNameLengthLimit=100; #当设定了 $msf_new_showNAMEhash 的时候， 可以控制这个shown name的长度，不超过 $shownNameLengthLimit
  
  #读入clustalw输出文件
  my $in  = Bio::AlignIO->new(-file   => $clustalwOUTfile ,
                              -format => 'msf'            );


  my $tempMarkPngFile="${outPng}.mark.png";
  my $tempTreePngFile="${outPng}.tree.png";
  my $tempDomnPngFile="${outPng}.doms.png";
  
  my $tempMgPng      ="${outPng}.tpMG.png";
  
  
  my $out_1; my $out_2; my $out_3;
  open ($out_1,">$tempMarkPngFile") or die "cannot create \$tempMarkPngFile=$tempMarkPngFile : $!";
  open ($out_2,">$tempTreePngFile") or die "cannot create \$tempTreePngFile=$tempTreePngFile : $!";
  open ($out_3,">$tempDomnPngFile") or die "cannot create \$tempDomnPngFile=$tempDomnPngFile : $!";
  #open (OUT,">$outPng") or die "cannot create \$outPng=$outPng : $!";
  my $trackHash; #建立一个hash，用来将需要绘图的Track装进去，以便后面的 绘图安排
  
 
  
  #读取clustalw 各个比对 的数据流，一个数据流 就是一个aln
  while (my $aln=$in->next_aln){
    print "1: ", $aln->length, "\n";   
    print "2: ", $aln->num_residues, "\n";  
    print "3: ", $aln->is_flush, "\n";  
    print "4: ", $aln->num_sequences, "\n";  
    print "5: ", $aln->percentage_identity, "\n";  
    print "6: ", $aln->consensus_string(50), "\n";
    
    
  
    my $DomainIDHash;              #这个Hash是用来 将所有的Domain ID，当做Key放入一个Hash中。  用于后面不同颜色标注的操作。
    
    my $DS_DomDadSonRelationHash;     #这个Hash是用来 将 有 Daddy-son relation的Domain 记录下来。
    my $DS_AllDadSonHash;             #这个Hash是用来 将 所有 ID，不管是dad还是son都记下来，和上面这个hash一起 ，建立关系树
    
    my $Bigest_showNameLength=-9999999999999999;  #这个是用来把所有名字对齐的
    
    #下面画每个 对齐对象的 track
    my $alnSeqIdx=0;
    foreach my $seq ($aln->each_seq) {          #warn "id: ",$seq->id(),"\n";     warn "Seq: ",$seq->seq(),"\n";    warn "Start: ",$seq->start(),"\n";    warn "End: ",$seq->end(),"\n";    warn "Number of gaps 1 .: ", $seq->num_gaps('.'),"\n";    warn "Number of gaps 2 -: ", $seq->num_gaps('-'),"\n";    warn "Number of gaps 2  : ", $seq->num_gaps(),"\n";            #warn "Col Number of rsd 20  : ", $seq->column_from_residue_number(20),"\n";    #warn "location_from_column 82  : ", $seq->location_from_column(82),"\n";       #warn "\n";    #my @tempMap=$seq->mapping();    #my $idxNb=0;foreach my $temp(@tempMap){warn "$idxNb: \$temp=$temp\t\t";$idxNb++}    #warn "\nAfter Print \@tempMap\n\n";        #my %hashTp=$seq->frameshifts();    #foreach my $keyH (keys (%hashTp)){warn "\$hashTp{\$keyH}=\$hashTp{$keyH}=$hashTp{$keyH}\t\t";} warn "\nAfter Print \%hashTp\n\n\n\n";        #my $loc82loc=$seq->location_from_column(82);    #if (defined ($loc82loc)){ warn "Loc82 LocationType: ",$loc82loc->location_type," Start: ",$loc82loc->start,"\tEnd: ",$loc82loc->end,"\nAfter Print \$loc82loc\n\n\n\n\n";}    
      
      
      my $segName=$seq->id();   
      my $segLength=$seq->end();       warn "\n\n\n\n1(1)-1-1-1\n\$segName=$segName\t\$segLength=$segLength\n\n"; print "\n\n\n\n1(1)-1-1-1\n\$segName=$segName\t\$segLength=$segLength\n\n";   
       
      my $showName=$segName; print "\n201812241707-0 \$segName=$segName\n";
      if ($segName=~m/^(\d+)_/){ 
      	my $NbKey=$1;  print "\n201812241707-1 \$segName=$segName \$NbKey=$NbKey \$msf_new_showNAMEhash->{$NbKey}=$msf_new_showNAMEhash->{$NbKey}\n";
      	if (   (  defined ( $msf_new_showNAMEhash )  ) && (  ref ( $msf_new_showNAMEhash ) eq 'HASH'  ) && (  defined ( $msf_new_showNAMEhash->{$NbKey} )  ) && ( $msf_new_showNAMEhash->{$NbKey}=~m/\S+/ )   ){
      		$showName=$msf_new_showNAMEhash->{$NbKey};  print "\n201812241707-2 \$showName=\$msf_new_showNAMEhash->{$NbKey}=$msf_new_showNAMEhash->{$NbKey}\n";
      	}
      	else{
      		#$die_MsgHead.
      		my $dieMsg_1="\$msf_new_showNAMEhash->{$NbKey}=$msf_new_showNAMEhash->{$NbKey} is not right!!\n\n\n"; print $dieMsg_1; warn $dieMsg_1; 
      	}
      }
      else{
      	my $dieMsg_2=$die_MsgHead."\$segName=$segName maybe not right!!\n\n\n"; #print $dieMsg_2; die $dieMsg_2;
      }
      
      #然后画各个比对结果
           
      my $clustalwSegMentList=&GetSeqPosList_2($seq, 1 , $segLength);
      #PrintDumper('tttt.txt',$clustalwSegMentList )   
      
      #$panel=&NEWAddFeatureTrack3Panel($clustalwSegMentList,$panel,'','',"$segName  $segLength");
      
      $trackHash->{'_alnSeqArray'}->[$alnSeqIdx]->{'_SegmentList'}   =$clustalwSegMentList;
      $trackHash->{'_alnSeqArray'}->[$alnSeqIdx]->{'_SegmentKeyName'}="$segName  $segLength";
      $trackHash->{'_alnSeqArray'}->[$alnSeqIdx]->{'_Segmentcolor'}  ='';
      $trackHash->{'_alnSeqArray'}->[$alnSeqIdx]->{'_SegmentProID'}  =$segName;   print "\n1(2)-1-1-1\n\$trackHash->{'_alnSeqArray'}->[$alnSeqIdx]->{'_SegmentProID'}  =\$segName=$segName\n"; #warn "\$trackHash->{'_alnSeqArray'}->[$alnSeqIdx]->{'_SegmentProID'}  =\$segName=$segName\n";
      $trackHash->{'_alnSeqArray'}->[$alnSeqIdx]->{'_SegmentShwID'}  =$showName;  print "\n1(2)-1-1-1\n\$trackHash->{'_alnSeqArray'}->[$alnSeqIdx]->{'_SegmentShwID'}  =\$showName=$showName\n"; #warn "\$trackHash->{'_alnSeqArray'}->[$alnSeqIdx]->{'_SegmentProID'}  =\$segName=$segName\n";
      my $showNmLength=length ( $showName );
      if ( $Bigest_showNameLength < $showNmLength ) { $Bigest_showNameLength=$showNmLength;}
      
      warn "\n1(2)-1-1-1\n\$trackHash->{'_alnSeqArray'}->[$alnSeqIdx]->{'_SegmentProID'}  =\$segName=$segName\n"; warn "1(3)-1-1-1\n\$trackHash->{'_alnSeqArray'}->[$alnSeqIdx]->{'_SegmentProID'}  =\$segName=$segName\n";
      
      
      
      #然后画出所有 domin
     
      #找到 html显示的domain文件，再进行解析等后续工作
       
      #print "1111111111110000000000000\$Whole_htmlParseOut->{$segName}=$Whole_htmlParseOut->{$segName}\n";
      
      if (  ( ref ( $Whole_htmlParseOut->{$segName}->{'Domains and repeats'}->{'DomainLineNBhash'} ) ) eq 'HASH' ){
      	my $DomainLineNBhash=$Whole_htmlParseOut->{$segName}->{'Domains and repeats'}->{'DomainLineNBhash'};    print "\n\n1-if2-1-1\$DomainLineNBhash=\$Whole_htmlParseOut->{$segName}->{'Domains and repeats'}->{'DomainLineNBhash'}=$DomainLineNBhash\n";  warn "\n\n1-if2-1-1\$DomainLineNBhash=\$Whole_htmlParseOut->{$segName}->{'Domains and repeats'}->{'DomainLineNBhash'}=$DomainLineNBhash\n";  
        foreach my $DomainLineNb (    sort {$a <=> $b} (   keys (  %{ $DomainLineNBhash }  )   )    ){
          
          
          if (  ( ref ( $DomainLineNBhash->{$DomainLineNb}->{'CavReg2DomainIDHash'} ) ) eq 'HASH' ){
          	my $CavReg2DomainIDHash=$DomainLineNBhash->{$DomainLineNb}->{'CavReg2DomainIDHash'};  print "\n\n\n1-2-foreachIf3-1\$CavReg2DomainIDHash=\$DomainLineNBhash->{$DomainLineNb}->{'CavReg2DomainIDHash'}=$CavReg2DomainIDHash\n"; warn "\n\n\n1-2-foreachIf3-1\$CavReg2DomainIDHash=\$DomainLineNBhash->{$DomainLineNb}->{'CavReg2DomainIDHash'}=$CavReg2DomainIDHash\n";
          	
          	my $dmListHere; my $dlIdx=0; 
          	my $domianSegPosHashList;
          	
            foreach my $domainRegion (    sort { $CavReg2DomainIDHash->{$a}->{'RegionStart'} <=> $CavReg2DomainIDHash->{$b}->{'RegionStart'} } (   keys (  %{ $CavReg2DomainIDHash }  )   )    ){
              
              
              
              my $DomainStt=$CavReg2DomainIDHash->{$domainRegion}->{'RegionStart'};
              my $DomainEnd=$CavReg2DomainIDHash->{$domainRegion}->{'RegionEnd'};
              
              my $DomainInD=$CavReg2DomainIDHash->{$domainRegion}->{'DomainID'};  
              my $DomainUrl=$DomainLineNBhash->{$DomainLineNb}->{'DomainIDhash'}->{$DomainInD}->{'Url'       };
              my $DomainDef=$DomainLineNBhash->{$DomainLineNb}->{'DomainIDhash'}->{$DomainInD}->{'Definition'}; print "1-2-3-foreach4(1)\$DomainDef=$DomainDef\t\t\$DomainInD=\$CavReg2DomainIDHash->{$domainRegion}->{'DomainID'}=$CavReg2DomainIDHash->{$domainRegion}->{'DomainID'}\n"; warn  "1-2-3-foreach4(1)\$DomainDef=$DomainDef\t\t\$DomainInD=\$CavReg2DomainIDHash->{$domainRegion}->{'DomainID'}=$CavReg2DomainIDHash->{$domainRegion}->{'DomainID'}\n";
              
              
              
               
               
                                                                                                                                                                                                            #$dmListHere->[$dlIdx]->[0]=     $DomainStt;                #$dmListHere->[$dlIdx]->[1]=     $DomainEnd;                 #$dmListHere->[$dlIdx]->[2]=     $DomainInD;  #家族id               #$dmListHere->[$dlIdx]->[3]=     $DomainStt;  #家族描述
              if ($DomainEnd>$segLength){$DomainEnd=$segLength; }  #这句语句的含义，主要是 当Coil等 patScan的结果 有可能比蛋白序列本身更长，这种情况，会导致 解析 msf的模块出错。故在这里进行强行的 将相关片段的长度进行减少。
              my $TempList=&GetSeqPosList_3($seq,$DomainStt,$DomainEnd,'_INCDS');                                                                                                                                                          #print "Print Temp seq List: ", &PrintList( $TempList ), "\n\n\n\n\n";                #$dmListHere->[$dlIdx]->[4]=     $TempList;          	            	                #$dmListHere->[$dlIdx]->[5]=     $DomainUrl;
              
              $domianSegPosHashList->[$dlIdx]=     $TempList;
              
              
              $trackHash->{'_alnSeqArray'}->[$alnSeqIdx]->{'_DomainInformHas'}->{$DomainLineNb}->{'domainRegionArray'}->[$dlIdx]->{'_SegmentList'}     =$TempList;
              $trackHash->{'_alnSeqArray'}->[$alnSeqIdx]->{'_DomainInformHas'}->{$DomainLineNb}->{'domainRegionArray'}->[$dlIdx]->{'_SegmentDmID'}     =$DomainInD;
              $trackHash->{'_alnSeqArray'}->[$alnSeqIdx]->{'_DomainInformHas'}->{$DomainLineNb}->{'domainRegionArray'}->[$dlIdx]->{'_Segmentcolor'}    ='red';
              $trackHash->{'_alnSeqArray'}->[$alnSeqIdx]->{'_DomainInformHas'}->{$DomainLineNb}->{'domainRegionArray'}->[$dlIdx]->{'_SegmentDmDf'}     =$DomainDef;
              $trackHash->{'_alnSeqArray'}->[$alnSeqIdx]->{'_DomainInformHas'}->{$DomainLineNb}->{'domainRegionArray'}->[$dlIdx]->{'_ProteinSqNm'}     =$segName;
              
              
              
              $DomainIDHash->{$DomainInD}->{'Definition' }=$DomainDef;                      #####################
              $DomainIDHash->{$DomainInD}->{'Url'        }=$DomainUrl;                      ##################### 
              $DomainIDHash->{$DomainInD}->{'seqNameHash'}->{$segName}=1;                  #####################
              
              
              
              $DS_AllDadSonHash->{$DomainInD}=1;    print "\n\n11\$DS_AllDadSonHash->{$DomainInD}=$DS_AllDadSonHash->{$DomainInD}\n";
              if (  ( ref ( $DomainLineNBhash->{$DomainLineNb}->{'DomainIDhash'}->{$DomainInD}->{'subDomainsHash'} ) ) eq 'HASH' ){
              	my $subDomHash =$DomainLineNBhash->{$DomainLineNb}->{'DomainIDhash'}->{$DomainInD}->{'subDomainsHash'};
                foreach my $subDmID (    sort { $a cmp $b } (   keys (  %{ $subDomHash }  )   )    ){
                	
                	
                  
                  my $subDef=$subDomHash->{$subDmID}->{'Definition'};
                  my $subUrl=$subDomHash->{$subDmID}->{'Url'};
                  
                  $DS_DomDadSonRelationHash->{$DomainInD}->{$subDmID}=1;  print "\n\$DS_DomDadSonRelationHash->{$DomainInD}->{$subDmID}=$DS_DomDadSonRelationHash->{$DomainInD}->{$subDmID}\n";
                	$DS_AllDadSonHash->{$subDmID}=1;  print "22\$DS_AllDadSonHash->{$subDmID}=$DS_AllDadSonHash->{$subDmID}\n\n\n";
                  
                  
                  $trackHash->{'_alnSeqArray'}->[$alnSeqIdx]->{'_DomainInformHas'}->{$DomainLineNb}->{'domainRegionArray'}->[$dlIdx]->{'_SegmentDmID'}     =$subDmID;
                  $trackHash->{'_alnSeqArray'}->[$alnSeqIdx]->{'_DomainInformHas'}->{$DomainLineNb}->{'domainRegionArray'}->[$dlIdx]->{'_SegmentDmDf'}     =$subDef;
                  
                  $DomainIDHash->{$subDmID}->{'Definition' }=$subDef;                      #####################
                  $DomainIDHash->{$subDmID}->{'Url'        }=$subUrl;                      ##################### 
                  $DomainIDHash->{$subDmID}->{'seqNameHash'}->{$segName}=1;               #####################
                  
                }
              }
              else{
                
                $DomainIDHash->{$DomainInD}->{'Definition' }=$DomainDef;                      #####################
                $DomainIDHash->{$DomainInD}->{'Url'        }=$DomainUrl;                      ##################### 
                $DomainIDHash->{$DomainInD}->{'seqNameHash'}->{$segName}=1;                  #####################
                
              }
          
              
              
              
              print "1-2-3-foreach4(2)\$trackHash->{'_alnSeqArray'}->[$alnSeqIdx]->{'_DomainInformHas'}->{$DomainLineNb}->{'domainRegionArray'}->[$dlIdx]->{'_SegmentDmDf'}=\$DomainDef=$DomainDef\n";
              warn "1-2-3-foreach4(2)\$trackHash->{'_alnSeqArray'}->[$alnSeqIdx]->{'_DomainInformHas'}->{$DomainLineNb}->{'domainRegionArray'}->[$dlIdx]->{'_SegmentDmDf'}=\$DomainDef=$DomainDef\n";
              
              $dlIdx++;
              
            }
            
            
          }
        }
      }
        
    
      $alnSeqIdx++;
    }
    
    if (  ( ref ( $DomainIDHash ) ) eq 'HASH' ){
      foreach my $DomIDe (   keys (  %{ $DomainIDHash }  )   ){
        $DomainIDHash->{$DomIDe}->{'times'}=(   keys (  %{ $DomainIDHash->{$DomIDe}->{'seqNameHash'} }  )   );
      }
    }
    
    #将所有的颜色字符串，写入到 $DomainIDHash这个Ｈａｓｈ中
    my @colorArray=('yellow','red','blue','mistyrose','darkcyan','orange','purple','greenyellow','hotpink','cyan','brown','fuchsia','pink','grey','rosybrown'); #my $colIdx=0;
    if (  ( ref ( $DomainIDHash ) ) eq 'HASH' ){
    	my $DmIDidx=0;
      foreach my $DmID (    sort { $DomainIDHash->{$b}->{'times'} <=> $DomainIDHash->{$a}->{'times'} } (   keys (  %{ $DomainIDHash }  )   )    ){
        my $colorNB=$DmIDidx%@colorArray;  print "\$DomainIDHash->{$DmID}->{'times'}=$DomainIDHash->{$DmID}->{'times'}\t\$DmIDidx=$DmIDidx\t\$colorNB=$colorNB\n";
        $DomainIDHash->{$DmID}->{'Domaincolor'}=$colorArray[$colorNB];
        $DmIDidx++;
      }
    }
    
    
    #建立 Domain之间的 关系树形结构
    my $domainTreeHash=&MakeTreeFrom_DadSonHashs_onlyMultiLevesTree($DS_AllDadSonHash, $DS_DomDadSonRelationHash);
    
    my $outTreeWord; $outTreeWord=&PrintTreePoints ($outTreeWord, $domainTreeHash, $DomainIDHash);
    #print "\n\n\$outTreeWord=$outTreeWord\n\n\n";
   
    #首先画 panel_1
    my $panel_1 = Bio::Graphics::Panel->new(
                                               -length    => $aln->length,        #图的全长宽度能表示的碱基数量总长度
                                               -width     => 400,                #图的实际宽度
                                               -key_style => 'right',
                                               -start     => 1,
                                               -pad_left  => 100,
                                               -pad_right => 100,
                                               
                                         );  
    
                                         
    #画图示                                      
    if (  ( ref ( $DomainIDHash ) ) eq 'HASH' ){    	
      foreach my $DmID_2 (    sort { $DomainIDHash->{$b}->{'times'} <=> $DomainIDHash->{$a}->{'times'} } (   keys (  %{ $DomainIDHash }  )   )    ){
        
        $panel_1->add_track ( generic => Bio::SeqFeature::Generic->new(-start=> ($aln->length)*5/100,
                                                                       -end  => ($aln->length)*95/100,
                                                                      #-id   => $DomainIDHash->{$DmID_2}->{'Definition'},
                                                                      #-name => $DomainIDHash->{$DmID_2}->{'Definition'}
                                                                      ),
                            -glyph  => 'generic',
                            -bgcolor=> $DomainIDHash->{$DmID_2}->{'Domaincolor'},
                            -key    => "$DomainIDHash->{$DmID_2}->{'times'}($DomainIDHash->{$DmID_2}->{'Domaincolor'}) $DmID_2: $DomainIDHash->{$DmID_2}->{'Definition'}",
                            -label  => 1,
                            
                            
                      );
        
        
      }
    }
    
    #首先画 panel_2
    my $panel_2 = Bio::Graphics::Panel->new(
                                               -length    => $aln->length,        #图的全长宽度能表示的碱基数量总长度
                                               -width     => 200,                  #图的实际宽度
                                               -key_style => 'right',
                                               -start     => 1,
                                               -pad_left  => 100,
                                               -pad_right => 100,
                                               
                                         );  
    
    $panel_2=&treePanelWork($panel_2, $domainTreeHash, $DomainIDHash, $aln->length);
    
    #首先画 panel_3
    my $panel_3 = Bio::Graphics::Panel->new(
                                               -length    => $aln->length,        #图的全长宽度能表示的碱基数量总长度
                                               -width     => 1000,                #图的实际宽度
                                               -key_style => 'left',
                                               -start     => 1,
                                               -pad_top  => 20,
                                               -pad_left  => 100,
                                               -pad_right => 100,
                                               
                                         );  
                                         
    
   #print "\n20190104-0-0 \$uPosClumNub=$uPosClumNub\n";                   
   if (   (  defined ( $uPosClumNub )  ) && ( $uPosClumNub=~m/^\d+$/) && ( $uPosClumNub >= 0 ) && ( $uPosClumNub <= $aln->length )   ){   print "\n20190104-0-1 \$uPosClumNub=$uPosClumNub\n";
   	 my $U_pos_show = Bio::SeqFeature::Generic->new(
                                                         -start => $uPosClumNub, 
                                                         -end   => $uPosClumNub
                                                    );
     $panel_3->add_track ( $U_pos_show,
                                              -glyph   => 'diamond',
                                              #-tick    => 2,
                                              -fgcolor => 'green',
                                              #-double  => 1
                         );                              
                                                    
                                                    
                                                    
   }
    
    #然后画 标线
    my $full_length = Bio::SeqFeature::Generic->new(
                                                         -start => 1, 
                                                         -end   => $aln->length
                                                    );
                                                    
                                                    
                                                    
    # 将标尺的 track加入到panel中。
    $panel_3->add_track ( $full_length,
                                           -glyph   => 'arrow',
                                           -tick    => 2,
                                           -fgcolor => 'black',
                                           -double  => 1
                      );
   
   
    
    
    
    my $colorChangeNB=0;
    
    #进行实际的绘画，画到具体的图上  $trackHash->{'_alnSeqArray'}->[$alnSeqIdx]->{'_DomainInformHas'}->{$DomainLineNb}->{'domainRegionArray'}->[$dlIdx]->{'_SegmentList'}   $trackHash->{'_alnSeqArray'}->[$alnSeqIdx]->{'_SegmentList'}
    if (  ( ref ( $trackHash->{'_alnSeqArray'} ) ) eq 'ARRAY' ){
      for( my $alnSeqIdx_2=0; $alnSeqIdx_2 < @{ $trackHash->{'_alnSeqArray'} }; $alnSeqIdx_2++ ){  #print "22222222221111111\$trackHash->{'_alnSeqArray'}->[$alnSeqIdx_2]->{'_DomainInformHas'}=$trackHash->{'_alnSeqArray'}->[$alnSeqIdx_2]->{'_DomainInformHas'}\n";
        
        my $tkColore='white'; if (($alnSeqIdx_2%2) == 0){$tkColore='whitesmoke';} #'whitesmoke', #$colorIdx, #'#FFF8DC',#'lightyellow',#'cyan',#'gray'
        my $ExistSubLev=0;  
        my $eachGeneTrack;
        
        #my $segueName =$trackHash->{'_alnSeqArray'}->[$alnSeqIdx_2]->{'_SegmentProID'};  #warn "\$segueName=\$trackHash->{'_alnSeqArray'}->[$alnSeqIdx_2]->{'_SegmentProID'}=$segueName}\n";   print "4444444444444444444444\$segueName=\$trackHash->{'_alnSeqArray'}->[$alnSeqIdx_2]->{'_SegmentProID'}=$segueName\n";
        my $segueName =$trackHash->{'_alnSeqArray'}->[$alnSeqIdx_2]->{'_SegmentShwID'};  #warn "\$segueName=\$trackHash->{'_alnSeqArray'}->[$alnSeqIdx_2]->{'_SegmentProID'}=$segueName}\n";   print "4444444444444444444444\$segueName=\$trackHash->{'_alnSeqArray'}->[$alnSeqIdx_2]->{'_SegmentProID'}=$segueName\n";
        if ( $Bigest_showNameLength > $shownNameLengthLimit){
        	$Bigest_showNameLength=$shownNameLengthLimit;
        	$segueName=substr ($segueName, 0, $Bigest_showNameLength);
        }
        $segueName=sprintf ("%-${Bigest_showNameLength}s",$segueName);  print "201812241739-0 \$segueName=$segueName=\$segueName\n";
        
        my $msfSeqList=$trackHash->{'_alnSeqArray'}->[$alnSeqIdx_2]->{'_SegmentList'};
        ($panel_3, $eachGeneTrack)=@{ &AddTracks2Panel($panel_3, $segueName, $tkColore) };
        ($panel_3, $eachGeneTrack)=@{ &AddFeature2PanelAndTrack($msfSeqList,$panel_3,$eachGeneTrack, 'white') };
        
        
        if (  ( ref ( $trackHash->{'_alnSeqArray'}->[$alnSeqIdx_2]->{'_DomainInformHas'} ) ) eq 'HASH' ){
          foreach my $DomainLineNb_2 (    sort {$a <=> $b} (   keys (  %{ $trackHash->{'_alnSeqArray'}->[$alnSeqIdx_2]->{'_DomainInformHas'} }  )   )    ){    print "22222222222223333333\$trackHash->{'_alnSeqArray'}->[$alnSeqIdx_2]->{'_DomainInformHas'}->{$DomainLineNb_2}->{'domainRegionArray'}=$trackHash->{'_alnSeqArray'}->[$alnSeqIdx_2]->{'_DomainInformHas'}->{$DomainLineNb_2}->{'domainRegionArray'}\n";
             
            if (  ( ref ( $trackHash->{'_alnSeqArray'}->[$alnSeqIdx_2]->{'_DomainInformHas'}->{$DomainLineNb_2}->{'domainRegionArray'} ) ) eq 'ARRAY' ){
              for ( my $dlIdx_2=0; $dlIdx_2 < @{ $trackHash->{'_alnSeqArray'}->[$alnSeqIdx_2]->{'_DomainInformHas'}->{$DomainLineNb_2}->{'domainRegionArray'} }; $dlIdx_2++ ){   print "2222222224444444\$dlIdx_2=$dlIdx_2\n";
                
                my $SegmeList=$trackHash->{'_alnSeqArray'}->[$alnSeqIdx_2]->{'_DomainInformHas'}->{$DomainLineNb_2}->{'domainRegionArray'}->[$dlIdx_2]->{'_SegmentList'};
                my $DomainIDe=$trackHash->{'_alnSeqArray'}->[$alnSeqIdx_2]->{'_DomainInformHas'}->{$DomainLineNb_2}->{'domainRegionArray'}->[$dlIdx_2]->{'_SegmentDmID'};
                my $DomainDef=$trackHash->{'_alnSeqArray'}->[$alnSeqIdx_2]->{'_DomainInformHas'}->{$DomainLineNb_2}->{'domainRegionArray'}->[$dlIdx_2]->{'_SegmentDmDf'};
                my $SegmeColo=$DomainIDHash->{$DomainIDe}->{'Domaincolor'};
                #my $segueName=$trackHash->{'_alnSeqArray'}->[$alnSeqIdx_2]->{'_DomainInformHas'}->{$DomainLineNb_2}->{'domainRegionArray'}->[$dlIdx_2]->{'_ProteinSqNm'};
                
                #print "2222222222225555555555555555\$DomainIDe=\$trackHash->{'_alnSeqArray'}->[$alnSeqIdx_2]->{'_DomainInformHas'}->{$DomainLineNb_2}->{'domainRegionArray'}->[$dlIdx_2]->{'_SegmentDmID'}=$DomainIDe=$trackHash->{'_alnSeqArray'}->[$alnSeqIdx_2]->{'_DomainInformHas'}->{$DomainLineNb_2}->{'domainRegionArray'}->[$dlIdx_2]->{'_SegmentDmID'}\n";
                #sprint "22222222227777\$alnSeqIdx_2=$alnSeqIdx_2\t\$DomainIDe=$DomainIDe\t\$tkColore=$tkColore\t\$segueName=$segueName\n";
                #$panel=&NEWAddFeatureTrack3Panel( $SegmeList,
                #                                  $panel,
                #                                  "$DomainIDe  $DomainDef",
                #                                 
                #                                  $SegmeColo, 
                #                                  $segueName,
                #                                  $tkColore
                #                                  );
                                                  
                #if ($ExistSubLev==0) {                	
                  #($panel, $eachGeneTrack)=@{ &AddTracks2Panel($panel, $segueName, $tkColore) };
                #}                                  
                ($panel_3, $eachGeneTrack)=@{ &AddFeature2PanelAndTrack($SegmeList, $panel_3, $eachGeneTrack, $SegmeColo, "$DomainIDe  $DomainDef") };
                
                                                  
                $ExistSubLev=1;

              }
            }
            
          }
        } 
        if ($ExistSubLev==1){$colorChangeNB++;}
        
      }
    }
    
    print $out_1  $panel_1->png;   #  $tempMarkPngFile
    print $out_2  $panel_2->png;   #  $tempTreePngFile
    print $out_3  $panel_3->png;   #  $tempDomnPngFile
    
    close ($out_1);
    close ($out_2);
    close ($out_3);
   
    #&MergeImage ("$tempMarkPngFile", "$tempDomnPngFile", $outPng);
    #my $tempMgPng;
    &MergeImage ($tempMarkPngFile, $tempTreePngFile, $tempMgPng);
    &MergeImage ($tempMgPng,       $tempDomnPngFile, $outPng);
    if (  -e ($tempMarkPngFile)  ) {system ("rm -rf $tempMarkPngFile"); }
    if (  -e ($tempTreePngFile)  ) {system ("rm -rf $tempTreePngFile"); }
    if (  -e ($tempDomnPngFile)  ) {system ("rm -rf $tempDomnPngFile"); }
    if (  -e ($tempMgPng      )  ) {system ("rm -rf $tempMgPng"      ); }
    
    
    #&MergeImage ($panel_1->gd(), $panel_3->gd(), $outPng); 
  
  }
    
}

sub MakeTreeFrom_DadSonHashs {          #利用父子关系建立树形结构           #($DS_AllDadSonHash, $DS_DomDadSonRelationHash)
  my ($AllPointKeyHash, $DadSonRelationHash)=@_;
  my $onlyDadPointHash=$AllPointKeyHash; foreach my $k (keys (%$onlyDadPointHash)){print "temptemp$k=$k\n";}
  if (  ( ref ( $DadSonRelationHash ) ) eq 'HASH' ){
    foreach my $DadKey (    sort {$a cmp $b} (   keys (  %{ $DadSonRelationHash }  )   )    ){  
    	
    	if (  ( ref ( $DadSonRelationHash->{$DadKey} ) ) eq 'HASH' ){
        foreach my $SonKey (    sort {$a cmp $b} (   keys (  %{ $DadSonRelationHash->{$DadKey} }  )   )    ){  
          warn "111 \$onlyDadPointHash->{$SonKey}=$onlyDadPointHash->{$SonKey}\n\n";  print  "111 \$onlyDadPointHash->{$SonKey}=$onlyDadPointHash->{$SonKey}\n\n";
          delete $onlyDadPointHash->{$SonKey}; 
          #$onlyDadPointHash->{$SonKey}='!S!o!n!'; 
          warn "222 \$onlyDadPointHash->{$SonKey}=$onlyDadPointHash->{$SonKey}\n\n";  print  "222 \$onlyDadPointHash->{$SonKey}=$onlyDadPointHash->{$SonKey}\n\n";
       
        }
      }
    	
    	
    }
  }
  
  my $finalTree;
  #$finalTree=&DiGuiDadSon($finalTree, $onlyDadPointHash, $DadSonRelationHash);
  if (  ( ref ( $onlyDadPointHash ) ) eq 'HASH' ){
    foreach my $DadKey (    sort {$a cmp $b} (   keys (  %{ $onlyDadPointHash }  )   )    ){   print "\$DadKey=$DadKey\n\n";
    	my $eachTree;
    	$eachTree=&DiGuiDadSon($eachTree, $DadKey, $DadSonRelationHash);
    	$finalTree->{$DadKey}=$eachTree->{$DadKey};
    }
  }
  
  return $finalTree;
}


sub MakeTreeFrom_DadSonHashs_onlyMultiLevesTree {    #利用父子关系建立树形结构                 #($DS_AllDadSonHash, $DS_DomDadSonRelationHash)
  my ($AllPointKeyHash, $DadSonRelationHash)=@_;
  my $onlyDadPointHash;#=$AllPointKeyHash; foreach my $k (keys (%$onlyDadPointHash)){print "temptemp$k=$k\n";}
  my $onlySonPointHash;
  if (  ( ref ( $DadSonRelationHash ) ) eq 'HASH' ){
    foreach my $DadKey (    sort {$a cmp $b} (   keys (  %{ $DadSonRelationHash }  )   )    ){  
    	if (  ( ref ( $DadSonRelationHash->{$DadKey} ) ) eq 'HASH' ){
    		$onlyDadPointHash->{$DadKey}=1;
    		
        foreach my $SonKey (    sort {$a cmp $b} (   keys (  %{ $DadSonRelationHash->{$DadKey} }  )   )    ){ 
        	$onlySonPointHash->{$SonKey}=1;                   
        }
      }    	
    }
  }
  
  my $rootPointHASH;
  if (  ( ref ( $onlyDadPointHash ) ) eq 'HASH' ){
    foreach my $DadKey (    sort {$a cmp $b} (   keys (  %{ $onlyDadPointHash }  )   )    ){  
    	if (   (  defined ( $onlySonPointHash->{$DadKey} )  ) && ( $onlySonPointHash->{$DadKey}==1 )   ){
    		
    	}
    	else {
    		$rootPointHASH->{$DadKey}=1;
    	}
    }
  }
  
  my $finalTree;
  #$finalTree=&DiGuiDadSon($finalTree, $onlyDadPointHash, $DadSonRelationHash);
  if (  ( ref ( $rootPointHASH ) ) eq 'HASH' ){
    foreach my $DadKey (    sort {$a cmp $b} (   keys (  %{ $rootPointHASH }  )   )    ){   print "\$DadKey=$DadKey\n\n";
    	my $eachTree;
    	$eachTree=&DiGuiDadSon($eachTree, $DadKey, $DadSonRelationHash,1);
    	$finalTree->{$DadKey}=$eachTree->{$DadKey};
    }
  }
  
  return $finalTree;
}

sub DiGuiDadSon{  #利用递归法 建立树
  my ($buildTree,  $workPoint, $DadSonTree, $nb)=@_;   #print "Insub DiGuiDadSon: \$buildTree=$buildTree,  \$workPoint=$workPoint, \$DadSonTree=$DadSonTree\n";
                                                #&PrintTreePoints($buildTree);
  
  if (  ( ref ( $DadSonTree->{$workPoint} ) ) eq 'HASH' ){
  	my $forNB=0;
    foreach my $SonKey (    sort {$a cmp $b} (   keys (  %{ $DadSonTree->{$workPoint} }  )   )    ){   #print "\$buildTree->{$workPoint}=$buildTree->{$workPoint} \n"; 
    	print "\$nb=$nb\t \$forNB=$forNB\t\$SonKey=$SonKey \$DadSonTree->{$workPoint}=$DadSonTree->{$workPoint}\n\n";
      $forNB++;  	
      $buildTree->{$workPoint}= &DiGuiDadSon($buildTree->{$workPoint}, $SonKey, $DadSonTree, $nb++);         	
        	
    }
  }
  
  else {
    $buildTree->{$workPoint}={'-!End-+-Point!-'};
  }
  
  return $buildTree;
  
}

sub OLD_DiGuiDadSon{  #利用递归法 建立树
  my ($buildTree,  $workPoint, $DadSonTree, $nb)=@_;   #print "Insub DiGuiDadSon: \$buildTree=$buildTree,  \$workPoint=$workPoint, \$DadSonTree=$DadSonTree\n";
                                                #&PrintTreePoints($buildTree);
  
  if (  ( ref ( $DadSonTree->{$workPoint} ) ) eq 'HASH' ){
    foreach my $SonKey (    sort {$a cmp $b} (   keys (  %{ $DadSonTree->{$workPoint} }  )   )    ){   print "\$buildTree->{$workPoint}=$buildTree->{$workPoint} \n"; 
    	print "\$nb=$nb\t \$SonKey=$SonKey \$DadSonTree=$DadSonTree\n\n";
        	
      $buildTree->{$workPoint}= &DiGuiDadSon($buildTree->{$workPoint}, $SonKey, $DadSonTree, $nb++);         	
        	
    }
  }
  
  else {
    $buildTree->{$workPoint}={'-!End-+-Point!-'};
  }
  
  return $buildTree;
  
}

sub PrintTreePoints{  #打印出树形结构
  my ($OutLine, $InHash, $domDefHash, $levNb)=@_;
  if (defined ($levNb)){$levNb++;}else {$levNb=0;print "\n";}
  
  if (  ( ref ( $InHash ) ) eq 'HASH' ){
    foreach my $ecKey (    sort {$a cmp $b} (   keys (  %{ $InHash }  )   )    ){ 
      if ($ecKey eq '-!End-+-Point!-'){}
      else { 
      	for (my $i=0; $i<$levNb; $i++){ print "    ";$OutLine.="    "; } 
      	print     "+--- $ecKey:$domDefHash->{$ecKey}->{'Definition' }\n";  
      	$OutLine.="+--- $ecKey:$domDefHash->{$ecKey}->{'Definition' }\n";
      }   
      $OutLine=&PrintTreePoints( $OutLine, $InHash->{$ecKey}, $domDefHash, $levNb );
      
    }
  }
  else {
    #for (my $i=0; $i<$levNb; $i++){ print " "; } #print "\n"; #print "\$levNb=$levNb \$InHash=$InHash\n";
  }
  return $OutLine;
}

sub treePanelWork{  #画出树的png图
  my ($panelHERE, $domTreHsh, $DomIDHash, $alnlength, $levNb)=@_;
  if (defined ($levNb)){$levNb++;}
  else {$levNb=0;print "\n";
  	$panelHERE->add_track (
                             -glyph  => 'line',
                          
                          );
  }
  
  if (  ( ref ( $domTreHsh ) ) eq 'HASH' ){
    foreach my $ecKey (    sort {  $DomIDHash->{$b}->{'times'} <=> $DomIDHash->{$a}->{'times'}  } (   keys (  %{ $domTreHsh }  )   )    ){ 
      if ($ecKey eq '-!End-+-Point!-'){}
      else { 
      	my $startLine=''; for (my $i=0; $i<$levNb; $i++){ print " ";$startLine.=" "; } 
      	$panelHERE->add_track ( generic => Bio::SeqFeature::Generic->new(-start=> $alnlength*75/100,
                                                                         -end  => $alnlength*95/100,
                                                                      
                                                                      ),
                            -glyph  => 'generic',
                            -bgcolor=> $DomIDHash->{$ecKey}->{'Domaincolor'},
                            -key    => "$startLine+---- $ecKey:$DomIDHash->{$ecKey}->{'Definition' }",
                            -label  => 1,
                            
                            
                      );
      	
        $panelHERE=&treePanelWork( $panelHERE, $domTreHsh->{$ecKey}, $DomIDHash, $alnlength, $levNb );
      }
      
    }
  }
  return $panelHERE;
}



sub MergeImage{  #将png文件$image_1和$image_2,合并成$image_out, 应该是上下合并
  my ($image_1, $image_2, $image_out)=@_;
  
  #warn "\$image_1=$image_1\t\$image_2=$image_2\t\$image_out=$image_out\n";
  
  my $fx=0; my $fy=0; my $dist=10;  my $eOn;
  my $im1 ; my( $x1,  $y1  ); if ( -e($image_1) ){ $eOn->{1 }=1;$im1 = GD::Image->newFromPng( $image_1 ) or die"die \$image_1=$image_1 $!\n"; 
  	( $x1,  $y1  ) = $im1->getBounds() ; if ($fx<=$x1 ){$fx=$x1 ;}$fy+=($y1 +$dist); }
  my $im2 ; my( $x2,  $y2  ); if ( -e($image_2) ){ $eOn->{2 }=1;$im2 = GD::Image->newFromPng( $image_2 ) or die"$!\n"; 
  	( $x2,  $y2  ) = $im2->getBounds() ; if ($fx<=$x2 ){$fx=$x2 ;}$fy+=($y2 +$dist); }
  #warn "\$fx=$fx; \$fy=$fy\n";
  
  my $imF = GD::Image->new( $fx, $fy, 0 );       my $startY=0;         my $dist0=1;               
  if (  ( defined($eOn->{1 }) ) && ($eOn->{1 }==1)  ){ $imF->copy( $im1 , 0,$startY, 0, 0, $x1 , $y1 );  $startY+=($y1 +$dist0);}
  if (  ( defined($eOn->{2 }) ) && ($eOn->{2 }==1)  ){ $imF->copy( $im2 , 0,$startY, 0, 0, $x2 , $y2 );  $startY+=($y2 +$dist);}
  #warn "\$startY=$startY\n";
                            
  open ALLPNG, ">$image_out" or die "cannot create \$image_out=$image_out\n";
  print ALLPNG $imF->png;
  close ALLPNG;
                          
}
  





sub AddTracks2Panel{  #将各个含有feature的Track加到Panel中
  my ($tpPanel,  $keyHere, $tkBgCol)=@_;  #warn "\$tpPanel=$tpPanel, \$tkBgCol=$tkBgCol, \$keyHere=$keyHere\n";  print   "\$tpPanel=$tpPanel,  \$tkBgCol=$tkBgCol, \$keyHere=$keyHere\n";
    
  #先建立各个分段的feature
  #$e1 = Bio::SeqFeature::Lite->new(-start=>1,-stop=>100,-type=>'exon');
  
  
  my $tpTrack=$tpPanel->add_track(   

                                           -glyph       => 'segments',
                                           -connector   => 'dashed',
                                           #-tkcolor     =>  $tkBgCol, #'whitesmoke', #$colorIdx, #'#FFF8DC',#'lightyellow',#'cyan',#'gray'
                                           -bgcolor     =>  #$colorHere,
                                                            sub {  #'green',
                                                                  my $feature = shift;  
                                                                  #warn "InFeature :", $feature->start,"..",$feature->end,"\n"; 
                                                                  my $returnColor= 'green'; 
                                                                  my $segType=$feature->seqname;  #warn "\$segType=$segType\n\$colorHere=$colorHere\n\n";
                                                                  
                                                                  #if    ($colorHere=~m/\S+/  )  {$returnColor=$colorHere;}
                                                                  if    ($segType eq '_INCDS')  {$returnColor='lightgreen';}
                                                                  elsif ($segType eq '_BtCDS')  {$returnColor='blue';}
                                                                  elsif ($segType eq '_UinCDS') {$returnColor='red';}
                                                                  elsif ($segType eq '_UbtCDS') {$returnColor='orange';}
                                                                  else  {$returnColor=$segType;}
                                                                  #warn "\$segType=$segType, \$returnColor=$returnColor \n";
                                                                  $returnColor;
                                                                },  
                                           
                                           
                                           
                                           -fgcolor     => 'black',
                                           -font2color  => 'red',
                                           -key         => $keyHere,
                                           -spacing     =>  0,
                                           -bump        =>  0,
                                           -height      =>  10,
                                           -label       =>  0,  
                                           #-part_labels =>  1,             #显示外显子编号
                                           -description =>  1
                                                          #sub {
                                                          #           my $feature = shift;
                                                          #           return unless $feature->has_tag('description');
                                                          #           my ($description) = $feature->each_tag_value('description');
                                                          #           "\$description=$description";
                                                          #      },
                              
                     );
  #print "33333333333333333",$tpTrack->box, "\n\n";
  
  
 
  
  return [$tpPanel,$tpTrack];                 
  


}

sub AddFeature2PanelAndTrack{ #将各个含有feature的Track加到Panel中
  my ($AAposList, $tpPanel, $tpTrack,  $colorHere)=@_;  #warn "\$AAposList=$AAposList, \$tpPanel=$tpPanel, \$colorHere=$colorHere, \$tpTrack=$tpTrack\n"; 
    
  #先建立各个分段的feature
  #$e1 = Bio::SeqFeature::Lite->new(-start=>1,-stop=>100,-type=>'exon');
  my $allSubFeatures;
  for (my $i=0; $i< @{ $AAposList }; $i++){  
    $allSubFeatures->[$i]=Bio::SeqFeature::Lite->new( -start=>$AAposList->[$i]->{'_tailPos'},-stop=>$AAposList->[$i]->{'_headPos'}, -name => $colorHere); #$AAposList->[$i]->{'_SegmentType'} );
  }
  
  #然后建立 总Feature
  my $AlignFeature = Bio::SeqFeature::Lite->new (-segments=>$allSubFeatures,  -strand       => 1);
  #然后画各个比对结果
  
  

  
  $tpTrack->add_feature($AlignFeature,  -bgcolor     =>$colorHere); 
  #$tpTrack->add_feature($AlignFeature_2);
  
  
  return [$tpPanel,$tpTrack];                 
  
  
    


}


sub GetSeqPosList{  #利用clustalw的结果 及column_from_residue_number 等方法，获得经过比对对齐后，重新排布的各个片段位置形成的数列
  my ($inSeqOgj, $inSeqStart, $inSeqEnd)=@_;
  my $outList;  my $ListIdx=0;
  for (my $i=$inSeqStart; $i<=$inSeqEnd; $i++){
    # 判断任何一个字符是不是 开头和结尾，需要进行四个步骤的判断：如开头判断：1看是不是 整个序列开头， 2 看是不是 单个片段开头， 3 看是不是特定字符 ，如U，   4 看是不是 特定字符后的 片段开头
    #首先判断开头
    #1先判断这个字符是不是在一个片段中是开头，就是让这个字符和前一个字符的位置数相减，看是否大于1，单由于如果是第一个片段，则不存在前一个字符，所以，先判断是否是第一个字符。
    if ($i==$inSeqStart){                                                                                                print "\n\n1.0 1.0 If  \$i==$i==\$inSeqStart\n";
    	$outList->[$ListIdx]->[0]=$inSeqOgj->column_from_residue_number($i);                                               print "\t\t\$outList->[\$ListIdx]->[0]=\$outList->[$ListIdx]->[0]=\$inSeqOgj->column_from_residue_number(\$i)=\$inSeqOgj->column_from_residue_number($i)= ", $inSeqOgj->column_from_residue_number($i), "\n";
    } #2经上步判断，如不是第一个片段，则用elsif 来判断，该字符是不是 片段头的位置
    elsif ( ($inSeqOgj->column_from_residue_number($i)  -$inSeqOgj->column_from_residue_number($i-1)) > 1 ){             print "1.1 1.1 If    (\$inSeqOgj->column_from_residue_number(\$i)  -\$inSeqOgj->column_from_residue_number(\$i-1)) = (\$inSeqOgj->column_from_residue_number($i)  -\$inSeqOgj->column_from_residue_number($i-1))= (",$inSeqOgj->column_from_residue_number($i),"-",$inSeqOgj->column_from_residue_number($i-1),") =",($inSeqOgj->column_from_residue_number($i) - $inSeqOgj->column_from_residue_number($i-1))," > 1\n" ;
     	$outList->[$ListIdx]->[0]=$inSeqOgj->column_from_residue_number($i);                                               print "\t\t\$outList->[\$ListIdx]->[0]=\$outList->[$ListIdx]->[0]=\$inSeqOgj->column_from_residue_number(\$i)=\$inSeqOgj->column_from_residue_number($i)= ",$inSeqOgj->column_from_residue_number($i),"\n";
    } #3 假如该位置为U，则设置该位置为 片段开头
    elsif ( $inSeqOgj->subseq($inSeqOgj->column_from_residue_number($i),$inSeqOgj->column_from_residue_number($i)) =~m /^U$/i ){  print "1.2 1.2  该位置为U \n";
    	$outList->[$ListIdx]->[0]=$inSeqOgj->column_from_residue_number($i);                                               print "\t\t\$outList->[\$ListIdx]->[0]=\$outList->[$ListIdx]->[0]=\$inSeqOgj->column_from_residue_number(\$i)=\$inSeqOgj->column_from_residue_number($i)= ",$inSeqOgj->column_from_residue_number($i),"\n";  
    } #4 假如该位置的上一个字符为U，则设置该位置为 片段开头
    elsif ( $inSeqOgj->subseq($inSeqOgj->column_from_residue_number($i)-1,$inSeqOgj->column_from_residue_number($i)-1) =~m /^U$/i ){  print "1.3 1.3  该位置的上个位置为U \n";
    	$outList->[$ListIdx]->[0]=$inSeqOgj->column_from_residue_number($i);                                                print "\t\t\$outList->[\$ListIdx]->[0]=\$outList->[$ListIdx]->[0]=\$inSeqOgj->column_from_residue_number(\$i)=\$inSeqOgj->column_from_residue_number($i)= ",$inSeqOgj->column_from_residue_number($i),"\n";  
    }
      
    #然后，再判断该字符是不是出现在一个片段的结尾，因为有的片段有可能由一个字符所组成，所以这里的判断和上面的判断头字符的步骤 应该是平行的，所以不用elseif而用if
    #1先判断是否是最后一个片段的尾部
    if ($i==$inSeqEnd){                                                                                                  print "\n2.0 2.0 Elsif \$i==$i==\$inSeqEnd\n";
    	$outList->[$ListIdx]->[1]=$inSeqOgj->column_from_residue_number($i);                                               print "\t\t\$outList->[\$ListIdx]->[0]=\$outList->[$ListIdx]->[0]=\$inSeqOgj->column_from_residue_number(\$i)=\$inSeqOgj->column_from_residue_number($i)= ", $inSeqOgj->column_from_residue_number($i), "\n";
    	$outList->[$ListIdx]->[2]=$inSeqOgj->subseq($outList->[$ListIdx]->[0],$outList->[$ListIdx]->[1]);                  warn "222222  \$outList->[$ListIdx]->[2]=$outList->[$ListIdx]->[2]\n";
    	$ListIdx++;
    }  #2下面再判断是不是 中部片段的尾部
    elsif ( ($inSeqOgj->column_from_residue_number($i+1)-$inSeqOgj->column_from_residue_number($i)  ) > 1 ){             print "2.1 2.1 ElsIf (\$inSeqOgj->column_from_residue_number(\$i+1)-\$inSeqOgj->column_from_residue_number(\$i)  ) = (\$inSeqOgj->column_from_residue_number($i+1)-\$inSeqOgj->column_from_residue_number($i)  )= (",$inSeqOgj->column_from_residue_number($i+1),"-",$inSeqOgj->column_from_residue_number($i), ")=" ,($inSeqOgj->column_from_residue_number($i+1)- $inSeqOgj->column_from_residue_number($i)  )," > 1\n" ;
      $outList->[$ListIdx]->[1]=$inSeqOgj->column_from_residue_number($i);                                               print "\t\t\$outList->[\$ListIdx]->[1]=\$outList->[$ListIdx]->[1]=\$inSeqOgj->column_from_residue_number(\$i)=\$inSeqOgj->column_from_residue_number($i)= ",$inSeqOgj->column_from_residue_number($i),"\n";
      $outList->[$ListIdx]->[2]=$inSeqOgj->subseq($outList->[$ListIdx]->[0],$outList->[$ListIdx]->[1]);                  warn "2.1 2.1 \$outList->[$ListIdx]->[2]=$outList->[$ListIdx]->[2]\n";
      $ListIdx++;
    }  #3 假如该位置为U，则设置该位置为 片段结尾
    elsif ( $inSeqOgj->subseq($inSeqOgj->column_from_residue_number($i),$inSeqOgj->column_from_residue_number($i)) =~m /^U$/i ){  print "2.2 2.2  该位置为U \n";
    	$outList->[$ListIdx]->[1]=$inSeqOgj->column_from_residue_number($i);                                               print "\t\t\$outList->[\$ListIdx]->[1]=\$outList->[$ListIdx]->[1]=\$inSeqOgj->column_from_residue_number(\$i)=\$inSeqOgj->column_from_residue_number($i)= ",$inSeqOgj->column_from_residue_number($i),"\n";
      $outList->[$ListIdx]->[2]=$inSeqOgj->subseq($outList->[$ListIdx]->[0],$outList->[$ListIdx]->[1]);                  warn "2.2 2.2 \$outList->[$ListIdx]->[2]=$outList->[$ListIdx]->[2]\n";
      $ListIdx++;
    } #4 假如该位置的下一个字符为U，则设置该位置为 片段结尾
    elsif ( $inSeqOgj->subseq($inSeqOgj->column_from_residue_number($i)+1,$inSeqOgj->column_from_residue_number($i)+1) =~m /^U$/i ){  print "2.3 2.3  该位置的上个位置为U \n";
    	$outList->[$ListIdx]->[1]=$inSeqOgj->column_from_residue_number($i);                                               print "\t\t\$outList->[\$ListIdx]->[1]=\$outList->[$ListIdx]->[1]=\$inSeqOgj->column_from_residue_number(\$i)=\$inSeqOgj->column_from_residue_number($i)= ",$inSeqOgj->column_from_residue_number($i),"\n";
      $outList->[$ListIdx]->[2]=$inSeqOgj->subseq($outList->[$ListIdx]->[0],$outList->[$ListIdx]->[1]);                  warn "2.3 2.3 \$outList->[$ListIdx]->[2]=$outList->[$ListIdx]->[2]\n";
      $ListIdx++;
    }
      
    
  }
  return $outList;
}


sub GetSeqPosList_2 {  #比GetSeqPosList多了个Type的信息输入 #利用clustalw的结果 及column_from_residue_number 等方法，获得经过比对对齐后，重新排布的各个片段位置形成的数列
  my ($inSeqOgj, $inSeqStart, $inSeqEnd, $inSegType)=@_;  #$inSegType 这个参数是 输入的seq 所属的类型，通常是指 特定氨基酸 是否是U，是否是外显子内 还是 外显子末端等
  
  my %keysUsedHere=(
    '_headPos'        =>   "这个 Segment 氨基酸序列的 起点",
    '_tailPos'        =>   "这个 Segment 氨基酸序列的 终点",
    '_SegmentType'    =>   "某个氨基酸AA位置的 类型，或者某段氨基酸序列Segment位置数标明区段的 类型，具体类型，就是下面这4中类型"
  );
  
  my $outList;  my $ListIdx=0;
  for (my $i=$inSeqStart; $i<=$inSeqEnd; $i++){   
    # 判断任何一个字符是不是 开头和结尾，需要进行四个步骤的判断：如开头判断：1看是不是 整个序列开头， 2 看是不是 单个片段开头， 3 看是不是特定字符 ，如U，   4 看是不是 特定字符后的 片段开头
    #首先判断开头
    #1先判断这个字符是不是在一个片段中是开头，就是让这个字符和前一个字符的位置数相减，看是否大于1，单由于如果是第一个片段，则不存在前一个字符，所以，先判断是否是第一个字符。
    if ($i==$inSeqStart){                                                                                                         print "\n\n1.0 1.0 If  \$i==$i==\$inSeqStart\n";
    	$outList->[$ListIdx]->{'_headPos'}=$inSeqOgj->column_from_residue_number($i);                                               print "\t\t\$outList->[\$ListIdx]->{'_headPos'}=\$outList->[$ListIdx]->{'_headPos'}=\$inSeqOgj->column_from_residue_number(\$i)=\$inSeqOgj->column_from_residue_number($i)= ", $inSeqOgj->column_from_residue_number($i), "\n";
    } #2经上步判断，如不是第一个片段，则用elsif 来判断，该字符是不是 片段头的位置
    elsif ( ($inSeqOgj->column_from_residue_number($i)  -$inSeqOgj->column_from_residue_number($i-1)) > 1 ){                      print "1.1 1.1 If    (\$inSeqOgj->column_from_residue_number(\$i)  -\$inSeqOgj->column_from_residue_number(\$i-1)) = (\$inSeqOgj->column_from_residue_number($i)  -\$inSeqOgj->column_from_residue_number($i-1))= (",$inSeqOgj->column_from_residue_number($i),"-",$inSeqOgj->column_from_residue_number($i-1),") =",($inSeqOgj->column_from_residue_number($i) - $inSeqOgj->column_from_residue_number($i-1))," > 1\n" ;
     	$outList->[$ListIdx]->{'_headPos'}=$inSeqOgj->column_from_residue_number($i);                                               print "\t\t\$outList->[\$ListIdx]->{'_headPos'}=\$outList->[$ListIdx]->{'_headPos'}=\$inSeqOgj->column_from_residue_number(\$i)=\$inSeqOgj->column_from_residue_number($i)= ",$inSeqOgj->column_from_residue_number($i),"\n";
    }
      
    #然后，再判断该字符是不是出现在一个片段的结尾，因为有的片段有可能由一个字符所组成，所以这里的判断和上面的判断头字符的步骤 应该是平行的，所以不用elseif而用if
    #1先判断是否是最后一个片段的尾部
    if ($i==$inSeqEnd){                                                                                                           print "\n2.0 2.0 Elsif \$i==$i==\$inSeqEnd\n";
    	$outList->[$ListIdx]->{'_tailPos'}=$inSeqOgj->column_from_residue_number($i);                                               print "\t\t\$outList->[\$ListIdx]->{'_headPos'}=\$outList->[$ListIdx]->{'_headPos'}=\$inSeqOgj->column_from_residue_number(\$i)=\$inSeqOgj->column_from_residue_number($i)= ", $inSeqOgj->column_from_residue_number($i), "\n";
    	$outList->[$ListIdx]->{'_SegmentType'}=$inSegType;                                                                          #warn "222222  \$outList->[$ListIdx]->{'_SegmentType'}=$outList->[$ListIdx]->{'_SegmentType'}\n";
    	$ListIdx++;
    }  #2下面再判断是不是 中部片段的尾部
    elsif ( ($inSeqOgj->column_from_residue_number($i+1)-$inSeqOgj->column_from_residue_number($i)  ) > 1 ){                      print "2.1 2.1 ElsIf (\$inSeqOgj->column_from_residue_number(\$i+1)-\$inSeqOgj->column_from_residue_number(\$i)  ) = (\$inSeqOgj->column_from_residue_number($i+1)-\$inSeqOgj->column_from_residue_number($i)  )= (",$inSeqOgj->column_from_residue_number($i+1),"-",$inSeqOgj->column_from_residue_number($i), ")=" ,($inSeqOgj->column_from_residue_number($i+1)- $inSeqOgj->column_from_residue_number($i)  )," > 1\n" ;
      $outList->[$ListIdx]->{'_tailPos'}=$inSeqOgj->column_from_residue_number($i);                                               print "\t\t\$outList->[\$ListIdx]->{'_tailPos'}=\$outList->[$ListIdx]->{'_tailPos'}=\$inSeqOgj->column_from_residue_number(\$i)=\$inSeqOgj->column_from_residue_number($i)= ",$inSeqOgj->column_from_residue_number($i),"\n";
      $outList->[$ListIdx]->{'_SegmentType'}=$inSegType;                                                                          #warn "2.1 2.1 \$outList->[$ListIdx]->{'_SegmentType'}=$outList->[$ListIdx]->{'_SegmentType'}\n";
      $ListIdx++;
    }  
      
    
  }
  return $outList;
}



sub GetSeqPosList_3 {  #这个就把整个的domain画出来，中间被隔开的连接起来处理
  my ($inSeqOgj, $inSeqStart, $inSeqEnd, $inSegType)=@_;  #$inSegType 这个参数是 输入的seq 所属的类型，通常是指 特定氨基酸 是否是U，是否是外显子内 还是 外显子末端等
  
  my %keysUsedHere=(
    '_headPos'        =>   "这个 Segment 氨基酸序列的 起点",
    '_tailPos'        =>   "这个 Segment 氨基酸序列的 终点",
    '_SegmentType'    =>   "某个氨基酸AA位置的 类型，或者某段氨基酸序列Segment位置数标明区段的 类型，具体类型，就是下面这4中类型"
  );
  
  my $outList;  
  $outList->[0]->{'_headPos'}    =$inSeqOgj->column_from_residue_number($inSeqStart);
  $outList->[0]->{'_tailPos'}    =$inSeqOgj->column_from_residue_number($inSeqEnd  ); 
  $outList->[0]->{'_SegmentType'}=$inSegType;  
  
  return $outList;
}





sub DrawClustalw{  #画出clustlaw结果的 png图 #这个函数在这个程序中没有被使用到
	my ($clustalwOUTfile, $formatHere, $clustalPNGfile)=@_;

  
  #读入clustalw输出文件
  my $in  = Bio::AlignIO->new(-file   => $clustalwOUTfile ,
                              -format => $formatHere            
                              #-report_type => 'blastn'    
                              );

  open (OUT,">$clustalPNGfile") or die "cannot create $clustalPNGfile : $!";

  #my $testNumber=130;

  #读取clustalw 各个比对 的数据流，一个数据流 就是一个aln
  my $howManySeq=0;
  my $aln=$in->next_aln();
  #while (my $aln=$in->next_aln)
  if (1){
    warn "1: ", $aln->length, "\n";   warn "2: ", $aln->num_residues, "\n";  warn "3: ", $aln->is_flush, "\n";  warn "4: ", $aln->num_sequences, "\n";  warn "5: ", $aln->percentage_identity, "\n";  warn "6: ", $aln->consensus_string(50), "\n";
    
    #首先画 panel
    my $panel = Bio::Graphics::Panel->new(
                                               -length    => $aln->length,   #图的全长宽度能表示的碱基数量总长度
                                               -width     => 1000,                #图的实际宽度
                                               -key_style => 'between',
                                               -start     => 1,
                                               -pad_left  => 100,
                                               -pad_right => 100,
                                               
                                         );  
    
    #然后画 标线
    my $full_length = Bio::SeqFeature::Generic->new(
                                                         -start => 1, 
                                                         -end   => $aln->length
                                                        );
    # 将标尺的 track加入到panel中。
    $panel->add_track ( $full_length,
                                           -glyph   => 'arrow',
                                           -tick    => 2,
                                           -fgcolor => 'black',
                                           -double  => 1
                      );
                      
    #下面画每个 对齐对象的 track
    
    foreach my $seq ($aln->each_seq) {          warn "id: ",$seq->id(),"\n";     warn "Seq: ",$seq->seq(),"\n";    warn "Start: ",$seq->start(),"\n";    warn "End: ",$seq->end(),"\n";    warn "Number of gaps 1 .: ", $seq->num_gaps('.'),"\n";    warn "Number of gaps 2 -: ", $seq->num_gaps('-'),"\n";    warn "Number of gaps 2  : ", $seq->num_gaps(),"\n";            #warn "Col Number of rsd 20  : ", $seq->column_from_residue_number(20),"\n";    #warn "location_from_column 82  : ", $seq->location_from_column(82),"\n";       #warn "\n";    #my @tempMap=$seq->mapping();    #my $idxNb=0;foreach my $temp(@tempMap){warn "$idxNb: \$temp=$temp\t\t";$idxNb++}    #warn "\nAfter Print \@tempMap\n\n";        #my %hashTp=$seq->frameshifts();    #foreach my $keyH (keys (%hashTp)){warn "\$hashTp{\$keyH}=\$hashTp{$keyH}=$hashTp{$keyH}\t\t";} warn "\nAfter Print \%hashTp\n\n\n\n";        #my $loc82loc=$seq->location_from_column(82);    #if (defined ($loc82loc)){ warn "Loc82 LocationType: ",$loc82loc->location_type," Start: ",$loc82loc->start,"\tEnd: ",$loc82loc->end,"\nAfter Print \$loc82loc\n\n\n\n\n";}    
      
      my $outAAposList=&GetSeqPosList($seq,$seq->start(),$seq->end());        warn "Print AA seq List: ", &PrintList( $outAAposList ), "\n\n\n\n\n";
      my $segName=$seq->id(); my $segLength=$seq->end();        
      #然后画各个比对结果
      $panel=&NEWAddFeatureTrack2Panel($outAAposList,$panel,'','green',"$segName  $segLength");
      $howManySeq++;
    }
    print OUT  $panel->png; 
  }
  close (OUT);  
     
  return $howManySeq;
	
	
}
sub PrintList{  #这个函数在这个程序中没有被使用到 #获得用于打印出来的 外显子位点的 数列如： 116..156,158..181,184..193,195..202,206..263,265..313,315..330
  my ($inList)=@_;
  my @outList;
  for (my $j=0; $j<@{ $inList }; $j++){
    push @outList, "$inList->[$j]->[0]..$inList->[$j]->[1]";
  }
  my $outListLine=join (",", @outList);
  return $outListLine;
}



sub NEWAddFeatureTrack2Panel{  #给Panel加Track
  my ($AAposList, $tpPanel, $idKeyHere, $colorHere, $keyHere)=@_;
    
  #先建立各个分段的feature
  #$e1 = Bio::SeqFeature::Lite->new(-start=>1,-stop=>100,-type=>'exon');
  my $allSubFeatures;
  for (my $i=0; $i< @{ $AAposList }; $i++){  
    $allSubFeatures->[$i]=Bio::SeqFeature::Lite->new( -start=>$AAposList->[$i]->[0],-stop=>$AAposList->[$i]->[1],-seq  =>  $AAposList->[$i]->[2], -id => "onYeah $i" );
  }
  
  #然后建立 总Feature
  my $AlignFeature = Bio::SeqFeature::Lite->new (-segments=>$allSubFeatures, -id =>$idKeyHere, -strand       => 1);
  #然后画各个比对结果
  
  $tpPanel->add_track(   $AlignFeature ,

                                           -glyph       => 'segments',
                                           -connector   => 'dashed',
                                           -bgcolor     =>  #$colorHere,
                                                            sub {  #'green',
                                                                  my $feature = shift;  &print_all_sub_array($feature);
                                                                  warn "InFeature :", $feature->start,"..",$feature->end,"\n"; 
                                                                  my $returnColor= 'green'; 
                                                                  my $sementSeq = $feature->{'seq'}; warn "\$sementSeq=",$sementSeq,"\n";
                                                                  my $idHere=$feature->{'id'};  warn "\$idHere=",$idHere,"\n";
                                                                  if ($sementSeq =~m/^U$/i  ){ $returnColor='red';}
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




###################################
#

#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#









1;