#!/usr/bin/perl -w
BEGIN{  push (@INC, '/home/fredjiang/fredProgs/PMs')  };


use strict;
use warnings;



############################################################################################################
#
#
#          sub clustalwForLongPathFile{  #��Ϊclusatlw�����е�ʱ���������ļ���·���������ͻ�������Ա༭���������
#          sub clustalwForLongPathFileOutOrgOrder{  #��Ϊclusatlw�����е�ʱ���������ļ���·���������ͻ�������Ա༭���������,��clustalwForLongPathFile���������ڣ�����ĳ������������Fasta������˳�������е�
#          sub readInSeq{ #��ȡһ�������е��ļ��������ļ����hash��hash����������������ÿ�����и������Ը�����primary_idΪKey�ĸ�����Ϣ
#          sub readInSeqAndUseClustalwOrder{  #��ȡһ�������е��ļ��������ļ����hash��hash����������������ÿ�����и������Ը�����primary_idΪKey�ĸ�����Ϣ,�������������һ�����ܣ����ø�ʽ��$MsfFormat��alignemnt�ļ�$orgMsfFile�е�˳�����Զ������ļ���������
#          sub runingClustalw{  #ֱ������Clustalw
#          sub runingClustalwInputOrder{  #ֱ������Clustalw,���ȶԺ� ���˳����Ȼ��Ҫ����� fasta���е�˳�����
#          sub GetClustalwOrderHash{   #����һ��algnment�ļ������һ��Hash������������ID��˳��֮��� �����Ӧ��ϵ��hash
#          sub ReWriteClustalw{  #��msf��ʽ��clusatlw����ļ�$orgMsfFile�����и�ʽת����ת��Ϊ$outFormat��ʽ��$outMutiAlgnFile�ļ�
#          sub RewriteClustalw20180408{     #ĳЩ���룬��Ҫ���� 2018.04.08.genomePostionAnalysis.pl �е��м����ݽṹ #����һ�� �����бȶԵĽ��������б�ŵ��½����ͬʱ���mega�����������һ��hash�����м�¼�� id�ͱ��֮��Ĺ�ϵ
#          sub sortFastaFile{ #��fastafile��������ʹ�ó����Ƚ�������
#          sub PrintPngForClustalw{  #����clustalw�����msf�ļ�$clustalwOUTfile,��interproscanģ��Ľ���ļ�$Whole_htmlParseOut,����ͼƬ�������$outPng��
#          sub MakeTreeFrom_DadSonHashs {          #���ø��ӹ�ϵ�������νṹ           #($DS_AllDadSonHash, $DS_DomDadSonRelationHash)
#          sub MakeTreeFrom_DadSonHashs_onlyMultiLevesTree {    #���ø��ӹ�ϵ�������νṹ                 #($DS_AllDadSonHash, $DS_DomDadSonRelationHash)
#          sub DiGuiDadSon{  ���õݹ鷨 ������
#          sub PrintTreePoints{  #��ӡ�����νṹ
#          sub treePanelWork{  #��������pngͼ
#          sub MergeImage{  #��png�ļ�$image_1��$image_2,�ϲ���$image_out, Ӧ�������ºϲ�
#          sub AddTracks2Panel{  #����������feature��Track�ӵ�Panel��
#          sub AddFeature2PanelAndTrack{ #����������feature��Track�ӵ�Panel��
#          sub GetSeqPosList{  #����clustalw�Ľ�� ��column_from_residue_number �ȷ�������þ����ȶԶ���������Ų��ĸ���Ƭ��λ���γɵ�����
#          sub GetSeqPosList_2 {  #��GetSeqPosList���˸�Type����Ϣ���� #����clustalw�Ľ�� ��column_from_residue_number �ȷ�������þ����ȶԶ���������Ų��ĸ���Ƭ��λ���γɵ�����
#          sub DrawClustalw{  #����clustlaw����� pngͼ #������������������û�б�ʹ�õ�
#          sub PrintList{  #������������������û�б�ʹ�õ� #������ڴ�ӡ������ ������λ��� �����磺 116..156,158..181,184..193,195..202,206..263,265..313,315..330
#          sub NEWAddFeatureTrack2Panel{  #��Panel��Track
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



sub clustalwForLongPathFile{  #��Ϊclusatlw�����е�ʱ���������ļ���·���������ͻ�������Ա༭���������
  my ($inPathFile, $outPathFile, $outFormat)=@_;
  my $inFile  =File::Basename::basename $inPathFile;
  my $inPath  =File::Basename::dirname ($inPathFile);   
  my $outPath =File::Basename::dirname $outPathFile;
  my $outFile =File::Basename::basename $outPathFile;
  my $InTempFl=TimeWork::GetTimeDirOrFileName()."in.txt";  
  my $otTempFl=TimeWork::GetTimeDirOrFileName()."ot.msf";
  my $otTdndFl=TimeWork::GetTimeDirOrFileName()."in.dnd";   #�����Clustalw�Զ����ɵ�һ���ļ� 
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


sub clustalwForLongPathFileOutOrgOrder{  #��Ϊclusatlw�����е�ʱ���������ļ���·���������ͻ�������Ա༭���������,��clustalwForLongPathFile���������ڣ�����ĳ������������Fasta������˳�������е�
  my ($inPathFile, $outPathFile, $outFormat)=@_;
  my $inFile  =File::Basename::basename $inPathFile;
  my $inPath  =File::Basename::dirname ($inPathFile);   
  my $outPath =File::Basename::dirname $outPathFile;
  my $outFile =File::Basename::basename $outPathFile;
  my $InTempFl=TimeWork::GetTimeDirOrFileName()."in.txt";  
  my $otTempFl=TimeWork::GetTimeDirOrFileName()."ot.msf";
  my $otTdndFl=TimeWork::GetTimeDirOrFileName()."in.dnd";   #�����Clustalw�Զ����ɵ�һ���ļ� 
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
                                -format => $Align_format   );   #������� aln�ļ��������� ����
  
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

sub readInSeq{   #  ClustalwRun::readInSeq($inseqFile, $inFormat);   #��ȡһ�������е��ļ��������ļ����hash��hash����������������ÿ�����и������Ը�����primary_idΪKey�ĸ�����Ϣ
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

sub NewReadInSeq_into_A_HASH{   #  ClustalwRun::NewReadInSeq_into_A_HASH($inseqFile, $inFormat);   #��ȡһ�������е��ļ��������ļ����hash��hash����������������ÿ�����и������Ը�����primary_idΪKey�ĸ�����Ϣ
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
                                -format => $org_format   );   #������� aln�ļ��������� ����
  
  	my $idRelationHash;
  
    my $howManyAln=0;                            
    while ( my $aln = $in->next_aln() ) {
    	
    	#��ʱ�ļ� �����м�fas��ʽ���ļ�
    	my $tempFile=TimeWork::GetTimeDirOrFileName().".fas";
    	
    	#������ʱ�ļ���������� ����
    	my $tempFast    = Bio::SeqIO->new(-file   => ">$tempFile" ,          #
                                        -format => 'fasta');               #��ʱ�ļ��Ĵ򿪣�����д��fasta��ʽ���ļ�, For Msf
    	
    	#ѭ����ȡ msf������Ϣ�� д����ʱ�ļ�  	
    	my $seqIdx=1;
      foreach my $seq ($aln->each_seq) {
        #��ȡ id �� msf������
        my $idHere=$seq->id(); my $seqHere=$seq->seq();                                   # ����õ�ԭ����id�ͱȶԵ�����
        
        
        
        
        my $outIdx=$seqIdx."-".$idHere;
        #For Mega �����е�Uת��ΪC
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
        
         
        #�����µ� SeqObj
        my $sortFasSeqObj=Bio::Seq->new( -display_id => $outIdx,   #$seqIdx,                          #��� 
                                                -seq => $seqHere                ); 
        #print "201811201803-1 \$cutStt=$cutStt\t$cutEnd=$cutEnd\n"; 
        if (  (  defined ( $cutStt )  ) && ( $cutStt=~m/\d+/) && ( $cutStt>0 ) && (  defined ( $cutEnd )  ) && ( $cutEnd=~m/\d+/) && ( $cutEnd>0 )   ){
        	my $truncateSeq=$sortFasSeqObj->trunc($cutStt, $cutEnd);  print "201811201803-2 \$truncateSeq=$truncateSeq\n";
        	#$sortFasSeqObj->seq()=$truncateSeq; 
        	$tempFast->write_seq($truncateSeq);   
        }
        
        else {
        	#д��
          $tempFast->write_seq($sortFasSeqObj);      
        }
                                                
                                                       
            
              
        $seqIdx++;
      }
      
      #������ȡ��ʱ�ļ��� IO����
      my $fastain= Bio::AlignIO->new(-file   => $tempFile ,              ##  #For Msf #��fasta����ʱ�ļ���Ϊ ����
                                     -format => 'fasta'             );   ##
      
      
      #���� д�� ��Msf�����¸�ʽ���ļ��ġ��ɣ϶���
      
      my $out    = Bio::AlignIO->new(-file   => ">$out_AlignFIle" ,     #   #For Msf #�����µ��������
                                     -format => $out_format          );   #   
      
      ##For Msf #������������ �ȶ��ļ��� д��
      while   ( my $FstAln = $fastain->next_aln() ) {
      	$out->write_aln($FstAln);
      }                             
      system ("rm -f $tempFile");    ##ɾ����ʱ�ļ�
      
      
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
                                -format => $org_format   );   #������� aln�ļ��������� ����
  
  	my $idRelationHash;
  
    my $howManyAln=0;                            
    while ( my $aln = $in->next_aln() ) {
    	
    	#��ʱ�ļ� �����м�fas��ʽ���ļ�
    	my $tempFile=TimeWork::GetTimeDirOrFileName().".fas";
    	
    	#������ʱ�ļ���������� ����
    	my $tempFast    = Bio::SeqIO->new(-file   => ">$tempFile" ,          #
                                        -format => 'fasta');               #��ʱ�ļ��Ĵ򿪣�����д��fasta��ʽ���ļ�, For Msf
    	
    	#ѭ����ȡ msf������Ϣ�� д����ʱ�ļ�  	
    	my $seqIdx=1;
      foreach my $seq ($aln->each_seq) {
        #��ȡ id �� msf������
        my $idHere=$seq->id(); my $seqHere=$seq->seq();                                   # ����õ�ԭ����id�ͱȶԵ�����
        warn "1,".$idHere."\n";
        
        
        
        my $outIdx=$seqIdx."-".$idHere;
        #For Mega �����е�Uת��ΪC
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
        
         
        #�����µ� SeqObj
        my $sortFasSeqObj=Bio::Seq->new( -display_id => $idHere,   #$seqIdx,                          #��� 
                                                -seq => $seqHere                ); 
        #print "201811201803-1 \$cutStt=$cutStt\t$cutEnd=$cutEnd\n"; 
        if (  (  defined ( $cutStt )  ) && ( $cutStt=~m/\d+/) && ( $cutStt>0 ) && (  defined ( $cutEnd )  ) && ( $cutEnd=~m/\d+/) && ( $cutEnd>0 )   ){
        	my $truncateSeq=$sortFasSeqObj->trunc($cutStt, $cutEnd);  print "201811201803-2 \$truncateSeq=$truncateSeq\n";
        	#$sortFasSeqObj->seq()=$truncateSeq; 
        	$tempFast->write_seq($truncateSeq);   
        }
        
        else {
        	#д��
          $tempFast->write_seq($sortFasSeqObj);      
        }
                                                
                                                       
            
              
        $seqIdx++;
      }
      
      #������ȡ��ʱ�ļ��� IO����
      my $fastain= Bio::AlignIO->new(-file   => $tempFile ,              ##  #For Msf #��fasta����ʱ�ļ���Ϊ ����
                                     -format => 'fasta'             );   ##
      
      
      #���� д�� ��Msf�����¸�ʽ���ļ��ġ��ɣ϶���
      
      my $out    = Bio::AlignIO->new(-file   => ">$out_AlignFIle" ,     #   #For Msf #�����µ��������
                                     -format => $out_format          );   #   
      
      ##For Msf #������������ �ȶ��ļ��� д��
      while   ( my $FstAln = $fastain->next_aln() ) {
      	$out->write_aln($FstAln);
      }                             
      system ("rm -f $tempFile");    ##ɾ����ʱ�ļ�
      
      
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
                                -format => $org_format   );   #������� aln�ļ��������� ����
  
  	my $idRelationHash;
  
    my $howManyAln=0;                            
    while ( my $aln = $in->next_aln() ) {
    	
    	#��ʱ�ļ� �����м�fas��ʽ���ļ�
    	my $tempFile=TimeWork::GetTimeDirOrFileName().".fas";
    	
    	#������ʱ�ļ���������� ����
    	my $tempFast    = Bio::SeqIO->new(-file   => ">$tempFile" ,          #
                                        -format => 'fasta');               #��ʱ�ļ��Ĵ򿪣�����д��fasta��ʽ���ļ�, For Msf
    	
    	#ѭ����ȡ msf������Ϣ�� д����ʱ�ļ�  	
    	my $seqIdx=1;
      foreach my $seq ($aln->each_seq) {
        #��ȡ id �� msf������
        my $idHere=$seq->id(); my $seqHere=$seq->seq();                                   # ����õ�ԭ����id�ͱȶԵ�����
        warn "1,".$idHere."\n";
        
        
        
        my $outIdx=$seqIdx."-".$idHere;
        #For Mega �����е�Uת��ΪC
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
        
         
        #�����µ� SeqObj
        my $sortFasSeqObj=Bio::Seq->new( -display_id => $idHere,   #$seqIdx,                          #��� 
                                                -seq => $seqHere                ); 
        #print "201811201803-1 \$cutStt=$cutStt\t$cutEnd=$cutEnd\n"; 
        if (  (  defined ( $cutStt )  ) && ( $cutStt=~m/\d+/) && ( $cutStt>0 ) && (  defined ( $cutEnd )  ) && ( $cutEnd=~m/\d+/) && ( $cutEnd>0 )   ){
        	my $truncateSeq=$sortFasSeqObj->trunc($cutStt, $cutEnd);  print "201811201803-2 \$truncateSeq=$truncateSeq\n";
        	#$sortFasSeqObj->seq()=$truncateSeq; 
        	$tempFast->write_seq($truncateSeq);   
        }
        
        else {
        	#д��
          $tempFast->write_seq($sortFasSeqObj);      
        }
                                                
                                                       
            
              
        $seqIdx++;
      }
      
      #������ȡ��ʱ�ļ��� IO����
      my $fastain= Bio::AlignIO->new(-file   => $tempFile ,              ##  #For Msf #��fasta����ʱ�ļ���Ϊ ����
                                     -format => 'fasta'             );   ##
      
      
      #���� д�� ��Msf�����¸�ʽ���ļ��ġ��ɣ϶���
      
      my $out    = Bio::AlignIO->new(-file   => ">$out_AlignFIle" ,     #   #For Msf #�����µ��������
                                     -format => $out_format          );   #   
      
      ##For Msf #������������ �ȶ��ļ��� д��
      while   ( my $FstAln = $fastain->next_aln() ) {
      	$out->write_aln($FstAln);
      }                             
      system ("rm -f $tempFile");    ##ɾ����ʱ�ļ�
      
      
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



sub readInSeqAndUseClustalwOrder{  #��ȡһ�������е��ļ��������ļ����hash��hash����������������ÿ�����и������Ը�����primary_idΪKey�ĸ�����Ϣ,�������������һ�����ܣ����ø�ʽ��$MsfFormat��alignemnt�ļ�$orgMsfFile�е�˳�����Զ������ļ���������
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

sub runingClustalw{  #ֱ������Clustalw
  my ($inseqFile, $outSeqFile)=@_;
  my $cmd="clustalw $inseqFile -OUTPUT=GCG -OUTFILE=$outSeqFile";  
  warn $cmd, "\n\n"; print  $cmd, "\n\n";
  system ("$cmd");
  
}

sub runingClustalwInputOrder{  #ֱ������Clustalw,���ȶԺ� ���˳����Ȼ��Ҫ����� fasta���е�˳�����
  my ($inseqFile, $outSeqFile)=@_;
  my $cmd="clustalw $inseqFile -OUTPUT=GCG -OUTFILE=$outSeqFile -OUTORDER=INPUT";  
  warn $cmd, "\n\n"; print  $cmd, "\n\n";
  system ("$cmd");
  
}


sub GetClustalwOrderHash{   #����һ��algnment�ļ������һ��Hash������������ID��˳��֮��� �����Ӧ��ϵ��hash
  my ($orgMsfFile, $InFormat)=@_;
  
  my $in  = Bio::AlignIO->new(-file   => $orgMsfFile ,     #
                              -format => $InFormat     );   #������� aln�ļ��������� ����
  
  my $idRelationHash;
  
  my $howManyAln=0;                            
  while ( my $aln = $in->next_aln() ) {
  	  	
  	#ѭ����ȡ ������Ϣ 	
  	my $seqIdx=1;
    foreach my $seq ($aln->each_seq) {
      #��ȡ id �� msf������
      my $idHere=$seq->id(); my $seqHere=$seq->seq();                                   # ����õ�ԭ����id�ͱȶԵ�����
      
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

sub clustalBackToFastaFile {  #��$inFormat��ʽ��clusatlw����ļ�$orgMsfFile�����и�ʽת����ת��ΪFasta�ļ�
  my ($orgMsfFile, $inFormat, $outFastaFile )=@_;
  
  my $in  = Bio::AlignIO->new(-file   => $orgMsfFile ,     #
                              -format => $inFormat    );   #������� aln�ļ��������� ����
  
  my $howManyAln=0;                            
  while ( my $aln = $in->next_aln() ) {
  	
  	#������� ����
  	my $tempFast    = Bio::SeqIO->new( -file   => ">$outFastaFile" ,          #
                                       -format => 'fasta'           );               #��ʱ�ļ��Ĵ򿪣�����д��fasta��ʽ���ļ�, For Msf
  	
  	#ѭ����ȡ msf������Ϣ�� д����ʱ�ļ�  	
  	my $seqIdx=1;
    foreach my $seq ($aln->each_seq) {
      #��ȡ id �� msf������
      my $idHere=$seq->id(); my $seqHere=$seq->seq();                                   # ����õ�ԭ����id�ͱȶԵ�����
      $seqHere=~s/\s//g; $seqHere=~s/\n//g; 
      $seqHere=~s/\.//g;  #print "\$seqHere=$seqHere\n";
       
      #�����µ� SeqObj
      my $sortFasSeqObj=Bio::Seq->new( -display_id => $idHere,                          #��� 
                                              -seq => $seqHere                );        
      #д��
      $tempFast->write_seq($sortFasSeqObj);          
            
      $seqIdx++;
    }
    
    $howManyAln++;
  }                         
  if ($howManyAln>1){
    die "\n\$orgMsfFile=$orgMsfFile\t\t\$howManyAln=$howManyAln\n\n\n";
  }                    
  
}


sub clustalBackToFastaHash {  #��$inFormat��ʽ��clusatlw����ļ�$orgMsfFile�����и�ʽת����ת��Ϊ��������Fasta��hash��key��id��ֵ������
  my ($orgMsfFile, $inFormat)=@_;
  
  my $in  = Bio::AlignIO->new(-file   => $orgMsfFile ,     #
                              -format => $inFormat    );   #������� aln�ļ��������� ����
  
  my $outHash;
  
  my $howManyAln=0;                            
  while ( my $aln = $in->next_aln() ) {
  	
  	#ѭ����ȡ msf������Ϣ�� д����ʱ�ļ�  	
  	my $seqIdx=1;
    foreach my $seq ($aln->each_seq) {
      #��ȡ id �� msf������
      my $idHere=$seq->id(); my $seqHere=$seq->seq();                                   # ����õ�ԭ����id�ͱȶԵ�����
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

sub ReWriteClustalw{  #��msf��ʽ��clusatlw����ļ�$orgMsfFile�����и�ʽת����ת��Ϊ$outFormat��ʽ��$outMutiAlgnFile�ļ�
  my ($orgMsfFile, $outMutiAlgnFile, $outFormat)=@_;
  
  my $in  = Bio::AlignIO->new(-file   => $orgMsfFile ,     #
                              -format => 'msf'        );   #������� aln�ļ��������� ����
  
  my $idRelationHash;
  
  my $howManyAln=0;                            
  while ( my $aln = $in->next_aln() ) {
  	
  	#��ʱ�ļ� �����м�fas��ʽ���ļ�
  	my $tempFile=TimeWork::GetTimeDirOrFileName()."fas";
  	
  	#������ʱ�ļ���������� ����
  	my $tempFast    = Bio::SeqIO->new(-file   => ">$tempFile" ,          #
                                      -format => 'fasta');               #��ʱ�ļ��Ĵ򿪣�����д��fasta��ʽ���ļ�, For Msf
  	
  	#ѭ����ȡ msf������Ϣ�� д����ʱ�ļ�  	
  	my $seqIdx=1;
    foreach my $seq ($aln->each_seq) {
      #��ȡ id �� msf������
      my $idHere=$seq->id(); my $seqHere=$seq->seq();                                   # ����õ�ԭ����id�ͱȶԵ�����
      
      
      
      #my $segName=$seq->id(); my $segLength=$seq->end();        
      #my $clustalwSegMentList=&GetSeqPosList_2($seq, 1 , $segLength);
      #PrintDumper('tttt.txt',$clustalwSegMentList );
      
      
      
      
      
      #For Mega �����е�Uת��ΪC
      if ($outFormat eq 'mega'){
        $seqHere=~s/U/C/g;   $seqHere=~s/\./-/g;                  
      } 
      
      
      
      $idRelationHash->{'Org2Chg'}->{$idHere}=$seqIdx;
      $idRelationHash->{'Chg2Org'}->{$seqIdx}=$idHere;
      
       
      #�����µ� SeqObj
      my $sortFasSeqObj=Bio::Seq->new( -display_id => $seqIdx,                          #��� 
                                              -seq => $seqHere                );        
      #д��
      $tempFast->write_seq($sortFasSeqObj);          
            
      $seqIdx++;
    }
    
    #������ȡ��ʱ�ļ��� IO����
    my $fastain= Bio::AlignIO->new(-file   => $tempFile ,              ##  #For Msf #��fasta����ʱ�ļ���Ϊ ����
                                   -format => 'fasta'             );   ##
    
    
    #���� д�� ��Msf�����¸�ʽ���ļ��ġ��ɣ϶���
    
    my $out    = Bio::AlignIO->new(-file   => ">$outMutiAlgnFile" ,     #   #For Msf #�����µ��������
                                   -format => $outFormat          );   #   
    
    ##For Msf #������������ �ȶ��ļ��� д��
    while   ( my $FstAln = $fastain->next_aln() ) {
    	$out->write_aln($FstAln);
    }                             
    system ("rm -f $tempFile");    ##ɾ����ʱ�ļ�
    
    
    $howManyAln++;
  }                         
  if ($howManyAln>1){
    die "\n\$orgMsfFile=$orgMsfFile\t\t\$howManyAln=$howManyAln\n\n\n";
  }                        
  return  $idRelationHash;
  
}
#my $Whole_htmlParseOut=Interproscan::SeqIn_PrasedHashOut($clustalwInSeqFile);


sub RewriteClustalw20180408{     #ĳЩ���룬��Ҫ���� 2018.04.08.genomePostionAnalysis.pl �е��м����ݽṹ #����һ�� �����бȶԵĽ��������б�ŵ��½����ͬʱ���mega�����������һ��hash�����м�¼�� id�ͱ��֮��Ĺ�ϵ
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
   
  my $outPutHash;             #�����hash�����м�¼�� id�ͱ��֮��Ĺ�ϵ
  
  if ($outFormat=~m/\S+/){} else {$outFormat=$inFormat;}   #��������ĸ�����û������Ļ����򲻸ı��ļ��ĸ�ʽ
  
  
  my $in  = Bio::AlignIO->new(-file   => $clutalwInfile ,     #
                              -format => $inFormat        );  #������� aln�ļ��������� ����
  my $timeNow=time();
  my $tempFile     = "${clutalwOutfile}.${timeNow}.temp.txt";                #��1�������м������ fasta��ʽ�� ��Ϊ���Ĵ��ں�ʹ�ã���������������ܽ��в��д��� for Msf
  my $tempFileIdx  = "${clutalwOutfile}.${timeNow}.Idx.temp.txt";            #��2�������м������ fasta��ʽ�� ��Ϊ���Ĵ��ں�ʹ�ã���������������ܽ��в��д��� for Idx Msf
  my $tempFileSort = "${clutalwOutfile}.${timeNow}.sot.temp.txt";            #��3�������м������ fasta��ʽ�� ��Ϊ���Ĵ��ں�ʹ�ã���������������ܽ��в��д��� for sort Msf  �����û��sort ����sort���ļ�
  my $sortedTempFil= "${clutalwOutfile}.${timeNow}.sorted.temp.txt";         #��4�������м������ fasta��ʽ�� ��Ϊ���Ĵ��ں�ʹ�ã���������������ܽ��в��д��� for sort Msf  ������Ѿ�sort���ļ�
  my $tempFileTwo  = "${clutalwOutfile}.${timeNow}.2temp2.txt";              #��5�������м������ fasta��ʽ�� ��Ϊ���Ĵ��ں�ʹ�ã���������������ܽ��в��д��� For Mega

  ################################################################# for Msf#############################################                           
  #���濪ʼ���� aln������������,
  my $idx_to_newId_hash;
  my $testIdx=0;
  while ( my $aln = $in->next_aln() ) {
    
  	my $tempFast    = Bio::SeqIO->new(-file   => ">$tempFile" ,          #
                                      -format => 'fasta');               #��ʱ�ļ��Ĵ򿪣�����д��fasta��ʽ���ļ�, For Msf
                                      
    my $tempFast_Idx= Bio::SeqIO->new(-file   => ">$tempFileIdx" ,          #
                                      -format => 'fasta');               #��ʱ�ļ��Ĵ򿪣�����д��fasta��ʽ���ļ�, For Msf
                                      

  	my $tempSortFast= Bio::SeqIO->new(-file   => ">$tempFileSort" ,         #
                                   -format => 'fasta');                  #��ʱ�ļ��Ĵ򿪣�����д��fasta��ʽ���ļ�, For sort Msf

  	my $tempFastTwo = Bio::SeqIO->new(-file   => ">$tempFileTwo" ,       #
                                      -format => 'fasta');               #��ʱ�ļ��Ĵ򿪣�����д��fasta��ʽ���ļ�, For Mega
  	
  	my $testIdx2=1;  #���ڱ�ʾ������ȶ��ļ��е�����       
  	
  	#���濪ʼһ�� aln�е� ����������
  	#my $sortHash;  #��msf�е��������� ��Seq�����������Ž��µ�hash��
  	
  	foreach my $seq ($aln->each_seq) {
      my $idHere=$seq->id(); my $seqHere=$seq->seq();                                   ##For sort Msf #�ص�     ����õ�ԭ����id�ͱȶԵ�����
      my $sortFasSeqObj=Bio::Seq->new( -display_id => $idHere,                          ##For sort Msf #�ص�     
                                              -seq => $seqHere                );        ##For sort Msf #�ص�     
      $tempSortFast->write_seq($sortFasSeqObj);                                         ##For sort Msf # д����ʱ fasta�ļ�
       
      my $msfId=  $idx_into_protId_Hash->{ $idHere }->{'ShortSP'}."-"
                 .$idx_into_protId_Hash->{ $idHere }->{'SmallFml'}."-"
                 .$idx_into_protId_Hash->{ $idHere }->{'smFmlSpUCinIdx'}."-"
                 .$idx_into_protId_Hash->{ $idHere }->{'CU'};
      $idx_to_newId_hash->{$idHere}=$msfId;

      my $seqobj = Bio::Seq->new( -display_id => $msfId,                            ##For Msf #�ص�     �������fasta��ʽ
                                         -seq => $seqHere                );         ##For Msf #�ص�     ��������Խ�һ�� �������еĸ���
      $tempFast->write_seq($seqobj);                                                 #For Msf # д����ʱ fasta�ļ�  
      
      
      my $idxMsfId=$idHere.$msfId;           
      my $seqobj_Idx = Bio::Seq->new( -display_id => $idxMsfId,                             ##For _Idx Msf #�ص�     �������fasta��ʽ
                                             -seq => $seqHere                );             ##For _Idx Msf #�ص�     ��������Խ�һ�� �������еĸ���
      $tempFast_Idx->write_seq($seqobj_Idx);                                                #For _IdxMsf # д����ʱ fasta�ļ�  
      
                                       
      my $U2CSeq=$seqHere;  $U2CSeq=~s/U/C/g;   $U2CSeq=~s/\./-/g;                 ##For Mega#�ص�     �����е�Uת��ΪC
      my $FmlMark;
      if ($idx_into_protId_Hash->{ $idHere }->{'SmallFml'} eq $idx_into_protId_Hash->{ $idHere }->{'BigFml'}){}
      else{ $FmlMark=$idx_into_protId_Hash->{ $idHere }->{'SmallFml'};$FmlMark=~s/^.*(\w)$/$1/; }
      my $megaId= $idx_into_protId_Hash->{ $idHere }->{'ShortSP'}
                 .$idx_into_protId_Hash->{ $idHere }->{'CU'}
                 .$idx_into_protId_Hash->{ $idHere }->{'smFmlSpUCinIdx'}
                 .$FmlMark;                                                        ##For Mega#�ص�     ��id����ת��                                                
      my $Megaseqobj = Bio::Seq->new( -display_id => $megaId,                      ##For Mega#�ص�     �������fasta��ʽ
                                             -seq => $U2CSeq                );     ##For Mega#�ص�     ��������Խ�һ�� �������еĸ���
                                   
                                         
      
      
      $tempFastTwo->write_seq($Megaseqobj);   #For Mega# д����ʱ fasta�ļ�
     
      $outPutHash->{$idHere}=$testIdx2;       #��������� id�� ��ŵ� ���  
      #print "$testIdx\t$testIdx2\t\t\$idHere=$idHere\n\$seqHere=$seqHere\n\n";   #\$seq->id()=$seq->id()\n";

      $testIdx2++;
    }
    $testIdx++;
    

    
  }

    
  if ($testIdx>0){
  	
  	
  	
  	
  	
    my $fastain= Bio::AlignIO->new(-file   => $tempFile ,              ##  #For Msf #��fasta����ʱ�ļ���Ϊ ����
                                   -format => 'fasta'             );   ##
    my $out    = Bio::AlignIO->new(-file   => ">$clutalwOutfile" ,     #   #For Msf #�����µ��������
                                   -format => $outFormat          );   #   
    
    ##For Msf #������������ �ȶ��ļ��� д��
    while   ( my $FstAln = $fastain->next_aln() ) {
    	$out->write_aln($FstAln);
    }                             
    system ("rm -f $tempFile");    ##ɾ����ʱ�ļ�
    
    
    my $fastain_Idx= Bio::AlignIO->new(-file   => $tempFileIdx ,                      ##  #For Msf_Idx #��fasta����ʱ�ļ���Ϊ ����
                                       -format => 'fasta'             );              ##
    my $out_Idx    = Bio::AlignIO->new(-file   => ">$In_fmlCluOutGdIdxNameFile" ,     #   #For Msf_Idx #�����µ��������
                                       -format => $outFormat          );              #   
    
    ##For Msf_Idx #������������ �ȶ��ļ��� д��
    while   ( my $FstAln_Idx = $fastain_Idx->next_aln() ) {
    	$out_Idx->write_aln($FstAln_Idx);
    }                             
    system ("rm -f $tempFileIdx");    ##ɾ����ʱ�ļ�
    
    &sortFastaFile ($tempFileSort, $sortedTempFil, $idx_to_newId_hash );
    my $fastainSort= Bio::AlignIO->new(-file   => $sortedTempFil,         ##  #For sort Msf #��fasta����ʱ�ļ���Ϊ ����
                                       -format => 'fasta'             );   ##
    my $outSort    = Bio::AlignIO->new(-file   => ">$clutalwOutSortfile" , #   #For sort Msf #�����µ��������
                                       -format => $outFormat          );   #
    
    
    ##For Msf #������������ �ȶ��ļ��� д��
    while   ( my $SortFstAln = $fastainSort->next_aln() ) {
    	$outSort->write_aln($SortFstAln);
    }                             
    system ("rm -f $tempFileSort");    ##ɾ����ʱ�ļ�
    system ("rm -f $sortedTempFil");    ##ɾ����ʱ�ļ�
    
    my $fastainTwo= Bio::AlignIO->new(-file   => $tempFileTwo ,           ##  #For Mega#��fasta����ʱ�ļ���Ϊ ����
                                      -format => 'fasta'             );   ##
    my $outTwo    = Bio::AlignIO->new(-file   => ">$megaOutFile" ,        #   #For Mega#�����µ��������
                                      -format => $megaFormat          );  #                                  
    
    ##For Mega#������������ �ȶ��ļ��� д��
    while   ( my $FstAlnTwo = $fastainTwo->next_aln() ) {
    	$outTwo->write_aln($FstAlnTwo);
    }                             
    system ("rm -f $tempFileTwo");    ##ɾ����ʱ�ļ�
    
  }
  ####################^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^################ for Msf#############################################
  
  
  return $outPutHash;                              
  
}

sub sortFastaFile{ #��fastafile��������ʹ�ó����Ƚ�������
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

sub PrintPngForClustalw{   #   ClustalwRun::PrintPngForClustalw  ($clustalwOUTfile, $Whole_htmlParseOut, $outPng, $msf_new_showNAMEhash);    #����clustalw�����msf�ļ�$clustalwOUTfile,��interproscanģ��Ľ���ļ�$Whole_htmlParseOut,����ͼƬ�������$outPng��
  
  my ($clustalwOUTfile, $Whole_htmlParseOut, $outPng, $msf_new_showNAMEhash, $uPosClumNub )=@_;
  
  
  my $warnMsgBody="\nIn package  ClustalwRun,\tIn sub PrintPngForClustalw,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $caller_inform=DirFileHandle::print_SubCallerInform; print "\n201812241644-1 \n";
  
  
  my $shownNameLengthLimit=100; #���趨�� $msf_new_showNAMEhash ��ʱ�� ���Կ������shown name�ĳ��ȣ������� $shownNameLengthLimit
  
  #����clustalw����ļ�
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
  my $trackHash; #����һ��hash����������Ҫ��ͼ��Trackװ��ȥ���Ա����� ��ͼ����
  
 
  
  #��ȡclustalw �����ȶ� ����������һ�������� ����һ��aln
  while (my $aln=$in->next_aln){
    print "1: ", $aln->length, "\n";   
    print "2: ", $aln->num_residues, "\n";  
    print "3: ", $aln->is_flush, "\n";  
    print "4: ", $aln->num_sequences, "\n";  
    print "5: ", $aln->percentage_identity, "\n";  
    print "6: ", $aln->consensus_string(50), "\n";
    
    
  
    my $DomainIDHash;              #���Hash������ �����е�Domain ID������Key����һ��Hash�С�  ���ں��治ͬ��ɫ��ע�Ĳ�����
    
    my $DS_DomDadSonRelationHash;     #���Hash������ �� �� Daddy-son relation��Domain ��¼������
    my $DS_AllDadSonHash;             #���Hash������ �� ���� ID��������dad����son�������������������hashһ�� ��������ϵ��
    
    my $Bigest_showNameLength=-9999999999999999;  #������������������ֶ����
    
    #���滭ÿ�� �������� track
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
      
      #Ȼ�󻭸����ȶԽ��
           
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
      
      
      
      #Ȼ�󻭳����� domin
     
      #�ҵ� html��ʾ��domain�ļ����ٽ��н����Ⱥ�������
       
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
              
              
              
               
               
                                                                                                                                                                                                            #$dmListHere->[$dlIdx]->[0]=     $DomainStt;                #$dmListHere->[$dlIdx]->[1]=     $DomainEnd;                 #$dmListHere->[$dlIdx]->[2]=     $DomainInD;  #����id               #$dmListHere->[$dlIdx]->[3]=     $DomainStt;  #��������
              if ($DomainEnd>$segLength){$DomainEnd=$segLength; }  #������ĺ��壬��Ҫ�� ��Coil�� patScan�Ľ�� �п��ܱȵ������б������������������ᵼ�� ���� msf��ģ����������������ǿ�е� �����Ƭ�εĳ��Ƚ��м��١�
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
    
    #�����е���ɫ�ַ�����д�뵽 $DomainIDHash����ȣ�����
    my @colorArray=('yellow','red','blue','mistyrose','darkcyan','orange','purple','greenyellow','hotpink','cyan','brown','fuchsia','pink','grey','rosybrown'); #my $colIdx=0;
    if (  ( ref ( $DomainIDHash ) ) eq 'HASH' ){
    	my $DmIDidx=0;
      foreach my $DmID (    sort { $DomainIDHash->{$b}->{'times'} <=> $DomainIDHash->{$a}->{'times'} } (   keys (  %{ $DomainIDHash }  )   )    ){
        my $colorNB=$DmIDidx%@colorArray;  print "\$DomainIDHash->{$DmID}->{'times'}=$DomainIDHash->{$DmID}->{'times'}\t\$DmIDidx=$DmIDidx\t\$colorNB=$colorNB\n";
        $DomainIDHash->{$DmID}->{'Domaincolor'}=$colorArray[$colorNB];
        $DmIDidx++;
      }
    }
    
    
    #���� Domain֮��� ��ϵ���νṹ
    my $domainTreeHash=&MakeTreeFrom_DadSonHashs_onlyMultiLevesTree($DS_AllDadSonHash, $DS_DomDadSonRelationHash);
    
    my $outTreeWord; $outTreeWord=&PrintTreePoints ($outTreeWord, $domainTreeHash, $DomainIDHash);
    #print "\n\n\$outTreeWord=$outTreeWord\n\n\n";
   
    #���Ȼ� panel_1
    my $panel_1 = Bio::Graphics::Panel->new(
                                               -length    => $aln->length,        #ͼ��ȫ������ܱ�ʾ�ļ�������ܳ���
                                               -width     => 400,                #ͼ��ʵ�ʿ��
                                               -key_style => 'right',
                                               -start     => 1,
                                               -pad_left  => 100,
                                               -pad_right => 100,
                                               
                                         );  
    
                                         
    #��ͼʾ                                      
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
    
    #���Ȼ� panel_2
    my $panel_2 = Bio::Graphics::Panel->new(
                                               -length    => $aln->length,        #ͼ��ȫ������ܱ�ʾ�ļ�������ܳ���
                                               -width     => 200,                  #ͼ��ʵ�ʿ��
                                               -key_style => 'right',
                                               -start     => 1,
                                               -pad_left  => 100,
                                               -pad_right => 100,
                                               
                                         );  
    
    $panel_2=&treePanelWork($panel_2, $domainTreeHash, $DomainIDHash, $aln->length);
    
    #���Ȼ� panel_3
    my $panel_3 = Bio::Graphics::Panel->new(
                                               -length    => $aln->length,        #ͼ��ȫ������ܱ�ʾ�ļ�������ܳ���
                                               -width     => 1000,                #ͼ��ʵ�ʿ��
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
    
    #Ȼ�� ����
    my $full_length = Bio::SeqFeature::Generic->new(
                                                         -start => 1, 
                                                         -end   => $aln->length
                                                    );
                                                    
                                                    
                                                    
    # ����ߵ� track���뵽panel�С�
    $panel_3->add_track ( $full_length,
                                           -glyph   => 'arrow',
                                           -tick    => 2,
                                           -fgcolor => 'black',
                                           -double  => 1
                      );
   
   
    
    
    
    my $colorChangeNB=0;
    
    #����ʵ�ʵĻ滭�����������ͼ��  $trackHash->{'_alnSeqArray'}->[$alnSeqIdx]->{'_DomainInformHas'}->{$DomainLineNb}->{'domainRegionArray'}->[$dlIdx]->{'_SegmentList'}   $trackHash->{'_alnSeqArray'}->[$alnSeqIdx]->{'_SegmentList'}
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

sub MakeTreeFrom_DadSonHashs {          #���ø��ӹ�ϵ�������νṹ           #($DS_AllDadSonHash, $DS_DomDadSonRelationHash)
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


sub MakeTreeFrom_DadSonHashs_onlyMultiLevesTree {    #���ø��ӹ�ϵ�������νṹ                 #($DS_AllDadSonHash, $DS_DomDadSonRelationHash)
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

sub DiGuiDadSon{  #���õݹ鷨 ������
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

sub OLD_DiGuiDadSon{  #���õݹ鷨 ������
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

sub PrintTreePoints{  #��ӡ�����νṹ
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

sub treePanelWork{  #��������pngͼ
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



sub MergeImage{  #��png�ļ�$image_1��$image_2,�ϲ���$image_out, Ӧ�������ºϲ�
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
  





sub AddTracks2Panel{  #����������feature��Track�ӵ�Panel��
  my ($tpPanel,  $keyHere, $tkBgCol)=@_;  #warn "\$tpPanel=$tpPanel, \$tkBgCol=$tkBgCol, \$keyHere=$keyHere\n";  print   "\$tpPanel=$tpPanel,  \$tkBgCol=$tkBgCol, \$keyHere=$keyHere\n";
    
  #�Ƚ��������ֶε�feature
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
                                           #-part_labels =>  1,             #��ʾ�����ӱ��
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

sub AddFeature2PanelAndTrack{ #����������feature��Track�ӵ�Panel��
  my ($AAposList, $tpPanel, $tpTrack,  $colorHere)=@_;  #warn "\$AAposList=$AAposList, \$tpPanel=$tpPanel, \$colorHere=$colorHere, \$tpTrack=$tpTrack\n"; 
    
  #�Ƚ��������ֶε�feature
  #$e1 = Bio::SeqFeature::Lite->new(-start=>1,-stop=>100,-type=>'exon');
  my $allSubFeatures;
  for (my $i=0; $i< @{ $AAposList }; $i++){  
    $allSubFeatures->[$i]=Bio::SeqFeature::Lite->new( -start=>$AAposList->[$i]->{'_tailPos'},-stop=>$AAposList->[$i]->{'_headPos'}, -name => $colorHere); #$AAposList->[$i]->{'_SegmentType'} );
  }
  
  #Ȼ���� ��Feature
  my $AlignFeature = Bio::SeqFeature::Lite->new (-segments=>$allSubFeatures,  -strand       => 1);
  #Ȼ�󻭸����ȶԽ��
  
  

  
  $tpTrack->add_feature($AlignFeature,  -bgcolor     =>$colorHere); 
  #$tpTrack->add_feature($AlignFeature_2);
  
  
  return [$tpPanel,$tpTrack];                 
  
  
    


}


sub GetSeqPosList{  #����clustalw�Ľ�� ��column_from_residue_number �ȷ�������þ����ȶԶ���������Ų��ĸ���Ƭ��λ���γɵ�����
  my ($inSeqOgj, $inSeqStart, $inSeqEnd)=@_;
  my $outList;  my $ListIdx=0;
  for (my $i=$inSeqStart; $i<=$inSeqEnd; $i++){
    # �ж��κ�һ���ַ��ǲ��� ��ͷ�ͽ�β����Ҫ�����ĸ�������жϣ��翪ͷ�жϣ�1���ǲ��� �������п�ͷ�� 2 ���ǲ��� ����Ƭ�ο�ͷ�� 3 ���ǲ����ض��ַ� ����U��   4 ���ǲ��� �ض��ַ���� Ƭ�ο�ͷ
    #�����жϿ�ͷ
    #1���ж�����ַ��ǲ�����һ��Ƭ�����ǿ�ͷ������������ַ���ǰһ���ַ���λ������������Ƿ����1������������ǵ�һ��Ƭ�Σ��򲻴���ǰһ���ַ������ԣ����ж��Ƿ��ǵ�һ���ַ���
    if ($i==$inSeqStart){                                                                                                print "\n\n1.0 1.0 If  \$i==$i==\$inSeqStart\n";
    	$outList->[$ListIdx]->[0]=$inSeqOgj->column_from_residue_number($i);                                               print "\t\t\$outList->[\$ListIdx]->[0]=\$outList->[$ListIdx]->[0]=\$inSeqOgj->column_from_residue_number(\$i)=\$inSeqOgj->column_from_residue_number($i)= ", $inSeqOgj->column_from_residue_number($i), "\n";
    } #2���ϲ��жϣ��粻�ǵ�һ��Ƭ�Σ�����elsif ���жϣ����ַ��ǲ��� Ƭ��ͷ��λ��
    elsif ( ($inSeqOgj->column_from_residue_number($i)  -$inSeqOgj->column_from_residue_number($i-1)) > 1 ){             print "1.1 1.1 If    (\$inSeqOgj->column_from_residue_number(\$i)  -\$inSeqOgj->column_from_residue_number(\$i-1)) = (\$inSeqOgj->column_from_residue_number($i)  -\$inSeqOgj->column_from_residue_number($i-1))= (",$inSeqOgj->column_from_residue_number($i),"-",$inSeqOgj->column_from_residue_number($i-1),") =",($inSeqOgj->column_from_residue_number($i) - $inSeqOgj->column_from_residue_number($i-1))," > 1\n" ;
     	$outList->[$ListIdx]->[0]=$inSeqOgj->column_from_residue_number($i);                                               print "\t\t\$outList->[\$ListIdx]->[0]=\$outList->[$ListIdx]->[0]=\$inSeqOgj->column_from_residue_number(\$i)=\$inSeqOgj->column_from_residue_number($i)= ",$inSeqOgj->column_from_residue_number($i),"\n";
    } #3 �����λ��ΪU�������ø�λ��Ϊ Ƭ�ο�ͷ
    elsif ( $inSeqOgj->subseq($inSeqOgj->column_from_residue_number($i),$inSeqOgj->column_from_residue_number($i)) =~m /^U$/i ){  print "1.2 1.2  ��λ��ΪU \n";
    	$outList->[$ListIdx]->[0]=$inSeqOgj->column_from_residue_number($i);                                               print "\t\t\$outList->[\$ListIdx]->[0]=\$outList->[$ListIdx]->[0]=\$inSeqOgj->column_from_residue_number(\$i)=\$inSeqOgj->column_from_residue_number($i)= ",$inSeqOgj->column_from_residue_number($i),"\n";  
    } #4 �����λ�õ���һ���ַ�ΪU�������ø�λ��Ϊ Ƭ�ο�ͷ
    elsif ( $inSeqOgj->subseq($inSeqOgj->column_from_residue_number($i)-1,$inSeqOgj->column_from_residue_number($i)-1) =~m /^U$/i ){  print "1.3 1.3  ��λ�õ��ϸ�λ��ΪU \n";
    	$outList->[$ListIdx]->[0]=$inSeqOgj->column_from_residue_number($i);                                                print "\t\t\$outList->[\$ListIdx]->[0]=\$outList->[$ListIdx]->[0]=\$inSeqOgj->column_from_residue_number(\$i)=\$inSeqOgj->column_from_residue_number($i)= ",$inSeqOgj->column_from_residue_number($i),"\n";  
    }
      
    #Ȼ�����жϸ��ַ��ǲ��ǳ�����һ��Ƭ�εĽ�β����Ϊ�е�Ƭ���п�����һ���ַ�����ɣ�����������жϺ�������ж�ͷ�ַ��Ĳ��� Ӧ����ƽ�еģ����Բ���elseif����if
    #1���ж��Ƿ������һ��Ƭ�ε�β��
    if ($i==$inSeqEnd){                                                                                                  print "\n2.0 2.0 Elsif \$i==$i==\$inSeqEnd\n";
    	$outList->[$ListIdx]->[1]=$inSeqOgj->column_from_residue_number($i);                                               print "\t\t\$outList->[\$ListIdx]->[0]=\$outList->[$ListIdx]->[0]=\$inSeqOgj->column_from_residue_number(\$i)=\$inSeqOgj->column_from_residue_number($i)= ", $inSeqOgj->column_from_residue_number($i), "\n";
    	$outList->[$ListIdx]->[2]=$inSeqOgj->subseq($outList->[$ListIdx]->[0],$outList->[$ListIdx]->[1]);                  warn "222222  \$outList->[$ListIdx]->[2]=$outList->[$ListIdx]->[2]\n";
    	$ListIdx++;
    }  #2�������ж��ǲ��� �в�Ƭ�ε�β��
    elsif ( ($inSeqOgj->column_from_residue_number($i+1)-$inSeqOgj->column_from_residue_number($i)  ) > 1 ){             print "2.1 2.1 ElsIf (\$inSeqOgj->column_from_residue_number(\$i+1)-\$inSeqOgj->column_from_residue_number(\$i)  ) = (\$inSeqOgj->column_from_residue_number($i+1)-\$inSeqOgj->column_from_residue_number($i)  )= (",$inSeqOgj->column_from_residue_number($i+1),"-",$inSeqOgj->column_from_residue_number($i), ")=" ,($inSeqOgj->column_from_residue_number($i+1)- $inSeqOgj->column_from_residue_number($i)  )," > 1\n" ;
      $outList->[$ListIdx]->[1]=$inSeqOgj->column_from_residue_number($i);                                               print "\t\t\$outList->[\$ListIdx]->[1]=\$outList->[$ListIdx]->[1]=\$inSeqOgj->column_from_residue_number(\$i)=\$inSeqOgj->column_from_residue_number($i)= ",$inSeqOgj->column_from_residue_number($i),"\n";
      $outList->[$ListIdx]->[2]=$inSeqOgj->subseq($outList->[$ListIdx]->[0],$outList->[$ListIdx]->[1]);                  warn "2.1 2.1 \$outList->[$ListIdx]->[2]=$outList->[$ListIdx]->[2]\n";
      $ListIdx++;
    }  #3 �����λ��ΪU�������ø�λ��Ϊ Ƭ�ν�β
    elsif ( $inSeqOgj->subseq($inSeqOgj->column_from_residue_number($i),$inSeqOgj->column_from_residue_number($i)) =~m /^U$/i ){  print "2.2 2.2  ��λ��ΪU \n";
    	$outList->[$ListIdx]->[1]=$inSeqOgj->column_from_residue_number($i);                                               print "\t\t\$outList->[\$ListIdx]->[1]=\$outList->[$ListIdx]->[1]=\$inSeqOgj->column_from_residue_number(\$i)=\$inSeqOgj->column_from_residue_number($i)= ",$inSeqOgj->column_from_residue_number($i),"\n";
      $outList->[$ListIdx]->[2]=$inSeqOgj->subseq($outList->[$ListIdx]->[0],$outList->[$ListIdx]->[1]);                  warn "2.2 2.2 \$outList->[$ListIdx]->[2]=$outList->[$ListIdx]->[2]\n";
      $ListIdx++;
    } #4 �����λ�õ���һ���ַ�ΪU�������ø�λ��Ϊ Ƭ�ν�β
    elsif ( $inSeqOgj->subseq($inSeqOgj->column_from_residue_number($i)+1,$inSeqOgj->column_from_residue_number($i)+1) =~m /^U$/i ){  print "2.3 2.3  ��λ�õ��ϸ�λ��ΪU \n";
    	$outList->[$ListIdx]->[1]=$inSeqOgj->column_from_residue_number($i);                                               print "\t\t\$outList->[\$ListIdx]->[1]=\$outList->[$ListIdx]->[1]=\$inSeqOgj->column_from_residue_number(\$i)=\$inSeqOgj->column_from_residue_number($i)= ",$inSeqOgj->column_from_residue_number($i),"\n";
      $outList->[$ListIdx]->[2]=$inSeqOgj->subseq($outList->[$ListIdx]->[0],$outList->[$ListIdx]->[1]);                  warn "2.3 2.3 \$outList->[$ListIdx]->[2]=$outList->[$ListIdx]->[2]\n";
      $ListIdx++;
    }
      
    
  }
  return $outList;
}


sub GetSeqPosList_2 {  #��GetSeqPosList���˸�Type����Ϣ���� #����clustalw�Ľ�� ��column_from_residue_number �ȷ�������þ����ȶԶ���������Ų��ĸ���Ƭ��λ���γɵ�����
  my ($inSeqOgj, $inSeqStart, $inSeqEnd, $inSegType)=@_;  #$inSegType ��������� �����seq ���������ͣ�ͨ����ָ �ض������� �Ƿ���U���Ƿ����������� ���� ������ĩ�˵�
  
  my %keysUsedHere=(
    '_headPos'        =>   "��� Segment ���������е� ���",
    '_tailPos'        =>   "��� Segment ���������е� �յ�",
    '_SegmentType'    =>   "ĳ��������AAλ�õ� ���ͣ�����ĳ�ΰ���������Segmentλ�����������ε� ���ͣ��������ͣ�����������4������"
  );
  
  my $outList;  my $ListIdx=0;
  for (my $i=$inSeqStart; $i<=$inSeqEnd; $i++){   
    # �ж��κ�һ���ַ��ǲ��� ��ͷ�ͽ�β����Ҫ�����ĸ�������жϣ��翪ͷ�жϣ�1���ǲ��� �������п�ͷ�� 2 ���ǲ��� ����Ƭ�ο�ͷ�� 3 ���ǲ����ض��ַ� ����U��   4 ���ǲ��� �ض��ַ���� Ƭ�ο�ͷ
    #�����жϿ�ͷ
    #1���ж�����ַ��ǲ�����һ��Ƭ�����ǿ�ͷ������������ַ���ǰһ���ַ���λ������������Ƿ����1������������ǵ�һ��Ƭ�Σ��򲻴���ǰһ���ַ������ԣ����ж��Ƿ��ǵ�һ���ַ���
    if ($i==$inSeqStart){                                                                                                         print "\n\n1.0 1.0 If  \$i==$i==\$inSeqStart\n";
    	$outList->[$ListIdx]->{'_headPos'}=$inSeqOgj->column_from_residue_number($i);                                               print "\t\t\$outList->[\$ListIdx]->{'_headPos'}=\$outList->[$ListIdx]->{'_headPos'}=\$inSeqOgj->column_from_residue_number(\$i)=\$inSeqOgj->column_from_residue_number($i)= ", $inSeqOgj->column_from_residue_number($i), "\n";
    } #2���ϲ��жϣ��粻�ǵ�һ��Ƭ�Σ�����elsif ���жϣ����ַ��ǲ��� Ƭ��ͷ��λ��
    elsif ( ($inSeqOgj->column_from_residue_number($i)  -$inSeqOgj->column_from_residue_number($i-1)) > 1 ){                      print "1.1 1.1 If    (\$inSeqOgj->column_from_residue_number(\$i)  -\$inSeqOgj->column_from_residue_number(\$i-1)) = (\$inSeqOgj->column_from_residue_number($i)  -\$inSeqOgj->column_from_residue_number($i-1))= (",$inSeqOgj->column_from_residue_number($i),"-",$inSeqOgj->column_from_residue_number($i-1),") =",($inSeqOgj->column_from_residue_number($i) - $inSeqOgj->column_from_residue_number($i-1))," > 1\n" ;
     	$outList->[$ListIdx]->{'_headPos'}=$inSeqOgj->column_from_residue_number($i);                                               print "\t\t\$outList->[\$ListIdx]->{'_headPos'}=\$outList->[$ListIdx]->{'_headPos'}=\$inSeqOgj->column_from_residue_number(\$i)=\$inSeqOgj->column_from_residue_number($i)= ",$inSeqOgj->column_from_residue_number($i),"\n";
    }
      
    #Ȼ�����жϸ��ַ��ǲ��ǳ�����һ��Ƭ�εĽ�β����Ϊ�е�Ƭ���п�����һ���ַ�����ɣ�����������жϺ�������ж�ͷ�ַ��Ĳ��� Ӧ����ƽ�еģ����Բ���elseif����if
    #1���ж��Ƿ������һ��Ƭ�ε�β��
    if ($i==$inSeqEnd){                                                                                                           print "\n2.0 2.0 Elsif \$i==$i==\$inSeqEnd\n";
    	$outList->[$ListIdx]->{'_tailPos'}=$inSeqOgj->column_from_residue_number($i);                                               print "\t\t\$outList->[\$ListIdx]->{'_headPos'}=\$outList->[$ListIdx]->{'_headPos'}=\$inSeqOgj->column_from_residue_number(\$i)=\$inSeqOgj->column_from_residue_number($i)= ", $inSeqOgj->column_from_residue_number($i), "\n";
    	$outList->[$ListIdx]->{'_SegmentType'}=$inSegType;                                                                          #warn "222222  \$outList->[$ListIdx]->{'_SegmentType'}=$outList->[$ListIdx]->{'_SegmentType'}\n";
    	$ListIdx++;
    }  #2�������ж��ǲ��� �в�Ƭ�ε�β��
    elsif ( ($inSeqOgj->column_from_residue_number($i+1)-$inSeqOgj->column_from_residue_number($i)  ) > 1 ){                      print "2.1 2.1 ElsIf (\$inSeqOgj->column_from_residue_number(\$i+1)-\$inSeqOgj->column_from_residue_number(\$i)  ) = (\$inSeqOgj->column_from_residue_number($i+1)-\$inSeqOgj->column_from_residue_number($i)  )= (",$inSeqOgj->column_from_residue_number($i+1),"-",$inSeqOgj->column_from_residue_number($i), ")=" ,($inSeqOgj->column_from_residue_number($i+1)- $inSeqOgj->column_from_residue_number($i)  )," > 1\n" ;
      $outList->[$ListIdx]->{'_tailPos'}=$inSeqOgj->column_from_residue_number($i);                                               print "\t\t\$outList->[\$ListIdx]->{'_tailPos'}=\$outList->[$ListIdx]->{'_tailPos'}=\$inSeqOgj->column_from_residue_number(\$i)=\$inSeqOgj->column_from_residue_number($i)= ",$inSeqOgj->column_from_residue_number($i),"\n";
      $outList->[$ListIdx]->{'_SegmentType'}=$inSegType;                                                                          #warn "2.1 2.1 \$outList->[$ListIdx]->{'_SegmentType'}=$outList->[$ListIdx]->{'_SegmentType'}\n";
      $ListIdx++;
    }  
      
    
  }
  return $outList;
}



sub GetSeqPosList_3 {  #����Ͱ�������domain���������м䱻������������������
  my ($inSeqOgj, $inSeqStart, $inSeqEnd, $inSegType)=@_;  #$inSegType ��������� �����seq ���������ͣ�ͨ����ָ �ض������� �Ƿ���U���Ƿ����������� ���� ������ĩ�˵�
  
  my %keysUsedHere=(
    '_headPos'        =>   "��� Segment ���������е� ���",
    '_tailPos'        =>   "��� Segment ���������е� �յ�",
    '_SegmentType'    =>   "ĳ��������AAλ�õ� ���ͣ�����ĳ�ΰ���������Segmentλ�����������ε� ���ͣ��������ͣ�����������4������"
  );
  
  my $outList;  
  $outList->[0]->{'_headPos'}    =$inSeqOgj->column_from_residue_number($inSeqStart);
  $outList->[0]->{'_tailPos'}    =$inSeqOgj->column_from_residue_number($inSeqEnd  ); 
  $outList->[0]->{'_SegmentType'}=$inSegType;  
  
  return $outList;
}





sub DrawClustalw{  #����clustlaw����� pngͼ #������������������û�б�ʹ�õ�
	my ($clustalwOUTfile, $formatHere, $clustalPNGfile)=@_;

  
  #����clustalw����ļ�
  my $in  = Bio::AlignIO->new(-file   => $clustalwOUTfile ,
                              -format => $formatHere            
                              #-report_type => 'blastn'    
                              );

  open (OUT,">$clustalPNGfile") or die "cannot create $clustalPNGfile : $!";

  #my $testNumber=130;

  #��ȡclustalw �����ȶ� ����������һ�������� ����һ��aln
  my $howManySeq=0;
  my $aln=$in->next_aln();
  #while (my $aln=$in->next_aln)
  if (1){
    warn "1: ", $aln->length, "\n";   warn "2: ", $aln->num_residues, "\n";  warn "3: ", $aln->is_flush, "\n";  warn "4: ", $aln->num_sequences, "\n";  warn "5: ", $aln->percentage_identity, "\n";  warn "6: ", $aln->consensus_string(50), "\n";
    
    #���Ȼ� panel
    my $panel = Bio::Graphics::Panel->new(
                                               -length    => $aln->length,   #ͼ��ȫ������ܱ�ʾ�ļ�������ܳ���
                                               -width     => 1000,                #ͼ��ʵ�ʿ��
                                               -key_style => 'between',
                                               -start     => 1,
                                               -pad_left  => 100,
                                               -pad_right => 100,
                                               
                                         );  
    
    #Ȼ�� ����
    my $full_length = Bio::SeqFeature::Generic->new(
                                                         -start => 1, 
                                                         -end   => $aln->length
                                                        );
    # ����ߵ� track���뵽panel�С�
    $panel->add_track ( $full_length,
                                           -glyph   => 'arrow',
                                           -tick    => 2,
                                           -fgcolor => 'black',
                                           -double  => 1
                      );
                      
    #���滭ÿ�� �������� track
    
    foreach my $seq ($aln->each_seq) {          warn "id: ",$seq->id(),"\n";     warn "Seq: ",$seq->seq(),"\n";    warn "Start: ",$seq->start(),"\n";    warn "End: ",$seq->end(),"\n";    warn "Number of gaps 1 .: ", $seq->num_gaps('.'),"\n";    warn "Number of gaps 2 -: ", $seq->num_gaps('-'),"\n";    warn "Number of gaps 2  : ", $seq->num_gaps(),"\n";            #warn "Col Number of rsd 20  : ", $seq->column_from_residue_number(20),"\n";    #warn "location_from_column 82  : ", $seq->location_from_column(82),"\n";       #warn "\n";    #my @tempMap=$seq->mapping();    #my $idxNb=0;foreach my $temp(@tempMap){warn "$idxNb: \$temp=$temp\t\t";$idxNb++}    #warn "\nAfter Print \@tempMap\n\n";        #my %hashTp=$seq->frameshifts();    #foreach my $keyH (keys (%hashTp)){warn "\$hashTp{\$keyH}=\$hashTp{$keyH}=$hashTp{$keyH}\t\t";} warn "\nAfter Print \%hashTp\n\n\n\n";        #my $loc82loc=$seq->location_from_column(82);    #if (defined ($loc82loc)){ warn "Loc82 LocationType: ",$loc82loc->location_type," Start: ",$loc82loc->start,"\tEnd: ",$loc82loc->end,"\nAfter Print \$loc82loc\n\n\n\n\n";}    
      
      my $outAAposList=&GetSeqPosList($seq,$seq->start(),$seq->end());        warn "Print AA seq List: ", &PrintList( $outAAposList ), "\n\n\n\n\n";
      my $segName=$seq->id(); my $segLength=$seq->end();        
      #Ȼ�󻭸����ȶԽ��
      $panel=&NEWAddFeatureTrack2Panel($outAAposList,$panel,'','green',"$segName  $segLength");
      $howManySeq++;
    }
    print OUT  $panel->png; 
  }
  close (OUT);  
     
  return $howManySeq;
	
	
}
sub PrintList{  #������������������û�б�ʹ�õ� #������ڴ�ӡ������ ������λ��� �����磺 116..156,158..181,184..193,195..202,206..263,265..313,315..330
  my ($inList)=@_;
  my @outList;
  for (my $j=0; $j<@{ $inList }; $j++){
    push @outList, "$inList->[$j]->[0]..$inList->[$j]->[1]";
  }
  my $outListLine=join (",", @outList);
  return $outListLine;
}



sub NEWAddFeatureTrack2Panel{  #��Panel��Track
  my ($AAposList, $tpPanel, $idKeyHere, $colorHere, $keyHere)=@_;
    
  #�Ƚ��������ֶε�feature
  #$e1 = Bio::SeqFeature::Lite->new(-start=>1,-stop=>100,-type=>'exon');
  my $allSubFeatures;
  for (my $i=0; $i< @{ $AAposList }; $i++){  
    $allSubFeatures->[$i]=Bio::SeqFeature::Lite->new( -start=>$AAposList->[$i]->[0],-stop=>$AAposList->[$i]->[1],-seq  =>  $AAposList->[$i]->[2], -id => "onYeah $i" );
  }
  
  #Ȼ���� ��Feature
  my $AlignFeature = Bio::SeqFeature::Lite->new (-segments=>$allSubFeatures, -id =>$idKeyHere, -strand       => 1);
  #Ȼ�󻭸����ȶԽ��
  
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
                                           #-part_labels =>  1,             #��ʾ�����ӱ��
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