
#!/usr/bin/perl -w

BEGIN{  push (@INC, '/home/fredjiang/fredProgs/PMs')  };

use strict;
use warnings;

use Storable ;
use Storable qw(nstore store_fd nstore_fd freeze thaw dclone);

use PrintSubArrayHash;
use FastaFileHandle;
use SeqSegmentsTools;
use BlastHandle;
use DirFileHandle;
use DieWork;
use InFileHandle;
use TaxonomyWork_NEW;

package  GenomeDATAfeach;

#NCBI 有关 基因组数据的 下载地址等信息 都在 ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/


#NCBI 有关 基因组数据的 下载地址等信息 都在 ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/

my $assembly_gbk_url="ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt";
my $assembly_ref_url="ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt";

#旧服务器上的 上层文件夹地址
my $upLevelDIR                                          ="/home/fredjiang/OneTSSD/fredjiang20190321";

#新服务器上的 下载数据的地址
#各种 解析文件夹 都在该 文件夹下
my $JL_GenomeWork_Dir                                   ="/mnt/md1/datahub";
#my $JL_GenomeWork_Dir                                   ="/home/fredjiang/md1/fredjiang20190916/GenomeData20190319";

#   reprot解析 文件夹

my $JL_all_assemblyDownload_DIR                         =$upLevelDIR."/GENOME_REPORTS";
my $JL_difVersion_DIR                                   =$JL_all_assemblyDownload_DIR."/G_REPORTS_20191227";  # 每次更新 更改这个文件夹的名字

my $JL_GenomeWork_assembly_report_Dir                   =$JL_difVersion_DIR."/ASSEMBLY_REPORTS";



#各种 从 assembly文件 解析出来的信息，这些hash文件的内容，由输入的assembly文件（下载于ncbi）所决定的。
#下载操作时也是依据 这些文件来进行下载和更新的
my $JL_GenomeWork_assembly_GenBankSummary               =$JL_GenomeWork_assembly_report_Dir."/1m0m_assembly_summary_genbank.txt";                             #原始文件，ncbi下载  未注释基因组信息
my $JL_GenomeWork_assembly_GenBankSummary_HASH          =$JL_GenomeWork_assembly_report_Dir."/1m1m_assembly_summary_genbank.txt.hash";
my $JL_GenomeWork_assembly_GenBankSummary_TaxidKey_HASH =$JL_GenomeWork_assembly_report_Dir."/1m2m_assembly_summary_genbank.txt.TaxidKey.hash";
my $JL_GenomeWork_assembly_GenBank_TaxidKey_Linage_HASH =$JL_GenomeWork_assembly_report_Dir."/1m3m_assembly_summary_genbank.txt.TaxidKey.lineage.hash";

my $JL_GenomeWork_assembly_refSeqSummary                =$JL_GenomeWork_assembly_report_Dir."/2m0m_assembly_summary_refseq.txt";                              #原始文件，ncbi下载  已注释基因组信息
my $JL_GenomeWork_assembly_refSeqSummary_HASH           =$JL_GenomeWork_assembly_report_Dir."/2m1m_assembly_summary_refseq.txt.hash";
my $JL_GenomeWork_assembly_refSeqSummary_TaxidKey_HASH  =$JL_GenomeWork_assembly_report_Dir."/2m2m_assembly_summary_refseq.txt.TaxidKey.hash";
my $JL_GenomeWork_assembly_refSeq_TaxidKey_Linage_HASH  =$JL_GenomeWork_assembly_report_Dir."/2m3m_assembly_summary_refseq.txt.TaxidKey.lineage.hash";


#下面是各个 包含下载genome的文件的地址的 hash，这些文件是随着下载过程 生成的。
my $JL_GenomeWork_assembly_download_Eukayotic_HASH      =$JL_GenomeWork_assembly_report_Dir."/3m1m_assembly_summary_eukaryotes_data.txt";
my $JL_GenomeWork_assembly_download_Bacterial_HASH      =$JL_GenomeWork_assembly_report_Dir."/3m2m_assembly_summary_Bacteria___data.txt";


######下面的reports 输出，由于没有ftp信息， 使用频率并不高
my $JL_GenomeWork_reports_Dir                           =$JL_difVersion_DIR."/GENOME_REPORTS";
my $JL_GenomeWork_reports_euakryotes                    =$JL_GenomeWork_reports_Dir."/eukaryotes.txt";
my $JL_GenomeWork_reports_euakryotes_HASH               =$JL_GenomeWork_reports_Dir."/eukaryotes.txt.hash";



#############################################################
#程序运行示例：

  #第0步，更改 文件夹$JL_difVersion_DIR =$JL_all_assemblyDownload_DIR."/GENOME_REPORTS_20191227"; 中双引号内的名字
  #然后利用下面的函数下载两个assembly report文件
#GenomeDATAfeach::wget_assembly_information( );
  #再初略比较是否需要更新

  #第1步，将assembly的summary的文本文件 转化为 以 assembly id，（如 GCF_000001215.4 GCA_000001215.4）为key的 hash
#GenomeDATAfeach::BuildSummaryHASH();
  #第2步，将assembly的summary的文本文件 转化为 以  taxnomy ID （24, 7953） 为key 的hash
#GenomeDATAfeach::GetTxidHash_for_Genebank_and_refSeq_File();
  #第3步，将所有assembly summary的文件中涉及到的 taxnomy id，都到 ncbi上获取xml文件，更新服务器上的 物种信息数据库
#GenomeDATAfeach::CollectAllTaxInform_for_Genebank_and_refSeq_File();
  #第4步，获得所有物种的 genbank和rsf类型 的linage信息
#GenomeDATAfeach::Build_all_taxid_to_linage_hash_for_Genebank_and_refSeq_File();
 
 #第5步，下载并生成带有物种保存地址信息的hash。
   #GenomeDATAfeach::RSYNC_all_genomeData_for_a_eukaryotes();  #下载真核生物
   #GenomeDATAfeach::RSYNC_all_genomeData_for_Bacteria();      #下载原核生物
#############################################################


sub wget_assembly_information{	 #GenomeDATAfeach::wget_assembly_information( );
	
	my ($warnMsgBody, $warnMsgHead, $die_MsgHead, $caller_inform)=@{ DieWork::BuildWarnDieHeadInform( 'GenomeDATAfeach', 'wget_assembly_information' ) };
	
	#DirFileHandle::BuildWarnDieHeadInform ($workingDIRpath);
	#system ("mkdir -p $JL_GenomeWork_assembly_report_Dir ");
	my $mkDir_wkingDirPath="mkdir -p $JL_GenomeWork_assembly_report_Dir";  	DieWork::Print_and_warn( "\n$warnMsgBody$mkDir_wkingDirPath\n\n");
	system ( "$mkDir_wkingDirPath");
	
	&wgetSub($assembly_gbk_url, $JL_GenomeWork_assembly_GenBankSummary);  DieWork::Check_FileDirExist_or_DIE   ( $JL_GenomeWork_assembly_GenBankSummary,     "\$JL_GenomeWork_assembly_GenBankSummary",    $die_MsgHead, $caller_inform  );
	&wgetSub($assembly_ref_url, $JL_GenomeWork_assembly_refSeqSummary);   DieWork::Check_FileDirExist_or_DIE   ( $JL_GenomeWork_assembly_refSeqSummary,      "\$JL_GenomeWork_assembly_refSeqSummary",     $die_MsgHead, $caller_inform  );
	
	
}

sub wgetSub{
	my ($urlHere, $localPath)=@_;
	$urlHere=~s/^ftp:\/\///;
	#rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/102/435/GCA_002102435.1_Ovir.te_1.0	           genomes/all/GCA/002/102/435
	my $rsyncCMD="rsync --copy-links --recursive --times --verbose rsync://$urlHere     $localPath";   warn "$rsyncCMD\n";  print "$rsyncCMD\n";
	system ( "$rsyncCMD");
	#my $wgetCMD="wget -o $localPath $urlHere  ";
	#DieWork::Print_and_warn( "\n\$wgetCMD=$wgetCMD\n\n");
	#system ( "$wgetCMD");
	
}


##################################
#                                #
# RSYNC 下载更新 函数组          # 
#                                #
##################################

##################################
# RSYNC 下载更新 函数组          # 1
################################## 
#从NCBI下载所有的 Bacteria 的基因组数据，如重复使用该函数用于已下载的数据，则运行效果是 仅仅更新数据，并不重新下载
#没有函数输入，有全局变量 'Eukaryota',  和 $JL_GenomeWork_assembly_download_Eukayotic_HASH
sub RSYNC_all_genomeData_for_a_eukaryotes{  #GenomeDATAfeach::RSYNC_all_genomeData_for_a_eukaryotes();
	GenomeDATAfeach::RSYNC_all_genomeData_for_a_phylaName ( 'Eukaryota', $JL_GenomeWork_assembly_download_Eukayotic_HASH );
}

##################################
# RSYNC 下载更新 函数组          # 2
################################## 
#从NCBI下载所有的 Bacteria 的基因组数据，如重复使用该函数用于已下载的数据，则运行效果是 仅仅更新数据，并不重新下载
#没有函数输入，有全局变量 'Eukaryota',  和 $JL_GenomeWork_assembly_download_Eukayotic_HASH
sub RSYNC_all_genomeData_for_Bacteria{  #GenomeDATAfeach::RSYNC_all_genomeData_for_Bacteria();
	GenomeDATAfeach::RSYNC_all_genomeData_for_a_phylaName ( 'Bacteria', $JL_GenomeWork_assembly_download_Bacterial_HASH );
}

##################################
# RSYNC 下载更新 函数组          # 3
################################## 

#  输入1 进化分支名，如 Cervidae， 和输入2 期望的输出文件名，获得所有该进化分支下的基因组，以及生成 数据情况HASH
#还有 隐藏的 全局变量
# 全局变量 1 $JL_GenomeWork_assembly_GenBankSummary_TaxidKey_HASH
# 全局变量 2 $JL_GenomeWork_assembly_GenBank_TaxidKey_Linage_HASH
# 全局变量 3 $JL_GenomeWork_assembly_refSeqSummary_TaxidKey_HASH
# 全局变量 4 $JL_GenomeWork_assembly_refSeq_TaxidKey_Linage_HASH
#输出是 包含新的 本地 存储地址的hash
#调用了 GenomeDATAfeach::FindAllTxid_for_phlyaName_fromHASHfile 用于寻找 phyla_Name下的所有物种
#调用了 ArrayHashChange::Change_Hash_to_Array 将hash变为Array
#调用了 GenomeDATAfeach::RSYNC_remote_GENOME_data_from_Taxid 来进行 下载和更新
sub RSYNC_all_genomeData_for_a_phylaName{  #  my $GBK_and_RSF_outHASH=GenomeDATAfeach::RSYNC_all_genomeData_for_a_phylaName ( $phlya_Name, $outHASH_withLocal_Name );
	my ( $phlya_Name, $outHASH_withLocal_Name )=@_;
	
	my $warnMsgBody="\nIn package  GenomeDATAfeach,\tIn sub RSYNC_all_genomeData_for_a_phylaName,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody"; 
  my $subCallereIfm=DirFileHandle::print_SubCallerInform; #warn $subCallereIfm;
	#warn "20190403-0-0-0 $phlya_Name, $JL_GenomeWork_assembly_GenBankSummary_TaxidKey_HASH, $JL_GenomeWork_assembly_GenBank_TaxidKey_Linage_HASH\n";
	my $GeneBank_outTaxid_to_AssembyHASH=GenomeDATAfeach::FindAllTxid_for_phlyaName_fromHASHfile($phlya_Name, $JL_GenomeWork_assembly_GenBankSummary_TaxidKey_HASH, $JL_GenomeWork_assembly_GenBank_TaxidKey_Linage_HASH);
	my $RefSeq___outTaxid_to_AssembyHASH=GenomeDATAfeach::FindAllTxid_for_phlyaName_fromHASHfile($phlya_Name, $JL_GenomeWork_assembly_refSeqSummary_TaxidKey_HASH,  $JL_GenomeWork_assembly_refSeq_TaxidKey_Linage_HASH);
	my $GeneBank_txid_Array=ArrayHashChange::Change_Hash_to_Array ($GeneBank_outTaxid_to_AssembyHASH);
	my $RefSeq___txid_Array=ArrayHashChange::Change_Hash_to_Array ($RefSeq___outTaxid_to_AssembyHASH);
	
	my $GeneBank_OutHASH_LocalPath_file;   $GeneBank_OutHASH_LocalPath_file=$outHASH_withLocal_Name.".GBK.HASH" if (   (  defined ( $outHASH_withLocal_Name )  ) && ( $outHASH_withLocal_Name =~ m/\S+/)   );
	my $RefSeq___OutHASH_LocalPath_file;   $RefSeq___OutHASH_LocalPath_file=$outHASH_withLocal_Name.".RSF.HASH" if (   (  defined ( $outHASH_withLocal_Name )  ) && ( $outHASH_withLocal_Name =~ m/\S+/)   );
	
	my $GeneBank_OutHASH_with_LocalPath=GenomeDATAfeach::RSYNC_remote_GENOME_data_from_Taxid ($GeneBank_txid_Array, $GeneBank_outTaxid_to_AssembyHASH, $GeneBank_OutHASH_LocalPath_file);
	my $RefSeq___OutHASH_with_LocalPath=GenomeDATAfeach::RSYNC_remote_GENOME_data_from_Taxid ($RefSeq___txid_Array, $RefSeq___outTaxid_to_AssembyHASH, $RefSeq___OutHASH_LocalPath_file);
	
	my $GBK_and_RSF_outHASH;
	$GBK_and_RSF_outHASH->{'0_0_0_GenBank'}=$GeneBank_OutHASH_with_LocalPath;
	$GBK_and_RSF_outHASH->{'0_0_1__RefSeq'}=$RefSeq___OutHASH_with_LocalPath;
	
	if (   (  defined ( $outHASH_withLocal_Name )  ) && ( $outHASH_withLocal_Name =~ m/\S+/)   ){
		if  (   (  defined ( $GBK_and_RSF_outHASH )  ) && (  ref ( $GBK_and_RSF_outHASH ) eq 'HASH'  )   ){
			DirFileHandle::PrintDumper( $outHASH_withLocal_Name, $GBK_and_RSF_outHASH);
		}
	}
	
	return $GBK_and_RSF_outHASH;
	
}

##################################
# RSYNC 下载更新 函数组          # 4
################################## 
#      核心 运算 函数            #  
##################################

#本函数的功能，是利用 taxnomy id为key的ARRAY及HASH中的ftp地址，到NCBI上用rsync命令获取相应物种的基因组数据，下载到本地的特定地址，并将该地址保存于输出的HASH
#输入1    所有需要 下载的taxnomy id为元素的ARRAY（这个ARRAY是由 输入2 的HASH，进行ArrayHashChange::Change_Hash_to_Array转换而来）， 
#输入2    所有需要 下载的taxnomy id为key 的HASH
#输入3    最后输出的 hash的dump文件的文件路径
#输出     包含了 基因组保存地址的 新的hash文件
sub RSYNC_remote_GENOME_data_from_Taxid{  #  GenomeDATAfeach::RSYNC_remote_GENOME_data_from_Taxid ($inTaxidArray, $inGeneBankTaxIDHASH, $OutHASH_FILE_withLocalPath);
	my ($inTaxidArray, $inGeneBankTaxIDHASH, $OutHASH_FILE_withLocalPath,  $only_obtain_theInformHASH)=@_;
	
	my $warnMsgBody="\nIn package  GenomeDATAfeach,\tIn sub RSYNC_remote_GENOME_data_from_Taxid,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody"; 
  my $subCallereIfm=DirFileHandle::print_SubCallerInform; #warn $subCallereIfm;
  
  if  (   (  defined ( $inTaxidArray )  ) && (  ref ( $inTaxidArray ) eq 'ARRAY'  )   ){}
  else{  	DieWork::Just_dieWork( $die_MsgHead."\n \$inTaxidArray=$inTaxidArray\n should be a ARRAY ref  $!\n\n\n".$subCallereIfm ); 	  }
  	
  if  (   (  defined ( $inGeneBankTaxIDHASH )  ) && (  ref ( $inGeneBankTaxIDHASH ) eq 'HASH'  )   ){}
  else {  	DieWork::Just_dieWork( $die_MsgHead."\n \$inGeneBankTaxIDHASH=$inGeneBankTaxIDHASH\n should be a Hash ref  $!\n\n\n".$subCallereIfm ); 	  }
  
  my $howManyTaxid=@ { $inTaxidArray };
  my $taxid_count=0;
  
  my $out_HASH_with_DataPath; 
	foreach my $eachSPecNM_id (  @ { $inTaxidArray }  ){ 
		my $realCount=$taxid_count+1;
		warn "\n\n\n\n\t\t\t$howManyTaxid\tTaxnomy id to search\n\t\tNow the $realCount is working!!!\n\n\n";
		if  (   (  defined ( $inGeneBankTaxIDHASH->{$eachSPecNM_id} )  ) && (  ref ( $inGeneBankTaxIDHASH->{$eachSPecNM_id} ) eq 'ARRAY'  )   ){
			
			my $howManyAssembly_in_this_taxid=@ { $inGeneBankTaxIDHASH->{$eachSPecNM_id} };
      my $Assembly_count=0;
			foreach my $eachAssemblyHASH (  @ { $inGeneBankTaxIDHASH->{$eachSPecNM_id} }  ){ 
			  my $realAssembly_count=$Assembly_count+1;
		    warn "\n\n\n\n\t\t\t$howManyAssembly_in_this_taxid\tAssembly to search for TaxID:$eachSPecNM_id($realCount/$howManyTaxid) \n\t\tNow the $realAssembly_count Assembly is working!!!\n\n\n";
			  if  (   (  defined ( $eachAssemblyHASH )  ) && (  ref ( $eachAssemblyHASH ) eq 'HASH'  ) && (  defined ( $eachAssemblyHASH->{'ftp_path'} )  ) && ( $eachAssemblyHASH->{'ftp_path'} =~ m/\S+/) && (  defined ( $eachAssemblyHASH->{'assembly_accession'} )  ) && ( $eachAssemblyHASH->{'assembly_accession'} =~ m/\S+/)   ){
			  	my $ftp_path=$eachAssemblyHASH->{'ftp_path'};
			  	
			  	                      #ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/789/105/GCF_003789105.1_ASM378910v1
			  	if ( $ftp_path=~m/^na$/ ){ 
			  	  
			  	}			  	
			  	elsif ( $ftp_path=~m/^\s*ftp:\/\/(ftp\.ncbi\.nlm\.nih\.gov\/(genomes\/all\/\w{3}\/\d+\/\d+\/\d+)\/(\S+))\s*$/ ){  
			  		my $rsync_wPath=$1;
			  		my $genome_Path=$2;
			  		my $assemblyDIR=$3;
			  		if (   (  defined ( $rsync_wPath )  ) && ( $rsync_wPath =~ m/\S+/) && (  defined ( $genome_Path )  ) && ( $genome_Path =~ m/\S+/) && (  defined ( $assemblyDIR )  ) && ( $assemblyDIR =~ m/\S+/)   ){
			  			#mkdir -p /home/fredjiang/EightT2/fredjiang20190319/GenomeData20190319/genomes/all/GCA/003/697/985
			  		  #rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/102/435/GCA_002102435.1_Ovir.te_1.0	           genomes/all/GCA/002/102/435
			  	    my $mkdirCMD="mkdir -p $JL_GenomeWork_Dir/$genome_Path"; warn "$mkdirCMD\n";  print "$mkdirCMD\n";
			  	    system ( "$mkdirCMD" );
			  	    my $rsyncCMD="rsync --copy-links --recursive --times --verbose rsync://$rsync_wPath     $JL_GenomeWork_Dir/$genome_Path";   warn "$rsyncCMD\n";  print "$rsyncCMD\n";
			  	    system ( "$rsyncCMD" );
			  	    my $TempHASH=Storable::dclone( $eachAssemblyHASH );
			  	    $TempHASH->{'0_0_0_DataPath'}="$JL_GenomeWork_Dir/$genome_Path/$assemblyDIR";
			  	    
			  	    $out_HASH_with_DataPath->{$eachSPecNM_id}->{ $eachAssemblyHASH->{'assembly_accession'} }=$TempHASH;
			  		}
			  		else{
			  			DieWork::Just_dieWork( $die_MsgHead."\n \$rsync_wPath=$rsync_wPath \$genome_Path=$genome_Path \$assemblyDIR=$assemblyDIR\n should all be a defined and notNull string  $!\n\n\n".$subCallereIfm ); 
			  		}
			  	}
			  	else{
			  		DieWork::Just_dieWork( $die_MsgHead."\n \$ftp_path=$ftp_path\n didnot fit the regular expression  $!\n\n\n".$subCallereIfm ); 
			  	}
			  	
			  }
			  $Assembly_count++;
			}
		}
		
		$taxid_count++;
	}
	
	if (   (  defined ( $OutHASH_FILE_withLocalPath )  ) && ( $OutHASH_FILE_withLocalPath =~ m/\S+/)   ){
		if  (   (  defined ( $out_HASH_with_DataPath )  ) && (  ref ( $out_HASH_with_DataPath ) eq 'HASH'  )   ){
			DirFileHandle::PrintDumper( $OutHASH_FILE_withLocalPath, $out_HASH_with_DataPath);
		}
	}
	return $out_HASH_with_DataPath;
	
}



##################################
#                                #
# 找寻phyla下所有 测序物种的函数 # 
#                                #
##################################

##################################
# 找寻phyla下所有 测序物种的函数 #1
################################## 

#这个函数和下一个函数FindAllTxid_for_phlyaName功能是一致的，只是 这个函数输入的是 hash的dump文件，下一个函数输入的就是hash
#输入一个 phyla name，可以理解为 taxnomy分类学的 任一个 节点的名， 如 脊椎动物veterbrate，或真核生物eukaryotics等等。
#还要输入 一个包含了所有taxnomy id的hash，从这个hash中找所有的taxnomy id，用来进行检测，看看该taxnomy id，是不是改phyla的下游，如“鸭嘴兽”的taxnomy id，是不是veterbrate的下游
#还要输入 一个包含了所有 taxnomy id的linage的信息的hash，从该hash中找到所有的 linage信息
#输出的是 该phyla name下的所有物种的taxnomy id形成的 hash
sub FindAllTxid_for_phlyaName_fromHASHfile{  #my $outTaxid_to_AssembyHASH=GenomeDATAfeach::FindAllTxid_for_phlyaName_fromHASHfile($phlya_Name, $GenBk_or_refSeq_SummaryTxidHASH_file, $Lineage_HASH_file);
	my ($phlya_Name, $GenBk_or_refSeq_SummaryTxidHASH_file, $Lineage_HASH_file)=@_;
	
	warn "\$phlya_Name=$phlya_Name, \$GenBk_or_refSeq_SummaryTxidHASH_file=$GenBk_or_refSeq_SummaryTxidHASH_file, \$Lineage_HASH_file=$Lineage_HASH_file\n";
	
	my $warnMsgBody="\nIn package  GenomeDATAfeach,\tIn sub FindAllTxid_for_phlyaName_fromHASHfile,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $subCallereIfm=DirFileHandle::print_SubCallerInform;
	
	my $GenBk_or_refSeq_SummaryTxidHASH;	
	if (   (  defined ( $GenBk_or_refSeq_SummaryTxidHASH_file )  ) && ( $GenBk_or_refSeq_SummaryTxidHASH_file=~m/\S+/ ) && (  -e ( $GenBk_or_refSeq_SummaryTxidHASH_file )  )   ){
		$GenBk_or_refSeq_SummaryTxidHASH=Storable::retrieve( $GenBk_or_refSeq_SummaryTxidHASH_file );				
	}
	else{
		DieWork::Just_dieWork( $die_MsgHead."\n \$GenBk_or_refSeq_SummaryTxidHASH_file=$GenBk_or_refSeq_SummaryTxidHASH_file\n should be a existed defined file  $!\n\n\n".$subCallereIfm ); 	
	}
  
  my $Lineage_HASH;
  if (   (  defined ( $Lineage_HASH_file )  ) && ( $Lineage_HASH_file=~m/\S+/ ) && (  -e ( $Lineage_HASH_file )  )   ){
		$Lineage_HASH=Storable::retrieve( $Lineage_HASH_file );				
	}
	else{
		DieWork::Just_dieWork( $die_MsgHead."\n \$Lineage_HASH_file=$Lineage_HASH_file\n should be a existed defined file  $!\n\n\n".$subCallereIfm ); 	
	}
	
	my $outTaxid_to_AssembyHASH=GenomeDATAfeach::FindAllTxid_for_phlyaName($phlya_Name, $GenBk_or_refSeq_SummaryTxidHASH, $Lineage_HASH);	
	return $outTaxid_to_AssembyHASH;
	
}


##################################
# 找寻phyla下所有 测序物种的函数 #2
################################## 
#      核心 运算 函数            #  
##################################

#输入一个 phyla name，可以理解为 taxnomy分类学的 任一个 节点的名， 如 脊椎动物veterbrate，或真核生物eukaryotics等等。
#还要输入 一个包含了所有taxnomy id的hash，从这个hash中找所有的taxnomy id，用来进行检测，看看该taxnomy id，是不是改phyla的下游，如“鸭嘴兽”的taxnomy id，是不是veterbrate的下游
#还要输入 一个包含了所有 taxnomy id的linage的信息的hash，从该hash中找到所有的 linage信息
#输出的是 该phyla name下的所有物种的taxnomy id形成的 hash
sub FindAllTxid_for_phlyaName{  #my $outTaxid_to_AssembyHASH=GenomeDATAfeach::FindAllTxid_for_phlyaName($phlya_Name, $GenBk_or_refSeq_SummaryTxidHASH , $Lineage_HASH);
	my ($phlya_Name, $GenBk_or_refSeq_SummaryTxidHASH, $Lineage_HASH)=@_;
	
	my $warnMsgBody="\nIn package  GenomeDATAfeach,\tIn sub FindAllTxid_for_phlyaName,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $subCallereIfm=DirFileHandle::print_SubCallerInform;
	
	my $outTaxid_to_AssembyHASH;
	if (   (  defined ( $GenBk_or_refSeq_SummaryTxidHASH )  ) && (  ref ( $GenBk_or_refSeq_SummaryTxidHASH ) eq 'HASH'  ) && (  defined ( $Lineage_HASH )  ) && (  ref ( $Lineage_HASH ) eq 'HASH'  )   ){
		
		foreach my $Txid (    sort  (   keys (  %{ $GenBk_or_refSeq_SummaryTxidHASH } )   )    ){
		  if (   (  defined ( $Lineage_HASH->{$Txid}->{$phlya_Name} )  ) && ( $Lineage_HASH->{$Txid}->{$phlya_Name} =~ m/\d+/ )   ){
			  $outTaxid_to_AssembyHASH->{$Txid}=$GenBk_or_refSeq_SummaryTxidHASH->{$Txid};
		  }	
		}	
	}
	else{
		DieWork::Just_dieWork( $die_MsgHead."\n \$GenBk_or_refSeq_SummaryTxidHASH=$GenBk_or_refSeq_SummaryTxidHASH and \$Lineage_HASH=$Lineage_HASH\n should all be a Hash ref  $!\n\n\n".$subCallereIfm ); 	
	}
	return $outTaxid_to_AssembyHASH;
	
}


##################################
#                                #
# 融合 GenBank 和 RefSeq 函数组  # 
#          没有写完              #
#                                #
##################################

##################################
# 融合 GenBank 和 RefSeq 函数组  # 1
################################## 

#这个没有写完
sub Merge_GenBank_RefSeq_TaxID_HASH{
	my ($In_GBK_TaxIDKEY_HASH, $In_RSF_TaxIDKEY_HASH)=@_;
	
	#my $In_GBK_AsmbIDKEY_HASH=GenomeDATAfeach::Change_TxidKeyHASH_backTO_AssmIDkeyHASH ($In_GBK_TaxIDKEY_HASH);
	#my $In_RSF_AsmbIDKEY_HASH=GenomeDATAfeach::Change_TxidKeyHASH_backTO_AssmIDkeyHASH ($In_RSF_TaxIDKEY_HASH);
	
	
	
}

##################################
# 融合 GenBank 和 RefSeq 函数组  # 2
################################## 

#这个没有写完
sub Change_TxidKeyHASH_backTO_AssmIDkeyHASH{  #  my $AssemblyID_hash=GenomeDATAfeach::Change_TxidKeyHASH_backTO_AssmIDkeyHASH ($inTxidKeyHASH);
	#将以Taxid为key的hash，变回以 assembly ID为Key的Hash。
  my ($inTxidKeyHASH)=@_;
  
  my $warnMsgBody="\nIn package  GenomeDATAfeach,\tIn sub Change_TxidKeyHASH_backTO_AssmIDkeyHASH,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody"; 
  my $subCallereIfm=DirFileHandle::print_SubCallerInform; #warn $subCallereIfm;
  	
	if  (   (  defined ( $inTxidKeyHASH )  ) && (  ref ( $inTxidKeyHASH ) eq 'HASH'  )   ){}
  else{  	DieWork::Just_dieWork( $die_MsgHead."\n \$inTxidKeyHASH=$inTxidKeyHASH\n should be a HASH ref  $!\n\n\n".$subCallereIfm ); 	  }
  
  my $AssemblyID_hash;
  foreach my $eachTxid (    sort  (   keys (  %{ $inTxidKeyHASH } )   )    ){
		if  (   (  defined ( $inTxidKeyHASH->{$eachTxid} )  ) && (  ref ( $inTxidKeyHASH->{$eachTxid} ) eq 'ARRAY'  )   ){
			foreach my $eachAssemblyHASH (  @ { $inTxidKeyHASH->{$eachTxid} }  ){ 
				if  (   (  defined ( $eachAssemblyHASH )  ) && (  ref ( $eachAssemblyHASH ) eq 'HASH'  ) && (  defined ( $eachAssemblyHASH->{'assembly_accession'} )  ) && ( $eachAssemblyHASH->{'assembly_accession'} =~ m/\S+/)   ){
				  $AssemblyID_hash->{ $eachAssemblyHASH->{'assembly_accession'} }=Storable::dclone( $eachAssemblyHASH );
				}
				else{
					DieWork::Just_dieWork( $die_MsgHead."\n \$eachAssemblyHASH=$eachAssemblyHASH \$eachAssemblyHASH->{'assembly_accession'}=$eachAssemblyHASH->{'assembly_accession'}\n should be all defined and as HASH ref and right asssemby id  $!\n\n\n".$subCallereIfm );
				}
			}
		}
	}
  
  return $AssemblyID_hash;
}


##################################
# 融合 GenBank 和 RefSeq 函数组  # 3
################################## 

#####这个没有写完
sub Add_gbrs_paired_inform{
	my ($in_out_AddedHASH, $in_TxidKeyHASH, $GBK_or_RSF, $in_Parid_assemb_HASH)=@_;
	
	my $warnMsgBody="\nIn package  GenomeDATAfeach,\tIn sub Add_gbrs_paired_inform,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody"; 
  my $subCallereIfm=DirFileHandle::print_SubCallerInform; #warn $subCallereIfm;
  
	if  (   (  defined ( $in_TxidKeyHASH )  ) && (  ref ( $in_TxidKeyHASH ) eq 'HASH'  )   ){}
  else{  	DieWork::Just_dieWork( $die_MsgHead."\n \$in_TxidKeyHASH=$in_TxidKeyHASH\n should be a HASH ref  $!\n\n\n".$subCallereIfm ); 	  }
  if  (   (  defined ( $in_Parid_assemb_HASH )  ) && (  ref ( $in_Parid_assemb_HASH ) eq 'HASH'  )   ){}
  else{  	DieWork::Just_dieWork( $die_MsgHead."\n \$in_Parid_assemb_HASH=$in_Parid_assemb_HASH\n should be a HASH ref  $!\n\n\n".$subCallereIfm ); 	  }
	
	my $GBK_RSF_KEY_1='0_0_0_GenBank';
	my $GBK_RSF_KEY_2='0_0_1__RefSeq';
	if (   (  defined ( $GBK_or_RSF )  ) && ( $GBK_or_RSF=~m/\S+/ ) && (  ( $GBK_or_RSF eq 'GBK' ) || ( $GBK_or_RSF eq 'RSF' )  )   ){
		if ( $GBK_or_RSF eq 'GBK' ){					}
		else{
			$GBK_RSF_KEY_1='0_0_1__RefSeq';
			$GBK_RSF_KEY_2='0_0_0_GenBank';
		}
	}
	else{  	DieWork::Just_dieWork( $die_MsgHead."\n \$GBK_or_RSF=$GBK_or_RSF\n should be defined and as GBK or RSF  $!\n\n\n".$subCallereIfm ); 	  }
	
	my $AssemblyID_hash;	
  foreach my $eachTxid (    sort  (   keys (  %{ $in_TxidKeyHASH } )   )    ){
		if  (   (  defined ( $in_TxidKeyHASH->{$eachTxid} )  ) && (  ref ( $in_TxidKeyHASH->{$eachTxid} ) eq 'ARRAY'  )   ){
			foreach my $eachAssemblyHASH (  @ { $in_TxidKeyHASH->{$eachTxid} }  ){ 
				if  (   (  defined ( $eachAssemblyHASH )  ) && (  ref ( $eachAssemblyHASH ) eq 'HASH'  )   ){
					#####################
					#
				  if (   (  defined ( $eachAssemblyHASH->{'gbrs_paired_asm'} )  ) && ( $eachAssemblyHASH->{'gbrs_paired_asm'} =~ m/\S+/)   ){
				    $AssemblyID_hash->{ $eachAssemblyHASH->{'assembly_accession'} }=Storable::dclone( $eachAssemblyHASH );
				  }
				}
				else{
					DieWork::Just_dieWork( $die_MsgHead."\n \$eachAssemblyHASH=$eachAssemblyHASH should be  defined and as HASH ref  $!\n\n\n".$subCallereIfm );
				}
			}
		}
	}
  
  return $AssemblyID_hash;
	
	
}






##################################
#                                #
#   NCBI taxnomy信息处理 函数    # 
#                                #
##################################

##################################
#   NCBI taxnomy信息处理 函数    #1
################################## 

#调用下面的CollectAllTaxnomyInform_from_NCBI两次，对GCA和GCF的taxnomy id为key的hash都进行 linage的hash的建立
#没有函数输入输出，但有4个 全局变量: 
#第一对全局变量，是 Genbank（未注释）的assembly的hash（key是taxonomy的id） 
#   $JL_GenomeWork_assembly_GenBankSummary_TaxidKey_HASH, $JL_GenomeWork_assembly_GenBank_TaxidKey_Linage_HASH 
#第二对全局变量，是 Ref  （已经注释）的assembly的hash（key是taxonomy的id） 
#   $JL_GenomeWork_assembly_refSeqSummary_TaxidKey_HASH,  $JL_GenomeWork_assembly_refSeq_TaxidKey_Linage_HASH 
sub Build_all_taxid_to_linage_hash_for_Genebank_and_refSeq_File{  #GenomeDATAfeach::Build_all_taxid_to_linage_hash_for_Genebank_and_refSeq_File();
	GenomeDATAfeach::Build_all_taxid_to_linage_hash( $JL_GenomeWork_assembly_GenBankSummary_TaxidKey_HASH, $JL_GenomeWork_assembly_GenBank_TaxidKey_Linage_HASH);
	GenomeDATAfeach::Build_all_taxid_to_linage_hash( $JL_GenomeWork_assembly_refSeqSummary_TaxidKey_HASH,  $JL_GenomeWork_assembly_refSeq_TaxidKey_Linage_HASH );
}

##################################
#   NCBI taxnomy信息处理 函数    #2
################################## 
#      核心 运算 函数            #1
##################################

##输入的是 以 taxonomy id为key的hash的dump文件，然后调用ArrayHashChange::Change_Hash_to_Array 将其转换为ARRAY， 
#调用 TaxonomyWork_NEW::Build_LineageHASH_for_quick_search获得这个ARRAY中所有taxanomy id对应的linage信息
#输出 以 taxonomy id为key的包含linage的hash的dump文件，
sub Build_all_taxid_to_linage_hash{  #GenomeDATAfeach::Build_all_taxid_to_linage_hash( $in_Summary_Taxid_Key_hash_file, $out_Summary_Taxid_Key_linage_hash_file );
	my ( $in_Summary_Taxid_Key_hash_file, $out_Summary_Taxid_Key_linage_hash_file )=@_;
	
	my $warnMsgBody="\nIn package  GenomeDATAfeach,\tIn sub Build_all_taxid_to_linage_hash,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $subCallereIfm=DirFileHandle::print_SubCallerInform;
	
	my $bigSearchHASH;
	if (   (  defined ( $in_Summary_Taxid_Key_hash_file )  ) && ( $in_Summary_Taxid_Key_hash_file=~m/\S+/ ) && (  -e ( $in_Summary_Taxid_Key_hash_file )  )   ){
		my $in_Summary_Taxid_Key_hash=Storable::retrieve( $in_Summary_Taxid_Key_hash_file );
		if (   (  defined ( $in_Summary_Taxid_Key_hash  )  ) && (  ref ( $in_Summary_Taxid_Key_hash  ) eq 'HASH'  )   ){
			my $taxid_search_array=ArrayHashChange::Change_Hash_to_Array ($in_Summary_Taxid_Key_hash);
      $bigSearchHASH=TaxonomyWork_NEW::Build_LineageHASH_for_quick_search( $taxid_search_array );
      if  (   (  defined ( $bigSearchHASH  )  ) && (  ref ( $bigSearchHASH  ) eq 'HASH'  )   ){
      	if (   (  defined ( $out_Summary_Taxid_Key_linage_hash_file )  ) && ( $out_Summary_Taxid_Key_linage_hash_file=~m/\S+/ )   ){
      		DirFileHandle::PrintDumper( $out_Summary_Taxid_Key_linage_hash_file,  $bigSearchHASH   ) ;  
      	}
      }
      else{
      	DieWork::Just_dieWork( $die_MsgHead."\n \$bigSearchHASH=$bigSearchHASH\n should be a Hash ref  $!\n\n\n".$subCallereIfm ); 	
      }
      
		}
		else{
		  DieWork::Just_dieWork( $die_MsgHead."\n \$in_Summary_Taxid_Key_hash=$in_Summary_Taxid_Key_hash\n should be a Hash ref  $!\n\n\n".$subCallereIfm ); 	
	  }
	}
	else{
		DieWork::Just_dieWork( $die_MsgHead."\n \$in_Summary_Taxid_Key_hash_file=$in_Summary_Taxid_Key_hash_file\n should be a readable file  $!\n\n\n".$subCallereIfm ); 	
	}
	return $bigSearchHASH;
}


##################################
#   NCBI taxnomy信息处理 函数    #3
################################## 


#调用下面的CollectAllTaxnomyInform_from_NCBI两次，对GCA和GCF的taxnomy id为key的hash都进行 ncbi信息的搜索，并更新服务器的taxnomy hash
#没有函数输入输出，但有2个 全局变量: 
#第一个全局变量，是 Genbank（未注释）的assembly的hash（key是taxonomy id） 
#   $JL_GenomeWork_assembly_GenBankSummary_TaxidKey_HASH
#第二个全局变量，是 Ref  （已经注释）的assembly的hash（key是taxonomy id） 
#   $JL_GenomeWork_assembly_refSeqSummary_TaxidKey_HASH
sub CollectAllTaxInform_for_Genebank_and_refSeq_File{  #GenomeDATAfeach::CollectAllTaxInform_for_Genebank_and_refSeq_File();
	GenomeDATAfeach::CollectAllTaxnomyInform_from_NCBI( $JL_GenomeWork_assembly_GenBankSummary_TaxidKey_HASH);
	GenomeDATAfeach::CollectAllTaxnomyInform_from_NCBI( $JL_GenomeWork_assembly_refSeqSummary_TaxidKey_HASH );
}

##################################
#   NCBI taxnomy信息处理 函数    #4
################################## 
#      核心 运算 函数            #2  
##################################

#输入的是 以 taxonomy id为key的hash的dump文件，然后调用ArrayHashChange::Change_Hash_to_Array 将其转换为ARRAY， 
#再调用TaxonomyWork_NEW::build_JL_taxonomy_HASH_from_TaxidArray获得这个ARRAY中所有taxanomy id对应的信息信息，集成到 服务器中的taxnomy信息hash中
#看起来没有输出，实际上 已经完成了 assembly中所有 taxnomy id的详细信息的获取
sub CollectAllTaxnomyInform_from_NCBI{  #GenomeDATAfeach::CollectAllTaxnomyInform_from_NCBI( $in_Summary_Taxid_Key_hash_file );
	my ( $in_Summary_Taxid_Key_hash_file )=@_;
	
	my $warnMsgBody="\nIn package  GenomeDATAfeach,\tIn sub CollectAllTaxnomyInform_from_NCBI,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $subCallereIfm=DirFileHandle::print_SubCallerInform;
	
	
	if (   (  defined ( $in_Summary_Taxid_Key_hash_file )  ) && ( $in_Summary_Taxid_Key_hash_file=~m/\S+/ ) && (  -e ( $in_Summary_Taxid_Key_hash_file )  )   ){
		my $in_Summary_Taxid_Key_hash=Storable::retrieve( $in_Summary_Taxid_Key_hash_file );
		if (   (  defined ( $in_Summary_Taxid_Key_hash  )  ) && (  ref ( $in_Summary_Taxid_Key_hash  ) eq 'HASH'  )   ){
			my $taxid_search_array=ArrayHashChange::Change_Hash_to_Array ($in_Summary_Taxid_Key_hash);
      TaxonomyWork_NEW::build_JL_taxonomy_HASH_from_TaxidArray( $taxid_search_array );
		}
		else{
		  DieWork::Just_dieWork( $die_MsgHead."\n \$in_Summary_Taxid_Key_hash=$in_Summary_Taxid_Key_hash\n should be a Hash ref  $!\n\n\n".$subCallereIfm ); 	
	  }
	}
	else{
		DieWork::Just_dieWork( $die_MsgHead."\n \$in_Summary_Taxid_Key_hash_file=$in_Summary_Taxid_Key_hash_file\n should be a readable file  $!\n\n\n".$subCallereIfm ); 	
	}
	
}





##################################
#                                #
# summary文件 转换为hash的函数组 # 
#                                #
##################################

##################################
# summary文件 转换为hash的函数组 # 1
################################## 
#下面这个函数，是调用BuildHashFileForSummaryFile，对rsf（已注释）和genbank（未注释）基因组的summary文件分别进行 hash的转化操作
sub BuildSummaryHASH{  #GenomeDATAfeach::BuildSummaryHASH();
	GenomeDATAfeach::BuildHashFileForSummaryFile($JL_GenomeWork_assembly_refSeqSummary,   $JL_GenomeWork_assembly_refSeqSummary_HASH );	 warn "20190401-0-1-0 \$JL_GenomeWork_assembly_refSeqSummary=$JL_GenomeWork_assembly_refSeqSummary \$JL_GenomeWork_assembly_refSeqSummary_HASH=$JL_GenomeWork_assembly_refSeqSummary_HASH\n";
	GenomeDATAfeach::BuildHashFileForSummaryFile($JL_GenomeWork_assembly_GenBankSummary,  $JL_GenomeWork_assembly_GenBankSummary_HASH);	 warn "20190401-0-1-1 \$JL_GenomeWork_assembly_GenBankSummary=$JL_GenomeWork_assembly_GenBankSummary \$JL_GenomeWork_assembly_GenBankSummary_HASH=$JL_GenomeWork_assembly_GenBankSummary_HASH\n";
}


##################################
# summary文件 转换为hash的函数组 # 2
##################################
#输入 文本的通常是csv或其它表格形式的文件，如assembly_summary_refseq.txt，转换为hash文件，key是 #assembly_accession，如 GCF_000001215.4 GCA_000001215.4
#本函数的另一个输入是 输出的hash文件 ，调用了 GetHash_fromSummaryFile 文件
sub BuildHashFileForSummaryFile{ #GenomeDATAfeach::BuildHashFileForSummaryFile($inSummaryFile,  $inSummaryHASH_File);	
	my ($inSummaryFile,  $inSummaryHASH_File)=@_;
	
	my $warnMsgBody="\nIn package  GenomeDATAfeach,\tIn sub BuildHashFileForSummaryFile,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $subCallereIfm=DirFileHandle::print_SubCallerInform;
	
	my $SummaryHash=GenomeDATAfeach::GetHash_fromSummaryFile($inSummaryFile);	  warn "20190401-0-0-0 \$inSummaryHASH_File=$inSummaryHASH_File \$SummaryHash=$SummaryHash\n";
	
	if (   (  defined ( $SummaryHash  )  ) && (  ref ( $SummaryHash  ) eq 'HASH'  ) && (  defined ( $inSummaryHASH_File )  ) && ( $inSummaryHASH_File=~m/\S+/ )    ){
		DirFileHandle::PrintDumper( $inSummaryHASH_File,  $SummaryHash   ) ;  warn "20190401-0-0-1 \$inSummaryHASH_File=$inSummaryHASH_File \$SummaryHash=$SummaryHash\n";
	}   
	else{
		DieWork::Just_dieWork( $die_MsgHead."\n \$inSummaryHASH_File=$inSummaryHASH_File \$SummaryHash=$SummaryHash\n Something wrong here  $!\n\n\n".$subCallereIfm ); 	
	}
	
}


##################################
# summary文件 转换为hash的函数组 # 3
##################################
# 输入 summary文件，生成hash, 调用了 GenomeAssembly_report_to_hash函数，并以0列 作为summury文件的hash的key的所在列
sub GetHash_fromSummaryFile{  #$SummaryHash=GenomeDATAfeach::GetHash_fromSummaryFile($inSummaryFile);	
	my ($inSummaryFile)=@_;
  
  my $warnMsgBody="\nIn package  GenomeDATAfeach,\tIn sub GetHash_fromSummaryFile,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $subCallereIfm=DirFileHandle::print_SubCallerInform;
	
	my $SummaryHash;
	if (   (  defined ( $inSummaryFile )  ) && ( $inSummaryFile=~m/\S+/ ) && (  -e ( $inSummaryFile )  )   ){	  
	  $SummaryHash=GenomeDATAfeach::GenomeAssembly_report_to_hash($inSummaryFile, 0);	
	}
	else{
		DieWork::Just_dieWork( $die_MsgHead."\n \$inSummaryFile=$inSummaryFile\n should  be defined and exist  $!\n\n\n".$subCallereIfm ); 	
	}	
	return $SummaryHash;	
}

##################################
# summary文件 转换为hash的函数组 # 4
##################################
#      核心 运算 函数            #  
##################################
# 输入 summary文件，以及特定的列列号，生成hash
sub GenomeAssembly_report_to_hash{   # my $outHash=GenomeDATAfeach::GenomeAssembly_report_to_hash($inFile, $keyColNb);  #输入 特定格式的 table文件，解析变成2维hash，1D的key是指定输入列的值，2D的key是各个关键字字段名
	my ($inFile, $keyColNb)=@_;
  
  my $warnMsgBody="\nIn package  GenomeDATAfeach,\tIn sub GenomeAssembly_report_to_hash,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $subCallereIfm=DirFileHandle::print_SubCallerInform;
  
  
  
  open (IN,$inFile)or DieWork::Just_dieWork( $die_MsgHead."\n\n\cannot open \$inFile=$inFile : $!\n\n\n".$subCallereIfm );
  	         
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
  	  	  	  #die "\n\nIn package InFileHandle, in sub FormTableToHash\n\$inFile=$inFile, The key line and explanning line should be on the top, please check the file: \n$inFile\n\n\n";  	  	  	
  	  	      DieWork::Just_dieWork( $die_MsgHead."\n\n\$inFile=$inFile, The key line and explanning line should be on the top, please check the file: \n$inFile\n\n\n".$subCallereIfm );
	              
  	  	  	}
  	  	  	if ($eachCell=~/^# \S+/){    #Regular NO2 !
  	  	  	  #print "\n\nIn package InFileHandle, in sub FormTableToHash\nRegular NO2 !\t\t\$eachCell=$eachCell\n\$lineHere=$lineHere\n\n";
  	  	  	  $isThisLineKeyWordLine=1;         #说明此行 是关键字行
  	  	  	  $HowManyKeyWorkLineFound++;
  	  	  	  $eachCell=~s/^#\s*//;
  	  	    }
  	  	    else {
  	  	      $isThisLineISnotCount=1;          #说明此行 是解释性语句
  	  	    }
  	  	  }
  	  	  else {
  	  	  	#因为关键字行 或 解释性语句 应该放在最上面，所以如果输入的文件 不符合这个格式的要求，也就是说 语句运行到这里的时候，$HowManyKeyWorkLineFound==0 则出错
  	  	  	if ($HowManyKeyWorkLineFound==0){
  	  	  	  #die "\n\nIn package InFileHandle, in sub FormTableToHash\n\$inFile=$inFile, The key line and explanning line should be on the top, please check the file: \n$inFile\n\n\n";  	  	  	
  	  	  	  DieWork::Just_dieWork( $die_MsgHead."\n\n\$inFile=$inFile, The key line and explanning line should be on the top, please check the file: \n$inFile\n\n\n".$subCallereIfm );
	            
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
  	  		  #die "\n\nIn package InFileHandle, in sub FormTableToHash\n\$inFile=$inFile,\n\$keyColNb=$keyColNb\nID:        \$tempAr[\$keyColNb]=\$tempAr[$keyColNb]=$tempAr[$keyColNb] is not unique!!  or KeyWord:\$keyhash->{\$ColNb}=\$keyhash->{$ColNb}=$keyhash->{$ColNb} is not unique!!\n\nPlease check the file: $inFile\n\n\n";
  	  	    DieWork::Just_dieWork( $die_MsgHead."\n\n\$inFile=$inFile,\n\$keyColNb=$keyColNb\nID:        \$tempAr[\$keyColNb]=\$tempAr[$keyColNb]=$tempAr[$keyColNb] is not unique!!  or KeyWord:\$keyhash->{\$ColNb}=\$keyhash->{$ColNb}=$keyhash->{$ColNb} is not unique!!\n\nPlease check the file: $inFile\n\n\n".$subCallereIfm );
	            
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


##################################
# summary文件 转换为hash的函数组 # 5
################################## 

#本函数，获得GCA GCF的以taxnomy id为key的两个hash
#没有函数输入输出，但有4个 全局变量: 
#第一对全局变量，是 Genbank（未注释）的assembly的hash（key是GCA形式的id） 
#   $JL_GenomeWork_assembly_GenBankSummary_HASH, $JL_GenomeWork_assembly_GenBankSummary_TaxidKey_HASH 
#第二对全局变量，是 Ref  （已经注释）的assembly的hash（key是GCF形式的id） 
#   $JL_GenomeWork_assembly_GenBankSummary_HASH, $JL_GenomeWork_assembly_GenBankSummary_TaxidKey_HASH 
sub GetTxidHash_for_Genebank_and_refSeq_File{  #GenomeDATAfeach::GetTxidHash_for_Genebank_and_refSeq_File();
	GenomeDATAfeach::GetTaxidHASH( $JL_GenomeWork_assembly_GenBankSummary_HASH, $JL_GenomeWork_assembly_GenBankSummary_TaxidKey_HASH);
	GenomeDATAfeach::GetTaxidHASH( $JL_GenomeWork_assembly_refSeqSummary_HASH,  $JL_GenomeWork_assembly_refSeqSummary_TaxidKey_HASH );
}

##################################
# summary文件 转换为hash的函数组 # 6
################################## 

#输入 转换为hash的assembly的summary文件，key是 #assembly_accession，如 GCF_000001215.4 GCA_000001215.4
#本函数  进一步 将以上hash转换为 taxnomy ID （24, 7953）等数字式编号 表示的 物种或分类学id 为key的 hash，其第二层是个ARRAY，因为同一个id可以对应于 多个 assembly accession（如多个同物种的不同 株）
#调用了 ArrayHashChange::Build_HashArray_from_reversed_HASH
#本函数的另一个输入是输出的hash的dump文件名
sub GetTaxidHASH{  # GenomeDATAfeach::GetTaxidHASH( $in_Summary_AssemblyID_KEY_hash_file, $out_Summary_Taxid_Key_hash_file);
	my ( $in_Summary_AssemblyID_KEY_hash_file, $out_Summary_Taxid_Key_hash_file)=@_;
	
	my $warnMsgBody="\nIn package  GenomeDATAfeach,\tIn sub GetTaxidHASH,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $subCallereIfm=DirFileHandle::print_SubCallerInform;
	
	if (   (  defined ( $in_Summary_AssemblyID_KEY_hash_file )  ) && ( $in_Summary_AssemblyID_KEY_hash_file=~m/\S+/ ) && (  -e ( $in_Summary_AssemblyID_KEY_hash_file )  )   ){
		my $in_Summary_AssemblyID_KEY_hash=Storable::retrieve( $in_Summary_AssemblyID_KEY_hash_file );
		if (   (  defined ( $in_Summary_AssemblyID_KEY_hash  )  ) && (  ref ( $in_Summary_AssemblyID_KEY_hash  ) eq 'HASH'  )   ){
			my $out_Summary_Taxid_Key_hash=ArrayHashChange::Build_HashArray_from_reversed_HASH( $in_Summary_AssemblyID_KEY_hash, 'taxid' );
			if (   (  defined ( $out_Summary_Taxid_Key_hash  )  ) && (  ref ( $out_Summary_Taxid_Key_hash  ) eq 'HASH'  )   ){
			  DirFileHandle::PrintDumper ( $out_Summary_Taxid_Key_hash_file,  $out_Summary_Taxid_Key_hash );	
			}
			else{
		    DieWork::Just_dieWork( $die_MsgHead."\n \$out_Summary_Taxid_Key_hash=$out_Summary_Taxid_Key_hash\n should be good Hash ref  $!\n\n\n".$subCallereIfm ); 	
	    }
		}
		else{
		  DieWork::Just_dieWork( $die_MsgHead."\n \$in_Summary_AssemblyID_KEY_hash=$in_Summary_AssemblyID_KEY_hash\n should be a Hash ref  $!\n\n\n".$subCallereIfm ); 	
	  }
	}
	else{
		DieWork::Just_dieWork( $die_MsgHead."\n \$in_Summary_AssemblyID_KEY_hash_file=$in_Summary_AssemblyID_KEY_hash_file\n should be a readable file  $!\n\n\n".$subCallereIfm ); 	
	}
	
}








##################################
#                                #
# report文件 转换为hash的函数组  # 
# 没有ftp地址，故使用较少        #
#                                #
##################################

##################################
# report文件 转换为hash的函数组  # 1
################################## 
#专门针对 eukarotic 的report文件转hash
#不用输入，但潜在输入是 本模块中的全局变量 $JL_GenomeWork_reports_euakryotes
#调用了 Extract_genome_hash_from_reportFIle
#有输出， 输出是解析后的hash，key是 物种名，但因为一个物种名可能对应多个 taxid，所以第二层是一个ARRAY
sub Extract_eukarotic_GenomeInfomHash{   #  my $GenomeInfomrHash= GenomeDATAfeach::Extract_genome_hash_from_reportFIle();
	#生成 真核生物基因组的hash
	
	my $warnMsgBody="\nIn package  GenomeDATAfeach,\tIn sub Extract_eukarotic_GenomeInfomHash,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $subCallereIfm=DirFileHandle::print_SubCallerInform;
	
	my $GenomeInfomrHash=GenomeDATAfeach::Extract_genome_hash_from_reportFIle($JL_GenomeWork_reports_euakryotes);
	
	return $GenomeInfomrHash;
	
}



##################################
# report文件 转换为hash的函数组  # 2
################################## 
#专门针对 eukarotic 的report文件转hash
#不用输入，但潜在输入是 本模块中的全局变量 $JL_GenomeWork_reports_euakryotes
#调用了 Extract_genome_hash_from_reportFIle
#没有输出，但潜在的输出是 本模块中的全局变量 JL_GenomeWork_reports_euakryotes_HASH
sub Extract_eukarotic_GenomeInfomHash_and_dumperIT{   #  GenomeDATAfeach::Extract_eukarotic_GenomeInfomHash_and_dumperIT();
	#在特定位置，生成 真核生物基因组的hash
	my $warnMsgBody="\nIn package  GenomeDATAfeach,\tIn sub Extract_eukarotic_GenomeInfomHash_and_dumperIT,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $subCallereIfm=DirFileHandle::print_SubCallerInform;
	
	my $GenomeInfomrHash=GenomeDATAfeach::Extract_genome_hash_from_reportFIle( $JL_GenomeWork_reports_euakryotes );
	
	if (   (  defined ( $GenomeInfomrHash  )  ) && (  ref ( $GenomeInfomrHash  ) eq 'HASH'  )    ){
		DirFileHandle::PrintDumper( $JL_GenomeWork_reports_euakryotes_HASH,  $GenomeInfomrHash   ) ;
	}   
	else{
		DieWork::Just_dieWork( $die_MsgHead."\n \$JL_GenomeWork_reports_euakryotes_HASH=$JL_GenomeWork_reports_euakryotes_HASH \$GenomeInfomrHash=$GenomeInfomrHash\n Something wrong here  $!\n\n\n".$subCallereIfm ); 	
	}
	
}

##################################
# report文件 转换为hash的函数组  # 3
################################## 
#      核心 运算 函数            #  
##################################
#这个是从report 文件，中获取相关信息，并转换为hash的，但report文件 没有基因组在ncbi上的存储位置信息，不能用于后续的下载操作
#输入是任意 report文件，输出是解析后的hash，key是 物种名，但因为一个物种名可能对应多个 taxid，所以第二层是一个ARRAY
sub Extract_genome_hash_from_reportFIle{   #  my $GenomeInfomrHash=GenomeDATAfeach::Extract_genome_hash_from_reportFIle($inReportFile);
	#  从某个genome的report文件中 获取全信息HASH
	my ($inReportFile)=@_;
	
	my $warnMsgBody="\nIn package  GenomeDATAfeach,\tIn sub Extract_genome_hash_from_reportFIle,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $subCallereIfm=DirFileHandle::print_SubCallerInform;
	
	my $GenomeInfomrHash;
	if (   (  defined ( $inReportFile )  ) && ( $inReportFile=~m/\S+/ ) && (  -e ( $inReportFile )  )   ){
	  
	  $GenomeInfomrHash=InFileHandle::ChangeIntoHashForTxt($inReportFile);
	}
	else{
		DieWork::Just_dieWork( $die_MsgHead."\n \$inReportFile=$inReportFile\n should  be defined and exist  $!\n\n\n".$subCallereIfm ); 	
	}
	
	return $GenomeInfomrHash;
	
}


sub GetGenomeSequence{
	my ($inTaxidHASH, $old_downloaded_HASH )=@_;
	
	my ($warnMsgBody, $warnMsgHead, $die_MsgHead, $caller_inform)=@{ DieWork::BuildWarnDieHeadInform( 'GenomeDATAfeach', 'GetGenomeSequence' ) };
	
	DieWork::Check_Hash_or_DIE   ( $inTaxidHASH,       "\$inTaxidHASH",       $die_MsgHead, $caller_inform  );  
	
	foreach my $taxID_here (    sort { $a <=> $b } (   keys (  %{ $inTaxidHASH } )   )    ){
		if ( DieWork::Check_Hash_or_NOT( $inTaxidHASH->{$taxID_here} ) ){
			
		}
	}
	
	
}

1;

##########################################################################################################################################
# 
