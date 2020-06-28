
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

#NCBI �й� ���������ݵ� ���ص�ַ����Ϣ ���� ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/


#NCBI �й� ���������ݵ� ���ص�ַ����Ϣ ���� ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/

my $assembly_gbk_url="ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt";
my $assembly_ref_url="ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt";

#�ɷ������ϵ� �ϲ��ļ��е�ַ
my $upLevelDIR                                          ="/home/fredjiang/OneTSSD/fredjiang20190321";

#�·������ϵ� �������ݵĵ�ַ
#���� �����ļ��� ���ڸ� �ļ�����
my $JL_GenomeWork_Dir                                   ="/mnt/md1/datahub";
#my $JL_GenomeWork_Dir                                   ="/home/fredjiang/md1/fredjiang20190916/GenomeData20190319";

#   reprot���� �ļ���

my $JL_all_assemblyDownload_DIR                         =$upLevelDIR."/GENOME_REPORTS";
my $JL_difVersion_DIR                                   =$JL_all_assemblyDownload_DIR."/G_REPORTS_20191227";  # ÿ�θ��� ��������ļ��е�����

my $JL_GenomeWork_assembly_report_Dir                   =$JL_difVersion_DIR."/ASSEMBLY_REPORTS";



#���� �� assembly�ļ� ������������Ϣ����Щhash�ļ������ݣ��������assembly�ļ���������ncbi���������ġ�
#���ز���ʱҲ������ ��Щ�ļ����������غ͸��µ�
my $JL_GenomeWork_assembly_GenBankSummary               =$JL_GenomeWork_assembly_report_Dir."/1m0m_assembly_summary_genbank.txt";                             #ԭʼ�ļ���ncbi����  δע�ͻ�������Ϣ
my $JL_GenomeWork_assembly_GenBankSummary_HASH          =$JL_GenomeWork_assembly_report_Dir."/1m1m_assembly_summary_genbank.txt.hash";
my $JL_GenomeWork_assembly_GenBankSummary_TaxidKey_HASH =$JL_GenomeWork_assembly_report_Dir."/1m2m_assembly_summary_genbank.txt.TaxidKey.hash";
my $JL_GenomeWork_assembly_GenBank_TaxidKey_Linage_HASH =$JL_GenomeWork_assembly_report_Dir."/1m3m_assembly_summary_genbank.txt.TaxidKey.lineage.hash";

my $JL_GenomeWork_assembly_refSeqSummary                =$JL_GenomeWork_assembly_report_Dir."/2m0m_assembly_summary_refseq.txt";                              #ԭʼ�ļ���ncbi����  ��ע�ͻ�������Ϣ
my $JL_GenomeWork_assembly_refSeqSummary_HASH           =$JL_GenomeWork_assembly_report_Dir."/2m1m_assembly_summary_refseq.txt.hash";
my $JL_GenomeWork_assembly_refSeqSummary_TaxidKey_HASH  =$JL_GenomeWork_assembly_report_Dir."/2m2m_assembly_summary_refseq.txt.TaxidKey.hash";
my $JL_GenomeWork_assembly_refSeq_TaxidKey_Linage_HASH  =$JL_GenomeWork_assembly_report_Dir."/2m3m_assembly_summary_refseq.txt.TaxidKey.lineage.hash";


#�����Ǹ��� ��������genome���ļ��ĵ�ַ�� hash����Щ�ļ����������ع��� ���ɵġ�
my $JL_GenomeWork_assembly_download_Eukayotic_HASH      =$JL_GenomeWork_assembly_report_Dir."/3m1m_assembly_summary_eukaryotes_data.txt";
my $JL_GenomeWork_assembly_download_Bacterial_HASH      =$JL_GenomeWork_assembly_report_Dir."/3m2m_assembly_summary_Bacteria___data.txt";


######�����reports ���������û��ftp��Ϣ�� ʹ��Ƶ�ʲ�����
my $JL_GenomeWork_reports_Dir                           =$JL_difVersion_DIR."/GENOME_REPORTS";
my $JL_GenomeWork_reports_euakryotes                    =$JL_GenomeWork_reports_Dir."/eukaryotes.txt";
my $JL_GenomeWork_reports_euakryotes_HASH               =$JL_GenomeWork_reports_Dir."/eukaryotes.txt.hash";



#############################################################
#��������ʾ����

  #��0�������� �ļ���$JL_difVersion_DIR =$JL_all_assemblyDownload_DIR."/GENOME_REPORTS_20191227"; ��˫�����ڵ�����
  #Ȼ����������ĺ�����������assembly report�ļ�
#GenomeDATAfeach::wget_assembly_information( );
  #�ٳ��ԱȽ��Ƿ���Ҫ����

  #��1������assembly��summary���ı��ļ� ת��Ϊ �� assembly id������ GCF_000001215.4 GCA_000001215.4��Ϊkey�� hash
#GenomeDATAfeach::BuildSummaryHASH();
  #��2������assembly��summary���ı��ļ� ת��Ϊ ��  taxnomy ID ��24, 7953�� Ϊkey ��hash
#GenomeDATAfeach::GetTxidHash_for_Genebank_and_refSeq_File();
  #��3����������assembly summary���ļ����漰���� taxnomy id������ ncbi�ϻ�ȡxml�ļ������·������ϵ� ������Ϣ���ݿ�
#GenomeDATAfeach::CollectAllTaxInform_for_Genebank_and_refSeq_File();
  #��4��������������ֵ� genbank��rsf���� ��linage��Ϣ
#GenomeDATAfeach::Build_all_taxid_to_linage_hash_for_Genebank_and_refSeq_File();
 
 #��5�������ز����ɴ������ֱ����ַ��Ϣ��hash��
   #GenomeDATAfeach::RSYNC_all_genomeData_for_a_eukaryotes();  #�����������
   #GenomeDATAfeach::RSYNC_all_genomeData_for_Bacteria();      #����ԭ������
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
# RSYNC ���ظ��� ������          # 
#                                #
##################################

##################################
# RSYNC ���ظ��� ������          # 1
################################## 
#��NCBI�������е� Bacteria �Ļ��������ݣ����ظ�ʹ�øú������������ص����ݣ�������Ч���� �����������ݣ�������������
#û�к������룬��ȫ�ֱ��� 'Eukaryota',  �� $JL_GenomeWork_assembly_download_Eukayotic_HASH
sub RSYNC_all_genomeData_for_a_eukaryotes{  #GenomeDATAfeach::RSYNC_all_genomeData_for_a_eukaryotes();
	GenomeDATAfeach::RSYNC_all_genomeData_for_a_phylaName ( 'Eukaryota', $JL_GenomeWork_assembly_download_Eukayotic_HASH );
}

##################################
# RSYNC ���ظ��� ������          # 2
################################## 
#��NCBI�������е� Bacteria �Ļ��������ݣ����ظ�ʹ�øú������������ص����ݣ�������Ч���� �����������ݣ�������������
#û�к������룬��ȫ�ֱ��� 'Eukaryota',  �� $JL_GenomeWork_assembly_download_Eukayotic_HASH
sub RSYNC_all_genomeData_for_Bacteria{  #GenomeDATAfeach::RSYNC_all_genomeData_for_Bacteria();
	GenomeDATAfeach::RSYNC_all_genomeData_for_a_phylaName ( 'Bacteria', $JL_GenomeWork_assembly_download_Bacterial_HASH );
}

##################################
# RSYNC ���ظ��� ������          # 3
################################## 

#  ����1 ������֧������ Cervidae�� ������2 ����������ļ�����������иý�����֧�µĻ����飬�Լ����� �������HASH
#���� ���ص� ȫ�ֱ���
# ȫ�ֱ��� 1 $JL_GenomeWork_assembly_GenBankSummary_TaxidKey_HASH
# ȫ�ֱ��� 2 $JL_GenomeWork_assembly_GenBank_TaxidKey_Linage_HASH
# ȫ�ֱ��� 3 $JL_GenomeWork_assembly_refSeqSummary_TaxidKey_HASH
# ȫ�ֱ��� 4 $JL_GenomeWork_assembly_refSeq_TaxidKey_Linage_HASH
#����� �����µ� ���� �洢��ַ��hash
#������ GenomeDATAfeach::FindAllTxid_for_phlyaName_fromHASHfile ����Ѱ�� phyla_Name�µ���������
#������ ArrayHashChange::Change_Hash_to_Array ��hash��ΪArray
#������ GenomeDATAfeach::RSYNC_remote_GENOME_data_from_Taxid ������ ���غ͸���
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
# RSYNC ���ظ��� ������          # 4
################################## 
#      ���� ���� ����            #  
##################################

#�������Ĺ��ܣ������� taxnomy idΪkey��ARRAY��HASH�е�ftp��ַ����NCBI����rsync�����ȡ��Ӧ���ֵĻ��������ݣ����ص����ص��ض���ַ�������õ�ַ�����������HASH
#����1    ������Ҫ ���ص�taxnomy idΪԪ�ص�ARRAY�����ARRAY���� ����2 ��HASH������ArrayHashChange::Change_Hash_to_Arrayת���������� 
#����2    ������Ҫ ���ص�taxnomy idΪkey ��HASH
#����3    �������� hash��dump�ļ����ļ�·��
#���     ������ �����鱣���ַ�� �µ�hash�ļ�
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
# ��Ѱphyla������ �������ֵĺ��� # 
#                                #
##################################

##################################
# ��Ѱphyla������ �������ֵĺ��� #1
################################## 

#�����������һ������FindAllTxid_for_phlyaName������һ�µģ�ֻ�� �������������� hash��dump�ļ�����һ����������ľ���hash
#����һ�� phyla name���������Ϊ taxnomy����ѧ�� ��һ�� �ڵ������ �� ��׵����veterbrate�����������eukaryotics�ȵȡ�
#��Ҫ���� һ������������taxnomy id��hash�������hash�������е�taxnomy id���������м�⣬������taxnomy id���ǲ��Ǹ�phyla�����Σ��硰Ѽ���ޡ���taxnomy id���ǲ���veterbrate������
#��Ҫ���� һ������������ taxnomy id��linage����Ϣ��hash���Ӹ�hash���ҵ����е� linage��Ϣ
#������� ��phyla name�µ��������ֵ�taxnomy id�γɵ� hash
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
# ��Ѱphyla������ �������ֵĺ��� #2
################################## 
#      ���� ���� ����            #  
##################################

#����һ�� phyla name���������Ϊ taxnomy����ѧ�� ��һ�� �ڵ������ �� ��׵����veterbrate�����������eukaryotics�ȵȡ�
#��Ҫ���� һ������������taxnomy id��hash�������hash�������е�taxnomy id���������м�⣬������taxnomy id���ǲ��Ǹ�phyla�����Σ��硰Ѽ���ޡ���taxnomy id���ǲ���veterbrate������
#��Ҫ���� һ������������ taxnomy id��linage����Ϣ��hash���Ӹ�hash���ҵ����е� linage��Ϣ
#������� ��phyla name�µ��������ֵ�taxnomy id�γɵ� hash
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
# �ں� GenBank �� RefSeq ������  # 
#          û��д��              #
#                                #
##################################

##################################
# �ں� GenBank �� RefSeq ������  # 1
################################## 

#���û��д��
sub Merge_GenBank_RefSeq_TaxID_HASH{
	my ($In_GBK_TaxIDKEY_HASH, $In_RSF_TaxIDKEY_HASH)=@_;
	
	#my $In_GBK_AsmbIDKEY_HASH=GenomeDATAfeach::Change_TxidKeyHASH_backTO_AssmIDkeyHASH ($In_GBK_TaxIDKEY_HASH);
	#my $In_RSF_AsmbIDKEY_HASH=GenomeDATAfeach::Change_TxidKeyHASH_backTO_AssmIDkeyHASH ($In_RSF_TaxIDKEY_HASH);
	
	
	
}

##################################
# �ں� GenBank �� RefSeq ������  # 2
################################## 

#���û��д��
sub Change_TxidKeyHASH_backTO_AssmIDkeyHASH{  #  my $AssemblyID_hash=GenomeDATAfeach::Change_TxidKeyHASH_backTO_AssmIDkeyHASH ($inTxidKeyHASH);
	#����TaxidΪkey��hash������� assembly IDΪKey��Hash��
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
# �ں� GenBank �� RefSeq ������  # 3
################################## 

#####���û��д��
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
#   NCBI taxnomy��Ϣ���� ����    # 
#                                #
##################################

##################################
#   NCBI taxnomy��Ϣ���� ����    #1
################################## 

#���������CollectAllTaxnomyInform_from_NCBI���Σ���GCA��GCF��taxnomy idΪkey��hash������ linage��hash�Ľ���
#û�к����������������4�� ȫ�ֱ���: 
#��һ��ȫ�ֱ������� Genbank��δע�ͣ���assembly��hash��key��taxonomy��id�� 
#   $JL_GenomeWork_assembly_GenBankSummary_TaxidKey_HASH, $JL_GenomeWork_assembly_GenBank_TaxidKey_Linage_HASH 
#�ڶ���ȫ�ֱ������� Ref  ���Ѿ�ע�ͣ���assembly��hash��key��taxonomy��id�� 
#   $JL_GenomeWork_assembly_refSeqSummary_TaxidKey_HASH,  $JL_GenomeWork_assembly_refSeq_TaxidKey_Linage_HASH 
sub Build_all_taxid_to_linage_hash_for_Genebank_and_refSeq_File{  #GenomeDATAfeach::Build_all_taxid_to_linage_hash_for_Genebank_and_refSeq_File();
	GenomeDATAfeach::Build_all_taxid_to_linage_hash( $JL_GenomeWork_assembly_GenBankSummary_TaxidKey_HASH, $JL_GenomeWork_assembly_GenBank_TaxidKey_Linage_HASH);
	GenomeDATAfeach::Build_all_taxid_to_linage_hash( $JL_GenomeWork_assembly_refSeqSummary_TaxidKey_HASH,  $JL_GenomeWork_assembly_refSeq_TaxidKey_Linage_HASH );
}

##################################
#   NCBI taxnomy��Ϣ���� ����    #2
################################## 
#      ���� ���� ����            #1
##################################

##������� �� taxonomy idΪkey��hash��dump�ļ���Ȼ�����ArrayHashChange::Change_Hash_to_Array ����ת��ΪARRAY�� 
#���� TaxonomyWork_NEW::Build_LineageHASH_for_quick_search������ARRAY������taxanomy id��Ӧ��linage��Ϣ
#��� �� taxonomy idΪkey�İ���linage��hash��dump�ļ���
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
#   NCBI taxnomy��Ϣ���� ����    #3
################################## 


#���������CollectAllTaxnomyInform_from_NCBI���Σ���GCA��GCF��taxnomy idΪkey��hash������ ncbi��Ϣ�������������·�������taxnomy hash
#û�к����������������2�� ȫ�ֱ���: 
#��һ��ȫ�ֱ������� Genbank��δע�ͣ���assembly��hash��key��taxonomy id�� 
#   $JL_GenomeWork_assembly_GenBankSummary_TaxidKey_HASH
#�ڶ���ȫ�ֱ������� Ref  ���Ѿ�ע�ͣ���assembly��hash��key��taxonomy id�� 
#   $JL_GenomeWork_assembly_refSeqSummary_TaxidKey_HASH
sub CollectAllTaxInform_for_Genebank_and_refSeq_File{  #GenomeDATAfeach::CollectAllTaxInform_for_Genebank_and_refSeq_File();
	GenomeDATAfeach::CollectAllTaxnomyInform_from_NCBI( $JL_GenomeWork_assembly_GenBankSummary_TaxidKey_HASH);
	GenomeDATAfeach::CollectAllTaxnomyInform_from_NCBI( $JL_GenomeWork_assembly_refSeqSummary_TaxidKey_HASH );
}

##################################
#   NCBI taxnomy��Ϣ���� ����    #4
################################## 
#      ���� ���� ����            #2  
##################################

#������� �� taxonomy idΪkey��hash��dump�ļ���Ȼ�����ArrayHashChange::Change_Hash_to_Array ����ת��ΪARRAY�� 
#�ٵ���TaxonomyWork_NEW::build_JL_taxonomy_HASH_from_TaxidArray������ARRAY������taxanomy id��Ӧ����Ϣ��Ϣ�����ɵ� �������е�taxnomy��Ϣhash��
#������û�������ʵ���� �Ѿ������ assembly������ taxnomy id����ϸ��Ϣ�Ļ�ȡ
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
# summary�ļ� ת��Ϊhash�ĺ����� # 
#                                #
##################################

##################################
# summary�ļ� ת��Ϊhash�ĺ����� # 1
################################## 
#��������������ǵ���BuildHashFileForSummaryFile����rsf����ע�ͣ���genbank��δע�ͣ��������summary�ļ��ֱ���� hash��ת������
sub BuildSummaryHASH{  #GenomeDATAfeach::BuildSummaryHASH();
	GenomeDATAfeach::BuildHashFileForSummaryFile($JL_GenomeWork_assembly_refSeqSummary,   $JL_GenomeWork_assembly_refSeqSummary_HASH );	 warn "20190401-0-1-0 \$JL_GenomeWork_assembly_refSeqSummary=$JL_GenomeWork_assembly_refSeqSummary \$JL_GenomeWork_assembly_refSeqSummary_HASH=$JL_GenomeWork_assembly_refSeqSummary_HASH\n";
	GenomeDATAfeach::BuildHashFileForSummaryFile($JL_GenomeWork_assembly_GenBankSummary,  $JL_GenomeWork_assembly_GenBankSummary_HASH);	 warn "20190401-0-1-1 \$JL_GenomeWork_assembly_GenBankSummary=$JL_GenomeWork_assembly_GenBankSummary \$JL_GenomeWork_assembly_GenBankSummary_HASH=$JL_GenomeWork_assembly_GenBankSummary_HASH\n";
}


##################################
# summary�ļ� ת��Ϊhash�ĺ����� # 2
##################################
#���� �ı���ͨ����csv�����������ʽ���ļ�����assembly_summary_refseq.txt��ת��Ϊhash�ļ���key�� #assembly_accession���� GCF_000001215.4 GCA_000001215.4
#����������һ�������� �����hash�ļ� �������� GetHash_fromSummaryFile �ļ�
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
# summary�ļ� ת��Ϊhash�ĺ����� # 3
##################################
# ���� summary�ļ�������hash, ������ GenomeAssembly_report_to_hash����������0�� ��Ϊsummury�ļ���hash��key��������
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
# summary�ļ� ת��Ϊhash�ĺ����� # 4
##################################
#      ���� ���� ����            #  
##################################
# ���� summary�ļ����Լ��ض������кţ�����hash
sub GenomeAssembly_report_to_hash{   # my $outHash=GenomeDATAfeach::GenomeAssembly_report_to_hash($inFile, $keyColNb);  #���� �ض���ʽ�� table�ļ����������2άhash��1D��key��ָ�������е�ֵ��2D��key�Ǹ����ؼ����ֶ���
	my ($inFile, $keyColNb)=@_;
  
  my $warnMsgBody="\nIn package  GenomeDATAfeach,\tIn sub GenomeAssembly_report_to_hash,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $subCallereIfm=DirFileHandle::print_SubCallerInform;
  
  
  
  open (IN,$inFile)or DieWork::Just_dieWork( $die_MsgHead."\n\n\cannot open \$inFile=$inFile : $!\n\n\n".$subCallereIfm );
  	         
  my $outHash;
  my $keyhash;
  
  #my $rowKeyword="#Organism\/Name";    #����� ������Ϊ Ψһ���ظ���key�� �ؼ���
 
  my $HowManyValueLineWereFound=0;    #��д��ֵ����  ����дΪ��0������ ��¼���ж��ٸ�������
  my $HowManyKeyWorkLineFound=0;      #�ǹؼ�����    ����дΪ��0������ ��¼���ж��ٸ�������
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
  	  
  	  
  	  my $isThisLineKeyWordLine=0;  #�ǹؼ�����    ����дΪ1
  	  my $isThisLineISnotCount=0;   #�ǽ��������  ����дΪ1
  	  my $isThisLineValueLine=0;    #��д��ֵ����  ����дΪ��0������
  	  
  	  
  	  my $ColNb=0;
  	  foreach my $eachCell (@tempAr){
  	  	
  	  	#��� ��һ��cell������һ��ʱ �ؼ����� ���ǽ�������䡣�ؼ�������Ҫȡ �ؼ��֣������������ ������һ�в�����
  	  	if ($ColNb==0){
  	  	  if ($eachCell=~/^#/){   #Regular NO1 ! #��ͷ��#���������� �ؼ����� �� ���������
  	  	  	#print "\n\nIn package InFileHandle, in sub FormTableToHash\n\Regular NO1 !\t\t$eachCell=$eachCell\n\$lineHere=$lineHere\n\n";
  	  	  	#��Ϊ�ؼ����� �� ��������� Ӧ�÷��������棬�������������ļ� �����������ʽ��Ҫ��Ҳ����˵ ������е������ʱ��$HowManyValueLineWereFound>0 �����
  	  	  	if ($HowManyValueLineWereFound>0){
  	  	  	  #die "\n\nIn package InFileHandle, in sub FormTableToHash\n\$inFile=$inFile, The key line and explanning line should be on the top, please check the file: \n$inFile\n\n\n";  	  	  	
  	  	      DieWork::Just_dieWork( $die_MsgHead."\n\n\$inFile=$inFile, The key line and explanning line should be on the top, please check the file: \n$inFile\n\n\n".$subCallereIfm );
	              
  	  	  	}
  	  	  	if ($eachCell=~/^# \S+/){    #Regular NO2 !
  	  	  	  #print "\n\nIn package InFileHandle, in sub FormTableToHash\nRegular NO2 !\t\t\$eachCell=$eachCell\n\$lineHere=$lineHere\n\n";
  	  	  	  $isThisLineKeyWordLine=1;         #˵������ �ǹؼ�����
  	  	  	  $HowManyKeyWorkLineFound++;
  	  	  	  $eachCell=~s/^#\s*//;
  	  	    }
  	  	    else {
  	  	      $isThisLineISnotCount=1;          #˵������ �ǽ��������
  	  	    }
  	  	  }
  	  	  else {
  	  	  	#��Ϊ�ؼ����� �� ��������� Ӧ�÷��������棬�������������ļ� �����������ʽ��Ҫ��Ҳ����˵ ������е������ʱ��$HowManyKeyWorkLineFound==0 �����
  	  	  	if ($HowManyKeyWorkLineFound==0){
  	  	  	  #die "\n\nIn package InFileHandle, in sub FormTableToHash\n\$inFile=$inFile, The key line and explanning line should be on the top, please check the file: \n$inFile\n\n\n";  	  	  	
  	  	  	  DieWork::Just_dieWork( $die_MsgHead."\n\n\$inFile=$inFile, The key line and explanning line should be on the top, please check the file: \n$inFile\n\n\n".$subCallereIfm );
	            
  	  	  	}
  	  	  	$isThisLineValueLine=1;             #˵������ ��д��ֵ����
  	  	  	$HowManyValueLineWereFound++;       #��¼���� �ж��ٸ� д��ֵ����
  	  	  	
  	  	  }
  	  	}
  	  	
  	  	#�������������еĵ�Ԫ��ֵ���в���
  	  	if ($isThisLineKeyWordLine==1){         #����ǹؼ�����, �򽫸��з���ؼ���hash�У���һ�еĲ���������ֵ�г���֮ǰ���
  	  		$keyhash->{$ColNb}=$eachCell;  print "\$keyhash->{$ColNb}=\$eachCell=$keyhash->{$ColNb}\n";
  	  	}
  	  	elsif ($isThisLineValueLine==1){        #�������ֵ�У�  �򽫸���ֵ���� ���$outHash�У�1������ $keyColNb�����ĵڼ��е� ֵ�� 2������ ���еĹؼ��� �� ֵ�Ǹ��и��е�ֵ
  	  		if (   defined (  $outHash->{ $tempAr[$keyColNb] }->{ $keyhash->{$ColNb} }  )   ){ #�����ֵ�Ѿ����������˵�� $tempAr[$keyColNb] ���� $keyhash->{$ColNb} �������ظ�
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
# summary�ļ� ת��Ϊhash�ĺ����� # 5
################################## 

#�����������GCA GCF����taxnomy idΪkey������hash
#û�к����������������4�� ȫ�ֱ���: 
#��һ��ȫ�ֱ������� Genbank��δע�ͣ���assembly��hash��key��GCA��ʽ��id�� 
#   $JL_GenomeWork_assembly_GenBankSummary_HASH, $JL_GenomeWork_assembly_GenBankSummary_TaxidKey_HASH 
#�ڶ���ȫ�ֱ������� Ref  ���Ѿ�ע�ͣ���assembly��hash��key��GCF��ʽ��id�� 
#   $JL_GenomeWork_assembly_GenBankSummary_HASH, $JL_GenomeWork_assembly_GenBankSummary_TaxidKey_HASH 
sub GetTxidHash_for_Genebank_and_refSeq_File{  #GenomeDATAfeach::GetTxidHash_for_Genebank_and_refSeq_File();
	GenomeDATAfeach::GetTaxidHASH( $JL_GenomeWork_assembly_GenBankSummary_HASH, $JL_GenomeWork_assembly_GenBankSummary_TaxidKey_HASH);
	GenomeDATAfeach::GetTaxidHASH( $JL_GenomeWork_assembly_refSeqSummary_HASH,  $JL_GenomeWork_assembly_refSeqSummary_TaxidKey_HASH );
}

##################################
# summary�ļ� ת��Ϊhash�ĺ����� # 6
################################## 

#���� ת��Ϊhash��assembly��summary�ļ���key�� #assembly_accession���� GCF_000001215.4 GCA_000001215.4
#������  ��һ�� ������hashת��Ϊ taxnomy ID ��24, 7953��������ʽ��� ��ʾ�� ���ֻ����ѧid Ϊkey�� hash����ڶ����Ǹ�ARRAY����Ϊͬһ��id���Զ�Ӧ�� ��� assembly accession������ͬ���ֵĲ�ͬ �꣩
#������ ArrayHashChange::Build_HashArray_from_reversed_HASH
#����������һ�������������hash��dump�ļ���
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
# report�ļ� ת��Ϊhash�ĺ�����  # 
# û��ftp��ַ����ʹ�ý���        #
#                                #
##################################

##################################
# report�ļ� ת��Ϊhash�ĺ�����  # 1
################################## 
#ר����� eukarotic ��report�ļ�תhash
#�������룬��Ǳ�������� ��ģ���е�ȫ�ֱ��� $JL_GenomeWork_reports_euakryotes
#������ Extract_genome_hash_from_reportFIle
#������� ����ǽ������hash��key�� ������������Ϊһ�����������ܶ�Ӧ��� taxid�����Եڶ�����һ��ARRAY
sub Extract_eukarotic_GenomeInfomHash{   #  my $GenomeInfomrHash= GenomeDATAfeach::Extract_genome_hash_from_reportFIle();
	#���� �������������hash
	
	my $warnMsgBody="\nIn package  GenomeDATAfeach,\tIn sub Extract_eukarotic_GenomeInfomHash,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $subCallereIfm=DirFileHandle::print_SubCallerInform;
	
	my $GenomeInfomrHash=GenomeDATAfeach::Extract_genome_hash_from_reportFIle($JL_GenomeWork_reports_euakryotes);
	
	return $GenomeInfomrHash;
	
}



##################################
# report�ļ� ת��Ϊhash�ĺ�����  # 2
################################## 
#ר����� eukarotic ��report�ļ�תhash
#�������룬��Ǳ�������� ��ģ���е�ȫ�ֱ��� $JL_GenomeWork_reports_euakryotes
#������ Extract_genome_hash_from_reportFIle
#û���������Ǳ�ڵ������ ��ģ���е�ȫ�ֱ��� JL_GenomeWork_reports_euakryotes_HASH
sub Extract_eukarotic_GenomeInfomHash_and_dumperIT{   #  GenomeDATAfeach::Extract_eukarotic_GenomeInfomHash_and_dumperIT();
	#���ض�λ�ã����� �������������hash
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
# report�ļ� ת��Ϊhash�ĺ�����  # 3
################################## 
#      ���� ���� ����            #  
##################################
#����Ǵ�report �ļ����л�ȡ�����Ϣ����ת��Ϊhash�ģ���report�ļ� û�л�������ncbi�ϵĴ洢λ����Ϣ���������ں��������ز���
#���������� report�ļ�������ǽ������hash��key�� ������������Ϊһ�����������ܶ�Ӧ��� taxid�����Եڶ�����һ��ARRAY
sub Extract_genome_hash_from_reportFIle{   #  my $GenomeInfomrHash=GenomeDATAfeach::Extract_genome_hash_from_reportFIle($inReportFile);
	#  ��ĳ��genome��report�ļ��� ��ȡȫ��ϢHASH
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
