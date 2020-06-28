#!/usr/bin/perl -w
BEGIN{  push (@INC, '/home/fredjiang/fredProgs/PMs')  };

use strict;
use warnings;
use HTML::TreeBuilder;
use XML::Simple;
use Data::Dumper;
$Data::Dumper::Purity=1;
use Storable ;
use Storable qw(nstore store_fd nstore_fd freeze thaw dclone);
use Excel::Writer::XLSX;
use Excel::Writer::XLSX::Utility;
use Bio::Seq;
use Bio::SeqIO;
use File::Basename;
use Bio::DB::EUtilities;
use Bio::DB::Taxonomy;
use DirFileHandle;
use TimeWork;
use ExcelHandle;
use ForeachHash;
use InFileHandle;
use OnlyOnePidWork;

                      
package TaxonomyWork_NEW;


##################################################################################################
#
#  使用以下工具可以通过 物种的taxid，获得物种的 taxnomy的详细数据
#
#
##################################################################################################


my $JL_taxonomy_Database_DIR   = "/home/fredjiang/OneTSSD/fredjiang20190321/Taxnomy_database.20190401";


my $JL_taxonomy_hashFILE       =$JL_taxonomy_Database_DIR."/0_0_taxonomy_hash.Hsh";
my $JL_taxonomy_Xml_holderDIR  =$JL_taxonomy_Database_DIR."/0_1_taxonomy_xmlf_DIR";
my $taxID_to_StNm_hashFILE     =$JL_taxonomy_Database_DIR."/1_0_taxID_to_StNm.Hsh";
my $taxID_to_Lnge_hashFILE     =$JL_taxonomy_Database_DIR."/1_1_taxID_to_Lnge.Hsh";
my $StNm_to_taxID_hashFILE     =$JL_taxonomy_Database_DIR."/1_2_StNm_to_taxID.Hsh";
my $StNm_to_Linge_hashFILE     =$JL_taxonomy_Database_DIR."/1_3_StNm_to_Linge.Hsh";
my $working_record_dirctry     =$JL_taxonomy_Database_DIR."/2_0_working_recod_DIR";

my $JL_taxnomy_working_flagFile=$JL_taxonomy_Database_DIR."/0_0_0_Working.flag.file";
my $Longest_waiting_time='1d';


my $JL_SubDIrSize=1000;  #

#my $InFl_SpeciseTreeDataInformationFile="/home/fredjiang/work/Algae/20150522NewDATA/2018.05.03.InPutDataForAlgaeAnalysis/specise.tree.inforrmation.txt";  
#if (   ref ($IdxHash_SpeciseTreeDataInformation)                  eq 'HASH' ) { DirFileHandle::PrintDumper ("IdxHash_SpeciseTreeDataInformation.hash", $IdxHash_SpeciseTreeDataInformation); }

sub Get_TaxID_from_StfNM{     #TaxonomyWork_NEW::Get_TaxID_from_StfNM($inSftName);
	my ($inSftName)=@_;
	
	my $warnMsgBody="\nIn package  TaxonomyWork_NEW,\tIn sub Get_TaxID_from_StfNM,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody"; #warn "$warnMsgHead\$inFastaFile=$inFastaFile\n";
  #my $subCallereIfm=DirFileHandle::print_SubCallerInform; #warn $subCallereIfm;
	
	
	my $JL_StNm_to_taxID_HASH; 
	if (  -e ( $StNm_to_taxID_hashFILE )  ) {  $JL_StNm_to_taxID_HASH=Storable::retrieve ( $StNm_to_taxID_hashFILE );  }
	else        { my $dieMSG=$die_MsgHead."Please check the \$StNm_to_taxID_hashFILE=$StNm_to_taxID_hashFILE, it is not right!!!\n\n\n"; 	    	  &Just_dieWork( $dieMSG );}
	
	my $outTaxID;
	if		(   (  defined ( $JL_StNm_to_taxID_HASH )  ) && (  ref ( $JL_StNm_to_taxID_HASH ) eq 'HASH'  )  && (  defined ( $JL_StNm_to_taxID_HASH->{$inSftName} )  ) && ( $JL_StNm_to_taxID_HASH->{$inSftName}=~m/\S+/ )    ){
		$outTaxID=$JL_StNm_to_taxID_HASH->{$inSftName};		
	}
	else {
		my $dieMSG=$die_MsgHead."Please check the \$StNm_to_taxID_hashFILE=$StNm_to_taxID_hashFILE, it is not right!!!\n\n\$StNm_to_taxID_hashFILE->{$inSftName}=$StNm_to_taxID_hashFILE->{$inSftName} is also not right!!!!\n\n\n"; 	    	  
		&Just_dieWork( $dieMSG );
	}
	
	return $outTaxID;
}


sub Get_Lineage_from_StfNM{     #TaxonomyWork_NEW::Get_Lineage_from_StfNM($inSftName);
	my ($inSftName)=@_;
	
	my $warnMsgBody="\nIn package  TaxonomyWork_NEW,\tIn sub Get_Lineage_from_StfNM,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody"; #warn "$warnMsgHead\$inFastaFile=$inFastaFile\n";
  #my $subCallereIfm=DirFileHandle::print_SubCallerInform; #warn $subCallereIfm;
	
	
	my $JL_StNm_to_Linge_HASH; 
	if (  -e ( $StNm_to_Linge_hashFILE )  ) {  $JL_StNm_to_Linge_HASH=Storable::retrieve ( $StNm_to_Linge_hashFILE );  }
	else        { my $dieMSG=$die_MsgHead."Please check the \$StNm_to_Linge_hashFILE=$StNm_to_Linge_hashFILE, it is not right!!!\n\n\n"; 	    	  &Just_dieWork( $dieMSG );}
	
	my $outLineage;
	if		(   (  defined ( $JL_StNm_to_Linge_HASH )  ) && (  ref ( $JL_StNm_to_Linge_HASH ) eq 'HASH'  )  && (  defined ( $JL_StNm_to_Linge_HASH->{$inSftName} )  ) && ( $JL_StNm_to_Linge_HASH->{$inSftName}=~m/\S+/ )    ){
		$outLineage=$JL_StNm_to_Linge_HASH->{$inSftName};		
	}
	else {
		my $dieMSG=$die_MsgHead."Please check the \$StNm_to_Linge_hashFILE=$StNm_to_Linge_hashFILE, it is not right!!!\n\n\$JL_StNm_to_Linge_HASH->{$inSftName}=$JL_StNm_to_Linge_HASH->{$inSftName} is also not right!!!!\n\n\n"; 	    	  
		&Just_dieWork( $dieMSG );
	}
	
	return $outLineage;
}


sub Change_LinageString_toARRAY{   #TaxonomyWork_NEW::Change_LinageString_toARRAY($inString, $inSpecesNm);
	
	my ($inString, $inSpecesNm, $addtionIfomr)=@_;
	
	my $warnMsgBody="\nIn package  TaxonomyWork_NEW,\tIn sub Change_LinageString_toARRAY,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody"; 
  #my $subCallereIfm=DirFileHandle::print_SubCallerInform; #warn $subCallereIfm;
	
	my $addtionMsg='';  $addtionMsg=$addtionIfomr if (   (  defined ( $addtionIfomr )  ) && ( $addtionIfomr=~m/\S+/ )   ) ;
	
	if (   (  defined ( $inString )  ) && ( $inString=~m/\S+/ )   ){		}
	else{
		my $msg="$addtionMsg \n\$inString=$inString \$inSpecesNm=$inSpecesNm\n the \$inString=$inString is not right\n\n\n";
		warn $msg; print $msg; #sleep(3);
	}
	my @stigArry=split (';',$inString);
	foreach my $eachStr ( @stigArry ){
		$eachStr=~s/^\s+//;
		$eachStr=~s/\s+$//;
	}
	
	if (   (  defined ( $inSpecesNm )  ) && ( $inSpecesNm=~m/\S+/ )   ){
		my $lastIdx=@stigArry;
		if ( $stigArry[$lastIdx-1] eq $inSpecesNm ){			
		}
		else {
			push @stigArry, $inSpecesNm;
		}	
	}
	
	if ( $stigArry[0] eq  'cellular organisms' ){
		shift @stigArry;
	}
	
	
	return [@stigArry];
	
}

sub Change_LinageString_toHASH_notTOaddSpeciseNAME{   #my $outHASH=TaxonomyWork_NEW::Change_LinageString_toHASH_notTOaddSpeciseNAME($inString);
	# 将 物种 lineage信息，变成HASH
	my ($inString)=@_;
	
	my $warnMsgBody="\nIn package  TaxonomyWork_NEW,\tIn sub Change_LinageString_toHASH_notTOaddSpeciseNAME,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody"; 
  my $subCallereIfm=DirFileHandle::print_SubCallerInform; #warn $subCallereIfm;
	
	
	
	if (   (  defined ( $inString )  ) && ( $inString=~m/\S+/ )   ){		}
	else{
		DieWork::Just_dieWork( $die_MsgHead."\$inString=$inString should be a defined string with containt!!\n\n" );
	}
	my @stigArry=split (';',$inString);
	
	my $outHASH;
	my $OrderNB=0;
	foreach my $eachStr ( @stigArry ){
		$eachStr=~s/^\s+//;
		$eachStr=~s/\s+$//;
		$outHASH->{$eachStr}=$OrderNB;
		$OrderNB++;
	}
	
	return $outHASH;
	
}


sub Change_LineageARRAY_to_String{   #TaxonomyWork_NEW::Change_LineageARRAY_to_String($inArray);
	
	my ($inArray)=@_;
	
	my $warnMsgBody="\nIn package  TaxonomyWork_NEW,\tIn sub Change_LineageARRAY_to_String,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody"; 
  my $subCallereIfm=DirFileHandle::print_SubCallerInform; #warn $subCallereIfm;
	
	
	if (   (  defined ( $inArray )  ) && (  ref ( $inArray )  eq 'ARRAY'  )   ){
		foreach my $eachStr ( @{ $inArray } ){
		  $eachStr=~s/^\s+//;
		  $eachStr=~s/\s+$//;
		  
	  }
	}
	
	my $outString=join ( ',', @{ $inArray });
	return $outString;	
	
	
}


sub Build_LineageHASH_for_quick_search{  #my $bigSearchHASH=TaxonomyWork_NEW::Build_LineageHASH_for_quick_search($inTaxidArray);
	my ($inTaxidArray)=@_;
	
	#针对 一个写有物种id的数组($inTaxidArray) ，建立其所有id对应的 lineage的hash，这个hash用来，搜索这个id数组的子集中，包含特定lineage名的成员，如搜索所有的eukaryotes。
	
	my $warnMsgBody="\nIn package  TaxonomyWork_NEW,\tIn sub Build_LineageHASH_for_quick_search,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody"; #warn "$warnMsgHead\$inFastaFile=$inFastaFile\n";
  my $subCallereIfm=DirFileHandle::print_SubCallerInform; #warn $subCallereIfm;
  
  
  my $taxID_to_Lnge_hash=Storable::retrieve( $taxID_to_Lnge_hashFILE )  if (  -e ( $taxID_to_Lnge_hashFILE )  ) ;
  
  my $bigSearchHASH;
  foreach my $eachSPecNM_id (  @ { $inTaxidArray }  ){   print "\$eachSPecNM_id=$eachSPecNM_id\n";
    if (   (  defined ( $taxID_to_Lnge_hash )  ) && (  ref ( $taxID_to_Lnge_hash ) eq 'HASH'  ) && (  defined ( $taxID_to_Lnge_hash->{$eachSPecNM_id} )  ) && ( $taxID_to_Lnge_hash->{$eachSPecNM_id}=~m/\S+/)   ){
		  my $LngeString=$taxID_to_Lnge_hash->{$eachSPecNM_id};	
		  my $outHASH=TaxonomyWork_NEW::Change_LinageString_toHASH_notTOaddSpeciseNAME($LngeString);
		  $bigSearchHASH->{$eachSPecNM_id}=$outHASH;
		}
		else{
			DieWork::Just_dieWork( $die_MsgHead."\$taxID_to_Lnge_hash=$taxID_to_Lnge_hash \$taxID_to_Lnge_hash->{$eachSPecNM_id}=$taxID_to_Lnge_hash->{$eachSPecNM_id} not right!!\n\n" );
			
		}
  }
	return $bigSearchHASH;
	
}

sub Update_JL_taxonomy_HASH_from_Algae_data{  #  TaxonomyWork_NEW::Update_JL_taxonomy_HASH_from_Algae_data( $inArray_or_In_SpeciesNAME )
	my ( $inArray_or_In_SpeciesNAME )=@_;
	my $warnMsgBody="\nIn package  TaxonomyWork_NEW,\tIn sub Update_JL_taxonomy_HASH_from_Algae_data,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody"; #warn "$warnMsgHead\$inFastaFile=$inFastaFile\n";
  #my $subCallereIfm=DirFileHandle::print_SubCallerInform; #warn $subCallereIfm;
	
	my $work_stfcNM_hash;
	if		(   (  defined ( $inArray_or_In_SpeciesNAME )  ) && (  ref ( $inArray_or_In_SpeciesNAME ) eq 'ARRAY'  )    ){
		foreach my $eachKey (  @{ $inArray_or_In_SpeciesNAME }  ){
			$work_stfcNM_hash->{$eachKey}=1;
		}
	}
	elsif (   (  defined ( $inArray_or_In_SpeciesNAME )  ) && ( $inArray_or_In_SpeciesNAME=~m/\S+/ )    ){
		$work_stfcNM_hash->{$inArray_or_In_SpeciesNAME}=1;
	} 
	else{
		DieWork::Just_dieWork( $die_MsgHead."\$inArray_or_In_SpeciesNAME=$inArray_or_In_SpeciesNAME is not right!!\n\n" );
	}  
	
	my $JL_taxonomy_HASH;    	 if (  -e ( $JL_taxonomy_hashFILE    )  ) {  $JL_taxonomy_HASH     =Storable::retrieve ($JL_taxonomy_hashFILE); 	}   
	
  my $JL_taxID_to_StNm_HASH; if (  -e ( $taxID_to_StNm_hashFILE )  ) {  $JL_taxID_to_StNm_HASH=Storable::retrieve ( $taxID_to_StNm_hashFILE );  }
  my $JL_taxID_to_Lnge_HASH; if (  -e ( $taxID_to_Lnge_hashFILE )  ) {  $JL_taxID_to_Lnge_HASH=Storable::retrieve ( $taxID_to_Lnge_hashFILE );  }
  my $JL_StNm_to_taxID_HASH; if (  -e ( $StNm_to_taxID_hashFILE )  ) {  $JL_StNm_to_taxID_HASH=Storable::retrieve ( $StNm_to_taxID_hashFILE );  }
  my $JL_StNm_to_Linge_HASH; if (  -e ( $StNm_to_Linge_hashFILE )  ) {  $JL_StNm_to_Linge_HASH=Storable::retrieve ( $StNm_to_Linge_hashFILE );  }
	
	
	my $entrezTaxnomy_db = Bio::DB::Taxonomy->new(-source => 'entrez');
	my $SubDIrSize=$JL_SubDIrSize;
	
	my $workIDhash;
	if		(   (  defined ( $work_stfcNM_hash )  ) && (  ref ( $work_stfcNM_hash ) eq 'HASH'  )    ){
		foreach my $eachSPecNM (    sort { $a cmp $b } (   keys (  %{ $work_stfcNM_hash }  )   )    ){  warn "\$eachSPecNM=$eachSPecNM\n";
		  my $eachSPecNM_id = $entrezTaxnomy_db->get_taxonid( $eachSPecNM );  warn "\$eachSPecNM=$eachSPecNM \$eachSPecNM_id=$eachSPecNM_id\n";
		  $workIDhash->{$eachSPecNM_id}=$eachSPecNM;
		}
	}
	
  
		
  my $isThere_anything_new=0;
  if		(   (  defined ( $workIDhash )  ) && (  ref ( $workIDhash ) eq 'HASH'  )    ){
		foreach my $eachSPecNM_id (    sort { $a <=> $b } (   keys (  %{ $workIDhash }  )   )    ){ 
		  
		  
      
      if (      (  defined ( $JL_taxonomy_HASH )  ) && (  ref ( $JL_taxonomy_HASH ) eq 'HASH'  ) 
	           && (  defined ( $JL_taxonomy_HASH->{  '0_1_0_taxonomy_ID_HASH'} )  ) && (  ref ( $JL_taxonomy_HASH->{'0_1_0_taxonomy_ID_HASH'} ) eq 'HASH'  ) 
	           && (  defined ( $JL_taxonomy_HASH->{  '0_1_0_taxonomy_ID_HASH'}->{$eachSPecNM_id} )  )    
	       )
	    {
	      #do nothing	
	    }
	    else{
	      $isThere_anything_new=1;
	    	
	    	my $efetchFactory= Bio::DB::EUtilities->new(  -eutil    =>   'efetch',
                                                      -email    =>   'fredjiang240@126.com',
                                                      -db       =>   'taxonomy',
                                                      -id       =>   $eachSPecNM_id, 
                                                      #-rettype =>   'taxid'
                                                                                            );
	    	
	    	    	  
	    	my $New_SubDir;
	    	my $New_DrNumb;
	    	my $New_InsubN; 	
	    	my $New_totNmb;
        if (      (  defined ( $JL_taxonomy_HASH )  ) && (  ref ( $JL_taxonomy_HASH ) eq 'HASH'  ) 
	    	       && (  defined ( $JL_taxonomy_HASH->{'0_0_0_LastSubDir'} )  ) && ( $JL_taxonomy_HASH->{'0_0_0_LastSubDir'}=~m/\d+/ )
	    	       && (  defined ( $JL_taxonomy_HASH->{'0_0_1_LsSbDrNumb'} )  ) && ( $JL_taxonomy_HASH->{'0_0_1_LsSbDrNumb'}=~m/\d+/ )
	    	       && (  defined ( $JL_taxonomy_HASH->{'0_0_2_LastInsubN'} )  ) && ( $JL_taxonomy_HASH->{'0_0_2_LastInsubN'}=~m/\d+/ )
	    	       && (  defined ( $JL_taxonomy_HASH->{'0_0_3_total_Numb'} )  ) && ( $JL_taxonomy_HASH->{'0_0_3_total_Numb'}=~m/\d+/ )
	    	   )
	    	{
	    	  	    	    		  
	    	  $New_SubDir=$JL_taxonomy_HASH->{'0_0_0_LastSubDir'};
	    	  $New_DrNumb=$JL_taxonomy_HASH->{'0_0_1_LsSbDrNumb'};
	    	  $New_InsubN=$JL_taxonomy_HASH->{'0_0_2_LastInsubN'}+1;
	    	  $New_totNmb=$JL_taxonomy_HASH->{'0_0_3_total_Numb'}+1;
	    	  
	    	  if ( $New_InsubN > $SubDIrSize ){
	    	  	$New_InsubN=1;
	    	  	$New_DrNumb=$New_DrNumb+1;
	    	  	$New_SubDir=sprintf ("%05d",$New_DrNumb);
	    	  }
	    	  
	    	  $JL_taxonomy_HASH->{'0_0_0_LastSubDir'}=$New_SubDir;
	    	  $JL_taxonomy_HASH->{'0_0_1_LsSbDrNumb'}=$New_DrNumb;
	    	  $JL_taxonomy_HASH->{'0_0_2_LastInsubN'}=$New_InsubN; 
	    	  $JL_taxonomy_HASH->{'0_0_3_total_Numb'}=$New_totNmb;  
	    	}
	    	else{
	    		$New_SubDir=1;     
	    	  $New_DrNumb=1;    $New_SubDir=sprintf ("%05d",$New_DrNumb);
	    	  $New_InsubN=1;
	    		
	    		$JL_taxonomy_HASH->{'0_0_0_LastSubDir'}=$New_SubDir;
	    	  $JL_taxonomy_HASH->{'0_0_1_LsSbDrNumb'}=$New_DrNumb;
	    	  $JL_taxonomy_HASH->{'0_0_2_LastInsubN'}=$New_InsubN;
	    	  $JL_taxonomy_HASH->{'0_0_3_total_Numb'}=1;
	    	}  
	    	
	    	my $sftcNMAE=  $workIDhash->{$eachSPecNM_id} ;
                            $JL_taxonomy_HASH->{ '0_2_0_StfNm_to_ID_HASH' }->{ $sftcNMAE }         =$eachSPecNM_id;   
                            
                            $JL_taxID_to_StNm_HASH->{$eachSPecNM_id}=$sftcNMAE;                            
                            $JL_StNm_to_taxID_HASH->{$sftcNMAE     }=$eachSPecNM_id;
                            
                            
                        
        my $txanomy_eacDIR                                                                                             =$JL_taxonomy_Xml_holderDIR."/".$New_SubDir; system ( "mkdir -p $txanomy_eacDIR " ) ;                                                                                    
        my $txanomyXmlFile= $JL_taxonomy_HASH->{ '0_1_0_taxonomy_ID_HASH' }->{$eachSPecNM_id}->{'0_1_0_txanomyXmlFile'}=$txanomy_eacDIR."/".$eachSPecNM_id.".xml.txt";                                                                                     
        my $XmlPrasedHshFl= $JL_taxonomy_HASH->{ '0_1_0_taxonomy_ID_HASH' }->{$eachSPecNM_id}->{'0_1_1_XmlPrasedHshFl'}=$txanomyXmlFile.".xmlhsh";                                                                                     
        #my $XmlPrasedHshFl= $JL_taxonomy_HASH->{  '0_1_0_taxonomy_ID_HASH'}->{$eachSPecNM_id}->{'0_1_1_XmlPrasedHshFl'}=$txanomyXmlFile.".xmlhsh";                                                                                     
                                                                                                   	    
 	      $efetchFactory->get_Response(-file => $txanomyXmlFile);
 	
 	      
 	
 	      my $Xml_Prased_HASH = XML::Simple::XMLin($txanomyXmlFile);
 	      DirFileHandle::PrintDumper($XmlPrasedHshFl, $Xml_Prased_HASH) if  (   (  defined ( $Xml_Prased_HASH )  ) && (  ref ( $Xml_Prased_HASH ) eq 'HASH'  )   );
 	      
 	      if (      (  defined ( $Xml_Prased_HASH )  ) && (  ref ( $Xml_Prased_HASH ) eq 'HASH'  ) && (  defined ( $Xml_Prased_HASH->{'Taxon'} )  ) && (  ref ( $Xml_Prased_HASH->{'Taxon'} ) eq 'HASH'  )     ){
 	      	                  $JL_taxonomy_HASH->{  '0_1_0_taxonomy_ID_HASH'}->{$eachSPecNM_id}->{'0_2_0_ScientificName'}=$Xml_Prased_HASH->{'Taxon'}->{'ScientificName'} if (      (  defined ( $Xml_Prased_HASH->{'Taxon'}->{'ScientificName'} )  ) && (  $Xml_Prased_HASH->{'Taxon'}->{'ScientificName'}=~m/\S+/  )    );
 	      	  my $LineageNme =$JL_taxonomy_HASH->{  '0_1_0_taxonomy_ID_HASH'}->{$eachSPecNM_id}->{'0_2_1________Lineage'}=$Xml_Prased_HASH->{'Taxon'}->{'Lineage'       } if (      (  defined ( $Xml_Prased_HASH->{'Taxon'}->{'Lineage'       } )  ) && (  $Xml_Prased_HASH->{'Taxon'}->{'Lineage'       }=~m/\S+/  )    );
 	      	                  
 	      	                  $JL_taxID_to_Lnge_HASH->{$eachSPecNM_id} =$LineageNme  if  (   (  defined ( $LineageNme )  ) && ( $LineageNme=~/\S+/ )   ); 
 	      	                  $JL_StNm_to_Linge_HASH->{$sftcNMAE     } =$LineageNme  if  (   (  defined ( $LineageNme )  ) && ( $LineageNme=~/\S+/ )   );        
 	      	                  
 	      	                  
 	      	                  $JL_taxonomy_HASH->{  '0_1_0_taxonomy_ID_HASH'}->{$eachSPecNM_id}->{'0_3_0_all_Other_HASH'}=Storable::dclone ( $Xml_Prased_HASH->{'Taxon'}  );
 	      }
 	      else {
 	    	  my $dieMSG=$die_MsgHead."Please check the \$txanomyXmlFile=$txanomyXmlFile, it is not right!!!\n\n\n";
 	    	  &Just_dieWork( $dieMSG );
 	      }
 	     
		  }
		  
		}
	}
	
	if ( $isThere_anything_new==1){
		DirFileHandle::PrintDumper($JL_taxonomy_hashFILE, $JL_taxonomy_HASH) if  (   (  defined ( $JL_taxonomy_HASH )  ) && (  ref ( $JL_taxonomy_HASH ) eq 'HASH'  )   );  
		
		DirFileHandle::PrintDumper($taxID_to_StNm_hashFILE, $JL_taxID_to_StNm_HASH) if  (   (  defined ( $JL_taxID_to_StNm_HASH )  ) && (  ref ( $JL_taxID_to_StNm_HASH ) eq 'HASH'  )   );  
		DirFileHandle::PrintDumper($taxID_to_Lnge_hashFILE, $JL_taxID_to_Lnge_HASH) if  (   (  defined ( $JL_taxID_to_Lnge_HASH )  ) && (  ref ( $JL_taxID_to_Lnge_HASH ) eq 'HASH'  )   );  
		DirFileHandle::PrintDumper($StNm_to_taxID_hashFILE, $JL_StNm_to_taxID_HASH) if  (   (  defined ( $JL_StNm_to_taxID_HASH )  ) && (  ref ( $JL_StNm_to_taxID_HASH ) eq 'HASH'  )   );  
		DirFileHandle::PrintDumper($StNm_to_Linge_hashFILE, $JL_StNm_to_Linge_HASH) if  (   (  defined ( $JL_StNm_to_Linge_HASH )  ) && (  ref ( $JL_StNm_to_Linge_HASH ) eq 'HASH'  )   );  
		

		
	}
	
}





sub build_JL_taxonomy_HASH_from_TaxidArray{   #TaxonomyWork_NEW::build_JL_taxonomy_HASH_from_TaxidArray( $inTaxidArray );
	# 在工作程序 前后加上判定函数，这个工作程序 在服务器上，一次只能有一个在运行
	my ($inTaxidArray)=@_;
	
	my $warnMsgBody="\nIn package  TaxonomyWork_NEW,\tIn sub build_JL_taxonomy_HASH_from_TaxidArray,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody"; #warn "$warnMsgHead\$inFastaFile=$inFastaFile\n";
  my $subCallereIfm=DirFileHandle::print_SubCallerInform; #warn $subCallereIfm;
  
  
  my $nowWorkTimePidInformMark=OnlyOnePidWork::CheckAndStart_step_forOnlyOnePidWork($JL_taxnomy_working_flagFile, $Longest_waiting_time);
  if (   ( defined ( $nowWorkTimePidInformMark )  ) && ( $nowWorkTimePidInformMark=~/\S+/ )   ){
  	
  	TaxonomyWork_NEW::build_JL_taxonomy_HASH_from_TaxidArray_inOnlyOnePidwork( $inTaxidArray );
  	
  	OnlyOnePidWork::CheckAnd_Done_step_forOnlyOnePidWork( $JL_taxnomy_working_flagFile, $nowWorkTimePidInformMark );
  	
  }
  
}


sub build_JL_taxonomy_HASH_from_TaxidArray_inOnlyOnePidwork{   #TaxonomyWork_NEW::build_JL_taxonomy_HASH_from_TaxidArray_inOnlyOnePidwork( $inTaxidArray );
	my ($inTaxidArray)=@_;
	
	my $warnMsgBody="\nIn package  TaxonomyWork_NEW,\tIn sub build_JL_taxonomy_HASH_from_TaxidArray,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  my $subCallereIfm=DirFileHandle::print_SubCallerInform; #warn $subCallereIfm;
	
	my $group_to_eFatch_size=1000;
	my $reTryMaxTime=10;
	my $sleepTime=5;
	
	
	my $JL_taxonomy_HASH;    	 if (  -e ( $JL_taxonomy_hashFILE    )  ) {  $JL_taxonomy_HASH    =Storable::retrieve ($JL_taxonomy_hashFILE); 	}   
	
  my $JL_taxID_to_StNm_HASH; if (  -e ( $taxID_to_StNm_hashFILE )  ) {  $JL_taxID_to_StNm_HASH=Storable::retrieve ( $taxID_to_StNm_hashFILE );  }
  my $JL_taxID_to_Lnge_HASH; if (  -e ( $taxID_to_Lnge_hashFILE )  ) {  $JL_taxID_to_Lnge_HASH=Storable::retrieve ( $taxID_to_Lnge_hashFILE );  }
  my $JL_StNm_to_taxID_HASH; if (  -e ( $StNm_to_taxID_hashFILE )  ) {  $JL_StNm_to_taxID_HASH=Storable::retrieve ( $StNm_to_taxID_hashFILE );  }
  my $JL_StNm_to_Linge_HASH; if (  -e ( $StNm_to_Linge_hashFILE )  ) {  $JL_StNm_to_Linge_HASH=Storable::retrieve ( $StNm_to_Linge_hashFILE );  }
	
	
	
	my $SubDIrSize=$JL_SubDIrSize;
	
	my $DeRepeatIDhash;
  if		(   (  defined ( $inTaxidArray )  ) && (  ref ( $inTaxidArray ) eq 'ARRAY'  )    ){
		foreach my $eachSPecNM_id (  @ { $inTaxidArray }  ){   print "\$eachSPecNM_id=$eachSPecNM_id\n";		  
		  $DeRepeatIDhash->{$eachSPecNM_id}=1;
		}
	}
	
	
	my $grouped_taxID_array;  my $group_IDX_lv_0=0; my $group_IDX_lv_1=0;
	my $isThere_anything_new   =0;   #检查Hash，看有没有新的taxid的信息要更新  
  if		(   (  defined ( $DeRepeatIDhash )  ) && (  ref ( $DeRepeatIDhash ) eq 'HASH'  )    ){
		foreach my $eachSPecNM_id (    sort { $a <=> $b } (   keys (  %{ $DeRepeatIDhash }  )   )    ){ 
		  
		  
      
      if (      (  defined ( $JL_taxonomy_HASH )  ) && (  ref ( $JL_taxonomy_HASH ) eq 'HASH'  ) 
	           && (  defined ( $JL_taxonomy_HASH->{  '0_1_0_taxonomy_ID_HASH'} )  ) && (  ref ( $JL_taxonomy_HASH->{'0_1_0_taxonomy_ID_HASH'} ) eq 'HASH'  ) 
	           && (  defined ( $JL_taxonomy_HASH->{  '0_1_0_taxonomy_ID_HASH'}->{$eachSPecNM_id} )  )  && (  ref ( $JL_taxonomy_HASH->{'0_1_0_taxonomy_ID_HASH'}->{$eachSPecNM_id} ) eq 'HASH'  ) 
	           && (  defined ( $JL_taxonomy_HASH->{  '0_1_0_taxonomy_ID_HASH'}->{$eachSPecNM_id}->{'0_2_0_ScientificName'} )  )  && (  $JL_taxonomy_HASH->{  '0_1_0_taxonomy_ID_HASH'}->{$eachSPecNM_id}->{'0_2_0_ScientificName'}=~m/\S+/  )
	           
	           && (  defined ( $JL_taxID_to_StNm_HASH )  ) && (  ref ( $JL_taxID_to_StNm_HASH ) eq 'HASH'  )
	           && (  defined ( $JL_taxID_to_StNm_HASH->{$eachSPecNM_id} )  ) && (  $JL_taxID_to_StNm_HASH->{$eachSPecNM_id}=~m/\S+/  )
	           
	           && (  defined ( $JL_taxID_to_Lnge_HASH )  ) && (  ref ( $JL_taxID_to_Lnge_HASH ) eq 'HASH'  )
	           && (  defined ( $JL_taxID_to_Lnge_HASH->{$eachSPecNM_id} )  ) && (  $JL_taxID_to_Lnge_HASH->{$eachSPecNM_id}=~m/\S+/  )
	           
	           && (  defined ( $JL_StNm_to_taxID_HASH )  ) && (  ref ( $JL_StNm_to_taxID_HASH ) eq 'HASH'  )
	           && (  defined ( $JL_StNm_to_taxID_HASH->{$JL_taxonomy_HASH->{  '0_1_0_taxonomy_ID_HASH'}->{$eachSPecNM_id}->{'0_2_0_ScientificName'}} )  ) && (  $JL_StNm_to_taxID_HASH->{$JL_taxonomy_HASH->{  '0_1_0_taxonomy_ID_HASH'}->{$eachSPecNM_id}->{'0_2_0_ScientificName'}}=~m/\S+/  )
	           
	           && (  defined ( $JL_StNm_to_Linge_HASH )  ) && (  ref ( $JL_StNm_to_Linge_HASH ) eq 'HASH'  )
	           && (  defined ( $JL_StNm_to_Linge_HASH->{$JL_taxonomy_HASH->{  '0_1_0_taxonomy_ID_HASH'}->{$eachSPecNM_id}->{'0_2_0_ScientificName'}} )  ) && (  $JL_StNm_to_Linge_HASH->{$JL_taxonomy_HASH->{  '0_1_0_taxonomy_ID_HASH'}->{$eachSPecNM_id}->{'0_2_0_ScientificName'}}=~m/\S+/  )
	           
	              
	       )
	    {
	      #do nothing	
	    }
	    else{
	      $isThere_anything_new=1;  
	    	
	    	$grouped_taxID_array->[$group_IDX_lv_0]->[$group_IDX_lv_1]=$eachSPecNM_id;
	    	if    (  $group_IDX_lv_1 >= ( $group_to_eFatch_size -1 )  ){
	    		$group_IDX_lv_0++;
	    		$group_IDX_lv_1=0;
	    	}
	    	else{
	    		$group_IDX_lv_1++;
	    	}
	    	
	    	
	    }
	  }
	}
	
	
	#############
	if (   (  defined ( $grouped_taxID_array )  ) && (  ref( $grouped_taxID_array ) eq 'ARRAY' )   ){
		system ( "mkdir -p $working_record_dirctry" );
	  my $NtPid=TimeWork::GetNowTimePid_microSecond ();
	  my $record_file="$working_record_dirctry/$NtPid\.txt";
    my $warnMsg= DirFileHandle::ReturnDumperInform ($grouped_taxID_array);	 
	  InFileHandle::PrintStringIntoFile($record_file, $warnMsg);  #print "$ptOutWrn"    ; 
	}
	
	
	#############
	
	my $New_SubDir;
	my $New_DrNumb;
	my $New_InsubN; 	
	my $New_totNmb;
	my $New_TxIDnb;
	if		(   (  defined ( $grouped_taxID_array )  ) && (  ref ( $grouped_taxID_array ) eq 'ARRAY'  )    ){
	  for ( my $i=0; $i<@{ $grouped_taxID_array }; $i++ ){
	  	if		(   (  defined ( $grouped_taxID_array->[$i] )  ) && (  ref ( $grouped_taxID_array->[$i] ) eq 'ARRAY'  )    ){
	  	  my $taxIDnumber=@{ $grouped_taxID_array->[$i] };
	  	  
	  	  #### 更新大Hash的数据统计项   # START ###########################################################################################################################################
	  	  
	  	  
	  	  if (      (  defined ( $JL_taxonomy_HASH )  ) && (  ref ( $JL_taxonomy_HASH ) eq 'HASH'  ) 
	    	       && (  defined ( $JL_taxonomy_HASH->{'0_0_0_LastSubDir'} )  ) && ( $JL_taxonomy_HASH->{'0_0_0_LastSubDir'}=~m/\d+/ )
	    	       && (  defined ( $JL_taxonomy_HASH->{'0_0_1_LsSbDrNumb'} )  ) && ( $JL_taxonomy_HASH->{'0_0_1_LsSbDrNumb'}=~m/\d+/ )
	    	       && (  defined ( $JL_taxonomy_HASH->{'0_0_2_LastInsubN'} )  ) && ( $JL_taxonomy_HASH->{'0_0_2_LastInsubN'}=~m/\d+/ )
	    	       && (  defined ( $JL_taxonomy_HASH->{'0_0_3_total_Numb'} )  ) && ( $JL_taxonomy_HASH->{'0_0_3_total_Numb'}=~m/\d+/ )
	    	       && (  defined ( $JL_taxonomy_HASH->{'0_0_4_taxID_Numb'} )  ) && ( $JL_taxonomy_HASH->{'0_0_4_taxID_Numb'}=~m/\d+/ )
	    	   )
	    	{
	    	  	    	    		  
	    	  $New_SubDir=$JL_taxonomy_HASH->{'0_0_0_LastSubDir'};
	    	  $New_DrNumb=$JL_taxonomy_HASH->{'0_0_1_LsSbDrNumb'};
	    	  $New_InsubN=$JL_taxonomy_HASH->{'0_0_2_LastInsubN'}+1;
	    	  $New_totNmb=$JL_taxonomy_HASH->{'0_0_3_total_Numb'}+1;
	    	  $New_TxIDnb=$JL_taxonomy_HASH->{'0_0_4_taxID_Numb'}+$taxIDnumber;
	    	  
	    	  if ( $New_InsubN > $SubDIrSize ){
	    	  	$New_InsubN=1;
	    	  	$New_DrNumb=$New_DrNumb+1;
	    	  	$New_SubDir=sprintf ("%05d",$New_DrNumb);
	    	  }
	    	  
	    	  $JL_taxonomy_HASH->{'0_0_0_LastSubDir'}=$New_SubDir;
	    	  $JL_taxonomy_HASH->{'0_0_1_LsSbDrNumb'}=$New_DrNumb;
	    	  $JL_taxonomy_HASH->{'0_0_2_LastInsubN'}=$New_InsubN; 
	    	  $JL_taxonomy_HASH->{'0_0_3_total_Numb'}=$New_totNmb;  
	    	  $JL_taxonomy_HASH->{'0_0_4_taxID_Numb'}=$New_TxIDnb;  
	    	}
	    	else{
	    		$New_SubDir=1;     
	    	  $New_DrNumb=1;    $New_SubDir=sprintf ("%05d",$New_DrNumb);
	    	  $New_InsubN=1;
	    	  $New_totNmb=1;
	    	  $New_TxIDnb=$taxIDnumber;
	    		
	    		$JL_taxonomy_HASH->{'0_0_0_LastSubDir'}=$New_SubDir;
	    	  $JL_taxonomy_HASH->{'0_0_1_LsSbDrNumb'}=$New_DrNumb;
	    	  $JL_taxonomy_HASH->{'0_0_2_LastInsubN'}=$New_InsubN;
	    	  $JL_taxonomy_HASH->{'0_0_3_total_Numb'}=$New_totNmb;
	    	  $JL_taxonomy_HASH->{'0_0_4_taxID_Numb'}=$New_TxIDnb;  
	    	}  
	    	#### 更新大Hash的数据统计项 #  END  ###########################################################################################################################################
	  	  
	  	  
	  	  my $txanomy_eacDIR= $JL_taxonomy_Xml_holderDIR."/".$New_SubDir; system ( "mkdir -p $txanomy_eacDIR " ) ;                                                                                    
        my $txanomyXmlFile= $txanomy_eacDIR."/".$New_totNmb.".xml.txt";                                                                                     
        my $XmlPrasedHshFl= $txanomyXmlFile.".xmlhsh";                                                                                     
        
        
        #### 网上获取数据 # START ############
        my $efetchFactory= Bio::DB::EUtilities->new(  -eutil    =>   'efetch',
                                                      -email    =>   'fredjiang240@126.com',
                                                      -db       =>   'taxonomy',
                                                      -id       =>    $grouped_taxID_array->[$i], 
                                                       #-rettype =>   'taxid'
                                                                                                  );                             
          
        if ( -e ( $txanomyXmlFile ) ){ system ("rm -f $txanomyXmlFile"); }
        
        
        my $Work_done_flag=0;
        FORMK: for ( my $i=0; $i<$reTryMaxTime; $i++ ) {                                  
          my $pidKey=fork();    
	        if (  !defined (  $pidKey )  ) {                  DieWork::Just_dieWork( $die_MsgHead."\n Error in fork: $!".$subCallereIfm );       }
              
          if ($pidKey == 0) {    warn "$warnMsgBody\n  Child   fork: My pid = $$\t\t \$txanomyXmlFile=$txanomyXmlFile\n \$i=$i\n\n";  print "$warnMsgBody\n Child   fork: My pid = $$\t\t \$txanomyXmlFile=$txanomyXmlFile\n\n";
                                              
            $efetchFactory->get_Response(-file => $txanomyXmlFile);
            exit 0;
          } 
          
          waitpid($pidKey, 0);
          
          if ( -e ( $txanomyXmlFile ) ){
          	$Work_done_flag=1;
          	last FORMK;    
          }
          else{
          	my $ciShu=$i+1;
          	my $sleep_Seconds=$sleepTime*$ciShu;
          	my $warnMsg="\n\n\$txanomyXmlFile=$txanomyXmlFile\n Try to efatch the $ciShu time!!\nSleep for $sleep_Seconds seconds!!!\n\n\n";
          	warn  $warnMsg;
          	print $warnMsg;
          	sleep( $sleep_Seconds );
          }
        } 
        if ( $Work_done_flag==0 ){
        	DieWork::Just_dieWork( $die_MsgHead."\n\$txanomyXmlFile=$txanomyXmlFile maybe cannot  connect to eutils.ncbi.nlm.nih.gov !!!: $!".$subCallereIfm );
        }   	
        #### 网上获取数据 #  END  ############
                                                                                                 	    
 	     
 	      #### 解析 xml 文件， 对输出hash进行 填写 ，并打出来 # START #####################################
 	      my $org_Xml_Prased_HASH = XML::Simple::XMLin($txanomyXmlFile);
 	      my $Xml_Prased_HASH;
 	      if (      (  defined ( $org_Xml_Prased_HASH )  ) && (  ref ( $org_Xml_Prased_HASH ) eq 'HASH'  ) && (  defined ( $org_Xml_Prased_HASH->{'Taxon'} )  ) && (  ref ( $org_Xml_Prased_HASH->{'Taxon'} ) eq 'HASH'  )     ){
 	      	$Xml_Prased_HASH->{'Taxon'}->[0]=Storable::dclone( $org_Xml_Prased_HASH->{'Taxon'} );
 	      }
 	      else{
 	      	$Xml_Prased_HASH=Storable::dclone ( $org_Xml_Prased_HASH );
 	      }
 	      
 	      DirFileHandle::PrintDumper($XmlPrasedHshFl, $Xml_Prased_HASH) if  (   (  defined ( $Xml_Prased_HASH )  ) && (  ref ( $Xml_Prased_HASH ) eq 'HASH'  )   );
        
        if (      (  defined ( $Xml_Prased_HASH )  ) && (  ref ( $Xml_Prased_HASH ) eq 'HASH'  ) && (  defined ( $Xml_Prased_HASH->{'Taxon'} )  ) && (  ref ( $Xml_Prased_HASH->{'Taxon'} ) eq 'ARRAY'  )     ){
 	      	for (  my $j=0; $j< @{ $Xml_Prased_HASH->{'Taxon'} }; $j++  ){
 	      	  my $taxid=$Xml_Prased_HASH->{'Taxon'}->[$j]->{'TaxId'};
 	      	  
 	      	                
 	      	  
 	      	  
 	      	  if (   (  defined ( $taxid )  ) && ( $taxid=~/\S+/ ) && (  defined ( $grouped_taxID_array->[$i]->[$j] )  ) && ( $grouped_taxID_array->[$i]->[$j]=~/\S+/ )   ){
 	      	  	if ( $taxid=~m/^\s*-?\d+(\.\d+)?\s*$/ ){
 	      	  		if ( $grouped_taxID_array->[$i]->[$j] == $taxid ){ }
 	      	  		else{ 
 	      	  			if (   (  defined ( $Xml_Prased_HASH->{'Taxon'}->[$j]->{'AkaTaxIds'}->{'TaxId'} )  ) && ( $grouped_taxID_array->[$i]->[$j] == $Xml_Prased_HASH->{'Taxon'}->[$j]->{'AkaTaxIds'}->{'TaxId'} )  ){
 	      	  				$taxid=$Xml_Prased_HASH->{'Taxon'}->[$j]->{'AkaTaxIds'}->{'TaxId'};
 	      	  			}
 	      	  			else{
 	      	  				DieWork::Just_dieWork( $die_MsgHead."\n\$grouped_taxID_array->[$i]->[$j]=$grouped_taxID_array->[$i]->[$j] should == $Xml_Prased_HASH->{'Taxon'}->[$j]->{'AkaTaxIds'}->{'TaxId'}=\$Xml_Prased_HASH->{'Taxon'}->[$j]->{'AkaTaxIds'}->{'TaxId'}\n    : $!".$subCallereIfm ); 
 	      	  			}
 	      	  			
 	      	  		}
 	      	  	}
 	      	  	else{
 	      	  		if ( $grouped_taxID_array->[$i]->[$j] eq $taxid ){ }
 	      	  		else{ 
 	      	  			if (   (  defined ( $Xml_Prased_HASH->{'Taxon'}->[$j]->{'AkaTaxIds'}->{'TaxId'} )  ) && ( $grouped_taxID_array->[$i]->[$j] eq $Xml_Prased_HASH->{'Taxon'}->[$j]->{'AkaTaxIds'}->{'TaxId'} )  ){
 	      	  				$taxid=$Xml_Prased_HASH->{'Taxon'}->[$j]->{'AkaTaxIds'}->{'TaxId'};
 	      	  			}
 	      	  			else{
 	      	  				DieWork::Just_dieWork( $die_MsgHead."\n\$grouped_taxID_array->[$i]->[$j]=$grouped_taxID_array->[$i]->[$j] should eq $Xml_Prased_HASH->{'Taxon'}->[$j]->{'AkaTaxIds'}->{'TaxId'}=\$Xml_Prased_HASH->{'Taxon'}->[$j]->{'AkaTaxIds'}->{'TaxId'}\n    : $!".$subCallereIfm ); 
 	      	  			}
 	      	  		}
 	      	  	}
 	      	  }
 	      	  else{ DieWork::Just_dieWork( $die_MsgHead."\n\$grouped_taxID_array->[$i]->[$j]=$grouped_taxID_array->[$i]->[$j]  and \$taxid=$taxid should be all defined and not as NULL\n    : $!".$subCallereIfm );  }
            
                          $JL_taxonomy_HASH->{ '0_1_0_taxonomy_ID_HASH' }->{$taxid}->{'0_1_0_txanomyXmlFile'}=$txanomyXmlFile;                                                                                     
                          $JL_taxonomy_HASH->{ '0_1_0_taxonomy_ID_HASH' }->{$taxid}->{'0_1_1_XmlPrasedHshFl'}=$XmlPrasedHshFl; 
                          
 	      	                $JL_taxonomy_HASH->{ '0_1_0_taxonomy_ID_HASH' }->{$taxid}->{'0_1_2_ThisID_HASHidx'}=$j;
 	      	                $JL_taxonomy_HASH->{ '0_1_0_taxonomy_ID_HASH' }->{$taxid}->{'0_1_3_wayToGet__HASH'}="HASH=\$Xml_Prased_HASH->{'Taxon'}->[\$j]=\$Xml_Prased_HASH->{'Taxon'}->[$j]";
 	      	  
 	      	  my $sftcNMAE= $JL_taxonomy_HASH->{ '0_1_0_taxonomy_ID_HASH' }->{$taxid}->{'0_2_0_ScientificName'}=$Xml_Prased_HASH->{'Taxon'}->[$j]->{'ScientificName'} if (      (  defined ( $Xml_Prased_HASH->{'Taxon'}->[$j]->{'ScientificName'} )  ) && (  $Xml_Prased_HASH->{'Taxon'}->[$j]->{'ScientificName'}=~m/\S+/  )    );
 	      	  
 	      	  if (  defined ( $sftcNMAE )  ) {           
                          $JL_taxonomy_HASH->{ '0_2_0_StfNm_to_ID_HASH' }->{ $sftcNMAE }         =$taxid;   
                          
                          $JL_taxID_to_StNm_HASH->{$taxid        }=$sftcNMAE;                            
                          $JL_StNm_to_taxID_HASH->{$sftcNMAE     }=$taxid;
 	      	  }  
 	      	  
 	      	  my $LineageNme =$JL_taxonomy_HASH->{  '0_1_0_taxonomy_ID_HASH'}->{$taxid}->{'0_2_1________Lineage'}=$Xml_Prased_HASH->{'Taxon'}->[$j]->{'Lineage'       } if (      (  defined ( $Xml_Prased_HASH->{'Taxon'}->[$j]->{'Lineage'       } )  ) && (  $Xml_Prased_HASH->{'Taxon'}->[$j]->{'Lineage'       }=~m/\S+/  )    );
 	      	                
 	      	                $JL_taxID_to_Lnge_HASH->{$taxid        } =$LineageNme  if  (   (  defined ( $LineageNme )  ) && ( $LineageNme=~/\S+/ )   ); 
 	      	                $JL_StNm_to_Linge_HASH->{$sftcNMAE     } =$LineageNme  if  (   (  defined ( $LineageNme )  ) && ( $LineageNme=~/\S+/ )   );        
 	      	                
 	      	                
 	      	                $JL_taxonomy_HASH->{  '0_1_0_taxonomy_ID_HASH'}->{$taxid}->{'0_3_0_all_Other_HASH'}=Storable::dclone ( $Xml_Prased_HASH->{'Taxon'}->[$j]  );
 	      	  
 	      	}
 	      	  
 	      }
 	      else {
 	    	  DieWork::Just_dieWork( $die_MsgHead."\n \$Xml_Prased_HASH=Xml_Prased_HASH  should be a HASH ref!!\n\$Xml_Prased_HASH->{'Taxon'}=$Xml_Prased_HASH->{'Taxon'} should not a ARRAY ref!!  $!\n\n\n".$subCallereIfm ); 
 	      }
        
	  	  DirFileHandle::PrintDumper($JL_taxonomy_hashFILE  , $JL_taxonomy_HASH     ) if  (   (  defined ( $JL_taxonomy_HASH )  ) && (  ref ( $JL_taxonomy_HASH ) eq 'HASH'  )   );  
		
		    DirFileHandle::PrintDumper($taxID_to_StNm_hashFILE, $JL_taxID_to_StNm_HASH) if  (   (  defined ( $JL_taxID_to_StNm_HASH )  ) && (  ref ( $JL_taxID_to_StNm_HASH ) eq 'HASH'  )   );  
		    DirFileHandle::PrintDumper($taxID_to_Lnge_hashFILE, $JL_taxID_to_Lnge_HASH) if  (   (  defined ( $JL_taxID_to_Lnge_HASH )  ) && (  ref ( $JL_taxID_to_Lnge_HASH ) eq 'HASH'  )   );  
		    DirFileHandle::PrintDumper($StNm_to_taxID_hashFILE, $JL_StNm_to_taxID_HASH) if  (   (  defined ( $JL_StNm_to_taxID_HASH )  ) && (  ref ( $JL_StNm_to_taxID_HASH ) eq 'HASH'  )   );  
		    DirFileHandle::PrintDumper($StNm_to_Linge_hashFILE, $JL_StNm_to_Linge_HASH) if  (   (  defined ( $JL_StNm_to_Linge_HASH )  ) && (  ref ( $JL_StNm_to_Linge_HASH ) eq 'HASH'  )   ); 
	  	  
	  	  #### 解析 xml 文件， 对输出hash进行 填写 ，并打出来 #  END  #####################################
	  	  
	  	}
	  }
	}
}	

	    	


sub build_JL_taxonomy_HASH_from_TaxidArray_inOnlyOnePidwork_oldVersion{   #TaxonomyWork_NEW::build_JL_taxonomy_HASH_from_TaxidArray_inOnlyOnePidwork_oldVersion( $inTaxidArray );
	my ($inTaxidArray)=@_;
	
	my $warnMsgBody="\nIn package  TaxonomyWork_NEW,\tIn sub build_JL_taxonomy_HASH_from_TaxidArray_inOnlyOnePidwork_oldVersion,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody"; #warn "$warnMsgHead\$inFastaFile=$inFastaFile\n";
  my $subCallereIfm=DirFileHandle::print_SubCallerInform; #warn $subCallereIfm;
	
	my $reTryMaxTime=10;
	my $sleepTime=5;
	
	
	my $JL_taxonomy_HASH;    	 if (  -e ( $JL_taxonomy_hashFILE    )  ) {  $JL_taxonomy_HASH    =Storable::retrieve ($JL_taxonomy_hashFILE); 	}   
	
  my $JL_taxID_to_StNm_HASH; if (  -e ( $taxID_to_StNm_hashFILE )  ) {  $JL_taxID_to_StNm_HASH=Storable::retrieve ( $taxID_to_StNm_hashFILE );  }
  my $JL_taxID_to_Lnge_HASH; if (  -e ( $taxID_to_Lnge_hashFILE )  ) {  $JL_taxID_to_Lnge_HASH=Storable::retrieve ( $taxID_to_Lnge_hashFILE );  }
  my $JL_StNm_to_taxID_HASH; if (  -e ( $StNm_to_taxID_hashFILE )  ) {  $JL_StNm_to_taxID_HASH=Storable::retrieve ( $StNm_to_taxID_hashFILE );  }
  my $JL_StNm_to_Linge_HASH; if (  -e ( $StNm_to_Linge_hashFILE )  ) {  $JL_StNm_to_Linge_HASH=Storable::retrieve ( $StNm_to_Linge_hashFILE );  }
	
	
	#my $entrezTaxnomy_db = Bio::DB::Taxonomy->new(-source => 'entrez');
	my $SubDIrSize=$JL_SubDIrSize;
	
	#my $IdxHash_SpeciseTreeDataInformation=InFileHandle::FormTableToHash( $InFl_SpeciseTreeDataInformationFile, 0);   
	
	my $workIDhash;
  if		(   (  defined ( $inTaxidArray )  ) && (  ref ( $inTaxidArray ) eq 'ARRAY'  )    ){
		foreach my $eachSPecNM_id (  @ { $inTaxidArray }  ){   print "\$eachSPecNM_id=$eachSPecNM_id\n";		  
		  $workIDhash->{$eachSPecNM_id}=1;
		}
	}
	
	
	
  my $isThere_anything_new   =0;   #检查Hash，看有没有新的taxid的信息要更新
  my $new_sum_nb_form_0      =0;   #检查Hash，看有没有新的taxid的信息要更新,并计数
  my $new_sum_nb_startToWrite=100; #检查Hash，看有没有新的taxid的信息要更新,并计数, 当达到 $new_sum_nb_startToWrite=100 个时，进行更新，写入各个Hash。
  if		(   (  defined ( $workIDhash )  ) && (  ref ( $workIDhash ) eq 'HASH'  )    ){
		foreach my $eachSPecNM_id (    sort { $a <=> $b } (   keys (  %{ $workIDhash }  )   )    ){ 
		  
		  
      
      if (      (  defined ( $JL_taxonomy_HASH )  ) && (  ref ( $JL_taxonomy_HASH ) eq 'HASH'  ) 
	           && (  defined ( $JL_taxonomy_HASH->{  '0_1_0_taxonomy_ID_HASH'} )  ) && (  ref ( $JL_taxonomy_HASH->{'0_1_0_taxonomy_ID_HASH'} ) eq 'HASH'  ) 
	           && (  defined ( $JL_taxonomy_HASH->{  '0_1_0_taxonomy_ID_HASH'}->{$eachSPecNM_id} )  )    
	       )
	    {
	      #do nothing	
	    }
	    else{
	      $isThere_anything_new=1;  $new_sum_nb_form_0++;  warn "\n     \$eachSPecNM_id=$eachSPecNM_id \$new_sum_nb_form_0=$new_sum_nb_form_0 \n\n";
	    	    	  
	    	my $New_SubDir;
	    	my $New_DrNumb;
	    	my $New_InsubN; 	
	    	my $New_totNmb;
        if (      (  defined ( $JL_taxonomy_HASH )  ) && (  ref ( $JL_taxonomy_HASH ) eq 'HASH'  ) 
	    	       && (  defined ( $JL_taxonomy_HASH->{'0_0_0_LastSubDir'} )  ) && ( $JL_taxonomy_HASH->{'0_0_0_LastSubDir'}=~m/\d+/ )
	    	       && (  defined ( $JL_taxonomy_HASH->{'0_0_1_LsSbDrNumb'} )  ) && ( $JL_taxonomy_HASH->{'0_0_1_LsSbDrNumb'}=~m/\d+/ )
	    	       && (  defined ( $JL_taxonomy_HASH->{'0_0_2_LastInsubN'} )  ) && ( $JL_taxonomy_HASH->{'0_0_2_LastInsubN'}=~m/\d+/ )
	    	       && (  defined ( $JL_taxonomy_HASH->{'0_0_3_total_Numb'} )  ) && ( $JL_taxonomy_HASH->{'0_0_3_total_Numb'}=~m/\d+/ )
	    	   )
	    	{
	    	  	    	    		  
	    	  $New_SubDir=$JL_taxonomy_HASH->{'0_0_0_LastSubDir'};
	    	  $New_DrNumb=$JL_taxonomy_HASH->{'0_0_1_LsSbDrNumb'};
	    	  $New_InsubN=$JL_taxonomy_HASH->{'0_0_2_LastInsubN'}+1;
	    	  $New_totNmb=$JL_taxonomy_HASH->{'0_0_3_total_Numb'}+1;
	    	  
	    	  if ( $New_InsubN > $SubDIrSize ){
	    	  	$New_InsubN=1;
	    	  	$New_DrNumb=$New_DrNumb+1;
	    	  	$New_SubDir=sprintf ("%05d",$New_DrNumb);
	    	  }
	    	  
	    	  $JL_taxonomy_HASH->{'0_0_0_LastSubDir'}=$New_SubDir;
	    	  $JL_taxonomy_HASH->{'0_0_1_LsSbDrNumb'}=$New_DrNumb;
	    	  $JL_taxonomy_HASH->{'0_0_2_LastInsubN'}=$New_InsubN; 
	    	  $JL_taxonomy_HASH->{'0_0_3_total_Numb'}=$New_totNmb;  
	    	}
	    	else{
	    		$New_SubDir=1;     
	    	  $New_DrNumb=1;    $New_SubDir=sprintf ("%05d",$New_DrNumb);
	    	  $New_InsubN=1;
	    		
	    		$JL_taxonomy_HASH->{'0_0_0_LastSubDir'}=$New_SubDir;
	    	  $JL_taxonomy_HASH->{'0_0_1_LsSbDrNumb'}=$New_DrNumb;
	    	  $JL_taxonomy_HASH->{'0_0_2_LastInsubN'}=$New_InsubN;
	    	  $JL_taxonomy_HASH->{'0_0_3_total_Numb'}=1;
	    	}  
	    	
	    	
                            
                            
                        
        my $txanomy_eacDIR                                                                                             =$JL_taxonomy_Xml_holderDIR."/".$New_SubDir; system ( "mkdir -p $txanomy_eacDIR " ) ;                                                                                    
        my $txanomyXmlFile= $JL_taxonomy_HASH->{ '0_1_0_taxonomy_ID_HASH' }->{$eachSPecNM_id}->{'0_1_0_txanomyXmlFile'}=$txanomy_eacDIR."/".$eachSPecNM_id.".xml.txt";                                                                                     
        my $XmlPrasedHshFl= $JL_taxonomy_HASH->{ '0_1_0_taxonomy_ID_HASH' }->{$eachSPecNM_id}->{'0_1_1_XmlPrasedHshFl'}=$txanomyXmlFile.".xmlhsh";                                                                                     
        #my $XmlPrasedHshFl= $JL_taxonomy_HASH->{  '0_1_0_taxonomy_ID_HASH'}->{$eachSPecNM_id}->{'0_1_1_XmlPrasedHshFl'}=$txanomyXmlFile.".xmlhsh";                                                                                     
                             
        
        my $Xml_Prased_HASH;    my $theTaxIdHashBuild_before=0;                 
        if (   (  defined ( $XmlPrasedHshFl )  ) && (  -e ( $XmlPrasedHshFl )  )   ){
        	$Xml_Prased_HASH=Storable::retrieve( $XmlPrasedHshFl );  #'TaxId'
        	if (      (  defined ( $Xml_Prased_HASH )  ) && (  ref ( $Xml_Prased_HASH ) eq 'HASH'  ) && (  defined ( $Xml_Prased_HASH->{'Taxon'} )  ) && (  ref ( $Xml_Prased_HASH->{'Taxon'} ) eq 'HASH'  )     ){
        		
        		if (   (  defined ( $Xml_Prased_HASH->{'Taxon'}->{'TaxId'} )  ) && ( $Xml_Prased_HASH->{'Taxon'}->{'TaxId'} == $eachSPecNM_id)  ){
        			if (  ( $eachSPecNM_id=~m/^-?\d+\.?(\d+)?$/ )  && ( $Xml_Prased_HASH->{'Taxon'}->{'TaxId'} == $eachSPecNM_id )  ){
        				if   ( $Xml_Prased_HASH->{'Taxon'}->{'TaxId'} == $eachSPecNM_id ){
        					$theTaxIdHashBuild_before=1;
        				}
        			}
        			else{
        				if  ( $Xml_Prased_HASH->{'Taxon'}->{'TaxId'} eq $eachSPecNM_id ){
        					$theTaxIdHashBuild_before=1;
        				}
        			}       		        			
        		}
        	}
        }
        
        if ( $theTaxIdHashBuild_before == 1){
        	
        }
        else{
        	my $efetchFactory= Bio::DB::EUtilities->new(  -eutil    =>   'efetch',
                                                      -email    =>   'fredjiang240@126.com',
                                                      -db       =>   'taxonomy',
                                                      -id       =>   $eachSPecNM_id, 
                                                       #-rettype =>   'taxid'
                                                                                            );                             
          
          if ( -e ( $txanomyXmlFile ) ){ system ("rm -f $txanomyXmlFile"); }
          
          
          my $Work_done_flag=0;
          FORMK: for ( my $i=0; $i<$reTryMaxTime; $i++ ) {                                  
            my $pidKey=fork();    
	          if (  !defined (  $pidKey )  ) {                  DieWork::Just_dieWork( $die_MsgHead."\n Error in fork: $!".$subCallereIfm );       }
                
            if ($pidKey == 0) {    warn "$warnMsgBody\n  Child   fork: My pid = $$\t\t \$txanomyXmlFile=$txanomyXmlFile\n \$i=$i\n\n";  print "$warnMsgBody\n Child   fork: My pid = $$\t\t \$txanomyXmlFile=$txanomyXmlFile\n\n";
                                                
              $efetchFactory->get_Response(-file => $txanomyXmlFile);
              exit 0;
            } 
            
            waitpid($pidKey, 0);
            
            if ( -e ( $txanomyXmlFile ) ){
            	$Work_done_flag=1;
            	last FORMK;    
            }
            else{
            	my $ciShu=$i+1;
            	my $sleep_Seconds=$sleepTime*$ciShu;
            	my $warnMsg="\n\n\$txanomyXmlFile=$txanomyXmlFile\n Try to efatch the $ciShu time!!\nSleep for $sleep_Seconds seconds!!!\n\n\n";
            	warn  $warnMsg;
            	print $warnMsg;
            	sleep( $sleep_Seconds );
            }
          } 
          if ( $Work_done_flag==0 ){
          	DieWork::Just_dieWork( $die_MsgHead."\n\$txanomyXmlFile=$txanomyXmlFile maybe cannot  connect to eutils.ncbi.nlm.nih.gov !!!: $!".$subCallereIfm );
          }   	
                                                                                                   	    
 	        #$efetchFactory->get_Response(-file => $txanomyXmlFile);
 	        
 	        
 	        
 	        $Xml_Prased_HASH = XML::Simple::XMLin($txanomyXmlFile);
 	        DirFileHandle::PrintDumper($XmlPrasedHshFl, $Xml_Prased_HASH) if  (   (  defined ( $Xml_Prased_HASH )  ) && (  ref ( $Xml_Prased_HASH ) eq 'HASH'  )   );
        }     
 	
 	      
 	
 	      
 	      
 	      
 	      if (      (  defined ( $Xml_Prased_HASH )  ) && (  ref ( $Xml_Prased_HASH ) eq 'HASH'  ) && (  defined ( $Xml_Prased_HASH->{'Taxon'} )  ) && (  ref ( $Xml_Prased_HASH->{'Taxon'} ) eq 'HASH'  )     ){
 	      	  my $sftcNMAE=   $JL_taxonomy_HASH->{  '0_1_0_taxonomy_ID_HASH'}->{$eachSPecNM_id}->{'0_2_0_ScientificName'}=$Xml_Prased_HASH->{'Taxon'}->{'ScientificName'} if (      (  defined ( $Xml_Prased_HASH->{'Taxon'}->{'ScientificName'} )  ) && (  $Xml_Prased_HASH->{'Taxon'}->{'ScientificName'}=~m/\S+/  )    );
 	      	  
 	      	  if (  defined ( $sftcNMAE )  ) {           
                            $JL_taxonomy_HASH->{ '0_2_0_StfNm_to_ID_HASH' }->{ $sftcNMAE }         =$eachSPecNM_id;   
                            
                            $JL_taxID_to_StNm_HASH->{$eachSPecNM_id}=$sftcNMAE;                            
                            $JL_StNm_to_taxID_HASH->{$sftcNMAE     }=$eachSPecNM_id;
 	      	  }  
 	      	  
 	      	  my $LineageNme =$JL_taxonomy_HASH->{  '0_1_0_taxonomy_ID_HASH'}->{$eachSPecNM_id}->{'0_2_1________Lineage'}=$Xml_Prased_HASH->{'Taxon'}->{'Lineage'       } if (      (  defined ( $Xml_Prased_HASH->{'Taxon'}->{'Lineage'       } )  ) && (  $Xml_Prased_HASH->{'Taxon'}->{'Lineage'       }=~m/\S+/  )    );
 	      	                  
 	      	                  $JL_taxID_to_Lnge_HASH->{$eachSPecNM_id} =$LineageNme  if  (   (  defined ( $LineageNme )  ) && ( $LineageNme=~/\S+/ )   ); 
 	      	                  $JL_StNm_to_Linge_HASH->{$sftcNMAE     } =$LineageNme  if  (   (  defined ( $LineageNme )  ) && ( $LineageNme=~/\S+/ )   );        
 	      	                  
 	      	                  
 	      	                  $JL_taxonomy_HASH->{  '0_1_0_taxonomy_ID_HASH'}->{$eachSPecNM_id}->{'0_3_0_all_Other_HASH'}=Storable::dclone ( $Xml_Prased_HASH->{'Taxon'}  );
 	      }
 	      else {
 	    	  my $dieMSG=$die_MsgHead."Please check the \$txanomyXmlFile=$txanomyXmlFile, it is not right!!!\n\n\n";
 	    	  &Just_dieWork( $dieMSG );
 	      }
 	     
		  }
		  
		  
		  if ( $new_sum_nb_form_0 >= $new_sum_nb_startToWrite){  #每有 $new_sum_nb_startToWrite （这里设定为10个或100个都可以）个新的结果，则更新一遍所有的Hash
		  	
		  	warn "\n\$new_sum_nb_form_0=$new_sum_nb_form_0 \$new_sum_nb_startToWrite=$new_sum_nb_startToWrite Now do the HASH upgrade writing work\n\n\n";
		  	$new_sum_nb_form_0=0;
		  	
		  	DirFileHandle::PrintDumper($JL_taxonomy_hashFILE, $JL_taxonomy_HASH) if  (   (  defined ( $JL_taxonomy_HASH )  ) && (  ref ( $JL_taxonomy_HASH ) eq 'HASH'  )   );  
		
		    DirFileHandle::PrintDumper($taxID_to_StNm_hashFILE, $JL_taxID_to_StNm_HASH) if  (   (  defined ( $JL_taxID_to_StNm_HASH )  ) && (  ref ( $JL_taxID_to_StNm_HASH ) eq 'HASH'  )   );  
		    DirFileHandle::PrintDumper($taxID_to_Lnge_hashFILE, $JL_taxID_to_Lnge_HASH) if  (   (  defined ( $JL_taxID_to_Lnge_HASH )  ) && (  ref ( $JL_taxID_to_Lnge_HASH ) eq 'HASH'  )   );  
		    DirFileHandle::PrintDumper($StNm_to_taxID_hashFILE, $JL_StNm_to_taxID_HASH) if  (   (  defined ( $JL_StNm_to_taxID_HASH )  ) && (  ref ( $JL_StNm_to_taxID_HASH ) eq 'HASH'  )   );  
		    DirFileHandle::PrintDumper($StNm_to_Linge_hashFILE, $JL_StNm_to_Linge_HASH) if  (   (  defined ( $JL_StNm_to_Linge_HASH )  ) && (  ref ( $JL_StNm_to_Linge_HASH ) eq 'HASH'  )   ); 
		  	
		  	
		  }
		  
		}
	}
	
	if ( $isThere_anything_new==1){   #最后完成后再讲所有Hash写一遍
		 
		DirFileHandle::PrintDumper($JL_taxonomy_hashFILE, $JL_taxonomy_HASH) if  (   (  defined ( $JL_taxonomy_HASH )  ) && (  ref ( $JL_taxonomy_HASH ) eq 'HASH'  )   );  
		
		DirFileHandle::PrintDumper($taxID_to_StNm_hashFILE, $JL_taxID_to_StNm_HASH) if  (   (  defined ( $JL_taxID_to_StNm_HASH )  ) && (  ref ( $JL_taxID_to_StNm_HASH ) eq 'HASH'  )   );  
		DirFileHandle::PrintDumper($taxID_to_Lnge_hashFILE, $JL_taxID_to_Lnge_HASH) if  (   (  defined ( $JL_taxID_to_Lnge_HASH )  ) && (  ref ( $JL_taxID_to_Lnge_HASH ) eq 'HASH'  )   );  
		DirFileHandle::PrintDumper($StNm_to_taxID_hashFILE, $JL_StNm_to_taxID_HASH) if  (   (  defined ( $JL_StNm_to_taxID_HASH )  ) && (  ref ( $JL_StNm_to_taxID_HASH ) eq 'HASH'  )   );  
		DirFileHandle::PrintDumper($StNm_to_Linge_hashFILE, $JL_StNm_to_Linge_HASH) if  (   (  defined ( $JL_StNm_to_Linge_HASH )  ) && (  ref ( $JL_StNm_to_Linge_HASH ) eq 'HASH'  )   );
		

		
	}
	
}




sub Just_dieWork{
	my ( $dieWord )=@_;
	
	print $dieWord;
	die   $dieWord;
}

1;