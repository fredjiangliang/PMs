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

                      
package TaxonomyWork;

my $JL_taxonomy_Database_DIR   = "/home/fredjiang/EightT/fredjiang.2018.04.02/Taxnomy_database.20181220";
my $JL_taxonomy_hashFILE       =$JL_taxonomy_Database_DIR."/0_0_taxonomy_hash.Hsh";
my $JL_taxonomy_Xml_holderDIR  =$JL_taxonomy_Database_DIR."/0_1_taxonomy_xmlf_DIR";
my $taxID_to_StNm_hashFILE    =$JL_taxonomy_Database_DIR."/1_0_taxID_to_StNm.Hsh";
my $taxID_to_Lnge_hashFILE    =$JL_taxonomy_Database_DIR."/1_1_taxID_to_Lnge.Hsh";
my $StNm_to_taxID_hashFILE    =$JL_taxonomy_Database_DIR."/1_2_StNm_to_taxID.Hsh";
my $StNm_to_Linge_hashFILE    =$JL_taxonomy_Database_DIR."/1_3_StNm_to_Linge.Hsh";

my $JL_SubDIrSize=1000;

my $InFl_SpeciseTreeDataInformationFile="/home/fredjiang/work/Algae/20150522NewDATA/2018.05.03.InPutDataForAlgaeAnalysis/specise.tree.inforrmation.txt";    

#if (   ref ($IdxHash_SpeciseTreeDataInformation)                  eq 'HASH' ) { DirFileHandle::PrintDumper ("IdxHash_SpeciseTreeDataInformation.hash", $IdxHash_SpeciseTreeDataInformation); }

sub Get_TaxID_from_StfNM{     #TaxonomyWork::Get_TaxID_from_StfNM($inSftName);
	my ($inSftName)=@_;
	
	my $warnMsgBody="\nIn package  TaxonomyWork,\tIn sub Get_TaxID_from_StfNM,\n\n";	
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


sub Get_Lineage_from_StfNM{     #TaxonomyWork::Get_Lineage_from_StfNM($inSftName);
	my ($inSftName)=@_;
	
	my $warnMsgBody="\nIn package  TaxonomyWork,\tIn sub Get_Lineage_from_StfNM,\n\n";	
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


sub Change_LinageString_toARRAY{   #TaxonomyWork::Change_LinageString_toARRAY($inString, $inSpecesNm);
	
	my ($inString, $inSpecesNm, $addtionIfomr)=@_;
	
	my $warnMsgBody="\nIn package  TaxonomyWork,\tIn sub Change_LinageString_toARRAY,\n\n";	
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

sub Change_LineageARRAY_to_String{   #TaxonomyWork::Change_LineageARRAY_to_String($inArray);
	
	my ($inArray)=@_;
	
	my $warnMsgBody="\nIn package  TaxonomyWork,\tIn sub Change_LineageARRAY_to_String,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody"; 
  #my $subCallereIfm=DirFileHandle::print_SubCallerInform; #warn $subCallereIfm;
	
	
	if (   (  defined ( $inArray )  ) && (  ref ( $inArray )  eq 'ARRAY'  )   ){
		foreach my $eachStr ( @{ $inArray } ){
		  $eachStr=~s/^\s+//;
		  $eachStr=~s/\s+$//;
		  
	  }
	}
	
	my $outString=join ( ',', @{ $inArray });
	return $outString;	
	
	
}

sub Update_JL_taxonomy_HASH_from_Algae_data{  #  TaxonomyWork::Update_JL_taxonomy_HASH_from_Algae_data( $inArray_or_In_SpeciesNAME )
	my ( $inArray_or_In_SpeciesNAME )=@_;
	my $warnMsgBody="\nIn package  TaxonomyWork,\tIn sub Update_JL_taxonomy_HASH_from_Algae_data,\n\n";	
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


sub build_JL_taxonomy_HASH_from_Algae_data{   #TaxonomyWork::build_JL_taxonomy_HASH_from_Algae_data();
	
	my $warnMsgBody="\nIn package  TaxonomyWork,\tIn sub build_JL_taxonomy_HASH_from_Algae_data,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody"; #warn "$warnMsgHead\$inFastaFile=$inFastaFile\n";
  #my $subCallereIfm=DirFileHandle::print_SubCallerInform; #warn $subCallereIfm;
	
	
	my $JL_taxonomy_HASH;    	 if (  -e ( $JL_taxonomy_hashFILE    )  ) {  $JL_taxonomy_HASH     =Storable::retrieve ($JL_taxonomy_hashFILE); 	}   
	
  my $JL_taxID_to_StNm_HASH; if (  -e ( $taxID_to_StNm_hashFILE )  ) {  $JL_taxID_to_StNm_HASH=Storable::retrieve ( $taxID_to_StNm_hashFILE );  }
  my $JL_taxID_to_Lnge_HASH; if (  -e ( $taxID_to_Lnge_hashFILE )  ) {  $JL_taxID_to_Lnge_HASH=Storable::retrieve ( $taxID_to_Lnge_hashFILE );  }
  my $JL_StNm_to_taxID_HASH; if (  -e ( $StNm_to_taxID_hashFILE )  ) {  $JL_StNm_to_taxID_HASH=Storable::retrieve ( $StNm_to_taxID_hashFILE );  }
  my $JL_StNm_to_Linge_HASH; if (  -e ( $StNm_to_Linge_hashFILE )  ) {  $JL_StNm_to_Linge_HASH=Storable::retrieve ( $StNm_to_Linge_hashFILE );  }
	
	
	my $entrezTaxnomy_db = Bio::DB::Taxonomy->new(-source => 'entrez');
	my $SubDIrSize=$JL_SubDIrSize;
	
	my $IdxHash_SpeciseTreeDataInformation=InFileHandle::FormTableToHash( $InFl_SpeciseTreeDataInformationFile, 0);   
	
	my $workIDhash;
  if		(   (  defined ( $IdxHash_SpeciseTreeDataInformation )  ) && (  ref ( $IdxHash_SpeciseTreeDataInformation ) eq 'HASH'  )    ){
		foreach my $eachSPecNM (    sort { $a cmp $b } (   keys (  %{ $IdxHash_SpeciseTreeDataInformation }  )   )    ){ 
		  my $eachSPecNM_id = $entrezTaxnomy_db->get_taxonid( $eachSPecNM );  print "\$eachSPecNM=$eachSPecNM \$eachSPecNM_id=$eachSPecNM_id\n";
		  $workIDhash->{$eachSPecNM_id}=$eachSPecNM;
		}
	}
	
	
	
  my $isThere_anything_new=0;
  if		(   (  defined ( $workIDhash )  ) && (  ref ( $workIDhash ) eq 'HASH'  )    ){
		foreach my $eachSPecNM_id (    sort { $a <=> $b } (   keys (  %{ $workIDhash }  )   )    ){ 
		  
		  my $efetchFactory= Bio::DB::EUtilities->new(  -eutil    =>   'efetch',
                                                    -email    =>   'fredjiang240@126.com',
                                                    -db       =>   'taxonomy',
                                                    -id       =>   $eachSPecNM_id, 
                                                    #-rettype =>   'taxid'
                                                                                            );
      
      if (      (  defined ( $JL_taxonomy_HASH )  ) && (  ref ( $JL_taxonomy_HASH ) eq 'HASH'  ) 
	           && (  defined ( $JL_taxonomy_HASH->{  '0_1_0_taxonomy_ID_HASH'} )  ) && (  ref ( $JL_taxonomy_HASH->{'0_1_0_taxonomy_ID_HASH'} ) eq 'HASH'  ) 
	           && (  defined ( $JL_taxonomy_HASH->{  '0_1_0_taxonomy_ID_HASH'}->{$eachSPecNM_id} )  )    
	       )
	    {
	      #do nothing	
	    }
	    else{
	      $isThere_anything_new=1;
	    	    	  
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




sub build_JL_taxonomy_HASH_from_TaxidArray{   #TaxonomyWork::build_JL_taxonomy_HASH_from_TaxidArray();
	my ($inTaxidArray)=@_;
	
	my $warnMsgBody="\nIn package  TaxonomyWork,\tIn sub build_JL_taxonomy_HASH_from_TaxidArray,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody"; #warn "$warnMsgHead\$inFastaFile=$inFastaFile\n";
  #my $subCallereIfm=DirFileHandle::print_SubCallerInform; #warn $subCallereIfm;
	
	
	my $JL_taxonomy_HASH;    	 if (  -e ( $JL_taxonomy_hashFILE    )  ) {  $JL_taxonomy_HASH     =Storable::retrieve ($JL_taxonomy_hashFILE); 	}   
	
  my $JL_taxID_to_StNm_HASH; if (  -e ( $taxID_to_StNm_hashFILE )  ) {  $JL_taxID_to_StNm_HASH=Storable::retrieve ( $taxID_to_StNm_hashFILE );  }
  my $JL_taxID_to_Lnge_HASH; if (  -e ( $taxID_to_Lnge_hashFILE )  ) {  $JL_taxID_to_Lnge_HASH=Storable::retrieve ( $taxID_to_Lnge_hashFILE );  }
  my $JL_StNm_to_taxID_HASH; if (  -e ( $StNm_to_taxID_hashFILE )  ) {  $JL_StNm_to_taxID_HASH=Storable::retrieve ( $StNm_to_taxID_hashFILE );  }
  my $JL_StNm_to_Linge_HASH; if (  -e ( $StNm_to_Linge_hashFILE )  ) {  $JL_StNm_to_Linge_HASH=Storable::retrieve ( $StNm_to_Linge_hashFILE );  }
	
	
	my $entrezTaxnomy_db = Bio::DB::Taxonomy->new(-source => 'entrez');
	my $SubDIrSize=$JL_SubDIrSize;
	
	my $IdxHash_SpeciseTreeDataInformation=InFileHandle::FormTableToHash( $InFl_SpeciseTreeDataInformationFile, 0);   
	
	my $workIDhash;
  if		(   (  defined ( $IdxHash_SpeciseTreeDataInformation )  ) && (  ref ( $IdxHash_SpeciseTreeDataInformation ) eq 'HASH'  )    ){
		foreach my $eachSPecNM (    sort { $a cmp $b } (   keys (  %{ $IdxHash_SpeciseTreeDataInformation }  )   )    ){ 
		  my $eachSPecNM_id = $entrezTaxnomy_db->get_taxonid( $eachSPecNM );  print "\$eachSPecNM=$eachSPecNM \$eachSPecNM_id=$eachSPecNM_id\n";
		  $workIDhash->{$eachSPecNM_id}=$eachSPecNM;
		}
	}
	
	
	
  my $isThere_anything_new=0;
  if		(   (  defined ( $workIDhash )  ) && (  ref ( $workIDhash ) eq 'HASH'  )    ){
		foreach my $eachSPecNM_id (    sort { $a <=> $b } (   keys (  %{ $workIDhash }  )   )    ){ 
		  
		  my $efetchFactory= Bio::DB::EUtilities->new(  -eutil    =>   'efetch',
                                                    -email    =>   'fredjiang240@126.com',
                                                    -db       =>   'taxonomy',
                                                    -id       =>   $eachSPecNM_id, 
                                                    #-rettype =>   'taxid'
                                                                                            );
      
      if (      (  defined ( $JL_taxonomy_HASH )  ) && (  ref ( $JL_taxonomy_HASH ) eq 'HASH'  ) 
	           && (  defined ( $JL_taxonomy_HASH->{  '0_1_0_taxonomy_ID_HASH'} )  ) && (  ref ( $JL_taxonomy_HASH->{'0_1_0_taxonomy_ID_HASH'} ) eq 'HASH'  ) 
	           && (  defined ( $JL_taxonomy_HASH->{  '0_1_0_taxonomy_ID_HASH'}->{$eachSPecNM_id} )  )    
	       )
	    {
	      #do nothing	
	    }
	    else{
	      $isThere_anything_new=1;
	    	    	  
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




sub Just_dieWork{
	my ( $dieWord )=@_;
	
	print $dieWord;
	die   $dieWord;
}

1;