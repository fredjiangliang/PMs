
#!/usr/bin/perl -w

BEGIN{  push (@INC, '/home/fredjiang/fredProgs/PMs')  };

use strict;
use warnings;


use Bio::Graphics;
use Bio::SearchIO;
use Bio::SearchIO::Writer::HTMLResultWriter;
use Storable ;
use Storable qw(nstore store_fd nstore_fd freeze thaw dclone);

use SeqSegmentsTools;
use ArrayHashChange;
use FastaFileHandle;
use MatrixCsvChange;
use BLOSUMalignSCORE;
use String_Work;

package  BlastPraser;




#新的blast praser，不排序 不求最佳 仅仅解析得到数组
sub BioPerlBlastPraser20191021_justPrase{       #    my $prasedOUThash=BlastPraser::BioPerlBlastPraser20191021_justPrase ($blastfile, $report_type, $idxFnm[ $totalQureyNumber, $numHits] )   #解析Blast文件
  my ($blastfile, $report_type, $idxFnm, $totalQureyNumber, $numHits )=@_;   
  
  
  my $warnMsgBody="\nIn package  BlastPraser,\tIn sub BioPerlBlastPraser20191021_justPrase,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $subCalInfom=DirFileHandle::print_SubCallerInform;
  
  
  if (   (  defined ( $idxFnm )  ) && ( $idxFnm=~m/\S+/ ) && (  -e ( $idxFnm )  )    ){ 	  }
  else{		DieWork::Just_dieWork( $die_MsgHead."\n \$idxFnm=$idxFnm should be a defined index file  !!  $!\n\n\n".$subCalInfom ); 	}
  
  my $Show_totalQureyNumber='';
  if (   (  defined ( $totalQureyNumber )  ) && ( $totalQureyNumber=~m/\S+/ ) && ( $totalQureyNumber > 0  )    ){ 
  	$Show_totalQureyNumber="/$totalQureyNumber";	  
  }
  
  
  #warn  "\n\nNow in Sub &BioPerlpraseBlast!\n";      warn  "Input 1:\$blastfile=$blastfile\t" if (defined ($blastfile)); warn  "Input 2:\$outBestBlastFile=$outBestBlastFile\t" if (defined ($outBestBlastFile)); warn  "Input 3:\$numHits=$numHits\t" if (defined ($numHits));  warn  "Input 4:\$report_type=$report_type\t" if (defined ($report_type));  warn  "\n\n";
  print "\n\cl\nNow in Sub BlastPraser::BioPerlBlastPraser20191021_justPrase!\n:";  print "Input 1:\$blastfile=$blastfile\t" if (defined ($blastfile)); 
  
  print "Input 2:\$numHits=$numHits\t" if (defined ($numHits));  print "Input 3:\$report_type=$report_type\t" if (defined ($report_type));  print "\n\cl\n";
  
  my $in;
  if (   ( defined ($report_type) )  && ($report_type=~m/^\s*tblastn\s*$/i)   ){ 
  	$report_type='tblastn';
    $in = Bio::SearchIO->new(-format => 'blast', -file => $blastfile, -report_type => $report_type); warn "1 \$report_type=$report_type\n";
  }
  elsif (   ( defined ($report_type) )  && ($report_type=~m/^\s*blastx\s*$/i)   ){ 
  	$report_type='blastx';
    $in = Bio::SearchIO->new(-format => 'blast', -file => $blastfile, -report_type => $report_type); warn "2 \$report_type=$report_type\n";
  }
  elsif (   ( defined ($report_type) )  && ($report_type=~m/^\s*blastp\s*$/i)   ){ 
  	$report_type='blastp';
    $in = Bio::SearchIO->new(-format => 'blast', -file => $blastfile, -report_type => $report_type); warn "3 \$report_type=$report_type\n";
  }
  else {
  	$in = Bio::SearchIO->new(-format => 'blast', -file => $blastfile                              ); warn "4 \$report_type=$report_type\n";
  }
  #sleep(3);
  
  my $prasedOUThash;
  my $resultArrayIdx=0;    
  my $maxBestinformToget=5;  #这个数字是用来 限制， $summuryHere这个变量的长度的
  
  
  my $QueryHit_excel_HASH;
  
  my $totalHit=0;
  
  my $Check_finish_HASH;  
  
  my $query_Conut=0;
  
  # extraction of information for each result recursively
  while ( my $result = $in->next_result ) {
	  # the name of the query sequence    	
	  print  $result->query_name . "\t";
	  my $queryMinStart =  99999999999999999999;
    my $queryMaxEnd   = -99999999999999999999;
	  my $queryName=$result->query_name;    DieWork::Print_and_warn( "\n\$queryName=$queryName ($resultArrayIdx$Show_totalQureyNumber)\n\n");
	  
	  my $onlyFindingAA_orNot_inQury=0; 
	  
	  $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_0_0_query_name'}=$queryName; 
   	
   	$Check_finish_HASH->{$queryName}=1;  #这个hash，记录blast结果中所有的query，最终用来判断 输入的fasta文件是否被blast完整计算
   	
   	my $databaseFind=$result->database_name;                                                                                                   #print  "\n20190322-0-1 database=".$databaseFind."\n";
   	#if (   (  defined ( $databaseFind )  ) && ( $databaseFind=~m/\S+/)   ){  #&& (  -e ( $databaseFind )  )    )  
   	#	                                                                                                                                       #print  "\n20190322-0-2 database=".$databaseFind."\n"; 		
   	#} 
   	
   	
   	
   	$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_0_1_query_fastaFile'}=$databaseFind;                                                    print  "\n20190322-0-5 database=".$databaseFind."\n";
  
	  $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_0_2_query_indexFile'}=$idxFnm;      
	 
	  
   
    my $query_length=$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_1_query_length'}=$result->query_length;
    
    $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_2_0_num_hits'}=$result->num_hits;  print  $result->num_hits;
    # output "no hits found" if there is no hits
    $totalHit+=$result->num_hits;  print "\n\$totalHit=$totalHit\t\$result->num_hits=$result->num_hits\n";
    
    my $summuryHere='';  my $summuryHereLength=0;
    if ( $result->num_hits == 0 ) {
		  print  "\tNo hits found\n";
    } 
    else {
		  my $count = 0;   
      # process each hit recursively
		  while (my $hit = $result->next_hit) {
			  print  "\t" if ($count > 0);
                          # get the accession numbers of the hits
			  print  "\t" . $hit->accession . "\t";
                          # get the lengths of the hit sequences
        print  $hit->length . "\t";
                          # get the description of the hit sequences
			  print  $hit->description . "\t";
                          # get the E value of the hit
			  print  $hit->significance . "\t";
                          #get the bit score of the hit
			  print  $hit->bits . "\t";
        
        my $HitName=$hit->accession;
        
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_0_0_locus'}   =$hit->locus();        
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_0_1__name'}   =$hit->name();  
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_0_2_giNub'}   =BlastHandle::GetGiNumber( $hit->name() );  
           
        #$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_0_3_accessionNB'}  =$hit->accession;
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_0_3_accessionNB'}   =$HitName;
        
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_0_4_hitSeqLength'}  =$hit->length;
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_0_5_hitDescription'}=$hit->description;
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_0_6_hitEvalue'}     =$hit->significance;
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_0_7_hitScore'}      =$hit->bits;
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_0_8_algorithm'}      =$hit->algorithm;
        
        
        my $accIDinShort=$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_0_3_accessionNB'};
        $accIDinShort=~s/^(\S+)\s*.*$/$1/;  #214g2702C.paradoxa20150522;thioredoxin_reductase2;Contig53901_-3

        
        my $hspcount = 0; 
        my $totalIdentical=0;
        my $totalConserved=0;
        
        my $totConsePercent=0;  #
        
        
        my $queryCoverHash;  
        my $hitCoverHash;
        
        my $queryCoverRate=0;  
        my $qureyCovRatePct=0;
        my $total_qry_covrAte=0;
        my $total_qry_covrAte_pct=0;
        my $hitCoverRate=0;
        my $hitCovRatePct=0;
        
        my $onlyFindingAA_orNot_inHit=0;  # 用来判断在一个Hit中是否 找到要找的特定的AA的char
        
                          # process the top HSP for the top hit
			  while (my $hsp = $hit->next_hsp) {
			  	
			  	
			  	
          print  "\t\t\t\t\t\t\t", if ($hspcount > 0);
                          	# get the frame of the query sequence
			  	print  $hsp->query->frame . "\t";
                                  # get the start and the end of the query sequence in the alignment
			  	print  $hsp->start('query') . "\t" . $hsp->end('query'). "\t";
                                  # get the start and the end of the hit sequence in the alignment
			  	print  $hsp->start('hit') . "\t" . $hsp->end('hit') . "\t";
                                  # get the similarity value
			  	printf  "%.1f" , ($hsp->frac_conserved * 100);
			  	print  "%\t";
                                  # get the identity value
			  	printf  "%.1f" , ($hsp->frac_identical * 100);
			  	
			  	#$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'_HspQurFrame'}        =$hsp->frame('query');  #oldVersion   #bioperl 得到的blast的hsp的frame是不对的
			  	#$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'_HspHitFrame'}        =$hsp->frame('hit');    #oldVersion   #bioperl 得到的blast的hsp的frame是不对的
			  	
			  	
			  	
			  	#$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'_Hsp_1stfm_3frm_2qur'}       =$hsp->frame('query');        my $QrAbsFm=$hsp->frame('query')+1;  #bioperl 得到的blast的hsp的frame是不对的
			  	#$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'_Hsp_1stfm_4frm_2hit'}       =$hsp->frame('hit');          my $HtAbsFm=$hsp->frame('hit')+1;   #bioperl 得到的blast的hsp的frame是不对的
			  	$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_0_0_Hsp_1stfm_1stn_1qur'}       =$hsp->strand('query');       my $QrZF; if ($hsp->strand('query') >=0){$QrZF='+';}else {$QrZF='-';} 
			  	$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_0_1_Hsp_1stfm_2stn_1hit'}       =$hsp->strand('hit');         my $HtZF; if ($hsp->strand('hit')   >=0){$HtZF='+';}else {$HtZF='-';} 
			  	#$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'_Hsp_1stfm_1qur'}       =$QrZF.$QrAbsFm;   #bioperl 得到的blast的hsp的frame是不对的
			  	#$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'_Hsp_1stfm_2hit'}       =$HtZF.$HtAbsFm;   #bioperl 得到的blast的hsp的frame是不对的
			  	$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_0_2_Hsp_1stfm_3ZoF_1qur'}       =$QrZF;       
			  	$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_0_3_Hsp_1stfm_4ZoF_1hit'}       =$HtZF;         
			  	
			  	$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_1_0_HspQueryStart'}      =$hsp->start('query');  #oldVersion
			  	$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_1_1_HspQueryEnd'}        =$hsp->end('query');    #oldVersion
			  	$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_1_2_HspHitStart'}        =$hsp->start('hit');    #oldVersion
			  	$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_1_3_HspHitEnd'}          =$hsp->end('hit');      #oldVersion
			  	
			  	my $QryHspStt=			  	
			  	$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_2_0_Hsp_2StEd_1Qr_1St'}        =$hsp->start('query');       $queryCoverHash->[$hspcount]->{'_start'}=$hsp->start('query'); 
			  	my $QryHspEnd=
			  	$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_2_1_Hsp_2StEd_1Qr_2Ed'}        =$hsp->end('query');         $queryCoverHash->[$hspcount]->{'_end'}  =$hsp->end('query');   
			  	$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_2_2_Hsp_2StEd_2Ht_1St'}        =$hsp->start('hit');         $hitCoverHash->[$hspcount]->{'_start'}  =$hsp->start('hit');     
			  	$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_2_3_Hsp_2StEd_2Ht_2Ed'}        =$hsp->end('hit');           $hitCoverHash->[$hspcount]->{'_end'}    =$hsp->end('hit');       
			  	

			  	
			  	$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_3_0_HspConserverd'}      =$hsp->frac_conserved ;
			  	$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_3_1_HspIdentical'}       =$hsp->frac_identical ;
			  	



			  	
			  	#$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'_Hsp_1Qry_string'}      =$hsp->query_string ;    #oldVersion
			  	#$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'_Hsp_3Hit_string'}      =$hsp->hit_string ;      #oldVersion
			  	#$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'_Hsp_2Hom_string'}      =$hsp->homology_string ; #oldVersion
			  	
			  	my $Qr_aln_string=
			  	$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_4_1_Hsp_3WhSting_1Qr'}      =$hsp->query_string ;
			  	my $mc_aln_string=
			  	$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_4_2_Hsp_3WhSting_2Ho'}      =$hsp->homology_string ;
			  	my $Ht_aln_string=
			  	$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'4_0_0_hspArray'}->[$hspcount]->{'5_4_3_Hsp_3WhSting_3Ht'}      =$hsp->hit_string ;
			  	
			  	
			  	
			  	$totalIdentical+=$hsp->frac_identical ;
			  	$totalConserved+=$hsp->frac_conserved ;
			  	
			  	my $smallQNB=$hsp->start('query'); my $bigQNB=$hsp->end('query');
			  	if ($smallQNB > $bigQNB ){ $bigQNB=$hsp->start('query'); my $smallQNB=$hsp->end('query'); }
			  	
			  	if ( $smallQNB < $queryMinStart ){  $queryMinStart = $smallQNB;   }
			  	if ( $bigQNB   > $queryMaxEnd   ){  $queryMaxEnd   = $bigQNB;     }
			  	  
		      print  "%\n";
          $hspcount++;
        }
        
        $queryCoverHash=BlastHandle::overLayerHadle( $queryCoverHash );  
        $hitCoverHash  =BlastHandle::overLayerHadle( $hitCoverHash ); 
        
        my $queryCoverLength=BlastHandle::GetCuttedLength($queryCoverHash);
        my $hitCoverLength=BlastHandle::GetCuttedLength($hitCoverHash);
        
        my ($hit_queryStart, $hit_queryEnd)=@{ BlastHandle::getWholeMatchedStartEnd($queryCoverHash) };
        my ($hit_HitStart,   $hit_HitEnd)  =@{ BlastHandle::getWholeMatchedStartEnd($hitCoverHash)   };
        
        
        #$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_1_0_queryCoverHash'}  =$queryCoverHash;        
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_1_1_queryCoverLength'}=$queryCoverLength;
        
        #$prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_1_3_hitCoverHash'}    =$hitCoverHash;
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_1_4_hitCoverLength'}  =$hitCoverLength;
       
        
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_1_5_hit_queryStart'}      =$hit_queryStart;
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_1_6_hit_queryEnd'}        =$hit_queryEnd;
        my $hit_queryRangeLength =(  (abs ($hit_queryEnd-$hit_queryStart)  ) + 1   );
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_1_7_hit_queryRangeLength'}=$hit_queryRangeLength;
        
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_1_8_hit_HitStart'}      =$hit_HitStart;
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_1_9_hit_HitEnd'}        =$hit_HitEnd;
        my $hit_HitRangeLength   =(  (abs ($hit_HitEnd-$hit_HitStart)  ) + 1   );
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_1_a_hit_HitRangeLength'}=$hit_HitRangeLength;
        
        
        
        if ($hspcount > 0){
          $totalIdentical=$totalIdentical/$hspcount;
          $totalConserved=$totalConserved/$hspcount;
          $totConsePercent=100*$totalConserved;  $totConsePercent=sprintf "%.2f",$totConsePercent; $totConsePercent="$totConsePercent%";
          
          $queryCoverRate=$queryCoverLength/$hit_queryRangeLength;  
          $qureyCovRatePct   =100*$queryCoverRate;  $qureyCovRatePct=sprintf "%.2f",$qureyCovRatePct; $qureyCovRatePct="$qureyCovRatePct%";
          
          $total_qry_covrAte=$queryCoverLength/$query_length;
          $total_qry_covrAte_pct=100*$total_qry_covrAte;  $total_qry_covrAte_pct=sprintf "%.2f",$total_qry_covrAte_pct; $total_qry_covrAte_pct="$total_qry_covrAte_pct%";
          
        
          $hitCoverRate=$hitCoverLength/$hit_HitRangeLength;
          $hitCovRatePct     =100*$hitCoverRate;  $hitCovRatePct=sprintf "%.2f",$hitCovRatePct; $hitCovRatePct="$hitCovRatePct%";
          
          
        }
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_2_0_hitTotalIdentical'}=$totalIdentical;
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_2_1_hitTotalConserved'}=$totalConserved;
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_2_2_totConsePercent'   }=$totConsePercent;
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_2_3_queryCoverRate'    }=$queryCoverRate;
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_2_4_qureyCovRatePct'   }=$qureyCovRatePct;
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_2_5_total_qry_covrAte' }=$total_qry_covrAte;
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_2_6_total_qry_covrAte_pct'}=$total_qry_covrAte_pct;
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_2_7_hitCoverRate'      }=$hitCoverRate;
        $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_2_8_hitCovRatePct'     }=$hitCovRatePct;
        
        
        if ($summuryHereLength<$maxBestinformToget){
        	my $singleSum="${accIDinShort}($prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_Z_hitArray'}->[$count]->{'3_0_6_hitEvalue'}|$totConsePercent)";
          $summuryHere.=$singleSum;
          $summuryHereLength++;
          $QueryHit_excel_HASH->{$queryName}->{$HitName}=$singleSum;
        }
        
       
        
			  
			  
			  $count++;
			  
			  
			  
			  
        
        # flow control for the number of hits needed
        if (  (defined ($numHits) ) && ($numHits=~m/\d+/) && ($numHits>0)   ){
			    last if ($count == $numHits);
			  }
		  }
    }                                                       
                                                         
    $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_3_QueryMinStart'}=$queryMinStart;
    $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_4_QueryMaxEnd'}  =$queryMaxEnd;
    $prasedOUThash->{'0_0_1_Query_array'}->[$resultArrayIdx]->{'2_0_5_Summury'}      =$summuryHere;
    #warn "\$resultArrayIdx=$resultArrayIdx\n";
    $resultArrayIdx++;	
    
    
    
  }
  

  

  
  $prasedOUThash->{'0_0_0_0_totalQry'}=$resultArrayIdx;
  $prasedOUThash->{'0_0_0_1_totalHit'}=$totalHit;
  $prasedOUThash->{'0_1_0_CheckFinishHASH'}=$Check_finish_HASH;
  
  
 
  return $prasedOUThash;
 
 
  

}



sub New_GetAllHit_HASH_fromBlastRstArrayHash{   # my $AllHit_Array=BlastPraser::New_GetAllHit_HASH_fromBlastRstArrayHash( $blastRstARRAYhash );
	my ($blastRstARRAYhash)=@_;
	
	my $warnMsgBody="\nIn package  BlastPraser,\tIn sub New_GetAllHit_HASH_fromBlastRstArrayHash,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  my $subCallereIfm=DirFileHandle::print_SubCallerInform; #warn $subCallereIfm;
  
  my $AllHit_Hash; my $AllHit_Array;
  if (   (  defined ( $blastRstARRAYhash )  ) && (  ref ( $blastRstARRAYhash ) eq 'HASH'  ) && (  defined ( $blastRstARRAYhash->{'0_0_1_Query_array'} )  ) && (  ref ( $blastRstARRAYhash->{'0_0_1_Query_array'} ) eq 'ARRAY'  )   ){
  	my $QueryIdx=0;                                                                                                                                                                                                                                #print "20190325-0-0-0-1 \$QueryIdx=$QueryIdx\n";
  	foreach my $eachQueryHASH (  @{ $blastRstARRAYhash->{'0_0_1_Query_array'} }  ){                                                                                                                                                                #print "20190325-0-0-0-2 \$QueryIdx=$QueryIdx\n";
  		if (   (  defined ( $eachQueryHASH )  ) && (  ref ( $eachQueryHASH ) eq 'HASH'  ) && (  defined ($eachQueryHASH->{'2_0_0_0_query_name'}) ) && ($eachQueryHASH->{'2_0_0_0_query_name'}=~m/\S+/)   ){
  		  my $queryName=$eachQueryHASH->{'2_0_0_0_query_name'};                                                                                                                                                                                      #print "20190325-0-0-0-3 \$queryName=\$eachQueryHASH->{'2_0_0_0_query_name'}=$eachQueryHASH->{'2_0_0_0_query_name'}\n";
  		  if (   (  defined ( $eachQueryHASH->{'2_0_Z_hitArray'} )  ) && (  ref ( $eachQueryHASH->{'2_0_Z_hitArray'} ) eq 'ARRAY'  )   ){
  		  	my $HitIdx=0;                                                                                                                                                                                                                            #print "20190325-0-0-0-4 \$HitIdx=$HitIdx\n";
  		  	foreach my $eachHitHASH (  @{ $eachQueryHASH->{'2_0_Z_hitArray'} }  ){                                                                                                                                                                   #print "20190325-0-0-0-5 \$HitIdx=$HitIdx\n";
  		  	  if (   (  defined ( $eachHitHASH )  ) && (  ref ( $eachHitHASH ) eq 'HASH'  ) && (  defined ($eachHitHASH->{'3_0_3_accessionNB'}) ) && ($eachHitHASH->{'3_0_3_accessionNB'}=~m/\S+/)   ){
  		  	  	my $hitName=$eachHitHASH->{'3_0_3_accessionNB'};                                                                                                                                                                                     #print "20190325-0-0-0-6 \$hitName=$hitName=\$eachHitHASH->{'3_0_3_accessionNB'}=$eachHitHASH->{'3_0_3_accessionNB'}\n";
  		  	  	my $gidName=$eachHitHASH->{'3_0_2_giNub'}; 
  		  	  	
  		  	  	if  (  defined ( $AllHit_Hash->{$hitName} )  ) {  		  	  		
  		  	  	}
  		  	  	else{
  		  	  		push @{ $AllHit_Array }, $hitName;                              #print "20190325-0-0-0-7 \$hitName=$hitName\n";
  		  	  	  $AllHit_Hash->{$hitName}=$gidName;                              #print "20190325-0-0-0-8 \$AllHit_Hash->{$hitName}=$AllHit_Hash->{$hitName}\n"; 
  		  	  	}
  		  	  	
  		  	  	
  		  	  }
  		  	  $HitIdx++;
  		  	}
  		  }
  		  
  		}
  		$QueryIdx++;
  	}
  }
	return $AllHit_Hash;
}

#Get the hash (in which the AccessionNB as KEY, Gi number as value) from the HASH (in which only the best hit was recorded for each query)
sub Get_Acc_to_Gi_HASH_from_Qry_bestHit_HASH{   # my $AllHit_Array=BlastPraser::Get_Acc_to_Gi_HASH_from_Qry_bestHit_HASH( $Qry_bestHit_HASH );
	my ($Qry_bestHit_HASH)=@_;
	
	my $warnMsgBody="\nIn package  BlastPraser,\tIn sub Get_Acc_to_Gi_HASH_from_Qry_bestHit_HASH,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  my $subCallereIfm=DirFileHandle::print_SubCallerInform; #warn $subCallereIfm;
  
  my $AllHit_Hash; my $AllHit_Array;
  DieWork::Check_Hash_or_DIE( $Qry_bestHit_HASH, "\$Qry_bestHit_HASH", $die_MsgHead, $subCallereIfm  );
	my $QueryIdx=0;  
	foreach my $eachQueryKEY (    sort { $a cmp $b } (   keys (  %{ $Qry_bestHit_HASH } )   )    ){                                                                                                                                                                                                                              #print "20190325-0-0-0-1 \$QueryIdx=$QueryIdx\n";  
    DieWork::Check_Hash_or_DIE( $Qry_bestHit_HASH->{$eachQueryKEY}, "\$Qry_bestHit_HASH->{$eachQueryKEY}", $die_MsgHead, $subCallereIfm  );
    my $QryHASH=$Qry_bestHit_HASH->{$eachQueryKEY};
    DieWork::Check_DfdNoEmptString_or_DIE( $QryHASH->{'3_0_3_accessionNB'}, "\$QryHASH->{'3_0_3_accessionNB'}", $die_MsgHead, $subCallereIfm  );
    DieWork::Check_DfdNoEmptString_or_DIE( $QryHASH->{'3_0_2_giNub'},       "\$QryHASH->{'3_0_2_giNub'}",       $die_MsgHead, $subCallereIfm  );
    
    my $hitName=$QryHASH->{'3_0_3_accessionNB'};                                                                                                                                                                                     #print "20190325-0-0-0-6 \$hitName=$hitName=\$eachHitHASH->{'3_0_3_accessionNB'}=$eachHitHASH->{'3_0_3_accessionNB'}\n";
  	my $gidName=$QryHASH->{'3_0_2_giNub'}; 
  	
  	if  (  defined ( $AllHit_Hash->{$hitName} )  ) {
  	}
  	else{
  		push @{ $AllHit_Array }, $hitName;                              #print "20190325-0-0-0-7 \$hitName=$hitName\n";
  	  $AllHit_Hash->{$hitName}=$gidName;                              #print "20190325-0-0-0-8 \$AllHit_Hash->{$hitName}=$AllHit_Hash->{$hitName}\n"; 
  	}
    
  }       
  
	return $AllHit_Hash;
}


sub New_GetAllHit_Array_fromBlastRstArrayHash{   # my $AllHit_Array=BlastPraser::New_GetAllHit_Array_fromBlastRstArrayHash( $blastRstARRAYhash );
	my ($blastRstARRAYhash)=@_;
	
	my $warnMsgBody="\nIn package  BlastPraser,\tIn sub New_GetAllHit_Array_fromBlastRstArrayHash,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  my $subCallereIfm=DirFileHandle::print_SubCallerInform; #warn $subCallereIfm;
  
  my $AllHit_Hash; my $AllHit_Array;
  if (   (  defined ( $blastRstARRAYhash )  ) && (  ref ( $blastRstARRAYhash ) eq 'HASH'  ) && (  defined ( $blastRstARRAYhash->{'0_0_1_Query_array'} )  ) && (  ref ( $blastRstARRAYhash->{'0_0_1_Query_array'} ) eq 'ARRAY'  )   ){
  	my $QueryIdx=0;                                                                                                                                                                                                                           #print "20190325-0-0-0-1 \$QueryIdx=$QueryIdx\n";
  	foreach my $eachQueryHASH (  @{ $blastRstARRAYhash->{'0_0_1_Query_array'} }  ){                                                                                                                                                                #print "20190325-0-0-0-2 \$QueryIdx=$QueryIdx\n";
  		if (   (  defined ( $eachQueryHASH )  ) && (  ref ( $eachQueryHASH ) eq 'HASH'  ) && (  defined ($eachQueryHASH->{'2_0_0_0_query_name'}) ) && ($eachQueryHASH->{'2_0_0_0_query_name'}=~m/\S+/)   ){
  		  my $queryName=$eachQueryHASH->{'2_0_0_0_query_name'};                                                                                                                                                                                        #print "20190325-0-0-0-3 \$queryName=\$eachQueryHASH->{'2_0_0_0_query_name'}=$eachQueryHASH->{'2_0_0_0_query_name'}\n";
  		  if (   (  defined ( $eachQueryHASH->{'2_0_Z_hitArray'} )  ) && (  ref ( $eachQueryHASH->{'2_0_Z_hitArray'} ) eq 'ARRAY'  )   ){
  		  	my $HitIdx=0;                                                                                                                                                                                                                       #print "20190325-0-0-0-4 \$HitIdx=$HitIdx\n";
  		  	foreach my $eachHitHASH (  @{ $eachQueryHASH->{'2_0_Z_hitArray'} }  ){                                                                                                                                                                   #print "20190325-0-0-0-5 \$HitIdx=$HitIdx\n";
  		  	  if (   (  defined ( $eachHitHASH )  ) && (  ref ( $eachHitHASH ) eq 'HASH'  ) && (  defined ($eachHitHASH->{'3_0_3_accessionNB'}) ) && ($eachHitHASH->{'3_0_3_accessionNB'}=~m/\S+/)   ){
  		  	  	my $hitName=$eachHitHASH->{'3_0_3_accessionNB'};                                                                                                                                                                                     print "20190325-0-0-0-6 \$hitName=$hitName=\$eachHitHASH->{'3_0_3_accessionNB'}=$eachHitHASH->{'3_0_3_accessionNB'}\n";
  		  	  	my $gidName=$eachHitHASH->{'3_0_2_giNub'}; 
  		  	  	
  		  	  	if (   (  defined ( $AllHit_Hash->{$hitName} )  ) && ( $AllHit_Hash->{$hitName} == 1 )   ){  		  	  		
  		  	  	}
  		  	  	else{
  		  	  		push @{ $AllHit_Array }, $hitName;                       print "20190325-0-0-0-7 \$hitName=$hitName\n";
  		  	  	  $AllHit_Hash->{$hitName}=$gidName;                              print "20190325-0-0-0-8 \$AllHit_Hash->{$hitName}=$AllHit_Hash->{$hitName}\n"; 
  		  	  	}
  		  	  	
  		  	  	
  		  	  }
  		  	  $HitIdx++;
  		  	}
  		  }
  		  
  		}
  		$QueryIdx++;
  	}
  }
	return $AllHit_Array;
}

sub New_ChangeBlastResult_into_Qry_Hit_HASH{ # my $blastQueryHit_hash=BlastPraser::New_ChangeBlastResult_into_Qry_Hit_HASH( $blastRstARRAYhash );
	my ($blastRstARRAYhash)=@_;
	
	my $warnMsgBody="\nIn package  BlastPraser,\tIn sub New_ChangeBlastResult_into_Qry_Hit_HASH,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  my $subCallereIfm=DirFileHandle::print_SubCallerInform; #warn $subCallereIfm;
                  
  my $blastQueryHit_hash;                                                                                                                                                                                                                     #print "20190124-0-0-0-0 \$blastQueryHit_hash=$blastQueryHit_hash\n";
  if (   (  defined ( $blastRstARRAYhash )  ) && (  ref ( $blastRstARRAYhash ) eq 'HASH'  ) && (  defined ( $blastRstARRAYhash->{'0_0_1_Query_array'} )  ) && (  ref ( $blastRstARRAYhash->{'0_0_1_Query_array'} ) eq 'ARRAY'  )   ){
  	my $QueryIdx=0;                                                                                                                                                                                                                           #print "20190124-0-0-0-1 \$QueryIdx=$QueryIdx\n";
  	foreach my $eachQueryHASH (  @{ $blastRstARRAYhash->{'0_0_1_Query_array'} }  ){                                                                                                                                                                #print "20190124-0-0-0-2 \$QueryIdx=$QueryIdx\n";
  		if (   (  defined ( $eachQueryHASH )  ) && (  ref ( $eachQueryHASH ) eq 'HASH'  ) && (  defined ($eachQueryHASH->{'2_0_0_0_query_name'}) ) && ($eachQueryHASH->{'2_0_0_0_query_name'}=~m/\S+/)   ){
  		  my $queryName=$eachQueryHASH->{'2_0_0_0_query_name'};                                                                                                                                                                                        #print "20190124-0-0-0-3 \$queryName=\$eachQueryHASH->{'2_0_0_0_query_name'}=$eachQueryHASH->{'2_0_0_0_query_name'}\n";
  		  if (   (  defined ( $eachQueryHASH->{'2_0_Z_hitArray'} )  ) && (  ref ( $eachQueryHASH->{'2_0_Z_hitArray'} ) eq 'ARRAY'  )   ){
  		  	my $HitIdx=0;                                                                                                                                                                                                                       #print "20190124-0-0-0-4 \$HitIdx=$HitIdx\n";
  		  	foreach my $eachHitHASH (  @{ $eachQueryHASH->{'2_0_Z_hitArray'} }  ){                                                                                                                                                                   #print "20190124-0-0-0-5 \$HitIdx=$HitIdx\n";
  		  	  if (   (  defined ( $eachHitHASH )  ) && (  ref ( $eachHitHASH ) eq 'HASH'  ) && (  defined ($eachHitHASH->{'3_0_3_accessionNB'}) ) && ($eachHitHASH->{'3_0_3_accessionNB'}=~m/\S+/)   ){
  		  	  	my $hitName=$eachHitHASH->{'3_0_3_accessionNB'};                                                                                                                                                                                     print "20190124-0-0-0-6 \$hitName=$hitName=\$eachHitHASH->{'3_0_3_accessionNB'}=$eachHitHASH->{'3_0_3_accessionNB'}\n";
  		  	  	$blastQueryHit_hash->{$queryName}->{$hitName}=Storable::dclone( $eachHitHASH );
  		  	  	
  		  	  	$blastQueryHit_hash->{$queryName}->{$hitName}->{'0_0_0_0_0_0_Query___name'}=$queryName;
  		  	  	$blastQueryHit_hash->{$queryName}->{$hitName}->{'0_0_0_0_0_1_Hit_____name'}=$hitName;
  		  	  	$blastQueryHit_hash->{$queryName}->{$hitName}->{'0_0_0_0_1_Query___length'}=$eachQueryHASH->{'2_0_1_query_length'};
  		  	  	$blastQueryHit_hash->{$queryName}->{$hitName}->{'0_0_0_0_2_Query__Min_Stt'}=$eachQueryHASH->{'2_0_3_QueryMinStart'};
  		  	  	$blastQueryHit_hash->{$queryName}->{$hitName}->{'0_0_0_0_3_Query__Max_End'}=$eachQueryHASH->{'2_0_4_QueryMaxEnd'};
  		  	  	$blastQueryHit_hash->{$queryName}->{$hitName}->{'0_0_0_0_4_0_QryFastaFile'}=$eachQueryHASH->{'2_0_0_1_query_fastaFile'};
  		  	  	$blastQueryHit_hash->{$queryName}->{$hitName}->{'0_0_0_0_4_1_QryIndexFile'}=$eachQueryHASH->{'2_0_0_2_query_indexFile'};
  		  	  	$blastQueryHit_hash->{$queryName}->{$hitName}->{'0_0_0_0_5_Query__Summury'}=$eachQueryHASH->{'2_0_5_Summury'};
  		  	  	$blastQueryHit_hash->{$queryName}->{$hitName}->{'0_0_0_0_6_Query_hits_num'}=$eachQueryHASH->{'2_0_2_0_num_hits'};
  		  	  	
  		  	  	$blastQueryHit_hash->{$queryName}->{$hitName}->{'0_0_0_1_0_Query_arreyIdx'}=$QueryIdx;
  		  	  	$blastQueryHit_hash->{$queryName}->{$hitName}->{'0_0_0_1_1_Hit___arreyIdx'}=$HitIdx;
  		  	  	
  		  	  }
  		  	  $HitIdx++;
  		  	}
  		  }
  		  
  		}
  		$QueryIdx++;
  	}
  }
	return $blastQueryHit_hash;
}


#my $blastQueryHit_hash=BlastPraser::GetTheBestHit_for_eachQuery_basedON( $in_QueryHit_HASH, $sort_Hit_Key );
#GetTheBestHit_for_eachQuery_basedON    '3_0_6_hitEvalue' => '4.0e-256', '3_0_7_hitScore' => '892.9','3_2_0_hitTotalIdentical' => 1, '3_2_1_hitTotalConserved' => 1,'3_2_5_total_qry_covrAte' => 1,'3_2_7_hitCoverRate' => 1,
sub GetTheBestHit_for_eachQuery_basedON{ 
	my ($in_QueryHit_HASH, $sort_Hit_Key)=@_;
	
	my $warnMsgBody="\nIn package  BlastPraser,\tIn sub GetTheBestHit_for_eachQuery_basedON,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  my $subCallereIfm=DirFileHandle::print_SubCallerInform; #warn $subCallereIfm;
	
	if  (   (  defined ( $in_QueryHit_HASH )  ) && (  ref ( $in_QueryHit_HASH ) eq 'HASH'  )   ){}
  else{		DieWork::Just_dieWork( $die_MsgHead."\n \$in_QueryHit_HASH=$in_QueryHit_HASH should be a HASH  !!  $!\n\n\n".$subCallereIfm ); 	}
	
	my $sort_reverse_or_not=0;
	if   (      (  defined ( $sort_Hit_Key )  ) && ( $sort_Hit_Key=~m/\S+/ ) 
	         && (      ( $sort_Hit_Key eq '3_0_6_hitEvalue'         ) || ( $sort_Hit_Key eq '3_0_7_hitScore'          ) 
	               ||  ( $sort_Hit_Key eq '3_2_0_hitTotalIdentical' ) || ( $sort_Hit_Key eq '3_2_1_hitTotalConserved' ) 
	               ||  ( $sort_Hit_Key eq '3_2_5_total_qry_covrAte' ) || ( $sort_Hit_Key eq '3_2_7_hitCoverRate' )
	            )
	     )
	{  
		if( $sort_Hit_Key eq '3_0_6_hitEvalue'         ){ $sort_reverse_or_not=1; }
	}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$sort_Hit_Key=$sort_Hit_Key should be a defined not empty file  !!  $!\n\n\n".$subCallereIfm ); 	}
	
	my $blastQueryHit_hash;                                                                                                                                                                                                                     #print "20190124-0-0-0-0 \$blastQueryHit_hash=$blastQueryHit_hash\n";
  foreach my $eachQueryKey (    sort { $a cmp $b } (   keys (  %{ $in_QueryHit_HASH } )   )    ){  
  	
  	if  (   (  defined ( $in_QueryHit_HASH->{$eachQueryKey} )  ) && (  ref ( $in_QueryHit_HASH->{$eachQueryKey} ) eq 'HASH'  )   ){}
    else{		DieWork::Just_dieWork( $die_MsgHead."\n \$in_QueryHit_HASH->{$eachQueryKey}=$in_QueryHit_HASH->{$eachQueryKey} should be a HASH  !!  $!\n\n\n".$subCallereIfm ); 	}
  	
  	my @HitKeyQrray=(   keys (  %{ $in_QueryHit_HASH->{$eachQueryKey} } )   );
  	my @sorted_hit_key_array;
  	if    ( $sort_reverse_or_not==0 ){  @sorted_hit_key_array=sort { $in_QueryHit_HASH->{$eachQueryKey}->{$b}->{$sort_Hit_Key} <=> $in_QueryHit_HASH->{$eachQueryKey}->{$a}->{$sort_Hit_Key} } @HitKeyQrray;  }
  	elsif ( $sort_reverse_or_not==1 ){  @sorted_hit_key_array=sort { $in_QueryHit_HASH->{$eachQueryKey}->{$a}->{$sort_Hit_Key} <=> $in_QueryHit_HASH->{$eachQueryKey}->{$b}->{$sort_Hit_Key} } @HitKeyQrray;  }
  	
  	my $stepNB=0;
  	FOREACHMARK: foreach my $each_Hit_Key ( @sorted_hit_key_array ){
      $blastQueryHit_hash->{$eachQueryKey}=Storable::dclone( $in_QueryHit_HASH->{$eachQueryKey}->{$each_Hit_Key} );
      $stepNB++;
      if ($stepNB>0){
      	last FOREACHMARK;
      }
  	}
  	
  }
	return $blastQueryHit_hash;
}

  
sub New_ChangeBlastResult_into_Hit_Qry_HASH{ # my $blastQueryHit_hash=BlastPraser::New_ChangeBlastResult_into_Hit_Qry_HASH( $blastRstARRAYhash );
	my ($blastRstARRAYhash)=@_;
	
	my $warnMsgBody="\nIn package  BlastPraser,\tIn sub New_ChangeBlastResult_into_Hit_Qry_HASH,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  my $subCallereIfm=DirFileHandle::print_SubCallerInform; #warn $subCallereIfm;
  
  my $blastQueryHit_hash=BlastPraser::New_ChangeBlastResult_into_Qry_Hit_HASH( $blastRstARRAYhash );
  my $blastHitQuery_hash;
  if (   (  defined ( $blastQueryHit_hash )  ) && (  ref ( $blastQueryHit_hash ) eq 'HASH'  )   ){
  	$blastHitQuery_hash=ArrayHashChange::Reverse_level1KEY_Level2Key_HASH($blastQueryHit_hash);
  }
  return $blastHitQuery_hash;
}





1;

##########################################################################################################################################
# 


