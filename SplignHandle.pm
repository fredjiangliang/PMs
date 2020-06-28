
#!/usr/bin/perl -w

BEGIN{  push (@INC, '/home/fredjiang/fredProgs/PMs')  };

use strict;
use warnings;

use PrintSubArrayHash;
use FastaFileHandle;
use SeqSegmentsTools;
use BlastHandle;
use DirFileHandle;
use Storable ;
use Storable qw(nstore store_fd nstore_fd freeze thaw dclone);

package  SplignHandle;



my ($ctg_stt_key,     $ctg_end_key,    $ctg_len_key,      $old_stt_key,     $old_end_key,     $old_DNA_key,     $old_PEP_key,   $new_stt_key,     $new_end_key,     $new_DNA_key,     $new_PEP_key,     $newDNA1lkey,             $newPEP1lkey            )
  =('7_0_ctg_0_stt',  '7_0_ctg_1_end', '7_0_ctg_2_len',   '7_1_old_0_stt',  '7_1_old_1_end', '7_1_old_2_DNA',  '7_1_old_3_PEP', '7_1_new_0_stt',  '7_1_new_1_end',  '7_1_new_2_DNA',  '7_1_new_3_PEP',  '7_1_new_4_oneLine_DNA',  '7_1_new_5_oneLine_PEP' );
#   ctgStt             ctgEnd          ctgLeth            old_stt            old_end                                             new_stt           new_end           newString 
  
my ($ups_stt_key,    $ups_end_key,     $ups_DNA_key,      $ups_PEP_key,     $dws_stt_key,    $dws_end_key,     $dws_DNA_key,      $dws_PEP_key   )
  =('8_0_ups_0_stt', '8_0_ups_1_end',  '8_0_ups_2_DNA',  '8_0_ups_3_PEP',   '8_1_dws_0_stt', '8_1_dws_1_end',  '8_1_dws_2_DNA',  '8_1_dws_3_PEP' );
#                                                                     
 
 my ($utr5_st_key,    $utr5_ed_key,     $utr5_DNA_key,       $utr3_st_key,   $utr3_ed_key,     $utr3_DNA_key,     $utr3_len_key,    )
  =('9_0_utr5_0_st', '9_0_utr5_1_ed',  '9_0_utr5_2_DNA',    '9_1_utr3_0_st', '9_1_utr3_1_ed',  '9_1_utr3_2_DNA', '9_1_utr3_3_len'   );
#                                                            utrStt           utrEnd            utrSequence        putative3UTRLength      
  
my ($upsATG_fKey,         $upsATGfMakK,          $upsSTP_fKey,         $upsSTPfMakK,          $dwsSTP_fKey,         $dwsSTPfMakK         )
  =('A_0_upsATG_0_found', 'A_0_upsATG_1_fdMak',  'A_1_upsSTP_0_found', 'A_1_upsSTP_1_fdMak',  'A_2_dwsSTP_0_found', 'A_2_dwsSTP_1_fdMak',);
#    upS__ATG_found         upS__ATG_found_mark   upS_stop_found         upS_stop_found_mark   dwS_stop_found        dwS_stop_found_mark



sub Splign_to_check_Geno_Est_MATCH{  #  SplignHandle::Splign_to_check_Geno_Est_MATCH
	my ($GenoDatabase, $GenoCtgName,                                     $Gno_proSplign_HASH, 
	                                                                                          $Est_Database, $Est_CtgName, $Est_proSplign_HASH, $work_DIR)=@_;
	#my ($GenoDatabase, $GenoCtgName, $GenoStt, $GenoEnd, $Gno_gene_ZoF, $Gno_proSplign_HASH, $Est_Database, $Est_CtgName, $Est_gene_ZoF, $work_DIR)=@_;
	
	
	my $warnMsgBody="\nIn package  SplignHandle,\tIn sub Splign_to_check_Geno_Est_MATCH,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $caller_inform=DirFileHandle::print_SubCallerInform;  
  
  
  
  my $out_HASH;
  
  my $ExpandLength=5000;
  
  #my $Gno_proSplign_HASH=retrieve ($Gno_proSplign_File);
  #my $Est_proSplign_HASH=retrieve ($Est_proSplign_File);
  
  if (      (  defined ( $Gno_proSplign_HASH                        )  ) && (  ref ( $Gno_proSplign_HASH                        ) eq 'HASH'   ) 
      	 && (  defined ( $Gno_proSplign_HASH->{'0_MatchArray'}      )  ) && (  ref ( $Gno_proSplign_HASH->{'0_MatchArray'}      ) eq 'ARRAY'  ) 
      	 && (  defined ( $Gno_proSplign_HASH->{'0_MatchArray'}->[0] )  ) && (  ref ( $Gno_proSplign_HASH->{'0_MatchArray'}->[0] ) eq 'HASH'   ) 
      	 && (  defined ( $Est_proSplign_HASH                        )  ) && (  ref ( $Est_proSplign_HASH                        ) eq 'HASH'   ) 
      	 && (  defined ( $Est_proSplign_HASH->{'0_MatchArray'}      )  ) && (  ref ( $Est_proSplign_HASH->{'0_MatchArray'}      ) eq 'ARRAY'  ) 
      	 && (  defined ( $Est_proSplign_HASH->{'0_MatchArray'}->[0] )  ) && (  ref ( $Est_proSplign_HASH->{'0_MatchArray'}->[0] ) eq 'HASH'   ) 
     )
  {
    my $org_GenoStt =$Gno_proSplign_HASH->{'0_MatchArray'}->[0]->{'0_1_aln_0_4_0_Dps_1stf1'}; 
    my $org_GenoEnd =$Gno_proSplign_HASH->{'0_MatchArray'}->[0]->{'0_1_aln_0_5_0_Dps_Lstf3'}; 
    
    my $GenoStt=SeqSegmentsTools::getSmallOne( $org_GenoStt, $org_GenoEnd); my $expandGnoStt=$GenoStt-$ExpandLength;
    my $GenoEnd=SeqSegmentsTools::getbigOne  ( $org_GenoStt, $org_GenoEnd); my $expandGnoEnd=$GenoEnd+$ExpandLength;
    
    my $Gno_gene_ZoF=$Gno_proSplign_HASH->{'0_MatchArray'}->[0]->{'2_DNA_ZF'};    
    my $Est_gene_ZoF=$Est_proSplign_HASH->{'0_MatchArray'}->[0]->{'2_DNA_ZF'};  
    
    my $GenoFatchARRAY = FastaFileHandle::FactchSeqBy_FastaCMD($GenoDatabase, $GenoCtgName, '+', $expandGnoStt, $expandGnoEnd); 
    my $Est_FatchARRAY = FastaFileHandle::FactchSeqBy_FastaCMD($Est_Database, $Est_CtgName                                   ) ; 
    if  (   (  defined ( $GenoFatchARRAY )  ) && (  ref ( $GenoFatchARRAY ) eq 'ARRAY'  ) && (  defined ( $Est_FatchARRAY )  ) && (  ref ( $Est_FatchARRAY ) eq 'ARRAY'  )   ){
    	my ( $GenoExtractSFasta, $Geno_Beyond_orNot, $Geno_Exr_real_ZoF, $Geno_Exr_real_Stt, $Geno_Exr_real_end, $Geno_Exr_Seqlenth ) = @{ $GenoFatchARRAY };
    	my ( $Est_ExtractSFasta, $Est__Beyond_orNot, $Est__Exr_real_ZoF, $Est__Exr_real_Stt, $Est__Exr_real_end, $Est__Exr_Seqlenth ) = @{ $Est_FatchARRAY }; 
    	system  ("mkdir -p $work_DIR");
    	
    	my $Geno_Seq_file_0_0_0=$work_DIR."/0_0_0_Geno_Seq_file.txt";                                      my $Est__Seq_file_0_0_1=$work_DIR."/0_0_1_Est__Seq_file.txt";  
      FastaFileHandle::BuildFastaFile_withFastaString($Geno_Seq_file_0_0_0, $GenoExtractSFasta );        FastaFileHandle::BuildFastaFile_withFastaString($Est__Seq_file_0_0_1, $Est_ExtractSFasta);
      $out_HASH->{'0_0_0_Geno_Seq_file'}=$Geno_Seq_file_0_0_0;                                           $out_HASH->{'0_0_1_Est__Seq_file'}=$Est__Seq_file_0_0_1;
      
      my $SplignAlnOut_file_0_0_2=$work_DIR."/0_0_2_SplignAlnOut_file.txt";  
      $out_HASH->{'0_0_2_SplignAlnOut_file'}=$SplignAlnOut_file_0_0_2; 
      SplignHandle::runSplign($Est__Seq_file_0_0_1, $Geno_Seq_file_0_0_0, $SplignAlnOut_file_0_0_2);  print "201811221547-1\n";
    	
    	my $SplignOUTarray=SplignHandle::insplignAlnOutPraser($SplignAlnOut_file_0_0_2, $Geno_Exr_real_Stt, $Geno_Exr_real_end, $Geno_Exr_real_ZoF); print "201811221547-2\n";
      if ( ref($SplignOUTarray) eq 'ARRAY' ){    	 print "201811221547-3\n";
        $SplignOUTarray=SplignHandle::AddQueryLength_get_covIdt($SplignOUTarray, $Est__Exr_Seqlenth);    print "201811221547-4\n";
        	if ( ref($SplignOUTarray) eq 'ARRAY' ){    print "201811221547-5\n";
        		if (    
        		       (  defined ( $Gno_proSplign_HASH->{'0_MatchArray'}->[0]->{'7_0_represtCDSpos'} )  ) && (  ref ( $Gno_proSplign_HASH->{'0_MatchArray'}->[0]->{'7_0_represtCDSpos'} ) eq 'ARRAY'  ) 
        		   )
        	  {   print "201811221547-6\n";
        		  my $Gno_CDSreprestPosList=$Gno_proSplign_HASH->{'0_MatchArray'}->[0]->{'7_0_represtCDSpos'};
              $SplignOUTarray=SplignHandle::Add_KeyPosList_and_Sort_By_covIDt($SplignOUTarray, $Gno_CDSreprestPosList, $Est_gene_ZoF, $Gno_gene_ZoF);  
              if ( ref($SplignOUTarray) eq 'ARRAY' ){   
              	my $SplignOUTArrayfil_0_0_3=$work_DIR."/0_0_3_SplignOUTArrayfil.ary";    DirFileHandle::PrintDumper($SplignOUTArrayfil_0_0_3, $SplignOUTarray);
              	$out_HASH->{'0_0_3_SplignOUTArrayfil'}=$SplignOUTArrayfil_0_0_3; 
              	    
                #my $Simp_splignARRAY=SplignHandle::simplify_splignPraser_Out_Hash($SplignOUTarray);
                #if ( ref($Simp_splignARRAY) eq 'ARRAY' ){ 
              	#  my $SimpSplArray_file_0_0_4=$work_DIR."/0_0_4_SimpSplArray_file.ary";    DirFileHandle::PrintDumper($SimpSplArray_file_0_0_4, $Simp_splignARRAY);
              	#  $out_HASH->{'0_0_4_SimpSplArray_file'}=$SimpSplArray_file_0_0_4; 
                #}
              } 
            }  
        	}
      	
      	
      }  
    }
    else {
    	my $dieMsg1=$die_MsgHead.$caller_inform."\$GenoFatchARRAY=$GenoFatchARRAY or \$Est_FatchARRAY=$Est_FatchARRAY is not right!!! :$!\n\n\n\n";
    	print $dieMsg1;
    	die $dieMsg1;
    }
  }
  else {
   	my $dieMsg2=$die_MsgHead.$caller_inform."\$Gno_proSplign_HASH->{'0_MatchArray'}->[0]=$Gno_proSplign_HASH->{'0_MatchArray'}->[0] or \$Est_proSplign_HASH->{'0_MatchArray'}->[0]=$Est_proSplign_HASH->{'0_MatchArray'}->[0]  is not right!!! :$!\n\n\n\n";
   	print $dieMsg2;
  	die $dieMsg2;
  }
    
  return $out_HASH;
  
  

}



sub runSplign{
	my $warnMsgHead="\n\n\n   In package  SplignHandle,\tIn sub runSplign,\n\n";	
	my ($quryFile, $subjFile, $outAlnFile)=@_;
	if (     ( defined ($quryFile) ) && ( defined ($subjFile) ) && ( defined ($outAlnFile) ) 
          && ( $quryFile=~m/^\S+$/ ) && ( $subjFile=~m/^\S+$/ ) && ( $outAlnFile=~m/^\S+$/ )
       ){
  }
  else {
    my $dieMsg="\n\n\nDIE!!!!!\n$warnMsgHead\nThe input files are not right!!!!\n\$quryFile=$quryFile, \$subjFile=$subjFile, \$outAlnFile=$outAlnFile\n\n";
    print $dieMsg; die $dieMsg;
  }
	
	my $command="splign -subj $subjFile -query $quryFile -aln $outAlnFile";
	print $warnMsgHead, $command, "\n\n\n";
	warn  $warnMsgHead, $command, "\n\n\n";
	system ("$command");
}



sub insplignAlnOutPraser{      #解析 splign 输出aln文件 Prase the aln output file of splign
	
  my $warnMsgHead="\n\n\n   In package  SplignHandle,\tIn sub insplignAlnOutPraser,\n\n";	
	
  my ($inSplignOutFile, $sbjSeg_stt, $sbjSeg_end, $sbjSeg_ZoF, $qurSeg_stt, $qurSeg_end, $qurSeg_ZoF)=@_;  print "\$inSplignOutFile=$inSplignOutFile\n";  
  
  
  my $splGno_0='1_0splGno';                            my $sj_char_4='4sj_char';               #'3_sj_line'                         
	my $splExN_1='1_1splExN';                            my $Ex_numb_6='6Ex_numb';               #'4_ex_numb'                         
	my $spl_CS_2='1_2spl_CS';                            my $CS_char_3='3CS_char';               #'2_cs_line'                         
	my $splCtg_3='1_3splCtg';                            my $qr_char_2='2qr_char';               #'1_qr_line'                         
	my $spl_AA_4='1_4spl_AA';                            my $AA_char_1='1AA_char';               #'0_AA_line'                         
	my $ctgPos_a='Not_include_spl_ctg_posi';             my $qr_posi_0='0qr_posi';               #'Not_include_0qr_posi'
	my $gnoPos_b='Not_include_spl_gno_posi';             my $sj_posi_5='5sj_posi';               #'Not_include_5sj_posi'   
	my $ctgZoF_a='Not_include_spl_ctg_ZorF';             my $qr_ZorF_0='0qr_ZorF';               #  
	my $gnoZoF_b='Not_include_spl_gno_ZorF';             my $sj_ZorF_5='5sj_ZorF';               #
	
  
  
  my $qur_ChgBack=0;
  if (  ( defined ($qurSeg_stt) ) || ( defined ($qurSeg_end) ) || ( defined ($qurSeg_ZoF) )  ){
    if (     ( defined ($qurSeg_stt) ) && ( defined ($qurSeg_end) ) && ( defined ($qurSeg_ZoF) ) 
          && ( $qurSeg_stt=~m/^\d+$/ ) && ( $qurSeg_end=~m/^\d+$/ ) && ( $qurSeg_ZoF=~m/^(\+|-)$/ )
       ){
      $qur_ChgBack=1;
    }
    else {
    	my $dieMsg="\n\n\nDIE!!!!!\n$warnMsgHead\nThe qurSegmentInformation is not right!!!!\n\$qurSeg_stt=$qurSeg_stt, \$qurSeg_end=$qurSeg_end, \$qurSeg_ZoF=$qurSeg_ZoF\n\n";
    	print $dieMsg; die $dieMsg;
    }
  }
  
  my $sbj_ChgBack=0;
  if (  ( defined ($sbjSeg_stt) ) || ( defined ($sbjSeg_end) ) || ( defined ($sbjSeg_ZoF) )  ){
    if (     ( defined ($sbjSeg_stt) ) && ( defined ($sbjSeg_end) ) && ( defined ($sbjSeg_ZoF) ) 
          && ( $sbjSeg_stt=~m/^\d+$/ ) && ( $sbjSeg_end=~m/^\d+$/ ) && ( $sbjSeg_ZoF=~m/^(\+|-)$/ )
       ){
      $sbj_ChgBack=1;
    }
    else {
    	my $dieMsg="\n\n\nDIE!!!!!\n$warnMsgHead\nThe sbjSegmentInformation is not right!!!!\n\$sbjSeg_stt=$sbjSeg_stt, \$sbjSeg_end=$sbjSeg_end, \$sbjSeg_ZoF=$sbjSeg_ZoF\n\n";
    	print $dieMsg; die $dieMsg;
    }
  }
  
  my $returnedArray;  #用于最终输出的 array
 
    
  open (INSPLIGNFILE,$inSplignOutFile) or die "In Sub &praseCapOut, Cannot open \$inSplignOutFile=$inSplignOutFile : $!\n";
  my @inSplignOut=(<INSPLIGNFILE>);
  my $wholeInSplignOut=join ('',@inSplignOut);
  
  $wholeInSplignOut=$wholeInSplignOut."\n";
  
  
  #>+1	1sig2388(+)	DNAforGeneWISE(+)
  my $gene_head_Regular='>(?:\+|-)\d+\s+[^\(]+\((?:\+|-)\)\s+[^\(]+\((?:\+|-)\)\s*\n';
  #Exon 1 (1-280,8614-8887) Len = 285 Identity = 0.898
  my $exon_head_Regular='.*Exon\s+\d+\s+\(\d+-\d+,\d+-\d+\)\s+Len\s+=\s+\d+\s+Identity\s+=\s+\S+\n';
  my $mainBody_AAm_line='.*\n';
  my $mainBody_Qur_line='[ ]*\d+\s+\S+\n';
  my $mainBody_Con_line='.*\n';
  my $mainBody_Sbj_line='[ ]*\d+\s+\S+\n';
  
  #print "test\n";
  #if($wholeInSplignOut=~m/$gene_head_Regular\n($exon_head_Regular)($mainBody_AAm_line)($mainBody_Qur_line)($mainBody_Con_line)/) {print "1 \$1=$1((((((((((1)))))))\n\$2=$2(((((22(2)))))\n\$3=$3(((((33))))))))\n\$4=$4(((((((4))))n";}
  if ($wholeInSplignOut=~/\S/){
  	if ($wholeInSplignOut=~m/
                                
                              (?:
                                $gene_head_Regular  
                                \n 
                                  
                                (?:
                                 $exon_head_Regular
                                  
                                  (?:
                                     $mainBody_AAm_line
                                     $mainBody_Qur_line
                                     $mainBody_Con_line
                                     $mainBody_Sbj_line
                                     
                                     \n+                                  
                                  )+
                                   
                                )+
                                
                              )+
                                
                            /x
         )
    {                                                                                                             	#print "\$1=$1\n\$2=$2\n\$3=$3\n\$4=$4\n$5=$5\n\n";  
    	
    	
    	my @GeneRstArray=( $wholeInSplignOut=~m/
                             
                                               (
                                                 $gene_head_Regular  
                                                 \n 
                                                   
                                                 (?:
                                                   $exon_head_Regular
                                                   
                                                   (?:
                                                      $mainBody_AAm_line
                                                      $mainBody_Qur_line
                                                      $mainBody_Con_line
                                                      $mainBody_Sbj_line
                                                      
                                                      \n+                                  
                                                   )+
                                                    
                                                 )+
                                                 
                                               )
                                                 
                                             /xg 
    	                         );  	                         
    	my $i=0; #create a val to show cicle number 
    	foreach my $eachGeneRst (@GeneRstArray){   		
    		
    		
    		#print "\$i=$i\t\$eachGeneRst=$eachGeneRst\n\n"; 
    		#gene_head_Regular='\s*>(?:\+|-)\d+\s+[^\(]+\((?:\+|-)\)\s+[^\(]+\((?:\+|-)\)\s*\n';
    		
    		
    		#                        $1    $2      $3        $4         $5        $6  
    		if ($eachGeneRst=~m/>(\+|-)(\d+)\s+([^\(]+)\((\+|-)\)\s+([^\(]+)\((\+|-)\)\s*\n/){
        	my ($gene_ZoF, $gene_num, $Qury_nam, $Qury_ZoF, $Subj_nam, $Subj_ZoF)
        	=  ($1,        $2,        $3,        $4,        $5,        $6       );
        	
        	
          if ($qur_ChgBack==1){ $Qury_ZoF=SeqSegmentsTools::Chang_ZoF_BackToOrgContig_ZoF($Qury_ZoF, $qurSeg_ZoF); }
        	if ($sbj_ChgBack==1){ $Subj_ZoF=SeqSegmentsTools::Chang_ZoF_BackToOrgContig_ZoF($Subj_ZoF, $sbjSeg_ZoF); }
        	
        	
        	$returnedArray->[$i]->{'0_gene_num'}=$gene_num;     $returnedArray->[$i]->{'1_gene_ZoF'}=$gene_ZoF;      	
        	$returnedArray->[$i]->{'2_Qury_nam'}=$Qury_nam;     $returnedArray->[$i]->{'3_Qury_ZoF'}=$Qury_ZoF;
        	$returnedArray->[$i]->{'4_Subj_nam'}=$Subj_nam;     $returnedArray->[$i]->{'5_Subj_ZoF'}=$Subj_ZoF;
        	
        	
        	my $Qury_ZoF_one=1; if ($Qury_ZoF eq '-'){$Qury_ZoF_one=-1;}
        	my $Subj_ZoF_one=1; if ($Subj_ZoF eq '-'){$Subj_ZoF_one=-1;}
        	
        	my @ExonArray=(  $eachGeneRst=~m/
    		                                    (
                                              $exon_head_Regular
                                              
                                              (?:
                                                 $mainBody_AAm_line
                                                 $mainBody_Qur_line
                                                 $mainBody_Con_line
                                                 $mainBody_Sbj_line
                                                 
                                                 \n+                                  
                                              )+
                                               
                                            )
    		                                  /xmg
    		               
    		                );
    		    
    		  my $ex_whole_line_idx=0;              
    		  foreach my $eachExonRst (@ExonArray){  print "\$eachExonRst=$eachExonRst((((((\n";
    		  	
    		  	#Exon 1 (1-280,8614-8887) Len = 285 Identity = 0.898
            #my $exon_head_Regular='Exon\s+\d+\s+\(\d+-\d+,\d+-\d+\)\s+Len\s+=\s+\d+\s+Identity\s+=\s+\S+\n';
            #                           $1        $2    $3    $4    $5                  $6                     $7  
    		    if ($eachExonRst=~m/Exon\s+(\d+)\s+\((\d+)-(\d+),(\d+)-(\d+)\)\s+Len\s+=\s+(\d+)\s+Identity\s+=\s+(\S+)\n(.*)$/s){   
    		    	my ($ex_number, $ex_qr_stt, $ex_qr_end, $ex_sj_stt, $ex_sj_end, $ex_al_len, $ex_al_idt, $ex_al_aln)
        	    =  ($1,         $2,         $3,         $4,         $5,         $6,         $7,         $8        ); #warn "\$ex_al_aln=$ex_al_aln\n"; sleep(3);
        	    
        	    
        	    if ($qur_ChgBack == 1){
                $ex_qr_stt=SeqSegmentsTools::Chang_pos_BackToOrgContigPositionk($ex_qr_stt, $qurSeg_stt, $qurSeg_end, $qurSeg_ZoF);
        	      $ex_qr_end=SeqSegmentsTools::Chang_pos_BackToOrgContigPositionk($ex_qr_end, $qurSeg_stt, $qurSeg_end, $qurSeg_ZoF);
        	    }
        	    if ($sbj_ChgBack == 1){
                $ex_sj_stt=SeqSegmentsTools::Chang_pos_BackToOrgContigPositionk($ex_sj_stt, $sbjSeg_stt, $sbjSeg_end, $sbjSeg_ZoF);
        	      $ex_sj_end=SeqSegmentsTools::Chang_pos_BackToOrgContigPositionk($ex_sj_end, $sbjSeg_stt, $sbjSeg_end, $sbjSeg_ZoF);
        	    }
        	    
        	    $returnedArray->[$i]->{'6_exon_Hsh'}->{$ex_number}->{'0_ex_number'}=$ex_number;     
        	    $returnedArray->[$i]->{'6_exon_Hsh'}->{$ex_number}->{'1_ex_qr_stt'}=$ex_qr_stt;     
        	    $returnedArray->[$i]->{'6_exon_Hsh'}->{$ex_number}->{'2_ex_qr_end'}=$ex_qr_end;     
        	    $returnedArray->[$i]->{'6_exon_Hsh'}->{$ex_number}->{'3_ex_sj_stt'}=$ex_sj_stt;     
        	    $returnedArray->[$i]->{'6_exon_Hsh'}->{$ex_number}->{'4_ex_sj_end'}=$ex_sj_end;     
        	    $returnedArray->[$i]->{'6_exon_Hsh'}->{$ex_number}->{'5_ex_al_len'}=$ex_al_len;     
        	    $returnedArray->[$i]->{'6_exon_Hsh'}->{$ex_number}->{'6_ex_al_idt'}=$ex_al_idt;     
        	    
        	    
        	    
        	
    		    	my @AlignArray=(  $ex_al_aln=~m/
    		                                                                                
                                                (
                                                   $mainBody_AAm_line
                                                   $mainBody_Qur_line
                                                   $mainBody_Con_line
                                                   $mainBody_Sbj_line
                                                   
                                                   \n+                                  
                                                )                                            
                                              
    		                                      /xg
    		               
    		                    );
    		  
    		    	my $alinIdx=0;   my $last_qr_pos; my $last_sj_pos;
    		    	foreach my $eachAlignRst (@AlignArray){  #print "\$AlignArray=$AlignArray\n";
    		    	  #my $mainBody_AAm_line='\s+\S*(?:\s+\S+)*\S*\n';
                #my $mainBody_Qur_line='\s*\d+\s+\S+\n';
                #my $mainBody_Con_line='\s+\S+(?:\s+\S+)*\S*\n';
                #my $mainBody_Sbj_line='\s*\d+\s+\S+\n';
    		    	  #                    #$1 
    		    	  if ($eachAlignRst=~m/^(.*)\n
    		    	  #                   #$2  $3       #4
    		    	                       ([ ]*(\d+)\s+)(\S+)\n
    		    	  #                   #$5      
    		    	                       (.*)\n
    		    	  #                   #$6  $7       $8
    		    	                       ([ ]*(\d+)\s+)(\S+)\n    
    		    	                      /x
    		    	     )
    		    	  {  
    		    	    my ($aln_AA_WLine, $aln_qr_hdWid, $aln_qr_hdNub, $aln_qr_ACTGl, $aln_cs_Wline, $aln_sj_hdWid, $aln_sj_hdNub, $aln_sj_ACTGl)
        	        =  ($1,           $2,             $3,            $4,            $5,            $6,            $7,            $8           );
        	        
        	        #$qur_ChgBack $sbj_ChgBack $qurSeg_stt, $qurSeg_end, $qurSeg_ZoF, $sbjSeg_stt, $sbjSeg_end, $sbjSeg_ZoF
                  #($inPosition, $segmentCtgStt, $segmentCtgEnd, $segmentCtgZoF)
                  if ($qur_ChgBack == 1){
                    $aln_qr_hdNub=SeqSegmentsTools::Chang_pos_BackToOrgContigPositionk($aln_qr_hdNub, $qurSeg_stt, $qurSeg_end, $qurSeg_ZoF);      	          
        	        }
        	        if ($sbj_ChgBack == 1){
                    $aln_sj_hdNub=SeqSegmentsTools::Chang_pos_BackToOrgContigPositionk($aln_sj_hdNub, $sbjSeg_stt, $sbjSeg_end, $sbjSeg_ZoF);      	          
        	        }
        	        
        	        my $aln_qr_hdwLength=length ($aln_qr_hdWid); my $aln_sj_hdwLength=length ($aln_sj_hdWid);
        	        if ($aln_qr_hdwLength != $aln_sj_hdwLength){ my $dieMsg="\n\nDIE!!!!\n$warnMsgHead\n\$eachAlignRst=$eachAlignRst\n\$aln_qr_hdWid=$aln_qr_hdWid\n\$aln_sj_hdWid=$aln_sj_hdWid\n\$aln_qr_hdwLength=$aln_qr_hdwLength should == $aln_sj_hdwLength=\$aln_sj_hdwLength\n\n\n";print $dieMsg; die $dieMsg;}
        	        my $aln_AA_rlLin=substr($aln_AA_WLine,$aln_qr_hdwLength ); 
        	        my $aln_CS_rlLin=substr($aln_cs_Wline,$aln_qr_hdwLength ); 
                  
        	        my $Qrlength=length ($aln_qr_ACTGl);  my $WLlength=$Qrlength;
        	        my $SJlength=length ($aln_sj_ACTGl);  if ($SJlength != $WLlength){ my $dieMsg="\n\nDIE!!!!\n$warnMsgHead\n\$SJlength=$SJlength should == $WLlength=\$WLlength\n\n\n";print $dieMsg; die $dieMsg;}
        	        my $AAlength=length ($aln_AA_rlLin);  if ($AAlength != $WLlength){ my $dieMsg="\n\nDIE!!!!\n$warnMsgHead\n\$AAlength=$AAlength should == $WLlength=\$WLlength\n\n\n";print $dieMsg; die $dieMsg;}
        	        my $CSlength=length ($aln_CS_rlLin);  if ($CSlength != $WLlength){ my $dieMsg="\n\nDIE!!!!\n$warnMsgHead\n\$CSlength=$CSlength should == $WLlength=\$WLlength\n\n\n";print $dieMsg; die $dieMsg;}
        	        
        	        
        	        
        	        my $aln_qr_hd_pointLine='';   my $aln_qr_tl_pointLine='';
        	        if ($aln_qr_ACTGl=~m/^(\.+)$/){ $aln_qr_hd_pointLine=$1;      } 
        	        else{
        	        	if ($aln_qr_ACTGl=~m/^(\.+)[^\.]*/){ $aln_qr_hd_pointLine=$1;      } 
        	          if ($aln_qr_ACTGl=~m/[^\.]*(\.+)$/){ $aln_qr_tl_pointLine=$1;      } 
        	        }
        	        my $aln_qr_hdPoitLen=length ($aln_qr_hd_pointLine);
        	        my $aln_qr_tlPoitLen=length ($aln_qr_tl_pointLine);
        	        my $aln_qr_middle_Len=$Qrlength-$aln_qr_hdPoitLen-$aln_qr_tlPoitLen;
        	        if ($aln_qr_middle_Len<0){
        	        	my $dieMsg="\n\nDIE!!!!\n$warnMsgHead\n\$aln_qr_middle_Len=\$Qrlength-\$aln_qr_hdPoitLen-\$aln_qr_tlPoitLen}=$Qrlength-$aln_qr_hdPoitLen-$aln_qr_tlPoitLen=$aln_qr_middle_Len Should be >=0!!!!\n\n\n";print $dieMsg; die $dieMsg;
        	        }
        	        
        	        
        	        my @qrWoLineAr=split '',$aln_qr_ACTGl;
        	        my @sjWoLineAr=split '',$aln_sj_ACTGl;
        	        my @AAWoLineAr=split '',$aln_AA_rlLin;
        	        my @csWoLineAr=split '',$aln_CS_rlLin;
        	        
        	        my $wholeLine_Hash; $wholeLine_Hash=$returnedArray->[$i]->{'7_whoL_Hsh'} if (  defined ( $returnedArray->[$i]->{'7_whoL_Hsh'} )  );
        	        my $query_____Hash; $query_____Hash=$returnedArray->[$i]->{'8_qurM_Hsh'} if (  defined ( $returnedArray->[$i]->{'8_qurM_Hsh'} )  );
        	        my $subject___Hash; $subject___Hash=$returnedArray->[$i]->{'9_sbjM_Hsh'} if (  defined ( $returnedArray->[$i]->{'9_sbjM_Hsh'} )  );
        	        
        	        my $qr_pos_nb;               my $sj_pos_nb;
        	        my $qr_pos_add_nb=0;         my $sj_pos_add_nb=0;
        	        for (my $wlIdx=0; $wlIdx<$WLlength; $wlIdx++){
        	          my $qr_char=$wholeLine_Hash->{$ex_whole_line_idx}->{$qr_char_2}=$qrWoLineAr[$wlIdx];
        	          my $sj_char=$wholeLine_Hash->{$ex_whole_line_idx}->{$sj_char_4}=$sjWoLineAr[$wlIdx];
        	          my $AA_char=$wholeLine_Hash->{$ex_whole_line_idx}->{$AA_char_1}=$AAWoLineAr[$wlIdx];
        	          my $CS_char=$wholeLine_Hash->{$ex_whole_line_idx}->{$CS_char_3}=$csWoLineAr[$wlIdx];
        	          
        	          $wholeLine_Hash->{$ex_whole_line_idx}->{$Ex_numb_6}=$ex_number;  
        	          
     	          
        	          
        	          if    ($qr_char eq '.'){
        	          }
        	          elsif ( ($qr_char eq '-') || ($qr_char=~m/^[a-zA-Z]$/) ){
        	          	if ($qr_char=~m/^[a-zA-Z]$/){
        	          		$qr_pos_nb=$aln_qr_hdNub+$qr_pos_add_nb*$Qury_ZoF_one;
        	          		$qr_pos_add_nb++;
        	          	}  
        	          	if (  ( defined ($qr_pos_nb) ) && ( $qr_pos_nb=~m/\d+/ )   ){
        	          		
        	          	  $query_____Hash->{$qr_pos_nb}->{'0_self_qr_char'}=$qr_char;
        	          	  $query_____Hash->{$qr_pos_nb}->{'0_self_qr_ZorF'}=$Qury_ZoF;
        	          	  $query_____Hash->{$qr_pos_nb}->{'1_mach_sj_char'}=$sj_char;
        	          	  $query_____Hash->{$qr_pos_nb}->{'7_whoL_positio'}=$ex_whole_line_idx;     
        	          	  $query_____Hash->{$qr_pos_nb}->{'8_exon__number'}=$ex_number;
        	          	  
        	          	  $wholeLine_Hash->{$ex_whole_line_idx}->{$qr_posi_0}=$qr_pos_nb;  	
        	          	  $wholeLine_Hash->{$ex_whole_line_idx}->{$qr_ZorF_0}=$Qury_ZoF;	
        	          	  
        	          	  #$qr_ZorF_0  $sj_ZorF_5
        	            }          	
        	          }
        	          else {
        	            my $dieMsg="\n\nDIE!!!!\n\n$warnMsgHead\n\$qr_char=\$qr_char Should be . or - or a-z or A-Z!!!!\n\n\n";print $dieMsg; die $dieMsg;      	        
        	          }
        	          
        	          if    ( ($sj_char eq '-') || ($sj_char=~m/^[a-zA-Z]$/) ){
        	            if ($sj_char=~m/^[a-zA-Z]$/){
        	              $sj_pos_nb=$aln_sj_hdNub+$sj_pos_add_nb*$Subj_ZoF_one;
        	          	  $sj_pos_add_nb++;
        	            }
        	          	if (    ( ($qr_char eq '-') || ($qr_char=~m/^[a-zA-Z]$/) ) && (  ( defined ($sj_pos_nb) ) && ( $sj_pos_nb=~m/\d+/ )   )    ){      	          	
        	          	  $subject___Hash->{$sj_pos_nb}->{'1_mach_qr_char'}=$qr_char;
        	          	  $subject___Hash->{$sj_pos_nb}->{'0_self_sj_char'}=$sj_char;
        	          	  $subject___Hash->{$sj_pos_nb}->{'0_self_sj_ZorF'}=$Subj_ZoF;
        	          	  $subject___Hash->{$sj_pos_nb}->{'7_whoL_positio'}=$ex_whole_line_idx;
        	          	  $subject___Hash->{$sj_pos_nb}->{'8_exon__number'}=$ex_number;
        	          	       
        	          	  $wholeLine_Hash->{$ex_whole_line_idx}->{$sj_posi_5}=$sj_pos_nb;  
        	          	  $wholeLine_Hash->{$ex_whole_line_idx}->{$sj_ZorF_5}=$Subj_ZoF;	
        	          	  
        	          	  #$qr_ZorF_0  $sj_ZorF_5
        	            }
        	          }
        	          else {
        	            my $dieMsg="\n\nDIE!!!!\n\n$warnMsgHead\n\$sj_char=\$sj_char Should be  - or a-z or A-Z!!!!\n\n\n";print $dieMsg; die $dieMsg;      	        
        	          }      	          
        	          
        	          if    ( ($qr_char=~m/^[a-zA-Z]$/) && ($sj_char=~m/^[a-zA-Z]$/) ){
        	            $query_____Hash->{$qr_pos_nb}->{'2_mach_sj_posi'}=$sj_pos_nb;
        	            $subject___Hash->{$sj_pos_nb}->{'2_mach_qr_posi'}=$qr_pos_nb;
        	          }
        	          
        	          if ($CS_char eq '|'){
        	          	$query_____Hash->{$qr_pos_nb}->{'3_mach_identity'}=1;
        	          	$subject___Hash->{$sj_pos_nb}->{'3_mach_identity'}=1;
        	          }
        	          
        	          if ($AA_char=~/\S/){
        	          	$query_____Hash->{$qr_pos_nb}->{'4_self_transAAm'}=$AA_char;   
        	          	$query_____Hash->{$qr_pos_nb}->{'5_self_trsAAfrm'}=2;
        	          }
        	          
        	          if ( ($alinIdx>0) && ($qr_pos_add_nb==1) ){
        	            if ( $qr_pos_nb != ($last_qr_pos+$Qury_ZoF_one) ){
        	              my $dieMsg="\n\nDIE!!!!\n\n$warnMsgHead\n\$qr_pos_nb=\$qr_pos_nb Should be equal to $last_qr_pos+$Qury_ZoF_one=\$last_qr_pos+\$Qury_ZoF_one!!!!\n\n\n";print $dieMsg; die $dieMsg;  
        	            }
        	          }
        	          
        	          if ( ($alinIdx>0) && ($sj_pos_add_nb==1) ){
        	            if ( $sj_pos_nb != ($last_sj_pos+$Subj_ZoF_one) ){
        	              my $dieMsg="\n\nDIE!!!!\n\n$warnMsgHead\n\$sj_pos_nb=\$sj_pos_nb Should be equal to $last_sj_pos+$Subj_ZoF_one=\$last_sj_pos+\$Subj_ZoF_one!!!!\n\n\n";print $dieMsg; die $dieMsg;  
        	            }
        	          }
        	          
        	          
        	          
        	          $ex_whole_line_idx++;
        	        }
        	        
        	        $last_qr_pos=$qr_pos_nb;
        	        $last_sj_pos=$sj_pos_nb;
        	        
        	        
        	        $returnedArray->[$i]->{'7_whoL_Hsh'}=$wholeLine_Hash;
        	        $returnedArray->[$i]->{'8_qurM_Hsh'}=$query_____Hash;
        	        $returnedArray->[$i]->{'9_sbjM_Hsh'}=$subject___Hash;
        	        
        	        
        	
    		    	  }
    		    	  
    		    	  $alinIdx++;
    		    	  
    		    	}
    		    	
    		    	
    		    	
    		    	
    		    }
    		    else {
    		    	my $dieMsg="\n\nDIE!!!!\n\n$warnMsgHead\n\$eachExonRst=$eachExonRst\n\nThe \$eachExonRst is not in correct format:!!!!\n\n\n";print $dieMsg; die $dieMsg;  
    		    }
    		  }
    		  
    		  ##
    		  #  $returnedArray->[$i]->{'7_whoL_Hsh'}=$wholeLine_Hash;
    		  if (  defined ( $returnedArray->[$i]->{'7_whoL_Hsh'} )  ){
    		    my $whoL_Hsh=$returnedArray->[$i]->{'7_whoL_Hsh'};
    		    if ( ref ($whoL_Hsh) eq 'HASH' ){
    		      foreach my $woLpos (    sort {$a <=> $b} (   keys (  %{ $whoL_Hsh }  )   )    ){
    		        $returnedArray->[$i]->{'6_0_alnShw'}->{$spl_AA_4}.=$whoL_Hsh->{$woLpos}->{$AA_char_1};
    		        $returnedArray->[$i]->{'6_0_alnShw'}->{$splCtg_3}.=$whoL_Hsh->{$woLpos}->{$qr_char_2};
    		        $returnedArray->[$i]->{'6_0_alnShw'}->{$spl_CS_2}.=$whoL_Hsh->{$woLpos}->{$CS_char_3};
    		        $returnedArray->[$i]->{'6_0_alnShw'}->{$splGno_0}.=$whoL_Hsh->{$woLpos}->{$sj_char_4};
    		        my $exNB=$whoL_Hsh->{$woLpos}->{$Ex_numb_6}; $exNB=~s/^.*(\S)$/$1/;
    		        $returnedArray->[$i]->{'6_0_alnShw'}->{$splExN_1}.=$exNB;
    		      }
    		      
    		      if (  defined ( $returnedArray->[$i]->{'8_qurM_Hsh'} )  ){
    		        my $qurM_Hsh=$returnedArray->[$i]->{'8_qurM_Hsh'};
    		        if ( ref ($qurM_Hsh) eq 'HASH' ){
    		        	my @qurM_hash_key_array=(   keys (  %{ $qurM_Hsh }  )   );
    		        	my $qrCDSlength=@qurM_hash_key_array;
    		        	$returnedArray->[$i]->{'3_Qury_0_0mLt'}=$qrCDSlength;
    		        	
    		        	my @Sorted_qurM_hash_key_array=sort { $qurM_Hsh->{$a}->{'7_whoL_positio'} <=> $qurM_Hsh->{$b}->{'7_whoL_positio'} } @qurM_hash_key_array;
     		          
    		        	$returnedArray->[$i]->{'3_Qury_2_0Stt'}=$Sorted_qurM_hash_key_array[0];
    		        	$returnedArray->[$i]->{'3_Qury_2_1end'}=$Sorted_qurM_hash_key_array[$qrCDSlength-1];
    		        	
    		        	for (my $qrPsIdx=0; $qrPsIdx<$qrCDSlength; $qrPsIdx++){
     		            my $qurPos=$Sorted_qurM_hash_key_array[$qrPsIdx];
    		            if (   (  defined ( $qurM_Hsh->{$qurPos}->{'3_mach_identity'} )  ) && ( $qurM_Hsh->{$qurPos}->{'3_mach_identity'}==1 )   ){
    		              $returnedArray->[$i]->{'3_Qury_1_0idL'}++;
    		            }
                    
                    
                    if (   (  defined ( $qurM_Hsh->{$qurPos}->{'5_self_trsAAfrm'} )  ) && ( $qurM_Hsh->{$qurPos}->{'5_self_trsAAfrm'}==2 )   ){
                      my $transAAm=$qurM_Hsh->{$qurPos}->{'4_self_transAAm'};  #print "\$transAAm=\$transAAm=\$qurM_Hsh->{$qurPos}->{'4_self_transAAm'}=$qurM_Hsh->{$qurPos}->{'4_self_transAAm'}\n";
                      my $threeCodon=$qurM_Hsh->{$qurPos}->{'0_self_qr_char'};
                      my @posArray=($qurPos);
                      if ( $qrPsIdx > 0 ){
                      	my $prePos=$Sorted_qurM_hash_key_array[$qrPsIdx-1];  push @posArray, $prePos;
                      	my $preChr=$qurM_Hsh->{$prePos}->{'0_self_qr_char'};
                      	$returnedArray->[$i]->{'8_qurM_Hsh'}->{$prePos}->{'5_self_trsAAfrm'}=1;
                      	$threeCodon=$preChr.$threeCodon;
                      }
                      if ( $qrPsIdx < ($qrCDSlength-1) ){
                      	my $aftPos=$Sorted_qurM_hash_key_array[$qrPsIdx+1];  push @posArray, $aftPos;
                      	my $aftChr=$qurM_Hsh->{$aftPos}->{'0_self_qr_char'};
                      	$returnedArray->[$i]->{'8_qurM_Hsh'}->{$aftPos}->{'5_self_trsAAfrm'}=3;
                      	$threeCodon=$threeCodon.$aftChr;
                      }            
                      
                      
                      
                      foreach my $psHere ( @posArray ){
                        $returnedArray->[$i]->{'8_qurM_Hsh'}->{$psHere}->{'4_self_transAAm'}=$transAAm;
                        $returnedArray->[$i]->{'8_qurM_Hsh'}->{$psHere}->{'6_self_threeCod'}=$threeCodon;
                        if(  ( length ($threeCodon) ) == 3  ){
                          $returnedArray->[$i]->{'8_qurM_Hsh'}->{$psHere}->{'7_self_fredJtrs'}=FastaFileHandle::CodonTable_TGA_U($threeCodon);
                          
                        }
                      }
                    
                    }
                    
                    
    
        	        
    		            
    		          }
    		        }
    		      }
    		      
    		      if (  defined ( $returnedArray->[$i]->{'9_sbjM_Hsh'} )  ){
    		        my $sbjM_Hsh=$returnedArray->[$i]->{'9_sbjM_Hsh'};
    		        if ( ref ($sbjM_Hsh) eq 'HASH' ){
    		        	my @sbjM_hash_key_array=(   keys (  %{ $sbjM_Hsh }  )   );
    		        	my $sjCDSlength=@sbjM_hash_key_array;
    		        	my @Sorted_sbjM_hash_key_array=sort { $sbjM_Hsh->{$a}->{'7_whoL_positio'} <=> $sbjM_Hsh->{$b}->{'7_whoL_positio'} } @sbjM_hash_key_array;
     		          
    		        	$returnedArray->[$i]->{'5_Subj_0_0mLt'}=$sjCDSlength;
    		        	$returnedArray->[$i]->{'3_Subj_2_0Stt'}=$Sorted_sbjM_hash_key_array[0];
    		        	$returnedArray->[$i]->{'3_Subj_2_1end'}=$Sorted_sbjM_hash_key_array[$sjCDSlength-1];
    		        	
    		        	foreach my $sbjPos ( @Sorted_sbjM_hash_key_array ){
    		            if (   (  defined ( $sbjM_Hsh->{$sbjPos}->{'3_mach_identity'} )  ) && ( $sbjM_Hsh->{$sbjPos}->{'3_mach_identity'}==1 )   ){
    		              $returnedArray->[$i]->{'5_Subj_1_0idL'}++;
    		            }
    		          }
    		        }
    		      }
    		      
    		      
    		    }
    		  }
    		  
    		  #  $returnedArray->[$i]->{'8_qurM_Hsh'}=$query_____Hash;
        	
    		  
    		  
        	#  $returnedArray->[$i]->{'9_sbjM_Hsh'}=$subject___Hash;
        	##        
    		  
    		  
    		  
    		}
    		else {
    			my $dieMsg="\n\nDIE!!!!\n\n$warnMsgHead\n\$eachGeneRst=$eachGeneRst\n\nThe \$eachGeneRst is not in correct format:!!!!\n\n\n";print $dieMsg; die $dieMsg; 
    		}
    		$i++;   	
    	} 
    
     
    } 
    else {
      my $dieMsg="\n\nDIE!!!!\n\n$warnMsgHead\n\$inSplignOutFile=$inSplignOutFile\n\$wholeInSplignOut=$wholeInSplignOut\n\nThe \$inSplignOutFile is not in correct format:!!!!\n\n\n";print $dieMsg; die $dieMsg; 
    }
  return $returnedArray;
  }
 
 
  
  
}

sub AddQueryLength_get_covIdt{
	my $warnMsgHead="\n\n\n   In package  SplignHandle,\tIn sub AddQueryLength_get_covIdt,\n\n";	
	my ($inArray, $queryLength)=@_;
	if ( ref ($inArray) eq 'ARRAY' ){
		for (  my $i=0; $i<@{ $inArray }; $i++  ){
			
			$inArray->[$i]->{'3_Qury_0_1mlR'}=$inArray->[$i]->{'3_Qury_0_0mLt'}/$queryLength;
			$inArray->[$i]->{'3_Qury_0_2mlP'}=BlastHandle::ChangeTo100PercentNB ( $inArray->[$i]->{'3_Qury_0_1mlR'} );
			
			$inArray->[$i]->{'3_Qury_1_1idR'}=$inArray->[$i]->{'3_Qury_1_0idL'}/$queryLength;
			$inArray->[$i]->{'3_Qury_1_2idP'}=BlastHandle::ChangeTo100PercentNB ( $inArray->[$i]->{'3_Qury_1_1idR'} );
			
			$inArray->[$i]->{'3_Qury_Idt_rt'}=$inArray->[$i]->{'3_Qury_1_0idL'}/$inArray->[$i]->{'3_Qury_0_0mLt'};
			$inArray->[$i]->{'3_Qury_Idt_pt'}=BlastHandle::ChangeTo100PercentNB ( $inArray->[$i]->{'3_Qury_Idt_rt'} );
			
		}
	}
	my $outPutArray;
	
	@{ $outPutArray }=sort { $b->{'3_Qury_1_1idR'} <=> $a->{'3_Qury_1_1idR'} } @{ $inArray };
	
	return $outPutArray;
	
}

sub Add_KeyPos_and_Sort_By_covIDt{ #use core position such as TGA-sec postiont to sort splign result
  my $warnMsgHead="\n\n\n   In package  SplignHandle,\tIn sub Add_KeyPos_and_Sort_By_covIDt,\n\n";	
	my ($inArray, $sbj_key_pos, $InEstZoF, $InGnoZoF)=@_;
	if ( ref ($inArray) eq 'ARRAY' ){
		for (  my $i=0; $i<@{ $inArray }; $i++  ){
			
			if   (  defined  ( $inArray->[$i]->{'9_sbjM_Hsh'}->{$sbj_key_pos} )  ){
			  $inArray->[$i]->{'3_Subj_0kyp'}=1;
			}
			else {$inArray->[$i]->{'3_Subj_0kyp'}=0;}
			
			$inArray->[$i]->{'3_Qury_ZoF_orgPredict'}=$InEstZoF;
			$inArray->[$i]->{'5_Subj_ZoF_orgPredict'}=$InGnoZoF;
			my $directionScore=0; 
			if    (  ( $inArray->[$i]->{'3_Qury_ZoF'} eq $InEstZoF ) && ( $inArray->[$i]->{'5_Subj_ZoF'} eq $InEstZoF )  ){ $directionScore=3; }
			elsif (  ( $inArray->[$i]->{'3_Qury_ZoF'} eq $InEstZoF ) && ( $inArray->[$i]->{'5_Subj_ZoF'} ne $InEstZoF )  ){ $directionScore=2; }
			elsif (  ( $inArray->[$i]->{'3_Qury_ZoF'} ne $InEstZoF ) && ( $inArray->[$i]->{'5_Subj_ZoF'} eq $InEstZoF )  ){ $directionScore=1; }
			$inArray->[$i]->{'5_z_QrySbj_ZoF_IdScor'}=$directionScore;
		}
	}
	
	my $outPutArray;
	
	@{ $outPutArray }=sort {  ( $b->{'3_Subj_0kyp'} <=> $a->{'3_Subj_0kyp'} ) or ( $b->{'5_z_QrySbj_ZoF_IdScor'} <=> $a->{'5_z_QrySbj_ZoF_IdScor'} ) or( $b->{'3_Qury_1_1idR'} <=> $a->{'3_Qury_1_1idR'} )  } @{ $inArray };
	
	return $outPutArray;
}

sub Add_KeyPosList_and_Sort_By_covIDt{ #use represented CDS position to sort splign result
  my $warnMsgHead="\n\n\n   In package  SplignHandle,\tIn sub Add_KeyPosList_and_Sort_By_covIDt,\n\n";	
	my ($inArray, $sbj_key_pos_list, $InEstZoF, $InGnoZoF)=@_;

	if ( ref ($inArray) eq 'ARRAY' ){
		for (  my $i=0; $i<@{ $inArray }; $i++  ){
			$inArray->[$i]->{'3_Subj_0kyp'}=0;
			if (   (  defined ( $sbj_key_pos_list )  ) && (  ref ( $sbj_key_pos_list ) eq 'ARRAY'  )   ){
				my $AllKeyPosNB=0; my $realMatchNB=0;
				for (  my $j=0; $j<@{ $sbj_key_pos_list }; $j++  ){
					my $keyPos=$sbj_key_pos_list->[$j];
					if (   (  defined ( $keyPos )  ) && ( $keyPos=~m/\d+/ )   ){
						
					  if   (  defined  ( $inArray->[$i]->{'9_sbjM_Hsh'}->{$keyPos} )  ){
			      
			        $realMatchNB++;
			      }
			      
			      $AllKeyPosNB++;
			    }
				}
				$inArray->[$i]->{'3_Subj_0kyp'}=$realMatchNB/$AllKeyPosNB if ($AllKeyPosNB>0);
			  
			}
			
			$inArray->[$i]->{'3_Qury_ZoF_orgPredict'}=$InEstZoF;
			$inArray->[$i]->{'5_Subj_ZoF_orgPredict'}=$InGnoZoF;
			my $directionScore=0; 
			if    (  ( $inArray->[$i]->{'3_Qury_ZoF'} eq $InEstZoF ) && ( $inArray->[$i]->{'5_Subj_ZoF'} eq $InEstZoF )  ){ $directionScore=3; }
			elsif (  ( $inArray->[$i]->{'3_Qury_ZoF'} eq $InEstZoF ) && ( $inArray->[$i]->{'5_Subj_ZoF'} ne $InEstZoF )  ){ $directionScore=2; }
			elsif (  ( $inArray->[$i]->{'3_Qury_ZoF'} ne $InEstZoF ) && ( $inArray->[$i]->{'5_Subj_ZoF'} eq $InEstZoF )  ){ $directionScore=1; }
			$inArray->[$i]->{'5_z_QrySbj_ZoF_IdScor'}=$directionScore;
		}
	}
	
	my $outPutArray;
	
	@{ $outPutArray }=sort {  ( $b->{'3_Subj_0kyp'} <=> $a->{'3_Subj_0kyp'} ) or ( $b->{'5_z_QrySbj_ZoF_IdScor'} <=> $a->{'5_z_QrySbj_ZoF_IdScor'} ) or( $b->{'3_Qury_1_1idR'} <=> $a->{'3_Qury_1_1idR'} )  } @{ $inArray };
	
	return $outPutArray;
}

sub simplify_splignPraser_Out_Hash{
	my ($inArray)=@_;                          #warn ref( $inArray ) , "  111 \n";
	my $outArray=Storable::dclone ($inArray);  #warn ref( $outArray ) ,   "  222 \n";
	
	if (  ref( $outArray ) eq 'ARRAY'  ){
		for (  my $i=0; $i<@{ $outArray }; $i++  ){   #warn "\$i=$i\n";
			
			
			if (  ref( $outArray->[$i]->{'7_whoL_Hsh'} ) eq 'HASH'  ){  #warn "7_whoL_Hsh\n";
        delete ( $outArray->[$i]->{'7_whoL_Hsh'} );
      }
       
      if (  ref( $outArray->[$i]->{'8_qurM_Hsh'} ) eq 'HASH'  ){ #warn "8_qurM_Hsh\n";
        delete ( $outArray->[$i]->{'8_qurM_Hsh'} );
      }
       
      if (  ref( $outArray->[$i]->{'9_sbjM_Hsh'} ) eq 'HASH'  ){ #warn "9_sbjM_Hsh\n";
        delete ( $outArray->[$i]->{'9_sbjM_Hsh'} );
      }
       
			
		}
	}
	
	return $outArray;
}

sub build_EstMatchReportHash_with_only_Cap3Out{
  my ($cap3Hash, $bestProSplHash)=@_;
  
  my $warnMsgBody="\nIn package  SplignHandle,\tIn sub build_EstMatchReportHash_with_only_Cap3Out,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	
  my $splGno_0='1_0splGno';                            my $sj_char_4='4sj_char';               #'3_sj_line'                         
	my $splExN_1='1_1splExN';                            my $Ex_numb_6='6Ex_numb';               #'4_ex_numb'                         
	my $spl_CS_2='1_2spl_CS';                            my $CS_char_3='3CS_char';               #'2_cs_line'                         
	my $splCtg_3='1_3splCtg';                            my $qr_char_2='2qr_char';               #'1_qr_line'                         
	my $spl_AA_4='1_4spl_AA';                            my $AA_char_1='1AA_char';               #'0_AA_line'                         
	my $ctgPos_a='Not_include_spl_ctg_posi';             my $qr_posi_0='0qr_posi';               #'Not_include_0qr_posi'
	my $gnoPos_b='Not_include_spl_gno_posi';             my $sj_posi_5='5sj_posi';               #'Not_include_5sj_posi'   
	my $ctgZoF_a='Not_include_spl_ctg_ZorF';             my $qr_ZorF_0='0qr_ZorF';               #  
	my $gnoZoF_b='Not_include_spl_gno_ZorF';             my $sj_ZorF_5='5sj_ZorF';               #

	
	my $ProSplignHash_DNA_ZoF;
	if (   ( ref ($bestProSplHash) eq 'HASH' ) && (  defined ( $bestProSplHash->{'2_DNA_ZF'} )  ) && ( $bestProSplHash->{'2_DNA_ZF'}=~/(-|\+)/)   ){
	  $ProSplignHash_DNA_ZoF=$bestProSplHash->{'2_DNA_ZF'};	
	}
	else {
		my $dieMsgB=$die_MsgHead."\$bestProSplHash->{'2_DNA_ZF'}=$bestProSplHash->{'2_DNA_ZF'} should be defined and should be + or -  !!!!!\n\n\n";
	}
	
  my $reverse_on_or_off=0; # The Contig ZorF was set to + in the original version of this method, but in new version we need set it the same as the prosplign DNA ZorF. So that, if the $ProSplignHash_DNA_ZoF eq '-', we will reverse the ACTG in the sub code programs
  if ( $ProSplignHash_DNA_ZoF eq '-'){ $reverse_on_or_off=1; }
	
	
	
	my $bestCap3Hash=$cap3Hash;
	
   
	my $betCap3_Wps2DifPos=$bestCap3Hash->{'_Consensus'}->{'_Wpos2DifTypePos'};
	my $betCap3_Wps2Difwrd=$bestCap3Hash->{'_Consensus'}->{'_WlPos2Word'};
	
	my $betCap3_CtgPos2Wps=$bestCap3Hash->{'_Consensus'}->{'_DifTypePos2Wpos'}->{'_NactgType'}; 
  
	my $Cap3EstPrtKeyHash=&BuildCap3OrderKeyHash($bestCap3Hash);  #DirFileHandle::PrintAndWarnDumper ($Cap3EstPrtKeyHash);
	
	
	my $new_WhoLine_Hash; my $new_WhoLine_idex=0;
	
	###############################
	#$qr_posi_0 => 1,
  #                                   $AA_char_1 => ' ',
  #                                   $qr_char_2 => 'A',
	#warn "\$betCap3_Wps2Difwrd=$betCap3_Wps2Difwrd\n";
	if (  ref( $betCap3_Wps2Difwrd ) eq 'HASH'  ) {  
		foreach my $woLinPos (    sort {$a <=> $b} (   keys (  %{ $betCap3_Wps2Difwrd }  )   )    ){
			if ($betCap3_Wps2DifPos->{$woLinPos}->[0] eq '_NactgType'){
			  $new_WhoLine_Hash->{$woLinPos}->{$qr_posi_0}=$betCap3_Wps2DifPos->{$woLinPos}->[1]+1;
			  $new_WhoLine_Hash->{$woLinPos}->{$AA_char_1}=' ';
			  $new_WhoLine_Hash->{$woLinPos}->{$qr_char_2}=$betCap3_Wps2Difwrd->{$woLinPos};    if ($reverse_on_or_off==1) { $new_WhoLine_Hash->{$woLinPos}->{$qr_char_2}=FastaFileHandle::ReverseComplementBase( $new_WhoLine_Hash->{$woLinPos}->{$qr_char_2} ); }
			}
			else {
				$new_WhoLine_Hash->{$woLinPos}->{$qr_posi_0}='';
			  $new_WhoLine_Hash->{$woLinPos}->{$AA_char_1}=' ';
			  $new_WhoLine_Hash->{$woLinPos}->{$qr_char_2}=$betCap3_Wps2Difwrd->{$woLinPos};    if ($reverse_on_or_off==1) { $new_WhoLine_Hash->{$woLinPos}->{$qr_char_2}=FastaFileHandle::ReverseComplementBase( $new_WhoLine_Hash->{$woLinPos}->{$qr_char_2} ); }
			}
			
			
			if (   (  ref( $bestCap3Hash->{'_Align'} ) eq 'HASH'  ) && (  ref( $bestCap3Hash->{'1_OrderHash'} ) eq 'HASH'  )   ) {
		    foreach my $ordKey (    sort {$a <=> $b} (   keys (  %{ $bestCap3Hash->{'1_OrderHash'} }  )   )    ){
		      my ($EstID, $EstZoF)=@{ $bestCap3Hash->{'1_OrderHash'}->{$ordKey} };
		      
		      $new_WhoLine_Hash->{$woLinPos}->{ $Cap3EstPrtKeyHash->{$ordKey} }=$bestCap3Hash->{'_Align'}->{$EstID}->{$EstZoF}->{'_WlPos2Word'}->{$woLinPos};  if ($reverse_on_or_off==1) { $new_WhoLine_Hash->{$woLinPos}->{ $Cap3EstPrtKeyHash->{$ordKey} }=FastaFileHandle::ReverseComplementBase( $new_WhoLine_Hash->{$woLinPos}->{ $Cap3EstPrtKeyHash->{$ordKey} } ); }
		    }
	    }
		 
		  
		}
	}
		
	#DirFileHandle::PrintAndWarnDumper ($new_WhoLine_Hash);	
	

	
	my $outPutHash;
	if ( ref ($new_WhoLine_Hash) eq 'HASH' ){
		foreach my $newWLposKey (    sort {$a<=>$b} (   keys (  %{ $new_WhoLine_Hash }  )   )    ) {# $ctgZoF_a
			
						                                                                                                                                   $outPutHash->{'5_alnShwHs'}->{$newWLposKey}->{$ctgPos_a}=$new_WhoLine_Hash->{$newWLposKey}->{$qr_posi_0};
						                                                                                                                                   $outPutHash->{'5_alnShwHs'}->{$newWLposKey}->{$ctgZoF_a}='+' if (   (  defined ( $outPutHash->{'5_alnShwHs'}->{$newWLposKey}->{$ctgPos_a} )  ) && ( $outPutHash->{'5_alnShwHs'}->{$newWLposKey}->{$ctgPos_a}=~m/\d+/ )   );
			                                                                                                                                         $outPutHash->{'5_alnShwHs'}->{$newWLposKey}->{$ctgZoF_a}=&reverseZoF  ( $outPutHash->{'5_alnShwHs'}->{$newWLposKey}->{$ctgZoF_a} );
			
			$outPutHash->{'0___alnShw'}->{$spl_AA_4}.=$new_WhoLine_Hash->{$newWLposKey}->{$AA_char_1};                                               $outPutHash->{'5_alnShwHs'}->{$newWLposKey}->{$spl_AA_4}=$new_WhoLine_Hash->{$newWLposKey}->{$AA_char_1};
    	$outPutHash->{'0___alnShw'}->{$splCtg_3}.=$new_WhoLine_Hash->{$newWLposKey}->{$qr_char_2};                                               $outPutHash->{'5_alnShwHs'}->{$newWLposKey}->{$splCtg_3}=$new_WhoLine_Hash->{$newWLposKey}->{$qr_char_2};
    	if (   (  ref( $bestCap3Hash->{'_Align'} ) eq 'HASH'  ) && (  ref( $bestCap3Hash->{'1_OrderHash'} ) eq 'HASH'  )   ) {
		    foreach my $ordKey (    sort {$a <=> $b} (   keys (  %{ $bestCap3Hash->{'1_OrderHash'} }  )   )    ){
		    	$outPutHash->{'0___alnShw'}->{ $Cap3EstPrtKeyHash->{$ordKey} }.=$new_WhoLine_Hash->{$newWLposKey}->{ $Cap3EstPrtKeyHash->{$ordKey} };  $outPutHash->{'5_alnShwHs'}->{$newWLposKey}->{ $Cap3EstPrtKeyHash->{$ordKey} }.=$new_WhoLine_Hash->{$newWLposKey}->{ $Cap3EstPrtKeyHash->{$ordKey} };
		    }
	    }
		}
	}          
	
	my $EstNameHash;
	if (   (  ref( $bestCap3Hash->{'_Align'} ) eq 'HASH'  ) && (  ref( $bestCap3Hash->{'1_OrderHash'} ) eq 'HASH'  )   ) {
		foreach my $ordKey (    sort {$a <=> $b} (   keys (  %{ $bestCap3Hash->{'1_OrderHash'} }  )   )    ){
		  my ($EstID, $EstZoF)=@{ $bestCap3Hash->{'1_OrderHash'}->{$ordKey} };
		  $EstNameHash->{ $Cap3EstPrtKeyHash->{$ordKey} }=[$EstID, $EstZoF];
	  }
	}
	
	$outPutHash->{'2_EstIdHsh'}=$EstNameHash;
	$outPutHash->{'4_whoL_Hsh'}=$new_WhoLine_Hash;
	
	
	if ( ref ($new_WhoLine_Hash) eq 'HASH' ){
		my $outSegMetHash; 
		my $lastStateHash; my $frcIdx=0;
	  foreach my $whoLinPos (    sort {$a<=>$b} (   keys (  %{ $new_WhoLine_Hash }  )   )    ) {
	  	#my $genoChar=$bstSpln_subjG_Hash->{$genoPosK}->{'0_self_sj_char'};   if ($splnCtgQuryZF eq '-') { $genoChar=FastaFileHandle::ReverseComplementBase($genoChar); }
	  	my $CotgChar=$new_WhoLine_Hash->{$whoLinPos}->{$qr_char_2};   
	  	my $CotgPosi=$new_WhoLine_Hash->{$whoLinPos}->{$qr_posi_0};               print "\$CotgPosi=\$new_WhoLine_Hash->{$whoLinPos}=$CotgPosi\n\n";
	 
	  	
	  	
	  	if (    (  defined ( $bestCap3Hash )  ) && (  ref ( $bestCap3Hash ) eq 'HASH'  ) && (  defined ( $bestCap3Hash->{'1_OrderHash'} )  ) && (  ref ( $bestCap3Hash->{'1_OrderHash'} ) eq 'HASH'  )    ){
	  	  foreach my $ordKey (    sort {$a <=> $b} (   keys (  %{ $bestCap3Hash->{'1_OrderHash'} }  )   )    ){
	  	  	my ($EstID, $EstZoF)=@{ $bestCap3Hash->{'1_OrderHash'}->{$ordKey} };
	  	  	my $EstIdxID=$Cap3EstPrtKeyHash->{$ordKey};     
	  	  	
	  	  	my $Cap3ResposWhoPos; $Cap3ResposWhoPos=$betCap3_CtgPos2Wps->{$CotgPosi-1} if (   ( defined ($CotgPosi) ) &&  ($CotgPosi=~m/\d+/) && (  defined ( $betCap3_CtgPos2Wps->{$CotgPosi-1} )  ) && ( $betCap3_CtgPos2Wps->{$CotgPosi-1}=~m/\d+/ )    );
	  	  	my $Est_Char;  $Est_Char=$bestCap3Hash->{'_Align'}->{$EstID}->{$EstZoF}->{'_WlPos2Word'}->{$Cap3ResposWhoPos}           if (      ( defined ($Cap3ResposWhoPos) ) &&  ($Cap3ResposWhoPos=~m/\d+/) && (  defined ( $bestCap3Hash->{'_Align'}->{$EstID}->{$EstZoF}->{'_WlPos2Word'}->{$Cap3ResposWhoPos}           )  )   );
	  	  	if ($reverse_on_or_off==1) { $Est_Char=FastaFileHandle::ReverseComplementBase( $Est_Char ); }
	  	  	my $Est_posi;  $Est_posi=$bestCap3Hash->{'_Align'}->{$EstID}->{$EstZoF}->{'_Wpos2DifTypePos'}->{$Cap3ResposWhoPos}->[1] if (      ( defined ($Cap3ResposWhoPos) ) &&  ($Cap3ResposWhoPos=~m/\d+/) && (  defined ( $bestCap3Hash->{'_Align'}->{$EstID}->{$EstZoF}->{'_Wpos2DifTypePos'}->{$Cap3ResposWhoPos}->[1] )  )   );
	  	  	
	  	  	my $EstStMet=&EstVSCtgSttMentDefine($CotgChar, $Est_Char);
	  	  	 
	  	  	###############################
	  	  	
	  	  	my $EstNewSegSt=0;
	  	    if ($EstStMet ne '5_OTHER_TYP' ){
	  	    	if    ($frcIdx==0)                                                                                                                                                                 { $EstNewSegSt=1; }
	  	    	elsif (   (  defined ( $lastStateHash->{$EstIdxID}->{'stageMet'} )  ) && ($lastStateHash->{$EstIdxID}->{'stageMet'} ne $EstStMet)   )                                              { $EstNewSegSt=1; }
	  	    	elsif (    (   (  abs ( $CotgPosi-$lastStateHash->{$EstIdxID}->{'ctg_PosK'} )  ) > 1   ) ||  (   (  abs ( $Est_posi-$lastStateHash->{$EstIdxID}->{'Est_posi'} )  ) > 1    )     )  { $EstNewSegSt=1; }
	  	    	else     {	  			$EstNewSegSt=0;	  		}
	  	    }
	  	    if ($EstNewSegSt){  
	  	      $lastStateHash->{$EstIdxID}->{'startPos'}=$CotgPosi;
	  	    	$outSegMetHash->{$EstIdxID}->{ $lastStateHash->{$EstIdxID}->{'startPos'} }->{$EstStMet}=$CotgPosi;  
	  	    }
	  	    		
	  	    
	  	    my $EstCotinues=0;	  	
	  	    if ( 
	  	         ( $EstStMet ne '5_OTHER_TYP' ) 
	  	         &&
	  	         (   (  defined ( $lastStateHash->{$EstIdxID}->{'stageMet'} )  ) && ($lastStateHash->{$EstIdxID}->{'stageMet'} eq $EstStMet)   )   
	  	         &&
	  	         (    (   (  abs ( $CotgPosi-$lastStateHash->{$EstIdxID}->{'ctg_PosK'} )  ) == 1    ) && (   (  abs ( $Est_posi-$lastStateHash->{$EstIdxID}->{'Est_posi'} )  ) == 1    )     )
	  	    )  	{ $EstCotinues=1; }	  		
	  	    if ($EstCotinues){  
	  	    	$outSegMetHash->{$EstIdxID}->{ $lastStateHash->{$EstIdxID}->{'startPos'} }->{$EstStMet}=$CotgPosi;  
	  	    }    
	  	    
	  	    $lastStateHash->{$EstIdxID}->{'stageMet'}=$EstStMet;
	  	    $lastStateHash->{$EstIdxID}->{'ctg_PosK'}=$CotgPosi;
	  	    $lastStateHash->{$EstIdxID}->{'Est_posi'}=$Est_posi;
	  	  	
	  	  	###############################
	  	  	
	  	  }
	  	}
	  	$frcIdx++;
	  }
	  
	  
	  if ( ref ($outSegMetHash) eq 'HASH'){
	  	my $EstMatchGeneStructure_Hash;
	  	foreach my $keyHere (    sort {$a cmp $b} (   keys (  %{ $outSegMetHash }  )   )    ){
	  	  if (  ref ( $outSegMetHash->{$keyHere} ) eq 'HASH'  ){                                                               
          foreach my $sgmHead (    sort {$a <=> $b} (   keys (  %{ $outSegMetHash->{$keyHere} }  )   )    ){  
          	my ($sgmType)=(   keys (  %{ $outSegMetHash->{$keyHere}->{$sgmHead} }  )   );
          	my $sgmTail=$outSegMetHash->{$keyHere}->{$sgmHead}->{$sgmType};
          	my $tempHash;
          	$tempHash->{'0_SgmtHead'}=$sgmHead;
          	$tempHash->{'1_SgmtTail'}=$sgmTail;
          	$tempHash->{'3_SgmtType'}=$sgmType;
          	
          	push @{ $EstMatchGeneStructure_Hash->{$keyHere} }, $tempHash;
          }
        }          
	  	}
	  	$outPutHash->{'3_GStr_Hsh'}=$EstMatchGeneStructure_Hash;
	  }	  
	  
	}
	
	if ($reverse_on_or_off==1){
		#0___alnShw  4_whoL_Hsh  5_alnShwHs
		$outPutHash->{'0___alnShw'} =&ReversePtShwHash( $outPutHash->{'0___alnShw'} );
		$outPutHash->{'4_whoL_Hsh'} =&ReverseNBkeyHash( $outPutHash->{'4_whoL_Hsh'} );
		$outPutHash->{'5_alnShwHs'} =&ReverseNBkeyHash( $outPutHash->{'5_alnShwHs'} );
		
	}
	
	return $outPutHash;
	
}

sub ReversePtShwHash{
	my ($inHash)=@_;
	
	my $warnMsgBody="\nIn package  SplignHandle,\tIn sub ReversePtShwHash,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $outHash;
	if ( ref ($inHash) eq 'HASH' ){
		foreach my $ecPtKey (    sort { $a cmp $b } (   keys (  %{ $inHash }  )   )    ){
			$outHash->{$ecPtKey} = reverse $inHash->{$ecPtKey};
		}
		return $outHash;
	}
	else {
		my $dieMsg=$die_MsgHead."\$inHash=$inHash should be a HASH ref!!!!\n\n\n\n";
	} 
}

sub ReverseNBkeyHash{
	my ($inHash)=@_;
	
	my $warnMsgBody="\nIn package  SplignHandle,\tIn sub ReverseNBkeyHash,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $outHash; my $outHashIdxKey=0;
	if ( ref ($inHash) eq 'HASH' ){
		foreach my $eachKey (    sort { $b <=> $a } (   keys (  %{ $inHash }  )   )    ){
			$outHash->{$outHashIdxKey} = $inHash->{$eachKey};
			$outHashIdxKey++;
		}
		return $outHash;
	}
	else {
		my $dieMsg=$die_MsgHead."\$inHash=$inHash should be a HASH ref!!!!\n\n\n\n";
	} 
}

sub build_EstGenoMatchReportHash_withSpilgnANDCap3Out{
  my ($splignArray, $cap3Hash, $bestProSplHash)=@_;
  
  my $warnMsgBody="\nIn package  SplignHandle,\tIn sub buildSplignOutHash_withCpp3OutHash,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	
	
  my $splGno_0='1_0splGno';                            my $sj_char_4='4sj_char';               #'3_sj_line'                         
	my $splExN_1='1_1splExN';                            my $Ex_numb_6='6Ex_numb';               #'4_ex_numb'                         
	my $spl_CS_2='1_2spl_CS';                            my $CS_char_3='3CS_char';               #'2_cs_line'                         
	my $splCtg_3='1_3splCtg';                            my $qr_char_2='2qr_char';               #'1_qr_line'                            
	my $spl_AA_4='1_4spl_AA';                            my $AA_char_1='1AA_char';               #'0_AA_line'                         
	my $ctgPos_a='Not_include_spl_ctg_posi';             my $qr_posi_0='0qr_posi';               #'Not_include_0qr_posi'
	my $gnoPos_b='Not_include_spl_gno_posi';             my $sj_posi_5='5sj_posi';               #'Not_include_5sj_posi'   
	my $ctgZoF_a='Not_include_spl_ctg_ZorF';             my $qr_ZorF_0='0qr_ZorF';               #  
	my $gnoZoF_b='Not_include_spl_gno_ZorF';             my $sj_ZorF_5='5sj_ZorF';               #
	
	my $cap3Wl_c='Not_include_cap3_wl_posi';
	
	my $ProSplignHash_DNA_ZoF;
	if (   ( ref ($bestProSplHash) eq 'HASH' ) && (  defined ( $bestProSplHash->{'2_DNA_ZF'} )  ) && ( $bestProSplHash->{'2_DNA_ZF'}=~/(-|\+)/)   ){
	  $ProSplignHash_DNA_ZoF=$bestProSplHash->{'2_DNA_ZF'};	
	}
	else {
		my $dieMsgB=$die_MsgHead."\$bestProSplHash->{'2_DNA_ZF'}=$bestProSplHash->{'2_DNA_ZF'} should be defined and should be + or -  !!!!!\n\n\n";
	}
	
  my $reverse_on_or_off=0; # The Contig ZorF was set to + in the original version of this method, but in new version we need set it the same as the prosplign DNA ZorF. So that, if the $ProSplignHash_DNA_ZoF eq '-', we will reverse the ACTG in the sub code programs
  if ( $ProSplignHash_DNA_ZoF eq '-'){ $reverse_on_or_off=1; }
	
	my $bestSplignHash;
	if (   ( ref ($splignArray) eq 'ARRAY' ) && (  ref ( $splignArray->[0] ) eq 'HASH'  )   ){ $bestSplignHash=$splignArray->[0];	}
	else {my $splignArrayDump=DirFileHandle::ReturnDumperInform($splignArray);my $dieMsg="$die_MsgHead\nThe \$splignArray=$splignArray should be a correct array ref\n\$splignArrayDump=$splignArrayDump\n\n\n"; print $dieMsg; die $dieMsg;} 
	
	#my $bestCap3Hash;
	#if (   ( ref ($cap3Hash) eq 'HASH' ) && (  ref ( $cap3Hash->{1} ) eq 'HASH'  )   ){ $bestCap3Hash=$cap3Hash->{1};	}  #注意这里的->{1}，不确定是不是正确的
	#else {my $cap3HashDump=DirFileHandle::ReturnDumperInform($cap3Hash);my $dieMsg="$die_MsgHead\nThe \$cap3Hash=$cap3Hash should be a correct array ref\n\$cap3HashDump=$cap3HashDump\n\n\n"; print $dieMsg; die $dieMsg;} 
	
	my $bestCap3Hash=$cap3Hash;
	
  #my $CtgLength=$bestCap3Hash->{'_Consensus'}->{'_NactgTypeLengthKw'};
	my $splnCtgQuryZF=$bestSplignHash->{'3_Qury_ZoF'};
	my $splnGnoSubjZF=$bestSplignHash->{'5_Subj_ZoF'};                  #my $realCtgZorF=$splnCtgQuryZF; if ($splnGnoSubjZF eq '-'){$realCtgZorF=&reverseZoF($realCtgZorF);}
  
  my $realSplnGnoZF;
  if    (   (  ( $splnCtgQuryZF eq '+' ) && ($splnCtgQuryZF eq '-')  ) || (  ( $splnCtgQuryZF eq '-' ) && ($splnCtgQuryZF eq '+')  )   ) { $realSplnGnoZF='-'; }
  elsif (   (  ( $splnCtgQuryZF eq '+' ) && ($splnCtgQuryZF eq '+')  ) || (  ( $splnCtgQuryZF eq '-' ) && ($splnCtgQuryZF eq '-')  )   ) { $realSplnGnoZF='+'; }
  
  #my $realSplnGnoZF=$splnGnoSubjZF;   if ($splnCtgQuryZF eq '-'){$realSplnGnoZF=&reverseZoF($realSplnGnoZF);}
  
	my $bstSpln_whoL__Hash=$bestSplignHash->{'7_whoL_Hsh'};
	my $betCap3_Wps2DifPos=$bestCap3Hash->{'_Consensus'}->{'_Wpos2DifTypePos'};
	my $betCap3_Wps2Difwrd=$bestCap3Hash->{'_Consensus'}->{'_WlPos2Word'};
	
	my $bstSpln_qury__Hash=$bestSplignHash->{'8_qurM_Hsh'};
	my $betCap3_CtgPos2Wps=$bestCap3Hash->{'_Consensus'}->{'_DifTypePos2Wpos'}->{'_NactgType'}; 
  
  my $bstSpln_subjG_Hash=$bestSplignHash->{'9_sbjM_Hsh'};
	
  
	my $Cap3EstPrtKeyHash=&BuildCap3OrderKeyHash($bestCap3Hash);
	
	
	my $new_WhoLine_Hash; my $new_WhoLine_idex=0;
	if ( ref ($bstSpln_whoL__Hash) eq 'HASH' ){
	  my @bstSpln_whoLhsKey=keys (  %{ $bstSpln_whoL__Hash }  ); 
	  my @sort_bstSpln_whoLhsKey=sort { $a <=> $b } @bstSpln_whoLhsKey;
	  if ($splnCtgQuryZF eq '-'){@sort_bstSpln_whoLhsKey=sort { $b <=> $a } @bstSpln_whoLhsKey;}
	  foreach my $bstSpln_whoL_key (@sort_bstSpln_whoLhsKey){
	  	
	  	
	  	
	  	
	  	my $SplgnQryCtgPos=$bstSpln_whoL__Hash->{$bstSpln_whoL_key}->{$qr_posi_0};
	  	#
	  	my $SplgnQryCtgChr=$bstSpln_whoL__Hash->{$bstSpln_whoL_key}->{$qr_char_2}; if ($splnCtgQuryZF eq '-') { $SplgnQryCtgChr=FastaFileHandle::ReverseComplementBase($SplgnQryCtgChr); }  if ($reverse_on_or_off==1) { $SplgnQryCtgChr=FastaFileHandle::ReverseComplementBase($SplgnQryCtgChr); }
	    my $SplgnSbjGnoChr=$bstSpln_whoL__Hash->{$bstSpln_whoL_key}->{$sj_char_4}; if ($splnCtgQuryZF eq '-') { $SplgnSbjGnoChr=FastaFileHandle::ReverseComplementBase($SplgnSbjGnoChr); }  if ($reverse_on_or_off==1) { $SplgnSbjGnoChr=FastaFileHandle::ReverseComplementBase($SplgnSbjGnoChr); }
	    
	    $new_WhoLine_Hash->{$new_WhoLine_idex}->{$qr_posi_0}=$SplgnQryCtgPos;
      $new_WhoLine_Hash->{$new_WhoLine_idex}->{$qr_ZorF_0}=$bstSpln_whoL__Hash->{$bstSpln_whoL_key}->{$qr_ZorF_0};     if ($reverse_on_or_off==1) { $new_WhoLine_Hash->{$new_WhoLine_idex}->{$qr_ZorF_0}=&reverseZoF( $new_WhoLine_Hash->{$new_WhoLine_idex}->{$qr_ZorF_0} ); }
      $new_WhoLine_Hash->{$new_WhoLine_idex}->{$AA_char_1}=$bstSpln_whoL__Hash->{$bstSpln_whoL_key}->{$AA_char_1};
      $new_WhoLine_Hash->{$new_WhoLine_idex}->{$qr_char_2}=$SplgnQryCtgChr;
      $new_WhoLine_Hash->{$new_WhoLine_idex}->{$CS_char_3}=$bstSpln_whoL__Hash->{$bstSpln_whoL_key}->{$CS_char_3};
      $new_WhoLine_Hash->{$new_WhoLine_idex}->{$sj_char_4}=$SplgnSbjGnoChr;
      $new_WhoLine_Hash->{$new_WhoLine_idex}->{$sj_posi_5}=$bstSpln_whoL__Hash->{$bstSpln_whoL_key}->{$sj_posi_5}   if ( $SplgnSbjGnoChr=~/\w/); 
      $new_WhoLine_Hash->{$new_WhoLine_idex}->{$sj_ZorF_5}=$bstSpln_whoL__Hash->{$bstSpln_whoL_key}->{$sj_ZorF_5};    if ($reverse_on_or_off==1) { $new_WhoLine_Hash->{$new_WhoLine_idex}->{$sj_ZorF_5}=&reverseZoF( $new_WhoLine_Hash->{$new_WhoLine_idex}->{$sj_ZorF_5} ); }
      $new_WhoLine_Hash->{$new_WhoLine_idex}->{$Ex_numb_6}=$bstSpln_whoL__Hash->{$bstSpln_whoL_key}->{$Ex_numb_6}; 
      #warn "\$SplgnQryCtgPos=$SplgnQryCtgPos\n";
      if ( (defined $SplgnQryCtgPos) && (defined $SplgnQryCtgChr) && ($SplgnQryCtgPos=~m/^\d+$/) && ($SplgnQryCtgChr=~m/^[a-zA-Z]$/) ) {}
      else {
        if (   (  ref( $bestCap3Hash->{'_Align'} ) eq 'HASH'  ) && (  ref( $bestCap3Hash->{'1_OrderHash'} ) eq 'HASH'  )   ) {
		  	  foreach my $ordKey (    sort {$a <=> $b} (   keys (  %{ $bestCap3Hash->{'1_OrderHash'} }  )   )    ){
		  	    my ($EstID, $EstZoF)=@{ $bestCap3Hash->{'1_OrderHash'}->{$ordKey} };
		  	    
		  	    $new_WhoLine_Hash->{$new_WhoLine_idex}->{ $Cap3EstPrtKeyHash->{$ordKey} }=' ';
		  	  }
		    }
      
      	$new_WhoLine_idex++;
      }
      	
      
      
      # if the position here is a contig ACGT, then do the follow stuff
      if ( (defined $SplgnQryCtgPos) && (defined $SplgnQryCtgChr) && ($SplgnQryCtgPos=~m/^\d+$/) && ($SplgnQryCtgChr=~m/^[a-zA-Z]$/) )  {
      	my $cap3_ctgPos=$SplgnQryCtgPos-1; #in cap3 result, the first position is 0, but in splign result, the first position is 1     
      	if (   (  defined ( $betCap3_CtgPos2Wps )  ) && (  defined ( $betCap3_CtgPos2Wps->{$cap3_ctgPos} )  )   ){
      	  my $ResposWhoPos=$betCap3_CtgPos2Wps->{$cap3_ctgPos};
      	  
      	  
      	  $new_WhoLine_Hash->{$new_WhoLine_idex}->{$cap3Wl_c}=$ResposWhoPos;  #print "111 20180925 \$new_WhoLine_Hash->{$new_WhoLine_idex}->{$cap3Wl_c}=\$ResposWhoPos=$ResposWhoPos\n";
      	  if (   (  ref( $bestCap3Hash->{'_Align'} ) eq 'HASH'  ) && (  ref( $bestCap3Hash->{'1_OrderHash'} ) eq 'HASH'  )   ) {
		  	    foreach my $ordKey (    sort {$a <=> $b} (   keys (  %{ $bestCap3Hash->{'1_OrderHash'} }  )   )    ){
		  	      my ($EstID, $EstZoF)=@{ $bestCap3Hash->{'1_OrderHash'}->{$ordKey} };
		  	      my $Est_Each_Char=$bestCap3Hash->{'_Align'}->{$EstID}->{$EstZoF}->{'_WlPos2Word'}->{$ResposWhoPos};    if ($reverse_on_or_off==1) { $Est_Each_Char=FastaFileHandle::ReverseComplementBase($Est_Each_Char); }
		  	      $new_WhoLine_Hash->{$new_WhoLine_idex}->{ $Cap3EstPrtKeyHash->{$ordKey} }=$Est_Each_Char;               
		  	    }
		  	    
		      }
		  	  $new_WhoLine_idex++;  
      	  
      	  
      	  
      	  
      	  my $nextWhoLiPos=$ResposWhoPos+1;
      	  while (   (  defined ( $betCap3_Wps2DifPos->{$nextWhoLiPos} )  ) && ( $betCap3_Wps2DifPos->{$nextWhoLiPos}->[0] ne '_NactgType' )   ){
      	  	
      	  	
            $new_WhoLine_Hash->{$new_WhoLine_idex}->{$AA_char_1}=' ';
            $new_WhoLine_Hash->{$new_WhoLine_idex}->{$qr_char_2}=$betCap3_Wps2Difwrd->{$nextWhoLiPos};
            $new_WhoLine_Hash->{$new_WhoLine_idex}->{$CS_char_3}=' ';
            $new_WhoLine_Hash->{$new_WhoLine_idex}->{$sj_char_4}=' ';
            
            $new_WhoLine_Hash->{$new_WhoLine_idex}->{$Ex_numb_6}='x'; 
            $new_WhoLine_Hash->{$new_WhoLine_idex}->{$cap3Wl_c}=$nextWhoLiPos;   #print "222 20180925 \$new_WhoLine_Hash->{$new_WhoLine_idex}->{$cap3Wl_c}=\$nextWhoLiPos=$nextWhoLiPos\n";
            
            if (   (  ref( $bestCap3Hash->{'_Align'} ) eq 'HASH'  ) && (  ref( $bestCap3Hash->{'1_OrderHash'} ) eq 'HASH'  )   ) {
		  	      foreach my $ordKey (    sort {$a <=> $b} (   keys (  %{ $bestCap3Hash->{'1_OrderHash'} }  )   )    ){
		  	        my ($EstID, $EstZoF)=@{ $bestCap3Hash->{'1_OrderHash'}->{$ordKey} };
		  	        $new_WhoLine_Hash->{$new_WhoLine_idex}->{ $Cap3EstPrtKeyHash->{$ordKey} }=$bestCap3Hash->{'_Align'}->{$EstID}->{$EstZoF}->{'_WlPos2Word'}->{$nextWhoLiPos};
		  	      }
		        }
            
            $new_WhoLine_idex++;
            $nextWhoLiPos++;
      	  }
      	}  
      	else {
      		$new_WhoLine_idex++; 
      	}
      }
	  }
	}
	#print "333333333333333333333333333333  20180925\n";
	
	my @sortedBstCp3WpsPosArray= sort {$a<=>$b} (   keys (  %{ $betCap3_Wps2Difwrd }  )   );
	my $add_cap3_miss_pos_hash; my $newIdxKey=0; my $capWlPosStt=0;
	foreach my $newWLposKey (    sort {$a<=>$b} (   keys (  %{ $new_WhoLine_Hash }  )   )    ) {  #print "20180925 \$newWLposKey=$newWLposKey\n";
		
		if (   (  defined  (  $new_WhoLine_Hash->{$newWLposKey}->{$cap3Wl_c} )  ) && (  $new_WhoLine_Hash->{$newWLposKey}->{$cap3Wl_c}=~m/\d+/ )   ){    #print "20180925 \$new_WhoLine_Hash->{$newWLposKey}->{$cap3Wl_c}=$new_WhoLine_Hash->{$newWLposKey}->{$cap3Wl_c}\n";
		  my $cap3WlPos=$new_WhoLine_Hash->{$newWLposKey}->{$cap3Wl_c};                                             #print "20180925 \$cap3WlPos=\$new_WhoLine_Hash->{$newWLposKey}->{$cap3Wl_c}=$new_WhoLine_Hash->{$newWLposKey}->{$cap3Wl_c}\n";
			for ( my $i=$capWlPosStt; $i<$cap3WlPos; $i++ ){                                                          #print "aaaa20180925 \$capWlPosStt=$capWlPosStt,          \$i=$i,         <             \$cap3WlPos=$cap3WlPos\n";
				if (  defined ( $betCap3_Wps2Difwrd->{$i} )  ) {                                                        #print "bbbbb20180925 \$betCap3_Wps2Difwrd->{$i}=$betCap3_Wps2Difwrd->{$i}\n";
					$add_cap3_miss_pos_hash->{$newIdxKey}->{$cap3Wl_c }=$i ;
					$add_cap3_miss_pos_hash->{$newIdxKey}->{$AA_char_1}=' ';
					$add_cap3_miss_pos_hash->{$newIdxKey}->{$qr_char_2}=$betCap3_Wps2Difwrd->{$i};
					#########
					
					$add_cap3_miss_pos_hash->{$newIdxKey}->{$qr_posi_0}=$betCap3_Wps2DifPos->{$i}->[1]+1 if (   (  defined ( $betCap3_Wps2DifPos->{$i}->[0] )  )  && ( $betCap3_Wps2DifPos->{$i}->[0] eq '_NactgType' )    );
         # $add_cap3_miss_pos_hash->{$newIdxKey}->{$qr_ZorF_0}=$bstSpln_whoL__Hash->{$bstSpln_whoL_key}->{$qr_ZorF_0};     if ($reverse_on_or_off==1) { $add_cap3_miss_pos_hash->{$newIdxKey}->{$qr_ZorF_0}=&reverseZoF( $add_cap3_miss_pos_hash->{$newIdxKey}->{$qr_ZorF_0} ); }

					#########
					$add_cap3_miss_pos_hash->{$newIdxKey}->{$CS_char_3}=' ';
					$add_cap3_miss_pos_hash->{$newIdxKey}->{$sj_char_4}=' ';
					$add_cap3_miss_pos_hash->{$newIdxKey}->{$Ex_numb_6}=' ';
          if (   (  ref( $bestCap3Hash->{'_Align'} ) eq 'HASH'  ) && (  ref( $bestCap3Hash->{'1_OrderHash'} ) eq 'HASH'  )   ) {      #print "cccc20180925 \$bestCap3Hash->{'1_OrderHash'}=$bestCap3Hash->{'1_OrderHash'}\n";
		  	    foreach my $ordKey (    sort {$a <=> $b} (   keys (  %{ $bestCap3Hash->{'1_OrderHash'} }  )   )    ){                     #print "ddddd20180925 \$ordKey=$ordKey\n";
		  	      my ($EstID, $EstZoF)=@{ $bestCap3Hash->{'1_OrderHash'}->{$ordKey} };                                                    #print "eeeee20180925 \$EstID=$EstID  \$EstZoF=$EstZoF\n";
		  	      $add_cap3_miss_pos_hash->{$newIdxKey}->{ $Cap3EstPrtKeyHash->{$ordKey} }=$bestCap3Hash->{'_Align'}->{$EstID}->{$EstZoF}->{'_WlPos2Word'}->{$i};   #print "ffffff20180925 \$add_cap3_miss_pos_hash->{$newIdxKey}->{ \$Cap3EstPrtKeyHash->{$ordKey} }=\$add_cap3_miss_pos_hash->{$newIdxKey}->{ $Cap3EstPrtKeyHash->{$ordKey} }=\$bestCap3Hash->{'_Align'}->{$EstID}->{$EstZoF}->{'_WlPos2Word'}->{$i}=$bestCap3Hash->{'_Align'}->{$EstID}->{$EstZoF}->{'_WlPos2Word'}->{$i}\n";
		  	    }
		      }
            
          
            
					$newIdxKey++;
				}
			}		
			$capWlPosStt=$cap3WlPos+1;
		}
		else {			
		}
		#print "111fff20180925 \$add_cap3_miss_pos_hash->{$newIdxKey}=$add_cap3_miss_pos_hash->{$newIdxKey}\n";
		$add_cap3_miss_pos_hash->{$newIdxKey}=Storable::dclone  ( $new_WhoLine_Hash->{ $newWLposKey } );                       #print "222fff20180925 \$add_cap3_miss_pos_hash->{$newIdxKey}=$add_cap3_miss_pos_hash->{$newIdxKey}\n";
		$newIdxKey++;
		
	} 
	
	
	for ( my $i=$capWlPosStt; $i< @sortedBstCp3WpsPosArray; $i++){
		if (  defined ( $betCap3_Wps2Difwrd->{$i} )  ) {
			$add_cap3_miss_pos_hash->{$newIdxKey}->{$cap3Wl_c }=$i;
			$add_cap3_miss_pos_hash->{$newIdxKey}->{$AA_char_1}=' ';
			$add_cap3_miss_pos_hash->{$newIdxKey}->{$qr_char_2}=$betCap3_Wps2Difwrd->{$i};  if ($reverse_on_or_off==1) { $add_cap3_miss_pos_hash->{$newIdxKey}->{$qr_char_2}=FastaFileHandle::ReverseComplementBase( $add_cap3_miss_pos_hash->{$newIdxKey}->{$qr_char_2} ); }
			
			$add_cap3_miss_pos_hash->{$newIdxKey}->{$qr_posi_0}=$betCap3_Wps2DifPos->{$i}->[1]+1 if (   (  defined ( $betCap3_Wps2DifPos->{$i}->[0] )  )  && ( $betCap3_Wps2DifPos->{$i}->[0] eq '_NactgType'  )    );
      #$add_cap3_miss_pos_hash->{$newIdxKey}->{$qr_ZorF_0}=$bstSpln_whoL__Hash->{$bstSpln_whoL_key}->{$qr_ZorF_0};     if ($reverse_on_or_off==1) { $add_cap3_miss_pos_hash->{$newIdxKey}->{$qr_ZorF_0}=&reverseZoF( $add_cap3_miss_pos_hash->{$newIdxKey}->{$qr_ZorF_0} ); }

			
			$add_cap3_miss_pos_hash->{$newIdxKey}->{$CS_char_3}=' ';
			$add_cap3_miss_pos_hash->{$newIdxKey}->{$sj_char_4}=' ';
			$add_cap3_miss_pos_hash->{$newIdxKey}->{$Ex_numb_6}=' ';
      if (   (  ref( $bestCap3Hash->{'_Align'} ) eq 'HASH'  ) && (  ref( $bestCap3Hash->{'1_OrderHash'} ) eq 'HASH'  )   ) {
		    foreach my $ordKey (    sort {$a <=> $b} (   keys (  %{ $bestCap3Hash->{'1_OrderHash'} }  )   )    ){
		      my ($EstID, $EstZoF)=@{ $bestCap3Hash->{'1_OrderHash'}->{$ordKey} };
		      $add_cap3_miss_pos_hash->{$newIdxKey}->{ $Cap3EstPrtKeyHash->{$ordKey} }=$bestCap3Hash->{'_Align'}->{$EstID}->{$EstZoF}->{'_WlPos2Word'}->{$i};   #print "3333ffff20180925 \$add_cap3_miss_pos_hash->{$newIdxKey}->{ \$Cap3EstPrtKeyHash->{$ordKey} }=\$add_cap3_miss_pos_hash->{$newIdxKey}->{ $Cap3EstPrtKeyHash->{$ordKey} }=\$bestCap3Hash->{'_Align'}->{$EstID}->{$EstZoF}->{'_WlPos2Word'}->{$i}=$bestCap3Hash->{'_Align'}->{$EstID}->{$EstZoF}->{'_WlPos2Word'}->{$i}\n";
		      if ($reverse_on_or_off==1) { $add_cap3_miss_pos_hash->{$newIdxKey}->{ $Cap3EstPrtKeyHash->{$ordKey} }=FastaFileHandle::ReverseComplementBase( $add_cap3_miss_pos_hash->{$newIdxKey}->{ $Cap3EstPrtKeyHash->{$ordKey} } ); }
		    }
	    } 
			$newIdxKey++;
		}
	}	
	
	my $outPutHash;
	if ( ref ($add_cap3_miss_pos_hash) eq 'HASH' ){
		foreach my $newWLposKey (    sort {$a<=>$b} (   keys (  %{ $add_cap3_miss_pos_hash }  )   )    ) {
			                                                                                                                                         $outPutHash->{'5_alnShwHs'}->{$newWLposKey}->{$ctgPos_a}=$add_cap3_miss_pos_hash->{$newWLposKey}->{$qr_posi_0} if (   (  defined ( $add_cap3_miss_pos_hash->{$newWLposKey}->{$qr_char_2} )  ) && ( $add_cap3_miss_pos_hash->{$newWLposKey}->{$qr_char_2}=~/\w/ )    );
			                                                                                                                                         $outPutHash->{'5_alnShwHs'}->{$newWLposKey}->{$gnoPos_b}=$add_cap3_miss_pos_hash->{$newWLposKey}->{$sj_posi_5} if (   (  defined ( $add_cap3_miss_pos_hash->{$newWLposKey}->{$sj_char_4} )  ) && ( $add_cap3_miss_pos_hash->{$newWLposKey}->{$sj_char_4}=~/\w/ )    );    	
			                                                                                                                                         $outPutHash->{'5_alnShwHs'}->{$newWLposKey}->{$ctgZoF_a}=$add_cap3_miss_pos_hash->{$newWLposKey}->{$qr_ZorF_0};
			                                                                                                                                         $outPutHash->{'5_alnShwHs'}->{$newWLposKey}->{$gnoZoF_b}=$add_cap3_miss_pos_hash->{$newWLposKey}->{$sj_ZorF_5};    	

			
			$outPutHash->{'0___alnShw'}->{$spl_AA_4}.=$add_cap3_miss_pos_hash->{$newWLposKey}->{$AA_char_1};                                               $outPutHash->{'5_alnShwHs'}->{$newWLposKey}->{$spl_AA_4}=$add_cap3_miss_pos_hash->{$newWLposKey}->{$AA_char_1};
    	$outPutHash->{'0___alnShw'}->{$splCtg_3}.=$add_cap3_miss_pos_hash->{$newWLposKey}->{$qr_char_2};                                               $outPutHash->{'5_alnShwHs'}->{$newWLposKey}->{$splCtg_3}=$add_cap3_miss_pos_hash->{$newWLposKey}->{$qr_char_2};
    	$outPutHash->{'0___alnShw'}->{$spl_CS_2}.=$add_cap3_miss_pos_hash->{$newWLposKey}->{$CS_char_3};                                               $outPutHash->{'5_alnShwHs'}->{$newWLposKey}->{$spl_CS_2}=$add_cap3_miss_pos_hash->{$newWLposKey}->{$CS_char_3};
    	$outPutHash->{'0___alnShw'}->{$splGno_0}.=$add_cap3_miss_pos_hash->{$newWLposKey}->{$sj_char_4};                                               $outPutHash->{'5_alnShwHs'}->{$newWLposKey}->{$splGno_0}=$add_cap3_miss_pos_hash->{$newWLposKey}->{$sj_char_4};
    	$outPutHash->{'0___alnShw'}->{$splExN_1}.=$add_cap3_miss_pos_hash->{$newWLposKey}->{$Ex_numb_6};                                               $outPutHash->{'5_alnShwHs'}->{$newWLposKey}->{$splExN_1}=$add_cap3_miss_pos_hash->{$newWLposKey}->{$Ex_numb_6};
    	if (   (  ref( $bestCap3Hash->{'_Align'} ) eq 'HASH'  ) && (  ref( $bestCap3Hash->{'1_OrderHash'} ) eq 'HASH'  )   ) {
		    foreach my $ordKey (    sort {$a <=> $b} (   keys (  %{ $bestCap3Hash->{'1_OrderHash'} }  )   )    ){
		    	$outPutHash->{'0___alnShw'}->{ $Cap3EstPrtKeyHash->{$ordKey} }.=$add_cap3_miss_pos_hash->{$newWLposKey}->{ $Cap3EstPrtKeyHash->{$ordKey} };  $outPutHash->{'5_alnShwHs'}->{$newWLposKey}->{ $Cap3EstPrtKeyHash->{$ordKey} }=$add_cap3_miss_pos_hash->{$newWLposKey}->{ $Cap3EstPrtKeyHash->{$ordKey} };
		    }
	    }
		}
	}          
	
	my $EstNameHash;
	if (   (  ref( $bestCap3Hash->{'_Align'} ) eq 'HASH'  ) && (  ref( $bestCap3Hash->{'1_OrderHash'} ) eq 'HASH'  )   ) {
		foreach my $ordKey (    sort {$a <=> $b} (   keys (  %{ $bestCap3Hash->{'1_OrderHash'} }  )   )    ){
		  my ($EstID, $EstZoF)=@{ $bestCap3Hash->{'1_OrderHash'}->{$ordKey} };
		  $EstNameHash->{ $Cap3EstPrtKeyHash->{$ordKey} }=[$EstID, $EstZoF];
	  }
	}
	
	$outPutHash->{'1_Geno_ZoF'}=$realSplnGnoZF;
	$outPutHash->{'2_EstIdHsh'}=$EstNameHash;
	$outPutHash->{'4_whoL_Hsh'}=$add_cap3_miss_pos_hash;
	
	
	if ( ref ($bstSpln_subjG_Hash) eq 'HASH' ){
		my $outSegMetHash; 
		my $lastStateHash; my $frcIdx=0;
	  foreach my $genoPosK (    sort {$a<=>$b} (   keys (  %{ $bstSpln_subjG_Hash }  )   )    ) {
	  	my $genoChar=$bstSpln_subjG_Hash->{$genoPosK}->{'0_self_sj_char'};   if ($splnCtgQuryZF eq '-') { $genoChar=FastaFileHandle::ReverseComplementBase($genoChar); }  if ($reverse_on_or_off==1) { $genoChar=FastaFileHandle::ReverseComplementBase($genoChar); }
	  	my $CotgChar=$bstSpln_subjG_Hash->{$genoPosK}->{'1_mach_qr_char'};   if ($splnCtgQuryZF eq '-') { $CotgChar=FastaFileHandle::ReverseComplementBase($CotgChar); }  if ($reverse_on_or_off==1) { $CotgChar=FastaFileHandle::ReverseComplementBase($CotgChar); }
	  	my $CotgPosi; $CotgPosi=$bstSpln_subjG_Hash->{$genoPosK}->{'2_mach_qr_posi'} if (  defined ( $bstSpln_subjG_Hash->{$genoPosK}->{'2_mach_qr_posi'} )  );               #print "\$CotgPosi=\$bstSpln_subjG_Hash->{$genoPosK}->{'2_mach_qr_posi'}=$CotgPosi\n\n";
	  	my $stageMet=&CtgVSGenoSttMentDefine($genoChar, $CotgChar);
	  	
	  	my $newSegSt=0;
	  	if ($stageMet ne '5_OTHER_TYP' ){
	  		if    ($frcIdx==0)                                                                                                                                                                 { $newSegSt=1; }
	  		elsif (   (  defined ( $lastStateHash->{'6ctgVSgn'}->{'stageMet'} )  ) && ($lastStateHash->{'6ctgVSgn'}->{'stageMet'} ne $stageMet)   )                                            { $newSegSt=1; }
	  		elsif (    (   (  abs ( $genoPosK-$lastStateHash->{'6ctgVSgn'}->{'genoPosK'} )  ) > 1   ) ||  (   (  abs ( $CotgPosi-$lastStateHash->{'6ctgVSgn'}->{'CotgPosi'} )  ) > 1    )     ){ $newSegSt=1; }
	  		else     {	  			$newSegSt=0;	  		}
	  	}
	  	if ($newSegSt){  
	  	  $lastStateHash->{'6ctgVSgn'}->{'startPos'}=$genoPosK;
	  		$outSegMetHash->{'6ctgVSgn'}->{ $lastStateHash->{'6ctgVSgn'}->{'startPos'} }->{$stageMet}=$genoPosK;  
	  	}
	  			
	  	
	  	my $cotinues=0;	  	
	  	if ( 
	  	     ( $stageMet ne '5_OTHER_TYP' ) 
	  	     &&
	  	     (   (  defined ( $lastStateHash->{'6ctgVSgn'}->{'stageMet'} )  ) && ($lastStateHash->{'6ctgVSgn'}->{'stageMet'} eq $stageMet)   )   
	  	     &&
	  	     (    (   (  abs ( $genoPosK-$lastStateHash->{'6ctgVSgn'}->{'genoPosK'} )  ) == 1    ) && (   (  abs ( $CotgPosi-$lastStateHash->{'6ctgVSgn'}->{'CotgPosi'} )  ) == 1    )     )
	  	)  	{ $cotinues=1; }	  		
	  	if ($cotinues){  
	  		$outSegMetHash->{'6ctgVSgn'}->{ $lastStateHash->{'6ctgVSgn'}->{'startPos'} }->{$stageMet}=$genoPosK;  
	  	}    
	  	
	  	$lastStateHash->{'6ctgVSgn'}->{'stageMet'}=$stageMet;
	  	$lastStateHash->{'6ctgVSgn'}->{'genoPosK'}=$genoPosK;
	  	$lastStateHash->{'6ctgVSgn'}->{'CotgPosi'}=$CotgPosi;
	  	
	  	
	  	if (    (  defined ( $bestCap3Hash )  ) && (  ref ( $bestCap3Hash ) eq 'HASH'  ) && (  defined ( $bestCap3Hash->{'1_OrderHash'} )  ) && (  ref ( $bestCap3Hash->{'1_OrderHash'} ) eq 'HASH'  )    ){
	  	  foreach my $ordKey (    sort {$a <=> $b} (   keys (  %{ $bestCap3Hash->{'1_OrderHash'} }  )   )    ){
	  	  	my ($EstID, $EstZoF)=@{ $bestCap3Hash->{'1_OrderHash'}->{$ordKey} };
	  	  	my $EstIdxID=$Cap3EstPrtKeyHash->{$ordKey};     
	  	  	
	  	  	#my $Cap3ResposWhoPos=$betCap3_CtgPos2Wps->{$CotgPosi-1};
	  	  	#my $Est_Char=$bestCap3Hash->{'_Align'}->{$EstID}->{$EstZoF}->{'_WlPos2Word'}->{$Cap3ResposWhoPos};
	  	  	#my $Est_posi=$bestCap3Hash->{'_Align'}->{$EstID}->{$EstZoF}->{'_Wpos2DifTypePos'}->{$Cap3ResposWhoPos}->[1];
	  	  	
	  	  	############
	  	  	
	  	  	my $Cap3ResposWhoPos; $Cap3ResposWhoPos=$betCap3_CtgPos2Wps->{$CotgPosi-1} if (   ( defined ($CotgPosi) ) &&  ($CotgPosi=~m/\d+/) && (  defined ( $betCap3_CtgPos2Wps->{$CotgPosi-1} )  ) && ( $betCap3_CtgPos2Wps->{$CotgPosi-1}=~m/\d+/ )    );
	  	  	my $Est_Char;  $Est_Char=$bestCap3Hash->{'_Align'}->{$EstID}->{$EstZoF}->{'_WlPos2Word'}->{$Cap3ResposWhoPos}           if (   ( defined ($Cap3ResposWhoPos) ) &&  ($Cap3ResposWhoPos=~m/\d+/) && (  defined ( $bestCap3Hash->{'_Align'}->{$EstID}->{$EstZoF}->{'_WlPos2Word'}->{$Cap3ResposWhoPos}           )  )   );
	  	  	if ($reverse_on_or_off==1) { $Est_Char=FastaFileHandle::ReverseComplementBase($Est_Char); }
	  	  	my $Est_posi;  $Est_posi=$bestCap3Hash->{'_Align'}->{$EstID}->{$EstZoF}->{'_Wpos2DifTypePos'}->{$Cap3ResposWhoPos}->[1] if (   ( defined ($Cap3ResposWhoPos) ) &&  ($Cap3ResposWhoPos=~m/\d+/) && (  defined ( $bestCap3Hash->{'_Align'}->{$EstID}->{$EstZoF}->{'_Wpos2DifTypePos'}->{$Cap3ResposWhoPos}->[1] )  )   );
	  	  	
	  	  	
	  	  	############
	  	  	
	  	  	
	  	  	my $EstStMet=&EstVSCtgVSGenoSttMentDefine($genoChar, $CotgChar,$Est_Char);
	  	  	
	  	  	###############################
	  	  	
	  	  	my $EstNewSegSt=0;
	  	    if ($EstStMet ne '5_OTHER_TYP' ){
	  	    	if    ($frcIdx==0)                                                                                                                                                                 { $EstNewSegSt=1; }
	  	    	elsif (   (  defined ( $lastStateHash->{$EstIdxID}->{'stageMet'} )  ) && ($lastStateHash->{$EstIdxID}->{'stageMet'} ne $EstStMet)   )                                              { $EstNewSegSt=1; }
	  	    	elsif (    (   (  abs ( $genoPosK-$lastStateHash->{$EstIdxID}->{'genoPosK'} )  ) > 1   ) ||  (   (  abs ( $Est_posi-$lastStateHash->{$EstIdxID}->{'Est_posi'} )  ) > 1    )     )  { $EstNewSegSt=1; }
	  	    	else     {	  			$EstNewSegSt=0;	  		}
	  	    }
	  	    if ($EstNewSegSt){  
	  	      $lastStateHash->{$EstIdxID}->{'startPos'}=$genoPosK;
	  	    	$outSegMetHash->{$EstIdxID}->{ $lastStateHash->{$EstIdxID}->{'startPos'} }->{$EstStMet}=$genoPosK;  
	  	    }
	  	    		
	  	    
	  	    my $EstCotinues=0;	  	
	  	    if ( 
	  	         ( $EstStMet ne '5_OTHER_TYP' ) 
	  	         &&
	  	         (   (  defined ( $lastStateHash->{$EstIdxID}->{'stageMet'} )  ) && ($lastStateHash->{$EstIdxID}->{'stageMet'} eq $EstStMet)   )   
	  	         &&
	  	         (    (   (  abs ( $genoPosK-$lastStateHash->{$EstIdxID}->{'genoPosK'} )  ) == 1    ) && (   (  abs ( $Est_posi-$lastStateHash->{$EstIdxID}->{'Est_posi'} )  ) == 1    )     )
	  	    )  	{ $EstCotinues=1; }	  		
	  	    if ($EstCotinues){  
	  	    	$outSegMetHash->{$EstIdxID}->{ $lastStateHash->{$EstIdxID}->{'startPos'} }->{$EstStMet}=$genoPosK;  
	  	    }    
	  	    
	  	    $lastStateHash->{$EstIdxID}->{'stageMet'}=$EstStMet;
	  	    $lastStateHash->{$EstIdxID}->{'genoPosK'}=$genoPosK;
	  	    $lastStateHash->{$EstIdxID}->{'Est_posi'}=$Est_posi;
	  	  	
	  	  	###############################
	  	  	
	  	  }
	  	}
	  	$frcIdx++;
	  }
	  
	  
	  if ( ref ($outSegMetHash) eq 'HASH'){
	  	my $EstMatchGeneStructure_Hash;
	  	foreach my $keyHere (    sort {$a cmp $b} (   keys (  %{ $outSegMetHash }  )   )    ){
	  	  if (  ref ( $outSegMetHash->{$keyHere} ) eq 'HASH'  ){                                                               
          foreach my $sgmHead (    sort {$a <=> $b} (   keys (  %{ $outSegMetHash->{$keyHere} }  )   )    ){  
          	my ($sgmType)=(   keys (  %{ $outSegMetHash->{$keyHere}->{$sgmHead} }  )   );
          	my $sgmTail=$outSegMetHash->{$keyHere}->{$sgmHead}->{$sgmType};
          	my $tempHash;
          	$tempHash->{'0_SgmtHead'}=$sgmHead;
          	$tempHash->{'1_SgmtTail'}=$sgmTail;
          	$tempHash->{'3_SgmtType'}=$sgmType;
          	
          	push @{ $EstMatchGeneStructure_Hash->{$keyHere} }, $tempHash;
          }
        }          
	  	}
	  	$outPutHash->{'3_GStr_Hsh'}=$EstMatchGeneStructure_Hash;
	  }	  
	  
	}
	
	
	return $outPutHash;
	
}

sub BuildCap3OrderKeyHash{
  my ($inBestCap3Hash)=@_;
  
  my $warnMsgBody="\nIn package  SplignHandle,\tIn sub BuildCap3OrderKeyHash,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";

  
  my $outHASH; 
  if (   (  ref( $inBestCap3Hash->{'_Align'} ) eq 'HASH'  ) && (  ref( $inBestCap3Hash->{'1_OrderHash'} ) eq 'HASH'  )   ) {
	  foreach my $ordKey (    sort {$a <=> $b} (   keys (  %{ $inBestCap3Hash->{'1_OrderHash'} }  )   )    ){
	    my ($EstID, $EstZoF)=@{ $inBestCap3Hash->{'1_OrderHash'}->{$ordKey} };
	    my $outKeySpritf=sprintf ("%04d", $ordKey); $outKeySpritf="7est_$outKeySpritf";
	    $outHASH->{$ordKey}=$outKeySpritf;
	    
	  }
	}
	return $outHASH;
}

sub reverseZoF{
  my ($inZF)=@_;
  
  my $warnMsgBody="\nIn package  SplignHandle,\tIn sub reverseZoF,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	
  #my $outZF;
  if  (  ( defined ( $inZF ) )&& (  $inZF=~m(-|\+)  )  ){
  	if    ($inZF eq '-'){ $inZF='+'; } 
  	elsif ($inZF eq '+'){ $inZF='-'; } 
  }
  return $inZF;
  
  #else { my $dieMsg="$die_MsgHead\nThe \$inZF=$inZF should be + or -!!!\n"; print $dieMsg; die $dieMsg; };
  
  
}


sub CtgVSGenoSttMentDefine{
	my ($In_genoChar, $In_CotgChar)=@_;
	my $outState;
	$In_genoChar= uc $In_genoChar if (  ( defined ($In_genoChar) ) && ( $In_genoChar=~m/^[a-zA-Z]$/ )  ); 
	$In_CotgChar= uc $In_CotgChar if (  ( defined ($In_CotgChar) ) && ( $In_CotgChar=~m/^[a-zA-Z]$/ )  ); 
	if    (  ( defined ($In_genoChar) ) && ( defined ($In_CotgChar) ) && ($In_genoChar=~m/^[a-zA-Z]$/) && ($In_CotgChar=~m/^[a-zA-Z]$/) && ($In_genoChar eq $In_CotgChar) ){	$outState='0_Geo_E_Ctg';	}
	elsif (  ( defined ($In_CotgChar) ) && ( defined ($In_genoChar) ) && ($In_genoChar=~m/^[a-zA-Z]$/) && ($In_CotgChar=~m/^[a-zA-Z]$/) && ($In_genoChar ne $In_CotgChar) ){	$outState='1_Geo_N_Ctg';	}
	else                                                                                                      {	$outState='5_OTHER_TYP';	}
	return $outState;
}

sub EstVSCtgVSGenoSttMentDefine{
	my ($In_genoChar, $In_CotgChar, $In_Est_Char)=@_;
	my $outState;
	$In_genoChar= uc $In_genoChar if (  ( defined ($In_genoChar) ) && ( $In_genoChar=~m/^[a-zA-Z]$/ )  ); 
	$In_CotgChar= uc $In_CotgChar if (  ( defined ($In_CotgChar) ) && ( $In_CotgChar=~m/^[a-zA-Z]$/ )  ); 
	$In_Est_Char= uc $In_Est_Char if (  ( defined ($In_Est_Char) ) && ( $In_Est_Char=~m/^[a-zA-Z]$/ )  ); 
	if    ( ( defined ($In_genoChar) ) && ( defined ($In_CotgChar) ) && ( defined ($In_Est_Char) ) && ($In_genoChar=~m/^[a-zA-Z]$/) && ($In_CotgChar=~m/^[a-zA-Z]$/) && ($In_Est_Char=~m/^[a-zA-Z]$/) && ($In_genoChar eq $In_CotgChar) && ($In_genoChar eq $In_Est_Char) ){	$outState='0_Geo_E_Ctg_E_Est';	}
	elsif ( ( defined ($In_genoChar) ) && ( defined ($In_CotgChar) ) && ( defined ($In_Est_Char) ) && ($In_genoChar=~m/^[a-zA-Z]$/) && ($In_CotgChar=~m/^[a-zA-Z]$/) && ($In_Est_Char=~m/^[a-zA-Z]$/) && ($In_genoChar ne $In_CotgChar) && ($In_genoChar eq $In_Est_Char) ){	$outState='1_Geo_E_Est_N_Ctg';	}
	elsif ( ( defined ($In_genoChar) ) && ( defined ($In_CotgChar) ) && ( defined ($In_Est_Char) ) && ($In_genoChar=~m/^[a-zA-Z]$/) && ($In_CotgChar=~m/^[a-zA-Z]$/) && ($In_Est_Char=~m/^[a-zA-Z]$/) && ($In_genoChar ne $In_CotgChar) && ($In_CotgChar eq $In_Est_Char) ){	$outState='2_Geo_N_Ctg_E_Est';	}
	elsif ( ( defined ($In_genoChar) ) && ( defined ($In_CotgChar) ) && ( defined ($In_Est_Char) ) && ($In_genoChar=~m/^[a-zA-Z]$/) && ($In_CotgChar=~m/^[a-zA-Z]$/) && ($In_Est_Char=~m/^[a-zA-Z]$/) && ($In_genoChar eq $In_CotgChar) && ($In_CotgChar ne $In_Est_Char) ){	$outState='3_Geo_E_Ctg_N_Est';	}
	elsif ( ( defined ($In_genoChar) ) && ( defined ($In_CotgChar) ) && ( defined ($In_Est_Char) ) && ($In_genoChar=~m/^[a-zA-Z]$/) && ($In_CotgChar=~m/^[a-zA-Z]$/) && ($In_Est_Char=~m/^[a-zA-Z]$/) && ($In_genoChar ne $In_CotgChar) && ($In_CotgChar ne $In_Est_Char) ){	$outState='4_Geo_N_Ctg_N_Est';	}
	else                                                                                                                                                                         {	$outState='5_OTHER_TYP';	      }
	return $outState;
}


sub EstVSCtgSttMentDefine{
	my ($In_CotgChar, $In_Est_Char)=@_;
	my $outState;
	$In_CotgChar= uc $In_CotgChar if (  ( defined ($In_CotgChar) ) && ( $In_CotgChar=~m/^[a-zA-Z]$/ )  );
	$In_Est_Char= uc $In_Est_Char if (  ( defined ($In_Est_Char) ) && ( $In_Est_Char=~m/^[a-zA-Z]$/ )  );
	if    ( ( defined ($In_Est_Char) ) && ( defined ($In_CotgChar) ) && ($In_Est_Char=~m/^[a-zA-Z]$/) && ($In_CotgChar=~m/^[a-zA-Z]$/) && ($In_Est_Char eq $In_CotgChar) ){	$outState='0_Ctg_E_Est';	}
	elsif ( ( defined ($In_CotgChar) ) && ( defined ($In_Est_Char) ) && ($In_Est_Char=~m/^[a-zA-Z]$/) && ($In_CotgChar=~m/^[a-zA-Z]$/) && ($In_Est_Char ne $In_CotgChar) ){	$outState='1_Ctg_N_Est';	}
	else                                                                                                      {	$outState='5_OTHER_TYP';	}
	return $outState;
}

sub Map_est_sgmt_2_div_array_bk_to_geno{ #SplignHandle::Map_est_sgmt_2_div_array_bk_to_geno
	my ($InSplignOUTarray, $InGeneStructure_2d_Array)=@_;
	
	my $warnMsgBody="\nIn package  SplignHandle,\tIn sub Map_est_sgmt_2_div_array_bk_to_geno,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	
	my $outArray;
	if ( ref ($InGeneStructure_2d_Array) eq 'ARRAY' ){
		for (  my $i=0; $i < @{ $InGeneStructure_2d_Array }; $i++  ){
			if (  ref ( $InGeneStructure_2d_Array->[$i] ) eq 'ARRAY'  ){
			  my $mpBkArray= &Map_est_gene_strucure_back_to_geno($InSplignOUTarray, $InGeneStructure_2d_Array->[$i] );
			  $outArray->[$i]=$mpBkArray;
			}
			else {
				my $dumpHere=DirFileHandle::ReturnDumperInform ( $InGeneStructure_2d_Array->[$i] );
				my $dieMsg_1="$die_MsgHead\nThe \$InGeneStructure_2d_Array->[$i]=$InGeneStructure_2d_Array->[$i] should be a ARRAY ref !!!!\n\n$dumpHere\n\n";
  	    print $dieMsg_1; die $dieMsg_1;
			}
		}
	}
	else {
		my $dumpHere=DirFileHandle::ReturnDumperInform ( $InGeneStructure_2d_Array );
	  my $dieMsg2="$die_MsgHead\nThe \$InGeneStructure_2d_Array=$InGeneStructure_2d_Array should be a ARRAY ref !!!!\n\n$dumpHere\n\n";
    print $dieMsg2; die $dieMsg2;
	}
	return $outArray;
	
}

sub Map_est_gene_strucure_back_to_geno{ #SplignHandle::Map_est_gene_strucure_back_to_geno
  my ($InSplignOUTarray, $InGeneStructureArray)=@_;
  
  my $warnMsgBody="\nIn package  SplignHandle,\tIn sub Map_est_gene_strucure_back_to_geno,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	
	my $bestSplignHash;
	if (   ( ref ($InSplignOUTarray) eq 'ARRAY' ) && (  ref ( $InSplignOUTarray->[0] ) eq 'HASH'  )   ){ $bestSplignHash=$InSplignOUTarray->[0];	}
	else {my $splignArrayDump=DirFileHandle::ReturnDumperInform($InSplignOUTarray);my $dieMsg="$die_MsgHead\nThe \$InSplignOUTarray=$InSplignOUTarray should be a correct array ref\n\$splignArrayDump=$splignArrayDump\n\n\n"; print $dieMsg; die $dieMsg;} 

	my $splnCtgQuryZF=$bestSplignHash->{'3_Qury_ZoF'};
	my $splnGnoSubjZF=$bestSplignHash->{'5_Subj_ZoF'};                  
           
  
  my $realSplnGnoZF=$splnGnoSubjZF; if ($splnCtgQuryZF eq '-'){$realSplnGnoZF=&reverseZoF($realSplnGnoZF);}
  my $bstSpln_qury__Hash=$bestSplignHash->{'8_qurM_Hsh'};
  
  
  if ( ref ($InGeneStructureArray) eq 'ARRAY' ){
  	my $pointArray;
    my $pointIdx=0;
  	for (  my $i=0; $i<@{ $InGeneStructureArray }; $i++  ){
  		#$outPutArray->[$i]=$InGeneStructureArray->[$i];
  		
  		
  		my $realSgmZoF=$InGeneStructureArray->[$i]->{'2_SgmtZorF'}; 
  		if     (  ( ($splnCtgQuryZF eq '+') || ($splnCtgQuryZF eq '-') ) && ( ($splnGnoSubjZF eq '+') || ($splnGnoSubjZF eq '-') ) && ($splnCtgQuryZF eq $splnGnoSubjZF)  ){  			
  		}
  		elsif  (  ( ($splnCtgQuryZF eq '+') || ($splnCtgQuryZF eq '-') ) && ( ($splnGnoSubjZF eq '+') || ($splnGnoSubjZF eq '-') ) && ($splnCtgQuryZF ne $splnGnoSubjZF)  ){
  			$realSgmZoF=&reverseZoF($realSgmZoF);  			
  		}
  		else {
  			my $dieMsg="$die_MsgHead\$splnCtgQuryZF=$splnCtgQuryZF\t\$splnGnoSubjZF=$splnGnoSubjZF these two should be + or -!!\n\n\n"; print  $dieMsg; die $dieMsg;
  		}
  		
  		my $head= $InGeneStructureArray->[$i]->{'0_SgmtHead'}; print "a-20180906\$head= \$InGeneStructureArray->[\$i]->{'0_SgmtHead'}= \$InGeneStructureArray->[$i]->{'0_SgmtHead'}=$head\n";
  		my $tail= $InGeneStructureArray->[$i]->{'1_SgmtTail'}; print "b-20180906\$tail= \$InGeneStructureArray->[\$i]->{'1_SgmtTail'}= \$InGeneStructureArray->[$i]->{'1_SgmtTail'}=$tail\n";
  		my $dirc=-1; if ($tail>=$head){$dirc=1;}
  		for (   my $k=0; $k<=( abs ($tail-$head) ); $k++   ){
  			my $pointPosition=$head+$dirc*$k;
  			$pointArray->[$pointIdx]=Storable::dclone ($InGeneStructureArray->[$i]);
  			$pointArray->[$pointIdx]->{'0_SgmtHead'}=$InGeneStructureArray->[$i]->{'0_SgmtHead'};
  			$pointArray->[$pointIdx]->{'1_SgmtTail'}=$InGeneStructureArray->[$i]->{'1_SgmtTail'};
  			$pointArray->[$pointIdx]->{'2_SgmtZorF'}=$realSgmZoF;
  			$pointArray->[$pointIdx]->{'8_del_StPs'}=$pointPosition;   print "0-20180906\$pointArray->[\$pointIdx]->{'8_del_StPs'}=\$pointPosition=\$pointArray->[$pointIdx]->{'8_del_StPs'}=$pointPosition\n";
  			$pointArray->[$pointIdx]->{'9_del_SpPs'}=$bstSpln_qury__Hash->{ $pointPosition }->{'2_mach_sj_posi'} if (  defined ( $bstSpln_qury__Hash->{ $pointPosition }->{'2_mach_sj_posi'} )  );  #print "c-20180906\$pointArray->[\$pointIdx]->{'9_del_SpPs'}=\$bstSpln_qury__Hash->{ \$pointPosition }->{'2_mach_sj_posi'}=\$pointArray->[$pointIdx]->{'9_del_SpPs'}=\$bstSpln_qury__Hash->{ $pointPosition }->{'2_mach_sj_posi'}=$pointArray->[$pointIdx]->{'9_del_SpPs'}\n"; 
  			my $pointArray_idx_dump= DirFileHandle::ReturnDumperInform($pointArray->[$pointIdx]); print "d-20180906\$pointArray->[\$pointIdx]=\$pointArray->[$pointIdx]=$pointArray->[$pointIdx]:\n$pointArray_idx_dump\n";
  			$pointIdx++;
  		}
  		
  		my $pointArray_dump= DirFileHandle::ReturnDumperInform($pointArray); print "1-20180906$pointArray_dump\n";
  		#$outPutArray->[$i]->{'0_SgmtHead'}=$bstSpln_qury__Hash->{ $InGeneStructureArray->[$i]->{'0_SgmtHead'} }->{'2_mach_sj_posi'};
  		#$outPutArray->[$i]->{'1_SgmtTail'}=$bstSpln_qury__Hash->{ $InGeneStructureArray->[$i]->{'1_SgmtTail'} }->{'2_mach_sj_posi'};
  		#$outPutArray->[$i]->{'2_SgmtZorF'}=$realSgmZoF;
  		
  	}
  	
  	my $SgmTypeHash;
  	my $lastPos; my $last_SgmtType; my $last_start_idx;
  	for (  my $j=0; $j<@{ $pointArray }; $j++  ){
  		
  		if ($j==0){
  		  if (   (  defined ( $pointArray->[$j]->{'9_del_SpPs'} )  ) && ( $pointArray->[$j]->{'9_del_SpPs'}=~m/\d+/ ) && (  defined ( $pointArray->[$j]->{'3_SgmtType'} )  ) && ( $pointArray->[$j]->{'3_SgmtType'}=~m/\S+/ )    ){
  		  	$SgmTypeHash->{ $j }->{ $pointArray->[$j]->{'3_SgmtType'} }=$j;                       print "2-20180906\$SgmTypeHash->{ $j }->{ \$pointArray->[$j]->{'3_SgmtType'} }=\$SgmTypeHash->{ $j }->{ $pointArray->[$j]->{'3_SgmtType'} }=\$j=$SgmTypeHash->{ $j }->{ $pointArray->[$j]->{'3_SgmtType'} }\n";
  		  	$last_start_idx=$j;
  		  }
  		}
  		else {
  			if (   (  defined ( $pointArray->[$j]->{'9_del_SpPs'} )  ) && ( $pointArray->[$j]->{'9_del_SpPs'}=~m/\d+/ ) && (  defined ( $pointArray->[$j]->{'3_SgmtType'} )  ) && ( $pointArray->[$j]->{'3_SgmtType'}=~m/\S+/ )   ){
  				if (  ($lastPos=~m/\d+/) && ( abs ( $pointArray->[$j]->{'9_del_SpPs'} - $lastPos) == 1) && ( $pointArray->[$j]->{'3_SgmtType'} eq $last_SgmtType)  ){ #                           
  					$SgmTypeHash->{ $last_start_idx }->{ $pointArray->[$j]->{'3_SgmtType'} }=$j;        print "3-20180906\$SgmTypeHash->{ \$last_start_idx }->{ \$pointArray->[\$j]->{'3_SgmtType'} }=\$SgmTypeHash->{ $last_start_idx }->{ $pointArray->[$j]->{'3_SgmtType'} }=\$j=$SgmTypeHash->{ $last_start_idx }->{ $pointArray->[$j]->{'3_SgmtType'} }\n";
  				}
  				else {
  					$SgmTypeHash->{ $j }->{ $pointArray->[$j]->{'3_SgmtType'} }=$j;                  print "e-20180906\$SgmTypeHash->{ \$j }->{ \$pointArray->[$j]->{'3_SgmtType'} }=\$SgmTypeHash->{ $j }->{ $pointArray->[$j]->{'3_SgmtType'} }=\$j=$j\n";
  					$last_start_idx=$j;
  				}
  			}
  			
  			
  		}
  		$lastPos='';       $lastPos=      $pointArray->[$j]->{'9_del_SpPs'} if (  defined ( $pointArray->[$j]->{'9_del_SpPs'} )  );  #print "4-20180906\$lastPos=$lastPos\$pointArray->[$j]->{'9_del_SpPs'}=$pointArray->[$j]->{'9_del_SpPs'}\n";
  		$last_SgmtType=''; $last_SgmtType=$pointArray->[$j]->{'3_SgmtType'} if (  defined ( $pointArray->[$j]->{'3_SgmtType'} )  );  print "5-20180906\$last_SgmtType=$last_SgmtType\$pointArray->[$j]->{'3_SgmtType'}=$pointArray->[$j]->{'3_SgmtType'}\n";
  		
  		
  	}
  	
  	my $SgmTypeHash_dump= DirFileHandle::ReturnDumperInform($SgmTypeHash); print "f-20180906 dump \$SgmTypeHash=$SgmTypeHash:\n$SgmTypeHash_dump\n";
  	
  	my $outPutArray; my $outPutArray_idx=0; 
  	if (   (  defined ( $SgmTypeHash )  ) && ( ref ($SgmTypeHash) eq 'HASH' )   ){         print "g-20180906 \$SgmTypeHash=$SgmTypeHash\n";
  	  foreach my $idxKey (    sort {$a<=>$b} (   keys (  %{ $SgmTypeHash }  )   )    ){    print "6-20180906 \$idxKey=$idxKey\n";
  	  	
  		  if (   (  defined ( $SgmTypeHash->{$idxKey} )  ) && (  ref ( $SgmTypeHash->{$idxKey} ) eq 'HASH'  )   ){
  		    foreach my $sgmTpKy (    sort {$a cmp $b} (   keys (  %{ $SgmTypeHash->{$idxKey} }  )   )    ){ 		  
  		    	my $lastIdx=$SgmTypeHash->{$idxKey}->{$sgmTpKy};                         print "7-20180906 \$SgmTypeHash->{$idxKey}=$SgmTypeHash->{$idxKey}\t\$lastIdx=\$SgmTypeHash->{\$idxKey}->{\$sgmTpKy}=\$SgmTypeHash->{$idxKey}->{$sgmTpKy}=$lastIdx\n";
  		    	$outPutArray->[$outPutArray_idx]                =Storable::dclone ($pointArray->[$idxKey]);
  		    	$outPutArray->[$outPutArray_idx]->{'0_SgmtHead'}=$pointArray->[ $idxKey  ]->{'9_del_SpPs'};
  		    	$outPutArray->[$outPutArray_idx]->{'1_SgmtTail'}=$pointArray->[ $lastIdx ]->{'9_del_SpPs'};
  		    	delete $outPutArray->[$outPutArray_idx]->{'8_del_StPs'} if (  defined ( $outPutArray->[$outPutArray_idx]->{'8_del_StPs'} )  );
  		    	delete $outPutArray->[$outPutArray_idx]->{'9_del_SpPs'} if (  defined ( $outPutArray->[$outPutArray_idx]->{'9_del_SpPs'} )  );
  		    	my $outPutArray_idx_dump= DirFileHandle::ReturnDumperInform($outPutArray->[$outPutArray_idx]); print "9-20180906 dump \$outPutArray->[$outPutArray_idx]=$outPutArray->[$outPutArray_idx]:\n$outPutArray_idx_dump\n";
  		      $outPutArray_idx++;
  		    }
  		  }
  		  
  		}  	 	
  	}
  	my $outPutArray_dump= DirFileHandle::ReturnDumperInform($outPutArray); print "10-20180906 dump \$outPutArray=$outPutArray:\n$outPutArray_dump\n";
  	return $outPutArray;
  	
  }
  else {
  	my $dieMsg="$die_MsgHead\n\$InGeneStructureArray=$InGeneStructureArray should be a ARRAY ref!!!\n\n\n";
  	print $dieMsg; die $dieMsg;
  }
}


#Merge_Est_proSplign_Splign_cap3_Inform($spling_cap3_merge_Hash,$bestProSplgnAlnHash);#
sub Merge_proSplign_Splign_cap3_Inform{  #  SplignHandle::Merge_Est_proSplign_Splign_cap3_Inform
  my ($In_spling_cap3_merge_Hash, $In_bestProSplgnAlnHash, $EstOrGenoType)=@_;
  
  my $warnMsgBody="\nIn package  SplignHandle,\tIn sub Merge_Est_proSplign_Splign_cap3_Inform,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	
  my $splGno_0='1_0splGno';                            my $sj_char_4='4sj_char';               #'3_sj_line'                         
	my $splExN_1='1_1splExN';                            my $Ex_numb_6='6Ex_numb';               #'4_ex_numb'                         
	my $spl_CS_2='1_2spl_CS';                            my $CS_char_3='3CS_char';               #'2_cs_line'                         
	my $splCtg_3='1_3splCtg';                            my $qr_char_2='2qr_char';               #'1_qr_line'                         
	my $spl_AA_4='1_4spl_AA';                            my $AA_char_1='1AA_char';               #'0_AA_line'                         
	my $ctgPos_a='Not_include_spl_ctg_posi';             my $qr_posi_0='0qr_posi';               #'Not_include_0qr_posi'
	my $gnoPos_b='Not_include_spl_gno_posi';             my $sj_posi_5='5sj_posi';               #'Not_include_5sj_posi'   
	my $ctgZoF_a='Not_include_spl_ctg_ZorF';             my $qr_ZorF_0='0qr_ZorF';               #  
	my $gnoZoF_b='Not_include_spl_gno_ZorF';             my $sj_ZorF_5='5sj_ZorF';               #
	
	
	my $whoLine_hash_Key;
	my ($splPosiKey, $splCharKey);
	my ($alinAA_key,         $PepFrm_key,         $trasAA_key,         $DNAchr_key        );
	my ($realAln_alinAA_key, $realAln_PepFrm_key, $realAln_trasAA_key, $realAln_DNAchr_key);
	if    ($EstOrGenoType eq 'Est')         { $whoLine_hash_Key='4_whoL_Hsh'; ($splPosiKey, $splCharKey)=($qr_posi_0, $qr_char_2); ($alinAA_key, $PepFrm_key, $trasAA_key, $DNAchr_key)=('00_est0alinAA', '00_est1PepFrm', '00_est2trasAA', '00_est3DNAchr'); ($realAln_alinAA_key, $realAln_PepFrm_key, $realAln_trasAA_key, $realAln_DNAchr_key)=('2_4e_alAA', '2_3e_ppFm', '2_1e_trAA', '2_0e_ACTG'); }
	elsif ($EstOrGenoType eq 'OneFiledGeno'){ $whoLine_hash_Key='1_whoL_Hsh'; ($splPosiKey, $splCharKey)=($sj_posi_5, $sj_char_4); ($alinAA_key, $PepFrm_key, $trasAA_key, $DNAchr_key)=('77_gno3alinAA', '77_gno2PepFrm', '77_gno1trasAA', '77_gno0DNAchr'); ($realAln_alinAA_key, $realAln_PepFrm_key, $realAln_trasAA_key, $realAln_DNAchr_key)=('0_2g_alAA', '0_3g_ppFm', '0_5g_trAA', '0_6g_ACTG'); }
	else {
	  my $dieMsg="$die_MsgHead\nThe \$EstOrGenoType=$EstOrGenoType should be a Est or OneFiledGeno !!!!\n\n";
  	print $dieMsg; die $dieMsg;
	}
	
  my $newOutHash;
  if  (   ( ref ($In_spling_cap3_merge_Hash) eq 'HASH' ) && (  ref ( $In_spling_cap3_merge_Hash->{$whoLine_hash_Key} ) eq 'HASH'  )   ){
  	my $newWoLiHash;
  	foreach my $whoLKey (    sort {$a <=> $b} (   keys (  %{ $In_spling_cap3_merge_Hash->{$whoLine_hash_Key} }  )   )    ){
  		my $tempHASH= $In_spling_cap3_merge_Hash->{$whoLine_hash_Key}->{$whoLKey};
  		
  		my $spl_dna_pos= $In_spling_cap3_merge_Hash->{$whoLine_hash_Key}->{$whoLKey}->{$splPosiKey};   ####################
  		my $spl_dna_chr= $In_spling_cap3_merge_Hash->{$whoLine_hash_Key}->{$whoLKey}->{$splCharKey};   ####################
  		if ( (defined $spl_dna_pos) && (defined $spl_dna_chr) && ($spl_dna_pos=~m/^\d+$/) && ($spl_dna_chr=~m/^[a-zA-Z]$/) ) {
  			my $DNAchr=' '; $DNAchr=$In_bestProSplgnAlnHash->{'7_DNA_to_PEP_hash'}->{$spl_dna_pos}->{'DNAChar'} if (   (  defined ( $In_bestProSplgnAlnHash->{'7_DNA_to_PEP_hash'}->{$spl_dna_pos}->{'DNAChar'} )  ) && (  $In_bestProSplgnAlnHash->{'7_DNA_to_PEP_hash'}->{$spl_dna_pos}->{'DNAChar'} =~m/^\S+$/)   );
  			my $PepFrm=' '; $PepFrm=$In_bestProSplgnAlnHash->{'7_DNA_to_PEP_hash'}->{$spl_dna_pos}->{'PepFrm' } if (   (  defined ( $In_bestProSplgnAlnHash->{'7_DNA_to_PEP_hash'}->{$spl_dna_pos}->{'PepFrm' } )  ) && (  $In_bestProSplgnAlnHash->{'7_DNA_to_PEP_hash'}->{$spl_dna_pos}->{'PepFrm' } =~m/^\d+$/)   );
  			my $trasAA=' '; $trasAA=$In_bestProSplgnAlnHash->{'7_DNA_to_PEP_hash'}->{$spl_dna_pos}->{'trasAA' } if (   (  defined ( $In_bestProSplgnAlnHash->{'7_DNA_to_PEP_hash'}->{$spl_dna_pos}->{'trasAA' } )  ) && (  $In_bestProSplgnAlnHash->{'7_DNA_to_PEP_hash'}->{$spl_dna_pos}->{'trasAA' } =~m/^\S+$/)   );
  			my $alinAA=' '; $alinAA=$In_bestProSplgnAlnHash->{'7_PEP_to_DNA_hash'}->{ $In_bestProSplgnAlnHash->{'7_DNA_to_PEP_hash'}->{$spl_dna_pos}->{'PepCorPos'} }->{'PePChar'}  if (   (  defined ( $In_bestProSplgnAlnHash->{'7_DNA_to_PEP_hash'}->{$spl_dna_pos}->{'PepCorPos'}  )  ) && (  $In_bestProSplgnAlnHash->{'7_DNA_to_PEP_hash'}->{$spl_dna_pos}->{'PepCorPos'} =~m/^\d+$/)  && (  defined ( $In_bestProSplgnAlnHash->{'7_PEP_to_DNA_hash'}->{ $In_bestProSplgnAlnHash->{'7_DNA_to_PEP_hash'}->{$spl_dna_pos}->{'PepCorPos'} }->{'PePChar'}  )  ) && (  $In_bestProSplgnAlnHash->{'7_PEP_to_DNA_hash'}->{ $In_bestProSplgnAlnHash->{'7_DNA_to_PEP_hash'}->{$spl_dna_pos}->{'PepCorPos'} }->{'PePChar'} =~m/^\S+$/)  );
  			if (   ($PepFrm=~m/^\d$/) && (  ($PepFrm == 1) || ($PepFrm == 3)  )   ){ $trasAA=' '; $alinAA=' '; }  # $trasAA=lc $trasAA; $alinAA=lc $alinAA;  
  			if (   ($PepFrm=~m/^\d$/) && (  ($PepFrm == 2) || ($PepFrm == 3)  )   ){ $DNAchr=lc $DNAchr;       }  # $trasAA=lc $trasAA; $alinAA=lc $alinAA;  
  			$tempHASH->{$DNAchr_key}=$DNAchr;
  			$tempHASH->{$alinAA_key}=$alinAA;                                              ################
  			$tempHASH->{$PepFrm_key}=$PepFrm;
  			$tempHASH->{$trasAA_key}=$trasAA;
  			
  		}
  		else {
  			$tempHASH->{$DNAchr_key}=' ';
  			$tempHASH->{$alinAA_key}=' ';
  			$tempHASH->{$PepFrm_key}=' ';
  			$tempHASH->{$trasAA_key}=' ';
  			
  		}
  	  
  	  $newWoLiHash->{$whoLKey}=$tempHASH;
  	  
  	}
  	$newOutHash->{'1_whoL_Hsh'}=$newWoLiHash;
  	
  	my $alnShwHash; my $newWlHash;
  	if ( ref ($newWoLiHash) eq 'HASH' ){
		  foreach my $newWLposKey (    sort {$a<=>$b} (   keys (  %{ $newWoLiHash }  )   )    ) {
		  	if (  ref ( $newWoLiHash->{$newWLposKey} ) eq 'HASH'  ){
		      foreach my $wordKey (    sort {$a cmp $b} (   keys (  %{ $newWoLiHash->{$newWLposKey} }  )   )    ) {
		      	my $realWordKey;
		      	if ( ($wordKey eq $qr_posi_0) || ($wordKey eq $sj_posi_5) ){ $newWlHash->{$newWLposKey}->{"Not_include_$wordKey"}=$newWoLiHash->{$newWLposKey}->{$wordKey};	}  #this two value is not used in the alignment shown hash
		      	else {   
		      		
		      	  if    ($wordKey eq '00_est0alinAA') {    $realWordKey='2_4e_alAA';       }         #$alnShwHash->{'2_2e_trAA'}.=$newWoLiHash->{$newWLposKey}->{$wordKey};              }
		      	  elsif ($wordKey eq '00_est1PepFrm') {    $realWordKey='2_3e_ppFm';       }         #$alnShwHash->{'2_1e_ppFm'}.=$newWoLiHash->{$newWLposKey}->{$wordKey};              }
		      	  elsif ($wordKey eq '00_est2trasAA') {    $realWordKey='2_4e_alAA';       }         #$alnShwHash->{'2_0e_alAA'}.=$newWoLiHash->{$newWLposKey}->{$wordKey};              }
		      	  elsif ($wordKey eq '00_est3DNAchr') {    $realWordKey='2_0e_ACTG';       }         #$alnShwHash->{'2_0e_alAA'}.=$newWoLiHash->{$newWLposKey}->{$wordKey};              }
		      	  elsif ($wordKey eq '77_gno3alinAA') {    $realWordKey='0_2g_alAA';       }         #$alnShwHash->{'2_2e_trAA'}.=$newWoLiHash->{$newWLposKey}->{$wordKey};              }
		      	  elsif ($wordKey eq '77_gno2PepFrm') {    $realWordKey='0_3g_ppFm';       }         #$alnShwHash->{'2_2e_trAA'}.=$newWoLiHash->{$newWLposKey}->{$wordKey};              }
		      	  elsif ($wordKey eq '77_gno1trasAA') {    $realWordKey='0_5g_trAA';       }         #$alnShwHash->{'2_1e_ppFm'}.=$newWoLiHash->{$newWLposKey}->{$wordKey};              }
		      	  elsif ($wordKey eq '77_gno0DNAchr') {    $realWordKey='0_6g_ACTG';       }         #$alnShwHash->{'2_0e_alAA'}.=$newWoLiHash->{$newWLposKey}->{$wordKey};              }
		      	  elsif ($wordKey eq $Ex_numb_6     ) {    $realWordKey=$splExN_1  ;       }         #$alnShwHash->{'1_0splExN'}.=$newWoLiHash->{$newWLposKey}->{$wordKey};              }
		      	  elsif ($wordKey eq $sj_char_4     ) {    $realWordKey=$splGno_0  ;       }         #$alnShwHash->{'1_1splGno'}.=$newWoLiHash->{$newWLposKey}->{$wordKey};              }
		      	  elsif ($wordKey eq $CS_char_3     ) {    $realWordKey=$spl_CS_2  ;       }         #$alnShwHash->{'1_2spl_CS'}.=$newWoLiHash->{$newWLposKey}->{$wordKey};              }
		      	  elsif ($wordKey eq $qr_char_2     ) {    $realWordKey=$splCtg_3  ;       }         #$alnShwHash->{'1_3splCtg'}.=$newWoLiHash->{$newWLposKey}->{$wordKey};              }
		      	  elsif ($wordKey eq $AA_char_1     ) {    $realWordKey=$spl_AA_4  ;       }         #$alnShwHash->{'1_4spl_AA'}.=$newWoLiHash->{$newWLposKey}->{$wordKey};              }
		      	  else                                {    $realWordKey=$wordKey   ;       }
		      	  
		      	  $alnShwHash->{$realWordKey}.=$newWoLiHash->{$newWLposKey}->{$wordKey};
		      	  $newWlHash->{$newWLposKey}->{$realWordKey}=$newWoLiHash->{$newWLposKey}->{$wordKey};
		      	  
		      	  ###########
		      	  
		
		      	  
		      	  ###########
		      	  
		      	  
		      	  
		      	}
		      }
		    }
  	    
  	  }
  	  $newOutHash->{'0___alnShw'}=$alnShwHash; 
  	  $newOutHash->{'2_alnShwHs'}=$newWlHash; 
  	  
    }
  	
  	
  	#my $direcSortednewOutHash=$newOutHash;
  	
  	my $direcSortednewOutHash=&CheckDirection_and_add_UM_mark($newOutHash, $EstOrGenoType);
  	
  	return $direcSortednewOutHash;
  	
  }
  
  else {
  	my $In_spling_cap3_merge_HashDump=DirFileHandle::ReturnDumperInform($In_spling_cap3_merge_Hash);
  	my $dieMsg="$die_MsgHead\nThe \$In_spling_cap3_merge_Hash=$In_spling_cap3_merge_Hash should be a correct array ref\n\$In_spling_cap3_merge_HashDump=$In_spling_cap3_merge_HashDump\n\n\n";
  	print $dieMsg; die $dieMsg;
  }  
}

sub CheckDirection_and_add_UM_mark{  #  SplignHandle::CheckDirection_and_add_UM_mark
	my ($inHash, $genoOrEst)=@_;
	
  my $warnMsgBody="\nIn package  SplignHandle,\tIn sub CheckDirection_and_add_UM_mark,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";

  my $splGno_0='1_0splGno';                            my $sj_char_4='4sj_char';               #'3_sj_line'                         
	my $splExN_1='1_1splExN';                            my $Ex_numb_6='6Ex_numb';               #'4_ex_numb'                         
	my $spl_CS_2='1_2spl_CS';                            my $CS_char_3='3CS_char';               #'2_cs_line'                         
	my $splCtg_3='1_3splCtg';                            my $qr_char_2='2qr_char';               #'1_qr_line'                         
	my $spl_AA_4='1_4spl_AA';                            my $AA_char_1='1AA_char';               #'0_AA_line'                         
	my $ctgPos_a='Not_include_spl_ctg_posi';             my $qr_posi_0='0qr_posi';               #'Not_include_0qr_posi'
	my $gnoPos_b='Not_include_spl_gno_posi';             my $sj_posi_5='5sj_posi';               #'Not_include_5sj_posi'
	my $ctgZoF_a='Not_include_spl_ctg_ZorF';             my $qr_ZorF_0='0qr_ZorF';               #  
	my $gnoZoF_b='Not_include_spl_gno_ZorF';             my $sj_ZorF_5='5sj_ZorF';               #

	my $chckFrmKey;  #my $ppFm_3; if ($GenoOrEst eq 'Est') { $ppFm_3='2_3e_ppFm'; }                elsif ($GenoOrEst eq 'OneFiledGeno') { $ppFm_3='0_3g_ppFm'; }
	if    (  ( defined ($genoOrEst) ) && ($genoOrEst eq 'Est'         )  ) {      $chckFrmKey='2_3e_ppFm';                        }
  elsif (  ( defined ($genoOrEst) ) && ($genoOrEst eq 'OneFiledGeno')  ) {      $chckFrmKey='0_3g_ppFm';                        }	
	else {
		my $dieMsg="$die_MsgHead\$genoOrEst=$genoOrEst should be defined and it should be OneFiledGeno or Est";
		print $dieMsg; die $dieMsg;
	}
	
	

  my $outHash;
	if ( ref ($inHash) eq 'HASH'){
		if (  ref ( $inHash->{'2_alnShwHs'} ) eq 'HASH'  ){ 
			
			my @sorted_a_b_keys_array=sort { $a<=>$b } (   keys (  %{ $inHash->{'2_alnShwHs'} }  )   );
			my @sorted_b_a_keys_array=sort { $b<=>$a } (   keys (  %{ $inHash->{'2_alnShwHs'} }  )   );
			
			my $lastFrmNB;  my $lastDirect;  my $Z_nb=0; my $F_nb=0;
		  foreach my $wlK ( @sorted_a_b_keys_array ) {

		  	my $frmNB; my $direct;  
		  	
		  	$frmNB=$inHash->{'2_alnShwHs'}->{$wlK}->{$chckFrmKey} if (   (  defined ( $inHash->{'2_alnShwHs'}->{$wlK}->{$chckFrmKey} )  ) && ( $inHash->{'2_alnShwHs'}->{$wlK}->{$chckFrmKey}=~m/\d+/ )   ); #print "+00_20180911 \$frmNB=$frmNB \$inHash->{'2_alnShwHs'}->{$wlK}->{$chckFrmKey}=$inHash->{'2_alnShwHs'}->{$wlK}->{$chckFrmKey}\n";
		  	#$frmNB=$inHash->{'2_alnShwHs'}->{$wlK}->{'0_1g_ppFm'} if (   (  defined ( $inHash->{'2_alnShwHs'}->{$wlK}->{'0_1g_ppFm'} )  ) && ( $inHash->{'2_alnShwHs'}->{$wlK}->{'0_1g_ppFm'}=~m/\d+/ )   );
		  	if (    (   (  defined ( $lastFrmNB )  ) && ( $lastFrmNB=~m/\d+/ )   ) && (   (  defined ( $frmNB )  ) && ( $frmNB=~m/\d+/ )   )    ) {
		  		if    (  ( ($lastFrmNB==1) && ($frmNB==2) ) || ( ($lastFrmNB==2) && ($frmNB==3) ) || ( ($lastFrmNB==3) && ($frmNB==1) )  ){ $direct='+';  $Z_nb++;  print "+00_20180911 \$lastFrmNB=$lastFrmNB \$frmNB=$frmNB \$direct=$direct\n"; }
		  		elsif (  ( ($lastFrmNB==3) && ($frmNB==2) ) || ( ($lastFrmNB==2) && ($frmNB==1) ) || ( ($lastFrmNB==1) && ($frmNB==3) )  ){ $direct='-';  $F_nb++;  print "-00_20180911 \$lastFrmNB=$lastFrmNB \$frmNB=$frmNB \$direct=$direct\n"; }
		  		else  {
		  			print "$warnMsgHead\n\$lastFrmNB=$lastFrmNB\t\$frmNB=$frmNB\tTHese two should be  12 23 31 or 32 21 13!!!\n\$wlK=$wlK\n\n\n";
		  			#my $dieMsg="$die_MsgHead\n\$lastFrmNB=$lastFrmNB\t\$frmNB=$frmNB\tTHese two should be  12 23 31 or 32 21 13!!!\n\$wlK=$wlK\n\n\n";
		  			#print $dieMsg; die $dieMsg;
		  		}
		  	}
		  	
		  	#if (    (   (  defined ( $direct )  ) && (  ( $direct eq '+' )  || ( $direct eq '-' )   )   ) && (   (  defined ( $lastDirect )  ) && (  ( $lastDirect eq '+' )  || ( $lastDirect eq '-' )   )   )    ) {
		  	#	if ($lastDirect ne $direct){
		  	#	  my $dieMsg="$die_MsgHead\n\$lastDirect=$lastDirect\t\$direct=$direct\tTHese two should be  ++ or --!!!\n\$wlK=$wlK\n\n\n";
		  	#		print $dieMsg; die $dieMsg;	
		  	#	}
		  	#}
		  	
		  	#print "00_20180911\$direct=$direct\t\$frmNB=\$inHash->{'2_alnShwHs'}->{$wlK}->{$chckFrmKey}=$frmNB\n";
		  	
		  	
		  	$lastFrmNB =$frmNB                                        if (   (  defined ( $frmNB )  ) && ( $frmNB=~m/\d+/ )   );
		  	#$lastDirect=$direct                                       if (   (  defined ( $direct )  ) && (  ( $direct eq '+' )  || ( $direct eq '-' )   )   );
		  }
		  print "0a_20180911\$Z_nb=$Z_nb\t\$F_nb=$F_nb\t\$lastFrmNB=$lastFrmNB\n";
		  
		  if    ( ( $Z_nb >0 ) && ( $F_nb==0 ) ) { $lastDirect = '+'; }
		  elsif ( ( $F_nb> 0 ) && ( $Z_nb==0 ) ) { $lastDirect = '-'; }
		  elsif ( ( $F_nb> 0 ) && ( $Z_nb> 0 ) ) { 
		    if ( $F_nb > $Z_nb ){
		    	if ($Z_nb <= 30){ $lastDirect = '-'; print "$warnMsgBody \$lastDirect=$lastDirect\t\$F_nb=$F_nb\t\$Z_nb=$Z_nb\n\n";              }
		    	else            { my $dieMsg="$die_MsgHead 11111 \$F_nb=$F_nb\t\$Z_nb=$Z_nb\n\n"; print $dieMsg; die $dieMsg;                     } 
		    }
		    if ( $F_nb < $Z_nb ){
		    	if ($F_nb <= 30){ $lastDirect = '+'; print "$warnMsgBody \$lastDirect=$lastDirect\t\$F_nb=$F_nb\t\$Z_nb=$Z_nb\n\n";              }
		    	else            { my $dieMsg="$die_MsgHead 22222 \$F_nb=$F_nb\t\$Z_nb=$Z_nb\n\n"; print $dieMsg; die $dieMsg;                     } 
		    }
		  }
		  elsif ( ( $F_nb==0 ) && ( $Z_nb==0 ) ) {
		    my $dieMsg="$die_MsgHead \$F_nb=$F_nb\t\$Z_nb=$Z_nb SHould not both be 0 \n\n"; print $dieMsg; die $dieMsg;
		  }
		  
		  print "0b_20180911\$lastDirect=$lastDirect\n";
		  
		  if (   (  defined ( $lastDirect )  ) && (  ( $lastDirect eq '+' )  || ( $lastDirect eq '-' )   )   ){
		  	my @Pep_a_b_keys_array;
		  	if    ( $lastDirect eq '-' ){ @Pep_a_b_keys_array=@sorted_b_a_keys_array; }
		  	elsif ( $lastDirect eq '+' ){ @Pep_a_b_keys_array=@sorted_a_b_keys_array; } 
		  	
		  	my $keyNb=0;	
		      foreach my $wlK ( @Pep_a_b_keys_array ){                                                                                                                                                                                            #print "10a_20180911\$wlK=$wlK\n";
		    	$outHash->{'1_whoL_Hsh'}->{$keyNb}=Storable::dclone( $inHash->{'1_whoL_Hsh'}->{$wlK} ) if (   (  defined ( $inHash->{'1_whoL_Hsh'}->{$wlK} )  ) && (  ref ( $inHash->{'1_whoL_Hsh'}->{$wlK} ) eq 'HASH'  )   ); 
		    	$outHash->{'2_alnShwHs'}->{$keyNb}=Storable::dclone( $inHash->{'2_alnShwHs'}->{$wlK} ) if (   (  defined ( $inHash->{'2_alnShwHs'}->{$wlK} )  ) && (  ref ( $inHash->{'2_alnShwHs'}->{$wlK} ) eq 'HASH'  )   ); 		
		    	#$outHash->{'2_alnShwHs'}->{$keyNb}=&add_U_M_mark($outHash->{'2_alnShwHs'}->{$keyNb}, '0_0g_alAA', '0_0g_0_UM', 'U' );          #print "11_20180911\$outHash->{'2_alnShwHs'}->{$keyNb}->{'0_0g_0_UM'}=$outHash->{'2_alnShwHs'}->{$keyNb}->{'0_0g_0_UM'}\n" if (  defined ( $outHash->{'2_alnShwHs'}->{$keyNb}->{'0_0g_0_UM'} )  );
		    	#$outHash->{'2_alnShwHs'}->{$keyNb}=&add_U_M_mark($outHash->{'2_alnShwHs'}->{$keyNb}, '0_2g_trAA', '0_2g_UMtA', 'U', 'M'  );    ##print "12_20180911\$outHash->{'2_alnShwHs'}->{$keyNb}->{'0_2g_UMtA'}=$outHash->{'2_alnShwHs'}->{$keyNb}->{'0_2g_UMtA'}\n" if (  defined ( $outHash->{'2_alnShwHs'}->{$keyNb}->{'0_2g_UMtA'} )  );
		    	#$outHash->{'2_alnShwHs'}->{$keyNb}=&add_U_M_mark($outHash->{'2_alnShwHs'}->{$keyNb}, '2_1e_trAA', '2_1e_umTA', 'U', 'M'  );    #print "13_20180911\$outHash->{'2_alnShwHs'}->{$keyNb}->{'2_1e_umTA'}=$outHash->{'2_alnShwHs'}->{$keyNb}->{'2_1e_umTA'}\n" if (  defined ( $outHash->{'2_alnShwHs'}->{$keyNb}->{'2_1e_umTA'} )  );
		    	#$outHash->{'2_alnShwHs'}->{$keyNb}=&add_U_M_mark($outHash->{'2_alnShwHs'}->{$keyNb}, '2_3e_alAA', '2_3e_uAlA', 'U'  );         #print "14_20180911\$outHash->{'2_alnShwHs'}->{$keyNb}->{'2_3e_uAlA'}=$outHash->{'2_alnShwHs'}->{$keyNb}->{'2_3e_uAlA'}\n" if (  defined ( $outHash->{'2_alnShwHs'}->{$keyNb}->{'2_3e_uAlA'} )  );
		    	$keyNb++;  	
		  	}
		  		
		  	
		  	if (   (  defined ( $inHash->{'0___alnShw'} )  ) && (  ref ( $inHash->{'0___alnShw'} ) eq 'HASH'  )   ){
		  		if (   (  defined ( $outHash->{'2_alnShwHs'} )  ) && (  ref ( $outHash->{'2_alnShwHs'} ) eq 'HASH'  )   ){
		  			foreach my $WlKei (    sort { $a <=>$b } (   keys (  %{ $outHash->{'2_alnShwHs'} }  )   )    ){
		  				if  (   (  defined ( $outHash->{'2_alnShwHs'}->{$WlKei} )  ) && (  ref ( $outHash->{'2_alnShwHs'}->{$WlKei} ) eq 'HASH'  )   ){
		  				  foreach my $prtKey (    sort { $a cmp $b } (   keys (  %{ $outHash->{'2_alnShwHs'}->{$WlKei} }  )   )    ){
		  				  	#if    ( ($prtKey eq $ctgPos_a) || ($prtKey eq $gnoPos_b) || ($prtKey eq 'Not_include_0Dn_posi')|| ($prtKey eq 'Not_include_5Pe_posi') ){	  	}
		  				  	if    (  $prtKey=~m/Not_include/   ){	  	}
		  				  	else {
		  				  		$outHash->{'0___alnShw'}->{$prtKey}.=$outHash->{'2_alnShwHs'}->{$WlKei}->{$prtKey};
		  				  	}
		  					}
		  				} 
		  			}
		  		}
		  		   
		  	}
		  	
		  	
		  }
		  
		  
		}
	}
	
	if ( ref ($outHash) eq 'HASH' ){
		return $outHash;
	}
	else {
		my $dieMsg="$die_MsgHead\n\$outHash=$outHash\tIt should be a HASH ref!!\n\$inHash=$inHash\n\n\n";
		print $dieMsg; die $dieMsg;	
	}
	
}


sub add_U_M_mark{  #  SplignHandle::add_U_M_mark   #($outHash->{'2_alnShwHs'}->{$keyNb}, '0_0g_alAA','0_0g_0_UM', 'U', 'M' )
	my ($inHash, $AA_wdKey, $ShowAAky, $AA_toSw1, $AA_toSw2)=@_;
	
	my $warnMsgBody="\nIn package  SplignHandle,\tIn sub add_U_M_mark,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";

  if (   (  ( defined ($AA_toSw1) ) && ($AA_toSw1=~m/^\S$/)  ) || (  ( defined ($AA_toSw2) ) && ($AA_toSw2=~m/^\S$/)  )   ){  #print "6 2018.09.11 \$AA_toSw2=$AA_toSw2\n";                                                  
  	if (  defined ( $inHash->{$AA_wdKey} )  ){            #print "7_2018.09.11 \$inHash->{$AA_wdKey}=$inHash->{$AA_wdKey}\n";                        
  		                         
  		if    (   (  ( defined ($AA_toSw1) ) && ($AA_toSw1=~m/^\S$/)  ) && ( $inHash->{$AA_wdKey} eq $AA_toSw1 )   ) {           #print "9_2018.09.11 \$inHash->{$AA_wdKey}=$inHash->{$AA_wdKey} eq \$AA_toSw2=$AA_toSw2\n";
  			$inHash->{$ShowAAky}=$AA_toSw1;                     #print "10_2018.09.11 \$inHash->{$ShowAAky}=$inHash->{$ShowAAky}=$AA_toSw2=\$AA_toSw2\n";       
  		}
  		elsif (   (  ( defined ($AA_toSw2) ) && ($AA_toSw2=~m/^\S$/)  ) && ( $inHash->{$AA_wdKey} eq $AA_toSw2 )   ) {           #print "9_2018.09.11 \$inHash->{$AA_wdKey}=$inHash->{$AA_wdKey} eq \$AA_toSw2=$AA_toSw2\n";
  			$inHash->{$ShowAAky}=$AA_toSw2;                     #print "10_2018.09.11 \$inHash->{$ShowAAky}=$inHash->{$ShowAAky}=$AA_toSw2=\$AA_toSw2\n";       
  		}
  		else {
  			$inHash->{$ShowAAky}=' ';                             #print "8_2018.09.11 \$inHash->{$ShowAAky}=$inHash->{$ShowAAky}\n";    
  		}
  	}
  }  
	return $inHash;                 #print "10_2018.09.11 \$inHash->{$ShowAAky}=$inHash->{$ShowAAky}=$AA_toSw2=\$AA_toSw2\n";
}

sub Merge_ExpandCtg { #  SplignHandle::Merge_ExpandCtg
	my ($inHash, $ExpHash)=@_;
	
  my $warnMsgBody="\nIn package  SplignHandle,\tIn sub Merge_ExpandCtg,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	
  my $splGno_0='1_0splGno';                            my $sj_char_4='4sj_char';               #'3_sj_line'                         
	my $splExN_1='1_1splExN';                            my $Ex_numb_6='6Ex_numb';               #'4_ex_numb'                         
	my $spl_CS_2='1_2spl_CS';                            my $CS_char_3='3CS_char';               #'2_cs_line'                         
	my $splCtg_3='1_3splCtg';                            my $qr_char_2='2qr_char';               #'1_qr_line'                         
	my $spl_AA_4='1_4spl_AA';                            my $AA_char_1='1AA_char';               #'0_AA_line'                         
	my $ctgPos_a='Not_include_spl_ctg_posi';             my $qr_posi_0='0qr_posi';               #'Not_include_0qr_posi'
	my $gnoPos_b='Not_include_spl_gno_posi';             my $sj_posi_5='5sj_posi';               #'Not_include_5sj_posi'
  my $ctgZoF_a='Not_include_spl_ctg_ZorF';             my $qr_ZorF_0='0qr_ZorF';               #  
	my $gnoZoF_b='Not_include_spl_gno_ZorF';             my $sj_ZorF_5='5sj_ZorF';               #


	
	my $outHash;
	if (  ( ref ($inHash) eq 'HASH') && ( ref ($ExpHash) eq 'HASH')  ){
		if (     (  ref ( $inHash->{'2_alnShwHs'} ) eq 'HASH'  ) 
		      && (  defined ( $ExpHash->{$new_stt_key} )  ) &&  ( $ExpHash->{$new_stt_key}=~m/\d+/ ) 
		      && (  defined ( $ExpHash->{$new_end_key} )  ) &&  ( $ExpHash->{$new_end_key}=~m/\d+/ ) 
		      && (  defined ( $ExpHash->{$old_stt_key} )  ) &&  ( $ExpHash->{$old_stt_key}=~m/\d+/ ) 
		      && (  defined ( $ExpHash->{$old_end_key} )  ) &&  ( $ExpHash->{$old_end_key}=~m/\d+/ )   ){ 
			
			my @sorted_a_b_keys_array=sort { $a<=>$b } (   keys (  %{ $inHash->{'2_alnShwHs'} }  )   );
			#my @sorted_b_a_keys_array=sort { $b<=>$a } (   keys (  %{ $inHash->{'2_alnShwHs'} }  )   );
			
			
			my $pos_to_char_hash;
			foreach my $wlK ( @sorted_a_b_keys_array ) {
			  if (    
			           (   (  defined ( $inHash->{'2_alnShwHs'}->{$wlK}->{$ctgPos_a} )  ) && ( $inHash->{'2_alnShwHs'}->{$wlK}->{$ctgPos_a}=~m/\d+/ )   )
			       &&  (   (  defined ( $inHash->{'2_alnShwHs'}->{$wlK}->{$splCtg_3} )  ) && ( $inHash->{'2_alnShwHs'}->{$wlK}->{$splCtg_3}=~m/\S+/ )   )
			     ){
			  	$pos_to_char_hash->{ $inHash->{'2_alnShwHs'}->{$wlK}->{$ctgPos_a} }=$inHash->{'2_alnShwHs'}->{$wlK}->{$splCtg_3};
			  }
			}
			
			
			my $exp_stt=$ExpHash->{$new_stt_key};
			my $exp_end=$ExpHash->{$new_end_key};
			my $old_stt=$ExpHash->{$old_stt_key};
			my $old_end=$ExpHash->{$old_end_key};
			my $exp_ZoF=$ExpHash->{'2_SgmtZorF'};
			
			
			my @nbArray;			
			my ($Up___strm_hd, $Up___strm_tl, $Down_strm_hd, $Down_strm_tl); 
			
			if ($exp_ZoF eq '+'){ 
				if     ($exp_stt == $old_stt){ }
				elsif  ($exp_stt <  $old_stt){   ($Up___strm_hd, $Up___strm_tl)=($exp_stt, $old_stt-1);    }
				elsif  ($exp_stt >  $old_stt){my $dieMsg="$die_MsgHead\n\111 $exp_ZoF=$exp_ZoF\t\$exp_stt=$exp_stt should <= \$old_stt=$old_stt!!\n\n"; print $dieMsg; die $dieMsg; }
				else                         {my $dieMsg="$die_MsgHead\n\222 $exp_ZoF=$exp_ZoF\t\$exp_stt=$exp_stt should <= \$old_stt=$old_stt!!\n\n"; print $dieMsg; die $dieMsg; }
				
				if     ($exp_end == $old_end){ }
				elsif  ($exp_end >  $old_end){   ($Down_strm_hd, $Down_strm_tl)=($old_end+1, $exp_end);    }
				elsif  ($exp_end <  $old_end){my $dieMsg="$die_MsgHead\n\111 $exp_ZoF=$exp_ZoF\t\$exp_end=$exp_end should >= \$old_end=$old_end!!\n\n"; print $dieMsg; die $dieMsg; }
				else                         {my $dieMsg="$die_MsgHead\n\222 $exp_ZoF=$exp_ZoF\t\$exp_end=$exp_end should >= \$old_end=$old_end!!\n\n"; print $dieMsg; die $dieMsg; }
				
				if (   (  ( defined ($Up___strm_hd) ) && ($Up___strm_hd=~m/\d+/)  )  && (  ( defined ($Up___strm_tl) ) && ($Up___strm_tl=~m/\d+/)  )   ){
					for (my $i=$Up___strm_hd; $i<=$Up___strm_tl; $i++){
						push @nbArray, $i;                                                  #print "01-20180913$inHash, $ExpHash   {$exp_ZoF [$exp_stt ($old_stt, $old_end) $exp_end] } \$i=$i\t\@nbArray=@nbArray\n";
					} 
				}
				if (   (  ( defined ($Down_strm_hd) ) && ($Down_strm_hd=~m/\d+/)  )  && (  ( defined ($Down_strm_tl) ) && ($Down_strm_tl=~m/\d+/)  )   ){
					for (my $i=$Down_strm_hd; $i<=$Down_strm_tl; $i++){
						push @nbArray, $i;                                                  #print "02-20180913$inHash, $ExpHash   {$exp_ZoF [$exp_stt ($old_stt, $old_end) $exp_end] } \$i=$i\t\@nbArray=@nbArray\n";
					}
				}
				
				
				
			}	
			
			if ($exp_ZoF eq '-'){ 
				if     ($exp_stt == $old_stt){ }
				elsif  ($exp_stt >  $old_stt){   ($Up___strm_hd, $Up___strm_tl)=($exp_stt, $old_stt+1);    }
				elsif  ($exp_stt <  $old_stt){my $dieMsg="$die_MsgHead\n\333 $exp_ZoF=$exp_ZoF\t\$exp_stt=$exp_stt should >= \$old_stt=$old_stt!!\n\n"; print $dieMsg; die $dieMsg; }
				else                         {my $dieMsg="$die_MsgHead\n\444 $exp_ZoF=$exp_ZoF\t\$exp_stt=$exp_stt should >= \$old_stt=$old_stt!!\n\n"; print $dieMsg; die $dieMsg; }
				
				if     ($exp_end == $old_end){ }
				elsif  ($exp_end <  $old_end){   ($Down_strm_hd, $Down_strm_tl)=($old_end-1, $exp_end);    }
				elsif  ($exp_end >  $old_end){my $dieMsg="$die_MsgHead\n\333 $exp_ZoF=$exp_ZoF\t\$exp_end=$exp_end should <= \$old_end=$old_end!!\n\n"; print $dieMsg; die $dieMsg; }
				else                         {my $dieMsg="$die_MsgHead\n\444 $exp_ZoF=$exp_ZoF\t\$exp_end=$exp_end should <= \$old_end=$old_end!!\n\n"; print $dieMsg; die $dieMsg; }
				
				
				if (   (  ( defined ($Up___strm_hd) ) && ($Up___strm_hd=~m/\d+/)  )  && (  ( defined ($Up___strm_tl) ) && ($Up___strm_tl=~m/\d+/)  )   ){
					for (my $i=$Up___strm_hd; $i>=$Up___strm_tl; $i--){
						push @nbArray, $i;                                                  #print "03-20180913$inHash, $ExpHash   {$exp_ZoF [$exp_stt ($old_stt, $old_end) $exp_end] } \$i=$i\t\@nbArray=@nbArray\n";                                                          
					}
				}
				if (   (  ( defined ($Down_strm_hd) ) && ($Down_strm_hd=~m/\d+/)  )  && (  ( defined ($Down_strm_tl) ) && ($Down_strm_tl=~m/\d+/)  )   ){
					for (my $i=$Down_strm_hd; $i>=$Down_strm_tl; $i--){
						push @nbArray, $i;                                                  #print "04-20180913$inHash, $ExpHash   {$exp_ZoF [$exp_stt ($old_stt, $old_end) $exp_end] } \$i=$i\t\@nbArray=@nbArray\n";
					}
				}		
				
			}	
			
			my $exp_hash;
			my @three_nb_array; 
			for (my $j=0; $j<@nbArray; $j++){ my $stpNB=$j+1;                            #print "05-20180913$inHash, $ExpHash    \$j=$j\t\$stpNB=$stpNB\n";
				my $step_3_2_1 = $stpNB % 3; if ( $step_3_2_1 == 0 ){ $step_3_2_1=3; }     #print "06-20180913$inHash, $ExpHash    \$step_3_2_1=$step_3_2_1\n";
				my $posNB=$nbArray[$j];                                                    #print "07-20180913$inHash, $ExpHash    \$posNB=\$nbArray[$j]=$posNB\n";       
				
				$exp_hash->{$posNB}->{'4_0x_ACTG'}=$pos_to_char_hash->{$posNB} if (   (  defined ( $pos_to_char_hash->{$posNB} )  ) && ( $pos_to_char_hash->{$posNB} =~m/\S+/)   );   #print "08-20180913$inHash, $ExpHash    \$exp_hash->{$posNB}->{'4_0x_ACTG'}=$exp_hash->{$posNB}->{'4_0x_ACTG'}\n";   
				$exp_hash->{$posNB}->{'4_3x_ppFm'}=$step_3_2_1;	                                                                                                                      #print "09-20180913$inHash, $ExpHash    \$exp_hash->{$posNB}->{'4_3x_ppFm'}=$exp_hash->{$posNB}->{'4_3x_ppFm'}\n";  
				if    ($step_3_2_1 == 1){ $exp_hash->{$posNB}->{'4_0x_ACTG'}=uc ( $exp_hash->{$posNB}->{'4_0x_ACTG'} )  if (   (  defined ( $exp_hash->{$posNB}->{'4_0x_ACTG'} )  ) && ( $exp_hash->{$posNB}->{'4_0x_ACTG'}=~m/\S+/)   ); }
				elsif ($step_3_2_1 == 2){ $exp_hash->{$posNB}->{'4_0x_ACTG'}=lc ( $exp_hash->{$posNB}->{'4_0x_ACTG'} )  if (   (  defined ( $exp_hash->{$posNB}->{'4_0x_ACTG'} )  ) && ( $exp_hash->{$posNB}->{'4_0x_ACTG'}=~m/\S+/)   ); }
				elsif ($step_3_2_1 == 3){ $exp_hash->{$posNB}->{'4_0x_ACTG'}=lc ( $exp_hash->{$posNB}->{'4_0x_ACTG'} )  if (   (  defined ( $exp_hash->{$posNB}->{'4_0x_ACTG'} )  ) && ( $exp_hash->{$posNB}->{'4_0x_ACTG'}=~m/\S+/)   ); }
				
				if    ($step_3_2_1==1) { @three_nb_array=($posNB); }					
				else                   { push @three_nb_array, $posNB; }
				
				if (   ($step_3_2_1==3) || (  ( $j == (@nbArray-1) ) && ( ($step_3_2_1==1) || ($step_3_2_1==2) )  )   ){ 
					my $three_actg;
					foreach my $ecPnb ( @three_nb_array )  { 							$three_actg.=$pos_to_char_hash->{$ecPnb}		if (   (  defined ( $pos_to_char_hash->{$ecPnb} )  ) && ( $pos_to_char_hash->{$ecPnb}=~m/\S+/ )   );	 		}  #warn 	"\$pos_to_char_hash->{$ecPnb}=$pos_to_char_hash->{$ecPnb}\t\$three_actg.=\$pos_to_char_hash->{$ecPnb}=$three_actg\n";
					my $Pep___char=' '; if ($step_3_2_1==3){              $Pep___char=FastaFileHandle::CodonTable_standard ($three_actg); }
					my $nb=1;
					foreach my $ecPnb ( @three_nb_array )  { 							
						if ($nb==2) {
							
							$exp_hash->{$ecPnb}->{'4_1x_trAA'}=$Pep___char if (  ( defined ($Pep___char) ) && ($Pep___char=~m/\S+/)  ) ;                       #print "10-20180913$inHash, $ExpHash    \$exp_hash->{$posNB}->{'4_1x_trAA'}=$exp_hash->{$posNB}->{'4_1x_trAA'}\n";   
							$exp_hash->{$ecPnb}->{'4_2x_AAum'}=$Pep___char if (  ( defined ($Pep___char) ) && ($Pep___char=~m/(U|M|\*)/)  ) ;                     #print "11-20180913$inHash, $ExpHash    \$exp_hash->{$posNB}->{'4_2x_AAum'}=$exp_hash->{$posNB}->{'4_2x_AAum'}\n";   
						}						
						$nb++;
					}
					
				}
				#$exp_hash->{$posNB}->{'4_0x_ACTG'}=' ';
				#$exp_hash->{$posNB}->{'4_1x_trAA'}=' ';
				#$exp_hash->{$posNB}->{'4_2x_AAum'}=' ';
				#$exp_hash->{$posNB}->{'4_3x_ppFm'}
				
			}
			
			foreach my $wlK ( @sorted_a_b_keys_array ){
				$inHash->{'2_alnShwHs'}->{$wlK}->{'4_0x_ACTG'}=' ';
				$inHash->{'2_alnShwHs'}->{$wlK}->{'4_1x_trAA'}=' ';
				$inHash->{'2_alnShwHs'}->{$wlK}->{'4_2x_AAum'}=' ';
				$inHash->{'2_alnShwHs'}->{$wlK}->{'4_3x_ppFm'}=' ';
				
				if  (   (  defined ( $inHash->{'2_alnShwHs'}->{$wlK}->{$ctgPos_a} )  ) && ( $inHash->{'2_alnShwHs'}->{$wlK}->{$ctgPos_a}=~m/\d+/ )   ){          #print "a12-20180913$inHash, $ExpHash    \$inHash->{'2_alnShwHs'}->{\$wlK}->{$ctgPos_a}=\$inHash->{'2_alnShwHs'}->{$wlK}->{$ctgPos_a}=$inHash->{'2_alnShwHs'}->{$wlK}->{$ctgPos_a}\n";   
					my $ctgPos=$inHash->{'2_alnShwHs'}->{$wlK}->{$ctgPos_a};                                                                                                  #print "b12-20180913$inHash, $ExpHash    \$ctgPos=\$inHash->{'2_alnShwHs'}->{\$wlK}->{$ctgPos_a}=\$inHash->{'2_alnShwHs'}->{$wlK}->{$ctgPos_a}=$ctgPos\n";   
					$inHash->{'2_alnShwHs'}->{$wlK}->{'4_0x_ACTG'}=$exp_hash->{$ctgPos}->{'4_0x_ACTG'} if (   (  defined ( $exp_hash->{$ctgPos}->{'4_0x_ACTG'} )  ) && ( $exp_hash->{$ctgPos}->{'4_0x_ACTG'}=~m/\S+/ )   );  #print "13-20180913$inHash, $ExpHash    \$inHash->{'2_alnShwHs'}->{\$wlK}->{'4_0x_ACTG'}=\$exp_hash->{\$ctgPos}->{'4_0x_ACTG'}=\$inHash->{'2_alnShwHs'}->{$wlK}->{'4_0x_ACTG'}=\$exp_hash->{$ctgPos}->{'4_0x_ACTG'}=$inHash->{'2_alnShwHs'}->{$wlK}->{'4_0x_ACTG'}\n"; 
					$inHash->{'2_alnShwHs'}->{$wlK}->{'4_1x_trAA'}=$exp_hash->{$ctgPos}->{'4_1x_trAA'} if (   (  defined ( $exp_hash->{$ctgPos}->{'4_1x_trAA'} )  ) && ( $exp_hash->{$ctgPos}->{'4_1x_trAA'}=~m/\S+/ )   );  #print "14-20180913$inHash, $ExpHash    \$inHash->{'2_alnShwHs'}->{\$wlK}->{'4_1x_trAA'}=\$exp_hash->{\$ctgPos}->{'4_1x_trAA'}=\$inHash->{'2_alnShwHs'}->{$wlK}->{'4_1x_trAA'}=\$exp_hash->{$ctgPos}->{'4_1x_trAA'}=$inHash->{'2_alnShwHs'}->{$wlK}->{'4_1x_trAA'}\n"; 
					$inHash->{'2_alnShwHs'}->{$wlK}->{'4_2x_AAum'}=$exp_hash->{$ctgPos}->{'4_2x_AAum'} if (   (  defined ( $exp_hash->{$ctgPos}->{'4_2x_AAum'} )  ) && ( $exp_hash->{$ctgPos}->{'4_2x_AAum'}=~m/\S+/ )   );  #print "15-20180913$inHash, $ExpHash    \$inHash->{'2_alnShwHs'}->{\$wlK}->{'4_2x_AAum'}=\$exp_hash->{\$ctgPos}->{'4_2x_AAum'}=\$inHash->{'2_alnShwHs'}->{$wlK}->{'4_2x_AAum'}=\$exp_hash->{$ctgPos}->{'4_2x_AAum'}=$inHash->{'2_alnShwHs'}->{$wlK}->{'4_2x_AAum'}\n"; 
					$inHash->{'2_alnShwHs'}->{$wlK}->{'4_3x_ppFm'}=$exp_hash->{$ctgPos}->{'4_3x_ppFm'} if (   (  defined ( $exp_hash->{$ctgPos}->{'4_3x_ppFm'} )  ) && ( $exp_hash->{$ctgPos}->{'4_3x_ppFm'}=~m/\S+/ )   );  #print "16-20180913$inHash, $ExpHash    \$inHash->{'2_alnShwHs'}->{\$wlK}->{'4_3x_ppFm'}=\$exp_hash->{\$ctgPos}->{'4_3x_ppFm'}=\$inHash->{'2_alnShwHs'}->{$wlK}->{'4_3x_ppFm'}=\$exp_hash->{$ctgPos}->{'4_3x_ppFm'}=$inHash->{'2_alnShwHs'}->{$wlK}->{'4_3x_ppFm'}\n"; 
				}				  
			}
	    
	    
	    if (  defined ( $inHash->{'0___alnShw'} )  ) {
	      foreach my $prtKey (    sort { $a cmp $b } (   keys (  %{ $inHash->{'0___alnShw'} }  )   )    ){		    
		    	$inHash->{'0___alnShw'}->{$prtKey}='';
		    }
		  }				
	    
	    if (  defined ( $inHash->{'2_alnShwHs'} )  ) {
		  	if (   (  defined ( $inHash->{'2_alnShwHs'} )  ) && (  ref ( $inHash->{'2_alnShwHs'} ) eq 'HASH'  )   ){
		  		foreach my $WlKei (    sort { $a <=>$b } (   keys (  %{ $inHash->{'2_alnShwHs'} }  )   )    ){
		  			if  (   (  defined ( $inHash->{'2_alnShwHs'}->{$WlKei} )  ) && (  ref ( $inHash->{'2_alnShwHs'}->{$WlKei} ) eq 'HASH'  )   ){
		  			  foreach my $prtKey (    sort { $a cmp $b } (   keys (  %{ $inHash->{'2_alnShwHs'}->{$WlKei} }  )   )    ){
		  			  	#if    ( ($prtKey eq $ctgPos_a) || ($prtKey eq $gnoPos_b) || ($prtKey eq 'Not_include_0Dn_posi')|| ($prtKey eq 'Not_include_5Pe_posi') ){	  	}
		  			  	if    (  $prtKey=~m/Not_include/   ){	  	}
		  			  	else {
		  			  		$inHash->{'0___alnShw'}->{$prtKey}.=$inHash->{'2_alnShwHs'}->{$WlKei}->{$prtKey};
		  			  	}
		  				}
		  			} 
		  		}
		  	}		  	   
		  }
	    

		  
		  
		}
	}
	
	
  return $inHash;
	

	
}

sub Better_Merge_proSplign_Splign_cap3_Inform{     #  SplignHandle::Better_Merge_proSplign_Splign_cap3_Inform  #($proSplign_Splign_cap3_merge_hash, $bestGenoProsplignHash, 'OneFiledGeno');
	my ($in_aln_hash_Splig, $in_aln_Hash_ProSP, $Geno_or_Est)=@_;
	
	
  my $warnMsgBody="\nIn package  SplignHandle,\tIn sub Better_Merge_proSplign_Splign_cap3_Inform,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	
	#my $proDnaZ_or_F;
	#if (   ( defined ($in_aln_Hash_ProSP) ) && ( ref ($in_aln_Hash_ProSP) eq 'HASH' ) &&  (  defined ( $in_aln_Hash_ProSP->{'2_DNA_ZF'} )  ) && ( $in_aln_Hash_ProSP->{'2_DNA_ZF'}=~m/(-|\+)/ )   ){
	#	$proDnaZ_or_F=$in_aln_Hash_ProSP->{'2_DNA_ZF'};
	#}
	#else {
	#	my $dieMsg="$die_MsgHead\$in_aln_Hash_ProSP=$in_aln_Hash_ProSP should be defined and it should be a HASH ref\t\$in_aln_Hash_ProSP->{'2_DNA_ZF'}=$in_aln_Hash_ProSP->{'2_DNA_ZF'} should be defined and it should be + or - !!!!\n\n\n";
	#	print $dieMsg; die $dieMsg;
	#}
  

  
	
	my $ProSp_est_dna_posi='Not_include_e_Dna_posi';        my $ProSp_gno_dna_posi='Not_include_g_Dna_posi'; 
	my $ProSp_est_dna_ZorF='Not_include_e_Dna_ZorF';        my $ProSp_gno_dna_ZorF='Not_include_g_Dna_ZorF';	
	#my $ProSp_est_pep_posi='Not_include_e_Pep_posi';        my $ProSp_gno_pep_posi='Not_include_g_Pep_posi';
	#my $ProSp_est_pep_ZorF='Not_include_e_Pep_ZorF';       #my $ProSp_gno_pep_ZorF='Not_include_g_Pep_ZorF';
	
	
	my $Splig_est_ctg_posi='Not_include_spl_ctg_posi';      my $Splig_genome__posi='Not_include_spl_gno_posi';   
	my $Splig_est_ctg_ZorF='Not_include_spl_ctg_ZorF';      my $Splig_genome__ZorF='Not_include_spl_gno_ZorF';   
	  
	
	
	
	my $ProSp_ruler_posi;
	my $ProSp_ruler_ZorF;
	
	my $Splig_ruler_posi;
	my $Splig_ruler_ZorF;
	
	
	if    (  ( defined ($Geno_or_Est) ) && ($Geno_or_Est eq 'Est'         )  ) { $ProSp_ruler_posi=$ProSp_est_dna_posi; $ProSp_ruler_ZorF=$ProSp_est_dna_ZorF; $Splig_ruler_posi=$Splig_est_ctg_posi; $Splig_ruler_ZorF=$Splig_est_ctg_ZorF; }
  elsif (  ( defined ($Geno_or_Est) ) && ($Geno_or_Est eq 'OneFiledGeno')  ) { $ProSp_ruler_posi=$ProSp_gno_dna_posi; $ProSp_ruler_ZorF=$ProSp_gno_dna_ZorF; $Splig_ruler_posi=$Splig_genome__posi; $Splig_ruler_ZorF=$Splig_genome__ZorF; }	
	else {
		my $dieMsg="$die_MsgHead\$Geno_or_Est=$Geno_or_Est should be defined and it should be OneFiledGeno or Est\n\n\n";
		print $dieMsg; die $dieMsg;
	}
	
	my $outHash=&Merge_two_aln_hash($in_aln_hash_Splig, $Splig_ruler_posi, $Splig_ruler_ZorF, $in_aln_Hash_ProSP, $ProSp_ruler_posi, $ProSp_ruler_ZorF);
	
	my $direcSortednewOutHash=&CheckDirection_and_add_UM_mark($outHash, $Geno_or_Est);
  	
  return $direcSortednewOutHash;
	
}

sub Merge_two_aln_hash{   #SplignHandle::Merge_two_aln_hash  
	my ($in_aln_Hash_1, $ruler_posi_1, $ruler_ZorF_1, $in_aln_Hash_2, $ruler_posi_2, $ruler_ZorF_2)=@_;
	
	
  my $warnMsgBody="\nIn package  SplignHandle,\tIn sub Merge_two_aln_hash,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	
	my $outHash;
	
	if     (  ( ref ($in_aln_Hash_1) eq 'HASH' ) && ( ref ($in_aln_Hash_2) ne 'HASH' )  ){
	  $outHash=$in_aln_Hash_1;
	}	
	elsif  (  ( ref ($in_aln_Hash_1) ne 'HASH' ) && ( ref ($in_aln_Hash_2) eq 'HASH' )  ){    
	  $outHash=$in_aln_Hash_2;	
	}
	elsif  (  ( ref ($in_aln_Hash_1) eq 'HASH' ) && ( ref ($in_aln_Hash_2) eq 'HASH' )  ){    
	  
	  #my $Dierct_1=&CheckDirection($in_aln_Hash_1, $ruler_posi_1);  
	  #my $Dierct_2=&CheckDirection($in_aln_Hash_2, $ruler_posi_2);  
	  
	  #my @sortKeyArr_1=sort { $a <=> $b } (   keys (  %{ $in_aln_Hash_1->{$alnShwHs_key_1} }  )   ); if ($Dierct_1 eq '-') { @sortKeyArr_1= reverse @sortKeyArr_1; }
	  #my @sortKeyArr_2=sort { $a <=> $b } (   keys (  %{ $in_aln_Hash_2->{$alnShwHs_key_2} }  )   ); if ($Dierct_2 eq '-') { @sortKeyArr_2= reverse @sortKeyArr_2; }
	  
	  my $All_befor_part_KeyHASH_1=&GetAllRuler_befor_part_KeyHASH($in_aln_Hash_1, $ruler_posi_1);
	  my $All_befor_part_KeyHASH_2=&GetAllRuler_befor_part_KeyHASH($in_aln_Hash_2, $ruler_posi_2);
	  
	  my $All_after_part_KeyHASH_1=&GetAllRuler_after_part_KeyHASH($in_aln_Hash_1, $ruler_posi_1);
	  my $All_after_part_KeyHASH_2=&GetAllRuler_after_part_KeyHASH($in_aln_Hash_2, $ruler_posi_2);
	  
	  
	  
	  my $alnShwHs_key_1=&Get_alnShwHs_key($in_aln_Hash_1);
	  my $alnShwHs_key_2=&Get_alnShwHs_key($in_aln_Hash_2);
	  
	  my $All_key_hash_1=&GetAllKey_from_hash( $in_aln_Hash_1->{$alnShwHs_key_1} );  my @allKeyArray_1=keys ( %{ $All_key_hash_1 } );
	  my $All_key_hash_2=&GetAllKey_from_hash( $in_aln_Hash_2->{$alnShwHs_key_2} );  my @allKeyArray_2=keys ( %{ $All_key_hash_2 } );
	  my $All_Key_type; 
	  foreach my $key1 (  keys ( %{ $All_key_hash_1 } )  ){ $All_Key_type->{$key1}=1; }
	  foreach my $key2 (  keys ( %{ $All_key_hash_2 } )  ){ $All_Key_type->{$key2}=1; }
	  
	  my $allRulerKeyHash;
	  $allRulerKeyHash=&GetAllRulerKeyHASH($allRulerKeyHash, $in_aln_Hash_1, $ruler_posi_1); my $allRulerKeyHASH_only_1; $allRulerKeyHASH_only_1=&GetAllRulerKeyHASH($allRulerKeyHASH_only_1, $in_aln_Hash_1, $ruler_posi_1);
	  $allRulerKeyHash=&GetAllRulerKeyHASH($allRulerKeyHash, $in_aln_Hash_2, $ruler_posi_2); my $allRulerKeyHASH_only_2; $allRulerKeyHASH_only_2=&GetAllRulerKeyHASH($allRulerKeyHASH_only_2, $in_aln_Hash_2, $ruler_posi_2);
	  
	  
	  my @allRulerKeyArray=sort { $a <=> $b } (   keys (  %{ $allRulerKeyHash }  )   ); 
	  my @allRulerKeyHASH_only_1_Array=sort { $b <=> $a } (   keys (  %{ $allRulerKeyHASH_only_1 }  )   ); my $hash_1_the_last_key=$allRulerKeyHASH_only_1_Array[0];
	  my @allRulerKeyHASH_only_2_Array=sort { $b <=> $a } (   keys (  %{ $allRulerKeyHASH_only_2 }  )   ); my $hash_2_the_last_key=$allRulerKeyHASH_only_2_Array[0]; 
	  
	  
    my $stepHer=0;
	  foreach my $rulerKeyPos ( @allRulerKeyArray ){
	  	if  (   (  defined ( $All_befor_part_KeyHASH_1->{$rulerKeyPos} )  ) &&  (  ref ( $All_befor_part_KeyHASH_1->{$rulerKeyPos} ) eq 'HASH'  ) && (  defined ( $All_befor_part_KeyHASH_1->{$rulerKeyPos}->{'last_no_rular_wl_key_arry'} )  ) &&  (  ref ( $All_befor_part_KeyHASH_1->{$rulerKeyPos}->{'last_no_rular_wl_key_arry'} ) eq'ARRAY'  )  ){
        foreach my $wlK	 (  @{ $All_befor_part_KeyHASH_1->{$rulerKeyPos}->{'last_no_rular_wl_key_arry'} }  ){
          foreach my $allPtKey (   keys (  %{ $All_Key_type }  )   ){        	
        	  #if (   (  defined ( $outHash->{'2_alnShwHs'}->{$stepHer}->{$allPtKey} )  ) && (  ( $outHash->{'2_alnShwHs'}->{$stepHer}->{$allPtKey} eq ' ' ) || ( $outHash->{'2_alnShwHs'}->{$stepHer}->{$allPtKey}=~/\S+/ )  )   ) { 
        	  	
        	  #}
        	  #else {
        	  	$outHash->{'2_alnShwHs'}->{$stepHer}->{$allPtKey}=' ';  
        	    $outHash->{'2_alnShwHs'}->{$stepHer}->{$allPtKey}=$in_aln_Hash_1->{$alnShwHs_key_1}->{$wlK}->{$allPtKey} if (   (  defined ( $in_aln_Hash_1->{$alnShwHs_key_1}->{$wlK}->{$allPtKey} )  ) && ( $in_aln_Hash_1->{$alnShwHs_key_1}->{$wlK}->{$allPtKey}=~m/\S+/ )   );
        	  
        	  #}
        	  
        	}
        	$stepHer++;
        }
        
	  	}
	  	
	  	if  (   (  defined ( $All_befor_part_KeyHASH_2->{$rulerKeyPos} )  ) &&  (  ref ( $All_befor_part_KeyHASH_2->{$rulerKeyPos} ) eq 'HASH'  ) && (  defined ( $All_befor_part_KeyHASH_2->{$rulerKeyPos}->{'last_no_rular_wl_key_arry'} )  ) &&  (  ref ( $All_befor_part_KeyHASH_2->{$rulerKeyPos}->{'last_no_rular_wl_key_arry'} ) eq'ARRAY'  )  ){
        foreach my $wlK	 (  @{ $All_befor_part_KeyHASH_2->{$rulerKeyPos}->{'last_no_rular_wl_key_arry'} }  ){
          foreach my $allPtKey (   keys (  %{ $All_Key_type }  )   ){        	
        	  #if (   (  defined ( $outHash->{'2_alnShwHs'}->{$stepHer}->{$allPtKey} )  ) && (  ( $outHash->{'2_alnShwHs'}->{$stepHer}->{$allPtKey} eq ' ' ) || ( $outHash->{'2_alnShwHs'}->{$stepHer}->{$allPtKey}=~/\S+/ )  )   ) {
        	  	
        	  #}
        	  #else {
        	  	$outHash->{'2_alnShwHs'}->{$stepHer}->{$allPtKey}=' ';  
        	    $outHash->{'2_alnShwHs'}->{$stepHer}->{$allPtKey}=$in_aln_Hash_2->{$alnShwHs_key_2}->{$wlK}->{$allPtKey} if (   (  defined ( $in_aln_Hash_2->{$alnShwHs_key_2}->{$wlK}->{$allPtKey} )  ) && ( $in_aln_Hash_2->{$alnShwHs_key_2}->{$wlK}->{$allPtKey}=~m/\S+/ )   );
        	  #}
        	}
        	$stepHer++;
        }
        
	  	}
	  	
	  	if (
	  	         (   (  defined ( $All_befor_part_KeyHASH_1->{$rulerKeyPos} )  ) &&  (  ref ( $All_befor_part_KeyHASH_1->{$rulerKeyPos} ) eq 'HASH'  ) && (  defined ( $All_befor_part_KeyHASH_1->{$rulerKeyPos}->{'itSelf'} )  ) &&  ( $All_befor_part_KeyHASH_1->{$rulerKeyPos}->{'itSelf'}=~m/\d+/ )  )
	  	      || (   (  defined ( $All_befor_part_KeyHASH_2->{$rulerKeyPos} )  ) &&  (  ref ( $All_befor_part_KeyHASH_2->{$rulerKeyPos} ) eq 'HASH'  ) && (  defined ( $All_befor_part_KeyHASH_2->{$rulerKeyPos}->{'itSelf'} )  ) &&  ( $All_befor_part_KeyHASH_2->{$rulerKeyPos}->{'itSelf'}=~m/\d+/ )  )
	  	      
	  	   ) {
	  	   	
	  	  
	  	  my $allPtKey_space_wright=0;  
	  	  my $CHar_ZorF_1; my $CHar_ZorF_2; my $CHar_1;
	  	  if(   (  defined ( $All_befor_part_KeyHASH_1->{$rulerKeyPos} )  ) &&  (  ref ( $All_befor_part_KeyHASH_1->{$rulerKeyPos} ) eq 'HASH'  ) && (  defined ( $All_befor_part_KeyHASH_1->{$rulerKeyPos}->{'itSelf'} )  ) &&  ( $All_befor_part_KeyHASH_1->{$rulerKeyPos}->{'itSelf'}=~m/\d+/ )  ){
	  	  	my $ptKeyHere=$All_befor_part_KeyHASH_1->{$rulerKeyPos}->{'itSelf'};
	  	  	$CHar_ZorF_1=$in_aln_Hash_1->{$alnShwHs_key_1}->{$ptKeyHere}->{$ruler_ZorF_1} if (   (  defined ( $in_aln_Hash_1->{$alnShwHs_key_1}->{$ptKeyHere}->{$ruler_ZorF_1} )  ) && ( $in_aln_Hash_1->{$alnShwHs_key_1}->{$ptKeyHere}->{$ruler_ZorF_1}=~m/\S+/ )   );
	  	  	#$CHar_1=     $in_aln_Hash_1->{$alnShwHs_key_1}->{$ptKeyHere}->{  };
	  	  	foreach my $allPtKey (   keys (  %{ $All_Key_type }  )   ){        	
        	  #if (   (  defined ( $outHash->{'2_alnShwHs'}->{$stepHer}->{$allPtKey} )  ) && (  ( $outHash->{'2_alnShwHs'}->{$stepHer}->{$allPtKey} eq ' ' ) || ( $outHash->{'2_alnShwHs'}->{$stepHer}->{$allPtKey}=~/\S+/ )  )   ){}
        	  #else {
        	    $outHash->{'2_alnShwHs'}->{$stepHer}->{$allPtKey}=' ';  
        	    $outHash->{'2_alnShwHs'}->{$stepHer}->{$allPtKey}=$in_aln_Hash_1->{$alnShwHs_key_1}->{$ptKeyHere}->{$allPtKey} if (   (  defined ( $in_aln_Hash_1->{$alnShwHs_key_1}->{$ptKeyHere}->{$allPtKey} )  ) && ( $in_aln_Hash_1->{$alnShwHs_key_1}->{$ptKeyHere}->{$allPtKey}=~m/\S+/ )   );
        	  	
        	  #}
        	  
        	} 
        	#$allPtKey_space_wright=1;
	  	  } 
	  	  if(   (  defined ( $All_befor_part_KeyHASH_2->{$rulerKeyPos} )  ) &&  (  ref ( $All_befor_part_KeyHASH_2->{$rulerKeyPos} ) eq 'HASH'  ) && (  defined ( $All_befor_part_KeyHASH_2->{$rulerKeyPos}->{'itSelf'} )  ) &&  ( $All_befor_part_KeyHASH_2->{$rulerKeyPos}->{'itSelf'}=~m/\d+/ )  ){
	  	  	my $ptKeyHere=$All_befor_part_KeyHASH_2->{$rulerKeyPos}->{'itSelf'};
	  	  	$CHar_ZorF_2=$in_aln_Hash_2->{$alnShwHs_key_2}->{$ptKeyHere}->{$ruler_ZorF_2} if (   (  defined ( $in_aln_Hash_2->{$alnShwHs_key_2}->{$ptKeyHere}->{$ruler_ZorF_2} )  ) && ( $in_aln_Hash_2->{$alnShwHs_key_2}->{$ptKeyHere}->{$ruler_ZorF_2}=~m/\S+/ )   );
	  	  	foreach my $allPtKey (   keys (  %{ $All_Key_type }  )   ){        	
        	  #if (   (  defined ( $outHash->{'2_alnShwHs'}->{$stepHer}->{$allPtKey} )  ) && (  ( $outHash->{'2_alnShwHs'}->{$stepHer}->{$allPtKey} eq ' ' ) || ( $outHash->{'2_alnShwHs'}->{$stepHer}->{$allPtKey}=~/\S+/ )  )   ) {}
        	  #else { 
        	  if ($allPtKey_space_wright==0){
        	  	if (   (  defined ( $outHash->{'2_alnShwHs'}->{$stepHer}->{$allPtKey} )  ) && (  ( $outHash->{'2_alnShwHs'}->{$stepHer}->{$allPtKey} eq ' ' ) || ( $outHash->{'2_alnShwHs'}->{$stepHer}->{$allPtKey}=~/\S+/ )  )   ) {}
        	  	else {
        	  		$outHash->{'2_alnShwHs'}->{$stepHer}->{$allPtKey}=' ';   
        	  	} 
        	  }       	  	
        	  $outHash->{'2_alnShwHs'}->{$stepHer}->{$allPtKey}=$in_aln_Hash_2->{$alnShwHs_key_2}->{$ptKeyHere}->{$allPtKey} if (   (  defined ( $in_aln_Hash_2->{$alnShwHs_key_2}->{$ptKeyHere}->{$allPtKey} )  ) && ( $in_aln_Hash_2->{$alnShwHs_key_2}->{$ptKeyHere}->{$allPtKey}=~m/\S+/ )   );
        	  
        	} 
	  	  } 
	  	  #if      (  ( ($CHar_ZorF_1 eq '+') && ($CHar_ZorF_2 eq '+') ) || ( ($CHar_ZorF_1 eq '-') && ($CHar_ZorF_2 eq '-') )   ){	  	  		  	  }
	  	  #elsif   (  ( ($CHar_ZorF_1 eq '+') && ($CHar_ZorF_2 eq '-') ) || ( ($CHar_ZorF_1 eq '-') && ($CHar_ZorF_2 eq '+') )   ){	  	  		  	  
	  	  #  if ()	
	  	  #}
	  	  $stepHer++;
	  	}
	  	
	  	################ when the position is the last one ,then just add the after array 
	  	if  (   ( defined ($hash_1_the_last_key) ) && ($rulerKeyPos == $hash_1_the_last_key) && (  defined ( $All_after_part_KeyHASH_1->{$hash_1_the_last_key} )  ) &&  (  ref ( $All_after_part_KeyHASH_1->{$hash_1_the_last_key} ) eq 'HASH'  ) && (  defined ( $All_after_part_KeyHASH_1->{$hash_1_the_last_key}->{'last_no_rular_wl_key_arry'} )  ) &&  (  ref ( $All_after_part_KeyHASH_1->{$hash_1_the_last_key}->{'last_no_rular_wl_key_arry'} ) eq'ARRAY'  )  ){
        foreach my $wlK	 (  @{ $All_after_part_KeyHASH_1->{$hash_1_the_last_key}->{'last_no_rular_wl_key_arry'} }  ){
          foreach my $allPtKey (   keys (  %{ $All_Key_type }  )   ){        	
        	  #if (   (  defined ( $outHash->{'2_alnShwHs'}->{$stepHer}->{$allPtKey} )  ) && (  ( $outHash->{'2_alnShwHs'}->{$stepHer}->{$allPtKey} eq ' ' ) || ( $outHash->{'2_alnShwHs'}->{$stepHer}->{$allPtKey}=~/\S+/ )  )   ) {}
        	  #else {
        	    $outHash->{'2_alnShwHs'}->{$stepHer}->{$allPtKey}=' '; 
        	    $outHash->{'2_alnShwHs'}->{$stepHer}->{$allPtKey}=$in_aln_Hash_1->{$alnShwHs_key_1}->{$wlK}->{$allPtKey} if (   (  defined ( $in_aln_Hash_1->{$alnShwHs_key_1}->{$wlK}->{$allPtKey} )  ) && ( $in_aln_Hash_1->{$alnShwHs_key_1}->{$wlK}->{$allPtKey}=~m/\S+/ )   );
        	  #}
        	}
        	$stepHer++;
        }
        
	  	}
	  	
	  	if  (   ( defined ($hash_2_the_last_key) ) && ($rulerKeyPos == $hash_2_the_last_key) && (  defined ( $All_after_part_KeyHASH_2->{$hash_2_the_last_key} )  ) &&  (  ref ( $All_after_part_KeyHASH_2->{$hash_2_the_last_key} ) eq 'HASH'  ) && (  defined ( $All_after_part_KeyHASH_2->{$hash_2_the_last_key}->{'last_no_rular_wl_key_arry'} )  ) &&  (  ref ( $All_after_part_KeyHASH_2->{$hash_2_the_last_key}->{'last_no_rular_wl_key_arry'} ) eq'ARRAY'  )  ){
        foreach my $wlK	 (  @{ $All_after_part_KeyHASH_2->{$hash_2_the_last_key}->{'last_no_rular_wl_key_arry'} }  ){
          foreach my $allPtKey (   keys (  %{ $All_Key_type }  )   ){        	
        	  #if (   (  defined ( $outHash->{'2_alnShwHs'}->{$stepHer}->{$allPtKey} )  ) && (  ( $outHash->{'2_alnShwHs'}->{$stepHer}->{$allPtKey} eq ' ' ) || ( $outHash->{'2_alnShwHs'}->{$stepHer}->{$allPtKey}=~/\S+/ )  )   ) {}
        	  #else {
        	  	$outHash->{'2_alnShwHs'}->{$stepHer}->{$allPtKey}=' ';  
        	  	$outHash->{'2_alnShwHs'}->{$stepHer}->{$allPtKey}=$in_aln_Hash_2->{$alnShwHs_key_2}->{$wlK}->{$allPtKey} if (   (  defined ( $in_aln_Hash_2->{$alnShwHs_key_2}->{$wlK}->{$allPtKey} )  ) && ( $in_aln_Hash_2->{$alnShwHs_key_2}->{$wlK}->{$allPtKey}=~m/\S+/ )   );
        	  #}
        	}
        	$stepHer++;
        }
        
	  	}
      #################
	  	
	  }
	  
	}	
	
	return $outHash;
	
}

sub GetAllRuler_befor_part_KeyHASH{  #SplignHandle::GetAllRuler_befor_part_KeyHASH  
	my ( $inHash, $rulerPosKey)=@_;
	
	
  my $warnMsgBody="\nIn package  SplignHandle,\tIn sub GetAllRuler_befor_part_KeyHASH,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	
	
	my $returnHash;
	
	if (  ref ($inHash) eq 'HASH' ){
		
		my $alnShwHs_key=&Get_alnShwHs_key($inHash);
		my $diercti=&CheckDirection($inHash, $rulerPosKey);
	
	   my @sortKeyArr=sort { $a <=> $b } (   keys (  %{ $inHash->{$alnShwHs_key} }  )   ); 
	   if ($diercti eq '-') { 
	   	@sortKeyArr= reverse @sortKeyArr; 
	    
	  }
		
		
		if (  ref ( $inHash->{$alnShwHs_key} ) eq 'HASH' ){
			
			my @last_no_rular_wl_key_arry=(); 
     
			foreach my $alKey ( @sortKeyArr ){
				#push @last_no_rular_wl_key_arry, $alKey;
				if (   (  defined ( $inHash->{$alnShwHs_key}->{$alKey}->{$rulerPosKey} )  ) && ( $inHash->{$alnShwHs_key}->{$alKey}->{$rulerPosKey}=~m/\d+/ )   ){
					$returnHash->{ $inHash->{$alnShwHs_key}->{$alKey}->{$rulerPosKey } }->{'last_no_rular_wl_key_arry'}=[ @last_no_rular_wl_key_arry ] ;	@last_no_rular_wl_key_arry=();	
					$returnHash->{ $inHash->{$alnShwHs_key}->{$alKey}->{$rulerPosKey } }->{'itSelf'}=	$alKey;		
				}
				else {
					push @last_no_rular_wl_key_arry, $alKey;					
				}
				
			}	
			
		}
		else {
		  my $dieMsg1="$die_MsgHead\$inHash->{$alnShwHs_key}=$inHash->{$alnShwHs_key} should be a Hash Ref !!!!!!!!\n\n\n";
		  print $dieMsg1; die $dieMsg1;
	  }
				
	}
	else {
		my $dieMsg="$die_MsgHead\$inHash=$inHash should be a Hash Ref !!!!!!!!\n\n\n";
		print $dieMsg; die $dieMsg;
	}
	
	return $returnHash;
	
}

sub GetAllRuler_after_part_KeyHASH{  #SplignHandle::GetAllRuler_befor_part_KeyHASH  
	my ( $inHash, $rulerPosKey)=@_;
	
	
  my $warnMsgBody="\nIn package  SplignHandle,\tIn sub GetAllRuler_befor_part_KeyHASH,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	
	
	my $returnHash;
	
	if (  ref ($inHash) eq 'HASH' ){
		
		my $alnShwHs_key=&Get_alnShwHs_key($inHash);
		my $diercti=&CheckDirection($inHash, $rulerPosKey);
	
	   my @sortKeyArr=sort { $b <=> $a } (   keys (  %{ $inHash->{$alnShwHs_key} }  )   ); 
	   if ($diercti eq '-') { 
	   	@sortKeyArr= reverse @sortKeyArr; 
	    
	  }
		
		
		if (  ref ( $inHash->{$alnShwHs_key} ) eq 'HASH' ){
			
			my @last_no_rular_wl_key_arry=();      
			foreach my $alKey ( @sortKeyArr ){
				#push @last_no_rular_wl_key_arry, $alKey;
				if (   (  defined ( $inHash->{$alnShwHs_key}->{$alKey}->{$rulerPosKey} )  ) && ( $inHash->{$alnShwHs_key}->{$alKey}->{$rulerPosKey}=~m/\d+/ )   ){
					
					$returnHash->{ $inHash->{$alnShwHs_key}->{$alKey}->{$rulerPosKey } }->{'last_no_rular_wl_key_arry'}=[ @last_no_rular_wl_key_arry ] ;	@last_no_rular_wl_key_arry=();	
					$returnHash->{ $inHash->{$alnShwHs_key}->{$alKey}->{$rulerPosKey } }->{'itSelf'}=	$alKey;		
				}
				else {
					
					unshift @last_no_rular_wl_key_arry, $alKey;					
				}				
			}	
			
		}
		else {
		  my $dieMsg1="$die_MsgHead\$inHash->{$alnShwHs_key}=$inHash->{$alnShwHs_key} should be a Hash Ref !!!!!!!!\n\n\n";
		  print $dieMsg1; die $dieMsg1;
	  }
				
	}
	else {
		my $dieMsg="$die_MsgHead\$inHash=$inHash should be a Hash Ref !!!!!!!!\n\n\n";
		print $dieMsg; die $dieMsg;
	}
	
	return $returnHash;
	
}


sub GetAllRulerKeyHASH{  #SplignHandle::GetAllRulerKeyHASH  
	my ($returnHash, $inHash, $rulerPosKey)=@_;
	
	
  my $warnMsgBody="\nIn package  SplignHandle,\tIn sub GetAllRulerKeyHASH,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	
	
	
	
	if (  ref ($inHash) eq 'HASH' ){
		
		my $alnShwHs_key=&Get_alnShwHs_key($inHash);
		my $diercti=&CheckDirection($inHash, $rulerPosKey);
	
	  my @sortKeyArr=sort { $a <=> $b } (   keys (  %{ $inHash->{$alnShwHs_key} }  )   ); if ($diercti eq '-') { @sortKeyArr= reverse @sortKeyArr; }
		
		
		if (  ref ( $inHash->{$alnShwHs_key} ) eq 'HASH' ){
			

			foreach my $alKey ( @sortKeyArr ){
				if (   (  defined ( $inHash->{$alnShwHs_key}->{$alKey}->{$rulerPosKey} )  ) && ( $inHash->{$alnShwHs_key}->{$alKey}->{$rulerPosKey}=~m/\d+/ )   ){
					$returnHash->{ $inHash->{$alnShwHs_key}->{$alKey}->{$rulerPosKey } }=1;				
				}				
			}	
			
		}
		else {
		  my $dieMsg1="$die_MsgHead\$inHash->{$alnShwHs_key}=$inHash->{$alnShwHs_key} should be a Hash Ref !!!!!!!!\n\n\n";
		  print $dieMsg1; die $dieMsg1;
	  }
				
	}
	else {
		my $dieMsg="$die_MsgHead\$inHash=$inHash should be a Hash Ref !!!!!!!!\n\n\n";
		print $dieMsg; die $dieMsg;
	}
	
	return $returnHash;
	
}


sub CheckDirection{ #SplignHandle::CheckDirection  
	my ($inHash, $rulerPosKey)=@_;
	
	
  my $warnMsgBody="\nIn package  SplignHandle,\tIn sub CheckDirection,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	
	if (  ref ($inHash) eq 'HASH' ){
		
		my $alnShwHs_key=&Get_alnShwHs_key($inHash);
		if (  ref ( $inHash->{$alnShwHs_key} ) eq 'HASH' ){
			
			my $diret; my $lastNb; my $now_nb; my $bigNb=0; my $smlNb=0; my $equNb=0;
			
			my $step=0;
			foreach my $alKey (    sort { $a <=> $b}  (   keys (  %{ $inHash->{$alnShwHs_key} }  )   )    ){
				if (   (  defined ( $inHash->{$alnShwHs_key}->{$alKey}->{$rulerPosKey} )  ) && ( $inHash->{$alnShwHs_key}->{$alKey}->{$rulerPosKey}=~m/\d+/ )   ){
					$now_nb=$inHash->{$alnShwHs_key}->{$alKey}->{$rulerPosKey}; #warn "\$now_nb=\$inHash->{$alnShwHs_key}->{$alKey}->{$rulerPosKey}=$now_nb \$lastNb=$lastNb\n";
					if ($step>0){
						if    ($now_nb> $lastNb){					$bigNb++;							}
						elsif ($now_nb< $lastNb){					$smlNb++;							}
						elsif ($now_nb==$lastNb){					$equNb++;							}
					}
					#warn "           \$bigNb=$bigNb\t\$smlNb=$smlNb\n";
					
					$lastNb=$now_nb;
				}
				$step++;
			}
			
			#warn "2222           \$bigNb=$bigNb\t\$smlNb=$smlNb\n";
			if    ( $bigNb > $smlNb ) { $diret='+'; } # warn "3333           \$bigNb=$bigNb\t\$smlNb=$smlNb\n";}
			elsif ( $bigNb < $smlNb ) { $diret='-'; } # warn "4444           \$bigNb=$bigNb\t\$smlNb=$smlNb\n";}
			else {
				my $dieMsg2="$die_MsgHead\$rulerPosKey=$rulerPosKey\$bigNb=$bigNb\t\$smlNb=$smlNb!!!!!!!!\n\n\n";  
		    print $dieMsg2;   DirFileHandle::PrintAndWarnDumper ( $inHash->{$alnShwHs_key} );   print $dieMsg2;
		    die $dieMsg2;
			}
			return $diret;
			
		}
		else {
		  my $dieMsg1="$die_MsgHead\$inHash->{$alnShwHs_key}=$inHash->{$alnShwHs_key} should be a Hash Ref !!!!!!!!\n\n\n"; DirFileHandle::PrintAndWarnDumper ($inHash->{$alnShwHs_key});
		  print $dieMsg1; die $dieMsg1;
	  }
		
		
	}
	else {
		my $dieMsg="$die_MsgHead\$inHash=$inHash should be a Hash Ref !!!!!!!!\n\n\n";
		print $dieMsg; die $dieMsg;
	}
	
	
}

sub GetAllKey_from_hash{  #SplignHandle::GetAllKey_from_hash  
	
	my ($inHash)=@_;
	
	
  my $warnMsgBody="\nIn package  SplignHandle,\tIn sub GetAllKey_from_hash,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	
	my $outHash;
	if (  ref ($inHash) eq 'HASH' ){
	  
	  foreach my $keyH (    sort (   keys (  %{ $inHash }  )   )   ){
	  	if (  ref ( $inHash->{$keyH} ) eq 'HASH'  ){
	  		
	  		foreach my $keyH (    sort (   keys (  %{ $inHash->{$keyH} }  )   )   ){
	  			
	  			$outHash->{$keyH}=1;
	  			
	  		}
	  		
	  	}
	  }
	  
	}
	else {
	  my $dieMsg="$die_MsgHead\$inHash=$inHash should be a HASH ref\n\n\n";
		print $dieMsg; die $dieMsg;
	}
  
  return $outHash

	
	
}


sub Get_alnShwHs_key{ #SplignHandle::Get_alnShwHs_key  
	
	my ($inHash)=@_;
	
	
  my $warnMsgBody="\nIn package  SplignHandle,\tIn sub Get_alnShwHs_key,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $alnShwHs_key;
  if     (  ref ($inHash) eq 'HASH' ){
		
		my $alnShwHs_key; my $found=0; my $diePtKeyString;
		foreach my $alKey (    sort (   keys (  %{ $inHash }  )   )    ){
			$diePtKeyString.="($alKey)";
			if ($alKey=~m/alnShwHs/){
				$alnShwHs_key=$alKey;  $found++;
				
			}
		}
		if ($found==1){
			#warn "$alnShwHs_key  $diePtKeyString、\n";
			return $alnShwHs_key;
		}
		else {
			my $dieMsg="$die_MsgHead\$found=$found\t\$alnShwHs_key=$alnShwHs_key should be a alnShwHs_key !!!!!!!!All key: $diePtKeyString\n\n\n";
		  print $dieMsg; die $dieMsg;
		}
		
		
		
	}
	else {
		my $dieMsg="$die_MsgHead\$inHash=$inHash should be a Hash Ref !!!!!!!!\n\n\n";
		print $dieMsg; die $dieMsg;
	}
  
  
	
	
	
}














1;

##########################################################################################################################################
# 
