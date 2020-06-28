
#!/usr/bin/perl -w

BEGIN{  push (@INC, '/home/fredjiang/fredProgs/PMs')  };

use strict;
use warnings;

use PrintSubArrayHash;
use FastaFileHandle;
use SeqSegmentsTools;
use BlastHandle;
use DirFileHandle;

package  SECISwork;


sub ChangeTheAllSECISinformHashStructure{  #  SECISwork::ChangeTheAllSECISinformHashStructure ($inSECISinformHash);
	my ($inSECISinformHash)=@_;
	
	my $warnMsgBody="\nIn package  SECISwork,\tIn sub ChangeTheAllSECISinformHashStructure,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	
	
	#'Thalassiosira_oceanica20150522' => {
  #                                              'Est' => {
  #                                                         '1sig1170' => {
  #                                                                         '+1' => {
  #                                                                                   '118' => [
  #                                                                                              {
  #                                                                                                '_CoveScoKwd' => '0',
  #                                                                                                '_SECIS_5_partKwd' => 'GACACGAGAGGGUGACCUCAAGUC UUCAACUCCUG',
  #                                                                                                '_SamDifPlsMins_DisMK' => 219,
  #                                                                                                '_SECIS_4_r2__part' => 'AA ',
  #                                                                                                '_rlTGA_SECIS_dis' => 219,
  #                                                                                                '_totalHitsKwd' => '10',
  #                                                                                                '_SECISsechSeqDircK' => '+',
  #                                                                                                '_upstemEnKwd' => '-8.00',
  #                                                                                                '_fullstructEnKwd' => '-26.40',
  #                                                                                                '_stdNUMBER' => 1,
  #                                                                                                '_HitOrderNbKwd' => '10',
  #                                                                                                '_SECIS_3_partKwd' => 'UGACGGAUGGA ',
  #                                                                                                '_secisTLK' => '447',
  #                                                                                                '_inFastaNBkwd' => '39',
  #                                                                                                '_SECIS_7_partKwd' => 'AUGAAUG UUGAUUU UAUUCCCAUC',
  #                                                                                                '_pngPathFileKwd' => '/home/fredjiang/work/Algae/20150522NewDATA/20150611_Out/SECISworkTest/1023/out/std1/3072_files/188569.png',
  #                                                                                                '_secisHDK' => '337',
  #                                                                                                '_secisDrctK' => '+',
  #                                                                                                '_SameDiffPlusMinusMK' => 0,
  #                                                                                                '_SECIS_2_r1__part' => 'GUGAA ',
  #                                                                                                '_SECIS_1_partKwd' => 'GGGUAGGGAA UGAAAAG CAGAUCGGCAUUU ',
  #                                                                                                '_SECIS_6_r3__part' => 'AGAA '
  #                                                                                              }
  #                                                                                            ]
  #                                                                                 }
  #                                                                       },
  #                                                       },
  #                                     }
  
  my $otHash;                                       
	my $inHash=$inSECISinformHash;
	if (  ( ref ($inHash) ) eq 'HASH' ){
    foreach my $keyLev_0 (    sort { $a cmp $b } (   keys (  %{ $inHash } )   )    ){   #print "    \$keyLev_0=$keyLev_0\n";
    	my $valLev_0=$inHash->{$keyLev_0};                                                #print "    \$valLev_0=$valLev_0\n";
    	my $refLev_0=ref ($valLev_0);                                                     #print "    \$refLev_0=$refLev_0\n\n";  
    	
      if ( $refLev_0 eq 'HASH' ){
        foreach my $keyLev_1 (    sort { $a cmp $b } (   keys (  %{ $valLev_0 } )   )    ){   #print "      \$keyLev_1=$keyLev_1\n";
        	my $valLev_1=$valLev_0->{$keyLev_1};                                                #print "      \$valLev_1=$valLev_1\n";                
        	my $refLev_1=ref ($valLev_1);                                                       #print "      \$refLev_1=$refLev_1\n\n";
        	
      	  if ( $refLev_1 eq 'HASH' ){
            foreach my $keyLev_2 (    sort { $a cmp $b } (   keys (  %{ $valLev_1 } )   )    ){   #print "        \$keyLev_2=$keyLev_2\n";
            	my $valLev_2=$valLev_1->{$keyLev_2};                                                #print "        \$valLev_2=$valLev_2\n";
            	my $refLev_2=ref ($valLev_2);                                                       #print "        \$refLev_2=$refLev_2\n\n";
     
            	if ( $refLev_2 eq 'HASH' ){
                foreach my $keyLev_3 (    sort { $a cmp $b } (   keys (  %{ $valLev_2 } )   )    ){   #print "          \$keyLev_3=$keyLev_3\n";
                	my $valLev_3=$valLev_2->{$keyLev_3};                                                #print "          \$valLev_3=$valLev_3\n";
                	my $refLev_3=ref ($valLev_3);                                                       #print "          \$refLev_3=$refLev_3\n\n";
                	
                	if ( $refLev_3 eq 'HASH' ){
                    foreach my $keyLev_4 (    sort { $a <=> $b } (   keys (  %{ $valLev_3 } )   )    ){   #print "            \$keyLev_4=$keyLev_4\n";
                    	my $valLev_4=$valLev_3->{$keyLev_4};                                                #print "            \$valLev_4=$valLev_4\n";
                    	my $refLev_4=ref ($valLev_4);                                                       #print "            \$refLev_4=$refLev_4\n\n";
       
            	        if ( $refLev_4 eq 'ARRAY' ){
            	        	for (  my $idxLev_5=0; $idxLev_5 < @{ $valLev_4 }; $idxLev_5++  ){                     #print "              \$idxLev_5=$idxLev_5\n";
                        #foreach my $keyLev_5 (    sort { $a cmp $b } (   keys (  %{ $valLev_4 } )   )    ){   print "              \$keyLev_5=$keyLev_5\n";
                        	my $valLev_5=$valLev_4->[$idxLev_5];                                                #print "              \$valLev_5=$valLev_5\n";
                        	my $refLev_5=ref ($valLev_5);                                                       #print "              \$refLev_5=$refLev_5\n\n";
                        	
                        	if ( $refLev_5 eq 'HASH' ){                        		
                        		
                        		my $tempHash;
                        		$tempHash->{ '30_tempTGAframe'	          }=$keyLev_3;
                        		$tempHash->{ '31_tempTGAPosit'	          }=$keyLev_4;
                        		
                            foreach my $keyLev_6 (    sort { $a cmp $b } (   keys (  %{ $valLev_5 } )   )    ){               #print "              \$keyLev_6=$keyLev_6\n";
                        	    my $valLev_6=$valLev_5->{$keyLev_6} if (  defined  ( $valLev_5->{$keyLev_6} )  );       #print "              \$valLev_6=$valLev_6\n";
                        	    #my $refLev_6=ref ($valLev_6);                                                                   print "              \$refLev_6=$refLev_6\n\n\n\n";                            
                              
                              if    ( $keyLev_6 eq '_CoveScoKwd'           ) {     $tempHash->{ '00_CoveScoKwd'	          }=$valLev_6;  }
                              elsif ( $keyLev_6 eq '_stdNUMBER'            ) {     $tempHash->{ '01_stdNUMBER'	          }=$valLev_6;  }
                              elsif ( $keyLev_6 eq '_secisDrctK'           ) {     $tempHash->{ '02_secisDrctK'	          }=$valLev_6;  }
                              elsif ( $keyLev_6 eq '_secisHDK'             ) {     $tempHash->{ '03_secisHDK'	            }=$valLev_6;  }
                              elsif ( $keyLev_6 eq '_secisTLK'             ) {     $tempHash->{ '04_secisTLK'	            }=$valLev_6;  }
                              elsif ( $keyLev_6 eq '_fullstructEnKwd'      ) {     $tempHash->{ '05_fullstructEnKwd'	    }=$valLev_6;  }
                              elsif ( $keyLev_6 eq '_upstemEnKwd'          ) {     $tempHash->{ '06_upstemEnKwd'	        }=$valLev_6;  }
                              elsif ( $keyLev_6 eq '_SECIS_1_partKwd'      ) {     $tempHash->{ '11_SECIS_1_partKwd'	    }=$valLev_6;  }
                              elsif ( $keyLev_6 eq '_SECIS_2_r1__part'     ) {     $tempHash->{ '12_SECIS_2_r1__part'	    }=$valLev_6;  }
                              elsif ( $keyLev_6 eq '_SECIS_3_partKwd'      ) {     $tempHash->{ '13_SECIS_3_partKwd'	    }=$valLev_6;  }
                              elsif ( $keyLev_6 eq '_SECIS_4_r2__part'     ) {     $tempHash->{ '14_SECIS_4_r2__part'	    }=$valLev_6;  }
                              elsif ( $keyLev_6 eq '_SECIS_5_partKwd'      ) {     $tempHash->{ '15_SECIS_5_partKwd'	    }=$valLev_6;  }
                              elsif ( $keyLev_6 eq '_SECIS_6_r3__part'     ) {     $tempHash->{ '16_SECIS_6_r3__part'	    }=$valLev_6;  }
                              elsif ( $keyLev_6 eq '_SECIS_7_partKwd'      ) {     $tempHash->{ '17_SECIS_7_partKwd'	    }=$valLev_6;  }
                              elsif ( $keyLev_6 eq '_pngPathFileKwd'       ) {     $tempHash->{ '20_pngPathFileKwd'	      }=$valLev_6;  }
                              elsif ( $keyLev_6 eq '_SECISsechSeqDircK'    ) {     $tempHash->{ '40_SECISsechSeqDircK'	  }=$valLev_6;  }
                              elsif ( $keyLev_6 eq '_SameDiffPlusMinusMK'	 ) {     $tempHash->{ '41_SameDiffPlusMinusMK'  }=$valLev_6;  }
                              elsif ( $keyLev_6 eq '_SamDifPlsMins_DisMK'	 ) {     $tempHash->{ '42_SamDifPlsMins_DisMK'  }=$valLev_6;  }
                              elsif ( $keyLev_6 eq '_rlTGA_SECIS_dis'      ) {     $tempHash->{ '43_rlTGA_SECIS_dis'	    }=$valLev_6;  }
                              elsif ( $keyLev_6 eq '_inFastaNBkwd'         ) {     $tempHash->{ '50_inFastaNBkwd'	        }=$valLev_6;  }
                              elsif ( $keyLev_6 eq '_totalHitsKwd'         ) {     $tempHash->{ '51_totalHitsKwd'	        }=$valLev_6;  }
                              elsif ( $keyLev_6 eq '_HitOrderNbKwd'        ) {     $tempHash->{ '52_HitOrderNbKwd'	      }=$valLev_6;  }                              
                              else                                           {     $tempHash->{ "60_$keyLev_6"	          }=$valLev_6;  }                            
                            
                            }
                            
                            if (        (    (   defined (  $tempHash->{ '11_SECIS_1_partKwd'	      }  )   ) && (  $tempHash->{ '11_SECIS_1_partKwd'	    }=~m/\S+/  )    )
                                    ||  (    (   defined (  $tempHash->{ '12_SECIS_2_r1__part'	    }  )   ) && (  $tempHash->{ '12_SECIS_2_r1__part'	    }=~m/\S+/  )    )
                                    ||  (    (   defined (  $tempHash->{ '13_SECIS_3_partKwd'	      }  )   ) && (  $tempHash->{ '13_SECIS_3_partKwd'	    }=~m/\S+/  )    )
                                    ||  (    (   defined (  $tempHash->{ '14_SECIS_4_r2__part'	    }  )   ) && (  $tempHash->{ '14_SECIS_4_r2__part'	    }=~m/\S+/  )    )
                                    ||  (    (   defined (  $tempHash->{ '15_SECIS_5_partKwd'	      }  )   ) && (  $tempHash->{ '15_SECIS_5_partKwd'	    }=~m/\S+/  )    )
                                    ||  (    (   defined (  $tempHash->{ '16_SECIS_6_r3__part'	    }  )   ) && (  $tempHash->{ '16_SECIS_6_r3__part'	    }=~m/\S+/  )    )
                                    ||  (    (   defined (  $tempHash->{ '17_SECIS_7_partKwd'	      }  )   ) && (  $tempHash->{ '17_SECIS_7_partKwd'	    }=~m/\S+/  )    )
                                       
                                ){
                                $tempHash->{ '18_SECIS_Sequence'	    }  =      $tempHash->{ '11_SECIS_1_partKwd'	      } if  (    (   defined (  $tempHash->{ '11_SECIS_1_partKwd'	    }  )   ) && (  $tempHash->{ '11_SECIS_1_partKwd'	      }=~m/\S+/  )    );
                                $tempHash->{ '18_SECIS_Sequence'	    } .= "  ".$tempHash->{ '12_SECIS_2_r1__part'	    } if  (    (   defined (  $tempHash->{ '12_SECIS_2_r1__part'	  }  )   ) && (  $tempHash->{ '12_SECIS_2_r1__part'	      }=~m/\S+/  )    );
                                $tempHash->{ '18_SECIS_Sequence'	    } .= "  ".$tempHash->{ '13_SECIS_3_partKwd'	      } if  (    (   defined (  $tempHash->{ '13_SECIS_3_partKwd'	    }  )   ) && (  $tempHash->{ '13_SECIS_3_partKwd'		    }=~m/\S+/  )    );
                                $tempHash->{ '18_SECIS_Sequence'	    } .= "  ".$tempHash->{ '14_SECIS_4_r2__part'      } if  (    (   defined (  $tempHash->{ '14_SECIS_4_r2__part'    }  )   ) && (  $tempHash->{ '14_SECIS_4_r2__part' 	    }=~m/\S+/  )    );
                                $tempHash->{ '18_SECIS_Sequence'	    } .= "  ".$tempHash->{ '15_SECIS_5_partKwd'	      } if  (    (   defined (  $tempHash->{ '15_SECIS_5_partKwd'	    }  )   ) && (  $tempHash->{ '15_SECIS_5_partKwd'		    }=~m/\S+/  )    );
                                $tempHash->{ '18_SECIS_Sequence'	    } .= "  ".$tempHash->{ '16_SECIS_6_r3__part'      } if  (    (   defined (  $tempHash->{ '16_SECIS_6_r3__part'    }  )   ) && (  $tempHash->{ '16_SECIS_6_r3__part' 	    }=~m/\S+/  )    );
                                $tempHash->{ '18_SECIS_Sequence'	    } .= "  ".$tempHash->{ '17_SECIS_7_partKwd'	      } if  (    (   defined (  $tempHash->{ '17_SECIS_7_partKwd'	    }  )   ) && (  $tempHash->{ '17_SECIS_7_partKwd'		    }=~m/\S+/  )    );
                                	
                            }
                            
                            if (   (  defined ( $otHash->{$keyLev_0}->{$keyLev_1}->{$keyLev_2}->{$keyLev_4} )  ) && (  ref ( $otHash->{$keyLev_0}->{$keyLev_1}->{$keyLev_2}->{$keyLev_4} ) eq 'ARRAY'  )   ){
                        			my $findTheSameHash_YoN=0;
                        			for (  my $i=0; $i < @{ $otHash->{$keyLev_0}->{$keyLev_1}->{$keyLev_2}->{$keyLev_4} }; $i++  ){
                        				my $outHash_lev5_hash=$otHash->{$keyLev_0}->{$keyLev_1}->{$keyLev_2}->{$keyLev_4}->[$i];
                        				$findTheSameHash_YoN=&CheckTwo1DivHashTheSameOrNot( $tempHash, $outHash_lev5_hash );
                        			}
                        			if ($findTheSameHash_YoN==0){
                            	  push ( @{ $otHash->{$keyLev_0}->{$keyLev_1}->{$keyLev_2}->{$keyLev_4} }, $tempHash );
                              }
                        		}
                        		else {
                        			$otHash->{$keyLev_0}->{$keyLev_1}->{$keyLev_2}->{$keyLev_4}->[0]=$tempHash;
                        		}
                            
                            
                            
                          }  	
                        	
                        }
                      }
                      
                    }
                  }
                  
                }
              }
              
            }
          }
          
        }
      }
      
    }
  }
  
  return $otHash;
	
}


sub ChangeTheAllSECISinformHashStructure_WithOutTGAIform{  #  SECISwork::ChangeTheAllSECISinformHashStructure
	my ($inSECISinformHash, $genome_key_changeHash)=@_;
	
	my $warnMsgBody="\nIn package  SECISwork,\tIn sub ChangeTheAllSECISinformHashStructure,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $subCallereIfm=DirFileHandle::print_SubCallerInform;
	
	#'Thalassiosira_oceanica20150522' => {
  #                                              'Est' => {
  #                                                         '1sig1170' => {
  #                                                                         '+1' => {
  #                                                                                   '118' => [
  #                                                                                              {
  #                                                                                                '_CoveScoKwd' => '0',
  #                                                                                                '_SECIS_5_partKwd' => 'GACACGAGAGGGUGACCUCAAGUC UUCAACUCCUG',
  #                                                                                                '_SamDifPlsMins_DisMK' => 219,
  #                                                                                                '_SECIS_4_r2__part' => 'AA ',
  #                                                                                                '_rlTGA_SECIS_dis' => 219,
  #                                                                                                '_totalHitsKwd' => '10',
  #                                                                                                '_SECISsechSeqDircK' => '+',
  #                                                                                                '_upstemEnKwd' => '-8.00',
  #                                                                                                '_fullstructEnKwd' => '-26.40',
  #                                                                                                '_stdNUMBER' => 1,
  #                                                                                                '_HitOrderNbKwd' => '10',
  #                                                                                                '_SECIS_3_partKwd' => 'UGACGGAUGGA ',
  #                                                                                                '_secisTLK' => '447',
  #                                                                                                '_inFastaNBkwd' => '39',
  #                                                                                                '_SECIS_7_partKwd' => 'AUGAAUG UUGAUUU UAUUCCCAUC',
  #                                                                                                '_pngPathFileKwd' => '/home/fredjiang/work/Algae/20150522NewDATA/20150611_Out/SECISworkTest/1023/out/std1/3072_files/188569.png',
  #                                                                                                '_secisHDK' => '337',
  #                                                                                                '_secisDrctK' => '+',
  #                                                                                                '_SameDiffPlusMinusMK' => 0,
  #                                                                                                '_SECIS_2_r1__part' => 'GUGAA ',
  #                                                                                                '_SECIS_1_partKwd' => 'GGGUAGGGAA UGAAAAG CAGAUCGGCAUUU ',
  #                                                                                                '_SECIS_6_r3__part' => 'AGAA '
  #                                                                                              }
  #                                                                                            ]
  #                                                                                 }
  #                                                                       },
  #                                                       },
  #                                     }
  
  my $otHash;                                       
	my $inHash=$inSECISinformHash;
	if (  ( ref ($inHash) ) eq 'HASH' ){
    foreach my $keyLev_0 (    sort { $a cmp $b } (   keys (  %{ $inHash } )   )    ){   #print "    \$keyLev_0=$keyLev_0\n";
    	my $valLev_0=$inHash->{$keyLev_0};                                                #print "    \$valLev_0=$valLev_0\n";
    	my $refLev_0=ref ($valLev_0);                                                     #print "    \$refLev_0=$refLev_0\n\n";  
    	
      if ( $refLev_0 eq 'HASH' ){
        foreach my $keyLev_1 (    sort { $a cmp $b } (   keys (  %{ $valLev_0 } )   )    ){   #print "      \$keyLev_1=$keyLev_1\n";
        	my $valLev_1=$valLev_0->{$keyLev_1};                                                #print "      \$valLev_1=$valLev_1\n";                
        	my $refLev_1=ref ($valLev_1);                                                       #print "      \$refLev_1=$refLev_1\n\n";
        	
      	  if ( $refLev_1 eq 'HASH' ){
            foreach my $keyLev_2 (    sort { $a cmp $b } (   keys (  %{ $valLev_1 } )   )    ){   #print "        \$keyLev_2=$keyLev_2\n";
            	my $valLev_2=$valLev_1->{$keyLev_2};                                                #print "        \$valLev_2=$valLev_2\n";
            	my $refLev_2=ref ($valLev_2);                                                       #print "        \$refLev_2=$refLev_2\n\n";
     
            	if ( $refLev_2 eq 'HASH' ){
                foreach my $keyLev_3 (    sort { $a cmp $b } (   keys (  %{ $valLev_2 } )   )    ){   #print "          \$keyLev_3=$keyLev_3\n";
                	my $valLev_3=$valLev_2->{$keyLev_3};                                                #print "          \$valLev_3=$valLev_3\n";
                	my $refLev_3=ref ($valLev_3);                                                       #print "          \$refLev_3=$refLev_3\n\n";
                	
                	if ( $refLev_3 eq 'HASH' ){
                    foreach my $keyLev_4 (    sort { $a <=> $b } (   keys (  %{ $valLev_3 } )   )    ){   #print "            \$keyLev_4=$keyLev_4\n";
                    	my $valLev_4=$valLev_3->{$keyLev_4};                                                #print "            \$valLev_4=$valLev_4\n";
                    	my $refLev_4=ref ($valLev_4);                                                       #print "            \$refLev_4=$refLev_4\n\n";
       
            	        if ( $refLev_4 eq 'ARRAY' ){
            	        	for (  my $idxLev_5=0; $idxLev_5 < @{ $valLev_4 }; $idxLev_5++  ){                     #print "              \$idxLev_5=$idxLev_5\n";
                        #foreach my $keyLev_5 (    sort { $a cmp $b } (   keys (  %{ $valLev_4 } )   )    ){   print "              \$keyLev_5=$keyLev_5\n";
                        	my $valLev_5=$valLev_4->[$idxLev_5];                                                #print "              \$valLev_5=$valLev_5\n";
                        	my $refLev_5=ref ($valLev_5);                                                       #print "              \$refLev_5=$refLev_5\n\n";
                        	
                        	if ( $refLev_5 eq 'HASH' ){                        		
                        		
                        		my $tempHash;
                        		#$tempHash->{ '30_tempTGAframe'	          }=$keyLev_3;
                        		#$tempHash->{ '31_tempTGAPosit'	          }=$keyLev_4;
                        		
                        		my $secisDrct;
                            foreach my $keyLev_6 (    sort { $a cmp $b } (   keys (  %{ $valLev_5 } )   )    ){               #print "              \$keyLev_6=$keyLev_6\n";
                        	    my $valLev_6=$valLev_5->{$keyLev_6} if (  defined  ( $valLev_5->{$keyLev_6} )  );       #print "              \$valLev_6=$valLev_6\n";
                        	    #my $refLev_6=ref ($valLev_6);                                                                   print "              \$refLev_6=$refLev_6\n\n\n\n";                            
                              
                              if    ( $keyLev_6 eq '_CoveScoKwd'           ) {     $tempHash->{ '00_CoveScoKwd'	          }=$valLev_6;  }
                              elsif ( $keyLev_6 eq '_stdNUMBER'            ) {     $tempHash->{ '01_stdNUMBER'	          }=$valLev_6;  }
                              elsif ( $keyLev_6 eq '_secisDrctK'           ) {     $tempHash->{ '02_secisDrctK'	          }=$valLev_6;       $secisDrct=$valLev_6; }
                              elsif ( $keyLev_6 eq '_secisHDK'             ) {     $tempHash->{ '03_secisHDK'	            }=$valLev_6;  }
                              elsif ( $keyLev_6 eq '_secisTLK'             ) {     $tempHash->{ '04_secisTLK'	            }=$valLev_6;  }
                              elsif ( $keyLev_6 eq '_fullstructEnKwd'      ) {     $tempHash->{ '05_fullstructEnKwd'	    }=$valLev_6;  }
                              elsif ( $keyLev_6 eq '_upstemEnKwd'          ) {     $tempHash->{ '06_upstemEnKwd'	        }=$valLev_6;  }
                              elsif ( $keyLev_6 eq '_SECIS_1_partKwd'      ) {     $tempHash->{ '11_SECIS_1_partKwd'	    }=$valLev_6;  }
                              elsif ( $keyLev_6 eq '_SECIS_2_r1__part'     ) {     $tempHash->{ '12_SECIS_2_r1__part'	    }=$valLev_6;  }
                              elsif ( $keyLev_6 eq '_SECIS_3_partKwd'      ) {     $tempHash->{ '13_SECIS_3_partKwd'	    }=$valLev_6;  }
                              elsif ( $keyLev_6 eq '_SECIS_4_r2__part'     ) {     $tempHash->{ '14_SECIS_4_r2__part'	    }=$valLev_6;  }
                              elsif ( $keyLev_6 eq '_SECIS_5_partKwd'      ) {     $tempHash->{ '15_SECIS_5_partKwd'	    }=$valLev_6;  }
                              elsif ( $keyLev_6 eq '_SECIS_6_r3__part'     ) {     $tempHash->{ '16_SECIS_6_r3__part'	    }=$valLev_6;  }
                              elsif ( $keyLev_6 eq '_SECIS_7_partKwd'      ) {     $tempHash->{ '17_SECIS_7_partKwd'	    }=$valLev_6;  }
                              elsif ( $keyLev_6 eq '_pngPathFileKwd'       ) {     $tempHash->{ '20_pngPathFileKwd'	      }=$valLev_6;  }
                              #elsif ( $keyLev_6 eq '_SECISsechSeqDircK'    ) {     $tempHash->{ '40_SECISsechSeqDircK'	  }=$valLev_6;  }
                              #elsif ( $keyLev_6 eq '_SameDiffPlusMinusMK'	 ) {     $tempHash->{ '41_SameDiffPlusMinusMK'  }=$valLev_6;  }
                              #elsif ( $keyLev_6 eq '_SamDifPlsMins_DisMK'	 ) {     $tempHash->{ '42_SamDifPlsMins_DisMK'  }=$valLev_6;  }
                              #elsif ( $keyLev_6 eq '_rlTGA_SECIS_dis'      ) {     $tempHash->{ '43_rlTGA_SECIS_dis'	    }=$valLev_6;  }
                              #elsif ( $keyLev_6 eq '_inFastaNBkwd'         ) {     $tempHash->{ '50_inFastaNBkwd'	        }=$valLev_6;  }
                              #elsif ( $keyLev_6 eq '_totalHitsKwd'         ) {     $tempHash->{ '51_totalHitsKwd'	        }=$valLev_6;  }
                              #elsif ( $keyLev_6 eq '_HitOrderNbKwd'        ) {     $tempHash->{ '52_HitOrderNbKwd'	      }=$valLev_6;  }                              
                              #else                                           {     $tempHash->{ "60_$keyLev_6"	          }=$valLev_6;  }                            
                            
                            }
                            
                            
                            if (        (    (   defined (  $tempHash->{ '11_SECIS_1_partKwd'	      }  )   ) && (  $tempHash->{ '11_SECIS_1_partKwd'	    }=~m/\S+/  )    )
                                    ||  (    (   defined (  $tempHash->{ '12_SECIS_2_r1__part'	    }  )   ) && (  $tempHash->{ '12_SECIS_2_r1__part'	    }=~m/\S+/  )    )
                                    ||  (    (   defined (  $tempHash->{ '13_SECIS_3_partKwd'	      }  )   ) && (  $tempHash->{ '13_SECIS_3_partKwd'	    }=~m/\S+/  )    )
                                    ||  (    (   defined (  $tempHash->{ '14_SECIS_4_r2__part'	    }  )   ) && (  $tempHash->{ '14_SECIS_4_r2__part'	    }=~m/\S+/  )    )
                                    ||  (    (   defined (  $tempHash->{ '15_SECIS_5_partKwd'	      }  )   ) && (  $tempHash->{ '15_SECIS_5_partKwd'	    }=~m/\S+/  )    )
                                    ||  (    (   defined (  $tempHash->{ '16_SECIS_6_r3__part'	    }  )   ) && (  $tempHash->{ '16_SECIS_6_r3__part'	    }=~m/\S+/  )    )
                                    ||  (    (   defined (  $tempHash->{ '17_SECIS_7_partKwd'	      }  )   ) && (  $tempHash->{ '17_SECIS_7_partKwd'	    }=~m/\S+/  )    )
                                       
                                ){
                                $tempHash->{ '18_SECIS_Sequence'	    }  =      $tempHash->{ '11_SECIS_1_partKwd'	      } if  (    (   defined (  $tempHash->{ '11_SECIS_1_partKwd'	    }  )   ) && (  $tempHash->{ '11_SECIS_1_partKwd'	      }=~m/\S+/  )    );
                                $tempHash->{ '18_SECIS_Sequence'	    } .= "  ".$tempHash->{ '12_SECIS_2_r1__part'	    } if  (    (   defined (  $tempHash->{ '12_SECIS_2_r1__part'	  }  )   ) && (  $tempHash->{ '12_SECIS_2_r1__part'	      }=~m/\S+/  )    );
                                $tempHash->{ '18_SECIS_Sequence'	    } .= "  ".$tempHash->{ '13_SECIS_3_partKwd'	      } if  (    (   defined (  $tempHash->{ '13_SECIS_3_partKwd'	    }  )   ) && (  $tempHash->{ '13_SECIS_3_partKwd'		    }=~m/\S+/  )    );
                                $tempHash->{ '18_SECIS_Sequence'	    } .= "  ".$tempHash->{ '14_SECIS_4_r2__part'      } if  (    (   defined (  $tempHash->{ '14_SECIS_4_r2__part'    }  )   ) && (  $tempHash->{ '14_SECIS_4_r2__part' 	    }=~m/\S+/  )    );
                                $tempHash->{ '18_SECIS_Sequence'	    } .= "  ".$tempHash->{ '15_SECIS_5_partKwd'	      } if  (    (   defined (  $tempHash->{ '15_SECIS_5_partKwd'	    }  )   ) && (  $tempHash->{ '15_SECIS_5_partKwd'		    }=~m/\S+/  )    );
                                $tempHash->{ '18_SECIS_Sequence'	    } .= "  ".$tempHash->{ '16_SECIS_6_r3__part'      } if  (    (   defined (  $tempHash->{ '16_SECIS_6_r3__part'    }  )   ) && (  $tempHash->{ '16_SECIS_6_r3__part' 	    }=~m/\S+/  )    );
                                $tempHash->{ '18_SECIS_Sequence'	    } .= "  ".$tempHash->{ '17_SECIS_7_partKwd'	      } if  (    (   defined (  $tempHash->{ '17_SECIS_7_partKwd'	    }  )   ) && (  $tempHash->{ '17_SECIS_7_partKwd'		    }=~m/\S+/  )    );
                                	
                            }
                            
                            my $shortGenomeID=$keyLev_2;
                            if  (   ( $keyLev_1 eq 'OneFiledGeno' ) && (  defined ( $genome_key_changeHash )  ) && (  ref ( $genome_key_changeHash ) eq 'HASH'  )   ){
                            	if (    (  defined ( $genome_key_changeHash->{$keyLev_2} )  ) && ( $genome_key_changeHash->{$keyLev_2}=~m/\S+/ )   ) {
                            		$shortGenomeID=$genome_key_changeHash->{$keyLev_2};                            		
                            	}
                            	else{
                            		DieWork::Just_dieWork( $die_MsgHead."\n\$keyLev_1=$keyLev_1 \$keyLev_2=$keyLev_2 is not found in \$genome_key_changeHash=$genome_key_changeHash right!!! \n\n$! \n$subCallereIfm\n\n" );
                            	}
                            }
                            
                            if (   (  defined ( $otHash->{$keyLev_0}->{$keyLev_1}->{$shortGenomeID}->{$secisDrct} )  ) && (  ref ( $otHash->{$keyLev_0}->{$keyLev_1}->{$shortGenomeID}->{$secisDrct} ) eq 'ARRAY'  )   ){
                        			my $findTheSameHash_YoN=0;
                        			for (  my $i=0; $i < @{ $otHash->{$keyLev_0}->{$keyLev_1}->{$shortGenomeID}->{$secisDrct} }; $i++  ){
                        				my $outHash_lev5_hash=$otHash->{$keyLev_0}->{$keyLev_1}->{$shortGenomeID}->{$secisDrct}->[$i];
                        				$findTheSameHash_YoN=&CheckTwo1DivHashTheSameOrNot( $tempHash, $outHash_lev5_hash );
                        			}
                        			if ($findTheSameHash_YoN==0){
                            	  push ( @{ $otHash->{$keyLev_0}->{$keyLev_1}->{$shortGenomeID}->{$secisDrct} }, $tempHash );
                              }
                        		}
                        		else {
                        			$otHash->{$keyLev_0}->{$keyLev_1}->{$shortGenomeID}->{$secisDrct}->[0]=$tempHash;
                        		}
                            
                            
                            
                          }  	
                        	
                        }
                      }
                      
                    }
                  }
                  
                }
              }
              
            }
          }
          
        }
      }
      
    }
  }
  
  return $otHash;
	
}

sub CheckTwo1DivHashTheSameOrNot{
	my ($hash1, $hash2)=@_;
	
	my $warnMsgBody="\nIn package  SECISwork,\tIn sub CheckTwo1DivHashTheSameOrNot,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	
	my $finalReturn=0;
	if (  ( ref ($hash1) eq 'HASH' ) && ( ref ($hash2) eq 'HASH' )  ){
		
		my $foundTheSame_hash2val_for_all_hash1val=1;
    foreach my $hash1key (   keys (  %{ $hash1 } )   ){                                             #print "    \$hash1key=$hash1key\n";
    	my $hash1val=$hash1->{$hash1key} if (  defined  ( $hash1->{$hash1key} )  );                   #print "    \$hash1val=$hash1val\n";
    	
    	my $foundTheSame_hash2val_for_one_hash1val=0;
    	foreach my $hash2key (   keys (  %{ $hash2 } )   ){                                             #print "    \$hash2key=$hash2key\n";
    	  my $hash2val=$hash2->{$hash2key} if (  defined  ( $hash2->{$hash2key} )  );                   #print "    \$hash2val=$hash2val\n";
    	  
    	  my $defined_hash1val=defined ($hash1val);
    	  my $defined_hash2val=defined ($hash2val);
    	  
    	  if (  ( ($defined_hash1val==0) && ($defined_hash2val==0) ) || ( ($defined_hash1val==1) && ($defined_hash2val==1) && ($hash2val eq $hash1val) )  ){
    	  	$foundTheSame_hash2val_for_one_hash1val=1;
    	  }
    	  
    	}
    	
    	if ($foundTheSame_hash2val_for_one_hash1val==0){
    		$foundTheSame_hash2val_for_all_hash1val=0;
    	}
    	
    }
    
    my $foundTheSame_hash1val_for_all_hash2val=1; 
    foreach my $hash2key (   keys (  %{ $hash2 } )   ){                                             #print "    \$hash2key=$hash2key\n";
      my $hash2val=$hash2->{$hash2key} if (  defined  ( $hash2->{$hash2key} )  );                   #print "    \$hash2val=$hash2val\n";
      
      
      my $foundTheSame_hash1val_for_one_hash2val=0;
      foreach my $hash1key (   keys (  %{ $hash1 } )   ){                                             #print "    \$hash1key=$hash1key\n";
        my $hash1val=$hash1->{$hash1key} if (  defined  ( $hash1->{$hash1key} )  );                   #print "    \$hash1val=$hash1val\n";
        
        my $defined_hash1val=defined ($hash1val);
    	  my $defined_hash2val=defined ($hash2val);
        if (  ( ($defined_hash1val==0) && ($defined_hash2val==0) ) || ( ($defined_hash1val==1) && ($defined_hash2val==1) && ($hash2val eq $hash1val) )  ){
      	  $foundTheSame_hash1val_for_one_hash2val=1;
        }
        
      }
      
      if ($foundTheSame_hash1val_for_one_hash2val==0){
    	  $foundTheSame_hash1val_for_all_hash2val=0;
      }
    
    }
    
    if ( ($foundTheSame_hash2val_for_all_hash1val==1) && ($foundTheSame_hash1val_for_all_hash2val==1) ) {
    	$finalReturn=1;
    }	
    
    	
    
  }
  
  return $finalReturn;
	
}

sub GetBestSECISarray{  # SECISwork::GetBestSECISarray
	my ($tga_pos, $CDS_Z_or_F, $SECISarray, $CDSlastPos, $ecSelPro)=@_;
	
	my $warnMsgBody="\nIn package  SECISwork,\tIn sub GetBestSECISarray,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	
	my $sortArray=&SortSECISarray($tga_pos, $CDS_Z_or_F, $SECISarray, $CDSlastPos, $ecSelPro);
	
	my $bestArray;
	if ( ref ($sortArray) eq 'ARRAY' ){
		
		my $j=0;
    FORMK: for (  my $i=0; $i < @{ $sortArray }; $i++  ){ 
    	if ($j>2){
    		last FORMK;
    	}
    	if ($i==0){
    		$bestArray->[$j]=$sortArray->[$i];  print "0-20180908$ecSelPro\$bestArray->[$j]=\$sortArray->[$i]=$bestArray->[$j]\n";
    		$j++;
    	}
    	else {
    		if (   (  defined ( $sortArray->[$i] )  ) && (  defined (  $sortArray->[$i]->{'00_CoveScoKwd'} )  ) && ( $sortArray->[$i]->{'00_CoveScoKwd'} >0 ) && ($j<=2)  ){  
    			my $foundSameSECIS=0;                                                                             
    			if (  (  defined ( $bestArray )  )  && ( ref ($bestArray) eq 'ARRAY' )  ){
    				for (  my $k=0; $k<@{ $bestArray}; $k++  ){                                                     print "1-20180908$ecSelPro\$sortArray->[$i]->{'02_secisDrctK'}=$sortArray->[$i]->{'02_secisDrctK'}\t\$sortArray->[$i]->{'03_secisHDK'}=$sortArray->[$i]->{'03_secisHDK'}\t\$sortArray->[$i]->{'04_secisTLK'}=$sortArray->[$i]->{'04_secisTLK'}\t\n";
    					                                                                                              print "2-20180908$ecSelPro\$bestArray->[$k]->{'02_secisDrctK'}=$bestArray->[$k]->{'02_secisDrctK'}\t\$bestArray->[$k]->{'03_secisHDK'}=$bestArray->[$k]->{'03_secisHDK'}\t\$bestArray->[$k]->{'04_secisTLK'}=$bestArray->[$k]->{'04_secisTLK'}\t\n";
    					my $theSameOrNot=&CheckTheSECISTheSameOrNot( $sortArray->[$i], $bestArray->[$k], $ecSelPro );            print "3-20180908$ecSelPro\$theSameOrNot=$theSameOrNot\n";
    					if ($theSameOrNot ==1){$foundSameSECIS=1;}                                                    print "4-20180908$ecSelPro\$foundSameSECIS=$foundSameSECIS\n";
    				}
    			}                                                                                                 print "5-20180908$ecSelPro\$foundSameSECIS=$foundSameSECIS\n";
    			if ($foundSameSECIS==0){
    				$bestArray->[$j]=$sortArray->[$i];                                                              print "6-20180908$ecSelPro\$sortArray->[$i]->{'02_secisDrctK'}=$sortArray->[$i]->{'02_secisDrctK'}\t\$sortArray->[$i]->{'03_secisHDK'}=$sortArray->[$i]->{'03_secisHDK'}\t\$sortArray->[$i]->{'04_secisTLK'}=$sortArray->[$i]->{'04_secisTLK'}\t\n";
    					                                                                                              print "7-20180908$ecSelPro\$bestArray->[$j]->{'02_secisDrctK'}=$bestArray->[$j]->{'02_secisDrctK'}\t\$bestArray->[$j]->{'03_secisHDK'}=$bestArray->[$j]->{'03_secisHDK'}\t\$bestArray->[$j]->{'04_secisTLK'}=$bestArray->[$j]->{'04_secisTLK'}\t\n";
    		    $j++;
    			}
    			
    		}
    	}
    	
    }
  }
  	
	return $bestArray;
}

sub SortSECISarray{  # SECISwork::SortSECISarray
	my ($tga_pos, $CDS_Z_or_F, $SECISarray, $CDSlastPos, $ecSelPro)=@_;
	
	print "a2-20180906$ecSelPro SortSECISarrayINPUT: $tga_pos, $CDS_Z_or_F, $SECISarray, $CDSlastPos, $ecSelPro\n";
	
	my $warnMsgBody="\nIn package  SECISwork,\tIn sub SortSECISarray,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	
	#The SECIS must in the downstream region of TGA
	#The SECIS must in the same direction with the TGA
	#the sort order: COVE SCORE > top stem energy > full energy > std > distance
	if ( ref ($SECISarray) eq 'ARRAY' ){
		my $newOutHash;  my $j=0;
    FORNEXTMK: for (  my $i=0; $i < @{ $SECISarray }; $i++  ){   print "a2-20180906$ecSelPro\$SECISarray->[$i]->{'02_secisDrctK'}=$SECISarray->[$i]->{'02_secisDrctK'}\t\$SECISarray->[$i]->{'03_secisHDK'}=$SECISarray->[$i]->{'03_secisHDK'}\n";
    	if (    ( $SECISarray->[$i]->{'02_secisDrctK'} eq $CDS_Z_or_F ) 
    	     && ( 
    	             (  ( $CDS_Z_or_F eq '+') && ($CDS_Z_or_F eq $SECISarray->[$i]->{ '02_secisDrctK'} ) && (    $tga_pos <  $SECISarray->[$i]->{'03_secisHDK'} )  )  
    	          || (  ( $CDS_Z_or_F eq '-') && ($CDS_Z_or_F eq $SECISarray->[$i]->{ '02_secisDrctK'} ) && (    $tga_pos >  $SECISarray->[$i]->{'03_secisHDK'} )  )  
    	         )    
    	    ){
    	  
    	  if (   (  defined (  $SECISarray->[$i]->{ '60__orgEstNmKwd' } )  ) && (  $SECISarray->[$i]->{ '60__orgEstNmKwd' } =~m/\S+/  )   ){
    	  	next FORNEXTMK;
    	  }
    	  
        $newOutHash->[$j]->{ '00_CoveScoKwd'	          }=$SECISarray->[$i]->{ '00_CoveScoKwd'	          } if (   defined (  $SECISarray->[$i]->{ '00_CoveScoKwd'	      } )  ); 
        #$newOutHash->[$j]->{ '0a_CoveScoBig'	          }=0; $newOutHash->[$j]->{ '0a_CoveScoBig'      }=1 if (   (  defined ( $SECISarray->[$i]->{ '00_CoveScoKwd'} )  ) && ( $SECISarray->[$i]->{ '00_CoveScoKwd'}>0 )   );
        $newOutHash->[$j]->{ '01_stdNUMBER'	            }=$SECISarray->[$i]->{ '01_stdNUMBER'	            } if (   defined (  $SECISarray->[$i]->{ '01_stdNUMBER'	        } )  );  
        $newOutHash->[$j]->{ '02_secisDrctK'	          }=$SECISarray->[$i]->{ '02_secisDrctK'	          } if (   defined (  $SECISarray->[$i]->{ '02_secisDrctK'	      } )  );
        $newOutHash->[$j]->{ '03_secisHDK'	            }=$SECISarray->[$i]->{ '03_secisHDK'	            } if (   defined (  $SECISarray->[$i]->{ '03_secisHDK'	        } )  );
        $newOutHash->[$j]->{ '04_secisTLK'	            }=$SECISarray->[$i]->{ '04_secisTLK'	            } if (   defined (  $SECISarray->[$i]->{ '04_secisTLK'	        } )  );
        $newOutHash->[$j]->{ '05_fullstructEnKwd'	      }=0; $newOutHash->[$j]->{ '05_fullstructEnKwd' }=$SECISarray->[$i]->{ '05_fullstructEnKwd'	      } if (   defined (  $SECISarray->[$i]->{ '05_fullstructEnKwd'	      } )  );  
        $newOutHash->[$j]->{ '06_upstemEnKwd'	          }=0; $newOutHash->[$j]->{ '06_upstemEnKwd'     }=$SECISarray->[$i]->{ '06_upstemEnKwd'	          } if (   defined (  $SECISarray->[$i]->{ '06_upstemEnKwd'	          } )  );
        $newOutHash->[$j]->{ '11_SECIS_1_partKwd'	      }=$SECISarray->[$i]->{ '11_SECIS_1_partKwd'	      } if (   defined (  $SECISarray->[$i]->{ '11_SECIS_1_partKwd'	  } )  );  
        $newOutHash->[$j]->{ '12_SECIS_2_r1__part'	    }=$SECISarray->[$i]->{ '12_SECIS_2_r1__part'	    } if (   defined (  $SECISarray->[$i]->{ '12_SECIS_2_r1__part'  } )  );
        $newOutHash->[$j]->{ '13_SECIS_3_partKwd'	      }=$SECISarray->[$i]->{ '13_SECIS_3_partKwd'	      } if (   defined (  $SECISarray->[$i]->{ '13_SECIS_3_partKwd'	  } )  );  
        $newOutHash->[$j]->{ '14_SECIS_4_r2__part'	    }=$SECISarray->[$i]->{ '14_SECIS_4_r2__part'	    } if (   defined (  $SECISarray->[$i]->{ '14_SECIS_4_r2__part'  } )  );
        $newOutHash->[$j]->{ '15_SECIS_5_partKwd'	      }=$SECISarray->[$i]->{ '15_SECIS_5_partKwd'	      } if (   defined (  $SECISarray->[$i]->{ '15_SECIS_5_partKwd'	  } )  );  
        $newOutHash->[$j]->{ '16_SECIS_6_r3__part'	    }=$SECISarray->[$i]->{ '16_SECIS_6_r3__part'	    } if (   defined (  $SECISarray->[$i]->{ '16_SECIS_6_r3__part'  } )  );
        $newOutHash->[$j]->{ '17_SECIS_7_partKwd'	      }=$SECISarray->[$i]->{ '17_SECIS_7_partKwd'	      } if (   defined (  $SECISarray->[$i]->{ '17_SECIS_7_partKwd'	  } )  );  
        $newOutHash->[$j]->{ '18_SECIS_Sequence'	      }=$SECISarray->[$i]->{ '18_SECIS_Sequence'	      } if (   defined (  $SECISarray->[$i]->{ '18_SECIS_Sequence'	  } )  );  
        $newOutHash->[$j]->{ '20_pngPathFileKwd'	      }=$SECISarray->[$i]->{ '20_pngPathFileKwd'	      } if (   defined (  $SECISarray->[$i]->{ '20_pngPathFileKwd'	  } )  );
        $newOutHash->[$j]->{ '30_tempTGAframe'	        }=$SECISarray->[$i]->{ '30_tempTGAframe'	        } if (   defined (  $SECISarray->[$i]->{ '30_tempTGAframe'	    } )  );  
        $newOutHash->[$j]->{ '31_tempTGAPosit'	        }=$SECISarray->[$i]->{ '31_tempTGAPosit'	        } if (   defined (  $SECISarray->[$i]->{ '31_tempTGAPosit'	    } )  );
        
        #$newOutHash->[$j]->{ '40_SECISsechSeqDircK'	   }=$SECISarray->[$i]->{ '40_SECISsechSeqDircK'	    };  
        #$newOutHash->[$j]->{ '41_SameDiffPlusMinusMK'   }=$SECISarray->[$i]->{ '41_SameDiffPlusMinusMK'   }; 
        #$newOutHash->[$j]->{ '42_SamDifPlsMins_DisMK'   }=$SECISarray->[$i]->{ '42_SamDifPlsMins_DisMK'   }; 
        #$newOutHash->[$j]->{ '43_TGA_SECIS_dis'	        }=abs ($tga_pos    -  $SECISarray->[$i]->{'03_secisHDK'} ); 
        #$newOutHash->[$j]->{ '44_End_SECIS_dis'	        }=abs ($CDSlastPos -  $SECISarray->[$i]->{'03_secisHDK'} );  #print "a2-20180906$ecSelPro\$newOutHash->[$j]->{ '44_End_SECIS_dis'	        }=$newOutHash->[$j]->{ '44_End_SECIS_dis'	        }\n";
        
        if    ( $CDS_Z_or_F eq '+'){
        	$newOutHash->[$j]->{ '43_TGA_SECIS_dis'       }= $SECISarray->[$i]->{'03_secisHDK'} - $tga_pos       if (   (  defined ( $tga_pos    )  ) && ( $tga_pos=~m/\d+/    ) && ( $tga_pos > 0    )   );  
        	$newOutHash->[$j]->{ '44_End_SECIS_dis'       }= $SECISarray->[$i]->{'03_secisHDK'} - $CDSlastPos    if (   (  defined ( $CDSlastPos )  ) && ( $CDSlastPos=~m/\d+/ ) && ( $CDSlastPos > 0 )   );          	 
        }
        elsif ( $CDS_Z_or_F eq '-'){
        	$newOutHash->[$j]->{ '43_TGA_SECIS_dis'       }= $tga_pos    -  $SECISarray->[$i]->{'03_secisHDK'}   if (   (  defined ( $tga_pos    )  ) && ( $tga_pos=~m/\d+/    ) && ( $tga_pos > 0    )   );  
        	$newOutHash->[$j]->{ '44_End_SECIS_dis'       }= $CDSlastPos -  $SECISarray->[$i]->{'03_secisHDK'}   if (   (  defined ( $CDSlastPos )  ) && ( $CDSlastPos=~m/\d+/ ) && ( $CDSlastPos > 0 )   );   
        }
        
        
        $newOutHash->[$j]->{ '50_inFastaNBkwd'	        }=$SECISarray->[$i]->{ '50_inFastaNBkwd'	        } if (   defined (  $SECISarray->[$i]->{ '50_inFastaNBkwd'		  } )  );
        $newOutHash->[$j]->{ '51_totalHitsKwd'	        }=$SECISarray->[$i]->{ '51_totalHitsKwd'	        } if (   defined (  $SECISarray->[$i]->{ '51_totalHitsKwd'		  } )  );
        $newOutHash->[$j]->{ '52_HitOrderNbKwd'	        }=$SECISarray->[$i]->{ '52_HitOrderNbKwd'	        } if (   defined (  $SECISarray->[$i]->{ '52_HitOrderNbKwd'	    } )  );  
        $j++;
      }
      
    }
    if ( ref ($newOutHash) eq 'ARRAY' ){
      my @sortArray=@{ $newOutHash };   #the sort order: COVE SCORE > top stem energy > full energy > std > distance
      my @sortArray_1; my @sortArray_2;  my @sortArray_3;
      foreach my $eachUnit (@sortArray){
      	#if ($eachUnit->{'00_CoveScoKwd'}>0){ #print "a2-20180906$ecSelPro 1\$eachUnit->{'00_CoveScoKwd'}=$eachUnit->{'00_CoveScoKwd'}\n";
      	if ($eachUnit->{'00_CoveScoKwd'}>=5){
      		push @sortArray_1, $eachUnit;
      	}
      	elsif ($eachUnit->{'44_End_SECIS_dis'}>=0){
      		push @sortArray_2, $eachUnit;
      	}
      	else {                               #print "a2-20180906$ecSelPro 2\$eachUnit->{'00_CoveScoKwd'}=$eachUnit->{'00_CoveScoKwd'}\n";
      		push @sortArray_3, $eachUnit;
      	}
      }
      #print "a2-20180906$ecSelPro 00 print \@sortArray_1 dump\n"; DirFileHandle::PrintAndWarnDumper(  [ @sortArray_1 ]  );
      #print "a2-20180906$ecSelPro 00 print \@sortArray_2 dump\n"; DirFileHandle::PrintAndWarnDumper(  [ @sortArray_2 ]  );
      @sortArray_1=sort {      $b->{'00_CoveScoKwd'     } <=> $a->{'00_CoveScoKwd'     }  
      	                   ||  $a->{'06_upstemEnKwd'    } <=> $b->{'06_upstemEnKwd'    }  
      	                   ||  $a->{'05_fullstructEnKwd'} <=> $b->{'05_fullstructEnKwd'}  
      	                   ||  $b->{'01_stdNUMBER'      } <=> $a->{'01_stdNUMBER'      }  
      	                   ||  $a->{'44_End_SECIS_dis'  } <=> $b->{'44_End_SECIS_dis'  }  
      	                                                            	                     } @sortArray_1;
      @sortArray_2=sort {       
      	                       $b->{'00_CoveScoKwd'     } <=> $a->{'00_CoveScoKwd'     }
      	                   ||  $a->{'06_upstemEnKwd'    } <=> $b->{'06_upstemEnKwd'    }  
      	                   ||  $a->{'05_fullstructEnKwd'} <=> $b->{'05_fullstructEnKwd'}  
      	                   ||  $b->{'01_stdNUMBER'      } <=> $a->{'01_stdNUMBER'      }  
      	                   ||  $a->{'44_End_SECIS_dis'  } <=> $b->{'44_End_SECIS_dis'  }   
      	                                                            	                    } @sortArray_2;  
      @sortArray_3=sort {       
      	                       $b->{'00_CoveScoKwd'     } <=> $a->{'00_CoveScoKwd'     }
      	                   ||  $a->{'06_upstemEnKwd'    } <=> $b->{'06_upstemEnKwd'    }  
      	                   ||  $a->{'05_fullstructEnKwd'} <=> $b->{'05_fullstructEnKwd'}  
      	                   ||  $b->{'01_stdNUMBER'      } <=> $a->{'01_stdNUMBER'      }  
      	                   ||  $b->{'44_End_SECIS_dis'  } <=> $a->{'44_End_SECIS_dis'  }   
      	                                                            	                    } @sortArray_3;   
      my @finalSortArray; 
      for ( my $i=0; $i<@sortArray_1; $i++ ){    	       
      	if ($i<=1){
      		push @finalSortArray, $sortArray_1[$i];
      	}
      }
      for ( my $i=0; $i<@sortArray_2; $i++ ){    	       
      	if ($i<=0){
      		push @finalSortArray, $sortArray_2[$i];
      	}
      }
      for ( my $i=0; $i<@sortArray_1; $i++ ){    	       
      	if ($i>1){
      		push @finalSortArray, $sortArray_1[$i];
      	}
      }
      for ( my $i=0; $i<@sortArray_2; $i++ ){    	       
      	if ($i>0){
      		push @finalSortArray, $sortArray_2[$i];
      	}
      }
      
                    
      
      #print "a2-20180906$ecSelPro print \@sortArray_1 dump\n"; DirFileHandle::PrintAndWarnDumper(  [ @sortArray_1 ]  );
      #print "a2-20180906$ecSelPro print \@sortArray_2 dump\n"; DirFileHandle::PrintAndWarnDumper(  [ @sortArray_2 ]  );
      #print "a2-20180906$ecSelPro print \@finalSortArray dump\n"; DirFileHandle::PrintAndWarnDumper(  [ @finalSortArray ]  );
                                                  	                      	                                                            	                  
      my $newSortEdOutHash=[@finalSortArray];
      return $newSortEdOutHash;
    }
  }
  else {
  	my $dieMsg="$die_MsgHead\nThe \$SECISarray=$SECISarray should be a ARRAY ref!\n\n\n"; 
  	print $dieMsg; die $dieMsg;
  	my $caller_inform=DirFileHandle::print_SubCallerInform;
  	DieWork::Just_dieWork( "$dieMsg\n$! \n$caller_inform\n\n" );
  }

}

sub CheckTheSECISTheSameOrNot{  #  SECISwork::CheckTheSECISTheSameOrNot
	my ($SECIS_Hash1, $SECIS_Hash2, $ecSelPro)=@_;
	
	my $warnMsgBody="\nIn package  SECISwork,\tIn sub CheckTheSECISTheSameOrNot,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	
	my $finalReturn=0;
	if (  ( ref ($SECIS_Hash1) eq 'HASH' ) && ( ref ($SECIS_Hash2) eq 'HASH' )  ){                                               #print "8-1-20180908$ecSelPro\n";
		if  (   (  defined ( $SECIS_Hash1->{'02_secisDrctK'} )  ) && (  defined ( $SECIS_Hash2->{'02_secisDrctK'} )  ) && ( $SECIS_Hash1->{'02_secisDrctK'} eq $SECIS_Hash2->{'02_secisDrctK'} )   ){    #print "8-2-20180908$ecSelPro\n";
			if   (   (  defined ( $SECIS_Hash1->{'03_secisHDK'} )  ) && (  defined ( $SECIS_Hash2->{'04_secisTLK'} )  )  && (  defined ( $SECIS_Hash2->{'03_secisHDK'} )  ) && (  defined ( $SECIS_Hash2->{'04_secisTLK'} )  )   ){    #print "8-3-20180908$ecSelPro\n";
			  my $SECIS_1_Hd=$SECIS_Hash1->{'03_secisHDK'};        print "8-4-20180908$ecSelPro\$SECIS_1_Hd=\$SECIS_Hash1->{'03_secisHDK'}=$SECIS_1_Hd\n";
		  	my $SECIS_1_Tl=$SECIS_Hash1->{'04_secisTLK'};        print "8-5-20180908$ecSelPro\$SECIS_1_Tl=\$SECIS_Hash1->{'04_secisTLK'}=$SECIS_1_Tl\n";
		  	my $SECIS_1_Lt=abs ($SECIS_1_Tl-$SECIS_1_Hd)+1;
		  	my $SECIS_2_Hd=$SECIS_Hash2->{'03_secisHDK'};        print "8-6-20180908$ecSelPro\$SECIS_2_Hd=\$SECIS_Hash2->{'03_secisHDK'}=$SECIS_2_Hd\n";
		  	my $SECIS_2_Tl=$SECIS_Hash2->{'04_secisTLK'};        print "8-7-20180908$ecSelPro\$SECIS_2_Tl=\$SECIS_Hash2->{'04_secisTLK'}=$SECIS_2_Tl\n";
		  	my $SECIS_2_Lt=abs ($SECIS_2_Tl-$SECIS_2_Hd)+1;
			  
			  my ($SECIS1_ovly_rate, $SECIS2_ovly_rate)= @{ SeqSegmentsTools::Get_Overlay_rate($SECIS_1_Hd, $SECIS_1_Tl, $SECIS_2_Hd, $SECIS_2_Tl, $ecSelPro) }; print "8-8-20180908\$SECIS1_ovly_rate=\$SECIS1_ovly_rate\t\$SECIS2_ovly_rate=\$SECIS2_ovly_rate\n";
			  if (  ($SECIS1_ovly_rate >= 0.9) || ($SECIS2_ovly_rate >= 0.9)  ){
			  	$finalReturn=1;
			  }
			}
			 
		}
		  
	}
	return $finalReturn;	
}

sub GetAllSECISsgmentARRAY{  # SECISwork::GetAllSECISsgmentARRAY
	my ($SECISarray)=@_;
	
	my $warnMsgBody="\nIn package  SECISwork,\tIn sub GetAllSECISsgmentARRAY,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	
	my $outArray;
	if ( ref ($SECISarray) eq 'ARRAY' ){
    for (  my $i=0; $i < @{ $SECISarray }; $i++  ){                     
      $outArray->[$i]->[0]->{'0_SgmtHead'}=$SECISarray->[$i]->{'03_secisHDK'};
      $outArray->[$i]->[0]->{'1_SgmtTail'}=$SECISarray->[$i]->{'04_secisTLK'};
      $outArray->[$i]->[0]->{'2_SgmtZorF'}=$SECISarray->[$i]->{'02_secisDrctK'};
      $outArray->[$i]->[0]->{'3_SgmtType'}='SECIS';
      my $SgmtInfm="Cove$SECISarray->[$i]->{'00_CoveScoKwd'} $SECISarray->[$i]->{'12_SECIS_2_r1__part'}_$SECISarray->[$i]->{'14_SECIS_4_r2__part'}"; 
      
      $outArray->[$i]->[0]->{'5_SgmtInfm'}=$SgmtInfm;
    }
  }
  else {
  	my $dieMsg="$die_MsgHead\nThe \$SECISarray=$SECISarray should be a ARRAY ref!\n\n\n"; 
  	print $dieMsg; die $dieMsg;
  }
  return $outArray;
}


#my $cuttedSECISarray=SECISwork::Get_around_SECISarray( $SECISarray, $wholGeneStartPos, $wholGeneEnd__Pos, $SECIS_distance_cut_off);
sub Get_around_SECISarray{
	my ($SECISarray, $wholGeneStartPos, $wholGeneEnd__Pos, $SECIS_distance_cut_off )=@_;
	
	my $warnMsgBody="\nIn package  SECISwork,\tIn sub GetAllSECISsgmentARRAY,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $subCallereIfm=DirFileHandle::print_SubCallerInform;
	
	my $expandLength=0; $expandLength=$SECIS_distance_cut_off if (   (  defined ( $SECIS_distance_cut_off )  ) && ( $SECIS_distance_cut_off=~m/\d+/ ) && ( $SECIS_distance_cut_off > 0 )   );
	my $in_nb____big=SeqSegmentsTools::getbigOne    ($wholGeneStartPos, $wholGeneEnd__Pos); $in_nb____big+=$expandLength;
  my $in_nb__small=SeqSegmentsTools::getSmallOne  ($wholGeneStartPos, $wholGeneEnd__Pos); $in_nb__small-=$expandLength;
	
	my $outArray; my $outArrayIdx=0;
	if ( ref ($SECISarray) eq 'ARRAY' ){
    for (  my $i=0; $i < @{ $SECISarray }; $i++  ){                     
      my $secisHD=$SECISarray->[$i]->{'03_secisHDK'};
      my $secisTL=$SECISarray->[$i]->{'04_secisTLK'};
      
      my $SECIS___big=SeqSegmentsTools::getbigOne    ($secisHD, $secisTL);
      my $SECIS_small=SeqSegmentsTools::getSmallOne  ($secisHD, $secisTL);
      
      my $yes_or_no=SeqSegmentsTools::a_b_including_c_d ($in_nb__small, $in_nb____big, $SECIS_small, $SECIS___big);
      if ( $yes_or_no == 1 ){
      	$outArray->[$outArrayIdx]=$SECISarray->[$i];
      	$outArrayIdx++;
      }
      
    }
  }
  else {  
  	DieWork::Just_dieWork( $die_MsgHead."\nThe \$SECISarray=$SECISarray should be a ARRAY ref!!! \n\n$! \n$subCallereIfm\n\n" );
  }
  return $outArray;
	
	
	
}








1;

##########################################################################################################################################
# 
#my @call = caller(0);
#    print $call[3];