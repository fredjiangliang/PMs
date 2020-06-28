#!/usr/bin/perl -w

BEGIN{  push (@INC, '/home/fredjiang/fredProgs/PMs')  };

use strict;
use warnings;

package ForeachHash;




sub foreachHash{
  
  my ($inHash)=@_;

  if (  ( ref ($inHash) ) eq 'HASH' ){
    foreach my $keyLev_0 (    sort { $a cmp $b } (   keys (  %{ $inHash } )   )    ){   print "    \$keyLev_0=$keyLev_0\n";
    	my $valLev_0=$inHash->{$keyLev_0};                                                print "    \$valLev_0=$valLev_0\n";
    	my $refLev_0=ref ($valLev_0);                                                     print "    \$refLev_0=$refLev_0\n\n";  
    	
      if ( $refLev_0 eq 'HASH' ){
        foreach my $keyLev_1 (    sort { $a cmp $b } (   keys (  %{ $valLev_0 } )   )    ){   print "      \$keyLev_1=$keyLev_1\n";
        	my $valLev_1=$valLev_0->{$keyLev_1};                                                print "      \$valLev_1=$valLev_1\n";                
        	my $refLev_1=ref ($valLev_1);                                                       print "      \$refLev_1=$refLev_1\n\n";
        	
      	  if ( $refLev_1 eq 'HASH' ){
            foreach my $keyLev_2 (    sort { $a cmp $b } (   keys (  %{ $valLev_1 } )   )    ){   print "        \$keyLev_2=$keyLev_2\n";
            	my $valLev_2=$valLev_1->{$keyLev_2};                                                print "        \$valLev_2=$valLev_2\n";
            	my $refLev_2=ref ($valLev_2);                                                       print "        \$refLev_2=$refLev_2\n\n";
     
            	if ( $refLev_2 eq 'HASH' ){
                foreach my $keyLev_3 (    sort { $a cmp $b } (   keys (  %{ $valLev_2 } )   )    ){   print "          \$keyLev_3=$keyLev_3\n";
                	my $valLev_3=$valLev_2->{$keyLev_3};                                                print "          \$valLev_3=$valLev_3\n";
                	my $refLev_3=ref ($valLev_3);                                                       print "          \$refLev_3=$refLev_3\n\n";
       
            	    if ( $refLev_3 eq 'HASH' ){
                    foreach my $keyLev_4 (    sort { $a cmp $b } (   keys (  %{ $valLev_3 } )   )    ){   print "            \$keyLev_4=$keyLev_4\n";
                    	my $valLev_4=$valLev_3->{$keyLev_4};                                                print "            \$valLev_4=$valLev_4\n";
                    	my $refLev_4=ref ($valLev_4);                                                       print "            \$refLev_4=$refLev_4\n\n";
       
            	        if ( $refLev_4 eq 'HASH' ){
                        foreach my $keyLev_5 (    sort { $a cmp $b } (   keys (  %{ $valLev_4 } )   )    ){   print "              \$keyLev_5=$keyLev_5\n";
                        	my $valLev_5=$valLev_4->{$keyLev_5};                                                print "              \$valLev_5=$valLev_5\n";
                        	my $refLev_5=ref ($valLev_5);                                                       print "              \$refLev_5=$refLev_5\n\n";
        
                          if ( $refLev_5 eq 'HASH' ){
                            foreach my $keyLev_6 (    sort { $a cmp $b } (   keys (  %{ $valLev_5 } )   )    ){   print "              \$keyLev_6=$keyLev_6\n";
                        	    my $valLev_6=$valLev_5->{$keyLev_6};                                                print "              \$valLev_6=$valLev_6\n";
                        	    my $refLev_6=ref ($valLev_6);                                                       print "              \$refLev_6=$refLev_6\n\n\n\n";
                            
                              if ( $refLev_6 eq 'HASH' ){
                                foreach my $keyLev_7 (    sort { $a cmp $b } (   keys (  %{ $valLev_6 } )   )    ){   print "              \$keyLev_7=$keyLev_7\n";
                        	        my $valLev_7=$valLev_6->{$keyLev_7};                                                print "              \$valLev_7=$valLev_7\n";
                        	        my $refLev_7=ref ($valLev_7);                                                       print "              \$refLev_7=$refLev_7\n\n\n\n";
                        	     
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
    }
  }

}


sub buildOrderHash{
  my ($orgHash, $orderType)=@_;
  my $outHash;  my $outIdx=0;
  if ($orderType eq 'number'){
  	if (  ( ref ( $orgHash ) ) eq 'HASH' ){
  		foreach my $keyHere (    sort { $a <=> $b } (   keys (  %{ $orgHash }  )   )    ){
  			$outHash->{$keyHere}=$outIdx;
  			$outIdx++;
  		}
  	}
  }
  elsif ($orderType eq 'string'){
  	if (  ( ref ( $orgHash ) ) eq 'HASH' ){
  		foreach my $keyHere (    sort { $a cmp $b } (   keys (  %{ $orgHash }  )   )    ){
  			$outHash->{$keyHere}=$outIdx;
  			$outIdx++;
  		}
  	}
  }
  return $outHash;
}


1;