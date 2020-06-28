#!/usr/bin/perl -w
BEGIN{  push (@INC, 'home/fredjiang/fredProgs/PMs')  };


use strict;
use warnings;



package BuildSelCisElmentsNameIdex;

####################################################################
#
# 这个是用在为硒蛋白的合成体系的分析时，建立了各个作用因子的名字和相关顺序、简写等信息
#
#
####################################################################


sub getAllCisElmentIdxHash{
  my $CisElmHash;  #"20170218_EFsec_SelB", "20170218_PSTK", "20170218_SBP2", "20170218_SPS2", "20170218_SecS_SelA", "20170218_YbbB", "20170218_ribosomalProteinL30", "20170218_secp43", "20170218_tRNAsec"
  $CisElmHash->{'20170218_EFsec_SelB'         }->{'TreeIdx'}=0;  $CisElmHash->{'20170218_EFsec_SelB'         }->{'RealName'}='EFsec_SelB';             $CisElmHash->{'20170218_EFsec_SelB'         }->{'ShortNm'}='EFse';         
  $CisElmHash->{'20170218_PSTK'               }->{'TreeIdx'}=1;  $CisElmHash->{'20170218_PSTK'               }->{'RealName'}='PSTK';                   $CisElmHash->{'20170218_PSTK'               }->{'ShortNm'}='PSTK';               
  $CisElmHash->{'20170218_SBP2'               }->{'TreeIdx'}=2;  $CisElmHash->{'20170218_SBP2'               }->{'RealName'}='SBP2';                   $CisElmHash->{'20170218_SBP2'               }->{'ShortNm'}='SBP2';               
  $CisElmHash->{'20170218_SecS_SelA'          }->{'TreeIdx'}=3;  $CisElmHash->{'20170218_SecS_SelA'          }->{'RealName'}='SecS_SelA';              $CisElmHash->{'20170218_SecS_SelA'          }->{'ShortNm'}='SecS';          
  $CisElmHash->{'20170218_ribosomalProteinL30'}->{'TreeIdx'}=4;  $CisElmHash->{'20170218_ribosomalProteinL30'}->{'RealName'}='ribosomalProteinL30';    $CisElmHash->{'20170218_ribosomalProteinL30'}->{'ShortNm'}='rP30';
  $CisElmHash->{'20170218_secp43'             }->{'TreeIdx'}=5;  $CisElmHash->{'20170218_secp43'             }->{'RealName'}='secp43';                 $CisElmHash->{'20170218_secp43'             }->{'ShortNm'}='sc43';             
  $CisElmHash->{'20170218_SPS2'               }->{'TreeIdx'}=6;  $CisElmHash->{'20170218_SPS2'               }->{'RealName'}='SPS';                    $CisElmHash->{'20170218_SPS2'               }->{'ShortNm'}='SPS2';                
  $CisElmHash->{'20170218_YbbB'               }->{'TreeIdx'}=7;  $CisElmHash->{'20170218_YbbB'               }->{'RealName'}='YbbB';                   $CisElmHash->{'20170218_YbbB'               }->{'ShortNm'}='YbbB';               
  $CisElmHash->{'20170218_tRNAsec'            }->{'TreeIdx'}=8;  $CisElmHash->{'20170218_tRNAsec'            }->{'RealName'}='tRNAsec';                $CisElmHash->{'20170218_tRNAsec'            }->{'ShortNm'}='tRse';            
  $CisElmHash->{'20180620__SELENOF'           }->{'TreeIdx'}=9;  $CisElmHash->{'20180620__SELENOF'           }->{'RealName'}='SELENOF';                $CisElmHash->{'20180620__SELENOF'           }->{'ShortNm'}='SelF';            
  $CisElmHash->{'20180702__SELENOM'           }->{'TreeIdx'}=9;  $CisElmHash->{'20180702__SELENOM'           }->{'RealName'}='SELENOM';                $CisElmHash->{'20180702__SELENOM'           }->{'ShortNm'}='SelM';            
  
  return  $CisElmHash;
}

sub BuildCisElmentArray{
  my $CisElmHash;
  $CisElmHash->[0]='20170218_EFsec_SelB'         ;
  $CisElmHash->[1]='20170218_PSTK'               ;
  $CisElmHash->[2]='20170218_SBP2'               ;
  $CisElmHash->[3]='20170218_SecS_SelA'          ;
  $CisElmHash->[4]='20170218_ribosomalProteinL30';
  $CisElmHash->[5]='20170218_secp43'             ;
  $CisElmHash->[6]='20170218_SPS2'               ;
  $CisElmHash->[7]='20170218_YbbB'               ;
  #$CisElmHash->[8]='20170218_tRNAsec'            ;
  
                                                 
  return $CisElmHash;
}

sub BuildCisElmOrderHash{
  my $CisElmHash;  
  $CisElmHash->{'20170218_EFsec_SelB'         }=0;         
  $CisElmHash->{'20170218_PSTK'               }=1;               
  $CisElmHash->{'20170218_SBP2'               }=2;               
  $CisElmHash->{'20170218_SecS_SelA'          }=3;          
  $CisElmHash->{'20170218_ribosomalProteinL30'}=4;  
  $CisElmHash->{'20170218_secp43'             }=5;             
  $CisElmHash->{'20170218_SPS2'               }=6;                
  $CisElmHash->{'20170218_YbbB'               }=7; 
  $CisElmHash->{'20170218_tRNAsec'            }=8;              
  
  return  $CisElmHash;
}

sub BuildCisElmShortNameHash{
  my $CisElmHash;  
  $CisElmHash->{'20170218_EFsec_SelB'         }='EFse';         
  $CisElmHash->{'20170218_PSTK'               }='PSTK';               
  $CisElmHash->{'20170218_SBP2'               }='SBP2';               
  $CisElmHash->{'20170218_SecS_SelA'          }='SecS';          
  $CisElmHash->{'20170218_ribosomalProteinL30'}='rP30';  
  $CisElmHash->{'20170218_secp43'             }='sc43';             
  $CisElmHash->{'20170218_SPS2'               }='SPS2';                
  $CisElmHash->{'20170218_tRNAsec'            }='tRse';            
  $CisElmHash->{'20170218_YbbB'               }='YbbB';               
  
  return  $CisElmHash;
}

sub BuildCisElmRealNameHash{
  my $CisElmHash;  
  $CisElmHash->{'20170218_EFsec_SelB'         }='EFsec_SelB';         
  $CisElmHash->{'20170218_PSTK'               }='PSTK';               
  $CisElmHash->{'20170218_SBP2'               }='SBP2';               
  $CisElmHash->{'20170218_SecS_SelA'          }='SecS_SelA';          
  $CisElmHash->{'20170218_ribosomalProteinL30'}='ribosomalProteinL30';
  $CisElmHash->{'20170218_secp43'             }='secp43';             
  $CisElmHash->{'20170218_SPS2'               }='SPS2';                
  $CisElmHash->{'20170218_tRNAsec'            }='tRNAsec';            
  $CisElmHash->{'20170218_YbbB'               }='YbbB';               
  
  return  $CisElmHash;
}

sub BuildCisElmDomianIDHash{
  my $CisElmHash;  
  $CisElmHash->{'20170218_EFsec_SelB'         }->{'IPR000795'}=1;  $CisElmHash->{'20170218_EFsec_SelB'         }->{'IPR009000'}=1;  $CisElmHash->{'20170218_EFsec_SelB'         }->{'IPR027417'}=1;
  $CisElmHash->{'20170218_PSTK'               }->{'IPR027417'}=1;                                 #IPR027417
  
  $CisElmHash->{'20170218_SBP2'               }->{'IPR004038'}=1;                                 #IPR004038
  $CisElmHash->{'20170218_ribosomalProteinL30'}->{'IPR004038'}=1;                                 #IPR004038
  $CisElmHash->{'20170218_SecS_SelA'          }->{'IPR015421'}=1;  $CisElmHash->{'20170218_SecS_SelA'          }->{'IPR015424'}=1; 
  $CisElmHash->{'20170218_secp43'             }->{'IPR000504'}=1;  $CisElmHash->{'20170218_secp43'             }->{'IPR012677'}=1;              
  $CisElmHash->{'20170218_SPS2'               }->{'IPR010918'}=1;  $CisElmHash->{'20170218_SPS2'               }->{'IPR016188'}=1;                
           
  $CisElmHash->{'20170218_YbbB'               }->{'IPR001763'}=1;               
  
  return  $CisElmHash;
}

1;

