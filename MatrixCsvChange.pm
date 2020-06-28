#!/usr/bin/perl -w
BEGIN{  push (@INC, '/home/fredjiang/fredProgs/PMs')  };

use strict;
use warnings;
use Excel::Writer::XLSX;
use Excel::Writer::XLSX::Utility;
use DieWork;
use InFileHandle;



use DirFileHandle;


package MatrixCsvChange;

sub Change_IP_for_share_in_hyperLink_in_ptout_excelTable{  #MatrixCsvChange::Change_IP_for_share_in_hyperLink_in_ptout_excelTable ($inFile, $org_IP, $new_IP);
	my ($inFile, $org_IP, $new_IP)=@_;
	
	my $warnMsgBody="\nIn package  MatrixCsvChange,\tIn sub Change_IP_for_share_in_hyperLink_in_ptout_excelTable,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	DirFileHandle::print_SubCallerInform($org_IP);
	DirFileHandle::print_SubCallerInform($new_IP);
  $org_IP=~s/\./\\\./g; #warn $org_IP;
  #$new_IP=~s/\./\\\./g; warn $new_IP;
	if (   (  defined ( $inFile )  ) && ( $inFile=~m/\S+/ )   ){
		my $outString=InFileHandle::readAllfileIntoAstring($inFile);
		$outString=~s/$org_IP/$new_IP/g;
		InFileHandle::PrintStringIntoFile($inFile, $outString);
	}
	else{
		DieWork::Just_dieWork( $die_MsgHead."\n\$inFile=$inFile is not right!! $!\n\n\n".$caller_inform );
	}
	
	
}

sub CheckIP_format{   #   my $checkEd=DirFileHandle::print_SubCallerInform($org_IP);
	
	my ($org_IP)=@_;
	
	my $warnMsgBody="\nIn package  MatrixCsvChange,\tIn sub CheckIP_format,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	my $checkEd=0;
	if ($org_IP=~m/^(\d+)\.(\d+)\.(\d+)\.(\d+)$/){ 
	  if (  ( (0<=$1) && ($1<=255) ) && ( (0<=$2) && ($2<=255) ) && ( (0<=$3) && ($3<=255) ) && ( (0<=$4) && ($4<=255) )  ){	
	  	$checkEd=1;  	
	  }
	  else{
	  	DieWork::Just_dieWork( $die_MsgHead."\n\$org_IP=$org_IP is not right!! $!\n\n\n".$caller_inform );
	  }
	}
	else{
    DieWork::Just_dieWork( $die_MsgHead."\n\$org_IP=$org_IP is not right!! $!\n\n\n".$caller_inform );
  }
	return $checkEd;
}

sub Change00Number_toBack{  #my $outNB= MatrixCsvChange::Change00Number_toBack( $inNB );
	my ($inNB)=@_;
	
	my $warnMsgBody="\nIn package  MatrixCsvChange,\tIn sub Change00Number_toBack,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	my $outNB;
	
	if ( $inNB=~m/^0*([123456789]*\d*\d)$/){
		$outNB=$1;		
	}
	else{
		DieWork::Just_dieWork( $die_MsgHead."\n\$inNB=$inNB is not right!! $!\n\n\n".$caller_inform );
	}
	
	
}

sub GetBigestNb_from_Array{   # my $biggestNb=MatrixCsvChange::GetBigestNb_from_Array( $inArray );
	my ($inArray)=@_;
	
	my $warnMsgBody="\nIn package  MatrixCsvChange,\tIn sub GetBigestNb_from_Array,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	my $biggestNb=-999999999999999999999999999999999;
	if  (   (  defined ( $inArray )  ) && (  ref ( $inArray ) eq 'ARRAY'  )   ){
	  
	  my $checkTime=0;
	  foreach my $eachNB (  @{ $inArray }  ){
	  	if  (   (  defined ( $eachNB )  ) && (  $eachNB=~m/(-)?\d+(\.\d+)?/  )   ){
	  	  if ( $eachNB > $biggestNb ){
	  		  $biggestNb=$eachNB;
	  		  $checkTime++;
	  	  }
	  	}
	  	else{
	  		DieWork::Just_dieWork( $die_MsgHead."\n\$eachNB=$eachNB Should be a number!! $!\n\n\n".$caller_inform );
	  	}
	  }
	  if ( $checkTime==0 ){
	  	DieWork::Just_dieWork( $die_MsgHead."\n\$checkTime=$checkTime Should be a bumber > 0, means there is at lest one number bigger than original $biggestNb !! $!\n\n\n".$caller_inform );
	  }
	}
	else{
		DieWork::Just_dieWork( $die_MsgHead."\n\$inArray=$inArray Should be a ARRAY ref!! $!\n\n\n".$caller_inform );
	}
	return $biggestNb;
}

# Get the digits of a input number
sub Get_Weishu{   # my $realWeishu=MatrixCsvChange::Get_Weishu($inputNb);
	my ($inputNb)=@_;
	my ($warnMsgBody, $warnMsgHead, $die_MsgHead, $caller_inform)=@{ DieWork::BuildWarnDieHeadInform( 'MatrixCsvChange', 'Get_Weishu' ) };
  my $weiSHu=&log10 ( $inputNb   );      
	my $weiShuZheng=int ( $weiSHu );     
	my $realWeishu=$weiShuZheng+1;      
  return $realWeishu;
}

sub SprintfKeyHead{    #my $showKey= MatrixCsvChange::SprintfKeyHead( $bigestNb, $prNB );
	my ($bigestNb, $prNB)=@_;
	
	my $warnMsgBody="\nIn package  MatrixCsvChange,\tIn sub SprintfKeyHead,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	#my $weiSHu=&log10 ( $bigestNb   );       	my $weiShuZheng=int ( $weiSHu );     	my $realWeishu=$weiShuZheng+1;    
	my $realWeishu=MatrixCsvChange::Get_Weishu($bigestNb);        
	my $showKey=sprintf ("%0${realWeishu}d", $prNB); 
	#return $realWeishu;
	return $showKey;
	
}

sub SprintfKeyHead_inPut_weishu{    #my $showKey= MatrixCsvChange::SprintfKeyHead_inPut_weishu( $weiSHu, $prNB );
	my ($weiSHu, $prNB)=@_;
	
	my $warnMsgBody="\nIn package  MatrixCsvChange,\tIn sub SprintfKeyHead_inPut_weishu,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	#my $weiSHu=&log10 ( $bigestNb   );      
	my $weiShuZheng=int ( $weiSHu );     
	my $realWeishu=$weiShuZheng+1;            
	my $showKey=sprintf ("%0${realWeishu}d", $prNB); 
	#return $realWeishu;
	return $showKey;
	
}


sub log10{
	my ($inNB)=@_;
	my $outNB=log($inNB)/log(10);
	return $outNB;
}

sub Change2dHashIntoCsvString{
  my ($inPutHash, $reverse)=@_;
  my $outString;
  
  
  my $allKey0Hash;
  my $allKey1Hash;
  
  my $inHash=$inPutHash;  #DirFileHandle::PrintDumper("inHash.hash", $inHash)  if ( ref ( $inHash ) eq 'HASH');;

  if (  ( ref ($inHash) ) eq 'HASH' ){
    foreach my $keyLev_0 (    sort { $a cmp $b } (   keys (  %{ $inHash } )   )    ){   print "In package MatrixCsvChange sub Change2dHashIntoCsvString   \$keyLev_0=$keyLev_0\n";
    	my $valLev_0=$inHash->{$keyLev_0};                                                print "In package MatrixCsvChange sub Change2dHashIntoCsvString    \$valLev_0=$valLev_0\n";
    	my $refLev_0=ref ($valLev_0);                                                     print "In package MatrixCsvChange sub Change2dHashIntoCsvString    \$refLev_0=$refLev_0\n\n";  
    	
    	$allKey0Hash->{$keyLev_0}=1;
    	
      if ( $refLev_0 eq 'HASH' ){
        foreach my $keyLev_1 (    sort { $a cmp $b } (   keys (  %{ $valLev_0 } )   )    ){   print "In package MatrixCsvChange sub Change2dHashIntoCsvString      \$keyLev_1=$keyLev_1\n";
        	my $valLev_1=$valLev_0->{$keyLev_1};                                                print "In package MatrixCsvChange sub Change2dHashIntoCsvString      \$valLev_1=$valLev_1\n";                
        	
        	$allKey1Hash->{$keyLev_1}=1;
        	
          
        }
      }
      
    }
  }
  
   my $colKeyHash=$allKey0Hash;
   my $rowKeyHash=$allKey1Hash;
 
  
  if ( defined ($reverse) ){
    if ($reverse == 1){
    	$colKeyHash=$allKey1Hash;
      $rowKeyHash=$allKey0Hash;
    }
  }
  
 
  #写下第一行  的字段名
  $outString.="TopLeft\t";   
  if (  ( ref ($rowKeyHash) ) eq 'HASH' ){
    foreach my $keyLev_row (    sort { $a cmp $b } (   keys (  %{ $rowKeyHash } )   )    ){   print "In package MatrixCsvChange sub Change2dHashIntoCsvString    \$keyLev_row=$keyLev_row\n";    	
    	$outString.="$keyLev_row\t"; 
    }
  }
  $outString.="\n"; 
  
  #开始写行标识  及 每行的内容
  if (  ( ref ($allKey1Hash) ) eq 'HASH' ){
    foreach my $keyLev_col (    sort { $a cmp $b } (   keys (  %{ $colKeyHash } )   )    ){   print "In package MatrixCsvChange sub Change2dHashIntoCsvString    \$keyLev_col=$keyLev_col\n";
    	
    	$outString.="$keyLev_col\t";  #写下每行标识
    	
      if  (  ( ref ($rowKeyHash) ) eq 'HASH' ){
        foreach my $keyLev_row (    sort { $a cmp $b } (   keys (  %{ $rowKeyHash } )   )    ){   print "In package MatrixCsvChange sub Change2dHashIntoCsvString      \$keyLev_row=$keyLev_row\n";
        	
        	my $Key_0=$keyLev_col; 
        	my $key_1=$keyLev_row;
        	if ( defined ($reverse) ){
            if ($reverse == 1){
            	$Key_0=$keyLev_row;
              $key_1=$keyLev_col;
            }
          }
        	
        	
        	
        if ( defined ($inPutHash->{$Key_0}->{$key_1}) ){
          $outString.="$inPutHash->{$Key_0}->{$key_1}\t"; 
        }
        else {
          $outString.="0\t";
        }
        	
        	
        	
        	
          
        }
      }
      $outString.="\n"; 
      
    }
  }

  return $outString;

}

sub PrintExcelTableFileFromHASH{  #              my $outString= MatrixCsvChange::PrintExcelTableFileFromHASH( $outPrintFile, $inPutHash, $reverse, $sortHashForKey0, $sortKeyForKey0, $sortHashForKey1, $sortKeyForKey1, $donot_show_fisrtKeyCloumn );
	my ($outPrintFile, $inPutHash, $reverse, $sortHashForKey0, $sortKeyForKey0, $sortHashForKey1, $sortKeyForKey1, $donot_show_fisrtKeyCloumn)=@_;
	   #1               2          3          4                5                6                  7               8 
	my $warnMsgBody="\nIn package  MatrixCsvChange,\tIn sub PrintExcelTableFileFromHASH,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	my $outString=MatrixCsvChange::Change2dHashIntoCsvString_withSortFunc( $inPutHash, $reverse, $sortHashForKey0, $sortKeyForKey0, $sortHashForKey1, $sortKeyForKey1, $donot_show_fisrtKeyCloumn );
	
	open  ( OUTFILE, ">$outPrintFile" ) or die "cannot create \$outPrintFile=$outPrintFile :$!\n\n\n";
  print   OUTFILE $outString ;
  close ( OUTFILE );

	return $outString;
}


sub Change2dHashIntoCsvString_withSortFunc{  #  my $outString=MatrixCsvChange::Change2dHashIntoCsvString_withSortFunc( $inPutHash, $reverse, $sortHashForKey0, $sortKeyForKey0, $sortHashForKey1, $sortKeyForKey1, $donot_show_fisrtKeyCloumn );
  my ($inPutHash, $reverse, $sortHashForKey0, $sortKeyForKey0, $sortHashForKey1, $sortKeyForKey1, $donot_show_fisrtKeyCloumn)=@_;  
  
  my $warnMsgBody="\nIn package  MatrixCsvChange,\tIn sub Change2dHashIntoCsvString_withSortFunc,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
	my $caller_inform=DirFileHandle::print_SubCallerInform;
  
  #print " In package MatrixCsvChange sub Change2dHashIntoCsvString_withSortFunc   ";
  if ( defined ($inPutHash                 ) ){print "\$inPutHash=$inPutHash\t"                       ;}
  if ( defined ($reverse                   ) ){print "\$reverse=$reverse\t"                           ;}
  if ( defined ($sortHashForKey0           ) ){print "\$sortHashForKey0=$sortHashForKey0\t"           ;}
  if ( defined ($sortKeyForKey0            ) ){print "\$sortKeyForKey0=$reverse\t"                    ;}
  if ( defined ($sortHashForKey1           ) ){print "\$sortHashForKey1=$sortHashForKey1\t"           ;}
  if ( defined ($sortKeyForKey1            ) ){print "\$sortKeyForKey1=$sortKeyForKey1\t"             ;}
  if ( defined ($donot_show_fisrtKeyCloumn ) ){print "\$donot_show_fisrtKeyCloumn=$sortKeyForKey1\t"  ;}
  print "\n\n";
  
  my $outString;
  
  
  my $allKey0Hash;
  my $allKey1Hash;
  
  my $inHash=$inPutHash;  #DirFileHandle::PrintDumper("inHash.hash", $inHash)  if ( ref ( $inHash ) eq 'HASH');

  if (  ( ref ($inHash) ) eq 'HASH' ){
    foreach my $keyLev_0 (    sort { $a cmp $b } (   keys (  %{ $inHash } )   )    ){   #print "In package MatrixCsvChange sub Change2dHashIntoCsvString_withSortFunc   \$keyLev_0=$keyLev_0\n";
    	my $valLev_0=$inHash->{$keyLev_0};                                                #print "In package MatrixCsvChange sub Change2dHashIntoCsvString_withSortFunc    \$valLev_0=$valLev_0\n";
    	my $refLev_0=ref ($valLev_0);                                                     #print "In package MatrixCsvChange sub Change2dHashIntoCsvString_withSortFunc    \$refLev_0=$refLev_0\n\n";  
    	
    	$allKey0Hash->{$keyLev_0}=1;
    	
      if ( $refLev_0 eq 'HASH' ){
        foreach my $keyLev_1 (    sort { $a cmp $b } (   keys (  %{ $valLev_0 } )   )    ){   #print "In package MatrixCsvChange sub Change2dHashIntoCsvString_withSortFunc      \$keyLev_1=$keyLev_1\n";
        	my $valLev_1=$valLev_0->{$keyLev_1};                                                #print "In package MatrixCsvChange sub Change2dHashIntoCsvString_withSortFunc      \$valLev_1=$valLev_1\n";                
        	
        	$allKey1Hash->{$keyLev_1}=1;
        	
          
        }
      }
      
    }
  }
  
  my $sortHashForKDef0=0;  
  my $sortHashForKDef1=0;
  if ( ( defined ($sortHashForKey0) ) && ( defined ($sortKeyForKey0) ) ){                        print "0000001 In package MatrixCsvChange sub Change2dHashIntoCsvString_withSortFunc   \$sortHashForKey0=$sortHashForKey0\n";
    if (   (  ( ref ($sortHashForKey0) ) eq 'HASH'  ) && ($sortKeyForKey0=~m/\S+/)    ){         print "0000002 In package MatrixCsvChange sub Change2dHashIntoCsvString_withSortFunc   \$sortKeyForKey0=$sortKeyForKey0\n";
    	$sortHashForKDef0=1;
    }
  }
  if ( ( defined ($sortHashForKey1) ) && ( defined ($sortKeyForKey1) ) ){                        print "0000003 In package MatrixCsvChange sub Change2dHashIntoCsvString_withSortFunc   \$sortHashForKey1=$sortHashForKey1\n";
    if (   (  ( ref ($sortHashForKey1) ) eq 'HASH'  ) && ($sortKeyForKey1=~m/\S+/)    ){         print "0000004 In package MatrixCsvChange sub Change2dHashIntoCsvString_withSortFunc   \$sortKeyForKey1=$sortKeyForKey1\n";
    	$sortHashForKDef1=1;
    }
  }
  
  
  
  
  my $colKeyHash;  my $colSortDef=0;  my $colSortHash; my $colSortKey;
  my $rowKeyHash;  my $rowSortDef=0;  my $rowSortHash; my $rowSortKey;
  my $realReverse=0;
  if ( defined ($reverse) ){
    if ($reverse == 1){
    	$realReverse=1;                                                                           print "0000005 In package MatrixCsvChange sub Change2dHashIntoCsvString_withSortFunc   \$realReverse=$realReverse\n";
    	
    } 
  }
  else {
    print "0000006 In package MatrixCsvChange sub Change2dHashIntoCsvString_withSortFunc   \$realReverse=$realReverse\n";
  }
  
  if ($realReverse){  #这是反过来的 第一个Key作 行， 第二个key作列
    $colKeyHash=$allKey1Hash;    if ($sortHashForKDef1) { $colSortDef=1; $colSortHash=$sortHashForKey1; $colSortKey=$sortKeyForKey1; print "11In package MatrixCsvChange sub Change2dHashIntoCsvString_withSortFunc \realReverse=realReverse\t\$colSortDef=$colSortDef\t\$colSortKey=$colSortKey\n";}
    $rowKeyHash=$allKey0Hash;    if ($sortHashForKDef0) { $rowSortDef=1; $rowSortHash=$sortHashForKey0; $rowSortKey=$sortKeyForKey0; print "22In package MatrixCsvChange sub Change2dHashIntoCsvString_withSortFunc \realReverse=realReverse\t\$colSortDef=$colSortDef\t\$colSortKey=$colSortKey\n";}
  }
  else{                #这是正的，第一个Key作 列， 第二个key作行
  	$colKeyHash=$allKey0Hash;    if ($sortHashForKDef0) { $colSortDef=1; $colSortHash=$sortHashForKey0; $colSortKey=$sortKeyForKey0; print "33In package MatrixCsvChange sub Change2dHashIntoCsvString_withSortFunc \realReverse=realReverse\t\$colSortDef=$colSortDef\t\$colSortKey=$colSortKey\n";}
    $rowKeyHash=$allKey1Hash;    if ($sortHashForKDef1) { $rowSortDef=1; $rowSortHash=$sortHashForKey1; $rowSortKey=$sortKeyForKey1; print "44In package MatrixCsvChange sub Change2dHashIntoCsvString_withSortFunc \realReverse=realReverse\t\$colSortDef=$colSortDef\t\$colSortKey=$colSortKey\n";}
  }
  
  
 
  my @colKeyAr;  if (  ( ref ($colKeyHash) ) eq 'HASH' ){   	@colKeyAr=(   keys (  %{ $colKeyHash } )   );    }
  my @rowKeyAr;  if (  ( ref ($rowKeyHash) ) eq 'HASH' ){   	@rowKeyAr=(   keys (  %{ $rowKeyHash } )   );    }
  my @colSortKeyAr;  	
  if ( $colSortDef ) { @colSortKeyAr = sort { $colSortHash->{$a}->{$colSortKey} <=> $colSortHash->{$b}->{$colSortKey} } @colKeyAr;} 
  else               { 
  	if   (  ( defined ($colKeyAr[0]) ) && ($colKeyAr[0]=~m/^\s*(\+|-)?\d+\.?\d*\s*$/)  ){@colSortKeyAr = sort { $a <=> $b } @colKeyAr;}
  	else                                          {@colSortKeyAr = sort { $a cmp $b } @colKeyAr;}
  }
  	
  my @rowSortKeyAr;  	
  if ( $rowSortDef ) { @rowSortKeyAr = sort { $rowSortHash->{$a}->{$rowSortKey} <=> $rowSortHash->{$b}->{$rowSortKey} } @rowKeyAr;} 
  else               { 
  	if   (  ( defined ($colKeyAr[0]) ) && ($rowKeyAr[0]=~m/^\s*(\+|-)?\d+\.?\d*\s*$/)  ){@rowSortKeyAr = sort { $a <=> $b } @rowKeyAr;}
  	else                                          {@rowSortKeyAr = sort { $a cmp $b } @rowKeyAr;}
  }
  
  #写下第一行  的字段名
  if (   (  defined ( $donot_show_fisrtKeyCloumn )  ) && ( $donot_show_fisrtKeyCloumn=~m/\d+/ ) && ( $donot_show_fisrtKeyCloumn == 1 )   ){ 	  }
  else                  {   $outString.="TopLeft\t";   }
  
  if (  ( ref ($rowKeyHash) ) eq 'HASH' ){
    my $idx=0; my $lastIdx=@rowSortKeyAr;
    foreach my $keyLev_row (@rowSortKeyAr){   print "In package MatrixCsvChange sub Change2dHashIntoCsvString_withSortFunc    \$keyLev_row=$keyLev_row\n";    	
    	if (  ($idx+1) == $lastIdx   ){
    		$outString.="$keyLev_row"; 
    	}
    	else{
    		$outString.="$keyLev_row\t"; 
    	}
    	$idx++;
    }
  }
  $outString.="\n"; 
  
  
  #开始写行标识  及 每行的内容
  if (  ( ref ($colKeyHash) ) eq 'HASH' ){ 	
    
    foreach my $keyLev_col (@colSortKeyAr){   print "In package MatrixCsvChange sub Change2dHashIntoCsvString_withSortFunc    \$keyLev_col=$keyLev_col\n";
    	
    	#写下每行标识
    	if (   (  defined ( $donot_show_fisrtKeyCloumn )  ) && ( $donot_show_fisrtKeyCloumn=~m/\d+/ ) && ( $donot_show_fisrtKeyCloumn == 1 )   ){ 	  }
    	else{     $outString.="$keyLev_col\t";     }   
    	
    	
      if  (  ( ref ($rowKeyHash) ) eq 'HASH' ){
      	my $idx=0; my $lastIdx=@rowSortKeyAr;
      	foreach my $keyLev_row (@rowSortKeyAr){   
      		#print "In package MatrixCsvChange sub Change2dHashIntoCsvString_withSortFunc      \$keyLev_row=$keyLev_row\t\$rowSortHash->{$keyLev_row}->{$rowSortKey}=$rowSortHash->{$keyLev_row}->{$rowSortKey}\n";
        	
        	my $Key_0=$keyLev_col; 
        	my $key_1=$keyLev_row;
        	if ( defined ($reverse) ){
            if ($reverse == 1){
            	$Key_0=$keyLev_row;
              $key_1=$keyLev_col;
            }
          }
        	
        	
          if (  ($idx+1) == $lastIdx   ){	
            if ( defined ($inPutHash->{$Key_0}->{$key_1}) ){
            	my $tempString=MatrixCsvChange::Cut_long_String_to_fit_excel_cell_maximal_length_limit( $inPutHash->{$Key_0}->{$key_1} );
              $outString.="$tempString"; 
            }
            else {
              $outString.="";
            }
        	}
        	else{	
            if ( defined ($inPutHash->{$Key_0}->{$key_1}) ){
            	my $tempString=MatrixCsvChange::Cut_long_String_to_fit_excel_cell_maximal_length_limit( $inPutHash->{$Key_0}->{$key_1} );
              $outString.="$tempString\t"; 
            }
            else {
              $outString.="\t";
            }
        	}
        	
        	$idx++;
          
        }
      }
      $outString.="\n"; 
      
    }
  }

  return $outString;

}

#my $outString=MatrixCsvChange::Cut_long_String_to_fit_excel_cell_maximal_length_limit($inString);
sub Cut_long_String_to_fit_excel_cell_maximal_length_limit{
	my ($inString)=@_;
	
	#my $warnMsgBody="\nIn package  MatrixCsvChange,\tIn sub Cut_long_String_to_fit_excel_cell_maximal_length_limit,\n\n";	
  #my $warnMsgHead="\n\n\n$warnMsgBody";
	#my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
	#my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	my $outString=$inString;
	my $maximal_cell_limit=32767;
	my $cut_subStrinLength=$maximal_cell_limit-1; if ( $cut_subStrinLength < 0) { $cut_subStrinLength=0; }
	
	my $StringLength=length( $inString );
	if ( $StringLength > $maximal_cell_limit ){
		$outString=substr( $inString, 0, $cut_subStrinLength );
	}
	return $outString;
	
	
}




1;