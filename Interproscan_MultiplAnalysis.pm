#!/usr/bin/perl -w
BEGIN{  push (@INC, '/home/fredjiang/fredProgs/PMs')  };

use strict;
use warnings;
use Interproscan;
use DieWork;
use Storable ;
use Storable qw(nstore store_fd nstore_fd freeze thaw dclone);


package Interproscan_MultiplAnalysis;

my $ipsTypes; my $ipsTypeNB;    #��interProScan�а����ĸ�����ͬ�� �ؼ��֣��������黯������ų� hash�� key��val
@{$ipsTypes} =('Family',    'Domain',    'SITE',     'Repeat',     'SUPERFAMILY',       'GENE3D',      'Pfam',      'PANTHER',      'PRINTS',      'PROSITE profiles',      'PROSITE patterns',       'PIRSF',       'UNKNOWN',      'SMART',        'HAMAP',       'COILS',       'Biological Process',      'Molecular Function',       'Celluar Component'    );
%{$ipsTypeNB}=('Family'=>0, 'Domain'=>1, 'SITE'=> 2, 'Repeat'=> 3, 'SUPERFAMILY' => 4,  'GENE3D' => 5, 'Pfam' => 6, 'PANTHER' => 7, 'PRINTS' => 8, 'PROSITE profiles' => 9, 'PROSITE patterns' => 10, 'PIRSF' => 11, 'UNKNOWN' => 12,'SMART' => 13,  'HAMAP' => 14, 'COILS' => 15, 'Biological Process' => 16,'Molecular Function' => 17, 'Celluar Component'=> 18);


#Interproscan_MultiplAnalysis::Build_excel_for_mutipl_IPShash($inHash, $inFastaHASH, $excel_file_path, $sheet_1, $sheet_2, $middle_clustalw_or_not, $middle_clustalRSTdir, $middle_hperLkCLUSTEDrstDIR);
sub Build_excel_for_mutipl_IPShash{
	my ($inHash, $inFastaHASH, $excel_file_path, $sheet_1, $sheet_2, $middle_clustalw_or_not, $middle_clustalRSTdir, $middle_hperLkCLUSTEDrstDIR)=@_;
	
	my ($warnMsgBody, $warnMsgHead, $die_MsgHead, $caller_inform)=@{ DieWork::BuildWarnDieHeadInform( 'Interproscan_MultiplAnalysis', 'Build_excel_for_mutipl_IPShash' ) };

  DieWork::Check_Hash_or_DIE( $inHash,                      "\$inHash",         $die_MsgHead, $caller_inform  );
  DieWork::Check_Hash_or_DIE( $inFastaHASH,                 "\$inFastaHASH",     $die_MsgHead, $caller_inform  );
  DieWork::Check_DfdNoEmptString_or_DIE( $excel_file_path, "\$excel_file_path", $die_MsgHead, $caller_inform  );
  
  my $real_sheet_1='st1_no_order';
  my $real_sheet_2='st2_in_order';
  if ( DieWork::Check_DfdNoEmptString_or_NOT( $sheet_1 )  ){ $real_sheet_1=$sheet_1;  }
  if ( DieWork::Check_DfdNoEmptString_or_NOT( $sheet_2 )  ){ $real_sheet_2=$sheet_2;  }
  
    #���ڿ�ʼ�����µ�sheet��Ӧ��hash
  my $newSheetHash_KeyedByFMT;        #�����µ�hash�����������ڴ�ӡ����sheet�еģ���fmtΪkey����hash
  my $inPtFmtHash;                    #format key hash
  my $middle_multCol_startAtNB=5;     #���������Ϊ�ˣ����½��� html���� sheet�У��������뼸�У�����һЩ��Ϣ�ġ� �������5����ʾ����5�С�
  
  
  my $colWidthHash;                   #����һ��hash������ȷ�������еĿ��
  my $Out_3Darray;                    #����һ����άhash���飬��һ��key�����кţ��ڶ���key�����кţ�������key��IPScan�������id��ָ�����html�ļ�������� ��ά�������ݽṹ
  my $Out_ClassedExcelHash;           #����һ����άhash���飬��һ��key�ǣ�family��domain��site�ȵ����������ͣ����ڶ���key��IPScan��id��ֵ���Ƕ����Ϣ��ɵ�����ĵ�ַ
  
  my $IPS_idx=0;
  ($Out_3Darray,$Out_ClassedExcelHash)=@{ Interproscan_MultiplAnalysis::getIformFromHtmlHash(
  	                                                              $inHash, 
    	                                                            #$oldRow, 
    	                                                            #$gmeColForType, 
    	                                                            #$rowColHash->{$oldRow}->{0}->[0], 
    	                                                            $IPS_idx,
    	                                                            $Out_3Darray,
    	                                                            #$Out_ClassedExcelHash,
    	                                                            $ipsTypeNB, 
    	                                                            $middle_multCol_startAtNB 
    	                                                            ) 
    	                                  };
  
  
  
  
  
  #DirFileHandle::PrintDumper_plReadFile('Out_3Darray.hsh', $Out_3Darray)               if ( ref ($Out_3Darray         ) eq 'HASH' ) ;
  #DirFileHandle::PrintDumper_plReadFile('Out_ClassedExcelHash.hsh', $Out_ClassedExcelHash)      if ( ref ($Out_ClassedExcelHash) eq 'HASH'  ) ;
  
  my ($out_clust3Darray, $out_subTYPEhash, $out_subCOLhash); #�½����������������沽���У���һ���� html�ĸ�����Ϣ����������� ���������ĺ��壬����sub makePThashFromHtmls�в鿴
  #���������ϲ����ɵ� �������ݽṹ������ html���������Ϣ ��ɵ� sheet�� hash����hash��format��Ϊkey
  ( $newSheetHash_KeyedByFMT,
    $inPtFmtHash,
    $colWidthHash, 
    $out_clust3Darray, 
    $out_subTYPEhash, 
    $out_subCOLhash )
  = @{ &makePThashFromHtmls( $Out_3Darray,
  	                         $Out_ClassedExcelHash,
  	                         $ipsTypes, 
  	                         $newSheetHash_KeyedByFMT,
  	                         #$inPtFmtHash,
  	                         $colWidthHash, 
  	                         $middle_multCol_startAtNB 
  	                       ) 
  	 };
  
  #DirFileHandle::PrintDumper_plReadFile('newSheetHash_KeyedByFMT.hsh', $newSheetHash_KeyedByFMT);      #DirFileHandle::PrintDumper('inPtFmtHash.hsh', $inPtFmtHash)      if ( ref ($inPtFmtHash) eq 'HASH'  ) ;  #DirFileHandle::PrintDumper_plReadFile('colWidthHash.hsh', $colWidthHash);       #DirFileHandle::PrintDumper_plReadFile('out_clust3Darray.hsh', $out_clust3Darray);       #DirFileHandle::PrintDumper_plReadFile('out_subTYPEhash.hsh', $out_subTYPEhash);        #DirFileHandle::PrintDumper_plReadFile('out_subCOLhash.hsh', $out_subCOLhash) ;     
  
  my $outXlsBook= Excel::Writer::XLSX->new($excel_file_path);                      #����ļ��Ķ�����
  my $inPPSheet_1=$outXlsBook->add_worksheet($real_sheet_1);                                                 
  #&PtExcelBasedOnFmtHash($ptOutBook, $inPPSheet_1, $newSheetHash_KeyedByFMT ,$inPtFmtHash );           
  
  ExcelHandle::PtExcelBasedOnFmtHash($outXlsBook, $inPPSheet_1, $newSheetHash_KeyedByFMT, $inPtFmtHash);
                                                                                                                                                                                     #print "\n\n333333333333333333333333333333 print_all_sub_array (\$colWidthHash)\n\cl\cl";&print_all_sub_array ($colWidthHash);  print "\cl\cl\n\n";
  setColumnFormat_2 ($inPPSheet_1,$colWidthHash);
  
  DieWork::Check_DfdNoEmptString_or_DIE( $middle_clustalRSTdir,       "\$middle_clustalRSTdir", $die_MsgHead, $caller_inform  );
  DieWork::Check_DfdNoEmptString_or_DIE( $middle_hperLkCLUSTEDrstDIR, "\$middle_hperLkCLUSTEDrstDIR", $die_MsgHead, $caller_inform  );
  
  if ( DieWork::Check_INTNB_equal_to_aNUMBER_or_NOT( $middle_clustalw_or_not, 1 ) ){  	
  }
  else{
  	#$middle_clustalRSTdir      ='middle_clustalRSTdir';
  	#$middle_hperLkCLUSTEDrstDIR='middle_hperLkCLUSTEDrstDIR';
  }
  
  
    #��������html��ͬ���������򣬶Բ�ͬ������clustalw���������������ս��
  my $SortedHTMLhash=&sortClustalByIPScanResult ( 
                                                  $newSheetHash_KeyedByFMT,  
                                                  $out_clust3Darray, 
                                                  $out_subTYPEhash, 
                                                  $out_subCOLhash, 
                                                  $inFastaHASH,
                                                  $middle_clustalRSTdir,
                                                  $middle_hperLkCLUSTEDrstDIR, 
                                                  $middle_clustalw_or_not, 
                                                  $ipsTypeNB, 
                                                  $middle_multCol_startAtNB
                                                );
  
                                                                                                                                                                                     #warn "Finished the &sortClustalByIPScanResult(\$newSheetHash_KeyedByFMT=$newSheetHash_KeyedByFMT,  \$out_clust3Darray=$out_clust3Darray, \$out_subTYPEhash=$out_subTYPEhash, \$out_subCOLhash=$out_subCOLhash, \$middle_clustalRSTdir=$middle_clustalRSTdir,\$middle_hperLkCLUSTEDrstDIR=$middle_hperLkCLUSTEDrstDIR, \$middle_clustalw_or_not=$middle_clustalw_or_not);      \n";
  my $inPPSheet_2=$outXlsBook->add_worksheet($real_sheet_2); 
  ExcelHandle::PtExcelBasedOnFmtHash($outXlsBook, $inPPSheet_2, $SortedHTMLhash ,$inPtFmtHash ); 
  &setColumnFormat_2 ($inPPSheet_2,$colWidthHash);

  
  
}






##sub8.4          getIformFromHtmlHash{   #��ȡhtml�ļ��еĸ�����Ϣ��������� ��Ҫ�� ���ݽṹ
sub getIformFromHtmlHash{   #��ȡhtml�ļ��еĸ�����Ϣ��������� ��Ҫ�� ���ݽṹ
	my ($eachHtmlHash, 
	    #$orgRowNb, 
	    #$g_m_e_typeNb, 
	    #$specNm, 
	    $RowStartNB,
	    $the3Darray, 
	    #$ClassedExcelHash, 
	    $inTypeArrayAddr, 
	    $multCol_startAtNB_1 )=@_;
	
	
	my ($warnMsgBody, $warnMsgHead, $die_MsgHead, $caller_inform)=@{ DieWork::BuildWarnDieHeadInform( 'Interproscan_MultiplAnalysis', 'getIformFromHtmlHash' ) };
  
  DieWork::Check_Hash_or_DIE( $eachHtmlHash,         "\$eachHtmlHash",          $die_MsgHead, $caller_inform  );
  #DieWork::Check_Hash_or_DIE( $ClassedExcelHash,     "\$ClassedExcelHash",      $die_MsgHead, $caller_inform  );
  DieWork::Check_Hash_or_DIE( $inTypeArrayAddr,      "\$inTypeArrayAddr",       $die_MsgHead, $caller_inform  );
  
  DieWork::Check_DfdNoEmptNUMBER_or_DIE( $RowStartNB,          "\$RowStartNB",          $die_MsgHead, $caller_inform  );
  DieWork::Check_DfdNoEmptNUMBER_or_DIE( $multCol_startAtNB_1, "\$multCol_startAtNB_1", $die_MsgHead, $caller_inform  );
  
	#my $multCol_startAtNB_1=5;   #��������֣������������ɵ� װ��html�������ݵ� sheet�У�ǰ�����ж��������������������ݵġ���������5�С�
	# =('Family'=>0, 'Domain'=>1, 'SITE'=> 2, 'Repeat'=> 3, 'SUPERFAMILY' => 4,  'GENE3D' => 5, 'Pfam' => 6, 'PANTHER' => 7, 'PRINTS' => 8, 'PROSITE profiles' => 9, 'PROSITE patterns' => 10, 'PIRSF' => 11, 'UNKNOWN' => 12,'SMART' => 13,  'HAMAP' => 14, 'COILS' => 15, 'Biological Process' => 16,'Molecular Function' => 17, 'Celluar Component'=> 18);
  my %typeColNb=%{$inTypeArrayAddr};   
  foreach my $fmdomKey (keys %typeColNb){  $typeColNb{$fmdomKey}+=$multCol_startAtNB_1;}  #����װ�и�����ͬ�������͵���ŵ�hash��key�����飬��ÿ��key����Ԥ���������������������������$multCol_startAtNB_1=5;

  my $ClassedExcelHash; #����һ����άhash���飬��һ��key�ǣ�family��domain��site�ȵ����������ͣ����ڶ���key��IPScan��id��ֵ���Ƕ����Ϣ��ɵ�����ĵ�ַ

	my $newRowNb=$RowStartNB;
  foreach my $htmlIDkey (    sort { $a cmp $b } (  keys ( %{ $eachHtmlHash } )  )    ){    #$htmlIDkey ���ﻹ��html�Ǹ��ļ���protein ��name
    
    if (  DieWork::Check_Hash_or_NOT( $eachHtmlHash->{$htmlIDkey} )  ){
      foreach my $subTpKey (  keys ( %{ $eachHtmlHash->{$htmlIDkey} } )  ){    #$subTpKey ������ length, Go term prediction, Domain, no IPR, Family, Repeat, SITE
      	printf "\$subTpKey=%50s\n", $subTpKey;
      
        ################################################################################################################################################################################################################################################################################################################################################
      	if    ($subTpKey =~ m/^(Family|Domain|SITE|Repeat)$/){
      		foreach my $IPRidKey (    sort { $a cmp $b } (  keys ( %{ $eachHtmlHash->{$htmlIDkey}->{$subTpKey} } )  )   ){  #$IPRidKey ������ ÿ��PIR��id ��IPR001650��IPR011545��IPR014001��IPR027417        
      			#print "\$eachHtmlHash->{$htmlIDkey}->{$subTpKey}->{$IPRidKey}\n"; $eachHtmlHash->{$htmlIDkey}->{$subTpKey}->{$IPRidKey} &print_all_sub_array( $eachHtmlHash->{$htmlIDkey}->{$subTpKey}->{$IPRidKey} ); print "\cl\n\n";
            
            my $IPR_defin    =                   $eachHtmlHash->{$htmlIDkey}->{$subTpKey}->{$IPRidKey}->{'Definition'} ;
            my $IPR_webURL   =                   $eachHtmlHash->{$htmlIDkey}->{$subTpKey}->{$IPRidKey}->{'Url'} ;
            DieWork::Print_and_warn( "20191104-0-0-0 \$eachHtmlHash->{$htmlIDkey}->{$subTpKey}->{$IPRidKey}->{'subArray'}=$eachHtmlHash->{$htmlIDkey}->{$subTpKey}->{$IPRidKey}->{'subArray'}\n" );
            my $IPR_subArray = Storable::dclone( $eachHtmlHash->{$htmlIDkey}->{$subTpKey}->{$IPRidKey}->{'subArray'} );
            
            my $typeColNbHere=$typeColNb{$subTpKey};
            #my $eachCellInfom=[ $orgRowNb, $g_m_e_typeNb, $specNm, $htmlIDkey, $IPRidKey, $IPR_defin, $IPR_webURL, $newRowNb, $typeColNbHere ];
            my $eachCellInfom=[ $htmlIDkey, $IPRidKey, $IPR_defin, $IPR_webURL, $newRowNb, $typeColNbHere ];
            $the3Darray->{$newRowNb}->{$typeColNbHere}->{$IPRidKey}=$htmlIDkey;    print "111 \$the3Darray->{$newRowNb}->{$typeColNbHere}->{$IPRidKey}=\$htmlIDkey=$htmlIDkey\n";
            
                  
      		  if (  defined ( $ClassedExcelHash->{$subTpKey}->{$IPRidKey} )  ){
      		    push @{ $ClassedExcelHash->{$subTpKey}->{$IPRidKey} }, $eachCellInfom;   print "\@\$eachCellInfom below\n"; map {printf "%10s\t",$_; } @{$eachCellInfom}; print "\n", "\$ClassedExcelHash->{$subTpKey}->{$IPRidKey}=$ClassedExcelHash->{$subTpKey}->{$IPRidKey}\n\n"
      		  }
      		  else { $ClassedExcelHash->{$subTpKey}->{$IPRidKey}=[$eachCellInfom]; }
      		  
      		  ######################################################################################
      		  #########################������domain��family���������ڲ�#############################
      		  #print "\n\n\cl\&print_all_sub_array(\$IPR_subArray)\n";&print_all_sub_array($IPR_subArray); print "\cl\n\n";
      		  
      		  foreach my $WithIPR_idKey_Hash ( @{ $IPR_subArray }  ){          # print "\$WithIPR_idKey_Hash=$WithIPR_idKey_Hash\n"; print "\n\n\cl\&print_all_sub_array(\$WithIPR_idKey_Hash)\n";&print_all_sub_array($WithIPR_idKey_Hash); print "\cl\n\n";
            
              foreach my $WithIPR_idKey (    sort { $a cmp $b } (   keys %{ $WithIPR_idKey_Hash } )   ){
              
                #my ( $with_IPR_webURL, $with_IPR_defin, $with_IPR_subArray)=@{ $WithIPR_idKey_Hash->{$WithIPR_idKey} }; $with_IPR_defin=~s/^\s*\S*\s*\(//; $with_IPR_defin=~s/\)\s*$//;
                
                my $with_IPR_defin    = $WithIPR_idKey_Hash->{$WithIPR_idKey}->{'Definition'} ;
                my $with_IPR_webURL   = $WithIPR_idKey_Hash->{$WithIPR_idKey}->{'Url'} ;
                DieWork::Print_and_warn( "20191104-0-0-1 \$WithIPR_idKey_Hash->{$WithIPR_idKey}->{'subSubArray'}=$WithIPR_idKey_Hash->{$WithIPR_idKey}->{'subSubArray'}\n" );
            
                my $with_IPR_subArray = Storable::dclone( $WithIPR_idKey_Hash->{$WithIPR_idKey}->{'subSubArray'} );
            
                #DirFileHandle::PrintAndWarnDumper ($with_IPR_subArray, "\n\n20191104-0-0-1  \$with_IPR_subArray=$with_IPR_subArray:\n");
                
                my $with_IPR_type=$with_IPR_subArray->[0]->{'DomainDataBase'};                                                                             print "\$with_IPR_type=\$with_IPR_subArray->[0]->{'DomainDataBase'}=$with_IPR_subArray->[0]->{'DomainDataBase'}\n";
                #$familyType{$with_IPR_type}=1;
                my $with_IPR_typeColNbHere=$typeColNb{$with_IPR_type};
                
                #my $with_IPR_eachCellInfom=[ $orgRowNb, $g_m_e_typeNb, $specNm, $htmlIDkey, $WithIPR_idKey, $with_IPR_defin, $with_IPR_webURL, $newRowNb, $with_IPR_typeColNbHere ];  #print "\n\n\cl\&print_all_sub_array(\$with_IPR_eachCellInfom)\n";&print_all_sub_array($with_IPR_eachCellInfom); print "\cl\n\n";
                my $with_IPR_eachCellInfom=[ $htmlIDkey, $WithIPR_idKey, $with_IPR_defin, $with_IPR_webURL, $newRowNb, $with_IPR_typeColNbHere ];  #print "\n\n\cl\&print_all_sub_array(\$with_IPR_eachCellInfom)\n";&print_all_sub_array($with_IPR_eachCellInfom); print "\cl\n\n";
                $the3Darray->{$newRowNb}->{$with_IPR_typeColNbHere}->{$WithIPR_idKey}=$htmlIDkey;  print "222 \$the3Darray->{$newRowNb}->{$with_IPR_typeColNbHere}->{$WithIPR_idKey}=\$htmlIDkey=$htmlIDkey\n";
                
      		      if (  defined ( $ClassedExcelHash->{$with_IPR_type}->{$WithIPR_idKey} )  ){
      		        push @{ $ClassedExcelHash->{$with_IPR_type}->{$WithIPR_idKey} }, $with_IPR_eachCellInfom;                                       print "\@\$with_IPR_eachCellInfom below\n"; map {printf "%10s\n",$_; } @{$with_IPR_eachCellInfom}; print "\n", "\$ClassedExcelHash->{$with_IPR_type}->{$WithIPR_idKey}=$ClassedExcelHash->{$with_IPR_type}->{$WithIPR_idKey}\n\n";
      		      }
      		      else { $ClassedExcelHash->{$with_IPR_type}->{$WithIPR_idKey}=[$with_IPR_eachCellInfom];                                          print "\@\$with_IPR_eachCellInfom below\n"; map {printf "%10s\n",$_; } @{$with_IPR_eachCellInfom}; print "\n", "\$ClassedExcelHash->{$with_IPR_type}->{$WithIPR_idKey}=$ClassedExcelHash->{$with_IPR_type}->{$WithIPR_idKey}\n\n";}
      		    }
      		  }
      		  ########################������domain��family���������ڲ�##############################
      		  ######################################################################################
      		  
      		  
      		}
      	}
      	elsif ($subTpKey =~ m/^no IPR$/){
      		
      		foreach my $noIPRidKey (    sort { $a cmp $b } (  keys ( %{ $eachHtmlHash->{$htmlIDkey}->{$subTpKey} } )  )    ){  #          print "\$eachHtmlHash->{$htmlIDkey}->{$subTpKey}->{$noIPRidKey}= $eachHtmlHash->{$htmlIDkey}->{$subTpKey}->{$noIPRidKey}\n"; &print_all_sub_array( $eachHtmlHash->{$htmlIDkey}->{$subTpKey}->{$noIPRidKey} ); print "\cl\n\n";
            
            #my ($noIPR_defin, $noIPR_webURL, $noIPR_subSubArray)=@{ $eachHtmlHash->{$htmlIDkey}->{$subTpKey}->{$noIPRidKey} }; 
            
            my $noIPR_defin        = $eachHtmlHash->{$htmlIDkey}->{$subTpKey}->{$noIPRidKey}->{'Definition'} ;
            my $noIPR_webURL       = $eachHtmlHash->{$htmlIDkey}->{$subTpKey}->{$noIPRidKey}->{'Url'} ;
            my $noIPR_Db_type      = $eachHtmlHash->{$htmlIDkey}->{$subTpKey}->{$noIPRidKey}->{'DomainDataBase'} ;
            DieWork::Print_and_warn( "20191104-0-0-2 \$eachHtmlHash->{$htmlIDkey}->{$subTpKey}->{$noIPRidKey}->{'subNoPirArray'}=$eachHtmlHash->{$htmlIDkey}->{$subTpKey}->{$noIPRidKey}->{'subNoPirArray'}\n" );
            my $noIPR_subSubArray  = Storable::dclone( $eachHtmlHash->{$htmlIDkey}->{$subTpKey}->{$noIPRidKey}->{'subNoPirArray'} );
            
            #DirFileHandle::PrintAndWarnDumper ($noIPR_subSubArray, "\n\n20191104-0-02  \$noIPR_subSubArray=$noIPR_subSubArray:\n");
                
            
            $noIPR_defin=~s/^\s*\(//; $noIPR_defin=~s/\)\s*$//;
            
            my $noIPRtype=$noIPR_subSubArray->[0]->{'DomainDataBase'};                                                                             print "\$noIPRtype=\$noIPR_subSubArray->[0]->{'DomainDataBase'}=$noIPRtype=$noIPR_subSubArray->[0]->{'DomainDataBase'}\n";
            #$familyType{$noIPRtype}=1;
            my $typeColNbHere=$typeColNb{$noIPRtype};
            
            #my $eachCellInfom=[ $orgRowNb, $g_m_e_typeNb, $specNm, $htmlIDkey, $noIPRidKey, $noIPR_defin, $noIPR_webURL, $newRowNb, $typeColNbHere ];
            my $eachCellInfom=[  $htmlIDkey, $noIPRidKey, $noIPR_defin, $noIPR_webURL, $newRowNb, $typeColNbHere ];
            $the3Darray->{$newRowNb}->{$typeColNbHere}->{$noIPRidKey}=$htmlIDkey;  print "222 \$the3Darray->{$newRowNb}->{$typeColNbHere}->{$noIPRidKey}=\$htmlIDkey=$htmlIDkey\n";
      
      		  if (  defined ( $ClassedExcelHash->{$noIPRtype}->{$noIPRidKey} )  ){
      		    push @{ $ClassedExcelHash->{$noIPRtype}->{$noIPRidKey} }, $eachCellInfom;                                     print "\@\$eachCellInfom below\n"; map {printf "%10s\n",$_; } @{$eachCellInfom}; print "\n", "\$ClassedExcelHash->{$noIPRtype}->{$noIPRidKey}=$ClassedExcelHash->{$noIPRtype}->{$noIPRidKey}\n\n";
      		  }
      		  else { $ClassedExcelHash->{$noIPRtype}->{$noIPRidKey}=[$eachCellInfom];                                          print "\@\$eachCellInfom below\n"; map {printf "%10s\n",$_; } @{$eachCellInfom}; print "\n", "\$ClassedExcelHash->{$noIPRtype}->{$noIPRidKey}=$ClassedExcelHash->{$noIPRtype}->{$noIPRidKey}\n\n";}
      		  
      		}
      		
      		
      		
      		
      	}         #GO term prediction
      	elsif ($subTpKey =~ m/GO term prediction/){ print "GOhere:\$eachHtmlHash->{$htmlIDkey}->{$subTpKey}=$eachHtmlHash->{$htmlIDkey}->{$subTpKey}\n";
      		
      		foreach my $subGOkey (    sort { $a cmp $b } (  keys ( %{ $eachHtmlHash->{$htmlIDkey}->{$subTpKey} } )  )    ){  #$subGOkey ������ GO����������          print "GOhere:\$eachHtmlHash->{$htmlIDkey}->{$subTpKey}->{$subGOkey}=$eachHtmlHash->{$htmlIDkey}->{$subTpKey}->{$subGOkey}\n";  &print_all_sub_array( $eachHtmlHash->{$htmlIDkey}->{$subTpKey}->{$subGOkey} ); print "\cl\n\n";
            
            foreach my $GoIdKey (    sort { $a cmp $b } ( keys ( %{ $eachHtmlHash->{$htmlIDkey}->{$subTpKey}->{$subGOkey} } )  )    ){
      
              #my ($GO_defin, $GO_webURL)=@{ $eachHtmlHash->{$htmlIDkey}->{$subTpKey}->{$subGOkey}->{$GoIdKey} };                       printf "GO_define:%30s\tURL:%50s\n", $GO_defin, $GO_webURL;
              my $GO_defin        = $eachHtmlHash->{$htmlIDkey}->{$subTpKey}->{$subGOkey}->{$GoIdKey}->{'Definition'} ;
              my $GO_webURL       = $eachHtmlHash->{$htmlIDkey}->{$subTpKey}->{$subGOkey}->{$GoIdKey}->{'Url'} ;
            
              
              my $typeColNbHere=$typeColNb{$subGOkey};                                                                                 print "\$typeColNbHere=\$typeColNb{$subGOkey}=$typeColNbHere\n";
              #my $eachCellInfom=[ $orgRowNb, $g_m_e_typeNb, $specNm, $htmlIDkey, $GoIdKey, $GO_defin, $GO_webURL, $newRowNb, $typeColNbHere ]; 
              my $eachCellInfom=[  $htmlIDkey, $GoIdKey, $GO_defin, $GO_webURL, $newRowNb, $typeColNbHere ]; 
              $the3Darray->{$newRowNb}->{$typeColNbHere}->{$GoIdKey}=$htmlIDkey;   print "333 \$the3Darray->{$newRowNb}->{$typeColNbHere}->{$GoIdKey}=\$htmlIDkey=$htmlIDkey\n";
              
      		    if (  defined ( $ClassedExcelHash->{$subGOkey}->{$GoIdKey} )  ){
      		      push @{ $ClassedExcelHash->{$subGOkey}->{$GoIdKey} }, $eachCellInfom;                                                print "\@\$eachCellInfom below\n"; map {printf "%10s\t",$_; } @{$eachCellInfom}; print "\n", "\$ClassedExcelHash->{$subGOkey}->{$GoIdKey}=$ClassedExcelHash->{$subGOkey}->{$GoIdKey}\n\n";
      		    }
      		    else { $ClassedExcelHash->{$subGOkey}->{$GoIdKey}=[$eachCellInfom];                                                    print "\@\$eachCellInfom below\n"; map {printf "%10s\t",$_; } @{$eachCellInfom}; print "\n", "\$ClassedExcelHash->{$subGOkey}->{$GoIdKey}=$ClassedExcelHash->{$subGOkey}->{$GoIdKey}\n\n"; }
      		  
      		  }
      		
      	  }
      	}
      	
      ################################################################################################################################################################################################################################################################################################################################################	
      	
      } 
      
    }
    $newRowNb++;
  }
  my $optHere=[$the3Darray, $ClassedExcelHash];
  return $optHere;
  
}

##sub8.5          subGroupSort3Darray {   
sub subGroupSort3Darray {  
	print "\n\n\n\cl\n\ninside of \&subGroupSort3Darray \n ";
  my ($in3Darray)=@_;
  
  #DirFileHandle::PrintAndWarnDumper ($in3Darray, "\n\n201911041431-0-0-0  \$in3Darray=$in3Darray:\n");
  my $bigestRowNb=0; my $biggestColNb=0;  
  foreach my $eachRowKey (keys (%{$in3Darray})){
  	if ($bigestRowNb < $eachRowKey){$bigestRowNb=$eachRowKey;}
  	foreach my $eachColKey (keys (%{$in3Darray->{$eachRowKey}})){
  	  if (   (  DieWork::Check_DfdNoEmptString_or_NOT( $biggestColNb )  ) && ($biggestColNb < $eachColKey)   ){$biggestColNb=$eachColKey;}
  	}    
  }                                                              print "\n\$bigestRowNb= $bigestRowNb\n"; print  "\n\$biggestColNb=$biggestColNb\n\n";
  
  my $subGroupEdArray;
  for (my $i=0; $i<=$biggestColNb; $i++){           #����                 
  	my $mergHash;   printf "colNbs \$i=\t%10d\n", $i;
    for (my $j=0; $j<=$bigestRowNb; $j++){              #����                          #my @idGroup=(keys(%{$in3Darray->{$j}->{$i}})); #���2d�ĵ������ж��ٸ���ͬ��id�������һ��@idGroup
    	my $idSubGroupHash;  my $idSubGroupRowNbHash;  my $sigOrGrpEd=0;   my $gotVal=0;    printf "rowNbs \$j=\t%20d\n", $j;
    	if (  (defined($mergHash->{$i}->{$j})) && ( $mergHash->{$i}->{$j} == 1 )  ){   #��������������Ѿ����ֵ�ĳһ�� �������ˣ�����Ҫ�����ͱ����Ŀ�Ƚ���
    	  print "\$mergHash->{$i}->{$j}=\$mergHash->{$i}->{\$j}=$mergHash->{$i}->{$j}\n";
    	}
    	else {

      	for (my $k=($j+1); $k<=$bigestRowNb; $k++){                          #����,����� һ�������� �к��еıȽ��е� ���Ƚϵ��Ǹ���
      		printf "Iner rowNbs \$k=\t%20d\n", $k;
      	  if (  (defined($mergHash->{$i}->{$k})) && ( $mergHash->{$i}->{$k} == 1 )  ){ #��������������Ѿ����ֵ�ĳһ�� �������ˣ�����Ҫ�ٺ����Ƚ���
      	    #print "\$mergHash->{$i}->{\$k}=\$mergHash->{$i}->{$k}=$mergHash->{$i}->{$k}\n";
      	    MKONE: foreach my $idKeyHere ( keys ( %{ $in3Darray->{$j}->{$i} } )  ){      printf "\n\$idKeyHere in \$in3Darray->{$j}->{$i}:\t%10s\n", $idKeyHere;# warn "\$idKeyHere=$idKeyHere\n";# print "\n\n\&print_all_sub_array(\$in3Darray->{$j}->{$i})\n";&print_all_sub_array($in3Darray->{$j}->{$i}); print "\n\n";
      	  		$gotVal=1;  last MKONE;
      	  	}
      	  }
      	  else {     
      	  	#my @idKeyArr=keys %{ $in3Darray->{$j}->{$i} }; if (@idKeyArr>0){  print "\n\nShowing \@idKeyArr :\cl\n"; map {printf "%5s\t", $_} @idKeyArr; print "\n\cl\@idKeyArr end!\n\n"; }
      	  	foreach my $idKeyHere (    sort { $a cmp $b } ( keys ( %{ $in3Darray->{$j}->{$i} } )   )     ){      printf "\n\$idKeyHere in \$in3Darray->{$j}->{$i}:\t%10s\n", $idKeyHere;# warn "\$idKeyHere=$idKeyHere\n";# print "\n\n\&print_all_sub_array(\$in3Darray->{$j}->{$i})\n";&print_all_sub_array($in3Darray->{$j}->{$i}); print "\n\n";
      	  		$gotVal=1;
      	  		if (  (defined ( $in3Darray->{$k}->{$i}->{$idKeyHere} ) ) && ($in3Darray->{$k}->{$i}->{$idKeyHere}=~m/\S+/)  ){ printf "Found this \$idKeyHere in \$in3Darray->{$k}->{$i}:\t%10s\n\n", $idKeyHere;
      	  		  foreach my $in2idKey ( keys( %{ $in3Darray->{$j}->{$i} } ) ){ $idSubGroupHash->{$in2idKey}=1; }
      	  		  foreach my $in3idKey ( keys( %{ $in3Darray->{$k}->{$i} } ) ){ $idSubGroupHash->{$in3idKey}=1; }
      	  		  $idSubGroupRowNbHash->{$j}=1; $idSubGroupRowNbHash->{$k}=1;
      	  		  $mergHash->{$i}->{$k}=1;  $mergHash->{$i}->{$j}=1;  print "Defined Here:\$mergHash->{$i}->{\$k}=\$mergHash->{$i}->{$k}=\$mergHash->{$i}->{\$j}=\$mergHash->{$i}->{$j}=1";
      	  		  $sigOrGrpEd=1;
      	  		  print "\$sigOrGrpEd=$sigOrGrpEd\t\$gotVal=$gotVal\n";
      	  		}
      	  	}
      	  }
      	} print "\$bigestRowNb=$bigestRowNb\n\n";
      	
      	if (  ($sigOrGrpEd==0) && ( $gotVal==1 )  ){
      		foreach my $in2idKey ( keys( %{ $in3Darray->{$j}->{$i} } ) ){ $idSubGroupHash->{$in2idKey}=1; }
      		$idSubGroupRowNbHash->{$j}=1; 
      		$mergHash->{$i}->{$j}=1;  print "Defined Here:\$mergHash->{$i}->{\$j}=\$mergHash->{$i}->{$j}=1";
      	  #print "\$gotVal=$gotVal\n";  print "\n\n\&print_all_sub_array(\$in3Darray->{$j}->{$i})\n";&print_all_sub_array($in3Darray->{$j}->{$i}); print "\n\n";
      	}
      	
      	if ($gotVal==1){
      	  #print "\n\n\cl\&print_all_sub_array(\$idSubGroupHash)\n";&print_all_sub_array($idSubGroupHash); print "\n\n";
      	  #print "\n\n\&print_all_sub_array(\$idSubGroupRowNbHash)\n";&print_all_sub_array($idSubGroupRowNbHash); print "\cl\n\n";
       	  push @{$subGroupEdArray->{$i}}, [ $idSubGroupHash , $idSubGroupRowNbHash  ]; #�ѷֳɵ�С�����Ϣ����������ЩС���row���Ž�һ����ַ
      	}
  
    	}
    }
  }
  #print "\n\n\cl\cl\&print_all_sub_array(\$subGroupEdArray)\n";&print_all_sub_array($subGroupEdArray); print "\cl\cl\n\n";
  
  #DirFileHandle::PrintAndWarnDumper ($subGroupEdArray, "\n\n201911041431-0-0-1  \$subGroupEdArray=$subGroupEdArray:\n");
  
  my $new3Darray_include_subs;  my $bigestColNbHash;       print "\n\cl\n\n\nStartOf buid \$new3Darray_include_subs and \$bigestColNbHash\n\n\n";
  foreach my $colNb (  sort {$a<=>$b} (keys (%{ $subGroupEdArray }) )  ) {      printf "\$colNbs in \$subGroupEdArray:\t%30d\n", $colNb; #    $colNb����������
  	my $bigestColNb=0;
    for (my $i=0; $i< @{ $subGroupEdArray->{$colNb} }; $i++){    printf "%-50s\t%15d\n", "Find Big Group: \$i in \$subGroupEdArray->{$colNb}:", $i;  #  $i������ ����С���Ԫ������
      my $eachIdGroupSize=(  keys %{ $subGroupEdArray->{$colNb}->[$i]->[0] }  ); map {printf "%10s\t", $_;} (  keys %{ $subGroupEdArray->{$colNb}->[$i]->[0] }  ); print "\n";
      if ($bigestColNb<$eachIdGroupSize){$bigestColNb=$eachIdGroupSize; printf "New \$bigestColNb=\t%40d\n", $bigestColNb;}
    }
    $bigestColNbHash->{$colNb}=$bigestColNb;  print "\n\$bigestColNbHash->{$colNb}=\$bigestColNb=$bigestColNb\n\n\n";
    
    for (my $i=0; $i< @{ $subGroupEdArray->{$colNb} }; $i++){ printf "%-50s\t%15d\n", "\$i in \$subGroupEdArray->{$colNb}:", $i;  #  $i������ ����С���Ԫ������      
      my $IdWholeGroupAddr=[ (  keys %{ $subGroupEdArray->{$colNb}->[$i]->[0] }  ) ];  #print "\n\n\&print_all_sub_array(\$subGroupEdArray->{$colNb}->[$i])\n";&print_all_sub_array($subGroupEdArray->{$colNb}->[$i]); print "\n\n";
      #print "\n\n\&print_all_sub_array(\$IdWholeGroupAddr)\n"; &print_all_sub_array($IdWholeGroupAddr); print "\n\n";      
      foreach my $rowNb (   sort {$a <=> $b} ( keys %{ $subGroupEdArray->{$colNb}->[$i]->[1] } )   ){    printf "\$rowNb in \$subGroupEdArray->{$colNb}->[$i]->[1]:\t%30d\n", $rowNb; # 
        foreach my $typeID (  keys %{ $in3Darray->{$rowNb}->{$colNb} }  ){  printf "\$typeID in \$in3Darray->{$rowNb}->{$colNb}:\t%35s\n", $typeID;
        	
        	my $IDinOne3dUnitAddr=[ (  keys %{ $in3Darray->{$rowNb}->{$colNb} }  ) ];      	
        	my $TargetPoint_and_ExtendPoints=&GetTargetPoint_and_ExtendPoints($bigestColNb, $IdWholeGroupAddr, $IDinOne3dUnitAddr,$typeID);
          $new3Darray_include_subs->{$rowNb}->{$colNb}->{$typeID}=[ $bigestColNb,  $TargetPoint_and_ExtendPoints->[0], $TargetPoint_and_ExtendPoints->[1] ];
          print "\$new3Darray_include_subs->{$rowNb}->{$colNb}->{$typeID}=$new3Darray_include_subs->{$rowNb}->{$colNb}->{$typeID}\n";
          #print "\n\n\cl\cl\n\n\&print_all_sub_array(\$new3Darray_include_subs->{$rowNb}->{$colNb}->{$typeID})\n"; &print_all_sub_array($new3Darray_include_subs->{$rowNb}->{$colNb}->{$typeID}); print "\n\n\cl\cl\n\n";      
        }
      }
    }
    
  }
  return ([ $new3Darray_include_subs, $bigestColNbHash ]);
  
##sub8.5.1          GetTargetPoint_and_ExtendPoints{
  sub GetTargetPoint_and_ExtendPoints{
    my ($ExtTotalColNb, $WholeGroupIds, $subGroupIds, $findPointID)=@_;  #printf "\n\n\n\cl\nInside of \&GetTargetPoint_and_ExtendPoints\n\n%20s\t%20s\t%20s\t%20s\n\n",$ExtTotalColNb, $WholeGroupIds, $subGroupIds, $findPointID;
    #print "\n\n\&print_all_sub_array(\$WholeGroupIds)\n";&print_all_sub_array($WholeGroupIds); print "\n\n";
    #print "\n\n\&print_all_sub_array(\$subGroupIds)\n";&print_all_sub_array($subGroupIds); print "\n\n";
    my @SbWholeIdGroup=@{ $WholeGroupIds };    @SbWholeIdGroup=sort {$a cmp $b} @SbWholeIdGroup;    #print "\n\n Showing the \@SbWholeIdGroup:\n\n"; map {printf "%5s\t", $_; } @SbWholeIdGroup; print "\nThe end of \@SbWholeIdGroup\n\n\n"; 
    my $sortWholeIdGpHash; for (my $j=0; $j<@SbWholeIdGroup; $j++){$sortWholeIdGpHash->{ $SbWholeIdGroup[$j] }=$j; }

    my @SbSubIdGroup  =@{ $subGroupIds };      @SbSubIdGroup  =sort {$sortWholeIdGpHash->{$a} <=> $sortWholeIdGpHash->{$b}} @SbSubIdGroup;   #print "\n\n Showing the \@SbSubIdGroup:\n\n"; map {printf "%5s\t", $_; } @SbSubIdGroup; print "\nThe end of \@SbSubIdGroup\n\n\n"; 
    my $sortSubIdGpHash; for (my $k=0; $k<@SbSubIdGroup; $k++)    {$sortSubIdGpHash->{ $SbSubIdGroup[$k] }=$k;     }  
    
    #print "\n\n\cl\&print_all_sub_array(\$sortWholeIdGpHash)\n";&print_all_sub_array($sortWholeIdGpHash); print "\n\n";
    #print "\n\n\&print_all_sub_array(\$sortSubIdGpHash)\n";&print_all_sub_array($sortSubIdGpHash); print "\cl\n\n";
    
    #print "\$findPointID=$findPointID\n";
    my $fPID_inwhole_idx= $sortWholeIdGpHash->{$findPointID};     #print "\$fPID_inwhole_idx=\$sortWholeIdGpHash->{$findPointID}=$fPID_inwhole_idx\n";
    my $fPID_inSub_idx=   $sortSubIdGpHash->{$findPointID};       #print "\$fPID_inSub_idx=\$sortSubIdGpHash->{$findPointID}=$fPID_inSub_idx\n";
    
    my @extendPoints_idxAr;
    
    
    
    #���õ���SubGroup�У�����  ����û�е�
    if ( ($fPID_inSub_idx-1) < 0 ){                                             #����inSub����  û�е�
      
      #for (my $i=0; $i<$fPID_inwhole_idx; $i++){push @extendPoints_idxAr,$i;}
      
      #print "a___a (\$fPID_inSub_idx=$fPID_inSub_idx -1) < 0 :\t";printf "%5d < 0\n",($fPID_inSub_idx-1);
      if ( ($fPID_inwhole_idx-1) < 0 ){                                                                                         #����inwhole����  û�е�
        #print "b___b (\$fPID_inwhole_idx=$fPID_inwhole_idx -1) < 0 :\t";printf "%5d < 0\n",($fPID_inwhole_idx-1);
      }
      else{                                                                                                                     #����inwhole����  �е�
      	push @extendPoints_idxAr, (0 .. ($fPID_inwhole_idx-1));  #printf "c___c Else 1, in (0 .. (\$fPID_inwhole_idx-1)), it is: (0 .. %-2d)\n", ($fPID_inwhole_idx-1) ;
      }
    }
    else {                                                                      #����inSub����  �е�
      #my $upStreamPoint_idx=($sortWholeIdGpHash->{$SbSubIdGroup[($fPID_inSub_idx-1)]} )
      if (  ( ($sortWholeIdGpHash->{$SbSubIdGroup[($fPID_inSub_idx-1)]} ) +1 ) > ($fPID_inwhole_idx-1)  ){       #����inwhole ���ε�� Ŀ���֮�� û��Ԫ��
        #printf "d___d (%10d > %-10d)\n\n",( ( ($sortWholeIdGpHash->{$SbSubIdGroup[($fPID_inSub_idx-1)]} ) ) + 1 ), ($fPID_inwhole_idx-1);
      }
      else {                                                                                                                    #����inwhole ���ε�� Ŀ���֮��   ��Ԫ��
        push @extendPoints_idxAr, (  ( ( ($sortWholeIdGpHash->{$SbSubIdGroup[($fPID_inSub_idx-1)]} ) ) + 1 ) .. ($fPID_inwhole_idx-1)  );
        #print "e___e Else 2, (  ( ( [(\$sortWholeIdGpHash->{\$SbSubIdGroup[($fPID_inSub_idx-1)]} ) ) + 1 ) .. (\$fPID_inwhole_idx-1)  )\n";
        #print "(  ( ( (\$sortWholeIdGpHash->{$SbSubIdGroup[($fPID_inSub_idx-1)]} ) ) + 1 ) .. ($fPID_inwhole_idx-1)  )\n";
        #print "(  ( ( ($sortWholeIdGpHash->{$SbSubIdGroup[($fPID_inSub_idx-1)]} ) ) + 1 ) .. ($fPID_inwhole_idx-1)  )\n";
        #print "(  ( ( ($sortWholeIdGpHash->{$SbSubIdGroup[($fPID_inSub_idx-1)]} ) ) + 1 ) .. ($fPID_inwhole_idx-1)  )\n";
        #printf "(%10d .. %-10d)\n\n",( ( $SbWholeIdGroup[($sortWholeIdGpHash->{$SbSubIdGroup[($fPID_inSub_idx-1)]} )] ) + 1 ), ($fPID_inwhole_idx-1);
        #print "\n\n1111     Showing the \@extendPoints_idxAr:\n\n"; map {printf "%5d\t", $_; } @extendPoints_idxAr; print "\n\n";
      }
    }
    #�ټ��õ���SubGroup�У�����  ����û�е�
    if ( ($fPID_inSub_idx+1) > (@SbSubIdGroup-1) ){                             #����inSub����  û�е�
    	#print "f___f (\$fPID_inSub_idx=$fPID_inSub_idx + 1) > (\@SbSubIdGroup-1) :\t";printf "%5d < %-5d\n",($fPID_inSub_idx-1),(@SbSubIdGroup-1) ;
    	if ( ($fPID_inwhole_idx+1) > ($ExtTotalColNb-1) ){                                                                       #���� �Ѿ�   �������group�����ұߵ�λ����
    	  #print "g___g (\$fPID_inSub_idx=$fPID_inSub_idx + 1) > (\$ExtTotalColNb=$ExtTotalColNb-1) :\t";printf "%5d > %-5d\n",($fPID_inSub_idx+1),($ExtTotalColNb-1) ;
    	}
    	else{                                                                                                                     #����     û �������group�����ұߵ�λ����
    		push @extendPoints_idxAr, (  ($fPID_inwhole_idx+1) .. ($ExtTotalColNb-1)  );
    		#print "h___h else 3: (\$fPID_inSub_idx=$fPID_inSub_idx + 1) < (\$ExtTotalColNb=$ExtTotalColNb-1) :\t";printf "%5d < %-5d\n",($fPID_inSub_idx+1),($ExtTotalColNb-1) ;
    		#printf "%5d .. %-5d\n\n",($fPID_inSub_idx+1),($ExtTotalColNb-1) ; 
    		#print "\n\n2222     Showing the \@extendPoints_idxAr:\n\n"; map {printf "%5d\t", $_; } @extendPoints_idxAr; print "\n\n";
    	}
    }
    else {                                                                      #����inSub����  �е�
    	#print "i___i \n";
    } 
    
    #print "\n\$fPID_inwhole_idx=$fPID_inwhole_idx\n121212 Showing the \@extendPoints_idxAr:\n"; map {printf "%5d\t", $_; } @extendPoints_idxAr; print "\n\n";
    
    my $extendPoints_idxAr_addr=[@extendPoints_idxAr];
    my $finalOut=[$fPID_inwhole_idx, $extendPoints_idxAr_addr];
    return $finalOut;
  }

}


##sub8.6          makePThashFromHtmls{               #�����ϲ����ɵ� �������ݽṹ������ html���������Ϣ ��ɵ� sheet�� hash����hash��format��Ϊkey
sub makePThashFromHtmls{               #�����ϲ����ɵ� �������ݽṹ������ html���������Ϣ ��ɵ� sheet�� hash����hash��format��Ϊkey
	
	my ($inpt3Darr,
	    $ClassfiedExcelHash, 
	    $inptIPStypes, 
	    $outPtSheetHash, 
	    #$lastStepFmtHash, 
	    $cWidthHash_toUpgrade, 
	    $multCol_startAtNB_2)=@_;
	
	#my $multCol_startAtNB_2=5;   ������������� html������ sheet����������ɵļ��и�����Ϣ�������ġ�
	my $homManySubCol=2;                                  #�������ָ��һ���м�����������һ��domain��Ϣ

	
	my @classifiedTypesInOrder=@{$inptIPStypes};  #('Family',    'Domain',    'SITE',     'Repeat',     'SUPERFAMILY',       'GENE3D',      'Pfam',      'PANTHER',      'PRINTS',      'PROSITE profiles',      'PROSITE patterns',       'PIRSF',       'UNKNOWN',      'SMART',        'HAMAP',       'COILS',       'Biological Process',      'Molecular Function',       'Celluar Component'    );
	my ($subSortedGrouped3Darray, $bigestCOlnbForEachOrgCol)=@{ (&subGroupSort3Darray( $inpt3Darr )) };     #��html���������Ϣ�����ɵĶ�ά���飬���з��������Ȳ���
  #DirFileHandle::PrintDumper_plReadFile ("subSortedGrouped3Darray.ary", $subSortedGrouped3Darray)  if (    ( defined ( $subSortedGrouped3Darray ) ) && (   (  ref ( $subSortedGrouped3Darray ) eq 'HASH'  ) || (  ref ( $subSortedGrouped3Darray ) eq 'ARRAY'  )   )    );
  #DirFileHandle::PrintDumper_plReadFile ("bigestCOlnbForEachOrgCol.hsh", $bigestCOlnbForEachOrgCol)  if (    ( defined ( $bigestCOlnbForEachOrgCol ) ) && (   (  ref ( $bigestCOlnbForEachOrgCol ) eq 'HASH'  ) || (  ref ( $bigestCOlnbForEachOrgCol ) eq 'ARRAY'  )   )    );
  #����ļ��д��룬��Ϊ��������Ӧ�� �е������ı仯
  my $addNB_Hash;       #�½���һ��hash��
  $addNB_Hash->{0}=0;   #��ʼ����hash��һ�� key �� val
  my $add_nb=0;         #�������
  my @sortEdBig2Small_colArray= sort {$b<=>$a} (keys (%{ $bigestCOlnbForEachOrgCol }) );   #����
  print "\$sortEdBig2Small_colArray[0]=$sortEdBig2Small_colArray[0]\n";
  for (my $i=0; $i<$sortEdBig2Small_colArray[0]; $i++){ 
  	print "\$i=$i\n";
    my $eachAddNb=0; if (   (  defined ( $bigestCOlnbForEachOrgCol->{$i} )  ) && ( $bigestCOlnbForEachOrgCol->{$i}>1 )   ){ $eachAddNb=$bigestCOlnbForEachOrgCol->{$i}-1; }
    $add_nb+=$eachAddNb;
    $addNB_Hash->{($i+1)}=$add_nb; 
  }
  
  
  #$insub_clust3Darray���������װ����һ����ά���ݽṹ�ĵ�ַ�������ݽṹ�Բ�ͬ��domain��family�������ڵ���Ϊ ��һ��key����ÿ�еĲ�ͬ�� IPScan��id ��Ϊ�ڶ���key����ײ���ָ���ֵ���� һ������ĵ�ַ����������� װ���ǣ�����domain��family���ڵ��к�IPScan��id��ͬ���У�����Ϊ���У�����ɵ�����
  my $insub_clust3Darray;  
  # $insub_subTYPEhash ���������װ����һ��hash�ĵ�ַ����hash��key�� ��ͬfamily��domain���ͣ��������sheet�еĵ�һ�е������� ����hash��ֵ����ָ�� �����%typeColNb��"��ͬ��family��domain" �Ȳ�ͬ����ָ������֡�
  my $insub_subTYPEhash; 
  # $insub_subCOLhash ���������װ����һ��hash�ĵ�ַ����hash��key�� ��ͬfamily��domain���ͣ��������sheet�еĵ�һ�е������� ����hash��ֵ����ָ�� �����ɵ�sheet�У���������ռ�Ķ���е���ֵ��ɵ����� �ĵ�ַ��
  my $insub_subCOLhash;

  
  
  
  my $ipsFormatKeyHash;   #����һ��hash�����е�key���� ����IPRid����Щid����ʾ�˲�ͬ�� domain ��family ��������������Ĵ����ʡ� ���hash��װ��key��������format��hash�У�Ҳ��key��
  
  foreach my $clTpKey (@classifiedTypesInOrder){                                                   print "warn\$clTpKey =$clTpKey \n";  #warn "\$clTpKey =$clTpKey \n";
    foreach my $IPRidK (    sort { $a cmp $b } ( keys ( %{ $ClassfiedExcelHash->{$clTpKey} } ) )    ){
      $ipsFormatKeyHash->{$IPRidK}=1; #�� IPRid����װ�����hash��key��
      
      for (my $i=0; $i<@{ $ClassfiedExcelHash->{$clTpKey}->{$IPRidK} }; $i++){
        #my ($orgRow, $gme_tpNb, $specName, $htmlID, $IPRid,  $IPRdf, $IPRweb, $newRow, $colNb)=@{ $ClassfiedExcelHash->{$clTpKey}->{$IPRidK}->[$i] };  #print "(\$orgRow=$orgRow, \$gme_tpNb=$gme_tpNb, \$specName=$specName, \$htmlID=$htmlID, \$IPRid=$IPRid,  \$IPRdf=$IPRdf, \$IPRweb=$IPRweb, \$newRow=$newRow, \$colNb=$colNb)=\@{ \$ClassfiedExcelHash->{\$clTpKey=$clTpKey}->{\$IPRidK=$IPRidK}->[\$i=$i] }\n";
        my (                                $htmlID, $IPRid,  $IPRdf, $IPRweb, $newRow, $colNb)=@{ $ClassfiedExcelHash->{$clTpKey}->{$IPRidK}->[$i] };   print "(\$htmlID=$htmlID, \$IPRid=$IPRid,  \$IPRdf=$IPRdf, \$IPRweb=$IPRweb, \$newRow=$newRow, \$colNb=$colNb)=\@{ \$ClassfiedExcelHash->{\$clTpKey=$clTpKey}->{\$IPRidK=$IPRidK}->[\$i=$i] }\n";
                                                                                                                                                                    #print "\n\n\cl\&print_all_sub_array(\$the3Darray->{$newRow}->{$colNb})\n";&print_all_sub_array($the3Darray->{$newRow}->{$colNb}); print "\cl\n\n";      #map {print "\$_=$_";} ( keys ( %{ $the3Darray->{$newRow}->{$colNb} } ) );
        push @{ $outPtSheetHash->{ 'PepID' } }, [$newRow, 0, $htmlID ];
        
        my $addNbHere=$addNB_Hash->{$colNb};                        print "\n\n\$clTpKey=$clTpKey\t\$IPRidK=$IPRidK\t\$addNbHere=\$addNB_Hash->{\$colNb}=$addNbHere=\$addNB_Hash->{$colNb}\n"; #��������ǰ������
        
        my $innerColNb=$subSortedGrouped3Darray->{$newRow}->{$colNb}->{$IPRidK}->[1];
        my @otherSameColNb=@{ $subSortedGrouped3Darray->{$newRow}->{$colNb}->{$IPRidK}->[2] };
        my @allInerColNbAr=@otherSameColNb; push  @allInerColNbAr, $innerColNb;
        
        foreach my $inerColNb (@allInerColNbAr){
        	
        	my $secLvCol=$addNbHere+$colNb+$inerColNb;  #my $outV=ExcelHandle::makeHyperLk($InAdr, $InVal);
        	my $hpLKwithIPRdf =ExcelHandle::makeHyperLk( $IPRweb, $IPRdf  );  #���� hyperlink�� ָ����ַ�� ��ʾ IPRdf
        	my $hpLKwithIPRidK=ExcelHandle::makeHyperLk( $IPRweb, $IPRidK );  #���� hyperlink�� ָ����ַ�� ��ʾ IPRidK
          
          my $NbAfterStart=$secLvCol-$multCol_startAtNB_2;
          
          my $muti_col_0=$multCol_startAtNB_2+$NbAfterStart*$homManySubCol-1+1;
          my $muti_col_1=$multCol_startAtNB_2+$NbAfterStart*$homManySubCol-1+2; 
          #my $muti_col_2=$multCol_startAtNB_2+$NbAfterStart*$homManySubCol-1+3;
          #my $muti_col_3=$multCol_startAtNB_2+$NbAfterStart*$homManySubCol-1+4;
          #my $muti_col_4=$multCol_startAtNB_2+$NbAfterStart*$homManySubCol-1+5;
  
          #��������������ʾ �����ID��
          push @{ $outPtSheetHash->{ $IPRidK } }, [$newRow, $muti_col_0, $hpLKwithIPRidK ];  printf "%-120s\t", "\$outPtSheetHash->{ $IPRidK }\t\$newRow, \$muti_col_0, \$hpLKwithIPRidK,  \$IPRidK, \$IPRweb:";  printf "%20s\t%20s\t%20s\t%20s\t%20s\n\n",$newRow, $muti_col_0, $hpLKwithIPRidK, $IPRidK, $IPRweb;
          #������������ʾ ע�ͻ���˵���ܽṹ����
          push @{ $outPtSheetHash->{ $IPRidK } }, [$newRow, $muti_col_1, $hpLKwithIPRdf];    printf "%-120s\t", "\$outPtSheetHash->{ $IPRidK }\t\$newRow, \$muti_col_1, \$hpLKwithIPRdf,   \$IPRdf,  \$IPRweb:";  printf "%20s\t%20s\t%20s\t%20s\t%20s\n\n",$newRow, $muti_col_1, $hpLKwithIPRdf,  $IPRdf, $IPRweb;
          
          
          #push @{ $specpt2sheet_hash_3->{ $FDS_Format_prop_string } }, [$newRow, $muti_col_0, 'Clustalw',   $IPRweb];  print "00 \$specpt2sheet_hash_3->{ $FDS_Format_prop_string }\t";  printf "%20s\t%20s\t%20s\t%20s\n\n",$newRow, $muti_col_0, $IPRidK, $IPRweb;
          #push @{ $specpt2sheet_hash_3->{ $FDS_Format_prop_string } }, [$newRow, $muti_col_1, 'Tree',       $IPRweb];  print "11 \$specpt2sheet_hash_3->{ $FDS_Format_prop_string }\t";  printf "%20s\t%20s\t%20s\t%20s\n\n",$newRow, $muti_col_1, $IPRidK, $IPRweb;
          #push @{ $specpt2sheet_hash_3->{ $FDS_Format_prop_string } }, [$newRow, $muti_col_2, 'BootedTree', $IPRweb];  print "22 \$specpt2sheet_hash_3->{ $FDS_Format_prop_string }\t";  printf "%20s\t%20s\t%20s\t%20s\n\n",$newRow, $muti_col_2, $IPRidK, $IPRweb;
          #push @{ $specpt2sheet_hash_3->{ $FDS_Format_prop_string } }, [$newRow, $muti_col_3, $IPRidK,      $IPRweb];  print "33 \$specpt2sheet_hash_3->{ $FDS_Format_prop_string }\t";  printf "%20s\t%20s\t%20s\t%20s\n\n",$newRow, $muti_col_3, $IPRidK, $IPRweb;
          #push @{ $specpt2sheet_hash_3->{ $FDS_Format_prop_string } }, [$newRow, $muti_col_4, $IPRdf,       $IPRweb];  print "44 \$specpt2sheet_hash_3->{ $FDS_Format_prop_string }\t";  printf "%20s\t%20s\t%20s\t%20s\n\n",$newRow, $muti_col_4, $IPRdf,  $IPRweb;
                  
          $cWidthHash_toUpgrade->{$muti_col_0}=1;   print "\$cWidthHash_toUpgrade->{\$muti_col_0}=\$cWidthHash_toUpgrade->{$muti_col_0}=$cWidthHash_toUpgrade->{$muti_col_0}\n";    #�п��������� ����$muti_col_0һ��
          $cWidthHash_toUpgrade->{$muti_col_1}=4;   print "\$cWidthHash_toUpgrade->{\$muti_col_1}=\$cWidthHash_toUpgrade->{$muti_col_1}=$cWidthHash_toUpgrade->{$muti_col_1}\n";    #�п��������� ����$muti_col_1һ��
          
          
          
          push @{ $insub_clust3Darray->{$secLvCol}->{$IPRidK} }, $newRow;
          $insub_subTYPEhash->{$secLvCol} =  $colNb;
          $insub_subCOLhash->{$secLvCol}=[$muti_col_0, $muti_col_1]; 
          #$insub_subCOLhash->{$secLvCol}=[$muti_col_0, $muti_col_1, $muti_col_2, $muti_col_3];
          
        }                                  
      }                                   
                                      
     
    }                                   
  }                      
   
  my $lastStepFmtHash; $lastStepFmtHash=Interproscan_MultiplAnalysis::MakeFormatHash_IPSid_Version($lastStepFmtHash, $ipsFormatKeyHash);
  
  return 
  [ 
    $outPtSheetHash, 
    $lastStepFmtHash, 
    $cWidthHash_toUpgrade, 
    $insub_clust3Darray, 
    $insub_subTYPEhash, 
    $insub_subCOLhash
  ];
  
}




#my $tempFmtHash=Interproscan_MultiplAnalysis::MakeFormatHash_IPSid_Version($tempFmtHash, $interproscanIDhash);
##sub8.7          MakeFormatHash_3 {   #������� format_Hash����װ�ض��� ��ʽ�ĺ����� ���е���Ҫ�ֶ�����ĸ�ʽ��Ϣ����������װ��
sub MakeFormatHash_IPSid_Version {   #������� format_Hash����װ�ض��� ��ʽ�ĺ����� ���е���Ҫ�ֶ�����ĸ�ʽ��Ϣ����������װ��
  my ($tempFmtHash, $interproscanIDhash)=@_;
  
  #interproscan ��id Ϊkey������hash
  my $BgrColNb=11;
  foreach my $IPSid (keys ( %{$interproscanIDhash})){
    $tempFmtHash->{$IPSid}=     { -bold => 1, -bg_color =>  $BgrColNb, -size => 8  };
    $BgrColNb++;
    if ($BgrColNb > 63) {$BgrColNb=11;}      
  }
   
  #�������ֶ���д�ĸ�ʽ
  $tempFmtHash->{'normal'}= {                                          };
  $tempFmtHash->{'PepID' }= {-bold => 1                                };
  
  return $tempFmtHash;
}

##sub8.9          setColumnFormat_2 {    #�����sheet�����п�����
sub setColumnFormat_2 {    #����������п����� 
  my ($InXlsSheet,$inColWidth_hash)=@_;                #����$inColWidth_hash�������п�
  
  foreach my $coKey (   sort {$a <=> $b} (  keys %{ $inColWidth_hash }  )   ){  #warn "\$coKey=$coKey\n";
 
     $InXlsSheet->set_column ( $coKey,  $coKey,  $inColWidth_hash->{$coKey} );  print "\$inColWidth_hash->{\$coKey}=\$inColWidth_hash->{$coKey}=$inColWidth_hash->{$coKey}\n";
    
  }   
  
}



##sub8.10          sortClustalByIPScanResult{      #�������ɵ�html������sheet��������clustalw����������clustawl����������µ�sheet��  
sub sortClustalByIPScanResult{    #�������ɵ�html������sheet��������clustalw����������clustawl����������µ�sheet��
	my ($inpt_htmlSheetHash, 
	    $clust_3d_arry, 
	    $sub_type_hash, 
	    $subSub_col_hash,
	    $inFastaHASH,
	    $inSub_clustalRSTdir,
	    $inSub_hperLkCLUSTEDrstDIR, 
	    $insub_clustalw_or_not, 
	    $typeColNbHASHaddress, 
	    $multCol_startAtNB_3)=@_;   #($newSheetHash_KeyedByFMT,  $out_clust3Darray, $out_subTYPEhash, $out_subCOLhash, , $middle_clustalRSTdir,$middle_hperLkCLUSTEDrstDIR, $middle_clustalw_or_not, $ipsTypeNB, $middle_multCol_startAtNB);
	#$clust_3d_arry���������װ����һ����ά���ݽṹ�ĵ�ַ�������ݽṹ�Բ�ͬ��domain��family�������ڵ���Ϊ ��һ��key����ÿ�еĲ�ͬ�� IPScan��id ��Ϊ�ڶ���key����ײ���ָ���ֵ���� һ������ĵ�ַ����������� װ���ǣ�����domain��family���ڵ��к�IPScan��id��ͬ���У�����Ϊ���У�����ɵ�����
  #$sub_type_hash ���������װ����һ��hash�ĵ�ַ����hash��key�� ��ͬfamily��domain���ͣ��������sheet�еĵ�һ�е������� ����hash��ֵ����ָ�� �����%typeColNb��"��ͬ��family��domain" �Ȳ�ͬ����ָ������֡�
  #$subSub_col_hash ���������װ����һ��hash�ĵ�ַ����hash��key�� ��ͬfamily��domain���ͣ��������sheet�еĵ�һ�е������� ����hash��ֵ����ָ�� �����ɵ�sheet�У���������ռ�Ķ���е���ֵ��ɵ����� �ĵ�ַ��

  my ($warnMsgBody, $warnMsgHead, $die_MsgHead, $caller_inform)=@{ DieWork::BuildWarnDieHeadInform( 'Interproscan_MultiplAnalysis', 'sortClustalByIPScanResult' ) };

  #�������䣬�Ǹ�������ͬ��domain���ͣ����������Ź�������������� ��Ҫ�� %typeColNb���hash��  ����� $type_order���hash�� ��ַ
  #my %typeColNb=        
  #('Family'=>0, 'Domain'=>1, 'SITE'=> 2, 'Repeat'=> 3, 'SUPERFAMILY' => 4,  'GENE3D' => 5, 'Pfam' => 6, 'PANTHER' => 7, 'PRINTS' => 8, 'PROSITE profiles' => 9, 'PROSITE patterns' => 10, 'PIRSF' => 11, 'UNKNOWN' => 12,'SMART' => 13,  'HAMAP' => 14, 'COILS' => 15, 'Biological Process' => 16,'Molecular Function' => 17, 'Celluar Component'=> 18);
  #('Family'=>2, 'Domain'=>3, 'SITE'=> 4, 'Repeat'=> 5, 'SUPERFAMILY' => 6,  'GENE3D' => 7, 'Pfam' => 8, 'PANTHER' => 9, 'PRINTS' => 10, 'PROSITE profiles' => 11, 'PROSITE patterns' => 12, 'PIRSF' => 13, 'UNKNOWN' => 14,'SMART' => 15,  'HAMAP' => 16, 'COILS' => 17, 'Biological Process' => 18,'Molecular Function' => 19, 'Celluar Component'=> 20);
  my %typeColNb=%{ $typeColNbHASHaddress };
  foreach my $fmdomKey (keys %typeColNb){  $typeColNb{$fmdomKey}+=$multCol_startAtNB_3;}  #����װ�и�����ͬ�������͵���ŵ�hash��key�����飬��ÿ��key����Ԥ���������������������������$multCol_startAtNB_1=5;
  my $type_order;  my @tempArray;                             
  @tempArray=( $typeColNb{'Family'}                                   );                   map{$type_order->{$_}=0;} @tempArray;
  @tempArray=( $typeColNb{'Domain'}                                   );                   map{$type_order->{$_}=1;} @tempArray;
  @tempArray=( $typeColNb{'SITE'}                                     );                   map{$type_order->{$_}=2;} @tempArray;#@tempArray=(6,7,8,9,10,11,12,13,14,15,16);map{$type_order->{$_}=3;} @tempArray;
  @tempArray=( $typeColNb{'SUPERFAMILY'}                              );                   map{$type_order->{$_}=3;} @tempArray;
  @tempArray=( $typeColNb{'GENE3D'}                                   );                   map{$type_order->{$_}=4;} @tempArray;
  @tempArray=( $typeColNb{'Pfam'}                                     );                   map{$type_order->{$_}=5;} @tempArray;
  @tempArray=( $typeColNb{'PANTHER'}                                  );                   map{$type_order->{$_}=6;} @tempArray;
  @tempArray=( $typeColNb{'PRINTS'}                                   );                   map{$type_order->{$_}=7;} @tempArray;
  @tempArray=( $typeColNb{'PROSITE profiles'}                         );                   map{$type_order->{$_}=8;} @tempArray;
  @tempArray=( $typeColNb{'PROSITE patterns'}                         );                   map{$type_order->{$_}=9;} @tempArray;
  @tempArray=( $typeColNb{'PIRSF'}                                    );                   map{$type_order->{$_}=10;}@tempArray;
  @tempArray=( $typeColNb{'UNKNOWN'}                                  );                   map{$type_order->{$_}=11;}@tempArray;
  @tempArray=( $typeColNb{'SMART'}                                    );                   map{$type_order->{$_}=12;}@tempArray;
  @tempArray=( $typeColNb{'HAMAP'}                                    );                   map{$type_order->{$_}=13;}@tempArray;#@tempArray=(18,19,20  );                  map{$type_order->{$_}=4;}@tempArray;
  @tempArray=( $typeColNb{'Biological Process'}                       );                   map{$type_order->{$_}=14;}@tempArray;
  @tempArray=( $typeColNb{'Molecular Function'}                       );                   map{$type_order->{$_}=15;}@tempArray;
  @tempArray=( $typeColNb{'Celluar Component'}                        );                   map{$type_order->{$_}=16;}@tempArray;
  @tempArray=( $typeColNb{'Repeat'},            $typeColNb{'COILS'}   );                   map{$type_order->{$_}=17;}@tempArray;
 
   
  #������Щ��������Ҫ���ڽ�cluastalw�Ľ����sheet�е���Ӧ��Ϣ֮���������
  my $clu_pep_hash;                                                              #���װÿ������ͬ�����ڵĵ��������ļ������������бȶ�  
  my $msf_file_hash;                                                             #����� �����бȶԽ�� �� ��Ӧ�� ������Ӧ����
  my $nwk_file_hash;                                                             #����� �����бȶԵĽ����� ��� ��Ӧ����
  my $cus_file_hash;
  my $fas_file_hash;
  my $link_RowCol_2_clu;                                                         #���������ÿ���㣬�����Ӧ��clustalw��������е������ַ�����Ӧ����
  #my $link_RowCol_2_msf;                                                        #���������ÿ���㣬�����Ӧ��clustalw��msf�����������
  my $big_link_RowCol_2_url;                                                     #���������ÿ����Ͷ����б��еĸ��������������� 

  #�����µ�hash�ĵ�ַ����hash��key��$type_order������潨����hash�ĸ�����ͬ���Ͷ�Ӧ������ֵ����hash��val���� �ɡ���ͬ����IPSid��Ӧ������ɵ�����Ĵ�С����ͬ���͵���������IPS��id���Լ�$clust_3d_arry���ݸ�����key�����Ķ��������ɵ�����ĵ�ַ������ɵ�����ĵ�ַ��
  my $longestColArray; 
  foreach my $colNbKeyHr ( sort {$a<=>$b}( keys %{ $clust_3d_arry } ) ){
  	foreach my $IPRidKKeyHr ( sort {$a cmp $b}( keys %{ $clust_3d_arry->{$colNbKeyHr} } ) ){    #printf  "%-20s\t\t", "Wo kao"; map {printf "%10s\t", $_;} @{ $clust_3d_arry->{$colNbKeyHr}->{$IPRidKKeyHr} }; print  "\n"; 
      my @sortedRow=@{ $clust_3d_arry->{$colNbKeyHr}->{$IPRidKKeyHr} };     #print "111 Before link_RowCol_2_clu define!\t";map {printf "%10s\t", $_;} @sortedRow; print  "\n";  
      @sortedRow= sort {$a<=>$b}  @sortedRow;                               #print "222 Before link_RowCol_2_clu define!\t";map {printf "%10s\t", $_;} @sortedRow; print  "\n"; 
      my $sotRowIn2List=join ',-123+,',@sortedRow;                             #print "333 Before link_RowCol_2_clu define!\t";map {printf "%10s\t", $_;} @sortedRow; print  "\n"; 
      $clu_pep_hash->{$sotRowIn2List}=[@sortedRow];                          #���װÿ������ͬ�����ڵĵ��������ļ������������бȶ� 
      #print "444 Before link_RowCol_2_clu define!\t";map {printf "%10s\t", $_;} @sortedRow; print  "\n"; 
      map {$link_RowCol_2_clu->{$_}->{$colNbKeyHr}=$sotRowIn2List;   
      	    #printf "%10s\t%10s\t%30s\t","CkMark",$colNbKeyHr,$IPRidKKeyHr; printf "%-50s\t%-30s\t","\$link_RowCol_2_clu->{$_}->{$colNbKeyHr}=\$sotRowIn2List=", $link_RowCol_2_clu->{$_}->{$colNbKeyHr};
      	  } @sortedRow; #print  "\n";  #���������ÿ���㣬�����Ӧ��clustalw�Ľ����������
      
      my $arSize=@{ $clust_3d_arry->{$colNbKeyHr}->{$IPRidKKeyHr} };  #printf "%10s\t%10s\t%30s\t","CkMark",$colNbKeyHr,$IPRidKKeyHr; printf  "%-20s\t\t", "\$arSize=$arSize"; map {printf "%10s\t", $_;} @{ $clust_3d_arry->{$colNbKeyHr}->{$IPRidKKeyHr} }; print  "\n"; 
     
      #                            �������        ��һ����            ʵ�ʵ���
      push @{ $longestColArray->{ $type_order->{ $sub_type_hash->{ $colNbKeyHr } } } }, [$arSize,$colNbKeyHr,$IPRidKKeyHr,$clust_3d_arry->{$colNbKeyHr}->{$IPRidKKeyHr}];  
      printf "%-95s\t%50s\t%50s\t%35s\t","push sth into \@{ \$longestColArray->{ \$type_order->{ \$sub_type_hash->{ \$colNbKeyHr } } } }=","\@{ \$longestColArray->{ \$type_order->{ \$sub_type_hash->{ $colNbKeyHr } } } }=","\@{ \$longestColArray->{ \$type_order->{ $sub_type_hash->{ $colNbKeyHr } } } }=","\@{ \$longestColArray->{ $type_order->{ $sub_type_hash->{ $colNbKeyHr } } } }:=";
      printf "\$arSize=%5s\t\$colNbKeyHr=%4s\t\$IPRidKKeyHr=%-30s\t\$clust_3d_arry->{\$colNbKeyHr}->{\$IPRidKKeyHr}=\t",$arSize,$colNbKeyHr,$IPRidKKeyHr,$clust_3d_arry->{$colNbKeyHr}->{$IPRidKKeyHr}; map {printf "%5d ", $_;} @{ $clust_3d_arry->{$colNbKeyHr}->{$IPRidKKeyHr} }; print "\n";
    }
  }
    
 
  ################################################################################################################################################################################################################################################################################################################################################################################################################################################
  #######  ���漸����������clustaw��� ###########################################################################################################################################################################################################################################################################################################################################################################################################
  
  my $rowCOLhtmlHASH=ExcelHandle::FmtHashChgIn2RowColFmtHash($inpt_htmlSheetHash); #��html������g  m e��һ�����к���µ�sheet��hash��ת��Ϊ row col ��Ϊkey��hash��
  #DirFileHandle::PrintAndWarnDumper ($rowCOLhtmlHASH, "\n 201911141608-0-0-0 \$rowCOLhtmlHASH=$rowCOLhtmlHASH");
  if (  -d ($inSub_clustalRSTdir)  ) {}
  else {  system ("mkdir -p $inSub_clustalRSTdir");  }
  my $cluPepNameNb=1;
  foreach my $clu_key (sort ( keys ( %{ $clu_pep_hash } ) ) ){
  	
  	my $CluAbsPathFile="$inSub_clustalRSTdir/${cluPepNameNb}.txt";  print "\$CluAbsPathFile=$CluAbsPathFile\n";
  	open (CLUINPEP,">$CluAbsPathFile") or die "cannot create \$CluAbsPathFile=$CluAbsPathFile $!\n";
  	my $howNamySeqs=0;
    foreach my $row_4_clu ( @{ $clu_pep_hash->{$clu_key} } ){   print "\$row_4_clu=$row_4_clu\n";
    	#������ ���е�fasta��ʽ���ַ�����Ϣ
    	DieWork::Print_and_warn( "\n 201911141608-0-0-1 $rowCOLhtmlHASH, $row_4_clu, 0, 'V' \n\n");
    	my $fastaHead=ExcelHandle::getValHplinkFmtFrom_RowColHash( $rowCOLhtmlHASH, $row_4_clu, 0, 'V'); #�� row col���� html����sheet��hash�У�ȡ����1��0��ͷ������1���е�value����ʵ����pep��ID
    	
    	DieWork::Check_DfdNoEmptString_or_DIE( $inFastaHASH->{$fastaHead}->{'0_0_6_____sequence'}, "\$inFastaHASH->{\$fastaHead}->{'0_0_6_____sequence'}=$inFastaHASH->{$fastaHead}->{'0_0_6_____sequence'}", $die_MsgHead, $caller_inform  );
    	my $fastaSeq =$inFastaHASH->{$fastaHead}->{'0_0_6_____sequence'};
    	#ExcelHandle::getValHplinkFmtFrom_RowColHash( $rowCOLhtmlHASH, $row_4_clu, 2, 'L'); #�� row col���� html����sheet��hash�У�ȡ����3��0��ͷ������2���е�label����ʵ����pep������
    	
      my $fastaSeqEachRow=">$fastaHead\n$fastaSeq\n";     print "\$fastaSeqEachRow=\n>\$fastaHead\n\$fastaSeq\n=\n$fastaSeqEachRow\n\n";
      print CLUINPEP $fastaSeqEachRow;
      $howNamySeqs++;
    }
    close (CLUINPEP);
    #my $maoFileNoBoot="noBoot.mao";                �⼸������������ megaʹ�õ��������Ϣ��
    #my $maoFIleWtBoot="withBoot.mao";
    #my $finalMaoFile=$maoFIleWtBoot; if ($howNamySeqs<=3) {$finalMaoFile=$maoFileNoBoot;}
    #my $inFastaFile="${cluPepNameNb}.fas"; 
    #printf MEGACOMAND "M6CC.exe -a %-20s -d %-20s -o %-20s\n\r", $finalMaoFile, $inFastaFile, $cluPepNameNb; 
      
    my $CluPathFileHpLK="$inSub_hperLkCLUSTEDrstDIR/${cluPepNameNb}.txt"; $CluPathFileHpLK=ExcelHandle::Lnx2Dos($CluPathFileHpLK); $CluPathFileHpLK=".\\$CluPathFileHpLK"; #HyperLinkд�� ������URL��ʽ #$CluPathFileHpLK="external:.\\$CluPathFileHpLK";
    print "\$CluPathFileHpLK=$CluPathFileHpLK\n";
    
    
    my $curDir=File::Spec->curdir();
    my $rltPath=File::Spec->abs2rel ($CluAbsPathFile, $curDir); print "ABS_Path=$CluAbsPathFile\nRELATIVE_Path=$rltPath\n";
 
    if ($insub_clustalw_or_not == 1) {
    	FREECPUCHECKPOINT: 
    	my $freeCpuNb=&GetFreeCpuNb;
    	if ($freeCpuNb>=4){
    	
    	my $clustalw_CMD="clustalw $rltPath -QUIET -OUTPUT=gcg  &";  DieWork::Print_and_warn( $clustalw_CMD );
      system ( "$clustalw_CMD" );      
      #warn   "clustalw $rltPath -QUIET -OUTPUT=fasta\n";
      #system ("clustalw $rltPath -QUIET -OUTPUT=fasta &");
      }
      else {goto FREECPUCHECKPOINT;}
    }
    
    open (FASTAFILE,"$inSub_clustalRSTdir/${cluPepNameNb}.txt") or warn "\$inSub_clustalRSTdir/${cluPepNameNb}.txt=$inSub_clustalRSTdir/${cluPepNameNb}.txt cannot be opened :$! \n\n";
    open (NEWFASTAFILE,">$inSub_clustalRSTdir/${cluPepNameNb}.fas") or die "\$inSub_clustalRSTdir/${cluPepNameNb}.fas=$inSub_clustalRSTdir/${cluPepNameNb}.fas cannot be created :$! \n\n";
    
    my $tpMk=$/;
    $/=">";
    while (<FASTAFILE>){
    	s/>$//; s/^>//;
    	if (m/\S+/){
    	  my $tpHr=">$_"; $tpHr=~s/\s+$//;
    	  if ($tpHr=~m/>(\S+.*)\n((:?.*\S+.*)(:?\n.*\S+.*)*)/){
    	    my $seqN=$1;
    	    my $seqSeq=$2;
    	    $seqSeq=~tr/Uu/Cc/;
    	    print NEWFASTAFILE ">$seqN\n$seqSeq\n";
    	  }
    	  else {
    	    die "\$tpHr=$tpHr\n\n";
    	  }
    	}
    }
    $/=$tpMk;
    close (FASTAFILE);
    close (NEWFASTAFILE);
    
    my $MsfPathFileHpLk="$inSub_hperLkCLUSTEDrstDIR/${cluPepNameNb}.msf";             $MsfPathFileHpLk=ExcelHandle::Lnx2Dos($MsfPathFileHpLk); $MsfPathFileHpLk=".\\$MsfPathFileHpLk"; #HyperLinkд�� ������URL��ʽ #$MsfPathFileHpLk="external:.\\$MsfPathFileHpLk";
    my $nwkPathFileHpLk="$inSub_hperLkCLUSTEDrstDIR/${cluPepNameNb}.nwk";             $nwkPathFileHpLk=ExcelHandle::Lnx2Dos($nwkPathFileHpLk); $nwkPathFileHpLk=".\\$nwkPathFileHpLk"; #HyperLinkд�� ������URL��ʽ #$nwkPathFileHpLk="external:.\\$nwkPathFileHpLk";
    my $cusPathFileHpLk="$inSub_hperLkCLUSTEDrstDIR/${cluPepNameNb}_consensus.nwk";   $cusPathFileHpLk=ExcelHandle::Lnx2Dos($cusPathFileHpLk); $cusPathFileHpLk=".\\$cusPathFileHpLk"; #HyperLinkд�� ������URL��ʽ #$cusPathFileHpLk="external:.\\$cusPathFileHpLk";
    my $FasPathFileHpLk="$inSub_hperLkCLUSTEDrstDIR/${cluPepNameNb}.txt";             $FasPathFileHpLk=ExcelHandle::Lnx2Dos($FasPathFileHpLk); $FasPathFileHpLk=".\\$FasPathFileHpLk"; #HyperLinkд�� ������URL��ʽ #$FasPathFileHpLk="external:.\\$FasPathFileHpLk";
    
    #$MsfPathFileHpLk=&Lnx2Dos ($MsfPathFileHpLk);
    $msf_file_hash->{$clu_key}=$MsfPathFileHpLk;  printf "%-70s\t%-20s\t%-200s\n", "\$msf_file_hash->{\$clu_key}=\$msf_file_hash->{$clu_key}=", "\$MsfPathFileHpLk=", $msf_file_hash->{$clu_key};
    $nwk_file_hash->{$clu_key}=$nwkPathFileHpLk;  printf "%-70s\t%-20s\t%-200s\n", "\$nwk_file_hash->{\$clu_key}=\$nwk_file_hash->{$clu_key}=", "\$nwkPathFileHpLk=", $nwk_file_hash->{$clu_key};
    $cus_file_hash->{$clu_key}=$cusPathFileHpLk;  printf "%-70s\t%-20s\t%-200s\n", "\$cus_file_hash->{\$clu_key}=\$clu_file_hash->{$clu_key}=", "\$cusPathFileHpLk=", $cus_file_hash->{$clu_key};
    $fas_file_hash->{$clu_key}=$FasPathFileHpLk;  printf "%-70s\t%-20s\t%-200s\n", "\$fas_file_hash->{\$clu_key}=\$fas_file_hash->{$clu_key}=", "\$FasPathFileHpLk=", $fas_file_hash->{$clu_key};
    $cluPepNameNb++;
  }
  
  foreach my $rKhr ( sort {$a<=>$b} ( keys ( %{ $link_RowCol_2_clu } ) ) ){
    foreach my $cKhr ( sort {$a<=>$b} ( keys ( %{ $link_RowCol_2_clu->{$rKhr} } ) ) ){
      #$link_RowCol_2_msf->{$rKhr}->{$cKhr}=$msf_file_hash->{ $link_RowCol_2_clu->{$rKhr}->{$cKhr} };   printf "%-50s\t%-50s\t%-50s\t%-80s\t%-100s\n",$link_RowCol_2_msf->{$rKhr}->{$cKhr},"=\$link_RowCol_2_msf->{$rKhr}->{$cKhr}=", $msf_file_hash->{ $link_RowCol_2_clu->{$rKhr}->{$cKhr} }, "\$msf_file_hash->{ \$link_RowCol_2_clu->{$rKhr}->{$cKhr} }=","\$msf_file_hash->{ $link_RowCol_2_clu->{$rKhr}->{$cKhr} }";
      
      my $msfCol =$subSub_col_hash->{$cKhr}->[0];
      #my $fasCol =$subSub_col_hash->{$cKhr}->[1];
      #my $nwkCol =$subSub_col_hash->{$cKhr}->[1];
      #my $cusCol =$subSub_col_hash->{$cKhr}->[2];
      #my $fasCol =$subSub_col_hash->{$cKhr}->[3];
      
      #if ($getUrl_or_not==1){
      $big_link_RowCol_2_url->{$rKhr}->{$msfCol}=$msf_file_hash->{ $link_RowCol_2_clu->{$rKhr}->{$cKhr} };    printf "%-50s\t%-50s\t%-50s\t%-80s\t%-100s\n",$big_link_RowCol_2_url->{$rKhr}->{$msfCol},"=\$big_link_RowCol_2_url->{$rKhr}->{$msfCol}=", $msf_file_hash->{ $link_RowCol_2_clu->{$rKhr}->{$cKhr} }, "\$msf_file_hash->{ \$link_RowCol_2_clu->{$rKhr}->{$cKhr} }=","\$msf_file_hash->{ $link_RowCol_2_clu->{$rKhr}->{$cKhr} }";
      #}
      #$big_link_RowCol_2_url->{$rKhr}->{$nwkCol}=$nwk_file_hash->{ $link_RowCol_2_clu->{$rKhr}->{$cKhr} };   printf "%-50s\t%-50s\t%-50s\t%-80s\t%-100s\n",$big_link_RowCol_2_url->{$rKhr}->{$nwkCol},"=\$big_link_RowCol_2_url->{$rKhr}->{$nwkCol}=", $nwk_file_hash->{ $link_RowCol_2_clu->{$rKhr}->{$cKhr} }, "\$nwk_file_hash->{ \$link_RowCol_2_clu->{$rKhr}->{$cKhr} }=","\$nwk_file_hash->{ $link_RowCol_2_clu->{$rKhr}->{$cKhr} }";
      #$big_link_RowCol_2_url->{$rKhr}->{$cusCol}=$cus_file_hash->{ $link_RowCol_2_clu->{$rKhr}->{$cKhr} };   printf "%-50s\t%-50s\t%-50s\t%-80s\t%-100s\n",$big_link_RowCol_2_url->{$rKhr}->{$cusCol},"=\$big_link_RowCol_2_url->{$rKhr}->{$cusCol}=", $nwk_file_hash->{ $link_RowCol_2_clu->{$rKhr}->{$cKhr} }, "\$cus_file_hash->{ \$link_RowCol_2_clu->{$rKhr}->{$cKhr} }=","\$cus_file_hash->{ $link_RowCol_2_clu->{$rKhr}->{$cKhr} }";    
      #$big_link_RowCol_2_url->{$rKhr}->{$fasCol}=$fas_file_hash->{ $link_RowCol_2_clu->{$rKhr}->{$cKhr} };   printf "%-50s\t%-50s\t%-50s\t%-80s\t%-100s\n",$big_link_RowCol_2_url->{$rKhr}->{$fasCol},"=\$big_link_RowCol_2_url->{$rKhr}->{$fasCol}=", $nwk_file_hash->{ $link_RowCol_2_clu->{$rKhr}->{$cKhr} }, "\$fas_file_hash->{ \$link_RowCol_2_clu->{$rKhr}->{$cKhr} }=","\$fas_file_hash->{ $link_RowCol_2_clu->{$rKhr}->{$cKhr} }";    
    }
  }
  #######  ���漸����������clustaw��� ###########################################################################################################################################################################################################################################################################################################################################################################################################################################################
  ################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
 
  #��������������ÿһ�д�֣���ֹ����ǣ� 
  
  my $TotalRowNb=(keys ( $rowCOLhtmlHASH) );  print "\$TotalRowNb=$TotalRowNb\n\n";  #���html�������sheet�������� 
  my $bigScore_forEach_row;  # ���潨��һ��hash��װ��ÿһ�еķ���
  foreach my $eachOrededGroup_Order ( sort { $a<=>$b } ( keys ( %{ $longestColArray } ) ) ){
  	my @iner_sizeArray=@{ $longestColArray->{ $eachOrededGroup_Order } };
  	  	
  	my $SizeScore=1;
    printf "Before Sort: \$eachOrededGroup_Order=%3d\t", $eachOrededGroup_Order; my $stNb=0; map {my $xy=$_; printf "%5s:%-5d ", "($stNb)", $xy->[0]; $stNb++;} @iner_sizeArray;  print "\n";
    @iner_sizeArray=sort { $a->[0] <=> $b->[0] } @iner_sizeArray;   #���������⣡����������������������������
    printf "Sorted:    : \$eachOrededGroup_Order=%3d\t", $eachOrededGroup_Order;    $stNb=0; map {my $xy=$_; printf "%5s:%-5d ", "($stNb)", $xy->[0]; $stNb++;} @iner_sizeArray;  print "\n";
   
    my $lastSizeNb=0;
    for (my $i=0; $i<@iner_sizeArray; $i++){    ########�����ÿһ��Ԫ�أ�ָ���ǰ�size��С����󣬵�һ���� ����� ��ַ��������� size�������� ����id�� ����ĵ�ַ
        	
    	my @tempSgArray=@{ $iner_sizeArray[$i] }; 
    	
    	my ($sizeNb, $colHhH, $idHrH, $rowAddHrH)=@tempSgArray;  print "\n\n\n\cl\cl\n"; printf "\$sizeNb=%-20s\t\$colHhH=%-20s\t\$idHrH=%-20s\t\$rowAddHrH=%-20s\n", $sizeNb, $colHhH, $idHrH, $rowAddHrH;
    	#print "\n\n\&print_all_sub_array(\$rowAddHrH)\n"; &print_all_sub_array($rowAddHrH); print "\n\n";
                 
                    #####@{ $tempSgArray[3] } �����⣡����������������������������
      foreach my $ecRowHr (  @{ $tempSgArray[3] }  ){ printf "\$ecRowHr=%-20s\t\$eachOrededGroup_Order=%-20d\n", $ecRowHr, $eachOrededGroup_Order;
      	if (defined ($bigScore_forEach_row->{$ecRowHr}->{$eachOrededGroup_Order}) ){
      	  $bigScore_forEach_row->{$ecRowHr}->{$eachOrededGroup_Order}+=$SizeScore;   printf "\$SizeScore=%-40s\t%-80s\t%-50s\n", $SizeScore, "Add numbers   \$bigScore_forEach_row->{\$ecRowHr}->{\$eachOrededGroup_Order}+=\$SizeScore=\$bigScore_forEach_row->{$ecRowHr}->{$eachOrededGroup_Order}=",$bigScore_forEach_row->{$ecRowHr}->{$eachOrededGroup_Order};
        }
        else {
        	$bigScore_forEach_row->{$ecRowHr}->{$eachOrededGroup_Order}=$SizeScore;    printf "\$SizeScore=%-40s\t%-80s\t%-50s\n", $SizeScore, "First Defined  \$bigScore_forEach_row->{\$ecRowHr}->{\$eachOrededGroup_Order} =\$SizeScore=\$bigScore_forEach_row->{$ecRowHr}->{$eachOrededGroup_Order}=",$bigScore_forEach_row->{$ecRowHr}->{$eachOrededGroup_Order};
        }
      }  unshift @tempSgArray, $SizeScore;  printf "%-20s\t\t", "\$i=$i"; map {printf "%-10s\t", $_;} @tempSgArray; print "\n";  print "\cl\n\n\n"; 
      if    ($sizeNb> $lastSizeNb){$SizeScore=$SizeScore*1.2;                               printf "2   \$SizeScore=\t%-50s\t%50s  > %-50s\n", $SizeScore, $sizeNb, $lastSizeNb;}
      elsif ($sizeNb==$lastSizeNb){$SizeScore=$SizeScore*1.1;                             printf "1.1 \$SizeScore=\t%-50s\t%50s == %-50s\n", $SizeScore, $sizeNb, $lastSizeNb;}
      else  {die "\$sizeNb=$sizeNb\n\$lastSizeNb=$lastSizeNb\n";}
      $lastSizeNb  = $sizeNb;
    } 
    
    for (my $i=0; $i<$TotalRowNb; $i++){
      if (defined ($bigScore_forEach_row->{$i}->{$eachOrededGroup_Order}) ){  printf "%-100s\t%-50s\n", "Defined    \$bigScore_forEach_row->{\$eachOrededGroup_Order}->{\$i}=\$bigScore_forEach_row->{$i}->{$eachOrededGroup_Order}=", $bigScore_forEach_row->{$i}->{$eachOrededGroup_Order}; }
      else {$bigScore_forEach_row->{$i}->{$eachOrededGroup_Order}=0;          printf "%-100s\t%-50s\n", "Not Defined\$bigScore_forEach_row->{\$eachOrededGroup_Order}->{\$i}=\$bigScore_forEach_row->{$i}->{$eachOrededGroup_Order}=", $bigScore_forEach_row->{$i}->{$eachOrededGroup_Order}; }
    }
  
  }
  
  #�������䣬����һ����������е��кŽ������У��������ݸ������塢domain�����������෽���������ݱȷ����У����û�бȷְ�ԭhtml������sheet��˳������ 
  my $finalChangeRowHash; my $finalNewRowNb=0;
  foreach my $roKey ( sort { 
  	                             ($bigScore_forEach_row->{$b}->{0}  <=> $bigScore_forEach_row->{$a}->{0})  
  	                         ||  ($bigScore_forEach_row->{$b}->{1}  <=> $bigScore_forEach_row->{$a}->{1})
  	                         ||  ($bigScore_forEach_row->{$b}->{2}  <=> $bigScore_forEach_row->{$a}->{2})
  	                         ||  ($bigScore_forEach_row->{$b}->{3}  <=> $bigScore_forEach_row->{$a}->{3})
  	                         ||  ($bigScore_forEach_row->{$b}->{4}  <=> $bigScore_forEach_row->{$a}->{4})
  	                         ||  ($bigScore_forEach_row->{$b}->{5}  <=> $bigScore_forEach_row->{$a}->{5})
  	                         ||  ($bigScore_forEach_row->{$b}->{6}  <=> $bigScore_forEach_row->{$a}->{6})
  	                         ||  ($bigScore_forEach_row->{$b}->{7}  <=> $bigScore_forEach_row->{$a}->{7})
  	                         ||  ($bigScore_forEach_row->{$b}->{8}  <=> $bigScore_forEach_row->{$a}->{8})
  	                         ||  ($bigScore_forEach_row->{$b}->{9}  <=> $bigScore_forEach_row->{$a}->{9})
  	                         ||  ($bigScore_forEach_row->{$b}->{10} <=> $bigScore_forEach_row->{$a}->{10})
  	                         ||  ($bigScore_forEach_row->{$b}->{11} <=> $bigScore_forEach_row->{$a}->{11})
  	                         ||  ($bigScore_forEach_row->{$b}->{12} <=> $bigScore_forEach_row->{$a}->{12})
  	                         ||  ($bigScore_forEach_row->{$b}->{13} <=> $bigScore_forEach_row->{$a}->{13})
  	                         ||  ($bigScore_forEach_row->{$b}->{14} <=> $bigScore_forEach_row->{$a}->{14})
  	                         ||  ($bigScore_forEach_row->{$b}->{15} <=> $bigScore_forEach_row->{$a}->{15})
  	                         ||  ($bigScore_forEach_row->{$b}->{16} <=> $bigScore_forEach_row->{$a}->{16})
  	                         ||  ($bigScore_forEach_row->{$b}->{17} <=> $bigScore_forEach_row->{$a}->{17})
  	                         ||  $a                                 <=> $b
  	                         
  	                       } ( keys ( %{ $bigScore_forEach_row } ) ) ){
    
    #print "\n\n\n\n\cl\cl";printf "\n\n\n%-20s\t%-20s\n","\$roKey=", $roKey; 
    #map {printf "%40s\t%-10.2f\n","\$bigScore_forEach_row->{$roKey}->{$_}=",$bigScore_forEach_row->{$roKey}->{$_};} (0 .. 17); print "\n\n";
    $finalChangeRowHash->{$roKey}=$finalNewRowNb;  #printf "%-70s\t%-50s\t%-70s\t%-50s\n", "\$finalChangeRowHash->{\$roKey}=\$finalChangeRowHash->{$roKey}=",$finalNewRowNb, "\$bigScore_forEach_row->{\$roKey}=\$bigScore_forEach_row->{$roKey}=",$bigScore_forEach_row->{$roKey};
    $finalNewRowNb++;
  }
  
  my $sortedHTMLhash=&regulateTheRowOFsheetHASH($inpt_htmlSheetHash, $finalChangeRowHash, $big_link_RowCol_2_url);
  
  return  $sortedHTMLhash;

}


###sub8.10.1          GetFreeCpuNb{   # ��ÿ���CPU�Ľ�����
sub GetFreeCpuNb{   # ��ÿ���CPU�Ľ�����
  my $mpTime=7;
  my $freeCPUnb=0;  my $idle=1;   my $cpuNb;
  my $cpuStat=` mpstat 1 $mpTime`;    #print "\$cpuStat=$cpuStat\n";
  if ($cpuStat=~m/\((\d+) CPU\).*Average:(?:\s+\S+)+\s+(\S+)\s?$/sm){
    $cpuNb=$1;            
   	my $idle=$2;
   	#my $present=$idle;
  	$freeCPUnb=$cpuNb*$idle/100;      
  	warn "\n\$cpuNb=$cpuNb\t       idle=     ${idle}%\n\n\$freeCPUnb=$freeCPUnb\n\n\n";
    #printf "\cl\n\n\n%-10s%-10s%-10s%-10s\n\n", "Node:", $noteAdd,"Total:",$cpuNb;
    #printf "%-10s%-10s\n%-10s%-10s\n\n\n",   "${present}%", "CPU LEFT\n", "$freeCPUnb" , "CPU FREE"; 
  }
  else {
    die "\n\nthe mpstat  is not correct,\n the command is \n mpstat 1 $mpTime the output is  \n$cpuStat\n\n\n";
  }
  
  return $freeCPUnb;
}
###sub8.10.2          regulateTheRowOFsheetHASH{    #���� �任������hash���Լ��滻cluastalw��hyperlink��hash����ԭ����hash���������޸�
sub regulateTheRowOFsheetHASH{    #���� �任������hash���Լ��滻cluastalw��hyperlink��hash����ԭ����hash���������޸�
  my ($inFMThash, $changROWhash, $changeCLUSTALWhyperLINKhash)=@_;
  
                                                                                                                                                                                                                                                                                                                                                                  #print "\n\cl\nprint_all_sub_array(\$changeCLUSTALWhyperLINKhash=$changeCLUSTALWhyperLINKhash)\n\cl\n";  print_all_sub_array($changeCLUSTALWhyperLINKhash) ; print "\cl\cl\cl\n\n\n";
  #��һ��������clustalw������� ��������
  my $newFMThash;#�����µ�hash������������
  foreach my $fmtKEY (   keys (  %{ $inFMThash }  )   ){
    foreach my $addressOFarray (  @{ $inFMThash->{$fmtKEY} }  ){
      my ($inROW, $inCOL, $inVAL)=@{ $addressOFarray };
      if (defined ( $changeCLUSTALWhyperLINKhash->{$inROW}->{$inCOL} )){
      	if ($changeCLUSTALWhyperLINKhash->{$inROW}->{$inCOL}=~m/\S+/){
      		my $labelINinVAL    =ExcelHandle::GetLableOfHpLk($inVAL);                                                                                                                                                                                                                                              #print "bbbbbbbbbbbb55555555555555\$labelINinVAL=$labelINinVAL\t\$inVAL=$inVAL\n";
      	  my $hyperlinkINinVAL=ExcelHandle::GetHpLk($inVAL);                                                                                                                                                                                                                                                     #print "bbbbbbbbbbbb66666666666666\$hyperlinkINinVAL=$hyperlinkINinVAL\t\$inVAL=$inVAL\n";  
      	  my $newHPLINK       =$changeCLUSTALWhyperLINKhash->{$inROW}->{$inCOL};  #˵��������棬�е���=hyplink����ʽ
      	  my $newVAL          =ExcelHandle::makeHyperLk($newHPLINK,$labelINinVAL);                                                                                                                                                                                                                               #print "bbbbbbbbbbbb77777777777777\$newVAL=$newVAL\t\$inVAL=$inVAL\n";  
      	  $newFMThash=ExcelHandle::PushCellIntoSheeetHash($newFMThash, $inROW, $inCOL, $newVAL,$fmtKEY);                                                                                                                                                                                                                                                  #print "aaaaaaaaaaaaa11111\ExcelHandle::PushCellIntoSheeetHash(\$newFMThash=$newFMThash, \$inROW=$inROW, \$inCOL=$inCOL, \$newVAL=$newVAL, \$fmtKEY=$fmtKEY)\n";
        }
      }
      else {
      	$newFMThash=ExcelHandle::PushCellIntoSheeetHash($newFMThash, $inROW, $inCOL, $inVAL,$fmtKEY);                                                                                                                                                                                                                                                     #print "aaaaaaaaaaaaa22222\ExcelHandle::PushCellIntoSheeetHash(\$newFMThash=$newFMThash, \$inROW=$inROW, \$inCOL=$inCOL, \$inVAL=$inVAL, \$fmtKEY=$fmtKEY)\n";
      }
    }
  }
  
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  #print "\n\cl\nprint_all_sub_array(\$changROWhash=$changROWhash)\n\cl\n";  print_all_sub_array($changROWhash) ; print "\cl\cl\cl\n\n\n";
  #�ڶ�������������˳�򣬱仯row��ֵ�������µ�hash�С�
  my $changeROWfmtHASH;   #�����µ�hash������������
  foreach my $fmt2KEY (   keys (  %{ $newFMThash }  )   ){
    foreach my $address2OFarray (  @{ $newFMThash->{$fmt2KEY} }  ){
      my ($in2ROW, $in2COL, $in2VAL)=@{ $address2OFarray };
      if (defined ( $changROWhash->{$in2ROW} )){
      	if ($changROWhash->{$in2ROW}=~m/\d+/){                                                                                                                                                                                                                                                                                                                                                                   #print "\$changROWhash->{\$in2ROW}=\$changROWhash->{$in2ROW}=$changROWhash->{$in2ROW}\n";
      		my $sortedNEWrow=$changROWhash->{$in2ROW};
      	  $changeROWfmtHASH=ExcelHandle::PushCellIntoSheeetHash($changeROWfmtHASH, $sortedNEWrow, $in2COL, $in2VAL,$fmt2KEY);                                                                                                                                                                                                 #print "aaaaaaaaaaaaa333333ExcelHandle::PushCellIntoSheeetHash(\$changeROWfmtHASH=$changeROWfmtHASH, \$sortedNEWrow=$sortedNEWrow, \$in2COL=$in2COL, \$in2VAL=$in2VAL, \$fmt2KEY=$fmt2KEY);\n";
      	}
      	
      	
      }
      else {
      	$changeROWfmtHASH=ExcelHandle::PushCellIntoSheeetHash($changeROWfmtHASH, $in2ROW, $in2COL, $in2VAL,$fmt2KEY);                                                                                                                                                                                                                                                 #print "aaaaaaaaaaaaa444444ExcelHandle::PushCellIntoSheeetHash(\$changeROWfmtHASH=$changeROWfmtHASH, \$in2ROW=$in2ROW, \$in2COL=$in2COL, \$in2VAL=$in2VAL, \$fmt2KEY=$fmt2KEY);\n";
      }
    }
  }
  
  return $changeROWfmtHASH;
}

1;
