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

use DirFileHandle;
use TimeWork;
use ExcelHandle;
use ForeachHash;


                               
                      
package ProteinFamilyWork;

sub BuildProteinFamilyHash_FromDir_with_ProtFmlFiles{   # my $outHash=ProteinFamilyWork::BuildProteinFamilyHash_FromDir_with_ProtFmlFiles($inDir); 
	my ($inDir)=@_;
	
	my $warnMsgBody="\nIn package  ProteinFamilyWork,\tIn sub BuildProteinFamilyHash_FromDir_with_ProtFmlFiles,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	my $outHash;
	if  (   (  defined ( $inDir )  ) && (  -d ( $inDir )  )   ){
		my $inDirArray=DirFileHandle::getDirArray ($inDir);
		if  (   (  defined ( $inDirArray )  ) && (  ( ref ($inDirArray) ) eq 'ARRAY'  )   ){
			foreach my $ecProtFmlyFile (  @{ $inDirArray }  ){ 
				my $wholeFilePath=$inDir."/$ecProtFmlyFile";
        my $SeqHash=FastaFileHandle::BuildHashFromFastaFile_seqID_as_Key( $wholeFilePath );    #DirFileHandle::PrintAndWarnDumper ($SeqHash, "\$wholeFilePath=$wholeFilePath\t\$SeqHash=$SeqHash\n");
        if (   (  defined ( $SeqHash )  ) && (  ( ref ($SeqHash) ) eq 'HASH'  )   ){
          $outHash->{$ecProtFmlyFile}=Storable::dclone( $SeqHash );
        }
        
      }
		}
		
	}
	return $outHash;
}


sub build_working_DirAndDatabase{  #my $outHash=ProteinFamilyWork::build_working_DirAndDatabase($in_protein_hash, $inDatabaseDir, $testProtein); 
	
	my ($in_protein_hash, $inDatabaseDir, $testProtein)=@_;
	
	my $warnMsgBody="\nIn package  ProteinFamilyWork,\tIn sub build_working_DirAndDatabase,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	
	my ( $ProteinDIRe_0_0_0_key, $PureProtDIR_0_0_1_key, $ProteinFile_0_0_2_Key )
	  =( '0_0_0_ProteinDIRe',    '0_0_1_ProteinFile',    '0_0_2_ProteinFile'    );
	                      
  my $inHash=$in_protein_hash;
  
  my $outHash;
  
  if (   (  defined ( $testProtein )  ) && ($testProtein=~m/\S+/)   ){ 
  	my @testPtArray=split ',',$testProtein;
  	if (  ( ref ($inHash) ) eq 'HASH' ){
  	  foreach my $eacPt  ( @testPtArray ){
  	    my $hashIdx=0;
        foreach my $keyLev_0 (    sort { $a cmp $b } (   keys (  %{ $inHash } )   )    ){     #print "    \$keyLev_0=$keyLev_0\n";
        	if ( ($keyLev_0=~/\S+/) && ($keyLev_0 eq $eacPt) ) {                       #warn "\n\n\n\n$keyLev_0=$keyLev_0\n\n\n";
        	  my $valLev_0=$inHash->{$keyLev_0};                                                #print "    \$valLev_0=$valLev_0\n";
        	  my $refLev_0=ref ($valLev_0);                                                     #print "    \$refLev_0=$refLev_0\n\n";  
        	  
        	  my $formatHashIdx=sprintf ("%03d",$hashIdx);
        	  my $purePro_dir=$formatHashIdx."__".$keyLev_0;
        	  my $protein_dir=$inDatabaseDir."/".$purePro_dir;    system ("mkdir -p $protein_dir"); 
        	  $outHash->{$keyLev_0}->{$ProteinDIRe_0_0_0_key}=$protein_dir;
        	  $outHash->{$keyLev_0}->{$PureProtDIR_0_0_1_key}=$purePro_dir;
        	  my $proteinFasta;
            if ( $refLev_0 eq 'HASH' ){
              foreach my $keyLev_1 (    sort { $a cmp $b } (   keys (  %{ $valLev_0 } )   )    ){   #print "      \$keyLev_1=$keyLev_1\n";
              	my $valLev_1=$valLev_0->{$keyLev_1};                                                #print "      \$valLev_1=$valLev_1\n";                
              	#my $refLev_1=ref ($valLev_1);                                                       print "      \$refLev_1=$refLev_1\n\n";
              	$proteinFasta.=">".$keyLev_1."\n".$valLev_1."\n\n";
              	$outHash->{$keyLev_0}->{'0_0_4_Name_seqHSH'}->{$keyLev_1}=$valLev_1;
              }
            }
            
            my $fastaFile=$protein_dir."/".$keyLev_0.".txt";
            $outHash->{$keyLev_0}->{$ProteinFile_0_0_2_Key}=$fastaFile;
            $outHash->{$keyLev_0}->{'0_0_3_fastaSeqIfm'}=$proteinFasta;
            
            FastaFileHandle::BuildFastaFile_withFastaString($fastaFile, $proteinFasta);
            my $formatdbCMD="formatdb -i $fastaFile -p T -o T";
            system ("$formatdbCMD");
          }
          $hashIdx++;
        }
      }
    }
  }
  else {
  	if (  ( ref ($inHash) ) eq 'HASH' ){
    	my $hashIdx=0;
      foreach my $keyLev_0 (    sort { $a cmp $b } (   keys (  %{ $inHash } )   )    ){     #print "    \$keyLev_0=$keyLev_0\n";
      	if ($keyLev_0=~m/\S+/) { 
      	  my $valLev_0=$inHash->{$keyLev_0};                                                #print "    \$valLev_0=$valLev_0\n";
      	  my $refLev_0=ref ($valLev_0);                                                     #print "    \$refLev_0=$refLev_0\n\n";  
      	  
      	  my $formatHashIdx=sprintf ("%03d",$hashIdx);
      	  my $purePro_dir=$formatHashIdx."__".$keyLev_0;
      	  my $protein_dir=$inDatabaseDir."/".$purePro_dir;    system ("mkdir -p $protein_dir"); 
      	  $outHash->{$keyLev_0}->{$ProteinDIRe_0_0_0_key}=$protein_dir;
      	  $outHash->{$keyLev_0}->{$PureProtDIR_0_0_1_key}=$purePro_dir;
      	  my $proteinFasta;
          if ( $refLev_0 eq 'HASH' ){
            foreach my $keyLev_1 (    sort { $a cmp $b } (   keys (  %{ $valLev_0 } )   )    ){   #print "      \$keyLev_1=$keyLev_1\n";
            	my $valLev_1=$valLev_0->{$keyLev_1};                                                #print "      \$valLev_1=$valLev_1\n";                
            	#my $refLev_1=ref ($valLev_1);                                                       print "      \$refLev_1=$refLev_1\n\n";
            	$proteinFasta.=">".$keyLev_1."\n".$valLev_1."\n\n";
            	$outHash->{$keyLev_0}->{'0_0_4_Name_seqHSH'}->{$keyLev_1}=$valLev_1;
            }
          }
          
          my $fastaFile=$protein_dir."/".$keyLev_0.".txt";
          $outHash->{$keyLev_0}->{$ProteinFile_0_0_2_Key}=$fastaFile;
          $outHash->{$keyLev_0}->{'0_0_3_fastaSeqIfm'}=$proteinFasta;
          
          FastaFileHandle::BuildFastaFile_withFastaString($fastaFile, $proteinFasta);
          my $formatdbCMD="formatdb -i $fastaFile -p T -o T";
          system ("$formatdbCMD");
          $hashIdx++;
        }
        
        
      }
      
      if (   (  defined ( $outHash  )  ) && (  ref ( $outHash  ) eq 'HASH'  )   ){                                                                      #print "20181204-m-1 \$outHash=$outHash\n";
      	foreach my $pt_key (    sort { $a cmp $b } (   keys (  %{ $outHash } )   )    ) {                                                               #print "20181204-m-2 \$pt_key=$pt_key\n";
      	}
      }
    
      foreach my $keyLev_0 (    sort { $a cmp $b } (   keys (  %{ $inHash } )   )    ){                                                                    #print "    \$keyLev_0=$keyLev_0\n";
      	if ($keyLev_0 eq '') { 
      	  my $valLev_0=$inHash->{$keyLev_0};                                                                                                               #print "    \$valLev_0=$valLev_0\n";
      	  my $refLev_0=ref ($valLev_0);                                                                                                                    #print "    \$refLev_0=$refLev_0\n\n";  
      	  
      	  if ( $refLev_0 eq 'HASH' ){
            foreach my $keyLev_1 (    sort { $a cmp $b } (   keys (  %{ $valLev_0 } )   )    ){                                                                  #print "      \$keyLev_1=$keyLev_1\n";
            	my $valLev_1=$valLev_0->{$keyLev_1};                                                                                                               #print "      \$valLev_1=$valLev_1\n";                
            	
            	if ($keyLev_1=~m/\S\S.\S\S-(\S+)-\d+-\S/){
            		my $ptNm=$1;
            		
            		my $proteinFasta=">".$keyLev_1."\n".$valLev_1."\n\n";                                                                                                                                #print "20181204-0 \$keyLev_1=$keyLev_1   \$ptNm=$ptNm\n";
            		
            		if (   (  defined ( $outHash->{$ptNm} )  ) && (  ref ( $outHash->{$ptNm} ) eq 'HASH'  )   ){                                                                #print "20181204-1 \$outHash->{$ptNm}=$outHash->{$ptNm}\n";
            			if (   (  defined ( $outHash->{$ptNm}->{'0_0_3_fastaSeqIfm'} )  ) && ( $outHash->{$ptNm}->{'0_0_3_fastaSeqIfm'}=~m/\S+/ )   ){
            			  $outHash->{$ptNm}->{'0_0_3_fastaSeqIfm'}.=$proteinFasta;                                                                              #print "20181204-2 \$outHash->{$ptNm}->{'0_0_3_fastaSeqIfm'}=$outHash->{$ptNm}->{'0_0_3_fastaSeqIfm'}\n";
            			}
            			else{
            				$outHash->{$ptNm}->{'0_0_3_fastaSeqIfm'}=$proteinFasta;                                                                 #print "20181204-3 \$outHash->{$ptNm}->{'0_0_3_fastaSeqIfm'}=$outHash->{$ptNm}->{'0_0_3_fastaSeqIfm'}\n";
            			}
            			$outHash->{$ptNm}->{'0_0_5_add_workTdo'}=1;                                                                 #print "20181204-3-2 \$outHash->{$ptNm}->{'0_0_5_add_workTdo'}=$outHash->{$ptNm}->{'0_0_5_add_workTdo'}\n";
            		}
            		else {
            			my $formatHashIdx=sprintf ("%03d",$hashIdx);  $hashIdx++;
      	          my $purePro_dir=$formatHashIdx."__".$ptNm;
      	          my $protein_dir=$inDatabaseDir."/".$purePro_dir;    system ("mkdir -p $protein_dir"); 
      	          $outHash->{$ptNm}->{$ProteinDIRe_0_0_0_key}=$protein_dir;                                                                  #print "20181204-4 \$outHash->{$ptNm}->{$ProteinDIRe_0_0_0_key}=$outHash->{$ptNm}->{$ProteinDIRe_0_0_0_key}\n";
          	      $outHash->{$ptNm}->{$PureProtDIR_0_0_1_key}=$purePro_dir;                                                                  #print "20181204-5 \$outHash->{$ptNm}->{$PureProtDIR_0_0_1_key}=$outHash->{$ptNm}->{$PureProtDIR_0_0_1_key}\n";
          	      
          	      if (   (  defined ( $outHash->{$ptNm}->{'0_0_3_fastaSeqIfm'} )  ) && ( $outHash->{$ptNm}->{'0_0_3_fastaSeqIfm'}=~m/\S+/ )   ){
            			  $outHash->{$ptNm}->{'0_0_3_fastaSeqIfm'}.=$proteinFasta;                                                                 #print "20181204-6 \$outHash->{$ptNm}->{'0_0_3_fastaSeqIfm'}=$outHash->{$ptNm}->{'0_0_3_fastaSeqIfm'}\n";
            			}
            			else{
            				$outHash->{$ptNm}->{'0_0_3_fastaSeqIfm'}=$proteinFasta;                                                                #print "20181204-7 \$outHash->{$ptNm}->{'0_0_3_fastaSeqIfm'}=$outHash->{$ptNm}->{'0_0_3_fastaSeqIfm'}\n";
            			}
          	      $outHash->{$ptNm}->{'0_0_5_add_workTdo'}=1;                                                                 #print "20181204-8 \$outHash->{$ptNm}->{'0_0_5_add_workTdo'}=$outHash->{$ptNm}->{'0_0_5_add_workTdo'}\n";
            		}
            		$outHash->{$ptNm}->{'0_0_4_Name_seqHSH'}->{$keyLev_1}=$valLev_1;
            	}
            }
          }
          
          
        }
        
        
      }
      
      if (   (  defined ( $outHash  )  ) && (  ref ( $outHash  ) eq 'HASH'  )   ){                                                                #print "20181204-9 \$outHash=$outHash\n";
      	foreach my $pt_key (    sort { $a cmp $b } (   keys (  %{ $outHash } )   )    ) {                                                                #print "20181204-a \$pt_key=$pt_key\n";
      		if (   (  defined ( $outHash->{$pt_key}->{'0_0_5_add_workTdo'} )  ) && ( $outHash->{$pt_key}->{'0_0_5_add_workTdo'}==1 )   ){                                                                #print "20181204-b \$outHash->{$pt_key}->{'0_0_5_add_workTdo'}=$outHash->{$pt_key}->{'0_0_5_add_workTdo'}\n";
      			my $fastaFile=$outHash->{$pt_key}->{$ProteinDIRe_0_0_0_key}."/".$pt_key.".txt";                                                                #print "20181204-c \$outHash->{$pt_key}->{$ProteinDIRe_0_0_0_key}=$outHash->{$pt_key}->{$ProteinDIRe_0_0_0_key}\n";
            $outHash->{$pt_key}->{$ProteinFile_0_0_2_Key}=$fastaFile;
            my $proteinFasta=$outHash->{$pt_key}->{'0_0_3_fastaSeqIfm'};                                                                 #print "20181204-d \$outHash->{$pt_key}->{'0_0_3_fastaSeqIfm'}=$outHash->{$pt_key}->{'0_0_3_fastaSeqIfm'}\n";
            FastaFileHandle::BuildFastaFile_withFastaString($fastaFile, $proteinFasta);  
            my $formatdbCMD="formatdb -i $fastaFile -p T -o T";
            system ("$formatdbCMD");
      		}
        }
      }
    }
  }
  
	
	return $outHash;
}



1;