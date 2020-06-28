
#!/usr/bin/perl -w

BEGIN{  push (@INC, '/home/fredjiang/fredProgs/PMs')  };

use strict;
use warnings;
use Data::Dumper;
$Data::Dumper::Purity=1;
use Storable ;
use Storable qw(nstore store_fd nstore_fd freeze thaw dclone);


package  AlgaeFinalProteinSeqHandle;







#这个用于从0_SelProteinDATABASEhash.txt 文件，建立一个以sequence作为Key的hash

sub UseSequenceAsKeyFromIDkeyHash{
  my ($inHashFile)=@_;
  my $inHash=Storable::retrieve ($inHashFile);
  my $outHash;
  if (  ( ref ($inHash) ) eq 'HASH' ){
    foreach my $keyLev_0 (    sort { $a cmp $b } (   keys (  %{ $inHash } )   )    ){   #print "    \$keyLev_0=$keyLev_0\n";
    	my $valLev_0=$inHash->{$keyLev_0};                                                #print "    \$valLev_0=$valLev_0\n";
    	my $sequnesHere=$inHash->{$keyLev_0}->{'_FastaPepSequ20180402'};
    	$sequnesHere=~s/\s//g; $sequnesHere=~s/\n//g;
    	$sequnesHere=~s/\*/X/g;
    	$sequnesHere=~s/x/X/g;
    	if (  defined ( $outHash->{$sequnesHere} )  ){
    	  push @{ $outHash->{$sequnesHere} }, $keyLev_0;
    	}
    	else {
    	  $outHash->{$sequnesHere}->[0]=$keyLev_0;
    	}
    
    	
    }
  }
  return $outHash;;
}


sub buildSortedHashFromOneColumnTxtFile{  #利用输入的只有1列的数据，建立 键为每行数据，值为行号的hash
  my ($inFile)=@_;
  open (IN, $inFile) or die "\n\nDIE:In Sub buildSortedHashFromOneColumnTxtFile,\nCannot open \$inFile=$inFile :$!\n\n";
  my $outHash;
  my $tpMakr=$/;
  $/="\n";  #my @tpAr=split ("\n",$inLines);
  my $idx=1;
  while (<IN>){
  my $tpVal=$_;
  	chomp $tpVal;
  	$tpVal=~s/\#.*$//;
    if ($tpVal=~m/\S+/){
      $tpVal=~s/^\s*//; $tpVal=~s/\s*$//;
      $outHash->{$tpVal}=$idx;
      $idx++;
    }
  }
  close(IN);
  $/=$tpMakr;
  return $outHash; 
}

sub GetProteinFmlNameIndex{  #这个函数 获得蛋白家族的名字 等信息，做成Hash
  my ($inFile)=@_;
  my $outHash;
  open (INF, "$inFile") or die "Cannot open \$inFile=$inFile :$!\n\n";
  my $tp=$/; $/="\n";
  while (<INF>){
  	$_=~s/^\s*//;
  	$_=~s/\s*$//;
  	if (m/^#/){}
  	elsif (m/\S+/){
	    my @tpAr=split("\t",$_);  #SelR	16	SelR
      my @handledAr; my $i=0;
      foreach my $eachCell (@tpAr){ chomp $eachCell; $eachCell=~s/^\s*//; $eachCell=~s/\s*$//; $handledAr[$i]=$eachCell; $i++;}
      
      $outHash->{ $handledAr[0] }->{'RealName'}  =$handledAr[4];         print "SSSSSSSSSSSSSS\$outHash->{ \$handledAr[0] }->{'RealName'}  =\$handledAr[4]=\$outHash->{ $handledAr[0] }->{'RealName'}  =$handledAr[4]\n";
      $outHash->{ $handledAr[0] }->{'TreeIdx'}   =$handledAr[1];
      $outHash->{ $handledAr[0] }->{'BigFml'}    =$handledAr[3];
      $outHash->{ $handledAr[0] }->{'SmallFml'}  =$handledAr[4];
      $outHash->{ $handledAr[0] }->{'longName'}  =$handledAr[2];     
                                                                           
        
    }
    
  }
  $/=$tp;
  close(INF);
  return $outHash;
}

sub hashHadlesub{  #解析特定的 物种相关文件，生成带顺序的hash
  my ($inFile)=@_;
  my $outHash;
  open (INF, "$inFile") or die "Cannot open \$inFile=$inFile :$!\n\n";
  my $specisesIDX=1;
  while (<INF>){
  	$_=~s/^\s*//;
  	$_=~s/\s*$//;
  	if (m/^#/){}
  	elsif (m/\S+/){
	    my @tpAr=split("\t",$_);  #Bigelowiella natans	Bigelowiella_natans20150522	3	93442513	93442.513	93.4	3740	Bigelowiella_natans20150522	2561149	2561.149	2.6	3464	Bigelowiella_natans20150522	1716130	1716.13	1.7	2402
      my @handledAr; my $i=0;
      foreach my $eachCell (@tpAr){ chomp $eachCell; $eachCell=~s/^\s*//; $eachCell=~s/\s*$//; $handledAr[$i]=$eachCell; $i++;}
      $outHash->{ $handledAr[1] }->{'TreeIdx'}        =$specisesIDX;
      $outHash->{ $handledAr[1] }->{'RealName'}       =$handledAr[0];  print "SSSSSSSSSSSSSS\$outHash->{ $handledAr[1] }->{'RealName'}=$outHash->{ $handledAr[1] }->{'RealName'}\n";
      $outHash->{ $handledAr[1] }->{'GenoFileSZ'}     =$handledAr[3];
      $outHash->{ $handledAr[1] }->{'GenoFileSZ_K'}   =$handledAr[4];
      $outHash->{ $handledAr[1] }->{'GenoFileSZ_M'}   =$handledAr[5];
      $outHash->{ $handledAr[1] }->{'GenoFileSeqNB'}  =$handledAr[6];
      $outHash->{ $handledAr[1] }->{'OrgEstFileSZ'}   =$handledAr[8];
      $outHash->{ $handledAr[1] }->{'OrgEstFileSZ_K'} =$handledAr[9];
      $outHash->{ $handledAr[1] }->{'OrgEstFileSZ_M'} =$handledAr[10];
      $outHash->{ $handledAr[1] }->{'OrgEstFileSeqNB'}=$handledAr[11]; 
      $outHash->{ $handledAr[1] }->{'AsmEstFileSZ'}   =$handledAr[13];      
      $outHash->{ $handledAr[1] }->{'AsmEstFileSZ_K'} =$handledAr[14];      
      $outHash->{ $handledAr[1] }->{'AsmEstFileSZ_M'} =$handledAr[15];      
      $outHash->{ $handledAr[1] }->{'AsmEstFileSeqNB'}=$handledAr[16];      
      $outHash->{ $handledAr[1] }->{'ShortNm'        }=$handledAr[17];  
                                                                           
      $specisesIDX++;
    }
  }
  close(INF);
  return $outHash;
}


1;