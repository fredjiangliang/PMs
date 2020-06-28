
#!/usr/bin/perl -w

BEGIN{  push (@INC, '/home/fredjiang/fredProgs/PMs')  };

use strict;
use warnings;

use Storable ;
use Storable qw(nstore store_fd nstore_fd freeze thaw dclone);
use DirFileHandle;
use DieWork;



package  ArrayHashChange;


sub Segment_bigArray_into_smallArrayes{    #   my $Segmented_array=ArrayHashChange::Segment_bigArray_into_smallArrayes($inArrayRef, $smallArraySize); 
	# 将大的array分成很多个小的array，比如原来的array有10000个元素，可以分成10个2维array，每个小array中有1000个元素.
	my ($inArrayRef, $smallArraySize)=@_;
	
	my $warnMsgBody="\nIn package  ArrayHashChange,\tIn sub Segment_bigArray_into_smallArrayes,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $caller_inform=DirFileHandle::print_SubCallerInform;
  
  my $Segmented_array;
  if (   (  defined ( $inArrayRef )  ) &&  (  ref ( $inArrayRef ) eq 'ARRAY'  ) && (  defined ( $smallArraySize )  ) &&  ( $smallArraySize=~m/\d+/) &&  ( $smallArraySize >= 1)   ){
  
    my $part_idx=0; my $inner_idx=0; 	  
	  for (  my $i=0; $i<@{ $inArrayRef }; $i++  ){
	  	
	  	$Segmented_array->[$part_idx]->[$inner_idx]=$inArrayRef->[$i];
	  	if (  $inner_idx >= ( $smallArraySize-1 )  ){
	  		
	  		$inner_idx=0;
	  		$part_idx++;
	  	}
	  	else{
	  		$inner_idx++;
	  	}
	  	
	  }
  	
  }
  else{
  	DieWork::Just_dieWork( $die_MsgHead."\n\$inArrayRef=$inArrayRef should be a ARRAY ref \$smallArraySize=$smallArraySize should be a number >= 1 !!!: $!".$caller_inform );
  }	
  return $Segmented_array;
	
}


sub PushSmallHash_into_bigHash{    #my $bigHash; $bigHash=ArrayHashChange::PushSmallHash_into_bigHash($bigHash, $smlHash);
	my ($bigHash, $smlHash)=@_;
	
	my $warnMsgBody="\nIn package  ArrayHashChange,\tIn sub PushSmallHash_into_bigHash,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	if (   (  defined ( $smlHash )  ) && (  ref ( $smlHash ) eq 'HASH'  )   ){
		foreach my $eachKey (   keys ( %{ $smlHash }  )   ){                   #DirFileHandle::PrintAndWarnDumper ($smlHash->{$eachKey}, "\$smlHash->{$eachKey}=$smlHash->{$eachKey}\n");
			$bigHash->{$eachKey}=$smlHash->{$eachKey};
		}
	}
	
	return $bigHash;
	
}


sub BuildReverseHash{  #  ArrayHashChange::BuildReverseHash($in_Hash);
	my ($in_Hash)=@_;
	my $outHash;
	if (  ( ref ($in_Hash) ) eq 'HASH' ){
	  foreach my $eachKey (   keys (  %{ $in_Hash }  )   ){
	  	$outHash->{  $in_Hash->{ $eachKey }  }=$eachKey;
	  }
	}
	return $outHash;
}

sub Reverse_level1KEY_Level2Key_HASH{  #  ArrayHashChange::Reverse_level1KEY_Level2Key_HASH($in_Hash);
	my ($in_Hash)=@_;
	
	my $warnMsgBody="\nIn package  TaxonomyWork,\tIn sub Get_TaxID_from_StfNM,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody"; 
  my $subCallereIfm=DirFileHandle::print_SubCallerInform; #warn $subCallereIfm;
  
	my $outHash; my $isThereReal_2_level_HashForWork=0;
	if (   (  defined ( $in_Hash )  ) && (  ref ( $in_Hash ) eq 'HASH'  )   ){	
	  foreach my $level_1_key (   keys (  %{ $in_Hash }  )   ){
	  	if (   (  defined ( $in_Hash->{$level_1_key} )  ) && (  ref ( $in_Hash->{$level_1_key} ) eq 'HASH'  )   ){	
	  	  foreach my $level_2_key (   keys (  %{ $in_Hash->{$level_1_key} }  )   ){
	  	  	if (  defined ( $in_Hash->{$level_1_key}->{$level_2_key} )  ){   #DieWork::Print_and_warn( "\$in_Hash->{$level_1_key}->{$level_2_key}=$in_Hash->{$level_1_key}->{$level_2_key}\n" );
	  	  		$isThereReal_2_level_HashForWork=1;
	  	  		if (   (  DieWork::Check_Hash_or_NOT( $in_Hash->{$level_1_key}->{$level_2_key} )  ) || (  DieWork::Check_Array_or_NOT( $in_Hash->{$level_1_key}->{$level_2_key} )  )   ){
	  	  			$outHash->{$level_2_key}->{$level_1_key}=Storable::dclone( $in_Hash->{$level_1_key}->{$level_2_key} );
	  	  		}
	  	  		else{
	  	  			$outHash->{$level_2_key}->{$level_1_key}= $in_Hash->{$level_1_key}->{$level_2_key} ;
	  	  		}
	  	  		
	  	  	}
	  	  }
	  	}
	  	
	  }
	}
	if ( $isThereReal_2_level_HashForWork==1 ){
	  return $outHash;	
	}
	else {
		my $dumpInHash=DirFileHandle::ReturnDumperInform ($in_Hash);
		DieWork::Just_dieWork( $die_MsgHead.$subCallereIfm."\$in_Hash=$in_Hash\n".$dumpInHash );
	}
}


sub ChangeArrayToHash { #把数组变成hash，key是数组中的元素， 值是序数，或者(有重复值的时候)多个序数形成的数组
	my ($inarray)=@_;
	my $outHash;
	if (  ( ref ($inarray) ) eq 'ARRAY' ){
		for ( my $i=0; $i< @{ $inarray }; $i++  ){
			if (   defined (  $outHash->{ $inarray->[$i] }  )   ){
    	  push @{  $outHash->{ $inarray->[$i] }  }, $i;
    	}
    	else {
    	  $outHash->{ $inarray->[$i] }->[0]=$i;
    	}
		}
	}
	return $outHash;
}

sub Change_Hash_to_Array {   #  my $outArray=ArrayHashChange::Change_Hash_to_Array ($inHash) ;#把hash变成数组，key是数组中的元素，
	my ($inHash)=@_;
	
	my $warnMsgBody="\nIn package  ArrayHashChange,\tIn Change_Hash_to_Array,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	if  (   (  defined (  $inHash  )  ) && (  ( ref ($inHash) ) eq 'HASH' )   ){
	}
	else{
		DieWork::Just_dieWork( $die_MsgHead."\n \$inHash=$inHash should be defined and a HASH ref $!\n\n\n".$caller_inform ); 
	}
	
	my @otAry=sort(  keys (  %{ $inHash } )  );
	
	my $outArray=[ @otAry ];
	
	
	return $outArray;
}


sub CHangeArrayToEasyHash{  #把数组变成hash，key是数组中的元素的序数， 值是数组中的元素
  my ($inarray)=@_;
	my $outHash;
	if (  ( ref ($inarray) ) eq 'ARRAY' ){
		for ( my $i=0; $i< @{ $inarray }; $i++  ){
			$outHash->{$i}=$inarray->[$i] ;
		}
	}
	return $outHash; 
}

#Get a sub hash form a big hash. Which keys will be hold in the new sub hash, is up to the $keyGroup( Hash or Array).
#my $outHASH=ArrayHashChange::GetSubHASH_from_bigHASH ($keyGroup, $bigHASH);
sub GetSubHASH_from_bigHASH{  #对一个大的hash，取其子集，子集的范围由$keyGroup（可以是Hash，也可以是Array）决定，
	my ($keyGroup, $bigHASH)=@_;
	
	my $warnMsgBody="\nIn package  ArrayHashChange,\tIn GetSubHASH_from_bigHASH,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	my $subKeyHash;
	if     (   (  defined (  $keyGroup  )  ) && (  ( ref ($keyGroup) ) eq 'HASH'  )   ){ $subKeyHash=$keyGroup;	                                    }
	elsif  (   (  defined (  $keyGroup  )  ) && (  ( ref ($keyGroup) ) eq 'ARRAY' )   ){ $subKeyHash=ArrayHashChange::ChangeArrayToHash($keyGroup);	}	
	else{		DieWork::Just_dieWork(    $die_MsgHead."\n \$keyGroup=$keyGroup should be defined and a HASH or ARRAY ref $!\n\n\n".$caller_inform ); 	}
	
	if  (   (  defined (  $bigHASH  )  ) && (  ( ref ($bigHASH) ) eq 'HASH' )   ){	}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$bigHASH=$bigHASH should be defined and a HASH ref $!\n\n\n".$caller_inform );  	}
	
	my $outHASH;
	foreach my $subKey (   keys (  %{ $subKeyHash }  )   ){
		if (  defined (  $bigHASH->{$subKey}  )  ){
			$outHASH->{$subKey}=Storable::dclone( $bigHASH->{$subKey} );
		}
		else{
			DieWork::Just_dieWork( $die_MsgHead."\n \$bigHASH->{$subKey}=$bigHASH->{$subKey} should be a defined value in \$bigHASH=$bigHASH $!\n\n\n".$caller_inform );
		}
	}
	
	return $outHASH;
	
}


#my $outPut3dHash=ArrayHashChange::Change_2d_hashINTO_3d_hash ($in2dHash, $subKeyForDived);
sub Change_2d_hashINTO_3d_hash{ #使用2dhash的下层特定的key，来讲2d hash变成3d，这时候，相当于 将2d hash，依据更下层的一个key 进行切分。
	my ($in2dHash, $subKeyForDived)=@_;
	
	my $warnMsgBody="\nIn package  ArrayHashChange,\tIn Change_2d_hashINTO_3d_hash,\n\n$subKeyForDived=$subKeyForDived\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	if  (   (  defined (  $in2dHash  )  ) && (  ( ref ($in2dHash) ) eq 'HASH' )   ){	}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$in2dHash=$in2dHash should be defined and a HASH ref $!\n\n\n".$caller_inform );  	}
	
	if  (   (  defined (  $subKeyForDived  )  ) && ( $subKeyForDived=~m/\S+/ )   ){	}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$subKeyForDived=$subKeyForDived should be a not empty String $!\n\n\n".$caller_inform );  	}
	
	my $outPut3dHash;
	foreach my $lv1_key (   keys (  %{ $in2dHash }  )   ){
		
		if  (    (  defined (  $in2dHash->{$lv1_key}  )  ) && (   (  ref ( $in2dHash->{$lv1_key} )  ) eq 'HASH'  ) && (  defined (  $in2dHash->{$lv1_key}->{$subKeyForDived}  )  ) && ( $in2dHash->{$lv1_key}->{$subKeyForDived}=~m/\S+/ )    ){	}
	  else{		DieWork::Just_dieWork( $die_MsgHead."\n \$in2dHash->{$lv1_key}=$in2dHash->{$lv1_key} should be defined and a HASH ref \n\$in2dHash->{$lv1_key}->{$subKeyForDived}=$in2dHash->{$lv1_key}->{$subKeyForDived} should be a not empty string $!\n\n\n".$caller_inform );  	}
		
		$outPut3dHash->{ $in2dHash->{$lv1_key}->{$subKeyForDived} }->{ $lv1_key }=Storable::dclone( $in2dHash->{$lv1_key} );
		
	}
	
	if  (   (  defined (  $outPut3dHash  )  ) && (  ( ref ($outPut3dHash) ) eq 'HASH' )   ){	}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$outPut3dHash=$outPut3dHash should be defined and a HASH ref $!\n\n\n".$caller_inform );  	}
	
	return $outPut3dHash;
}

sub CheckHashKey_exactly_theSAME {  #  ArrayHashChange::CheckHashKey_exactly_theSAME ($Hash_1, $Hash_2);   #检查两个hash中的key是否完全相同， check the keies of 2 hash, if all of them are the same, then return 1, else return 0
	my ($Hash_1, $Hash_2, $addtionIfomr)=@_;
	
	my $warnMsgBody="\nIn package  ArrayHashChange,\tIn subCheckHashKey_exactly_theSAME,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $subCallereIfm=DirFileHandle::print_SubCallerInform;
  
  my $addtionMsg=$warnMsgBody.$subCallereIfm;  $addtionMsg=$addtionMsg.$addtionIfomr if (   (  defined ( $addtionIfomr )  ) && ( $addtionIfomr=~m/\S+/ )   ) ;
 
	my $outReal=0;
	
	if (  ( ref($Hash_1) ne 'HASH' ) && ( ref($Hash_2) ne 'HASH' )  ){
		my $warnMsg= DirFileHandle::ReturnDumperInform ($Hash_1, "\$Hash_1=$Hash_1");	
	  $warnMsg.= DirFileHandle::ReturnDumperInform ($Hash_2, "\$Hash_2=$Hash_2");			
		DieWork::Just_dieWork( $die_MsgHead."20190105-4-1 $addtionMsg \$Hash_1=$Hash_1 \$Hash_2=$Hash_2 should be two HASH!!!\n$!\n\n\n".$warnMsg );
	}
	if (  ( ref($Hash_1) eq 'HASH' ) && ( ref($Hash_2) eq 'HASH' )  ){
		
		my $size_1=( keys %{ $Hash_1 } );		my $size_2=( keys %{ $Hash_2 } );
		my $Hash_a;                     		my $Hash_b;
		
		#print "\n20190105-4-0-0-0 \$size_1=$size_1 \$size_2=$size_2 !!!\n";
		
		if (  ($size_1 ==0) && ($size_2 ==0)   ){  #print "\n20190105-4-0-0-1 \$size_1=$size_1 \$size_2=$size_2 !!!\n";
			 DieWork::Just_dieWork( $die_MsgHead."20190105-4-1 $addtionMsg \$size_1=$size_1 \$size_2=$size_2 shold be 2 number > 0\n \$Hash_1=$Hash_1 \$Hash_2=$Hash_2 should be two HASH!!!\n$!\n\n\n" );
		}
		
		elsif (  ($size_1 >0) && ($size_2 >0)   ){ #print "\n20190105-4-0-0-2 \$size_1=$size_1 \$size_2=$size_2 !!!\n";
			$Hash_a=$Hash_1;			$Hash_b=$Hash_2;    #print "\n20190105-4-0-0-3 \$size_1=$size_1 \$size_2=$size_2 !!!\n";
			if ($size_1==$size_2){ 
				my $diffrentFoumd=0;
				FOREACHMAR1: foreach my $key_a ( keys %{ $Hash_a } ){  #print "\n20190105-4-0-0-0 \$key_a=$key_a\n"; # warn "\$key_a=$key_a\n"; print "\$key_a=$key_a\n";
			    if (  defined ( $Hash_b->{$key_a} )  ){		
			    	#print "\n20190105-4-0-0 \$Hash_b->{$key_a}=$Hash_b->{$key_a}\n";
			    }
			    else {
			    	$diffrentFoumd=1;
			    	print "\n20190105-4-0-1 \$Hash_a->{$key_a}=$Hash_a->{$key_a} didnot found in \$Hash_b=$Hash_b\n";
			    	last FOREACHMAR1;
			    }
			  }
			  if (  $diffrentFoumd==0 ){
			  	$outReal=1;
			  }
			}
			else{
				$outReal=0;
				my $warnMsg= DirFileHandle::ReturnDumperInform ($Hash_1, "\$Hash_1=$Hash_1");	
			  $warnMsg.= DirFileHandle::ReturnDumperInform ($Hash_2, "\$Hash_2=$Hash_2");			
		  	print "\n20190105-4-0-2  \$size_1=$size_1 \$size_2=$size_2 shold be equal \n \$Hash_1=$Hash_1 \$Hash_2=$Hash_2 should be two HASH!!!\n$addtionMsg \n$warnMsg\n";
			}
		}		
				
		else{
			$outReal=0;
			my $warnMsg= DirFileHandle::ReturnDumperInform ($Hash_1, "\$Hash_1=$Hash_1");	
			$warnMsg.= DirFileHandle::ReturnDumperInform ($Hash_2, "\$Hash_2=$Hash_2");			
			print "\n20190105-4-0-2  \$size_1=$size_1 \$size_2=$size_2 shold be 2 number > 0\n \$Hash_1=$Hash_1 \$Hash_2=$Hash_2 should be two HASH!!!\n$addtionMsg \n$warnMsg\n";
			#DieWork::Just_dieWork( $die_MsgHead."20190105-4-1 $addtionMsg \$size_1=$size_1 \$size_2=$size_2 shold be 2 number > 0\n \$Hash_1=$Hash_1 \$Hash_2=$Hash_2 should be two HASH!!!\n$!\n\n\n".$warnMsg );
		}
		  
		
	}	
	else {
		$outReal=0;
		my $warnMsg2= DirFileHandle::ReturnDumperInform ($Hash_1, "\$Hash_1=$Hash_1");
		$warnMsg2= DirFileHandle::ReturnDumperInform ($Hash_2, "\$Hash_2=$Hash_2");	
		print "\n20190105-4-0-3   \$Hash_1=$Hash_1 \$Hash_2=$Hash_2 should be two HASH!!!!!\n$addtionMsg $warnMsg2\n";
			#DieWork::Just_dieWork( $die_MsgHead."20190105-4-2 $addtionMsg \$Hash_1=$Hash_1 \$Hash_2=$Hash_2 should be two HASH!!!\n$!\n\n\n".$warnMsg2 );
	}
	
	return $outReal;
	
}

#  ArrayHashChange::CheckHashKey_includeRelation ($Hash_1, $Hash_2);   
sub CheckHashKey_includeRelation {  #检查两个hash，第一个的key是不是都被第二个所包含。
	my ($Hash_1, $Hash_2, $addtionIfomr)=@_;
	
	my $warnMsgBody="\nIn package  ArrayHashChange,\tIn subCheckHashKey_exactly_theSAME,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
  
  my $subCallereIfm=DirFileHandle::print_SubCallerInform;
  
  my $addtionMsg=$warnMsgBody.$subCallereIfm;  $addtionMsg=$addtionMsg.$addtionIfomr if (   (  defined ( $addtionIfomr )  ) && ( $addtionIfomr=~m/\S+/ )   ) ;
 
	my $outReal=0;
	
	if (  ( ref($Hash_1) ne 'HASH' ) && ( ref($Hash_2) ne 'HASH' )  ){
		my $warnMsg= DirFileHandle::ReturnDumperInform ($Hash_1, "\$Hash_1=$Hash_1");	
	  $warnMsg.= DirFileHandle::ReturnDumperInform ($Hash_2, "\$Hash_2=$Hash_2");			
		DieWork::Just_dieWork( $die_MsgHead."20190105-4-1 $addtionMsg \$Hash_1=$Hash_1 \$Hash_2=$Hash_2 should be two HASH!!!\n$!\n\n\n".$warnMsg );
	}
	if (  ( ref($Hash_1) eq 'HASH' ) && ( ref($Hash_2) eq 'HASH' )  ){
		
		my $size_1=( keys %{ $Hash_1 } );		my $size_2=( keys %{ $Hash_2 } );
		my $Hash_a;                     		my $Hash_b;
		
		print "\n20190105-4-0-0-0 \$size_1=$size_1 \$size_2=$size_2 !!!\n";
		
		if (  ($size_1 ==0) && ($size_2 ==0)   ){  #print "\n20190105-4-0-0-1 \$size_1=$size_1 \$size_2=$size_2 !!!\n";
			 DieWork::Just_dieWork( $die_MsgHead."20190105-4-1 $addtionMsg \$size_1=$size_1 \$size_2=$size_2 shold be 2 number > 0\n \$Hash_1=$Hash_1 \$Hash_2=$Hash_2 should be two HASH!!!\n$!\n\n\n" );
		}
		
		elsif (  ($size_1 >0) && ($size_2 >0)   ){ #print "\n20190105-4-0-0-2 \$size_1=$size_1 \$size_2=$size_2 !!!\n";
			$Hash_a=$Hash_1;			$Hash_b=$Hash_2;    #print "\n20190105-4-0-0-3 \$size_1=$size_1 \$size_2=$size_2 !!!\n";
			#if ($size_1==$size_2){ 
				my $diffrentFoumd=0;
				FOREACHMAR1: foreach my $key_a ( keys %{ $Hash_a } ){  #print "\n20190105-4-0-0-0 \$key_a=$key_a\n"; # warn "\$key_a=$key_a\n"; print "\$key_a=$key_a\n";
			    if (  defined ( $Hash_b->{$key_a} )  ){		
			    	#print "\n20190105-4-0-0 \$Hash_b->{$key_a}=$Hash_b->{$key_a}\n";
			    }
			    else {
			    	$diffrentFoumd=1;
			    	my $warnMsg= DirFileHandle::ReturnDumperInform ($Hash_1, "\$Hash_1=$Hash_1");	
			      $warnMsg.= DirFileHandle::ReturnDumperInform ($Hash_2, "\$Hash_2=$Hash_2");	
			    	print "\n20190105-4-0-1 \$Hash_a->{$key_a}=$Hash_a->{$key_a} didnot found in \$Hash_b=$Hash_b\n$warnMsg\n\n";
			    	last FOREACHMAR1;
			    }
			  }
			  if (  $diffrentFoumd==0 ){
			  	$outReal=1;
			  }
			#}
			#else{
			#	$outReal=0;
			#	my $warnMsg= DirFileHandle::ReturnDumperInform ($Hash_1, "\$Hash_1=$Hash_1");	
			#  $warnMsg.= DirFileHandle::ReturnDumperInform ($Hash_2, "\$Hash_2=$Hash_2");			
		  #	print "\n20190105-4-0-2  \$size_1=$size_1 \$size_2=$size_2 shold be equal \n \$Hash_1=$Hash_1 \$Hash_2=$Hash_2 should be two HASH!!!\n$addtionMsg \n$warnMsg\n";
			#}
		}		
				
		else{
			$outReal=0;
			my $warnMsg= DirFileHandle::ReturnDumperInform ($Hash_1, "\$Hash_1=$Hash_1");	
			$warnMsg.= DirFileHandle::ReturnDumperInform ($Hash_2, "\$Hash_2=$Hash_2");			
			print "\n20190105-4-0-2  \$size_1=$size_1 \$size_2=$size_2 shold be 2 number > 0\n \$Hash_1=$Hash_1 \$Hash_2=$Hash_2 should be two HASH!!!\n$addtionMsg \n$warnMsg\n";
			#DieWork::Just_dieWork( $die_MsgHead."20190105-4-1 $addtionMsg \$size_1=$size_1 \$size_2=$size_2 shold be 2 number > 0\n \$Hash_1=$Hash_1 \$Hash_2=$Hash_2 should be two HASH!!!\n$!\n\n\n".$warnMsg );
		}
		  
		
	}	
	else {
		$outReal=0;
		my $warnMsg2= DirFileHandle::ReturnDumperInform ($Hash_1, "\$Hash_1=$Hash_1");
		$warnMsg2= DirFileHandle::ReturnDumperInform ($Hash_2, "\$Hash_2=$Hash_2");	
		print "\n20190105-4-0-3   \$Hash_1=$Hash_1 \$Hash_2=$Hash_2 should be two HASH!!!!!\n$addtionMsg $warnMsg2\n";
			#DieWork::Just_dieWork( $die_MsgHead."20190105-4-2 $addtionMsg \$Hash_1=$Hash_1 \$Hash_2=$Hash_2 should be two HASH!!!\n$!\n\n\n".$warnMsg2 );
	}
	
	return $outReal;
	
}

#my $outReal=ArrayHashChange::CheckHashKeyPreat ($Hash_1, $Hash_2);
sub CheckHashKeyPreat{  #检查两个hash中的key是否有重复， check the keies of 2 hash, if the same kay found in two hash, then return 1, else return 0
	my ($Hash_1, $Hash_2)=@_;
	my $outReal=0;
	if (  ( ref($Hash_1) eq 'HASH' ) && ( ref($Hash_2) eq 'HASH' )  ){
		
		my $size_1=( keys %{ $Hash_1 } );		my $size_2=( keys %{ $Hash_2 } );
		my $Hash_a;                     		my $Hash_b;
		
		if ($size_1<=$size_2) {			$Hash_a=$Hash_1;			$Hash_b=$Hash_2;  }
		else                  {			$Hash_a=$Hash_2;			$Hash_b=$Hash_1;  }
			
		FOREACHMAR1: foreach my $key_a ( keys %{ $Hash_a } ){# warn "\$key_a=$key_a\n"; print "\$key_a=$key_a\n";
			if (  defined ( $Hash_b->{$key_a} )  ){
				$outReal=1;
		    last FOREACHMAR1;
			}
		}		  
		
	}	
	else {
		warn "\n\$Hash_1=$Hash_1 is not a HASH! :$!\n\n" unless   ( ref($Hash_1) eq 'HASH' ) ; 
		warn "\n\$Hash_2=$Hash_2 is not a HASH! :$!\n\n" unless   ( ref($Hash_2) eq 'HASH' ) ; 
	}
	
	return $outReal;
	
}


#ArrayHashChange::CheckHashKeyPreat_showRepeat ($Hash_1, $Hash_2);  
sub CheckHashKeyPreat_showRepeat{  #check to find the repeat key in two hashes, and show the repeat fisrt found key.
	my ($Hash_1, $Hash_2)=@_;
	my $outReal=0;
	if (  ( ref($Hash_1) eq 'HASH' ) && ( ref($Hash_2) eq 'HASH' )  ){
		
		my $size_1=( keys %{ $Hash_1 } );		my $size_2=( keys %{ $Hash_2 } );
		my $Hash_a;                     		my $Hash_b;
		
		if ($size_1<=$size_2) {			$Hash_a=$Hash_1;			$Hash_b=$Hash_2;  }
		else                  {			$Hash_a=$Hash_2;			$Hash_b=$Hash_1;  }
			
		FOREACHMAR1: foreach my $key_a ( keys %{ $Hash_a } ){
			if (  defined ( $Hash_b->{$key_a} )  ){
				$outReal=1;
				return "\n\n\$key_a=$key_a is repeated!!!\n\n";
		    last FOREACHMAR1;
			}
		}		  
		
	}	
	else {
		warn "\n\$Hash_1=$Hash_1 is not a HASH! :$!\n\n" unless   ( ref($Hash_1) eq 'HASH' ) ; 
		warn "\n\$Hash_2=$Hash_2 is not a HASH! :$!\n\n" unless   ( ref($Hash_2) eq 'HASH' ) ; 
	}
	

	
}


sub Check4HashForReapeatKey{  # ckeck 4 hashes to find is there repeat key in those hashes.
  my ($H1, $H2, $H3, $H4)=@_;
  my $outRl=0;
  my @HsAr= ($H1, $H2, $H3, $H4);
  FOREACHMK2: foreach my $EcH_a (@HsAr){
    foreach my $EcH_b (@HsAr){
      if ($EcH_a ne $EcH_b){
        if 	( &CheckHashKeyPreat($EcH_a, $EcH_b) ){
        	$outRl=1;
        	last FOREACHMK2;
        }
      }
    } 
  }
  return $outRl;
}

sub Check4HashForReapeatKey_showMSG{  # ckeck 4 hashes to find is there repeat key in those hashes. And show the first found repeat key
  my ($H1, $H2, $H3, $H4)=@_;
  my $outRl=0;
  my @HsAr= ($H1, $H2, $H3, $H4);
  FOREACHMK2: foreach my $EcH_a (@HsAr){
    foreach my $EcH_b (@HsAr){
      if ($EcH_a ne $EcH_b){
        if 	( &CheckHashKeyPreat($EcH_a, $EcH_b) ){
        	$outRl=1; 
        	my $repeatMSG=&CheckHashKeyPreat_showRepeat($EcH_a, $EcH_b);
        	return "\n\n\$H1=$H1, \$H2=$H2, \$H3=$H3, \$H4=$H4\nrepeat found in $EcH_a and $EcH_b\n$repeatMSG\n\n\n";
        	last FOREACHMK2;
        }
      }
    } 
  }
  
}


sub ChangArrayIntoSimpleString{  #just join array into string     # my $outString=ArrayHashChange::ChangArrayIntoSimpleString ($inArrayRef);
	my ($inArrayRef)=@_;
	my $warnMsgHead="\n\n\n		In package ArrayHashChange,		\n In sub ChangArrayIntoSimpleString,\n\n"; 
	my $outString;
	if (ref($inArrayRef) eq 'ARRAY'){
		$outString=join (',', @{$inArrayRef});
	}
	else {
		die "$warnMsgHead\nTHe input \$inArrayRef=$inArrayRef is not a ARRAY ref!!\n\n\n\n";
	}
	return $outString;
}

# my $outString=ArrayHashChange::ChangArrayIntoSimpleString_divToTable ($inArrayRef);
sub ChangArrayIntoSimpleString_divToTable{  #just join array into string  #将数组，变成字符串，数组间用 制表符 分开    
	my ($inArrayRef)=@_;
	my $warnMsgHead="\n\n\n		In package ArrayHashChange,		\n In sub ChangArrayIntoSimpleString,\n\n"; 
	my $outString;
	if (ref($inArrayRef) eq 'ARRAY'){
		$outString=join ("\t", @{$inArrayRef});
	}
	else {
		die "$warnMsgHead\nTHe input \$inArrayRef=$inArrayRef is not a ARRAY ref!!\n\n\n\n";
	}
	return $outString;
}


sub Build_HashArray_from_reversed_HASH{  # my $outHash=ArrayHashChange::Build_HashArray_from_reversed_HASH($inHash, $new_key);
	# 对于一个2d的hash，其第一层的key是唯一的，第二层是可能有重复的，所以 如挑选第二层的一个key作为 新的hash的一层key，则会变成一个 hash ->array的 形式
	
	my ($inHash, $new_key)=@_;
	
	my $warnMsgBody="\nIn package  ArrayHashChange,\tIn Build_HashArray_from_reversed_HASH,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
		
	
	if  (   (  defined (  $new_key  )  ) && ( $new_key=~/\S+/ )   ){		
	}   
	else{
		DieWork::Just_dieWork( $die_MsgHead."\n \$new_key=$new_key should be defined and not NULL $!\n\n\n".$caller_inform ); 	
	}
	
	my $outHash;
	if  (   (  defined (  $inHash  )  ) && (  ( ref ($inHash) ) eq 'HASH' )   ){
	  foreach my $key_1 (   keys (  %{ $inHash }  )   ){
	  	if (   (  defined (  $inHash->{$key_1}  )  ) && (  ( ref $inHash->{$key_1} ) eq 'HASH' )    ){
	  		if (   defined (  $inHash->{$key_1}->{$new_key}  )   ) {
	  		  if ( $inHash->{$key_1}->{$new_key}=~/\S+/ ){
	  		  	
	  		  	if (   defined (  $outHash->{ $inHash->{$key_1}->{$new_key} }  )   ){
    	        push @{  $outHash->{ $inHash->{$key_1}->{$new_key} }  }, $inHash->{$key_1};
    	      }
    	      else {
    	        $outHash->{ $inHash->{$key_1}->{$new_key} }->[0]=$inHash->{$key_1};
    	      }
	  		  	
	  		  }
	  		  else {
	  		  	
	  		  	if (   defined (  $outHash->{ '' }  )   ){
    	        push @{  $outHash->{ $inHash->{$key_1}->{$new_key} }  }, $inHash->{$key_1};
    	      }
    	      else {
    	        $outHash->{ '' }->[0]=$inHash->{$key_1};
    	      }
    	      
	  		  }
	  		}
	  		else {
	  			if (   defined (  $outHash->{ 'undefined' }  )   ){
    	      push @{  $outHash->{ $inHash->{$key_1}->{$new_key} }  }, $inHash->{$key_1};
    	    }
    	    else {
    	      $outHash->{ 'undefined' }->[0]=$inHash->{$key_1};
    	    }
	  		}
	  			
	  	}		
	  		
	  }
	}
	else{
		DieWork::Just_dieWork( $die_MsgHead."\n \$inHash=$inHash should be defined and a HASH ref $!\n\n\n".$caller_inform ); 
	}
		
		
		
	return $outHash;
	
}



# example $inHashFile_subKey_HASH
#                 my $inHashFile_subKey_HASH;
#                 
#                 $inHashFile_subKey_HASH->{$musFileHASH_file    }->{'2m3m_delnum'}=0;
#                 $inHashFile_subKey_HASH->{$itprscFileHASH_file }->{'3m2m_iPsRed'}=1;
#                 $inHashFile_subKey_HASH->{$domPNGFileHASH_file }->{'4m1m_domPng'}=2;
#                 $inHashFile_subKey_HASH->{$nexTreeFileHASH_file}->{'5m1m_MrBays'}=3;


# my $out_HASH=ArrayHashChange::GetBigHash_fromLotsOfSubHASH_FILE( $inHashFile_subKey_HASH ); 
sub GetBigHash_fromLotsOfSubHASH_FILE{  #input the hash ,in which the hashdumpFile as key, the inHash key as 2d key. then extract all of those 2d key=> val, and build up a new big hash
	my ($inHashFile_subKey_HASH)=@_;
	
	my ($warnMsgBody, $warnMsgHead, $die_MsgHead, $caller_inform)=@{ DieWork::BuildWarnDieHeadInform( 'ArrayHashChange', 'GetBigHash_fromLotsOfSubHASH_FILE' ) };

  DieWork::Check_Hash_or_DIE( $inHashFile_subKey_HASH, "\$inHashFile_subKey_HASH", $die_MsgHead, $caller_inform  );
  
  my $out_HASH;
	foreach my $inHashFile (    sort { $a cmp $b } (  keys ( %{ $inHashFile_subKey_HASH } )  )    ){ 
		DieWork::Check_FileDirExist_or_DIE( $inHashFile, "\$inHashFile", $die_MsgHead, $caller_inform  );
		my $inHashAdd=Storable::retrieve ( $inHashFile );
		DieWork::Check_Hash_or_DIE( $inHashAdd, "\$inHashAdd", $die_MsgHead, $caller_inform  );
  
		
		DieWork::Check_Hash_or_DIE( $inHashFile_subKey_HASH->{$inHashFile}, "\$inHashFile_subKey_HASH->{\$inHashFile}=\$inHashFile_subKey_HASH->{$inHashFile}", $die_MsgHead, $caller_inform  );
		
		
		foreach my $inHashKEY (    sort { $a cmp $b } (  keys ( %{ $inHashFile_subKey_HASH->{$inHashFile} } )  )    ){ 
			
			DieWork::Check_DfdNoEmptNUMBER_or_NOT( $inHashFile_subKey_HASH->{$inHashFile}->{$inHashKEY}, "\$inHashFile_subKey_HASH->{\$inHashFile}->{\$inHashKEY}=\$inHashFile_subKey_HASH->{$inHashFile}->{$inHashKEY}", $die_MsgHead, $caller_inform  );
		  my $colNUMber=$inHashFile_subKey_HASH->{$inHashFile}->{$inHashKEY};
		
		  foreach my $fastaPathKey (    sort { $a cmp $b } (  keys ( %{ $inHashAdd } )  )    ){ 
		    $out_HASH->{$fastaPathKey}->{$inHashKEY}=$inHashAdd->{$fastaPathKey}->{$inHashKEY};
		  }  	
		
		}	
		
		
	}
	
	
	
	return $out_HASH;
		
}








1;







#  sub fatchSomething_from_dumpedHashOrArray{
#  	my ($dumpFile, $arry_or_hash, 
#  	                              $lv_0_type, $lv_0_Key, 
#  	                              $lv_1_type, $lv_1_key, 
#  	                              $lv_2_type, $lv_2_Key, 
#  	                              $lv_3_type, $lv_3_key,
#  	                              $lv_4_type, $lv_4_key,
#  	                              $lv_5_type, $lv_5_key,
#  	                              $lv_6_type, $lv_6_key,
#  	                              $lv_7_type, $lv_7_key,
#  	                              $lv_8_type, $lv_8_key,
#  	                              $lv_9_type, $lv_9_key,
#  	                              
#  	                              )=@_;
#  	                              
#  	my $warnMsgBody="\nIn package  ArrayHashChange,\tIn sub fatchSomething_from_dumpedHashOrArray,\n\n";	
#    my $warnMsgHead="\n\n\n$warnMsgBody";
#  	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
#    
#    my $caller_inform=DirFileHandle::print_SubCallerInform;
#    
#    my $outPut;
#    
#    
#    
#    my $fatchLevel=-1;
#    if  (    (   defined ( $lv_0_Key )  ) && (  defined ( $lv_0_type )   ) && (   (  ( $lv_0_type eq 'A' ) && ( $lv_0_Key =~m/^\d+$/ )  ) || (  ( $lv_0_type eq 'H' )  )   )    ){ 
#    	$fatchLevel=0;   	
#    	if  (    (   defined ( $lv_1_Key )  ) && (  defined ( $lv_1_type )   ) && (   (  ( $lv_1_type eq 'A' ) && ( $lv_1_Key =~m/^\d+$/ )  ) || (  ( $lv_1_type eq 'H' )  )   )    ){
#    		$fatchLevel=1;
#    		if  (    (   defined ( $lv_2_Key )  ) && (  defined ( $lv_2_type )   ) && (   (  ( $lv_2_type eq 'A' ) && ( $lv_2_Key =~m/^\d+$/ )  ) || (  ( $lv_2_type eq 'H' )  )   )    ){
#    		$fatchLevel=2;
#    		  if  (    (   defined ( $lv_3_Key )  ) && (  defined ( $lv_3_type )   ) && (   (  ( $lv_3_type eq 'A' ) && ( $lv_3_Key =~m/^\d+$/ )  ) || (  ( $lv_3_type eq 'H' )  )   )    ){
#    		  $fatchLevel=3;
#    		    if  (    (   defined ( $lv_3_Key )  ) && (  defined ( $lv_3_type )   ) && (   (  ( $lv_3_type eq 'A' ) && ( $lv_3_Key =~m/^\d+$/ )  ) || (  ( $lv_3_type eq 'H' )  )   )    ){
#    		    $fatchLevel=4;
#    		      if  (    (   defined ( $lv_4_Key )  ) && (  defined ( $lv_4_type )   ) && (   (  ( $lv_4_type eq 'A' ) && ( $lv_4_Key =~m/^\d+$/ )  ) || (  ( $lv_4_type eq 'H' )  )   )    ){
#    		      $fatchLevel=5;
#    		        if  (    (   defined ( $lv_5_Key )  ) && (  defined ( $lv_5_type )   ) && (   (  ( $lv_5_type eq 'A' ) && ( $lv_5_Key =~m/^\d+$/ )  ) || (  ( $lv_5_type eq 'H' )  )   )    ){
#    		        $fatchLevel=5;
#    		          if  (    (   defined ( $lv_6_Key )  ) && (  defined ( $lv_6_type )   ) && (   (  ( $lv_6_type eq 'A' ) && ( $lv_6_Key =~m/^\d+$/ )  ) || (  ( $lv_6_type eq 'H' )  )   )    ){
#    		          $fatchLevel=6;
#    		            if  (    (   defined ( $lv_7_Key )  ) && (  defined ( $lv_7_type )   ) && (   (  ( $lv_7_type eq 'A' ) && ( $lv_7_Key =~m/^\d+$/ )  ) || (  ( $lv_7_type eq 'H' )  )   )    ){
#    		            $fatchLevel=7;
#    		              if  (    (   defined ( $lv_8_Key )  ) && (  defined ( $lv_8_type )   ) && (   (  ( $lv_8_type eq 'A' ) && ( $lv_8_Key =~m/^\d+$/ )  ) || (  ( $lv_8_type eq 'H' )  )   )    ){
#    		              $fatchLevel=8;
#    		                if  (    (   defined ( $lv_9_Key )  ) && (  defined ( $lv_9_type )   ) && (   (  ( $lv_9_type eq 'A' ) && ( $lv_9_Key =~m/^\d+$/ )  ) || (  ( $lv_9_type eq 'H' )  )   )    ){
#    		                $fatchLevel=9;
#    		                warn "There are 10 level of hash or array you want to fatch!!!!!";
#    		                }
#    		              }
#    		            }
#    		          }
#    		        }
#    		      }
#    		    }
#    		  }
#    		}
#    	}
#    }
#    
#    
#    
#    
#    my $lastLevel=-1;
#    
#    my $outLevelHASH;
#    
#    if ( -e ( $dumpFile ) ){
#    	if (   (  defined ( $arry_or_hash )  ) && (  ( $arry_or_hash eq 'ARRAY' ) || ( $arry_or_hash eq 'HASH' )  )   ){
#    		my $dump_HASH=Storable::retrieve( $dumpFile );
#    		
#    		if  (    (  defined ( $dump_HASH )  ) && (   (  ref ( $dump_HASH )  ) eq $arry_or_hash   )   ){
#    		 	
#    		 	if  (   (  defined ( $lv_0_type )  ) && ( $lv_0_type eq 'A' ) && (  defined ( $lv_0_Key )  ) && ( $lv_0_Key =~m/^\d+$/ )    ){
#    		 		
#    		 		if (   (  defined ( $dump_HASH->[$lv_0_type] )  ){
#    		 			$lastLevel=0;
#    		 			$outLevelHASH->{$lastLevel}=$dump_HASH->[$lv_0_type];
#    		 			
#    		 			
#    		 			
#    		 		}
#    		 		else {
#    		 			DieWork::Just_dieWork( $die_MsgHead."\$dump_HASH->[$lv_0_type]=$dump_HASH->[$lv_0_type] should be defined!!! : $! \n\n\n\n\n" );
#    		 		}
#    		 		
#    		 	}
#    		 	elsif (   (  defined ( $lv_0_type )  ) && ( $lv_0_type eq 'H' ) && (  defined ( $lv_0_Key )  )   ){
#    		 		if (   (  defined ( $dump_HASH->{$lv_0_type} )  ){
#    		 			$lastLevel=0;
#    		 			$outLevelHASH->{$lastLevel}=$dump_HASH->{$lv_0_type};
#    		 		}
#    		 		else {
#    		 			DieWork::Just_dieWork( $die_MsgHead."\$dump_HASH->[$lv_0_type]=$dump_HASH->[$lv_0_type] should be defined!!! : $! \n\n\n\n\n" );
#    		 		}
#    		 	}
#    		 	
#    		}
#    		else {
#    			DieWork::Just_dieWork( $die_MsgHead."\$dump_HASH=$dump_HASH should be \$arry_or_hash=$arry_or_hash !!! : $! \n\n\n\n\n" );
#    		}
#    	}
#    	else {
#    		DieWork::Just_dieWork( $die_MsgHead."\$arry_or_hash=$arry_or_hash should be ARRAY or HASH!!! : $! \n\n\n\n\n" );
#    	}
#    } 
#    else {
#    	DieWork::Just_dieWork( $die_MsgHead."\$dumpFile=$dumpFile is not exist or cannot be read!!! : $! \n\n\n\n\n" );
#    }                             
#  	                              
#  	                              
#  }