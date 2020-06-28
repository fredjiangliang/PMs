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


package DieWork;

# my ($warnMsgBody, $warnMsgHead, $die_MsgHead, $caller_inform)=@{ DieWork::BuildWarnDieHeadInform( $pm_name, $subMethod_name ) };
sub BuildWarnDieHeadInform{
	my ($pm_name, $subMethod_name)=@_;
	
	#DieWork::Check_DfdNoEmptString_or_DIE( $pm_name, "\$pm_name", $die_MsgHead, $caller_inform  );
	#DieWork::Check_DfdNoEmptString_or_DIE( $inVal, $inVal_look, $die_MsgHead, $caller_inform  );
	
	my $warnMsgBody="\nIn package  $pm_name,\tIn sub $subMethod_name,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	return [$warnMsgBody, $warnMsgHead, $die_MsgHead, $caller_inform];
	
}

sub Just_dieWork{  # DieWork::Just_dieWork( $dieWord );
	my ( $dieWord )=@_;
	
	print $dieWord;
	die   $dieWord;
}


sub Print_and_warn{  # DieWork::Print_and_warn( $WarnMsg );
	my ( $WarnMsg )=@_;
	
	print $WarnMsg."\n";
	warn   "\n".'nomarl warn:'.$WarnMsg."\n";
}

# DieWork::Check_Hash_or_DIE( $inVal, $inVal_look, $die_MsgHead, $caller_inform  );
sub Check_Hash_or_DIE{  # 检查相应的量是不是 HASH ref ，否则 DIE
	my ( $inVal, $inVal_look, $die_MsgHead, $caller_inform )=@_;
	
	my $hashOrNot=0;
	$hashOrNot=DieWork::Check_Hash_or_NOT ( $inVal );
	if (   $hashOrNot == 0   ){ 	
		DieWork::Just_dieWork( $die_MsgHead."\n $inVal_look=$inVal should be a HASH ref !!  $!\n\n\n".$caller_inform );	
	}	
}

# DieWork::Check_Hash_or_NOT( $inVal );
sub Check_Hash_or_NOT{  # 检查相应的量是不是 HASH ref，是返回1，否则返回0
	my ( $inVal )=@_;
  if (   (  defined ( $inVal )  ) && (  ref ( $inVal ) eq 'HASH'  )   ){ return 1;	}
  else{ return 0; }
}

# DieWork::Check_Array_or_DIE( $inVal, $inVal_look, $die_MsgHead, $caller_inform  );
sub Check_Array_or_DIE{  # 检查相应的量是不是 Array ref ，否则 DIE
	my ( $inVal, $inVal_look, $die_MsgHead, $caller_inform )=@_;
	
	my $ArrayOrNot=0;
	$ArrayOrNot=DieWork::Check_Array_or_NOT ( $inVal );
	if (   $ArrayOrNot == 0   ){ 	
		DieWork::Just_dieWork( $die_MsgHead."\n $inVal_look=$inVal should be a Array ref !!  $!\n\n\n".$caller_inform );	
	}	
}

# DieWork::Check_Array_or_NOT( $inVal );
sub Check_Array_or_NOT{  # 检查相应的量是不是 Array ref，是返回1，否则返回0
	my ( $inVal )=@_;
  if (   (  defined ( $inVal )  ) && (  ref ( $inVal ) eq 'ARRAY'  )   ){ return 1;	}
  else{ return 0; }
}


## DieWork::Check_FileDirExist_or_DIE( $inVal, $inVal_look, $die_MsgHead, $caller_inform  );
sub Check_FileDirExist_or_DIE{ # 检查相应的量是不是 已存在的 文件或文件夹 ，否则 DIE
	my ( $inVal, $inVal_look, $die_MsgHead, $caller_inform )=@_;
	
	my $yes_or_not=0;
	$yes_or_not=DieWork::Check_FileDirExist_or_NOT ( $inVal );
	if   (    $yes_or_not == 0    ){
		DieWork::Just_dieWork( $die_MsgHead."\n $inVal_look=$inVal should be a defined file or dir path string  !! !!  $!\n\n\n".$caller_inform );	
	}
	
}

# DieWork::Check_FileDirExist_or_NOT( $inVal );
sub Check_FileDirExist_or_NOT{  # 检查相应的量是不是 已存在的 文件或文件夹 ，是返回1，否则返回0
	my ( $inVal )=@_;
  if (   (  defined ( $inVal )  ) && ( $inVal=~m/\S+/ ) && (  ( -e $inVal )  )   ){ return 1;	}
  else{ return 0; }
}





## DieWork::Check_DfdNoEmptString_or_DIE( $inVal, $inVal_look, $die_MsgHead, $caller_inform  );
sub Check_DfdNoEmptString_or_DIE{ # 检查相应的量是不是 非空 ，否则 DIE
	my ( $inVal, $inVal_look, $die_MsgHead, $caller_inform )=@_;
	
	my $yes_or_not=0;
	$yes_or_not=DieWork::Check_DfdNoEmptString_or_NOT ( $inVal );
	if   (    $yes_or_not == 0    ){		
		DieWork::Just_dieWork( $die_MsgHead."\n $inVal_look=$inVal should be a defined not empty string  !! !!  $!\n\n\n".$caller_inform );	
	}
	
}
# DieWork::Check_DfdNoEmptString_or_NOT( $inVal );
sub Check_DfdNoEmptString_or_NOT{  # 检查相应的量是不是 非空 ，是返回1，否则返回0
	my ( $inVal )=@_;
  if (   (  defined ( $inVal )  ) && ( $inVal=~m/\S+/ )   ){ return 1;	}
  else{ return 0; }
}

## DieWork::Check_DfdNoEmptNUMBER_or_DIE( $inVal, $inVal_look, $die_MsgHead, $caller_inform  );
sub Check_DfdNoEmptNUMBER_or_DIE{ # 检查相应的量是不是 整数字 ，否则 DIE
	my ( $inVal, $inVal_look, $die_MsgHead, $caller_inform )=@_;
	
	my $yes_or_not=0;
	$yes_or_not=DieWork::Check_DfdNoEmptNUMBER_or_NOT ( $inVal );
	if   (    $yes_or_not == 0    ){		
		DieWork::Just_dieWork( $die_MsgHead."\n $inVal_look=$inVal should be a defined number  !! !!  $!\n\n\n".$caller_inform );	
	}
	
}


# DieWork::Check_DfdNoEmptNUMBER_or_NOT( $inVal );
sub Check_DfdNoEmptNUMBER_or_NOT{  # 检查相应的量是不是 整数字 ，是返回1，否则返回0
	my ( $inVal )=@_;
  if (   (  defined ( $inVal )  ) && ( $inVal=~m/\d+/ )   ){ return 1;	}
  else{ return 0; }
}


## DieWork::Check_INTNB_equal_to_aNUMBER_or_DIE( $inVal, $inNub, $inVal_look, $die_MsgHead, $caller_inform  );
sub Check_INTNB_equal_to_aNUMBER_or_DIE{ # 检查相应的量是不是 非空 ，否则 DIE
	my ( $inVal,  $inNub, $inVal_look, $die_MsgHead, $caller_inform )=@_;
	
	my $yes_or_not=0;
	$yes_or_not=DieWork::Check_INTNB_equal_to_aNUMBER_or_NOT ( $inVal, $inNub  );
	if   (    $yes_or_not == 0    ){		
		DieWork::Just_dieWork( $die_MsgHead."\n $inVal_look=$inVal != \$inNub=$inNub, These 2 should be a equal  !! !!  $!\n\n\n".$caller_inform );	
	}
	
}


# DieWork::Check_INTNB_equal_to_aNUMBER_or_NOT( $inVal, $inNub );
sub Check_INTNB_equal_to_aNUMBER_or_NOT{  # 检查相应的量是不是 一个包括0的确定的整数 ，是返回1，否则返回0
	my ( $inVal, $inNub )=@_;
  if (   (  defined ( $inVal )  ) && ( $inVal=~m/^\d+$/ ) && ( $inVal == $inNub )   ){ return 1;	}
  else{ return 0; }
}



## DieWork::Check_INTNB_equal_to_aString_or_DIE( $inVal, $inStr, $inVal_look, $die_MsgHead, $caller_inform  );
sub Check_INTNB_equal_to_aString_or_DIE{ # 检查相应的量是不是 非空 ，否则 DIE
	my ( $inVal,  $inStr, $inVal_look, $die_MsgHead, $caller_inform )=@_;
	
	my $yes_or_not=0;
	$yes_or_not=DieWork::Check_INTNB_equal_to_aString_or_NOT ( $inVal, $inStr  );
	if   (    $yes_or_not == 0    ){		
		DieWork::Just_dieWork( $die_MsgHead."\n $inVal_look=$inVal != \$inStr=$inStr, These 2 should be a equal  !! !!  $!\n\n\n".$caller_inform );	
	}	
}

# DieWork::Check_INTNB_equal_to_aString_or_NOT( $inVal, $inStr );
sub Check_INTNB_equal_to_aString_or_NOT{  # 检查相应的量是不是 一个和 $inStr 相等的 字符串
	my ( $inVal, $inStr )=@_;
  if (   (  defined ( $inVal )  ) && ( $inVal=~m/^\S+$/ ) && ( $inVal eq $inStr )   ){ return 1;	}
  else{ return 0; }
}










1;