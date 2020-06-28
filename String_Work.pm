#!/usr/bin/perl -w
BEGIN{  push (@INC, '/home/fredjiang/fredProgs/PMs')  };

use strict;
use warnings;

use DirFileHandle;

package String_Work;

sub Chang_String_to_array{  #将字符串转化为一个字符为一个元素的 数组
	my ($inString)=@_;
	
	my $warnMsgBody="\nIn package  String_Work,\tIn sub Chang_String_to_array,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	if   (   (  defined ( $inString )  ) && ( $inString=~m/\S+/ )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$inString=$inString should be a defined not empty string  !!  $!\n\n\n".$caller_inform ); 	}
	
	my @strArray=split ('',$inString);
	return [ @strArray ];
}


sub Chang_String_into_subStr_array{  #my $subStrAraay=String_Work::Chang_String_into_subStr_array($inString, $subLength);
	my ($inString, $subLength)=@_;
	
	my ($warnMsgBody, $warnMsgHead, $die_MsgHead, $caller_inform)=@{ DieWork::BuildWarnDieHeadInform( 'String_Work', 'Chang_String_into_subStr_array' ) };

	
	DieWork::Check_DfdNoEmptString_or_DIE( $inString,  "\$inString",  $die_MsgHead, $caller_inform  );
	DieWork::Check_DfdNoEmptNUMBER_or_DIE( $subLength, "\$subLength", $die_MsgHead, $caller_inform  );
	
	my $stringLength=length( $inString );
	my $remainder=$stringLength % $subLength;
	if ( $remainder == 0 ){			}
	else {
		DieWork::Just_dieWork( $die_MsgHead."\n \$remainder=$remainder=\$stringLength=$stringLength % \$subLength=$subLength \n  should be 0  !! !!  $!\n\n\n".$caller_inform );	
	}
	
	
	my @strArray;
	
	my $stNB=0; 
	while (  $stNB<$stringLength ){
		my $subString=substr ( $inString, $stNB, $subLength);
		push @strArray, $subString;
		$stNB+=$subLength;		
	}
	
	
	return [ @strArray ];
}


sub GetSubString_from_middlePoint{  #my ( $outString, $NewPos_or_org_middlePostion_in_subStr) = @{ String_Work::GetSubString_from_middlePoint($middle_postion, $flanking_length, $in_string) };
	my ($middle_postion, $flanking_length, $in_string)=@_;
	
	my $warnMsgBody="\nIn package  String_Work,\tIn sub GetSubString_from_middlePoint,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	if   (   (  defined ( $middle_postion )  ) && ( $middle_postion=~m/\d+/ )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$middle_postion=$middle_postion should be a number  !!  $!\n\n\n".$caller_inform ); 	}
	
	if   (   (  defined ( $flanking_length )  ) && ( $flanking_length=~m/\d+/ ) && ( $flanking_length > 0 )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$flanking_length=$flanking_length should be a a number > 0 !!  $!\n\n\n".$caller_inform ); 	}
	
	
	if   (   (  defined ( $in_string )  ) && ( $in_string=~m/\S+/ )   ){}
	else{		DieWork::Just_dieWork( $die_MsgHead."\n \$in_string=$in_string should be a defined not empty string  !!  $!\n\n\n".$caller_inform ); 	}
	
	my $string_length=length $in_string;
	
	my $realStartPos=$middle_postion-$flanking_length; 
	if ( $realStartPos < 0 ) { $realStartPos = 0 ; }
	
	my $realEnd__Pos=$middle_postion+$flanking_length; 
	if (  $realStartPos > ( $string_length - 1 )  ) { $realStartPos = ( $string_length - 1 ) ; }
	
	my $sgment_length=$realEnd__Pos-$realStartPos+1;
	
	my $outString=substr ($in_string, $realStartPos,  $sgment_length);
	
	my $NewPos_or_org_middlePostion_in_subStr=$middle_postion-$realStartPos;
	return [ $outString, $NewPos_or_org_middlePostion_in_subStr ];
	
}


1;