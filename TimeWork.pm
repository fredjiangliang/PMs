
#!/usr/bin/perl -w

BEGIN{  push (@INC, '/home/fredjiang/fredProgs/PMs')  };

use strict;
use warnings;
use Time::HiRes qw(sleep time);
use DirFileHandle;

package TimeWork;


#sub GetTimeDirOrFileName{
                           #输入: 1，无
                           #输出：1，当前时间 及进程名 等信息，所组成的字符串，该字符串可用来生成 临时文件或临时文件夹




sub GetNowTime{  # my $tmNow=TimeWork::GetNowTime;   #获得此时此刻的时间，没有微秒数据
	my $tmNow=time;
	return $tmNow;
}

sub GetHumanRead_Time{  # my $HmReadTimeNow=TimeWork::GetHumanRead_Time ($timeIn);   #获得此时此刻的时间，没有微秒数据
	my ($timeIn)=@_;
	my ($sec, $min, $hou, $mdy, $mon, $yer, $wdy, $ydy, $isd) = gmtime( $timeIn );  #print $sec,"\t", $min,"\t", $hou,"\t", $mdy,"\t", $mon,"\t", $yer,"\t", $wdy,"\t", $ydy,"\t", $isd, "\n";
  $yer=$yer+1900;  $mon=sprintf "%02d", $mon; $mdy=sprintf "%02d", $mdy; $hou=sprintf "%02d", $hou; $min=sprintf "%02d", $min; $sec=sprintf "%02d", $sec; 
  my $HmReadTimeNow="$yer-$mon-${mdy}-$hou:$min:${sec}";
  return $HmReadTimeNow;
}

sub GetHumanRead_NOW_Time{  #my $HmReadTimeNow=TimeWork::GetHumanRead_NOW_Time();
  my $tmNow        =TimeWork::GetNowTime;
  my $HmReadTimeNow=TimeWork::GetHumanRead_Time( $tmNow );
  return $HmReadTimeNow;
}




sub GetNowTimePid{
	
	my $tmNow =TimeWork::GetNowTime;	  
  my $tpDir =TimeWork::GetTimePid( $tmNow );  
  return $tpDir;
}


sub GetTimeDirOrFileName{			    
  my $tpDir =TimeWork::GetNowTimePid;
  my $AdNb=0;  while (  ( -d ($tpDir) ) || ( -e ($tpDir) )  ){   $tpDir="${tpDir}A$AdNb";   $AdNb++; }
  return $tpDir;
}

sub GetTimePid{	
	my ($tmNow)=@_;
	
  my ($sec, $min, $hou, $mdy, $mon, $yer, $wdy, $ydy, $isd) = gmtime( $tmNow );  #print $sec,"\t", $min,"\t", $hou,"\t", $mdy,"\t", $mon,"\t", $yer,"\t", $wdy,"\t", $ydy,"\t", $isd, "\n";
  $yer=$yer+1900;  $mon=sprintf "%02d", $mon; $mdy=sprintf "%02d", $mdy; $hou=sprintf "%02d", $hou; $min=sprintf "%02d", $min; $sec=sprintf "%02d", $sec; 
  my $pid=$$;
  my $tpDir="D$yer$mon${mdy}T$hou$min${sec}P$pid";
  
  #print "\n\$tpDir=$tpDir\n\n";
  return $tpDir;
}



sub GetNowTimeMicroSeconds{   #获得此时此刻的时间，有微秒数据
	my $tmNow=Time::HiRes::time;
	return $tmNow;
}

sub GetNowTimePid_microSecond{   # my $NtPid=TimeWork::GetNowTimePid_microSecond ();   #获得此时此刻的时间，有微秒数据
  my $tmNow=TimeWork::GetNowTimeMicroSeconds;  
  my $NtPid=TimeWork::GetTimePid_microSecond( $tmNow );
  return $NtPid;
}

sub GetTimePid_microSecond{	
  my ($tmNow)=@_;
  
  my ($sec, $min, $hou, $mdy, $mon, $yer, $wdy, $ydy, $isd) = gmtime( $tmNow );  #print $sec,"\t", $min,"\t", $hou,"\t", $mdy,"\t", $mon,"\t", $yer,"\t", $wdy,"\t", $ydy,"\t", $isd, "\n";
  $yer=$yer+1900;  $mon=sprintf "%02d", $mon; $mdy=sprintf "%02d", $mdy; $hou=sprintf "%02d", $hou; $min=sprintf "%02d", $min; 
  
               #print "\$tmNow=$tmNow\n"; 
  my $timeInt=int $tmNow;                  #print "\$timeInt=$timeInt\n"; 
  my $flotTime=$tmNow-$timeInt;            #print "\$flotTime=$flotTime\n"; 
  $flotTime=sprintf "%06f", $flotTime;     #print "2 \$flotTime=$flotTime\n";  
  $flotTime=~s/^0\./_/;                    #print "3 \$flotTime=$flotTime\n"; 
  
  $sec=sprintf "%02d", $sec; $sec.=$flotTime;
  my $pid=$$;
  my $tpDir="D$yer$mon${mdy}T$hou$min${sec}P$pid";
  #my $AdNb=0;  while (  ( -d ($tpDir) ) || ( -e ($tpDir) )  ){   $tpDir="${tpDir}A$AdNb";   $AdNb++; }
  #print "\n\$tpDir=$tpDir\n\n";
  return $tpDir;
}

sub GetLastTime{  #计算两个时间(格式：2012-11-15 09:35:59)之间的间隔时间长度，并用可读的形式输出, $time1 >  $time2
  my ($time1, $time2)=@_;  #print "\$time1=$time1\n\$time2=$time2\n\n";
  my @time1Ar; my @time2Ar;  # 2012-11-15 09:35:59
  if ($time1=~m/^(\d+)-(\d+)-(\d+)-(\d+):(\d+):(\d+)$/){
  	@time1Ar=($1, $2, $3, $4, $5, $6); #print "\@time1Ar=@time1Ar\n";
  }
  if ($time2=~m/^(\d+)-(\d+)-(\d+)-(\d+):(\d+):(\d+)$/){
  	@time2Ar=($1, $2, $3, $4, $5, $6); #print "\@time2Ar=@time2Ar\n";
  }
  my $year=   ($time1Ar[0]-$time2Ar[0]);
  my $month=  ($time1Ar[1]-$time2Ar[1]);
  my $day=    ($time1Ar[2]-$time2Ar[2]);
  my $hour=   ($time1Ar[3]-$time2Ar[3]);
  my $minute= ($time1Ar[4]-$time2Ar[4]);
  my $second= ($time1Ar[5]-$time2Ar[5]);
  
  my $lastTime=( (60*60*24*365*     ($year   ))
                +(60*60*24*30.4375* ($month  ))
                +(60*60*24*         ($day    ))
                +(60*60*            ($hour   ))
                +(60*               ($minute ))
                +(1*                ($second ))
               );
  return $lastTime;
}

sub HumanRead_Last_Time{  #把一个以秒为单位的数字，转化为年月日
  my ($seconds)=@_;
  my $number=0; my $unit="Second";  
  if    ($seconds>=(60*60*24*365))     {$number=$seconds/(60*60*24*365);     $unit="Year";  }
  elsif ($seconds>=(60*60*24*30.4375) ){$number=$seconds/(60*60*24*30.4375); $unit="Month"; }
  elsif ($seconds>=(60*60*24) )        {$number=$seconds/(60*60*24);         $unit="Day"; }
  elsif ($seconds>=(60*60) )           {$number=$seconds/(60*60);            $unit="Hour"; }
  elsif ($seconds>=(60) )              {$number=$seconds/(60);               $unit="Minute"; }
  elsif ($seconds>(0) )                {$number=$seconds;                    $unit="Second"; }
  else                                 {$number=$seconds;                    $unit="Second";  }
  
  my $outPut=sprintf "%4.2f\t%8s",$number,$unit; 
}

sub HumanReadTime{  #把一个以秒为单位的数字，转化为年月日
  my ($seconds)=@_;
  my $number=0; my $unit="Second";  
  if    ($seconds>=(60*60*24*365))     {$number=$seconds/(60*60*24*365);     $unit="Year";  }
  elsif ($seconds>=(60*60*24*30.4375) ){$number=$seconds/(60*60*24*30.4375); $unit="Month"; }
  elsif ($seconds>=(60*60*24) )        {$number=$seconds/(60*60*24);         $unit="Day"; }
  elsif ($seconds>=(60*60) )           {$number=$seconds/(60*60);            $unit="Hour"; }
  elsif ($seconds>=(60) )              {$number=$seconds/(60);               $unit="Minute"; }
  elsif ($seconds>(0) )                {$number=$seconds;                    $unit="Second"; }
  else                                 {$number=$seconds;                    $unit="Second";  }
  
  my $outPut=sprintf "%4.2f\t%8s",$number,$unit; 
}

sub ChangeHumanReadTime_into_Seconds{  #   my $seconds=TimeWork::ChangeHumanReadTime_into_Seconds($inString);  # 将 123y这种形式的时间变成秒
	my ($inString)=@_;
	
	my $warnMsgBody="\nIn package  TimeWork,\tIn sub ChangeHumanReadTime_into_Seconds,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	my $seconds;
	
	$inString=~s/^\s+//;
	$inString=~s/\s+$//;
	$inString=~s/\s+//g;
	$inString=lc $inString;
	if ( $inString=~m/(\d+)(\w)/ ){
		my $timeNb=$1;
		my $timeUt=$2;
		
		if    ( $timeUt eq 'y' )  { $seconds=$timeNb*(60*60*24*365    ); }
		elsif ( $timeUt eq 'm' )  { $seconds=$timeNb*(60*60*24*30.4375); }
		elsif ( $timeUt eq 'd' )  { $seconds=$timeNb*(60*60*24        ); }
		elsif ( $timeUt eq 'h' )  { $seconds=$timeNb*(60*60           ); }
		elsif ( $timeUt eq 'f' )  { $seconds=$timeNb*(60              ); }
		elsif ( $timeUt eq 's' )  { $seconds=$timeNb                   ; }
		else {
			DieWork::Just_dieWork( $die_MsgHead."\n \$inString=$inString should be 123s 123S 123f 123F 123m 123M 123d 123D 123y 123Y, 123 here could be any number!!\ns=second f=fenzhong(minute) h=hour d=day m=monty y=year\n  $!\n\n\n".$caller_inform ); 
		}
		
	}
	else {
		DieWork::Just_dieWork( $die_MsgHead."\n \$inString=$inString should be 123s 123S 123f 123F 123m 123M 123d 123D 123y 123Y, 123 here could be any number!!\ns=second f=fenzhong(minute) h=hour d=day m=monty y=year\n  $!\n\n\n".$caller_inform ); 
	}
	
	return $seconds;
}

1;