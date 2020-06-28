#!/usr/bin/perl -w
BEGIN{  push (@INC, '/home/fredjiang/fredProgs/PMs')  };


use strict;
use warnings;


use InFileHandle;
use TimeWork;

############################################################################################################
#

#         
############################################################################################################



                      
package OnlyOnePidWork;




#######################################################################################################################################################################################################################################################################################################################################
#
#  ʹ�÷�����
#
#    use OnlyOnePidWork;
#
#    my $flgFile="/home/fredjiang/xxxx/yyy/zzz/ddd.txt";                  #���������� flag���ļ�����Ҫ�༭����������ļ������⣬����ֱ��ɾ��
#    my $lastCheckTime='1h';                                              #��һ��ֻ�ܵ������еĴ����������е�ʱ�������������ȴ���ã� �����ʽ �� 123y, 123m ,123d, 123h, 123f, 123s �� 123
#
#    while/for/foreach{                                                   #���������ѭ�� 
#      my $nowWorkTimePidInformMark='';                                     #��ʼ��һ�� ��������װ�������е�ʱ�� ��timePidMark�Ĵ��룬�����ݸ����յ� ����״̬��д ����
#	     $nowWorkTimePidInformMark=OnlyOnePidWork::CheckAndStart_step_forOnlyOnePidWork($flgFile, $lastCheckTime);  # ��飬�������������е�ָʾ�źţ��� ������Ч�� timePidMark�Ĵ���
#      if (   ( defined ( $nowWorkTimePidInformMark )  ) && ( $nowWorkTimePidInformMark=~/\S+/ )   ){ #һ�� timePidMark�Ĵ��� ������ ����� ��������Ҳ���� ֻ�ܽ���һ�����е� ���� ���� �� ����
#      
#        code code code code code code code code ................... code; # һ��ֻ������һ�εĴ��� ����
#
#        OnlyOnePidWork::CheckAnd_Done_step_forOnlyOnePidWork( $flgFile, $nowWorkTimePidInformMark ); #���ϳ�����ɺ� ���� Flag�ļ���д�� Done��ʱ�䣬�Ա���������ʹ��
#      }
#    }
#
#
#
#
#######################################################################################################################################################################################################################################################################################################################################




sub CheckAndStart_step_forOnlyOnePidWork{  # my $nowWorkTimePidInformMark=OnlyOnePidWork::CheckAndStart_step_forOnlyOnePidWork( $FlagFile, $longestTime_set_to_Wait );  #���� FlagFile��·�� �� �������̵���ĵȴ�ʱ��
	my ($FlagFile, $longestTime_set_to_Wait)=@_;
	
	my $warnMsgBody="\nIn package  OnlyOnePidWork,\tIn sub CheckAndStart_step_forOnlyOnePidWork,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	                                                                                                                                       #warn "     201903281824-0-0-0-0 \$longestTime_set_to_Wait=$longestTime_set_to_Wait\n";
	if (   (  defined ( $longestTime_set_to_Wait )  )  && ( $longestTime_set_to_Wait=~m/\S+/ )   ){
		$longestTime_set_to_Wait=~s/\s+//g;                                                                                                  #warn "     201903281824-0-0-0-0 \$longestTime_set_to_Wait=$longestTime_set_to_Wait\n";
		if ( $longestTime_set_to_Wait=~m/^\d+$/ ){                                                                                           #warn "     201903281824-0-0-0-1 \$longestTime_set_to_Wait=$longestTime_set_to_Wait\n";
			$longestTime_set_to_Wait.='s';                                                                                                     #warn "     201903281824-0-0-0-2 \$longestTime_set_to_Wait=$longestTime_set_to_Wait\n";
		}
	}
	else {
		$longestTime_set_to_Wait='5h';                                                                                                       #warn "     201903281824-0-0-0-3 \$longestTime_set_to_Wait=$longestTime_set_to_Wait\n";
	}
	$longestTime_set_to_Wait=TimeWork::ChangeHumanReadTime_into_Seconds($longestTime_set_to_Wait);
	                                                                                                                                       #warn "     201903281824-0-0-0-4 \$longestTime_set_to_Wait=$longestTime_set_to_Wait\n";
	
	                                              my $nowTimeMicroSecond=TimeWork::GetNowTimeMicroSeconds; 
	my $TimePidMark      ='TimePidMark';          my $nowWorkTimePidInformMark=TimeWork::GetTimePid_microSecond( $nowTimeMicroSecond );
	my $ProcessID        ='ProcessID';            my $pid                     =$$;
	my $WorkIng_StartTime='WorkIng_StartTime';    my $HmRead_StartTime        =TimeWork::GetHumanRead_Time     ( $nowTimeMicroSecond );    #warn "     201903281824-0-0-1 \$pid=$pid \$HmRead_StartTime=$HmRead_StartTime\n"; sleep (2);
	my $WorkIng_Done_Time='WorkIng_Done_Time';  
	
	my $start_new_work=0;
	
	my $FlagFile_noDefined=OnlyOnePidWork::CheckTheFlagFile_noDefined($FlagFile);            #warn "               201903290915-0-0-0 \$FlagFile_noDefined=$FlagFile_noDefined\n";
	if ( $FlagFile_noDefined==1 ){
		$start_new_work=1;
		
		OnlyOnePidWork::Write_Psudo_WorkingDone__intoFlagFile($FlagFile, $nowWorkTimePidInformMark, $pid, $HmRead_StartTime);  # д��ٵ� �����Ϣ����������� �����ļ�����״̬����
		
	}
	else{
		
		#$nowWorkTimePidInformMark=OnlyOnePidWork::Digui_checkAndWright ($FlagFile, $pid, $HmRead_StartTime, $longestTime_set_to_Wait, 0, 5);
		#return $nowWorkTimePidInformMark;
		
		#######
		my $couldBeStart=OnlyOnePidWork::Waiting_to_start_new_work ($FlagFile, $pid, $HmRead_StartTime, $longestTime_set_to_Wait); 
		if ( $couldBeStart==1 ){
			$start_new_work=1;
	  }
	  #######
	}
	
	#######
	if   ( $start_new_work==1 ){
		
		$nowWorkTimePidInformMark=OnlyOnePidWork::Digui_checkAndWright ($FlagFile, $pid, $HmRead_StartTime, $longestTime_set_to_Wait, 0, 5);
		return $nowWorkTimePidInformMark;
		#$nowTimeMicroSecond      =TimeWork::GetNowTimeMicroSeconds; 		                        		                     #����ʱ��
		#$nowWorkTimePidInformMark=TimeWork::GetTimePid_microSecond    ( $nowTimeMicroSecond ); 		                       #����ʱ��
		#$HmRead_StartTime        =TimeWork::GetHumanRead_Time         ( $nowTimeMicroSecond) ; 		                       #����ʱ��
		#
		#
		#
		#
		#print "\n20190328-0-3-time:$nowWorkTimePidInformMark      \$nowTimeMicroSecond=$nowTimeMicroSecond     \$HmRead_StartTime=$HmRead_StartTime\n"    ;
		#OnlyOnePidWork::Write_WorkingStart_intoFlagFile($FlagFile, $nowWorkTimePidInformMark, $pid, $HmRead_StartTime);
		#
		#
		#sleep(2);
		#
		#my $reCheck_nowWorkTimePidInformMark=OnlyOnePidWork::Get_which_information_from_FlagFile ($FlagFile, $TimePidMark);
		#if ( $reCheck_nowWorkTimePidInformMark eq $nowWorkTimePidInformMark){
		#	return $nowWorkTimePidInformMark;
		#}
		#else{
		#	my $RE_couldBeStart=OnlyOnePidWork::Waiting_to_start_new_work ($FlagFile, $pid, $HmRead_StartTime, $longestTime_set_to_Wait); 
		#  if ( $RE_couldBeStart==1 ){
		#  	$nowTimeMicroSecond      =TimeWork::GetNowTimeMicroSeconds; 		                        		                     #����ʱ��
		#    $nowWorkTimePidInformMark=TimeWork::GetTimePid_microSecond    ( $nowTimeMicroSecond ); 		                       #����ʱ��
		#    $HmRead_StartTime        =TimeWork::GetHumanRead_Time         ( $nowTimeMicroSecond) ; 		                       #����ʱ��
		#    
		#    print "\n20190328-0-4-time:$nowWorkTimePidInformMark      \$nowTimeMicroSecond=$nowTimeMicroSecond     \$HmRead_StartTime=$HmRead_StartTime\n"    ;
		#    OnlyOnePidWork::Write_WorkingStart_intoFlagFile($FlagFile, $nowWorkTimePidInformMark, $pid, $HmRead_StartTime);
		#	  return $nowWorkTimePidInformMark;
	  #  }
	  #  else{
	  #  	DieWork::Just_dieWork( $die_MsgHead."\n \$RE_couldBeStart=$RE_couldBeStart\n the 2nd run not ok Still!!!!  $!\n\n\n".$caller_inform ); 
	  #  }
		#}
		
		
	}
	#########
	else {
		return '';
	}
	
}


sub Digui_checkAndWright{
	my ( $FlagFile, $pid, $HmRead_StartTime, $longestTime_set_to_Wait, $diguiNB, $diguiLimit)=@_;
	
	my $warnMsgBody="\nIn package  OnlyOnePidWork,\tIn sub Digui_checkAndWright,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	my $sleepTime=2;
	
	my $TimePidMark      ='TimePidMark';          
	my $ProcessID        ='ProcessID';            
	my $WorkIng_StartTime='WorkIng_StartTime';    
	my $WorkIng_Done_Time='WorkIng_Done_Time'; 
	
	my $couldBeStart=OnlyOnePidWork::Waiting_to_start_new_work ($FlagFile, $pid, $HmRead_StartTime, $longestTime_set_to_Wait); 
  if ( $couldBeStart==1 ){
		my $nowTimeMicroSecond      =TimeWork::GetNowTimeMicroSeconds; 		                        		                       #����ʱ��
		my $nowWorkTimePidInformMark=TimeWork::GetTimePid_microSecond    ( $nowTimeMicroSecond ); 		                       #����ʱ��
		$HmRead_StartTime           =TimeWork::GetHumanRead_Time         ( $nowTimeMicroSecond) ; 		                       #����ʱ��
		
		print "\n20190328-9-$diguiNB-time:$nowWorkTimePidInformMark      \$nowTimeMicroSecond=$nowTimeMicroSecond     \$HmRead_StartTime=$HmRead_StartTime\n"    ;
		OnlyOnePidWork::Write_WorkingStart_intoFlagFile($FlagFile, $nowWorkTimePidInformMark, $pid, $HmRead_StartTime);		
		
		sleep($sleepTime);
		
	  my $reCheck_nowWorkTimePidInformMark=OnlyOnePidWork::Get_which_information_from_FlagFile ($FlagFile, $TimePidMark);
	  if ( $reCheck_nowWorkTimePidInformMark eq $nowWorkTimePidInformMark){
	    return $nowWorkTimePidInformMark;
	  }
	  else {
	  	$diguiNB++;
	  	if ( $diguiNB > $diguiLimit ){
	  		DieWork::Just_dieWork( $die_MsgHead."\n \$diguiNB=$diguiNB  > $diguiLimit=\$diguiLimit die to save time ,and check your work!!!!  $!\n\n\n".$caller_inform ); 
	  	}
	  	else{
	  		$nowWorkTimePidInformMark=OnlyOnePidWork::Digui_checkAndWright ($FlagFile, $pid, $HmRead_StartTime, $longestTime_set_to_Wait, $diguiNB, $diguiLimit);
	  	}
	  }
	}
}

sub CheckAnd_Done_step_forOnlyOnePidWork{  # my $finished_or_not=OnlyOnePidWork::CheckAnd_Done_step_forOnlyOnePidWork( $FlagFile, $statTimPidMARK );  #���� FlagFile��·�� �� �Ѿ���ʼ�Ĺ��� TimePid����Ϣ��������� CheckAndStart_step_forOnlyOnePidWork���ɵģ�������״̬��һ��������д���µ�finished���ź�
	my ($FlagFile, $statTimPidMARK)=@_;
	
	my $warnMsgBody="\nIn package  OnlyOnePidWork,\tIn sub CheckAnd_Done_step_forOnlyOnePidWork,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	
	my $TimePidMark      ='TimePidMark';         
	my $ProcessID        ='ProcessID';           
	my $WorkIng_StartTime='WorkIng_StartTime'; 
	my $WorkIng_Done_Time='WorkIng_Done_Time';  
	
	my $finished_or_not=0;                                                                                                                                                                             print "201903281331-0-0-0-0 \$statTimPidMARK  =$statTimPidMARK\n";
	if (   (  defined ( $FlagFile )  ) && ( $FlagFile=~m/\S+/ ) && (  -e ( $FlagFile )  )   ){
    my $workIngStatuHash=OnlyOnePidWork::Read__WorkingStatu_fromFlagFile($FlagFile);                                                                                                                 print "201903281331-0-0-0-1 \$workIngStatuHash=$workIngStatuHash\n";    
		if  (   (  defined ( $workIngStatuHash )  ) && (  ref ( $workIngStatuHash ) eq 'HASH'  )   ){                                                                                                    print "201903281331-0-0-0-2 \$workIngStatuHash->{$TimePidMark}=$workIngStatuHash->{$TimePidMark}\n";    
			if (   (  defined ( $workIngStatuHash->{$TimePidMark} )  ) && ( $workIngStatuHash->{$TimePidMark}=~m/\S+/ ) && ( $workIngStatuHash->{$TimePidMark} eq $statTimPidMARK )   ){                   print "201903281331-0-0-0-3 \$workIngStatuHash->{$TimePidMark}=$workIngStatuHash->{$TimePidMark} $statTimPidMARK=\$statTimPidMARK\n";  
			  if (   (  defined ( $workIngStatuHash->{$WorkIng_Done_Time} )  ) && ( $workIngStatuHash->{$WorkIng_Done_Time}=~m/\S+/ )   ){                                                                 print "201903281331-0-0-0-4 \$workIngStatuHash->{$WorkIng_Done_Time}=$workIngStatuHash->{$WorkIng_Done_Time}\n";  
          DieWork::Just_dieWork( $die_MsgHead."\n \$workIngStatuHash->{$WorkIng_Done_Time}=$workIngStatuHash->{$WorkIng_Done_Time}\n should not be defined  $!\n\n\n".$caller_inform ); 		
			  }
			  else {			  	                                                              #warn "\$statTimPidMARK=$statTimPidMARK\n";
			  	OnlyOnePidWork::Write_WorkingDone__intoFlagFile($FlagFile, $statTimPidMARK);			 
			  	$finished_or_not=1; 	
			  }
			}
			else{		
				DieWork::Just_dieWork( $die_MsgHead."\n \$workIngStatuHash->{$TimePidMark}=$workIngStatuHash->{$TimePidMark}\n should be defined and equal \$statTimPidMARK=$statTimPidMARK  $!\n\n\n".$caller_inform ); 		
			}
		}
		else {
		  DieWork::Just_dieWork( $die_MsgHead."\n \$FlagFile=$FlagFile\n was empty or something wrong  $!\n\n\n".$caller_inform ); 
	  }
		
	}
	else {
		DieWork::Just_dieWork( $die_MsgHead."\n \$FlagFile=$FlagFile\n was not defined or cannot open  $!\n\n\n".$caller_inform ); 

	}	
	return $finished_or_not;
}

sub Get_which_information_from_FlagFile{  #  my $whichInform_val=OnlyOnePidWork::Get_which_information_from_FlagFile ($FlagFile, $whichInformKey); #�� $FlagFile �л�ȡ $whichInformKey��ֵ
	my ( $FlagFile, $whichInformKey )=@_;
	
	my $warnMsgBody="\nIn package  OnlyOnePidWork,\tIn sub Get_which_information_from_FlagFile,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	my $sleeplimt=10;
	my $sleepTime=0;      
	my $whichInform_val;
	WLMK: while (1){
		
	
	  if (   (  defined ( $FlagFile )  ) && ( $FlagFile=~m/\S+/ ) && (  -e ( $FlagFile )  )   ){
      my $FlagHASH=OnlyOnePidWork::Read__WorkingStatu_fromFlagFile($FlagFile);
	  	if  (   (  defined ( $FlagHASH )  ) && (  ref ( $FlagHASH ) eq 'HASH'  )  && (  defined ( $FlagHASH->{$whichInformKey} )  ) && (  $FlagHASH->{$whichInformKey}=~m/\S+/  )   ){
        $whichInform_val=$FlagHASH->{$whichInformKey};
        last WLMK;
	  	}
	  	else {
	  		sleep(1);     $sleepTime++;
	  		if ($sleepTime >= $sleeplimt){
	  			DieWork::Just_dieWork( $die_MsgHead."\n \$sleepTime=$sleepTime >= $sleeplimt=\$sleeplimt \n\$FlagFile=$FlagFile \$whichInformKey=$whichInformKey\n \$FlagHASH=$FlagHASH\n\$FlagHASH->{$whichInformKey}=$FlagHASH->{$whichInformKey} were not right!!!  : $!\n\n\n".$caller_inform ); 
	  		}
	  		#
	  	}
	  }
	  else{
	  	DieWork::Just_dieWork( $die_MsgHead."\n \$FlagFile=$FlagFile \$whichInformKey=$whichInformKey\n was not right!!!  : $!\n\n\n".$caller_inform );
	  }  
	  
	}
	
	return $whichInform_val;
	
}


sub Waiting_to_start_new_work{  #  my $couldBeStart=OnlyOnePidWork::Waiting_to_start_new_work ($FlagFile, $pid, $HmRead_StartTime, $longestTime_set_to_Wait); #������� $FlagFile��״̬��ֱ���ﵽ�ɹ���״̬���򷵻ؿɹ���״̬�ź�
	my ($FlagFile, $pid, $HmRead_StartTime, $longestTime_set_to_Wait)=@_;
	
	my $warnMsgBody="\nIn package  OnlyOnePidWork,\tIn sub Waiting_to_start_new_work,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	my $TimePidMark      ='TimePidMark';         
	my $ProcessID        ='ProcessID';           
	my $WorkIng_StartTime='WorkIng_StartTime'; 
	my $WorkIng_Done_Time='WorkIng_Done_Time';  
	
	my $org_TimePidMark=OnlyOnePidWork::Get_which_information_from_FlagFile ($FlagFile, $TimePidMark);
	my $flagFile_changed_by_other_proscess=0;
	
	my $couldBeStart=0;
	WHLMK: while (1){
	  
	 
	  my ( $noDefined_Empty_or_done, $DonePidMark );
	  my $tempDoneHash=OnlyOnePidWork::CheckTheFlagFile_done($FlagFile);
	  ( $noDefined_Empty_or_done, $DonePidMark )=@{ $tempDoneHash } if (   (  defined ( $tempDoneHash )  ) && (  ref ( $tempDoneHash ) eq 'ARRAY')   ); warn "\$FlagFile=$FlagFile \$noDefined_Empty_or_done=$noDefined_Empty_or_done \$DonePidMark=$DonePidMark\n";
	  
	  if    ( $noDefined_Empty_or_done==2 ){  #��flag�ļ�Ϊ�յ�ʱ�� ��1��
	  	sleep (1);
	  	next WHLMK;
	  }
	  elsif (  ( $noDefined_Empty_or_done==1 )   ){  #��(1)flag�ļ�Ϊdone��״̬��ʱ�� �� ��2��Flag�ļ�û�ж����ʱ�� ��ʾ״̬Ϊ���Կ�ʼ start   #
	  	#my $Real_CouldBe_Start=OnlyOnePidWork::Re_re_re_check_FlagFileDone($FlagFile, $pid, $DonePidMark);
	  	#if ($Real_CouldBe_Start==1){
	  		$couldBeStart=1;
	  		last WHLMK;
	  	#}
	  }
	  else{
	  	
	  	my $nowWtPointTime=TimeWork::GetHumanRead_NOW_Time (); 
	  	
	  	
	  	my $how_long_thisPid_waiting   =TimeWork::GetLastTime            ( $nowWtPointTime, $HmRead_StartTime );
	  	my $readable_last_pidWtgtime   =TimeWork::HumanRead_Last_Time    ( $how_long_thisPid_waiting          );
	  	
	  	
	  	my $now_flagFile_hash       =OnlyOnePidWork::Read__WorkingStatu_fromFlagFile($FlagFile);
	  	my $now_flagFile_TimePidMark=$now_flagFile_hash->{$TimePidMark      };
	  	my $now_flagFile_StartTime  =$now_flagFile_hash->{$WorkIng_StartTime};
	  		
	  	my $how_long_from_lastWorkStart=TimeWork::GetLastTime        ( $nowWtPointTime, $now_flagFile_StartTime );
	  	my $readable_last_time         =TimeWork::HumanRead_Last_Time( $how_long_from_lastWorkStart             );
	  	
				  
		  warn "\n\n\nNow, this \$pid=$pid have been Waiting for  $how_long_from_lastWorkStart seconds  =  $readable_last_time \n\n";
	  	
	  	if ( $how_long_from_lastWorkStart > $longestTime_set_to_Wait ){
	  	#if ( $how_long_thisPid_waiting > $longestTime_set_to_Wait ){
	  		
	  		
	  		
	  		#if ( $now_flagFile_TimePidMark eq $org_TimePidMark){  #��ʾ Ŀǰ�Ľ��� �������ֵȴ���ʱ�䳬ʱ�ˡ����ʱ�򣬱Ƚ��ʼ��¼�� flagfile��timePidMark�����һ���Ƿ� �� flagfile���������̸����ˡ�û�и��£���die����������ʾ��Ϣ����������ˣ������else�Ĳ�����
	  			my $old_work_statu_string      =InFileHandle::readAllfileIntoAstring($FlagFile);
	  	    my $warnMsg0="0  equal:\n\$now_flagFile_TimePidMark=$now_flagFile_TimePidMark eq \$org_TimePidMark=$org_TimePidMark\n";  
	  	    my $warnMsg1="1  Input:\n\$FlagFile=$FlagFile\n\$pid=$pid\n\$HmRead_StartTime=$HmRead_StartTime\n";
	  	    my $warnMsg2="2  The last work information:\n$old_work_statu_string\n\n";
	  	    my $warnMsg3="3  No change was made to FlagFile from the satrt!!!";
	  	    $warnMsg3="3 The FlagFile has been upgrade by other process $flagFile_changed_by_other_proscess times!!!" if ( $flagFile_changed_by_other_proscess > 0);	  	       
	  	    my $warnMsg4="4  Waiting for too long time, Check it out !!\n\n\$longestTime_set_to_Wait=$longestTime_set_to_Wait  <  $how_long_thisPid_waiting=\$how_long_thisPid_waiting\n";
	  	    my $warnMsg5="5  Waiting for too long time, Check it out !!\n\n\$longestTime_set_to_Wait=$longestTime_set_to_Wait  <  $how_long_from_lastWorkStart=\$how_long_from_lastWorkStart\n";
	  	    my $warnMsg6="6  Form that last work it past $how_long_from_lastWorkStart seconds = $readable_last_time  $!\n";
	  		  DieWork::Just_dieWork( $die_MsgHead."\n\n$warnMsg0\n$warnMsg1\n\n$warnMsg2\n$warnMsg3\n$warnMsg4\n$warnMsg5\n$warnMsg6\n\n\n".$caller_inform ); 
	  		#}
	  		#else{  #��ʾ flagfile���������̸����ˡ� �򽫼���ȴ�ʱ�䣬����Ϊ�µ�FlagFile�е�StartTime��ʱ��
	  			$HmRead_StartTime=$now_flagFile_StartTime;
	  			$flagFile_changed_by_other_proscess++;
	  		#}
	  		
	  	}  
	  	
	    sleep (1);	
	  }
		
	}
	
	return $couldBeStart;
}


#�������Re_re_re_check_FlagFileDone �ϳ����ˣ�û����
sub Re_re_re_check_FlagFileDone{ #  my $Real_CouldBe_Start=OnlyOnePidWork::Re_re_re_check_FlagFileDone(  $FlagFile, $pid  ); #����3�λ����Σ�����Ƿ����������ʼ�µ� ���й��̣����μ��ļ����0.3����� pid��100000��֮�롣
	my ($FlagFile, $pid, $DonePidMark)=@_;
	
	my $warnMsgBody="\nIn package  OnlyOnePidWork,\tIn sub Re_re_re_check_FlagFileDone,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	my $TimePidMark      ='TimePidMark';         
	my $ProcessID        ='ProcessID';           
	my $WorkIng_StartTime='WorkIng_StartTime'; 
	my $WorkIng_Done_Time='WorkIng_Done_Time';
	
	my $Real_CouldBe_Start=0;
	
	my $pidMicroSeconds=$pid/1000000;
	
	
	
	my $allDone=0; my $roundNb=10;
	FORMK: for (my $i=0; $i<$roundNb; $i++){
		
		
		my ( $noDefined_Empty_or_done, $re_re_re_DonePidMark );
	  my $tempDoneHash=OnlyOnePidWork::CheckTheFlagFile_done($FlagFile);
	  
	  ( $noDefined_Empty_or_done, $re_re_re_DonePidMark )=@{ $tempDoneHash } if (   (  defined ( $tempDoneHash )  ) && (  ref ( $tempDoneHash ) eq 'ARRAY'  )   );
	  if    ( $noDefined_Empty_or_done==2 ){
	  	sleep (1);
	  	next FORMK;
	  }
	  elsif ( $noDefined_Empty_or_done==1 ){
	  	if (   (  defined ( $DonePidMark )  ) && ( $DonePidMark=~m/\S+/ ) && (  defined ( $re_re_re_DonePidMark )  ) && ( $re_re_re_DonePidMark=~m/\S+/ )   ){  print "20190328-0-0-time:$DonePidMark $i ", TimeWork::GetNowTimeMicroSeconds (), ", \$DonePidMark=$DonePidMark eq ?  eq $re_re_re_DonePidMark=\$re_re_re_DonePidMark\n"    ;
	  		if  ( $DonePidMark eq $re_re_re_DonePidMark )  {                                                                                                      print "20190328-0-1-time:$DonePidMark $i ", TimeWork::GetNowTimeMicroSeconds (), ", \$DonePidMark=$DonePidMark     eq   $re_re_re_DonePidMark=\$re_re_re_DonePidMark\n"    ;
	  		  $allDone++;		                     
	  		}
	  	  else{
	  	  	                                                                                                                                                    print "20190328-0-2-time:$DonePidMark $i ", TimeWork::GetNowTimeMicroSeconds (), ", \$DonePidMark=$DonePidMark     ne   $re_re_re_DonePidMark=\$re_re_re_DonePidMark\n"    ;
	  	  }
	  	}
	  }
	  my $sleepTime=0.3+$pidMicroSeconds; print "20190328-1-0-time:$DonePidMark  $i   ", TimeWork::GetNowTimeMicroSeconds (), ", \$sleepTime=$sleepTime \$allDone=$allDone =?= $roundNb=\$roundNb $DonePidMark=\$DonePidMark\n"    ;
	                                      print "20190328-1-1-time:$DonePidMark  $i   ", TimeWork::GetNowTimeMicroSeconds (), ", \$sleepTime=$sleepTime \$allDone=$allDone =?= $roundNb=\$roundNb $re_re_re_DonePidMark=\$re_re_re_DonePidMark\n"    ;
		sleep ($sleepTime);
	}
	if ( $allDone == $roundNb ){
		$Real_CouldBe_Start=1;              print "20190328-2-0-time:$DonePidMark         ", TimeWork::GetNowTimeMicroSeconds (), ", \$Real_CouldBe_Start=$Real_CouldBe_Start\n"    ;   
	}
	return $Real_CouldBe_Start;
}

sub CheckTheFlagFile_noDefined_Empty{  # $noDefined_Empty_or_done=OnlyOnePidWork::CheckTheFlagFile_noDefined_Empty($FlagFile);  #����Ƿ�ﵽ �¹��� ���Կ�ʼ��״̬���� 1 �ϸ�������״̬��¼�ļ�û�У�2 Ϊ�� 
	my ($FlagFile)=@_;
	
	my $warnMsgBody="\nIn package  OnlyOnePidWork,\tIn sub CheckTheFlagFile_noDefined_Empty,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	my $TimePidMark      ='TimePidMark';         
	my $ProcessID        ='ProcessID';           
	my $WorkIng_StartTime='WorkIng_StartTime'; 
	my $WorkIng_Done_Time='WorkIng_Done_Time';  
	
	my $noDefined_Empty_or_done=0;                                                                                      #warn "                     201903290914-0-0-0 \$noDefined_Empty_or_done=$noDefined_Empty_or_done\n";
  if (   (  defined ( $FlagFile )  ) && ( $FlagFile=~m/\S+/ ) && (  -e ( $FlagFile )  )   ){  
    my $workIngStatuHash=OnlyOnePidWork::Read__WorkingStatu_fromFlagFile($FlagFile);                                  #warn "                     201903290914-0-0-1 \$workIngStatuHash=$workIngStatuHash\n";
		if  (   (  defined ( $workIngStatuHash )  ) && (  ref ( $workIngStatuHash ) eq 'HASH'  )   ){                     #warn "                     201903290914-0-0-2 \$workIngStatuHash=$workIngStatuHash\n"; DirFileHandle::PrintAndWarnDumper ($workIngStatuHash, '201903290914-0-0-3');			
		}
		else {
		  $noDefined_Empty_or_done=1;
	  }
		
	}
	else {
		$noDefined_Empty_or_done=1;
	}
	return $noDefined_Empty_or_done;
}


sub CheckTheFlagFile_noDefined_Empty_or_done{  # $noDefined_Empty_or_done=OnlyOnePidWork::CheckTheFlagFile_noDefined_Empty_or_done($FlagFile);  #����Ƿ�ﵽ �¹��� ���Կ�ʼ��״̬���� 1 �ϸ�������״̬��¼�ļ�û�У�2 Ϊ�� 3 �Ѿ� done
	my ($FlagFile)=@_;
	
	my $warnMsgBody="\nIn package  OnlyOnePidWork,\tIn sub CheckTheFlagFile_noDefined_Empty_or_done,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	my $TimePidMark      ='TimePidMark';         
	my $ProcessID        ='ProcessID';           
	my $WorkIng_StartTime='WorkIng_StartTime'; 
	my $WorkIng_Done_Time='WorkIng_Done_Time';  
	
	my $noDefined_Empty_or_done=0;                                                                                      #warn "                     201903290914-0-0-0 \$noDefined_Empty_or_done=$noDefined_Empty_or_done\n";
  if (   (  defined ( $FlagFile )  ) && ( $FlagFile=~m/\S+/ ) && (  -e ( $FlagFile )  )   ){  
    my $workIngStatuHash=OnlyOnePidWork::Read__WorkingStatu_fromFlagFile($FlagFile);                                  #warn "                     201903290914-0-0-1 \$workIngStatuHash=$workIngStatuHash\n";
		if  (   (  defined ( $workIngStatuHash )  ) && (  ref ( $workIngStatuHash ) eq 'HASH'  )   ){                     #warn "                     201903290914-0-0-2 \$workIngStatuHash=$workIngStatuHash\n"; DirFileHandle::PrintAndWarnDumper ($workIngStatuHash, '201903290914-0-0-3');
		}
		else {
		  $noDefined_Empty_or_done=1;
	  }
		
	}
	else {
		DieWork::Just_dieWork( $die_MsgHead."\n It is empty check!!\n\$FlagFile=$FlagFile should not be undefined !!  $!\n\n\n".$caller_inform ); 
	}
	return $noDefined_Empty_or_done;
}

sub CheckTheFlagFile_done{  # $noDefined_Empty_or_done=OnlyOnePidWork::CheckTheFlagFile_done($FlagFile);  #����Ƿ�ﵽ �¹��� ���Կ�ʼ��״̬�� �Ѿ� done
	my ($FlagFile)=@_;
	
	my $warnMsgBody="\nIn package  OnlyOnePidWork,\tIn sub CheckTheFlagFile_done,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	my $TimePidMark      ='TimePidMark';         
	my $ProcessID        ='ProcessID';           
	my $WorkIng_StartTime='WorkIng_StartTime'; 
	my $WorkIng_Done_Time='WorkIng_Done_Time';  
	
	my $done=0;          
	my $DonePidMark;                                                                            #warn "                     201903290914-0-0-0 \$done=$done\n";
  if (   (  defined ( $FlagFile )  ) && ( $FlagFile=~m/\S+/ ) && (  -e ( $FlagFile )  )   ){  
    my $workIngStatuHash=OnlyOnePidWork::Read__WorkingStatu_fromFlagFile($FlagFile);                                  #warn "                     201903290914-0-0-1 \$workIngStatuHash=$workIngStatuHash\n";
		if  (   (  defined ( $workIngStatuHash )  ) && (  ref ( $workIngStatuHash ) eq 'HASH'  )   ){                     #warn "                     201903290914-0-0-2 \$workIngStatuHash=$workIngStatuHash\n"; DirFileHandle::PrintAndWarnDumper ($workIngStatuHash, '201903290914-0-0-3');
			if (   (  defined ( $workIngStatuHash->{$WorkIng_Done_Time} )  ) && ( $workIngStatuHash->{$WorkIng_Done_Time}=~m/\S+/ )   ){  
        $done=1;     
        $DonePidMark            =$workIngStatuHash->{$TimePidMark}  ;
			}
			else{				
				
			}
		}
		else {
		  #my $warnMasg= $warnMsgHead."\n \$FlagFile=$FlagFile should not be empty !!  $!\n\n\n".$caller_inform ; 
		  #print $warnMasg;
		  #warn $warnMasg;
		  #DieWork::Just_dieWork( $die_MsgHead."\n \$FlagFile=$FlagFile should not be empty !!  $!\n\n\n".$caller_inform );
		  $done=2;
	  }
		
	}
	else {
		#$done=3;
		DieWork::Just_dieWork( $die_MsgHead."\n \$FlagFile=$FlagFile should not be undefined !!  $!\n\n\n".$caller_inform ); 
	}
	return [ $done, $DonePidMark ];
}


sub CheckTheFlagFile_Empty{  # $noDefined_Empty_or_done=OnlyOnePidWork::CheckTheFlagFile_Empty($FlagFile);  #����Ƿ�ﵽ �¹��� ���Կ�ʼ��״̬���� 2 Ϊ�� 
	my ($FlagFile)=@_;
	
	my $warnMsgBody="\nIn package  OnlyOnePidWork,\tIn sub CheckTheFlagFile_Empty,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	my $TimePidMark      ='TimePidMark';         
	my $ProcessID        ='ProcessID';           
	my $WorkIng_StartTime='WorkIng_StartTime'; 
	my $WorkIng_Done_Time='WorkIng_Done_Time';  
	
	my $Empty=0;                                                                                                        #warn "                     201903290914-0-0-0 \$Empty=$Empty\n";
  if (   (  defined ( $FlagFile )  ) && ( $FlagFile=~m/\S+/ ) && (  -e ( $FlagFile )  )   ){  
    my $workIngStatuHash=OnlyOnePidWork::Read__WorkingStatu_fromFlagFile($FlagFile);                                  #warn "                     201903290914-0-0-1 \$workIngStatuHash=$workIngStatuHash\n";
		if  (   (  defined ( $workIngStatuHash )  ) && (  ref ( $workIngStatuHash ) eq 'HASH'  )   ){                     #warn "                     201903290914-0-0-2 \$workIngStatuHash=$workIngStatuHash\n"; DirFileHandle::PrintAndWarnDumper ($workIngStatuHash, '201903290914-0-0-3');			
		}
		else {
		  $Empty=1;
	  }		
	}
	else {
		DieWork::Just_dieWork( $die_MsgHead."\n It is empty check!!\n\$FlagFile=$FlagFile should not be undefined !!  $!\n\n\n".$caller_inform ); 
	}
	return $Empty;
}

sub CheckTheFlagFile_noDefined{  # $noDefined_Empty_or_done=OnlyOnePidWork::CheckTheFlagFile_noDefined($FlagFile);  #����Ƿ�ﵽ �¹��� ���Կ�ʼ��״̬���� 1 �ϸ�������״̬��¼�ļ�û�У�
	my ($FlagFile)=@_;
	
	my $warnMsgBody="\nIn package  OnlyOnePidWork,\tIn sub CheckTheFlagFile_noDefined,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	my $TimePidMark      ='TimePidMark';         
	my $ProcessID        ='ProcessID';           
	my $WorkIng_StartTime='WorkIng_StartTime'; 
	my $WorkIng_Done_Time='WorkIng_Done_Time';  
	
	my $noDefined=0;                                                                                      #warn "                     201903290914-0-0-0 \$noDefined=$noDefined\n";
  if (   (  defined ( $FlagFile )  ) && ( $FlagFile=~m/\S+/ ) && (  -e ( $FlagFile )  )   ){      
	}
	else {
		$noDefined=1;
	}
	return $noDefined;
}

sub Write_WorkingStart_intoFlagFile{   # OnlyOnePidWork::Write_WorkingStart_intoFlagFile($FlagFile, $nowWorkTimePidInformMark, $pid, $HmRead_StartTime);  # �� ��ʼ����ʱ�����Ϣ д���ļ� $FlagFile ��
	my ($FlagFile, $nowWorkTimePidInformMark, $pid, $HmRead_StartTime)=@_; 
	
	my $warnMsgBody="\nIn package  OnlyOnePidWork,\tIn sub Write_WorkingStart_intoFlagFile,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";  
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	my $TimePidMark      ='TimePidMark';          #my $nowWorkTimePidInformMark=TimeWork::GetNowTimePid_microSecond ();    
	my $ProcessID        ='ProcessID';            #my $pid                     =$$;
	my $WorkIng_StartTime='WorkIng_StartTime';    #my $HmRead_StartTime        =TimeWork::GetHumanRead_Time         ($nowWorkTimePidInformMark); 
	
  my $workIngString=	"TimePidMark\t$nowWorkTimePidInformMark\nProcessID\t$pid\nWorkIng_StartTime\t$HmRead_StartTime\n";  my $ptOutWrn="20190328-3-0-time:". TimeWork::GetNowTimeMicroSeconds (). " ,$nowWorkTimePidInformMark, \n\$workIngString=$workIngString\n\n"    ; 
	InFileHandle::PrintStringIntoFile($FlagFile, $workIngString, $ptOutWrn);  #print "$ptOutWrn"    ;   
	
	return $workIngString;
}


sub Read__WorkingStatu_fromFlagFile{   #  my $outputHash=OnlyOnePidWork::Read__WorkingStatu_fromFlagFile($FlagFile);   # ��ȡ ״̬��ʾ�ļ� $FlagFile �е���Ϣ
	my ($FlagFile)=@_;
	
	my $warnMsgBody="\nIn package  OnlyOnePidWork,\tIn sub Read__WorkingStatu_fromFlagFile,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";  
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	my $TimePidMark      ='TimePidMark';         
	my $ProcessID        ='ProcessID';           
	my $WorkIng_StartTime='WorkIng_StartTime'; 
	my $WorkIng_Done_Time='WorkIng_Done_Time';  
	
	my $outputHash;
	my $workIngString=InFileHandle::readAllfileIntoAstring($FlagFile);  #warn "\$workIngString=$workIngString\n";
	if ( $workIngString=~/   
	                         $TimePidMark      \t(.*)\n
	                         $ProcessID        \t(.*)\n
	                         $WorkIng_StartTime\t(.*)\n
	                         (?:
	                           $WorkIng_Done_Time\t(.*)\n
	                         )?
	                         
	                     /x ){  
	  $outputHash->{$TimePidMark       }=$1;
	  $outputHash->{$ProcessID         }=$2;
	  $outputHash->{$WorkIng_StartTime }=$3;
	  my $four_here=$4;
	  
	  $outputHash->{$WorkIng_Done_Time }=$four_here if (   (  defined ( $four_here )  ) && ( $four_here=~m/\S+/ )   ) ;
			
  }
  else{
  	#DieWork::Just_dieWork( $die_MsgHead."\n \$workIngString=$workIngString didnot fit the regular expression !!  $!\n\n\n".$caller_inform ); 
  	my $warnMSG= $warnMsgHead."\n \$workIngString=$workIngString didnot fit the regular expression !!  $!\n\n\n".$caller_inform ; 
  	#print     $warnMSG ;
  	#warn      $warnMSG ;
  }
	
	return $outputHash;
}


sub Write_Psudo_WorkingDone__intoFlagFile{   # OnlyOnePidWork::Write_Psudo_WorkingDone__intoFlagFile($FlagFile, $nowWorkTimePidInformMark, $pid, $HmRead_StartTime);  
	#���ֻ�� ��Flagfileû�б������ʱ�򣬴�����д���һ���ٵ� ��ɵ� Flag����Ϊ�˺����� ���״̬�ļ�� �ܹ���ִ��
	my ($FlagFile, $nowWorkTimePidInformMark, $pid, $HmRead_StartTime)=@_; 
	
	my $warnMsgBody="\nIn package  OnlyOnePidWork,\tIn sub Write_Psudo_WorkingDone__intoFlagFile,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";  
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	my $TimePidMark      ='TimePidMark';          #my $nowWorkTimePidInformMark=TimeWork::GetNowTimePid_microSecond ();    
	my $ProcessID        ='ProcessID';            #my $pid                     =$$;
	my $WorkIng_StartTime='WorkIng_StartTime';    #my $HmRead_StartTime        =TimeWork::GetHumanRead_Time         ($nowWorkTimePidInformMark); 
	
	my $WorkIng_Done_Time='WorkIng_Done_Time';     my $HmRead_Done_Time =TimeWork::GetHumanRead_NOW_Time         ();     #warn "\$HmRead_Done_Time=$HmRead_Done_Time\n";  ####д����ɹ�����ʱ��
	
  my $workIngString=	"TimePidMark\t$nowWorkTimePidInformMark\nProcessID\t$pid\nWorkIng_StartTime\t$HmRead_StartTime\n$WorkIng_Done_Time\t$HmRead_Done_Time\n";  my $ptOutWrn="20190328-3-1-time:". TimeWork::GetNowTimeMicroSeconds (). " ,$nowWorkTimePidInformMark, \n\$workIngString=$workIngString\n\n"    ; 
	InFileHandle::PrintStringIntoFile($FlagFile, $workIngString, $ptOutWrn);  #print "$ptOutWrn"    ;   
	
	return $workIngString;
}

sub Write_WorkingDone__intoFlagFile{  #  OnlyOnePidWork::Write_WorkingDone__intoFlagFile($FlagFile, $statTimPidMARK);
  my ($FlagFile, $statTimPidMARK)=@_;
	
	my $warnMsgBody="\nIn package  OnlyOnePidWork,\tIn sub Write_WorkingDone__intoFlagFile,\n\n";	
  my $warnMsgHead="\n\n\n$warnMsgBody";
	my $die_MsgHead="\n\nDIE!!!!!\n$warnMsgBody";  
	my $caller_inform=DirFileHandle::print_SubCallerInform;
	
	my $TimePidMark      ='TimePidMark';         
	my $ProcessID        ='ProcessID';           
	my $WorkIng_StartTime='WorkIng_StartTime'; 
	                                                                
	my $WorkIng_Done_Time='WorkIng_Done_Time';     my $HmRead_Done_Time =TimeWork::GetHumanRead_NOW_Time         ();     #warn "\$HmRead_Done_Time=$HmRead_Done_Time\n";  ####д����ɹ�����ʱ��
	
	if (   (  defined ( $FlagFile )  ) && ( $FlagFile=~m/\S+/ ) && (  -e ( $FlagFile )  )   ){
		my $workIngStatuHash=OnlyOnePidWork::Read__WorkingStatu_fromFlagFile($FlagFile);                   #w#arn "\$workIngStatuHash=$workIngStatuHash\n";
		if  (   (  defined ( $workIngStatuHash )  ) && (  ref ( $workIngStatuHash ) eq 'HASH'  )   ){
			if (   (  defined ( $workIngStatuHash->{$WorkIng_Done_Time} )  ) && ( $workIngStatuHash->{$WorkIng_Done_Time}=~m/\S+/ )   ){
				DieWork::Just_dieWork( $die_MsgHead."\n \$workIngStatuHash->{$WorkIng_Done_Time}=$workIngStatuHash->{$WorkIng_Done_Time} should not be pre-defined, it is wrong !!  $!\n\n\n".$caller_inform ); 
			}
			else{
				my $workIngStatuString=InFileHandle::readAllfileIntoAstring($FlagFile);                        #warn "\$workIngStatuString=$workIngStatuString\n";
				$workIngStatuString.="$WorkIng_Done_Time\t$HmRead_Done_Time\n";   my $ptOutWrn="20190328-4-0-time:". TimeWork::GetNowTimeMicroSeconds (). " ,$statTimPidMARK, \n\$workIngStatuString=$workIngStatuString\n\n"    ; 
				InFileHandle::PrintStringIntoFile($FlagFile, $workIngStatuString, $ptOutWrn);
			}
		}
		else{
		  DieWork::Just_dieWork( $die_MsgHead."\n \$FlagFile=$FlagFile is not Defined or empty or no such file !!  $!\n\n\n".$caller_inform ); 
	  }
		
	}
	else{
		DieWork::Just_dieWork( $die_MsgHead."\n \$FlagFile=$FlagFile is not Defined or empty or no such file !!  $!\n\n\n".$caller_inform ); 
	}
  

}


1;