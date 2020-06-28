#!/usr/bin/perl -w
BEGIN{  push (@INC, '/home/fredjiang/fredProgs/PMs')  };


use strict;
use warnings;

use POSIX ":sys_wait_h";
use DieWork;

package MultiThreadingWork;


# for each mutiple thread working perl script, the method below was used to figure out howmany core to really use!!
# my $user_Set_Core=MultiThreadingWork::HowManyCore_to_USE ( $maxCoreLimit, $default_Core[, $user_Set_Core ] );
sub HowManyCore_to_USE{
	my  ( $maxCoreLimit, $default_Core, $user_Set_Core )=@_;
	
	my ($warnMsgBody, $warnMsgHead, $die_MsgHead, $caller_inform)=@{ DieWork::BuildWarnDieHeadInform( 'MultiThreadingWork', 'HowManyCore_to_USE' ) };
  
  DieWork::Check_DfdNoEmptNUMBER_or_DIE( $maxCoreLimit, "\$maxCoreLimit", $die_MsgHead, $caller_inform  );
  DieWork::Check_DfdNoEmptNUMBER_or_DIE( $default_Core, "\$default_Core", $die_MsgHead, $caller_inform  );
	if (  DieWork::Check_DfdNoEmptNUMBER_or_NOT ( $user_Set_Core )  ){} 
	else{ $user_Set_Core=$default_Core;   DieWork::Print_and_warn( $warnMsgHead."\n\$user_Set_Core=$user_Set_Core is not set as a intergar number, it changed into the default number \$default_Core=$default_Core\n\n\n" );  }
	
	return $user_Set_Core;
}



#sub MuiltiThreadWork{       #MultiThreadingWork::MuiltiThreadWork ($inArray, $inSub);
#	my ($inArray, $inSub, $coreNb2Use)=@_;
#	my $chdPidNb=0;
#  my $num_proc = 0;                     ## == number of proc ==
#  my $num_collect = 0;                  ## == number of collected ==
#  my $collect;
#  
#  $SIG{CHLD} = \&proNBminus;            ## == get the child signal ==    
#  sub proNBminus{ 
#  	$num_proc--;                      	#print "\naaaaa -- \$num_proc=$num_proc\n"; warn "\naaaaa -- \$num_proc=$num_proc\n";
#  } 
#	my $runningPidHash;
#  if (    (  defined ( $inArray )  ) && (  ref ( $inArray ) eq 'ARRAY'  )   ){
#  	foreach my $eachHash (  @{ $inArray }  ){
#  		$chdPidNb++;
#			my $pid=fork();    
#			if (  !defined( $pid )  )  {
#        print "\naaaaa Error in fork: $!";
#        exit 1;
#      }
#      if ($pid == 0) {
#        ## == child proc ==           #print "\naaaaa Child $chdPidNb : My pid = $$\n";  
#        warn "\nChild $chdPidNb : My pid = $$\n";
#        
#        &$inSub ($eachHash);
#        
#        
#        
#        exit 0;  
#      }
#  		
#  		
#  		$runningPidHash->{$pid}=1;
#			Marktwo: while(1) {
#			  my $runningPidNB=0; 
#			  my $updatePidHash;
#			  foreach my $pidKey (    sort {$a<=>$b}(   keys(  %{ $runningPidHash }  )   )    ){
#			  	if ( (waitpid($pidKey, WNOHANG)) > 0 ){my $tpVal=waitpid($pidKey, WNOHANG); }#  print "\naaaaaa\$pidKey=$pidKey is dead!!\t\$tpVal=$tpVal\n\n";}
#			  	else {$updatePidHash->{$pidKey}=1; $runningPidNB++;                         }# print "\naaaaaa\$runningPidNB child pid is runing!!\t$pidKey=$pidKey is still running!!! \n\n";}
#			  }
#			  warn "\n\n                                     [       $runningPidNB process running now !      ]      \n";   #print "\naaaaaa Totally \$runningPidNB=$runningPidNB is running!\n"; 
#			  $runningPidHash=$updatePidHash;
#			  if ($runningPidNB<$coreNb2Use){ warn"                   [       process ID: $pid\t $runningPidNB process running! <  We set $coreNb2Use CPU core to use !! Countinue to next fork!      ]\n\n";last Marktwo;    }
#			  else                          { warn"                   [       process ID: $pid\t $runningPidNB process running! >= We set $coreNb2Use CPU core to use !! So Just Sleep 1 seconds!      ]\n\n"; sleep(1);   }
#			}
#  		
#  		
#  	}
#  	
#  	my $time_1=`date`; chomp $time_1;
#    warn "\n\n             [       Start to wait !!      ]\n\n             [       ", $time_1, "      ]\n\n" ;
#    
#    Markthree: while(1) {
#      my $runningPidNB=0; 
#      my $updatePidHash;
#      foreach my $pidKey (    sort {$a<=>$b}(   keys(  %{ $runningPidHash }  )   )    ){
#      	if ( (waitpid($pidKey, WNOHANG)) > 0 ){my $tpVal=waitpid($pidKey, WNOHANG); my $ps=`ps $pidKey`; print "\n             [       \$pidKey=$pidKey is dead!!\t\$tpVal=$tpVal\t\$ps=$ps      ]\n\n"; warn "\na             [       \$pidKey=$pidKey is dead!!\t\$tpVal=$tpVal\n\$ps=$ps      ]\n\n";}
#      	else {
#      		my $ps=`ps $pidKey`; 
#      		if ($ps=~m/$pidKey/){
#      		  $updatePidHash->{$pidKey}=1; $runningPidNB++;  
#      		  print "\n             [       $runningPidNB child pid is runing!!\t",(waitpid($pidKey, WNOHANG)),"\t\$pidKey=$pidKey is still running!!!       ]\n\$ps=$ps\n\n"; 
#      		  warn  "\n             [       $runningPidNB child pid is runing!!\t",(waitpid($pidKey, WNOHANG)),"\t\$pidKey=$pidKey is still running!!!       ]\t\$ps=$ps\n\n";
#      		}else {warn "\n             [       Running Pid \$pidKey=$pidKey found\tWAITPID=",(waitpid($pidKey, WNOHANG)),"\t, but no ps was found!!\$ps=$ps      ]\n\n"; print "\n             [       Running Pid \$pidKey=$pidKey found\tWAITPID=",(waitpid($pidKey, WNOHANG)),"\t, but no ps was found!!\$ps=$ps      ]\n\n";}
#      	}
#      }
#      
#      warn "\n             [       Totally  \$runningPidNB=$runningPidNB      ]\n";   print "\n             [       Totally \$runningPidNB=$runningPidNB is running!\n"; 
#      $runningPidHash=$updatePidHash;
#      if ($runningPidNB>0){ print"             [       \$runningPidNB=$runningPidNB > 0 !! Sleep 1 seconds!\tCountinue to next Wait checking run!      ]\n\n"; sleep(1);   }   
#      else                 {print"             [       \$runningPidNB=$runningPidNB <= 0 !! all child pid dead      ]\n\n";  sleep(1);last Markthree;   }
#    }
#    
#    my $waitPid=wait(); my $time_2=`date`; chomp $time_2;  warn "             [       Start Wait AT:      ]\n             [       $time_1      ]\n\n             [       End of Wait, \$waitPid=$waitPid!!      ]\n\n             [       ", $time_2, "      ]\n\n" ;
#
#  	
#  }
#}
#
#
#
#
#my $runningPidHash;
#
#   $chdPidNb++;
#				  my $pid=fork();    
#				  if (!defined($pid)) {
#            print "\naaaaa Error in fork: $!";
#          exit 1;
#          }
#         
#          if ($pid == 0) {
#            ## == child proc ==           #print "\naaaaa Child $chdPidNb : My pid = $$\n";  
#            warn "\nChild $chdPidNb : My pid = $$\n";
#    
#    	    }
#    	    
#    	    $runningPidHash->{$pid}=1;
#				  Marktwo: while(1) {
#				    my $runningPidNB=0; 
#				    my $updatePidHash;
#				    foreach my $pidKey (    sort {$a<=>$b}(   keys(  %{ $runningPidHash }  )   )    ){
#				    	if ( (waitpid($pidKey, WNOHANG)) > 0 ){my $tpVal=waitpid($pidKey, WNOHANG); }#  print "\naaaaaa\$pidKey=$pidKey is dead!!\t\$tpVal=$tpVal\n\n";}
#				    	else {$updatePidHash->{$pidKey}=1; $runningPidNB++;                         }# print "\naaaaaa\$runningPidNB child pid is runing!!\t$pidKey=$pidKey is still running!!! \n\n";}
#				    }
#				    warn "\n\n                                     [       $runningPidNB process running now !      ]      \n";   #print "\naaaaaa Totally \$runningPidNB=$runningPidNB is running!\n"; 
#				    $runningPidHash=$updatePidHash;
#				    if ($runningPidNB<$coreNb2Use){ warn"                   [       process ID: $pid\t $runningPidNB process running! <  We set $coreNb2Use CPU core to use !! Countinue to next fork!      ]\n\n";last Marktwo;    }
#				    else                          { warn"                   [       process ID: $pid\t $runningPidNB process running! >= We set $coreNb2Use CPU core to use !! So Just Sleep 1 seconds!      ]\n\n"; sleep(1);   }
#				  }
#
#
#my $time_1=`date`; chomp $time_1;
#warn "\n\n             [       Start to wait !!      ]\n\n             [       ", $time_1, "      ]\n\n" ;
#
#Markthree: while(1) {
#  my $runningPidNB=0; 
#  my $updatePidHash;
#  foreach my $pidKey (    sort {$a<=>$b}(   keys(  %{ $runningPidHash }  )   )    ){
#  	if ( (waitpid($pidKey, WNOHANG)) > 0 ){my $tpVal=waitpid($pidKey, WNOHANG); my $ps=`ps $pidKey`; print "\n             [       \$pidKey=$pidKey is dead!!\t\$tpVal=$tpVal\t\$ps=$ps      ]\n\n"; warn "\na             [       \$pidKey=$pidKey is dead!!\t\$tpVal=$tpVal\n\$ps=$ps      ]\n\n";}
#  	else {
#  		my $ps=`ps $pidKey`; 
#  		if ($ps=~m/$pidKey/){
#  		  $updatePidHash->{$pidKey}=1; $runningPidNB++;  
#  		  print "\n             [       $runningPidNB child pid is runing!!\t",(waitpid($pidKey, WNOHANG)),"\t\$pidKey=$pidKey is still running!!!       ]\n\$ps=$ps\n\n"; 
#  		  warn  "\n             [       $runningPidNB child pid is runing!!\t",(waitpid($pidKey, WNOHANG)),"\t\$pidKey=$pidKey is still running!!!       ]\t\$ps=$ps\n\n";
#  		}else {warn "\n             [       Running Pid \$pidKey=$pidKey found\tWAITPID=",(waitpid($pidKey, WNOHANG)),"\t, but no ps was found!!\$ps=$ps      ]\n\n"; print "\n             [       Running Pid \$pidKey=$pidKey found\tWAITPID=",(waitpid($pidKey, WNOHANG)),"\t, but no ps was found!!\$ps=$ps      ]\n\n";}
#  	}
#  }
#  
#  warn "\n             [       Totally  \$runningPidNB=$runningPidNB      ]\n";   print "\n             [       Totally \$runningPidNB=$runningPidNB is running!\n"; 
#  $runningPidHash=$updatePidHash;
#  if ($runningPidNB>0){ print"             [       \$runningPidNB=$runningPidNB > 0 !! Sleep 1 seconds!\tCountinue to next Wait checking run!      ]\n\n"; sleep(1);   }   
#  else                 {print"             [       \$runningPidNB=$runningPidNB <= 0 !! all child pid dead      ]\n\n";  sleep(1);last Markthree;   }
#}
#
#my $waitPid=wait(); my $time_2=`date`; chomp $time_2;  warn "             [       Start Wait AT:      ]\n             [       $time_1      ]\n\n             [       End of Wait, \$waitPid=$waitPid!!      ]\n\n             [       ", $time_2, "      ]\n\n" ;








sub Step_0_MultipleThread{    # my $num_proc = 0;  $SIG{CHLD} = \MultiThreadingWork::Step_0_MultipleThread($num_proc);
  my ($in_num_proc)=@_;
	$in_num_proc--;                print "\naaaaa -- \$in_num_proc=$in_num_proc\n"; warn "\naaaaa -- \$in_num_proc=$in_num_proc\n";
}


sub Step_1_MultipleThread{  #my ($chdPidNb, $runningPidHash)=@{ MultiThreadingWork::Step_1_MultipleThread };
		
	my $chdPidNb=0;
  #my $num_proc = 0;                     ## == number of proc ==
  #my $num_collect = 0;                  ## == number of collected ==
  #my $collect;   
	my $runningPidHash;
	
	
	return [
	          $chdPidNb,
	#          $num_proc,
	#          $num_collect,
	 #         $collect,
	          $runningPidHash
	       ];
}

sub Step_2_MultipleThread{  #  my $pid; ($chdPidNb, $pid)=@{ MultiThreadingWork::Step_2_MultipleThread($chdPidNb) };
	my (
	      $chdPidNb,	      
	   )=@_;
	
	$chdPidNb++;
	my $pid=fork();    
  if ( !defined($pid) ) {
    print "\naaaaa Error in fork: $!";
    exit 1;
  }
  return [
           $chdPidNb,
           $pid
         ];
}


#sub Step_3_4_MultipleThread{
  #if ($pid==0){
  	
  	
  # 	exit 0;  
  #}
#}


sub Step5_MultipleThread { #   $runningPidHash=  MultiThreadingWork::Step5_MultipleThread($runningPidHash, $pid, $coreNb2Use);
	my (
	      $runningPidHash,	
	      $pid,
	      $coreNb2Use
	            
	   )=@_;
	   
	$runningPidHash->{$pid}=1; warn "\$runningPidHash->{$pid}=$runningPidHash->{$pid}\n";
	my $sigHere='WNOHANG';
	Marktwo: while(1) {
	  my $runningPidNB=0; 
	  my $updatePidHash;
	  foreach my $pidKey (    sort {$a<=>$b}(   keys(  %{ $runningPidHash }  )   )    ){
	  	if ( (waitpid($pidKey, $sigHere) ) > 0 ){my $tpVal=waitpid($pidKey, $sigHere);   warn "\naaaaaa\$pidKey=$pidKey is dead!!\t\$tpVal=$tpVal\n\n";}
	  	else {$updatePidHash->{$pidKey}=1; $runningPidNB++;                           warn "\naaaaaa\$runningPidNB child pid is runing!!\t$pidKey=$pidKey is still running!!! \n\n";}
	  }
	  warn "\n\n                                     [       $runningPidNB process running now !      ]      \n";   #print "\naaaaaa Totally \$runningPidNB=$runningPidNB is running!\n"; 
	  $runningPidHash=$updatePidHash; 
	  if ($runningPidNB<$coreNb2Use){ warn"                   [       process ID: $pid\t $runningPidNB process running! <  We set $coreNb2Use CPU core to use !! Countinue to next fork!      ]\n\n";last Marktwo;    }
	  else                          { warn"                   [       process ID: $pid\t $runningPidNB process running! >= We set $coreNb2Use CPU core to use !! So Just Sleep 1 seconds!      ]\n\n"; sleep(1);   }
	}
	return $runningPidHash;
}

sub Step6_MultipleThread {  #MultiThreadingWork::Step6_MultipleThread($runningPidHash);
	my (
	      $runningPidHash,	
	   )=@_;
	   
	
	my $time_1=`date`; chomp $time_1;
  warn "\n\n             [       Start to wait !!      ]\n\n             [       ", $time_1, "      ]\n\n" ;
  
  Markthree: while(1) {
    my $runningPidNB=0; 
    my $updatePidHash;
    foreach my $pidKey (    sort {$a<=>$b}(   keys(  %{ $runningPidHash }  )   )    ){
    	if ( (waitpid($pidKey, 'WNOHANG')) > 0 ){my $tpVal=waitpid($pidKey, 'WNOHANG'); my $ps=`ps $pidKey`; print "\n             [       \$pidKey=$pidKey is dead!!\t\$tpVal=$tpVal\t\$ps=$ps      ]\n\n"; warn "\na             [       \$pidKey=$pidKey is dead!!\t\$tpVal=$tpVal\n\$ps=$ps      ]\n\n";}
    	else {
    		my $ps=`ps $pidKey`; 
    		if ($ps=~m/$pidKey/){
    		  $updatePidHash->{$pidKey}=1; $runningPidNB++;  
    		  print "\n             [       $runningPidNB child pid is runing!!\t",(waitpid($pidKey, 'WNOHANG')),"\t\$pidKey=$pidKey is still running!!!       ]\n\$ps=$ps\n\n"; 
    		  warn  "\n             [       $runningPidNB child pid is runing!!\t",(waitpid($pidKey, 'WNOHANG')),"\t\$pidKey=$pidKey is still running!!!       ]\t\$ps=$ps\n\n";
    		}else {warn "\n             [       Running Pid \$pidKey=$pidKey found\tWAITPID=",(waitpid($pidKey, 'WNOHANG')),"\t, but no ps was found!!\$ps=$ps      ]\n\n"; print "\n             [       Running Pid \$pidKey=$pidKey found\tWAITPID=",(waitpid($pidKey, 'WNOHANG')),"\t, but no ps was found!!\$ps=$ps      ]\n\n";}
    	}
    }
    
    warn "\n             [       Totally  \$runningPidNB=$runningPidNB      ]\n";   print "\n             [       Totally \$runningPidNB=$runningPidNB is running!\n"; 
    $runningPidHash=$updatePidHash;
    if ($runningPidNB>0){ print"             [       \$runningPidNB=$runningPidNB > 0 !! Sleep 1 seconds!\tCountinue to next Wait checking run!      ]\n\n"; sleep(1);   }   
    else                 {print"             [       \$runningPidNB=$runningPidNB <= 0 !! all child pid dead      ]\n\n";  sleep(1);last Markthree;   }
  }

  my $waitPid=wait(); my $time_2=`date`; chomp $time_2;  warn "             [       Start Wait AT:      ]\n             [       $time_1      ]\n\n             [       End of Wait, \$waitPid=$waitPid!!      ]\n\n             [       ", $time_2, "      ]\n\n" ;

	
	return $runningPidHash;
}


1;
