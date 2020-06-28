#!/usr/bin/perl -w
BEGIN{  push ( @INC, ('/home/fredjiang/fredProgs/PMs') )  };


use strict;
use warnings;
use POSIX ":sys_wait_h";
use DieWork;
use MultiThreadingWork;
use FastaFileHandle;
use DirFileHandle;
use Storable ;
use Storable qw(nstore store_fd nstore_fd freeze thaw dclone);
use Interproscan;

use Getopt::Std;
use vars qw($opt_a $opt_b $opt_i $opt_r $opt_s $opt_t $opt_d $opt_e $opt_o $opt_w $opt_h $opt_f $opt_c $opt_m $opt_p $opt_j);
getopts"a:b:i:r:s:t:d:e:o:w:h:f:c:m:p:j:";


my ($warnMsgBody, $warnMsgHead, $die_MsgHead, $caller_inform)=@{ DieWork::BuildWarnDieHeadInform( 'Mutilple_InterproScan.pl', 'Mutilple_InterproScan.pl' ) };

my $welcomeMsg="
-h the input hash holding all the path of fasta files and their msf files for interproscan
-c howmany core to use 
-o out put hash file
";

warn "\$opt_h=$opt_h\n" if (defined ($opt_h));
warn "\$opt_f=$opt_f\n" if (defined ($opt_f));
warn "\$opt_b=$opt_b\n" if (defined ($opt_b));
warn "\$opt_i=$opt_i\n" if (defined ($opt_i));
warn "\$opt_j=$opt_j\n" if (defined ($opt_j));
warn "\$opt_c=$opt_c\n" if (defined ($opt_c));
warn "\$opt_o=$opt_o\n" if (defined ($opt_o));
warn "\$opt_t=$opt_t\n" if (defined ($opt_t));
warn "\$opt_e=$opt_e\n" if (defined ($opt_e));
warn "\$opt_p=$opt_p\n" if (defined ($opt_p));
warn "\$opt_d=$opt_d\n" if (defined ($opt_d));


my $maxCoreLimit=60; #max core to use , nomatter the -opt_c is , not any cores more will be used than this number.
my $default_Core=4;  #default core to use, if the user didnot set core numbers.
my $user_Set_Core=MultiThreadingWork::HowManyCore_to_USE ( $maxCoreLimit, $default_Core, $opt_c );


my $inFstFl_hash_file=$opt_h;                                           DieWork::Check_FileDirExist_or_DIE( $opt_h,        "\$opt_h",        $die_MsgHead, $caller_inform  );  
my $inFstFl_hash=Storable::retrieve ( $inFstFl_hash_file );             DieWork::Check_Hash_or_DIE        ( $inFstFl_hash, "\$inFstFl_hash", $die_MsgHead, $caller_inform  );
  
my $out_HASH_FILE=$inFstFl_hash_file."_3m_.intproscan.hsh";             if ( DieWork::Check_DfdNoEmptString_or_NOT( $opt_o )  ){ $out_HASH_FILE=$opt_o;  }       

############################### 1             multiple thread ########## TOTAL HEAD ###########
############################### 1             multiple thread ########## TOTAL HEAD #################################
############################### 1             multiple thread ########## TOTAL HEAD ##############################################################


######################### part 1 HEAD   ###########################  multiple thread ###################
my $chdPidNb=0;
my $num_proc = 0;                     ## == number of proc ==
my $num_collect = 0;                  ## == number of collected ==
my $collect;
my $runningPidHash;

$SIG{CHLD} = \&proNBminus;            ## == get the child signal ==    
sub proNBminus{ 
	$num_proc--;                      	#print "\naaaaa -- \$num_proc=$num_proc\n"; warn "\naaaaa -- \$num_proc=$num_proc\n";
} 
######################### part 1 TAIL   ###########################  multiple thread ###################


######

#########

my $outChunk_hash;

my @FastaFilePathArray=(    sort { $a cmp $b } (   keys (  %{ $inFstFl_hash }  )   )    )  ;
my $howManyGroup=@FastaFilePathArray;
my $working_group_idx=0;                


my $temp_file_holding_all_fasta=TimeWork::GetTimeDirOrFileName;
$temp_file_holding_all_fasta=$out_HASH_FILE.".$temp_file_holding_all_fasta";
	  	
foreach my $eachFastaFile (    @FastaFilePathArray    ){ DieWork::Print_and_warn( "\n 20191216-0-0-1 \$eachFastaFile=$eachFastaFile\n"."\n" );
	DieWork::Check_FileDirExist_or_DIE( $eachFastaFile,        "\$eachFastaFile",        $die_MsgHead, $caller_inform  );  
  system ( "cat $eachFastaFile >> $temp_file_holding_all_fasta");
}
Interproscan::Doing_all_the_Ipsc_Search_MutipleThread_fromFASTAfile($temp_file_holding_all_fasta, $user_Set_Core);
system ( "rm -rf $temp_file_holding_all_fasta") ; 

foreach my $eachFastaFile (    @FastaFilePathArray    ){   DieWork::Print_and_warn( "\n 20191216-0-0-2 \$eachFastaFile=$eachFastaFile\n"."\n" );
  
  
  
  my $interProScanHashFile       =$eachFastaFile."_3m1m_itPrSc.hsh";  
  my $interProScanHash_read_File =$eachFastaFile."_3m2m_iPsRed.read.plx";  
  
  
  $outChunk_hash->{$eachFastaFile}->{'3m1m_itPrSc'}=$interProScanHashFile  ;
  $outChunk_hash->{$eachFastaFile}->{'3m2m_iPsRed'}=$interProScanHash_read_File  ;
    
  
  ######################### part 2 HEAD   ###########################  multiple thread ###################
  $chdPidNb++;
	my $pid=fork();    
	if (!defined($pid)) {
		  warn "\     Error in fork: $!";
      print "\     Error in fork: $!";
      exit 1;
  }
  ######################### part 2 MIDDLE #############  multiple thread ##########  
  if ($pid == 0) { ###### !!!!!!!#######
    ## == child proc ==           #print "\naaaaa Child $chdPidNb : My pid = $$\n";  
    warn "\nRound x Child   $chdPidNb : My pid = $$\n\t\t\t\t\t There are\$howManyGroup=$howManyGroup group, This is the \$working_group_idx=$working_group_idx one!!\n\n";
    ######################### part 2 TAIL   ###########################  multiple thread ###################
    	 
    
    
    my $interProScanHash=Interproscan::Eatract_interProscanHash_from_JL_ipsDatabase_inputFastaFile ($eachFastaFile);
	  DirFileHandle::PrintDumper_plReadFile($interProScanHashFile, $interProScanHash, $interProScanHash_read_File) ; 
    
    
    ######################### part 3 HEAD   ###########################  multiple thread ###################
    exit 0; 
  } ###### !!!!!!!#######       
  ######################### part 3 MIDDLE #############  multiple thread ########## 
  $runningPidHash->{$pid}=1;
	Marktwo: while(1) {
	  my $runningPidNB=0; 
	  my $updatePidHash;
	  foreach my $pidKey (    sort {$a<=>$b}(   keys(  %{ $runningPidHash }  )   )    ){
	  	if ( (waitpid($pidKey, WNOHANG)) > 0 ){my $tpVal=waitpid($pidKey, WNOHANG); }#  print "\naaaaaa\$pidKey=$pidKey is dead!!\t\$tpVal=$tpVal\n\n";}
	  	else {$updatePidHash->{$pidKey}=1; $runningPidNB++;                         }#  print "\naaaaaa\$runningPidNB child pid is runing!!\t$pidKey=$pidKey is still running!!! \n\n";}
	  }
	  warn "\n\n                                     [       $runningPidNB process running now !   ($working_group_idx/$howManyGroup)    ]      \n";   #print "\naaaaaa Totally \$runningPidNB=$runningPidNB is running!\n"; 
	  $runningPidHash=$updatePidHash;
	  if ($runningPidNB<$user_Set_Core){ warn"                   [       process ID: $pid\t $runningPidNB process running! <  We set $user_Set_Core CPU core to use !! Countinue to next fork!       ]\n\n";last Marktwo;    }
	  else                             { warn"                   [       process ID: $pid\t $runningPidNB process running! >= We set $user_Set_Core CPU core to use !! So Just Sleep 1 seconds!      ]\n\n";  sleep(1);   }
	}
  ######################### part 3 TAIL   ###########################  multiple thread ###################
  
  $working_group_idx++;
  
}       

######################### part 4 HEAD   ###########################  multiple thread ###################                                                                                                                                                                                                                    
my $time_1=`date`; chomp $time_1;                                                                                                                                                                                                                                                                                           
warn "\n\n             [       Start to wait !!      ]\n\n             [       ", $time_1, "      ]\n\n" ;                                                                                                                                                                                                                  
Markthree: while(1) {                                                                                                                                                                                                                                                                                                       
  my $runningPidNB=0;                                                                                                                                                                                                                                                                                                       
  my $updatePidHash;                                                                                                                                                                                                                                                                                                        
  foreach my $pidKey (    sort {$a<=>$b}(   keys(  %{ $runningPidHash }  )   )    ){                                                                                                                                                                                                                                        
  	if ( (waitpid($pidKey, WNOHANG)) > 0 ){my $tpVal=waitpid($pidKey, WNOHANG); my $ps=`ps $pidKey`; print "\n             [       \$pidKey=$pidKey is dead!!\t\$tpVal=$tpVal\t\$ps=$ps      ]\n\n"; warn "\na             [       \$pidKey=$pidKey is dead!!\t\$tpVal=$tpVal\n\$ps=$ps      ]\n\n";}                       
  	else {                                                                                                                                                                                                                                                                                                                  
  		my $ps=`ps $pidKey`;                                                                                                                                                                                                                                                                                                  
  		if ($ps=~m/$pidKey/){                                                                                                                                                                                                                                                                                                 
  		  $updatePidHash->{$pidKey}=1; $runningPidNB++;                                                                                                                                                                                                                                                                       
  		  #print "\n             [       $runningPidNB child pid is runing!!\t",(waitpid($pidKey, WNOHANG)),"\t\$pidKey=$pidKey is still running!!!       ]\n\$ps=$ps\n\n";                                                                                                                                                    
  		  warn  "\n             [       $runningPidNB child pid is runing!!\t",(waitpid($pidKey, WNOHANG)),"\t\$pidKey=$pidKey is still running!!!       ]\t\$ps=$ps\n\n";                                                                                                                                                    
  		}
  		else {
  			warn  "\n             [       Running Pid \$pidKey=$pidKey found\tWAITPID=",(waitpid($pidKey, WNOHANG)),"\t, but no ps was found!!\$ps=$ps      ]\n\n"; 
  			#print "\n             [       Running Pid \$pidKey=$pidKey found\tWAITPID=",(waitpid($pidKey, WNOHANG)),"\t, but no ps was found!!\$ps=$ps      ]\n\n";
  		}
  	}                                                                                                                                                                                                                                                                                                                       
  }                                                                                                                                                                                                                                                                                                                         
                                                                                                                                                                                                                                                                                                                            
  warn "\n             [       Totally  \$runningPidNB=$runningPidNB      ]\n";   print "\n             [       Totally \$runningPidNB=$runningPidNB is running!\n";                                                                                                                                                        
  $runningPidHash=$updatePidHash;                                                                                                                                                                                                                                                                                           
  if ($runningPidNB>0){ print"             [       \$runningPidNB=$runningPidNB > 0 !! Sleep 1 seconds!\tCountinue to next Wait checking run!      ]\n\n"; sleep(1);   }                                                                                                                                                    
  else                 {print"             [       \$runningPidNB=$runningPidNB <= 0 !! all child pid dead      ]\n\n";  sleep(1);last Markthree;   }                                                                                                                                                                       
}                                                                                                                                                                                                                                                                                                                           
                                                                                                                                                                                                                                                                                                                            
my $waitPid=wait();                                                                                                                                                                                                                                                                                                         
my $time_2=`date`; chomp $time_2;  warn "             [       Start Wait AT:      ]\n             [       $time_1      ]\n\n             [       End of Wait, \$waitPid=$waitPid!!      ]\n\n             [       ", $time_2, "      ]\n\n" ;                                                                               
######################### part 4 TAIL   ###########################  multiple thread ###################                                                                                                                                                                                                                    
                            
DirFileHandle::PrintDumper_plReadFile ($out_HASH_FILE, $outChunk_hash);
