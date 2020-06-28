#!/usr/bin/perl -w
BEGIN{  push (@INC, '/home/fredjiang/fredProgs/PMs')  };

use strict;
use warnings;

use DieWork;
use InFileHandle;
use ClustalwRun;
use File::Basename;
use File::Spec;


                      
package LGP_run_mrbayes;


##################################################################################################
#
#   
#
#
##################################################################################################

my $mb_g_bigEstNb=50000000;
my $run_mrbayes_sh="/home/fredjiang/OneTSSD/tools20191105/run_mrbayes20191202/run_mrbayes.sh";
my $modify_nexus_aln_py="/home/fredjiang/OneTSSD/tools20191105/run_mrbayes20191202/modify_nexus_aln.py";
my $modify_nexus_tree_py="/home/fredjiang/OneTSSD/tools20191105/run_mrbayes20191202/modify_nexus_tree.py";
my $run_mrbayes_parallel_sh="/home/fredjiang/OneTSSD/tools20191105/run_mrbayes20191202/run_mrbayes_parallel.sh";
my $runSH_byJL20191203_pl="/home/fredjiang/OneTSSD/tools20191105/run_mrbayes20191202/runSH_byJL20191203.pl";

sub Run_LGP_MRbayes{  #my $outNEXfilePATH=LGP_run_mrbayes::Run_LGP_MRbayes($msfFILE, $out_dir);
	my ($msfFILE, $out_dir)=@_;
	
	my ($warnMsgBody, $warnMsgHead, $die_MsgHead, $caller_inform)=@{ DieWork::BuildWarnDieHeadInform( 'LGP_run_mrbayes', 'Run_LGP_MRbayes' ) };
	
	DieWork::Check_FileDirExist_or_DIE( $msfFILE, "\$msfFILE", $die_MsgHead, $caller_inform  );
  
  my $msfFILE_DirNAME=File::Basename::dirname  ($msfFILE);
  
  my $real_out_basedir=$msfFILE."_mrbys.dir";  if (  DieWork::Check_DfdNoEmptString_or_NOT ( $out_dir )  ) { $real_out_basedir=$out_dir;   }
  $real_out_basedir=File::Basename::basename  ($real_out_basedir);
  
  my $real_out_dir=$msfFILE_DirNAME."/".$real_out_basedir;
  
  my $mkdirCMD="mkdir -p $real_out_dir"; DieWork::Print_and_warn( $mkdirCMD."\n" ); 
  if (  -d ($real_out_dir)  ) {} else {    system ( "$mkdirCMD" ) ; }
  
  
  my $msfFILE_baseNAME='m.msf';  #               File::Basename::basename ( $msfFILE ); 
  my $realWorkingMSF="$real_out_dir/$msfFILE_baseNAME";
  
  #  
  
#  cp -f $run_mrbayes_sh          $real_out_dir/
#  cp -f $modify_nexus_aln_py     $real_out_dir/
#  cp -f $modify_nexus_tree_py    $real_out_dir/
#  cp -f $run_mrbayes_parallel_sh $real_out_dir/  
#  cp -f $runSH_byJL20191203_pl   $real_out_dir/  
my $cpCMD="  
    cp -f $msfFILE $realWorkingMSF
  ";
  DieWork::Print_and_warn( $cpCMD."\n" );  system ( "$cpCMD" ) ;
  
	my $Seq_Number=ClustalwRun::GetSeq_Number_in_aln_file( $msfFILE, 'msf' );     	
	
	
	
	
	if (   (  defined ( $Seq_Number )  ) && (  $Seq_Number>0 )   ){
    my $gNmuber=int ( $mb_g_bigEstNb/$Seq_Number );     	                                           	                    	                                                     
  
    #my $perlCMD="perl $real_out_dir/runSH_byJL20191203.pl -d $real_out_dir -g 100 "; #$gNmuber ";
    #DieWork::Print_and_warn( $perlCMD."\n" );    system ( "$perlCMD");
    
    
    
    #my $workingDIR = `pwd`;
    #my $cd_CMD_1="cd $real_out_dir";  DieWork::Print_and_warn( $cd_CMD_1."\n" );     `$cd_CMD_1`; #system ( "$cd_CMD_1");
    
    my $rltv_wkMsfPath=File::Spec->abs2rel( $realWorkingMSF );
    my $mbTreeCmd="bash $run_mrbayes_sh -f $rltv_wkMsfPath -g $gNmuber ";  
    #my $mbTreeCmd="bash run_mrbayes.sh -f $realWorkingMSF -g 100 "; #$gNmuber ";  
    DieWork::Print_and_warn( $mbTreeCmd."\n" );     system ( "$mbTreeCmd");
    
    #my $cd_CMD_2="cd $workingDIR";    DieWork::Print_and_warn( $cd_CMD_2."\n" );     `$cd_CMD_2`; #system ( "$cd_CMD_2");
    
  }
	
	my $outNEXfilePATH=$realWorkingMSF;  
	$outNEXfilePATH=~s/msf$/nex.con_tre.nex/;
	
	return $outNEXfilePATH;
}

sub Not_Run_JUST_gettheNEXpath{  #my $outNEXfilePATH=LGP_run_mrbayes::Not_Run_JUST_gettheNEXpath($msfFILE, $out_dir);
	my ($msfFILE, $out_dir)=@_;
	
	my ($warnMsgBody, $warnMsgHead, $die_MsgHead, $caller_inform)=@{ DieWork::BuildWarnDieHeadInform( 'LGP_run_mrbayes', 'Not_Run_JUST_gettheNEXpath' ) };
	
	DieWork::Check_FileDirExist_or_DIE( $msfFILE, "\$msfFILE", $die_MsgHead, $caller_inform  );
  
  my $msfFILE_DirNAME=File::Basename::dirname  ($msfFILE);
  
  my $real_out_basedir=$msfFILE."_mrbys.dir";  if (  DieWork::Check_DfdNoEmptString_or_NOT ( $out_dir )  ) { $real_out_basedir=$out_dir;   }
  $real_out_basedir=File::Basename::basename  ($real_out_basedir);
  my $real_out_dir=$msfFILE_DirNAME."/".$real_out_basedir;
  #my $mkdirCMD="mkdir -p $real_out_dir";
  #if (  -d ($real_out_dir)  ) {} else {   DieWork::Print_and_warn( $mkdirCMD."\n" );  system ( "$mkdirCMD" ) ; }
  
  
  #my $msfFILE_baseNAME=File::Basename::basename ( $msfFILE ); 
  my $msfFILE_baseNAME='m.msf';
  my $realWorkingMSF="$real_out_dir/$msfFILE_baseNAME";
  my $outNEXfilePATH=$realWorkingMSF;
	
	$outNEXfilePATH=~s/msf$/nex.con_tre.nex/;
	
	return $outNEXfilePATH;
}




1;