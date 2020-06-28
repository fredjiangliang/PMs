#!/usr/bin/perl -w
BEGIN{  push (@INC, 'home/fredjiang/fredProgs/PMs')  };


use strict;
use warnings;



############################################################################################################

#         
############################################################################################################



                      
package PrintSubArrayHash;





#sub999999999999999          print_all_sub_array {    #打印多维数组的内部情况
sub  print_all_sub_array {    #打印多维数组的内部情况
       my ($temp_address, $temp_level, $array_nb, $print_line )=@_;  #print "$temp_level, $array_nb, $print_line\n";
       if ( ( $temp_level eq '') && ( $array_nb eq '')  && ( $print_line eq '') ){
                ($temp_level, $array_nb, $print_line ) =(0, 0, "\$ref");                #print  "kkkkkkkkkkk  $temp_level, $array_nb, $print_line\n";
       }
       else {                   $print_line="$print_line\->[$array_nb]";       }
       my $next_temp_level=($temp_level+1);
       if ($temp_level >5 ){
       	  for (my $j=0; $j<=$temp_level; $j++){                        printf "%5s\t", "|" ;               }               printf "%-5s\t%-30s\t%-30s\n", "$temp_level>", "$print_line=","$temp_address";
       }
       else {
           if (   $temp_address =~m/HASH/ ){                 
                    for (my $j=0; $j<=$temp_level; $j++){                        printf "%5s\t", "|" ;               }               printf "%-5s\t%-5s\t%-30s\n", "$temp_level>", "[$array_nb]", "$temp_address";
                    my @temp_array=%$temp_address;
                #for (my $i=0; $i<@temp_array; $i++) {    print "$temp_array[$i] ne  $temp_address\n"; if ( $temp_array[$i] ne $temp_address ) {print_all_sub_array ($temp_array[$i] , $next_temp_level, $i, $print_line );}        }
                for (my $i=0; $i<@temp_array; $i++) {     if ( $temp_array[$i] ne $temp_address ) {print_all_sub_array ($temp_array[$i] , $next_temp_level, $i, $print_line );}        }
           }       
           elsif (  $temp_address =~m/ARRAY/ ){
                    for (my $j=0; $j<=$temp_level; $j++){                        printf "%5s\t", "|" ;               }               printf "%-5s\t%-5s\t%-30s\n", "$temp_level>", "[$array_nb]", "$temp_address";
                    my @temp_array=@$temp_address;
                #for (my $i=0; $i<@temp_array; $i++) {    print "$temp_array[$i] ne  $temp_address\n"; if ( $temp_array[$i] ne $temp_address ) {print_all_sub_array ( $temp_array[$i] , $next_temp_level, $i, $print_line );}       }
                for (my $i=0; $i<@temp_array; $i++) {     if ( $temp_array[$i] ne $temp_address ) {print_all_sub_array ( $temp_array[$i] , $next_temp_level, $i, $print_line );}       }
           }     
           else {                
                    for (my $j=0; $j<=$temp_level; $j++){                        printf "%5s\t", "|" ;               }               printf "%-5s\t%-30s\t%-30s\n", "$temp_level>", "$print_line=","$temp_address";
           }
        }
}  


sub change_array_into_a_printStyle_string{
	my ($inArray)=@_;
	my $msgHead="\n\n\nNow In package PrintSubArrayHash,\nIn Sub change_array_into_a_printStyle_string\n\n";	
	my $outString; 
	my $lineLen; my $HuanHangLimit=10000;
	if ( ref($inArray) eq 'ARRAY' ){
		for (my $i=0; $i<@{ $inArray }; $i++){
			my $addString="[$i]$inArray->[$i]";
			my $addLength=length $addString;
			$lineLen+=$addLength+1;
			if ($lineLen>$HuanHangLimit){
				$outString.="$addString\n";
				$lineLen=0;
			}
			else{
				$outString.="$addString\t";
			}
		}
	}
	else {
		die "$msgHead\n\$inArray=$inArray IS NOT A ARRAY REF!!!!\n\n\n\n";
	}
	return $outString;
}

sub change_hash_into_a_printStyle_string{
	my ($inHash)=@_;
	my $msgHead="\n\n\nNow In package PrintSubArrayHash,\nIn Sub change_hash_into_a_printStyle_string\n\n";;	
	my $outString; 
	my $lineLen; my $HuanHangLimit=10000;
	if ( ref($inHash) eq 'HASH' ){
		foreach my $ecKy (  @{ DirFileHandle::AutoSortKey( $inHash ) }  ){ #print "\$ecKy=$ecKy\n";
			#my $addString="{$ecKy}$inHash->{$ecKy}"; print "{$ecKy}$inHash->{$ecKy}\n";
			my $addString="$inHash->{$ecKy}";
			my $addLength=length $addString;
			$lineLen+=$addLength+1;
			if ($lineLen>$HuanHangLimit){
				$outString.="$addString\n";
				$lineLen=0;
			}
			else{
				$outString.="$addString\t";
			}
		}
	}
	else {
		die "$msgHead\n\$inHash=$inHash IS NOT A Hash REF!!!!\n\n\n\n";
	}
	return $outString;
}


sub  print_all_sub_array_into_string {    #打印多维数组的内部情况
       my ($temp_address,  $temp_level, $array_nb, $print_line )=@_;  #print "$temp_level, $array_nb, $print_line\n";
       
       my $out_string;
       
       if ( ( $temp_level eq '') && ( $array_nb eq '')  && ( $print_line eq '') ){
                ($temp_level, $array_nb, $print_line ) =(0, 0, "\$ref");                #print  "kkkkkkkkkkk  $temp_level, $array_nb, $print_line\n";
       }
       else {                   $print_line="$print_line\->[$array_nb]";       }
       my $next_temp_level=($temp_level+1);
       if ($temp_level >5 ){
       	  for (my $j=0; $j<=$temp_level; $j++){                     $out_string=    sprintf "%5s\t", "|" ;               }               $out_string=    sprintf "%-5s\t%-30s\t%-30s\n", "$temp_level>", "$print_line=","$temp_address";
       }
       else {
           if (   $temp_address =~m/HASH/ ){                 
                    for (my $j=0; $j<=$temp_level; $j++){                        $out_string=    sprintf "%5s\t", "|" ;               }               $out_string=    sprintf "%-5s\t%-5s\t%-30s\n", "$temp_level>", "[$array_nb]", "$temp_address";
                    my @temp_array=%$temp_address;
                #for (my $i=0; $i<@temp_array; $i++) {    print "$temp_array[$i] ne  $temp_address\n"; if ( $temp_array[$i] ne $temp_address ) {print_all_sub_array ($temp_array[$i] , $next_temp_level, $i, $print_line );}        }
                for (my $i=0; $i<@temp_array; $i++) {     if ( $temp_array[$i] ne $temp_address ) {print_all_sub_array ($temp_array[$i] , $next_temp_level, $i, $print_line );}        }
           }       
           elsif (  $temp_address =~m/ARRAY/ ){
                    for (my $j=0; $j<=$temp_level; $j++){                        $out_string=    sprintf "%5s\t", "|" ;               }               $out_string=    sprintf "%-5s\t%-5s\t%-30s\n", "$temp_level>", "[$array_nb]", "$temp_address";
                    my @temp_array=@$temp_address;
                #for (my $i=0; $i<@temp_array; $i++) {    print "$temp_array[$i] ne  $temp_address\n"; if ( $temp_array[$i] ne $temp_address ) {print_all_sub_array ( $temp_array[$i] , $next_temp_level, $i, $print_line );}       }
                for (my $i=0; $i<@temp_array; $i++) {     if ( $temp_array[$i] ne $temp_address ) {print_all_sub_array ( $temp_array[$i] , $next_temp_level, $i, $print_line );}       }
           }     
           else {                
                    for (my $j=0; $j<=$temp_level; $j++){                        $out_string=    sprintf "%5s\t", "|" ;               }               $out_string=    sprintf "%-5s\t%-30s\t%-30s\n", "$temp_level>", "$print_line=","$temp_address";
           }
       }
       return $out_string;
}  



1;