#!/usr/bin/perl

use strict;
use warnings;

# Author: Stefan Seidel (date:15-mar-2013)

# this program parses 3 strings from a file(multi fasta format) what is 
# given as 2nd command line argument and applies the
# Needleman Wunsch algoritm in "3D" to them.
# 
# file format: 
#	>seq 
#	AATTTGGGCCTTTGGGT
#	>seq
#	AATTGGCCTTGGTTTGG
#	>seq
#	AAATTCCCTTTTGGCCT
#
#
# for a global alignment. The sequences are printed
# to the screen.

# however the program has problems with long and "equal long" sequences


# costs:
# match = 2
# mismatch = -1
# gap = -4

my $fh_in;      # filehandle
my $filename;	# filename
$filename = $ARGV[0];
chomp $filename;


open $fh_in, "$filename" or die "ERROR: could not open file '$filename': $!";

my $r_array = [];   # right matrix (1st dimension)
my @r_array;
 push(@r_array, 'space');

my $d_array = [];   # down matrix (2nd dimension)
my @d_array;
 push(@d_array, 'space');

my $u_array = [];    # up matrix (3rd dimension)
my @u_array;
 push(@u_array, 'space');

my @temp; # temporary array used for sorting of length of the sequences

my $dia_array = []; # diagonal matrix

my $var1; # variables taking sequences as input from file
my $var2; # variables taking sequences as input from file
my $var3; # variables taking sequences as input from file




my @temp_array;		# temporary array if r_array longer than d_array

# scoring for match/missmatch insertion or deletion for the scoring matrices
my $match = 2;
my $mismatch = -1;
my $gap = -4;	# no string match
#my $gap2 = -2;	# 2 strings match
my $matchquality; # stores eiter match or missmatch in scoring matrix

# case variables for tracking
my $ac = '';
my $bc = '';
my $cc = '';
my $dc = '';
my $ec = '';
my $fc = '';
my $gc = '';
my $hc = '';

#my $result_matrix = [];


my $track = ''; 	#tracks the moves for the path back from the bottom right

 # if set to 1 2 or 3 indicates that one sequence string (1,2,3) is over and filled up with gaps;
my $gapswitch1 = 0;
my $gapswitch2 = 0;
my $gapswitch3 = 0;

# sequences for the alignment output
my $rev_seq1 = '';	# reverse aligned sequence1
my $rev_seq2 = '';	# reverse aligned sequence2
my $rev_seq3 = '';      # reverse aligned sequence3

my $seq1;	# aligned sequence1
my $seq2;	# aligned sequence2
my $seq3;	# aligned sequence3



my $count;	  # matrixcounter for 2nd dimension
my $bool_switch;  # switch for down/right/diagonal weighted edge array matrix
my $seq_switch = 0;   # switch for 1st 2nd or 3rd sequence

while (<$fh_in>)
	{ #while
	if ( ($_ =~ /^>/) && ($seq_switch == 0))	# selection for down weight matrix (1st sequence)	
		{
		$bool_switch = 0;
		$count = 0;
		$seq_switch = 1;
		next;
		}#if

	elsif (($_ =~ /^>/) && ($seq_switch == 1))	# selection for right weight matrix (2nd sequence)
		{
		$bool_switch = 1;
		$count = 0;
		$seq_switch = 2;
		next;
		}#elsif

	elsif (($_ =~ /^>/) && ($seq_switch == 2))	# selection for up weight matrix (3nd sequence)
		{
		$bool_switch = 2;
		$count = 0;
		next;
		}#elsif
####################################################################


if ($bool_switch == 0 )
	{#if
	if ( ($_ =~ /^>/) && ($seq_switch == 0))	# selection for down weight matrix (1st sequence)	
		{#if
		next;
		}#if end
	else
	{#else
		$var1.=$_; 
		chomp($var1);
#		print "$_"; # debugging purpose
		@r_array = split(//m, $var1);

	#debugging
	#	print "$r_array[$count] [scalar @{$r_array}] \n";
	#	print "before for\n";
	#	for ($i = 0 ; $i < scalar @r_array ; $i++)
	#	{#for
	#	print "$r_array[$i] \n";
	#	}# end for
	#	print "after for\n";
		$count++;	
		#referencing array
#		$d_array = \@d_array;
	} #else
	}#if
	

if ($bool_switch == 1 )
	{#if
	if ( ($_ =~ /^>/) && ($seq_switch == 1))	# selection for right weight matrix (2nd sequence)	
		{#if
		next;
		}#if end
	else
		{#else
		$var2.=$_;
		chomp($var2); 

		# splitting up string into d_array
		@d_array = split(//m, $var2);

		$count++;
		}# else



	}# if

if ($bool_switch == 2 )
	{#if
	if ( ($_ =~ /^>/) && ($seq_switch == 2))	# selection for up weight matrix (3nd sequence)	
		{#if
		next;
		}#if end
	else
		{#else
		$var3.=$_;
		chomp($var3); 
#		print "$_"; # debugging purpose

		# splitting up string into d_array
		@u_array = split(//m, $var3);


		}#else end

	}# if end
	}# while end


close $fh_in;


# assigning length of arrays to variables for sorting the sequences in descending length
my $d = scalar @d_array; # length of array in $d_array

my $r = scalar @r_array; # length of array in $r_array

my $u = scalar @u_array; # length of array in $u_array

# sort sequences in descending lengths
#=comment
if (( $u > $d) and ($u > $r))
	{
	@temp = @r_array;
	@r_array = @u_array;
	@u_array = @temp;

	# reassigning length of arrays to variables
	$r = scalar @r_array ;
	$u = scalar @u_array ;
	}
if(( $d > $u ) and ($d > $r))
	{
	@temp = @r_array;
	@r_array = @d_array;
	@d_array = @temp;

	# reassigning length of arrays to variables
	$r = scalar @r_array ;
	$d = scalar @d_array ;
	}
	
if (( $u > $d )) 
	{
	@temp = @d_array;
	@d_array = @u_array;
	@u_array = @temp;
	
	# reassigning length of arrays to variables
	$d = scalar @d_array ;
	$u = scalar @u_array ;
	}


# assigning length of arrays for scoring matrix
 $d = scalar @d_array; # length of array in $d_array

 $r = scalar @r_array; # length of array in $r_array

 $u = scalar @u_array; # length of array in $u_array


	# referencing of arrays
		$r_array = \@r_array;
		$d_array = \@d_array;
		$u_array = \@u_array;
=end comment
=cut


#print "\n length of sequence#1:$r \n length of sequence#2:$d\n length of sequence#3:$u\n";


my $result_matrix = [];

# assingning value to starting point of the 3D matrix to zero
$result_matrix -> [0][0][0] = 0;




#######################################################################################
#######################################################################################
#######################################################################################

# assigning values to first colum (adding up numbers)


for (my $j = 1; $j < ($d+1); $j++) #d_array
	{
			
#		$d_array->[0][$j][0] = $gap;
		
	$result_matrix->[0][$j][0] = $result_matrix->[0][$j-1][0] + $gap; #$d_array->[0][$j-1][0];
	#}#for end
###############################################################################

# assigning 0 to 0,0,0
#$result_matrix->[0][0][0]=0;

# debugging
#print "resultarray d(j=$j): $result_matrix->[0][$j][0]\n";  # line for debugging purposes
	}# for end


# assigning values to the first row (adding up numbers)
for (my $i = 1; $i < ($r+1); $i++) #r_array
	{
#		$r_array->[$i][0][0] = $gap;

	# generation of first row by adding up the numbers
	$result_matrix->[$i][0][0] = $result_matrix->[$i-1][0][0] + $gap; #$r_array->[$i-1][0][0];
#	}#for end

# debugging
#print "resultarray d(i=$i): $result_matrix->[$i][0][0]\n"; # line for debugging purposes
	}#for end

# assigning values to the first row (adding up numbers)

for (my $k = 1; $k < ($u+1); $k++) #r_array
	{
#		$u_array->[0][0][$k] = $gap;

	# generation of first row by adding up the numbers
	$result_matrix->[0][0][$k] = $result_matrix->[0][0][$k-1] + $gap; #$u_array->[0][0][$k-1];
#	}#for end

# debugging
#print "resultarray u(k=$k) : $result_matrix->[0][0][$k]\n"; # line for debugging purposes
	}#for end


######################################################################################################
#################### generation of scoring matrix #################################################### 
######################################################################################################
 
# initialising the matrix with '0'
for (my $i = 0; $i < ($r+1); $i++)
	{# for $i
	for ( my $j = 0; $j < ($d+1); $j++)
		{ #for $j 
		for ( my $k = 0; $k < ($u+1); $k++)
			{
			$result_matrix->[$i][$j][$k] = 0;
			} # for end $k
		} # for end $j
	}# for end $i


# initialising matrix right down

for (my $i = 1; $i < ($r+1); $i++)
	{# for $i
	for ( my $j = 1; $j < ($d+1); $j++)
		{ #for $j 
		
		my @value = ();
		$value[0] = $result_matrix->[$i][$j-1][0] + $gap ; 	#right weighted value: gap
		$value[1] = $result_matrix->[$i-1][$j][0] + $gap ; 	#right down weighted value: gap
=comment
		print "value1-gap: $result_matrix->[$i-1][$j][0]\n";
		print " i = $i , j = $j   i:1 to r+1   j:1 to d+1  \n";
		print "r_array[$i-1] = $r_array[$i-1]\n";
		print "d_array[$j-1] = $d_array[$j-1]\n";
=comment
=cut

		# literal evaluation if match or mismatch in the sequence is given 	
		if (( $r_array->[$i-1] eq $d_array->[$j-1] ) )
			{#if
			$value[2] = $result_matrix->[$i-1][$j-1][0] + $match ; 	#diagonal weighted value: match
			}#end if
		else
			{#else
			$value[2] = $result_matrix->[$i-1][$j-1][0] + $mismatch ; 	#diagonal weighted value: missmatch
			}# end else		

		# selection of highest value what is then added to matrix
		@value = sort { $b<=>$a } @value;
		$result_matrix->[$i][$j][0] = $value[0];
	}#for end $j

	} # for end $i


# initialising the matrix down up
for (my $j = 1; $j < ($d+1); $j++)
	{# for $i
	for ( my $k = 1; $k < ($u+1); $k++)
		{ #for $j 
		
		my @value = ();
		$value[0] = $result_matrix->[0][$j-1][$k] + $gap ; 	#up weighted value: gap
		$value[1] = $result_matrix->[0][$j][$k-1] + $gap ; 	#down weighted value: gap

		# literal evaluation if match or mismatch in the sequence is given 	
		if (( $d_array->[$j-1] eq $u_array->[$k-1] ) )
			{#if
			$value[2] = $result_matrix->[0][$j-1][$k-1] + $match ; 	#diagonal weighted value: match
			}#end if
		else
			{#else
			$value[2] = $result_matrix->[0][$j-1][$k-1] + $mismatch ; 	#diagonal weighted value: missmatch
			}# end else		

		# selection of highest value what is then added to matrix
		@value = sort { $b<=>$a } @value;
		$result_matrix->[0][$j][$k] = $value[0];
	}#for end $k

	} # for end $j


# initialsing the matrix right up:
for (my $i = 1; $i < ($r+1); $i++)
	{# for $i
	for ( my $k = 1; $k < ($u+1); $k++)
		{ #for $j 
		
		my @value = ();
		$value[0] = $result_matrix->[$i-1][0][$k] + $gap ; 	#up weighted value: gap
		$value[1] = $result_matrix->[$i][0][$k-1] + $gap ; 	#right weighted value: gap

		# literal evaluation if match or mismatch in the sequence is given 	
		if (( $r_array->[$i-1] eq $u_array->[$k-1] ) )
			{#if
			$value[2] = $result_matrix->[$i-1][0][$k-1] + $match ; 	#diagonal weighted value: match
			}#end if
		else
			{#else
			$value[2] = $result_matrix->[$i-1][0][$k-1] + $mismatch ; 	#diagonal weighted value: missmatch
			}# end else		

		# selection of highest value what is then added to matrix
		@value = sort { $b<=>$a } @value;
		$result_matrix->[$i][0][$k] = $value[0];
		}#for end $j
	}# for end $k	


#=comment
# assigning values to the rest of the scoring matrix (resultmatrix)

for (my $i = 1; $i < ($r+1); $i++)
	{#for $i
	for (my $j = 1; $j < ($d+1); $j++)
		{ #for $j
		for (my $k = 1; $k < ($u+1); $k++)
		
	{ #for $k
#		if ($i < 1 ) {$i = 0;}	
#		if ($j < 1 ) {$j = 0;}	
#		if ($k < 1 ) {$k = 0;}	
		
		my @value = ();
		$value[0] = $result_matrix->[$i][$j-1][$k-1] + $gap ; 	#right weighted value: gap
		$value[1] = $result_matrix->[$i][$j][$k-1] + $gap ; 	#right down weighted value: gap
		$value[2] = $result_matrix->[$i-1][$j-1][$k] + $gap ; 	#up weighted value: gap
		$value[3] = $result_matrix->[$i-1][$j][$k] + $gap ; 	#up down weighted value: gap
		$value[4] = $result_matrix->[$i][$j-1][$k] + $gap ; 	#right up weighted value: gap
		$value[5] = $result_matrix->[$i-1][$j][$k-1] + $gap ; 	#down weighted value: gap
	#	print "value5-gap = $result_matrix->[$i-1][$j][$k-1]\n";
	#	print "i= $i j=$j k=$k\n";

		# literal evaluation if match or mismatch in the sequence is given 
		# "$i-1" because $r_array starts with "0" 	
		if (( $r_array->[$i-1] eq $d_array->[$j-1] )
		 and ($r_array->[$i-1] eq $u_array->[$k-1])
		 and ($d_array->[$j-1] eq $u_array->[$k-1]))
			{#if
			$value[6] = $result_matrix->[$i-1][$j-1][$k-1] + $match ; 	#diagonal weighted value: match
			}#end if
		else
			{#else
			$value[6] = $result_matrix->[$i-1][$j-1][$k-1] + $mismatch ; 	#diagonal weighted value: missmatch
			}# end else		

		# selection of highest value what is then added to matrix
		@value = sort { $b<=>$a } @value;
		$result_matrix->[$i][$j][$k] = $value[0];
		
#	print "my value[$i][$j][$k] = @value[0]\n";	
			}# for end $k
		}#for $j end	
	}# for $i end

#print "end of assigning values to scoring matrix\n";

=printing for debugging
# printing resultmatrix to screen

for (my $i = 0; $i < ($r+1); $i++) 
	{#for
for (my $j = 0; $j < ($d+1); $j++)
	{
	for (my $k = 0; $k < ($u+1); $k++)
		{#for
		print "$result_matrix->[$i][$j][$k] ";
		}# for end
	print "\n";
	}# for end
}# for end
=comment
=cut



###############################################################################
# for backtracking ################## for backtracking ########################
###############################################################################


# setting the starting point coordinates for the back tracking
my $i = ($r+1);
my $j = ($d+1);
my $k = ($u+1);

while (($i > 0) or ($j > 0) or ($k > 0))
{ # while

	if (($i < 1))  { $gapswitch1 = 1;} 
	if (($j < 1))  { $gapswitch2 = 2;}
	if (($k < 1))  { $gapswitch3 = 3;} 
	
	#print "i = $i , j = $j , k = $k\n";
		# avoiding that $i $j and $k are "escaping"(negative values) the matrix 

#=comment
			
			if ( $k < 1) 
				{ 
				$k = 0;
#				print "in kk\n";
				}

			if ( $j < 1) 
				{ 
				$j = 0;
#				print "in jj\n";
				}
		#	if (($i < 1 ) && ( $j < 1 ) && ($k < 1)) { print "stop\n";}

=end comment
=cut
		# if not defined allocate -100000 as value
		my $not_def_insert = "-10000";
		
		my $max_score = $result_matrix->[$r][$d][$u];
		my $not_def_match = $max_score-1;


	# filling up "reverse aligned sequence1" and "reverse aligned sequence2"
	$ac = $result_matrix->[$i-1][$j-1][$k-1]; # match
	#	print "ac = $ac \n";

		# insert right ( = sequence1)
		$bc = $result_matrix->[$i][$j-1][$k-1];
		# insert down ( = sequence2)
		$cc = $result_matrix->[$i-1][$j][$k-1];
		# insert up ( = sequence3)
		$dc = $result_matrix->[$i-1][$j-1][$k];
		# insert right and down ( = sequence1 and 2)
		$ec = $result_matrix->[$i][$j][$k-1];
		# insert right and up ( = sequence1 and 3)
		$fc = $result_matrix->[$i][$j-1][$k];
		# insert down and up ( = sequence2 and 3)
		$gc = $result_matrix->[$i-1][$j][$k];
		# insert right and up ( = sequence1 and 3)
		$hc = $result_matrix->[$i][$j-1][$k];
		
		$ac //= $not_def_match;
		$bc //= $not_def_insert; 
		$cc //= $not_def_insert;
		$dc //= $not_def_insert; 
		$ec //= $not_def_insert; 
		$fc //= $not_def_insert; 
		$gc //= $not_def_insert; 
		$hc //= $not_def_insert; 

#		print "ac = $ac \n";
#		print "bc = $bc \n";
#		print "cc = $cc \n";
#		print "dc = $dc \n";
#		print "ec = $ec \n";
#		print "fc = $fc \n";
#		print "gc = $gc \n";

		# match/missmatch showed highest score in the matrix
		if (($ac ge $bc) or ($ac ge $cc) or ( $ac ge $dc) or ($ac ge $ec) or ($ac ge $fc) or ($ac ge $gc) or ($ac ge $hc))
#		if ($ac ge ( $bc or $cc or $dc or $ec or $fc or $gc or $hc))
			{# if
			if ($gapswitch1 == 1)			
				{
				$rev_seq1.='-';
				}
			else
			{#else
			if (defined $r_array[$i-1]){              # if array element not defined, don't try to append to sequence
				$rev_seq1.=$r_array[$i-1];
		#	print "$r_array[$i-1]\n";
				}

			# moving one step further in the matrix ($i-1) and ($j-1)
			$i=$i-1;

			if ($gapswitch2 == 2)
				{
				$rev_seq2.='-';
				}
			else
			{
			if (defined $d_array[$j-1]){	             # if array element not defined, don't try to append to sequence
			$rev_seq2.=$d_array[$j-1];
				}
			$j=$j-1;
		#	print "m i= $i, j=$j resultmatrix:$result_matrix->[$i+1][$j+1]\n";
			}#else

			if ($gapswitch3 == 3)
				{
				$rev_seq3.='-'
				}
			elsif ($gapswitch3 != 3)
			{
			if (defined $u_array[$k-1]){                       # if array element not defined, don't try to append to sequence
				$rev_seq3.=$u_array[$k-1];
				}
			$k=$k-1;
			$track.='m-';
			#print "in ac if\n";
			# leaving for loop
			}
			}#else
			}# if end

		# insert right ( =sequence1)
		if (($bc gt $ac) and ($bc gt $cc) and ( $bc gt $dc) and ($bc gt $ec) and ($bc gt $fc) and ($bc gt $gc) and ($bc gt $hc))
			{# elsif
			if ($gapswitch1 == 1)
				{
				$rev_seq1.='-';
				}
			else
			{
			$rev_seq1.='-'; #deletion
			$i=$i;
			}
			if ($gapswitch2 == 2)
				{
				$rev_seq2.='-';
				}
			else
			{ 
			$rev_seq2.=$d_array[$j-1];
			# moving one step further in $j axis
			$j=$j-1;
		#	print "r = del i=$i, j=$j resultmatrix:$result_matrix->[$i][$j+1]\n";
			}
			if ($gapswitch3 == 3)
				{
				$rev_seq3.='-';
				}
			elsif ($gapswitch3 != 3)
			{
			$rev_seq3.=$u_array[$k-1];
			$k=$k-1;
			$track.='ir-';
			#print "in bc if\n";
			# leaving for loop
			}
			}# elsif end
		
		#insert down (= sequence 2)
		if (($cc gt $ac) and ($cc ge $bc) and ( $cc gt $dc) and ($cc gt $ec) and ($cc gt $fc) and ($cc gt $gc) and ($cc gt $hc))
	 		{# else 
			
			if ($gapswitch1 == 1)
				{
				$rev_seq1.='-';
				}
			else
			{
			$rev_seq1.=$r_array[$i-1];  
			# moving one step furhter in $i axis
			$i=$i-1;
			}
			
			if ($gapswitch2 == 2)
				{
				$rev_seq2.='-';
				}
			else
			{
			$rev_seq2.='-';
			$j=$j;
		#	print "d = ins i=$i, j=$j resultmatrix:$result_matrix[$i+1][$j]\n";
			}
			if ($gapswitch3 == 3)
				{
				$rev_seq3.='-';
				}
			elsif ($gapswitch3 != 3)
			{
			$rev_seq3.=$u_array[$k-1];
			$k=$k-1;
			$track.='id-';
			#print "in cc if\n";
			# leaving for loop
			}
			}
		
		#insert up (= sequence 3)
		if (($dc gt $ac) and ($dc gt $bc) and ( $dc gt $cc) and ($dc gt $ec) and ($dc gt $fc) and ($dc gt $gc) and ($dc gt $hc))
	 		{# else 
			if ($gapswitch1 == 1)
				{
				$rev_seq1.='-';
				}
			else
			{
			$rev_seq1.=$r_array[$i-1];  
			# moving one step furhter in $i axis
			$i=$i-1;
			}
			
			if ($gapswitch2 == 2)
				{
				$rev_seq2.='-';
				}
			else
			{
			$rev_seq2.=$d_array[$j-1];
			$j=$j-1;
		#	print "d = ins i=$i, j=$j resultmatrix:$result_matrix[$i+1][$j]\n";
			}
			
			if ($gapswitch3 == 3)
				{
				$rev_seq3.='-';
				}
			elsif ($gapswitch3 != 3)
			{
			$rev_seq3.='-';
			$k=$k;
			$track.='iu-';
			#print "in dc if\n";
 			}# else end
			}

		#insert right down (= sequence 1 and sequence 2 )
		if (($ec gt $ac) and ($ec gt $bc) and ( $ec gt $cc) and ($ec gt $dc) and ($ec gt $fc) and ($ec gt $gc) and ($ec gt $hc))
	 		{# else 
			if ($gapswitch1 == 1)
				{
				$rev_seq1.='-';
				}
			else
			{
			$rev_seq1.='-';  
			# moving one step furhter in $i axis
			$i=$i;
			}
			if ($gapswitch2 == 2)
				{
				$rev_seq2.='-';
				}
			else
			{
			$rev_seq2.='-';
			$j=$j;
#			print "d = ins i=$i, j=$j resultmatrix:$result_matrix[$i+1][$j]\n";
			}
			
			if ($gapswitch3 == 3)
				{
				$rev_seq3.='-';
				}
			elsif ($gapswitch3 != 3)
			{
			$rev_seq3.=$u_array[$k-1];
			$k=$k-1;
			$track.='ird-';
			#print "in ec if\n";
 			}
			}# else end


		#insert right up (= sequence 1 and sequence 3 )
		if (($fc gt $ac) and ($fc gt $bc) and ( $fc gt $cc) and ($fc gt $dc) and ($fc gt $ec) and ($fc gt $gc) and ($fc gt $hc))
	 		{# else 
			if ($gapswitch1 == 1)
				{
				$rev_seq1.='-';
				}
			else
			{
			$rev_seq1.='-';  
			# moving one step furhter in $i axis
			$i=$i;
			}
			
			if ($gapswitch2 == 2)
				{
				$rev_seq2.='-';
				}
			else
			{
			$rev_seq2.=$d_array[$j-1];
			$j=$j-1;
		#	print "d = ins i=$i, j=$j resultmatrix:$result_matrix[$i+1][$j]\n";
			}
			
			if ($gapswitch3 == 3)
				{
				$rev_seq3.='-';
				}
			elsif ($gapswitch3 != 3)
			{
			$rev_seq3.='-';
			$k=$k;
			$track.='iru-';
			#print "in fc if\n";
			}
 			}# else end


		#insert down up (= sequence 2 and sequence 3 )
		if (($gc gt $ac) and ($gc gt $bc) and ( $gc gt $cc) and ($gc ge $dc) and ($gc ge $ec) and ($gc ge $fc) and ($gc ge $hc))
	 		{# else
			 
			if ($gapswitch1 == 1)
				{
				$rev_seq1.='-';
				}
			else
			{
			$rev_seq1.=$r_array[$i-1];  
			# moving one step furhter in $i axis
			$i=$i-1;
			}
			
			if ($gapswitch2 == 2)
				{
				$rev_seq2.='-';
				}
			else
			{
			$rev_seq2.='-';
			$j=$j;
		#	print "d = ins i=$i, j=$j resultmatrix:$result_matrix[$i+1][$j]\n";
			}
			
			if ($gapswitch3 == 3)
				{
				$rev_seq3.='-';
				}
			elsif ($gapswitch3 != 3)
			{
			$rev_seq3.='-';
			$k=$k;
			$track.='idu-';
			#print "in gc if\n";
			}
 			}# else end





			
#		} # for end
#	}#for end

	} # while end	
#print "end of comparsion of strings\n";

=print out for debugging
# debugging


print "revseq1:$rev_seq1\n";
for( $i=0 ; $i < ($d); $i++)
	{
	print "d_array element: $d_array[$i]\n";
	}
print "revseq2:$rev_seq2\n";

print "r_array: ";
for( $i=0 ; $i < ($r); $i++)
	{
	print "$r_array[$i] ";
	}
	print "\n";
=end comment
=cut


# inversing the "backtracked" sequences to get the aligned sequence in originally orientation
$seq1 = reverse $rev_seq1;
$seq2 = reverse $rev_seq2;
$seq3 = reverse $rev_seq3;

# output of the alignment to the screen
print "\n global Alignment\n";
print "seq(length:$r):  $seq1\n"; 
print "seq(length:$d):  $seq2\n"; 
print "seq(length:$u):  $seq3\n";


#debugging 
 my $rev_track = reverse $track;
# print "tracking = $rev_track \n";
print "maximum of scoring matrix: $result_matrix->[$r][$d][$u] \n";
