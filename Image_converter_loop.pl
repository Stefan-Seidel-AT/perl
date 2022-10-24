#!/usr/bin/perl


# This program expects a in image file (.jpg) as an argument and 
# converts bits of the image. It iterates through various settings and saves output (images) 


# usage:     ./Image_converter_loop "your_image_file.jpg"
# example:   ./Image_converter_loop IMAGE2.jpg

use warnings;
use strict;
use Imager;

# takes file as first commandline argument
my $file = shift;
my $img = Imager->new();

$img->open(file=>$file) or die $img->errstr();

#my $newimg = $img->convert(preset=>'grey');
#$newimg->write(file=>'grey.jpg');

#my $newimg = $img->convert(preset=>'red');
#$newimg->write(file=>'red.jpg');

# Swap red/green channel  
#my $new = $img->convert(matrix=>[ [ 0, 1, 0 ],
#                                  [ 1, 0, 0 ],
#                                [ 0, 0, 1 ] ]);

# Swap channel  
my $new = $img->convert(matrix=>[ [ 1, 0, 0 ],
                                  [ 1, 0, 0 ],
                                [ 1, 0, 0 ] ]);




#$new -> write(file=>'convert.jpg');

=begin comment

# if Gamma should be altered as well enable this block
my @i = (1..9);

for (@i) {
Apply a Gamma of 1.4
my $gamma = '1.'.$_;
my @map = map { int( 0.5 + 255*($_/255)**$gamma ) } 0..255;
$new->map(all=>\@map);  # inplace conversion
my $outputfile = 'raw_gamma.'.$gamma.'.jpg';
$new->write(file=>$outputfile);
print ("$gamma\n");
}
=cut


# create arrays for loops
my @a = ();
my  $y=1;
for ($y=1; $y<25; $y++ ) 
	 { push(@a,($y*10));
		#print ("$y","\n");
		}
=comment		
my @b = ();
$y=1;
for ($y=200; $y<256; $y=$y+10 ) 
	 { push(@b,($y));
#		print ("$y","\n");
		}
=cut
		
# loop throug arrays to create pictures with different thresholds and substitution colors
my $counter = 0;
my $bitcounter = 0;
for my $a(@a) 
	{ 
	#for my $b(@b)
		#{
	# obtain input file for every loop?
	#$img->open(file=>$file) or die $img->errstr();
	my $new = $img->convert(matrix=>[ [ 1, 0, 0 ],
                                      [ 0, 1 , 0 ],
                                      [ 0, 0, 1 ] ]);

	
	my $gamma = '1.2';
	my $threshold = $a;
	my $mask_color = '255';
	
	
	my @map = map { if (int($_) < $threshold) {($_ = $mask_color); ($bitcounter++)}} 0..255;
	#my $input2 = $_;
	#	@map = map { int( 0.5 + 255*($input2/255)**$gamma ) } 0..255;

	#for(@map) {print ("$_","\n");}
	$new->map(all=>\@map);  # inplace conversion
	my $outputfile = 'try4_mask.'.$gamma.'_trs_'.$threshold.'_mskcol_'.$mask_color.'_'.$counter.'.jpg';
	$new->write(file=>$outputfile);
	print ("write to file: ".$outputfile."  bitcounter = ".$bitcounter. "\n");
	$bitcounter = 0;
	$counter = $counter+1;
	@map = {};
		#}
	}
