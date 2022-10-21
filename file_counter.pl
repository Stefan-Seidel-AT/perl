#!/usr/bin/perl -w

# count all ".txt", ".doc", ".docx" , ".xls", ".xlsx", ".pl" , ".jpg" , ".ppt" , ".pptx" files in
# the current directory
my $dirname;
my $count;
my $count2;
my @fileends;
my $end;
my $input;
my $other_files;

$dirname = 'C:\Users\stefans\Desktop';

@fileends = ( "txt", "doc", "docx" , "xls", "xlsx", "pl" , "jpg" , "ppt" , "pptx" , ".*");

print " You have these numbers of different files in following directory: \n $dirname \n\n";


foreach $end (@fileends) { 

opendir(DIR, $dirname);

@txt_files = grep(/\.$end$/,readdir(DIR));
$count = @txt_files;

if ( grep( /\.*/, @txt_files) ) {
  #TODO:
 # print "@txt_files \n";
  $count2 = @txt_files
}


closedir(DIR);

# print all the filenames in our array
foreach my $file (@txt_files) {
   #print "$file\n";
   
}

# normalize the length for human readability
while (length($end) < 4) {
	$end.=" ";
	}
print "   .$end  = $count files found \n";
}
print " In total $count2 files and/or directorys were found. \n";

# to manually end the programm
print " \n press Enter to contine... ";

$input = <STDIN>;
