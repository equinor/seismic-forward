#!/usr/bin/perl -w

use strict;

#
# Reading directory
#
opendir MYDIR, ".";
my @content = readdir MYDIR;
closedir MYDIR;

print "Updating all cohiba.log files in answers from last run\n";

#
# Only keep test case directories
#
my @dirs = ();
foreach ( @content ) {
    if (/^\d+\_/) {
        push @dirs, $_;
    }
    @dirs = sort @dirs;
}

#
# Loop over test cases
#
my $count = 0;
foreach my $dir ( @dirs ) {
    print "cp $dir/output/Logfile.txt $dir/answer/\n";
    system("cp $dir/output/Logfile.txt $dir/answer/ >/dev/null 2>&1");
    $count++;
}

print "Successfully updated $count files.\n";
