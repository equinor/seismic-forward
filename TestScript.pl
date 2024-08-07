#!/usr/bin/env perl

# o==========================================================================o
# | AUTHOR : Pål Dahle, Norsk regnesentral                                   |
# | VERSION: 2017/09/28                                                      |
# | PURPOSE: Compare Seismic-Forward output result with correct answers. To include |
# |          new test cases the adopted file and directory structure must be |
# |          be followed.                                                    |
# |                                                                          |
# |          We compare all traces and all samples in the volume. For each   |
# |          sample a[i] in the answer volume and o[i] in the output volume  |
# |          we calculate their absolute differences, ie.                    |
# |                                                                          |
# |          abs_diff[i] = abs(abs(a[i]) - abs(o[i]))                        |
# |                                                                          |
# |          If the the maximum absolute relative difference is less than a  |
# |          given threshold, the test has passed                            |
# |                                                                          |
# o==========================================================================o

use strict;
use File::Path;
use File::Copy;
use Cwd;

#----------------------------------------------------------------
#                       SUB section
#----------------------------------------------------------------

#---------------------------
sub FindBaseAndExeDir( $$$ )
#---------------------------
{
    my $exedir  = shift;
    my $basedir = shift;
    my $debug   = shift;
    my $curdir  = Cwd::getcwd();

    if (!$exedir) {
        $exedir = $curdir;
    }
    else {
        $exedir = JoinPaths($curdir, $exedir);
    }
    if (File::Basename::basename($exedir) eq File::Basename::basename($basedir)) {
        print "Please specify the path to the executable as \'TestScript.pl dir=path-to-executable\'\n";
        exit;
    }
    $basedir = JoinPaths($curdir, $basedir);
    return ($exedir, $basedir);
}

#------------------
sub JoinPaths( $$ )
#------------------
{
    my $path1 = shift;
    my $path2 = shift;
    my $path  = File::Spec->catfile($path1, $path2);
    $path =~ s/\/[\.\_\-\w]+\/\.\.//;
    return $path;
}

#---------------------------
sub Initialize( $$$$$$$$$$ )
#---------------------------
{
    my $exedir        = shift;
    my $basedir       = shift;
    my $testdir       = shift;
    my $make_segy     = shift;
    my $make_storm    = shift;
    my $geo2seis      = shift;
    my $compare_segy  = shift;
    my $compare_storm = shift;
    my $dir           = shift;
    my $debug         = shift;

    if (! -e $geo2seis) {
        print "\nThe executable \'$geo2seis\' does not exist\n";
        print "\nQuitting...\n\n";
        exit;
    }

    if (! -e $compare_segy) {
        if (! -e $basedir) {
            print "\nCannot find source code directory \'$basedir\'\n";
            print "\nQuitting...\n\n";
            exit;
        }
        if (! -e $make_segy) {
            print "\nCannot find shell script for making SegY cube comparison function \'$make_segy\'\n";
            print "\nQuitting...\n\n";
            exit;
        }

        Cwd::chdir($exedir);
        print "\nMaking Segy comparison executable \'$compare_segy\'\n";
        open (MAKE, "$make_segy $dir |");
        while ( <MAKE> ) { print $_; }
        close(MAKE);
    }

    if (! -e $compare_segy) {
        print "\nThe shell script \'$make_segy\' did not produce the executable \'$compare_segy\'";
        print "\nQuitting...\n\n";
        exit;
    }

    if (! -e $compare_storm) {
        if (! -e $basedir) {
            print "\nCannot find source code directory \'$basedir\'\n";
            print "\nQuitting...\n\n";
            exit;
        }
        if (! -e $make_storm) {
            print "\nCannot find shell script for making SegY cube comparison function \'$make_storm\'\n";
            print "\nQuitting...\n\n";
            exit;
        }

        Cwd::chdir($exedir);
        print "\nMaking STORM comparison executable \'$compare_storm\'\n";
        open (MAKE, "$make_storm $dir |");
        while ( <MAKE> ) { print $_; }
        close(MAKE);
    }

    if (! -e $compare_storm) {
        print "\nThe shell script \'$make_storm\' did not produce the executable \'$compare_storm\'";
        print "\nQuitting...\n\n";
        exit;
    }

    if ($debug) {
        print "Changing directory to $testdir\n";
    }
    Cwd::chdir($testdir);
}

#----------------------
sub ExtractOptions( @ )
#----------------------
{
    my @ARG = @_;

    my $exedir;
    my $passive  = 0;
    my @cases    = ();

    foreach my $arg (@ARG) {
        if ($arg =~ /passive=yes/) {
            $passive = 1;
        }
        elsif ($arg =~ /dir=([\_\-\.\/\w]+)/) {
            $exedir = $1;
        }
        elsif ($arg =~ /case=/) {
            if ($arg =~ /(\d+)\-(\d+)/) {                     # Test the following cases (specified as range: 1-4)
                my $first = $1;
                my $last  = $2;
                for (my $i = $first ; $i < $last +1 ; $i++) {
                    push @cases, $i;
                }
            }
            elsif ($arg =~ /([\d+\,]+)/) {
                @cases = split(',', $1);                      # Test the following cases (specified as CSV: 1,2,3,...)
            }
            else {
                die "Usage: $0 dir=executable-dir [passive=yes] [comma-separated case numbers] [range-specified case numbers]\n";
            }
        }
        else {
            die "Usage: $0 dir=executable-dir [passive=yes] [comma-separated case numbers] [range-specified case numbers]\n";
        }
    }
    return ($exedir, $passive, @cases);
}

#-----------------------
sub RunGeo2Seis( $$$$$ )
#-----------------------
{
    my $geo2seis   = shift;
    my $testdir    = shift;
    my $modeldir   = shift;
    my $modelfile  = shift;
    my $debug      = shift;
    my $outputdir  = $modeldir."/output";

    #
    # Remove old results...
    #
    if (-e $outputdir) {
        if ($debug) {
            print "   Removing old result...\n";
        }
        File::Path::rmtree($outputdir,$debug,1);
    }
    File::Path::mkpath($outputdir,$debug);

    #
    # Copy Geo2Seis executable and model file to output directory...
    #
    if (-e $geo2seis) {
        File::Copy::copy($geo2seis, $outputdir);
    }
    else {
        print "   Cannot find the executable \'$geo2seis\'\n";
        exit;
    }

    #
    # Copy model file
    #
    if (-e "$modeldir/$modelfile") {
        File::Copy::copy("$modeldir/$modelfile", $outputdir);
    }
    else {
        print "   Cannot find model file \'$modelfile\' in directory \'$modeldir\'\n";
        print "   Current directory is \'".Cwd::getcwd()."\'\n";
        exit;
    }

    #
    # Change directory and run Geo2Seis...
    #
    Cwd::chdir($outputdir);
    my $executable = File::Basename::basename($geo2seis);
    chmod 0755, $executable;

    if (!$debug) {
        print "   Running Seismic Forward ...\n";
    }
    open (GEO2SEIS, "./$executable $modelfile |");
    while ( <GEO2SEIS> ) {
        if ($debug) {
            print $_;
        }
    }
    close(GEO2SEIS);
    Cwd::chdir($testdir);
}

#---------------------------
sub CheckSegYGrids( $$$$$$ )
#---------------------------
{
    my $modeldir   = shift;
    my $executable = shift;
    my $r_params   = shift;
    my $threshold  = shift;
    my $debug      = shift;
    my $showall    = shift;

    my $outputdir  = "$modeldir/output";
    my $answerdir  = "$modeldir/answer";

    my @parameters = @$r_params;

    #
    # What kind of SegY grids to look for
    #
    my @filenames = ();
    for my $i (0 .. $#parameters) {
        push @filenames, $parameters[$i].".segy";
    }

    my @messages    = ();
    my $active_test = 0;

    my $match = 1;

    for my $i (0 .. $#filenames) {
        my $name    = $filenames[$i];
        my $file1   = "$outputdir/$name";
        my $file2   = "$answerdir/$name";
        my $exists1 = -e $file1;
        my $exists2 = -e $file2;

        if ($exists2) { # Require answer volume before looking for output
            my ($grid_defs_are_equal,
                $max_amp,
                $max_diff,
                $avg_diff,
                $bias,
                $max_trace,
                $max_sample) = CompareSegYGrids($executable, $file1,$file2,$debug);

            if ($exists1) {
                if ($grid_defs_are_equal == -1) {
                    push @messages, sprintf("   Grid : %43s  =>  ERROR: The Segy-comparison program failed (run from command line to see error message).\n",$name);
                    $match = 0;
                }
                else {
                    if ($grid_defs_are_equal == 1) {
                        if ($max_amp == 0) {
                            push @messages, sprintf("   Grid : %43s  =>  ERROR: Maximum amplitude of answer cube is zero\n",$name);
                            $match = 0;
                        }
                        else {
                            my $relative_diff = $max_diff/$max_amp;
                            if ($relative_diff < $threshold) {
                                push @messages, sprintf("   Grid : %43s  =>  MATCH\n",$name);
                            }
                            else {
                                $bias = 0.0 if ($bias < -0.000000);
                                push @messages, sprintf("   Grid : %43s  =>  ERROR: MaxAmplitude=%.6f  MaxDiff=%.6f  AvgDiff=%.6f  Bias=%.6f  Worst(Trace,Sample)=(%d,%d)\n",
                                                        $name,$max_amp,$max_diff,$avg_diff,$bias,$max_trace,$max_sample);
                                $match = 0;
                            }
                        }
                    }
                    else {
                        push @messages, sprintf("   Grid : %43s  =>  ERROR: Grid definitions are not equal\n",$name);
                        $match = 0;
                    }
                }
            }
            else {
                push @messages, sprintf("   Grid : %43s  =>  ERROR: Got correct answer but no Seismic Forward output to compare with.\n",
                                        $name);
                $match = 0;
            }
            $active_test = 1;
        }
    }
    if ($match && $active_test && !$showall) {
        printf "   Grids :                                        All  =>  MATCH\n";
    }
    else {
        foreach ( @messages ) {
            print;
        }
    }
    return $match;
}

#---------------------------
sub CompareSegYGrids( $$$$ )
#---------------------------
{
    my $executable = shift;
    my $file1      = shift;
    my $file2      = shift;
    my $debug      = shift;

    if (!-e $executable) {
        print "Could not find SegY comparison executable \'$executable\'\n";
        print "Current directory is \'".Cwd::getcwd()."\'\n";
        exit;
    }
    if ($debug) {
        print "Comparing Segy cubes \'$file1\' and \'$file2\'\n";
    }
    open (COMPARE, "$executable $file1 $file2 |");
    while ( <COMPARE> ) {
        if ($debug) {
            print $_;
        }
    }
    close(COMPARE);

    my $diff_file = "segy_amplitude_difference.txt";

    open(IN, $diff_file) or die "$!\nCould not find file \'$diff_file\'\n";
    my $values = <IN>;
    close(IN);

    my @values = split(' ', $values);

    my $grid_defs_are_equal = $values[0];
    my $max_amp             = $values[1];
    my $max_diff            = $values[2];
    my $avg_diff            = $values[3];
    my $bias                = $values[4];
    my $max_trace           = $values[5];
    my $max_sample          = $values[6];

    return ($grid_defs_are_equal,
            $max_amp,
            $max_diff,
            $avg_diff,
            $bias,
            $max_trace,
            $max_sample);
}

#----------------------------
sub CheckStormGrids( $$$$$$ )
#----------------------------
{
    my $modeldir   = shift;
    my $executable = shift;
    my $r_params   = shift;
    my $threshold  = shift;
    my $debug      = shift;
    my $showall    = shift;

    my $outputdir  = "$modeldir/output";
    my $answerdir  = "$modeldir/answer";

    my @parameters = @$r_params;

    #
    # What kind of SegY grids to look for
    #
    my @filenames = ();
    for my $i (0 .. $#parameters) {
        push @filenames, $parameters[$i].".storm";
    }

    my @messages    = ();
    my $active_test = 0;

    my $match = 1;

    for my $i (0 .. $#filenames) {
        my $name    = $filenames[$i];
        my $file1   = "$outputdir/$name";
        my $file2   = "$answerdir/$name";
        my $exists1 = -e $file1;
        my $exists2 = -e $file2;

        if ($exists2) { # Require answer volume before looking for output
            my ($grid_defs_are_equal,
                $largest_difference,
                $mean_abs_difference,
                $mean_val_difference,
                $mean_val_both_cubes) = CompareStormGrids($executable, $file1,$file2,$debug);

            if ($exists1) {
                if ($grid_defs_are_equal == 1) {
                    my $scale_factor        = 1.0;
                    my $relative_difference = $largest_difference/$mean_val_both_cubes;
                    if ($relative_difference < 0.0) {
                        $relative_difference *= -1.0;
                    }
                    if ($relative_difference < $threshold*$scale_factor || $largest_difference < $threshold*$scale_factor) {
                        push @messages, sprintf("   Grid : %43s  =>  MATCH\n",$name);
                    }
                    else {
                        push @messages, sprintf("   Grid : %43s  =>  MaxDiff=%.6f  MeanAbsDiff=%.6f  Bias=%.6f\n",
                                                $name,$largest_difference,$mean_abs_difference,$mean_val_difference);
                        $match = 0;
                    }
                }
                else {
                    push @messages, sprintf("   Grid : %43s  =>  ERROR: Grid definitions are not equal\n",$name);
                    $match = 0;
                }
            }
            else {
                push @messages, sprintf("   Grid : %43s  =>  ERROR: Got correct answer but no Seismic Forward output to compare with.\n",
                                        $name);
                $match = 0;
            }
            $active_test = 1;
        }
    }
    if ($match && $active_test && !$showall) {
        printf "   Grids :                                        All  =>  MATCH\n";
    }
    else {
        foreach ( @messages ) {
            print;
        }
    }
    return $match;
}

#----------------------------
sub CompareStormGrids( $$$$ )
#----------------------------
{
    my $executable = shift;
    my $file1      = shift;
    my $file2      = shift;
    my $debug      = shift;

    if (!-e $executable) {
        print "Could not find STORM comparison executable \'$executable\'\n";
        print "Current directory is \'".Cwd::getcwd()."\'\n";
        exit;
    }
    if ($debug) {
        print "Comparing STORM volumes \'$file1\' and \'$file2\'\n";
    }
    open (COMPARE, "$executable $file1 $file2 |");
    while ( <COMPARE> ) {
        if ($debug) {
            print $_;
        }
    }
    close(COMPARE);

    my $diff_file = "storm_volume_difference.txt";

    open(IN, $diff_file) or die "$!\nCould not find file \'$diff_file\'\n";
    my $values = <IN>;
    close(IN);

    my @values = split(' ', $values);

    my $grid_defs_are_equal = $values[0];
    my $largest_difference  = $values[1];
    my $mean_abs_difference = $values[2];
    my $mean_val_difference = $values[3];
    my $mean_val_both_cubes = $values[4];

    return ($grid_defs_are_equal,
            $largest_difference,
            $mean_abs_difference,
            $mean_val_difference,
            $mean_val_both_cubes);
}


#----------------------------------------------------------------
#                       EDIT section
#----------------------------------------------------------------

my @modeldir    =  (
                    "01_nmo_pp",                           # test directories
                    "02_nmo_pp_noise",
                    "03_nmo_ps",
                    "04_off_pp",
                    "05_off_pp_wavelet_from_file",
                    "06_pp",
                    "07_ps",
                    "08_ps_noise",
                    "09_remove_negative_dz",
                    "10_keep_negative_dz",
                    "11_corner_point",
                    "12_Reek_center_point_interpolation",
                    "13_Reek_corner_point_interpolation",
                    "14_Drogon_center_point_interpolation_nmo",
                    "15_Drogon_corner_point_interpolation_nmo",
                    "16_Drogon_center_point_interpolation_nmo_white_noise"
                    );

my @parameters  =  ("seismic_depth",
                    "seismic_depth_stack",
                    "seismic_time",
                    "seismic_time_stack",
                    "seismic_timeshift",
                    "seismic_time_prenmo",
                    "vp",
                    "vs",
                    "rho",
                    "twt",
                    "zgrid",
                    "vrms",
                    "reflections_0");

my $threshold   =  1.0e-08;                                # Match if diff < threshold
#my $threshold   =  3.0e-05;                                # Match if diff < threshold

my $showall     =  1;                                      # 1=yes, 0=no (write result for each file tested)
my $debug       =  0;                                      # 1=yes, 0=no

#----------------------------------------------------------------
#                       MAIN section
#----------------------------------------------------------------

my $time0   = time();

my ($tmpdir,
    $passive,
    @cases)       = ExtractOptions(@ARGV);

$0 =~ m|^(.*/)|;
my $dir = $1;

my ($exedir,                                               # This is where the seismic forward and comparison executables are stored
    $basedir)     = FindBaseAndExeDir($tmpdir, $dir, $debug);

my $testdir       = $basedir."/test_suite";                  # Where tests are stored
my $make_segy     = $basedir."/make_compare_traces.sh";      # Script for making the comparison function
my $make_storm    = $basedir."/make_compare_storm_grids.sh"; # Script for making the comparison function
my $geo2seis      = $exedir."/seismic_forward";              # The geo2seis executable
my $compare_segy  = $exedir."/compare_traces";               # The Segy comparison executable
my $compare_storm = $exedir."/compare_storm_grids";          # The STORM comparison executable
my $modelfile     = "modelfile.xml";
my $os            = $^O;

print "\n*******************************************************************";
print "\n*****            Seismic Forward test suite                   *****";
print "\n*******************************************************************\n";

print "\nPlatform                    : $os";
print "\nGeo2Seis executable         : $geo2seis";
print "\nSegY comparison executable  : $compare_segy";
print "\nSTORM comparison executable : $compare_storm";
print "\nTest directory              : $testdir\n";
print "\nPassive mode                : $passive\n";
print "\nThreshold                   : $threshold\n";

Initialize($exedir,$basedir,$testdir,$make_segy,$make_storm,$geo2seis,$compare_segy,$compare_storm,$dir,$debug);

my $ok1;
my $ok2;
my $test_ok;
my $all_ok = 1;
my @ok = ();

for my $i (0 .. $#modeldir) {
    $ok[$i] = 1;
}

my @run_case = @ok;

if (@cases) {     # When only a subset of the cases are wanted.
    print "\nRunning test cases          : ";
    for (my $i = 0 ; $i < @modeldir ; $i++) {
        $run_case[$i] = 0;
    }
    foreach my $case (@cases) {
        if ($case > 0 && $case <= @modeldir) {
            print $modeldir[$case - 1]." ";
            $run_case[$case - 1] = 1;
        }
    }
    print "\n";
}

for my $i (0 .. $#modeldir) {
    if ($run_case[$i]) {
        my $modeldir = $modeldir[$i];

        print "\nTEST : $modeldir\n";

        # Run test
        # --------

        if (!$passive) {
            RunGeo2Seis($geo2seis,$testdir,$modeldir,$modelfile,$debug);
        }

        # Check seismic cubes
        # -------------------

        $ok1 = CheckSegYGrids($modeldir,$compare_segy,\@parameters,$threshold,$debug,$showall);

        # Check storm volumes
        # -------------------

        $ok2 = CheckStormGrids($modeldir,$compare_storm,\@parameters,$threshold,$debug,$showall);

        $test_ok = $ok1 && $ok2;
        $all_ok  = $all_ok && $test_ok;
        $ok[$i]  = $test_ok;
    }
}

my $time1 = time();

if ($all_ok) {
    print "\n*** FINISHED SUCCESSFULLY - NO WARNINGS/ERRORS ***\n\n";
    print "*** Time : ",$time1-$time0,"s\n";
    exit 0;
}
else {
    print "\n*** FINISHED WITH WARNINGS/ERRORS FOR MODELS:***\n\n";
    for my $i (0 .. $#modeldir) {
        print "   $modeldir[$i]\n" if ( !$ok[$i] );
    }
    print "\n";
    print "*** Time : ",$time1-$time0,"s\n";
    exit 1;
}
