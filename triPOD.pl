#!/usr/bin/perl

############################################################################### 
#
# The triPOD software detects chromosomal abnormalities using 
# the Parent-of-Origin-based Detection (POD) method.
# Author and maintainer: Joseph D. Baugher. 
# Email:jbaughe2(at)jhmi.edu.
#
# This file is part of the triPOD software.
# triPOD is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	(at your option) any later version.
# triPOD is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
#   along with triPOD.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

use strict;
use warnings;
use Algorithm::Cluster 'kcluster';
use Cwd 'getcwd';
use File::Path qw(make_path remove_tree);
use File::Temp qw(tempfile tempdir);
use Getopt::Long;
use POSIX qw(ceil floor);
use threads;
use Time::HiRes 'usleep';
use Tree::Interval;
use Version_info;

my $version_info = Version_info->new();
print "\n", $version_info->get_name, " ", $version_info->get_version, "\n\n";
my $CMD_LINE = join(" ", $0, @ARGV);

# Detect the total number of cores
my $MAX_CORES = 1;
my $os = $^O;
if ($os =~ /mswin/i){warn("triPOD has not been tested on Windows!")}
elsif($os =~ /darwin/i){$MAX_CORES= `sysctl -n hw.ncpu`}
elsif($os =~ /linux/i){$MAX_CORES= `cat /proc/cpuinfo | grep processor | wc -l`}

my $THE_TIME = cpu_time();

# Get process ID, current working directory, and current time
my $CWD = getcwd();
my $PID_VALUE = $$;

# Default Parameters
my $ALPHA       = 0.1;
my $BATCH;
my $BUILD_FILE  = "./genome_build/hg18_centromeres.txt";
my $CITE;
my $CORES       = $MAX_CORES - 1;
my $GENDER      = 0;
my $GRAPHICS    = "none";
my $HD          = "T";
my $HELP;
# Number of stdevs from the heterozygous BAF mean - used as boundaries 
# for detection of informative SNPs. (The intent is to ignore at least 50%
# of the members of the "normal" heterozygous SNP distribution.)
my $HET_SD      = 1.414213562373095; # sqrt(2)
# Number of stdevs from the homozygous BAF means - used as boundaries 
# for mBAF conversions. (The intent is to ignore nearly all 
# of the members of the "normal" homozygous SNP distribution.)
my $HOM_SD      = 4;
my $MI1         = "T";
my $NC_THRESH   = 0.03;
my $OUTPUT_DIR;
my $POD         = "T";
my $PODcr       = "T";
my $SIZE_SMALL_WINDOW = 100;
my $STATS;
my $verbose     = "T";

# Get optional input parameters (e.g. --cores=2 --alpha=0.05)
GetOptions( 'alpha=f'  => \$ALPHA,            'homSD=f'  => \$HOM_SD,
            'batch=s'  => \$BATCH,            'mi1!'     => \$MI1,
            'build=s'  => \$BUILD_FILE,       'nc=f'     => \$NC_THRESH,
            'cite!'    => \$CITE,             'out=s'    => \$OUTPUT_DIR,
            'cores=i'  => \$CORES,            'pod!'     => \$POD,
            'gender=s' => \$GENDER,           'podcr!'   => \$PODcr,
            'graph=s'  => \$GRAPHICS,         'stats!'   => \$STATS,
            'hd!'      => \$HD,               'verbose!' => \$verbose,
            'help!'    => \$HELP,             'win=i'    => \$SIZE_SMALL_WINDOW,
            'hetSD=f'  => \$HET_SD                
          );

help() if $HELP;
if ($CITE) {print $version_info->get_reference, "\n\n"; exit 0} 
my $SIZE_LARGE_WINDOW = $SIZE_SMALL_WINDOW * 5;
my $CORES_PER_SAMPLE = $CORES;
my $AMP = 0.1;
my ($HD_LOW, $HD_HI) = (-2, -1.5);
my ($DEL_UPPER, $DEL_LOWER)   = (-0.1, -2);

# Array indices and various constants
use constant {
    CHR         => 1,  POS     => 2, START       => 1, STOP    => 2,
    LRR         => 1,  BAF     => 2, MALE        => 1, FEMALE  => 2,
    # Detection method
    POD         => 1,  HD      => 2, MI1         => 3, PODcr   => 4,
    # ID and contribution
    FATHER      => 1,  MOTHER  => 2, CHILD       => 3, NONE    => 4,
    BOTH        => 5,  UNKNOWN => 6, UNKNOWN_MI1 => 7, OUTLIER => 8,

    CHR_X       => 23, NA      => -1000,
    PI          => 4 * atan2(1, 1)
};

# Determine if gender was assigned
if    ($GENDER =~ /^m/i) {$GENDER = MALE}
elsif ($GENDER =~ /^f/i) {$GENDER = FEMALE}
elsif ($GENDER =~ /^n/i) {$GENDER = 0}
elsif ($GENDER) {
    $GENDER = 0; 
    print "The gender specification cannot be recognized. Please use M or F.",
        "\nProceeding without analysis of Chromosome X.\n\n";
}

if ($GRAPHICS ne "none" && $GRAPHICS ne "NA") {
    unless ($GRAPHICS eq "png" || $GRAPHICS eq "pdf" || $GRAPHICS eq "both") { 
        if ($GRAPHICS eq "PNG") {$GRAPHICS = "png"}
        elsif ($GRAPHICS eq "PDF") {$GRAPHICS = "pdf"}
        elsif ($GRAPHICS eq "BOTH") {$GRAPHICS = "both"}
        else {
            print STDERR "WARNING: The graphics input was not recognized\n. ",
                "Acceptable graphics input options are --graph=none, ",
                "--graph=png, --graph=pdf, or --graph=both.\n",
                "The analysis will proceed without graphical output.\n";
            $GRAPHICS = "";
        }
    }
}
else {$GRAPHICS = ""}

# Check if batch or single mode is correctly set by user
my $check_batch = ($BATCH) ? $BATCH : $ARGV[-1];

if(!open(INPUT_FILE, "<", $check_batch)) { 
    print STDERR "ERROR: Could not open input file: $check_batch!\n";
    die;
}
my $line = <INPUT_FILE>;
chomp($line);
my @line =  split("\t", $line);
close INPUT_FILE;
if (scalar(@line) == 1 && -e $line[0] && !$BATCH) {
    $BATCH = $check_batch; 
    print STDERR "WARNING: The input appears to be a list of files. ",
        "Switching to BATCH mode.\n";
}
elsif (scalar(@line) == 12 && $BATCH) {
    push(@ARGV, $check_batch);
    $BATCH = 0;
    print STDERR "WARNING: The input does not appear to be a list of files. ",
        "Switching to single sample mode.\n";
}

# Open batch file if supplied
my (@files, @batch, @batches);
if ($BATCH) {
    $verbose = "";
    open(BATCH_FILE, "<", $BATCH) || die("Could not open batch file!\n");
    while (<BATCH_FILE>) {
        my $file = $_;
        chomp($file);
        push(@files, $file);
    }
    close BATCH_FILE;
    $CORES_PER_SAMPLE = 1;
}
else {
    if ($ARGV[-1]) {
        $batches[0] = [$ARGV[-1]];
        push(@files, $ARGV[-1]);
    }
    else {die("Please include an input file.\n")}
}

# Open genome build file
my @CENTRO_REFS;
open(BUILD_FILE, "<", $BUILD_FILE)
    || die("\nCould not open genome build file!\n",
    "Please provide a correct file path.\n",
    "(e.g. --build=./hg18_centromeres.txt)\n\n");
    
<BUILD_FILE> for 1..2;
while (<BUILD_FILE>) {
    my $line = $_;
    chomp($line);
    my @line_array = split(/\t/, $line);
    $line_array[1] = (split(/chr/, $line_array[1]))[1];
    
    # Check Chromosome column for anything other than acceptable
    # chromosomes
    if ($line_array[1] !~ /^[1-9]$|^[1][0-9]|^[2][0-2]$|x|y|xy|m|mt/i) {
        print "\nUnrecognized character in Chromosome column of genome ",
        "build file!\n";
        print_proper_format();
        die;
    }
    # Check Position columns for non-numeric characters and negative
    # numbers  
    elsif ($line_array[2] !~ /^-?\d/ || $line_array[2] < 0) {
        print "\nUnrecognized character in Position Start column of genome ",
        "build file!\n";
        print_proper_format();
        die();
    }
    elsif ($line_array[3] !~ /^-?\d/ || $line_array[3] < 0) {
        print "\nUnrecognized character in Position End column of genome ",
        "build file!\n";
        print_proper_format();
        die();
    }
    
    if ($line_array[CHR] !~ /^\d/) {
        ($line_array[CHR] eq "X" && $GENDER) ? $line_array[CHR] = CHR_X : next;
    }
    $CENTRO_REFS[$line_array[1]] = [@line_array[1..3]];
}

my $input_name;
if ($BATCH) {$input_name = (split(/[\/]/,(split(/.txt/, $BATCH))[0]))[-1]}
else {$input_name = (split(/[\/]/,(split(/.txt/, $ARGV[-1]))[0]))[-1]}
if (!$OUTPUT_DIR) {$OUTPUT_DIR = "$CWD/triPOD_Results/$input_name"}

# If user input contains a trailing backslash, remove it
else {$OUTPUT_DIR =~ s|/\z||}

make_path("$OUTPUT_DIR/log_files");
-e "$OUTPUT_DIR/log_files" || die "Could not create output directory: ",
    "$OUTPUT_DIR\n";
print "Output directory = $OUTPUT_DIR/\n";
print "Number of Samples = ", scalar(@files), "\n";

my $OUTPUT_FILE;
if ($BATCH){$OUTPUT_FILE = "$OUTPUT_DIR/$input_name\_triPOD_Batch_Results.txt"}
else {$OUTPUT_FILE = "$OUTPUT_DIR/$input_name\_triPOD_Results.txt"}
open(OUTPUT_FILE, '>', $OUTPUT_FILE) || die 
    "Could not create output file - $OUTPUT_FILE!\n";
print OUTPUT_FILE join("\n", $THE_TIME, $CMD_LINE), "\n";
print OUTPUT_FILE $version_info->get_name, "\t", $version_info->get_version,
	"\tAuthor: ", $version_info->get_author, "\t", 
	$version_info->get_author_email, "\n\n"; 
if ($BATCH) {
    print OUTPUT_FILE "Detected Chromosomal Abnormalities", "\n",
        join("\t", qw(Sample Chr Start Stop Type Parent_of_Origin 
            Inheritance Size(SNPs) Informative_SNPs Size(bp) 
            Detection Median_mBAF Median_LRR Father-Median_mBAF 
            Father-Median_LRR Mother-Median_mBAF Mother-Median_LRR)
        ),"\n";
}

my $BED_FILE = "$OUTPUT_DIR/$input_name\_triPOD_Results.bed";
open(BED_FILE,'>',$BED_FILE) || die "Could not create bed file - $BED_FILE!\n"; 
print BED_FILE "track name=\"triPOD Results\" type=bedDetail ",
    "description=\"triPOD : Chromosomal Abnormality Detection\" visibility=2 ",
    "itemRgb=\"On\"\n";

# Redirect STDERR to a log file
my $PERL_LOG = "$OUTPUT_DIR/log_files/$input_name\_triPOD_log.txt";
open(STDERR, '>', $PERL_LOG) || die 
    "Could not create log file - $PERL_LOG!\n"; 

# Create stats file if requested
my $STATS_FILE;
if ($STATS) {
    $STATS_FILE = "$OUTPUT_DIR/$input_name\_triPOD_stats.txt";
    open(STATS_FILE, '>', $STATS_FILE) || die 
        "Could not create stats file - $STATS_FILE!\n"; 
    print STATS_FILE "File\tChild\tFather\tMother\tChild_NC_rate\t",
    "Father_NC_rate\tMother_NC_rate\tChild_HD_rate\tChild_Min_HD\t",
    "Father_HD_rate\tFather_Min_HD\tMother_HD_rate\tMother_Min_HD\t",
    "MI1_rate\tMin_MI1\tMI1_Upper_Thresh\tMI1_Lower_Thresh\tMin_POD\t",
    "Corrected_alpha\talpha\tError_rate\tAcceptable_Errors\tAA_Bound\t",
    "BB_bound\tAB_upper_bound\tAB_lower_bound\tDetected_Regions\t",
    "POD_regions\tPODcr_regions\tMI1_regions\tHD_regions\n";
}

# Declare global variables
my ($AA_BOUND, $AB_LOWER_BOUND, $AB_UPPER_BOUND, $ACCEPTABLE_ERRORS, 
    $BB_BOUND, $BOUNDARY_EXTENSION, $CH_LRR_MED, $CH_mBAF_MED, $CH_NAME,
    $CURR_FILE, $ERROR_RATE, $FILENAME, $HET_MI1_CT, $HET_MI1_RATE, $hUPD_RATE,
    $INPUT_FILE,$LARGE_POD_ALPHA,$LARGE_WIN_ACCEPTABLE_ERRORS,$LOCAL_CH_BAF_MEAN, 
    $LOCAL_CH_mBAF_MED, $LOCAL_CH_LRR_MED, $LOCAL_P1_mBAF_MED,$LOCAL_P1_LRR_MED,
    $LOCAL_P2_mBAF_MED, $LOCAL_P2_LRR_MED, $MI1_LOWER_THRESH, $MI1_UPPER_THRESH, 
    $MIN_BAF_OUT, $MIN_CH_HD, $MIN_hUPD, $MIN_MI1, $MIN_LARGE_POD, $MIN_P1_HD, 
    $MIN_P2_HD, $MIN_POD, $P1_AA_BOUND, $P1_BB_BOUND, $P2_AA_BOUND,$P2_BB_BOUND,
    $P1_LRR_MED, $P1_mBAF_MED, $P1_NAME, $P2_LRR_MED, $P2_mBAF_MED, $P2_NAME,
    $PERL_TO_R_FILE, $PERL_TO_R_FILENAME, $POD_ALPHA, $PODCR_LOWER_THRESH, 
    $PODCR_UPPER_THRESH, $R_ATTEMPTS, $R_PID, $R_PID_FILE, $R_PID_FILENAME, 
    $R_TO_PERL_FILE, $R_TO_PERL_FILENAME, $REFINING);
my (@BAF_OUTLIERS, @CURR_CHR_REFS, @LRR_REFS, @P1_INF_SNPS, @P2_INF_SNPS);
my (%BAF_BY_POS, %BAF_OUTLIERS_BY_POS, %CH_HD_BY_POS, %LRR_BY_POS, %MI1_BY_POS,  
    %P1_HD_BY_POS, %P1_INF_SNP_BY_POS, %P2_HD_BY_POS, %P2_INF_SNP_BY_POS, 
    %SNP_BY_NUM, %SNP_BY_POS);
# Count of completed samples - [0] = with abn, [1] = no abn
my @completed = (0,0);  
my (@batch_threads, @output, @sorted_output, @stats, @trio_def);
my $next_file = 0;
my $main_results;
my @bed_array;
if ($BATCH) {
    for (my $i = 1; $i <= $CORES; $i++) {
        if ($next_file <= $#files) {
            $CURR_FILE = $files[$i - 1];
            my $batch_thread = threads->new(\&main_program);
            push(@batch_threads, $batch_thread);
            $next_file++;
        }
    }
}
else {
    $CURR_FILE = $files[0];
    $main_results = main_program();
    if ($main_results) {
        if ($main_results->[0]) {
            push(@output, @{$main_results->[0]});
            @sorted_output = sort {
                $$a[1] <=> $$b[1] || $$a[2] <=> $$b[2]} @output;
            if ($STATS) {
                push(@stats, @{$main_results->[1]});
                map {print STATS_FILE join("\t", @{$_}), "\n"} @stats;
            }
            if ($main_results->[3]) {$completed[0]++}
        }
        elsif($main_results->[3]) {$completed[1]++}
    }
    print OUTPUT_FILE join("\n", @trio_def), "\n\n";
    print OUTPUT_FILE "Detected Chromosomal Abnormalities", "\n",
        join("\t", qw(Sample Chr Start Stop Type Parent_of_Origin 
            Inheritance Size(SNPs) Informative_SNPs Size(bp) 
            Detection Median_mBAF Median_LRR Father-Median_mBAF 
            Father-Median_LRR Mother-Median_mBAF Mother-Median_LRR)
        ),"\n";
    map {print OUTPUT_FILE join("\t", @{$_}), "\n"} @sorted_output;
    to_bed(\@sorted_output);
}

if ($BATCH) {
    # Check if threads are finished
    while (threads->list(threads::running)) {
        sleep(1);
        for (0..$#batch_threads) { 
            if ($batch_threads[$_]->is_joinable()) {
                $main_results = $batch_threads[$_]->join();
                if ($main_results) {
                    if ($main_results->[0]) {
                        push(@output, @{$main_results->[0]});
                        @sorted_output = sort {
                            $$a[1] <=> $$b[1] || $$a[2] <=> $$b[2]} @output;
                        map {print OUTPUT_FILE join("\t", @{$_}), "\n"}
                            @sorted_output;
                        to_bed(\@sorted_output);
                        if ($STATS) {
                            push(@stats, @{$main_results->[1]});
                            map{print STATS_FILE join("\t", @{$_}),"\n"}@stats;
                        }
                        if ($main_results->[3]) {$completed[0]++}
                        (@output, @sorted_output, @stats) = (()) x 3;
                    }    
                    elsif($main_results->[3]) {$completed[1]++}
                }
                
                if ($next_file <= $#files) {
                    $CURR_FILE = $files[$next_file];
                    my $batch_thread = threads->new(\&main_program);
                    push(@batch_threads, $batch_thread);
                    $next_file++;
                }
            }
        }
    }
}

for (0..$#batch_threads) { 
    if ($batch_threads[$_]->is_joinable()) {
        $main_results = $batch_threads[$_]->join();
        if ($main_results) {
            if ($main_results->[0]) {
                push(@output, @{$main_results->[0]});
                @sorted_output = sort {
                    $$a[1] <=> $$b[1] || $$a[2] <=> $$b[2]} @output;
                map {print OUTPUT_FILE join("\t",@{$_}),"\n"}@sorted_output;
                to_bed(\@sorted_output);
                if ($STATS) {
                    push(@stats, @{$main_results->[1]});
                    map {print STATS_FILE join("\t", @{$_}), "\n"} @stats;
                }
                if ($main_results->[3]) {$completed[0]++}
                (@output, @sorted_output, @stats) = (()) x 3;
            }
            elsif($main_results->[3]) {$completed[1]++}
        }
        (@output, @sorted_output, @stats) = (()) x 3;
    }
}

if    ($GENDER == 1) {$GENDER = "male"}          
elsif ($GENDER == 2) {$GENDER = "female"}          
$GENDER ||= "NA"; 
my @tot_param;
push(@tot_param, "alpha=$ALPHA");
push(@tot_param, "batch=$BATCH") if ($BATCH);
push(@tot_param, "build=$BUILD_FILE", "cores=$CORES", "gender=$GENDER",
    ($GRAPHICS) ? "graph=$GRAPHICS" : "nograph", ($HD) ? "hd" : "nohd", 
    "hetSD=$HET_SD", "homSD=$HOM_SD", ($MI1) ? "mi1" : "nomi1", 
    "out=$OUTPUT_DIR/", ($POD) ? "pod" : "nopod", "nc=$NC_THRESH", 
    "win=$SIZE_SMALL_WINDOW");
print OUTPUT_FILE "\n\nPARAMETERS:\n--", join(" --", @tot_param);

close OUTPUT_FILE; 
close BED_FILE;    
close STATS_FILE if $STATS; 

my $temp = join("_", "$OUTPUT_DIR\/log_files\/temp.txt");
my $cleanup_script 
    = "grep -v 'Wanted\\|Parameter' $PERL_LOG > $temp\nmv $temp $PERL_LOG";
system($cleanup_script);
my $successful = $completed[0] + $completed[1];
print "Number of successful samples = $successful\n";
print "Number of failed samples = ", scalar(@files) - $successful, "\n"
    if (scalar(@files) - $successful);
print "Samples with no detectable abnormalities = $completed[1]\n"
    if ($completed[1]);
print STDERR "Number of successful samples = $successful\n";
print STDERR "Number of failed samples = ", 
    scalar(@files) - $successful, "\n";
print STDERR "Samples with no detectable abnormalities = $completed[1]\n";
if (!$successful && !$completed[1]) {
    print STDERR "ERROR: This analysis was unsuccessful!\n";
    die;
}
close STDERR; 
exit 0;



#THE END


################################# Sub Functions ################################

############################### Main Sub Function ##############################
sub main_program {
    my (@output, @stats, @trio_def);
    my $completed = 0;
    clear_global_variables();

    # Input file, tab delim, sorted by chromosome and position, 
    # format = SNP Name, Chromosome, Position, Father, Father BAF, Father LRR, 
    # Mother, Mother BAF, Mother LRR, Child, Child BAF, Child LRR
    $INPUT_FILE = $CURR_FILE;
    if(!open(INPUT_FILE, "<", $INPUT_FILE)) { 
        print "Could not open input file: $INPUT_FILE!\n";
        print STDERR "ERROR: Could not open input file: $INPUT_FILE!\n";            
        return();
    }
    
    # Get filename
    $FILENAME = (split(/[\/]/,(split(/.txt/, $CURR_FILE))[0]))[-1];

    # Get sample names
    my $HEADERS = <INPUT_FILE>;
    $P1_NAME = (split(/[.]/, (split(/\t/, $HEADERS))[3]))[0];
    $P2_NAME = (split(/[.]/, (split(/\t/, $HEADERS))[6]))[0];
    $CH_NAME = (split(/[.]/, (split(/\t/, $HEADERS))[9]))[0];

    # Create temp files for communication with R and call Rscript    
    if ($GRAPHICS) {
        ($R_PID_FILE, $R_PID_FILENAME)         = tempfile(UNLINK => 1);
        ($R_TO_PERL_FILE, $R_TO_PERL_FILENAME) = tempfile(UNLINK => 1);
        ($PERL_TO_R_FILE, $PERL_TO_R_FILENAME) = tempfile(UNLINK => 1);
        start_R();
        $R_ATTEMPTS = 0;
    }
    
    # Initial Calculations
    print "Performing Initial Calculations...\n\n" if $verbose;
    my $init_results = parse_file($INPUT_FILE, \&initial_calculations, 
        "check input");
    return() unless $init_results;
    my (@init_ch_stats, @init_p1_stats, @init_p2_stats);
    map {push(@init_ch_stats, $$_[0])} @$init_results;
    map {push(@init_p1_stats, $$_[1])} @$init_results;
    map {push(@init_p2_stats, $$_[2])} @$init_results;
    my @sorted_init_p1 = sort {
    $$a[0] <=> $$b[0] || $$a[1] <=> $$b[1]} @init_p1_stats;
    my @sorted_init_p2 = sort {
    $$a[0] <=> $$b[0] || $$a[1] <=> $$b[1]} @init_p2_stats;
    my @sorted_init_ch = sort {
    $$a[0] <=> $$b[0] || $$a[1] <=> $$b[1]} @init_ch_stats;

    (my $p1_normals, $P1_mBAF_MED, $P1_LRR_MED, my $p1_nc_rate, 
     $P1_AA_BOUND, $P1_BB_BOUND) = calc_init_stats(\@sorted_init_p1, 9);
    (my $p2_normals, $P2_mBAF_MED, $P2_LRR_MED, my $p2_nc_rate, 
     $P2_AA_BOUND, $P2_BB_BOUND) = calc_init_stats(\@sorted_init_p2, 9);
    (my $ch_normals, $CH_mBAF_MED, $CH_LRR_MED, my $ch_nc_rate, 
     $AA_BOUND, $BB_BOUND, $AB_UPPER_BOUND, $AB_LOWER_BOUND)
     = calc_init_stats(\@sorted_init_ch, 9);
        
    # Calculate MI1 thresholds
    my ($ab_sums, $ab_sum_sq, $ab_ct) = (0) x 3;
    foreach (@$ch_normals) {
        $ab_sums   += $$_[2]->[3];
        $ab_sum_sq += $$_[2]->[4];
        $ab_ct     += $$_[2]->[5];
    }    
    my @AB_stats = st_dev($ab_sums, $ab_sum_sq, $ab_ct);
    $MI1_UPPER_THRESH = $AB_stats[0] + ($AB_stats[1] * 5);
    $MI1_UPPER_THRESH = 0.95 if $MI1_UPPER_THRESH > 0.95;
    $MI1_LOWER_THRESH = $AB_stats[0] - ($AB_stats[1] * 5);
    $MI1_LOWER_THRESH = 0.05 if $MI1_LOWER_THRESH < 0.05;
    $PODCR_UPPER_THRESH = $AB_stats[0] + ($AB_stats[1] * 2);
    $PODCR_LOWER_THRESH = $AB_stats[0] - ($AB_stats[1] * 2);
   
    # Set initial estimate of acceptable errors
    $ACCEPTABLE_ERRORS = 3;

    # Determine largest number of informative SNPs in window for each arm
    my $max_windows = parse_file($INPUT_FILE, \&max_window);
    return() unless $max_windows;
    my @sorted_max_windows = sort {
        $$a[0] <=> $$b[0] || $$a[1] <=> $$b[1]} @$max_windows;
    map {$$_[2] ||= 5} @sorted_max_windows;
    for (0..$#sorted_init_ch) {
        push(@{$sorted_init_ch[$_]}, $sorted_max_windows[$_]->[2]);
    }

    ($ch_normals, $CH_mBAF_MED, $CH_LRR_MED, $ch_nc_rate, 
     $AA_BOUND, $BB_BOUND, $AB_UPPER_BOUND, $AB_LOWER_BOUND)
     = calc_init_stats(\@sorted_init_ch, -1);
    my $rounded_AA_bound     = sprintf("%.3f",$AA_BOUND);
    my $rounded_BB_bound     = sprintf("%.3f",$BB_BOUND);
    my $rounded_AB_up_bound  = sprintf("%.3f",$AB_UPPER_BOUND);
    my $rounded_AB_low_bound = sprintf("%.3f",$AB_LOWER_BOUND);

    # Recalculate MI1 thresholds
    ($ab_sums, $ab_sum_sq, $ab_ct) = (0) x 3;
    foreach (@$ch_normals) {
        $ab_sums   += $$_[2]->[3];
        $ab_sum_sq += $$_[2]->[4];
        $ab_ct     += $$_[2]->[5];
    }  
    
    @AB_stats = st_dev($ab_sums, $ab_sum_sq, $ab_ct);
    $MI1_UPPER_THRESH = $AB_stats[0] + ($AB_stats[1] * 5);
    $MI1_UPPER_THRESH = 0.95 if $MI1_UPPER_THRESH > 0.95;
    $MI1_LOWER_THRESH = $AB_stats[0] - ($AB_stats[1] * 5);
    $MI1_LOWER_THRESH = 0.05 if $MI1_LOWER_THRESH < 0.05;
    $PODCR_UPPER_THRESH = $AB_stats[0] + ($AB_stats[1] * 2);
    $PODCR_LOWER_THRESH = $AB_stats[0] - ($AB_stats[1] * 2);

    ###### Recalculate largest window sizes and get counts for HD, MI1, ########
    ###### and outliers, and pvalues for FDR calculations
    my $secondary_results = parse_file($INPUT_FILE, \&secondary_counts);
    return() unless $secondary_results;
    my @sorted_secondary = sort {
        $$a[0] <=> $$b[0] || $$a[1] <=> $$b[1]} @$secondary_results;
    for (0..$#sorted_secondary) {
        push(@{$sorted_secondary[$_]}, @{$sorted_init_ch[$_]}[6,7]);
    }

    my (@scan_stat_array, @normals);       
    # Make an array of refs for kmeans
    map {push(@scan_stat_array, [$$_[2]])} @sorted_secondary;
    # Cluster refined scan statistics
    my $y = 0.5;
    my ($num_clusters, $normal_cluster, $cluster_ref) 
        = cluster(\@scan_stat_array, $y);
    if ($num_clusters == 1) {@normals = @sorted_secondary}
    else {
        my @temp = @sorted_secondary;
        until ($num_clusters == 1) {
            @normals = @scan_stat_array = ();
            for (my $i = 0; $i < @temp; $i++) {
                if ($$cluster_ref[$i] == $normal_cluster) {
                    push(@normals, $temp[$i]);                
                }
            }
            map {push(@scan_stat_array, [$$_[2]])} @normals;
            @temp = @normals;
            #last if (@normals == 1); 
            $y -= 0.025;
            ($num_clusters, $normal_cluster, $cluster_ref) 
                = cluster(\@scan_stat_array, $y);
        }
    }
    # Remove any values in normal cluster > 4SDs of mean
    if (@normals > 1) {
        my ($stdev, $mean) = st_dev([map{$$_[2]} @normals]);
        my @adj_normals;
        for (0..$#normals) {
            my $curr_value = $normals[$_]->[2];
            unless (abs($curr_value) > $mean + $stdev * 4) {
                push (@adj_normals, $normals[$_]);
            }
        }
        @normals = @adj_normals;
    }

    my ($baf_snp_ct, $norm_baf_ct) = (0) x 2;
    $baf_snp_ct  += $_ for map {$$_[-2]} @sorted_secondary;
    $norm_baf_ct += $_ for map {$$_[-2]} @normals;
    my $adjusted_mi1_ct = 0;
    $adjusted_mi1_ct += $_ for map {$$_[3]} @normals;
    my $mi1_rate = $adjusted_mi1_ct / $norm_baf_ct;
    $MIN_MI1  = calc_min_adj_snps($norm_baf_ct, $mi1_rate, $ALPHA);
    
    my $total_mi1_ct;
    map {$total_mi1_ct += "$$_[3]\t"} @sorted_secondary;    
    my $total_mi1_rate = $total_mi1_ct / $baf_snp_ct;
    map {$HET_MI1_CT += "$$_[8]\t"} @normals; 
    $HET_MI1_RATE = $HET_MI1_CT / $norm_baf_ct;
    my $baf_outliers;
    $ab_ct = 0;
    map {$baf_outliers += "$$_[9]\t"} @normals;    
    map {$ab_ct += "$$_[10]\t"} @normals;    
    my $baf_outlier_rate = $baf_outliers / $ab_ct;
    $MIN_BAF_OUT = calc_min_adj_snps($ab_ct, $baf_outlier_rate, $ALPHA);

    my $next = check_quality($total_mi1_rate, $p1_nc_rate, $p2_nc_rate, 
        $ch_nc_rate); 
    return() if $next;

    # Calculate hd rates and minimum size of hd regions
    my ($p1_hd_ct, $p2_hd_ct, $ch_hd_ct, $norm_lrr_snp_ct) = (0) x 4;                
    $p1_hd_ct        += $_ for map {$$_[5]}  @normals;
    $p2_hd_ct        += $_ for map {$$_[6]}  @normals;
    $ch_hd_ct        += $_ for map {$$_[7]}  @normals;        
    $norm_lrr_snp_ct += $_ for map {$$_[-1]} @normals;
    my $p1_hd_rate  = $p1_hd_ct / $norm_lrr_snp_ct;
    my $p2_hd_rate  = $p2_hd_ct / $norm_lrr_snp_ct;
    my $ch_hd_rate  = $ch_hd_ct / $norm_lrr_snp_ct;
    $MIN_P1_HD      = calc_min_adj_snps($norm_lrr_snp_ct, $p1_hd_rate, $ALPHA);
    $MIN_P2_HD      = calc_min_adj_snps($norm_lrr_snp_ct, $p2_hd_rate, $ALPHA);
    $MIN_CH_HD      = calc_min_adj_snps($norm_lrr_snp_ct, $ch_hd_rate, $ALPHA);
    
    # Default assumption of error rate is 3.75 times the rate of single
    # Mendelian errors as determined by the ratio of all single 
    # errors derived from mutating all genotype combinations 
    # of normal inheritance to all detectable errors (24 detectable 
    # changes out of 90 possible changes). The ratio (90/24) is then used
    # to estimate the overall "error" rate (single Mendelian errors and 
    # single SNP biological abnormalities). In a previous step, single 
    # Mendelian errors were counted while discarding any which were 
    # adjacent to reduce influence of large biological abnormalities in
    # error rate estimations.
    $ERROR_RATE = $mi1_rate * 3.75;
    my $rounded_error_rate = sprintf("%.5f", $ERROR_RATE); 
    
    # Calculate minimum window size (informative SNPs) 
    my $inf_snp_ct;
    $inf_snp_ct += $_ for map {$$_[4]} @normals;
    my $inf_hUPD_ct;
    map {$inf_hUPD_ct += "$$_[11]\t"} @normals;
    unless (@normals) {print STDERR "cluster problem - no normals $FILENAME\n"}
    $hUPD_RATE = $inf_hUPD_ct / $inf_snp_ct;    
    my $inf_snp_rate = $inf_snp_ct / $norm_baf_ct;
    $inf_snp_ct = int(($inf_snp_rate * $baf_snp_ct) + 0.5);
    
    ($MIN_POD, $POD_ALPHA) = calculate_threshold($SIZE_SMALL_WINDOW, 6,
        $inf_snp_ct, $baf_snp_ct, $ALPHA);
    ($MIN_LARGE_POD, $LARGE_POD_ALPHA) = calculate_threshold(
        $SIZE_SMALL_WINDOW * 5, 12, $inf_snp_ct, $baf_snp_ct, $ALPHA);
        
    # If no possible POD regions given alpha, continue running
    # without POD analysis if other detection methods are selected
    if (!$MIN_POD) {
        if ($HD || $MI1 || $PODcr) {($POD, $MIN_POD) = ("",9)}
        else {
            print OUTPUT_FILE "\n\nNo regions of abnormal parental ",
                "contribution were detected at this alpha level.\n\n" 
                unless $BATCH;
            print STDERR "Sample $FILENAME : No regions of abnormal parental ",
                "contribution were detected at this alpha level.\n" if $BATCH;    
            print "\n\nNo regions of abnormal parental contribution were ",
                "detected at this alpha level.\n\n" if $verbose;
            return([0,0,0,1]);
        }
    }
  
    # Calculate Acceptable Errors 
    my ($stop, $errors) = (0) x 2;
    until ($stop) {
        my $pvalue = bin_test($errors, $SIZE_SMALL_WINDOW, $ERROR_RATE, "g");
        if ($pvalue <= 0.001) {
            $stop = 1;
            $ACCEPTABLE_ERRORS = $errors - 1;
            $ACCEPTABLE_ERRORS ||= 0;
        }
        $errors++;
    }
    
    ($stop, $errors) = (0) x 2;
    until ($stop) {
        my $pvalue = bin_test($errors, $SIZE_LARGE_WINDOW, $ERROR_RATE, "g");
        if ($pvalue <= 0.001) {
            $stop = 1;
            $LARGE_WIN_ACCEPTABLE_ERRORS = $errors - 1;
            $LARGE_WIN_ACCEPTABLE_ERRORS ||= 0;
        }
        $errors++;
    }

    if ($GRAPHICS) { 
        # Open file from R containing the current R PID
        get_R_PID();

        # Check if the R script is still running
        until(determine_R_status()) {
            if ($R_ATTEMPTS > 1) {
                print "A fatal error has occurred in communication with the R",
                    " software.\n" if $verbose;
                print STDERR "ERROR: A fatal error has occurred in ",
                    "communication with the R software.\n";  
                return();
            }        
        }
    }
    
    ################################## Analysis ################################
    my $analysis_results = parse_file($INPUT_FILE, \&manage_chromosome_arm,
        "main");
    my (@detected_regions, @local_meds, @local_meds_by_chr);
    for (@$analysis_results) {    
        push(@detected_regions, @{$$_[0]});
        # Store local medians
        if ($$_[1]) {
            push(@local_meds, $$_[2]);
            $local_meds_by_chr[$$_[1]] = $$_[2];
        }
    }
    
    # Refine autosomal medians to reflect only normal regions in child
    $P1_mBAF_MED = median([map {$$_[0]} @local_meds]);
    $P2_mBAF_MED = median([map {$$_[1]} @local_meds]);
    $CH_mBAF_MED = median([map {$$_[2]} @local_meds]);
    $P1_LRR_MED  = median([map {$$_[3]} @local_meds]);
    $P2_LRR_MED  = median([map {$$_[4]} @local_meds]);
    $CH_LRR_MED  = median([map {$$_[5]} @local_meds]);
    
    ####################### Descriptive Calculations #######################            
    # Check if regions were detected
    if ($#detected_regions < 0) {
        print OUTPUT_FILE "\n\nNo regions of abnormal parental contribution ",
            "were detected at this alpha level.\n\n" unless $BATCH;
        print STDERR "Sample $FILENAME : No regions of abnormal parental ",
            "contribution were detected at this alpha level.\n" if $BATCH;    
        print "\n\nNo regions of abnormal parental contribution were ",
            "detected at this alpha level.\n\n" if $verbose;
        return([0,0,0,1]);
    }
   
    # Sort output array by chromosome and position     
    my @sorted_detected_regions = sort {
        $$a[0] <=> $$b[0] || $$a[1] <=> $$b[1]} @detected_regions;

    ######################## Prepare results for output file ####################
            
    my ($count, $num_pod, $num_mi1, $num_hd, $num_podcr) = (0) x 5;    
    foreach my $ref (@sorted_detected_regions) {
        if    ($$ref[17] == POD)   {$num_pod++}
        elsif ($$ref[17] == MI1)   {$num_mi1++}
        elsif ($$ref[17] == HD)    {$num_hd++}
        elsif ($$ref[17] == PODcr) {$num_podcr++}            
        
        # If a region comprised most of a chromosome arm and underwent
        # a single round of analyses, adjust the regions medians by
        # the medians of the normal SNPs from the other arm, or if 
        # most of the chromosome is abnormal, adjust using the autosomal
        # medians.
        my ($p1_mBAF, $p2_mBAF, $ch_mBAF, $p1_LRR, $p2_LRR, $ch_LRR);
        if ($$ref[19] == 1) {
            if ($local_meds_by_chr[$$ref[0]]) {
                # local medians 
                ($p1_mBAF, $p2_mBAF, $ch_mBAF, $p1_LRR, $p2_LRR, $ch_LRR) = 
                    @{$local_meds_by_chr[$$ref[0]]};
            }
            else {
                # refined autosomal medians
                ($p1_mBAF, $p2_mBAF, $ch_mBAF, $p1_LRR, $p2_LRR, $ch_LRR) = 
                    ($P1_mBAF_MED, $P2_mBAF_MED, $CH_mBAF_MED,
                    $P1_LRR_MED, $P2_LRR_MED, $CH_LRR_MED);
            }            
            $$ref[10] = sprintf("%.4f", $$ref[10] - ($ch_mBAF - 0.5))
                unless $$ref[10] eq "NA";
            $$ref[11] = sprintf("%.4f", $$ref[11] - $ch_LRR)
                unless $$ref[11] eq "NA";
            $$ref[12] = sprintf("%.4f", $$ref[12] - ($p1_mBAF - 0.5))
                unless $$ref[12] eq "NA";
            $$ref[13] = sprintf("%.4f", $$ref[13] - $p1_LRR)
                unless $$ref[13] eq "NA";
            $$ref[14] = sprintf("%.4f", $$ref[14] - ($p2_mBAF - 0.5))
                unless $$ref[14] eq "NA";
            $$ref[15] = sprintf("%.4f", $$ref[15] - $p2_LRR)
                unless $$ref[15] eq "NA";
            if ($$ref[11] > $AMP + $ch_LRR) {$$ref[3] = "AMP"}
            elsif ($$ref[11] < $DEL_UPPER + $ch_LRR && 
                $$ref[11] > $HD_HI + $ch_LRR) {$$ref[3] = "DEL"}
            elsif ($$ref[11] <= $HD_HI + $ch_LRR) {$$ref[3] = "HD"}
            #else {$$ref[3] = " "}
            $$ref[20] = call_inheritance($ref);
            unless($$ref[10] eq "NA") {
                $$ref[10] = sprintf("%.4f", 0.5) if ($$ref[10] < 0.5)}
            unless($$ref[12] eq "NA") {
                $$ref[12] = sprintf("%.4f", 0.5) if ($$ref[12] < 0.5)}
            unless($$ref[14] eq "NA") {
                $$ref[14] = sprintf("%.4f", 0.5) if ($$ref[14] < 0.5)}
        }

        if    ($$ref[4] == FATHER)  {$$ref[4] = "Father"}
        elsif ($$ref[4] == MOTHER)  {$$ref[4] = "Mother"}
        elsif ($$ref[4] == NONE)    {$$ref[4] = "None"}
        elsif ($$ref[4] == BOTH)    {$$ref[4] = "Both"}        
        elsif ($$ref[4] == UNKNOWN) {$$ref[4] = "Unknown"}
        if ($$ref[3] eq " " && ($$ref[4] eq "Father" || $$ref[4] eq "Mother")) {
            $$ref[4] = $$ref[4] . "(C)";
        }
        
        my $detection;         
        if    ($$ref[17] == HD)    {$detection = "HD"}
        elsif ($$ref[17] == MI1)   {$detection = "MI1"}
        elsif ($$ref[17] == POD)   {$detection = "POD"}
        elsif ($$ref[17] == PODcr) {$detection = "PODcr"}            
        
        push(@output, [$CH_NAME, @$ref[0..4,20,5..7], $detection, 
            @$ref[10..15]]);                    
    }
    push(@trio_def, "Trio Definition Autosomal NoCall Rate",
        "Father =   $P1_NAME\t$p1_nc_rate",
        "Mother =   $P2_NAME\t$p2_nc_rate",
        "Child  =   $CH_NAME\t$ch_nc_rate") unless ($BATCH);        
            
    ############## Create file to pass results to R script ##################
    if ($GRAPHICS) { 
        print $PERL_TO_R_FILE "Chr\tStart\tStop\tType\tParent\tSNPs",
            "\tInformative SNPs\tSize of Region (bp)\tRadius of region for R",
            "\tbp of midpoint for R\tMedian BAF\tMedian LRR\tFather-Median ",
            "BAF\tFather-Median LRR\tMother-Median BAF\tMother-Median LRR\n";

        map {print $PERL_TO_R_FILE join("\t", @{$_}[0..15]), "\n"} 
            @sorted_detected_regions;
        print "Creating Graphics\n" if $verbose;
        
        until(determine_R_status()) {
            if ($R_ATTEMPTS > 1) {
                print "A fatal error has occurred in communication with the R",
                    " software.\n" if $verbose;
                print STDERR "ERROR: A fatal error has occurred in ",
                    "communication with the R software.\n";  
                return();
            }        
        }
        
        # Wait for R script to finish and return results
        my $waiting_for_R = 0;
        until ($waiting_for_R) {
            if (-s $R_TO_PERL_FILENAME) {$waiting_for_R = 1}
            else {
                until(determine_R_status()) {
                    if ($R_ATTEMPTS > 1) {
                        print "A fatal error has occurred in communication ",
                            "with the R software.\n" if $verbose;
                        print STDERR "A fatal error has occurred in ",
                            "communication with the R software.\n";  
                        return();
                    }        
                }
                usleep(500000);
            }
        }
    }  
    map{$$_[4] = "NA" if $$_[4] eq " "} @output;
    map{$$_[6] = "NA" if $$_[6] eq " "} @output;    
    if ($STATS) {
        push(@stats, [$FILENAME, $CH_NAME, $P1_NAME, $P2_NAME, $ch_nc_rate,
        $p1_nc_rate, $p2_nc_rate, $ch_hd_rate, $MIN_CH_HD, $p1_hd_rate, 
        $MIN_P1_HD, $p2_hd_rate, $MIN_P2_HD, $mi1_rate, $MIN_MI1, 
        $MI1_UPPER_THRESH, $MI1_LOWER_THRESH, $MIN_POD, $POD_ALPHA, $ALPHA, 
        $ERROR_RATE, $ACCEPTABLE_ERRORS, $AA_BOUND, $BB_BOUND, 
        $AB_UPPER_BOUND, $AB_LOWER_BOUND, $#sorted_detected_regions, 
        $num_pod, $num_podcr, $num_mi1, $num_hd]);
    }
    $completed = 1; 
    return([\@output, \@stats, \@trio_def, $completed]);
}

########################### Sub Functions ###########################

sub baf_stats_for_stdev {
    my ($baf, $lower, $upper, $baf_ref, $mbaf_ref) = @_;
    my ($aa_sum, $aa_sum_squares, $aa_count, $ab_sum, $ab_sum_squares,
        $ab_count, $bb_sum, $bb_sum_squares, $bb_count, $baf_sum, 
        $baf_ct);
        
    # $baf_ref array indices = aa_sum, aa_sum_squares, aa_count, ab_sum, 
    # ab_sum_squares, ab_count, bb_sum, bb_sum_squares, bb_count, 
    # baf sum for mbaf, baf ct for mbaf
    if ($baf < $upper) {
        if ($baf <= $lower) {
            # AA
            $$baf_ref[0] += $baf;
            $$baf_ref[1] += $baf**2;
            $$baf_ref[2]++;
        }
        else {
            # AB
            $$baf_ref[3] += $baf;
            $$baf_ref[4] += $baf**2;
            $$baf_ref[5]++;
            push(@$mbaf_ref, $baf);
            $$baf_ref[9] += $baf;
            $$baf_ref[10]++;
        }
    }
    else {
        # BB
        $$baf_ref[6] += $baf;
        $$baf_ref[7] += $baf**2;
        $$baf_ref[8]++;
    }
}

sub bin_coeff {
    # Estimates and returns the binomial coefficient for n and k.  
    my ($n, $k) = @_;
    return 1 if $k == 0;
    my $r = $n - $k + 1;
    for (2..$k) {$r *= (($n - $k)/$_ + 1)}
    return ($r);
}

sub bin_test {
    # Performs an approximation of the binomial test.
    # alternative hypothesis = "t" = two-tailed, "l" = one-tailed <= k,
    # "g" = one-tailed >= k
    my ($k, $n, $p, $tail) = @_;
    my ($pvalue, $x, $y) = (0) x 3;
    if ($tail eq "t") {
        if ($n / 2 == $k) {return(1)}
        else {($x, $y) = ($k <= $n / 2) ? (0, $k) : ($k, $n)}
    }
    elsif ($tail eq "l") {($x, $y) = (0, $k)}
    elsif ($tail eq "g") {($x, $y) = ($k, $n)}
    for ($x..$y) {$pvalue += dbinom($_, $n, $p)}
    $pvalue = 2 * $pvalue if ($tail eq "t");
    $pvalue = 1 if ($pvalue > 1);
    return($pvalue);
}

sub calc_init_stats {
    my ($init_stats, $index) = @_;
    # Identify and ignore abnormal chromosome arms
    my (@values, @normals) = ();
    map {push(@values, [$$_[$index]])} @$init_stats;
    my $y = 0.5;
    my ($num_clusters, $normal_cluster, $cluster_ref) 
        = cluster(\@values, $y);
    my $ct = 0;    
    if ($num_clusters == 1) {@normals = @$init_stats}
    else {
        my @temp = @$init_stats;
        until ($num_clusters == 1) {
            @normals = @values = ();
            for (my $i = 0; $i < @temp; $i++) {
                if ($$cluster_ref[$i] == $normal_cluster) {
                    push(@normals, $temp[$i]);                
                }
            }
            map {push(@values, [$$_[$index]])} @normals;
            @temp = @normals;
            #last if (@normals == 1); 
            $ct++;
            $y -= 0.025 * $ct;
            ($num_clusters, $normal_cluster, $cluster_ref) 
                = cluster(\@values, $y);
        }
    }
    # Remove any values in normal cluster > 4SDs of mean
    if (@normals > 1) {
        my ($stdev, $mean) = st_dev([map{$$_[$index]} @normals]);
        my @adj_normals;
        for (0..$#normals) {
            my $curr_value = $normals[$_]->[$index];
            unless (abs($curr_value) > $mean + $stdev * 4) {
                push (@adj_normals, $normals[$_]);
            }
        }
        @normals = @adj_normals;
    }

    # Calculate median of mBAF medians of "normal" chromosome arms
    my $mBAF_med    = median([map{$$_[3]} @normals]);
    my $lrr_med     = median([map{$$_[4]} @normals]);
    # Calculate NC rate of "normal" chromosome arms                
    my ($nc_ct, $nc_snp_ct, $nc_rate) = (0) x 3;
    $nc_ct      += $_ for map {$$_[5]} @normals;
    $nc_snp_ct  += $_ for map {$$_[8]} @normals; 
    $nc_rate = sprintf("%.4f", $nc_ct / $nc_snp_ct) if $nc_snp_ct; 
    $nc_rate ||= 0;

    # Calculate BAF thresholds
    my ($aa_sums, $aa_sum_sq, $aa_ct) = (0) x 3;
    foreach (@normals) {
        $aa_sums   += $$_[2]->[0];
        $aa_sum_sq += $$_[2]->[1];
        $aa_ct     += $$_[2]->[2];
    }    
    my @AA_stats = st_dev($aa_sums, $aa_sum_sq, $aa_ct);
    my ($ab_sums, $ab_sum_sq, $ab_ct) = (0) x 3;
    foreach (@normals) {
        $ab_sums   += $$_[2]->[3];
        $ab_sum_sq += $$_[2]->[4];
        $ab_ct     += $$_[2]->[5];
    }    
    my @AB_stats = st_dev($ab_sums, $ab_sum_sq, $ab_ct);
    my ($bb_sums, $bb_sum_sq, $bb_ct) = (0) x 3;
    foreach (@normals) {
        $bb_sums   += $$_[2]->[6];
        $bb_sum_sq += $$_[2]->[7];
        $bb_ct     += $$_[2]->[8];
    }    
    my @BB_stats = st_dev($bb_sums, $bb_sum_sq, $bb_ct);

    my $aa_bound = $AA_stats[0] + ($AA_stats[1] * $HOM_SD);
    my $bb_bound = $BB_stats[0] - ($BB_stats[1] * $HOM_SD);
    my $ab_upper = $AB_stats[0] + ($AB_stats[1] * $HET_SD);
    my $ab_lower = $AB_stats[0] - ($AB_stats[1] * $HET_SD);
    
    return(\@normals, $mBAF_med, $lrr_med, $nc_rate,
            $aa_bound, $bb_bound, $ab_upper, $ab_lower);
}

sub calc_min_adj_snps {
    # Calculate minimum number of SNPs, given a significance 
    # threshold and estimated rate.
    # Iterate until the number SNPs < the significance threshold. 
    my ($n, $p, $thresh) = @_;
    my $k = 2;
    my $min_SNPs;
    until ($min_SNPs) {
        my $pvalue = streak_pvalue($n, $k, $p);
        if ($pvalue <= $thresh) {$min_SNPs = $k}
        $k++;
    }
    return($min_SNPs);
}
  
sub calculate_expected_number {
    # E(R) = (N - k + 1) sum(i=k+1 to N) of b(i;N,w/M)
    # where k = informative SNPs in a window, w = window size  
    # N = total informative SNPs, M = total SNPs, $n = trials,
    # E(R) = the number of expected windows with k informative SNPs    
    # The binomial probability of occurrence of clump sizes have been 
    # previously calculated and stored in p_array
    my ($k, $N, $p_array_ref) = @_;
    my $sum = 0;
    map {$sum += $_} @$p_array_ref[$k - 1..@$p_array_ref - 1];
    my $expected = ($N - $k + 1) * $sum;
    return($expected);    
}
      
sub calculate_region_median {
    my ($region_ref, $input, $index) = @_;    
    my ($array_ref, $hash_ref, @median_array);
    my $count = 0;
    my ($median, $lowest_index, $highest_index, $lowest_value, $highest_value);
    my ($aa_bound, $bb_bound) = ($AA_BOUND, $BB_BOUND);
    
    if ($input == LRR) {($array_ref, $hash_ref) = (\@LRR_REFS, \%LRR_BY_POS)}
    elsif ($input == BAF) {
        ($array_ref, $hash_ref) = (\@CURR_CHR_REFS, \%BAF_BY_POS);
        if    ($index == 4) {($aa_bound,$bb_bound)=($P1_AA_BOUND,$P1_BB_BOUND)}
        elsif ($index == 7) {($aa_bound,$bb_bound)=($P2_AA_BOUND,$P2_BB_BOUND)}
    }
    
    my ($start_index, $stop_index, $present) 
        = find_closest_indices(@$region_ref[START,STOP], $hash_ref);
    if ($present) {
        for ($start_index..$stop_index) {
            if ($input == LRR || $$array_ref[$_]->[$index] >= $aa_bound &&
                $$array_ref[$_]->[$index] <= $bb_bound) {
                if ($_ != $start_index && $_ != $stop_index) {
                    if (!defined($lowest_value)) {
                        $lowest_value = $$array_ref[$_]->[$index];
                        $lowest_index = $_;
                    }
                    elsif ($$array_ref[$_]->[$index] < $lowest_value) {
                        $lowest_value = $$array_ref[$_]->[$index];
                        $lowest_index = $_;
                    }
                    if (!defined($highest_value)) {
                        $highest_value = $$array_ref[$_]->[$index]; 
                        $highest_index = $_;
                    }
                    elsif ($$array_ref[$_]->[$index] > $highest_value) {
                        $highest_value = $$array_ref[$_]->[$index];                
                        $highest_index = $_;
                    }
                }
                if ($input == BAF && $$array_ref[$_]->[$index] < 0.5) {
                    push(@median_array, (1 - $$array_ref[$_]->[$index]));
                }
                else {push(@median_array, $$array_ref[$_]->[$index]);
                }
            }
        }
        if (@median_array >= 2) {$median = median(\@median_array)}
        elsif ($input == BAF) {$median = 1}   
        else {$median = NA}
    }
    else {$median = NA}
    return([$median, $start_index, $stop_index, $lowest_index, $highest_index,
    $lowest_value, $highest_value]);
}
   
sub calculate_region_stats {
    # [0..10] = Chr, Start, Stop, Type, Parent, #SNPs, #Informative SNPs, 
    # Size(Mb), Radius, Midpoint 
    my $region_ref = $_[0];
    my @stats = @$region_ref;
    # Calculate number of SNPs in region    
    $stats[5] = count_SNPs(@stats[START,STOP], \%SNP_BY_POS);
    # Calculate size(bp) of region
    $stats[7] = $stats[STOP] - $stats[START];
    # Calculate radius for graphics
    $stats[8] = $stats[7] / 2;
    # Calculate midpoint for graphics
    $stats[9] = $stats[STOP] - $stats[8];
             
    # Calculate parental and progeny BAF and LRR medians for a 
    # detected region.            
    my $p1_lrr_median  = calculate_region_median($region_ref, LRR, 3)->[0];
    my $p2_lrr_median  = calculate_region_median($region_ref, LRR, 4)->[0];
    my $ch_lrr_median  = calculate_region_median($region_ref, LRR, 5)->[0];
    my $p1_mbaf_median = calculate_region_median($region_ref, BAF, 4)->[0];
    my $p2_mbaf_median = calculate_region_median($region_ref, BAF, 7)->[0];
    my $ch_mbaf_median = calculate_region_median($region_ref, BAF, 10)->[0];

    if($REFINING) {
        if ($ch_mbaf_median == 1){$stats[10] = sprintf("%.4g",$ch_mbaf_median)}
        else {$stats[10] = sprintf("%.4g", $ch_mbaf_median 
            - ($LOCAL_CH_mBAF_MED - 0.5))}
        $stats[11] = sprintf("%.4g", $ch_lrr_median - $LOCAL_CH_LRR_MED);
        if ($p1_mbaf_median == 1){$stats[12] = sprintf("%.4g",$p1_mbaf_median)}
        else {$stats[12] = sprintf("%.4g", $p1_mbaf_median 
            - ($LOCAL_P1_mBAF_MED - 0.5))}
        $stats[13] = sprintf("%.4g", $p1_lrr_median - $LOCAL_P1_LRR_MED);
        if ($p2_mbaf_median == 1){$stats[14] = sprintf("%.4g",$p2_mbaf_median)}
        else {$stats[14] = sprintf("%.4g", $p2_mbaf_median 
            - ($LOCAL_P2_mBAF_MED - 0.5))}
        $stats[15] = sprintf("%.4g", $p2_lrr_median - $LOCAL_P2_LRR_MED);
    }
    else {
        $stats[10] = sprintf("%.4g", $ch_mbaf_median);
        $stats[11] = sprintf("%.4g", $ch_lrr_median);
        $stats[12] = sprintf("%.4g", $p1_mbaf_median);
        $stats[13] = sprintf("%.4g", $p1_lrr_median);
        $stats[14] = sprintf("%.4g", $p2_mbaf_median);
        $stats[15] = sprintf("%.4g", $p2_lrr_median);
    }
    
    if ($stats[10] != NA && $stats[10] < 0.5) {$stats[10] = sprintf("%.4f", 0.5)}
    if ($stats[12] != NA && $stats[12] < 0.5) {$stats[12] = sprintf("%.4f", 0.5)}
    if ($stats[14] != NA && $stats[14] < 0.5) {$stats[14] = sprintf("%.4f", 0.5)}
    if ($ch_lrr_median < $HD_LOW || $ch_mbaf_median == NA) {$stats[10] = "NA"}
    if ($p1_lrr_median < $HD_LOW || $p1_mbaf_median == NA) {$stats[12] = "NA"}
    if ($p2_lrr_median < $HD_LOW || $p2_mbaf_median == NA) {$stats[14] = "NA"}
    
    if ($stats[11] >= $AMP) {$stats[3] = "AMP"}
    elsif ($stats[11] <= $DEL_UPPER && $stats[11] > $HD_HI) {$stats[3] = "DEL"}
    elsif ($stats[11] <= $HD_HI && $stats[11] != NA) {$stats[3] = "HD"}
    else {$stats[3] = " "}
    if ($ch_lrr_median == NA) {$stats[11] = "NA"}
    if ($p1_lrr_median == NA) {$stats[13] = "NA"}
    if ($p2_lrr_median == NA) {$stats[15] = "NA"}
    @$region_ref = @stats;
}

sub calculate_threshold {
    # Determine minimum POD region size and the adjusted alpha
    my ($w, $k, $N, $M, $alpha) = @_;
    my ($pvalue, $prob, $done) = (0) x 3;
    my @p_array = (0) x 6; 
    for ($k - 1..$N) {
        my $bin_coeff = bin_coeff($N, $_);
        $p_array[$_] = $bin_coeff * (($w / $M)**$_) 
            * ((1 - ($w / $M))**($N - $_));
        last if ($p_array[$_] == 0);
    }
    until ($done || $k > $w) {
        my $expected = int(calculate_expected_number($k, $N, \@p_array)+0.5);
        my $beta = bin_test($k, $k, 0.5, "t");
        my $gamma_k;
        for ($k / 2..$k) {
            $pvalue = bin_test($_, $k, 0.5, "t");
            $gamma_k = $pvalue if $pvalue <= $beta;
            last if $pvalue < $beta;
        }
        $prob = bin_test(1, $expected, $gamma_k, "g");
        $k++;
        next unless ($prob < $alpha);
        for (my $i = $k; $i <= $w; $i++) {
            my $expected = calculate_expected_number($i, $N, \@p_array);
            last if $expected < 1;
            $expected = int($expected + 0.5);
            my $gamma_k;
            for ($i / 2..$i) {
                $pvalue = bin_test($_, $i, 0.5, "t");
                $gamma_k = $pvalue if $pvalue < $beta;
                last if $pvalue < $beta;
            }
            my $test = bin_test(1, $expected, $gamma_k, "g");
            $prob += $test;
            last unless ($prob < $alpha);    
        }
        $done = 1 if ($prob < $alpha);
    }
    return($k - 1, $pvalue);
}  

sub call_inheritance {
    # Determines whether a likely inheritance pattern can be called.
    # Calls are made for two categories : inherited (INH) and 
    # inherited with novel copy number state (INH_CN). The absence of a called
    # inheritance pattern does not imply a de novo abnormality, but 
    # includes de novo, UPD, low-level mosaicism, and regions which
    # the artificial thresholds have failed to determine the true pattern.    
    my $region = $_[0];
    my $ch_AMP = $$region[3] eq "AMP";
    my $ch_DEL = $$region[3] eq "DEL";
    my $ch_HD  = $$region[3] eq "HD";
    
    if ($$region[3] eq "DEL" && !$$region[19]) {
        if ($$region[4] == FATHER) {$$region[4] = MOTHER}
        elsif ($$region[4] == MOTHER) {$$region[4] = FATHER}
    }
    my $p1_AMP = $$region[13] > $AMP unless $$region[13] eq "NA";
    my $p2_AMP = $$region[15] > $AMP unless $$region[15] eq "NA";
    my $p1_HD  = $$region[13] <= $HD_HI unless $$region[13] eq "NA";
    my $p2_HD  = $$region[15] <= $HD_HI unless $$region[15] eq "NA";
    my $p1_DEL = $$region[13] < $DEL_UPPER && $$region[13] > $HD_HI
        unless $$region[13] eq "NA";
    my $p2_DEL = $$region[15] < $DEL_UPPER && $$region[15] > $HD_HI
        unless $$region[15] eq "NA";
    if ($$region[4] == UNKNOWN) {
        if ($ch_AMP) {
            if ($p1_AMP && !$p2_AMP)    {$$region[4] = FATHER}
            elsif ($p2_AMP && !$p1_AMP) {$$region[4] = MOTHER}
        }
        if ($ch_DEL) {
            if ($p1_DEL && $p2_DEL)    {$$region[4] = FATHER}
            elsif ($p2_DEL && $p1_DEL) {$$region[4] = MOTHER}
        }
    }
    my $P1   = $$region[4] == FATHER;
    my $P2   = $$region[4] == MOTHER;
    my $NONE = $$region[4] == NONE;
    
    # Likely inherited 
    my $inh    = $ch_AMP && ($p1_AMP && $P1 || $p2_AMP && $P2)
              || $ch_DEL && (($p1_DEL || $p1_HD) && $P1
                  || ($p2_DEL || $p2_HD) && $P2)    
              || $ch_HD && (($p1_HD && ($p2_HD || $p2_DEL))  
                  || ($p2_HD && ($p1_HD || $p1_DEL)));
    
    # Likely inherited with unique copy number state
    my $inh_cn = $ch_DEL && (!$p1_DEL && !$p1_HD && $p2_HD && $P2
                          || !$p2_DEL && !$p2_HD && $p1_HD && $P1)
              || $ch_HD && $p1_DEL && $p2_DEL;
              
    if ($ch_HD && (($p1_DEL || $p1_HD)
        && ($p2_DEL || $p2_HD))) {$$region[4] = BOTH}
    if ($$region[4] == NONE) {$$region[4] = UNKNOWN}
    
    my $state = " ";
    $state = "INH" if ($inh);
    $state = "INH-CN" if ($inh_cn);
    return($state);                        
}              

sub check_input_format {
    # Perform checks for input formatting errors
    my $array_ref = $_[0];
    my $prev_position = 0;
    for (@$array_ref) {
        my @line_array = @$_;
        # Check Chromosome column for anything other than acceptable
        # chromosomes
        if ($line_array[CHR] 
            !~ /^[1-9]$|^[1][0-9]|^[2][0-2]$|x|y|xy|m|mt/i) {
            print "\nSample $FILENAME : FAILED - Unrecognized ",
                "character ($line_array[CHR]) in Chromosome column!\n";
            print STDERR "ERROR: Sample $FILENAME : FAILED - Unrecognized ",
                "character ($line_array[CHR]) in Chromosome column!\n";
            print_proper_format();
            return();
        }
        # Check Position column for non-numeric characters and negative
        # numbers  
        elsif ($line_array[POS] !~ /^-?\d/ || $line_array[POS] < 0) {
            print "\nSample $FILENAME : FAILED - Unrecognized number ",
                "($line_array[POS]) in Position column!\n";
            print STDERR "ERROR: Sample $FILENAME : FAILED - Unrecognized ",
                "number ($line_array[POS]) in Position column!\n";                    
            print_proper_format();
            return();
        }
        # Check Genotype columns for proper genotype format
        elsif ($line_array[3] !~ /AA|AB|BB|^N/ || 
               $line_array[6] !~ /AA|AB|BB|^N/ || 
               $line_array[9] !~ /AA|AB|BB|^N/) {
            print "\nSample $FILENAME : FAILED - Unrecognized ",
                "character (@line_array[3,6,9]) in a Genotype column!\n";
            print STDERR "ERROR: Sample $FILENAME : FAILED - Unrecognized ",
                "character (@line_array[3,6,9]) in a Genotype ",
                "column!\n";                    
            print_proper_format();
            return();
        }
        # Check BAF columns for numbers < 0 or > 1 or 2 for cnv SNP
        elsif ((($line_array[4]  < 0 || ($line_array[4] > 1
               && $line_array[4] != 2) || ($line_array[4] == 2
               && substr($line_array[0], 0, 2) !~ /cn/i)) 
               && $line_array[4] ne "NaN") ||
               (($line_array[7]  < 0 || ($line_array[7] > 1
               && $line_array[7] != 2) || ($line_array[7] == 2
               && substr($line_array[0], 0, 2) !~ /cn/i)) 
               && $line_array[7] ne "NaN") ||
               (($line_array[10]  < 0 || ($line_array[10] > 1
               && $line_array[10] != 2) || ($line_array[10] == 2
               && substr($line_array[0], 0, 2) !~ /cn/i)) 
               && $line_array[10] ne "NaN")) {
            print "\nSample $FILENAME : FAILED - B allele frequency ",
                "is out of range (@line_array[4,7,10]) !\n";
            print STDERR "ERROR: Sample $FILENAME : FAILED - B allele ",
                "frequency is out of range (@line_array[4,7,10]) !\n";                    
            print_proper_format();
            return();
        }
        # Check for input in last expected column
        elsif (!defined($line_array[11])) {
            print "\nSample $FILENAME : FAILED - Too few columns.\n";
            print STDERR "ERROR: Sample $FILENAME : FAILED - Too few ",
                "columns.\n";                    
            print_proper_format();
            return();
        }
        # Check if Position column is properly sorted
        elsif ($prev_position > $line_array[POS]) {
            print "\nSample $FILENAME : FAILED - The input file ",
                "is not properly sorted!\n";
            print STDERR "ERROR: Sample $FILENAME : FAILED - The input ",
                "file is not properly sorted!\n";                    
            print_proper_format();
            return();
        }
        $prev_position = $line_array[POS];
    }
    return(1);
}

sub check_quality {
    # Ends the current analysis if quality checks are not met.
    my ($mi1_rate, $p1_nc_rate, $p2_nc_rate, $ch_nc_rate) = @_;    
    my ($next, $alerts) = (0) x 2;
    # Analysis cannot continue if estimated non-adjacent single 
    # Mendelian error rate is above 2%.
    if ($mi1_rate >= 0.02) {
        $next = 1;
        print "FAILED - The Mendelian error rate (", 
            sprintf("%.3f",$mi1_rate), " %) is above acceptable levels,\n",
            "indicating a relationship annotation or quality problem.",
            "\nThe analysis cannot be completed!\n" if $verbose;
        print STDERR "ERROR: Sample $FILENAME : FAILED - The Mendelian ",
            "error rate (", sprintf("%.3f",$mi1_rate), " %) is above ",
            "acceptable levels, indicating a relationship annotation or ",
            "quality problem.\n";
        clean_up() unless $BATCH;
    }
  
    # Next file if NC rate is above the NC threshold       
    if ($p1_nc_rate >= $NC_THRESH) {
        $next = 1;
        print "FAILED - The Father's NoCall rate ($p1_nc_rate) is greater ",
            "than the NoCall threshold ($NC_THRESH),\nindicating a quality ",
            "problem. The analysis cannot be completed!\nRaising the NoCall ",
            "threshold may result in increased false positive calls!\n\n" 
            if $verbose;
        print STDERR "ERROR: Sample $FILENAME : FAILED - The Father's NoCall ",
            "rate ($p1_nc_rate) is greater than the NoCall threshold ",
            "($NC_THRESH), indicating a quality problem!\n",
            "The analysis cannot be completed!\nRaising the NoCall threshold ",
            "may result in increased false positive calls!\n\n";
        clean_up() unless $BATCH;            
    }   
    if ($p2_nc_rate >= $NC_THRESH) {
        $next = 1;
        print "FAILED - The Mother's NoCall rate ($p2_nc_rate) is greater ",
            "than the NoCall threshold ($NC_THRESH),\nindicating a quality ",
            "problem. The analysis cannot be completed!\nRaising the NoCall ",
            "threshold may result in increased false positive calls!\n\n" 
            if $verbose;
        print STDERR "ERROR: Sample $FILENAME : FAILED - The Mother's NoCall ",
            "rate ($p2_nc_rate) is greater than the NoCall threshold ",
            "($NC_THRESH), indicating a quality problem!\n",
            "The analysis cannot be completed!\nRaising the NoCall threshold ",
            "may result in increased false positive calls!\n\n";
        clean_up() unless $BATCH;
    }   
    if ($ch_nc_rate >= $NC_THRESH) {
        $next = 1;
        print "FAILED - The Child's NoCall rate ($ch_nc_rate) is greater ",
            "than the NoCall threshold ($NC_THRESH),\nindicating a quality ",
            "problem. The analysis cannot be completed!\nRaising the NoCall ",
            "threshold may result in increased false positive calls!\n\n" 
            if $verbose;
        print STDERR "ERROR: Sample $FILENAME : FAILED - The Child's NoCall ",
            "rate ($ch_nc_rate) is greater than the NoCall threshold ",
            "($NC_THRESH), indicating a quality problem!\n",
            "The analysis cannot be completed!\nRaising the NoCall threshold ",
            "may result in increased false positive calls!\n\n";
        clean_up() unless $BATCH;
    }   
    return($next); 
}

sub clean_up {
    # Removes files if analysis was not completed due to quality checks
    my $cmd = "rm $OUTPUT_FILE $STATS_FILE";
    system($cmd);
}

sub clear_global_variables {
    ($AA_BOUND, $AB_LOWER_BOUND, $AB_UPPER_BOUND, $ACCEPTABLE_ERRORS, 
     $BB_BOUND, $BOUNDARY_EXTENSION, $CH_LRR_MED, $CH_mBAF_MED, $CH_NAME, 
     $ERROR_RATE, $FILENAME,$HET_MI1_CT,$HET_MI1_RATE,$hUPD_RATE,$INPUT_FILE,
     $LARGE_POD_ALPHA, $LARGE_WIN_ACCEPTABLE_ERRORS, $LOCAL_CH_BAF_MEAN,
     $LOCAL_CH_mBAF_MED,$LOCAL_CH_LRR_MED,$LOCAL_P1_mBAF_MED,$LOCAL_P1_LRR_MED,
     $LOCAL_P2_mBAF_MED, $LOCAL_P2_LRR_MED,$MI1_LOWER_THRESH,$MI1_UPPER_THRESH, 
     $MIN_BAF_OUT, $MIN_CH_HD, $MIN_hUPD, $MIN_LARGE_POD, $MIN_MI1, $MIN_P1_HD,$MIN_P2_HD, 
     $MIN_POD, $P1_AA_BOUND,$P1_BB_BOUND,$P2_AA_BOUND,$P2_BB_BOUND,$P1_LRR_MED, 
     $P1_mBAF_MED,$P1_NAME,$P2_LRR_MED, $P2_mBAF_MED,$P2_NAME,$PERL_TO_R_FILE, 
     $PERL_TO_R_FILENAME, $POD_ALPHA, $PODCR_LOWER_THRESH, $PODCR_UPPER_THRESH, 
     $R_ATTEMPTS, $R_PID, $R_PID_FILE, $R_PID_FILENAME, $R_TO_PERL_FILE, 
     $R_TO_PERL_FILENAME, $REFINING) = (0) x 56;
}

sub cluster {
    # Returns optimal number of clusters (1 or 2) and array of cluster results
    my ($data_ref, $y) = @_; 
    if (@$data_ref < 2) {return(1, 1, [1])}
    my %param = (
        nclusters => 2, data => $data_ref, mask => '', weight => '',
        transpose => 0, npass => 10, method => 'a', dist => 'e',
        initialid => []
    );
    my @clusters = Algorithm::Cluster::kcluster(%param);
    my (@cluster1, @cluster2);

    # The jump method is used to determine the optimal number of clusters
    # between 1 and 2.
    for (my $i = 0; $i < @{$clusters[0]}; $i++) {
        if ($clusters[0]->[$i] == 0) {push(@cluster1, $$data_ref[$i]->[0])} 
        if ($clusters[0]->[$i] == 1) {push(@cluster2, $$data_ref[$i]->[0])} 
    }
        
    my @clustered = ([@cluster1, @cluster2], \@cluster1, \@cluster2);
    my @sum_dists;
    
    foreach my $cluster (@clustered) {
        my ($sum, $sum_dist) = (0) x 2;
        $sum += $_ for @$cluster;
        my $mean = $sum / @$cluster;
        $sum_dist += ($_ - $mean)**2 for @$cluster;
        push(@sum_dists, $sum_dist);
    }
    
    my @distortions;
    push(@distortions, $sum_dists[0] / @$data_ref); 
    push(@distortions, ($sum_dists[1] + $sum_dists[2]) / @$data_ref);
    
    # Compute transformed distortions
    my @trans_dist;
    for (@distortions) {push(@trans_dist, $_**-$y)} 
    # Compute max jump
    my $optimal_num_clusters = ($trans_dist[0] >= $trans_dist[1] 
        - $trans_dist[0]) ? 1 : 2;
    my $normal_cluster = 0;
    if ($optimal_num_clusters == 2) {
        my ($a_sum, $a_mean, $b_sum, $b_mean) = (0) x 4;
        map {$a_sum += $_} @cluster1;
        $a_mean = $a_sum / @cluster1;
        map {$b_sum += $_} @cluster2;
        $b_mean = $b_sum / @cluster2;
        
        if ($a_mean > $b_mean) {$normal_cluster = 1}
    }
    return($optimal_num_clusters, $normal_cluster, $clusters[0]);
}

sub collapse_region_list {
    my $regions_ref  = $_[0];    
    my $done = 0;
    return unless ($regions_ref && @$regions_ref && scalar(@$regions_ref) > 1);
    my @revised;
    until ($done) {
        $done = 1 if !@$regions_ref;
        OUTER: for (my $i = 0; $i <= $#$regions_ref; $i++) {
            my $reg1 = $$regions_ref[$i];
            for (my $j = 0; $j <= $#$regions_ref; $j++) {
                if ($i == $j) {next}
                my $reg2 = $$regions_ref[$j];
                if ($$reg1[17] == HD) {next unless ($$reg2[17] == HD)}
                if (overlap(@$reg1[1,2], @$reg2[1,2])->[0]) {
                    for (my $k = 0; $k < @$regions_ref; $k++) {
                        push(@revised, $$regions_ref[$k]) 
                            if ($k != $i && $k != $j);
                    }
                    @$regions_ref = @revised; 
                    @revised = ();
                    my @overlap = @{splice_overlap_regions($reg1, $reg2)};
                    map {$$_[5] = count_SNPs(@$_[START,STOP], \%SNP_BY_POS)} 
                        @overlap;
                    push(@$regions_ref, @overlap);
                    last OUTER;
                }
            }
            $done = 1 if (!@$regions_ref || $i == $#$regions_ref);
        }
        map{for my $x (4..15,18){$$_[$x] ||= 0}}@$regions_ref if @$regions_ref;
    }
}

sub count_SNPs {
    my ($position1, $position2, $hash_ref) = @_;
    my $num_snps  = 0;
    if (!$position1 || !$position2) {return($num_snps)}
    my ($value1, $value2) = ($$hash_ref{$position1}, $$hash_ref{$position2});
    if (!$value1 || !$value2) {
        foreach my $key (keys %$hash_ref) {
            $num_snps++ if ($key >= $position1 && $key <= $position2);
        }
    }
    else {$num_snps = $value2 - $value1 + 1}
    return($num_snps);
}    

sub count_std_mi1_overlap {
    my $regions_ref  = $_[0];
    map{$$_[18] = 0} @$regions_ref;
    for (my $i = 0; $i < @$regions_ref; $i++) {
        my $reg1 = $$regions_ref[$i];
        if ($$reg1[17] == POD) {
            for (my $j = 0; $j < @$regions_ref; $j++) {
                next if ($i == $j);
                my $reg2 = $$regions_ref[$j];
                if ($$reg2[17] == MI1) {
                    if (overlap($$reg1[1], $$reg1[2], $$reg2[1], 
                        $$reg2[2])->[0]) {$$reg1[18]++}
                }
            }
        }
    }
}
        
sub cpu_time {
    # Retrieves and returns the time and date as set on the current machine
    my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
    my @week_days = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
    my ($second, $minute, $hour, $day_of_month, $month, $year_offset, 
        $day_of_week, $day_of_year, $daylight_savings) = localtime();
    my $year = 1900 + $year_offset;
    $day_of_year = $daylight_savings = 0; # to make -w happy
    my $THE_TIME = join("", $hour, ":", $minute,":", $second, " ",
        $week_days[$day_of_week], " ", $months[$month], $day_of_month,
        " ", $year);
    return($THE_TIME);
}

sub cusum {
    # Use cumulative sums of LRR difference from a threshold to detect 
    # critical points indicating start and stop of cnv.
    # For efficiency, the LRR medians and thresholds for deletions are
    # inverted to allow for detection of maximum critical points.
    my ($region_ref, $array_ref, $reg_index, $median, $start_index,$stop_index,
        $critical_index) = @_;
    my @data_refs = @$array_ref; 
    my $data;
    my ($prior_median, $post_median, $begin_cusum, $end_cusum, $min,
        $max_index, $min_index, $found_start_index, $found_end, $cusum, $hd, 
        $del, $amp, $other, $thresh, $local_thresh, $start_crit_index, 
        $end_crit_index, $adj_abn, $father, $mother) = (0) x 21;
    my $max = undef;  

    # For region of 2 snps begin cusum from both ends
    if ($$region_ref[5] == 2 || !$critical_index) {
        ($start_crit_index, $end_crit_index)  = ($stop_index, $start_index);
    }
    else {$start_crit_index = $end_crit_index = $critical_index}
    
    my (@prior, @post);
    my $count = 1;
    my $extend = $$region_ref[5] * 5;
    $extend = 25 if $extend > 25;
    if    ($$region_ref[4] == FATHER) {$father = 1}
    elsif ($$region_ref[4] == MOTHER) {$mother = 1}
        
    if ($REFINING) {
        if ($reg_index == 10)   {$local_thresh = $LOCAL_CH_mBAF_MED}
        elsif ($reg_index == 3) {$local_thresh = -$LOCAL_P1_LRR_MED}        
        elsif ($reg_index == 4) {$local_thresh = -$LOCAL_P2_LRR_MED}
        else                    {$local_thresh = -$LOCAL_CH_LRR_MED}
    }
    else {
        if ($reg_index == 10)   {$local_thresh = $CH_mBAF_MED}
        elsif ($reg_index == 3) {$local_thresh = -$P1_LRR_MED}        
        elsif ($reg_index == 4) {$local_thresh = -$P2_LRR_MED}
        else                    {$local_thresh = -$CH_LRR_MED}
    }
    
    if    ($$region_ref[3] eq "HD")   {$hd  = 1}
    elsif ($$region_ref[3] eq "DEL")  {$del = 1}
    elsif ($$region_ref[3] eq "AMP")  {$amp = 1}
    else {$other = 1}
    $median = -$median if ($hd || $del);
    
    if ($median > -$HD_HI * 2) {$thresh = -$HD_HI}
    else {$thresh = $median - (($median - $local_thresh) / 2)}
    
    # If a non-hd region is large, start cusum $MIN_POD number of 
    # informative SNPs from each end.
    if ($$region_ref[6] > $MIN_POD * 3) {
        if ($$region_ref[17] == POD) {
            my $inf_array_ref;
            my ($start, $stop) 
                = find_closest_indices(@$region_ref[START,STOP,16]);
            if    ($father) {$inf_array_ref = \@P1_INF_SNPS}
            elsif ($mother) {$inf_array_ref = \@P2_INF_SNPS}
            my $hash_ref = ($other) ? \%BAF_BY_POS : \%LRR_BY_POS;
            my $start_pos = $$inf_array_ref[$start + $MIN_POD - 1]->[POS];
            my $end_pos   = $$inf_array_ref[$stop  - $MIN_POD]->[POS];    
            ($start_crit_index, $end_crit_index) 
                = find_closest_indices($start_pos, $end_pos, $hash_ref);
        }
    }

    $count = 0;
    if ($start_index - $extend >= 0) {
        for (my $x = 1; @prior < $extend; $x++) {
            $count++;
            my $data = $data_refs[$start_index - $x]->[$reg_index];
            if ($other) {
                next if ($data < $AA_BOUND || $data > $BB_BOUND);
                $data = 1 - $data if ($data < $LOCAL_CH_BAF_MEAN);
            }
            push(@prior, $data);
        }
        $prior_median = ($hd || $del) ? -median(\@prior) : median(\@prior);
        # Determine if adjacent region is a unique abnormal region
        # with a greater median
        $adj_abn = 1 if ($prior_median >= $median * 2);
        $begin_cusum = $start_index - $count;
        $count  = 0;
        FIND_PRIOR_NORM:
        until ($prior_median < $thresh) {
            for (my $x = $count; $x <= $extend + $count; $x++) {
                if ($start_index - ($extend + $x) < 0) {
                   $$region_ref[START] = $data_refs[0]->[POS];
                   $found_start_index = 1;
                   last FIND_PRIOR_NORM;
                }
                shift(@prior);
                $data = $data_refs[$start_index -($extend + $x)]->[$reg_index];
                if ($other) {
                    next if ($data < $AA_BOUND || $data > $BB_BOUND);
                    $data = 1 - $data if ($data < $LOCAL_CH_BAF_MEAN);
                }
                push(@prior, $data);
            }
            $prior_median = ($hd || $del) ? -median(\@prior) : median(\@prior);
            $begin_cusum = $start_index - ($extend + $count);
            $count += $extend;
        }    
    }
    else {
        for (my $x = 1; $start_index - $x >= 0; $x++) {
            $data = $data_refs[$start_index - $x]->[$reg_index];
            if ($other) {
                next if ($data < $AA_BOUND || $data > $BB_BOUND);
                $data = 1 - $data if ($data < $LOCAL_CH_BAF_MEAN);
            }
            push(@prior, $data);
        }
        $prior_median = median(\@prior);
        $prior_median = -$prior_median if ($prior_median!=NA && ($hd || $del));       
        if ($prior_median > $thresh) {
            $$region_ref[START] = $data_refs[0]->[POS];
            $found_start_index = 1;
        }
        else {$begin_cusum = 0}
    }
    
    unless ($found_start_index) {
        for (my $x = $start_crit_index; $x >= $begin_cusum; $x--) {
            my $value = $data_refs[$x]->[$reg_index];
            # Minimize affect of large outliers in del region
            next if $del && $median < $HD_HI && $value > $HD_HI; 
            if ($other) {
                next if ($value < $AA_BOUND || $value > $BB_BOUND);
                $value = 1 - $value if ($value < $LOCAL_CH_BAF_MEAN);
            }
            # Minimize affect of large outliers in upd region
            next if ($other && $value > $median + 0.1); 
            $value = -$value if ($hd || $del);
            $cusum += $value - $thresh;
            if (!defined($max)) {
                $max = $min = $cusum;
                $max_index = $min_index = $x;
            }
            if ($cusum > $max) {($max, $max_index) = ($cusum, $x)}
            elsif ($cusum < $min) {($min, $min_index) = ($cusum, $x)}
        }
        # If adjacent SNPs are likely to be in another abnormal region keep  
        # initial boundary unless there is a detectable cusum minimum
        if ($adj_abn) {
            if ($min_index != $start_crit_index 
                && $min_index < $begin_cusum + ($extend / 5)) { 
                $$region_ref[START] = $data_refs[$min_index]->[POS];
            }
        }
        elsif ($max_index != $start_crit_index) { 
            $$region_ref[START] = $data_refs[$max_index]->[POS];
        }
        if ($$region_ref[17] == HD && ($$region_ref[4] == FATHER 
            || $$region_ref[4] == MOTHER) && $reg_index != 3 && $reg_index!=4){
            if ($$region_ref[START] < $data_refs[$start_index]->[POS]) {
                $$region_ref[START] = $data_refs[$start_index]->[POS];
            }
        }
    }  
    ($max_index, $min, $min_index, $cusum, $adj_abn) = (0) x 5;
    ($max, $count) = (undef, 0);

    if ($stop_index + $extend <= $#data_refs) { 
        for (my $x = 1; @post < $extend; $x++) {
            $count++;    
            $data = $data_refs[$stop_index + $x]->[$reg_index];
            if ($other) {
                last if ($stop_index + $x >= $#data_refs);
                next if ($data < $AA_BOUND || $data > $BB_BOUND);
                $data = 1 - $data if ($data < $LOCAL_CH_BAF_MEAN);
            }
            push(@post, $data); 
            last if ($stop_index + $x >= $#data_refs);
        }
        $post_median = ($hd || $del) ? -median(\@post) : median(\@post);
        # Determine if adjacent region is a unique abnormal region
        # with a greater median
        $adj_abn = 1 if ($post_median >= $median * 2);
        $end_cusum = $stop_index + $count;
        FIND_POST_NORM:
        until ($post_median < $thresh) {
            for (my $x = $count; $x <= $extend + $count; $x++) {
                if ($stop_index + $extend + $x > $#data_refs) {
                   $$region_ref[STOP] = $data_refs[-1]->[POS];
                   $found_end = 1;
                   last FIND_POST_NORM;
                }
                shift(@post);
                $data = $data_refs[$stop_index + ($extend + $x)]->[$reg_index];
                if ($other) {
                    next if ($data < $AA_BOUND || $data > $BB_BOUND);
                    $data = 1 - $data if ($data < $LOCAL_CH_BAF_MEAN);
                }
                push(@post, $data);
            }
            $post_median = ($hd || $del) ? -median(\@post) : median(\@post);
            $end_cusum = $stop_index + ($extend + $count);                            
            $count += $extend;
        }
    }
    else {
        for (my $x = 1; $x <= $#data_refs - $stop_index; $x++) {
            $data = $data_refs[$stop_index + $x]->[$reg_index];
            if ($other) {
                next if ($data < $AA_BOUND || $data > $BB_BOUND);
                $data = 1 - $data if ($data < $LOCAL_CH_BAF_MEAN);
            }
            push(@post, $data);                    
        }
        $post_median = median(\@post);
        $post_median = -$post_median if ($prior_median != NA && ($hd || $del));               
        if ($post_median > $thresh) {
            $$region_ref[STOP] = $data_refs[-1]->[POS];
            $found_end = 1;
        }
        else {$end_cusum = $#data_refs}
    }

    unless ($found_end) {
        for (my $x = $end_crit_index; $x <= $end_cusum; $x++) {
            my $value = $data_refs[$x]->[$reg_index];
            # Minimize affect of large outliers in del region
            next if $del && $median < $HD_HI && $value > $HD_HI; 
            if ($other) {
                next if ($value < $AA_BOUND || $value > $BB_BOUND);
                $value = 1 - $value if ($value < $LOCAL_CH_BAF_MEAN);
            }
            # Minimize affect of large outliers in upd region
            next if ($other && $value > $median + 0.1);  
            $value = -$value if ($hd || $del);
            $cusum += $value - $thresh;
            if (!defined($max)) {
                $max = $min = $cusum;
                $max_index = $min_index = $x;
            }
            if ($cusum > $max) {($max, $max_index) = ($cusum, $x)}
            elsif ($cusum < $min) {($min, $min_index) = ($cusum, $x)}
        }
        # If adjacent SNPs are likely to be in another abnormal region keep  
        # initial boundary unless there is a detectable cusum minimum
        if ($adj_abn) {
            if ($min_index != $end_crit_index 
                && $min_index < $end_cusum - ($extend / 5)) { 
                $$region_ref[STOP] = $data_refs[$min_index]->[POS];
            }
            elsif ($max_index != $end_crit_index) {
                $$region_ref[STOP] = $data_refs[$max_index]->[POS];
            }
        }
        elsif ($max_index != $end_crit_index) { 
            $$region_ref[STOP] = $data_refs[$max_index]->[POS];
        }
        if ($$region_ref[17] == HD && ($$region_ref[4] == FATHER 
            || $$region_ref[4] == MOTHER) && $reg_index != 3 && $reg_index!=4){
            if ($$region_ref[STOP] > $data_refs[$stop_index]->[POS]) {
                $$region_ref[STOP] = $data_refs[$stop_index]->[POS];
            }
        }
    }
    return($region_ref);
}
    
sub define_mi1_regions {
    my $mi1_ref = $_[0];
    map {@{$_}[16,17] = (\%MI1_BY_POS, MI1)} @$mi1_ref;
    map {for my $x (7..15,18) {$$_[$x] ||= 0}} @$mi1_ref;
    map {$_ = evaluate_region($_)} @$mi1_ref;
    collapse_region_list(\@$mi1_ref);
    map {$_ = evaluate_region($_) unless $$_[6]} @$mi1_ref;
    my @mi1_regions;
    map {push(@mi1_regions, $_) if ($$_[6] >= $MIN_MI1)} @$mi1_ref;
    map {calculate_region_stats($_)} @mi1_regions;
    return (\@mi1_regions);
}

sub define_POD_region {
    # Combines and calculates descriptive info for a detected POD region.
    my ($ref, $region1, $last_sig_SNP) = @_;
    my (@abn_region, $array_ref);
    for (my $i = 0; $i < @$ref; $i++) {
        if ($i == 0) {@abn_region[0,START,STOP] = @{$$ref[$i]}[0,START,STOP]}
        elsif ($i == $last_sig_SNP) {$abn_region[STOP] = $$ref[$i]->[STOP]}
        if ($i == @$ref - 1) {
            if ($region1) {
                @abn_region[3,4,16,17] = (" ",FATHER,\%P1_INF_SNP_BY_POS, POD);
                $array_ref = \@P1_INF_SNPS;
            } 
            else {
                @abn_region[3,4,16,17] = (" ",MOTHER,\%P2_INF_SNP_BY_POS, POD);
                $array_ref = \@P2_INF_SNPS;
            }
        }
    }   
    for (4..17) {$abn_region[$_] ||= 0}
    my ($closest_start_ind, $closest_stop_ind) 
        = find_closest_indices(@abn_region[START,STOP,16]);
    $abn_region[START] = $$array_ref[$closest_start_ind]->[POS];
    $abn_region[STOP]  = $$array_ref[$closest_stop_ind]->[POS];
    $abn_region[6] = count_SNPs(@abn_region[START,STOP,16]);
    @abn_region = @{evaluate_region(\@abn_region)};
    $abn_region[6] = count_SNPs(@abn_region[START,STOP,16]);
    return(\@abn_region);   
}                

sub define_PODcr_regions {
    my $outlier_regs = $_[0];
    my $array_ref;
    map {@{$_}[16,17] = (\%BAF_OUTLIERS_BY_POS, PODcr)}@$outlier_regs;
    map {for my $x (7..15,18) {$$_[$x] ||= 0}} @$outlier_regs;
    collapse_region_list(\@$outlier_regs);
    for my $ref (@$outlier_regs) {
        my ($closest_start_ind, $closest_stop_ind) 
            = find_closest_indices(@$ref[START,STOP,16]);
        $$ref[START] = $BAF_OUTLIERS[$closest_start_ind]->[POS];
        $$ref[STOP]  = $BAF_OUTLIERS[$closest_stop_ind]->[POS];
        $$ref[6] = count_SNPs(@$ref[START,STOP,16]);
    }
    map {$_ = evaluate_region($_)} @$outlier_regs;
    collapse_region_list(\@$outlier_regs);
    map {$_ = evaluate_region($_) unless $$_[6]} @$outlier_regs;
    map {calculate_region_stats($_)} @$outlier_regs;
    
    # To avoid regions due to BAF genomic waves, determine if the distribution
    # of AB calls on is as expected above the upper thresh and below the lower
    my @results;
    foreach my $ref (@$outlier_regs) {
        my ($upper, $lower) = (0) x 2;
        my ($reg_start, $reg_stop, $success) 
            = find_closest_indices(@$ref[START,STOP], \%BAF_OUTLIERS_BY_POS);
        if ($success) {
            for (my $j = $reg_start; $j <= $reg_stop; $j++) {
                if ($BAF_OUTLIERS[$j]->[10] > $AB_UPPER_BOUND) {$upper++}
                elsif ($BAF_OUTLIERS[$j]->[10] < $AB_LOWER_BOUND) {$lower++}
            }
            my $pvalue = bin_test($lower, $lower + $upper, 0.5, "t");
            my $sidak_alpha = 1 - (1 - $ALPHA)**(1 / scalar(@$outlier_regs));
            if ($pvalue > $sidak_alpha) {push(@results, $ref)}
        }
    }
    return (\@results);
}

sub detect_hom_del { 
    my (@p1_hd_indices, @p2_hd_indices, @child_hd_indices); 
    my ($prev_snp, $next_snp, $size);
    # Locate SNPs indicating a homozygous deletion
    # Will be labeled as parent lacking information content and 
    # switched to parent contributing before printing
    for (my $i = 0; $i < @LRR_REFS; $i++) {
        my $line_ref  = $LRR_REFS[$i];
        my ($p1_LRR, $p2_LRR, $child_LRR) = @$line_ref[3..5];
        if ($child_LRR <= $HD_HI) {
            $CH_HD_BY_POS{$$line_ref[POS]} = keys %CH_HD_BY_POS;
            push(@child_hd_indices, $i);
        }
        if ($p1_LRR <= $HD_HI) {
            $P2_HD_BY_POS{$$line_ref[POS]} = keys %P2_HD_BY_POS;
            push(@p2_hd_indices, $i);
        }
        if ($p2_LRR <= $HD_HI) {
            $P1_HD_BY_POS{$$line_ref[POS]} = keys %P1_HD_BY_POS;
            push(@p1_hd_indices, $i);
        }
    }
     
    # Extend each SNP into regions
    my @p1_hd = @{extend_hd(\@p1_hd_indices, 3)};
    my @p2_hd = @{extend_hd(\@p2_hd_indices, 4)};
    my @ch_hd = @{extend_hd(\@child_hd_indices, 5)};
    
    if (@p1_hd) {
        map {for my $x (4..18) {$$_[$x] ||= 0}} @p1_hd;
        map {$_ = evaluate_region($_, FATHER)} @p1_hd;
        map {$_ = evaluate_region($_)} @p1_hd;
        collapse_region_list(\@p1_hd);
        map {@{$_}[16,17] = (\%P1_HD_BY_POS, HD)} @p1_hd;
        map {$_ = evaluate_region($_) unless $$_[6]} @p1_hd;
    }
    if (@p2_hd) {
        map {for my $x (4..18) {$$_[$x] ||= 0}} @p2_hd;
        map {$_ = evaluate_region($_, MOTHER)} @p2_hd;
        map {$_ = evaluate_region($_)} @p2_hd;    
        collapse_region_list(\@p2_hd);
        map {@{$_}[16,17] = (\%P2_HD_BY_POS, HD)} @p2_hd;
        map {$_ = evaluate_region($_) unless $$_[6]} @p2_hd;
    }
        
    if (@ch_hd) {
        map {for my $x (4..18) {$$_[$x] ||= 0}} @ch_hd;
        map {$_ = evaluate_region($_)} @ch_hd;
        collapse_region_list(\@ch_hd);
        map {@{$_}[16,17] = (\%CH_HD_BY_POS, HD)} @ch_hd;
        map {$_ = evaluate_region($_) unless $$_[6]} @ch_hd;
    }
    # Correct any overlap between trio members
    my (@combine_regions, @return_regions);
    push(@combine_regions, @p1_hd, @p2_hd, @ch_hd);
    if (@combine_regions) { 
        collapse_region_list(\@combine_regions);
        map {$_ = evaluate_region($_) unless $$_[6]} @combine_regions;
        
        # Count SNPs and informative SNPS
        my $hd_thresh;
        foreach my $ref (@combine_regions) {
            $$ref[5] = count_SNPs(@$ref[1,2], \%SNP_BY_POS);
            if ($$ref[4] == FATHER) {
                $$ref[16] = \%P1_HD_BY_POS;
                $$ref[6] = count_SNPs(@$ref[1,2,16]);
                $hd_thresh = $MIN_P2_HD;
            }
            elsif ($$ref[4] == MOTHER) {
                $$ref[16] = \%P2_HD_BY_POS;
                $$ref[6] = count_SNPs(@$ref[1,2,16]);
                $hd_thresh = $MIN_P1_HD;
            }
            elsif ($$ref[4] == NONE) {
                $$ref[16] = \%CH_HD_BY_POS;
                $$ref[6] = count_SNPs(@$ref[1,2,16]);
                $hd_thresh = $MIN_CH_HD;
            }
            elsif ($$ref[4] == BOTH) {
                if ($MIN_P1_HD <= $MIN_P2_HD) {
                    $$ref[16] = \%P2_HD_BY_POS;
                    $hd_thresh = $MIN_P1_HD;
                }
                else {
                    $$ref[16] = \%P1_HD_BY_POS;
                    $hd_thresh = $MIN_P2_HD;
                }
            }
            else {$$ref[6] = NA}
            if ($$ref[6] >= $hd_thresh) {push(@return_regions, $ref)}
        }        
    }
    map {@{$_}[17] = HD} @return_regions if @return_regions;
    map {calculate_region_stats($_)} @return_regions if @return_regions;
    return(\@return_regions);
}            
  
sub determine_R_status {
    # Determines if the R script is still running and makes one attempt to 
    # restart the Rscript.
    my $R_status = `ps -p $R_PID -o pid=`;
    if (!$R_status && !$R_ATTEMPTS) { 
        if (!$R_ATTEMPTS) {
            $R_ATTEMPTS++;
            start_R();
            sleep(5);
            ($R_PID_FILE, $R_PID_FILENAME) = tempfile(UNLINK => 1);
            if(!get_R_PID()) {$R_ATTEMPTS++}
        }            
    }
    return($R_status);
} 

sub evaluate_region {
    my ($region_ref, $contrib) = @_;
    my $count = 1;
    my ($index, $local_lrr);
    $contrib ||= CHILD;
    my @region_array = @$region_ref;
    my @array = ();
    return() if (!@$region_ref);
    if ($contrib == FATHER) {
        $local_lrr = ($REFINING) ? $LOCAL_P2_LRR_MED : $P2_LRR_MED;
        $index = 4;
    }
    elsif ($contrib == MOTHER) {
        $local_lrr = ($REFINING) ? $LOCAL_P1_LRR_MED : $P1_LRR_MED;
        $index = 3;
    }
    else {
        $local_lrr = ($REFINING) ? $LOCAL_CH_LRR_MED : $CH_LRR_MED;
        $index = 5;
    }
    $region_array[5] = count_SNPs(@region_array[START,STOP], \%SNP_BY_POS);
    my ($median, $start_index, $stop_index, $lowest_index, 
        $highest_index) = @{calculate_region_median($region_ref, LRR, $index)};
    if ($median < $DEL_UPPER + $local_lrr) {
        $region_array[3] = "DEL" if ($region_array[3] ne "HD");
        @region_array = @{cusum(\@region_array, \@LRR_REFS, $index, $median,
            $start_index, $stop_index, $lowest_index)};
    }
    elsif ($median >= $AMP + $local_lrr) {
        $region_array[3] = "AMP";
        @region_array = @{cusum(\@region_array, \@LRR_REFS, $index, $median,
            $start_index, $stop_index, $highest_index)};
    }
    else {
        $index = 10;
        $region_array[3] = " ";
        my $median_index;
        ($median, $start_index, $stop_index, $lowest_index, $highest_index) 
            = @{calculate_region_median($region_ref, BAF, $index)};
        if ($median == NA) {return($region_ref)}
        elsif ($median > $LOCAL_CH_mBAF_MED && $median < $BB_BOUND) {
            my (@inf, $median_index);
            my $hash_ref = $region_array[16];
            foreach my $key (keys %$hash_ref) {
                push(@inf, $key) if ($key >= $region_array[1] 
                    && $key <= $region_array[2]);
            }
            my @sorted_inf = sort { $a <=> $b } @inf;
            my $median_inf = $sorted_inf[int(@sorted_inf / 2)];
            if ($median_inf) {
                for (my $i = $start_index; $i <= $stop_index; $i++) {
                    $median_index = $i if ($CURR_CHR_REFS[$i]->[POS] 
                        == $median_inf)
                }
                @region_array = @{cusum(\@region_array, \@CURR_CHR_REFS,$index,
                    $median, $start_index, $stop_index, $median_index)};
            }
        }
        else {$region_array[3] = " "}    
    }
    if ($contrib != CHILD) {$region_array[3] = " "}
    if ($region_array[START] >= $region_array[STOP]) {
        $region_array[START] = $$region_ref[START];
        $region_array[STOP]  = $$region_ref[STOP];
    }
    $region_array[5] = count_SNPs(@region_array[START,STOP], \%SNP_BY_POS);
    $region_array[6] = count_SNPs(@region_array[START,STOP,16]);
    return(\@region_array);  
}   
  
sub extend_hd {
    my ($indices_ref, $lrr_index) = @_;
    my $tree = Tree::Interval->new();
    # Extension array - first index = 0 or 1 has been extended, 
    # second index = fallback position, third index = count
    my @regions;
    foreach my $index_ref (@$indices_ref) {
        my ($up_index, $low_index) =  (1) x 2;
        my @reg_array = ();
        my $line_ref = $LRR_REFS[$index_ref];
        @reg_array[0..3] = (@$line_ref[CHR,POS], $$line_ref[POS], " ");
        if (@regions) {next if ($tree->find($$line_ref[POS]))}

        my ($upper_end, $lower_end) = (0) x 2;
        until ($upper_end) {
            $upper_end = 1;
            if ($LRR_REFS[$index_ref + $up_index]) {
                my $curr_line = $LRR_REFS[$index_ref + $up_index];
                if (!$tree->find($$curr_line[POS])) {
                    if ($$curr_line[$lrr_index] < $HD_HI) {
                        $reg_array[STOP] = $$curr_line[POS];
                        $upper_end = 0;
                    }
                }
            }
            $up_index++;
        }
        until ($lower_end) {
            $lower_end = 1;
            if ($LRR_REFS[$index_ref - $low_index]) {
                my $curr_line = $LRR_REFS[$index_ref - $low_index];
                if (!$tree->find($$curr_line[POS])) {                        
                    if ($$curr_line[$lrr_index] < $HD_HI) {
                        $reg_array[START] = $$curr_line[POS];
                        $lower_end = 0;
                    }
                }
            }
            $low_index++;
        }
        
        my $hd_thresh;
        if ($lrr_index == 3)    {
            @reg_array[4,16] = (FATHER, \%P1_HD_BY_POS);
            $hd_thresh = $MIN_P2_HD;
        }
        elsif ($lrr_index == 4) {            
            @reg_array[4,16] = (MOTHER, \%P2_HD_BY_POS);
            $hd_thresh = $MIN_P1_HD;
        }
        elsif ($lrr_index == 5)   { 
            @reg_array[3,4,16] = ("HD", NONE, \%CH_HD_BY_POS);        
            $hd_thresh = $MIN_CH_HD;
        }
        @reg_array[6,17] = (count_SNPs(@reg_array[START,STOP,16]), HD);

        if ($reg_array[6] >= $hd_thresh) {
            for (4..15,18) {$reg_array[$_] ||= 0}
            push(@regions, \@reg_array);
            $tree->insert(@reg_array[START,STOP,START]);
        }
    }
    return(\@regions);
}                            

sub find_closest_indices {
    my ($start, $stop, $hash_ref) = @_;
    return unless ($start && $stop && $hash_ref);
    my ($closest_lower, $closest_higher, $present) = (0) x 3;
    my @indices;
    my ($value1, $value2) = ($$hash_ref{$start}, $$hash_ref{$stop});
    if (!$value1 || !$value2) {
        foreach my $key (keys %$hash_ref) {
            push(@indices, $$hash_ref{$key}) if($key>=$start && $key<=$stop);
        }
        if (@indices) {
            my @sorted_indices = sort { $a <=> $b } @indices;
            ($closest_lower, $closest_higher) = @sorted_indices[0,-1];
            $present = 1;
        }
    }
    else {
        ($closest_lower, $closest_higher) = ($value1, $value2);
        $present = 1;
    }
    return($closest_lower, $closest_higher, $present);
}

sub find_informative_snps {
    # Analyzes each SNP for parental contribution. 
    # Results are added to the input array in binary form.
    # Returns a number of variables for genomewide calculations.
    my $array_ref = $_[0];
    my ($p1, $p2, $p1_mi1, $p2_mi1, $un_mi1, $outlier, $ch_upper_MI1, 
        $ch_lower_MI1, $outlier_ct, $iso_or_het, $inf_hUPD);
        
    if ($$array_ref[10] > $MI1_UPPER_THRESH) {$ch_upper_MI1 = 1} 
    elsif ($$array_ref[10] < $MI1_LOWER_THRESH) {$ch_lower_MI1 = 1} 
    
    if ($$array_ref[3] eq "BB") { 
        if ($$array_ref[6] eq "BB" && $$array_ref[9] eq "AB") {$un_mi1 = 1}
        elsif ($$array_ref[6] eq "AA") {
            if ($$array_ref[10] > $AB_UPPER_BOUND) {
                $p1 = $iso_or_het = 1;
                $p1_mi1 = 1 if ($ch_upper_MI1);
            }
            elsif ($$array_ref[10] < $AB_LOWER_BOUND) {
                $p2 = $iso_or_het = 1;
                $p2_mi1 = 1 if ($ch_lower_MI1);
            }
            $inf_hUPD = 1 if ($$array_ref[9] eq "BB" || $$array_ref[9] eq "AA");
        }
        elsif ($$array_ref[6] eq "AB") {
            if ($$array_ref[10] < $AB_LOWER_BOUND) {
                $p2 = 1;
                $p2_mi1 = 1 if ($ch_lower_MI1);
            }
            $inf_hUPD = 1 if ($$array_ref[9] eq "BB" || $$array_ref[9] eq "AB");
        }
        elsif (substr($$array_ref[6], 0, 1) eq "N") {
            if (($$array_ref[10] < $AB_LOWER_BOUND || $$array_ref[9] eq "AA") 
                && $$array_ref[8] > $HD_LOW) {$p2 = 1}
        }
    }
    elsif ($$array_ref[3] eq "AA") {    
        if ($$array_ref[6] eq "BB") {
            if ($$array_ref[10] > $AB_UPPER_BOUND) {
                $p2 = $iso_or_het = 1;
                $p2_mi1 = 1 if ($ch_upper_MI1);
            }
            elsif ($$array_ref[10] < $AB_LOWER_BOUND) {
                $p1 = $iso_or_het = 1;
                $p1_mi1 = 1 if ($ch_lower_MI1);
            }
            $inf_hUPD = 1 if ($$array_ref[9] eq "BB" || $$array_ref[9] eq "AA");
        }
        elsif ($$array_ref[6] eq "AB") {
            if ($$array_ref[10] > $AB_UPPER_BOUND) {
                $p2 = 1;
                $p2_mi1 = 1 if ($ch_upper_MI1);
            }
            elsif ($$array_ref[10] < $AB_LOWER_BOUND 
                && $$array_ref[9] eq "AB") {$outlier = 1}
            $inf_hUPD = 1 if ($$array_ref[9] eq "AA" || $$array_ref[9] eq "AB");
        }
        elsif (substr($$array_ref[6], 0, 1) eq "N") {
            if (($$array_ref[10] > $AB_UPPER_BOUND || $$array_ref[9] eq "BB") 
                && $$array_ref[8] > $HD_LOW) {$p2 = 1}
        }
        elsif ($$array_ref[6] eq "AA" && $$array_ref[9] eq "AB") {$un_mi1 = 1}
    }
    elsif ($$array_ref[3] eq "AB") {
        if ($$array_ref[6] eq "BB") {
            if ($$array_ref[10] < $AB_LOWER_BOUND) {
                $p1 = 1;
                $p1_mi1 = 1 if ($ch_lower_MI1);
            }
            $inf_hUPD = 1 if ($$array_ref[9] eq "BB" || $$array_ref[9] eq "AB");
        }
        elsif ($$array_ref[6] eq "AA") {
            if ($$array_ref[10] > $AB_UPPER_BOUND) {
                $p1 = 1;
                $p1_mi1 = 1 if ($ch_upper_MI1);
            }
            $inf_hUPD = 1 if ($$array_ref[9] eq "AA" || $$array_ref[9] eq "AB");
        }
    }
    elsif (substr($$array_ref[3], 0, 1) eq "N") {
        if ($$array_ref[6] eq "AA") {
            if (($$array_ref[10] > $AB_UPPER_BOUND || $$array_ref[9] eq "BB") 
                && $$array_ref[5] > $HD_LOW) {$p1 = 1}
        }                            
        elsif ($$array_ref[6] eq "BB") {
            if (($$array_ref[10] < $AB_LOWER_BOUND || $$array_ref[9] eq "AA") 
                && $$array_ref[5] > $HD_LOW) {$p1 = 1}
        }                    
    }
            
    if ($$array_ref[9] eq "AB" && ($$array_ref[10] > $PODCR_UPPER_THRESH
        || $$array_ref[10] < $PODCR_LOWER_THRESH)) {$outlier_ct = 1}
            
    my ($inf_snp, $mi1) = (0) x 2;
    if ($$array_ref[5] < $HD_LOW || $$array_ref[8] < $HD_LOW 
        || $$array_ref[11] < $HD_LOW) {push(@$array_ref, UNKNOWN)}   
    elsif ($p1) {
        push(@$array_ref, FATHER);
        $inf_snp = 1 unless ($$array_ref[CHR] == CHR_X);
    }
    elsif ($p2) {
        push(@$array_ref, MOTHER);
        $inf_snp = 1 unless ($$array_ref[CHR] == CHR_X);
    }
    elsif ($un_mi1) {push(@$array_ref, UNKNOWN_MI1)}
    else {push(@$array_ref, UNKNOWN)}
    my $type = 0;
    if ($p1_mi1) {($mi1, $type) = ($p1_mi1, FATHER)}
    elsif ($p2_mi1) {($mi1, $type) = ($p2_mi1, MOTHER)}
    if ($un_mi1) {($mi1, $type) = ($un_mi1, NONE)}
    return([$p1, $p2, $inf_snp, $mi1, $type, $outlier_ct, $iso_or_het, 
    $inf_hUPD]);
}

sub find_POD_regions {
    # Calculates pvalues for windows with adequate 
    # informative SNPs and combines overlapping abnormal windows.
    # Calls define_POD_region and returns a reference to the refined regions.
    my ($array_ref, $window_size, $acceptable_errors) = @_; 
    my (@abnormal_regions, @abn_region, @temp_region_refs, @return_regions);
    my ($flag1, $flag2, $in_region1, $in_region2, $last_sig_SNP, $p1_sig, 
        $p2_sig, $chr, $end_of_region, $abn_region, $region_size) = (0) x 11;
    my $alpha = $POD_ALPHA; 
    my $min_size = $MIN_POD;
    return() if (scalar(@$array_ref) < $window_size);
 
    my $windows = overlapping_windows($array_ref, $window_size);
    if ($window_size > $SIZE_SMALL_WINDOW) {
        ($alpha, $min_size) = ($LARGE_POD_ALPHA, $MIN_LARGE_POD);
    }

    for (my $i = 0; $i < @$windows; $i++) {
        my @curr_window = @{$$windows[$i]};
        my ($p1_ct, $p2_ct) = @curr_window[3,4];
        my $not_in_region = (!$in_region1 && !$in_region2);
        $chr = $curr_window[0];  
        # If the last window overlaps a previous abnormal region, it may 
        # contain another small nonoverlapping region. In order to detect
        # the additional region, the window boundaries are adjusted to be
        # nonoverlapping and the informative SNPs are recalculated.
        if ($curr_window[START] <= $end_of_region && $not_in_region) { 
            if ($i == @$windows - 1) { 
                $curr_window[START] 
                    = $SNP_BY_NUM{$SNP_BY_POS{$end_of_region} + 1};
                $p1_ct = count_SNPs(@curr_window[START,STOP],
                    \%P1_INF_SNP_BY_POS);
                $p2_ct = count_SNPs(@curr_window[START,STOP], 
                    \%P2_INF_SNP_BY_POS);
            }
            else {$p1_sig = $p2_sig = 0; next}
        }
        # Ignore windows with number of informative SNPs less than
        # the previously calculated minimum number of SNPs for 
        # which the probability can possibly be below the significance
        # threshold. 
        my $abnormal = my $next = 0;
        if ($p1_ct + $p2_ct < $min_size) { 
            if ($not_in_region) {$p1_sig = $p2_sig = 0; next}
        }
        else {
            my $pvalue = bin_test($p1_ct, $p1_ct + $p2_ct, 0.5, "t");

            # If the probability is significantly abnormal, determine
            # which parent was the major contributor 
            if ($pvalue <= $alpha) {
                $abnormal = 1;
                ($p1_sig, $p2_sig) = ($p1_ct > $p2_ct) ? (1,0) : (0,1);
            }
            elsif ($not_in_region) {$p1_sig = $p2_sig = 0; next}
            else {$p1_sig = $p2_sig = 0}
        }
         
        if ($in_region1) {
            if (!$p1_sig && $p2_ct > $ACCEPTABLE_ERRORS) {$flag2 = 1}
            elsif (!$p1_sig && $curr_window[5]) {
                my $pvalue = bin_test($p2_ct, $window_size, $ERROR_RATE, "g");
            }
        }
        elsif ($in_region2) {
            if (!$p2_sig && $p1_ct > $ACCEPTABLE_ERRORS) {$flag1 = 1}
            elsif (!$p1_sig && $curr_window[5]) {
                my $pvalue = bin_test($p1_ct, $window_size, $ERROR_RATE, "g");
            }
        }
            
        # Define abnormal regions
        if ($not_in_region) {
            if ($p1_sig) {
                push(@temp_region_refs, \@curr_window);
                ($last_sig_SNP, $in_region1) = (0,1);
                next;
            }
            elsif ($p2_sig) {
                push(@temp_region_refs, \@curr_window);
                ($last_sig_SNP, $in_region2) = (0,1);
                next;
            }
        }
        elsif ($in_region1) {
            # Store current line or finalize calculations for region
            if (!$flag2) {
                push(@temp_region_refs, \@curr_window);
                $last_sig_SNP = $#temp_region_refs if ($p1_sig);           
            }
            else {
                $abn_region = define_POD_region(\@temp_region_refs,
                    $in_region1, $last_sig_SNP); 
                push(@abnormal_regions, $abn_region);
                $end_of_region = $$abn_region[2];
                @temp_region_refs = ();
                $in_region1 = $flag2 = $last_sig_SNP = 0;            
                if ($p2_sig) { 
                    push(@temp_region_refs, \@curr_window);
                    ($last_sig_SNP, $in_region2) = (0,1);
                }
            }      
        }
        elsif ($in_region2) {
            # Store current line or finalize calculations for region
            if (!$flag1) {
                push(@temp_region_refs, \@curr_window);            
                $last_sig_SNP = $#temp_region_refs if ($p2_sig);
            }
            else { 
                unless ($chr == CHR_X && $GENDER == MALE) {
                    $abn_region = define_POD_region(\@temp_region_refs, 
                        $in_region1, $last_sig_SNP); 
                    push(@abnormal_regions, $abn_region);
                    $end_of_region = $$abn_region[2];                                    
                }
                $end_of_region = $$abn_region[2];                
                @temp_region_refs = ();
                $in_region2 = $flag1 = $last_sig_SNP = 0;
                
                if ($p1_sig) { 
                    push(@temp_region_refs, \@curr_window);
                    ($last_sig_SNP, $in_region1) = (0,1);
                }
            }
        }   
        $p1_sig = $p2_sig = 0;    
    }

    if ($temp_region_refs[-1]) { 
        unless ($in_region2 && $chr == CHR_X && $GENDER == MALE) {
            $abn_region = define_POD_region(\@temp_region_refs, $in_region1,
                $last_sig_SNP); 
            push(@abnormal_regions, $abn_region);
        }
    }
    if (@abnormal_regions) {
        collapse_region_list(\@abnormal_regions);
        map {$_ =  evaluate_region($_)} @abnormal_regions;
        map {calculate_region_stats($_)} @abnormal_regions;
    }
    return(\@abnormal_regions);    
}

sub get_R_PID {
    $R_PID = <$R_PID_FILE>;
    if (!$R_PID) {
        print STDERR "An error has occurred in communication with the R ", 
        "software.\n"
    }
    return($R_PID);
}

sub help {
    print <<END;

Usage: perl $0 [--options] [INPUT_FILE]
  --alpha   A threshold for which a p-value <= alpha is considered significant.
            Default = --alpha=0.1
  --batch   Submit a file containing a list of file names to be analyzed
  --build   The path to a file containing the UCSC Genome assembly 
            for centromere locations.
            Default = --build=./genome_build/hg18_centromeres.txt
  --cite    Prints reference info for citations
  --cores   Number of CPU cores to employ
            Default = maximum cores - 1 (e.g. --cores=8)
  --gender  Gender designation for sample (M or F) (--gender=F)
            Currently, chromosome X is analyzed only if a female gender is 
            specified. Default = NA
  --graph   Creates graphic output of results in PNG or PDF format
            (--graph=none, --graph=png, --graph=pdf, or --graph=both)
            Default = --graph=none
  --hd      Abnormality detection by homozygous deletion analysis (PODhd)
            (--hd or --nohd) Default = --hd 
  --help    Prints a help message and exits
  --hetSD   Heterozygous SNP threshold. StDevs from mean BAF
            Default = --hetSD=1.414213562373095 (sqrt of 2)
  --homSD   Homozygous SNP threshold. StDevs from mean BAF
            Default = --homSD=4
  --mi1     Abnormality detection by Mendelian "error" analysis (PODmi1)
            (--mi1 or --nomi1) Default = --mi1 
  --nc      Samples must have a No Call rate below this threshold.
            Default = --nc=0.03             
  --out     Specify an output directory (e.g. --out=results)
            Default = --out=triPOD_Results  
  --pod     Abnormality detection by standard POD algorithm 
            (--pod or --nopod) Default = --pod    
  --podcr   Abnormality detection by cryptic POD algorithm (PODcr) 
            (--podcr or --nopodcr) Default = --podcr     
  --stats   Create a file including calculated parameters, etc.
            (--stats) Not default.
  --verbose Prints progress info to screen. Negated using --noverbose.
            --batch mode will run in --noverbose mode.
            Default = --verbose
  --win     Number of SNPs per window for sliding window analysis
            Default = --win=100

END
    print_proper_format();
    exit 0;
}

sub initial_calculations {
    my (@ch_stats, @p1_stats, @p2_stats);
    @p1_stats = @p2_stats = map {$ch_stats[$_] = 0} (0..10);
    # The following are *_stats array indices
    my ($CHR, $START, $BAF_REF, $MBAF_MED, $LRR_MED, $NC_CT,
        $BAF_SNP_CT, $LRR_SNP_CT, $NC_SNP_CT, $STDEV) = (0..9);
    my ($ch_baf_clust_sum, $ch_baf_clust_sum_sq, $ch_baf_clust_ct, 
        $p1_baf_clust_sum, $p1_baf_clust_sum_sq, $p1_baf_clust_ct,
        $p2_baf_clust_sum, $p2_baf_clust_sum_sq, $p2_baf_clust_ct);        
    my (@ch_baf, @ch_mbaf, @ch_lrr, @mi1_count, @p1_baf, @p1_mbaf, @p1_lrr,
        @p2_baf, @p2_mbaf, @p2_lrr);
    for (0..10) {$ch_baf[$_] = 0}
    @p1_baf = @p2_baf = @ch_baf;
    my $snp_ct = 0;
    foreach my $line_ref (@CURR_CHR_REFS) {
        $snp_ct++;
        my @line_array = @$line_ref;
        next if (!$line_array[3] || !$line_array[6] || !$line_array[9]);
        $p1_stats[$NC_CT]++ if substr($line_array[3], 0, 1) eq "N";
        $p2_stats[$NC_CT]++ if substr($line_array[6], 0, 1) eq "N";
        $ch_stats[$NC_CT]++ if substr($line_array[9], 0, 1) eq "N";
        $p1_stats[$NC_SNP_CT]++;
        next if (substr($line_array[4], 0, 1) eq "N" 
            ||  substr($line_array[7],  0, 1) eq "N" 
            ||  substr($line_array[10], 0, 1) eq "N");
        push(@p1_lrr, $line_array[5]);
        push(@p2_lrr, $line_array[8]);
        push(@ch_lrr, $line_array[11]);
        $p1_stats[$LRR_SNP_CT]++;
        
        if (substr($line_array[0], 0, 2) =~ /cn/i) {
            $p1_stats[$NC_CT]-- if substr($line_array[3], 0, 1) eq "N";
            $p2_stats[$NC_CT]-- if substr($line_array[6], 0, 1) eq "N";
            $ch_stats[$NC_CT]-- if substr($line_array[9], 0, 1) eq "N";
            next;
        }
        $p1_stats[$BAF_SNP_CT]++;
        
        # Avoid 1% of SNPs near centromere and telomere 
        # for global BAF calculations
        if ($snp_ct > $#CURR_CHR_REFS * 0.01 
            && $snp_ct < $#CURR_CHR_REFS * 0.99) {
            # Store values for BAF StDev calculations
            baf_stats_for_stdev($line_array[4], 0.25,0.75,\@p1_baf,\@p1_mbaf); 
            baf_stats_for_stdev($line_array[7], 0.25,0.75,\@p2_baf,\@p2_mbaf); 
            baf_stats_for_stdev($line_array[10],0.25,0.75,\@ch_baf,\@ch_mbaf); 
                     
            # Store data for stdev clustering
            if ($line_array[10] > 0.02 && $line_array[10] < 0.98) {
                $ch_baf_clust_sum += $line_array[10];
                $ch_baf_clust_sum_sq += $line_array[10]**2;
                $ch_baf_clust_ct++;
            }
            if ($line_array[4] > 0.02 && $line_array[4] < 0.98) {
                $p1_baf_clust_sum += $line_array[4];
                $p1_baf_clust_sum_sq += $line_array[4]**2;
                $p1_baf_clust_ct++;
            }
            if ($line_array[7] > 0.02 && $line_array[7] < 0.98) {
                $p2_baf_clust_sum += $line_array[7];
                $p2_baf_clust_sum_sq += $line_array[7]**2;
                $p2_baf_clust_ct++;
            }      
        }                
    }
    return() unless ($p1_stats[$BAF_SNP_CT] >= $SIZE_SMALL_WINDOW);
    
    (@ch_stats[$CHR,$START], @p2_stats[$CHR,$START], @p1_stats[$CHR,$START]) 
        = (@{$CURR_CHR_REFS[0]}[1,2]) x 3;
    $ch_stats[$BAF_SNP_CT] = $p2_stats[$BAF_SNP_CT] = $p1_stats[$BAF_SNP_CT];
    $ch_stats[$LRR_SNP_CT] = $p2_stats[$LRR_SNP_CT] = $p1_stats[$LRR_SNP_CT];
    $ch_stats[$NC_SNP_CT]  = $p2_stats[$NC_SNP_CT]  = $p1_stats[$NC_SNP_CT];
    $ch_stats[$BAF_REF] = \@ch_baf;
    $p1_stats[$BAF_REF] = \@p1_baf;
    $p2_stats[$BAF_REF] = \@p2_baf;
    
    my $p1_baf_mean = $p1_baf[9] / $p1_baf[10];
    my $p2_baf_mean = $p2_baf[9] / $p2_baf[10];
    my $ch_baf_mean = $ch_baf[9] / $ch_baf[10];
    map {$_ = 1 - $_ if ($_ < $p1_baf_mean)} @p1_mbaf;    
    map {$_ = 1 - $_ if ($_ < $p2_baf_mean)} @p2_mbaf; 
    map {$_ = 1 - $_ if ($_ < $ch_baf_mean)} @ch_mbaf;

    $p1_stats[$MBAF_MED] = median(\@p1_mbaf);
    $p2_stats[$MBAF_MED] = median(\@p2_mbaf);
    $ch_stats[$MBAF_MED] = median(\@ch_mbaf);
    $p1_stats[$LRR_MED]  = median(\@p1_lrr);
    $p2_stats[$LRR_MED]  = median(\@p2_lrr);
    $ch_stats[$LRR_MED]  = median(\@ch_lrr);    
    $ch_stats[$STDEV] = st_dev($ch_baf_clust_sum, $ch_baf_clust_sum_sq,
                                $ch_baf_clust_ct);
    $p1_stats[$STDEV] = st_dev($p1_baf_clust_sum, $p1_baf_clust_sum_sq,
                                $p1_baf_clust_ct);
    $p2_stats[$STDEV] = st_dev($p2_baf_clust_sum, $p2_baf_clust_sum_sq,
                                $p2_baf_clust_ct);
    return([\@ch_stats, \@p1_stats, \@p2_stats]);
}

sub manage_chromosome_arm {
    # This is the hub for the analysis. Each thread uses this function to 
    # start a cascade of function calls and to return the detected regions to 
    # the main program. 
    # Returns a reference to the final detected regions and a number of 
    # descriptive variables.
    ($LOCAL_CH_BAF_MEAN, $LOCAL_CH_mBAF_MED, $LOCAL_P1_mBAF_MED, 
     $LOCAL_P2_mBAF_MED, $LOCAL_CH_LRR_MED, $LOCAL_P1_LRR_MED, 
     $LOCAL_P2_LRR_MED, $BOUNDARY_EXTENSION) = (0) x 8;
    my (@ch_mbaf_array, @ch_lrr_array, @p1_mbaf_array, @p1_lrr_array, 
        @p2_mbaf_array, @p2_lrr_array, @remove);
    my ($ch_baf_sum, $ch_baf_ct, $p1_baf_sum, $p1_baf_ct, $p2_baf_sum, 
        $p2_baf_ct) = (0) x 6;

    $REFINING = 0;
    for (my $i = 0; $i < @CURR_CHR_REFS; $i++) {
        if (!$CURR_CHR_REFS[$i]->[3] || !$CURR_CHR_REFS[$i]->[6] 
            || !$CURR_CHR_REFS[$i]->[9]) {
            push(@remove, $i);
            next;
        }
        # Create hashes for SNP number calculations  
        $SNP_BY_NUM{keys %SNP_BY_POS} = $CURR_CHR_REFS[$i]->[POS];        
        $SNP_BY_POS{$CURR_CHR_REFS[$i]->[POS]} = keys %SNP_BY_POS; 

        if (substr($CURR_CHR_REFS[$i]->[4], 0, 1) eq "N" 
            || substr($CURR_CHR_REFS[$i]->[7], 0, 1) eq "N"
            || substr($CURR_CHR_REFS[$i]->[10],0, 1) eq "N") {
            push(@remove, $i);
            next;
        }
        # Create chromosome array including cnv markers for LRR calculations
        my $lrr_ref = [@{$CURR_CHR_REFS[$i]}[0..2,5,8,11]];
        push(@LRR_REFS, $lrr_ref);
        $LRR_BY_POS{$$lrr_ref[POS]} = $#LRR_REFS;    
                
        # Ignore intensity only markers for baf analyses
        if (substr($CURR_CHR_REFS[$i]->[0], 0, 2) =~ /cn/i) {
            push(@remove, $i);
        }
        else {
            my $ch_baf = $CURR_CHR_REFS[$i]->[10];
            if ($ch_baf > $AA_BOUND && $ch_baf < $BB_BOUND) {
                $ch_baf_sum += $ch_baf;
                $ch_baf_ct++;
            }
        }
    }
   
    delete @CURR_CHR_REFS[@remove] if @remove;
    my @temp;
    map {push(@temp, $_) if $_} @CURR_CHR_REFS;
    @CURR_CHR_REFS = @temp;
    return() unless ($#CURR_CHR_REFS >= $SIZE_SMALL_WINDOW);
    $LOCAL_CH_BAF_MEAN = $ch_baf_sum / $ch_baf_ct;
    (@ch_mbaf_array, @ch_lrr_array, @p1_mbaf_array, @p1_lrr_array, 
     @p1_mbaf_array, @p1_lrr_array, @remove) = (()) x 7;
    $ch_baf_sum = $ch_baf_ct = 0;

    ############################## Process chromosome ##########################
    my @final_regions = @{process_chromosome_arm()};
    my (@baf_normal, @lrr_normal);
    # Adjusted for 2% near centromere and telomere
    my $total_snps = $#CURR_CHR_REFS * 0.98;
    my $abn_snps = 0;
    if (@final_regions) {
        foreach my $ref (@final_regions) {$abn_snps += $$ref[5]}
        my (@baf_abn_indices, @lrr_abn_indices);
        # 75% threshold adjusted for 2% near centromere and telomere
        if ($abn_snps < $total_snps * 0.75) {
            foreach my $ref (@final_regions) {
                my ($baf_start_ind, $baf_stop_ind, $baf_success) 
                    = find_closest_indices(@$ref[START,STOP], \%BAF_BY_POS);
                if ($baf_success) {
                    push(@baf_abn_indices, $baf_start_ind..$baf_stop_ind);
                }
                my ($lrr_start_ind, $lrr_stop_ind, $lrr_success) 
                    = find_closest_indices(@$ref[START,STOP], \%LRR_BY_POS);
                if ($lrr_success) {
                    push(@lrr_abn_indices, $lrr_start_ind..$lrr_stop_ind);        
                }
            }
        }
        else {
            # If most of the chromosome arm is abnormal analyses are not 
            # repeated using local medians.
            # Indicate single pass through analyses and return results.
            map {$$_[19] = 1} @final_regions;
            return([\@final_regions]);
        }

        if (@baf_abn_indices) {
            # Remove 1% of SNPs near ends
            my $start = ($#CURR_CHR_REFS * 0.01) + 1;
            my $end = ($#CURR_CHR_REFS * 0.99) - 1;
            my @baf_temp = @CURR_CHR_REFS[$start..$end];
            delete @baf_temp[@baf_abn_indices];
            map {push(@baf_normal, $_) if $_} @baf_temp;
            @baf_temp = (); 
        }
        else {@baf_normal = @CURR_CHR_REFS}
                    
        my @lrr_temp = @LRR_REFS;
        delete @lrr_temp[@lrr_abn_indices];
        map {push(@lrr_normal, $_) if $_} @lrr_temp;
        @lrr_temp = (); 
    }
    else {
        @baf_normal = @CURR_CHR_REFS;
        @lrr_normal = @LRR_REFS;
    }       
    my ($aa_sum, $aa_sum_squares, $aa_count, $ab_sum, $ab_sum_squares,
        $ab_count, $bb_sum, $bb_sum_squares, $bb_count);
    foreach my $ref (@baf_normal) {
        # Store values for BAF StDev calculations
        my $ch_baf = $$ref[10];   
        if ($ch_baf < 0.75) {
            if ($ch_baf <= 0.25) {
                $aa_sum += $ch_baf;
                $aa_sum_squares += $ch_baf**2;
                $aa_count++;
            }
            else {
                $ab_sum += $ch_baf;
                $ab_sum_squares += $ch_baf**2;
                $ab_count++;
            }
            if ($ch_baf > $AA_BOUND) {
                push(@ch_mbaf_array, $ch_baf);
                $ch_baf_sum += $ch_baf;
                $ch_baf_ct++;
            }
        }
        else {
            $bb_sum += $ch_baf;
            $bb_sum_squares += $ch_baf**2;
            $bb_count++;
            if ($ch_baf < $BB_BOUND) {
                push(@ch_mbaf_array, $ch_baf) ;
                $ch_baf_sum += $ch_baf;
                $ch_baf_ct++;
            }                
        }
        my $baf = $$ref[4];
        if ($baf > $P1_AA_BOUND && $baf < $P1_BB_BOUND) {
            push(@p1_mbaf_array, $baf); 
            $p1_baf_sum += $baf;
            $p1_baf_ct++;            
        }
        $baf = $$ref[7];
        if ($baf > $P2_AA_BOUND && $baf < $P2_BB_BOUND) {
            push(@p2_mbaf_array, $baf);                
            $p2_baf_sum += $baf;
            $p2_baf_ct++;            
        }
    }
    
    #Store initial thresholds
    my @temp_bounds = ($AA_BOUND, $BB_BOUND, $AB_UPPER_BOUND, $AB_LOWER_BOUND,
        $MI1_UPPER_THRESH, $MI1_LOWER_THRESH, $ACCEPTABLE_ERRORS);
    # Calculate standard deviation for BAFs and BAF thresholds for analysis

    my @AA_stats = st_dev($aa_sum, $aa_sum_squares, $aa_count);
    my @BB_stats = st_dev($bb_sum, $bb_sum_squares, $bb_count);
    my @AB_stats = st_dev($ab_sum, $ab_sum_squares, $ab_count);
    $AA_BOUND         = $AA_stats[0] + ($AA_stats[1] * $HOM_SD);
    $BB_BOUND         = $BB_stats[0] - ($BB_stats[1] * $HOM_SD);
    $AB_UPPER_BOUND   = $AB_stats[0] + ($AB_stats[1] * $HET_SD);
    $AB_LOWER_BOUND   = $AB_stats[0] - ($AB_stats[1] * $HET_SD);
    $MI1_UPPER_THRESH = $AB_stats[0] + ($AB_stats[1] * 5);
    $MI1_LOWER_THRESH = $AB_stats[0] - ($AB_stats[1] * 5);
    @AA_stats = @AB_stats = @BB_stats = ();
    foreach my $lrr_ref (@lrr_normal) {
        push(@p1_lrr_array, $$lrr_ref[3]);
        push(@p2_lrr_array, $$lrr_ref[4]);
        push(@ch_lrr_array, $$lrr_ref[5]);
    }
    
    $LOCAL_P1_LRR_MED     = median(\@p1_lrr_array);
    $LOCAL_P2_LRR_MED     = median(\@p2_lrr_array);
    $LOCAL_CH_LRR_MED     = median(\@ch_lrr_array);    
    $LOCAL_CH_BAF_MEAN    = $ch_baf_sum / $ch_baf_ct;
    my $local_p1_baf_mean = $p1_baf_sum / $p1_baf_ct;
    my $local_p2_baf_mean = $p2_baf_sum / $p2_baf_ct;
    map {$_ = 1 - $_ if ($_ < $LOCAL_CH_BAF_MEAN)} @ch_mbaf_array;
    map {$_ = 1 - $_ if ($_ < $local_p1_baf_mean)} @p1_mbaf_array;    
    map {$_ = 1 - $_ if ($_ < $local_p2_baf_mean)} @p2_mbaf_array; 
    $LOCAL_CH_mBAF_MED = median(\@ch_mbaf_array);
    $LOCAL_P1_mBAF_MED = median(\@p1_mbaf_array);
    $LOCAL_P2_mBAF_MED = median(\@p2_mbaf_array);
    
    my $chr = $baf_normal[0]->[CHR];
    (@baf_normal, @lrr_normal, @p1_lrr_array, @p2_lrr_array, @ch_lrr_array,
     @ch_mbaf_array) = (()) x 6; 

    $REFINING = 1;
    # Reprocess chromosome
    my $return = process_chromosome_arm();
    map {$$_[19] = 2} @$return;
    
    # Restore initial thresholds     
    ($AA_BOUND, $BB_BOUND, $AB_UPPER_BOUND, $AB_LOWER_BOUND, 
     $MI1_UPPER_THRESH, $MI1_LOWER_THRESH, $ACCEPTABLE_ERRORS) = @temp_bounds;
    (%LRR_BY_POS, %SNP_BY_POS, %SNP_BY_NUM) = (()) x 3; 
    (@CURR_CHR_REFS, @LRR_REFS) = (()) x 2;
    return([$return, $chr, [($LOCAL_P1_mBAF_MED, $LOCAL_P2_mBAF_MED,
        $LOCAL_CH_mBAF_MED, $LOCAL_P1_LRR_MED, $LOCAL_P2_LRR_MED, 
        $LOCAL_CH_LRR_MED)]]);
}

sub max_window {
    my @max_window = (@{$CURR_CHR_REFS[0]}[CHR,POS], 0);
    my ($p1_ct, $p2_ct, $snp_ct) = (0) x 3;
    my @bafsnp_array;
    for (my $i = 0; $i < @CURR_CHR_REFS; $i++) {
        next if (!$CURR_CHR_REFS[$i]->[3] || !$CURR_CHR_REFS[$i]->[6] 
            || !$CURR_CHR_REFS[$i]->[9]);
        next if ((substr($CURR_CHR_REFS[$i]->[0], 0, 2) =~ /cn/i)
            || substr($CURR_CHR_REFS[$i]->[4],  0, 1) eq "N" 
            || substr($CURR_CHR_REFS[$i]->[7],  0, 1) eq "N" 
            || substr($CURR_CHR_REFS[$i]->[10], 0, 1) eq "N");
        $snp_ct++;
        my $results = find_informative_snps($CURR_CHR_REFS[$i]);
        if ($snp_ct < $SIZE_SMALL_WINDOW) {
            if ($CURR_CHR_REFS[$i]->[-1] == FATHER) {$p1_ct++} 
            elsif ($CURR_CHR_REFS[$i]->[-1] == MOTHER) {$p2_ct++}
            push(@bafsnp_array, $CURR_CHR_REFS[$i]);
        }
        else {
            push(@bafsnp_array, $CURR_CHR_REFS[$i]);        
            if ($bafsnp_array[-1]->[-1] == FATHER) {$p1_ct++} 
            elsif ($bafsnp_array[-1]->[-1] == MOTHER) {$p2_ct++}
            $max_window[2] = $p1_ct + $p2_ct if($p1_ct+$p2_ct>$max_window[2]);
            if ($bafsnp_array[0]->[-1] == FATHER) {$p1_ct--}
            if ($bafsnp_array[0]->[-1] == MOTHER) {$p2_ct--}
            shift(@bafsnp_array);
        }            
    }
    return(\@max_window);
}

sub median {
    my $array_ref = $_[0];
    my $median = 0;
    if (!$$array_ref[0] || @$array_ref == 1) {$median = NA}
    elsif (@$array_ref == 2) {$median = ($$array_ref[0] + $$array_ref[1]) / 2}
    else {
        my @sorted_array = sort { $a <=> $b } @$array_ref;
        my $mid_index = @sorted_array / 2;
        if (@sorted_array % 2) {
            $median = $sorted_array[int(($mid_index) - 0.5)];
        }
        else {
            $median = ($sorted_array[$mid_index] 
                + $sorted_array[$mid_index - 1]) / 2;
        }
    }
    return($median);
}

sub overlap {
    my ($start1, $stop1, $start2, $stop2) = @_;
    my ($overlap_start, $overlap_stop, $end1_start, $end1_stop, 
        $end2_start, $end2_stop, $prev_snp, $next_snp, $type) = (0) x 9;
    # -----
    # -----       
    if ($start1 == $start2 && $stop1 == $stop2) {
        ($overlap_start, $overlap_stop)  = ($start1, $stop1);
        $type = 1;
    }
    #   -----
    # -----
    elsif ($start1 > $start2 && $stop2 >= $start1 && $stop2 < $stop1) {
        $prev_snp = $SNP_BY_NUM{$SNP_BY_POS{$start1} - 1};
        $next_snp = $SNP_BY_NUM{$SNP_BY_POS{$stop2} + 1};
        ($end1_start, $end1_stop)       = ($start2, $prev_snp);
        ($overlap_start, $overlap_stop) = ($start1, $stop2);
        ($end2_start, $end2_stop)       = ($next_snp, $stop1);
        $type = 2;
    }
    # -----
    # ---
    elsif ($start1 == $start2 && $stop2 > $start1 && $stop2 < $stop1) {
        $next_snp = $SNP_BY_NUM{$SNP_BY_POS{$stop2} + 1};
        ($overlap_start, $overlap_stop) = ($start2, $stop2);
        ($end2_start, $end2_stop)       = ($next_snp, $stop1);
        $type = 3;
    }
    #   ---
    # -----
    elsif ($start1 > $start2 && $stop2 > $start1 && $stop2 == $stop1) {
        $prev_snp = $SNP_BY_NUM{$SNP_BY_POS{$start1} - 1};
        ($end1_start, $end1_stop)       = ($start2, $prev_snp);
        ($overlap_start, $overlap_stop) = ($start1, $stop1);
        $type = 4;
    }
    #  ---
    # -----
    elsif ($start1 > $start2 && $stop1 < $stop2) { 
        $prev_snp = $SNP_BY_NUM{$SNP_BY_POS{$start1} - 1};
        $next_snp = $SNP_BY_NUM{$SNP_BY_POS{$stop1}  + 1};        
        ($end1_start, $end1_stop)       = ($start2, $prev_snp);
        ($overlap_start, $overlap_stop) = ($start1, $stop1);
        ($end2_start, $end2_stop)       = ($next_snp, $stop2);
        $type = 5;        
    }
    # -----
    #   -----
    elsif ($start2 > $start1 && $stop1 >= $start2 && $stop1 < $stop2) {
        $prev_snp = $SNP_BY_NUM{$SNP_BY_POS{$start2} - 1};
        $next_snp = $SNP_BY_NUM{$SNP_BY_POS{$stop1}  + 1};        
        ($end1_start, $end1_stop)       = ($start1, $prev_snp);
        ($overlap_start, $overlap_stop) = ($start2, $stop1);
        ($end2_start, $end2_stop)       = ($next_snp, $stop2);
        $type = 6;        
    }
    # ---
    # -----
    elsif ($start2 == $start1 && $stop1 > $start2 && $stop1 < $stop2) {
        $next_snp = $SNP_BY_NUM{$SNP_BY_POS{$stop1}  + 1};        
        ($overlap_start, $overlap_stop) = ($start1, $stop1);
        ($end2_start, $end2_stop)       = ($next_snp, $stop2);
        $type = 7;        
    }
    # -----
    #   ---
    elsif ($start2 > $start1 && $stop1 > $start2 && $stop1 == $stop2) {
        $prev_snp = $SNP_BY_NUM{$SNP_BY_POS{$start2} - 1};
        ($end1_start, $end1_stop)       = ($start1, $prev_snp);
        ($overlap_start, $overlap_stop) = ($start2, $stop2);
        $type = 8;        
    }
    # -----
    #  ---
    elsif ($start2 > $start1 && $stop2 < $stop1) {
        $prev_snp = $SNP_BY_NUM{$SNP_BY_POS{$start2} - 1};
        $next_snp = $SNP_BY_NUM{$SNP_BY_POS{$stop2}  + 1};        
        ($end1_start, $end1_stop)       = ($start1, $prev_snp);
        ($overlap_start, $overlap_stop) = ($start2, $stop2);
        ($end2_start, $end2_stop)       = ($next_snp, $stop1);
        $type = 9;
    }
    return([$type, $end1_start, $end1_stop, $overlap_start, $overlap_stop,
          $end2_start, $end2_stop]);
}

sub overlapping_windows {
    # Counts and stores parental contribution for each overlapping window
    my ($array_ref, $window_size) = @_;
    my ($p1_ct, $p2_ct, $index_count) = (0) x 3;
    my $new_info = 1;
    my $chr = $$array_ref[0]->[CHR];
    my $array_size = scalar(@$array_ref);
    my @return_refs;
    for (my $i = 0; $i < $window_size; $i++) {   
        my $curr_ref = $$array_ref[$i];
        if ($$curr_ref[-1] == FATHER) {$p1_ct++} 
        elsif ($$curr_ref[-1] == MOTHER) {$p2_ct++}
    }
    
    while ($index_count + $window_size < $array_size) {     
        my $start_ref = $$array_ref[$index_count];
        my $end_ref = $$array_ref[$index_count + $window_size];
        if ($$end_ref[-1] == FATHER) {$p1_ct++; $new_info = 1} 
        elsif ($$end_ref[-1] == MOTHER) {$p2_ct++; $new_info = 1} 
        push(@return_refs, [$chr, $$start_ref[POS], $$end_ref[POS], 
            $p1_ct, $p2_ct, $new_info]);
        $new_info = 0;
        if ($$start_ref[-1] == FATHER) {$p1_ct--; $new_info = 1} 
        elsif ($$start_ref[-1] == MOTHER) {$p2_ct--; $new_info = 1}         
        $index_count++;  
    }
    return(\@return_refs)
}

sub parse_file {
    my ($INPUT_FILE, $function, $info) = @_;
    my $check_input = 1 if ($info && $info eq "check input");
    my $main = 1 if ($info && $info eq "main"); 
    my ($curr_chr, $centro_start, $centro_stop, $gap, $thread_list);
    my (@results, @threads);
    if(!open(INPUT_FILE, "<", $INPUT_FILE)) { 
        print "Could not open input file: $INPUT_FILE!\n";
        print STDERR "Could not open input file: $INPUT_FILE!\n";            
        return();
    }

    while (<INPUT_FILE>) {
        my $line = $_;
        chomp($line);
        my @line_array = split(/\t/, $line);
        next if (!$line_array[0]);
        #Setup
        if ($. > 1) {
            # Convert chromosome identifier to numeric 
            if ($line_array[CHR] !~ /^\d/) {
                if ($main && $GENDER && $line_array[CHR] eq "X"){
                    $line_array[CHR] = CHR_X;
                }
                else {next}
            }
            $curr_chr = $line_array[CHR];
            if ($main && $verbose) { 
                print "Analyzing Chromosome ", $curr_chr, "\n";
            }
            # Store centromere info
            ($centro_start, $centro_stop) = @{$CENTRO_REFS[$curr_chr]}[1,2];
            $gap = 1 if ($line_array[POS] < $centro_stop);
            push(@CURR_CHR_REFS, \@line_array); 
            last;
        }
    }        

    while (<INPUT_FILE>) {
        my $line = $_;
        chomp($line);
        my @line_array = split(/\t/, $line);
        next if (!$line_array[0]);
        
        # Convert chromosome identifier to numeric 
        if ($line_array[CHR] !~ /^\d/) {
            if ($main && $GENDER && $line_array[CHR] eq "X"){
                $line_array[CHR] = CHR_X;
            }
            else {next}
        }

        # Start new thread when centromere or new chromosome is encountered                
        if ($curr_chr != $line_array[CHR]
            || ($gap && $line_array[POS] > $centro_start)) { 
            $gap = 0;
            
            if (@CURR_CHR_REFS >= $SIZE_SMALL_WINDOW) {
                # Create threads only if >1 processor is available                
                if ($CORES_PER_SAMPLE == 1) {
                    my $results = $function->();
                    push(@results, $results) if $results;
                }
                else {            
                    while (threads->list(threads::running)
                        >= $CORES_PER_SAMPLE) {usleep(100000)}
                    for (0..$#threads) { 
                        if ($threads[$_]->is_joinable()) {
                            my $results = $threads[$_]->join();
                            if ($results) {push(@results, $results)}
                        }               
                    }
                    $thread_list = threads->new($function);
                    push(@threads, $thread_list);
                } 
            }
            if ($curr_chr != $line_array[CHR]) {
                $curr_chr = $line_array[CHR];    
                if ($main && $verbose) { 
                    print "Analyzing Chromosome ", $curr_chr, "\n";
                }
                # Store centromere info
                ($centro_start, $centro_stop) =
                    @{$CENTRO_REFS[$curr_chr]}[1,2];
                $gap = 1 if ($line_array[POS] < $centro_stop);                    
            }
            @CURR_CHR_REFS = ();
        }
        push(@CURR_CHR_REFS, \@line_array); 
        if ($check_input && $. == 1000) {
            my $return = check_input_format(\@CURR_CHR_REFS);
            return() unless $return;
        }
    }
    close(INPUT_FILE);
    if (@CURR_CHR_REFS >= $SIZE_SMALL_WINDOW) {
        if ($CORES_PER_SAMPLE == 1) {
            my $results = $function->();
            push(@results, $results) if $results;
        }
        else {            
            while (threads->list(threads::running) 
                >= $CORES_PER_SAMPLE) {usleep(100000)}
            for (0..$#threads) { 
                if ($threads[$_]->is_joinable()) {
                    my $results = $threads[$_]->join();
                    if ($results) {push(@results, $results)}
                }               
            }
            $thread_list = threads->new($function);
            push(@threads, $thread_list);
        }
    }
    if ($CORES_PER_SAMPLE > 1) {
        # Check if threads are finished
        while (threads->list(threads::running)) {usleep(100000)}
        for (0..$#threads) { 
            if ($threads[$_]->is_joinable()) {
                my $results = $threads[$_]->join();
                if ($results) {push(@results, $results)}
            }               
        }
    }
    @CURR_CHR_REFS = ();
    return(\@results);
}

sub print_proper_format {
    # Prints a message describing the proper format of the input file.
    print <<END;
    
The input file must be tab delimited, sorted by chromosome and position, 
  and in the following order: SNP Name, Chromosome, Position, 
  Father GType, Father BAF, Father LRR, Mother GType, Mother BAF, Mother LRR, 
  Child GType, Child BAF, Child LRR.
The genotypes must be AA, AB, BB, NC or NoCall.
B allele frequencies must be >= 0 and <= 1 for polymorphic markers.
A header line is expected but is not used to determine column identity.
    
END
}    

sub process_chromosome_arm {
    my ($adjusted_mi1_ct, $chr, $end, $hom_reg_ref, $inf_SNP_count, $mi1, 
        $mi1_contrib, $mi1_ct, $mi1_reg_ref, $POD_ref, $PODcr_ref, $podcr_regs,
        $start, $curr_mi1_het_ct, $out_start, $out_end, $outlier_ct) = (0) x 17;
    my (@abnormal_regions, @adj_mi1_regions, @return, @large_regions, 
        @outlier_regions);
    (@BAF_OUTLIERS, @P1_INF_SNPS, @P2_INF_SNPS) = (()) x 3;
    (%BAF_BY_POS, %BAF_OUTLIERS_BY_POS,%CH_HD_BY_POS,%MI1_BY_POS,%P1_HD_BY_POS, 
     %P1_INF_SNP_BY_POS, %P2_HD_BY_POS, %P2_INF_SNP_BY_POS) = (()) x 8;
    $chr = $CURR_CHR_REFS[0]->[CHR];

    for (my $i = 0; $i < @CURR_CHR_REFS; $i++) {
        $BAF_BY_POS{$CURR_CHR_REFS[$i]->[POS]} = $i;
        my $results = find_informative_snps($CURR_CHR_REFS[$i]);
        if ($$results[0]) {
            push(@P1_INF_SNPS, $CURR_CHR_REFS[$i]);
            $P1_INF_SNP_BY_POS{$CURR_CHR_REFS[$i]->[POS]} 
                = keys %P1_INF_SNP_BY_POS;
        }
        elsif ($$results[1]) {
            push(@P2_INF_SNPS, $CURR_CHR_REFS[$i]);
            $P2_INF_SNP_BY_POS{$CURR_CHR_REFS[$i]->[POS]} 
                = keys %P2_INF_SNP_BY_POS;
        }
        $inf_SNP_count += $$results[2];
        unless ($chr == CHR_X && $GENDER == MALE) {
            if ($$results[3] && ($$results[4] == FATHER 
                || $$results[4] == MOTHER)) {
                if ($mi1_ct && $mi1_contrib != $$results[4]) {
                    if ($mi1_ct >= $MIN_MI1) {
                        push(@adj_mi1_regions,[$chr, $start, 
                            $CURR_CHR_REFS[$i]->[POS], " ",$mi1_contrib,$mi1_ct, 
                            $mi1_ct]);
                    }
                    $mi1_ct = $start = 0;
                }
                $mi1_ct++;
                $mi1_contrib = $$results[4];
                $MI1_BY_POS{$CURR_CHR_REFS[$i]->[POS]} = keys %MI1_BY_POS;
                if ($mi1_ct == 1) {$start = $end = $CURR_CHR_REFS[$i]->[POS]}
                else {$end = $CURR_CHR_REFS[$i]->[POS]}
                $adjusted_mi1_ct++ if !$mi1_ct;
            }
            else {
                if ($mi1_ct >= $MIN_MI1) {
                    push(@adj_mi1_regions, [$chr, $start, $end," ",$mi1_contrib, 
                        $mi1_ct, $mi1_ct]);
                }
                $mi1_ct = $start = $end = $mi1_contrib = 0;
            }
        }
        if ($$results[5]) {
            $outlier_ct++;
            $BAF_OUTLIERS_BY_POS{$CURR_CHR_REFS[$i]->[POS]} 
                = keys %BAF_OUTLIERS_BY_POS;
            push(@BAF_OUTLIERS, $CURR_CHR_REFS[$i]);    
            if ($outlier_ct == 1) {$out_start = $out_end 
                = $CURR_CHR_REFS[$i]->[POS]}
            else {$out_end = $CURR_CHR_REFS[$i]->[POS]}
        }
        elsif ($CURR_CHR_REFS[$i]->[9] eq "AB") {
            if ($outlier_ct >= $MIN_BAF_OUT) {
                push(@outlier_regions, [$chr, $out_start, $out_end, 
                    " ", UNKNOWN, $outlier_ct, $outlier_ct]);
            }
            $outlier_ct = $out_start = $out_end = 0;
        }
    }    

    unless ($chr == CHR_X && $GENDER == MALE) {
        if ($mi1_ct >= $MIN_MI1) {
            push(@adj_mi1_regions, [$chr, $start, $end, " ", $mi1_contrib, 
                $mi1_ct, $mi1_ct]);
        }
    }
    if ($outlier_ct >= $MIN_BAF_OUT) {
        push(@outlier_regions, [$chr, $out_start, $out_end, 
            " ", UNKNOWN, $outlier_ct, $outlier_ct]);
    }
    if ($POD) {
        $POD_ref = find_POD_regions(\@CURR_CHR_REFS, $SIZE_SMALL_WINDOW, 
            $ACCEPTABLE_ERRORS);
        my @indices;
        if (@$POD_ref) {
            push(@abnormal_regions, @$POD_ref);
            my @sorted_regions = sort {$$a[1] <=> $$b[1]} @abnormal_regions;
            for (@abnormal_regions) {
                push(@indices, [find_closest_indices(@$_[START,STOP], 
                    \%BAF_BY_POS)]);
            }
        }
        unless ($chr == CHR_X) {
            # Send normal regions to large window POD detection
            if (@indices) {
                my @sorted_indices = sort {$$a[0] <=> $$b[0]} @indices;
                my ($start_norm, $end_norm) = (0) x 2;
                for (my $j = 0; $j <= $#sorted_indices; $j++) {
                    my ($start_abn, $end_abn) = @{$sorted_indices[$j]}[0,1];
                    $end_norm = $start_abn - 1;
                    $end_norm ||= 0;

                    if ($end_norm - $start_norm >= $SIZE_LARGE_WINDOW) {
                        my $large_POD_ref = find_POD_regions(
                            [@CURR_CHR_REFS[$start_norm..$end_norm]], 
                            $SIZE_LARGE_WINDOW, $LARGE_WIN_ACCEPTABLE_ERRORS);
                        if (@$large_POD_ref) {    
                            push(@abnormal_regions, @$large_POD_ref);
                        }
                    }
                    $start_norm = $end_abn + 1;
                    if ($start_norm > $#CURR_CHR_REFS) {
                        $start_norm = $#CURR_CHR_REFS;
                    }
                }
                $end_norm = $#CURR_CHR_REFS;
                if ($end_norm - $start_norm >= $SIZE_LARGE_WINDOW) {
                    my $large_POD_ref = find_POD_regions(
                        [@CURR_CHR_REFS[$start_norm..$end_norm]], 
                        $SIZE_LARGE_WINDOW, $LARGE_WIN_ACCEPTABLE_ERRORS);
                    if (@$large_POD_ref) {        
                        push(@abnormal_regions, @$large_POD_ref);
                    }
                }
            }
            else {
                my $large_POD_ref = find_POD_regions(\@CURR_CHR_REFS, 
                    $SIZE_LARGE_WINDOW, $LARGE_WIN_ACCEPTABLE_ERRORS);
                push(@abnormal_regions, @$large_POD_ref) if (@$large_POD_ref);
            }
            collapse_region_list(\@abnormal_regions);
        }
    }
    if ($PODcr && @outlier_regions) {
        $podcr_regs = define_PODcr_regions(\@outlier_regions);
        push(@abnormal_regions, @$podcr_regs) if @$podcr_regs;
    }      
    if ($MI1 && @adj_mi1_regions && ($chr != CHR_X && $GENDER != MALE)) {
        $mi1_reg_ref = define_mi1_regions(\@adj_mi1_regions);
        push(@abnormal_regions, @$mi1_reg_ref) if @$mi1_reg_ref;
    }  
    if ($HD) {
        $hom_reg_ref = detect_hom_del();
        push(@abnormal_regions, @$hom_reg_ref) if @$hom_reg_ref;
    }
    
    if (@abnormal_regions) {        
        count_std_mi1_overlap(\@abnormal_regions);
        collapse_region_list(\@abnormal_regions);
        map {$_ = evaluate_region($_) unless $$_[6]} @abnormal_regions;
        map {calculate_region_stats($_) unless $$_[7]} @abnormal_regions;
        count_std_mi1_overlap(\@abnormal_regions);    
        collapse_region_list(\@abnormal_regions);
        map{$$_[6] = count_SNPs(@$_[START,STOP,16]) unless $$_[6]} 
            @abnormal_regions;
        map {calculate_region_stats($_) unless $$_[7]} @abnormal_regions;
        my $different;
        for (my $i = 0; $i <= $#abnormal_regions; $i++) {
            my $ref = $abnormal_regions[$i];
            next if ($$ref[0] == CHR_X && $GENDER == MALE && $$ref[3] eq "DEL");
            my $curr_mi1_het_ct = 0;
            my $pc_index;
            if ($$ref[4] == FATHER) {$pc_index = 3}
            elsif ($$ref[4] == MOTHER) {$pc_index = 6}
            if ($pc_index && $$ref[3] eq " ") {
                my ($cusum, $max, $min, $snp_ct, $match, $homz, $inf_ct, $hUPD) = (0) x 8;
                my @cusum;
                my ($reg_start, $reg_stop, $success) 
                    = find_closest_indices(@$ref[START,STOP] ,\%BAF_BY_POS);

                if ($success && $reg_stop - $reg_start > 0) {
                    for ($reg_start..$reg_stop) {
                        my $curr = $CURR_CHR_REFS[$_];    
                        my @res = @{find_informative_snps($curr)};
                        if ($res[6]) {
                            $hUPD++; 
                            $cusum++;
                            $inf_ct++;
                            push(@cusum, [$_,$cusum]);
                        }
                        elsif ($res[2]) {
                            $cusum--; 
                            $inf_ct++;
                            push(@cusum, [$_,$cusum]);
                        }
                        $max = $cusum if ($cusum > $max);
                        $min = $cusum if ($cusum < $min);
                        if (substr($$curr[$pc_index], 0, 1) ne "N"
                            && substr($$curr[9], 0, 1) ne "N") {
                            $snp_ct++;
                            $match++ if ($$curr[$pc_index] eq $$curr[9]);
                            $homz++ if ($$curr[9] eq "AA" || $$curr[9] eq "BB");
                        }
                    }
                    my $curr_hUPD_rate = 0;
                    $curr_hUPD_rate = $hUPD / $inf_ct if $inf_ct;
                    if ($curr_hUPD_rate > $hUPD_RATE) {
                        my $prob = bin_test($hUPD, $inf_ct, $hUPD_RATE, "g");
                        if ($prob <= 0.001) {
                            my $peak = $min > -($$ref[6] * 0.025) && $min < ($$ref[6] * 0.025)
                                && $max > $cusum + ($$ref[6] * 0.025);
                            my $valley = $max > $cusum - ($$ref[6] * 0.025) && $max < $cusum 
                                + ($$ref[6] * 0.025) && $min < -($$ref[6] * 0.025);
                            my $change_pt;
                            if ($peak) {$change_pt = $max}
                            elsif ($valley) {$change_pt = $min}
                            
                            if ($peak || $valley) {
                                my @reg1 = my @reg2 = @$ref;
                                my ($index, $min_size);
                                for (@cusum) {if ($$_[1] == $change_pt) {$index = $$_[0]}}
                                $reg1[STOP] = $CURR_CHR_REFS[$index]->[POS];
                                $reg2[START] = $CURR_CHR_REFS[$index+1]->[POS];    
                                calculate_region_stats(\@reg1);
                                calculate_region_stats(\@reg2);
                                $reg1[6] = count_SNPs(@reg1[START,STOP,16]);
                                $reg2[6] = count_SNPs(@reg2[START,STOP,16]);                                
                                if    ($$ref[17] == POD) {$min_size = $MIN_POD}
                                elsif ($$ref[17] == MI1) {$min_size = $MIN_MI1}
                                push(@abnormal_regions, \@reg1) if $reg1[6] >= $min_size;
                                push(@abnormal_regions, \@reg2) if $reg2[6] >= $min_size;
                                next;
                            }
                            else {
                                if ($match >= $snp_ct - ceil($snp_ct * $ERROR_RATE)) {
                                    $$ref[3] = "UPhD";
                                }
                            }
                        }
                    }
                    if ($$ref[3] eq " ") {
                        if ($snp_ct && $homz / $snp_ct >= 0.9 
                        || ($$ref[10] =~ /^[+-]?(?=\.?\d)\d*\.?\d*(?:e[+-]?\d+)?\z/i 
                        && $$ref[10] > 0.55 
                        && $$ref[11] =~ /^[+-]?(?=\.?\d)\d*\.?\d*(?:e[+-]?\d+)?\z/i 
                        && $$ref[11] < 0.05 && $$ref[11] > -0.05)) {$$ref[3] = "UPiD"}
                    }
                }
            }
            # Check for MI1s indicative of abn region in parent instead of child
            my ($reg_start, $reg_stop, $success) 
                = find_closest_indices(@$ref[START,STOP], \%BAF_BY_POS);
            if ($success && $reg_stop - $reg_start > 0) {
                for (my $j = $reg_start; $j <= $reg_stop; $j++) {
                    $curr_mi1_het_ct++ if ($CURR_CHR_REFS[$j]->[-1] 
                        == UNKNOWN_MI1);
                }
                if($curr_mi1_het_ct / ($reg_stop - $reg_start) > $HET_MI1_RATE){
                    my($n, $k, $p, $q)=($reg_stop-$reg_start,$curr_mi1_het_ct, 
                        $HET_MI1_RATE, 1 - $HET_MI1_RATE);
                    my $prob = bin_test($k, $n, $p, "g");
                    next if ($prob < 0.001);
                }
            }

            if ($$ref[17] == HD) {
                my $hd_thresh;
                if    ($$ref[4] == FATHER) {$hd_thresh = $MIN_P2_HD}
                elsif ($$ref[4] == MOTHER) {$hd_thresh = $MIN_P1_HD}
                elsif ($$ref[4] == NONE)   {$hd_thresh = $MIN_CH_HD}
                elsif ($$ref[4] == BOTH)   {$hd_thresh = 
                    ($MIN_P1_HD <= $MIN_P2_HD) ? $MIN_P1_HD : $MIN_P2_HD}
                push(@return, $ref) if ($$ref[6] >= $hd_thresh);
            }
            elsif ($$ref[17]==MI1) {push(@return,$ref) if ($$ref[6]>=$MIN_MI1)}
            elsif ($$ref[17] == POD || $$ref[17] == PODcr) {push(@return,$ref)}
        }
        map {$$_[20] = call_inheritance($_)} @return;
    }
    my $results = \@return;
    return($results);
}

sub secondary_counts {
    # Returns various counts and the scan statistic.
    my ($snp_ct, $prev_p1_hd, $prev_p2_hd, $prev_ch_hd, $prev_MI1, 
        $prev_outlier, $outlier_ct, $p1_ct, $p2_ct) = (0) x 9;
    my (@stats, @bafsnp_array);
    # Initialize array
    # The following are *_stats array indices
    my ($CHR, $START, $SCAN_STAT, $MI1_CT, $INF_SNP_CT, $P1_HD_CT, 
        $P2_HD_CT, $CH_HD_CT, $HET_MI1_CT, $BAF_OUTLIERS, $AB_CT,
        $ISO_HET, $INF_hUPD) = (0..12);
    @stats[$CHR,$START] = @{$CURR_CHR_REFS[0]}[1,2];
    $stats[$_] = 0 for (2..12); 
        
    for (my $i = 0; $i < @CURR_CHR_REFS; $i++) {
        next if (!$CURR_CHR_REFS[$i]->[3] || !$CURR_CHR_REFS[$i]->[6] 
            || !$CURR_CHR_REFS[$i]->[9]);
        next if ((substr($CURR_CHR_REFS[$i]->[0], 0, 2) =~ /cn/i)
              || substr($CURR_CHR_REFS[$i]->[4],  0, 1) eq "N" 
              || substr($CURR_CHR_REFS[$i]->[7],  0, 1) eq "N" 
              || substr($CURR_CHR_REFS[$i]->[10], 0, 1) eq "N");
        $snp_ct++;
        my $results = find_informative_snps($CURR_CHR_REFS[$i]);
        $stats[$INF_SNP_CT] += $$results[2];
        $stats[$ISO_HET]++  if $$results[6];
        $stats[$INF_hUPD]++ if $$results[7];
        
        # Count nonconsecutive SNPs with MI1 (or imputed)
        if ($$results[3]) {
            $stats[$HET_MI1_CT]++ if ($$results[4] == NONE);
            $stats[$MI1_CT]++ unless $prev_MI1;
            $stats[$MI1_CT]-- if ($prev_MI1 == 1);
            $prev_MI1++;
        }
        else {$prev_MI1 = 0}
        if ($$results[5]) {
            $stats[$BAF_OUTLIERS]++ unless $prev_outlier;
            $stats[$BAF_OUTLIERS]-- if ($prev_outlier == 1);
            $prev_outlier++;
        }
        else {$prev_outlier = 0 if ($CURR_CHR_REFS[$i]->[9] eq "AB")}

        $stats[$AB_CT]++ if ($CURR_CHR_REFS[$i]->[9] eq "AB");

        # Count nonconsecutive SNPs with LRR < HD threshold        
        if ($CURR_CHR_REFS[$i]->[5] < $HD_HI) {
            $stats[$P1_HD_CT]++ unless ($prev_p1_hd);
            $stats[$P1_HD_CT]-- if ($prev_p1_hd == 1);
            $prev_p1_hd++;
        }
        else {$prev_p1_hd = 0} 
        if ($CURR_CHR_REFS[$i]->[8] < $HD_HI) {
            $stats[$P2_HD_CT]++ unless ($prev_p2_hd);
            $stats[$P2_HD_CT]-- if ($prev_p2_hd == 1);
            $prev_p2_hd++;
        }
        else {$prev_p2_hd = 0} 
        if ($CURR_CHR_REFS[$i]->[11] < $HD_HI) {
            $stats[$CH_HD_CT]++ unless ($prev_ch_hd);
            $stats[$CH_HD_CT]-- if ($prev_ch_hd == 1);
            $prev_ch_hd++;
        }
        else {$prev_ch_hd = 0} 
        my $total = 0;
        
        if ($snp_ct < $SIZE_SMALL_WINDOW) {
            if ($CURR_CHR_REFS[$i]->[-1] == FATHER) {$p1_ct++} 
            elsif ($CURR_CHR_REFS[$i]->[-1] == MOTHER) {$p2_ct++}
            push(@bafsnp_array, $CURR_CHR_REFS[$i]);
        }
        else {  
            push(@bafsnp_array, $CURR_CHR_REFS[$i]);        
            if ($bafsnp_array[-1]->[-1] == FATHER) {$p1_ct++} 
            elsif ($bafsnp_array[-1]->[-1] == MOTHER) {$p2_ct++} 
            if ($p1_ct + $p2_ct > $stats[$SCAN_STAT]) {
                $stats[$SCAN_STAT] = $p1_ct + $p2_ct;  
            }                         
            if ($bafsnp_array[0]->[-1] == FATHER) {$p1_ct--}
            elsif ($bafsnp_array[0]->[-1] == MOTHER) {$p2_ct--} 
            shift(@bafsnp_array);
        } 
    }
    return() unless ($snp_ct >= $SIZE_SMALL_WINDOW); 
    return(\@stats);
}

sub splice_overlap_regions {
    # Splices two overlapping regions depending on type of overlap 
    my ($reg1, $reg2) = @_;
    my (@results, @overlap);
    if (!@$reg1 && !@$reg2) {return(\@results)}
    elsif (!@$reg1) {return($reg2)}
    elsif (!@$reg2) {return($reg1)}
    my ($chr, $start1, $stop1, $type1, $contrib1, $size1, $inf1, $detect1
        ) = @$reg1[0..6,17];
    my ($start2, $stop2, $type2, $contrib2, $size2, $inf2, $detect2
        ) = @$reg2[1..6,17]; 
    my ($orient, $end1_start, $end1_stop, $overlap_start, $overlap_stop,
        $end2_start, $end2_stop, $no_split) 
        = @{overlap($start1, $stop1, $start2, $stop2)};
    # Combinations of contribution and type 
    my $eq = ($contrib1 == $contrib2 && $size1 == $size2 
        && $inf1 == $inf2 && $detect1 == $detect2);
   # my $mi1_1_trumps = $detect2 == POD && $detect1 == MI1 
   #     && $$reg2[18] <= 2 && $inf1 >= 5 && (($inf2 - $inf1) < $MIN_POD);
   # my $mi1_2_trumps = $detect1 == POD && $detect2 == MI1 
    #    && $$reg1[18] <= 2 && $inf2 >= 5 && (($inf1 - $inf2) < $MIN_POD);
    my $mi1_1_trumps = $detect2 == POD && $detect1 == MI1 
        && $$reg2[18] <= 2 && (($inf2 - $inf1) < $MIN_POD);
    my $mi1_2_trumps = $detect1 == POD && $detect2 == MI1 
        && $$reg1[18] <= 2 && (($inf1 - $inf2) < $MIN_POD);
    my $combine1 = !$mi1_2_trumps && $contrib1 == $contrib2;
    my $combine2 = !$mi1_1_trumps && $contrib1 == $contrib2;            
    my $reg1_trumps = ($detect1 == POD && !$mi1_2_trumps) 
        || ($contrib1 == $contrib2 && $detect1 == $detect2 
            && ((!$inf2 || $size1 > $size2) 
            || ($size1 == $size2 && $inf1 > $inf2)))
        || $contrib1 == NONE && $contrib2 != NONE 
        || $contrib2 == UNKNOWN
        || $mi1_1_trumps;
    my $reg2_trumps = ($detect2 == POD && !$mi1_1_trumps) 
        || ($contrib1 == $contrib2 && $detect1 == $detect2
            && ((!$inf1 || $size2 > $size1) 
            || ($size1 == $size2 && $inf2 > $inf1)))
        || $contrib2 == NONE && $contrib1 != NONE 
        || $contrib1 == UNKNOWN
        || $mi1_2_trumps;    
    
    # Output combinations
    my $combo111 = [$chr, $end1_start, $end1_stop, $type1, $contrib1,
                    (0) x 11, $$reg1[16], $detect1];
    my $combo112 = [$chr, $end1_start, $end1_stop, $type2, $contrib2,
                    (0) x 11, $$reg2[16], $detect2];
    my $combo121 = [$chr, $end1_start, $end2_stop, $type1, $contrib1,
                    (0) x 11, $$reg1[16], $detect1];
    my $combo122 = [$chr, $end1_start, $end2_stop, $type2, $contrib2,
                    (0) x 11, $$reg2[16], $detect2];
    my $combo221 = [$chr, $end2_start, $end2_stop, $type1, $contrib1,
                    (0) x 11, $$reg1[16], $detect1];
    my $combo222 = [$chr, $end2_start, $end2_stop, $type2, $contrib2,
                    (0) x 11, $$reg2[16], $detect2];
    my @both_contrib;
    if ($detect2 == POD || $inf2 > $inf1){@both_contrib =($$reg2[16],$detect2)}
    else {@both_contrib = ($$reg1[16], $detect1)}
    my $comboOOB = [$chr, $overlap_start, $overlap_stop, " ", BOTH,
                    (0) x 11, @both_contrib[0,1]];

    # -----
    # -----
    if ($orient == 1) {
        $no_split = 1;
        if ($eq || $reg1_trumps) {push(@results, $reg1)}
        elsif ($reg2_trumps)     {push(@results, $reg2)}                
        else {push(@results, $comboOOB)}
    }
    #   -----
    # -----
    elsif ($orient == 2) {
        if ($eq || $combine1) {push(@results, $combo121); $no_split = 1}
        elsif ($combine2)     {push(@results, $combo122); $no_split = 1}
        elsif ($reg1_trumps)  {push(@results, $reg1, $combo112)}
        elsif ($reg2_trumps)  {push(@results, $reg2, $combo221)}
        else {push(@results, $combo112, $comboOOB, $combo221)} 
    }
    # -----
    # ---
    elsif ($orient == 3) {
        if ($eq || $reg1_trumps) {push(@results, $reg1); $no_split = 1}
        elsif ($reg2_trumps) {push(@results, $reg2, $combo221)}
        else {push(@results, $comboOOB, $combo221)}
    }
    #   ---
    # -----
    elsif ($orient == 4) {
        if ($eq || $reg2_trumps) {push(@results, $reg2); $no_split = 1}
        elsif ($reg1_trumps) {push(@results, $reg1, $combo112)}
        else {push(@results, $combo112, $comboOOB)}
    }
    #  ---
    # -----
    elsif ($orient == 5) { 
        if ($eq || $reg2_trumps) {push(@results, $reg2); $no_split = 1}
        elsif ($reg1_trumps) {push(@results, $reg1, $combo112, $combo222)}
        else {push(@results, $combo112, $comboOOB, $combo222)}
    }
    # -----
    #   -----
    elsif ($orient == 6) {
        if ($eq || $combine1) {push(@results, $combo121); $no_split = 1}
        elsif ($combine2)     {push(@results, $combo122); $no_split = 1}
        elsif ($reg1_trumps)  {push(@results, $reg1, $combo222)}
        elsif ($reg2_trumps)  {push(@results, $reg2, $combo111)}
        else {push(@results, $combo111, $comboOOB, $combo222)}
    }
    # ---
    # -----
    elsif ($orient == 7) {
        if ($eq || $reg2_trumps) {push(@results, $reg2); $no_split = 1}
        elsif ($reg1_trumps) {push(@results, $reg1, $combo222)}
        else {push(@results, $comboOOB, $combo222)}
    }
    # -----
    #   ---
    elsif ($orient == 8) {
        if ($eq || $reg1_trumps) {push(@results, $reg1); $no_split = 1}
        elsif ($reg2_trumps) {push(@results, $reg2, $combo111)}
        else {push(@results, $combo111, $comboOOB)}
    }
    # -----
    #  ---
    elsif ($orient == 9) {
        if ($eq || $reg1_trumps) {push(@results, $reg1); $no_split = 1}
        elsif ($reg2_trumps) {push(@results, $reg2, $combo111, $combo221)}
        else {push(@results, $combo111, $comboOOB, $combo221)}
    }
    else {push(@results, $reg1, $reg2); $no_split = 1}
    
    map {$$_[6] = count_SNPs(@$_[START,STOP,16])} @results;
    my @return;
    if ($no_split) {@return = @results}
    else {
        for my $ref (@results) {
            if ($$ref[17] == HD) {
                my $hd_thresh;
                if    ($$ref[4] == FATHER) {$hd_thresh = $MIN_P2_HD}
                elsif ($$ref[4] == MOTHER) {$hd_thresh = $MIN_P1_HD}
                elsif ($$ref[4] == NONE)   {$hd_thresh = $MIN_CH_HD}
                elsif ($$ref[4] == BOTH)   {
                    $hd_thresh = ($MIN_P1_HD <= $MIN_P2_HD) ? $MIN_P1_HD 
                        : $MIN_P2_HD;
                }
                push(@return, $ref) if ($$ref[6] >= $hd_thresh);
            }
            elsif ($$ref[17] == MI1) {
                push(@return, $ref) if ($$ref[6] >= $MIN_MI1);
            }
            elsif ($$ref[17] == POD || $$ref[17] == PODcr) {
                push(@return, $ref) if ($$ref[6] >= $MIN_POD);
            }
            else {push(@return, $ref)}
        }
    }
    return(\@return);
}

sub st_dev {
    # Calculates standard deviation.
    # Takes either an array ref or sum, sum_squares, count
    # Returns the mean, standard deviation, variance, and COV 
    my ($in1, $in2, $in3) = @_;
    my ($sum, $sum_squares, $ct, $mean, $variance, $coef_of_var, 
        $st_dev) = (0) x 7;   
    if (!$in1) {return(0,0,0,0)}        
    elsif ($in2) {($sum, $sum_squares, $ct) = @_}
    else {
        $sum += $_ for @$in1;
        $sum_squares += $_**2 for @$in1;
        $ct = @$in1;
    }
    if ($ct > 1) { 
        $mean = $sum / $ct;
        if ($mean) {
            $variance = ($sum_squares - ($sum**2 / $ct)) / ($ct - 1);      
            if ($variance < 0 && $variance > -0.000001) {
                $variance = $st_dev = $coef_of_var = 0; 
            }
            else {
                $st_dev = sqrt($variance); 
                $coef_of_var = $st_dev / $mean;
            }
        }
    } 
    return($mean, $st_dev, $variance, $coef_of_var);
}

sub start_R {
    my $graphics_file = "./triPOD_graphics.R";
    unless (-e $graphics_file) { 
        print STDERR "The file entitled $graphics_file cannot be ",
        "located.\nPlease place in the same directory as the Perl script.\n",
        "The analysis will continue without graphical output.\n" if $BATCH;    
        print "\n\nThe file entitled $graphics_file cannot be ",
        "located.\nPlease place in the same directory as the Perl script.\n",
        "The analysis will continue without graphical output.\n" if $verbose;
        $GRAPHICS = "";
    }
    # Makes a system call to kick off the Rscript.
    my $r_script = join(" ", "Rscript --slave --no-save --no-restore",
        "--no-environ --silent $graphics_file",
        $OUTPUT_DIR, $INPUT_FILE, $FILENAME, $P1_NAME, 
        $P2_NAME, $CH_NAME, $PID_VALUE, $R_PID_FILENAME,
        $PERL_TO_R_FILENAME, $R_TO_PERL_FILENAME, $GRAPHICS, "&> /dev/null &");
    system($r_script);
}

sub streak_pvalue {
    # Approximate the probability Q of a streak size >= m in n trials
    # Q = Q2(Q2 / Q3)^((n/m)/2)
    # where Q2 = 1 - p^m(1 = mq)
    # Q3 = 1 - p^m(1 = 2mq) + 0.5p^(2m) * (2mq + m(m-1) * q^2) 
    my($n, $m, $p) = @_;
    my $q = 1 - $p;
    my $result;
    if (!$n || $n < $m) {return(1)}
    my $Q2 = 1 - (($p**$m) * (1 + $m * $q));
    my $Q3 = 1 - (($p**$m) * (1 + 2 * $m * $q))
        + 0.5 * $p**(2*$m) * (2 * $m * $q + $m * ($m - 1) * $q**2);
    $result = 1 - ($Q2 * ($Q3 / $Q2)**(($n / $m) - 2));    
    return($result);
}

sub to_bed {
    my $region = $_[0];
    for (@$region) {
        my ($chr, $start, $end) = ("chr$$_[1]", $$_[2] - 1, $$_[3] + 1);
        my $par_contrib = substr($$_[5], 0, 1);
        $par_contrib ||= "U";
        my $inheritance;
        if ($$_[6] eq "INH") {$inheritance = "IN"}
        elsif ($$_[6] eq "INH-CN") {$inheritance = "IC"}
        else {$inheritance = "NA"}
        my $sample = substr($$_[0], 0, 11);
        my $name = "$sample($par_contrib-$inheritance)";
        my $color;    
        if    ($$_[4] eq "DEL") {$color = "255,0,0"}     #red
        elsif ($$_[4] eq "AMP") {$color = "0,0,255"}     #blue
        elsif ($$_[4] eq "HD")  {$color = "190,190,190"} #gray
        elsif ($$_[4] eq "UPhD"){$color = "160,32,240"}  #purple
        elsif ($$_[4] eq "UPiD"){$color = "0,255,0"}     #green        
        else                    {$color = "0,0,0"}       #black       
        
        print BED_FILE join("\t", $chr, $start, $end, $name, 0, "+", 
            $start, $end, $color, $name);
            
        print BED_FILE "\t<H2>Description</H2><P>The triPOD software is a ",
        "trio-based chromosomal abnormality detection program based on the ",
        "Parent of Origin-based Detection (POD) method.</P><H2>Methods</H2>",
        "<P>Chromosomal abnormalities are reported, along with parental ",
        "origin, inheritance state, and the type of abnormality.<BR> The name ",
        "field includes the sample name, the parent of origin or parental ",
        "contributor, and the inheritance state (e.g. Sample1(F-IN)).",
        "<H3>Parental Origin</H3><ul><li>Father (F)</li><li>Mother (M)</li>",
        "<li>Both (B)</li><li>Unknown (U)</li></ul><H3>Inheritance States</H3>",
        "<ul><li>Likely inherited (IN)</li><li>Likely inherited with a unique ",
        "copy number state (IC)</li><li>Undetermined inheritance patterns (NA)",
        "</li></ul><H3>The regions are colored:</H3><ul><li>",
        "<B><FONT COLOR = RED>Red</FONT></B> for deletions (normalized median ",
        "LRR >-1.5 and <-0.1)</li><li><B><FONT COLOR = BLUE>Blue</FONT></B> for ",
        "amplifications (normalized median LRR >0.1)</li><li>",
        "<B><FONT COLOR = GRAY>Gray</FONT></B> for homozygous deletions ",
        "(normalized median LRR <-1.5)</li><li><B><FONT COLOR = ",
        "GREEN>Green</FONT></B> for uniparental isodisomy </li><li><B><FONT ",
        "COLOR = PURPLE>Purple</FONT></B> for uniparental heterodisomy</li><li>",
        "<B><FONT COLOR = BLACK>Black</FONT></B> for all other options, ",
        "including low-level mosaic events (normalized median LRR ",
        ">-0.1 and <0.1)</li></ul></P>",
        "<H2>Credits</H2><P></P><H2>References</H2><P></P><H2>Contact</H2>\n";
    }
}
 
sub stirlerr {
    # stirlerr(n) = log(n!) - log( sqrt(2*pi*n)*(n/e)^n ) 
    # Called by dbinom to estimate the binomial probability
    # Translated to perl from C functions developed by 
    # Catherine Loader in "Fast and Accurate Computation of 
    # Binomial Probabilities (2000)"
    # This method is also used by R and Matlab.
    my $n = $_[0];
    my $s0 = 0.083333333333333333333;        # 1/12
    my $s1 = 0.00277777777777777777778;      # 1/360 
    my $s2 = 0.00079365079365079365079365;   # 1/1260 
    my $s3 = 0.000595238095238095238095238;  # 1/1680 
    my $s4 = 0.0008417508417508417508417508; # 1/1188
    my @sfe = (0, 0.081061466795327258219670264,
        0.041340695955409294093822081, 0.0276779256849983391487892927,
        0.020790672103765093111522771, 0.0166446911898211921631948653,
        0.013876128823070747998745727, 0.0118967099458917700950557241,
        0.010411265261972096497478567, 0.0092554621827127329177286366,
        0.008330563433362871256469318, 0.0075736754879518407949720242,
        0.006942840107209529865664152, 0.0064089941880042070684396310,
        0.005951370112758847735624416, 0.0055547335519628013710386899
        );
    if ($n < 16) {return($sfe[int($n)])}
    my $nn = $n;
    $nn = $nn*$nn;
    if ($n>500) {return(($s0-$s1/$nn)/$n)}
    elsif ($n>80) {return(($s0-($s1-$s2/$nn)/$nn)/$n)}
    elsif ($n>35) {return(($s0-($s1-($s2-$s3/$nn)/$nn)/$nn)/$n)}
    return(($s0-($s1-($s2-($s3-$s4/$nn)/$nn)/$nn)/$nn)/$n);
}
    
sub bd0 {
    # Evaluate the deviance term bd0(x,np) = x log(x/np) + np - x
    # Called by dbinom to estimate the binomial probability
    # Translated to perl from C functions developed by 
    # Catherine Loader in "Fast and Accurate Computation of 
    # Binomial Probabilities (2000)"
    # This method is also used by R and Matlab.
    my ($x, $np) = @_;
    my ($ej, $s, $s1, $v, $j);
    if (abs($x-$np)<0.1*($x+$np))
    {    $s = ($x-$np)*($x-$np)/($x+$np);
        $v = ($x-$np)/($x+$np);
        $ej = 2*$x*$v;
        for ($j=1; ;$j++)
        {    $ej *= $v*$v;
            $s1 = $s+$ej/(2*$j+1);
            if ($s1==$s) {return($s1)}
            $s = $s1;
        }    
    } 
    return($x*log($x/$np)+$np-$x);
}

sub dbinom {
    # Estimates the binomial probability
    # Translated to perl from C functions developed by 
    # Catherine Loader in "Fast and Accurate Computation of 
    # Binomial Probabilities (2000)"
    # This method is also used by R and Matlab.
    my ($x, $n, $p) = @_;
    my $PI2 = 6.283185307179586476925286;
    my $lc;
    if ($p==0.0) {return( ($x==0) ? 1.0 : 0.0)}
    if ($p==1.0) {return( ($x==$n) ? 1.0 : 0.0)}
    if ($x==0) {return(exp($n*log(1-$p)))}
    if ($x==$n) {return(exp($n*log($p)))}
    $lc = stirlerr($n) - stirlerr($x) - stirlerr($n-$x)
        - bd0($x,$n*$p) - bd0($n-$x,$n*(1.0-$p));
    return(exp($lc)*sqrt($n/($PI2*$x*($n-$x))));
}
