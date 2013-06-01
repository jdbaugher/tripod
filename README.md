triPOD
======

This software program detects chromosomal abnormalities in parent-child
trio-based microarray data. It has been shown to provide improved 
sensitivity for abnormality detection, with substantial improvement in 
detection of low-level mosaic abnormalities, as described by 
[Joseph D. Baugher, Benjamin D. Baugher, Matthew D. Shirley, and Jonathan Pevsner. Sensitive and specific detection of mosaic chromosomal abnormalities using the Parent-of-Origin-based Detection (POD) method. BMC Genomics 2013.](http://www.biomedcentral.com/1471-2164/14/367/)

## Requirements
triPOD is designed for use on Unix-based operating systems.  It was
developed and tested on Red Hat Enterprise Linux 6 using Perl v5.10.1
and R version 2.14.1.

#### Perl modules:
    Algorithm::Cluster
    threads
    Tree::Interval
#### R libraries:
    shape
    TTR

## Usage
    perl triPOD.pl [--options] [INPUT_FILE]

    Optional Arguments:
      --alpha   The threshold under which the familywise error rate is 
                controlled. Default = --alpha=0.1
      --batch   Submit a file containing a list of file names to be run 
                sequentially.
      --build   The path to a file containing centromere locations from
                the appropriate UCSC assembly.
                Default = --build=./genome_build/hg18_centromeres.txt
      --cite    Prints reference info for citations
      --cores   Number of CPU cores to employ 
                Default = maximum cores - 1 (e.g. --cores=8)
      --gender  Gender designation for sample (M or F). 
                Currently, chromosome X is analyzed only if a female 
                gender is specified. Default = NA
      --graph   Creates graphic output of results in PNG or PDF format 
                (--graph=none, --graph=png, --graph=pdf, or --graph=both)
                Default = --graph=none
      --hd      Abnormality detection by homozygous deletion analysis (PODhd)
                (--hd or --nohd). Default = --hd
      --help    Prints a message describing optional arguments and formatting
      --hetSD   Heterozygous SNP threshold as standard deviations from mean BAF
                Default = --hetSD=1.414213562373095 (sqrt of 2)
      --homSD   Homozygous SNP threshold as standard deviations from mean BAF 
                Default = --homSD=4
      --mi1     Abnormality detection by Mendelian "error" analysis (PODmi1) 
                (--mi1 or --nomi1). Default = --mi1
      --nc      Samples must have a No Call rate below this threshold to be 
                analyzed. Default = --nc=0.03
      --out     Specify an output directory (e.g. --out=results) 
                Default = --out=triPOD_Results
      --pod     Abnormality detection by standard POD algorithm 
                (--pod or --nopod). Default = --pod
      --podcr   Abnormality detection by cryptic POD algorithm (PODcr) 
                (--podcr or --nopodcr). Default = --podcr
      --stats   Create a file including calculated parameters, etc. (--stats) 
                Not created by default.
      --verbose Prints progress info to screen. Negated using --noverbose. 
                Default = --verbose. 
                Note that --batch mode will run in --noverbose mode.
      --win     Number of SNPs per window for sliding window analysis
                Default = --win=100

**Note** 

  - It is highly recommended that parameters such as --alpha, --hetSD, --homSD,
  and --win are altered only by experienced users.

## Input Formatting

  The input file must be **tab delimited**, **sorted** by *chromosome* and 
  *position*, and in the following order (columns): 
    
    SNP_Name Chromosome Position Father.GType Father.BAF Father.LRR
     Mother.GType Mother.BAF Mother.LRR Child.GType Child.BAF Child.LRR

  A header line is expected and is used to extract sample names 
  (i.e. Sample1.GType = Sample1), but not to determine column identity.

**Genotype Annotation**

  The genotypes must be AA, AB, BB, NC (or NoCall).  
	triPOD has been developed using genotypes annotated by the Illumina method 
  described in ["TOP/BOT" Strand and "A/B" Allele.](www.illumina.com/documents/products/technotes/technote_topbot.pdf) 
  If converting from HapMap format (ATCG) to Illumina format (AB), 
  use simple replacement as follows: 
    
    AA = AA, TT; BB = CC, GG; AB = AC, AG, TC, TG; -- = NC. 
  
  Any markers with alternative genotyping combinations (e.g. CG or TA) should 
  be discarded unless the user performs an additional analysis of the 
  surrounding sequence.

  B allele frequencies must be >= 0 and <= 1 for polymorphic markers. 

## Output 

Several types of output are produced by triPOD.
	
  - *triPOD_Results.txt
  
  A tabular listing of detected abnormalities including the following information:
	  Sample Name, Chromosome, Abnormal region start position, Abnormal
	  region end position, Type of Abnormality, Parental Origin,
	  Inheritance Pattern, Size of abnormal region (SNPs), Number of
	  Informative SNPs, Size of abnormal region (bps), Detection
	  Method, 
    Descriptive statistics for the abnormal region: 
	  Child's normalized median mirrored BAF (mBAF) and LRR, 
	  Father's normalized median mirrored BAF (mBAF) and LRR,
	  Mother's normalized median mirrored BAF (mBAF) and LRR.
	
	The results file also provides the date and time, triPOD version, 
  the command line information, and the parameters employed.
	
  - *triPOD_Results.bed
  
	A .bed file is provided for visualization with genome browsers.  
	
  - *triPOD_stats.txt
  
	An optional file which contains sample-specific calculations and
	parameters.  See [Baugher JD, et al.](http://www.biomedcentral.com/1471-2164/14/367/) 
  for explanations of rate and probability calculations. The 
  following information is included for each trio: 
  file name, trio member names, NoCall rates, 
  nonadjacent homozygous deletion (HD) and single
	Mendelian error (MI1) rates, minimum region sizes, calculated
	MI1 BAF thresholds, initial and corrected alpha values for FDR
	control, estimated error rate, acceptable errors for POD
	region extension, global BAF thresholds for detecting
	informative SNPs, counts of regions detected by each detection
	method.  
	
  - ./log_files/*triPOD_log.txt 
	
  STDERR is redirected to a log file.  
	
  - *.png and/or *.pdf
  
	Optional graphical output is generated for each chromosome 
	harboring an abnormality. The upper panel is a plot of the log 
	R ratio (LRR) values with the LRR moving average plotted in green
	and a black horizontal line at y=0 to assist with visualization of
	copy number changes. The center panel is a plot of the B allele 
	frequency (BAF), in which abnormalities can frequently be 
	visualized as splitting of the heterozygous (center) BAF band. The 
	lower panel is a plot of the abnormalities detected by triPOD. 
	Abnormalities with known parental origin or contribution are 
	plotted adjacent to the Father and Mother axis labels. 
	Abnormalities with unknown contribution or abnormal contribution 
	from both parents are plotted along the y=0.5 axis, between the 
	Father and Mother labels. The types of abnormalities are color 
	coded as follows: red = deletion, blue = amplification, gray = 
	homozygous deletion, purple = uniparental heterodisomy, green = 
	uniparental isodisomy, and black = all other abnormalities for 
	which a type cannot be assigned (including low-level mosaic 
	events of all types). 

### Additional Considerations

  - Care should be taken when interpreting large abnormalities in 
  commonly variable regions due to a tendency to combine small 
  adjacent abnormalities of the same type.
  
  - Interesting findings should be graphically investigated until 
  the user has gained expertise with the strengths and weaknesses 
  of the various detection and annotation methods employed by triPOD. 

## Author

Joseph D. Baugher, Ph.D.  -  jbaughe2(at)jhmi.edu
