### Overview
Mendel Gene Dropping is a component of the umbrella [OpenMendel](https://openmendel.github.io) project. This analysis option “drops” genes into a pedigree through the founders and propagates them through the pedigrees. In many circumstances, it is useful to simulate genetic data consistent with a postulated map. For example, you might want to generate p-values empirically or to estimate the power of a collection of pedigrees to detect linkage. Gene dropping randomly fills in genotypes subject to prescribed allele frequencies, a given genetic map, and Hardy-Weinberg and linkage equilibrium. After performing gene dropping the package generates a new pedigree file for further analysis. The missing data pattern in the new pedigrees can mimic the input pedigrees, or the user can specify a new pattern.

### Appropriate Problems and Data Sets
The raw material for gene dropping consists of sets of pedigrees and loci. People within the pedigrees must be assigned either blank phenotypes or Mendelian consistent phenotypes. Gene dropping is carried out independently of observed phenotypes at those loci common to the definition and map files. By varying the content of the map file, you can choose exactly which loci to subject to gene dropping. Phenotypes at the remaining loci of the definition file are left untouched. Simulated genotypes rather than simulated phenotypes are reported. There are no limits on the complexity of the pedigrees or the number of loci. You can use founders from different populations, provided these populations are defined. Each founder should be assigned to a population; any unassigned founders are assumed to come from the first population.

### Installation
*Note: The three OpenMendel packages (1) [SnpArrays](https://openmendel.github.io/SnpArrays.jl/latest/), (2) [MendelSearch](https://openmendel.github.io/MendelSearch.jl), and (3) [MendelBase](https://openmendel.github.io/MendelBase.jl) must be installed before any other OpenMendel package will run. It is easiest if these three packages are installed in the above order and before any other OpenMendel package.*

Within Julia, use the package manager to install Mendel		:

    pkg> add https://github.com/OpenMendel/MendelGeneDropping.jl.git

This package supports Julia v1.0+

### Input Files
The MendelGeneDropping analysis package uses the following input files. Example input files can be found in the [data](https://github.com/OpenMendel/MendelGeneDropping.jl/tree/master/data) subfolder of the MendelGeneDropping project. (An analysis won't always need every file type below.)

* [Control File](#control-file): Specifies the names of your data input and output files and any optional parameters (*keywords*) for the analysis. (For a list of common keywords, see [Keywords Table](https://openmendel.github.io/MendelBase.jl/#keywords-table)).
* [Locus File](https://openmendel.github.io/MendelBase.jl/#locus-file): Names and describes the genetic loci in your data.
* [Pedigree File](https://openmendel.github.io/MendelBase.jl/#pedigree-file): Gives information about your individuals, such as name, sex, family structure, and ancestry.
* [Phenotype File](https://openmendel.github.io/MendelBase.jl/#phenotype-file): Lists the available phenotypes.
* [SNP Definition File](https://openmendel.github.io/MendelBase.jl/#snp-definition-file): Defines your SNPs with information such as SNP name, chromosome, position, allele names, allele frequencies.
* [SNP Data File](https://openmendel.github.io/MendelBase.jl/#snp-data-file): Holds the genotypes for your data set. Must be a standard binary PLINK BED file in SNP major format. If you have a SNP data file you must have a SNP definition file.

<a id="control-file"></a>
### Control file
The Control file is a text file consisting of keywords and their assigned values. The format of the Control file is:

	Keyword = Keyword_Value(s)

Below is an example of a simple Control file to run Gene Dropping:

	#
	# Input and Output files.
	#
	locus_file = genedropping LocusFrame.txt
	pedigree_file = genedropping PedigreeFrame.txt
	phenotype_file = genedropping PhenotypeFrame.txt
	new_pedigree_file = genedropping NewPedigreeFrame.txt
	output_file = genedropping Output.txt
	#
	# Analysis parameters for GeneDropping option.
	#
	repetitions = 2
	populations = European, African, Chinese

In the example above, there are seven keywords. Three keywords specify the input files: *genedropping LocusFrame.txt*, *genedropping PedigreeFrame.txt*, and *genedropping PhenotypeFrame.txt*. Two keywords specify the output files: *genedropping Output.txt* is the results file and *genedropping NewPedigreeFrame.txt* is the new pedigree OpenMendel generates, with the added simulated genotypes. The last two keywords specify analysis parameters: *repetitions* (2), and *populations* (European, African, and Chinese). The text after the '=' are the keyword values.

<a id="keywords-table"></a>
### Keywords
This is a list of OpenMendel keywords specific to Gene Dropping. A list of OpenMendel keywords common to most analysis package can be found [here](https://openmendel.github.io/MendelBase.jl/#keywords-table). The names of keywords are *not* case sensitive. (The keyword values *may* be case sensitive.)

Keyword          |   Default Value    | Allowed Values |  Short Description       
----------------      |  ----------------       |  ----------------      |  ----------------
gene_drop_output  | Unordered | Unordered, Ordered, Sourced, or Population |   Output format style for gene drop pedigrees 
interleaved          | FALSE |  TRUE, FALSE  |  Whether or not genes are dropped into pedigrees sequentially or interleaved
keep_founder_genotypes           | TRUE  |  TRUE, FALSE  |  Whether or not founder genotypes are retained
missing_rate   | 0.0 |   Real     |       
repetitions    |   1   |   integer     |       Repetitions for sharing statistics

### Data Files
Gene Dropping requires a [Control file](https://openmendel.github.io/MendelBase.jl/#control-file), and a [Pedigree file](https://openmendel.github.io/MendelBase.jl/#pedigree-file). Genotype data can be included in the Pedigree file, in which case a [Locus file](https://openmendel.github.io/MendelBase.jl/#locus-file) is required. Alternatively, genotype data can be provided in a [SNP data file](https://openmendel.github.io/MendelBase.jl/#snp-data-file), in which case a [SNP Definition File](https://openmendel.github.io/MendelBase.jl/#snp-definition-file) is required. OpenMendel will also accept [PLINK format](http://zzz.bwh.harvard.edu/plink) FAM and BIM files. Details on the format and contents of the Control and data files can be found on the [MendelBase](https://openmendel.github.io/MendelBase.jl) documentation page. There are example data files in the Gene Dropping [data](https://github.com/OpenMendel/MendelGeneDropping.jl/tree/master/data) folder.

### Running the Analysis

To run this analysis package, first launch Julia. Then load the package with the command:

     julia> using MendelGeneDropping

Next, if necessary, change to the directory containing your files, for example,

     julia> cd("~/path/to/data/files/")

Finally, to run the analysis using the parameters in your Control file, for example, Control_file.txt, use the command:

     julia> GeneDropping("Control_file.txt")

*Note: The package is called* MendelGeneDropping *but the analysis function is called simply* GeneDropping.

<!--- ### Interpreting the results
 ... --->

### Citation

If you use this analysis package in your research, please cite the following reference in the resulting publications:

*Lange K, Papp JC, Sinsheimer JS, Sripracha R, Zhou H, Sobel EM (2013) Mendel: The Swiss army knife of genetic analysis programs. Bioinformatics 29:1568-1570.*

<!--- ### Contributing
We welcome contributions to this Open Source project. To contribute, follow this procedure ... --->

### Acknowledgments

This project is supported by the National Institutes of Health under NIGMS awards R01GM053275 and R25GM103774 and NHGRI award R01HG006139.
