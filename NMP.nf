#!/usr/bin/env nextflow
	
/**
	Nsilico Metagenomic Pipeline (NMP)
	Copyright (C) 2018 	Dr Bruno Andrade 	      
	
*/

version='1.0'
timestamp='20181903'

/**
	Prints version when asked for
*/
if (params.version) {
	System.out.println("")
	System.out.println("Nsilico metagenomic pipeline (NMP) - Version: $version ($timestamp)")
	exit 1
}

/**
	Prints help when asked for
*/

if (params.help) {
	System.out.println("")
	System.out.println("Nsilico metagenomic pipeline (NMP) - Version: $version ($timestamp)")
	System.out.println("")
	System.out.println("Usage")
	System.out.println("   nextflow run -c nextflow.config NMP.nf --reads_R1 R1 --reads_R2 R2 --date date --prefix mysample --outdir path")
	System.out.println("                [options] [-with-docker]")
	System.out.println("")
	System.out.println("Mandatory arguments:")
	System.out.println("    --reads_R1   R1      Path for forward reads (library = paired-end) or for all reads (library = single-end")
	System.out.println("    [--reads_R2] R2      Path for reverse reads (library = paired-end)")
	System.out.println("    --date date      Date of collection of an individual microbiome sample")
	System.out.println("    --prefix   prefix  Prefix used to name the result files")
	System.out.println("    --outdir   path    Output directory (will be outdir/prefix/date)")
	System.out.println("Options:")
	System.out.println("    --library <single-end|paired-end>")
	System.out.println("")
	System.out.println("Container:")
	System.out.println("    Docker image to use with docker is")
	System.out.println("    'docker://bgabriel/nmp'")
    exit 1
}

	
/**

	Checking user-defined parameters

*/


//--reads_R2 can be omitted if the library layout is "single-end"

if (params.library != "single-end" && (params.reads_R2 == "null") ) {
	exit 1, "If dealing with paired-end reads, please set the reads_R2 arguments and if dealing with single-end reads, please set the library argument to 'single-end'"
}


//Creates working dir
workingpath = params.outdir + "/" + params.prefix
workingdir = file(workingpath)
if( !workingdir.exists() ) {
    if( !workingdir.mkdirs() ) 	{
        exit 1, "Cannot create working directory: $workingpath"
    } 
}

//Creating other folders
matrixdir = file(params.outdir + "/" + params.prefix + "/Matrix")
plotdir = file(params.outdir + "/" + params.prefix + "/Plots")
QCdir = file(params.outdir + "/" + params.prefix + "/QC")

if( !matrixdir.exists() ) {
    if( !matrixdir.mkdirs() ) 	{
        exit 1, "Cannot create Matrix directory"
    } 
}

if( !plotdir.exists() ) {
    if( !plotdir.mkdirs() ) {
	exit 1, "Cannot create directory for plots"
    }
}

if( !QCdir.exists() ) {
    if( !QCdir.mkdirs() ) {
        exit 1, "Cannot create directory for Quality check files"
    }
}



/**
	Quality Assessment 

*/


if (params.library == "paired-end") {
	rawreads = Channel.from( ['1', file(params.reads_R1), '_R1', '_rawreads'], ['1', file(params.reads_R2), '_R2', '_rawreads'] )
}
else {
	rawreads = Channel.value( ['1', file(params.reads_R1), '', '_rawreads'] )
}


// ------------------------------------------------------------------------------   
// 				Quality control (QC)				
// ------------------------------------------------------------------------------   

/**
	Quality control - STEP 1. Trimming of low quality bases and of adapter sequences. Short reads
	are discarded. A decontamination of synthetic sequences is also peformed.
	If dealing with paired-end reads, when either forward or reverse of a paired-read
	are discarded, the surviving read is saved on a file of singleton reads.

*/




if (params.library == "paired-end") {
	Trimming = Channel.value( [file(params.reads_R1), file(params.reads_R2)] )
} 
else {
	Trimming = Channel.value( [file(params.reads_R1), "null"] )
}


//When single-end reads are used, the input tuple (singleton) will not match input set 
//cardinality declared by 'trim' process (pair), so two mock files were pushed in the channel,
//and then the first two files were taked
mocktrim = Channel.from("null")
process trim {

	input:
   	set file(reads1), file(reads2) from Trimming.concat(mocktrim).flatMap().take(2).buffer(size : 2)
	file(adapters) from Channel.from( file(params.adapters) )
	file(artifacts) from Channel.from( file(params.artifacts) )
	file(phix174ill) from Channel.from( file(params.phix174ill) )

	output:

	file("${params.prefix}*.fq") into todecontaminate
	file("${params.prefix}*.fq") into trimmedreads


   	script:
	"""	

	#Sets the maximum memory to the value requested in the config file
	maxmem=\$(echo ${task.memory} | sed 's/ //g' | sed 's/B//g')

	#Defines command for trimming of adapters and low quality bases
	if [ \"$params.library\" = \"paired-end\" ]; then
		CMD=\"bbduk.sh -Xmx\"\$maxmem\" in=$reads1 in2=$reads2 out=${params.prefix}_trimmed_R1_tmp.fq out2=${params.prefix}_trimmed_R2_tmp.fq outs=${params.prefix}_trimmed_singletons_tmp.fq ktrim=r k=$params.kcontaminants mink=$params.mink hdist=$params.hdist qtrim=rl trimq=$params.Quality  minlength=$params.minlength ref=$adapters qin=$params.Pcoding threads=${task.cpus} tbo tpe ow\"
	else
		CMD=\"bbduk.sh -Xmx\"\$maxmem\" in=$reads1 out=${params.prefix}_trimmed_tmp.fq ktrim=r k=$params.kcontaminants mink=$params.mink hdist=$params.hdist qtrim=rl trimq=$params.Quality  minlength=$params.minlength ref=$adapters qin=$params.Pcoding threads=${task.cpus} tbo tpe ow\"
	fi

	#Trims adapters and low quality bases	
	exec \$CMD 2>&1
	
	#Removing synthetic contaminants
	if [ \"$params.library\" = \"paired-end\" ]; then
		bbduk.sh -Xmx\"\$maxmem\" in=${params.prefix}_trimmed_R1_tmp.fq in2=${params.prefix}_trimmed_R2_tmp.fq out=${params.prefix}_trimmed_R1.fq out2=${params.prefix}_trimmed_R2.fq k=31 ref=$phix174ill,$artifacts qin=$params.Pcoding threads=${task.cpus} ow
	else
		bbduk.sh -Xmx\"\$maxmem\" in=${params.prefix}_trimmed_tmp.fq out=${params.prefix}_trimmed.fq k=31 ref=$phix174ill,$artifacts qin=$params.Pcoding threads=${task.cpus} ow
	fi


	if [ \"$params.library\" = \"paired-end\" ]; then
	bbduk.sh -Xmx\"\$maxmem\" in=${params.prefix}_trimmed_singletons_tmp.fq out=${params.prefix}_trimmed_singletons.fq k=31 ref=$phix174ill,$artifacts qin=$params.Pcoding threads=${task.cpus} ow
	
	fi
	#Removes tmp files. This avoids adding them to the output channels
	rm -rf ${params.prefix}_trimmed*_tmp.fq 
	"""
}



/**
	Decontamination step.
*/


mockdecontaminate = Channel.from("null", "null")

process decontaminate {

	publishDir  QCdir, mode: 'copy', pattern: "*_clean.fq"
		
	input:
	set file(infile1), file(infile2), file(infile12) from todecontaminate.concat(mockdecontaminate).flatMap().take(3).buffer(size : 3)
	file(RefGenome) from Channel.from( file(params.RefGenome, type: 'dir') )
	
	output:

	file "${params.prefix}*_cont.fq" into decontaminatedreads
	file "${params.prefix}*_clean.fq" into toprofiletaxa

	script:
	"""

	#Sets the maximum memory to the value requested in the config file
	maxmem=\$(echo ${task.memory} | sed 's/ //g' | sed 's/B//g')
	
	if [ \"$params.library\" = \"paired-end\" ]; then
		CMD=\"bbwrap.sh  -Xmx\"\$maxmem\" mapper=bbmap append=t in1=$infile1,$infile12 in2=$infile2,null outu=${params.prefix}_clean.fq outm=${params.prefix}_cont.fq minid=$params.mind maxindel=$params.maxindel bwr=$params.bwr bw=12 minhits=2 qtrim=rl trimq=$params.Quality path=$RefGenome qin=$params.Pcoding threads=${task.cpus} untrim quickmatch fast ow\"
	else
		CMD=\"bbwrap.sh  -Xmx\"\$maxmem\" mapper=bbmap append=t in1=$infile1 outu=${params.prefix}_clean.fq outm=${params.prefix}_cont.fq minid=$params.mind maxindel=$params.maxindel bwr=$params.bwr bw=12 minhits=2 qtrim=rl trimq=$params.Quality path=$RefGenome qin=$params.Pcoding threads=${task.cpus} untrim quickmatch fast ow\"
	fi

	#Execute decontamination
	exec \$CMD 2>&1
	
	"""
}

// ------------------------------------------------------------------------------
//                             Microbiome profiling                            
// ------------------------------------------------------------------------------

/**
        Microbiome profiling.

        Metaphlan2 will be used to anotate shotgun reads in taxonomic levels, this step will generate 3 outputs, a matrix containing the taxonomic levels, relative abundaces, read counts for each clade, sample size and average genome size in the working file and two others containing only the relative abundance and the read_counts in Matrix/Abundance and Matrix/Read_counts, respectively.
*/


process profileTaxa {

	publishDir  matrixdir, mode: 'copy', pattern: "*.txt"
	
	input:
	file(infile) from toprofiletaxa
	file(mpa_pkl) from Channel.from( file(params.mpa_pkl) )


    output:

	file "${params.date}*.txt" into toprofilefunctionbugs


	script:
	"""
	#If a file with the same name is already present, Metaphlan2 will crash
	rm -rf ${params.prefix}_bt2out.txt
	
	#Defines command for estimating abundances

	CMD=\"python /opt/biobakery-metaphlan2-*/metaphlan2.py --input_type fastq --bowtie2out=${params.prefix}.bt2 -t rel_ab --bt2_ps $params.bt2options --nproc ${task.cpus} $infile -o ${params.date}.txt \"



	#Estimates microbial abundances
	exec \$CMD 2>&1


	"""
}

// ------------------------------------------------------------------------------
//                     Quality assessment and visualization                    
// ------------------------------------------------------------------------------



if (params.library == "paired-end") {
	trimmedreads2qc = Channel.from('4').combine(trimmedreads.flatMap().merge( Channel.from( ['_R1', '_R2'] ) ){ a, b -> [a, b] }).combine(Channel.from('_trimmedreads'))
} 
else {
	trimmedreads2qc = Channel.from('4').combine(trimmedreads.flatMap()).combine( Channel.from( '' ) ).combine(Channel.from('_trimmedreads'))
}
decontaminatedreads2qc = Channel.from('6').combine(decontaminatedreads).combine( Channel.from( '' ) ).combine(Channel.from('_decontaminatedreads'))

//Creates the channel which performs the QC
toQC = rawreads.mix(trimmedreads2qc, decontaminatedreads2qc) 

//Process performing all the Quality Assessment
process qualityAssessment {
	
	publishDir  QCdir, mode: 'copy', pattern: "*.{html,txt}"
	  	
	input:
   	set val(step), file(reads), val(label), val(stem) from toQC

	output:
	file "*.html" 

   	script:
	"""	
	
	fastqc --quiet --noextract --format fastq --outdir=. --threads ${task.cpus} $reads
	

	"""	
}
