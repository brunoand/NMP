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
	System.out.println("   nextflow run -c nextflow.config NMP.nf --reads_R1 R1 --reads_R2 R2 --date date --prefix mysample --outdir path --mode MODE  ")
	System.out.println("                [options] [-with-docker]")
	System.out.println("")
	System.out.println("Mandatory arguments:")
	System.out.println("    --reads_R1   R1      Path for forward reads (library = paired-end) or for all reads (library = single-end")
	System.out.println("    [--reads_R2] R2      Path for reverse reads (library = paired-end)")
	System.out.println("    --date date      Date of collection of an individual microbiome sample")
	System.out.println("    --prefix   prefix  Prefix used to name the result files")
	System.out.println("    --outdir   path    Output directory (will be outdir/prefix/date)")
	System.out.println("    --mode     <QC|complete>")
	System.out.println("Options:")
	System.out.println("    --library <single-end|paired-end>")
	System.out.println("    --dedup         <true|false>   whether to perform de-duplication")
	System.out.println("")
	System.out.println("Container:")
	System.out.println("    Docker image to use with docker is")
	System.out.println("    'docker://bgabriel/nmp'")
    exit 1
}

	
/**
	STEP 0. 
	
	Checks input parameters and (if it does not exists) creates the directory 
	where the results will be stored. 
	Initialises the log file.
	
	The working directory is named after the prefix and the correspondent sample date
        located in the outdir folder. The log file, that will save summary statistics,
        execution time,	and warnings generated during the pipeline execution, will be 
        saved in the working directory as "prefix.log".
*/


//Checking user-defined parameters	
if (params.mode != "QC" && params.mode != "complete") {
	exit 1, "Mode not available. Choose any of <QC, characterisation, complete>"
}	

if (params.library != "paired-end" && params.library != "single-end") { 
	exit 1, "Library layout not available. Choose any of <single-end, paired-end>" 
}   

if (params.Pcoding != 33 && params.Pcoding != 64) { 
	exit 1, "Input quality offset (Pcoding) not available. Choose either 33 (ASCII+33) or 64 (ASCII+64)" 
}   

//--reads_R2 can be omitted when the library layout is "single-end"

if (params.library != "single-end" && (params.reads_R2 == "null") ) {
	exit 1, "If dealing with paired-end reads, please set the reads_R2 arguments and if dealing with single-end reads, please set the library argument to 'single-end'"
}



if (params.mode != "characterisation" && ( (params.library == "paired-end" && (params.reads_R1 == "null" || params.reads_R2 == "null")) ||  params.library == "single-end" && params.reads_R1 == "null") ) {
	exit 1, "Please set the reads_R1 and/or reads_R2 parameters"
}


//Creates working dir
workingpath = params.outdir + "/" + params.prefix + "/" + params.date
workingdir = file(workingpath)
if( !workingdir.exists() ) {
    if( !workingdir.mkdirs() ) 	{
        exit 1, "Cannot create working directory: $workingpath"
    } 
}

//Creating other folders
matrixpath = params.outdir + "/" + params.prefix + "/Matrix"
abundancepath = params.outdir + "/" + params.prefix + "/Matrix/Abundance"
read_count_path = params.outdir + "/" + params.prefix + "/Matrix/Read_count"
plotpath = params.outdir + "/" + params.prefix + "/Plots"
matrixdir = file(matrixpath)
abundancedir = file(abundancepath)
read_count_dir = file(read_count_path)
plotdir = file(plotpath)


if( !matrixdir.exists() ) {
    if( !matrixdir.mkdirs() ) 	{
        exit 1, "Cannot create working directory: $matrixpath"
    } 
}
if( !abundancedir.exists() ) {
    if( !abundancedir.mkdirs() ) 	{
        exit 1, "Cannot create working directory: $abundancepath"
    } 
}
if( !read_count_dir.exists() ) {
    if( !read_count_dir.mkdirs() ) 	{
        exit 1, "Cannot create working directory: $read_count_path"
    } 
}
if( !plotdir.exists() ) {
    if( !plotdir.mkdirs() ) {
	exit 1, "Cannot creat working directory: $plotpath"
    }
}


//Creates main log file
mylog = file(params.outdir + "/" + params.prefix + "/" + params.date + "/" + params.prefix + ".log")

//Logs headers
mylog <<  """---------------------------------------------
|Nsilico Metagenomic Pipeline (NMP) - Version: $version ($timestamp)|
----------------------------------------------

"""
	   
//Fetches information on OS and java versions, including user name
osname = System.getProperty("os.name") //Operating system name
osarch = System.getProperty("os.arch") //Operating system architecture
osversion = System.getProperty("os.version") //Operating system version
osuser = System.getProperty("user.name") //User's account name

javaversion = System.getProperty("java.version") //Java Runtime Environment version
javaVMname = System.getProperty("java.vm.name") //Java Virtual Machine implementation name
javaVMVersion = System.getProperty("java.vm.version") //Java Virtual Machine implementation version

//Gets starting time		
sysdate = new java.util.Date() 
		
//Logs starting time and other information about the run		
mylog << """ 
Analysis starting at $sysdate by user: $osuser
Analysed sample(s): $params.reads_R1 and $params.reads_R2
Results will be saved at $workingdir
New files will be saved using the '$params.prefix' prefix

Analysis mode? $params.mode
Library layout? $params.library
Performing de-duplication? $params.dedup	

---------------------------------------------------------------------------------------
									    
	Analysis introspection:            				    
									    
	Operating System:						    
		name:         $osname					    
		architecture: $osarch					    
		version:      $osversion				    
									    
	Java:								    
		version: $javaversion					    
		Java Virtual Machine: $javaVMname ; version: $javaVMVersion 
									    
	nextflow:							    
		version:   $nextflow.version				    
		build:     $nextflow.build				    	
		timestamp: $nextflow.timestamp				    	
									    
	Container:							    
		Docker image: $workflow.container			    
									    
	Analysis environment:						    
									    
		projectDir: $workflow.projectDir			    
		launchDir:  $workflow.launchDir				    
		workingDir: $workflow.workDir				    
									    
		command line: $workflow.commandLine			    
									    
		Run name:   $workflow.runName				    
		Session ID: $workflow.sessionId				    
		profile:    $workflow.profile				    
									    
---------------------------------------------------------------------------------------   
	   
""" 

/**
	Quality Assessment - STEP 1. Read quality assessment using the FastQC software. 
	Multiple  plots are generated to show average phred quality scores and other metrics.
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
	Quality Control - STEP 1. De-duplication. Only exact duplicates are removed.

	This step is OPTIONAL. De-duplication should be carried on iff you are
        using PCR amplification (in this case identical reads are technical artefacts)
	but not otherwise (identical reads will identify natural duplicates).
*/

// Defines channel with <readfile1, readfile2> as input for de-duplicates
// When layout is single, the params.reads_R2 is not used.

if (params.library == "paired-end") {
	deduplicate = Channel.value( [file(params.reads_R1), file(params.reads_R2)] )
} 
else {
	deduplicate = Channel.value( [file(params.reads_R1), "null"] )
}

process Deduplicate {
	
	input:
	set file(in1), file(in2) from deduplicate

	output:

	file("${params.prefix}_dedupe*.fq.gz") into Trimming

	when:
	(params.mode == "QC" || params.mode == "complete") && params.dedup

	script:
	"""
	#Measures execution time
	sysdate=\$(date)
	starttime=\$(date +%s.%N)
	echo \"Performing Quality Control. STEP 1 [De-duplication] at \$sysdate\" >>  $mylog
	echo \" \" >>  $mylog
	
	#Sets the maximum memory to the value requested in the config file
	maxmem=\$(echo \"$task.memory\" | sed 's/ //g' | sed 's/B//g')
	
	#Defines command for de-duplication
	if [ \"$params.library\" = \"paired\" ]; then
		CMD=\"clumpify.sh -Xmx\"\$maxmem\" in1=$in1 in2=$in2 out1=${params.prefix}_dedupe_R1.fq.gz out2=${params.prefix}_dedupe_R2.fq.gz qin=$params.Pcoding dedupe subs=0 threads=${task.cpus}\"
	else
		CMD=\"clumpify.sh -Xmx\"\$maxmem\" in=$in1 out=${params.prefix}_dedupe.fq.gz qin=$params.Pcoding dedupe subs=0 threads=${task.cpus}\"
	fi

	#De-duplicates
	exec \$CMD 2>&1 | tee tmp.log

	#Logs some figures about sequences passing de-duplication
	echo  \"Clumpify's de-duplication stats: \" >> $mylog
	echo \" \" >>  $mylog
	sed -n '/Reads In:/,/Duplicates Found:/p' tmp.log >>  $mylog
	echo \" \" >>  $mylog
	totR=\$(grep \"Reads In:\" tmp.log | cut -f 1 | cut -d: -f 2 | sed 's/ //g')
	remR=\$(grep \"Duplicates Found:\" tmp.log | cut -f 1 | cut -d: -f 2 | sed 's/ //g')
	survivedR=\$((\$totR-\$remR))
	percentage=\$(echo \$survivedR \$totR | awk '{print \$1/\$2*100}' )
	echo \"\$survivedR out of \$totR paired reads survived de-duplication (\$percentage%, \$remR reads removed)\" >>  $mylog
	echo \" \" >>  $mylog
		
	#Logs version of the software and executed command (BBmap prints on stderr)
	version=\$(clumpify.sh --version 2>&1 >/dev/null | grep \"BBMap version\") 
	echo \"Using clumpify.sh in \$version \" >>  $mylog
	echo \" \" >>  $mylog

	"""

}



/**
	Quality control - STEP 2. Trimming of low quality bases and of adapter sequences. Short reads
a	are discarded. A decontamination of synthetic sequences is also peformed.
	If dealing with paired-end reads, when either forward or reverse of a paired-read
	are discarded, the surviving read is saved on a file of singleton reads.

	If layout is "paired-end", three compressed FASTQ file (forward/reverse paired-end 
	and singleton reads) are outputed, if "single-end", only one compressed FASTQ file
	is returned. Several information are logged.
*/


//When the de-suplication is not done, the raw file should be pushed in the corret channel
if (!params.dedup) {
	if (params.library == "paired-end") {
		Trimming = Channel.value( [file(params.reads_R1), file(params.reads_R2)] )
	} 
	else {
		Trimming = Channel.value( [file(params.reads_R1), "null"] )
	}
}

//When single-end reads are used, the input tuple (singleton) will not match input set 
//cardinality declared by 'trim' process (pair), so I push two mock files in the channel,
//and then I take only the first two files
mocktrim = Channel.from("null")
process trim {

	input:
   	set file(reads1), file(reads2) from Trimming.concat(mocktrim).flatMap().take(2).buffer(size : 2)
	file(adapters) from Channel.from( file(params.adapters) )
	file(artifacts) from Channel.from( file(params.artifacts) )
	file(phix174ill) from Channel.from( file(params.phix174ill) )

	output:

	file("${params.prefix}_trimmed*.fq") into todecontaminate
	file("${params.prefix}_trimmed*.fq") into trimmedreads

	when:
	params.mode == "QC" || params.mode == "complete"

   	script:
	"""	
	#Measures execution time
	sysdate=\$(date)
	starttime=\$(date +%s.%N)
	echo \"Performing Quality Control. STEP 2 [Trimming] at \$sysdate\" >>  $mylog
	echo \" \" >>  $mylog

	#Sets the maximum memory to the value requested in the config file
	maxmem=\$(echo ${task.memory} | sed 's/ //g' | sed 's/B//g')

	#Defines command for trimming of adapters and low quality bases
	if [ \"$params.library\" = \"paired\" ]; then
		CMD=\"bbduk.sh -Xmx\"\$maxmem\" in=$reads1 in2=$reads2 out=${params.prefix}_trimmed_R1_tmp.fq out2=${params.prefix}_trimmed_R2_tmp.fq outs=${params.prefix}_trimmed_singletons_tmp.fq ktrim=r k=$params.kcontaminants mink=$params.mink hdist=$params.hdist qtrim=rl trimq=$params.Quality  minlength=$params.minlength ref=$adapters qin=$params.Pcoding threads=${task.cpus} tbo tpe ow\"
	else
		CMD=\"bbduk.sh -Xmx\"\$maxmem\" in=$reads1 out=${params.prefix}_trimmed_tmp.fq ktrim=r k=$params.kcontaminants mink=$params.mink hdist=$params.hdist qtrim=rl trimq=$params.Quality  minlength=$params.minlength ref=$adapters qin=$params.Pcoding threads=${task.cpus} tbo tpe ow\"
	fi

	#Trims adapters and low quality bases	
	exec \$CMD 2>&1 | tee tmp.log
	
	#Logs some figures about sequences passing trimming
	echo  \"BBduk's trimming stats (trimming adapters and low quality reads): \" >>  $mylog
	sed -n '/Input:/,/Result:/p' tmp.log >>  $mylog
	echo \" \" >>  $mylog			
	if [ \"$params.library\" = \"paired\" ]; then
		unpairedR=\$(wc -l ${params.prefix}_trimmed_singletons_tmp.fq | cut -d\" \" -f 1)
		unpairedR=\$((\$unpairedR/4))
		echo  \"\$unpairedR singleton reads whose mate was trimmed shorter preserved\" >>  $mylog
		echo \" \" >>  $mylog
	fi

	
	#Logs version of the software and executed command (BBMap prints on stderr)
	version=\$(bbduk.sh --version 2>&1 >/dev/null | grep \"BBMap version\") 
	echo \"Using bbduk.sh in \$version \" >>  $mylog
	echo \"Using adapters in $params.adapters \" >>  $mylog
	echo \"Using synthetic contaminants in $params.phix174ill and in $params.artifacts \" >>  $mylog

	
			
	if [ \"$params.library\" = \"paired\" ]; then
		unpairedR=\$(wc -l ${params.prefix}_trimmed_singletons_tmp.fq | cut -d\" \" -f 1)
		unpairedR=\$((\$unpairedR/4))
		echo  \"\$unpairedR singleton reads whose mate was trimmed shorter preserved\" >>  $mylog
		echo \" \" >>  $mylog
	fi

	#Defines command for removing synthetic contaminants
	if [ \"$params.library\" = \"paired\" ]; then
		bbduk.sh -Xmx\"\$maxmem\" in=${params.prefix}_trimmed_R1_tmp.fq in2=${params.prefix}_trimmed_R2_tmp.fq out=${params.prefix}_trimmed_R1.fq out2=${params.prefix}_trimmed_R2.fq k=31 ref=$phix174ill,$artifacts qin=$params.Pcoding threads=${task.cpus} ow
	else
		bbduk.sh -Xmx\"\$maxmem\" in=${params.prefix}_trimmed_tmp.fq out=${params.prefix}_trimmed.fq k=31 ref=$phix174ill,$artifacts qin=$params.Pcoding threads=${task.cpus} ow
	fi


	#Removes synthetic contaminants and logs some figures (singleton read file, 
	#that exists iif the library layout was 'paired')
	if [ \"$params.library\" = \"paired\" ]; then
	bbduk.sh -Xmx\"\$maxmem\" in=${params.prefix}_trimmed_singletons_tmp.fq out=${params.prefix}_trimmed_singletons.fq k=31 ref=$phix174ill,$artifacts qin=$params.Pcoding threads=${task.cpus} ow
		
		
	fi
	
	#Removes tmp files. This avoids adding them to the output channels
	rm -rf ${params.prefix}_trimmed*_tmp.fq 
	"""
}



/**
	Quality control - STEP 3. Decontamination. Removes external organisms' contamination, 
	using a previously created index. When paired-end reads are used, decontamination is 
	carried on idependently on paired reads and on singleton reads thanks to BBwrap, 
	that calls BBmap once on the paired reads and once on the singleton ones, merging
	 the results on a single output file.

	Two files are outputted: the FASTQ of the decontaminated reads (including both
	paired-reads and singletons) and that of the contaminating reads (that can be
	used for refinements/checks).
	Please note that if keepQCtmpfile is set to false, the file of the contaminating 
	reads is discarded

	TODO: use BBsplit for multiple organisms decontamination, or fix the ref to a 
	FASTA file pangenome
*/

//When single-end reads are used, the input tuple (singleton) will not match input set 
//cardinality declared by 'trim' process (triplet), so I push two mock files in the channel,
//and then I take only the first three files.
mockdecontaminate = Channel.from("null", "null")

process decontaminate {

	publishDir  workingdir, mode: 'copy', pattern: "*_clean.fq"
		
	input:
	set file(infile1), file(infile2), file(infile12) from todecontaminate.concat(mockdecontaminate).flatMap().take(3).buffer(size : 3)
	file(RefGenome) from Channel.from( file(params.RefGenome, type: 'dir') )
	
	output:

	file "${params.prefix}_clean.fq" into decontaminatedreads
	file "${params.prefix}_clean.fq" into toprofiletaxa

	when:
	params.mode == "QC" || params.mode == "complete"

	script:
	"""
	#Measures execution time
	sysdate=\$(date)
	starttime=\$(date +%s.%N)
	echo \"Performing Quality Control. STEP 3 [Decontamination] at \$sysdate\" >>  $mylog
	echo \" \" >>  $mylog

	#Sets the maximum memory to the value requested in the config file
	maxmem=\$(echo ${task.memory} | sed 's/ //g' | sed 's/B//g')
	
	if [ \"$params.library\" = \"paired\" ]; then
		CMD=\"bbwrap.sh  -Xmx\"\$maxmem\" mapper=bbmap append=t in1=$infile1,$infile12 in2=$infile2,null outu=${params.prefix}_clean.fq outm=${params.prefix}_cont.fq minid=$params.mind maxindel=$params.maxindel bwr=$params.bwr bw=12 minhits=2 qtrim=rl trimq=$params.Quality path=$RefGenome qin=$params.Pcoding threads=${task.cpus} untrim quickmatch fast ow\"
	else
		CMD=\"bbwrap.sh  -Xmx\"\$maxmem\" mapper=bbmap append=t in1=$infile1 outu=${params.prefix}_clean.fq outm=${params.prefix}_cont.fq minid=$params.mind maxindel=$params.maxindel bwr=$params.bwr bw=12 minhits=2 qtrim=rl trimq=$params.Quality path=$RefGenome qin=$params.Pcoding threads=${task.cpus} untrim quickmatch fast ow\"
	fi

	#Decontaminates
	exec \$CMD 2>&1 | tee tmp.log
	
	

	#Informations about decontaminated/contaminated reads
	echo  \"BBwrap's human decontamination stats (paired reads): \" >>  $mylog
	sed -n '/Read 1 data:/,/N Rate:/p' tmp.log | head -17 >>  $mylog
	echo \" \" >>  $mylog
	sed -n '/Read 2 data:/,/N Rate:/p' tmp.log >>  $mylog
	echo \" \" >>  $mylog
	
	#Information of the software and executed command
	version=\$(bbwrap.sh --version 2>&1 >/dev/null | grep \"BBMap version\") 
	echo \"Using bbwrap.sh in \$version \" >>  $mylog
	echo \"Using contaminant (pan)genome indexed in $params.RefGenome \" >>  $mylog
	echo \" \" >>  $mylog


	nClean=\$(wc -l ${params.prefix}_clean.fq | cut -d\" \" -f 1)
	nClean=\$((\$nClean/4))
	nCont=\$(wc -l ${params.prefix}_cont.fq | cut -d\" \" -f 1)
	nCont=\$((\$nCont/4))
	echo \"\$nClean reads survived decontamination (\$nCont reads removed)\" >>  $mylog
	echo \" \" >>  $mylog

	#Measures and log execution time
	endtime=\$(date +%s.%N)
	exectime=\$(echo \"\$endtime \$starttime\" | awk '{print \$1-\$2}')
	sysdate=\$(date)
	echo \"STEP 3 (Quality Control) terminated at \$sysdate (\$exectime seconds)\" >>  $mylog
	echo \" \" >>  $mylog
	echo \"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\" >>  $mylog
	echo \"\" >>  $mylog
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

	publishDir  workingdir, mode: 'copy', pattern: "*.tsv"
	
	input:
	file(infile) from toprofiletaxa
	file(mpa_pkl) from Channel.from( file(params.mpa_pkl) )


    output:

	file "${params.prefix}_${params.date}*.txt" into toprofilefunctionbugs

	when:
	params.mode == "complete"

	script:
	"""
	#Measures execution time
	sysdate=\$(date)
	starttime=\$(date +%s.%N)
	echo \"Performing Community Characterisation. STEP 1 [Taxonomic binning and profiling] at \$sysdate\" >> $mylog
	echo \" \" >> $mylog
	
	#If a file with the same name is already present, Metaphlan2 will crash
	rm -rf ${params.prefix}_bt2out.txt
	
	#Defines command for estimating abundances
	CMD=\"metaphlan2.py --input_type fastq --bowtie2out=${params.prefix}.bt2 -t rel_ab_w_read_stats --bt2_ps $params.bt2options --nproc ${task.cpus} $infile -o ${params.prefix}_${params.date}_matrix.txt \"


	#Estimates microbial abundances
	exec \$CMD 2>&1 | tee tmp.log


	#Generating Relative abundance and read count outputs for every taxonomic level
	cp ${params.prefix}_${params.date}_matrix.txt $matrixdir

	cut -f1,2 ${params.prefix}_${params.date}_matrix.txt > $matrixdir/Abundance/${params.prefix}_${params.date}_abundance.txt
	cut -f1,5 ${params.prefix}_${params.date}_matrix.txt > $matrixdir/Read_count/${params.prefix}_${params.date}_counts.txt


	#Measures and log execution time
	endtime=\$(date +%s.%N)
	exectime=\$(echo \"\$endtime \$starttime\" | awk '{print \$1-\$2}')
	sysdate=\$(date)
	echo \"\" >> $mylog
	echo \"STEP 1 (Community Characterisation) terminated at \$sysdate (\$exectime seconds)\" >> $mylog
	echo \" \" >> $mylog
	echo \"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\" >> $mylog
	echo \"\" >> $mylog
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
	
	publishDir  workingdir, mode: 'copy', pattern: "*.{html,txt}"
	  	
	input:
   	set val(step), file(reads), val(label), val(stem) from toQC

	output:
	file "${params.prefix}*_fastqc.html" 

	when:
	params.mode == "QC" | params.mode == "complete"

   	script:
	"""	
	
	#Logs version of the software and executed command
	version=\$(fastqc --version) 
	fastqc --quiet --noextract --format fastq --outdir=. --threads ${task.cpus} $reads
	

	"""	
}
