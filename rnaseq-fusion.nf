#! /usr/bin/env nextflow

// Copyright (C) 2017 IARC/WHO

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


params.CTAT_folder = '.'
params.input_folder = '.'
params.output_folder= "results_fusion"
params.mem  = 2
params.cpu  = 2
params.fastq_ext = "fq.gz"
params.suffix1 = "_1"
params.suffix2 = "_2"
params.junction_suffix = "Chimeric.SJ.out.junction"

params.twopass  = null

params.help = null

log.info ""
log.info "--------------------------------------------------------"
log.info "  <PROGRAM_NAME> <VERSION>: <SHORT DESCRIPTION>         "
log.info "--------------------------------------------------------"
log.info "Copyright (C) IARC/WHO"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "--------------------------------------------------------"
log.info ""

if (params.help) {
    log.info "--------------------------------------------------------"
    log.info "  USAGE                                                 "
    log.info "--------------------------------------------------------"
    log.info ""
    log.info "nextflow run iarcbioinfo/rnaseq-transcript-nf [-with-docker] [OPTIONS]"
    log.info ""
    log.info "Mandatory arguments:"
    log.info '    --input_folder   FOLDER                  Folder containing fastq files and STAR junction files.'
    log.info '    --CTAT_folder     FOLDER                  Folder with STAR-Fusion bundle (CTAT).'
    log.info ""
    log.info "Optional arguments:"
    log.info '    --output_folder     STRING                Output folder (default: results_fusion).'
    log.info '    --cpu          INTEGER                 Number of cpu used by bwa mem and sambamba (default: 2).'
    log.info '    --mem          INTEGER                 Size of memory used for mapping (in GB) (default: 2).' 
    log.info ""
    log.info "Flags:"
    log.info "--<FLAG>                                                    <DESCRIPTION>"
    log.info ""
    exit 0
} else {
/* Software information */
   log.info "input_folder = ${params.input_folder}"
   log.info "cpu          = ${params.cpu}"
   log.info "mem          = ${params.mem}"
   log.info "output_folder= ${params.output_folder}"
   log.info "help:                               ${params.help}"
}

file(params.input_folder).listFiles().println()

if ( file(params.input_folder).listFiles().findAll { it.name ==~ /.*junction/ }.size() > 0){
       println "Junction files found, proceed with fusion genes discovery"

//test1 = Channel
//    .fromPath( params.input_folder )
//    .println()


// Gather files ending with _1 suffix
   reads1 = Channel
    .fromPath( params.input_folder+'/*'+params.suffix1+'.'+params.fastq_ext )
    .map {  path -> [ path.name.replace("${params.suffix1}.${params.fastq_ext}",""), path ] }


// Gather files ending with _2 suffix
   reads2 = Channel
    .fromPath( params.input_folder+'/*'+params.suffix2+'.'+params.fastq_ext )
    .map {  path -> [ path.name.replace("${params.suffix2}.${params.fastq_ext}",""), path ] }


//   Channel
//	.fromPath(params.input_folder+'/*junction' )
//	.set{ test}
//   test.println()
// Gather files ending with junction
   junctions = Channel
    .fromPath( params.input_folder+'/*junction')
    .map {  path -> [ path.name.replace("STAR.","").replace(".Chimeric.SJ.out.junction",""), path ] }
//    .view()
//    .println()
//    .map {  path -> [ path.name.replace("STAR.",""), path ] }
//    .println()

// Match the pairs on two channels having the same 'key' (name) and emit a new pair containing the expected files
   input_triplet = reads1
                     .phase(reads2)
		     .map { pair1, pair2 -> [pair1[0], pair1[1], pair2[1] ] }
 //   		     .view()


   input_triplet = input_triplet.phase(junctions)
//				.view()
//		     .map { pairs,  -> [ pair1[1], pair2[1] ] }
   input_triplet = input_triplet
		     .map { pairs, junction -> [ pairs[1], pairs[2], junction[1] ] }
//		     .view()
//input_triplet = input_triplet.view()


//       fastq_files      = Channel.fromFilePairs( params.input_folder+'/*{1.fq.gz,2.fq.gz,junction}')
// 				 .println()
//	keys1 = file(params.input_folder).listFiles().findAll { it.name ==~ /.*${params.suffix1}.${params.fastq_ext}/ }.collect { it.getName() }
//		   			.collect { it.replace("${params.suffix1}.${params.fastq_ext}",'') }
//
//	junction_files   = Channel.fromFilePairs( params.input_folder+'/*{1.fq.gz,junction}')
//				 .println()
}else{
       println "ERROR: input folder contains no fastq files"; System.exit(1)
}

println input_triplet

process STAR_Fusion {
	cpus params.cpu
	memory params.mem+'G'
	tag { file_tag }

	input:
	set pair1, pair2 , junction from input_triplet
	
	output:
	file("*") into outputs

	publishDir "${params.output_folder}/${file_tag}", mode: 'copy'	

	shell:
	file_tag=pair1[0].baseName

    	'''
	STAR-Fusion --genome_lib_dir !{params.CTAT_folder} -J !{junction} --left_fq !{pair1} --right_fq !{pair2} --output_dir . --FusionInspector validate --denovo_reconstruct
    	'''
}

