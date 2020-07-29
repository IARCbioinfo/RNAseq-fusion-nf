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
params.junctions = null

params.help = null

log.info ""
log.info "--------------------------------------------------------"
log.info "  rnaseq-fusion-nf v1.0: nextflow pipeline to run STAR-fusion "
log.info "--------------------------------------------------------"
log.info "Copyright (C) IARC/WHO"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "--------------------------------------------------------"
log.info ""

if (params.help) {
    log.info "--------------------------------------------------------"
    log.info "  USAGE nextflow run rnaseq-fusion-nf --input_folder fastq/ --CTAT_folder GRCh38_CTAT/                                                 "
    log.info "--------------------------------------------------------"
    log.info ""
    log.info "nextflow run iarcbioinfo/rnaseq-transcript-nf [-with-docker] [OPTIONS]"
    log.info ""
    log.info "Mandatory arguments:"
    log.info '    --input_folder    FOLDER                  Folder containing fastq files and STAR junction files.'
    log.info '    --CTAT_folder     FOLDER                  Folder with STAR-Fusion bundle (CTAT).'
    log.info ""
    log.info "Optional arguments:"
    log.info '    --output_folder   STRING                  Output folder (default: results_fusion).'
    log.info '    --fastq_ext       STRING                  Extension of fastq files (default: fq.gz).'
    log.info '    --suffix1         STRING                  Suffix of 1st element of fastq files pair (default: _1).'
    log.info '    --suffix2         STRING                  Suffix of 2nd element of fastq files pair (default: _2).'
    log.info '    --junction_suffix STRING                  Suffix of STAR chimeric junction files (default: Chimeric.SJ.out.junction).'
    log.info '    --cpu             INTEGER                 Number of cpu used by bwa mem and sambamba (default: 2).'
    log.info '    --mem             INTEGER                 Size of memory used for mapping (in GB) (default: 2).' 
    log.info ""
    log.info "Flags:"
    log.info "--junctions          Option to use STAR junction files (default: null)."
    log.info ""
    exit 0
} else {
/* Software information */
   log.info "input_folder    = ${params.input_folder}"
   log.info "cpu             = ${params.cpu}"
   log.info "mem             = ${params.mem}"
   log.info "output_folder   = ${params.output_folder}"
   log.info "CTAT_folder     = ${params.CTAT_folder}"
   log.info "fastq_ext       = ${params.fastq_ext}"
   log.info "suffix1         = ${params.suffix1}"
   log.info "suffix2         = ${params.suffix2}"
   log.info "junction_suffix = ${params.junction_suffix}"
   log.info "junctions       = ${params.junctions}"
   log.info "help:             ${params.help}"
}


// Gather paired fastq files
   readPairs = Channel.fromFilePairs(params.input_folder +"/*{${params.suffix1},${params.suffix2}}" +'.'+ params.fastq_ext)
                      .map {  row -> [ row[0], row[1][0], row[1][1] ] }

   if(params.junctions){
   println "Gather STAR junction files"
   if ( file(params.input_folder).listFiles().findAll { it.name ==~ /.*junction/ }.size() > 0){
       println "Junction files found, proceed with fusion genes discovery"
   }else{
	println "ERROR: input folder contains no junction files"; System.exit(1)
   }
   junctions = Channel.fromPath( params.input_folder+'/*' +params.junction_suffix)
    .map {  path -> [ path.name.replace("STAR.","").replace(".${params.junction_suffix}",""), path ] }

// Match the pairs on two channels having the same 'key' (name) and emit a new pair containing the expected files
   input_triplet = readPairs.phase(junctions)
   input_triplet = input_triplet
		     .map { pairs, junction -> [ pairs[0],pairs[1], pairs[2], junction[1] ] }
   }else{
	println "Do not gather STAR junction files; STAR will be used for alignment"
 	input_triplet = readPairs.map { pairs -> [ pairs[0],pairs[1], pairs[2], 'NO_FILE' ] }
}


process STAR_Fusion {
	cpus params.cpu
	memory params.mem+'G'
	tag { file_tag }

	input:
	set file_tag, file(pair1), file(pair2) , file(junction) from input_triplet
	file(CTAT_folder) from file(params.CTAT_folder)	

	output:
	file("FusionInspector*") into FIoutputs
	file("star-fusion*") into SFoutputs

	publishDir "${params.output_folder}/${file_tag}", mode: 'copy'	

	shell:
	if(params.junctions){
		SF_junction="-J "+junction+" "
	}else{
	        SF_junction=" "
	}
    '''
	STAR-Fusion --genome_lib_dir $PWD/!{CTAT_folder} !{SF_junction} --left_fq !{pair1} --right_fq !{pair2} --output_dir . --FusionInspector validate --denovo_reconstruct --examine_coding_effect --CPU !{params.cpu}
    '''
}
