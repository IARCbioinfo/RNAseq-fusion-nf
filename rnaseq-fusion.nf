#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

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


/* ============================
 * PARAMETERS
 * ============================ */

params.CTAT_folder      = '.'
params.input_folder     = '.'
params.input_file       = null
params.output_folder    = 'results_fusion'
params.mem              = 2
params.cpu              = 2
params.fastq_ext        = 'fq.gz'
params.suffix1          = '_1'
params.suffix2          = '_2'
params.junction_suffix  = 'Chimeric.SJ.out.junction'
params.junctions        = null
params.starfusion_path  = '/usr/local/src/STAR-Fusion/STAR-Fusion'
params.help             = null

log.info ""
log.info "--------------------------------------------------------"
log.info "  rnaseq-fusion-nf v1.1: nextflow pipeline to run STAR-fusion "
log.info "--------------------------------------------------------"
log.info "Copyright (C) IARC/WHO"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "--------------------------------------------------------"
log.info ""

if (params.help) {
    log.info "--------------------------------------------------------"
    log.info "  USAGE            "
    log.info "--------------------------------------------------------"
    log.info ""
    log.info "nextflow run IARCbioinfo/RNAseq-fusion-nf --input_folder fastq/ --CTAT_folder GRCh38_CTAT/  [-with-docker] [OPTIONS]"
    log.info ""
    log.info "Mandatory arguments:"
    log.info '    --input_folder    FOLDER                  Folder containing fastq files and STAR junction files.'
    log.info '    --CTAT_folder     FOLDER                  Folder with STAR-Fusion bundle (CTAT).'
    log.info ""
    log.info "Optional arguments:"
    log.info '    --input_file      STRING                  Input file (comma-separated) with 4 columns:'
    log.info '                                              SM(sample name), pair1 (path to fastq pair 1), '
    log.info '                                              pair2 (path to fastq pair 2), and junction (path to junction file).'
    log.info '    --output_folder   STRING                  Output folder (default: results_fusion).'
    log.info '    --fastq_ext       STRING                  Extension of fastq files (default: fq.gz).'
    log.info '    --suffix1         STRING                  Suffix of 1st element of fastq files pair (default: _1).'
    log.info '    --suffix2         STRING                  Suffix of 2nd element of fastq files pair (default: _2).'
    log.info '    --junction_suffix STRING                  Suffix of STAR chimeric junction files (default: Chimeric.SJ.out.junction).'
    log.info '    --starfusion_path STRING                  Path to STAR-fusion executable (default: /usr/local/src/STAR-Fusion/STAR-Fusion for singularity image'
    log.info '    --cpu             INTEGER                 Number of cpu used by bwa mem and sambamba (default: 2).'
    log.info '    --mem             INTEGER                 Size of memory used for mapping (in GB) (default: 2).' 
    log.info ""
    log.info "Flags:"
    log.info "    --junctions                               Option to use STAR junction files (default: null)."
    log.info ""
    exit 0
} else {
/* Software information */
   log.info "input_folder    = ${params.input_folder}"
   log.info "input_file      = ${params.input_file}"
   log.info "cpu             = ${params.cpu}"
   log.info "mem             = ${params.mem}"
   log.info "output_folder   = ${params.output_folder}"
   log.info "CTAT_folder     = ${params.CTAT_folder}"
   log.info "fastq_ext       = ${params.fastq_ext}"
   log.info "suffix1         = ${params.suffix1}"
   log.info "suffix2         = ${params.suffix2}"
   log.info "junction_suffix = ${params.junction_suffix}"
   log.info "starfusion_path = ${params.starfusion_path}"
   log.info "junctions       = ${params.junctions}"
   log.info "help:             ${params.help}"
}

/* ============================
 * PROCESS
 * ============================ */

process STAR_FUSION {

    cpus params.cpu
    memory "${params.mem}G"
    tag { sample_id }

    input:
    tuple val(sample_id), path(pair1), path(pair2), path(junction)
    path CTAT_folder

    output:
    path "FusionInspector*", emit: fi
    path "star-fusion*", emit: sf

    publishDir "${params.output_folder}/${sample_id}", mode: 'copy'

    script:
    def sf_junction = params.junctions ? "-J ${junction}" : ""
    """
    input_txt="${sample_id}\t${pair1}\t${pair2}"

    echo -e "\$input_txt" > input.txt

    ${params.starfusion_path} \
        --genome_lib_dir \$PWD/${CTAT_folder} \
        ${sf_junction} \
        --samples_file input.txt \
        --output_dir . \
        --FusionInspector validate \
        --denovo_reconstruct \
        --examine_coding_effect \
        --CPU ${params.cpu}
    """
}

/* ============================
 * WORKFLOW
 * ============================ */

workflow {

    def input_triplet

    if (params.input_file) {

        input_triplet = Channel
            .fromPath(params.input_file)
            .splitCsv(header: true, sep: '\t', strip: true)
            .map { row ->
                tuple(
                    row.SM,
                    file(row.pair1),
                    file(row.pair2),
                    file(row.junction)
                )
            }

    } else {

        reads = Channel
            .fromFilePairs(
                "${params.input_folder}/*{${params.suffix1},${params.suffix2}}.${params.fastq_ext}"
            )
            .map { id, files ->
                tuple(id, files[0], files[1])
            }

        if (params.junctions) {

            junctions = Channel
                .fromPath("${params.input_folder}/*${params.junction_suffix}")
                .map { path ->
                    tuple(
                        path.name
                            .replace('STAR.', '')
                            .replace(".${params.junction_suffix}", ''),
                        path
                    )
                }

            input_triplet = reads
                .join(junctions)
                .map { id, p1, p2, junction ->
                    tuple(id, p1, p2, junction)
                }

        } else {

            input_triplet = reads
                .map { id, p1, p2 ->
                    tuple(id, p1, p2, file('NO_FILE'))
                }
        }
    }

    STAR_FUSION(
        input_triplet,
        Channel.value(file(params.CTAT_folder))
    )
}
