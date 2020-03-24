import "single_sample_kallisto.wdl" as single_sample_kallisto
import "multiqc.wdl" as multiqc
import "fastqc.wdl" as fastqc
import "report.wdl" as reporting


workflow KallistoQuantWorkflow{
    # This workflow is a 'super' workflow that parallelizes
    # RNA-seq analysis over multiple samples

    Array[File] r1_files
    Array[File] r2_files
    String genome
    File kallisto_index_path
    File transcript_to_gene_mapping
    String output_zip_name
    String git_repo_url
    String git_commit_hash
    Boolean is_pdx

    # if we are running PDx, we will also generate files for human and mouse separately.
    # these will be named similarly, e.g. 'A_vs_B.mouse.tpm.tsv' and 'A_vs_B.human.tpm.tsv'
    # and so on.
    String human_tag = "human"
    String mouse_tag = "mouse"
    String merged_quants_tag = "merged_tpm"

    Array[Pair[File, File]] fastq_pairs = zip(r1_files, r2_files)


    scatter(item in fastq_pairs){

        call fastqc.run_fastqc as fastqc_for_read1 {
            input:
                fastq = item.left
        }

        call fastqc.run_fastqc as fastqc_for_read2 {
            input:
                fastq = item.right
        }

        call single_sample_kallisto.SingleSampleKallistoWorkflow as single_sample_process{
            input:
                r1_fastq = item.left,
                r2_fastq = item.right,
                kallisto_index_path = kallisto_index_path
                        
            }

        call filter_for_pdx {
            input:
                is_pdx = is_pdx,
                abundance_tsv = single_sample_process.abundance_tsv,
                sample_name = single_sample_process.sample_name_output,
                human_tag = human_tag,
                mouse_tag = mouse_tag
        }
    }

    call multiqc.create_qc as experimental_qc {
        input:
            kallisto_stdout = single_sample_process.kallisto_stdout,
            r1_fastqc_zips = fastqc_for_read1.fastqc_zip,
            r2_fastqc_zips = fastqc_for_read2.fastqc_zip
    }

    # call reporting.generate_report as make_report {
    #     input:
    #         r1_files = r1_files,
    #         r2_files = r2_files,
    #         genome = genome,
    #         git_commit_hash = git_commit_hash,
    #         git_repo_url = git_repo_url
    # }

   call concatenate {
        input:
            main_abundance_files = single_sample_process.abundance_tsv,
            human_abundance_files = filter_for_pdx.human_tsv,
            mouse_abundance_files = filter_for_pdx.mouse_tsv,
            mouse_tag = mouse_tag,
            human_tag = human_tag,
            merged_quants_tag = merged_quants_tag,
            is_pdx = is_pdx
    }

    call zip_results {
        input:
            zip_name = output_zip_name,
            multiqc_report = experimental_qc.report,
            main_abundance_file = concatenate.main_tpm_matrix,
            human_abundance_file = concatenate.human_tpm_matrix,
            mouse_abundance_file = concatenate.mouse_tpm_matrix,
            #analysis_report = make_report.report,
            is_pdx = is_pdx
    }

    output {
        File zip_out = zip_results.zip_out
    }

    meta {
        workflow_title : "Kallisto quantification"
        workflow_short_description : "For determining transcript and gene abundances using Kallisto and tximport"
        workflow_long_description : "Use this workflow for performing pseudo-alignments with Kallisto, which produces estimated transcript abundances. Estimates of gene level quantification performed by tximport."
    }
}


task concatenate {
        # This concatenates the featureCounts count files into a 
        # raw count matrix.

        Array[File] main_abundance_files
        Array[File] human_abundance_files
        Array[File] mouse_abundance_files
        String human_tag
        String mouse_tag
        String merged_quants_tag
        Boolean is_pdx

        Int disk_size = 30

        command {
            if [ "${is_pdx}" = "true" ]
            then
                /usr/bin/python3 /opt/software/concat_tpm.py -o "${merged_quants_tag}.${mouse_tag}.tsv" -s "abundance.${mouse_tag}.tsv" ${sep=" " mouse_abundance_files}
                /usr/bin/python3 /opt/software/concat_tpm.py -o "${merged_quants_tag}.${human_tag}.tsv" -s "abundance.${human_tag}.tsv" ${sep=" " human_abundance_files}
                touch "${merged_quants_tag}.tsv"
            else
                /usr/bin/python3 /opt/software/concat_tpm.py -o "${merged_quants_tag}.tsv" -s "abundance.tsv" ${sep=" " main_abundance_files}
                touch "abundance.${human_tag}.tsv" ${sep=" " human_abundance_files}
                touch "abundance.${mouse_tag}.tsv" ${sep=" " mouse_abundance_files}
        fi
        }

        output {
            File main_tpm_matrix = "${merged_quants_tag}.tsv"
            File mouse_tpm_matrix = "${merged_quants_tag}.${mouse_tag}.tsv"
            File human_tpm_matrix = "${merged_quants_tag}.${human_tag}.tsv"
        }

        runtime {
            docker: "docker.io/blawney/kallisto:v0.0.2"
            cpu: 2
            memory: "4 G"
            disks: "local-disk " + disk_size + " HDD"
            preemptible: 0 
        }

    }

task filter_for_pdx {

    # Note that this assumes we have a human and mouse PDx situation
    # If this is not the case, then any files created here are irrelevant,
    # so the naming doesn't matter.

    Boolean is_pdx
    File abundance_tsv
    String sample_name
    String human_tag
    String mouse_tag

    Int disk_size = 30

    command {
        if [ "${is_pdx}" = "true" ]
        then
            head -1 ${abundance_tsv} > "${sample_name}.abundance.${human_tag}.tsv"
            grep -P "^ENST" ${abundance_tsv} >> "${sample_name}.abundance.${human_tag}.tsv"
            grep -vP "^ENST" ${abundance_tsv} >> "${sample_name}.abundance.${mouse_tag}.tsv"
        else
            touch "${sample_name}.abundance.${human_tag}.tsv"
            touch "${sample_name}.abundance.${mouse_tag}.tsv"
        fi
    }

    output {
        File human_tsv ="${sample_name}.abundance.${human_tag}.tsv"
        File mouse_tsv ="${sample_name}.abundance.${mouse_tag}.tsv"
    }

    runtime {
        docker: "docker.io/blawney/kallisto:v0.0.2"
        cpu: 2
        memory: "6 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}

task zip_results {

    String zip_name 
    File multiqc_report
    #File analysis_report
    #File gene_level_count_file
    #File transcript_level_abundance
    File main_abundance_file
    File human_abundance_file
    File mouse_abundance_file
    Boolean is_pdx

    Int disk_size = 100

    command {

        mkdir report
        mkdir report/qc
        mkdir report/quantifications

        mv ${multiqc_report} report/qc/

        if [ "${is_pdx}" = "true" ]
        then
            mv ${human_abundance_file} report/quantifications
            mv ${mouse_abundance_file} report/quantifications
        else
            mv ${main_abundance_file} report/quantifications
        fi

        zip -r "${zip_name}.zip" report
    }

    output {
        File zip_out = "${zip_name}.zip"
    }

    runtime {
        docker: "docker.io/blawney/kallisto:v0.0.2"
        cpu: 2
        memory: "6 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}
