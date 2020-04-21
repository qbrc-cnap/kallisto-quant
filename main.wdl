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
    }

    call multiqc.create_qc as experimental_qc {
        input:
            kallisto_stdout = single_sample_process.kallisto_stdout,
            r1_fastqc_zips = fastqc_for_read1.fastqc_zip,
            r2_fastqc_zips = fastqc_for_read2.fastqc_zip
    }

   call concatenate {
        input:
            main_abundance_files = single_sample_process.abundance_tsv,
            merged_quants_tag = merged_quants_tag,
    }

    call run_tximport {
        input:
            main_abundance_files = single_sample_process.abundance_tsv,
            transcript_to_gene_mapping = transcript_to_gene_mapping,
            human_tag = human_tag,
            mouse_tag = mouse_tag,
            is_pdx = is_pdx
    }

    call zip_results {
        input:
            zip_name = output_zip_name,
            multiqc_report = experimental_qc.report,
            main_tpm_file = concatenate.main_tpm_matrix,
            gene_level_counts = run_tximport.gene_level_counts
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


task run_tximport {
        Array[File] main_abundance_files
        File transcript_to_gene_mapping
        String human_tag
        String mouse_tag
        Boolean is_pdx

        String output_base = "gene_level_counts.tximport"

        Int disk_size = 50

        command {
            if [ "${is_pdx}" = "true" ]
            then
                Rscript /opt/software/run_tximport.R \
                    ${transcript_to_gene_mapping} \
                    abundance.tsv \
                    ${output_base} \
                    true \
                    ${human_tag} \
                    ${mouse_tag} \
                    ${sep=" " main_abundance_files}
            else
                Rscript /opt/software/run_tximport.R \
                    ${transcript_to_gene_mapping} \
                    abundance.tsv \
                    ${output_base} \
                    false \
                    ${human_tag} \
                    ${mouse_tag} \
                    ${sep=" " main_abundance_files}
            fi
        }

        output {
            Array[File] gene_level_counts = glob("${output_base}*tsv")
        }

        runtime {
            docker: "docker.io/blawney/kallisto-quant:v0.0.2"
            cpu: 2
            memory: "4 G"
            disks: "local-disk " + disk_size + " HDD"
            preemptible: 0 
        }
}


task concatenate {
        # This concatenates the abundance files into a TPM matrix

        Array[File] main_abundance_files
        String merged_quants_tag

        Int disk_size = 30

        command {
            /usr/bin/python3 /opt/software/concat_tpm.py -o "${merged_quants_tag}.tsv" -s "abundance.tsv" ${sep=" " main_abundance_files}
        }

        output {
            File main_tpm_matrix = "${merged_quants_tag}.tsv"
        }

        runtime {
            docker: "docker.io/blawney/kallisto-quant:v0.0.2"
            cpu: 2
            memory: "4 G"
            disks: "local-disk " + disk_size + " HDD"
            preemptible: 0 
        }

    }

task zip_results {

    String zip_name 
    File multiqc_report
    File main_tpm_file
    Array[File] gene_level_counts

    Int disk_size = 100

    command {

        mkdir report
        mkdir report/qc
        mkdir report/tpm_quantifications
        mkdir report/gene_level_quantifications

        mv ${multiqc_report} report/qc/
        mv ${main_tpm_file} report/tpm_quantifications
        mv ${sep=" " gene_level_counts} report/gene_level_quantifications

        zip -r "${zip_name}.zip" report
    }

    output {
        File zip_out = "${zip_name}.zip"
    }

    runtime {
        docker: "docker.io/blawney/kallisto-quant:v0.0.2"
        cpu: 2
        memory: "6 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}
