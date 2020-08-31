import "single_sample_kallisto.wdl" as single_sample_kallisto
import "multiqc.wdl" as multiqc
import "fastqc.wdl" as fastqc
import "report.wdl" as reporting


workflow KallistoQuantWorkflow{
    # This workflow is a 'super' workflow that parallelizes
    # RNA-seq analysis over multiple samples

    Array[File] r1_files
    Array[File] r2_files
    Int trim_length_bp
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

        call trim_reads {
            input:
            r1 = item.left,
            r2 = item.right,
            trim_length_bp = trim_length_bp
        }

        call fastqc.run_fastqc as fastqc_for_read1 {
            input:
                fastq = trim_reads.trimmed_r1
        }

        call fastqc.run_fastqc as fastqc_for_read2 {
            input:
                fastq = trim_reads.trimmed_r2
        }

        call single_sample_kallisto.SingleSampleKallistoWorkflow as single_sample_process{
            input:
                r1_fastq = trim_reads.trimmed_r1,
                r2_fastq = trim_reads.trimmed_r2,
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
            source activate r36

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
            docker: "docker.io/blawney/kallisto-quant:v0.0.3"
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
            docker: "docker.io/blawney/kallisto-quant:v0.0.3"
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
        docker: "docker.io/blawney/kallisto-quant:v0.0.3"
        cpu: 2
        memory: "6 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}

task trim_reads {

    File r1
    File r2
    Int trim_length_bp

    String suffix="_R1.fastq.gz"

    # Extract the samplename from the fastq filename
    String sample_name = basename(r1, suffix)

    Int disk_size = 200


    command {
        java -jar /opt/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
            -trimlog ${sample_name}.trim.log \
            -summary ${sample_name}.trim_summary.log \
            ${r1} ${r2} \
            -baseout ${sample_name}.trimmed.fastq.gz \
            CROP:${trim_length_bp}

        mv "${sample_name}.trimmed_1P.fastq.gz" "${sample_name}_R1.fastq.gz"
        mv "${sample_name}.trimmed_2P.fastq.gz" "${sample_name}_R2.fastq.gz"
    }

    output {
        File trimmed_r1 = "${sample_name}_R1.fastq.gz"
        File trimmed_r2 = "${sample_name}_R2.fastq.gz"
    }

    runtime {
        docker: "docker.io/blawney/kallisto-quant:v0.0.3"
        cpu: 4
        memory: "12 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}
