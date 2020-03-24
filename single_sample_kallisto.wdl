workflow SingleSampleKallistoWorkflow {

    File r1_fastq
    File r2_fastq
    File kallisto_index_path

    # Extract the samplename from the fastq filename
    String sample_name = basename(r1_fastq, "_R1.fastq.gz")

    Int disk_size = 100

    call KallistoQuantification {
        input:
            r1_fastq = r1_fastq,
            r2_fastq = r2_fastq,
            kallisto_index_path = kallisto_index_path,
            sample_name = sample_name
    }

    output {
        File kallisto_stdout = KallistoQuantification.kallisto_stdout
        File abundance_tsv = KallistoQuantification.abundance_tsv
        File run_info = KallistoQuantification.run_info
        String sample_name_output = "${sample_name}"
    }
}

task KallistoQuantification {
    File r1_fastq
    File r2_fastq
    File kallisto_index_path
    String sample_name

    Int threads = 2
    Int disk_size = 30

    String stdout_log = "stdout.log"

    command {
        source activate r36

        kallisto quant -i ${kallisto_index_path} \
            -o "kallisto_results_${sample_name}" \
            -t ${threads} \
            ${r1_fastq} \
            ${r2_fastq} >> ${stdout_log} 2>&1

        mv "kallisto_results_${sample_name}/abundance.tsv" "${sample_name}.abundance.tsv"
        mv "kallisto_results_${sample_name}/run_info.json" "${sample_name}.run_info.json"
    }

    output {
        File kallisto_stdout = "${stdout_log}"
        File abundance_tsv = "${sample_name}.abundance.tsv"
        File run_info = "${sample_name}.run_info.json"
    }

    runtime {
        docker: "docker.io/blawney/kallisto:v0.0.2"
        cpu: 4
        memory: "20 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}
