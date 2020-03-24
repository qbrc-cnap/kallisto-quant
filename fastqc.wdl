workflow test_fastqc {
    File fastq

    call run_fastqc{
        input:
            fastq=fastq
    }

    output {
        File fqc_zip = run_fastqc.fastqc_zip
    }
}

task run_fastqc {
    File fastq

    String bn = basename(fastq)
    String target_zip = sub(bn, "\\.fastq\\.gz", "_fastqc.zip")
    Int disk_size = 20

    command {
        fastqc ${fastq} -o .
    }

    output {
        File fastqc_zip = "${target_zip}"
    }

    runtime {
        zones: "us-east4-c"
        docker: "docker.io/blawney/kallisto:v0.0.2"
        cpu: 2
        memory: "4 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}
