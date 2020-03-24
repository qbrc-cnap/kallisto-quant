task create_qc {
    Array[File] kallisto_stdout
    Array[File] r1_fastqc_zips
    Array[File]? r2_fastqc_zips

    Int disk_size = 30

    command {
        multiqc .
    }

    output {
        File report = "multiqc_report.html"
    }
        
    runtime {
        zones: "us-east4-c"
        docker: "docker.io/blawney/kallisto:v0.0.2"
        cpu: 2
        memory: "3 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}
