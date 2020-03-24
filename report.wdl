task generate_report {

    Array[String] r1_files
    Array[String] r2_files
    File annotations
    String genome
    Array[File] sleuth_results
    String git_repo_url
    String git_commit_hash
    String normalized_counts_suffix
    String sleuth_output_suffix
    String versus_sep
    Float qval_threshold
    String pca_suffix
    String top_heatmap_suffix
    Int max_transcripts
    Int num_bootstraps

    Int disk_size = 15

    command <<<

        # make a json file with various parameters:
        echo "{" >> config.json
        echo '"genome": "${genome}",' >>config.json
        echo '"pca_suffix": "${pca_suffix}",' >>config.json
        echo '"qval": "${qval_threshold}",' >>config.json
        echo '"sleuth_output_suffix": "${sleuth_output_suffix}",' >>config.json
        echo '"normalized_counts_suffix": "${normalized_counts_suffix}",' >>config.json
        echo '"versus_sep": "${versus_sep}",' >>config.json
        echo '"git_repo": "${git_repo_url}",' >>config.json
        echo '"git_commit": "${git_commit_hash}",' >>config.json
        echo '"num_hits": "${max_transcripts}",' >>config.json
        echo '"num_bootstraps": "${num_bootstraps}",' >>config.json
        echo '"top_heatmap_suffix": "${top_heatmap_suffix}"}' >>config.json

        generate_report.py \
          -r1 ${sep=" " r1_files} \
          -r2 ${sep=" " r2_files} \
          -a ${annotations} \
          -d ${sep=" " sleuth_results} \
          -j config.json \
          -t /opt/report/report.md \
          -o completed_report.md

        pandoc -H /opt/report/report.css -s completed_report.md -o analysis_report.html
    >>>

    output {
        File report = "analysis_report.html"
    }

    runtime {
        docker: "docker.io/blawney/kallisto:v0.0.2"
        cpu: 2
        memory: "6 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}