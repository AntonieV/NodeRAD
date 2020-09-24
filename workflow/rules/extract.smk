rule edit_distance:  # calculates edit distances with minimap2
    input:
        config["resources"]["fastq-data"]
    output:
        "results/test.txt"
    params:
        ""
    log:
        "logs/test.log"
    conda:
        "../envs/mappy.yaml"
    script:
        "../scripts/edit_distance.py"

