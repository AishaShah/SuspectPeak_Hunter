rule plot_diffbind_heatmap_for_Validation_Samples:
    input:
        validation_samples="validation_samples.tsv",
        peaks=rules.Overlap_SLs_and_Validation_Samples.input.peaks,
        bams=rules.Overlap_SLs_and_Validation_Samples.input.bams,
    output:
        diffbind_array_noSL="diffbind_array.noSL.validation_samples.tsv",
        diffbind_array="diffbind_array.validation_samples.tsv",
        diffbind_plot="diffbind_heatmaps.pdf",
        diffbind_plot_noSL="diffbind_heatmaps.noSL.pdf"
    script:
        "scripts/diffbind_analysis.R"
    params:
        results_dir=config["results_dir"],
        num_reads=config["DownSample"]["Reads"],
        seacr_threshold=config["peak_calling"]["threshold"],
        min_depth=config["peak_calling"]["min_depth"]
    shell:
        """
        # Add new header to the diffbind array files
        echo -e "SampleID\tSampleID2\tStage\tTarget\tPeakCaller\tFactor\tReplicate\tbamReads\tPeaks" > {output.diffbind_array}
        echo -e "SampleID\tSampleID2\tStage\tTarget\tPeakCaller\tFactor\tReplicate\tbamReads\tPeaks" > {output.diffbind_array_noSL}

        # Process the file, remove the first line, and append transformed data
        tail -n +2 {input.validation_samples} | \
        awk -v bam_path={params.results_dir}/05.bootstrapping/round_Validation/01.DownSampling.{params.num_reads}/ \
            -v peak_path={params.results_dir}/05.bootstrapping/round_Validation/04.Peaks.DS_{params.num_reads}/filtered_mindepth_{params.min_depth}.SEACR_{params.seacr_threshold}/ \
            '{{print $4 "\t" $5 "\t" $6 "\t" "bed" "\t" $7 "\t" $9 "\t" bam_path $4 ".no_MT.sorted.bam" "\t" peak_path "noSL." $4 ".stringent.bed"}}' >> {output.diffbind_array_noSL}

        tail -n +2 {input.validation_samples} | \
        awk -v bam_path={params.results_dir}/05.bootstrapping/round_Validation/01.DownSampling.{params.num_reads}/ \
            -v peak_path={params.results_dir}/05.bootstrapping/round_Validation/04.Peaks.DS_{params.num_reads}/filtered_mindepth_{params.min_depth}.SEACR_{params.seacr_threshold}/ \
            '{{print $4 "\t" $5 "\t" $6 "\t" "bed" "\t" $7 "\t" $9 "\t" bam_path $4 ".no_MT.sorted.bam" "\t" peak_path $4 ".stringent.bed"}}' >> {output.diffbind_array}
        # Run the external R script
        Rscript scripts/diffbind_analysis.R {output.diffbind_array} {output.diffbind_array_noSL} {output.diffbind_plot} {output.diffbind_plot_noSL}
        """