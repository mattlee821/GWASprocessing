nextflow.enable.dsl = 2

workflow {
    PREPARE_MANIFEST()

    rows = PREPARE_MANIFEST.out.manifest.splitCsv(header: true, sep: '\t').map { row ->
            tuple(
                row.GWAS_location,
                row.phenotype,
                row.ancestry,
                row.author,
                row.year,
                row.PMID,
                row.sex,
                row.input_build,
                row.population,
                row.delimiter,
                row.format,
                row.study_yaml,
                row.study_id,
                row.output_id,
                row.row_hash,
                row.source_type,
                row.raw_dir,
                row.standardise_dir,
                row.output_prefix
            )
    }

    row_records = STANDARDISE_ROW(rows)
    UPDATE_STATE(row_records.records.collect())
}

process PREPARE_MANIFEST {
    tag "prepare-manifest"
    publishDir "${params.log_root}", mode: 'copy', overwrite: true

    output:
    path 'run_manifest.tsv', emit: manifest
    path 'run_manifest_summary.json', emit: summary

    script:
    """
    python3 "${params.pipeline_root}/bin/prepare_manifest.py" \\
      --manifest "${params.manifest}" \\
      --repo-root "${params.repo_root}" \\
      --work-root "${params.work_root}" \\
      --state-file "${params.work_root}/.standardise_state.tsv" \\
      --output run_manifest.tsv \\
      --summary run_manifest_summary.json \\
      --only-study "${params.only_study}" \\
      --force "${params.force}" \\
      --qcplot "${params.qcplot}" \\
      --row-limit "${params.row_limit}"
    """
}

process STANDARDISE_ROW {
    tag "${study_id}:${output_id}"

    input:
    tuple val(gwas_location), val(phenotype), val(ancestry), val(author), val(year), val(pmid), val(sex), val(input_build), val(population), val(delimiter), val(gwas_format), val(row_study_yaml), val(study_id), val(output_id), val(row_hash), val(source_type), val(raw_dir), val(standardise_dir), val(output_prefix)

    output:
    path '*.row_record.json', emit: records

    script:
    """
    python3 "${params.pipeline_root}/bin/stage_source.py" \\
      --gwas-location "${gwas_location}" \\
      --phenotype "${phenotype}" \\
      --ancestry "${ancestry}" \\
      --author "${author}" \\
      --year "${year}" \\
      --pmid "${pmid}" \\
      --sex "${sex}" \\
      --input-build "${input_build}" \\
      --population "${population}" \\
      --delimiter "${delimiter}" \\
      --format "${gwas_format}" \\
      --manifest-study-yaml "${row_study_yaml}" \\
      --study-id "${study_id}" \\
      --output-id "${output_id}" \\
      --row-hash "${row_hash}" \\
      --source-type "${source_type}" \\
      --raw-dir "${raw_dir}" \\
      --standardise-dir "${standardise_dir}" \\
      --output-prefix "${output_prefix}" \\
      --config-root "${params.config_root}" \\
      --study-yaml "${params.study_yaml}" \\
      --reference-root "${params.reference_root}" \\
      --dry-run "${params.dry_run}" \\
      --output stage_record.json

    record_id=\$(python3 -c 'import json, sys; print(json.load(open(sys.argv[1]))["output_id"])' stage_record.json)
    python3 "${params.pipeline_root}/bin/standardise_gwas.py" \\
      --stage-record stage_record.json \\
      --config-root "${params.config_root}" \\
      --study-yaml "${params.study_yaml}" \\
      --reference-root "${params.reference_root}" \\
      --force "${params.force}" \\
      --dry-run "${params.dry_run}" \\
      --output "\${record_id}.standardise_record.json"

    if [[ "${params.qcplot}" == "true" ]]; then
      python3 "${params.pipeline_root}/bin/qc_plot.py" \\
        --standardise-record "\${record_id}.standardise_record.json" \\
        --force "${params.force}" \\
        --dry-run "${params.dry_run}" \\
        --output "\${record_id}.row_record.json"
    else
      cp "\${record_id}.standardise_record.json" "\${record_id}.row_record.json"
    fi
    """
}

process UPDATE_STATE {
    tag "update-state"
    publishDir "${params.log_root}", mode: 'copy', overwrite: true

    input:
    path records

    output:
    path 'state_update_summary.json', emit: summary

    script:
    """
    python3 "${params.pipeline_root}/bin/update_state.py" \\
      --state-file "${params.work_root}/.standardise_state.tsv" \\
      --records ${records.join(' ')} \\
      --qcplot "${params.qcplot}" \\
      --output state_update_summary.json
    """
}
