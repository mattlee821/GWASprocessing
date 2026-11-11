nextflow.enable.dsl = 2

// One manifest row = one STANDARDISE_ROW task. Rows travel as JSON files written by prepare_manifest.py,
// so manifest text (phenotype names, URLs) is never interpolated into a shell command.
workflow {
    PREPARE_MANIFEST()
    row_records = STANDARDISE_ROW(PREPARE_MANIFEST.out.rows.flatten())
    UPDATE_STATE(row_records.records.collect().ifEmpty([]))
}

process PREPARE_MANIFEST {
    tag "prepare-manifest"
    publishDir "${params.log_root}", mode: 'copy', overwrite: true, pattern: 'run_manifest*'

    output:
    path 'run_manifest.tsv', emit: manifest
    path 'run_manifest_summary.json', emit: summary
    path 'rows/*.json', emit: rows, optional: true

    script:
    """
    python3 "${params.pipeline_root}/bin/prepare_manifest.py" \\
      --manifest "${params.manifest}" \\
      --repo-root "${params.repo_root}" \\
      --work-root "${params.work_root}" \\
      --state-file "${params.work_root}/.standardise_state.tsv" \\
      --config-root "${params.config_root}" \\
      --reference-root "${params.reference_root}" \\
      --study-yaml "${params.study_yaml}" \\
      --analysis-qc "${params.analysis_qc}" \\
      --output run_manifest.tsv \\
      --summary run_manifest_summary.json \\
      --rows-dir rows \\
      --only-study "${params.only_study}" \\
      --force "${params.force}" \\
      --qcplot "${params.qcplot}" \\
      --row-limit "${params.row_limit}"
    """
}

process STANDARDISE_ROW {
    tag "${row_json.baseName}"

    input:
    path row_json

    output:
    path '*.row_record.json', emit: records

    script:
    """
    python3 "${params.pipeline_root}/bin/stage_source.py" \\
      --row-json "${row_json}" \\
      --config-root "${params.config_root}" \\
      --study-yaml "${params.study_yaml}" \\
      --reference-root "${params.reference_root}" \\
      --dry-run "${params.dry_run}" \\
      --output stage_record.json

    python3 "${params.pipeline_root}/bin/standardise_gwas.py" \\
      --stage-record stage_record.json \\
      --config-root "${params.config_root}" \\
      --study-yaml "${params.study_yaml}" \\
      --reference-root "${params.reference_root}" \\
      --force "${params.force}" \\
      --dry-run "${params.dry_run}" \\
      --analysis-qc "${params.analysis_qc}" \\
      --workers ${task.cpus} \\
      --output standardise_record.json

    if [[ "${params.qcplot}" == "true" ]]; then
      python3 "${params.pipeline_root}/bin/qc_plot.py" \\
        --standardise-record standardise_record.json \\
        --force "${params.force}" \\
        --dry-run "${params.dry_run}" \\
        --output "${row_json.baseName}.row_record.json"
    else
      cp standardise_record.json "${row_json.baseName}.row_record.json"
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
      --qcplot "${params.qcplot}" \\
      --output state_update_summary.json \\
      --records ${records}
    """
}
