nextflow.enable.dsl = 2


process SIMULATE {
    tag "run_${run_index}"
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(run_index), path(run_config)

    output:
    tuple val(run_index), path("run_${run_index}"), emit: runs

    script:
    """
    mkdir "run_${run_index}"
    cd "run_${run_index}"
    sim.py \
        --config "../${run_config}" \
        --run-index ${run_index}
    """
}


workflow {
    def run_config = file(params.runs, checkIfExists: true)
    def runs = new groovy.json.JsonSlurper().parse(run_config.toFile())

    SIMULATE(
        channel.fromList(runs.indices).map { index -> tuple(index, run_config) }
    )
}
