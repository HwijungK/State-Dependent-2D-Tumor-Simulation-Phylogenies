nextflow.enable.dsl = 2


process SIMULATE {
    tag "run_${run_index}"
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(run_index), path(run_config)
    path simulation_script

    output:
    tuple val(run_index), path("run_${run_index}"), emit: runs

    script:
    """
    mkdir "run_${run_index}"
    cd "run_${run_index}"
    "${params.python}" "../${simulation_script}" \
        --config "../${run_config}" \
        --run-index ${run_index}
    """
}


workflow {
    def run_config = file(params.runs, checkIfExists: true)
    def runs = new groovy.json.JsonSlurper().parse(run_config.toFile())
    if (!(runs instanceof List) || runs.isEmpty()) {
        error 'The runs file must contain a nonempty JSON array of parameter objects.'
    }

    def run_inputs = channel.fromList(runs.withIndex().collect { parameters, index ->
        if (!(parameters instanceof Map)) {
            error "Run ${index} must be a JSON object of simulation parameters."
        }
        tuple(index, run_config)
    })

    SIMULATE(run_inputs, file("${projectDir}/bin/sim.py", checkIfExists: true))
}
