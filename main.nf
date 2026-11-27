nextflow.enable.dsl = 2


process SIMULATE {
    tag "run_${run_index}"
    publishDir { "${params.outdir}/run_${run_index}/simulation" }, mode: 'copy'

    input:
    tuple val(run_index), path(run_config), val(log_file), val(update_log_file)

    output:
    tuple val(run_index), path(log_file), emit: log
    tuple val(run_index), path(update_log_file), emit: updates

    script:
    """
    sim.py --config "${run_config}" --run-index ${run_index}
    """
}


process LINEAGE_MEAN_BIRTHRATE {
    tag "run_${run_index}"
    publishDir { "${params.outdir}/run_${run_index}/rates" }, mode: 'copy'

    input:
    tuple val(run_index), path(sim_log)

    output:
    tuple val(run_index), path('survivor_lineage_mean_birthrate.tsv'), emit: rates

    script:
    """
    lineage_mean_birthrate.py "${sim_log}" --output survivor_lineage_mean_birthrate.tsv
    """
}


process GROUND_TRUTH_BIRTHRATE {
    tag "run_${run_index}"
    publishDir { "${params.outdir}/run_${run_index}/rates" }, mode: 'copy'

    input:
    tuple val(run_index), path(update_log)

    output:
    tuple val(run_index), path('ground_truth_birthrate.tsv'), emit: rates

    script:
    """
    calculate_birthrate.py "${update_log}" > ground_truth_birthrate.tsv
    """
}


process SUBSAMPLE {
    tag "run_${run_index}"
    publishDir { "${params.outdir}/run_${run_index}/subsample" }, mode: 'copy'

    input:
    tuple val(run_index), path(sim_log)

    output:
    tuple val(run_index), path('subsampled.tsv'), emit: ancestry
    tuple val(run_index), path('subsampled_survivors.tsv'), emit: survivors

    script:
    def seed_option = params.seed != null ? "--seed ${params.seed}" : ''
    """
    subsample.py --input "${sim_log}" --output subsampled.tsv \
        --endnodes subsampled_survivors.tsv --N ${params.n} ${seed_option}
    """
}


process LOG_TO_FASTA {
    tag "run_${run_index}"
    publishDir { "${params.outdir}/run_${run_index}/subsample" }, mode: 'copy'

    input:
    tuple val(run_index), path(survivors)

    output:
    tuple val(run_index), path('subsampled.fasta'), emit: fasta

    script:
    """
    log_to_fasta.py "${survivors}" --output subsampled.fasta
    """
}


process CREATE_BEAST_XML {
    tag "run_${run_index}"
    publishDir { "${params.outdir}/run_${run_index}/xml" }, mode: 'copy'

    input:
    tuple val(run_index), path(fasta), path(sim_log)
    path template

    output:
    tuple val(run_index), path('sdevo.xml'), emit: sdevo
    tuple val(run_index), path('strict_clock.xml'), emit: strict_clock

    script:
    """
    create_beast_xml.py "${fasta}" --simulation-log "${sim_log}" --template "${template}" \
        --clock-model sdevo --output sdevo.xml \
        --chain-length ${params.chain_length} --log-every ${params.log_every}
    create_beast_xml.py "${fasta}" --simulation-log "${sim_log}" --template "${template}" \
        --clock-model strict --output strict_clock.xml \
        --chain-length ${params.chain_length} --log-every ${params.log_every}
    """
}


process VISUALIZE {
    tag "run_${run_index}"
    publishDir { "${params.outdir}/run_${run_index}/visualizations" }, mode: 'copy'

    input:
    tuple val(run_index), path(rates), path(ancestry), path(sim_log)

    output:
    tuple val(run_index), path('*.{png,svg,nwk}'), emit: plots

    script:
    """
    export MPLCONFIGDIR="\$PWD/.matplotlib"
    visualize_sim.py --simulation-log "${sim_log}" --ancestry-log "${ancestry}" \
        --rates "${rates}" --output-dir .
    """
}


process RUN_BEAST {
    tag "run_${run_index}:${model}"
    publishDir { "${params.outdir}/run_${run_index}/beast/${model}" }, mode: 'copy'
    cpus 1

    input:
    tuple val(run_index), val(model), path(xml)

    output:
    tuple val(run_index), val(model), path("${model}.log"), path("${model}.typed.node.trees"), emit: inference
    tuple val(run_index), val(model), path("${model}.trees"), path("${model}.typed.trees"), emit: trees
    path '*.state', optional: true

    script:
    """
    "${params.beast}" -java -overwrite -threads ${task.cpus} \
        -seed ${params.beast_seed + run_index} "${xml}"
    """
}


workflow {
    def run_config = file(params.runs, checkIfExists: true)
    def runs = new groovy.json.JsonSlurper().parse(run_config.toFile())

    SIMULATE(
        channel.fromList(runs.withIndex().collect { parameters, index ->
            def update_log = parameters.update_log_file ?: parameters.log_file.replaceFirst(/\.[^\/.]+$/, '') + '_updates.tsv'
            tuple(index, run_config, parameters.log_file, update_log)
        })
    )
    LINEAGE_MEAN_BIRTHRATE(SIMULATE.out.log)
    GROUND_TRUTH_BIRTHRATE(SIMULATE.out.updates)
    SUBSAMPLE(SIMULATE.out.log)
    LOG_TO_FASTA(SUBSAMPLE.out.survivors)
    CREATE_BEAST_XML(
        LOG_TO_FASTA.out.fasta.join(SIMULATE.out.log),
        file(params.template, checkIfExists: true)
    )
    VISUALIZE(
        LINEAGE_MEAN_BIRTHRATE.out.rates.join(SUBSAMPLE.out.ancestry).join(SIMULATE.out.log)
    )
    RUN_BEAST(
        CREATE_BEAST_XML.out.sdevo.map { index, xml -> tuple(index, 'sdevo', xml) }
            .mix(CREATE_BEAST_XML.out.strict_clock.map { index, xml -> tuple(index, 'strict_clock', xml) })
    )
}
