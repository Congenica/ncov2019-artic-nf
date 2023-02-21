process articGuppyPlex {
    tag { params.prefix + "-" + fastqDir }

    label 'largemem'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${params.prefix}*.fastq", mode: "copy"

    input:
        path fastqDir

    output:
        path "${params.prefix}*.fastq", emit: fastq

    script:
        """
        artic guppyplex \
        --min-length ${params.min_length} \
        --max-length ${params.max_length} \
        --prefix ${params.prefix} \
        --directory ${fastqDir}
        """
}

process articMinIONMedaka {
    tag { sampleName }

    label 'largecpu'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}*", mode: "copy"

    input:
        tuple path(fastq), path(scheme)

    output:
        path "${sampleName}*"

        tuple val(sampleName), path("${sampleName}.primertrimmed.rg.sorted.bam"), emit: ptrim
        tuple val(sampleName), path("${sampleName}.sorted.bam"), emit: mapped
        tuple val(sampleName), path("${sampleName}.consensus.fasta"), emit: consensus_fasta
        tuple val(sampleName), path("${sampleName}.pass.vcf.gz"), emit: vcf

    script:
        // Make an identifier from the fastq filename
        sampleName = fastq.getBaseName().replaceAll(~/\.fastq.*$/, '')

        // Configure artic minion pipeline
        minionRunConfigBuilder = []

        if ( params.normalise ) {
        minionRunConfigBuilder.add("--normalise ${params.normalise}")
        }

        if ( params.bwa ) {
        minionRunConfigBuilder.add("--bwa")
        } else {
        minionRunConfigBuilder.add("--minimap2")
        }

        minionFinalConfig = minionRunConfigBuilder.join(" ")

        """
        artic minion --medaka \
        ${minionFinalConfig} \
        --medaka-model ${params.medakaModel} \
        --threads ${task.cpus} \
        --scheme-directory ${params.schemeDir}/${params.scheme} \
        --read-file ${fastq} \
        --scheme-version SARS-CoV-2/${params.schemeVersion} \
        SARS-CoV-2 \
        ${sampleName}

        """
}

process articRemoveUnmappedReads {
    tag { sampleName }

    cpus 1

    input:
        tuple val(sampleName), path(bamfile)

    output:
        tuple val(sampleName), path("${sampleName}.mapped.sorted.bam")

    script:
        """
        samtools view -F4 -o ${sampleName}.mapped.sorted.bam ${bamfile}
        """
}
