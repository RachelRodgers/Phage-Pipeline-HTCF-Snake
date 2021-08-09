#----- tBLASTx Contig Dictionary against the Gut Phage DB (GPD) -----#

# query_GPD.smk

#----- Snakemake Set Up -----#
configfile: "./config/phage_pipeline_config.yaml"

# Tools
MMSEQS = config["Tools"]["MMSeqs2"]

# Databases
GPD = config["DBs"]["GPD"]

rule mmseqs_GPD:
        """
        Run mmseqs2 easy-search of contig dictionary against GPD
        """
        input:
                queryFile = os.path.join("results", "contig_dictionary", "assembly_1kb.fasta"),
                targetFile = {GPD}
        output:
                tmpDir = directory(temp(os.path.join("results", "mmseqs_GPD_results_tmp"))),
                alignmentFile = os.path.join("results", "mmseqs_GPD_results", "mmseqs_GPD_results.m8")
        threads: 16
        resources:
                mem_mb = 64000
        shell:
                """
                module load {MMSEQS}
                mmseqs easy-search \
                        {input.queryFile} \
                        {input.targetFile} \
                        {output.alignmentFile} \
                        {output.tmpDir} \
                        --search-type 2 \
                        --threads {threads}
                """
