#----- tBLASTx Contig Dictionary against the Metagenomic Gut Virus Catalog (MGV) -----#

# query_MGV.smk

#----- Snakemake Set Up -----#
configfile: "./config/phage_pipeline_config.yaml"

# Tools
MMSEQS = config["Tools"]["MMSeqs2"]

# Databases
MGV = config["DBs"]["MGV"]

rule mmseqs_MGV:
        """
        Run mmseqs2 easy-search of contig dictionary against MGV
        """
        input:
                queryFile = os.path.join("results", "contig_dictionary", "assembly_1kb.fasta"),
                targetFile = {MGV}
        output:
                tmpDir = directory(temp(os.path.join("results", "mmseqs_MGV_results_tmp"))),
                alignmentFile = os.path.join("results", "mmseqs_MGV_results", "mmseqs_MGV_results.m8")
        threads: 16
        resources:
                mem_mb = 96000
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
