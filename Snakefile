# Snakefile

configfile: "config.yaml"

# Rule all: track only the main .h5ad files
rule all:
    input:
        config["OBJ"],
        config["FILTERED_OBJ"],
        config["ANALYSED_OBJ"],
        f"figures/{config['PREFIX']}_umap.png",
        f"figures/{config['PREFIX']}_spatial.png",
        f"figures/{config['PREFIX']}_neighbors.png",
        f"figures/{config['PREFIX']}_top_moranI.png",

# Generate dummy data
rule generate_dummy:
    output:
        config["OBJ"]
    shell:
        "python dummy.py --output {output}"

# Filter AnnData
rule filter:
    input:
        config["OBJ"]
    output:
        config["FILTERED_OBJ"]
    shell:
        """
        python filter.py --input {input} --output {output} \
        --prefix {config[PREFIX]} --min_genes {config[MIN_GENES]} \
        --min_cells {config[MIN_CELLS]} --mt {config[MT]}
        """

# Analyse AnnData
rule analyse:
    input:
        config["FILTERED_OBJ"]
    output:
        config["ANALYSED_OBJ"]
    shell:
        "python analyse.py --input {input} --output {output} --prefix {config[PREFIX]}"

# Generate plots
rule plot:
    input:
        config["ANALYSED_OBJ"]
    shell:
        "python plots.py --input {input} --prefix {config[PREFIX]} --markers {config[MAKER_GENES]}"
