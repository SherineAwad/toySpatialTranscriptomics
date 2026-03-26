# Snakef vile

configfile: "config.yaml"

# Rule all: track only the main .h5ad files
rule all:
    input:
        config["OBJ"],
        config["FILTERED_OBJ"],
        config["ANALYSED_OBJ"],
        expand("figures/{prefix}_spatial.png", prefix=config["PREFIX"]),
        expand("figures/{prefix}_neighbors.png", prefix=config["PREFIX"]),
        expand("figures/{prefix}_top_moranI.png", prefix=config["PREFIX"])

# Generate dummy data
rule generate_dummy:
    output:
        config["OBJ"]
    shell:
        "python src/dummy.py --output {output}"

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
        "python src/analyse.py --input {input} --output {output} --prefix {config[PREFIX]}"


# Generate plots
rule plot:
    input:
        config["ANALYSED_OBJ"]
    output:
        "figures/{prefix}_spatial.png",
        "figures/{prefix}_neighbors.png",
        "figures/{prefix}_top_moranI.png",
    params:
        markers=config["MAKER_GENES"]
    shell:
        "python src/plots.py --input {input} --prefix {wildcards.prefix} --markers {params.markers}"
