# Snakef vile

configfile: "config.yaml"

# Rule all: track only the main .h5ad files
rule all:
    input:
        config["OBJ"],
        config["FILTERED_OBJ"],
        config["ANALYSED_OBJ"],
        config["ANNOTATED_OBJ"],
        expand("figures/{prefix}_nhood_enrichment.png", prefix=config["PREFIX"]),
        config["SPATIAL_DE_OBJ"]

rule get_data:
    output:
        config["OBJ"]
    shell:
        "python src/get_visium.py --output {output}"

# Filter AnnData
rule filter:
    input:
        config["OBJ"]
    output:
        config["FILTERED_OBJ"]
    params: config['PREFIX']
    shell:
        """
        python src/filter.py --input {input} --output {output} \
        --prefix {params} --min_genes {config[MIN_GENES]} \
        --min_cells {config[MIN_CELLS]} --mt {config[MT]}
        """

# Analyse AnnData
rule analyse:
    input:
        config["FILTERED_OBJ"]
    output:
        config["ANALYSED_OBJ"]
    params: 
           config['PREFIX']
    shell:
        "python src/analyse.py --input {input} --output {output} --prefix {params}"


# Annotate AnnData
rule annotate:
    input:
        config["ANALYSED_OBJ"]
    output:
        config["ANNOTATED_OBJ"]
    params:
        annot_file=config["ANNOT_FILE"],
        prefix = config['PREFIX']
    shell:
        """
        python src/annotate.py \
        --input {input} \
        --output {output} \
        --markers {params.annot_file} --prefix {params.prefix}
        """

# Generate plots
rule plot:
    input:
        config["ANNOTATED_OBJ"]
    output:
        "figures/{prefix}_nhood_enrichment.png",
    params:
        prefix =config['PREFIX']
    shell:
        "python src/plots.py --input {input} --prefix {params.prefix}"



# Spatial differential Expression per tissue 

rule spatial_DE:
    input: 
       config["ANNOTATED_OBJ"]
    output: 
       config["SPATIAL_DE_OBJ"] 
    params: 
      config['PREFIX']
    shell: 
      """ 
      python src/spatialDE.py --input {input} --output {output} --prefix {params}  
      """ 
