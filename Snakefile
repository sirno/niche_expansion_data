import yaml

from misosoup.library.common import merge_yaml

configfile: "config.yaml"


rule all:
    input:
        expand("out/tolerance_{tolerance}.yaml", tolerance=config["tolerances"])

rule misosoup:
    input:
        config["media"]
    output:
        "out/final/{tolerance}/{strain}/{carbon_source}.yaml"
    params:
        models=config["models"],
        params=lambda wildcards: config["parameters"]["all"],
        cache="out/cache/{tolerance}/{strain}/{carbon_source}.yaml.cache"
    resources:
        time="05-00",
        mem_per_cpu="4G",
        cpus_per_task=2,
    shell:
        "misosoup {params.models} "
        "--media {input} "
        "--media-select \"{wildcards.carbon_source}\" "
        "--output \"{output}\" "
        "--strain \"{wildcards.strain}\" "
        "--tolerance {wildcards.tolerance} "
        "{params.params} "


rule gather:
    input:
        [
            f"out/final/{{tolerance}}/{strain}/{carbon_source}.yaml"
            for carbon_source in config["carbon_sources"]
            for strain in config["strains"]
        ],
    output:
        "out/tolerance_{tolerance}.yaml",
    resources:
        mem_mb_per_cpu="4096"
    run:
        solutions = merge_yaml(input)
        with open(output[0], "w") as fd:
            yaml.dump(solutions, fd, Dumper=yaml.CSafeDumper)

