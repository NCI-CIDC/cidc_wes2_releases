rule download_vep_plugins:
    output:
        directory(Path(PREDIR) / "resources/vep/plugins")
    params:
        release=config["vep"]["release"]
    wrapper:
        "v3.7.0/bio/vep/plugins"
