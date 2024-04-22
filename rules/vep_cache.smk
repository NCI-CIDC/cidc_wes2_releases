rule get_vep_cache:
    output:
        directory(Path(PREDIR) / "resources/vep/cache")
    params:
        species=config["vep"]["species"],
        build=config["vep"]["build"],
        release=config["vep"]["release"]
    log:
        "logs/vep/cache.log",
    cache: "omit-software"  # save space and time with between workflow caching (see docs)
    wrapper:
        "v3.7.0/bio/vep/cache"