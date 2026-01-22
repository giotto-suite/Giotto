giotto_env <- basilisk::BasiliskEnvironment(envname="giotto_env",
                              pkgname="Giotto",
                              packages=c("python==3.11",
                                         "zarr==3.1.5",
                                         "anndata==0.12.6")
)
