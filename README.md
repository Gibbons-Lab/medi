<img src="https://github.com/Gibbons-Lab/medi/blob/main/.github/logo.png" width="50%">

# Metagenomic Estimation of Dietary Intake (MEDI)

This repo contains the MEDI nextflow pipeline(s) for recovering food abundances and
nutrient composition from metagenomic shotgun sequencing data.

It contains capabilities for the following individual functionalities:

1. Mapping items in FOODB to all currently available items in NCBI NIH databases and
   prioritizing hits
2. Downloading, consolidating, ANI distance calculation using minhasing, for all full and
   partial assemblies and annotation of  NCBI Taxonomy IDs for use with Kraken2.
3. Building the Kraken2 and BRACKEN hashes/databases including decoy sequences.
4. Perform quantification of food DNA and per-portion nutrient composition for metagenomic
   samples, starting from raw FASTQ files.

## Installation

> [!WARNING]
> MEDI is built and tested on Linux only which should be the most common configuration
> due to the large resource requirements (500GB RAM).

You will need a working miniforge or miniconda to start. YOu can follow the [installation
instructions here](https://github.com/conda-forge/miniforge?tab=readme-ov-file#install). After
create an environment with the included conda environment file.

Start by cloning the repository and cd-ing into it.

```bash
git clone https://github.com/gibbons-lab/medi
cd medi
```

```bash
conda env create -n medi -f medi.yml
```

After that activate the environment.

```bash
conda activate medi
```

### Compile the report generator

> [!TIP]
> This will only be necessary if the provided binary does not work.
> You can test this by running `./bin/kraken2-report` if this returns
> "malformed taxonomy file" you are good. If you get errors about the
> ELF class or a missing libc version you will need to recompile.

Kraken2 does not support generating reports on filtered output files by default.
We provide a pre-compiled report generator from a [Kraken2 fork](https://github.com/daydream-boost/kraken2).
If this does not work you can compile it using the provided Makefile:

```bash
conda activate medi
make report
```

This will compile the report generator for your platform and replace the binary.

---

And you are done. If you are running this on a HPC cluster or a cloud provider, you
moght need to [adjust your nextflow settings](https://www.nextflow.io/docs/latest/config.html#config-scopes) for your setup.

All pipelines support a `--threads` parameter that defines the maximum number of threads
to use for any single process.

## Calling MEDI steps

If you want to run steps in the install directory no additional steps are needed.

After setting up the conda environment there are two ways to make MEDI steps run
in arbitrary data directories.

First you can create a variable holding the MEDI installation and use this.

```bash
export MEDI=/my/install/location

cd /my/data/directory
nextflow run $MEDI/quant.nf --db=/my/medi_db
```

Second, you can create symlinks to the step and and the `bin` directory.

```bash
cd /my/data/directory
ln -s /my/install/location/bin
ln -s /my/install/location/quant.nf

nextflow run $MEDI/quant.nf --db=/my/medi_db
```

Both will work the same.

# MEDI steps

Here are the full steps to build the database and run it on your data.

> [!NOTE]
> If you are asking yourself why we don't just provide the built database for download
> please [see the comments here](docs/db_download.md).

## (1-2) Matching and downloading

```bash
nextflow run -resume database.nf
```

This will boostrap the database from nothing, downloading all required files and performing
the matching against the current versions of NCBI Genbank and Nucleotide. You can speed up the querying by obtaining and NCBI API key and adding it to your `.Rprofile` with

```text
options(reutil.api.key="XXXXXX")
```

Please also the the [troubleshooting guide](docs/db_download.md) in case you encounter issues.

## (3) Building the hashes

After running the previous step continue with

```bash
nextflow run build_kraken.nf --max_db_size=500
```

Here `--max_db_size` denotes the maximum size of the database in GB. The default
will use no reduction but you can set this to a lower level which will create a
smaller but less accurate hash. Note that for good performance you will need as much
RAM as what you choose here.

Note that this step of the pipeline will not work with the `-resume` option. The `add_*` need
to finish completely or the pipeline needs to be restarted from the beginning. Should this
work and the later steps crash, you can trigger just the hash building using the
`--rebuild` option which will rebuild the database but not attempt to add sequences again.

> [!WARNING]
> Do not run this step with the `-resume` Nextflow option as this will result in a
> corrupted database. In case, all adding sequences worked and you only want to run the
> build step again, use the `--rebuild=true` option.

## (4) Quantification for MGS samples

For your own sequencing data create a directory and use either of the setups described above.

Then create a `data` directory and within that a `raw` folder containing your unprocessed
demultiplexed FASTQ files. So it should look like:

```text
|- quant.nf // optional
|- bin/     // optional
|- data/
   |-- raw/
```

After that you can run MEDI with

```bash
nextflow run quant.nf -resume --db=/path/to/medi_db
```

Where `/path/to/medi_db/` should be the output directory from step (3). Usually `medi/data/medi_db`.

> [!TIP]
> For HPC clusters it can be beneficial to tune the `--batchsize` parameter to obtain
> massive performance benefits from memory-mapping. [See here](docs/concepts.md) for
> more info. For single machines we recommend `--batchsize 1`

### Output files

Those are all found in the output directory, usually `data/` if not specified otherwise.

| filename           | description                                          |
| --------------     | -----------------------------------------------------|
| food_abundance.csv | Food abundances for each sample in the data set. Samples with only NA entries had no detected food reads. |
| food_content.csv   | Estimates for nutrient and compound abundances normalized to a representative portion of 100g. |
| multiqc_report.html | Quality control metrics for the sequencing data. |

### Bonus: Getting microbial abundances

By default MEDI also reports bacteria, archaea, and host abundances with the same
filtering and decoy steps as for the food quantification. Those can be found in the
output directory (usually `data/`) and are names `[RANK]_counts.csv` where `[RANK]` denotes
the taxonomic rank (S - species, G - genus, D - domain).

## TODO

- [] see if we can provide a reduced DB for download
- [x] make execution more flexible
- [x] add resource limits for individual steps for HPC clusters
