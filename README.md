<img src=".github/logo.png" width="50%">

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

And you are done. If you are running this on a HPC cluster or a cloud provider, you
moght need to [adjust your nextflow settings](https://www.nextflow.io/docs/latest/config.html#config-scopes) for your setup.

All pipelines support a `--max_threads` parameter that defines the maximum number of threads
to use for any single process.

## (1-2) Matching and downloading

```bash
nextflow run -resume database.nf
```

This will boostrap the database from nothing, downloading all required files and performing
the matching against the current versions of NCBI Genbank and Nucleotide. You can speed up the
querying by obtaining and NCBI API key and adding it to your `.Rprofile` with

```text
options(reutil.api.key="XXXXXX")
```

## (3) Building the hashes

After running the previous step continue with

```bash
nextflow run build_kraken.nf --max_db_size=500
```

Here `--max_db_size` denotes the maximum size of the database in GB. The default
will use no reduction but you can set this to a lower level which will create a
smaller but less accurate hash. Note that for good performance you will need more
RAM than what you choose here.

Note that this step of the pipeline will not work with the `-resume` option. The `add_*` need
to finish completely or the pipeline needs to be restarted from the beginning. Should this
work and the later steps crash, you can trigger just the hash building using the
`--rebuild` option which will rebuild the database but not attempt to add sequences again.

## (4) Quantification for MGS samples

For your own sequencing data create a directory and either copy or link the medi
pipeline there. You will need at least the `quant.nf` file and the `scripts` folder in there.
Than create a `data` directory and within that a `raw` folder containing your unprocessed
demultiplexed FASTQ files. So it should look like:

```text
|- quant.nf
|- scripts/
|- data/
   |-- raw/
```

After that you can run MEDI with

```bash
nextflow run quant.nf -resume --db=/path/to/medi_db
```

Where `/path/to/medi_db/` should be the output directory from step (3). Usually `medi/data/medi_db`.

## TODO

- see if we can provide a reduced DB for download
- make execution more flexible
- add resource limits for individual steps for Grid clusters
