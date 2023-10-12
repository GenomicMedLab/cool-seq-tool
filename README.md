# **C**ommon **O**perations **O**n **L**ots-of **Seq**uences Tool

The **cool-seq-tool** provides:

  - Transcript alignment data from the [UTA](https://github.com/biocommons/uta) database
  - Fast access to sequence data using [SeqRepo](https://github.com/biocommons/biocommons.seqrepo)
  - Liftover between assemblies (GRCh38 <--> GRCh37) from [PyLiftover](https://github.com/konstantint/pyliftover)
  - Lifting over to preferred [MANE](https://www.ncbi.nlm.nih.gov/refseq/MANE/) compatible transcript. See [here](docs/TranscriptSelectionPriority.md) for more information.

## Installation

### pip

```commandline
pip install cool-seq-tool[dev,tests]
```

### Development

Clone the repo:

```commandline
git clone https://github.com/GenomicMedLab/cool-seq-tool
cd cool_seq_tool
```

[Install Pipenv](https://pipenv-fork.readthedocs.io/en/latest/#install-pipenv-today) if necessary.

Install backend dependencies and enter Pipenv environment:

```commandline
pipenv shell
pipenv update
pipenv install --dev
```

### UTA Database Installation

`cool-seq-tool` uses intalls local UTA database. For other ways to install, visit [biocommons.uta](https://github.com/biocommons/uta).

#### Local Installation

_The following commands will likely need modification appropriate for the installation environment._
1. Install [PostgreSQL](https://www.postgresql.org/)
2. Create user and database.

    ```
    $ createuser -U postgres uta_admin
    $ createuser -U postgres anonymous
    $ createdb -U postgres -O uta_admin uta
    ```

3. To install locally, from the _cool_seq_tool/data_ directory:
```
export UTA_VERSION=uta_20210129.pgd.gz
curl -O http://dl.biocommons.org/uta/$UTA_VERSION
gzip -cdq ${UTA_VERSION} | grep -v "^REFRESH MATERIALIZED VIEW" | psql -h localhost -U uta_admin --echo-errors --single-transaction -v ON_ERROR_STOP=1 -d uta -p 5433
```

##### UTA Installation Issues
If you have trouble installing UTA, you can visit [these two READMEs](https://github.com/ga4gh/vrs-python/tree/main/docs/setup_help).

#### Connecting to the database

To connect to the UTA database, you can use the default url (`postgresql://uta_admin:uta@localhost:5433/uta/uta_20210129`).

If you do not wish to use the default, you must set the environment variable `UTA_DB_URL` which has the format of `driver://user:password@host:port/database/schema`.

### Data Downloads

#### SeqRepo
`cool-seq-tool` relies on [seqrepo](https://github.com/biocommons/biocommons.seqrepo), which you must download yourself.

Use the `SEQREPO_ROOT_DIR` environment variable to set the path of an already existing SeqRepo directory. The default is `/usr/local/share/seqrepo/latest`.

From the _root_ directory:
```
pip install seqrepo
sudo mkdir /usr/local/share/seqrepo
sudo chown $USER /usr/local/share/seqrepo
seqrepo pull -i 2021-01-29  # Replace with latest version using `seqrepo list-remote-instances` if outdated
```

If you get an error similar to the one below:
```
PermissionError: [Error 13] Permission denied: '/usr/local/share/seqrepo/2021-01-29._fkuefgd' -> '/usr/local/share/seqrepo/2021-01-29'
```

You will want to do the following:\
(*Might not be ._fkuefgd, so replace with your error message path*)
```console
sudo mv /usr/local/share/seqrepo/2021-01-29._fkuefgd /usr/local/share/seqrepo/2021-01-29
exit
```

#### LRG_RefSeqGene

`cool-seq-tool` fetches the latest version of `LRG_RefSeqGene` if the environment variable `LRG_REFSEQGENE_PATH` is not set. When `LRG_REFSEQGENE_PATH` is set, `cool-seq-tool` will look at this path and expect the LRG_RefSeqGene file. This file is found can be found [here](https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/RefSeqGene).

#### MANE Summary Data

`cool-seq-tool` fetches the latest version of `MANE.GRCh38.*.summary.txt.gz` if the environment variable `MANE_SUMMARY_PATH` is not set. When `MANE_SUMMARY_PATH` is set, `cool-seq-tool` will look at this path and expect the MANE Summary Data file. This file is found can be found [here](https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/).

#### transcript_mapping.tsv
`cool-seq-tool` is packaged with transcript mapping data acquired from [Ensembl BioMart](http://www.ensembl.org/biomart/martview). If the environment variable `TRANSCRIPT_MAPPINGS_PATH` is not set, `cool-seq-tool` will use the built-in file. When `TRANSCRIPT_MAPPINGS_PATH` is set, `cool_seq_tool` will look at this path and expect to find the transcript mapping TSV file.

To acquire this data manually from the [BioMart](https://www.ensembl.org/biomart/martview), select the `Human Genes (GRCh38.p13)` dataset and choose the following attributes:

* Gene stable ID
* Gene stable ID version
* Transcript stable ID
* Transcript stable ID version
* Protein stable ID
* Protein stable ID version
* RefSeq match transcript (MANE Select)
* Gene name

![image](biomart.png)

## Starting the UTA Tools Service Locally

To start the service, run the following:

```commandline
uvicorn cool_seq_tool.api:app --reload
```

Next, view the FastAPI on your local machine: http://127.0.0.1:8000/cool_seq_tool

## Init coding style tests

Code style is managed by [Ruff](https://github.com/astral-sh/ruff) and [Black](https://github.com/psf/black), and should be checked prior to commit.

We use [pre-commit](https://pre-commit.com/#usage) to run conformance tests.

This ensures:

* Check code style
* Check for added large files
* Detect AWS Credentials
* Detect Private Key

Before first commit run:

```
pre-commit install
```

## Testing
From the _root_ directory of the repository:
```
pytest
```
