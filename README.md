# **C**ommon **O**perations **O**n **L**ots-of **Seq**uences Tool

The **cool-seq-tool** provides:

  - Transcript alignment data from the [UTA](https://github.com/biocommons/uta) database
  - Fast access to sequence data using [SeqRepo](https://github.com/biocommons/biocommons.seqrepo)
  - Liftover between assemblies (GRCh38 <--> GRCh37) from [PyLiftover](https://github.com/konstantint/pyliftover)
  - Lifting over to preferred [MANE](https://www.ncbi.nlm.nih.gov/refseq/MANE/) compatible transcript. See [here](docs/TranscriptSelectionPriority.md) for more information.

## Installation

### pip

```commandline
pip install cool-seq-tool
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
pipenv lock && pipenv sync
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

To connect to the UTA database, you can use the default url (`postgresql://uta_admin@localhost:5433/uta/uta_20210129`). If you use the default url, you must either set the password using environment variable `UTA_PASSWORD` or setting the parameter `db_pwd` in the UTA class.

If you do not wish to use the default, you must set the environment variable `UTA_DB_URL` which has the format of `driver://user:pass@host:port/database/schema`.

### Data Downloads

#### SeqRepo
`cool-seq-tool` relies on [seqrepo](https://github.com/biocommons/biocommons.seqrepo), which you must download yourself.

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

#### transcript_mappings.tsv
`cool-seq-tool` uses [Ensembl BioMart](http://www.ensembl.org/biomart/martview) to retrieve `cool_seq_tool/data/transcript_mappings.tsv`. We currently use `Human Genes (GRCh38.p13)` for the dataset and the following attributes we use are: Gene stable ID, Gene stable ID version, Transcript stable ID, Transcript stable ID version, Protein stable ID, Protein stable ID version, RefSeq match transcript (MANE Select), Gene name.

![image](biomart.png)

#### LRG_RefSeqGene

`cool-seq-tool` fetches the latest version of `LRG_RefSeqGene`. This file is found can be found [here](https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/RefSeqGene).

#### MANE Summary Data

`cool-seq-tool` fetches the latest version of `MANE.GRCh38.*.summary.txt.gz`. This file is found can be found [here](https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/).

## Starting the UTA Tools Service Locally

To start the service, run the following:

```commandline
uvicorn cool_seq_tool.main:app --reload
```

Next, view the FastAPI on your local machine: http://127.0.0.1:8000/cool_seq_tool

## Init coding style tests
Code style is managed by [flake8](https://github.com/PyCQA/flake8) and checked prior to commit.

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
