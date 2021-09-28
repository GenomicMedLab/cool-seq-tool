# uta_tools

Service for querying the [UTA](https://github.com/biocommons/uta) database

## Development

### Installation

Clone the repo:

```commandline
git clone https://github.com/cancervariants/uta_tools
cd uta_tools
```

[Install Pipenv](https://pipenv-fork.readthedocs.io/en/latest/#install-pipenv-today) if necessary.

Install backend dependencies and enter Pipenv environment:

```commandline
pipenv shell
pipenv lock && pipenv sync
```

## UTA Database Installation

uta_tools uses intalls local UTA database. For other ways to install, visit [biocommons.uta](https://github.com/biocommons/uta).

### Local Installation

_The following commands will likely need modification appropriate for the installation environment._
1. Install [PostgreSQL](https://www.postgresql.org/)
2. Create user and database.

    ```
    $ createuser -U postgres uta_admin
    $ createuser -U postgres anonymous
    $ createdb -U postgres -O uta_admin uta
    ```

3. To install locally, from the _variation/data_ directory:
```
export UTA_VERSION=uta_20210129.pgd.gz
curl -O http://dl.biocommons.org/uta/$UTA_VERSION
gzip -cdq ${UTA_VERSION} | grep -v "^REFRESH MATERIALIZED VIEW" | psql -h localhost -U uta_admin --echo-errors --single-transaction -v ON_ERROR_STOP=1 -d uta -p 5433
```

### Connecting to the database

To connect to the UTA database, you can use the default url (`postgresql://uta_admin@localhost:5433/uta/uta_20210129`). If you use the default url, you must either set the password using environment variable `UTA_PASSWORD` or setting the parameter `db_pwd` in the UTA class.

If you do not wish to use the default, you must set the environment variable `UTA_DB_URL` which has the format of `driver://user:pass@host/database/schema`.

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