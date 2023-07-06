Installation Instructions
=========================

Requirements
------------

* A UNIX-like environment (e.g. MacOS, WSL, Ubuntu) with superuser permissions
* Python 3.8+
* A recent version of PostgreSQL (ideally at least 11+), or Docker
* An available Java runtime (version 8.x or newer)

Package installation
--------------------

We recommend installing from PyPI:

.. code-block:: shell

   pip install cool-seq-tool

Data source setup
-----------------

``cool-seq-tool`` makes use of data drawn from several parent libraries. Some additional setup is required to make this data available upon installation.

Universal Transcript Archive (UTA)
++++++++++++++++++++++++++++++++++

``cool-seq-tool`` requires access to the `Universal Transcript Archive <https://github.com/biocommons/uta>`_ PostgreSQL database. See the UTA readme for all setup options and more detailed instructions.

By default, ``cool-seq-tool`` will attempt to connect to the UTA database at ``postgresql://uta_admin@localhost:5433/uta/uta_20210129``. If you use the default url, you must either set the password using environment variable ``UTA_PASSWORD`` or setting the parameter ``db_pwd`` in the UTA class.

If you do not wish to use the default, you must set the environment variable ``UTA_DB_URL`` with a `valid PostgreSQL connection string <https://www.postgresql.org/docs/current/libpq-connect.html#LIBPQ-CONNSTRING>`_.

Local UTA setup
_______________

The most accessible way to set up UTA is from a local PostgreSQL instance. Depending on the environment, something like the following should be sufficient:

1. Install `PostgreSQL <https://www.postgresql.org/>`_
2. Create user and database:

.. code-block:: shell

    createuser -U postgres uta_admin
    createuser -U postgres anonymous
    createdb -U postgres -O uta_admin uta

3. Acquire and load UTA data:

.. code-block:: shell

    export UTA_VERSION=uta_20210129.pgd.gz
    curl -O http://dl.biocommons.org/uta/$UTA_VERSION
    gzip -cdq ${UTA_VERSION} | grep -v "^REFRESH MATERIALIZED VIEW" | psql -h localhost -U uta_admin --echo-errors --single-transaction -v ON_ERROR_STOP=1 -d uta -p 5433

For troubleshooting, see additional documentation `here <https://github.com/ga4gh/vrs-python/tree/main/docs/setup_help>`_.

SeqRepo
+++++++

``cool-seq-tool`` relies on `SeqRepo <https://github.com/biocommons/biocommons.seqrepo>`_ for sequence metadata.

The following commands will acquire the latest SeqRepo data:

.. code-block:: shell

   sudo mkdir /usr/local/share/seqrepo
   sudo chown $USER /usr/local/share/seqrepo
   seqrepo pull -i 2021-01-29

Troubleshooting SeqRepo setup
_____________________________

If you get a PermissionError similar to the following:

.. code-block::

    PermissionError: [Error 13] Permission denied: '/usr/local/share/seqrepo/2021-01-29._XXXXXX' -> '/usr/local/share/seqrepo/2021-01-29'

You will want to do the following:

.. code-block:: shell

    # replace XXXXXXX with correct file path
    sudo mv /usr/local/share/seqrepo/2021-01-29._XXXXXXX /usr/local/share/seqrepo/2021-01-29

Static data files
+++++++++++++++++

``cool-seq-tool`` makes use of some additional gene and transcript data acquired from Ensembl and NCBI. On startup, it should automatically acquire the most recent available versions. See :ref:`static-files` for more information and for configuration options.
