.. _installation:

Installation
============

Prerequisites
-------------

* Python 3.11+
* A UNIX-like environment (e.g. MacOS, WSL, Ubuntu)
* A recent version of PostgreSQL (ideally at least 11+)

Install Cool-Seq-Tool
---------------------

Install Cool-Seq-Tool from `PyPI <https://pypi.org/project/cool-seq-tool/>`_:

    pip install cool-seq-tool

.. _dependency-groups:

.. note::

   Cool-Seq-Tool provides extra dependency groups for development and testing purposes. Most users won't need to install them.

   * ``dev`` includes packages for linting and performing other kinds of quality checks
   * ``tests`` includes packages for running tests
   * ``docs`` includes packages for writing and building documentation

Set up SeqRepo
--------------

Cool-Seq-Tool requires access to `SeqRepo <https://github.com/biocommons/biocommons.seqrepo>`_ data. In general, we recommend the following for local setup:

.. long-term, it would be best to move this over to seqrepo to avoid duplication

.. code-block::

   pip install seqrepo
   export SEQREPO_VERSION=2024-12-20
   sudo mkdir /usr/local/share/seqrepo
   sudo chown $USER /usr/local/share/seqrepo
   seqrepo pull -i $SEQREPO_VERSION

.. note::

   Our lab typically uses the latest SeqRepo release, which is ``2024-12-20`` as of this commit. To check for the presence of newer snapshots, use the ``seqrepo list-remote-instances`` CLI command.

While this should no longer occur with the latest SeqRepo release, some users in the past have reported experiencing the following error:

.. code-block::

   PermissionError: [Error 13] Permission denied: '/usr/local/share/seqrepo/2024-12-20._fkuefgd' -> '/usr/local/share/seqrepo/2024-12-20'

Try moving data manually with ``sudo``:

.. code-block::

   sudo mv /usr/local/share/seqrepo/$SEQREPO_VERSION.* /usr/local/share/seqrepo/$SEQREPO_VERSION

See `mirroring documentation <https://github.com/biocommons/biocommons.seqrepo/blob/main/docs/mirror.rst>`_ on the SeqRepo GitHub repo for instructions and additional troubleshooting.

Set up using Docker
-------------------

Cool-Seq-Tool's dependencies can be installed using a Docker container. We only provide guidance on setting up external dependencies using Docker.

.. important::

   This section assumes you have a local
   `SeqRepo <https://github.com/biocommons/biocommons.seqrepo>`_
   installed at ``/usr/local/share/seqrepo/2024-12-20``.
   If you have it installed elsewhere, please add a
   ``SEQREPO_ROOT_DIR`` environment variable in ``.env.shared``.
   See the `SeqRepo setup section <#set-up-seqrepo>`_ for additional information.

   You must download `uta_20241220.pgd.gz` from
   <https://dl.biocommons.org/uta/> using a web browser and
   move it to the root of the repository.

   If you're using Docker Desktop, you must go to
   **Settings → Resources → File sharing** and add
   ``/usr/local/share/seqrepo`` under the *Virtual file shares*
   section. Otherwise, you will get the following error::

      OSError: Unable to open SeqRepo directory /usr/local/share/seqrepo/2024-12-20

To build, (re)create, and start containers:

.. code-block:: shell

   docker volume create uta_vol
   docker compose up

.. tip::

   If you want a clean slate, run ``docker compose down -v`` to remove
   containers and volumes, then run
   ``docker compose up --build`` to rebuild and start fresh containers.

In Docker Desktop, you should see the following for a successful setup:

.. figure:: ../../docker-desktop-container.png
   :alt: Docker Desktop Container
   :align: center


Check data availability
-----------------------

The :py:meth:`check_status <cool_seq_tool.resources.status.check_status>` method is provided to check that data dependencies like SeqRepo and UTA are available and configured correctly:

.. code-block:: pycon

   >>> from cool_seq_tool.resources.status import check_status
   >>> import asyncio
   >>> asyncio.run(check_status())
   ResourceStatus(uta=True, seqrepo=True, transcript_mappings=True, mane_summary=True, lrg_refseqgene=True, liftover=True)
