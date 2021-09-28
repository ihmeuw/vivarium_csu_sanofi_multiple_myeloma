===============================
vivarium_csu_sanofi_multiple_myeloma
===============================

Research repository for the vivarium_csu_sanofi_multiple_myeloma project.

.. contents::
   :depth: 1

Installation
------------

You will need ``git``, ``git-lfs`` and ``conda`` to get this repository
and install all of its requirements.  You should follow the instructions for
your operating system at the following places:

- `git <https://git-scm.com/downloads>`_
- `git-lfs <https://git-lfs.github.com/>`_
- `conda <https://docs.conda.io/en/latest/miniconda.html>`_

Once you have all three installed, you should open up your normal shell
(if you're on Linux or macOS) or the ``git bash`` shell if you're on Windows.
You'll then make an environment, clone this repository, then install
all necessary requirements as follows::

  :~$ conda create --name=vivarium_csu_sanofi_multiple_myeloma python=3.8
  ...conda will download python and base dependencies...
  :~$ conda activate vivarium_csu_sanofi_multiple_myeloma
  (vivarium_csu_sanofi_multiple_myeloma) :~$ git clone https://github.com/ihmeuw/vivarium_csu_sanofi_multiple_myeloma.git
  ...git will copy the repository from github and place it in your current directory...
  (vivarium_csu_sanofi_multiple_myeloma) :~$ cd vivarium_csu_sanofi_multiple_myeloma
  (vivarium_csu_sanofi_multiple_myeloma) :~$ pip install -e .
  ...pip will install vivarium and other requirements...


Note the ``-e`` flag that follows pip install. This will install the python
package in-place, which is important for making the model specifications later.

Cloning the repository should take a fair bit of time as git must fetch
the data artifact associated with the demo (several GB of data) from the
large file system storage (``git-lfs``). **If your clone works quickly,
you are likely only retrieving the checksum file that github holds onto,
and your simulations will fail.** If you are only retrieving checksum
files you can explicitly pull the data by executing ``git-lfs pull``.

Vivarium uses the Hierarchical Data Format (HDF) as the backing storage
for the data artifacts that supply data to the simulation. You may not have
the needed libraries on your system to interact with these files, and this is
not something that can be specified and installed with the rest of the package's
dependencies via ``pip``. If you encounter HDF5-related errors, you should
install hdf tooling from within your environment like so::

  (vivarium_csu_sanofi_multiple_myeloma) :~$ conda install hdf5

The ``(vivarium_csu_sanofi_multiple_myeloma)`` that precedes your shell prompt will probably show
up by default, though it may not.  It's just a visual reminder that you
are installing and running things in an isolated programming environment
so it doesn't conflict with other source code and libraries on your
system.


Usage
-----

You'll find six directories inside the main
``src/vivarium_csu_sanofi_multiple_myeloma`` package directory:

- ``artifacts``

  This directory contains all input data used to run the simulations.
  You can open these files and examine the input data using the vivarium
  artifact tools.  A tutorial can be found at https://vivarium.readthedocs.io/en/latest/tutorials/artifact.html#reading-data

- ``components``

  This directory is for Python modules containing custom components for
  the vivarium_csu_sanofi_multiple_myeloma project.

- ``data``

  This directory contains  **small scale** external data for use in the simulation.

- ``model_specifications``

  This directory holds all model specifications and branch files
  associated with the project.

- ``results_processing``

  The post-processing and analysis code is stored in this directory.

- ``tools``

  This directory holds Python files used to run scripts used to prepare input data or process outputs.


Running Simulations
-------------------

With your conda environment active, you can then run simulations by, e.g.::

   (vivarium_csu_sanofi_multiple_myeloma) :~$ simulate run -v /<REPO_INSTALLATION_DIRECTORY>/vivarium_csu_sanofi_multiple_myeloma/src/vivarium_csu_sanofi_multiple_myeloma/model_specifications/united_states_of_america.yaml

The ``-v`` flag will log verbosely, so you will get log messages every time
step. For more ways to run simulations, see the tutorials at
https://vivarium.readthedocs.io/en/latest/tutorials/running_a_simulation/index.html
and https://vivarium.readthedocs.io/en/latest/tutorials/exploration.html
