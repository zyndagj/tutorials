GPU Containers on Slurm
========================

This tutorial introduces software containers, how to build them, and how to run them on Slurm clusters using apptainer.
This is not meant to teach container mastery, but expose you to some best practices with containers on HPC systems.

Last updated: |today|

**Objectives**

* Running GPU containers
* Building and testing your first GPU container
* Best practices for building python-based containers
* Developing python scripts inside a running container on BCP
* Running multi-node containers

**Requirements**

* Container build system - Apptainer build **OR** Docker CLI Pre-authenticated with a container registry such as:

   * `Docker Hub <hub.docker.com>`_
   * `nvcr.io <https://org.ngc.nvidia.com/setup/api-key>`_

* Container runtime - Slurm GPU cluster with `apptainer <https://apptainer.org/>`_ **OR** `enroot <https://github.com/NVIDIA/enroot/>`_

.. _prep:

**Prepare your environment**

.. tabs::

    .. group-tab:: Apptainer

        .. literalinclude:: assets/prep_apptainer.sh
            :caption: :download:`prep_apptainer.sh <assets/prep_apptainer.sh>`

    .. group-tab:: Docker
    
        .. literalinclude:: assets/prep_docker.sh
            :caption: :download:`prep_docker.sh <assets/prep_docker.sh>`

Introduction to containers
---------------------------

Containers are a common method for building, distributing, and running applications, web services, development, and more.
Containers have gained popularity because they package up an application and all dependencies, provide isolation from the host environment, and allow for a consistent deployment across platforms.
While you may have also heard of virtual machines (VMs), containers are separate and rely on namespace virtualization without hardware emulation, so there are no performance losses.

The portability and reproducibility, without sacrificing performance, make containers ideal for scientific applications. Whole environments can be saved to ensure a published tool can always be used over time.

`Docker <https://www.docker.com/>`_ is the most common container runtime, but there are many that can consume Open Container Initiative (OCI) images .
This tutorial will be focusing on building containers with `apptainer <https://apptainer.org/>`_ or docker, and then running containers in a shared HPC environment with apptainer or `enroot <https://github.com/NVIDIA/enroot>`_.

Quickstart: Running your first GPU container
----------------------------------------------

Think of this as a quickstart to running GPU containers on HPC systems.

Pulling the container
###############################

Containers should only be run on a compute node inside a job, so I recommend starting a job on a GPU node.
If you're not already on one, take a look at how to `prepare your environment <prep_>`_ or the ``srun`` command below if you're using ``enroot``.

To run a container on HPC systems, we first need to pull the layers from the container registry and then convert them into a single image.

.. tabs::

    .. group-tab:: Apptainer

        .. literalinclude:: assets/pull_apptainer.sh
            :caption: :download:`pull_apptainer.sh <assets/pull_apptainer.sh>`

    .. group-tab:: Enroot
    
        .. literalinclude:: assets/pull_enroot.sh
            :caption: :download:`pull_enroot.sh <assets/pull_enroot.sh>`

.. note::

    Keep this job running for the rest of this tutorial

Examining the CUDA environment
###############################

To start off, take a look at the CUDA environment outside of the container by running the `NVIDIA System Management Interface <https://docs.nvidia.com/deploy/nvidia-smi/index.html>`_ program (``nvidia-smi``).

.. code-block:: shell

    $ nvidia-smi

    Tue Dec  3 18:46:40 2024       
    +-----------------------------------------------------------------------------------------+
    | NVIDIA-SMI 565.57.01              Driver Version: 565.57.01      CUDA Version: 12.7     |
    |-----------------------------------------+------------------------+----------------------+
    | GPU  Name                 Persistence-M | Bus-Id          Disp.A | Volatile Uncorr. ECC |
    | Fan  Temp   Perf          Pwr:Usage/Cap |           Memory-Usage | GPU-Util  Compute M. |
    |                                         |                        |               MIG M. |
    |=========================================+========================+======================|
    |   0  NVIDIA L40S                    Off |   00000000:41:00.0 Off |                    0 |
    | N/A   36C    P8             35W /  350W |       1MiB /  46068MiB |      0%      Default |
    |                                         |                        |                  N/A |
    +-----------------------------------------+------------------------+----------------------+
                                                                                            
    +-----------------------------------------------------------------------------------------+
    | Processes:                                                                              |
    |  GPU   GI   CI        PID   Type   Process name                              GPU Memory |
    |        ID   ID                                                               Usage      |
    |=========================================================================================|
    |  No running processes found                                                             |
    +-----------------------------------------------------------------------------------------+

Running ``nvidia-smi`` is the easiest way to see if there's a GPU on your system and what driver it is running.
In addition to using it on your host, it works inside GPU-capable containers.

.. tabs::

    .. group-tab:: Apptainer

        .. literalinclude:: assets/nvidia-smi_apptainer.sh
            :caption: :download:`nvidia-smi_apptainer.sh <assets/nvidia-smi_apptainer.sh>`

    .. group-tab:: Enroot
    
        .. literalinclude:: assets/nvidia-smi_enroot.sh
            :caption: :download:`nvidia-smi_enroot.sh <assets/nvidia-smi_enroot.sh>`

If you're running ``aptainer``, you'll notice that the CUDA version doesn't change with the ``--nv`` flag.
This will change if the ``--nvccli`` option (nvidia container cli) is enabled on your system.

Optional Exercises
##########################

* What happens if you exclude the ``--nv`` flag with ``apptainer``?
* What happens if you run on the container on a system without a GPU?

Building and testing your first GPU container
---------------------------------------------

In this section, we'll be building the `nbody sample benchmark <https://github.com/NVIDIA/cuda-samples/tree/master/Samples/5_Domain_Specific/nbody>`_ from https://github.com/NVIDIA/cuda-samples.

The nbody benchmark demonstrates efficient all-pairs simulation of a gravitational n-body simulation in CUDA and provides a GFLOP/s metric at the end.
While this GFLOP/s metric is not meant for true performance comparisons, this sample code supports multiple GPUs and is relatively easy to build.

Containers are built using recipe files like Docker's `Dockerfile <https://docs.docker.com/reference/dockerfile/>`_ or Apptainer's `Definition file <https://apptainer.org/docs/user/main/definition_files.html#>`_, which are essentially scripts for provisioning a linux environment.

Choosing a starting container
#############################

The first step to building any container is choosing an image to start from.
This starting image is often a clean OS like this `ubuntu image <https://hub.docker.com/_/ubuntu>`_, from which you can add any necessary dependencies to build/run your software. Alternatively, you can start from an image that already contains software they're pre-installed.

We're going to be building and running a GPU application, so I recommend starting from NVIDIA's `CUDA container <https://catalog.ngc.nvidia.com/orgs/nvidia/containers/cuda>`_ on NGC.
NGC is NVIDIA's container registry, where NVIDIA software, SDKs, and models are published in container format.
Not only are these meant to make your development easier, they're also serve as a common environment for NVIDIA to reproduce and troubleshoot any issues you might encounter through `enterprise support <https://enterprise-support.nvidia.com/>`_ with `NVAIE <https://www.nvidia.com/en-us/data-center/products/ai-enterprise/>`_.

Looking at the tags tab, you'll see many different containers.
To help you understand the naming convention, containers usually have a ``<project>/<name>:<tag>`` format.
If you browse through the available containers, you'll see that each container is named cuda, but tags have some common elements along with a CUDA version prefix:

* ``base``: Includes the CUDA runtime (cudart)
* ``runtime``: base + CUDA math libraries, and NCCL
* ``devel``: runtime + headers, development tools for compiling CUDA applications
* ``cudnn-``: (prefix) any of the above + cuDNN libraries

There are a ton of options, so here are some recommendations on choosing a container:

* Latest CUDA version (unless a specific one is needed)
   
   * Newer libraries work on older drivers

* ``base`` for simple CUDA applications
* ``devel`` for multi-staged builds
* Choose an OS with a package manager you're familiar with

.. note::

    We'll cover multi-staged builds in container optimization

In this case, we're going to start from the ``nvcr.io/nvidia/cuda:12.4.1-devel-ubuntu22.04`` container that we already pulled and cached during the quickstart.

Installing dependencies and building
####################################

Just like when trying to run an application, identifying and installing compatible dependencies is the hardest part of container development.
If you look at the `dependencies for nbody <https://github.com/NVIDIA/cuda-samples/tree/master/Samples/5_Domain_Specific/nbody#dependencies-needed-to-buildrun>`_, X11 and GL are required to build and run.
On an ubuntu system (notice container tag), we can install the development headers and libraries along with ``curl`` using:

.. code-block:: shell

    apt-get update && apt-get install -y --no-install-recommends \
		freeglut3-dev libgl1-mesa-dev libglu1-mesa-dev curl

These commands won't work for non-root users because they modify the host system.
If you're figuring out how to build a container, you can prototype commands in an interactive container:

.. tabs::

    .. group-tab:: Apptainer

        .. literalinclude:: assets/interactive_build_apptainer.sh
            :language: shell
            :caption: :download:`interactive_build_apptainer.sh <assets/interactive_build_apptainer.sh>`

    .. group-tab:: Docker
    
        .. literalinclude:: assets/interactive_build_docker.sh
            :caption: :download:`interactive_build_docker.sh <assets/interactive_build_docker.sh>`

Once the dependencies are installed, you can download, build, and install the nbody application with the following commands:

.. code-block:: shell

    # Grab the sample code
    curl -sL https://github.com/NVIDIA/cuda-samples/archive/refs/tags/v12.4.1.tar.gz -o v12.4.1.tar.gz

    # Unpack the tarball
    tar -xzf v12.4.1.tar.gz

    # Build the nbody executable
    cd cuda-samples-12.4.1/Samples/5_Domain_Specific/nbody \
	    && make && mv nbody /usr/local/bin

Wrapping it all up and building the container
##############################################

Your desired starting container and installation commands can be wrapped up into a single file.
Apptainer uses `Definition files <https://apptainer.org/docs/user/main/definition_files.html>`_ and Docker uses `Dockerfiles <https://docs.docker.com/reference/dockerfile/>`_.

``exit`` your interactive container instance and ``wget`` your corresponding build file.

.. tabs::

    .. group-tab:: Apptainer

        .. literalinclude:: assets/Definition.nbody
            :caption: :download:`Definition.nbody <assets/Definition.nbody>`

    .. group-tab:: Docker
    
        .. literalinclude:: assets/Dockerfile.nbody
            :caption: :download:`Dockerfile.nbody <assets/Dockerfile.nbody>`

.. note::

    You can either download this file directly or copy and paste into your favorite text editor

You can then build a container named **nbody** from your build script as follows:

.. tabs::

    .. group-tab:: Apptainer

        .. literalinclude:: assets/build_nbody_apptainer.sh
            :caption: :download:`build_nbody_apptainer.sh <assets/build_nbody_apptainer.sh>`

    .. group-tab:: Docker
    
        .. literalinclude:: assets/build_nbody_docker.sh
            :caption: :download:`build_nbody_docker.sh <assets/build_nbody_docker.sh>`

This is a relatively large image, so not only does it take up a lot of space on the filesystem, but it also would take a while to upload to a remote registry for sharing or archive.
Lets instead figure out how to make our final image more space efficient.

Making your container more space efficient
###########################################

We can make this much smaller using the following techniques:

#. Use a multi-staged build - Building in one container and copying build binaries to a runtime container

   * `Docker multi-staged build documentation <https://docs.docker.com/build/building/multi-stage/>`_
   * `Apptainer multi-staged build documentation <https://apptainer.org/docs/user/main/definition_files.html#multi-stage-builds>`_

#. Only install runtime libraries in the final container

   #. Using the base container instead of devel
   #. Not installing ``*-devel`` packages from apt

#. Copy the finished binary instead of the full source repo

.. tabs::

    .. group-tab:: Apptainer

        .. literalinclude:: assets/Definition.nbody-efficient
            :caption: :download:`Definition.nbody-efficient <assets/Definition.nbody-efficient>`

    .. group-tab:: Docker

        .. literalinclude:: assets/Dockerfile.nbody-efficient
            :caption: :download:`Dockerfile.nbody-efficient <assets/Dockerfile.nbody-efficient>`

Make sure to change the name or tag of the container when building it.

.. tabs::

    .. group-tab:: Apptainer

        .. literalinclude:: assets/build_nbody-efficient_apptainer.sh
            :caption: :download:`build_nbody-efficient_apptainer.sh <assets/build_nbody-efficient_apptainer.sh>`

    .. group-tab:: Docker
    
        .. literalinclude:: assets/build_nbody-efficient_docker.sh
            :caption: :download:`build_nbody-efficient_docker.sh <assets/build_nbody-efficient_docker.sh>`

Once again, lets look at the final size of the containers we built.

.. code-block:: shell

    $ ls -lh nbody*sif

    -rwxr-xr-x 1 greg.zynda greg.zynda.grp 147M Dec  3 20:34 nbody-efficient.sif
    -rwxr-xr-x 1 greg.zynda greg.zynda.grp 4.2G Dec  3 08:51 nbody.sif

In the case of these apptainer ``.sif`` images built by ``apptainer``, you'll notice that the efficient build is much smaller: 147MB vs 4.2GB!
Not only will this take up less space on your filesystem, but it's also easier to archive with a publication.

Running the nbody sample benchmark
###################################

You should already be inside a job with an allocated GPU, so you can run the benchmark with the following:

.. tabs::

    .. group-tab:: Apptainer

        .. literalinclude:: assets/run_nbody_apptainer.sh
            :caption: :download:`run_nbody_apptainer.sh <assets/run_nbody_apptainer.sh>`

    .. group-tab:: Enroot
    
        .. literalinclude:: assets/run_nbody_enroot.sh
            :caption: :download:`run_nbody_enroot.sh <run_nbody_enroot.sh>`

When your job is done, you should see output similar to the following:

.. code-block:: shell

    > Windowed mode
    > Simulation data stored in video memory
    > Single precision floating point simulation
    > 1 Devices used for simulation
    GPU Device 0: "Ada" with compute capability 8.9

    > Compute 8.9 CUDA device: [NVIDIA L40S]
    Warning: "number of bodies" specified 2000000 is not a multiple of 256.
    Rounding up to the nearest multiple: 2000128.
    2000128 bodies, total time for 10 iterations: 21772.984 ms
    = 1837.374 billion interactions per second
    = 36747.484 single-precision GFLOP/s at 20 flops per interaction

.. note::

    These performance results will change based on the GPU type your were allocated.

Optional Exercises
##########################

* Looking at the help text, try using a different number of GPUs (requires new job)
* Try increasing the number of bodies in the simulation
* Try using double precision

Best practices for building python-based containers
---------------------------------------------------

One of the most common things I encounter when folks use containers with pre-existing python packages and libraries, is accidentally replacing or overwriting them with ``conda`` or ``pip``.
NVIDIA's NGC containers have patched version of PyTorch and supporting libraries that shouldn't be altered if you're looking for optimal and verified performance.

This section will focus on how to install python packages in a way that will prevent changes to the pre-installed packages.

To illustrate this, try installing pytorch from the base ``pytorch:24.03-py3`` container.

.. tabs::

    .. group-tab:: Apptainer

        .. literalinclude:: assets/pip-install_apptainer.sh
            :language: shell
            :caption: :download:`pip-install_apptainer.sh <assets/pip-install_apptainer.sh>`

    .. group-tab:: Docker
    
        .. literalinclude:: assets/pip-install_docker.sh
            :caption: :download:`pip-install_docker.sh <assets/pip-install_docker.sh>`

You'll notice that installing these packages changes the toch package and installs a bunch of CUDA libraries even though both already exist.
As you learned with our efficient builds, this greatly increases the size of the container layers while also potentially breaking any applications linked against these libaries and the "known working state".

Lets exit this container create a fresh overlay.

.. code-block:: shell

    # Be sure to exit your interactive container session
    exit

Luckily, you can lock the versions by creating a package `constraints file <https://pip.pypa.io/en/stable/user_guide/#constraints-files>`_, which has the same format as a requirements file.

.. tabs::

    .. group-tab:: Apptainer

        .. literalinclude:: assets/pip-constraints_apptainer.sh
            :language: shell
            :caption: :download:`pip-constraints_apptainer.sh <assets/pip-constraints_apptainer.sh>`

    .. group-tab:: Docker
    
        .. literalinclude:: assets/pip-constraints_docker.sh
            :caption: :download:`pip-constraints_docker.sh <assets/pip-constraints_docker.sh>`

This install should now fail because the pre-built torchaudio wheels can't be installed with the NVIDIA patched versions of torch.

.. note::

    If you actually want to install torchaudio into the Pytorch NGC container, take a look at `this recipe <https://github.com/NVIDIA/NeMo/blob/main/scripts/installers/install_torchaudio_latest.sh#L97>`_.

Lets practice using this constraint method by building a new container with the `PyTorch Lightning <https://lightning.ai/>`_ framework starting FROM the ``pytorch:24.03-py3`` container.

.. tabs::

    .. group-tab:: Apptainer

        .. literalinclude:: assets/Definition.lightning
            :language: shell
            :caption: :download:`Definition.lightning <assets/Definition.lightning>`

    .. group-tab:: Docker
    
        .. literalinclude:: assets/Dockerfile.lightning
            :caption: :download:`Dockerfile.lightning <assets/Dockerfile.lightning>`

After download the corresponding build script, the container can be built with the following commands.

.. tabs::

    .. group-tab:: Apptainer

        .. literalinclude:: assets/build_lightning_apptainer.sh
            :language: shell
            :caption: :download:`build_lightning_apptainer.sh <assets/build_lightning_apptainer.sh>`

    .. group-tab:: Docker
    
        .. literalinclude:: assets/build_lightning_docker.sh
            :caption: :download:`build_lightning_docker.sh <assets/build_lightning_docker.sh>`

Unlike the ``torchaudio`` install, this went fine, and no existing packages changed.
If a package or its dependencies require a different version of PyTorch, you can either change the container version based on the `NVIDIA support matrix <https://docs.nvidia.com/deeplearning/frameworks/support-matrix/index.html>`_ to match the required version or determine if the `package's dependencies can be relaxed <https://pip.pypa.io/en/stable/topics/dependency-resolution/#loosen-your-top-level-requirements>`_ to match the package version in the container.

Optional Exercises
####################

* Try installing another python package

Developing python scripts inside a running container
-----------------------------------------------------------

`Containers <https://www.docker.com/resources/what-container/>`_ are meant to be static, reproducible checkpoints for your code that can always be started in the same way.
This makes them ideal for porting software to different systems, reproducing results, archiving software, and more.
However, since containers *shouldn't* change once they're built (because that would break reproducibility), developing software in them is not always intuitive.

If you try to incorporate all your code in the container and rebuilding as it evolves, this can get tedious - especially if you're pushing and pulling these containers between a registry.
Instead, I recommend making a container with most or all of your dependencies, and mounting your code into the container at runtime.

To explore these concepts, lets launch an interactive environment with our lightning container.

.. tabs::

    .. group-tab:: Apptainer

        .. literalinclude:: assets/lightning_interactive-apptainer.sh
            :language: shell
            :caption: :download:`lightning_interactive-apptainer.sh <assets/lightning_interactive-apptainer.sh>`

    .. group-tab:: Enroot
    
        .. literalinclude:: assets/lightning_interactive-enroot.sh
            :caption: :download:`lightning_interactive-enroot.sh <assets/lightning_interactive-enroot.sh>`

First, lets open another terminal to the cluster.
That could be another tmux pane or a whole new terminal connection from your local system.
Once you have that open, lets look around in the running container.

.. list-table:: Exploring Environment
    :widths: 40 30 30
    :header-rows: 1
    
    * - 
      - Container shell
      - Second shell
    * - Who are you running as?
      - whoami
      - whoami
    * - Where are you running from?
      - pwd
      - cd $MYDATA/containers
    * - Do files match?
      - ls -lh
      - ls -lh
    * - Do changes propogate?
      - echo "hello" > container.txt
      - cat container.txt
    * - What else is in the container by default?
      - ls -lh $HOME; ls -lh /tmp
      - ls -lh $HOME; ls -lh /tmp
    * - What if you create a file somewhere else?
      - touch /workspace/test
      - ls /workspace
    * - Should you be able to create files?
      - ls -lhd /workspace
      - ls -lhd /workspace

.. note::

    ``$MYDATA/containers`` was available in the container because the container mounts our current working directory at runtime.
    If you need additional locations available in the container, you can make them available with (similar to Docker's ``-v``):

    * `apptainer - bind paths <https://apptainer.org/docs/user/main/bind_paths_and_mounts.html>`_ (-B)
    * `enroot - mount <https://github.com/NVIDIA/enroot/blob/master/doc/cmd/start.md>`_ (-m)

Running external scripts
##############################

As you experienced when trying to create a test file in ``/workspace``, which is open for writing, you discovered that the container has a read-only filesystem.
This means, that you can't make any changes without an overlay.
This might be tedious for prototyping, but it's good if you're sharing a container with colleagues on a project, or if you just want to make sure you can't accidentally make changes.

First, download :download:`python_dev.tar.gz <assets/python_dev.tar.gz>` to your current working directory with ``wget`` (may need to use ``--no-check-certificate``).
After downloading, unpack the tarball with ``tar``.

.. code-block:: shell

    # Unpack
    tar -xzf python_dev.tar.gz

    cd python_dev

    ls *

The script ``self_contained.py`` doesn't require any extra python modules other than PyTorch, which exists in the container, and can be run directly.
Try running it.

.. code-block:: shell

    python self_contained.py

Not only can you run scripts from inside the container, you can interact with them outside the container too.
If you have your other terminal still open, find these files and open the script in your favorite editor.
Not only can you open the files, you can edit them too - all while beign able to run them inside a container.

Developing packages from inside a container
#############################################

If you're developing a whole package that needs to be updated, you either have to rely on relative imports or install the package.
Relative imports often work, but may not depending on the complexity of the package.
In our example python code, there's a ``pt_bench`` python module that gets loaded and used by ``bench.py``.

.. code-block:: shell

    # Prints where pt_bench was loaded from
    python bench.py

    # Change directories
    cd ..

    # Copy bench.py to break relative imports
    cp python_dev/bench.py .
    python bench.py

You can see that it's easy to go wrong with relative imports, so I often recommend fully installing the package.

We already know that the container can't be modified.
Luckily, python can install packages in a user directory, that defaults to ``$HOME/.local``, using the ``--user`` flag.

.. code-block:: shell

    # Install pt_bench using our constraint file
    pip install -c /root/base_constraints.txt --user python_dev/

    # Try running bench.py again
    python bench.py

You should see that ``pt_bench`` is being loaded from ``$HOME/.local``, which is where user packages are installed.
While this works, this location is universally shared by all python packages, which will lead to collisions between containers.
I recommend launching the container with ``-c``, which will not mount any external locations, and ``-B`` to mount the current working directory.
Since many things require a valid ``$HOME`` for writing files, ``apptainer`` creates a temporary filesystem (tmpfs) for ``/home``.
You'll be able to make changes, like installing a small package, and it won't affect the container or bleed into other python environments.

First, lets clean our environment

.. code-block:: shell

    # remove pt_bench
    pip uninstall -y pt_bench

    # Exit the container
    exit

and then relaunch.

.. tabs::

    .. group-tab:: Apptainer

        .. literalinclude:: assets/lightning_interactive_nohome-apptainer.sh
            :language: shell
            :caption: :download:`lightning_interactive_nohome-apptainer.sh <assets/lightning_interactive_nohome-apptainer.sh>`

    .. group-tab:: Enroot
    
        .. literalinclude:: assets/lightning_interactive_nohome-enroot.sh
            :caption: :download:`lightning_interactive_nohome-enroot.sh <assets/lightning_interactive_nohome-enroot.sh>`

.. code-block:: shell

    # Make sure home is empty
    ls $HOME

    # Change to container directory
    cd $MYDATA/containers

    # Try running bench.py
    python bench.py
    # Install wasn't found

    # Do a local install in $HOME tmpfs
    pip install -c /root/base_constraints.txt --user python_dev/

    # Run bench.py
    python bench.py

Lastly, if you're making changes to the package, you can do an `editable install <https://pip.pypa.io/en/stable/topics/local-project-installs/#editable-installs>`_ with `-e`.
This means that when the package is installed, it's really just linked to it's current location instead of copying files.

.. code-block:: shell

    # remove pt_bench
    pip uninstall -y pt_bench

    # Editable install (-e)
    pip install -c /root/base_constraints.txt --user -e python_dev/

    # Make a change to a package file
    echo "print('New Change')" >> python_dev/pt_bench/__init__.py

    # Run bench, and see if change works
    python bench.py

When you exit the container, make sure the pt_bench package no longer exists.

.. code-block:: shell

    # Exit the container
    exit

    # Make sure pt_bench doesn't exist
    find $HOME/.local/ | grep pt_bench

If you're done exploring the container, feel free to exit the job in preparation for the next section.

.. code-block:: shell

    # Exit the job
    exit

Running multi-node containers
--------------------------------

Multi-node, or distributed, computing is a model of computation that takes the tasks from an algorithm that can be run independently and executes them across multiple computers.
While it's easy to spawn threads and processes on system, distributed applications need to be launched across all nodes and told how to communicate with eachother.

Multi-node MPI NCCL Test
#########################

PyTorch containers from NGC ship with `NCCL tests <https://github.com/NVIDIA/nccl-tests>`_, which are useful for diagnosing MPI and bandwidth issues.
If I'm ever questioning the performance of the compute fabric between GPUs, this is the first thing I run.

These can be run as single-line jobs using ``srun`` to handle the allocation and process spawning.

.. tabs::

    .. group-tab:: Apptainer

        .. literalinclude:: assets/nccl-apptainer.sh
            :language: shell
            :caption: :download:`nccl-apptainer.sh <assets/nccl-apptainer.sh>`

    .. group-tab:: Enroot
    
        .. literalinclude:: assets/nccl-enroot.sh
            :caption: :download:`nccl-enroot.sh <assets/nccl-enroot.sh>`


Multi-node PyTorch
######################

Using ``wget``, download :download:`pt_ddp_example.py <assets/pt_ddp_example.py>`, which is a simple script to demonstrate strong scaling using `PyTorch DDP <https://pytorch.org/tutorials/intermediate/ddp_tutorial.html>`_.
We'll be skipping over PyTorch specifics to focus on how to launch multi-node PyTorch containers with Slurm.
Download the following ``sbatch`` script as well.

.. tabs::

    .. group-tab:: Apptainer

        .. literalinclude:: assets/pt_ddp_example.sbatch
            :language: shell
            :caption: :download:`pt_ddp_example.sbatch <assets/pt_ddp_example.sbatch>`

    .. group-tab:: Enroot
    
PyTorch needs the following variables set for multi-node runs:

* MASTER_ADDR - Address of the main node
* MASTER_PORT - Port of connect on
* WORLD_SIZE - Total number of workers/GPUs

While ``srun`` launches the initial process on each node, it calls ``torchrun``, which spawns additional processes based on the argument ``--nproc_per_node``.
Think of ``torchrun`` as a helper script that handles a lot of the global and local rank logic.

Optional variables:

* LOGLEVEL - pytorch log level
* NCCL_DEBUG - NCCL log level

Submit the script with ``sbatch``, which will generate a ``.out`` file with a number corresponding to the job with all output text.
You'll see that this runs a training job on 4 GPUs in total, distributed across 2 nodes.
If you increase the resources allocated by the ``SBATCH`` arguements, training will scale as well.

Multi-node Pytorch Lightning
#############################

This is the same task as the `Multi-node PyTorch`_ script, just adapted to PyTorch Lightning.
You'll notice that code is clearner because PyTorch Lightning does it's best to simplify common training tasks, including multi-GPU and multi-node trainging.

Download both the training script :download:`ptl_ddp_example.py <assets/ptl_ddp_example.py>` and the ``sbatch`` script below.

.. tabs::

    .. group-tab:: Apptainer

        .. literalinclude:: assets/ptl_ddp_example.sbatch
            :language: shell
            :caption: :download:`ptl_ddp_example.sbatch <assets/ptl_ddp_example.sbatch>`

    .. group-tab:: Enroot

With PyTorch Lightning, both the ``MASTER_ADDR`` and ``MASTER_PORT`` need to be set, but also the ``NODE_RANK``, which is the 0-based index of the node the process is on.
In this example, it's being set in a bash shell, with the ``$`` escaped so it's substituted after being launched on each node by ``srun``.
When it's running, you'll see that Lightning has nice logging about the process pool at the start, and produces nice output during the training progress.

Optional Exercises
##########################

* Try using more GPUs to see how the number of steps run by each GPU scales.
* Try comparing training and NCCL performance on different types of nodes.

Next Steps
----------

Apptainer/Singularity is a well known container runtime in the world of HPC, but NVIDIA `recommends using enroot <https://github.com/NVIDIA/enroot/issues/25>`_ as a container runtime for several reasons.
Enroot doesn't have a build functionality, but can consume OCI images built by Docker or buildah and can be combined with `pyxis <https://github.com/NVIDIA/pyxis>`_ for Slurm support.
I also highly recommend checking out Docker for building containers due to the size of the community and support availability.

NVIDIA Containers:

* `NGC Container Catalog <https://catalog.ngc.nvidia.com/containers>`_
* `Containers for DL Frameworks <https://docs.nvidia.com/deeplearning/frameworks/user-guide/index.html>`_
* `HPC with Containers DLI <https://learn.nvidia.com/courses/course-detail?course_id=course-v1:DLI+L-AC-25+V1>`_

Container workshops/tutorials:

* `Containers@TACC <https://containers-at-tacc.readthedocs.io/en/latest/>`_
* `Getting started with Docker <https://docs.docker.com/get-started/>`_

.. |br| raw:: html

    <br>