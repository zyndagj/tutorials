Containers on Slurm
=====================

This tutorial introduces software containers, how to build them, and how to run them on Slurm clusters using apptainer.

*Last updated 12/2/2024*

**Objectives**

* Running GPU containers
* Building and testing your first GPU container
* Best practices for building python-based containers
* Developing python scripts inside a running container on BCP
* Running multi-node containers

**Requirements**

* Slurm GPU cluster with `apptainer <https://apptainer.org/>`_ or `enroot <https://github.com/NVIDIA/enroot/>`_
* Docker/Apptainer CLI Pre-authenticated with a container registry such as:

   * `Docker Hub <hub.docker.com>`_
   * `nvcr.io <https://org.ngc.nvidia.com/setup/api-key>`_
   * https://cloud.sylabs.io/library

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

Containers are common method for applications, web services, development, and more.
Containers have gained popularity because they package up an application and all dependencies, provide isolation from the host environment, and allow for consistent deployment across platforms.
While you may have also heard of virtual machines (VMs), containers are separate and rely on namespace virtualization without emulating any hardware, so there are no performance losses.

The portability, reproducibility, and performance make them ideal packing scientific applications, so whole environments can be saved to ensure a tool can always be used over time.

You've probably heard of `Docker <https://www.docker.com/>`_, but there are many different container runtimes that can consume Open Container Initiative (OCI) images . This tutorial will be focusing on building containers with `apptainer <https://apptainer.org/>`_ or docker, and then running containers in a shared HPC environment with apptainer or `enroot <https://github.com/NVIDIA/enroot>`_.

Running your first GPU container
--------------------------------

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

To start off, take a look at the CUDA environment outside of the container

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

* What happens if you run on the container on a system without a GPU?

Building and testing your first GPU container
---------------------------------------------

In this section, we'll be building the `nbody sample benchmark <https://github.com/NVIDIA/cuda-samples/tree/master/Samples/5_Domain_Specific/nbody>`_ from https://github.com/NVIDIA/cuda-samples.
The nbody benchmark demonstrates efficient all-pairs simulation of a gravitational n-body simulation in CUDA and provides a GFLOP/s metric at the end.
While this GFLOP/s metric is not meant for performance comparisons, this sample code supports multiple GPUs and is relatively easy to build.

Containers are built using recipe files like Docker's `Dockerfile <https://docs.docker.com/reference/dockerfile/>`_ or Apptainer's `Definition file <https://apptainer.org/docs/user/main/definition_files.html#>`_, which are essentially scripts for provisioning a linux environment.

Choosing a starting container
#############################

The first step to building any container is choosing an image to start from.
This is often a `base ubuntu image <https://hub.docker.com/_/ubuntu>`_, which is similar to a rootfs.
From there, you can ``apt-get`` any necessary dependencies and then add and build your software.

We're going to be building a GPU application from source, so I recommend starting from NVIDIA's `CUDA container <https://catalog.ngc.nvidia.com/orgs/nvidia/containers/cuda>`_ on NGC.
NGC is NVIDIA's container registry, where NVIDIA software, SDKs, and models are published in container format.

Container types:

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

Installing dependencies
############################

Just like when trying to run an application, identifying and installing compatible dependencies is the hardest part of container development.
If you look at the `dependencies for nbody <https://github.com/NVIDIA/cuda-samples/tree/master/Samples/5_Domain_Specific/nbody#dependencies-needed-to-buildrun>`_, X11 and GL are required to build and run.
On an ubuntu system, we can install the development headers and libraries along with ``curl`` using:

.. code-block:: shell

    apt-get update && apt-get install -y --no-install-recommends \
		freeglut3-dev libgl1-mesa-dev libglu1-mesa-dev curl

If you're figuring out how to build a container, you can prototype commands in an interactive container

.. tabs::

    .. group-tab:: Apptainer

        .. literalinclude:: assets/interactive_build_apptainer.sh
            :language: shell
            :caption: :download:`interactive_build_apptainer.sh <assets/interactive_build_apptainer.sh>`

    .. group-tab:: Docker
    
        .. literalinclude:: assets/interactive_build_docker.sh
            :caption: :download:`interactive_build_docker.sh <assets/interactive_build_docker.sh>`

Building and installing application
####################################

.. code-block:: shell

    # Grab the sample code
    curl -sL https://github.com/NVIDIA/cuda-samples/archive/refs/tags/v12.4.1.tar.gz -o v12.4.1.tar.gz

    # Unpack the tarball to /root
    tar -C /root -xzf /root/v12.4.1.tar.gz

    # Build the nbody executable
    cd /root/cuda-samples-12.4.1/Samples/5_Domain_Specific/nbody \
	    && make && mv nbody /usr/local/bin

Wrapping it all up
############################

.. tabs::

    .. group-tab:: Apptainer

        .. literalinclude:: assets/Definition.nbody
            :caption: :download:`Definition.nbody <assets/Definition.nbody>`

    .. group-tab:: Docker
    
        .. literalinclude:: assets/Dockerfile.nbody
            :caption: :download:`Dockerfile.nbody <assets/Dockerfile.nbody>`

.. note::

    You can either download this file directly or copy and paste into your favorite text editor

Building the container
###############################

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

#. Use a `multi-staged build <https://docs.docker.com/build/building/multi-stage/>`_
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

Make sure to change the name of the container when building it.

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

In the case of these apptainer ``.sif`` images, you'll notice that the efficient build is much smaller: 147MB vs 4.2GB!

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

Optional Exercises
##########################

* Looking at the help text, try using a different number of GPUs (requires new job)
* Try increasing the number of bodies in the simulation
* Try using double precision

Best practices for building python-based containers
---------------------------------------------------

Python packages can specify dependencies, and sometimes those dependencies can be strictly written where existing packages get changed.
Containers on NVIDIA's NGC contain patched versions of PyTorch and matching libraries that shouldn't be altered if you're looking for optimal performance.
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

You'll notice that installing these packages installs a bunch of CUDA libraries and downgrades torch.
This is not ideal since the pytorch containers on NGC already ship cuda libraries, so this not only breaks anything compiled against these libraries, but also needlessly increases the size of the container by replacing what already exists.

Luckily, you can lock the versions by creating a package `constraints file <https://pip.pypa.io/en/stable/user_guide/#constraints-files>`_, which is similar to a requirements file.

.. tabs::

    .. group-tab:: Apptainer

        .. literalinclude:: assets/pip-constraints_apptainer.sh
            :language: shell
            :caption: :download:`pip-constraints_apptainer.sh <assets/pip-constraints_apptainer.sh>`

    .. group-tab:: Docker
    
        .. literalinclude:: assets/pip-constraints_docker.sh
            :caption: :download:`pip-constraints_docker.sh <assets/pip-constraints_docker.sh>`

.. note::

    This should now fail because the pre-built torchaudio wheels can't be installed with NVIDIA patched versions.
    If you actually want to install torchaudio into the Pytorch NGC container, take a look at `this recipe <https://github.com/NVIDIA/NeMo/blob/main/scripts/installers/install_torchaudio_latest.sh#L97>`_.

Lets practice this constraint method by building a new container with the `PyTorch Lightning <https://lightning.ai/>`_ framework starting FROM the ``pytorch:24.03-py3`` container.

.. tabs::

    .. group-tab:: Apptainer

        .. literalinclude:: assets/Definition.lightning
            :language: shell
            :caption: :download:`Definition.lightning <assets/Definition.lightning>`

    .. group-tab:: Docker
    
        .. literalinclude:: assets/Dockerfile.lightning
            :caption: :download:`Dockerfile.lightning <assets/Dockerfile.lightning>`

And then build the container with the following commands

.. tabs::

    .. group-tab:: Apptainer

        .. literalinclude:: assets/build_lightning_apptainer.sh
            :language: shell
            :caption: :download:`build_lightning_apptainer.sh <assets/build_lightning_apptainer.sh>`

    .. group-tab:: Docker
    
        .. literalinclude:: assets/build_lightning_docker.sh
            :caption: :download:`build_lightning_docker.sh <assets/build_lightning_docker.sh>`

Unlike the ``torchaudio`` install, this went fine, and no existing packages changed.
If a package or its dependencies require a different version of PyTorch, you can either change container version based on the `NVIDIA support matrix <https://docs.nvidia.com/deeplearning/frameworks/support-matrix/index.html>`_ to match the required version, or determine if the package's dependencies can be relaxed to match the package version in the container.

Developing python scripts inside a running container on BCP
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

First, download :download:`python_dev.tar.gz <assets/python_dev.tar.gz>` to your current working directory with ``wget``. After downloading, unpack the tarball with ``tar``.

.. code-block:: shell

    # Unpack
    tar -xzf python_dev.tar.gz

    cd python_dev

    ls *

The script ``self_contained.py`` doesn't require any extra python modules other than pytorch, which exists in the container, and can be run directly.
If you take a look at the code, the ``timeit`` functions import functions from the script itself using ``__main__``.
Try running it.

.. code-block:: shell

    python self_contained.py

You'll notice that you can see the file in your other terminal outside the container, while you can run it inside the container.
This means you can be editing scripts in your favorite editor while running them in a container all at the same time.

Developing packages from inside a container
#############################################

If you're developing a whole package that needs to be updated, you either have to rely on relative imports or install the package.
Relative imports often work, but may not depending on the complexity of the package.
In our example python code, there's a ``pt_bench`` python module that we can test with the ``bench.py`` script that imports it.

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
Luckily, python can install packages in a user director that usually defaults to ``$HOME/.local`` using the ``--user`` flag.

.. code-block:: shell

    # Install pt_bench using our constraint file
    pip install -c /root/base_constraints.txt --user python_dev/

    # Try running bench.py again
    python bench.py

You should see that ``pt_bench`` is being loaded from ``$HOME/.local``, which is where the package was installed.
While this works, this location is universally shared by all python packages, which will lead to collisions between containers.
I recommend launching the container with ``--no-home``, which will launch the container with a tmpfs /home.
You'll be able to make changes, like installing a small package, it won't affect the container or bleed into other python environments.

First, lets clean our environment

.. code-block:: shell

    # remove pt_bench
    pip uninstall -y pt_bench

    # Exit the container
    exit

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

Running multi-node containers
--------------------------------

Multi-node NCCL Test
######################

PyTorch containers from NGC ship with `NCCL tests <https://github.com/NVIDIA/nccl-tests>`_, which are useful for diagnosing MPI and bandwidth issues.
These can be run as single-line jobs using ``srun`` to handle the allocation and process spawning.

srun -p gpu -N 2 -n 2 --gpus-per-node 1 --mpi=pmi2 apptainer exec --nv lightning.sif all_reduce_perf_mpi -b 1G -e 4G -f 2 -g 1

# If H100s are available
srun -p gpu --mem=32G -N 2 -n 2 --gpus-per-node h100:3 --mpi=pmi2 apptainer exec --nv lightning.sif all_reduce_perf_mpi -b 1G -e 4G -f 2 -g 3

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
    
Submit the script with ``sbatch``, which will generate a ``.out`` file with a number corresponding to the job with all output text.
You'll see that this runs a training job on 4 GPUs in total, distributed across 2 nodes.
If you increase the resources allocated by the ``SBATCH`` arguements, training will scale as well.

Multi-node Pytorch Lightning
#############################

This is the same task as the `Multi-node PyTorch`_ script, just adapted to PyTorch Lightning.
Download both the training script :download:`ptl_ddp_example.py <assets/ptl_ddp_example.py>` and the ``sbatch`` script below.

.. tabs::

    .. group-tab:: Apptainer

        .. literalinclude:: assets/ptl_ddp_example.sbatch
            :language: shell
            :caption: :download:`ptl_ddp_example.sbatch <assets/ptl_ddp_example.sbatch>`

    .. group-tab:: Enroot

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