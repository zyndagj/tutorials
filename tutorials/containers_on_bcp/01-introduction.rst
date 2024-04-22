Containers on BCP
=====================

This tutorial introduces how to locally build python-based containers and how to run them on `NVIDIA Base Command Platform <https://docs.nvidia.com/base-command-platform/user-guide/latest/index.html>`_ (BCP). This is not meant to teach container building mastery, just important topics and how to run on BCP.

*Last updated 4/21/2024*

Objectives
-----------------

* Building and testing your first GPU container
* Best practices for building python-based containers
* Developing python scripts inside a running container on BCP

Requirements
------------

* `DGX Cloud <https://www.nvidia.com/en-us/data-center/dgx-cloud/>`_ instance running BCP
* Pre-authenticated `NGC CLI <https://docs.nvidia.com/base-command-platform/user-guide/latest/index.html#introduction-to-the-ngc-cli>`_
* `Docker CLI Pre-authenticated with nvcr.io <https://org.ngc.nvidia.com/setup/api-key>`_

Building and testing your first GPU container
---------------------------------------------

In this section, we'll be building the `nbody sample benchmark <https://github.com/NVIDIA/cuda-samples/tree/master/Samples/5_Domain_Specific/nbody>`_ from https://github.com/NVIDIA/cuda-samples.
The nbody benchmark demonstrates efficient all-pairs simulation of a gravitational n-body simulation in CUDA and provides a GFLOP/s metric at the end.
While this GFLOP/s metric is not meant for performance comparisons, this sample code supports multiple GPUs and is relatively easy to build.

One of the most popular ways to build containers is with `Docker <https://www.docker.com/>`_, using a `Dockerfile <https://docs.docker.com/reference/dockerfile/>`_.
The Dockerfile is a text file that contains all the commands a user could call to to assemble an image, and feels like making a brand new linux environment ready for development.

A Dockerfile to build the nbody sample
############################

.. literalinclude:: assets/Dockerfile.nbody
    :caption: :download:`Dockerfile.nbody <assets/Dockerfile.nbody>`

Find your ORG and TEAM on BCP
#############################

When building a docker container, it's best practice to `tag it at the same time <https://docs.docker.com/reference/cli/docker/image/build/#tag>`_ with the ``-t`` argument so it's named in your image list. We'll be pushing to our private repositories on ``nvcr.io``, which have the following URL format:

.. code-block:: shell

    nvcr.io/<org>/<team>/<container name>:<tag>

To find your specific ``<org>`` and ``<team>``, print out your current NGC config using

.. code-block:: shell

    ngc config current

    +-------------+----------------------+---------------+
    | key         | value                | source        |
    +-------------+----------------------+---------------+
    | apikey      | ******************** | user settings |
    |             | ******************** |               |
    |             | ******************** |               |
    |             | ******************** |               |
    |             | ZWRk                 |               |
    | format_type | ascii                | user settings |
    | org         | SA-NVEX              | user settings |
    |             | (r2kuatviomfd)       |               |
    | team        | internal-sandbox     | user settings |
    | ace         | sa-nvex-scus-ace     | user settings |
    +-------------+----------------------+---------------+

Lets set some environment variables to help remember your ORG and TEAM values (and make this tutorial generally applicable).

.. code-block:: shell

    # This should be the hash value
    export ORG=r2kuatviomfd
    # Include a slash so we can handle cases without teams
    export TEAM=/internal-sandbox

Make sure to use your own values for these variables, and if you don't have a TEAM, just set it to be empty.

.. code-block:: shell

    export TEAM=

Building your first GPU container
#################################

.. code-block:: shell

    docker build -f Dockerfile.nbody -t nvcr.io/${ORG}${TEAM}/${USER}_nbody:simple .

Lets look at the size of the finished image.

.. code-block:: shell

    docker images

This is a relatively large image that would take a while to push from a local system.
Lets instead figure out how to make our final image more space efficient so it takes up less room and requires less time to push.

Making your container more space efficient
#################################

We can make this much smaller using the following techniques:

#. Use a `multi-staged build <https://docs.docker.com/build/building/multi-stage/>`_
#. Only install runtime libraries in the final container

   #. Using the base container instead of devel
   #. Not installing ``*-devel`` packages.

#. Copy the finished binary instead of the full source repo

.. literalinclude:: assets/Dockerfile.nbody-efficient
    :caption: :download:`Dockerfile.nbody-efficient <assets/Dockerfile.nbody-efficient>`

Download that Dockerfile and build it on your system.

.. code-block:: shell

    docker build -f Dockerfile.nbody-efficient -t nvcr.io/${ORG}${TEAM}/${USER}_nbody:efficient .

Once again, lets look at the final size of that image.

.. code-block:: shell

    docker images

    REPOSITORY                                            TAG        IMAGE ID       CREATED          SIZE
    nvcr.io/r2kuatviomfd/internal-sandbox/gzynda_nbody    simple     3765675fad29   17 seconds ago   8.26GB
    nvcr.io/r2kuatviomfd/internal-sandbox/gzynda_nbody    efficient  d8c9e8efef45   8 minutes ago    443MB


You'll notice it's much smaller: 443MB vs 8.26GB!

Lets push that container to NGC so we can test it.

.. code-block:: shell

    docker push nvcr.io/${ORG}${TEAM}/${USER}_nbody:efficient

Running the nbody sample benchmark
#################################

Now that it's pushed, lets run a sample job.
To give you an idea for how it should be run, you test it out locally with ``docker run``.
You'll notice that you'll get help text, but the actual benchmark won't run unless you have a GPU.

.. code-block:: shell

    docker run --rm -it nvcr.io/${ORG}${TEAM}/${USER}_nbody:efficient bash -l

    root@f5e06e291d60:/# nbody -h
    Run "nbody -benchmark [-numbodies=<numBodies>]" to measure performanc
    e.
            -fullscreen       (run n-body simulation in fullscreen mode)
            -fp64             (use double precision floating point values
    for simulation)
            -hostmem          (stores simulation data in host memory)
            -benchmark        (run benchmark to measure performance) 
            -numbodies=<N>    (number of bodies (>= 1) to run in simulati
    on) 
            -device=<d>       (where d=0,1,2.... for the CUDA device to u
    se)
            -numdevices=<i>   (where i=(number of CUDA devices > 0) to us
    e for simulation)
            -compare          (compares simulation results running once o
    n the default GPU and once on the CPU)
            -cpu              (run n-body simulation on the CPU)
            -tipsy=<file.bin> (load a tipsy model file for simulation)

    NOTE: The CUDA Samples are not meant for performance measurements. Re
    sults may vary when GPU Boost is enabled.

    Error: only 0 Devices available, 1 requested.  Exiting.

.. literalinclude:: assets/nbody_2M_1gpu.sh
    :caption: :download:`Dockerfile.nbody_2M_1gpu.sh <assets/nbody_2M_1gpu.sh>`

When your job is done, you should see output similar to the following:

.. code-block:: shell

    > Windowed mode
    > Simulation data stored in video memory
    > Single precision floating point simulation
    > 1 Devices used for simulation
    GPU Device 0: "Ampere" with compute capability 8.0

    > Compute 8.0 CUDA device: [NVIDIA A100-SXM4-80GB]                    
    Warning: "number of bodies" specified 2000000 is not a multiple of 256.
    Rounding up to the nearest multiple: 2000128.
    2000128 bodies, total time for 10 iterations: 57412.602 ms
    = 696.800 billion interactions per second
    = 13936.007 single-precision GFLOP/s at 20 flops per interaction

Optional Exercises
##########################

* Looking at the help text, try using a different number of GPUs
* Try increasing the number of bodies in the simulation
* Try using double precision

Best practices for building python-based containers
---------------------------------------------------

Python packages can specify dependencies, and sometimes those dependencies can be strictly written where existing packages get changed.
Containers on NVIDIA's NGC contain patched versions of PyTorch and matching libraries that shouldn't be altered if you're looking for optimal performance.
This section will focus on how to install python packages in a way that will prevent changes to the pre-installed packages.

To illustrate this, try installing pytorch from the base ``pytorch:24.03-py3`` container.

.. code-block:: shell

    # Launch pytorch container on your host
    docker run --rm -it nvcr.io/nvidia/pytorch:24.03-py3 bash -l

    # Pip install pytorch inside the running container
    pip install torch torchvision torchaudio

You'll notice that installing these packages installs a bunch of CUDA libraries and downgrades torch.
This is not ideal, since the pytorch containers on NGC already ship cuda libraries, so this not only breaks anything compiled against these libraries, but also needlessly increases the size of the container.

Luckily, you can lock the versions by creating a package `constraints file <https://pip.pypa.io/en/stable/user_guide/#constraints-files>`_, which is similar to a requirements file.

.. literalinclude:: assets/Dockerfile.lightning
    :caption: :download:`Dockerfile.lightning <assets/Dockerfile.lightning>`

.. code-block:: shell

    docker build -f Dockerfile.lightning -t nvcr.io/${ORG}${TEAM}/${USER}_lightning:latest .

If you try to do the torchvision/torchaudio install again using this recipe, it will fail, which lets you know an installation is trying to edit the base pytorch packages.
When this happens, you can either change container version based on the `NVIDIA support matrix <https://docs.nvidia.com/deeplearning/frameworks/support-matrix/index.html>`_ to match the required version, or determine if the package's dependencies can be relaxed to match the package version in the container.

Developing python scripts inside a running container on BCP
-----------------------------------------------------------

`Containers <https://www.docker.com/resources/what-container/>`_ are meant to be static, reproducible checkpoints for your code that can always be started in the same way.
This makes them ideal for porting software to different systems, reproducing results, archiving software, and more.
However, since containers *shouldn't* change once their built (because that would break reproducibility), developing software in them is not always intuitive.

If you try to incorporate all your code in the container as it evolves, this can get tedious - especially if you're pushing these containers to DGX Cloud.
Instead, I recommend making a container with most or all of your dependencies, and mounting your code into the container at runtime.

To explore these concepts, lets launch an interactive environment.

.. literalinclude:: assets/interactive.sh
    :caption: :download:`interactive.sh <assets/interactive.sh>`

Once that job is running, open up the jupyter environment and ``wget`` (or upload) the file :download:`python_dev.tar.gz <assets/python_dev.tar.gz>` to ``/mnt/workspace``.

Developing scripts from inside a container
##############################

Normally, whenever files are created or edited **inside a continer** on BCP, such as in ``/usr/local/lib`` where python installs packages, those files are lost between jobs because they didn't exist in the container that was uploaded to NGC.
If you want files to persist, they need to be written to a workspace outside of the container.
In the case of the job we launched, any changes you make to files in `/mnt/workspace` will persist between jobs because the ``${USER}_containers`` workspace was mounted to this path.

If you head over to the workspace and unpack ``python_dev.tar.gz``, we can explore some quirks.

.. code-block:: shell

    cd /mnt/workspace
    tar -xzf python_dev.tar.gz

First, you'll see the ``requirements.txt`` file.
If you have a script that depends on any python modules, you can list them out along with any versions (ideally not strictly) and install them at runtime.

.. code-block:: shell

    pip list | awk '{print$1"=="$2}' | tail +3 > /root/base_constraints.txt
    pip install -c /root/base_constraints.txt -r requirements.txt

After any dependencies are installed, you can develop and run your script as necessary in your job.
After you're done, you'll just need to download the files through jupyter or with a ``ngc workspace download``.

Developing packages from inside a container
##############################

Packages like ``pt_bench`` are normally installed with

.. code-block:: shell

    pip install -c /root/base_constraints.txt .

However, this copies the whole package to ``/usr/local/lib/python3.10/dist-packages/`` and you'll have to keep uninstalling and installing it when changing files.
To make this easier, you can install packages in *editable* mode with the ``-e`` argument.

.. code-block:: shell

    # Uninstall package
    pip uninstall -y py_bench
    # Install in editable mode
    pip install -c /root/base_constraints.txt -e .

Now, edit the ``__init__.py`` and try re-running ``bench.py`` after moving it to a location that can't locally load the ``pt_bench`` module to see if those changes were propogated.

.. code-block:: shell

    # Move bench.py
    cp bench.py /mnt/bench.py
    # Run from the new location
    cd /mnt && python bench.py

Once again, you can download any modified through jupyter or with a ``ngc workspace download``.

Hopefully this helps improve your development workflow!

Next Steps
----------

NVIDIA Containers:

* `NGC Container Catalog <https://catalog.ngc.nvidia.com/containers>`_
* `Containers for DL Frameworks <https://docs.nvidia.com/deeplearning/frameworks/user-guide/index.html>`_
* `HPC with Containers DLI <https://learn.nvidia.com/courses/course-detail?course_id=course-v1:DLI+L-AC-25+V1>`_

Container workshops/tutorials:

* `Containers@TACC <https://containers-at-tacc.readthedocs.io/en/latest/>`_
* `Getting started with Docker <https://docs.docker.com/get-started/>`_