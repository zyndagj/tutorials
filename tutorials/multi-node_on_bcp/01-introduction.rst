Multi-Node on BCP
=====================

This tutorial introduces ways to run multi-node applications on `NVIDIA Base Command Platform <https://docs.nvidia.com/base-command-platform/user-guide/latest/index.html>`_ (BCP).

Objectives
-----------------

Run multi-node applications using the following:

* Introduction to distributed environments
* MPI
* PyTorch DDP
* PyTorch Lightning
* PyTorch + Ray Train

Requirements
------------

* `DGX Cloud <https://www.nvidia.com/en-us/data-center/dgx-cloud/>`_ instance running BCP
* Pre-authenticated `NGC CLI <https://docs.nvidia.com/base-command-platform/user-guide/latest/index.html#introduction-to-the-ngc-cli>`_

Introduction to distributed environments
------------

Multi-node, or distributed, computing is a model of computation that takes the tasks from an algorithm that can be run independently and executes them across multiple computers.
In the deep learning world, the simplest version of this is the data-distributed parallel model, where a dataset is split across GPUs allocated to the job and a single model is trained in collaboration by all the GPUs.

While DGXs and DGX Cloud were designed for multi-gpu computing, a single process can't span across computers, so code still needs to be written to coordinate the processing on each system.

To learn some basics about working in distributed environments, lets start a 2-node job and explore

.. literalinclude:: assets/interactive_2node.sh
    :caption: :download:`interactive_2node.sh <assets/interactive_2node.sh>`

In this job, we'll be using ``bcprun`` to launch our distributed processes.

.. literalinclude:: assets/bcprun_help.txt
    :caption: bcprun help text

Once your job is running, connect to the jupyter portal or ``ngc batch exec`` into the job and work through the following:

.. literalinclude:: assets/explore_2node.sh
    :caption: :download:`explore_2node.sh <assets/explore_2node.sh>`

Traditionally, multi-node applications utilized Message Passing Interface (MPI) implementations like `MVAPICH <https://mvapich.cse.ohio-state.edu/>`_ for collective operations like broadcast, reduce, and allreduce. Many high performance computing applications still use it, so this tutorial will start with how to run applications with the MPI launcher, ``mpirun``.

MPI
---

MPI requires environment variables to be populated on each node. Schedulers like SLURM automatically populate these variables so MPI knows how many processes and threads to spawn, and how to communicate with other processes. On DGX Cloud with BCP, those variables get set by specifying the ``--array-type "MPI"`` argument when spawning a job.

Multi-node or distributed computing will enable huge speedups by allowing GPUs across multiple computers, or nodes, to cooperate when working when training deep learning models by allowing your program to 
While DGXs and DGX Cloud was designed for distributed computing, code still needs to be written to communicate

.. code-block:: shell

    mpirun --allow-run-as-root -x IBV_DRIVERS=/usr/lib/libibverbs/libmlx5 \
        -np ${NGC_ARRAY_SIZE} -npernode 1 \
        bash -c "all_reduce_perf_mpi -b 64M -e 4G -f 2 -c 0 -n 100 -g ${NGC_GPUS_PER_NODE}"

This can also be run with ``bcprun`` as follows:

.. code-block:: shell

    NGC_ARRAY_TYPE=MPIJob bcprun -no_redirect \
        --launcher 'mpirun --allow-run-as-root' \
        -c "all_reduce_perf_mpi -b 64M -e 4G -f 2 -c 0 -n 100 -g ${NGC_GPUS_PER_NODE}"

You'll notice that both the ``mpirun`` and ``bcprun`` commands are using two environment variables to help make these scripts generally applicable for jobs of various sizes.

.. code-block:: shell

    NGC_ARRAY_SIZE    - Number of nodes allocated to job
    NGC_GPUS_PER_NODE - Number of GPUs allocated per node

This was a 2 node job with 8 gpus per node, and we can double-check these values with

.. code-block:: shell

    $ env | grep -E "(NODE|SIZE)="
    NGC_ARRAY_SIZE=2
    NGC_GPUS_PER_NODE=8

Pytorch DDP
------------

PyTorch Lightning
-----------------


Pytorch + Ray Train
-------------------

`Ray <https://docs.ray.io/en/latest/index.html>`_ is an open source framework to build and scale ML and Python applications.
At it's core, it contains collective operations that can be run on Ray "clusters", or collections of worker processes.

To use Ray across multiple nodes, a Ray cluster needs to first be started across the nodes.
I recommend using the following helper script to start these processes in a BCP environment

.. literalinclude:: assets/start_ray_cluster.sh
    :caption: :download:`start_ray_cluster.sh <assets/start_ray_cluster.sh>`

and the ``bcprun`` command to run the script on each node:

.. code-block:: shell

    bcprun -no_redirect -c 'bash start_ray_cluster.sh'

This first starts the main process on the head node and then worker processes across all other nodes after sleeping for 10 seconds.
Once the cluster is started, you're able to submit jobs for execution on the cluster simply by utilizing the Ray library.

To illustrate this, the `Train a PyTorch Model on Fashion MNIST <https://docs.ray.io/en/latest/train/examples/pytorch/torch_fashion_mnist_example.html>`_ example was modified to accept an argument for the number of workers

.. literalinclude:: assets/ray+pt.py
    :caption: :download:`ray+pt.py <assets/ray+pt.py>`

The script can then be run across 2 nodes (16 GPUs) with the following command:

.. code-block:: shell

    python ray+pt.py --num-workers=16

This can also scale with your job size by using the ``NGC_*`` environment variables to calculate the number of GPUs in your job. 

.. code-block:: shell

    python ray+pt.py --num-workers=$(( $NGC_ARRAY_SIZE * $NGC_GPUS_PER_NODE ))


Once you're done, you can stop the ray cluster with

.. code-block:: shell

    ray stop -f

Next Steps
----------