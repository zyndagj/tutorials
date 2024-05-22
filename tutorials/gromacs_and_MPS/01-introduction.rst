GROMACS and MPS
=====================

This tutorial introduces `MPS <https://docs.nvidia.com/deploy/mps/index.html#topic_2>`_ and the benefits of using it with GROMACS.

Objectives
-----------------

* Introduction to GROMACS benchmark
* NVIDIA MPS and benefits
* FOR and MPS
* XARGS and MPS
* Summary
* Next Steps

Requirements
------------

To run all example scripts, you'll need

* `DGX Cloud <https://www.nvidia.com/en-us/data-center/dgx-cloud/>`_ instance running BCP
* Pre-authenticated `NGC CLI <https://docs.nvidia.com/base-command-platform/user-guide/latest/index.html#introduction-to-the-ngc-cli>`_
* 2 `Volta or newer NVIDIA GPUs <https://docs.nvidia.com/deploy/mps/index.html#topic_2_1_2>`_

Introduction to Gromacs benchmark
----------------------------------

This tutorial will be using the GROMACS container built for DGX on NGC

https://catalog.ngc.nvidia.com/orgs/hpc/containers/gromacs

GROMACS can also be built from source, but the pre-compiled container is convenient for this tutorial since it contains all necessary libraries.

First, lets kick off an interactive 2-GPU job on DGX Cloud with the ``gromacs:2023.2`` container:

.. literalinclude:: assets/interactive_job.sh
    :language: bash
    :caption: :download:`interactive_job.sh <assets/interactive_job.sh>`

You'll notice that this job script installs the following packages at runtime:

* **xterm** - used to resize terminal window
* **curl** - used to download benchmark data
* **wget** - used to download benchmark data
* **zip** - used to unzip benchmark data
* **vim-nox** - used for editing files
* **less** - used for viewing files

The GROMACS container is meant to just contain enough packages to run GROMACS, so many useful interactive programs are excluded.
This means the container is small and great for running batch jobs, but makes it sparse for interactive sessions.
Luckily, the container is built from an Ubuntu base, so we can install these extra packages easily.

After your job is running, connect to it and load the AVX2 GROMACS environment

.. code-block:: shell

    # From your localhost
    ngc batch exec <jobid>
    
    # Fix terminal size from inside job
    resize
    
    # Load GROMACS environment
    . /usr/local/gromacs/avx2_256/bin/GMXRC.bash
    
    # cd to the local DGX filesystem
    cd /raid

The benchMEM benchmark
##########################

This tutorial will be using input files from `bench <https://www.mpinat.mpg.de/grubmueller/bench>`__, a free GROMACS benchmark suite with a variety of sizes and run scripts.
Specifically, we will be working with the benchMEM benchmark, which features a protein in membrane system surrounded by water, comprised of 82,000 atoms, and employs a 2 fs time step.
This particular benchmark is well-suited for interactive exploration due to its relatively fast runtime.

Download and unpack it with the following:

.. literalinclude:: assets/downloading_data.sh
    :language: bash
    :caption: :download:`downloading_data.sh <assets/downloading_data.sh>`

Running the benchMEM benchmark
##############################

Now that the benchMEM benchmark is downloaded, we can run it as follows:

.. literalinclude:: assets/run_1gpu.sh
    :language: bash
    :caption: :download:`run_1gpu.sh <assets/run_1gpu.sh>`

The final performance is printed in nanoseconds per day (ns/day). If you think of the simulation as a movie, ns/day is the length and not the walltime of the job, so a higher ns/day is better.

Since the interactive job we submitted has two GPUs, which you can check with ``nvidia-smi``, we can scale this run across both with:

.. literalinclude:: assets/run_2gpu.sh
    :language: bash
    :caption: :download:`run_2gpu.sh <assets/run_2gpu.sh>`

Because this is a small simulation, you should notice that running benchMEM on two GPUs results in a lower ns/day.

benchMEM and GPU utilization
##############################

While these tasks are running, you can look at the telemetry for this BCP job to get an idea of how much of the GPU is being utilized.
If you're not using BCP, you can also use ``nvidia-smi`` as follows:

.. literalinclude:: assets/monitor_1gpu.sh
    :language: bash
    :caption: :download:`monitor_1gpu.sh <assets/monitor_1gpu.sh>`

This script starts by querying the GPU utilization and memory usage every 5 seconds and write it out to the file ``utilization.csv``.

Optional Exercises
##########################

* What happens to performance if you increase the number of steps?
* What happens to performance if you run on the CPU?
* Monitor GPU utilization when using two GPUs

NVIDIA MPS and benefits
------------------------

`NVIDIA Multi-Process Service <https://docs.nvidia.com/deploy/mps/index.html>`_ (MPS) is a feature that allows multiple CUDA applications to share the same GPU, improving system utilization and reducing the overhead of context switching between applications.
By providing a single, unified address space for all MPS clients, MPS enables efficient sharing of GPU resources, making it ideal for scenarios where multiple applications need to access the GPU simultaneously.

If you remember, the utilization of the benchMEM simulation wasn't very good.
Due to the small size of the benchMEM simulation (82k atoms) and the capability of modern GPU infrastructure, the A100 GPUs used in this tutorial did not have very good utilization.
This means that if you ran a bunch of benchMEM simulations on the GPU, the hardware would mostly be idle.

Usually, utilization is improved in the code at the application level by a developer.
Different (larger) inputs can also improve hardware utilization. 
If you're studying the simulation of small molecules, hardware utilization can be improved by running multiple simulations on the same GPU.
Multiple applications on the same GPU, but MPS improves performance by making context switching more efficient.

To illustrate this with an example, lets run four benchMEM simulations at a time on GPU 0 with and without MPS using the following code:

.. literalinclude:: assets/4sim_1gpu.sh
    :language: bash
    :caption: :download:`4sim_1gpu.sh <assets/4sim_1gpu.sh>`

To make it easier to compare throughput performance from the generated log files, use the :download:`calc_throughput.sh <assets/calc_throughput.sh>` script as follows:

.. code-block:: shell

    # Calculate throughput WITHOUT MPS
    $ bash calc_throughput.sh 1gpu_sim*-4.log

    1gpu_sim1-4.log 1gpu_sim2-4.log 1gpu_sim3-4.log 1gpu_sim4-4.log
    1gpu_sim1-4.log:Performance:    41.402
    1gpu_sim2-4.log:Performance:    39.433
    1gpu_sim3-4.log:Performance:    51.208
    1gpu_sim4-4.log:Performance:    41.574
    Total throughput: 173.617

    # Calculate throughput WITH MPS
    $ bash calc_throughput.sh 1gpu_sim*-4_mps.log

    1gpu_sim1-4_mps.log 1gpu_sim2-4_mps.log 1gpu_sim3-4_mps.log 1gpu_sim4-4_mps.log
    1gpu_sim1-4_mps.log:Performance:        45.262
    1gpu_sim2-4_mps.log:Performance:        44.495
    1gpu_sim3-4_mps.log:Performance:        53.269
    1gpu_sim4-4_mps.log:Performance:        44.811
    Total throughput: 187.837

You'll notice total throughput (ns/day) is higher when using MPS to share the GPU. In the next section, we'll figure out what the maximum throughput can be on the A100.

Optional Exercises
##########################

* How is GPU utilization when running these concurrent simulations?

FOR and MPS
------------

In the previous section, we ran each GROMACS processes individually as a separate line in the bash script.
This can also be done in a ``for`` loop if you're looping over files or a range:

.. literalinclude:: assets/Nsim_1gpu.sh
    :language: bash
    :caption: :download:`Nsim_1gpu.sh <assets/Nsim_1gpu.sh>`

By default, this script will run 8 concurrent GROMACS tasks at the same time, but that can also be controlled at runtime with an integer argument.
For example, you can run 5 tasks at at a time with

.. code-block:: shell

    bash Nsim_1gpu.sh 5

Using the :download:`calc_throughput.sh <assets/calc_throughput.sh>` script to determine the total throughput, take some time to identify the optimal number of processes to run at a time.

Optional Exercises
##########################

* Try visualizing the throughput results (N x throughput) in your favorite plotting program (excel counts)

XARGS and MPS
---------------

In the previous example, we used FOR loops to launch GROMACS processes in the background and waited for them to complete.
This worked great since we were just experimenting with 

The ``xargs`` program can be thought of as a "map" operation, where a list of inputs is given and a function is applied to each one.
This is often used to process file contents line by line, but it can also be used with a FOR loop in a script if the loop prints out the command.

.. literalinclude:: assets/xargs.sh
    :language: bash
    :lines: 1-4
    :caption: :download:`xargs.sh <assets/xargs.sh>`

``xargs`` can also call a function in a subshell

.. literalinclude:: assets/xargs.sh
    :language: bash
    :lines: 7-18
    :caption: :download:`xargs.sh <assets/xargs.sh>`

One of the cooler features, is that ``xargs`` can also run tasks in parallel with the ``-P`` argument

.. literalinclude:: assets/xargs.sh
    :language: bash
    :lines: 21-25
    :caption: :download:`xargs.sh <assets/xargs.sh>`

You should notice that this should run 4x faster than the version that called each function sequentially because it runs 4 tasks at the same time.

XARGS on 1 GPU
##########################

Now that we know how to run tasks in parallel in a function, we can apply that format to our GROMACS benchmark.
In the snippet below, you'll notice that we create variables ``NP`` for the number of tasks to run at a time (used by ``-P`` argument) and ``NT`` the number of tasks to generate.
We then create the ``run_gmx`` function to run the benchmark, which takes one argument, ``TID`` the task ID.
First, since this is a sub shell, we need to reload the GROMACS environment.
Then, we can run the GROMACS benchmark along with some help text to know what's running.

.. literalinclude:: assets/xargs_1gpu.sh
    :language: bash
    :caption: :download:`xargs_1gpu.sh <assets/xargs_1gpu.sh>`

Depending on the tasks you need to process, you may need to add additional arguments to this function to accept parameters or files you want to explore.

Optional Exercises
+++++++++++++++++++

* If you increase the number of tasks, does this solution scale nicely?
* If you have time, try looking at the utilization

XARGS on multiple GPUs
##########################

We can also do some math to calculate the GPU index based on the ``SLOT`` index.

.. literalinclude:: assets/xargs_Ngpu.sh
    :language: bash
    :caption: :download:`xargs_Ngpu.sh <assets/xargs_Ngpu.sh>`

If we restricted the number of CPUs to 1 core per task, this would be about double the single-GPU performance.
As an exercise, try changing the both scripts to allocate a single GPU per task to see if this is true.

Optional Exercises
##########################

* If you increase the number of tasks, does the ``xargs`` solution scale nicely?
* Try looking at utilization while running xargs to make sure both GPUs are actually being used.
* Try changing both xargs scripts to allocate a single CPU per task to see if 2-GPU throughput is 2x that of 1-GPU.

Summary
----------

After completing this tutorial, you should have learned the benefit of MPS when running multiple applications and how to efficiently process many tasks across multiple GPUs.


Next Steps
----------

NVIDIA Developer Blog Posts:

* `Maximizing GROMACS Throughput with Multiple Simulations per GPU Using MPS and MIG <https://developer.nvidia.com/blog/maximizing-gromacs-throughput-with-multiple-simulations-per-gpu-using-mps-and-mig/>`_
* `Massively Improved Multi-node NVIDIA GPU Scalability with GROMACS <https://developer.nvidia.com/blog/massively-improved-multi-node-nvidia-gpu-scalability-with-gromacs/>`_

NVIDIA NGC Container:

* `GROMACS Container <https://catalog.ngc.nvidia.com/orgs/hpc/containers/gromacs>`_

Documentation:

* `NVIDIA MPS <https://docs.nvidia.com/deploy/mps/index.html>`_
* `bench <https://www.mpinat.mpg.de/632182/bench.pdf>`_
* `GROMACS <https://manual.gromacs.org/>`_