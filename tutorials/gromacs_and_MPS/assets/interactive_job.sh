ngc batch run --name "gromacs-2gpu" \
	--total-runtime 7200s --instance dgxa100.80g.2.norm \
	--commandline "apt update && DEBIAN_FRONTEND=noninteractive apt install -y xterm curl wget zip vim-nox less && sleep 2h" \
	--result /results \
	--image "nvcr.io/hpc/gromacs:2023.2"
