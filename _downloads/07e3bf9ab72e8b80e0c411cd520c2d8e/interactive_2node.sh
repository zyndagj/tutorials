# Make a workspace
ngc workspace create --name ${USER}_tutorial

# Launch a 2-node job with the workspace mounted
ngc batch run --name "N2-test" --total-runtime 1h \
        --instance dgxa100.80g.8.norm \
        --commandline "jupyter lab --ip=0.0.0.0 --allow-root \
            --no-browser --NotebookApp.token='' \
            --notebook-dir=/ --NotebookApp.allow_origin='*'" \
        --result /results --array-type "PYTORCH" --replicas "2" \
        --workspace ${USER}_tutorial:/workspace:RW \
        --image "nvidia/pytorch:23.03-py3" --port 8888