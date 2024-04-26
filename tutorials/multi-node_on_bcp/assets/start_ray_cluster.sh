#!/bin/bash

if [ "$NGC_ARRAY_SIZE" -gt "1" ]; then
  export NGC_MASTER_ADDR=launcher-svc-${NGC_JOB_ID}
fi

export PORT=6379

if [ "$NGC_REPLICA_ID" -eq "0" ]; then
	ray start --head --node-ip-address=${NGC_MASTER_ADDR} --port=${PORT}
else
	sleep 10
	ray start --address=${NGC_MASTER_ADDR} --port=${PORT}
fi
