#!/bin/bash

module purge

GPU_DEVICES_COUNT=`lspci | grep -i nvidia | wc -l`
if [ $GPU_DEVICES_COUNT -gt 0 ]
then
  module load cuda/10.0 cudnn pgi trilinos/12.12.1
else
  module load intel/18.0.1.163 impi/2018.1.163 trilinos/12.12.1 phdf5/1.10.2
fi
