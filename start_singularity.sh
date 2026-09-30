#!/bin/bash

singularity shell \
  -B /pnfs:/pnfs,/cvmfs:/cvmfs,/exp/annie/data:/exp/annie/data,/exp/annie/app:/exp/annie/app \
  /cvmfs/singularity.opensciencegrid.org/anniesoft/toolanalysis\:latest/

