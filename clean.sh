#!/bin/bash

rm gama/output/bake/*.*

if [ $# -gt 0 ]; then
  if [ $1 == "all" ]; then
    rm -r gama/output
  fi
fi