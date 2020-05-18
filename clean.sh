#!/bin/bash

rm $1/output/bake/*.*

if [ $# -gt 1 ]; then
  if [ $2 == "all" ]; then
    rm -r $1/output
  fi
fi