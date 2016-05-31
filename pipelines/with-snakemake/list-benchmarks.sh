#!/bin/bash

for benchmark in `ls benchmarks`; do
  echo $benchmark && cat ./benchmarks/$benchmark
done
