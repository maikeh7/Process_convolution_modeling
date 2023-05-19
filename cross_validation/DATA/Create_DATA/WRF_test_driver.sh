#!/bin/bash

for i in {1..10}
  do
  echo "running file $i"
  nohup Rscript LME_normalized_wrftest.R $i > nohup-$i.out &
  echo "done file $i"
  done
