#!/bin/bash

readonly this_dir=`cd $(dirname $0); pwd`
cd $this_dir
source .cmake-format-common.sh

for f in ${files}; do
  cmake-format -i $f
done
