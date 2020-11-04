#!/bin/bash

readonly this_dir=`cd $(dirname $0); pwd`
cd $this_dir
source .cmake-format-common.sh

tmpfile=$(mktemp /tmp/cmake-format-check.XXXXXX)
out=0
for f in ${files}; do
  cmake-lint $f > $tmpfile
  if [ ! $? ]; then
    echo "Wrong formatting in $f"
    cat $tmpfile
    out=1
  fi
done
rm -f $tmpfile
if [[ $out -eq 1 ]]; then
  echo "You can run ./.cmake-format-fix.sh to fix the issues locally, then commit/push again"
fi
exit $out
