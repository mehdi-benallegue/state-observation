#!/bin/bash

export INSTALL_PREFIX=/c/devel/install
export CLONE_ONLY=false
export BUILD_TYPE=RelWithDebInfo
readonly this_dir=`cd $(dirname $0); pwd`
readonly SOURCE_DIR=`cd $this_dir/../..; pwd`
cd $this_dir
export BUILD_CORE=`nproc`
export CMAKE_BUILD_PARALLEL_LEVEL=${BUILD_CORE}


build_system_dependency()
{
  REPO=$1
  REF=$2
  SRC=$3
  if [ ! -f $SRC/build/install_manifest.txt ]
  then
    if [ ! -d $SRC ]
    then
      git clone --recursive https://github.com/$REPO $SRC
      cd $SRC
      git checkout $REF
    fi
    if $CLONE_ONLY
    then
      return
    fi
    mkdir -p $SRC/build
    cd $SRC/build
    cmake ../ -DCMAKE_BUILD_TYPE=$BUILD_TYPE         \
              -DCMAKE_INSTALL_PREFIX=$INSTALL_PREFIX \
              -DBUILD_TESTING:BOOL=OFF               \
              ${@:3}
    cmake --build . --config ${BUILD_TYPE}
    if [ $? -ne 0 ]
    then
      echo "Build failed for $1"
      exit 1
    fi
    ${SUDO_CMD} cmake --build . --target install --config ${BUILD_TYPE}
    if [ $? -ne 0 ]
    then
      echo "Installation failed for $1"
      exit 1
    fi
  fi
}


build_system_dependency eigenteam/eigen-git-mirror 3.3.7 "$SOURCE_DIR/eigen"
cd $this_dir/..
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$INSTALL_PREFIX ..
cmake --build . --config ${BUILD_TYPE}
cmake --install . --config ${BUILD_TYPE}
ctest -V ${BUILD_TYPE}
