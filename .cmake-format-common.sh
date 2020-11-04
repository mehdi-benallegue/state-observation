# This script is meant to be sourced from other scripts

# Check for cmake-format
if [[ ! -x "$(command -v cmake-format)" ]]; then
  echo "cmake-format must be installed (pip install cmakelang)"
  exit 1
fi

files=`find src include misc-tools unit-testings utils ./CMakeLists.txt -name '*.cmake' -o -name 'CMakeLists.txt'`
