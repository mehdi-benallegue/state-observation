State-Observation
===========

[![License](https://img.shields.io/badge/License-BSD%202--Clause-green.svg)](https://opensource.org/licenses/BSD-2-Clause)
[![Hosted By: Cloudsmith](https://img.shields.io/badge/OSS%20hosting%20by-cloudsmith-blue?logo=cloudsmith)](https://cloudsmith.com)
[![CI](https://github.com/jrl-umi3218/state-observation/workflows/CI%20of%20state-observation/badge.svg?branch=master)](https://github.com/jrl-umi3218/state-observation/actions?query=workflow%3A%22CI+of+state-observation%22)
[![Documentation](https://img.shields.io/badge/website-online-brightgreen?logo=read-the-docs&style=flat)](https://jrl-umi3218.github.io/state-observation/)


This software provides tools for running state observation. It has three levels of uses:
* The high level users can exploit the developped estimators.
* The intermediate level users can derive these estimators to fit their system.
* The low level users can develop their novel estimators relying on the available tools.


Installing
------

## Ubuntu LTS (16.04, 18.04, 20.04)

You must first setup our package mirror:

```
curl -1sLf \
  'https://dl.cloudsmith.io/public/mc-rtc/stable/setup.deb.sh' \
  | sudo -E bash
```

You can also choose the head mirror which will have the latest version of this package:

```
curl -1sLf \
  'https://dl.cloudsmith.io/public/mc-rtc/stable/setup.deb.sh' \
  | sudo -E bash
```

You can then install the package:

```bash
sudo apt install libstate-observation-dev
# Install documentation
sudo apt install libstate-observation-doc
```

## Manually build from source

To compile you need the following tools:

 * [Git]()
 * [CMake]() >= 2.8
 * [pkg-config]()
 * [doxygen]()
 * [g++]() >= 4.7 (for C++11 support)
 * [Boost](http://www.boost.org/doc/libs/1_58_0/more/getting_started/unix-variants.html) >= 1.49
 * [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) >= 3.2

### Building

```sh
git clone --recursive https://github.com/jrl-umi3218/state-observation
cd state-observation
mkdir build
cd build
cmake [options] ..
make && make intall
```
