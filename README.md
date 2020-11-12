State-Observation
===========

[![License](https://img.shields.io/badge/License-BSD%202--Clause-green.svg)](https://opensource.org/licenses/BSD-2-Clause)
[![CI](https://github.com/jrl-umi3218/state-observation/workflows/CI%20of%20state-observation/badge.svg?branch=master)](https://github.com/jrl-umi3218/state-observation/actions?query=workflow%3A%22CI+of+state-observation%22)
[![Documentation](https://img.shields.io/badge/website-online-brightgreen?logo=read-the-docs&style=flat)](https://jrl-umi3218.github.io/state-observation/)


This software provides tools for running state observation for humanoid robots.


Installing
------

## Ubuntu LTS (16.04, 18.04, 20.04)

```bash
# Make sure you have required tools
sudo apt install apt-transport-https lsb-release
# Add our key
sudo apt-key adv --keyserver 'hkp://keyserver.ubuntu.com:80' --recv-key 892EA6EE273707C6495A6FB6220D644C64666806
# Add our repository (stable versions)
sudo sh -c 'echo "deb https://dl.bintray.com/gergondet/multi-contact-release $(lsb_release -sc) main" | sudo tee /etc/apt/sources.list.d/multi-contact.list'
# Use this to setup the HEAD version
# sudo sh -c 'echo "deb https://dl.bintray.com/gergondet/multi-contact-release $(lsb_release -sc) main" | sudo tee /etc/apt/sources.list.d/multi-contact.list'
# Update packages list
sudo apt update
# Install state-observation packages
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
