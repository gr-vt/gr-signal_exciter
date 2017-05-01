# gr-signal_exciter

Website: https://github.com/gr-vt/gr-signal_exciter

Last GNU Radio update: 5/1/17

## Installation

Installation has been tested on Ubuntu 14.04 and 16.04 successfully.

### Through PyBOMBS

To install using pybombs you must of the gr-etcetra recipe directory available. As of writing (5/1/17) the updates to gr-etcetra recipes have not been merged into master, but are located within this repo in the directory 'recipes.'

With the recipes installed run (swap 'prefix' for your pybombs prefix)
  ```
  $ pybombs -p prefix install gr-signal_exciter
  ```
In the case that there are cloning issues, and the above fails
  ```
  $ cd path/to/prefix/src/
  $ git clone https://github.com/gr-vt/gr-signal_exciter.git
  $ pybombs -p prefix install gr-signal_exciter
  ```
### Manual Install

Dependencies are:
  ```
  1. cmake
  2. gnuradio
  3. boost
  4. libboost-random
  5. libjsoncpp
  ```
  ```
  $ git clone https://github.com/gr-vt/gr-signal_exciter.git
  $ cd gr-signal_exciter
  $ mkdir build && cd build
  $ cmake .. -DCMAKE_BUILD_TYPE=$cmakebuildtype -DCMAKE_INSTALL_PREFIX=$prefix -Wno-dev
  $ make
  $ make install
  ```

## Updates

If the previous version that made use of libconfig is already present, please do a clean install.
