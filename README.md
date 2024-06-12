# Tempus

Requires the [Boost](https://www.boost.org/) library to be installed.

Requires the [oneTBB](https://github.com/oneapi-src/oneTBB) library to be installed.

May require Linux to run properly.

This project uses [Dynamic Connectivity](https://github.com/tomtseng/dynamic-connectivity-hdt) as the underlying graph structure.

Might require a change to line 45 in `include/dynamic-connectivity-hdt/CMakeLists.txt` to the g++ version to the one installed on your machine.
