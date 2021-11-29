# ppmapm

Header only implementation of "Progressive Photon Mapping: A Probabilistic Approach"(PPMAPA) in C++.

In this reformulation of (stochastic) progressive photon mapping, not so many code changes are required from original photon mapping. You don't have to prepare hit points and local statistics of photons.

Due to this simplicity, most of the codes come from [the implementation of original photon mapping](https://github.com/yumcyaWiz/photon_mapping)

Algorithms are summarized as follows.

1. Photon tracing
2. Build photon map
3. Ray tracing from the eye and estimate the reflected radiance with photon map
4. Reduce the search radius for density estimation
5. go back to 1

Note that the procedure of 1, 2, 3 is the same as original photon mapping.

## Features

* Depth of field with thin lens camera
* Load obj model

## Requirements

* C++ (20>=)
* CMake (3.20>=)
* OpenMP
* [spdlog](https://github.com/gabime/spdlog)
* [Embree](https://github.com/embree/embree) (>=3)

## Build

|CMake option|Description|
|:--|:--|
|BUILD_TESTS|build tests|

```
git submodule update --init
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

## Externals

* [spdlog](https://github.com/gabime/spdlog)
* [Embree](https://github.com/embree/embree)
* [tinyobjloader](https://github.com/tinyobjloader/tinyobjloader)

## References

* [yumcyaWiz/photon_mapping](https://github.com/yumcyaWiz/photon_mapping)
* [Hachisuka, Toshiya, Shinji Ogaki, and Henrik Wann Jensen. "Progressive photon mapping." ACM SIGGRAPH Asia 2008 papers. 2008. 1-8.](https://doi.org/10.1145/1457515.1409083)
* [Hachisuka, Toshiya, and Henrik Wann Jensen. "Stochastic progressive photon mapping." ACM SIGGRAPH Asia 2009 papers. 2009. 1-8.](https://doi.org/10.1145/1661412.1618487)
* [Knaus, Claude, and Matthias Zwicker. "Progressive photon mapping: A probabilistic approach." ACM Transactions on Graphics (TOG) 30.3 (2011): 1-13.](https://doi.org/10.1145/1966394.1966404)