# ELMA-BraCat
ELMA - Extensions of the LeMonADE library
Repository with updates, analyzers, and projects for sharing BFM stuff related to various topics.
BraCat - Branched Catalysis

## Installation

* Clone and Install `git clone https://github.com/LeMonADE-project/LeMonADE.git`
* Install cmake (minimum version 2.8)
* Just do for standard compilation:
 
````sh
    # generates the projects
    mkdir build
    cd build
    cmake -DLEMONADE_INCLUDE_DIR=/path/to/LeMonADE-library/include/ -DLEMONADE_LIBRARY_DIR=/path/to/LeMonADE-library/lib/ ..
    make
````
or

````sh
    cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=Release -DLEMONADE_INCLUDE_DIR=/path/to/LeMonADE-library/include/ -DLEMONADE_LIBRARY_DIR=/path/to/LeMonADE-library/lib/
    cmake --build build/
````

## Environment
* tested with g++ version GNU 10.2.1
* [LeMonADE library](https://github.com/LeMonADE-project/LeMonADE/commit/0cd5499379a6abeab2f9e1c5f67b10edd288481b) version 2.2.1
* cmake version 3.18.4


## License

See the LICENSE in the root directory.
