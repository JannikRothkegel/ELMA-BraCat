# ELMA-BraCat
ELMA - Extensions of the LeMonADE library
BraCat - Branched Catalysis related to the [paper](https://pubs.acs.org/doi/abs/10.1021/jacs.9b06785)

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

## Run  

### Create branched structures

````sh
    ./CreatorChainWalkingTertiaryBondWalking -f filename.bfm -n 500 -p 0.5 -t 1.0 -b 256
````

## License

See the LICENSE in the root directory.

## How to cite

The development of the code is mainly funded by academic research grants. 
If you find the source code useful please cite the related paper:

[1] R. Dockhorn, L. Plüschke, M. Geisler, J. Zessin, P. Lindner, R. Mundil, J. Merna, J.-U. Sommer, A. Lederer,
    "Polyolefins Formed by Chain Walking Catalysis—A Matter of Branching Density Only?", [J. Am. Chem. Soc. 2019, 141 (31), 15586-15596](https://pubs.acs.org/doi/10.1021/jacs.9b06785)

