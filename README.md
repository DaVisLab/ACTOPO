# Accelerated Clustered TOPOlogical (*ACTOPO*) Data Structure

This repository contains the implementation of the ACTOPO data structure integrated into [Topology Toolkit](https://topology-tool-kit.github.io/) framework. 


## Organization

The repo contains the TTK framework with the necessary plugins that we used in the experiments. It has two branches: `sequential-algos`, which allows the topological data structures to run in parallel while the plugins run in sequential, and `parallel-algos`, which has both the data structure and plugins run in parallel. 

The implementation of ACTOPO can be found in `core/base/actopo`. The implementation of each evaluation plugins can also be found in the `core/base` folder with the following folder names.  
- `testTopoRelations`
- `scalarFieldCriticalPoints`
- `discreteGradient`
- `morseSmaleComplex3D` 


## How to Install

To install ACTOPO, you can simply follow the [official guide](https://topology-tool-kit.github.io/installation.html) to install TTK on your system. Here we provide a simplified walkthrough without the ParaView support on the Linux system (i.e., plugins can only run through the command line). 

### 1. Clone the repository
Simply run the following command to clone this repo to your local machine.  
```
git clone https://github.com/DaVisLab/ACTOPO.git
```

### 2. Install the dependencies 
Several dependencies need to be installed in order to compile TTK from source. Please enter the following commands in a termainal to install them:  
```
sudo apt-get install cmake-qt-gui libboost-system-dev libpython3.8-dev libxt-dev libxcursor-dev libopengl-dev
sudo apt-get install qt5-default qttools5-dev libqt5x11extras5-dev libqt5svg5-dev qtxmlpatterns5-dev-tools
```

### 3. Configure, build and install 
#### Configuration
To enter the configuration menu of TTK's build, enter the following commands:  
```
cd ACTOPO/
mkdir build-sequential
cd build-sequential
cmake-gui ../
```

The configuration window opens. Click on the `Configure` button to proceed. Since we do not need to activate ParaView support, please set the `TTK_BUILD_PARAVIEW_PLUGINS` to `OFF`, and also make sure `TTK_BUILD_STANDALONE_APPS` is `ON`. 

#### Build
Now you can start the compilation process by entering the following command, where `N` is the number of available cores on your system (this probably takes a **LONG** time):  
```
make -jN
```

#### Installation
Once the build is finished, enter the following command to install TTK on your system:  
```
sudo make install
```

To check if you successfully installed the TTK, enter the following commands:  
```
ls /usr/local/bin/*Cmd
ls /usr/local/bin/*Gui
```

### Switch branches 
To switch between sequential and parallel algorithms, we suggest that you create a separate build folder, e.g., `build-parallel`, and save the bulid configurations of `parallel-algos` branch. Every time you want to install a specific version on the system, just go to the corrseponding build folder and run `sudo make install` command.  


## How to Use 

ACTOPO is activated in the same way as `TTK CompactTriangulation`, i.e., when the input dataset contains the vertex-based clustering information as a scalar field named `_index` (defined in the file `core/base/common/Datatypes.h`). If you already have such datasets, you can skip the clustering step. 

### Cluster the dataset 
TTK provides a preconditioning plugin (`ttkCompactTriangulationPreconditioning`) that can be used to divide the input mesh based on the point region octree. To use it, enter the following command in the terminal (replace `input.vtu` and `output`): 
```
ttkCompactTriangulationPreconditioningCmd -i <input.vtu> -o <output>
```

To verify if the generated dataset contains the `_index` field, you can run any standalone plugin with `-l` option. For example, the followig command checks all data arrays in `output.vtu`: 
```
ttkTestTopoRelationsCmd -i <output.vtu> -l
``` 

### Run with ACTOPO
After you have the dataset with clustering information, you can run TTK standalone apps with ACTOPO. Some command-line options are dedicated to ACTOPO data structure: 
- `-t <int>` specifies the number of consumer threads (has no effects in `sequential-algo` build as algorithms run sequentially). 
- `-p <int>` specifies the number of producer threads (including one leader producer and (num-1) worker producers).
- `-b <float>` specifies the buffer size as a ratio of total number of blocks in the dataset.

The following command runs the `TestTopoRelations` plugin with 6 producer threads, and the buffer size is set to 20\% of the total number of blocks.
```
ttkTestTopoRelationsCmd -i <input.vtu> -a "field" -t 1 -p 6 -b 0.2
``` 

## License 
ACTOPO follows the BSD license used by TTK, see [LICENSE](./LICENSE) for more details. 

## Acknowledgements
[To be updated]
