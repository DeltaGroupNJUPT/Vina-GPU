


# Vina-GPU
A heterogeneous OpenCL implementation of AutoDock Vina

An CUDA version of Vina-GPU is avaliable at [here](https://github.com/Glinttsd/Vina-GPU-CUDA)
## Compiling and Running
**Note**: at least one GPU card is required and make sure the version of GPU driver is up to date
### Windows
#### Run on the executable file
1. For the first time to use Vina-GPU, please run `Vina-GPU-K.exe` with command `./Vina-GPU-K.exe --config=./input_file_example/2bm2_config.txt`
You are supposed to have the docking results `2bm2_out.pdbqt` of our example complex and a `Kernel2_Opt.bin` file
2. Once you have the `Kernel2_Opt.bin` file, you can run `Vina-GPU.exe` without compiling the kernel files (thus to save more runtime)
>When you run `Vina-GPU.exe`, please make sure `Kernel2_Opt.bin` file are in the same directory

For the usage and limitaiton of Vina-GPU, please check [Usage](#Usage) and [Limitation](#Limitation).
A graphic user interface (GUI) is also provided for Windows users, please check [GUI](#GUI)
#### Build from source file
>Visual Studio 2019 is recommended for build Vina-GPU from source
1. install [boost library](https://www.boost.org/) (current version is 1.77.0)
2. install [CUDA Toolkit](https://developer.nvidia.com/zh-cn/cuda-toolkit) (current version is v11.5) if you are using NVIDIA GPU cards

    Note: the OpenCL library can be found in CUDA installation path for NVIDIA or in the driver installation path for AMD

3. add `./lib` `./OpenCL/inc` `$(YOUR_BOOST_LIBRARY_PATH)/boost` `$(YOUR_CUDA_TOOLKIT_LIBRARY_PATH)/CUDA/v11.5/include` in the include directories
4. add `$(YOUR_BOOST_LIBRARY_PATH)/stage/lib` `$(YOUR_CUDA_TOOLKIT_PATH)/CUDA/lib/Win32`in the addtional library 
5. add `OpenCL.lib` in the additional dependencies 
6. add `--config=./input_file_example/2bm2_config.txt` in the command arguments
7.  add `WIN32` in the preprocessor definitions if necessary
8. if you want to compile the binary kernel file on the fly, add `BUILD_KERNEL_FROM_SOURCE` in the preprocessor definitions
9. build & run

Note: ensure the line ending are CLRF
### Linux
**Note**: At least 8M stack size is needed. To change the stack size, use `ulimit -s 8192`.
1. install [boost library](https://www.boost.org/) (current version is 1.77.0)
2. install [CUDA Toolkit](https://developer.nvidia.com/zh-cn/cuda-toolkit) (current version is 11.5) if you are using NVIDIA GPU cards

	note: OpenCL library can be usually in `/usr/local/cuda` (for NVIDIA GPU cards)
3. change the `BOOST_LIB_PATH` and `OPENCL_LIB_PATH` accordingly in `Makefile`
4. set GPU platform `GPU_PLATFORM` and OpenCL version `OPENCL_VERSION`in `Makefile`. some options are given below:

	**Note**: -DOPENCL_3_0 is highly recommended in Linux. To check the OpenCL version on a given platform, use `clinfo`.
	
|Macros|Options|Descriptions|
|--|--|--|	
|GPU_PLATFORM|-DNVIDIA_PLATFORM / -DAMD_PLATFORM|NVIDIA / AMD GPU platform
|  OPENCL_VERSION | -DOPENCL_3_0 / -DOPENCL_1_2|OpenCL version 3.0 / 1.2

6. type `make clean` and `make source` to build Vina-GPU that compile the kernel files on the fly (this would take some time at the first use)
7. after a successful compiling, `Vina-GPU` can be seen in the directory 
8. type `./Vina-GPU --config ./input_file_example/2bm2_config.txt` to run Vina-GPU
9. once you successfully run Vina-GPU, its runtime can be further reduced by typing `make clean` and `make` to build it without compiling kernel files (but make sure the `Kernel2_Opt.bin` file is **unchanged**)

10. other compile options: 
  
|Options| Description|
|--|--|
| -g | debug|
|-DDISPLAY_ADDITION_INFO|print addition information
### macOS
>Note: The running Vina-GPU on macOS is not recommended and has not fully tested yet
1. install [boost library](https://www.boost.org/) (current version is 1.77.0)
>modify `Makefile` as follows: 
2. annotate `OPENCL_LIB_PATH`, `OPENCL_INC_PATH` and `-L$(OPENCL_LIB_PATH)/lib64` 
3. add `-framework OpenCL` in `LIB3`
4. type `make` and run

## Usage
|Arguments| Description|Default value
|--|--|--|
|--config | the config file (in .txt format) that contains all the following arguments for the convenience of use| no default
| --receptor | the recrptor file (in .pdbqt format)| no default
|--ligand| the ligand file (in .pdbqt fotmat)| no default
|--thread| the scale of parallelism (docking lanes)|1000
|--search_depth| the number of searching iterations in each docking lane| heuristically determined
|--center_x/y/z|the center of searching box in the receptor|no default
|--size_x/y/z|the volume of the searching box|no default 

## Limitation
|Arguments| Description|Limitation
|--|--|--|
|--thread| the scale of parallelism (docking lanes)| preferably less than 10000
|--size_x/y/z|the volume of the searching box |less than 30/30/30

## Graphic User Interface (GUI)
A graphic user interface (GUI) is provided for users on **Windows** OS
1. first make sure that  `Vina-GPU.exe` can run on a terminal
2. put the `Vina-GPU.exe` and `Kernel2_Opt.bin` files in `./Vina-GPU/GUI/exec` and overwrite the original files
3. run the `Vina-GPU-GUI.exe`file within  `./Vina-GPU/GUI` to start up the Vina-GPU GUI
4. select the input and out output files
5. set the box center, the box size, thread and search_depth
6. click the `start` button to run Vina-GPU
## Citation
* Tang, S.; Chen, R.; Lin, M.; Lin, Q.; Zhu, Y.; Ding, J.; Hu, H.; Ling, M.; Wu, J. Accelerating AutoDock Vina with GPUs. Molecules 2022, 27, 3041. https://doi.org/10.3390/molecules27093041

* O. Trott, A. J. Olson, AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization and multithreading, Journal of Computational Chemistry 31 (2010) 455-461.
