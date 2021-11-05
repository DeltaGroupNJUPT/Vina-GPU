
# Vina-GPU
A heterogeneous OpenCL implementation of AutoDock Vina

## Compiling and Running
Notice that at least one GPU card is required and make sure the version of GPU driver is up to date
### Windows
#### run from executable file
Direcetly run Vina-GPU from executable file `Vina-GPU.exe`.
Take an example on PDBid:2bm2 (all example files in the input_file_example directory) :
1. Create the `2bm2_config.txt` file which includes the receptor and ligand files, center and volume of the search box, the number of threads and search depth
2. make sure that `Vina-GPU.exe`, `2bm2_config.txt` and `Kernel2_Opt.bin` are in the same directory
3. Type `Vina-GPU.exe --config 2bm2_config.txt`
4. Wait untill the docking process finishes and the output file will be `2bm2_out.pdbqt` 
#### build from source file
Visual Studio 2019 is recommended for build Vina-GPU from source
1. install [boost library](https://www.boost.org/) (current version is 1.77.0)
2. install [CUDA Toolkit](https://developer.nvidia.com/zh-cn/cuda-toolkit) (current version is v11.5) if you are using NVIDIA GPU cards
**note: OpenCL library can be found in CUDA installation path for NVIDIA or in the driver installation path for AMD**
3. add `./lib` `./OpenCL/inc` `$(YOUR_BOOST_LIBRARY_PATH)/boost` `$(YOUR_CUDA_TOOLKIT_LIBRARY_PATH)/CUDA/v11.5/include` in the include directories
4. add `$(YOUR_BOOST_LIBRARY_PATH)/stage/lib` `$(YOUR_CUDA_TOOLKIT_PATH)/CUDA/lib/Win32`in the addtional library 
5. add `OpenCL.lib` in the additional dependencies 
6. add `--config=./input_file_example/2bm2_config.txt` in the command arguments
7. if you want to compile the binary kernel file on the fly, add `BUILD_KERNEL_FROM_SOURCE` in the preprocessor definitions
8. build & run
### Linux
1. install [boost library](https://www.boost.org/) (current version is 1.77.0)
2. install [CUDA Toolkit](https://developer.nvidia.com/zh-cn/cuda-toolkit) (current version is 11.5) if you are using NVIDIA GPU cards
**note: OpenCL library can be usually in `/usr/local/cuda` (for NVIDIA GPU cards)**
3. change the `BOOST_LIB_PATH` and `OPENCL_LIB_PATH` accordingly
4. type `make clean` and `make`
5. other compile options: 
  
|Debug| Build kernels on the fly | Output additional info|
|--|--|--|
| -g | -DBUILD_KERNEL_FROM_SOURCE | -DDISPLAY_ADDITION_INFO|

## Citation
* Shidi, Tang, Chen Ruiqi, Lin Mengru, Lin Qingde, Zhu Yanxiang, Wu Jiansheng, Hu Haifeng, and Ling Ming. "Accelerating AutoDock VINA with GPUs." ChemRxiv (2021). Print.  

* O. Trott, A. J. Olson, AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization and multithreading, Journal of Computational Chemistry 31 (2010) 455-461.