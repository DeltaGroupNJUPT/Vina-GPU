
# Vina-GPU
A heterogeneous OpenCL implementation of AutoDock Vina

## Compiling and Running
Notice that at least one GPU card is required and make sure the version of GPU driver is up to date
### Windows
#### Run from executable file
Direcetly run Vina-GPU from executable file `Vina-GPU.exe`.
Take an example on PDBid:2bm2 (all example files in the input_file_example directory) :
1. Create the `2bm2_config.txt` file which includes the receptor and ligand files, center and volume of the search box, the number of threads and search depth
2. make sure that `Vina-GPU.exe`, `2bm2_config.txt` and `Kernel2_Opt.bin` are in the same directory
3. Type `Vina-GPU.exe --config 2bm2_config.txt`
4. Wait untill the docking process finishes and the output file will be `2bm2_out.pdbqt` 
#### Build from source file
Visual Studio 2019 is recommended for build Vina-GPU from source
1. install [boost library](https://www.boost.org/) (current version is 1.77.0)
2. install [CUDA Toolkit](https://developer.nvidia.com/zh-cn/cuda-toolkit) (current version is v11.5) if you are using NVIDIA GPU cards

    note: OpenCL library can be found in CUDA installation path for NVIDIA or in the driver installation path for AMD

3. add `./lib` `./OpenCL/inc` `$(YOUR_BOOST_LIBRARY_PATH)/boost` `$(YOUR_CUDA_TOOLKIT_LIBRARY_PATH)/CUDA/v11.5/include` in the include directories
4. add `$(YOUR_BOOST_LIBRARY_PATH)/stage/lib` `$(YOUR_CUDA_TOOLKIT_PATH)/CUDA/lib/Win32`in the addtional library 
5. add `OpenCL.lib` in the additional dependencies 
6. add `--config=./input_file_example/2bm2_config.txt` in the command arguments
7. if you want to compile the binary kernel file on the fly, add `BUILD_KERNEL_FROM_SOURCE` in the preprocessor definitions
8. build & run

note: ensure the line ending are CLRF
### Linux
1. install [boost library](https://www.boost.org/) (current version is 1.77.0)
2. install [CUDA Toolkit](https://developer.nvidia.com/zh-cn/cuda-toolkit) (current version is 11.5) if you are using NVIDIA GPU cards

	note: OpenCL library can be usually in `/usr/local/cuda` (for NVIDIA GPU cards)
3. change the `BOOST_LIB_PATH` and `OPENCL_LIB_PATH` accordingly in `Makefile`
4. set GPU platform `GPU_PLATFORM` and OpenCL version `OPENCL_VERSION`in `Makefile`. some options are given below:

	**note: -DOPENCL_3_0 is highly recommanded in Linux**
	
|GPU_PLATFORM| OPENCL_VERSION|
|--|--|
| -DNVIDIA_PLATFORM/-DAMD_PLATFORM | -DOPENCL_3_0/-DOPENCL_1_2|

6. type `make clean` and `make source` to build Vina-GPU that compile the kernel files on the fly (this could take some time the first time you build Vina-GPU)
7. after a successful compiling, `Vina-GPU` can be seen in the directory 
8. type `./Vina-GPU --config ./input_file_example/2bm2_config.txt` to run Vina-GPU
9. once you successfully ran Vina-GPU, its runtime can be further reduced by typing `make clean` and `make` to build it without compiling kernel files (but make sure `Kernel_Opt.bin` file is **unchanged**)

10. other compile options: 
  
|Options| Description|
|--|--|
| -g | debug|
|-DDISPLAY_ADDITION_INFO|print addition information

## Usage
|Arguments| Description|Default value
|--|--|--|
| --receptor | the recrptor file (in .pdbqt format)| no default
|--ligand| the ligand file (in .pdbqt fotmat)| no default
|--thread| the scale of parrallelism (docking lanes)|1000
|--search_depth| the number of searching iterations in each docking lane| heuristically determined (please check our paper)
|--center_x/y/z|the center of the searching box in the receptor|no default
|--size_x/y/z|the volumn of the searching box (cannot be too large)|no default 
## Citation
* Shidi, Tang, Chen Ruiqi, Lin Mengru, Lin Qingde, Zhu Yanxiang, Wu Jiansheng, Hu Haifeng, and Ling Ming. "Accelerating AutoDock VINA with GPUs." ChemRxiv (2021). Print.  

* O. Trott, A. J. Olson, AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization and multithreading, Journal of Computational Chemistry 31 (2010) 455-461.
