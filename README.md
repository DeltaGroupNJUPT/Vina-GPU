# VINA-GPU
A heterogeneous OpenCL implementation of AutoDock Vina

Notice that at least one GPU card is required and make sure the version of GPU driver is up to date
## Usage
### from executable file
Direcetly run VINA-GPU from executable file `VINA-GPU.exe`.
Take an example on PDBid:2bm2 (all example files in the input_file_example directory) :
1. Create the `2bm2_config.txt` file which includes the receptor and ligand files, center and volume of the search box, the number of threads and search steps 
2. make sure that `VINA-GPU.exe`, `2bm2_config.txt` and `Kernel2_Opt.bin` are in the same directory
3. Type `VINA-GPU.exe --config 2bm2_config.txt`
4. Wait untill the docking process finishes and the output file will be `2bm2_out.pdbqt` 
