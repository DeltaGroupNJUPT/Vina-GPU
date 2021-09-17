#include <wrapcl.h>
//#define DISPLAY_SUCCESS
//#define DISPLAY_ADDITION_INFO

void checkErr(cl_int err) {
    if (CL_SUCCESS != err) {
        //printf("OpenCL error(%d)", err);
        printf("Err(%d)", err);
    }
    else {
        //printf("Success!");
        //printf(" - ");
    }
    fflush(stdout);
}


void read_file(char** program_file, size_t* program_size, std::string file_path) {
    const char* file_path_trans = file_path.data();
    //FILE* program_handle = fopen(file_path_trans, "rb");
    //fseek(program_handle, 0, SEEK_END);
    //if (program_handle == NULL) {
    //    printf("Couldn't find the program file");
    //}
    ////确定文件大小
    //*program_size = (size_t)ftell(program_handle);
    //rewind(program_handle);
    //*program_file = (char*)malloc(*program_size + 1);
    ////读取文件内容，获得源码
    //fread((*program_file), sizeof(char), *program_size, program_handle);
    //(*program_file)[*program_size] = (const char)'\0';

    //fclose(program_handle);

    //std::fstream kernelFile(file_path_trans);
    //std::string content(
    //    (std::istreambuf_iterator<char>(kernelFile)),
    //    std::istreambuf_iterator<char>()
    //);

    //*program_file = new char[content.size()];
    //*program_file = content.c_str();
    //*program_size = content.size();
    size_t size;
    char* str;
    std::fstream f(file_path, (std::fstream::in | std::fstream::binary));

    if (f.is_open())
    {
        size_t fileSize;
        f.seekg(0, std::fstream::end);
        size = fileSize = (size_t)f.tellg();
        f.seekg(0, std::fstream::beg);
        str = new char[size + 1];
        if (!str)
        {
            f.close();
            return;
        }

        f.read(str, fileSize);
        f.close();
        str[size] = '\0';
        *program_file = str;
        //delete[] str;
        *program_size = size+1;
        return;
    }
   printf("Error: failed to open file\n: %s", file_path);
}


void read_n_file(char** program_file, size_t* program_size, std::string file_paths[], size_t num_file) {
    for (int i = 0; i < num_file; i++) {
        read_file(&(program_file[i]), &(program_size[i]), file_paths[i]);
    }
}


void SetupPlatform(cl_platform_id** platforms, cl_int* gpu_platform) {
    cl_int err;
    cl_uint num_platform;
    size_t size;
    std::string nvidia = "NVIDIA";
    std::string amd = "AMD";
    err = clGetPlatformIDs(0, NULL, &num_platform); checkErr(err);
    *platforms = (cl_platform_id*)malloc(sizeof(cl_platform_id) * (num_platform));
    err = clGetPlatformIDs(num_platform, *platforms, NULL); checkErr(err);
    //获取平台信息

    for (int i = 0; i < num_platform; i++) {
        err = clGetPlatformInfo((*platforms)[i], CL_PLATFORM_NAME, 0, NULL, &size); checkErr(err);
        char* platform_name = (char*)malloc(sizeof(char) * size);
        err = clGetPlatformInfo((*platforms)[i], CL_PLATFORM_NAME, size, platform_name, NULL); checkErr(err);
        std::string tmp = platform_name;
        if (tmp.find(nvidia) != std::string::npos || tmp.find(amd) != std::string::npos) {
            *gpu_platform = i;
        }
#ifdef DISPLAY_ADDITION_INFO
        printf("\nPlatform %d : %s\n", i, platform_name);
        free(platform_name);
        err = clGetPlatformInfo((*platforms)[i], CL_PLATFORM_VERSION, 0, NULL, &size); checkErr(err);
        char* platform_version = (char*)malloc(sizeof(char) * size);
        err = clGetPlatformInfo((*platforms)[i], CL_PLATFORM_VERSION, size, platform_version, NULL); checkErr(err);
        printf("\nPlatform %d version: %s\n", i, platform_version);
        free(platform_version);
#endif
    }

}


void SetupDevice(cl_platform_id* platforms, cl_device_id** devices, cl_int gpu_platform) {
    cl_uint num_device;
    cl_int err;
    size_t device_name_size;
    cl_ulong mem_size;
    cl_int N = gpu_platform;
    //初始化第N个平台的设备信息
    err = clGetDeviceIDs(platforms[N], CL_DEVICE_TYPE_GPU, 0, NULL, &num_device); checkErr(err);
    *devices = (cl_device_id*)malloc(sizeof(cl_device_id) * num_device);
    err = clGetDeviceIDs(platforms[N], CL_DEVICE_TYPE_GPU, num_device, *devices, NULL); checkErr(err);
#ifdef DISPLAY_ADDITION_INFO
    //Display devices info
    //for (int i = 0; i < num_device; i++) {
    //    err = clGetDeviceInfo((*devices)[i], CL_DEVICE_NAME, 0, NULL, &device_name_size); checkErr(err);
    //    char* device_name = (char*)malloc(sizeof(char) * device_name_size);
    //    err = clGetDeviceInfo((*devices)[i], CL_DEVICE_NAME, device_name_size, device_name, NULL);
    //    printf("\nPlatform %d device name:%s\n", N, device_name);

    //    err = clGetDeviceInfo((*devices)[i], CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), &mem_size, NULL);
    //    printf("\nPlatform %d global memory size:%f GB\n", N, (double)mem_size/1000000000);

    //    err = clGetDeviceInfo((*devices)[i], CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), &mem_size, NULL);
    //    printf("\nPlatform %d local memory size:%f KB\n", N, (double)mem_size / 1000);
    //}
    //初始化第二个平台的设备信息
    //err = clGetDeviceIDs(platforms[1], CL_DEVICE_TYPE_CPU, 0, NULL, &num_device);checkErr(err);
    //*devices = (cl_device_id*)malloc(sizeof(cl_device_id)*num_device);
    //err = clGetDeviceIDs(platforms[1], CL_DEVICE_TYPE_CPU, num_device, *devices, NULL);checkErr(err);
#endif
}


void SetupContext(cl_platform_id* platforms, cl_device_id* devices, cl_context* context, cl_uint num_device, cl_int gpu_platform_id) {
    cl_int err;
    //选择第一个平台
    cl_context_properties properties[] = { CL_CONTEXT_PLATFORM, (cl_context_properties)(platforms[gpu_platform_id]), 0 };
    cl_int num_device_in_context;
    *context = clCreateContext(properties, num_device, devices, NULL, NULL, &err);
    if (err == CL_SUCCESS) {
#ifdef DISPLAY_SUCCESS
        printf("Create context success!\n");
#endif
    }
    else {
        printf("Fail to create context! Err(%d)\n", err);
    }
    err = clGetContextInfo(*context, CL_CONTEXT_NUM_DEVICES, sizeof(cl_int), &num_device_in_context, NULL); checkErr(err);
    //printf("\nDevice number in Context: %d\n", num_device_in_context);
}


void SetupQueue(cl_command_queue* queue, cl_context context, cl_device_id* devices) {
    cl_int err = 0;
    cl_command_queue_properties props[] = { CL_QUEUE_PROFILING_ENABLE };
    //选择第0个设备
    *queue = clCreateCommandQueue(context, devices[0], *props, &err); checkErr(err);
     
    if (err == CL_SUCCESS) {
#ifdef DISPLAY_SUCCESS
        printf("Create queue success!\n");
#endif
    }
    else {
        printf("Fail to create queue! Err(%d)\n", err);
    }
    
}


void SetupBuildProgramWithSource(cl_program program_cl, cl_program program_head, cl_device_id* devices, std::string include_path, std::string addtion) {
    cl_int err;
    //const char* input_head_names[1] = { "defines.h" };
    //cl_program input_head[1] = { program_head };
    //  -cl-opt-disable -cl-single-precision-constant -cl-unsafe-math-optimizations -cl-finite-math-only -cl-no-signed-zeros -cl-strict-aliasing
    std::string option = " -Werror   -cl-std=CL3.0 -cl-mad-enable -cl-single-precision-constant";
    //std::string option = " -Werror -cl-std=CL3.0 ";
    std::string head_inc_path = "-I ";
    std::string full_path = head_inc_path + include_path + option + addtion;
    const char* options = full_path.data();

    
    //创建程序
    err = clBuildProgram(program_cl, 1, devices, options, NULL, NULL); checkErr(err);
    if (CL_SUCCESS != err) {
        printf("\nError: Failed to build program executable!");
        char* buffer;
        size_t logsize;
        //查看构建日志
        err = clGetProgramBuildInfo(program_cl, *devices, CL_PROGRAM_BUILD_LOG, 0, NULL, &logsize); checkErr(err);
        buffer = (char*)malloc(logsize * sizeof(char));
        err = clGetProgramBuildInfo(program_cl, *devices, CL_PROGRAM_BUILD_LOG, logsize, buffer, NULL); checkErr(err);
        printf("\nlog:%s", buffer);
        free(buffer);
    }
    else {
#ifdef  DISPLAY_SUCCESS
        printf("\nBuild program success!");
        size_t num_kernel;
        err = clGetProgramInfo(program_cl, CL_PROGRAM_NUM_KERNELS, sizeof(size_t), &num_kernel, NULL);
        printf("  Program kernel number: %d\n", (int)num_kernel);
#endif
    }
    /*
        //把头文件程序对象与内核程序对象一起编译
        err = clCompileProgram(program_cl, 0, NULL, 0, 1, input_head, input_head_names, NULL, NULL);
        //链接为程序
        program = clLinkProgram(context, 1, device, NULL, 1, &program_cl, NULL, NULL, &err);
    */
}

void SaveProgramToBinary(cl_program program_cl,const char* file_name) {
    cl_uint err;
    cl_uint numDevices = 0;
    err = clGetProgramInfo(program_cl, CL_PROGRAM_NUM_DEVICES, sizeof(cl_uint), &numDevices, NULL); checkErr(err);
    cl_device_id* devices = (cl_device_id*)malloc(sizeof(cl_device_id) * numDevices);
    err = clGetProgramInfo(program_cl, CL_PROGRAM_DEVICES, sizeof(cl_device_id) * numDevices, devices, NULL); checkErr(err);
    size_t* programBinarySize = (size_t*)malloc(sizeof(size_t) * numDevices);
    err = clGetProgramInfo(program_cl, CL_PROGRAM_BINARY_SIZES, sizeof(size_t) * numDevices, programBinarySize, NULL); checkErr(err);
    unsigned char** programBinaries = (unsigned char**)malloc(sizeof(unsigned char*) * numDevices);
    for (int i = 0; i < numDevices; i++) {
        programBinaries[i] = (unsigned char*)malloc(sizeof(unsigned char) * programBinarySize[i]);
    }
    err = clGetProgramInfo(program_cl, CL_PROGRAM_BINARIES, sizeof(unsigned char*) * numDevices, programBinaries, NULL); checkErr(err);
    FILE* fp = fopen(file_name, "w");
    for (int i = 0; i < numDevices; i++) {
        fwrite(programBinaries[i], 1, programBinarySize[i], fp);
    }
    fclose(fp);
}


cl_program SetupBuildProgramWithBinary(cl_context context, cl_device_id* devices, const char* binary_file_name) {
    cl_int err;
    cl_int binary_status;
    FILE* program_handle = fopen(binary_file_name, "r");
    fseek(program_handle, 0, SEEK_END);
    size_t program_size = ftell(program_handle);
    rewind(program_handle);
    char* binary_buffer = (char*)malloc(program_size + 1);
    fread(binary_buffer, sizeof(char), program_size, program_handle);
    binary_buffer[program_size] = '\0';
    fclose(program_handle);
    cl_program program_cl = clCreateProgramWithBinary(context, 1, devices, &program_size, (const unsigned char**)&binary_buffer, &binary_status, &err); checkErr(err);
    err = clBuildProgram(program_cl, 1, devices, NULL, NULL, NULL); checkErr(err);
    return program_cl;
}



void SetupKernel(cl_kernel* kernels, cl_program program, size_t num_kernel, const char kernel_name[][50]) {
    cl_int err;
    for (int i = 0; i < num_kernel; i++) {
        kernels[i] = clCreateKernel(program, kernel_name[i], &err); checkErr(err);
    }
}

void CreateDeviceBuffer(cl_mem* mem, cl_mem_flags flag, size_t size, cl_context context) {
    cl_int err;
    *mem = clCreateBuffer(context, flag, size, NULL, &err); checkErr(err);
}

void SetKernelArg(cl_kernel kernel, cl_uint num, size_t size, const void *ptr) {
    cl_int err;
    err = clSetKernelArg(kernel, num, size, ptr); checkErr(err);
}