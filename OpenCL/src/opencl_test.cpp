#include "commonMacros.h"
#include "wrapcl.h"

bool opencl_test() {
    cl_int err;
    size_t size;
    //初始化
    cl_platform_id* platforms;
    cl_device_id* devices;
    cl_context context;
    cl_command_queue queue;
    SetupPlatform(&platforms);
    SetupDevice(platforms, &devices);
    SetupContext(platforms, devices, &context, 1);
    SetupQueue(&queue, context, devices);

    //unsigned char* program_file;
    //cl_program program_cl;
    //cl_program program_head;
    //cl_program program;
    //size_t program_size;
    ////读取内核源码
    //read_file(&program_file, &program_size, "D:/Glinttsd/VSProjects/OpenCL_test/kernel1.cl");
    //program_cl = clCreateProgramWithSource(context, 1, (const char**)&program_file, &program_size, &err); checkErr(err);
    //read_file(&program_file, &program_size, "D:/Glinttsd/VSProjects/OpenCL_test/defines.h");
    //program_head = clCreateProgramWithSource(context, 1, (const char**)&program_file, &program_size, &err); checkErr(err);
    //SetupBuildProgram(program_cl, program_head, devices);
    //err = clUnloadPlatformCompiler(platforms[1]); checkErr(err);
    ////设置内核函数
    //cl_int num_kernel = 1; // user define
    //cl_kernel kernels[1];
    //char kernel_name[][50] = { "hello" };
    //SetupKernel(kernels, program_cl, num_kernel, kernel_name);
    //size_t global_size = 10;
    //err = clEnqueueNDRangeKernel(queue, kernels[0], 1, NULL, &global_size, NULL, 0, NULL, NULL); checkErr(err);
    printf("YES");
    return true;
}