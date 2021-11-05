LIB0=-I../boost_1_77_0/ -I../boost_1_77_0/boost -I./lib -I./OpenCL/inc -I/usr/local/cuda/include
LIB1=-lboost_program_options -lboost_system -lboost_filesystem -lOpenCL
LIB2=-lstdc++
LIB3=-lm -lpthread
LIB4=-fpermissive
LIB5=-L../boost_1_77_0/stage/lib -L/usr/local/cuda/lib64
SRC=./lib/*.cpp ./OpenCL/src/wrapcl.cpp ../boost_1_77_0/libs/thread/src/pthread/thread.cpp ../boost_1_77_0/libs/thread/src/pthread/once.cpp #../boost_1_77_0/boost/filesystem/path.hpp
MACRO=-DBUILD_KERNEL_FROM_SOURCE -DDISPLAY_SUCCESS -DDISPLAY_ADDITION_INFO
OPTION=-Wnarrowing -O3 #-g
all:out
out:./main/main.cpp
	gcc -o Vina-GPU $(LIB0) ./main/main.cpp $(SRC) $(LIB1) $(LIB2) $(LIB3) $(LIB4) $(LIB5) $(MACRO) $(OPTION)
clean:
	rm Vina-GPU
