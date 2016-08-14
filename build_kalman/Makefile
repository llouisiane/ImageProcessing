# Makefile for segmentation

PROG_NAME := kalman
OBJS += ./kalman.o ./kalman_functions.o 

LIBS := -lopencv_core -lopencv_imgproc -lopencv_highgui -lopencv_imgcodecs -lboost_system -lopencv_video

LIB_DIR := /usr/local/lib
OPENCV_DIR := /usr/local/include/opencv
OPENCV2_DIR := /usr/local/include/opencv2

RM := rm -rf

# All Target
all: $(PROG_NAME)

# Tool invocations
$(PROG_NAME): $(OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++ -L$(LIB_DIR) -o $(PROG_NAME) $(OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++14 -I$(OPENCV_DIR) -I$(OPENCV2_DIR) -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


# Other Targets
clean:
	-$(RM) *.d *.o $(PROG_NAME)
	-@echo ' '

