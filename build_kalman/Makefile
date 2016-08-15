# Makefile for kalman

PROG_NAME := kalman
OBJS += ./kalman.o ./kalman_functions.o 

LIBS := -lopencv_core -lopencv_imgproc -lopencv_highgui -lopencv_imgcodecs -lboost_system -lopencv_video


OPENCV_LIB_DIRS := -L/usr/local/lib
OPENCV_INCLUDES := -I/usr/local/include/opencv
OPENCV2_INCLUDES := -I/usr/local/include/opencv2

OPENCV_LIB_DIRS += $(shell pkg-config --libs-only-L opencv 2> /dev/null)
OPENCV_INCLUDES += $(shell pkg-config --cflags opencv 2> /dev/null)
OPENCV2_INCLUDES += $(shell echo $(OPENCV_INCLUDES) | sed -e 's/opencv/opencv2/g')
RM := rm -rf

# All Target
all: $(PROG_NAME)

# Tool invocations
$(PROG_NAME): $(OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++ $(OPENCV_LIB_DIRS) -o $(PROG_NAME) $(OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++14 $(OPENCV_INCLUDES) $(OPENCV2_INCLUDES) -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


# Other Targets
clean:
	-$(RM) *.d *.o $(PROG_NAME)
	-@echo ' '

