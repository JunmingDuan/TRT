LIBDIR:=              # 静态库目录
LIBS :=               # 静态库文件名
INCLUDES:=../include /home/jimmy/eigen        # 头文件目录

CC:=g++
CFLAGS := -g -Wall -O3 -std=c++11 -fopenmp
CPPFLAGS := $(CFLAGS)
CPPFLAGS += $(addprefix -I,$(INCLUDES))

RM-F := rm -f

SRCS := $(wildcard *.cpp)
OBJS := $(patsubst %.cpp,%.o,$(SRCS))

.PHONY : all clean info

all : $(OBJS)
	$(CC) -c -o $(OBJS) $(addprefix -L,$(LIBDIR)) $(addprefix -l,$(LIBS)) $(CPPFLAGS)

clean :
	@$(RM-F) $(OBJS)

info:
	@echo $(SRCS)
	@echo $(OBJS)

