CC := g++
CPPFLAGS = -O0 -g --std=c++17 -I ../src
# CPPFLAGS = $(CPPFLAGS) -D_GLIBCXX_DEBUG
MATLAB := matlab -nodesktop -nosplash -logfile Y.log -r

.PHONY: clean clear test run all doxygen
.SUFFIXES: .o .cpp .ex .hpp

# 正则表达式表示目录下所有.c文件，相当于：SRCS = main.c a.c b.c
SRCS = TestOrientedJordanCurve.cpp ../src/YinSet/OrientedJordanCurve.cpp ../src/Core/Curve.cpp ../src/Core/Config.cpp ../src/YinSet/SegmentedRealizableSpadjor.cpp ../src/YinSet/YinSet.cpp ../src/YinSet/SegmentsIntersector.cpp ../src/YinSet/PointsLocater.cpp ../src/YinSet/PastingMap.cpp

TARGET = TestOrientedJordanCurve.exe

all: $(TARGET)

$(TARGET): $(SRCS)
	$(CC) -o $@ $^ $(CPPFLAGS)

run: $(TARGET)
	./$(TARGET)

clean:
	-rm $(TARGET) *.o core* *.log *.d
