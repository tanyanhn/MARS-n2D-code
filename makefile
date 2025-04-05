# Makefile for the project
SHELL := /bin/bash   # 确保使用 Bash
.PHONY: clean run
.SUFFIXES: .o .cpp .ex .hpp

TARGET = build/test/TestMARSn2D

run: 
	. init.sh && . releaseconfigure.sh gcc && . compile.sh && ./$(TARGET)

clean:
	source clear.sh
