# Makefile for the project
SHELL := /bin/bash   # 确保使用 Bash
.PHONY: clean run
.SUFFIXES: .o .cpp .ex .hpp


$(TARGET): $(SRCS)
	$(CC) -o $@ $^ $(CPPFLAGS)

run: 
	. init.sh
	. releaseconfigure.sh gcc
	. compile.sh
	. parallelTest.sh

clean:
	source clear.sh
