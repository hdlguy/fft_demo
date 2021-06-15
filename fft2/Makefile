all: demo 

demo: demo.c fft.c fft.h
# $@: target file name, $^: all prerequisites
	gcc -g -o $@ $^ -lm
#	g++ -std=gnu++11 -o demo fft.cpp demo.cpp

clean: 
	rm demo
