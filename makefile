all: demo

run: demo
	./demo

demo: demo.cpp fft.cpp fft.hpp
	g++ -std=gnu++11 -o demo fft.cpp demo.cpp

clean: 
	rm demo
