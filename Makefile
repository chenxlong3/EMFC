all:
	g++  -std=c++14 -O3 subsim.cpp -o frim

clean:
	rm -f frim
