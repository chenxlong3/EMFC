all:
	g++  -std=c++14 -O3 subsim.cpp -o frim

frim:
	g++  -std=c++14 -O3 subsim.cpp -o frim

frim_2:
	g++  -std=c++14 -O3 subsim.cpp -o frim_2

EMFC:
	g++  -std=c++14 -O3 subsim.cpp -o EMFC

frim_hyp:
	g++  -std=c++14 -O3 subsim.cpp -o frim_hyp

frim_xi:
	g++  -std=c++14 -O3 subsim.cpp -o frim_xi

frim_orkut:
	g++  -std=c++14 -O3 subsim.cpp -o frim_orkut

frim_eps:
	g++  -std=c++14 -O3 subsim.cpp -o frim_eps

frim_lj:
	g++  -std=c++14 -O3 subsim.cpp -o frim_lj

frim_root_stat:
	g++  -std=c++14 -O3 subsim.cpp -o frim_root_stat

clean:
	rm -f frim frim_2 EMFC frim_hyp frim_xi frim_orkut frim_eps frim_lj frim_root_stat
