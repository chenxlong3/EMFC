all:
	g++  -std=c++14 -O3 EMFC.cpp -o frim

frim:
	g++  -std=c++14 -O3 EMFC.cpp -o frim

frim_2:
	g++  -std=c++14 -O3 EMFC.cpp -o frim_2

EMFC:
	g++  -std=c++14 -O3 EMFC.cpp -o EMFC

frim_hyp:
	g++  -std=c++14 -O3 EMFC.cpp -o frim_hyp

frim_xi:
	g++  -std=c++14 -O3 EMFC.cpp -o frim_xi

frim_orkut:
	g++  -std=c++14 -O3 EMFC.cpp -o frim_orkut

frim_eps:
	g++  -std=c++14 -O3 EMFC.cpp -o frim_eps

frim_lj:
	g++  -std=c++14 -O3 EMFC.cpp -o frim_lj

frim_root_stat:
	g++  -std=c++14 -O3 EMFC.cpp -o frim_root_stat

# oursP-style K_lb: B(u) && xi_gate>xi_lo && single-flip q-hit (see FRIM_K_LB_GATE_XILO).
frim_rr_graph_prune_xilo: EMFC.cpp alg.cpp
	g++ -std=c++14 -O3 -DFRIM_K_LB_GATE_XILO=1 EMFC.cpp -o frim_rr_graph_prune_xilo

clean:
	rm -f frim frim_2 EMFC frim_hyp frim_xi frim_orkut frim_eps frim_lj frim_root_stat frim_rr_graph_prune_xilo
