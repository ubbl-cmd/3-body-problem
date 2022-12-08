all:
	g++ -Wall src/methods/explicit_rk.hpp src/methods/explicit_rk.cpp src/methods/explicit_adams.cpp src/methods/explicit_adams.hpp src/main.cpp -o main
clear:
	rm -f trajectory barycenter_coorinates barycenter_speed energy potential_energy kinetic_energy main
