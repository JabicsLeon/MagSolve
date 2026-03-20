#include "MagSolve.hpp"

#include <iostream>

int main(int argc, char* argv[]) 
{
	if (magsolve::pipline(argc, argv)) return 1;

	std::cout << "Creat by LeoLib. MSU 2026\n";

	return 0;
}

