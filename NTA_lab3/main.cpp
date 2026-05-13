#include "IndexCalculus.h"

int main()
{
	uint64_t a = 70272181;
	uint64_t x = 1673757;
	uint64_t p = 179424671;

	uint64_t b = Gorn(a, x, p);
	uint64_t test = index_calculus(a, b, p);
	cout << test << " here \n";
	cout << Gorn(a, test, p) << " \n" << Gorn(a, x, p);
	return 0;
}
