#include "IndexCalculus.h"

uint64_t inv(uint64_t a, uint64_t n)
{
	int64_t v = 0;
	int64_t v_ = 1;
	int64_t r = n;
	int64_t r_ = a;

	while (r_ != 0)
	{
		int64_t q = r / r_;

		int64_t temp = v;
		v = v_;
		v_ = temp - q * v_;

		temp = r;
		r = r_;
		r_ = temp - q * r_;
	}

	if (r > 1)
	{
		return 0;
	}

	if (v < 0)
	{
		v += n;
	}

	return (uint64_t)v;
}

uint64_t gcd(uint64_t x, uint64_t p)
{
	while (x != 0)
	{
		p %= x;
		swap(x, p);
	}
	return p;
}

vector<uint64_t> rozvazok(uint64_t a, uint64_t b, uint64_t n) // ax = b mod n
{
	uint64_t d = gcd(a, n);
	vector<uint64_t> roz;

	if (b % d != 0)
	{
		return roz;
	}

	uint64_t a_ = a / d;
	uint64_t b_ = b / d;
	uint64_t n_ = n / d;

	uint64_t in = inv(a_, n_);
	uint64_t x_0 = (in * b_) % n_;

	for (int k = 0; k < d; k++)
	{
		roz.push_back(x_0 + k * n_);
	}

	return roz;
}

uint64_t Gorn(uint64_t a, uint64_t b, uint64_t n)
{
	uint64_t k = floor(log2(b));
	uint64_t y = 1;
	while (b > 0)
	{
		if (b & 1)
		{
			y = (y * a) % n;
		}

		a = (a * a) % n;
		b >>= 1;
	}
	return y;
}

uint64_t bf(uint64_t a, uint64_t b, uint64_t p, vector<uint64_t> x)
{
	int s = x.size();
	for (int i = 0; i < s; i++)
	{
		if (Gorn(a, x[i], p) == b)
		{
			return x[i];
		}
	}
	cerr << "Error, no element in vector is correct \n";
	return 0;
}

vector<int64_t> gladkist(vector<int64_t>& B, int64_t x, int64_t& B_smooth)
{
	vector<int64_t> rez(B.size(), 0);


	for (size_t i = 1; i < B.size(); i++)
	{
		int pow = 0;
		while (x % B[i] == 0)
		{
			pow++;
			x /= B[i];
		}
		rez[i] = pow;
	}

	B_smooth = (x == 1);
	return rez;
}

uint64_t index_calculus(uint64_t a, uint64_t b, uint64_t p)
{
	double L = exp(sqrt(0.5) * sqrt(log(p) * log(log(p))));
	ifstream file("D:/problems/term_6/NTA/lab_1/First10MillionPrimes.txt");

	vector<int64_t> S;

	uint64_t current_num;
	while ((file >> current_num) && current_num < L)
	{
		S.push_back(current_num);
	}
	file.close();

	return 1;
}
