#include "IndexCalculus.h"

uint64_t mult(uint64_t a, uint64_t b, uint64_t p)
{
	uint64_t res = 0;
	a %= p;

	while (b > 0)
	{
		if (b % 2 == 1)
		{
			if (res >= p - a)
			{
				res -= p;
			}
			res = (res + a) % p;
		}

		if (a >= p - a)
		{
			a = a * 2 - p;
		}
		else
		{
			a = a * 2;
		}

		b /= 2;
	}
	return res;
}

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

	v %= n;
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
	uint64_t x_0 = mult(in, b_, n_);

	for (int k = 0; k < d; k++)
	{
		roz.push_back(x_0 + k * n_);
	}

	return roz;
}

uint64_t Gorn(uint64_t a, uint64_t b, uint64_t n)
{
	uint64_t y = 1;
	while (b > 0)
	{
		if (b & 1)
		{
			y = mult(y, a, n);
		}

		a = mult(a, a, n);
		b >>= 1;
	}
	return y;
}

uint64_t bf(uint64_t a, uint64_t b, uint64_t p, vector<uint64_t>& x)
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

vector<int64_t> gladkist(vector<int64_t>& B, int64_t x, bool& B_smooth)
{
	vector<int64_t> rez(B.size(), 0);


	for (size_t i = 0; i < B.size(); i++)
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

vector<int64_t> gauss(vector<vector<int64_t>>& M, uint64_t mod)
{
	int n = M.size();
	int m = M[0].size() - 1;
	int row = 0;

	for (int k = 0; k < m && row < n; k++)
	{
		int pivot = row;
		while (pivot < n && M[pivot][k] == 0) pivot++;

		if (pivot == n) continue;
		swap(M[row], M[pivot]);

		for (int i = row + 1; i < n; i++)
		{
			while (M[i][k] != 0)
			{
				int64_t q = M[row][k] / M[i][k];
				if (q == 0)
				{
					swap(M[row], M[i]);
					continue;
				}
				for (int j = k; j <= m; j++)
				{
					int64_t sub = mult((uint64_t)q, (uint64_t)M[i][j], mod);
					M[row][j] = (M[row][j] - sub + mod) % mod;
				}
				swap(M[row], M[i]);
			}
		}

		uint64_t common = gcd(M[row][k], mod);
		if (common == 1)
		{
			uint64_t inv_val = inv(M[row][k], mod);
			for (int j = k; j <= m; j++)
			{
				M[row][j] = mult(M[row][j], inv_val, mod);
			}

			for (int i = 0; i < row; i++)
			{
				if (M[i][k] != 0)
				{
					uint64_t factor = M[i][k];
					for (int j = k; j <= m; j++)
					{
						int64_t sub = mult(factor, M[row][j], mod);
						M[i][j] = (M[i][j] - sub + mod) % mod;
					}
				}
			}
			row++;
		}
	}

	vector<int64_t> rez(m, 0);
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < m; j++)
		{
			if (M[i][j] == 1)
			{
				rez[j] = M[i][m];
				break;
			}
		}
	}
	return rez;
}

uint64_t get_random_R(uint64_t p)
{
	static std::random_device rd;
	static std::mt19937_64 gen(rd());

	std::uniform_int_distribution<uint64_t> dist(1, p - 2);

	return dist(gen);
}

vector<uint64_t> factor(uint64_t n)
{
	vector<uint64_t> primes;

	for (uint64_t d = 2; d * d <= n; d++)
	{
		if (n % d == 0)
		{
			primes.push_back(d);

			while (n % d == 0)
				n /= d;
		}
	}

	if (n > 1)
		primes.push_back(n);

	return primes;
}

bool is_generator(uint64_t a, uint64_t p)
{
	uint64_t phi = p - 1;

	vector<uint64_t> divisors = factor(phi);

	for (uint64_t q : divisors)
	{
		if (Gorn(a, phi / q, p) == 1)
			return false;
	}

	return true;
}

uint64_t index_calculus(uint64_t a, uint64_t b, uint64_t p)
{
	if (!is_generator(a, p))
	{
		cerr << "a is not generator\n";
		return 0;
	}

	double L = exp(sqrt(0.5) * sqrt(log(p) * log(log(p))));
	ifstream file("D:/problems/term_6/NTA/lab_1/First10MillionPrimes.txt");

	vector<int64_t> S;

	uint64_t current_num;
	while ((file >> current_num) && current_num < L)
	{
		S.push_back(current_num);
	}
	file.close();

	int p_ = S.size()+20;

	uint64_t k_ = 0;
	uint64_t k = 1;

	vector<vector<int64_t>> M;
	while (M.size() < p_)
	{
		uint64_t temp = Gorn(a, k, p);
		bool t = 0;
		vector<int64_t> rozklad = gladkist(S, temp, t);
		if (t)
		{
			cout << k_ << "\n";
			k_++;
			rozklad.push_back(k);
			M.push_back(rozklad);
		}
		k++;
	}

	vector<int64_t> logs = gauss(M, p - 1);

	cout << "-----" << "\n";
	for (int i = 0; i < logs.size(); i++)
	{
		cout << logs[i] << "\n";
	}
	cout << "-----";

	bool all_zero = true;

	for (auto l : logs)
	{
		if (l != 0)
		{
			all_zero = false;
		}
	}

	if (all_zero)
	{
		cerr << "Warning: Gaussian elimination returned empty/zero logs. Need more equations.\n";
	}


	bool is_smooth = false;
	while (true)
	{
		uint64_t R = get_random_R(p);

		uint64_t a_pow_R = Gorn(a, R, p);
		uint64_t b_prime = mult(b, a_pow_R, p);

		vector<int64_t> exponents = gladkist(S, b_prime, is_smooth);

		if (is_smooth)
		{
			int64_t sum_logs = 0;
			for (size_t i = 0; i < S.size(); i++)
			{
				if (exponents[i] > 0)
				{
					int64_t term = mult(exponents[i], logs[i], p - 1);
					sum_logs = (sum_logs + term) % (p - 1);
				}
			}

			int64_t x = (sum_logs - (int64_t)R) % (int64_t)(p - 1);
			if (x < 0)
			{
				x += p - 1;
			}

			if (Gorn(a, x, p) == b)
			{
				return x;
			}
		}
	}

}

uint64_t find_generator(uint64_t p)
{
    for (uint64_t g = 2; g < p; g++)
    {
        if (is_generator(g, p))
            return g;
    }

    return 0;
}

uint64_t index_calculus_patalel(uint64_t a, uint64_t b, uint64_t p)
{
	if (!is_generator(a, p))
	{
		cerr << "a is not generator\n";
		return 0;
	}

	double L = exp(sqrt(0.5) * sqrt(log(p) * log(log(p))));
	ifstream file("D:/problems/term_6/NTA/lab_1/First10MillionPrimes.txt");

	vector<int64_t> S;

	uint64_t current_num;
	while ((file >> current_num) && current_num < L)
	{
		S.push_back(current_num);
	}
	file.close();

	int p_ = S.size() + 20;

	vector<vector<int64_t>> M;

	std::mutex mtx;
	std::atomic<uint64_t> current_k{ 1 };

	auto worker = [&]()
		{
		while (true)
		{
			bool done = false;
			{
				std::lock_guard<std::mutex> lock(mtx);
				if (M.size() >= p_) done = true;
			}
			if (done)
			{
				break;
			}

			uint64_t k = current_k.fetch_add(1);

			uint64_t temp = Gorn(a, k, p);
			bool t = false;
			vector<int64_t> rozklad = gladkist(S, temp, t);

			if (t)
			{
				rozklad.push_back(k);

				std::lock_guard<std::mutex> lock(mtx);
				if (M.size() < p_)
				{
					cout << M.size() << "\n";
					M.push_back(rozklad);
				}
				else
				{
					break;
				}
			}
		}
	};

	unsigned int num_threads = std::thread::hardware_concurrency();
	if (num_threads == 0) num_threads = 4;

	vector<std::thread> threads;
	for (unsigned int i = 0; i < num_threads; ++i)
	{
		threads.emplace_back(worker);
	}

	for (auto& th : threads)
	{
		th.join();
	}

	vector<int64_t> logs = gauss(M, p - 1);

	cout << "-----" << "\n";
	for (int i = 0; i < logs.size(); i++)
	{
		cout << logs[i] << "\n";
	}
	cout << "-----";

	bool all_zero = true;

	for (auto l : logs)
	{
		if (l != 0)
		{
			all_zero = false;
		}
	}

	if (all_zero)
	{
		cerr << "Gauss elimination returned zero logs. Need more equations.\n";
	}

	bool is_smooth = false;
	while (true)
	{
		uint64_t R = get_random_R(p);

		uint64_t a_pow_R = Gorn(a, R, p);
		uint64_t b_prime = mult(b, a_pow_R, p);

		vector<int64_t> exponents = gladkist(S, b_prime, is_smooth);

		if (is_smooth)
		{
			int64_t sum_logs = 0;
			for (size_t i = 0; i < S.size(); i++)
			{
				if (exponents[i] > 0)
				{
					int64_t term = mult(exponents[i], logs[i], p - 1);
					sum_logs = (sum_logs + term) % (p - 1);
				}
			}

			int64_t x = (sum_logs - (int64_t)R) % (int64_t)(p - 1);
			if (x < 0)
			{
				x += p - 1;
			}

			if (Gorn(a, x, p) == b)
			{
				return x;
			}
		}
	}
}
