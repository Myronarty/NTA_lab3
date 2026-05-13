#pragma once
#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <thread>
#include <mutex>
#include <atomic>

using namespace std;

uint64_t bf(uint64_t a, uint64_t b, uint64_t p, vector<uint64_t>& x);

uint64_t Gorn(uint64_t a, uint64_t b, uint64_t n);

uint64_t index_calculus(uint64_t a, uint64_t b, uint64_t p);

uint64_t index_calculus_patalel(uint64_t a, uint64_t b, uint64_t p);

uint64_t find_generator(uint64_t p);
