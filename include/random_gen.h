#pragma once
#include <random>
// setup random generator
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> uniform(0.0, 1.0);


