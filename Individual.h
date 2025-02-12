#pragma once
#include <cmath>
#include <fplll.h>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include<random>
#include<memory>
#include "test_utils.h"
using namespace std;
using namespace fplll;
template<class ZT,class FT>
class Individual {
public:
    FT norm;
    // static  ZT* common_memory; // Initialization with nullptr
    // ZT* x;
    int dim;
    unique_ptr<ZT[]>x;
    unique_ptr<ZT[]>y;
    unique_ptr<bool[]>bitvec;
    FT get_norm(ZT* vect, int dim);
    unique_ptr<ZT[]>  matrix_multiply(ZT* vec, unique_ptr<ZT[]>* mat, int dim);
    unique_ptr<ZT[]>  YtoX(ZT* y,unique_ptr<FT[]>* mu,int dim);
    Individual(int dim, unique_ptr<FT[]>* mu, FT* alpha, unique_ptr<ZT[]>* B, ZT* length, ZT totLength);
    unique_ptr<bool[]>  encode(ZT* y, ZT* length, ZT totalLength);
    Individual();
    Individual(Individual &&t) noexcept;
    Individual& operator=(Individual &&t) noexcept;
    Individual(int dim);
    Individual(const Individual &t);
    Individual& operator=(const Individual &t);
    ~Individual();
};