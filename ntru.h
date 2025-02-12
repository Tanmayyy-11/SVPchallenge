#ifndef NTRU_H
#define NTRU_H

#pragma once
#include <cmath>
#include <fplll.h>
#include "test_utils.h"
#include <chrono>
#include <iostream>
#include<utility>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <fstream>
#include <memory>
#include<random>
#include "rng.h"

class Poly;
class NTRU;
int mpow(int base, int exp);
void printvec(vector<int> &v);
vector<int> vector_mul(vector<int> &a, vector<int> &b, int mod);
bool is_degree_one(vector<int>& v);
vector<int> multiplyByConstant(vector<int> &b, int constant, int mod);
int mod_inverse(int a, int m);
void add_vector(vector<int> &a, vector<int>&b, int mod);
void sub_vector(vector<int> &a, vector<int>&b, int mod);
void multiplyByX(vector<int>& poly);
void divideByX(vector<int>& poly);
int deg(vector<int> &v);
vector<int> roundMultiplication(vector<int> &a, int power, int mod);
vector<int> minus_roundMultiplication(vector<int> &a, int power, int mod);
bool is_one(vector<int> &f);
bool is_minusone(vector<int> &f, int mod);

// Poly class declaration
class Poly {
public:
    int n;               // Order
    vector<int> coeff;   // Coefficients

    Poly(int order); 
    void sample(vector<bool> &b, int d1, int d2);
    void multiply(Poly* a, Poly* b, int q);
    vector<int> inverse2();
    vector<int> inverse3();
    vector<int> inversep(int mod);
    vector<int> inverse_p_power_r(int p, int r);
};

// NTRU class declaration
class NTRU {
public:
    Poly *f;             // Polynomial f
    Poly *g;             // Polynomial g
    Poly *h;             // Polynomial h
    Poly *Fp;            // Polynomial Fp
    Poly *Fq;            // Polynomial Fq
    int p, q, n, d, r;   // Moduli, order, and other parameters
    unsigned long long sample_size; 
    ZZ_mat<mpz_t> basis; // Lattice basis matrix
    unsigned char *seed;
    vector<bool> b;

    NTRU(int order, int modulusP, int modulusQ, int dVal);

    void createkey();
    void generate_lattice();
    void writeBasisToFile(const std::string& filename);
};

#endif