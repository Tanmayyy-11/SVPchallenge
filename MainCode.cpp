#include <cstring>
#include <memory>
#include <fplll.h>
#include "test_utils.h"
#include "paper-1ls.h"
#include <iostream>
#include <fplll.h>
#include <cmath>
#include <fplll.h>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include<random>
using namespace std;
using namespace fplll;
using namespace std::chrono;

#include "Individual.h"


#ifndef TESTDATADIR
#define TESTDATADIR ".."
#endif

int main()
{
    //Paper1<Z_NR<mpz_t>,FP_NR<mpfr_t>>* p1 = new Paper1<Z_NR<mpz_t>,FP_NR<mpfr_t>>(n,p,q,d,BKZ_DEFAULT,GSO_DEFAULT,120,FT_DEFAULT);
    // cout<<"Program Started\n";
    // Paper1<Z_NR<mpz_t>,FP_NR<mpfr_t>>* p1 = new Paper1<Z_NR<mpz_t>,FP_NR<mpfr_t>>("ntrulattices/newbasis.txt",BKZ_DEFAULT,GSO_DEFAULT,120,FT_DEFAULT);
    Paper1<Z_NR<mpz_t>,FP_NR<mpfr_t>>* p1 = new Paper1<Z_NR<mpz_t>,FP_NR<mpfr_t>>("lattices/svpchallengedim50seed0.txt",BKZ_DEFAULT,GSO_DEFAULT,120,FT_DEFAULT);
    // Paper1_ls<Z_NR<mpz_t>,FP_NR<mpfr_t>>* p1 = new Paper1_ls<Z_NR<mpz_t>,FP_NR<mpfr_t>>("lattices/svpchallengedim50seed0.txt",BKZ_DEFAULT,GSO_DEFAULT,120,FT_DEFAULT);
    // Paper1<Z_NR<mpz_t>,FP_NR<mpfr_t>>* p1 = new Paper1<Z_NR<mpz_t>,FP_NR<mpfr_t>>("ntrulattices/newlattice_67.txt",67,3,256,22,BKZ_DEFAULT,GSO_DEFAULT,120,FT_DEFAULT);;

    Individual<Z_NR<mpz_t>,FP_NR<mpfr_t>>v0;
    // cout<<"v0 created\n";
    auto start = high_resolution_clock::now();

    unique_ptr<Z_NR<mpz_t>[]> ans=p1->runGA(1703,10,v0);
    cout<<"runGA complete\n";
    for(int i=0;i<60;i++)cout<<ans[i]<<" ";
    cout<<"\n";
    auto end = high_resolution_clock::now();
    // Calculate duration and print in milliseconds
    auto duration = duration_cast<milliseconds>(end - start);
    cout << "Time taken by function: " << duration.count() << " miliseconds" << endl;
    cout<<endl;
}
