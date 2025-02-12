#include "GA.h"
using namespace std;
using namespace fplll;
template<class ZT, class FT>
GA<ZT, FT>::GA(const char *input_filename,int flags_bkz, int flags_gso, int prec, FloatType float_type) {
    popObj = make_unique<repres<ZT, FT>>(input_filename,flags_bkz, flags_gso, prec, float_type);
}

template<class ZT, class FT>
GA<ZT, FT>::GA(const char* output_filename, int n, int p, int q, int d,int flags_bkz,int flags_gso,int prec,FloatType float_type) {
    popObj = make_unique<repres<ZT, FT>>(output_filename,n,p,q,d, flags_bkz, flags_gso, prec, float_type);
}

template class GA<Z_NR<mpz_t>,FP_NR<mpfr_t>>;
