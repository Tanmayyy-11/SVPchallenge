
#include "ntru.h"

using namespace std;

int mpow(int base, int exp) {

  int result = 1;
  while (exp > 0) {
    if (exp & 1) result = (result * base) ;
    base = (base * base) ;
    exp >>= 1;
  }
  return result;
}

void printvec(vector<int> &v){
    for(int i=0;i<v.size();i++){
        cout<<v[i]<<" ";
    }
    cout<<endl;
}

vector<int> vector_mul(vector<int> &a, vector<int> &b, int mod){
    int n = a.size();
    vector<int> res(n,0);
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            res[(i+j)%n] = (a[i]*b[j] + res[(i+j)%n] + mod*mod) % mod;
        }
    }
    return res;
}

bool is_degree_one(vector<int>& v){
    int n = v.size();
    for(int i=1;i<n;i++){
        if(v[i] != 0) return false;
    }
    return v[0] != 0;
}

vector<int> multiplyByConstant(vector<int> &b, int constant, int mod){
    int n = b.size();
    vector<int> res;
    for(int i=0;i<n;i++){
        res.push_back( (b[i]*constant + constant*mod)%mod);
    }
    return res;
}

int mod_inverse(int a, int m) {
    int m0 = m;
    int t, q;
    int x0 = 0, x1 = 1;
    if (m == 1) {
        return 0;
    }
    while (a > 1) {
        q = a / m;
        t = m;
        m = a % m;
        a = t;
        t = x0;
        x0 = x1 - q * x0;
        x1 = t;
    }
    if (x1 < 0) {
        x1 += m0;
    }
    return x1;
}

void add_vector(vector<int> &a, vector<int>&b, int mod){
    if(a.size()!= b.size()){
        cout<<"error in size"<<endl;
        return;
    }
    for(int i=0;i<a.size();i++){
        a[i] = (a[i]+b[i] + mod)%mod;
    }
}

void sub_vector(vector<int> &a, vector<int>&b, int mod){
    if(a.size()!= b.size()){
        cout<<"error in size"<<endl;
        return;
    }
    for(int i=0;i<a.size();i++){
        a[i] = (a[i]%mod-b[i]%mod + mod)%mod;
    }
}

void multiplyByX(vector<int>& poly) {
    poly.insert(poly.begin(), 0); // Shift all coefficients by 1 to multiply by X...lose the coeff for x^n if exists
    poly.pop_back();
}

void divideByX(vector<int>& poly) {//delete first element that was 0...maintain length by appending a 0
    if (!poly.empty()) {
        poly.erase(poly.begin());
        poly.push_back(0);
    }
}

int deg(vector<int> &v){
    for(int i= v.size()-1;i>-1;i--){
        if(v[i]!=0){
            return i;
        }
    }
}

vector<int> roundMultiplication(vector<int> &a, int power, int mod){
    
    int n = a.size() -1 ;
    vector<int> result (n+1,0);
    for(int i=0;i<n+1;i++){
        result[(i+power + n*n)%n] = (a[i] + result[(i+power + n*n)%n])% mod;
    }
    result[0] = (result[0]+result[n])%mod;
    result[n] = 0;
    return result;

}

vector<int> minus_roundMultiplication(vector<int> &a, int power, int mod){
    
    vector<int> result = roundMultiplication(a,power,mod);
    for(int i=0;i<result.size();i++){
        result[i] = (-result[i]%mod + mod)%mod;
    }
    return result;

}

bool is_one(vector<int> &f){//check f(x) = 1 
    int n = f.size();
    for(int i=1;i<n;i++){
        if(f[i]!=0)return false;
    }
    return (f[0] == 1);
}

bool is_minusone(vector<int> &f, int mod){//check f(x) = 1 
    int n = f.size();
    for(int i=1;i<n;i++){
        if(f[i]!=0)return false;
    }
    return (f[0]%mod == (-1 + mod)%mod);
}

Poly::Poly(int order) : n(order), coeff(order, 0) {}

void Poly::sample(vector<bool> &b, int d1, int d2){// d1 number of 1 and d2 number of -1
    for(int i=0;i<n;i++){
        for(int j =0 ;j<30;j++){
            coeff[i] += (1<<(2+j)) * b[30*i + j];
        }
        if(i<d1)coeff[i]++;
        else if(i<d1+d2)coeff[i]+=2;
    }

    sort(coeff.begin(),coeff.end());

    for(int i=0;i<n;i++){
        coeff[i]&=3; //mod4
        if(coeff[i]==2) coeff[i] = -1;
    }
}

void Poly::multiply(Poly* a, Poly*b, int q){
    for(int i=0;i<n;i++) coeff[i] =0;
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            coeff[(i+j)%n] = ((coeff[(i+j)%n] + a->coeff[i]*b->coeff[j])%q + q)%q;
        }
    }
}

vector<int> Poly::inverse2(){
    vector<int> a (n+1, 0); //initialize a vector of size n+1 is taking quotient x^n - 1 and copy the contents from coeff of 1 to x^n-1
    for(int i=0;i<n;i++){
        a[i] = coeff[i];
    }
    int k=0;
    vector<int> b (n+1, 0);//intialize b as 1... all other coeff as 0
    b[0] = 1;
    vector<int> c (n+1, 0);//intialize c as 0
    vector<int>& f = a ;//reference the same vector a
    vector<int> g (n+1, 0);// represent x*N -1 
    g[0] = 1;
    g[n] = 1;

    bool t = true;

    while(t){

        while(f[0] == 0){   //step3
            divideByX(f);   //step4
            multiplyByX(c);
            k++;
        }

        if(is_one(f)){
            vector<int> ans =  roundMultiplication(b,n-k,2);//step 5
            ans.pop_back();
            return ans;
        }

        if(deg(f) < deg(g)){// step 6
            swap(f,g);  //step 7
            swap(b,c);
        }

        add_vector(f,g,2);
        add_vector(b,c,2);
    }
    
}


vector<int> Poly::inverse3(){
    int mod = 3;
    vector<int> ans;
    vector<int> a (n+1, 0); //initialize a vector of size n+1 is taking quotient x^n - 1 and copy the contents from coeff of 1 to x^n-1
    for(int i=0;i<n;i++){
        a[i] = coeff[i];
    }
    int k=0;
    vector<int> b (n+1, 0);//intialize b as 1... all other coeff as 0
    b[0] = 1;
    vector<int> c (n+1, 0);//intialize c as 0
    vector<int>f (a);//copy
    vector<int> g (n+1, 0);// represent x*N -1 
    g[0] = 2;
    g[n] = 1;

    bool t = true;

    while(t){

        while(f[0] == 0){   //step3
            divideByX(f);   //step4
            multiplyByX(c);
            k++;
        }

        if(is_one(f)){
            vector<int> ans= roundMultiplication(b,n-k,mod);//step 5
            ans.pop_back();
            return ans;
        }
        else if(is_minusone(f,mod)){
            vector<int> ans =  minus_roundMultiplication(b,n-k,mod);
            ans.pop_back();
            return ans;
        }

        if(deg(f) < deg(g)){// step 6
            swap(f,g);  //step 7
            swap(b,c);
        }

        if(f[0] == g[0]){
            sub_vector(f,g,mod);
            sub_vector(b,c,mod);
        }
        else{
            add_vector(f,g,3);
            add_vector(b,c,3);
        }
    }
}

vector<int> Poly::inversep(int mod){
    vector<int> a (n+1, 0); //initialize a vector of size n+1 is taking quotient x^n - 1 and copy the contents from coeff of 1 to x^n-1
    for(int i=0;i<n;i++){
        a[i] = coeff[i];
    }
    int k=0;
    vector<int> b (n+1, 0);//intialize b as 1... all other coeff as 0
    b[0] = 1;
    vector<int> c (n+1, 0);//intialize c as 0
    int upperbound = n + deg(a);
    int counter = 0;

    vector<int>& f =a;
    vector<int> g (n+1, 0);// represent x*N -1 
    g[0] = 2;
    g[n] = 1;
    bool t = true;

    while(t){
        while(f[0] == 0){   //step3
            divideByX(f);   //step4
            counter++;
            if(counter == upperbound){
                cout<<"inverse doesnt exist"<<endl;
                return {};
            }
            multiplyByX(c);
            k++;
        }
        if(is_degree_one(f)){
            int f0_inv = mod_inverse(f[0],mod);
            b = multiplyByConstant(b,f0_inv,mod);
            vector<int> ans = roundMultiplication(b,n-k,mod);//step 5
            ans.pop_back();
            return ans;
        }
        if(deg(f) < deg(g)){// step 6
            swap(f,g);  //step 7
            swap(b,c);
        }
        int g0_inv = mod_inverse(g[0],mod);
        int u = f[0]*g0_inv % mod;
        vector<int> gTemp = multiplyByConstant(g,u,mod);
        sub_vector(f,gTemp,mod);
        vector<int> cTemp = multiplyByConstant(c,u,mod);
        sub_vector(b,cTemp,mod);
    }
}

vector<int> Poly::inverse_p_power_r(int p, int r){
    vector<int> a (n, 0); //initialize a vector of size n+1 is taking quotient x^n - 1 and copy the contents from coeff of 1 to x^n-1
    for(int i=0;i<n;i++){
        a[i] = coeff[i];
    }

    vector<int> b = this->inverse2();
    if(!b.size()) return b;
    int q = p;
    int bound  = mpow(p,r);
    while(q<bound){
        q = q*q;
        vector<int> t = vector_mul(a,b,q);
        t[0] = 2 - t[0]; //check
        for(int i=1;i<t.size();i++) t[i] = -t[i];
        b = vector_mul(b,t,q);
    }
    for(int i=0;i<b.size();i++){
        b[i] = (b[i] + bound)%bound;
    }
    return b;
}

NTRU::NTRU(int order, int modulusP, int modulusQ, int dVal): p(modulusP), q(modulusQ), n(order), d(dVal)
{
    f = new Poly(order);
    g = new Poly(order);
    h = new Poly(order);
    Fp = new Poly(order);
    Fq = new Poly(order);
    basis.resize(2 * n, 2 * n);

    unsigned char entropy_input[48];
    unsigned char personalization_string[48] = {0};
    srand(static_cast<unsigned int>(time(nullptr)));
    for (int i = 0; i < 48; i++) {
        entropy_input[i] = rand() & 0xFF;
    }
    randombytes_init(entropy_input, personalization_string, 256);

    sample_size = ceil((30 * n) / 8);
    seed = new unsigned char[sample_size];
    randombytes(seed,sample_size);
    for (size_t i = 0; i < sample_size; ++i) {
        unsigned char byte = seed[i];
        for (int bit = 7; bit >= 0; --bit) {
            b.push_back((byte >> bit) & 1);
        }
    }

}

void NTRU::createkey(){
    while(true){
        f->sample(b,d+1,d);
        // f->coeff = {0, 0, 1, 1, 0, 1, 1, 0, -1, -1, -1};
        Fp->coeff = f->inverse3();
        if(Fp->coeff.size()){
            Fq->coeff = f->inverse_p_power_r(2,floor(log2(q)));
            if(Fq->coeff.size()){
                break;
            }
        }
    }
    g->sample(b,d,d);
    // g->coeff = {0, 0, 0, -1, 0, -1, 1, 1, -1, 0, 1};
    h->multiply(Fq,g,q);
}

void NTRU::generate_lattice(){
    mpz_t q_mpz;
    mpz_init_set_si(q_mpz, q);

    for(int i=0;i<n;i++){
        basis(i,i) = q_mpz;
        basis(i+n,i+n) = 1;
    }
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            basis(i+n,j) = h->coeff[(i-j + n )%n];
        }
    }

}


void NTRU::writeBasisToFile(const std::string& filename) {
    std::ofstream outfile(filename);
    if (!outfile) {
        std::cerr << "Error: Could not open file " << filename << " for writing.\n";
        return;
    }
    outfile << "[";
    for (long i = 0; i < basis.get_rows(); i++) {
        outfile << "[";
        for (long j = 0; j < basis.get_cols(); j++) {
            outfile << basis(i, j);
            if (j < basis.get_cols() - 1) outfile << " ";
        }
        outfile << "]\n";
    }
    outfile << "]";
    
    outfile.close();
    std::cout << "Matrix written to " << filename << "\n";
}

