#include "Individual.h"
template<class ZT,class FT>
unique_ptr<ZT[]> Individual<ZT,FT>::matrix_multiply(ZT* vec, unique_ptr<ZT[]>* mat, int dim) {
    unique_ptr<ZT[]> result = make_unique<ZT[]>(dim);
    for (int i = 0; i < dim; ++i) {
        result[i] = 0;
        for (int j = 0; j < dim; ++j) {
            ZT temp; 
            temp.mul(vec[j], mat[j][i]);
            result[i].add(result[i], temp);
        }
    }
    return result;
}

template<class ZT,class FT>
unique_ptr<ZT[]> Individual<ZT,FT>::YtoX(ZT* y,unique_ptr<FT[]>* mu, int dim) {
    unique_ptr<ZT[]> x = make_unique<ZT[]>(dim);
    for (int i = dim - 1; i >= 0; i--) {
        FT t=0.0;
        for(int j=i+1;j<dim;j++)
        {
            FT temp;
            temp.set_z(x[j]);
            t+=mu[j][i]*temp;
        }
        t.rnd(t);
        ZT t_z;
        t_z.set_f(t);
        x[i].sub(y[i],t_z);
    }
    return x;
}
template<class ZT,class FT>
Individual<ZT,FT>::Individual(int dim, unique_ptr<FT[]>* mu, FT* alpha, unique_ptr<ZT[]>* B,ZT* length, ZT totLength) {

    this->x = make_unique<ZT[]>(dim);
    this->y = make_unique<ZT[]>(dim);
    this->bitvec = make_unique<bool[]>(totLength.get_si());
    this->totLength=totLength;
    this->dim = dim;
    unique_ptr<ZT[]> alpha_r=make_unique<ZT[]>(dim);
    for(int i=0;i<dim;i++)
    {
        if(i<dim/2) continue;
        FT temp=alpha[i];
        alpha[i].floor(alpha[i]);
        ZT use;
        use.set_f(alpha[i]);
        alpha_r[i].mul_ui(use,2);// alpha_r = 2floor(alpha)
        alpha_r[i].add_ui(alpha_r[i],1);// alpha_r = 2*floor[alpha] + 1 
        y[i].randm(alpha_r[i]);// y = [0...2floor(alpha)] -floor(alpha)  [-alpha....0...alpha]
        y[i].sub(y[i],use);//[0...2floor(alpha)] -floor(alpha)  [-alpha....0...alpha]
        alpha[i]=temp;
    }
    for (int i = dim - 1; i >= 0; i--) {

        FT t=0.0;
        for(int j=i+1;j<dim;j++)
        {
            FT temp;
            temp.set_z(x[j]);
            t+=mu[j][i]*temp;
        }
        t.rnd(t);
        ZT t_z;
        t_z.set_f(t);
        x[i].sub(y[i],t_z);    }
        // cout<<endl;
    unique_ptr<ZT[]> vect = (matrix_multiply(x.get(), B, dim));
    this->norm = get_norm(vect.get(), dim);
    bitvec = (encode(y.get(),length,totLength));
    cout<<"tot length "<<totLength.get_si()<<endl;
    for(int i=0;i<totLength.get_si();i++){
        cout<<this->bitvec[i]<<" ";
    }
    cout<<endl;

}
template<class ZT,class FT>
Individual<ZT,FT>::~Individual() 
{
}
template<class ZT,class FT>
Individual<ZT,FT>::Individual()
{
    // cout<<"Hello\n";
    this->x = nullptr;
    this->y = nullptr;
    this->norm = 0.0;
}
template<class ZT,class FT>
Individual<ZT,FT>::Individual(int dim)
{
    this->x = make_unique<ZT[]>(dim);
    this->y = make_unique<ZT[]>(dim);
    this->norm = 0.0;
}
// template<class ZT,class FT>
// Individual<ZT,FT>::Individual(const Individual &t)
// {
//     x = make_unique<ZT[]>(dim);
//     y = make_unique<ZT[]>(dim);
//     norm = t.norm;
//     for(int i = 0;i<dim;i++)
//     {
//         x[i] = t.x[i];
//         y[i] = t.y[i];

//     }
    
// }
template<class ZT, class FT>
Individual<ZT, FT>::Individual(const Individual &t) {
    // Copy the dimension and norm from the original object
    // cout<<"1\n";
    this->dim = t.dim;
    this->norm = t.norm;
    this->totLength = t.totLength;

    // Deep copy the unique_ptr<ZT[]> x
    if (t.x) {
        // Allocate new memory for x in the copied object
        // cout<<dim<<"\n";
        this->x = std::make_unique<ZT[]>(dim);
        // Copy the elements from the original object's x array
        // std::copy(t.x.get(), t.x.get() + dim, this->x.get());
        // cout<<"X copied\n";
        // cout<<t.x[0]<<"sd\n";
        
        for(int i = 0;i<dim;i++)
        {
            // cout<<t.x[i]<<" ";
            x[i] = t.x[i];
        }
        // cout<<"X copied\n";

        
    } else {
        // If the original x was null, keep the new one null as well
        this->x = nullptr;
    }

    // Deep copy the unique_ptr<ZT[]> y
    if (t.y) {
        // Allocate new memory for y in the copied object
        // cout<<dim<<"\n";
        this->y = std::make_unique<ZT[]>(dim);
        // Copy the elements from the original object's y array
        // std::copy(t.y.get(), t.y.get() + dim, this->y.get());
        for(int i = 0;i<dim;i++)y[i] = t.y[i];

    } else {
        // If the original y was null, keep the new one null as well
        this->y = nullptr;
    }


    // Deep copy the unique_ptr<ZT[]> y
    if (t.bitvec) {
        // Allocate new memory for y in the copied object
        // cout<<dim<<"\n";
        this->bitvec = std::make_unique<bool[]>(totLength.get_si());
        // Copy the elements from the original object's y array
        // std::copy(t.y.get(), t.y.get() + dim, this->y.get());
        for(int i = 0;i<totLength.get_si();i++)bitvec[i] = t.bitvec[i];

    } else {
        // If the original y was null, keep the new one null as well
        this->bitvec = nullptr;
    }

    
    // cout<<"Copy Done\n";
}
template<class ZT, class FT>
Individual<ZT, FT>::Individual(Individual &&t) noexcept : dim(t.dim), norm(t.norm),totLength(t.totLength), x(std::move(t.x)), y(std::move(t.y)), bitvec(std::move(t.bitvec)){
    // cout<<"2\n";
    // Leave t in a valid state
    t.dim = 0;
    t.norm = FT();  // Reset the norm
    t.totLength = ZT();
}
template<class ZT, class FT>
Individual<ZT, FT>& Individual<ZT, FT>::operator=(const Individual &t) {

    if (this != &t) {
        dim = t.dim;
        norm = t.norm;
        totLength = t.totLength;
        if (t.x) {
            x = std::make_unique<ZT[]>(dim);
            for(int i = 0;i<dim;i++)x[i] = t.x[i];
        }
        if (t.y) {
            y = std::make_unique<ZT[]>(dim);
            for(int i = 0;i<dim;i++)y[i] = t.y[i];
        }
        if (t.bitvec) {
            bitvec = std::make_unique<bool[]>(totLength.get_si());
            for(int i = 0;i<totLength.get_si();i++)bitvec[i] = t.bitvec[i];
        }
    }
    // cout<<"done\n";
    return *this;
}
template<class ZT, class FT>
Individual<ZT, FT>& Individual<ZT, FT>::operator=(Individual &&t) noexcept {

    if (this != &t) {
        dim = t.dim;
        norm = t.norm;
        x = std::move(t.x);
        y = std::move(t.y);
        bitvec = std::move(t.bitvec);
        totLength = std::move(t.totLength);

        
        // Leave t in a valid state
        t.dim = 0;
        t.norm = FT();
        t.totLength = ZT();
    }
    return *this;
}

template<class ZT,class FT>
unique_ptr<bool[]>  Individual<ZT,FT>::encode(ZT* y, ZT* length, ZT totalLength) {

    // for(int i=0;i<dim;i++){
    //     cout<<y[i]<<" ";
    // }
    // cout<<endl;
    // cout<<"totlengrh "<<totalLength<<endl;
    unique_ptr<bool[]>  chromosome = make_unique<bool[]>((int)totalLength.get_d());
    ZT tot;
    tot = (long)0;
    for (int i = 0; i < dim; i++) {
        bool sign = 0;
        ZT element = y[i];
        if (y[i] < (long)0) {
            sign = 1;
            element.neg(element);
        }
        else if(y[i] == 0){
            random_device rd;
            mt19937 gen(rd());
            uniform_int_distribution<int> distribution(0, 1);
            sign = distribution(gen);
        }
        ZT curEl = y[i];
        for (int j = tot.get_d() + length[i].get_d() - 1; j > tot.get_d(); j--) {
            ZT curBit, mod_;
            mod_=(long)2;
            curBit.mod(curEl, mod_);
            chromosome[j] = (int)curBit.get_d();
            curEl.div_2si(curEl, 1);//divide by 2^1
        }
        chromosome[(int)tot.get_d()] = sign;
        tot.add(tot, length[i]);
    }
    return chromosome;
}

template class Individual<Z_NR<mpz_t>,FP_NR<mpfr_t>>;

