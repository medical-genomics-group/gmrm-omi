#include <iostream>
#include <mpi.h>
#include <limits.h>
#include "utilities.hpp"

double round_dp(const double in) {
    return in;

    printf("in = %20.15f\n", in);
    double out = round(in * 1.0E12) / 1.0E12;
    printf("ou = %20.15f\n", out);
    return out;
}

void check_malloc(const void* ptr, const int linenumber, const char* filename) {
    if (ptr == NULL) {
        fprintf(stderr, "#FATAL#: malloc failed on line %d of %s\n", linenumber, filename);
        fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

void check_mpi(const int error, const int linenumber, const char* filename) {
    if (error != MPI_SUCCESS) {
        fprintf(stderr, "*FATAL*: MPI error %d at line %d of file %s\n", error, linenumber, filename);
        fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

// Check whether a size_t can be casted to int or would overflow
int check_int_overflow(const size_t n, const int linenumber, const char* filename) {

    if (n > INT_MAX) {
        fprintf(stderr, "FATAL  : integer overflow detected on line %d of %s. %lu does not fit in type int.\n", linenumber, filename, n);
        fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    return static_cast<int>(n);
}

double normal_cdf(double value)
{
   return 0.5 * erfc(-value * M_SQRT1_2);
}

double erfcx (double x)
{
    double a, d, e, m, p, q, r, s, t;

    a = fmax (x, 0.0 - x); // NaN preserving absolute value computation

    /* Compute q = (a-4)/(a+4) accurately. [0,INF) -> [-1,1] */
    m = a - 4.0;
    p = a + 4.0;
    r = 1.0 / p;
    q = m * r;
    t = fma (q + 1.0, -4.0, a); 
    e = fma (q, -a, t); 
    q = fma (r, e, q); 

    /* Approximate (1+2*a)*exp(a*a)*erfc(a) as p(q)+1 for q in [-1,1] */
    p =             0x1.edcad78fc8044p-31;  //  8.9820305531190140e-10
    p = fma (p, q,  0x1.b1548f14735d1p-30); //  1.5764464777959401e-09
    p = fma (p, q, -0x1.a1ad2e6c4a7a8p-27); // -1.2155985739342269e-08
    p = fma (p, q, -0x1.1985b48f08574p-26); // -1.6386753783877791e-08
    p = fma (p, q,  0x1.c6a8093ac4f83p-24); //  1.0585794011876720e-07
    p = fma (p, q,  0x1.31c2b2b44b731p-24); //  7.1190423171700940e-08
    p = fma (p, q, -0x1.b87373facb29fp-21); // -8.2040389712752056e-07
    p = fma (p, q,  0x1.3fef1358803b7p-22); //  2.9796165315625938e-07
    p = fma (p, q,  0x1.7eec072bb0be3p-18); //  5.7059822144459833e-06
    p = fma (p, q, -0x1.78a680a741c4ap-17); // -1.1225056665965572e-05
    p = fma (p, q, -0x1.9951f39295cf4p-16); // -2.4397380523258482e-05
    p = fma (p, q,  0x1.3be1255ce180bp-13); //  1.5062307184282616e-04
    p = fma (p, q, -0x1.a1df71176b791p-13); // -1.9925728768782324e-04
    p = fma (p, q, -0x1.8d4aaa0099bc8p-11); // -7.5777369791018515e-04
    p = fma (p, q,  0x1.49c673066c831p-8);  //  5.0319701025945277e-03
    p = fma (p, q, -0x1.0962386ea02b7p-6);  // -1.6197733983519948e-02
    p = fma (p, q,  0x1.3079edf465cc3p-5);  //  3.7167515521269866e-02
    p = fma (p, q, -0x1.0fb06dfedc4ccp-4);  // -6.6330365820039094e-02
    p = fma (p, q,  0x1.7fee004e266dfp-4);  //  9.3732834999538536e-02
    p = fma (p, q, -0x1.9ddb23c3e14d2p-4);  // -1.0103906603588378e-01
    p = fma (p, q,  0x1.16ecefcfa4865p-4);  //  6.8097054254651804e-02
    p = fma (p, q,  0x1.f7f5df66fc349p-7);  //  1.5379652102610957e-02
    p = fma (p, q, -0x1.1df1ad154a27fp-3);  // -1.3962111684056208e-01
    p = fma (p, q,  0x1.dd2c8b74febf6p-3);  //  2.3299511862555250e-01

    /* Divide (1+p) by (1+2*a) ==> exp(a*a)*erfc(a) */
    d = a + 0.5;
    r = 1.0 / d;
    r = r * 0.5;
    q = fma (p, r, r); // q = (p+1)/(1+2*a)
    t = q + q;
    e = (p - q) + fma (t, -a, 1.0); // residual: (p+1)-q*(1+2*a)
    r = fma (e, r, q);

    /* Handle argument of infinity */
    if (a > 0x1.fffffffffffffp1023) r = 0.0;

    /* Handle negative arguments: erfcx(x) = 2*exp(x*x) - erfcx(|x|) */
    if (x < 0.0) {
        s = x * x;
        d = fma (x, x, -s);
        e = exp (s);
        r = e - r;
        r = fma (e, d + d, r); 
        r = r + e;
        if (e > 0x1.fffffffffffffp1023) r = e; // avoid creating NaN
    }
    return r;
}

 // inner product implementation
 double inner_prod(std::vector<double> const& u, std::vector<double> const& v, int sync){

    double accum = 0;
    #ifdef _OPENMP
    #pragma omp parallel for reduction(+ : accum)
    #endif
    for (int i = 0; i < u.size(); i++) {
        accum += u[i] * v[i];
    }

    double accum_total = 0;
    if (sync == 1){
        //std::cout << "accum = " << accum << ", rank = " << rank << std::endl;
        MPI_Allreduce(&accum, &accum_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    } 
    else {
        accum_total = accum;    
    }

    return accum_total;
 }

 double l2_norm2(std::vector<double> const& u, int sync){
    return inner_prod(u,u, sync);
 }