#include <iostream>
#include <mpi.h>
#include "options.hpp"
#include "phenotype.hpp"
#include "bayes.hpp"
#include "dimensions.hpp"

int main(int argc, char *argv[]) {

    MPI_Init(NULL, NULL);

    const Options    opt(argc, argv);
    const Dimensions dims(opt);
    BayesRR brr(opt, dims);
    if (opt.infer()) {
        brr.infer();
    }
    if (opt.test()){
        brr.test();
    } 
    if (opt.predict()){
        brr.predict();
    }

    MPI_Finalize();

    return 0;
}
