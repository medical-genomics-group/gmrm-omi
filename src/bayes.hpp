#pragma once
#include <iostream>
#include <typeinfo>
#include <map>
#include "options.hpp"
#include "phenotype.hpp"
#include "dimensions.hpp"


class Bayes {

public:

    const double V0E  = 0.01;//0.0001;
    const double S02E = 0.01;//0.0001;
    const double V0G  = 0.01;//0.0001;
    const double S02G = 0.01;//0.0001;


    Bayes(const Options& opt, const Dimensions& dims) : opt(opt),
                                                        rank(dims.get_rank()),
                                                        nranks(dims.get_nranks()),
                                                        N(dims.get_nt()),
                                                        Mt(dims.get_mt()),
                                                        K(opt.get_nmixtures()),
                                                        G(opt.get_ngroups()),
                                                        cva(opt.get_cva()), 
                                                        cvai(opt.get_cvai()) {

        M = Mt;
        S = 0;
        Mm = Mt;

        if (!opt.predict()) {
            
            pi_prior.resize(G);
            mtotgrp.resize(G);

            for (int i=0 ; i<G; i++) {
                mtotgrp.at(i) = 0;
                double sum_cva = 0.0;
                for (int j=0; j<K-1; j++) {
                    sum_cva += cva[i][j+1];
                }
                pi_prior[i].resize(K, 0.5);
                for (int j=1; j<K; j++) {
                    pi_prior[i][j] =  pi_prior[i][0] * cva[i][j] / sum_cva;
                }
            }        
        }

        setup_processing();
    }

    ~Bayes() {
        //std::cout << "## calling Bayes dtor" << std::endl;
        if (meth_data != nullptr)  _mm_free(meth_data);
    }

    void predict();
    void process();

    double dot_product(const int mrki, double* __restrict__ phen, const double mu, const double sigma);
    void list_phen_files() const { opt.list_phen_files(); }
    int  get_N()  { return N;  } // Invariant over tasks
    int  get_M()  { return M;  } // Number of markers processed by task
    int  get_Mt() { return Mt; } // Total number of markers, sum over tasks
    int  get_Mm() { return Mm; } // Maximum number of markers per task (others may have M + 1)
    int  get_K()  { return K;  }
    void shuffle_markers();
    int  get_marker_group(const int mglob) { 
        return group_index.at(mglob);
    }
    void update_epsilon(const int* counts, const double* dbetas, const double* recv_meth);
    void check_openmp();
    void print_cva();
    void print_cvai();

private:
    const Options opt;
    PhenMgr pmgr;
    const int N = 0;
    const int Mt = 0;
    const int rank = 0;
    const int nranks = 0;
    double* meth_data = nullptr;
    const int K = 0;
    const int G = 0;
    
    std::vector<int> mtotgrp;
    std::vector<int> group_index;
    const std::vector<std::vector<double>> cva;
    const std::vector<std::vector<double>> cvai;
    std::vector<std::vector<double>> pi_prior;

    std::vector<std::string>   rsid;
    std::map<std::string, int> m_refrsid;

    std::vector<std::vector<double>> Z; //Covariates
    std::vector<double> delta; // Covariate effects

    int S = 0;              // task marker start 
    int M = 0;              // task marker length
    int Mm = 0;
    size_t mbytes = 0;

    void check_options();
    void setup_processing();

    void load_methylation_data();
    void load_covariates();
    void check_processing_setup();
    void read_group_index_file(const std::string& file);
};


class BayesRR : public Bayes {

public:
    
    BayesRR(const Options& opt, 
            const Dimensions& dims) : Bayes(opt,
                                            dims) { 
        //std::cout << "calling BayesRR constructor" << std::endl;
    }
};
