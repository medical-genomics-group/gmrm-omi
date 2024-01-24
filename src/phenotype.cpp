#include <iostream>
#include <fstream>
#include <regex>
#include <math.h>
#include <immintrin.h>
#include <mpi.h>
#include "utilities.hpp"
#include "phenotype.hpp"
#include <boost/range/algorithm.hpp>
#include <boost/random/uniform_int.hpp>
#include <filesystem>

namespace fs = std::filesystem;


Phenotype::Phenotype(std::string fp, const Options& opt, const int N, const int M) :
    filepath(fp),
    N(N),
    M(M),
    im4(N%4 == 0 ? N/4 : N/4+1),
    K(opt.get_nmixtures()),
    C(opt.get_cov_num()),
    G(opt.get_ngroups()) {

    deltas = std::vector<double> (C, 0.0);
    cov_denom = std::vector<double> (C, 0.0);

    mave = (double*) _mm_malloc(size_t(M) * sizeof(double), 32);
    check_malloc(mave, __LINE__, __FILE__);
    msig = (double*) _mm_malloc(size_t(M) * sizeof(double), 32);
    check_malloc(msig, __LINE__, __FILE__);
    epsilon_ = (double*) _mm_malloc(size_t(N) * sizeof(double), 32);
    check_malloc(epsilon_, __LINE__, __FILE__);
    if(opt.get_model() == "probit"){
        z_ = (double*) _mm_malloc(size_t(N) * sizeof(double), 32);
        check_malloc(z_, __LINE__, __FILE__);
        y_ = (double*) _mm_malloc(size_t(N) * sizeof(double), 32);
        check_malloc(y_, __LINE__, __FILE__);
    }

    if (!opt.predict()) {
        cass = (int*) _mm_malloc(G * K * sizeof(int), 32);
        check_malloc(cass, __LINE__, __FILE__);

        betas.assign(M, 0.0);
        acum.assign(M, 0.0);
        muk.assign(K, 0.0);
        denom.resize(K - 1);
        logl.resize(K);
        comp.assign(M, 0);
        beta_sqn.resize(G);
        m0.resize(G);
        sigmag.resize(G);

        dirich.clear();
        for (int i=0; i<K; i++) dirich.push_back(1.0);
    } else {
        set_prediction_filenames(opt.get_out_dir(), opt.get_infname_base());
    }

    set_output_filenames(opt.get_out_dir());
    read_file(opt);
}

//copy ctor
Phenotype::Phenotype(const Phenotype& rhs) :
    dist_m(rhs.dist_m),
    dist_d(rhs.dist_d),
    filepath(rhs.filepath),
    inbet_fp(rhs.inbet_fp),
    outbet_fp(rhs.outbet_fp),
    outcpn_fp(rhs.outcpn_fp),
    outcsv_fp(rhs.outcsv_fp),
    outmlma_fp(rhs.outmlma_fp),
    nonas(rhs.nonas),
    nas(rhs.nas),
    im4(rhs.im4),
    N(rhs.N),
    M(rhs.M),
    G(rhs.G),
    K(rhs.K),
    betas(rhs.betas),
    data(rhs.data),
    mask4(rhs.mask4),
    midx(rhs.midx),
    denom(rhs.denom),
    muk(rhs.muk),
    logl(rhs.logl),
    acum(rhs.acum),
    comp(rhs.comp),
    beta_sqn(rhs.beta_sqn),
    m0(rhs.m0),
    pi_est(rhs.pi_est),
    dirich(rhs.dirich),
    epssum(rhs.epssum),
    sigmae_(rhs.sigmae_),
    sigmag(rhs.sigmag),
    mu(rhs.mu) {
    epsilon_ = (double*) _mm_malloc(size_t(N) * sizeof(double), 32);
    check_malloc(epsilon_, __LINE__, __FILE__);
    for (int i=0; i<N; i++)
        epsilon_[i] = rhs.epsilon_[i];
    mave = (double*) _mm_malloc(size_t(M) * sizeof(double), 32);
    check_malloc(mave, __LINE__, __FILE__);
    msig = (double*) _mm_malloc(size_t(M) * sizeof(double), 32);
    check_malloc(msig, __LINE__, __FILE__);
    for (int i=0; i<M; i++) {
        mave[i] = rhs.mave[i];
        msig[i] = rhs.msig[i];
    }
    cass = (int*) _mm_malloc(size_t(K * G) * sizeof(int), 32);
    check_malloc(cass, __LINE__, __FILE__);
}

int Phenotype::get_m0_sum() {
    int sum = 0;
    for (int i=0; i<m0.size(); i++)
        sum += m0.at(i);
    return sum;
}

// Set input filenames based on input phen file (as per output)
void Phenotype::set_prediction_filenames(const std::string out_dir, const std::string in_fname_base) {
    fs::path pphen = filepath;
    fs::path base  = out_dir;
    base /= pphen.stem();

    inbet_fp = out_dir + "/" + in_fname_base + ".bet";
    std::cout << "inbet_fp = " << inbet_fp << std::endl;

    fs::path pmlma = base;
    pmlma += ".mlma";
    outmlma_fp = pmlma.string();

    fs::path pyest = base;
    pyest += ".yest";
    outyest_fp = pyest.string();
    std::cout << outyest_fp << std::endl;

    incov_fp = out_dir + "/" + in_fname_base + "_cov.csv";
    std::cout << "incov_fp = " << incov_fp << std::endl;
}

void Phenotype::set_output_filenames(const std::string out_dir) {
    fs::path pphen = filepath;
    fs::path base  = out_dir;
    base /= pphen.stem();
    fs::path pbet = base;
    pbet += ".bet";
    fs::path pcpn = base;
    pcpn += ".cpn";
    fs::path pcsv = base;
    pcsv += ".csv";
    fs::path pcov = base;
    pcov += "_cov.csv";

    outbet_fp = pbet.string();
    outcpn_fp = pcpn.string();
    outcsv_fp = pcsv.string();
    outcov_fp = pcov.string();
}

// Input and output
void Phenotype::open_prediction_files() {

    check_mpi(MPI_File_open(MPI_COMM_WORLD,
                            get_inbet_fp().c_str(),
                            MPI_MODE_RDONLY,
                            MPI_INFO_NULL,
                            get_inbet_fh()),
              __LINE__, __FILE__);

    check_mpi(MPI_File_open(MPI_COMM_WORLD,
                            get_outmlma_fp().c_str(),
                            MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL,
                            MPI_INFO_NULL,
                            get_outmlma_fh()),
              __LINE__, __FILE__);
}

void Phenotype::close_prediction_files() {
    check_mpi(MPI_File_close(get_inbet_fh()), __LINE__, __FILE__);
    check_mpi(MPI_File_close(get_outmlma_fh()), __LINE__, __FILE__);
}

void Phenotype::delete_output_prediction_files() {
    MPI_File_delete(get_outmlma_fp().c_str(), MPI_INFO_NULL);
}

void Phenotype::open_output_files() {
    check_mpi(MPI_File_open(MPI_COMM_WORLD,
                            get_outbet_fp().c_str(),
                            MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL,
                            MPI_INFO_NULL,
                            get_outbet_fh()),
              __LINE__, __FILE__);
    check_mpi(MPI_File_open(MPI_COMM_WORLD,
                            get_outcpn_fp().c_str(),
                            MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL,
                            MPI_INFO_NULL,
                            get_outcpn_fh()),
              __LINE__, __FILE__);
    check_mpi(MPI_File_open(MPI_COMM_WORLD,
                            get_outcsv_fp().c_str(),
                            MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL,
                            MPI_INFO_NULL,
                            get_outcsv_fh()),
              __LINE__, __FILE__);
    check_mpi(MPI_File_open(MPI_COMM_WORLD,
                            get_outcov_fp().c_str(),
                            MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL,
                            MPI_INFO_NULL,
                            get_outcov_fh()),
              __LINE__, __FILE__);
}

void Phenotype::close_output_files() {
    check_mpi(MPI_File_close(get_outbet_fh()), __LINE__, __FILE__);
    check_mpi(MPI_File_close(get_outcpn_fh()), __LINE__, __FILE__);
    check_mpi(MPI_File_close(get_outcsv_fh()), __LINE__, __FILE__);
    check_mpi(MPI_File_close(get_outcov_fh()), __LINE__, __FILE__);
}

void Phenotype::delete_output_files() {
    MPI_File_delete(get_outbet_fp().c_str(), MPI_INFO_NULL);
    MPI_File_delete(get_outcpn_fp().c_str(), MPI_INFO_NULL);
    MPI_File_delete(get_outcsv_fp().c_str(), MPI_INFO_NULL);
    MPI_File_delete(get_outcov_fp().c_str(), MPI_INFO_NULL);
}

void Phenotype::print_cass() {
    printf("INFO   : cass for phenotype [...]\n");
    for (int i=0; i<G; i++) {
        printf("         %2d : ", i);
        for (int j=0; j<K; j++) {
            printf("%7d", cass[i*K+j]);
        }
        printf("\n");
    }
}

void Phenotype::print_cass(const std::vector<int>& mtotgrp) {
    printf("INFO   : cass for phenotype [...]\n");
    for (int i=0; i<G; i++) {
        printf("         group %2d: %7d | cass: ", i, mtotgrp.at(i));
        for (int j=0; j<K; j++) {
            printf("%7d", cass[i*K + j]);
        }
        printf("\n");
    }
}

void Phenotype::update_pi_est_dirichlet(const int g) {
    std::vector<double> tmp;
    double sum = 0.0;
    for (int i=0; i<K; i++) {
        double val = dist_d.rgamma((double)cass[g * K + i] + dirich[i], 1.0);
        set_pi_est(g, i, val);
        sum += val;
    }
    for (int i=0; i<K; i++)
        set_pi_est(g, i, get_pi_est(g, i) / sum);
}

double Phenotype::epsilon_sum() {
    double* epsilon = get_epsilon();
    double sum = 0.0;
#ifdef _OPENMP
#pragma omp parallel for simd aligned(epsilon:32) reduction(+:sum)
#endif
    for (int i=0; i<N; i++) {
        sum += epsilon[i];
    }
    return sum;
}

double Phenotype::epsilon_sumsqr() {
    double* epsilon = get_epsilon();
    double sumsqr = 0.0;
#ifdef _OPENMP
#pragma omp parallel for simd aligned(epsilon:32) reduction(+:sumsqr)
#endif
    for (int i=0; i<N; i++) {
        sumsqr += epsilon[i] * epsilon[i];
    }
    return sumsqr;
}

void Phenotype::increment_beta_sqn(const int group, const double val) {
    if(val != val) std::cout << "NAN" << std::endl;
    beta_sqn.at(group) += val;
}

int Phenotype::get_marker_local_index(const int shuff_idx) {
    return midx[shuff_idx];
}

double Phenotype::sample_inv_scaled_chisq_rng(const double a, const double b) {
    return dist_d.inv_scaled_chisq_rng(a, b);
}

double Phenotype::sample_norm_rng(const double a, const double b) {
    return dist_d.norm_rng(a, b);
}

double Phenotype::sample_norm_rng() {
    //printf("sampling mu with epssum = %20.15f and sigmae = %20.15f; nonas = %d\n", epssum, sigmae, nonas);
    return dist_d.norm_rng(epssum / double(nonas), get_sigmae() / double(nonas));
}

double Phenotype::sample_trunc_norm_rng(const double a, const double b, const double c) {
    return dist_d.trunc_norm_rng(a, b, c);
}

double Phenotype::sample_beta_rng(const double a, const double b) {
    return dist_d.beta_rng(a, b);
}

double Phenotype::sample_unif_rng() {
    return dist_d.unif_rng();
}

void Phenotype::sample_for_free(const int n) {
    for (int i=0; i<n; i++)
        double fake = dist_d.unif_rng();
}

//void Phenotype::sample_sigmag_beta_rng() {
//    sigmag = dist.beta_rng(epssum / double(nonas), sigmae / double(nonas));
//}

void Phenotype::set_prng_m(const unsigned int seed) {
    dist_m.set_prng(seed);
}
void Phenotype::set_prng_d(const unsigned int seed) {
    dist_d.set_prng(seed);
}

void Phenotype::set_midx() {
    midx.clear();
    for (int i=0; i<M; ++i)
        midx.push_back(i);
}

void Phenotype::shuffle_midx(const bool mimic_hydra) {
    boost::uniform_int<> unii(0, M-1);
    if (mimic_hydra) {
        boost::variate_generator< boost::mt19937&, boost::uniform_int<> > generator(dist_d.get_rng(), unii);
        boost::range::random_shuffle(midx, generator);
    } else {
        boost::variate_generator< boost::mt19937&, boost::uniform_int<> > generator(dist_m.get_rng(), unii);
        boost::range::random_shuffle(midx, generator);
    }
}

// Set latent variable to 0 for all individuals
void Phenotype::init_latent(){
    double* z = get_z();

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i=0; i<N; i++) 
        z[i] = 0.0;
}

// Update latent variable based on current marker and effect
void Phenotype::update_latent(const int mloc, const double* meth) {
    double* z = get_z();
    double beta = get_marker_beta(mloc);
    double mave = get_marker_ave(mloc);
    double msig = get_marker_sig(mloc);

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i=0; i<N; i++) {
        double val = (meth[i] - mave) * msig;
        z[i] += val * beta;
    } 
}

void Phenotype::offset_latent(const double offset) {
    double* z = get_z();

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i=0; i<N; i++) {
        z[i] += offset;
    }
}

// Update latent variable based on current covariate effect delta
void Phenotype::update_latent_cov(const int covi, double delta) {
    double* z = get_z();
    for (int i=0; i<N; i++) {
        double offset = Z_[i][covi] * delta;
        z[i] += offset;
    }
}

// Dot product for covariates
double  Phenotype::dot_product_cov(int covi){
    double* epsilon = get_epsilon();
    double Ze = 0.0;
    for (int i=0; i<N; i++) {
        Ze += epsilon[i] * Z_[i][covi];
    }
    return Ze;
}

// Load covariate effects
void Phenotype::load_cov_deltas(){

    double ts = MPI_Wtime();

    std::ifstream covf(incov_fp);
    std::string line; 
    std::regex re("\\s+");

    int line_i = 0;

    while (std::getline(covf, line)) // read the current line
    {
        int Cobs = 0;
        std::vector<double> entries;
        std::sregex_token_iterator iter(line.begin(), line.end(), re, -1);
        std::sregex_token_iterator re_end;
        
        ++iter; // skip line number
        ++iter; // skip iteration
        ++iter; // skip number of covariates
        for ( ; iter != re_end; ++iter){
            entries.push_back(std::stod(*iter));
            Cobs++;
        }

        if (Cobs != C){
            printf("FATAL   : number of covariate deltas = %d does not match to the specified number of covariates = %d.\n", Cobs, C);
            exit(EXIT_FAILURE);
        }
        
        deltas_it.push_back(entries); 
        line_i++;   
    }

    printf("INFO   : Number of loaded lines from _cov.csv file = %d \n", line_i);

    double te = MPI_Wtime();

    printf("INFO   : time to load %d covariate effects = %.2f seconds.\n", C, te - ts);
}

// Sample artificial target from truncated normal and update residual 
void Phenotype::init_epsilon(){
    double* z = get_z();
    double* epsilon = get_epsilon();
    double* y = get_y();

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i=0; i<N; i++) {
        epsilon[i] = sample_trunc_norm_rng(z[i], 1.0, y[i]) - z[i];     
    }
}

// Add covariate contributions to base epsilon
void Phenotype::update_epsilon_cov(const int covi, double delta) {
    double* epsilon = get_epsilon();

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i=0; i<N; i++) {
        epsilon[i] += delta * Z_[i][covi];    
    }
}

// Add contributions to base epsilon
void Phenotype::update_epsilon(const double* dbeta, const double* meth) {

    double* epsilon = get_epsilon();
    const double db = dbeta[0];
    const double sig_inv = dbeta[2];
    const double ave = dbeta[1];
            
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i=0; i<N; i++){
        epsilon[i] += (meth[i] - ave) * db * sig_inv;
    }
}

void Phenotype::offset_epsilon(const double offset) {

    double* epsilon = get_epsilon();

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i=0; i<N; i++)
        epsilon[i] += offset;
}

/*
void Phenotype::update_epsilon_sum() {
    __m256d eps4, sig4, lutna;
    __m256d sums = _mm256_set1_pd(0.0);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(addpd4:sums)
#endif
    for (int i=0; i<im4; i++) {
        eps4  = _mm256_load_pd(&epsilon[i*4]);
        lutna = _mm256_load_pd(&na_lut[mask4[i] * 4]);
        eps4  = _mm256_mul_pd(eps4, lutna);
        sums  = _mm256_add_pd(sums, eps4);
    }
    epssum = sums[0] + sums[1] + sums[2] + sums[3];
}
*/

// Only depends on NAs
void Phenotype::update_epsilon_sigma() {
    double* epsilon = get_epsilon();

    double sigmae = 0.0;
#ifdef _OPENMP
#pragma omp parallel for simd reduction(+:sigmae)
#endif
    for (int i=0; i<N; i++) {

            sigmae += epsilon[i] * epsilon[i];
        
    }
    set_sigmae(sigmae / double(nonas) * 0.5);
}


// Compute mean and associated standard deviation for markers

void PhenMgr::compute_markers_statistics(const double* meth_data, const int N, const int M) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    Phenotype& phen = get_phens()[0];

    double* mave = phen.get_mave();
    double* msig = phen.get_msig();

    double start = MPI_Wtime();

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i=0; i<M; i++) {
        size_t methix = size_t(i) * size_t(N);
        const double* methm = &meth_data[methix];
        double suma = 0.0;

#ifdef _OPENMP
#pragma omp simd reduction(+:suma)
#endif
        for (int j = 0; j < N; j++) 
            suma += methm[j];

        //calculating vector of marker precision
        mave[i] = suma / double( N );
            
        double sumsqr = 0.0;

#ifdef _OPENMP
#pragma omp simd reduction(+:sumsqr)
#endif
        for (int j=0; j<N; j++) {
            double val = (methm[j] - mave[i]);
            sumsqr += val * val;
        }
            
        if (sumsqr != 0.0)
            msig[i] = 1.0 / sqrt(sumsqr / (double( N ) - 1.0));
        else 
            msig[i] = 1.0;
    }    
    
        double end = MPI_Wtime();

        std::cout << "rank = " << rank << ": statistics took " << end - start << " seconds to run." << std::endl;
        //std::cout << "msig[0] =" << msig[0] << ", msig[1] =" << msig[1] << ", msig[2] =" << msig[2] << std::endl;
        //std::cout << "mave[0] =" << mave[0] << ", mave[1] =" << mave[1] << ", mave[2] =" << mave[2] << std::endl;
        
}

void PhenMgr::display_markers_statistics(const int n) {
    for (auto& phen : get_phens()) {
        phen.print_info();
        for (int i=0; i<n; i++) {
            printf("avg for marker %d = %20.15f +/- %20.15f  1/sig = %20.15f\n", i, phen.get_mave()[i], phen.get_msig()[i], 1.0 / phen.get_msig()[i]);
        }
    }
}

void PhenMgr::print_info() {
    for (auto& phen : get_phens())
        phen.print_info();
}

void PhenMgr::read_phen_files(const Options& opt, const int N, const int M) {
    std::vector<std::string> phen_files = opt.get_phen_files();
    int pi = 0;
    for (auto fp = phen_files.begin(); fp != phen_files.end(); ++fp) {
        pi++;
        if (opt.verbosity_level(3))
            std::cout << "Reading phenotype file " << pi << ": " << *fp << std::endl;
        phens.emplace_back(*fp, opt, N, M);
        //std::cout << "emplace_back for pi " << pi << ": " << *fp << std::endl;
    }
}


// Read phenotype file assuming PLINK format:
// Family ID, Individual ID, Phenotype; One row per individual
void Phenotype::read_file(const Options& opt) {

    std::ifstream infile(filepath);
    std::string line;
    std::regex re("\\s+");

    double* epsilon = get_epsilon();
    double* y = get_y();
    double sum = 0.0;

    if (infile.is_open()) {
        int line_n = 0;
        nonas = 0, nas = 0;
        while (getline(infile, line)) {
            //std::cout << line << std::endl;
            int m4 = line_n % 4;
            if (m4 == 0)  mask4.push_back(0b00001111);

            std::sregex_token_iterator first{line.begin(), line.end(), re, -1}, last;
            std::vector<std::string> tokens{first, last};

            if (tokens[2] == "NA") {
                nas += 1;
                data.push_back(std::numeric_limits<double>::max());
                if (opt.verbosity_level(2)) {
                    std::cout << " ... found NA on line " << line_n << ", m4 = " << m4 << " on byte " << int(line_n / 4) << std::endl;
                    fflush(stdout);
                }
                mask4.at(int(line_n / 4)) &= ~(0b1 << m4);
            } else {
                nonas += 1;
                data.push_back(atof(tokens[2].c_str()));
                sum += atof(tokens[2].c_str());
            }

            line_n += 1;

            if (opt.verbosity_level(3)) {
                if (line_n % 4 == 0 && line_n > 3 && line_n < 30) {
                    std::cout << "mask4[" << int(line_n / 4) - 1 << "] = " << unsigned(mask4.at(int(line_n / 4) - 1)) << std::endl;
                }
            }
        }
        infile.close();

        assert(nas + nonas == N);

        if (opt.get_model()=="probit"){
            // Init epsilon as a sample from truncated normal
            for (int i=0; i<N; i++) {
                y[i] = double(data[i]);
                epsilon[i] = sample_trunc_norm_rng(0.0, 1.0, y[i]);
            }
        } else {
            // Center and scale
            double avg = sum / double(nonas);
            if (opt.verbosity_level(3))
                printf("phen avg = %20.15f\n", avg);

            double sqn = 0.0;
            for (int i=0; i<data.size(); i++) {
                if (opt.verbosity_level(3) && i < 10)
                    std::cout << data[i] - avg  << std::endl;
                if (data[i] == std::numeric_limits<double>::max()) {
                    epsilon[i] = 0.0;
                } else {
                    epsilon[i] = data[i] - avg;
                    sqn += epsilon[i] * epsilon[i];
                }
            }
            sqn = sqrt(double(nonas-1) / sqn);
            if (opt.verbosity_level(3))
                printf("phen sqn = %20.15f\n", sqn);
            for (int i=0; i<data.size(); i++)
                epsilon[i] *= sqn;
        }    
    } else {
        std::cout << "FATAL: could not open phenotype file: " << filepath << std::endl;
        exit(EXIT_FAILURE);
    }
}

void Phenotype::print_info() const {
    printf("INFO   : %s has %d NAs and %d non-NAs.\n", get_filepath().c_str(), nas, nonas);
}

void Phenotype::set_nas_to_zero(double* y, const int N) {

}
