#include <iostream>
#include <fstream>
#include <iterator>
#include <limits.h>
#include <cmath>
#include <immintrin.h>
#include <omp.h>
#include "bayes.hpp"
#include "utilities.hpp"
#include "xfiles.hpp"
#include <boost/math/special_functions/gamma.hpp>
#include <regex>

void Bayes::predict() {

    double ts = MPI_Wtime();

    check_openmp();

    MPI_Status status;
    MPI_Offset file_size = 0;
    MPI_File*  fh;

    int C = opt.get_cov_num();

    Phenotype& phen = pmgr.get_phens()[0];

    phen.delete_output_prediction_files();
    phen.open_prediction_files();

    fh = phen.get_inbet_fh();
    check_mpi(MPI_File_get_size(*fh, &file_size), __LINE__, __FILE__);
    printf("file_size = %u B\n", file_size);
    std::cout << "file name: " << phen.get_inbet_fp() << std::endl;

    // First element of the .bet is the total number of processed markers
    // Then: iteration (uint) beta (double) for all markers
    uint Mtot_ = 0;
    MPI_Offset betoff = size_t(0);
    check_mpi(MPI_File_read_at_all(*fh, betoff, &Mtot_, 1, MPI_UNSIGNED, &status), __LINE__, __FILE__);
    if (Mtot_ != Mt) {
        printf("Mismatch between expected and Mtot read from .bet file: %d vs %d\n", Mt, Mtot_);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    assert((file_size - sizeof(uint)) % (Mtot_ * sizeof(double) + sizeof(uint)) == 0);
    uint niter = (file_size - sizeof(uint)) / (Mtot_ * sizeof(double) + sizeof(uint));
    if (rank == 0)
        printf("INFO   : Number of recorded iterations in .bet file: %u\n", niter);

    double* beta_sum = (double*) _mm_malloc(size_t(Mtot_) * sizeof(double), 32);
    check_malloc(beta_sum, __LINE__, __FILE__);
    for (int i=0; i<Mtot_; i++) beta_sum[i] = 0.0;

    double* beta_it = (double*) _mm_malloc(size_t(Mtot_) * sizeof(double), 32);
    check_malloc(beta_it, __LINE__, __FILE__);

    uint start_iter = opt.get_burn_in();

    for (uint i=start_iter; i<niter; i++) {
        betoff
            = sizeof(uint) // Mtot
            + (sizeof(uint) + size_t(Mtot_) * sizeof(double)) * size_t(i)
            + sizeof(uint);
        check_mpi(MPI_File_read_at_all(*fh, betoff, beta_it, Mtot_, MPI_DOUBLE, &status), __LINE__, __FILE__);
        for (int j=0; j<Mtot_;j++)
            beta_sum[j] += beta_it[j];
    }

    for (int j=0; j<Mtot_;j++)
        beta_sum[j] /= double(niter);
        /*
        if(C > 0){
            phen.load_cov_deltas();
            phen.set_Z(Z);
            phen.avg_deltas_it(niter);
        }
        std::vector<double>* deltas = phen.get_deltas();

        double* c = (double*) _mm_malloc(size_t(N) * sizeof(double), 32);
        check_malloc(c, __LINE__, __FILE__);
        
        for (int i=0; i<N; i++){
            c[i] = 0.0;
            for(int covi = 0; covi<C; covi++){
                c[i] += deltas->at(covi) * Z[i][covi];
            }
        }*/
        
    fflush(stdout);

    double* g = (double*) _mm_malloc(size_t(N) * sizeof(double), 32);
    check_malloc(g, __LINE__, __FILE__);

    for (int i = 0; i < N; i++) g[i] = 0.0;

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int mrki = 0; mrki < M; mrki++) {

        size_t methix = size_t(mrki) * size_t(N);
        const double* methm = &meth_data[methix];

        double mave = phen.get_marker_ave(mrki);
        double msig = phen.get_marker_sig(mrki);

        for (int i = 0; i < N; i++) {
            double val = (methm[i] - mave);

#ifdef _OPENMP
#pragma omp atomic update
#endif
            g[i] += val * beta_sum[mrki];
        }
    }

    fflush(stdout);

    std::ofstream yest_stream;
    yest_stream.open(phen.get_outyest_fp());
    for (int i=0; i<N; i++){
        double z = g[i];// + c[i];
        yest_stream << z << std::endl;
    }
    yest_stream.close();

    double te = MPI_Wtime();

    if (rank == 0)
        printf("INFO   : Time to compute the predictions: %.2f seconds.\n", te - ts);
}

void Bayes::test() {

    double ts = MPI_Wtime();

    check_openmp();

    MPI_Status status;
    MPI_Offset file_size = 0;
    MPI_File*  fh;

    Phenotype& phen = pmgr.get_phens()[0];

    phen.delete_output_test_files();
    phen.open_test_files();

    fh = phen.get_inbet_fh();
    check_mpi(MPI_File_get_size(*fh, &file_size), __LINE__, __FILE__);
    printf("file_size = %u B\n", file_size);
    std::cout << "file name: " << phen.get_inbet_fp() << std::endl;

    // First element of the .bet is the total number of processed markers
    // Then: iteration (uint) beta (double) for all markers
    uint Mtot_ = 0;
    MPI_Offset betoff = size_t(0);
    check_mpi(MPI_File_read_at_all(*fh, betoff, &Mtot_, 1, MPI_UNSIGNED, &status), __LINE__, __FILE__);
    if (Mtot_ != Mt) {
        printf("Mismatch between expected and Mtot read from .bet file: %d vs %d\n", Mt, Mtot_);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    assert((file_size - sizeof(uint)) % (Mtot_ * sizeof(double) + sizeof(uint)) == 0);
    uint niter = (file_size - sizeof(uint)) / (Mtot_ * sizeof(double) + sizeof(uint));
    if (rank == 0)
        printf("INFO   : Number of recorded iterations in .bet file: %u\n", niter);

    double* beta_p = (double*) _mm_malloc(size_t(Mtot_) * sizeof(double), 32);
    check_malloc(beta_p, __LINE__, __FILE__);
    for (int i=0; i<Mtot_; i++) beta_p[i] = 0.0;

    double* beta_sum = (double*) _mm_malloc(size_t(Mtot_) * sizeof(double), 32);
    check_malloc(beta_sum, __LINE__, __FILE__);
    for (int i=0; i<Mtot_; i++) beta_sum[i] = 0.0;

    double* beta_it = (double*) _mm_malloc(size_t(Mtot_) * sizeof(double), 32);
    check_malloc(beta_it, __LINE__, __FILE__);

    uint start_iter = opt.get_burn_in();

    for (uint i=start_iter; i<niter; i++) {
        betoff
            = sizeof(uint) // Mtot
            + (sizeof(uint) + size_t(Mtot_) * sizeof(double)) * size_t(i)
            + sizeof(uint);
        check_mpi(MPI_File_read_at_all(*fh, betoff, beta_it, Mtot_, MPI_DOUBLE, &status), __LINE__, __FILE__);
        for (int j=0; j<Mtot_;j++){
            beta_sum[j] += beta_it[j];
            if(beta_it[j] != 0.0)  
                beta_p[j] += 1.0;
        }
    }

    MPI_File* mlma_fh = phen.get_outmlma_fh();

    for (int j=0; j<Mtot_;j++){
        beta_p[j] /= double(niter);
        beta_sum[j] /= double(niter);

        char buff[LENBUF];

        int cx = snprintf(buff, LENBUF, "%8d %20.15f %20.15f\n", j, beta_sum[j], beta_p[j]);
        assert(cx >= 0 && cx < LENBUF);

        MPI_Offset offset = size_t(j) * strlen(buff);
        check_mpi(MPI_File_write_at(*mlma_fh, offset, &buff, strlen(buff), MPI_CHAR, &status), __LINE__, __FILE__);
    }

    if (rank == 0)
        printf("INFO   : Association testing results soterd in .mlma file.\n"); 

    double* g = (double*) _mm_malloc(size_t(N) * sizeof(double), 32);
    check_malloc(g, __LINE__, __FILE__);
    for (int i = 0; i < N; i++) g[i] = 0.0;

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int mrki = 0; mrki < M; mrki++) {

        size_t methix = size_t(mrki) * size_t(N);
        const double* methm = &meth_data[methix];

        double mave = phen.get_marker_ave(mrki);
        double msig = phen.get_marker_sig(mrki);

        for (int i = 0; i < N; i++) {
            double val = (methm[i] - mave);

#ifdef _OPENMP
#pragma omp atomic update
#endif
            g[i] += val * beta_sum[mrki];
            
        }
    }

    std::ofstream yest_stream;
    yest_stream.open(phen.get_outyest_fp());
    for (int i=0; i<N; i++){
        yest_stream << g[i] << std::endl;
    }
    yest_stream.close();

    if (rank == 0)
        printf("INFO   : Predictions stored in .yest file.\n"); 

    double* y = phen.get_epsilon();
    double y_mean = 0.0;
    double SSres = 0.0;
    double SStot = 0.0;

    for (int i=0; i<N; i++){
        SSres += (y[i] - g[i]) * (y[i] - g[i]);
        y_mean += y[i];
    }
    y_mean /= double(N);
    for (int i=0; i<N; i++){
        SStot += (y[i] - y_mean) * (y[i] - y_mean);
    }
    double r2 = 1.0 - SSres / SStot;

    if (rank == 0)
        printf("INFO   : R2 = %.4f \n", r2);

    fflush(stdout);

    double te = MPI_Wtime();

    if (rank == 0)
        printf("INFO   : Time to compute the association testing: %.2f seconds.\n", te - ts);
}

void Bayes::process() {

    double ts_overall = MPI_Wtime();

    check_openmp();

    Phenotype& phen = pmgr.get_phens()[0];

    phen.delete_output_files();
    phen.open_output_files();
    phen.set_midx();
    for (int i=0; i<opt.get_ngroups(); i++) {
        phen.set_sigmag_for_group(i, phen.sample_beta_rng(1.0, 1.0));
        if (mtotgrp.at(i) == 0)
            phen.set_sigmag_for_group(i, 0.0);
    }

    phen.set_pi_est(pi_prior);

    // number of covariates
    int C = opt.get_cov_num();

    for (unsigned int it = 1; it <= opt.get_iterations(); it++) {

        double ts_it = MPI_Wtime();

        printf("\n\n@@@ ITERATION %5d\n", it);

        std::string model = opt.get_model();

        if (model == "probit"){
            // Init residual based on current latent variable 
            phen.init_epsilon();

            double* z_ = phen.get_z();

            double z_mean = 0.0;
            double z_std = 0.0;

            for (int i=0; i<N; i++) {
                z_mean += z_[i];
            }
            z_mean /= double(N);

            for (int i=0; i<N; i++) {
                z_std += (z_[i] - z_mean) * (z_[i] - z_mean);
            }
            z_std = sqrt(z_std / double(N));
            printf("z_mean = %0.4f, z_std = %0.4f\n", z_mean, z_std);
                
            
            // Init latent to 0
            phen.init_latent();
            
        }

        phen.offset_epsilon(phen.get_mu());

        if (it == 1){
            phen.update_epsilon_sigma();
            printf("epssum = %20.15f, sigmae = %20.15f\n", phen.get_epsilon_sum(), phen.get_sigmae());
        }
        phen.set_mu(phen.sample_norm_rng());
        phen.offset_epsilon(-phen.get_mu());

        if(model=="probit"){
            phen.offset_latent(phen.get_mu());
        }

        // Shuffling of the markers on its own PRNG (see README/wiki)
        if (opt.shuffle_markers()){
            phen.shuffle_midx(opt.mimic_hydra());
        }

        phen.reset_m0();
        phen.reset_cass();
        
        fflush(stdout);

        // Covariates
        if(C > 0){
            for(int covi = 0; covi < C; covi++){
                // Update epsilon with respect to previous covariate effect
                double delta = phen.get_cov_delta(covi);
                phen.update_epsilon_cov(covi, delta);

                double cov_num = phen.dot_product_cov(covi);
                double cov_denom = phen.get_cov_denom(covi);
                double delta_new = 0.0;

                if(cov_denom > 0){
                    // Sample new covariate effect delta
                    delta_new = phen.sample_norm_rng(cov_num / cov_denom, 1.0 / cov_denom);
                }

                phen.set_cov_delta(covi, delta_new);
                //printf("New delta[%d] = %0.6f\n", covi, delta_new);
                //fflush(stdout);

                // Update epsilon with respect to the new covariate effect
                phen.update_epsilon_cov(covi, -delta_new);
            }
        }
        //fflush(stdout);

        double dbetas[3]; // [ dbeta:mave:msig ]

        for (int mrki=0; mrki<Mm; mrki++) {

            bool share_mrk = false;

            for (int i=0; i<3; i++){
                dbetas[i] = 0.0;
            }

            const int mloc = pmgr.get_phens()[0].get_marker_local_index(mrki);
            const int mglo = mloc;
            const int mgrp = get_marker_group(mglo);
            //std::cout << "mloc = " << mloc << ", mglo = " << mglo << ", mgrp = " << mgrp << std::endl;

            if (phen.get_sigmag_for_group(mgrp) == 0.0) {
                phen.set_marker_acum(mloc, 1.0);
                phen.set_marker_beta(mloc, 0.0);
                continue;
            }

            double beta   = phen.get_marker_beta(mloc);
            double sige_g = phen.get_sigmae() / phen.get_sigmag_for_group(get_marker_group(mglo));
            double sigg_e = 1.0 / sige_g;
            double inv2sige = 1.0 / (2.0 * phen.get_sigmae());

            std::vector<double> denom = phen.get_denom();
            std::vector<double> muk   = phen.get_muk();
            std::vector<double> logl  = phen.get_logl();

            for (int i=1; i<=K-1; ++i) {
                denom.at(i-1) = (double)(N - 1) + sige_g * cvai[mgrp][i];
                //printf("it %d, rank %d, m %d: denom[%d] = %20.15f, cvai = %20.15f\n", it, rank, mloc, i-1, denom.at(i-1), cvai[mgrp][i]);
            }
            double ts_debug_dp = MPI_Wtime();
            double num = dot_product(mloc, phen.get_epsilon(), phen.get_marker_ave(mloc), phen.get_marker_sig(mloc));

            //printf("num = %20.15f\n", num);
            num += beta * double(phen.get_nonas() - 1);

            //printf("i:%d r:%d m:%d: num = %.17g, %20.15f, %20.15f\n", it, rank, mloc, num, phen.get_marker_ave(mloc), phen.get_marker_sig(mloc));

            for (int i=1; i<=K-1; ++i){
                muk.at(i) = num / denom.at(i - 1);
            }

            for (int i=0; i<K; i++){
                logl[i] = log(phen.get_pi_est(mgrp, i));
                if (i>0){
                    logl[i] += -0.5 * log(sigg_e * double(phen.get_nonas() - 1) * cva[mgrp][i] + 1.0) + muk[i] * num * inv2sige;
                }
                //printf("logl[%d] = %20.15f\n", i, logl[i]);
            }

            double prob = phen.sample_unif_rng();

            bool zero_acum = false;
            double tmp1 = 0.0;
            for (int i=0; i<K; i++) {
                if (abs(logl[i] - logl[0]) > 700.0){
                    zero_acum = true;
                }
                tmp1 += exp(logl[i] - logl[0]);
            }
            zero_acum ? tmp1 = 0.0 : tmp1 = 1.0 / tmp1;
            phen.set_marker_acum(mloc, tmp1);
            //printf("i:%d r:%d m:%d p:%d  num = %20.15f, acum = %20.15f, prob = %20.15f\n", it, rank, mloc, pheni, num, phen.get_marker_acum(mloc), prob);
            //fflush(stdout);
            double dbeta = phen.get_marker_beta(mloc);

            for (int i=0; i<K; i++) {
                if (prob <= phen.get_marker_acum(mloc) || i == K - 1) {
                    if (i == 0) {
                        phen.set_marker_beta(mloc, 0.0);
                        //printf("@0@ i:%4d r:%4d m:%4d: beta reset to 0.0\n", it, rank, mloc);
                    } else {
                        phen.set_marker_beta(mloc, phen.sample_norm_rng(muk[i], phen.get_sigmae() / denom[i-1]));
                        //printf("@B@ i:%4d r:%4d m:%4d:  dbetat = %20.15f, muk[%4d] = %15.10f with prob=%15.10f <= acum = %15.10f, denom = %15.10f, sigmaE = %15.10f: beta = %15.10f\n", it, rank, mloc, phen.get_marker_beta(mloc) - dbeta, i, muk[i], prob, phen.get_marker_acum(mloc), denom[i-1], phen.get_sigmae(), phen.get_marker_beta(mloc));
                        //fflush(stdout);
                    }
                    phen.increment_cass(mgrp, i, 1);
                    //std::cout << "cass " << mgrp << " " << i << " = " << phen.get_cass_for_group(mgrp,i) << std::endl;
                    phen.set_comp(mloc, i);
                    break;
                } else {
                    bool zero_inc = false;
                    for (int j=i+1; j<K; j++) {
                        if (abs(logl[j] - logl[i+1]) > 700.0){
                            zero_inc = true;
                        }
                    }
                    if (!zero_inc) {
                        double esum = 0.0;
                        for (int k=0; k<logl.size(); k++){
                            esum += exp(logl[k] - logl[i+1]);
                        }
                        phen.set_marker_acum(mloc, phen.get_marker_acum(mloc) + 1.0 / esum);
                    }
                }
            }

            dbeta -= phen.get_marker_beta(mloc);
            
            if(model=="probit"){
                if(phen.get_marker_beta(mloc) != 0.0){
                    // Add current marker effect to latent variable
                    size_t methix = size_t(mloc) * size_t(N);
                    double* meth = &meth_data[methix];
                    phen.update_latent(mloc, meth);
                }
            }

            //printf("iteration %3d, rank %3d, marker %5d: dbeta = %20.15f, pheni = %d, %20.15f %20.15f\n", it, rank, mloc, dbeta, pheni, dbetas[pheni + 1], dbetas[pheni + 2])
            if (abs(dbeta) > 0.0) {
                share_mrk = true;
                dbetas[0] = dbeta;
                dbetas[1] = phen.get_marker_ave(mloc);
                dbetas[2] = phen.get_marker_sig(mloc);
                //printf("@¢@ i:%3d r:%3d m:%5d: dbeta = %20.15f, pheni = %d, %20.15f %20.15f\n", it, rank, mloc, dbeta, pheni, dbetas[pheni * 3 + 1], dbetas[pheni * 3 + 2]);
            }
                
            size_t methix = size_t(mloc) * size_t(N);

            if (share_mrk){
                phen.update_epsilon(dbetas, &meth_data[methix]);
            }

            fflush(stdout);

        } // End marker loop

        phen.reset_beta_sqn_to_zero();
        for (int i=0; i<M; i++) {
            phen.increment_beta_sqn(get_marker_group(i), phen.get_marker_beta(i) * phen.get_marker_beta(i));
        }

        // Update global parameters
        //
        for (int i=0; i<G; i++) {

            // Skip empty groups
            if (mtotgrp.at(i) == 0)
                continue;

            phen.set_m0_for_group(i, mtotgrp.at(i) - phen.get_cass_for_group(i, 0));

            // Skip groups with m0 being null or empty cass (adaV in action)
            if (phen.get_m0_for_group(i) == 0 || phen.get_cass_sum_for_group(i) == 0) {
                phen.set_sigmag_for_group(i, 0.0);
                continue;
            }

            phen.set_sigmag_for_group(i, phen.sample_inv_scaled_chisq_rng(V0G + (double) phen.get_m0_for_group(i), (phen.get_beta_sqn_for_group(i) * (double) phen.get_m0_for_group(i) + V0G * S02G) / (V0G + (double) phen.get_m0_for_group(i))));
            phen.update_pi_est_dirichlet(i);
        }

        double e_sqn = phen.epsilon_sumsqr();

        phen.set_sigmae(phen.sample_inv_scaled_chisq_rng(V0E + (double)N, (e_sqn + V0E * S02E) / (V0E + (double)N)));

        double sigmae_r0 = phen.get_sigmae();

        phen.set_sigmae(sigmae_r0);

        printf("RESULT : i:%d r:%d  sum sigmaG = %20.15f  sigmaE = %20.15f\n", it, rank, phen.get_sigmag_sum(), phen.get_sigmae());        

        double te_it = MPI_Wtime();

        //printf("INFO : It %d  Global parameters time = %7.3f\n", it, te_it - ts_debug_gp);

        printf("RESULT : It %d  total proc time = %7.3f sec", it, te_it - ts_it);

        // Write output files
        if (it % opt.get_output_thin_rate() == 0) {
            const unsigned nthinned = it / opt.get_output_thin_rate() - 1;
            write_ofile_csv(*(phen.get_outcsv_fh()), it,  phen.get_sigmag(), phen.get_sigmae(), phen.get_m0_sum(), nthinned, phen.get_pi_est());
            if(C > 0){
                write_ofile_cov(*(phen.get_outcov_fh()), it,  phen.get_deltas(), nthinned);
            }
            
            write_ofile_h1(*(phen.get_outbet_fh()), rank, Mt, it, nthinned, S, M, phen.get_betas().data(), MPI_DOUBLE);
            write_ofile_h1(*(phen.get_outcpn_fh()), rank, Mt, it, nthinned, S, M, phen.get_comp().data(),  MPI_INTEGER); 
        }

    } // End iteration loop

    phen.close_output_files();

    double te_overall = MPI_Wtime();

    printf("OVERALL : Overall proc time = %7.3f sec\n",te_overall - ts_overall);
}

double Bayes::dot_product(const int mrki, double* __restrict__ phen, const double mu, const double sigma_inv) {
    size_t methix = size_t(mrki) * size_t(N);
    double* meth = &meth_data[methix];
    double dpa = 0.0;

#ifdef _OPENMP
#pragma omp parallel for simd schedule(static) reduction(+:dpa)
#endif
    for (int i=0; i<N; i++)
        dpa += (meth[i] - mu) * phen[i];

    return sigma_inv * dpa;
}


// Setup processing: load input files and define MPI task workload
void Bayes::setup_processing() {

    mbytes = N * sizeof(double);

    load_methylation_data();

    int C = opt.get_cov_num();
    if(C > 0)
        load_covariates();

    pmgr.read_phen_files(opt, get_N(), get_M());

    check_processing_setup();

    if (rank == 0)
        printf("INFO   : output directory: %s\n", opt.get_out_dir().c_str());

    MPI_Barrier(MPI_COMM_WORLD);
    double ts = MPI_Wtime();
    pmgr.compute_markers_statistics(meth_data, get_N(), get_M());
    MPI_Barrier(MPI_COMM_WORLD);
    double te = MPI_Wtime();
    if (rank == 0)
        printf("INFO   : Time to compute the markers' statistics: %.2f seconds.\n", te - ts);

    if (opt.predict() || opt.test()) return;


        Phenotype& phen = pmgr.get_phens()[0];
        phen.set_Z(Z); // set covariates 
        phen.set_prng_m((unsigned int)(opt.get_seed() + (rank + 0)));
        if (opt.mimic_hydra()) {
            phen.set_prng_d((unsigned int)(opt.get_seed() + (rank + 0) * 1000));
        } else {
            phen.set_prng_d((unsigned int)(opt.get_seed() + (rank + 1) * 1000));
        }
    

    read_group_index_file(opt.get_group_index_file());

    for (int i=0; i<Mt; i++) {
        mtotgrp.at(get_marker_group(i)) += 1;
    }
    //for (int i=0; i<G; i++)
    //    printf("mtotgrp at %d = %d\n", i, mtotgrp.at(i));
    //pmgr.display_markers_statistics(get_M());
}


void Bayes::check_openmp() {
#ifdef _OPENMP
#pragma omp parallel
    {
        if (rank == 0) {
            int nt = omp_get_num_threads();
            if (omp_get_thread_num() == 0)
                printf("INFO   : OMP parallel regions will use %d thread(s)\n", nt);
        }
    }
#else
    printf("WARNING: no OpenMP support!\n");
#endif
}

void Bayes::read_group_index_file(const std::string& file) {

    std::ifstream infile(file.c_str());
    if (! infile)
        throw ("Error: can not open the group file [" + file + "] to read. Use the --group-index-file option!");

    if (rank == 0)
        std::cout << "INFO   : Reading groups from " + file + "." << std::endl;

    std::string label;
    int group;

    group_index.clear();

    while (infile >> label >> group) {
        //std::cout << label << " - " << group << std::endl;
        if (group > G) {
            printf("FATAL  : group index file contains a value that exceeds the number of groups given in group mixture file.\n");
            printf("       : check the consistency between your group index and mixture input files.\n");
            exit(1);
        }
        group_index.push_back(group);
    }
}


void Bayes::check_processing_setup() {
    Phenotype& phen = pmgr.get_phens()[0];
        const int Np = phen.get_nas() + phen.get_nonas();
        if (Np != N) {
            std::cout << "Fatal: N = " << N << " while phen file " << phen.get_filepath() << " has " << Np << " individuals!" << std::endl;
            exit(1);
        }
    
}

void Bayes::load_methylation_data(){

    double ts = MPI_Wtime();

    MPI_File methfh;
    const std::string methfp = opt.get_bin_file();

    check_mpi(MPI_File_open(MPI_COMM_WORLD, methfp.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &methfh),  __LINE__, __FILE__);

    const size_t size_bytes = size_t(M) * size_t(N) * sizeof(double);

    meth_data = (double*)_mm_malloc(size_bytes, 32);

    printf("INFO  : rank %d has allocated %zu bytes (%.3f GB) for raw data.\n", rank, size_bytes, double(size_bytes) / 1.0E9);

    // Offset to section of bin file to be processed by task
    MPI_Offset offset = size_t(0) + size_t(S) * size_t(N) * sizeof(double);

    // Gather the sizes to determine common number of reads
    size_t max_size_bytes = 0;
    check_mpi(MPI_Allreduce(&size_bytes, &max_size_bytes, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD), __LINE__, __FILE__);

    const int NREADS = size_t( ceil(double(max_size_bytes)/double(INT_MAX/2)) );
    size_t bytes = 0;
    mpi_file_read_at_all <double*> (size_t(M) * size_t(N), offset, methfh, MPI_DOUBLE, NREADS, meth_data, bytes);

    MPI_File_close(&methfh);
    
    double te = MPI_Wtime();

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank >= 0)
        printf("INFO   : time to load methylation data = %.2f seconds.\n", te - ts);
}

// values should be separate with space delimiter
void Bayes::load_covariates(){ 

    std::string covfp = opt.get_cov_file();
    int C = opt.get_cov_num();
    if (C==0)
        return;

    double ts = MPI_Wtime();

    std::ifstream covf(covfp);
    std::string line; 
    std::regex re("\\s+");

    int line_i = 0;

    while (std::getline(covf, line)) // read the current line
    {
        if(line_i == 0) {
            line_i++;
            continue; //skip header
        } 

        int Cobs = 0;
        std::vector<double> entries;
        std::sregex_token_iterator iter(line.begin(), line.end(), re, -1);
        std::sregex_token_iterator re_end;
        
        ++iter; // skip individual ID
        ++iter; // skip family ID
        for ( ; iter != re_end; ++iter){
            entries.push_back(std::stod(*iter));
            Cobs++;
        }

        if (Cobs != C){
            printf("FATAL   : number of covariates = %d does not match to the specified number of covariates = %d.\n", Cobs, C);
            exit(EXIT_FAILURE);
        }   
        Z.push_back(entries);
        line_i++;       
    }

    if (rank == 0)
        printf("INFO   : Number of loaded covariate lines from .cov file = %d \n", line_i - 1);

    double te = MPI_Wtime();

    if (rank == 0)
        printf("INFO   : time to load covariates = %.2f seconds.\n", te - ts);
    
    // Normalize covariates
    for(int covi = 0; covi < C; covi++){
            
        long double cavg = 0.0;
        long double csig = 0.0;

        for (int i = 0; i < N; i++) {
            cavg += Z[i][covi];    
        }
        cavg = cavg / double(N);

        for (int i = 0; i < N; i++) {
            csig += ((Z[i][covi] - cavg) * (Z[i][covi] - cavg));            
        }
        csig = sqrt(csig / double(N));

        for (int i = 0; i < N; i++) {
            if(csig < 0.00000001)
                Z[i][covi] = 0;
            else 
                Z[i][covi] = (Z[i][covi] - cavg) / csig;
        }
    }
}

void Bayes::print_cva() {
    printf("INFO   : mixtures for all groups:\n");
    for (int i=0; i<G; i++) {
        printf("         grp %2d: ", i);
        for (int j=0; j<K; j++) {
            printf("%7.5f ", cva[i][j]);
        }
        printf("\n");
    }
}

void Bayes::print_cvai() {
    printf("INFO   : inverse mixtures for all groups:\n");
    for (int i=0; i<G; i++) {
        printf("         grp %2d: ", i);
        for (int j=0; j<K; j++) {
            printf("%10.3f ", cvai[i][j]);
        }
        printf("\n");
    }
}
