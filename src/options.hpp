#pragma once
#include <string>
#include <vector>

class Options {

public:
    Options() = default;
    Options(int argc, char** argv) {
        read_command_line_options(argc, argv);
        check_options();
        if (!predict() && !test())
            read_group_mixture_file();
    }
    void read_command_line_options(int argc, char** argv);
    std::string get_bin_file() const { return bin_file; }
    std::string get_dim_file() const { return dim_file; }
    std::string get_bim_file() const { return bim_file; }
    std::string get_ref_bim_file() const { return ref_bim_file; }
    std::string get_group_index_file() const { return group_index_file; }
    std::string get_group_mixture_file() const { return group_mixture_file; }
    std::string get_out_dir() const { return out_dir; }
    unsigned get_output_thin_rate() const { return thin; }
    const std::vector<std::string>& get_phen_files() const { return phen_files; }
    void list_phen_files() const;
    int  count_phen_files() const { return phen_files.size(); }
    int  get_verbosity() const { return verbosity; }
    bool verbosity_level(const int level) const { return level > get_verbosity() ? false : true; }
    bool shuffle_markers() const { return shuffle; }
    bool mimic_hydra() const { return mimic_hydra_; }
    bool predict() const { return predict_; }
    bool test() const { return test_; }
    unsigned int get_seed() const { return seed; }
    unsigned int get_iterations() const { return iterations; }
    unsigned int get_burn_in() const { return burn_in; }
    unsigned int get_truncm() const { return truncm; }
    const std::vector<double>& get_s() const { return S; }
    int get_ngroups()   const { return ngroups; }
    int get_nmixtures() const { return nmixtures; }
    const std::vector<std::vector<double>>& get_cva()  const { return cva;  }
    const std::vector<std::vector<double>>& get_cvai() const { return cvai; }
    std::string get_cov_file() const { return cov_file; };
    int get_cov_num() const { return cov_num; };
    std::string get_infname_base() const { return in_fname_base; };
    std::string get_model() const { return model; };

private:
    std::string bin_file = "";
    std::string dim_file = "";
    std::string bim_file = "";
    std::string ref_bim_file = "";
    std::string group_index_file   = "";
    std::string group_mixture_file = "";
    std::string out_dir = "";
    int verbosity = 0;
    bool shuffle = true;
    bool mimic_hydra_ = false;
    bool predict_ = false;
    bool test_ = false;
    unsigned int seed = 0;
    unsigned int iterations = 1;
    unsigned int truncm = 0;
    unsigned int thin = 1;
    unsigned int burn_in = 0;
    std::vector<std::string> phen_files;
    std::vector<double> S;
    int ngroups = 0;
    void _set_ngroups(const int i) { ngroups = i; }
    int nmixtures = 0;
    void _set_nmixtures(const int i) { nmixtures = i; }
    void check_options();
    void fail_if_last(char** argv, const int i);
    void read_group_mixture_file();
    std::vector<std::vector<double>> cva, cvai;
    std::string cov_file = "";
    int cov_num = 0;
    std::string in_fname_base = "";
    std::string model = "linear";
};
