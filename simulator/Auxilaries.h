#ifndef AUXILARIES___
#define AUXILARIES___

#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <cstring>      
#include <vector>
#include <ctime>
#include <map>
#include <algorithm>
#include <random>

using namespace std;

#define PRINT_MODE 0
#define DEBUG_MODE 0
#define PRINT_MUTATION_HISTORY 0
#define WRITE_MUTATION_HISTORY 1
#define READ_MUTATIONS 0

int write_ostringstream_to_file(ostringstream & out_stream, const string out_path);

int read_control_file(const string in_path, int& pop_size, int& num_of_pops, long& max_num_of_steps, long& super_step, size_t& locus_len, double& p_a, double& p_c, double& p_g, double& p_transition, double& mutation_rate_per_bacterium_per_step, double& recombination_rate_per_bacterium_per_step, double& between_pops_recombination_rate_per_bacterium_per_step, double& selection_intensity, int& selective_allele, string& sampled_sequences_out_prefix, string& done_path);

string remove_trailing_whitespaces(string str);

static default_random_engine _generator1;

#endif