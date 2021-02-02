
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <cstring>      
#include <vector>
#include <ctime>
#include <cmath>
#include <map>
#include <algorithm>
#include <fstream>
#include <sys/stat.h>

using namespace std;

#define DEBUG_MODE 0
#define PRINT_MUTATION_HISTORY 0
#define WRITE_MUTATION_HISTORY 1
#define TRACK_SELECTIVE_ALLELE 0
#define NOTIFY_BETWEEN_POPULATION_RECOMBINATION 1
#define WRITE_POPS_ALLELE_NUMBERS 0 
#define WRITE_POPS_SEQUENCES 1 
#define READ_MUTATIONS 0

double draw_from_zero_to_one_in_high_resolution_steps(int fraction);

int write_ostringstream_to_file(ostringstream & out_stream, const string out_path);

int read_control_file(const string in_path, int& pop_size, int& num_of_pops, long& num_of_generations, long& super_generation, size_t& num_of_loci, size_t& locus_len, size_t& num_of_non_recombining_loci, double& p_a, double& p_c, double& p_g, double& p_transition, double& mutation_rate_per_bacterium_per_generation, double& recombination_rate_per_bacterium_per_generation, double& between_pops_recombination_rate_per_bacterium_per_generation, double& selection_intensity, double& selective_allele_timing_proportion_with_respect_to_num_of_generations, string& output_dir);

string remove_trailing_whitespaces(string str);