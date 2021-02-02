#ifndef CLASS_SUPERPOPULATION
#define CLASS_SUPERPOPULATION

#include "Population.h"

class SuperPopulation{
public:
	vector<Population> _pops;
	vector<size_t> _num_of_different_alleles_in_locus; // how many different alleles at each locus (its size is num_of_loci)
	size_t _selective_allele;
	size_t _selective_locus;
	double _selection_intensity;

	SuperPopulation(size_t num_of_pops, vector<size_t> num_of_different_alleles, size_t pop_size, size_t selective_allele, size_t selective_locus, double selection_intensity);

	Population & get_pop(const size_t pop_index) { return _pops[pop_index];}

	size_t get_num_of_pops() const { return _pops.size();}
	
	size_t get_num_of_loci() const { return _num_of_different_alleles_in_locus.size();}
	
	size_t get_num_of_different_alleles_at_locus(const size_t locus_index) const { return _num_of_different_alleles_in_locus[locus_index];}
	
	vector<int> get_locus_i_superpop_wide(size_t locus_index, size_t num_of_bacteria, size_t pop_size);

	void check_superpop_validity(void); // verify that each allele in the super_populaiton is not "newer" (i.e. larger) the the "newest" (i.e. largest) allele at the locus it resides in


	
	void add_new_allele_to_locus(int locus_index) { _num_of_different_alleles_in_locus[locus_index]++; }

	size_t draw_pop_index_uniformly() const { return rand() % get_num_of_pops();}
	
	void set_pop(Population pop, const size_t pop_index)  { _pops[pop_index] = pop; }



	int write_pops_loci(string output_dir) ;

	int write_pops_sequences(const vector<LocusSequenceSimulator> & locus_sequences_simulators, string output_dir) ;

	int	sample_and_write_sequences(const vector<LocusSequenceSimulator> & locus_sequences_simulators, const string output_dir, size_t num_of_non_recombining_loci) ;

	void print_locus_i_superpop_wide_unique(size_t locus_index, size_t num_of_bacteria, size_t pop_size);
	
	void print_num_of_mutations(void);
};

#endif