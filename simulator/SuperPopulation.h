#ifndef CLASS_SUPERPOPULATION
#define CLASS_SUPERPOPULATION

#include "Population.h"

class SuperPopulation{
public:
	vector<Population> _pops;
	size_t _num_of_different_alleles; // how many different alleles at the locus
	size_t _selective_allele;
	double _selection_intensity;
	size_t _pop_size;

	explicit SuperPopulation(	size_t num_of_pops,
								size_t num_of_different_alleles,
								size_t pop_size,
								size_t selective_allele,
								double selection_intensity);

	Population & get_pop(const size_t pop_index) { return _pops[pop_index];}
	Population &operator[] (const size_t pop_index) { return _pops[pop_index]; }

	size_t get_num_of_pops() const { return _pops.size();}
	size_t get_pop_size() const { return _pop_size; }

	size_t get_num_of_different_alleles_at_locus() const { return _num_of_different_alleles;}
	
	vector<int> get_locus_superpop_wide(size_t total_num_of_bacteria, size_t pop_size);

	void check_superpop_validity(void); // verify that each allele in the super_populaiton is not "newer" (i.e. larger) the the "newest" (i.e. largest) allele at the locus it resides in


	
	void add_new_allele_to_locus() { _num_of_different_alleles++; }

	size_t draw_pop_index_uniformly() const { return rand() % get_num_of_pops();}
	
	void set_pop(Population pop, const size_t pop_index)  { _pops[pop_index] = pop; }



	int write_pops_loci(const string out_prefix, const string out_suffix) ;

	int write_pops_sequences(const LocusSequencesSimulator & locus_sequences_simulator, const string out_prefix, const string out_suffix) ;

	void sample_and_write_sequences(const LocusSequencesSimulator & locus_sequences_simulator, const string out_path, const string out_suffix) ;

	void print_locus_superpop_wide_unique(size_t total_num_of_bacteria, size_t pop_size);
	void print();
	void print_num_of_mutations(void);
};

#endif