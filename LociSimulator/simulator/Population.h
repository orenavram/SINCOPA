#ifndef CLASS_POPULATION
#define CLASS_POPULATION

#include "Bacterium.h"
#include "LocusSequencesSimulator.h"

class Population{
public:
	vector<Bacterium> _neutral_bacteria;
	vector<Bacterium> _selective_bacteria;
	size_t _selective_pop_size;
	size_t _neutral_pop_size;
	size_t _num_of_loci;
	size_t _selective_allele;
	size_t _selective_locus;
	double _selection_intensity;

	//constructors
	Population(vector<size_t> num_of_different_alleles, size_t pop_size, size_t selective_allele, size_t selective_locus, double selection_intensity);

	size_t get_total_size() const { return _neutral_pop_size + _selective_pop_size; }

	size_t get_neutral_pop_size() const { return _neutral_pop_size; }

	size_t get_selective_pop_size() const  { return _selective_pop_size; }
	
	size_t get_num_of_loci() const { return _num_of_loci; }

	//return a pointer to a bacterium with index bacterium_index
	Bacterium & get_bacterium(size_t bacterium_index) ;

	int get_allele_from_bacterium(size_t bacterium_index, size_t locus_index) ;
	
	double get_p_neutral() const { return get_neutral_pop_size() * 1.0 / get_population_fitness(); }

	double get_population_fitness() const {	return get_neutral_pop_size() + get_selective_pop_size() * _selection_intensity; }

	bool is_selective_bacterium(const Bacterium bacterium) const;

	bool is_selective_allele(size_t locus_index, size_t allele);

	void check_pop_validity(void); // verify that _neutral_bacteria hold only neutral bacteria and _selective_bacteria hold only selective bacteria



	void set_allele_to_bacterium(size_t bacterium_index, size_t locus_index, size_t old_allele, size_t new_allele);
	
	void add_bacterium(Bacterium & bacterium_to_add) ;



	size_t draw_bacterium_index_selectively() ;
	
	size_t draw_other_bacterium_index_uniformly(size_t bacterium_index) const ;

	int write_pop_loci(const string out_path) ;

	int write_pop_sequences(const vector<LocusSequenceSimulator> & locus_sequences_simulators, const string out_path, const size_t pop_index) ;
	
	void print_neutral_pop();
	
	void print_selective_pop();
	
	void print();
};

#endif