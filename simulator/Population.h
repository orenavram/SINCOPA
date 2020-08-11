#ifndef CLASS_POPULATION
#define CLASS_POPULATION

#include "Bacterium.h"
#include "LocusSequencesSimulator.h"
#include <random>
using namespace std;

class Population{
public:
	vector<Bacterium> _neutral_bacteria;
	vector<Bacterium> _selective_bacteria;
	size_t _selective_pop_size;
	size_t _neutral_pop_size;
	size_t _selective_allele;
	double _selection_intensity;

	
	


	//constructors
	Population(size_t num_of_different_alleles, size_t pop_size, size_t selective_allele, double selection_intensity);
	//Population(vector<size_t> num_of_different_alleles); // empty population

	size_t get_total_size() const { return _neutral_pop_size + _selective_pop_size; }

	size_t get_neutral_pop_size() const {
		return _neutral_pop_size; 
	}

	size_t get_selective_pop_size() const  { return _selective_pop_size; }
	
	//return a pointer to a bacterium with index bacterium_index
	Bacterium & get_bacterium(size_t bacterium_index) ;
	Bacterium & operator[](const size_t bacterium_index) {
		if (bacterium_index < get_neutral_pop_size()) {
			return _neutral_bacteria[bacterium_index];
		}
		else {
			return _selective_bacteria[bacterium_index - get_neutral_pop_size()];
		}
	}
	

	int get_allele_from_bacterium(size_t bacterium_index) ;
	
	double get_p_neutral() const { return get_neutral_pop_size() * 1.0 / get_population_fitness(); }

	double get_population_fitness() const {	return get_neutral_pop_size() + get_selective_pop_size() * _selection_intensity; }

	bool is_selective_bacterium(const Bacterium & bacterium) const;

	bool is_selective_allele(size_t allele);

	void check_pop_validity(void); // verify that _neutral_bacteria hold only neutral bacteria and _selective_bacteria hold only selective bacteria



	void set_allele_to_bacterium(size_t bacterium_index, size_t old_allele, size_t new_allele);
	
	void add_bacterium(Bacterium & bacterium_to_add) ;



	size_t draw_bacterium_index_selectively() ;
	
	size_t draw_other_bacterium_index_uniformly(size_t bacterium_index) const ;

	int write_pop_loci(const string out_path) ;

	int write_pop_sequences(const LocusSequencesSimulator & locus_sequences_simulators, const string out_path, const size_t pop_index) ;
	
	//void set_bacteria_loci(size_t bacterium_index, vector<int> loci);

	void print_neutral_pop();
	
	void print_selective_pop();
	
	void print();
};

#endif