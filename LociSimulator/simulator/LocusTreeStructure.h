#include "Auxilaries.h"

#ifndef CLASS_LOCUSTREESTRUCTURE
#define CLASS_LOCUSTREESTRUCTURE

class LocusTreeStructure{

public:
	map<int, int> _allelic_tree; //mapping from an allele to its least common ancestor (with respect to another allele)
	map<int, int> _num_of_mutations_tree; //mapping from an allele to the number of mutations occured from its least common ancestor

	//constructor
	LocusTreeStructure(const vector< vector<int> > & mutations_history, vector<int> current_alleles, size_t locus_index, size_t num_of_bacteria, size_t pop_size);

	vector<int> coalesce(vector<int> current_alleles, vector<int> & last_coalecsent, vector<int> & mutation_counters, int lca_allele, int num_of_occurences);

	bool there_is_more_than_one_allele(vector<int> & last_coalecsent);

	int get_oldest_allele(void);
};

#endif