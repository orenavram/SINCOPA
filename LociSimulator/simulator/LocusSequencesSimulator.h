#ifndef CLASS_LOCUSSEQUENCESSIMULATOR
#define CLASS_LOCUSSEQUENCESSIMULATOR

#include "LocusTreeStructure.h"

class LocusSequenceSimulator{

public:
	LocusTreeStructure _locus_tree_structure; //see LocusTreeStructure.h for documentation
	const size_t _locus_len;
	map<size_t, string> _sequences; //mapping from allele index to nucleotide sequence

	//constructor
	LocusSequenceSimulator(LocusTreeStructure & locus_tree_structure, const size_t locus_len);

	void simulate(const double p_a, const double p_c, const double p_g, const double p_transition);

	void mutate_sequence(const size_t allele, const double p_transition);

	void mutate_nucleotide(const size_t random_nucleotide_index, const int allele, const double p_transition);

	string draw_root_sequence(const double p_a, const double p_c, const double p_g) const;

	char apply_transition(const char nuc) const;

	char apply_transversion(const char nuc) const;

	void check_data_validity3(void);

};

#endif