#ifndef CLASS_BACTERIUM
#define CLASS_BACTERIUM

#include "Auxilaries.h"

class Bacterium{
public:
	vector<int> _loci; // a vector containing the loci of the bacterium (_loci[i] is the allele of the ith locus)

	//constructor
	Bacterium(vector<int> inputLociVector);
	Bacterium(const Bacterium & inputBacterium);
	

	//return the "genome" length
	size_t get_num_of_loci() const { return _loci.size(); }

	int get_allele(const size_t locus_index) const { return _loci[locus_index]; };

	void set_allele(const size_t locus_index, const size_t allele) { _loci[locus_index] = allele; };

	//return how many loci differ between "this" and "other"
	int diff(const Bacterium & other) const;

	//return all loci concatenated
	void print(const size_t bacterium_index) const { cout << bacterium_index << ":\t" << loci_to_string() << endl; };

	string loci_to_string(const string delimiter = "\t") const;
};

#endif
