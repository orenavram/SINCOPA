#ifndef CLASS_BACTERIUM
#define CLASS_BACTERIUM

#include "Auxilaries.h"

class Bacterium{
public:
	int _locus; // _locus is the allele of the locus.

	//constructor
	explicit Bacterium(int inputLocus) : _locus(inputLocus) {}
	explicit Bacterium(const Bacterium & inputBac) : _locus(inputBac._locus) {}
	
	int get_allele(void) const ;

	void set_allele(const size_t allele); 



	////return how many loci differ between "this" and "other"
	//int diff(const Bacterium & other) const;

	
	//return all loci concatenated
	void print(size_t ) const;

	string loci_to_string(const string delimiter = "\t") const;
};

#endif
