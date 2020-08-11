#include <sstream>
#include "Bacterium.h"





void Bacterium::print(const size_t bacterium_index) const {
	cout << bacterium_index << ":\t" << loci_to_string() << endl;
}

int Bacterium::get_allele(void) const {
	return _locus;
}

void Bacterium::set_allele(const size_t allele) {
	_locus = allele;
}

string Bacterium::loci_to_string(const string delimiter) const{
	ostringstream result;
	result << get_allele() << delimiter;
	return result.str();
}