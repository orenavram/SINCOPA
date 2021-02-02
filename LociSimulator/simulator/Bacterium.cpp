#include <sstream>
#include "Bacterium.h"

Bacterium::Bacterium(vector<int> loci) : _loci(loci) {
}

Bacterium::Bacterium(const Bacterium & input_bac) : _loci(input_bac._loci){
}


int Bacterium::diff(const Bacterium & other) const {
	int count = 0;
	for(size_t i=0; i < get_num_of_loci(); i++){
		if(get_allele(i) != other.get_allele(i)){
			count++;
		}
	}
	return count;
}

string Bacterium::loci_to_string(const string delimiter) const{
	ostringstream result;
	for(size_t i=0; i < get_num_of_loci(); i++){
		result << get_allele(i) << delimiter;
	}
	return result.str();
}