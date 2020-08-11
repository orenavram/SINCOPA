#include "Population.h"

Population::Population(	size_t num_of_different_alleles,
						size_t pop_size, 
						size_t selective_allele,
						double selection_intensity) :
	_selective_allele(selective_allele), _selection_intensity(selection_intensity){
	int locus = 0; // starting an empty bacteria
	_selective_pop_size = 0;
	_neutral_pop_size = 0;
	for(size_t i=0; i < pop_size; i++){
		_neutral_bacteria.push_back(Bacterium(locus));
		_selective_bacteria.push_back(Bacterium(locus));
		if(is_selective_allele(locus)){ // each pop counter determine whether the current bacteria is "relevant" or not
			_selective_pop_size++;
		}
		else{
			_neutral_pop_size++;
		}
	}
}

Bacterium & Population::get_bacterium(size_t bacterium_index) {
	if(bacterium_index < get_neutral_pop_size()){
		return _neutral_bacteria[bacterium_index];
	}
	else{
		return _selective_bacteria[bacterium_index-get_neutral_pop_size()];
	}
}

int Population::get_allele_from_bacterium(size_t bacterium_index) {
	return get_bacterium(bacterium_index).get_allele();
	/*if(bacterium_index < get_neutral_pop_size()){
		return _neutral_bacteria[bacterium_index].get_allele();
	}
	else{
		return _selective_bacteria[bacterium_index - get_neutral_pop_size()].get_allele();
	}*/
}

void Population::print_neutral_pop(){
	cout << "Neutral bacteria (" << _neutral_pop_size << " in total):" << endl;
	for (size_t i=0; i<get_neutral_pop_size(); ++i){
		_neutral_bacteria[i].print(i);
	}
	cout << endl;
}

void Population::print_selective_pop(){
	cout << "Selective bacteria (" << _selective_pop_size << " in total):" << endl;
	for (size_t i=0; i<get_selective_pop_size(); ++i){
		_selective_bacteria[i].print(i);
	}
	cout << endl;
}

void Population::print() {
	print_neutral_pop();
	print_selective_pop();
}

void Population::set_allele_to_bacterium(size_t bacterium_index, size_t old_allele, size_t new_allele){
	if(!is_selective_allele(old_allele) && is_selective_allele(new_allele)){ //adaptation!
		_neutral_bacteria[bacterium_index].set_allele(new_allele);
		_selective_bacteria.push_back(_neutral_bacteria[bacterium_index]); //set the bacterium in the right group
		_neutral_bacteria.erase(_neutral_bacteria.begin() + bacterium_index); //remove the bacterium from the old group
	}
	else if(is_selective_allele(old_allele) && !is_selective_allele(new_allele)){ //back to neutrality
		_selective_bacteria[bacterium_index].set_allele(new_allele);
		_neutral_bacteria.push_back(_selective_bacteria[bacterium_index]); //set the bacterium in the right group
		_selective_bacteria.erase(_selective_bacteria.begin() + bacterium_index); //remove the bacterium from the old group
	}
	else if(!is_selective_allele(old_allele) && !is_selective_allele(new_allele)){ //stay neutral
		_neutral_bacteria[bacterium_index].set_allele(new_allele);
	}
	else{ //i.e. (is_selective_allele(old_allele) && is_selective_allele(new_allele)) //stay selective
		_selective_bacteria[bacterium_index].set_allele(new_allele);
	}
}

bool Population::is_selective_allele(size_t allele){
	return allele == _selective_allele;
}

bool Population::is_selective_bacterium(const Bacterium & bacterium) const{
	return bacterium.get_allele() == static_cast<int>(_selective_allele);
}

//draw selectively bacterium index to be set as a parent
size_t Population::draw_bacterium_index_selectively() {
	uniform_real_distribution<double> _distribution1(0, get_population_fitness());
	double randomN = _distribution1(_generator1);
	//cout << "xxx " << randomN << " xxx" <<endl;
	if (randomN < get_neutral_pop_size()) {
		// neutral bacterium was drawn
		return static_cast<size_t>(randomN);
	}
	// else: selective bacterium was drawn
	double x = ((randomN - get_neutral_pop_size()) / _selection_intensity) + get_neutral_pop_size();
	return static_cast<size_t>(x);

}

//draw uniformly bacterium index (other than bac_ind) in order to make recombination with it
size_t Population::draw_other_bacterium_index_uniformly(size_t bacterium_index) const {
	size_t new_ind = rand() % (get_total_size() - 1);
	if(new_ind >= bacterium_index){
		new_ind++;
	}
	return new_ind;
}

void Population::add_bacterium(Bacterium & bacterium_to_add) {
	if(is_selective_bacterium(bacterium_to_add)){
		_selective_bacteria[_selective_pop_size] = bacterium_to_add;
		_selective_pop_size++;
	}
	else{
		_neutral_bacteria[_neutral_pop_size] = bacterium_to_add;
		_neutral_pop_size++;
	}
}

//write final_pop to an output file
int Population::write_pop_loci(const string out_path) {
	ostringstream out_stream;
	cout << endl << "Writing final population to:"<< endl<< out_path << endl;
	for(size_t i=0; i < get_total_size(); i++){
		out_stream << ">" << i << endl;
		out_stream << get_bacterium(i).loci_to_string() << endl;
	}
	
	write_ostringstream_to_file(out_stream, out_path);
	return 1;
}

//write sequences of all loci concatenated together to an output file
int Population::write_pop_sequences(const LocusSequencesSimulator & locus_sequences_simulators, const string out_path, const size_t pop_index) {
	ostringstream out_stream;
	size_t allele;
	cout << endl << "Writing sequences of pop " << pop_index << " to:" << endl << out_path <<endl;
	for(size_t i=0; i < get_total_size(); i++){
		/*if(i%10==0){
			cout << "Writing bac "<<i<<"-"<<min(i+9,pop.get_total_size()-1)<<" sequences..." <<endl;
		}*/
		if (PRINT_MODE && DEBUG_MODE) {
			cout << "Writing bac "<<i<<" sequences..." <<endl;
		}
		out_stream << ">" << pop_index << "_" << i << endl;
		allele = get_bacterium(i).get_allele();
		LocusSequencesSimulator locus_sequences_simulator = locus_sequences_simulators;
		out_stream <<  locus_sequences_simulator._sequences[allele];
		out_stream  << endl;
	/*
	for(size_t i=0; i < pop.get_total_size(); i++){
		string out_string("");
		stringstream ss;
		ss.str("");
		ss << i;
		cout << "Writing bac "<< i <<" sequences..." <<endl;
		out_string += ">" + ss.str() + "\n";
		allele = pop.get_bacterium(i)._loci[j];
		LocusSequencesSimulator locus_sequences_simulator = locus_sequences_simulators[j];
		out_string +=  locus_sequences_simulator._sequences[allele];
		out_string +=  "\n";
	*/
	}
	write_ostringstream_to_file(out_stream, out_path);

	return 1;
}

void Population::check_pop_validity(void){
	for (size_t i = 0; i < get_neutral_pop_size(); i++){
		if(_neutral_bacteria[i].get_allele() == static_cast<int>(_selective_allele)){
			cout << "Invariant error! bacterium " << i << " in _neutral_bacteria is SELECTIVE bacterium!" <<endl;
		}
	}
	for (size_t i = 0; i < get_selective_pop_size(); i++){
		if(_selective_bacteria[i].get_allele() != static_cast<int>(_selective_allele)){
			cout << "Invariant error! bacterium " << i << " in _selective_bacteria is NEUTRAL bacterium!" <<endl;
		}
	}
}

//void Population::set_bacteria_loci(size_t bacterium_index, vector<int> loci) {
//	if(bacterium_index < get_neutral_pop_size()){
//		_neutral_bacteria[bacterium_index]._loci = loci;
//	}
//	else{
//		_selective_bacteria[bacterium_index - get_neutral_pop_size()]._loci = loci;
//	}
//}