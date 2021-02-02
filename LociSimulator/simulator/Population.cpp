#include "Population.h"

Population::Population(vector<size_t> num_of_different_alleles, size_t pop_size, size_t selective_allele, size_t selective_locus, double selection_intensity) :
	_num_of_loci(num_of_different_alleles.size()), _selective_allele(selective_allele), _selective_locus(selective_locus), _selection_intensity(selection_intensity){
	vector<int> loci(_num_of_loci,0); // starting an empty bacteria
	_selective_pop_size = 0;
	_neutral_pop_size = 0;
	for(size_t i=0; i < pop_size; i++){
		_neutral_bacteria.push_back(Bacterium(loci));
		_selective_bacteria.push_back(Bacterium(loci));
		if(is_selective_allele(selective_locus, loci[selective_locus])){ // each population counter determines whether the current bacteria is "relevant" or not
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

int Population::get_allele_from_bacterium(size_t bacterium_index, size_t locus_index) {
	return get_bacterium(bacterium_index).get_allele(locus_index);
	/*if(bacterium_index < get_neutral_pop_size()){
		return _neutral_bacteria[bacterium_index].get_allele(locus_index);
	}
	else{
		return _selective_bacteria[bacterium_index - get_neutral_pop_size()].get_allele(locus_index);
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

void Population::set_allele_to_bacterium(size_t bacterium_index, size_t locus_index, size_t old_allele, size_t new_allele){
	if(!is_selective_allele(locus_index, old_allele) && is_selective_allele(locus_index, new_allele)){ //adaptation!
		_neutral_bacteria[bacterium_index].set_allele(locus_index,new_allele);
		_selective_bacteria.push_back(_neutral_bacteria[bacterium_index]); //set the bacterium in the right group
		_neutral_bacteria.erase(_neutral_bacteria.begin() + bacterium_index); //remove the bacterium from the old group
	}
	else if(is_selective_allele(locus_index, old_allele) && !is_selective_allele(locus_index, new_allele)){ //back to neutrality
		_selective_bacteria[bacterium_index].set_allele(locus_index, new_allele);
		_neutral_bacteria.push_back(_selective_bacteria[bacterium_index]); //set the bacterium in the right group
		_selective_bacteria.erase(_selective_bacteria.begin() + bacterium_index); //remove the bacterium from the old group
	}
	else if(!is_selective_allele(locus_index, old_allele) && !is_selective_allele(locus_index, new_allele)){ //stay neutral
		_neutral_bacteria[bacterium_index].set_allele(locus_index, new_allele);
	}
	else{ //i.e. (is_selective_allele(locus_index, old_allele) && is_selective_allele(locus_index, new_allele)) //stay selective
		_selective_bacteria[bacterium_index].set_allele(locus_index, new_allele);
	}
}

bool Population::is_selective_allele(size_t locus_index, size_t allele){
	return locus_index == _selective_locus && allele == _selective_allele;
}

bool Population::is_selective_bacterium(const Bacterium bacterium) const{
	return bacterium.get_allele(_selective_locus) == static_cast<int>(_selective_allele);
}

//draw selectively bacterium index to be set as a parent
size_t Population::draw_bacterium_index_selectively() {
	
	if (get_selective_pop_size() == 0){
		return rand() % get_neutral_pop_size();
	}

	if (get_neutral_pop_size() == 0){
		return rand() % get_selective_pop_size();
	}

	double random_event = rand()*1.0/RAND_MAX;
	if (random_event < get_p_neutral()){
		return rand() % get_neutral_pop_size();
	}
	return (rand() % get_selective_pop_size()) + get_neutral_pop_size(); //always return an index with respect to the total size
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
int Population::write_pop_sequences(const vector<LocusSequenceSimulator> & locus_sequences_simulators, const string out_path, const size_t pop_index) {
	ostringstream out_stream;
	size_t allele;
	cout << endl << "Writing sequences of pop " << pop_index << " to:" << endl << out_path <<endl;
	for(size_t bacterium_index=0; bacterium_index < get_total_size(); bacterium_index++){
		/*if(i%10==0){
			cout << "Writing bac "<<i<<"-"<<min(i+9,pop.get_total_size()-1)<<" sequences..." <<endl;
		}*/
		if (DEBUG_MODE) {
			cout << "Writing bac "<<bacterium_index<<" sequences..." <<endl;
		}

		if (get_bacterium(bacterium_index).get_allele(_selective_locus) == _selective_allele) {
			// note that the sample allele is the selective allele 
			out_stream << ">" << pop_index << "_" << bacterium_index << "_selective";
		}
		else {
			out_stream << ">" << pop_index << "_" << bacterium_index;
		}
		out_stream << endl;

		for(size_t j=0; j < get_num_of_loci(); j++){
			allele = get_bacterium(bacterium_index).get_allele(j);
			LocusSequenceSimulator locus_sequences_simulator = locus_sequences_simulators[j];
			out_stream <<  locus_sequences_simulator._sequences[allele];
		}
		out_stream << endl;
	}
	write_ostringstream_to_file(out_stream, out_path);

	return 1;
}

void Population::check_pop_validity(void){
	for (size_t i = 0; i < get_neutral_pop_size(); i++){
		if(_neutral_bacteria[i].get_allele(_selective_locus) == static_cast<int>(_selective_allele)){
			cout << "Invariant error! bacterium " << i << " in _neutral_bacteria is SELECTIVE bacterium!" <<endl;
		}
	}
	for (size_t i = 0; i < get_selective_pop_size(); i++){
		if(_selective_bacteria[i].get_allele(_selective_locus) != static_cast<int>(_selective_allele)){
			cout << "Invariant error! bacterium " << i << " in _selective_bacteria is NEUTRAL bacterium!" <<endl;
		}
	}
}
