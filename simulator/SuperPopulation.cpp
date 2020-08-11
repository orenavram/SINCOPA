#include "SuperPopulation.h"



SuperPopulation::SuperPopulation(size_t num_of_pops,
								size_t num_of_different_alleles,
								size_t pop_size,
								size_t selective_allele,
								double selection_intensity) :
	_num_of_different_alleles(num_of_different_alleles), _pop_size(pop_size),_selective_allele(selective_allele), _selection_intensity(selection_intensity){
	for(size_t i=0; i < num_of_pops; i++){
		_pops.push_back(Population(num_of_different_alleles, pop_size, selective_allele, selection_intensity));
	}
}

void SuperPopulation::print() {
	cout << "The super population has " << get_num_of_pops() << " populations" << endl;
	cout << "Each population has size: " << get_pop_size() << endl;
	for (size_t i = 0; i < get_num_of_pops(); ++i) {
		_pops[i].print();
	}
	cout << "##################################################" << endl;
}

int SuperPopulation::write_pops_loci(const string out_prefix, const string out_suffix) {
	ostringstream out_path;
	for (size_t i=0; i<get_num_of_pops(); i++){
		out_path << out_prefix << "simulatedLoci.pop." << i << out_suffix;
		get_pop(i).write_pop_loci(out_path.str());
		out_path.str(""); //clean stream
	}
	return 1;
}

int SuperPopulation::write_pops_sequences(const LocusSequencesSimulator & locus_sequences_simulator, const string out_prefix, const string out_suffix) {
	ostringstream out_path;
	for (size_t i=0; i<get_num_of_pops(); i++){
		out_path << out_prefix << "pop." << i << out_suffix;
		get_pop(i).write_pop_sequences(locus_sequences_simulator, out_path.str(), i);
		out_path.str(""); //clean stream
	}
	return 1;
}

//this function samples a bacterium from each population at the and of the evolution simulation and write the sampled bacteria sequences to all_loci_outpath.
//in addition it write to non_recombining_loci_outpath all the non recombining loci of the sampled bacteria (i.e. a SUBsequence of the full "genome"
void SuperPopulation::sample_and_write_sequences(const LocusSequencesSimulator & locus_sequences_simulator,
												const string out_prefix,
												const string out_suffix) {
	ostringstream all_out_stream;
	size_t allele;
	if (get_num_of_pops() == 1){
		cout << endl << "Num of populations is 1.  Writing the whole population instead of sampling." <<endl;
		
		Population pop = get_pop(0);
		
		for (size_t i=0; i<pop.get_total_size(); i++){
		
			if(i%25 == 0){
				cout << endl << "Writing sequence " << i << "." << endl;
			}
			if(PRINT_MODE){cout << "The sequence is:" <<endl;}
		
			all_out_stream << ">0_" << i << endl;
		
			allele = pop.get_bacterium(i).get_allele();
			
			if (locus_sequences_simulator._sequences.at(allele).size() != locus_sequences_simulator._locus_len && locus_sequences_simulator._sequences.at(allele).size() != 0){
				cout << endl << "error in locus !!! (single population run)"<< endl;
				cout << "locus_sequences_simulator._sequences[" << allele << "] is " << locus_sequences_simulator._sequences.at(allele) << endl;
				cout << "locus_sequences_simulator._sequences[" << allele << "].size() is " << locus_sequences_simulator._sequences.at(allele).length() << endl;
				cout << "locus_sequences_simulator._locus_len is " << locus_sequences_simulator._locus_len << endl;
			}
			if(PRINT_MODE){cout << locus_sequences_simulator._sequences.at(allele) << endl;}
			all_out_stream <<  locus_sequences_simulator._sequences.at(allele);
			all_out_stream  << endl;
		}
	}

	else{

		uniform_int_distribution<size_t> _distribution1(0, _pops[0].get_total_size()); // uniform distribution from zero to number of bacteria in population

		for(size_t i=0; i < get_num_of_pops(); i++){
			Population pop = get_pop(i);

			size_t random_bacterium_index = _distribution1(_generator1);
			//size_t random_bacterium_index = rand() % pop.get_total_size();
		
			cout << endl << "Sequence " << random_bacterium_index << " was drawn from pop " << i <<endl;
			if(PRINT_MODE){cout << "The sequence is:" <<endl;}
		
			all_out_stream << ">" << i << "_" << random_bacterium_index << endl;

			allele = pop.get_bacterium(random_bacterium_index).get_allele();
			
			if (locus_sequences_simulator._sequences.at(allele).size() != locus_sequences_simulator._locus_len && locus_sequences_simulator._sequences.at(allele).size() != 0){
				cout << endl << "error in locus !!!"<< endl;
				cout << "locus_sequences_simulator._sequences[" << allele << "] is " << locus_sequences_simulator._sequences.at(allele) << endl;
				cout << "locus_sequences_simulator._sequences[" << allele << "].size() is " << locus_sequences_simulator._sequences.at(allele).size() << endl;
				cout << "locus_sequences_simulator._locus_len is " << locus_sequences_simulator._locus_len << endl;
			}
			if(PRINT_MODE){cout << locus_sequences_simulator._sequences.at(allele) << endl;}
			all_out_stream <<  locus_sequences_simulator._sequences.at(allele);
			all_out_stream  << endl;
		}
	}

	cout << endl << "Writing all loci sequences to:" << endl << out_prefix << "all_loci" << out_suffix <<endl;
	write_ostringstream_to_file(all_out_stream, out_prefix + "all_loci" + out_suffix);
	
}

vector<int> SuperPopulation::get_locus_superpop_wide(size_t total_num_of_bacteria, size_t pop_size){
	vector<int> locus_of_bacteria;
	size_t bacterium_index = 0;
	size_t pop_index = 0;
	for(size_t i=0; i < total_num_of_bacteria; i++){
		pop_index = i/pop_size;
		bacterium_index = i%pop_size;
		size_t tmpAllele = _pops[pop_index][bacterium_index].get_allele();
		locus_of_bacteria.push_back(tmpAllele);
	}
	return locus_of_bacteria;
}

void SuperPopulation::print_locus_superpop_wide_unique(size_t total_num_of_bacteria, size_t pop_size){
	
	vector<int> locus_of_bacteria = get_locus_superpop_wide(total_num_of_bacteria, pop_size);
	
	if(PRINT_MODE){
		cout << "BEFORE UNIQUE: locus_of_bacteria.size() = " << locus_of_bacteria.size() << endl;
		cout << "Current alleles in locus are:" << endl;
		for(size_t i=0; i<locus_of_bacteria.size(); i++){
			cout << locus_of_bacteria[i] << "\t";
		}
		cout << endl;
	}

	sort(locus_of_bacteria.begin(), locus_of_bacteria.end());
	locus_of_bacteria.erase(unique(locus_of_bacteria.begin(), locus_of_bacteria.end() ), locus_of_bacteria.end() );
	
	if(PRINT_MODE){
		cout << "AFTER UNIQUE: locus_i_of_bacteria.size() = " << locus_of_bacteria.size() << endl;
		cout << "Current alleles in locus are:" << endl;
		for(size_t i=0; i<locus_of_bacteria.size(); i++){
			cout << locus_of_bacteria[i] << "\t";
		}
		cout << endl;
	}
}

void SuperPopulation::print_num_of_mutations(){
	cout << endl << "Number of mutations in locus is:"<<endl;
	cout << _num_of_different_alleles << "\t";
	cout << endl << endl;
	cout << "Max is:" <<endl;
	//cout << *max_element(_num_of_different_alleles.begin(), _num_of_different_alleles.end()) << endl << endl;
}

void SuperPopulation::check_superpop_validity(void){
	for (size_t i = 0; i < _pops.size(); i++){
		for (size_t j = 0; j < _pops[i].get_total_size(); j++){
			/*for (size_t k = 0; k < _num_of_different_alleles.size(); k++){
				Bacterium bac = _pops[i].get_bacterium(j);
				if(static_cast <int>(_num_of_different_alleles[k]) < bac.get_allele(k)){
					cout << "Invariant Error!"<< endl << "Bacterium " << j << " in pop " << i << " has:" << endl ;
					for (size_t l=0; l<_num_of_different_alleles.size(); l++){
						cout << bac.get_allele(l) << "\t";
					}
					cout << endl << "While the populations alleles are:" <<endl;
					for (size_t l=0; l<_num_of_different_alleles.size(); l++){
						cout << _num_of_different_alleles[l] << "\t";
					}
					cout << endl << "Allele " << k << " is illegal!!" << endl;
					return;
				}
			}*/
		}
	}
}