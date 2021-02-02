#include "SuperPopulation.h"



SuperPopulation::SuperPopulation(size_t num_of_pops, vector<size_t> num_of_different_alleles, size_t pop_size, size_t selective_allele, size_t selective_locus, double selection_intensity) :
	_num_of_different_alleles_in_locus(num_of_different_alleles), _selective_allele(selective_allele), _selective_locus(selective_locus), _selection_intensity(selection_intensity){
	for(size_t i=0; i < num_of_pops; i++){
		_pops.push_back(Population(num_of_different_alleles, pop_size, selective_allele, selective_locus, selection_intensity));
	}
}

int SuperPopulation::write_pops_loci(string output_dir) {
	ostringstream out_path;

	// create a dedicated folder
	output_dir = output_dir + "/pop_allele_numbers";
	mkdir(output_dir.c_str(), 0755);

	for (size_t i=0; i<get_num_of_pops(); i++){
		out_path << output_dir << "/simulatedLoci.pop." << i << ".txt";
		get_pop(i).write_pop_loci(out_path.str());
		out_path.str(""); //clean stream
	}
	return 1;
}

int SuperPopulation::write_pops_sequences(const vector<LocusSequenceSimulator> & locus_sequences_simulators, string output_dir) {
	ostringstream out_path;

	// create a dedicated folder
	output_dir = output_dir + "/populations_sequences";
	mkdir(output_dir.c_str(), 0755);

	for (size_t i=0; i<get_num_of_pops(); i++){
		out_path << output_dir << "/pop." << i << ".txt";
		get_pop(i).write_pop_sequences(locus_sequences_simulators, out_path.str(), i);
		out_path.str(""); //clean stream
	}
	return 1;
}

//this function samples a bacterium from each population at the and of the evolution simulation and write the sampled bacteria sequences to all_loci_outpath.
//in addition it write to non_recombining_loci_outpath all the non recombining loci of the sampled bacteria (i.e. a SUBsequence of the full "genome"
int SuperPopulation::sample_and_write_sequences(const vector<LocusSequenceSimulator> & locus_sequences_simulators, const string output_dir, size_t num_of_non_recombining_loci) {
	ostringstream allOutStream, nonRecombiningOutStream, selective_allele_stats_stream;
	size_t allele;
	size_t selectiveAlleleCount = 0;
	if (get_num_of_pops() == 1){
		cout << endl << "Num of populations is 1.  Writing the whole population instead of sampling." <<endl;
		
		Population pop = get_pop(0);
		
		for (size_t i=0; i<pop.get_total_size(); i++){
		
			if(i%25 == 0){
				cout << endl << "Writing sequence " << i << "." << endl;
			}
			if(DEBUG_MODE){cout << "The sequence is:" <<endl;}
		
			allOutStream << ">0_" << i << endl;
			nonRecombiningOutStream << ">0_" << i << endl;
		
			for(size_t j=0; j < get_num_of_loci(); j++){

				allele = pop.get_bacterium(i).get_allele(j);
				LocusSequenceSimulator locus_sequences_simulator = locus_sequences_simulators[j];
			
				if (locus_sequences_simulator._sequences[allele].size() != locus_sequences_simulator._locus_len && locus_sequences_simulator._sequences[allele].size() != 0){
					cout << endl << "error in locus " << j << " !!! (single population run)"<< endl;
					cout << "locus_sequences_simulator._sequences[" << allele << "] is " << locus_sequences_simulator._sequences[allele] << endl;
					cout << "locus_sequences_simulator._sequences[" << allele << "].size() is " << locus_sequences_simulator._sequences[allele].size() << endl;
					cout << "locus_sequences_simulator._locus_len is " << locus_sequences_simulator._locus_len << endl;
				}
				if(DEBUG_MODE){cout << locus_sequences_simulator._sequences[allele];}
				allOutStream << locus_sequences_simulator._sequences[allele];
				if (j<num_of_non_recombining_loci){
					nonRecombiningOutStream << locus_sequences_simulator._sequences[allele];
				}
			}
			allOutStream << endl;
			nonRecombiningOutStream << endl;
		}
	}

	else{
		for(size_t popIndex=0; popIndex < get_num_of_pops(); popIndex++){
			Population pop = get_pop(popIndex);
		
			size_t randomBacteriumIndex = rand() % pop.get_total_size();
		
			cout << endl << "Sequence " << randomBacteriumIndex << " was drawn from pop " << popIndex << endl;
			if(DEBUG_MODE){cout << "The sequence is:" <<endl;}
		
			// bacterium name
            allOutStream << ">" << popIndex << "_" << randomBacteriumIndex ;
            nonRecombiningOutStream << ">" << popIndex << "_" << randomBacteriumIndex ;

            allele = pop.get_bacterium(randomBacteriumIndex).get_allele(_selective_locus);
            if (allele == _selective_allele) {
				// note that the sample allele is the selective allele 
				allOutStream << "_selective" ;
                nonRecombiningOutStream << "_selective" ;
				selectiveAlleleCount ++;
			}

			for (size_t locusIndex = 0; locusIndex < get_num_of_loci(); locusIndex++) {
				allele = pop.get_bacterium(randomBacteriumIndex).get_allele(locusIndex);
				allOutStream << "_" << allele ;
                // headers of the non_recombining should be identical!
                // (although not whole sequence is written, only non-recombining alleles...)
                nonRecombiningOutStream << "_" << allele ;
			}

			allOutStream << endl ;
            nonRecombiningOutStream << endl ;

			for(size_t locusIndex=0; locusIndex < get_num_of_loci(); locusIndex++){

				allele = pop.get_bacterium(randomBacteriumIndex).get_allele(locusIndex);
				LocusSequenceSimulator locus_sequences_simulator = locus_sequences_simulators[locusIndex];
			
				if (locus_sequences_simulator._sequences[allele].size() != locus_sequences_simulator._locus_len && locus_sequences_simulator._sequences[allele].size() != 0){
					cout << endl << "error in locus " << locusIndex << " !!!" << endl;
					cout << "locus_sequences_simulator._sequences[" << allele << "] is " << locus_sequences_simulator._sequences[allele] << endl;
					cout << "locus_sequences_simulator._sequences[" << allele << "].size() is " << locus_sequences_simulator._sequences[allele].size() << endl;
					cout << "locus_sequences_simulator._locus_len is " << locus_sequences_simulator._locus_len << endl;
				}
				if (allele == _selective_allele) { cout << "Selective allele (" << _selective_allele << ") was sampled for population " << popIndex << " (bacteria " << randomBacteriumIndex << ")."; }
				if (DEBUG_MODE) { cout << locus_sequences_simulator._sequences[allele]; }
				
				// write sequence to all_loci file
				allOutStream << locus_sequences_simulator._sequences[allele];
				
				// write sequence to non_recombining_loci file if the locus is non_recombining (one of the first loci)
				if (locusIndex < num_of_non_recombining_loci){
					nonRecombiningOutStream << locus_sequences_simulator._sequences[allele];
				}
			}
			allOutStream << endl;
			nonRecombiningOutStream << endl;
		}

		cout << endl << "###############################################################################" << endl;
		cout << "@ " << 100.0 * selectiveAlleleCount / get_num_of_pops() << "% of the sampled populations contains the selective allele (" << selectiveAlleleCount << " out of " << get_num_of_pops() << ")";
		cout << endl << "###############################################################################" << endl << endl;

		string selective_allele_stats_path = output_dir + "/selective_allele_proportion.txt";
		selective_allele_stats_stream << 1.0 * selectiveAlleleCount / get_num_of_pops() << endl;
		write_ostringstream_to_file(selective_allele_stats_stream, selective_allele_stats_path);
	}

	string all_loci_path = output_dir + "/all_loci.txt";
	cout << endl << "Writing all loci sequences to:" << endl << all_loci_path <<endl;
	write_ostringstream_to_file(allOutStream, all_loci_path);

	string non_recombining_loci_path = output_dir + "/non_recombining_loci.txt";
	cout << endl << "Writing non_recombining loci sequences (for tree reconstruction) to:" << endl << non_recombining_loci_path << endl;
	write_ostringstream_to_file(nonRecombiningOutStream, non_recombining_loci_path);

	return 1;
}

vector<int> SuperPopulation::get_locus_i_superpop_wide(size_t locus_index, size_t num_of_bacteria, size_t pop_size){
	vector<int> locus_i_of_bacteria(num_of_bacteria, 0);
	size_t bacterium_index = 0;
	size_t pop_index = 0;
	for(size_t i=0; i<num_of_bacteria; i++){
		pop_index = i/pop_size;
		bacterium_index = i%pop_size;
		locus_i_of_bacteria[i] = get_pop(pop_index).get_bacterium(bacterium_index).get_allele(locus_index);
	}
	return locus_i_of_bacteria;
}

void SuperPopulation::print_locus_i_superpop_wide_unique(size_t locus_index, size_t num_of_bacteria, size_t pop_size){
	
	vector<int> locus_i_of_bacteria = get_locus_i_superpop_wide(locus_index, num_of_bacteria, pop_size);
	
	/*
	cout << "BEFORE UNIQUE: locus_i_of_bacteria.size() = " << locus_i_of_bacteria.size() << endl;
	cout << "Current alleles in locus " << locus_index << " are:" << endl;
	for(size_t i=0; i<locus_i_of_bacteria.size(); i++){
		cout << locus_i_of_bacteria[i] << "\t";
	}
	cout << endl;
	*/
	sort(locus_i_of_bacteria.begin(), locus_i_of_bacteria.end());
	locus_i_of_bacteria.erase(unique(locus_i_of_bacteria.begin(), locus_i_of_bacteria.end() ), locus_i_of_bacteria.end() );
	
	cout << "UNIQUE: locus_i_of_bacteria.size() = " << locus_i_of_bacteria.size() << endl;
	cout << "Current alleles in locus " << locus_index << " are:" << endl;
	for(size_t i=0; i<locus_i_of_bacteria.size(); i++){
		cout << locus_i_of_bacteria[i] << "\t";
	}
	cout << endl;
}

void SuperPopulation::print_num_of_mutations(void){
	cout << endl << "Number of mutations in each locus is: ";
	for(size_t i = 0; i < _num_of_different_alleles_in_locus.size(); i++){
		cout << _num_of_different_alleles_in_locus[i] << "\t";
	}
	cout << endl << "Max is: " << *max_element(_num_of_different_alleles_in_locus.begin(), _num_of_different_alleles_in_locus.end()) << endl << endl;
}

void SuperPopulation::check_superpop_validity(void){
	for (size_t i = 0; i < _pops.size(); i++){
		for (size_t j = 0; j < _pops[i].get_total_size(); j++){
			for (size_t k = 0; k < _num_of_different_alleles_in_locus.size(); k++){
				Bacterium bac = _pops[i].get_bacterium(j);
				if(static_cast <int>(_num_of_different_alleles_in_locus[k]) < bac.get_allele(k)){
					cout << "Invariant Error!"<< endl << "Bacterium " << j << " in pop " << i << " has:" << endl ;
					for (size_t l=0; l<_num_of_different_alleles_in_locus.size(); l++){
						cout << bac.get_allele(l) << "\t";
					}
					cout << endl << "While the populations alleles are:" <<endl;
					for (size_t l=0; l<_num_of_different_alleles_in_locus.size(); l++){
						cout << _num_of_different_alleles_in_locus[l] << "\t";
					}
					cout << endl << "Allele " << k << " is illegal!!" << endl;
					return;
				}
			}
		}
	}
}