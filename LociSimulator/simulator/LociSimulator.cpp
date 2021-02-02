

#include "SuperPopulation.h"
#include "LocusTreeStructure.h"

void mutate_bacterium(Bacterium & bacterium_to_mutate, SuperPopulation & new_pops, vector< vector<int> > & mutations_history, const size_t current_bacterium_index, const size_t pop_index);
void recombine_bacterium(Population & old_pop, Bacterium & bacterium_to_recombine, const size_t current_bacterium_index, const int num_of_non_recombining_loci, const size_t pop_index);
size_t rnd_locus_index_to_recombine(size_t num_of_loci, size_t num_of_non_recombining_loci);
void recombination_between_populations(SuperPopulation & old_pops, const size_t target_pop_index, Bacterium & target_bacterium, const size_t target_bacterium_index, const size_t num_of_non_recombining_loci);
void wright_fisher_generation(SuperPopulation & old_pops, SuperPopulation & new_pops, const size_t pop_index, vector< vector<int> > & mutations_history, const double mutation_rate_per_bacterium, const double recombination_rate_per_bacterium, const double between_pops_recombination_rate_per_bacterium, const int num_of_non_recombining_loci);
void wright_fisher_process(SuperPopulation & pops, long super_generation, long num_of_generations, double mutation_rate_per_bacterium_per_generation, double recombination_rate_per_bacterium_per_generation, double between_pops_recombination_rate_per_bacterium_per_generation, const int num_of_non_recombining_loci, vector< vector<int> > & mutations_history);
//void get_statistics(Population & pop, vector<int> & F);
void print_mutation_history(const vector< vector<int> > & mutations_history);
int write_mutation_history(SuperPopulation pops, const size_t num_of_loci, vector< vector<int> > & mutations_history, string out_path);
void print_F_and_G(long generation, vector<int> & g, vector<int> & f, Population & pop, double & number_of_F_computations);
void print_and_save_parameters_values(string output_dir, const int pop_size, const int num_of_pops, const long num_of_generations, const long super_generation, const size_t num_of_loci, const size_t locus_len,  const size_t num_of_non_recombining_loci, const double p_a, const double p_c, const double p_g, const double p_transition, const double mutation_rate_per_bacterium_per_generation, const double recombination_rate_per_bacterium_per_generation, const double between_pops_recombination_rate_per_bacterium_per_generation, const double selection_intensity, const size_t selective_locus, const double selective_allele_timing_proportion, const size_t selective_allele, const int num_of_bacteria, const int expected_num_of_mutations_per_locus, const int expectedNumMutationsPerLocusPlus5Percent);
void read_mutations(vector< vector<int> > & mutations_history, size_t num_of_loci);

void read_pops(SuperPopulation & pops);

int main(int argc, char** argv){
	
	cout << endl << endl << "###############################################################################" << endl << endl;

	int pop_size; // num of bacteria in each pop
	int num_of_pops;
	double selective_allele_timing_proportion;

	long num_of_generations; // WF process length
	long super_generation; 

	size_t num_of_loci;
	size_t locus_len;
	size_t num_of_non_recombining_loci; //to be able to reconstruct a tree (as if it is a 16s gene or something similar)

	double selection_intensity; // 1 for no selection
	double p_a;
	double p_c;
	double p_g;
	double p_transition;
	double mutation_rate_per_bacterium_per_generation;
	double recombination_rate_per_bacterium_per_generation;
	double between_pops_recombination_rate_per_bacterium_per_generation;

	string output_dir;

	if(argc > 1){
		read_control_file(argv[1], pop_size, num_of_pops, num_of_generations, super_generation, num_of_loci, locus_len, num_of_non_recombining_loci, p_a, p_c, p_g, p_transition, mutation_rate_per_bacterium_per_generation, recombination_rate_per_bacterium_per_generation, between_pops_recombination_rate_per_bacterium_per_generation, selection_intensity, selective_allele_timing_proportion, output_dir);
	}
	else{
	    string default_control_file = "./test_control_file.txt";
		cout << "Using DEFAULT control file from:" << endl << default_control_file <<endl;
        read_control_file(default_control_file, pop_size, num_of_pops, num_of_generations, super_generation, num_of_loci, locus_len, num_of_non_recombining_loci, p_a, p_c, p_g, p_transition, mutation_rate_per_bacterium_per_generation, recombination_rate_per_bacterium_per_generation, between_pops_recombination_rate_per_bacterium_per_generation, selection_intensity, selective_allele_timing_proportion, output_dir);
	}

	mkdir(output_dir.c_str(), 0755);

	size_t selective_locus = num_of_non_recombining_loci; // first loci are non_recombining then comes a selective locus and rest of loci are neutral

	int num_of_bacteria = num_of_pops * pop_size;
	
	int expected_num_of_mutations_per_locus = int(std::ceil(mutation_rate_per_bacterium_per_generation * num_of_bacteria * num_of_generations / num_of_loci));

	int expected_num_of_mutations_per_locus_plus_1_percent = static_cast<int>(expected_num_of_mutations_per_locus * 1.01);

	size_t selective_allele = static_cast<int>(expected_num_of_mutations_per_locus * selective_allele_timing_proportion); // for example, if the timing is 0.95, then the selective allele number will be: 0.95 * expected number of mutations (per locus)

	print_and_save_parameters_values(output_dir, pop_size, num_of_pops, num_of_generations, super_generation, num_of_loci, locus_len, num_of_non_recombining_loci, p_a, p_c, p_g, p_transition, mutation_rate_per_bacterium_per_generation, recombination_rate_per_bacterium_per_generation, between_pops_recombination_rate_per_bacterium_per_generation, selection_intensity, selective_locus, selective_allele_timing_proportion, selective_allele, num_of_bacteria, expected_num_of_mutations_per_locus, expected_num_of_mutations_per_locus_plus_1_percent);

	srand(static_cast<unsigned int>(time(NULL)));

	vector<size_t> num_of_different_alleles(num_of_loci,0);
	
	SuperPopulation pops(num_of_pops, num_of_different_alleles, pop_size, selective_allele, selective_locus, selection_intensity);
	
	cout << "Initial SuperPops were generated." << endl;

	vector< vector<int> > mutationHistory(num_of_loci); // in each of these vectors, the i'th cell contains the allele that (after mutation) led to the i'th allele. the vector is filled with -1's to avoid reallocation. thus, after the last mutation all cells will contain -1
	
	if(READ_MUTATIONS){
		read_mutations(mutationHistory, num_of_loci);
		read_pops(pops);
	}
	else{
		cout << "Trying to initialize mutations_history... (vector of size " << expected_num_of_mutations_per_locus_plus_1_percent << ")" << endl;
		for (size_t i=0; i<mutationHistory.size(); ++i) {
			mutationHistory[i].resize(expected_num_of_mutations_per_locus_plus_1_percent,-1);
		}
		cout << "Initialization succeeded." << endl;
	
		cout << "Executing WF process..." << endl;
		wright_fisher_process(pops, super_generation, num_of_generations, mutation_rate_per_bacterium_per_generation, recombination_rate_per_bacterium_per_generation, between_pops_recombination_rate_per_bacterium_per_generation, num_of_non_recombining_loci, mutationHistory);
		cout << "Done! Simulations are ready.\n\n" << endl;
		if(DEBUG_MODE || PRINT_MUTATION_HISTORY){
			print_mutation_history(mutationHistory);
		}
		if (DEBUG_MODE || WRITE_POPS_ALLELE_NUMBERS) {
			pops.write_pops_loci(output_dir);
		}
		if (DEBUG_MODE || WRITE_MUTATION_HISTORY) {
			write_mutation_history(pops, num_of_loci, mutationHistory, output_dir);
		}
	}

	pops.print_num_of_mutations();
	
	vector<LocusTreeStructure> locus_tree_structures; // containing a "local tree" for each locus

	for(size_t i=0; i < num_of_loci; i++){ //a sequence simulator for each locus
		cout << "Extracting superpop-wide alleles of locus " << i << endl;
		vector<int> current_alleles_at_locus_i = pops.get_locus_i_superpop_wide(i, num_of_bacteria, pop_size);

		if (DEBUG_MODE) {
			pops.print_locus_i_superpop_wide_unique(i, num_of_bacteria, pop_size);
		}

		cout << "Starting construct of tree struct for locus " << i << endl;
		LocusTreeStructure lts(mutationHistory, current_alleles_at_locus_i, i, num_of_bacteria, pop_size);
		cout << "lts " << i <<" was initialized." << endl;
		locus_tree_structures.push_back(lts);
		cout << "lts " << i <<" was pushed_back." << endl;
	}


	vector<LocusSequenceSimulator> locusSequenceSimulators;
	for(size_t i=0; i < num_of_loci; i++){ //a sequence simulator for each locus
		cout << "Starting simulate sequences for locus " << i << endl;
		LocusSequenceSimulator lss(locus_tree_structures[i], locus_len);
		//cout << "lss was initialized."<< endl;
		locusSequenceSimulators.push_back(lss);
		//cout << "lss was pushed_back."<< endl;
		locusSequenceSimulators[i].simulate(p_a, p_c, p_g, p_transition);
	}
	
	pops.sample_and_write_sequences(locusSequenceSimulators, output_dir, num_of_non_recombining_loci);
	
	if(WRITE_POPS_SEQUENCES || DEBUG_MODE){
		pops.write_pops_sequences(locusSequenceSimulators, output_dir);
	}

	return 0;
}

//mutate locus random_locus_index of bacterium bac and also update num_of_different_alleles in new_pop
void mutate_bacterium(Bacterium & bacterium_to_mutate, SuperPopulation & new_pops, vector< vector<int> > & mutations_history_of_locus, const size_t current_bacterium_index, const size_t pop_index){

	size_t random_locus_index = rand() % bacterium_to_mutate.get_num_of_loci(); //locus to mutate
	size_t old_allele = bacterium_to_mutate.get_allele(random_locus_index);
	
	if(DEBUG_MODE){
		cout<<"Mutating locus " << random_locus_index << "(MUTATION: num_of_different_alleles is: ";
		for(size_t i = 0; i < new_pops._num_of_different_alleles_in_locus.size()-1; i++){
			cout << new_pops._num_of_different_alleles_in_locus[i] << "\t";
		}
		cout << new_pops._num_of_different_alleles_in_locus[new_pops._num_of_different_alleles_in_locus.size()-1] << ")" << endl;
	}

	if(new_pops.get_num_of_different_alleles_at_locus(random_locus_index) < old_allele){
		cout << "new < old!!!!"<< endl;
	}
	
	if(DEBUG_MODE){cout << "locus " << random_locus_index << " has BEFORE " << new_pops.get_num_of_different_alleles_at_locus(random_locus_index) << " different alleles" << endl;}

	new_pops.add_new_allele_to_locus(random_locus_index);
	
	if(DEBUG_MODE){cout << "locus " << random_locus_index << " has AFTER " << new_pops.get_num_of_different_alleles_at_locus(random_locus_index) << " different alleles" << endl;}
	
	size_t new_allele = new_pops.get_num_of_different_alleles_at_locus(random_locus_index);
	
	if (random_locus_index == new_pops._selective_locus && new_allele == new_pops._selective_allele) { 
		cout << endl << endl << "SELECTIVE ALLELE in population " << pop_index << " was popped out !! " << endl << endl << endl;
	}

	if(DEBUG_MODE){cout << "MUTATION: Bacterium " << current_bacterium_index << " was mutated at locus " << random_locus_index << " from " << old_allele << " to " << new_allele << endl;}
	
	if (mutations_history_of_locus[random_locus_index].size() == new_allele){
		mutations_history_of_locus[random_locus_index].push_back(-1);  // mutation history is full. reallocation is needed...
	}
	mutations_history_of_locus[random_locus_index][new_allele] = old_allele; // documenting that old_allele is the ancestor allele of new_allele (for the relevant locus)

	bacterium_to_mutate.set_allele(random_locus_index, new_allele);//set mutation
}

//replace locus random_locus_index of bacterium bac1 with the corresponding locus of bac2
void recombine_bacterium(Population & old_pop, Bacterium & bacterium_to_recombine, const size_t bacterium_to_recombine_index, const int num_of_non_recombining_loci, const size_t pop_index){
	
	size_t random_locus_index = rnd_locus_index_to_recombine(bacterium_to_recombine.get_num_of_loci(), num_of_non_recombining_loci);
	
	size_t random_other_bacterium_index = old_pop.draw_other_bacterium_index_uniformly(bacterium_to_recombine_index);
	
	Bacterium other_bacterium = old_pop.get_bacterium(random_other_bacterium_index);
	
	int old_allele = bacterium_to_recombine.get_allele(random_locus_index);

	int new_allele = other_bacterium.get_allele(random_locus_index);

	if (random_locus_index == old_pop._selective_locus && new_allele == old_pop._selective_allele && (DEBUG_MODE || TRACK_SELECTIVE_ALLELE)) {
		// note that the sample allele is the selective allele  
		cout << "Selective allele was recombined within POP " << pop_index << " from Bacterium " << random_other_bacterium_index << " to Bacterium " << bacterium_to_recombine_index << endl;
	}

	if (DEBUG_MODE){cout << "RECOMBINATION: Bacterium " << bacterium_to_recombine_index << " from Bacterium " << random_other_bacterium_index << " at locus " << random_locus_index << " received allele number " << new_allele << " (had " << old_allele << ")" <<  endl; }
	
	bacterium_to_recombine.set_allele(random_locus_index, new_allele);//set recombination
}

size_t rnd_locus_index_to_recombine(size_t num_of_loci, size_t num_of_non_recombining_loci){
	return rand() % (num_of_loci - num_of_non_recombining_loci) + num_of_non_recombining_loci; //locus to recombine
}

void recombination_between_populations(SuperPopulation & old_pops, const size_t target_pop_index, Bacterium & target_bacterium, const size_t target_bacterium_index, const size_t num_of_non_recombining_loci){
	
	if(old_pops.get_num_of_pops()<2){ //superpop is degenerate. No "other" pop to get allele from..
		return;
	}
		
	size_t random_locus_index = rnd_locus_index_to_recombine(old_pops.get_num_of_loci(), num_of_non_recombining_loci);

	// donating population
	size_t random_source_population_index = rand() % (old_pops.get_num_of_pops() - 1); // (excluding the source_pop)
	if (random_source_population_index >= target_pop_index){
		random_source_population_index++;
	}

	// donating bacterium
	size_t random_source_bacterium_index = rand() % old_pops.get_pop(random_source_population_index).get_total_size(); // index can be the same as the source bacterium since they are from different pops
	Bacterium source_bacterium = old_pops.get_pop(random_source_population_index).get_bacterium(random_source_bacterium_index);

	int old_allele = target_bacterium.get_allele(random_locus_index);

	int new_allele = source_bacterium.get_allele(random_locus_index);

	if (random_locus_index == old_pops._selective_locus && new_allele == old_pops._selective_allele && (DEBUG_MODE || TRACK_SELECTIVE_ALLELE || NOTIFY_BETWEEN_POPULATION_RECOMBINATION)) {
		// note that the sample allele is the selective allele  
		cout << "Selective allele was recombined between POPs " << random_source_population_index << " and " << target_pop_index << " from Bacterium " << random_source_bacterium_index << " to Bacterium " << target_bacterium_index << " at locus " << random_locus_index << " received allele number " << new_allele << " (had " << old_allele << ")" << endl;
	}
	
	if (DEBUG_MODE) {
		cout << "BETWEEN POPULATION RECOMBINATION: POP " << target_pop_index << " Bacterium " << target_bacterium_index << " from POP " << random_source_population_index << " Bacterium " << random_source_bacterium_index << " at locus " << random_locus_index << " received allele number " << new_allele << " (had " << old_allele << ")" << endl; 
	}

	// finalize recombination
	target_bacterium.set_allele(random_locus_index, new_allele);

}

//execute a single wright_fisher step (generation) for a given the population
void wright_fisher_generation(SuperPopulation & old_pops, SuperPopulation & new_pops, const size_t pop_index, vector< vector<int> > & mutations_history, const double mutation_rate_per_bacterium, const double recombination_rate_per_bacterium, const double between_pops_recombination_rate_per_bacterium, const int num_of_non_recombining_loci){
	double mutate, recombine, recombine_between_pops;
	if(DEBUG_MODE){cout << "BEFORE: new_pops._num_of_different_alleles[0] == " << new_pops._num_of_different_alleles_in_locus[0] << endl;}
	new_pops._num_of_different_alleles_in_locus = old_pops._num_of_different_alleles_in_locus;
	if(DEBUG_MODE){cout << "AFTER: new_pops._num_of_different_alleles[0] == " << new_pops._num_of_different_alleles_in_locus[0] << endl;}
	new_pops.get_pop(pop_index)._selective_pop_size = 0;
	new_pops.get_pop(pop_index)._neutral_pop_size = 0;
	
	if(DEBUG_MODE){new_pops.check_superpop_validity();}

	//old_pops.get_pop(pop_index).check_pop_validity();

	for (size_t child_bacterium_index=0; child_bacterium_index < old_pops.get_pop(pop_index).get_total_size(); child_bacterium_index++){

		if (DEBUG_MODE) {
			cout << "______________________________________________" << endl;
			cout << "Generating bacteria " << child_bacterium_index << endl;
		}

		int randomParentIndex = old_pops.get_pop(pop_index).draw_bacterium_index_selectively();
		
		Bacterium child_bacterium = Bacterium(old_pops.get_pop(pop_index).get_bacterium(randomParentIndex)); // MUST be a copy of the old bacterium
		
		if(DEBUG_MODE){
			cout<<"Parent is " << randomParentIndex << endl; 
			if(static_cast <int> (new_pops.get_num_of_different_alleles_at_locus(0)) < child_bacterium.get_allele(0)){
				cout << "new < old!!!!"<< endl;
			}
		}

		//mutation
		mutate = draw_from_zero_to_one_in_high_resolution_steps(4);
		if(mutate < mutation_rate_per_bacterium){
			mutate_bacterium(child_bacterium, new_pops, mutations_history, child_bacterium_index, pop_index);
		}
		
		//recombination
		recombine = draw_from_zero_to_one_in_high_resolution_steps(4);
		if(recombine < recombination_rate_per_bacterium){
			recombine_bacterium(old_pops.get_pop(pop_index), child_bacterium, child_bacterium_index, num_of_non_recombining_loci, pop_index);
		}

		//between population recombination
		recombine_between_pops = draw_from_zero_to_one_in_high_resolution_steps(4);
		if(recombine_between_pops < between_pops_recombination_rate_per_bacterium){
			recombination_between_populations(old_pops, pop_index, child_bacterium, child_bacterium_index, num_of_non_recombining_loci);
		}
		
		new_pops.get_pop(pop_index).add_bacterium(child_bacterium);
		
		if(DEBUG_MODE){
			new_pops.check_superpop_validity();
			child_bacterium.print(child_bacterium_index);
		}
	}
}

//execute a wright_fisher process with num_of_generation generations
void wright_fisher_process(SuperPopulation & oldSuperPopulation, long super_generation, long num_of_generations, double mutation_rate_per_bacterium_per_generation, double recombination_rate_per_generation, double between_pops_recombination_rate_per_bacterium_per_generation, const int num_of_non_recombining_loci, vector< vector<int> > & mutations_history){
	SuperPopulation newSuperPopulation(oldSuperPopulation.get_num_of_pops(), oldSuperPopulation._num_of_different_alleles_in_locus,oldSuperPopulation.get_pop(0).get_total_size(), oldSuperPopulation._selective_allele, oldSuperPopulation._selective_locus, oldSuperPopulation._selection_intensity); // starting an empty new_populations vector
	long generation = 0;
	//vector<int> F(num_of_loci+1,0);
	//vector<int> G(F.size(), 0);
	//double number_of_F_computations = 0.0;
	while(generation < num_of_generations){
		generation += 2;
		//print_F_and_G(generation, F, G, pop, number_of_F_computations);
		
		if (generation % super_generation == 0){
			
			for (size_t i=0; i < oldSuperPopulation.get_num_of_pops(); i++){ // old_pops is more updated now. thus the super_generation is applied on old_pops
				
				size_t randomParentPopulationIndex = newSuperPopulation.draw_pop_index_uniformly();
		
				oldSuperPopulation.set_pop(newSuperPopulation.get_pop(randomParentPopulationIndex), i);

			}
			if(DEBUG_MODE){newSuperPopulation.check_superpop_validity();}
			//old_pops = new_pops;
			//old_pops._num_of_different_alleles = old_pops._num_of_different_alleles;
			//new_pops.check_superpop_validity();
		}

		// apply a Wright Fisher process over each population
		for(size_t i=0; i<oldSuperPopulation.get_num_of_pops(); i++){
			wright_fisher_generation(oldSuperPopulation, newSuperPopulation, i, mutations_history, mutation_rate_per_bacterium_per_generation, recombination_rate_per_generation, between_pops_recombination_rate_per_bacterium_per_generation, num_of_non_recombining_loci);
			if(DEBUG_MODE){
				cout << "generation = " << generation << endl;
				newSuperPopulation.get_pop(i).print();
				// print_mutation_history(mutations_history);
			}
			wright_fisher_generation(newSuperPopulation, oldSuperPopulation, i, mutations_history, mutation_rate_per_bacterium_per_generation, recombination_rate_per_generation, between_pops_recombination_rate_per_bacterium_per_generation, num_of_non_recombining_loci);
			// oldSuperPop is now the "new" superPop...
			if(DEBUG_MODE){
				cout << "generation = " << generation << endl;
				oldSuperPopulation.get_pop(i).print();
				// print_mutation_history(mutations_history);
			}
		}
		if (generation % 1000 == 0) {
			cout << "Wright Fisher generation number " << generation << ". Current alleles status is ";
			for (size_t i = 0; i < oldSuperPopulation.get_num_of_loci(); i++) {
				cout << oldSuperPopulation.get_num_of_different_alleles_at_locus(i) << " ";
			}
			cout << endl;
		}
	}
	//get_statistics(pop, F);
}

//calculate F_i_k of a given population
void get_statistics(Population & pop, vector<int> & F){
	for (size_t i=0; i < F.size(); ++i) {
		F[i]=0;
	}
	for(size_t i=0; i < pop.get_total_size() - 1; i++){
		for(size_t j=i+1; j < pop.get_total_size(); j++){
			int count = pop.get_bacterium(i).diff(pop.get_bacterium(j));
			F[count]++;
		}
	}
}

//print mutation_history to the screen
void print_mutation_history(const vector< vector<int> > & mutations_history){
	cout << endl <<  "Mutation history:" << endl;
	for(size_t i=0; i < mutations_history.size(); i++){
		cout << "Locus " << i << ":" << endl;
		for (size_t j=0; j < mutations_history[i].size(); j++){
			cout << mutations_history[i][j] << "\t";
		}
		cout << endl;
	}
}

//write mutation_history of each locus to its own output file
int write_mutation_history(const SuperPopulation pops, const size_t num_of_loci, vector< vector<int> > & mutations_history, string output_dir){
	//for each locus, mutation history is saved in the following manner:
	//the i'th row in the file contains the allele from which the i'th allele was derived.
	//if it has -1 it means that this allele was one of the alleles in the first generation.
	ofstream myFile;
	ostringstream out_stream_path, out_stream_string;

	// create a dedicated folder
	output_dir = output_dir + "/mutation_histories/";
	mkdir(output_dir.c_str(), 0755);

	int old_allele = -1;
	cout << "Writing mutation histories..." << endl;
	for(size_t i=0; i < num_of_loci; i++){
		out_stream_path << output_dir << "simulatedMutationsHistory.locus." << i << ".txt";
		cout << "Writing mutations history of locus " << i << " to:" <<endl << out_stream_path.str() << endl;
		for (size_t j=0; j<=pops._num_of_different_alleles_in_locus[i]; j++){
			old_allele = mutations_history[i][j];
			out_stream_string << old_allele << endl;
		}
		write_ostringstream_to_file(out_stream_string, out_stream_path.str());
		out_stream_path.str(""); //clean stream
	}
	cout << "Done writing mutations history.\n\n"<<endl;
	return 1;
}


void print_F_and_G(long generation, vector<int> & G, vector<int> & F, Population & pop, double & number_of_F_computations){
	if(!DEBUG_MODE && generation%100000==0 && generation>500000 ){
		cout << "generation = " << generation << endl;
		get_statistics(pop, F);
		number_of_F_computations++;

		cout<<"F: \t";
		for (size_t i=0; i < F.size(); ++i) {
			cout<<F[i]<<"\t";
		}
		cout<<endl;

		cout<<"G: \t";
		for (size_t i=0; i < F.size(); ++i) {
			G[i] += F[i];
			cout<<G[i]/number_of_F_computations<<"\t";
		}
		cout<<endl;
	}
}

void read_mutations(vector< vector<int> > & mutations_history, size_t num_of_loci){
	int value, number_of_lines, j;
    string line;
	ostringstream in_path;
	ifstream myfile;
	for (size_t i=0; i < num_of_loci; i++){
		number_of_lines = 0;
		j = 0;
		in_path << "D:\\GoogleDrive\\SelectiveSweepProject\\simulations\\FraserEtAl\\Simulations\\simulatedMutationsHistory.DEBUG.locus." << i << ".txt";
		cout << "Reading mutations history of locus " << i << endl;
		myfile.open(in_path.str());
		while (getline(myfile, line)){
			++number_of_lines;
		}
		cout << "Number of lines in " << in_path.str() << " is " << number_of_lines << endl;
		myfile.close();
		mutations_history[i].resize(number_of_lines,-1);
		myfile.open(in_path.str());
		while (myfile >> value ) {
			mutations_history[i][j] = value;
			j++;
		}
		in_path.str("");
		myfile.close();
	}
	cout << "Done!" << endl;
}


void read_pops(SuperPopulation & pops){
	int number_of_lines, j;
    string line;
	ostringstream in_path;
	ifstream myfile;
	for (size_t pop_index=0; pop_index < pops.get_num_of_pops(); pop_index++){
		Population pop = pops.get_pop(pop_index);
		number_of_lines = 0;
		j = 0;
		in_path << "D:\\GoogleDrive\\SelectiveSweepProject\\simulations\\FraserEtAl\\Simulations\\simulatedLoci.DEBUG.pop." << pop_index << ".txt";
		cout << "Reading pop " << pop_index << endl;
		myfile.open(in_path.str());
		while (getline(myfile, line)){
			++number_of_lines;
		}
		myfile.close();
		cout << "Number of lines in " << in_path.str() << " is " << number_of_lines << " (Number of bacteria in pop is " << pop.get_total_size() << ")" << endl;
		pop._neutral_pop_size = 0;
		pop._selective_pop_size = 0;
		myfile.open(in_path.str());
		while (myfile >> line ) {
			if(j % 2 ==1){
				vector<int> loci(pops.get_num_of_loci());
				string allele = "";
				int locus_index=0;
				for(size_t c=0; c<line.size(); c++){
					if (isdigit(line[c])){
						allele = allele + line[c];
					}
					else
					{
						allele = "";
						locus_index++;
					}
				}
				loci[locus_index] = stoi(allele);
				Bacterium bac(loci);
				pop.add_bacterium(bac);
			}
			j++;
		}
		in_path.str("");
		myfile.close();
		pops._pops[pop_index] = pop;
	}
	cout << "Done!" << endl;
	
}

void print_and_save_parameters_values(string output_dir, const int pop_size, const int num_of_pops, const long num_of_generations, const long super_generation, const size_t num_of_loci, const size_t locus_len,  const size_t num_of_non_recombining_loci, const double p_a, const double p_c, const double p_g, const double p_transition, const double mutation_rate_per_bacterium_per_generation, const double recombination_rate_per_bacterium_per_generation, const double between_pops_recombination_rate_per_bacterium_per_generation, const double selection_intensity, const size_t selective_locus, const double selective_allele_timing_proportion, const size_t selective_allele, const int num_of_bacteria, const int expected_num_of_mutations_per_locus, const int expected_num_of_mutations_per_locus_plus_1_percent){
	ostringstream out_stream_string;

	out_stream_string << "Parameters values are:" << endl;
	out_stream_string << "population size = " << pop_size << endl;
	out_stream_string << "number of populations = " << num_of_pops << endl;
	out_stream_string << "number of generations = " << num_of_generations << endl;
	out_stream_string << "super generation = " << super_generation << endl;
	out_stream_string << "number of loci = " << num_of_loci << endl; 
	out_stream_string << "locus length = " << locus_len << endl;
	out_stream_string << "number of non recombining loci = " << num_of_non_recombining_loci << endl;
	out_stream_string << "p_a = " << p_a << endl;
	out_stream_string << "p_c = " << p_c << endl;
	out_stream_string << "p_g = " << p_g << endl;
	out_stream_string << "p_transition = " << p_transition << endl;
	out_stream_string << "selective allele timing proportion allele = " << selective_allele_timing_proportion << endl;
	out_stream_string << "selective locus = " << selective_locus << endl;
	out_stream_string << "selection intensity = " << selection_intensity << endl;
	out_stream_string << "mutation rate per bacterium per generation = " << mutation_rate_per_bacterium_per_generation << endl;
	out_stream_string << "recombination rate per bacterium per generation = " << recombination_rate_per_bacterium_per_generation << endl;
	out_stream_string << "between pops recombination rate per bacterium per generation = " << between_pops_recombination_rate_per_bacterium_per_generation << endl;
	out_stream_string << "output directory = " << output_dir << endl;

	out_stream_string << "total number of simulated bacteria = " << num_of_bacteria << endl;
	out_stream_string << "expected num of mutations per locus = " << expected_num_of_mutations_per_locus << endl;
	out_stream_string << "selective allele number = " << selective_allele << endl;
	out_stream_string << "selective allele will pop out after " << selective_allele_timing_proportion * 100 << "% of the expected number of mutations" << endl;
	out_stream_string << "expected num of mutations per locus + 1% = " << expected_num_of_mutations_per_locus_plus_1_percent << endl;
	out_stream_string << "expected total number of mutations events = " << expected_num_of_mutations_per_locus*num_of_loci << endl;
	out_stream_string << "expected total number of recombinations events (within pops) = " << int(std::ceil(recombination_rate_per_bacterium_per_generation * num_of_bacteria * num_of_generations)) << endl;
	double expected_generation_in_which_the_selective_allele_will_pop_out = num_of_generations * selective_allele / expected_num_of_mutations_per_locus;
	out_stream_string << "expected generation in which the selective allele will pop out = " << expected_generation_in_which_the_selective_allele_will_pop_out << endl;
	int expected_total_num_of_recombinations_events_between_pops = int(std::ceil(between_pops_recombination_rate_per_bacterium_per_generation * num_of_bacteria * num_of_generations));
	out_stream_string << "expected total number of recombinations events (between pops) = " << expected_total_num_of_recombinations_events_between_pops << endl;
	int num_of_generations_after_selective_allele_was_popped_out = num_of_generations - expected_generation_in_which_the_selective_allele_will_pop_out;
	out_stream_string << "number of generations after selective allele was popped out = " << num_of_generations_after_selective_allele_was_popped_out << endl;
	out_stream_string << "expected total number of recombinations events (between pops) after selective allele was popped out = " << int(std::ceil(between_pops_recombination_rate_per_bacterium_per_generation * num_of_bacteria * num_of_generations_after_selective_allele_was_popped_out)) << endl;

	cout << out_stream_string.str();

	write_ostringstream_to_file(out_stream_string, output_dir + "/parameters.txt");
}
