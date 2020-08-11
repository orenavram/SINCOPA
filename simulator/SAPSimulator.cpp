#include "SuperPopulation.h"
#include "LocusTreeStructure.h"
#include "WrightFisherProcess.h"

void print_mutation_history(const vector<int> & mutations_history);
int write_mutation_history(const SuperPopulation & pops, vector<int> & mutations_history, const string out_path);
void print_F_and_G(long step, vector<int> & g, vector<int> & f, Population & pop, double & number_of_F_computations);
void print_parameters_values(const int pop_size, const int num_of_pops, const long max_num_of_steps, const long super_step, const size_t locus_len, const double p_a, const double p_c, const double p_g, const double p_transition, const double mutation_rate_per_bacterium_per_step, const double recombination_rate_per_bacterium_per_step, const double between_pops_recombination_rate_per_bacterium_per_step, const double selection_intensity, const size_t selective_allele, const int total_num_of_bacteria, const int expected_num_of_mutations_per_locus, const int expectedNumMutationsPerLocusPlus5Percent, const string sampled_sequences_out_prefix);
void write_parameters_values(string out_path, const int pop_size, const int num_of_pops, const long max_num_of_steps, const long super_step, const size_t locus_len, const double p_a, const double p_c, const double p_g, const double p_transition, const double mutation_rate_per_bacterium_per_step, const double recombination_rate_per_bacterium_per_step, const double between_pops_recombination_rate_per_bacterium_per_step, const double selection_intensity, const size_t selective_allele, const int total_num_of_bacteria, const int expected_num_of_mutations_per_locus, const int expectedNumMutationsPerLocusPlus5Percent, const string sampled_sequences_out_prefix);
void read_mutations(vector<int> & mutations_history);
void read_pops(SuperPopulation & pops);



int main(int argc, char** argv){

	cout << endl << endl << "###############################################################################" << endl << endl;

	int pop_size; // num of bacteria in each pop
	int num_of_pops;
	int selective_allele;
	long max_num_of_steps; // WF process length
	long super_step; 
	size_t locus_len;
	double selection_intensity; // 1 for no selection
	double p_a;
	double p_c;
	double p_g;
	double p_transition;
	double mutation_rate_per_bacterium_per_step;
	double recombination_rate_per_bacterium_per_step;
	double between_pops_recombination_rate_per_bacterium_per_step;
	string sampled_sequences_out_prefix;
	string out_prefix;
	string out_suffix;
	string done_path;
	
	_generator1.seed(0); // fixed seed 

	if(argc > 1){
		read_control_file(argv[1], pop_size, num_of_pops, max_num_of_steps, super_step, locus_len, p_a, p_c, p_g, p_transition, mutation_rate_per_bacterium_per_step, recombination_rate_per_bacterium_per_step, between_pops_recombination_rate_per_bacterium_per_step, selection_intensity, selective_allele, sampled_sequences_out_prefix, done_path);
		out_prefix = sampled_sequences_out_prefix;
		out_suffix = ".DEBUG.txt";
	}
	else{
		cout << "Reading DEFAULT control file..." <<endl;
		read_control_file(
			"C:\\Tal\\Dropbox\\workWithStudents\\Oren\\SAPSimulator\\control_file_new.txt",
			pop_size, 
			num_of_pops,
			max_num_of_steps, super_step, locus_len, p_a, p_c, p_g, p_transition, 
			mutation_rate_per_bacterium_per_step,
			recombination_rate_per_bacterium_per_step, between_pops_recombination_rate_per_bacterium_per_step, selection_intensity, selective_allele, sampled_sequences_out_prefix, done_path);
		//out_prefix = "D:\\Dropbox\\Sweeps Project\\SAPSimulator\\newVersion\\1.";
		out_prefix = "C:\\Tal\\Dropbox\\workWithStudents\\Oren\\SAPSimulator\\newVersion2\\1.";
		out_suffix = ".DEBUG.txt";
	}

	int total_num_of_bacteria = num_of_pops * pop_size;
	
	int expected_num_of_mutations_per_locus = int(std::ceil(mutation_rate_per_bacterium_per_step * total_num_of_bacteria * max_num_of_steps));

	int expectedNumMutationsPerLocusPlus5Percent = static_cast<int>(expected_num_of_mutations_per_locus * 1.05);

	print_parameters_values(pop_size, num_of_pops, max_num_of_steps, super_step, locus_len, p_a, p_c, p_g, p_transition, mutation_rate_per_bacterium_per_step, recombination_rate_per_bacterium_per_step, between_pops_recombination_rate_per_bacterium_per_step, selection_intensity, selective_allele, total_num_of_bacteria, expected_num_of_mutations_per_locus, expectedNumMutationsPerLocusPlus5Percent, sampled_sequences_out_prefix);
	write_parameters_values(out_prefix, pop_size, num_of_pops, max_num_of_steps, super_step, locus_len, p_a, p_c, p_g, p_transition, mutation_rate_per_bacterium_per_step, recombination_rate_per_bacterium_per_step, between_pops_recombination_rate_per_bacterium_per_step, selection_intensity, selective_allele, total_num_of_bacteria, expected_num_of_mutations_per_locus, expectedNumMutationsPerLocusPlus5Percent, sampled_sequences_out_prefix);
	
	srand(static_cast<unsigned int>(time(NULL)));

	size_t num_of_different_alleles(1);
	
	SuperPopulation pops(	num_of_pops,
							num_of_different_alleles,
							pop_size,
							selective_allele,
							selection_intensity);
	
	cout << "number of different alleles in pops = " << pops.get_num_of_different_alleles_at_locus() << endl;

	cout << "Initial SuperPops were generated." << endl;
	//pops.print();
	//vector<int> mutations_history; //the i'th cell contains the allele from which the i'th allele was derived. the vector is filled with -1's to avoid reallocation. thus, after the last mutation all cells will contain -1
	
	//if(READ_MUTATIONS){
	//	read_mutations(mutations_history);
	//	read_pops(pops);
	//}
	//else{
	//mutations_history.resize(expectedNumMutationsPerLocusPlus5Percent,-1);
	
	cout << "###############################################################################" << endl << endl; 
	cout << "Executing WF process" << endl;

	WrightFisherProcess wfp(pops,
							super_step,
							max_num_of_steps,
							mutation_rate_per_bacterium_per_step,
							recombination_rate_per_bacterium_per_step,
							between_pops_recombination_rate_per_bacterium_per_step);
	wfp.runWF();
	cout << "###############################################################################" << endl << endl;
	cout << "Done! Simulations are ready." << endl;
	if(PRINT_MODE || PRINT_MUTATION_HISTORY){
			print_mutation_history(wfp._mutations_history);
	}
	if(DEBUG_MODE || PRINT_MODE || WRITE_MUTATION_HISTORY){
		wfp._new_pops->write_pops_loci(out_prefix, out_suffix);
		write_mutation_history(*wfp._new_pops, wfp._mutations_history, out_prefix);
	}

	wfp._new_pops->print_num_of_mutations();
		
	vector<int> current_alleles = wfp._new_pops->get_locus_superpop_wide(total_num_of_bacteria, pop_size);

	if(DEBUG_MODE){
		wfp._new_pops->print_locus_superpop_wide_unique(total_num_of_bacteria, pop_size);
	}

	cout << "Starting construct of tree struct" << endl;
	LocusTreeStructure lts(wfp._mutations_history, current_alleles, total_num_of_bacteria, pop_size);
	cout << "lts was initialized." << endl;

	cout << "Starting simulate sequences for locus" << endl;
	LocusSequencesSimulator lss(lts, locus_len);
	lss.simulate(p_a, p_c, p_g, p_transition);
	
	wfp._new_pops->sample_and_write_sequences(lss, sampled_sequences_out_prefix, out_suffix);
	
	if(DEBUG_MODE){
		wfp._new_pops->write_pops_sequences(lss, out_prefix, out_suffix);
	}


	/*size_t num_of_simulations = 20;
	for(size_t i=0; i < num_of_simulations; i++){
	
		wright_fisher_process(pop, max_num_of_steps, mutation_rate_per_bacterium_per_step, recombination_rate_per_bacterium_per_step);

		ostringstream out_path;
		out_path << "D:\\GoogleDrive\\SelectiveSweepProject\\simulations\\FraserEtAl\\Simulations\\simulatedNeisseriaLoci\\
			simulatedNeisseriaLoci." << i << ".txt";
		write_pop(pop, out_path.str());

	}*/

	// generating a dummy file at the end of the process
	ostringstream out_stream_string;
	cout << "Touching..." << endl;
	out_stream_string << endl;
	write_ostringstream_to_file(out_stream_string, done_path);
	
	cout << "Done." << endl;

	return 0;
}

//mutate locus of bacterium bac and also update num_of_different_alleles in new_pop


//replace locus of bacterium bac1 with the corresponding locus of bac2









////calculate F_i_k of a given population
//void get_statistics(Population & pop, vector<int> & F){
//	for (size_t i=0; i < F.size(); ++i) {
//		F[i]=0;
//	}
//	for(size_t i=0; i < pop.get_total_size() - 1; i++){
//		for(size_t j=i+1; j < pop.get_total_size(); j++){
//			int count = pop.get_bacterium(i).diff(pop.get_bacterium(j));
//			F[count]++;
//		}
//	}
//}

//print mutation_history to the screen
void print_mutation_history(const vector<int> & mutations_history){
	cout << endl <<  "Mutation history:" << endl;
	for (size_t j=0; j < mutations_history.size(); j++){
		cout << mutations_history[j] << "\t";
	}
	cout << endl;
}

//write mutation_history of each locus to its own output file
int write_mutation_history(const SuperPopulation & pops, vector<int> & mutations_history, const string out_path){
	//mutation history is saved in the following manner:
	//the i'th row in the file contains the allele from which the i'th allele was derived.
	//if it has -1 it means that this allele was one of the alleles in the first generation.
	ofstream myFile;
	ostringstream out_stream_path, out_stream_string;
	int old_allele = -1;
	cout << "Starting to write mutations history..." << endl;
	out_stream_path << out_path << "simulatedMutationHistory.txt";
	cout << "Writing mutations history to:" <<endl << out_stream_path.str() << endl;
	for (size_t j=0; j<pops._num_of_different_alleles; j++){
		old_allele = mutations_history[j];
		out_stream_string << old_allele << endl;
	}
	write_ostringstream_to_file(out_stream_string, out_stream_path.str());
	out_stream_path.str(""); //clean stream
	cout << "Done writing mutations history."<<endl;
	return 1;
}


void print_F_and_G(long step, vector<int> & G, vector<int> & F, Population & pop, double & number_of_F_computations){
	if(!PRINT_MODE && step%100000==0 && step>500000 ){
		cout << "step = " << step << endl;
		//get_statistics(pop, F);
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

void read_mutations(vector<int> & mutations_history){
	int value, number_of_lines, j;
    string line;
	ostringstream in_path;
	ifstream myfile;
	number_of_lines = 0;
	j = 0;
	in_path << "D:\\GoogleDrive\\SelectiveSweepProject\\simulations\\FraserEtAl\\Simulations\\simulatedMutationsHistory.DEBUG.txt";
	cout << "Reading mutations history of locus" << endl;
	myfile.open(in_path.str());
	while (getline(myfile, line)){
		++number_of_lines;
	}
	cout << "Number of lines in " << in_path.str() << " is " << number_of_lines << endl;
	myfile.close();
	mutations_history.resize(number_of_lines,-1);
	myfile.open(in_path.str());
	while (myfile >> value ) {
		mutations_history[j] = value;
		j++;
	}
	in_path.str("");
	myfile.close();
	
	cout << "Done reading mutations!" << endl;
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
				int locus = 0;
				string allele = "";
				for(size_t c=0; c<line.size(); c++){
					if (isdigit(line[c])){
						allele = allele + line[c];
					}
					else
					{
						allele = "";
					}
				}
				locus = stoi(allele);
				Bacterium bac(locus);
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

void print_parameters_values(	const int pop_size,
								const int num_of_pops,
								const long max_num_of_steps, 
								const long super_step, 
								const size_t locus_len,
								const double p_a, 
								const double p_c,
								const double p_g,
								const double p_transition, 
								const double mutation_rate_per_bacterium_per_step,
								const double recombination_rate_per_bacterium_per_step,
								const double between_pops_recombination_rate_per_bacterium_per_step,
								const double selection_intensity, 
								const size_t selective_allele,
								const int total_num_of_bacteria,
								const int expected_num_of_mutations_per_locus, 
								const int expectedNumMutationsPerLocusPlus5Percent,
								const string sampled_sequences_out_prefix){
	cout << "Parameters values are:" << endl;
	cout << "pop_size = " << pop_size << endl;
	cout << "num_of_pops = " << num_of_pops << endl;
	cout << "max_num_of_steps = " << max_num_of_steps << endl;
	cout << "super_step = " << super_step << endl;
	cout << "locus_len = " << locus_len << endl;
	cout << "p_a = " << p_a << endl;
	cout << "p_c = " << p_c << endl;
	cout << "p_g = " << p_g << endl;
	cout << "p_transition = " << p_transition << endl;
	cout << "selective_allele = " << selective_allele << endl;
	cout << "selection_intensity = " << selection_intensity << endl;
	cout << "mutation_rate_per_bacterium_per_step = " << mutation_rate_per_bacterium_per_step << endl;
	cout << "recombination_rate_per_bacterium_per_step = " << recombination_rate_per_bacterium_per_step << endl;
	cout << "between_pops_recombination_rate_per_bacterium_per_step = " << between_pops_recombination_rate_per_bacterium_per_step << endl;
	cout << "total_num_of_bacteria = " << total_num_of_bacteria << endl;
	cout << "expected_num_of_mutations_per_locus = " << expected_num_of_mutations_per_locus << endl;
	cout << "expectedNumMutationsPerLocusPlus5Percent = " << expectedNumMutationsPerLocusPlus5Percent << endl;
	cout << "sampled_sequences_out_prefix = " << sampled_sequences_out_prefix << endl;
}

void write_parameters_values(	string out_path,
								const int pop_size, 
								const int num_of_pops,
								const long max_num_of_steps,
								const long super_step, 
								const size_t locus_len,
								const double p_a,
								const double p_c,
								const double p_g, 
								const double p_transition, 
								const double mutation_rate_per_bacterium_per_step,
								const double recombination_rate_per_bacterium_per_step,
								const double between_pops_recombination_rate_per_bacterium_per_step, 
								const double selection_intensity,
								const size_t selective_allele, 
								const int total_num_of_bacteria, 
								const int expected_num_of_mutations_per_locus,
								const int expectedNumMutationsPerLocusPlus5Percent, 
								const string sampled_sequences_out_prefix){
	ostringstream out_stream_string;

	out_path = out_path + "parameters";

	cout << "Writing parameters configuration to:" <<endl << out_path << endl;
	
	out_stream_string  << pop_size << "\t#pop_size"<< endl;
	out_stream_string << num_of_pops << "\t#num_of_pops" << endl;
	out_stream_string << max_num_of_steps << "\t#max_num_of_steps" << endl;
	out_stream_string << super_step << "\t#super_step" << endl;
	out_stream_string << locus_len << "\t#locus_len" << endl;
	out_stream_string << p_a << "\t#p_a" << endl;
	out_stream_string << p_c << "\t#p_c" << endl;
	out_stream_string << p_g << "\t#p_g" << endl;
	out_stream_string << p_transition << "\t#p_transition" << endl;
	out_stream_string << selective_allele << "\t#selective_allele" << endl;
	out_stream_string << selection_intensity << "\t#selection_intensity" << endl;
	out_stream_string << mutation_rate_per_bacterium_per_step << "\t#mutation_rate_per_bacterium_per_step" << endl;
	out_stream_string << recombination_rate_per_bacterium_per_step << "\t#recombination_rate_per_bacterium_per_step" << endl;
	out_stream_string << between_pops_recombination_rate_per_bacterium_per_step << "\t#between_pops_recombination_rate_per_bacterium_per_step" << endl;
	out_stream_string << total_num_of_bacteria << "\t#total_num_of_bacteria" << endl;
	out_stream_string << expected_num_of_mutations_per_locus << "\t#expected_num_of_mutations_per_locus" << endl;
	out_stream_string << expectedNumMutationsPerLocusPlus5Percent << "\t#expectedNumMutationsPerLocusPlus5Percent" << endl;
	out_stream_string << sampled_sequences_out_prefix << "\t#sampled_sequences_out_prefix" << endl;

	write_ostringstream_to_file(out_stream_string, out_path);
}
