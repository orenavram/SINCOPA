#include "Auxilaries.h"


double draw_from_zero_to_one_in_high_resolution_steps(int fraction)
{
	// draw a random number from 0 to 1 with steps of RAND_MAX^-fraction
	// if x ~ U[0,1] with steps of 1/k, then, x_1/k^1 + x_2/k^2 + ... + x_i/k^i ~ U[0,1] with steps of 1/k^i
	double result = 0;
	for (size_t i = 1; i <= fraction; i++) {
		result += rand()*1.0 / pow(RAND_MAX, i);
	}
	return result;
}


int write_ostringstream_to_file(ostringstream & out_stream, const string out_path){
	ofstream myFile;
	myFile.open(out_path);
	myFile << out_stream.str();
	myFile.close();
	out_stream.str(""); //clean stream
	cout << "Done."<<endl;
	return 1;
}


int read_control_file(const string in_path, int& pop_size, int& num_of_pops, long& number_of_generation, long& super_generation, size_t& num_of_loci, size_t& locus_len, size_t& num_of_non_recombining_loci, double& p_a, double& p_c, double& p_g, double& p_transition, double& mutation_rate_per_bacterium_per_generation, double& recombination_rate_per_bacterium_per_generation, double& between_pops_recombination_rate_per_bacterium_per_generation, double& selection_intensity, double& selective_allele_timing_proportion, string& output_dir){
    string line;
	ifstream myfile;	
	
	cout << "Openning control file from:" << endl << in_path << endl;
	myfile.open(in_path);
	cout << "Reading control file." << endl;
	
	getline(myfile, line);
	pop_size = stoi(line.substr(0, line.find("#")));

	getline(myfile, line);
	num_of_pops = stoi(line.substr(0, line.find("#")));

	getline(myfile, line);
	number_of_generation = stoi(line.substr(0, line.find("#")));

	getline(myfile, line);
	super_generation = stoi(line.substr(0, line.find("#")));

	getline(myfile, line);
	num_of_loci = stoi(line.substr(0, line.find("#")));

	getline(myfile, line);
	locus_len = stoi(line.substr(0, line.find("#")));

	getline(myfile, line);
	num_of_non_recombining_loci = stoi(line.substr(0, line.find("#")));

	getline(myfile, line);
	selective_allele_timing_proportion = stod(line.substr(0, line.find("#")));

	getline(myfile, line);
	selection_intensity = stod(line.substr(0, line.find("#")));

	getline(myfile, line);
	p_a = stod(line.substr(0, line.find("#")));
	
	getline(myfile, line);
	p_c = stod(line.substr(0, line.find("#")));

	getline(myfile, line);
	p_g = stod(line.substr(0, line.find("#")));

	getline(myfile, line);
	p_transition = stod(line.substr(0, line.find("#")));

	getline(myfile, line);
	size_t nucleotides_per_bacterium = locus_len * num_of_loci;
	double mutation_rate_per_nucleotide_per_generation = stod(line.substr(0, line.find("#")));
	mutation_rate_per_bacterium_per_generation = mutation_rate_per_nucleotide_per_generation * nucleotides_per_bacterium;
	if (mutation_rate_per_bacterium_per_generation > 1) {
		cout << "\nmutation_rate_per_bacterium_per_generation is " << mutation_rate_per_bacterium_per_generation << " (> 1). Setting it to 1.\n" << endl;
		mutation_rate_per_bacterium_per_generation = 1;
	}

	getline(myfile, line);
	double r_to_m_rate = stod(line.substr(0, line.find("#")));
	recombination_rate_per_bacterium_per_generation = r_to_m_rate * mutation_rate_per_bacterium_per_generation;

	getline(myfile, line);
	double between_pops_recombination_rate_to_recombination_rate_per_bacterium_per_generation = stod(line.substr(0, line.find("#")));
	between_pops_recombination_rate_per_bacterium_per_generation = between_pops_recombination_rate_to_recombination_rate_per_bacterium_per_generation * recombination_rate_per_bacterium_per_generation;

	getline(myfile, line);
	cout << line.substr(0, line.find("#")) << endl;
	output_dir = remove_trailing_whitespaces(line.substr(0, line.find("#")));

	cout << "Closing control file." << endl;
	myfile.close();
	cout << "Done!" << endl;
	
	return 1;
}

string remove_trailing_whitespaces(string str){
	size_t i = str.size();
	while (i>0){
		i--;
		if(str[i] != ' ' && str[i] != '\t' && str[i] != '\r' && str[i] != '\n'){
			break;
		}
	}
	return str.substr(0, i+1);
}