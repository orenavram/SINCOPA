#include "Auxilaries.h"

int write_ostringstream_to_file(ostringstream & out_stream, const string out_path){
	ofstream myFile;
	myFile.open(out_path);
	myFile << out_stream.str();
	myFile.close();
	out_stream.str(""); //clean stream
	cout << "Done."<<endl;
	return 1;
}


int read_control_file(const string in_path, int& pop_size, int& num_of_pops, long& max_num_of_steps, long& super_step, size_t& locus_len, double& p_a, double& p_c, double& p_g, double& p_transition, double& mutation_rate_per_bacterium_per_step, double& recombination_rate_per_bacterium_per_step, double& between_pops_recombination_rate_per_bacterium_per_step, double& selection_intensity, int& selective_allele, string& sampled_sequences_out_prefix, string& done_path){
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
	max_num_of_steps = stoi(line.substr(0, line.find("#")));

	getline(myfile, line);
	super_step = stoi(line.substr(0, line.find("#")));

	getline(myfile, line);
	locus_len = stoi(line.substr(0, line.find("#")));
	getline(myfile, line);
	selective_allele = stoi(line.substr(0, line.find("#")));

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
	int nucleotides_per_bacterium = locus_len;
	double mutation_rate_per_nucleotide_per_step = stod(line.substr(0, line.find("#")));
	mutation_rate_per_bacterium_per_step = mutation_rate_per_nucleotide_per_step * nucleotides_per_bacterium;

	getline(myfile, line);
	double r_to_m_rate = stod(line.substr(0, line.find("#")));
	recombination_rate_per_bacterium_per_step = r_to_m_rate * mutation_rate_per_bacterium_per_step;

	getline(myfile, line);
	double between_pops_recombination_rate_to_recombination_rate_per_bacterium_per_step = stod(line.substr(0, line.find("#")));
	between_pops_recombination_rate_per_bacterium_per_step = between_pops_recombination_rate_to_recombination_rate_per_bacterium_per_step * recombination_rate_per_bacterium_per_step;

	getline(myfile, line);
	cout << line.substr(0, line.find("#")) << endl;
	sampled_sequences_out_prefix = remove_trailing_whitespaces(line.substr(0, line.find("#")));

	getline(myfile, line);
	cout << line.substr(0, line.find("#")) << endl;
	done_path = remove_trailing_whitespaces(line.substr(0, line.find("#")));
	
	cout << "Closing control file." << endl;
	myfile.close();
	cout << "Done!" << endl;
	
	return 1;
}

string remove_trailing_whitespaces(string str){
	int i = str.size();
	while (i>0){
		i--;
		if(str[i] != ' ' && str[i] != '\t' && str[i] != '\r' && str[i] != '\n'){
			break;
		}
	}
	return str.substr(0, i+1);
}