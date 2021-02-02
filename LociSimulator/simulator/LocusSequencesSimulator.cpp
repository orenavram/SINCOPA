#include "LocusSequencesSimulator.h"

LocusSequenceSimulator::LocusSequenceSimulator(LocusTreeStructure & locus_tree_structure, size_t locus_len) : _locus_tree_structure(locus_tree_structure), _locus_len(locus_len){
}

void LocusSequenceSimulator::simulate(const double p_a, const double p_c, const double p_g, const double p_transition){
		
	cout << "Simulating..."<<endl;
	
	int lca_allele = _locus_tree_structure.get_oldest_allele();

	cout << "Least common ancestor is " << lca_allele <<endl;
	_sequences[lca_allele] = draw_root_sequence(p_a, p_c, p_g);
	
	if(DEBUG_MODE){cout << "with sequence: " << endl << _sequences[lca_allele] << endl;}

	if (_sequences[lca_allele].size() != _locus_len){
		cout << endl << "error!!!" << endl;
		cout << "_sequences[lca_allele] is " << _sequences[lca_allele] << endl;
		cout << "_sequences[lca_allele].size() is " << _sequences[lca_allele].size() << endl;
		cout << "_locus_len is " << _locus_len << endl;
	}

	vector<int> queue;
	queue.push_back(lca_allele);

	int ancestor, descendant;

	//TODO: get child and simulate sequence for it mutation as much as in mutation mapping. then remove it from the tree and mutation mapping
	while(queue.size() > 0){

		//pop first element
		lca_allele = queue[0];
		queue.erase(queue.begin());

		for(map<int, int>::iterator iterator = _locus_tree_structure._allelic_tree.begin(); iterator != _locus_tree_structure._allelic_tree.end(); iterator++) {
			descendant = iterator->first;
			ancestor = iterator->second;
			if(ancestor == lca_allele){
				if(DEBUG_MODE){cout << "Handling branch from " << ancestor << " to " << descendant << " with " << _locus_tree_structure._num_of_mutations_tree[descendant] << " mutations." <<endl;}
				
				// copy ancestor sequence and mutate it
				_sequences[descendant] = _sequences[ancestor];
				if (_sequences[descendant].size() != _locus_len){
					cout << "_sequences[" << descendant << "].size() != _locus_len" << endl;
					cout << _sequences[descendant] << endl;
					cout << _sequences[descendant].size() << " != " << _locus_len << endl;
				}
				
				
				if(DEBUG_MODE){cout << descendant << " was mutated from this: " << _sequences[descendant];}

				mutate_sequence(descendant, p_transition);
				
				if(DEBUG_MODE){cout << " to:" << _sequences[descendant] << endl;}


				if (_sequences[descendant].size() != _locus_len){
					cout << endl << "error!!!" << endl;
					cout << "_sequences[descendant] is " << _sequences[descendant] << endl;
					cout << "_sequences[descendant].size() is " << _sequences[descendant].size() << endl;
					cout << "_locus_len is " << _locus_len << endl;
				}
				if(DEBUG_MODE){cout << "Allele "<< descendant <<" was generated." <<endl;}
				
				queue.push_back(descendant); // enter to queue

				// _locus_tree_structure._allelic_tree.erase(descendant); // finished handling descendant sequence
				// ITERATOR SHOULD NOT BE PROMOTED. hold the iterator in place since the current element was removed from the tree map
			}
		}
	}

	for(map<size_t, string>::iterator iterator = _sequences.begin(); iterator != _sequences.end(); iterator++) {
		if(iterator->second.size() != _locus_len){
			cout << endl << "error!!!" << endl;
			cout << "iterator->second is " << iterator->second << endl;
			cout << "iterator->second.size() is " << iterator->second.size() << endl;
			cout << "_locus_len is " << _locus_len << endl;
		}
	}
	
	cout << "All mutations were handled."<<endl;
}

void LocusSequenceSimulator::mutate_sequence(const size_t allele, const double p_transition){

	int mutations_left = _locus_tree_structure._num_of_mutations_tree[allele];

	while(mutations_left){
		size_t random_nucleotide_index = rand() % _locus_len; //nucleotide to mutate
		mutate_nucleotide(random_nucleotide_index, allele, p_transition);
		mutations_left--;
	}
	
	_locus_tree_structure._num_of_mutations_tree.erase(allele); // no more mutations to handle -> no need to save the mapping anymore
}

void LocusSequenceSimulator::mutate_nucleotide(const size_t random_nucleotide_index, const int allele, const double p_transition){
	
	char new_nucleotide = '%';
	double random_event = rand()*1.0/RAND_MAX; //determine if transition or transversion
	//cout << "Mutating sequence " << child_index <<endl;
	if (random_event < p_transition){
		new_nucleotide = apply_transition(_sequences[allele][random_nucleotide_index]);
	}
	else{
		new_nucleotide = apply_transversion(_sequences[allele][random_nucleotide_index]);
	}
	
	_sequences[allele][random_nucleotide_index] = new_nucleotide; // set the new nucleotide
}

string LocusSequenceSimulator::draw_root_sequence(const double p_a, const double p_c, const double p_g) const {
	double random_event;
	string result(_locus_len, '#');
	for(size_t i=0; i<result.size(); i++){
		random_event = rand()*1.0/RAND_MAX;
		if(random_event < p_a){
			result[i] = 'a';
		}
		else if(random_event < p_a+p_c){
			result[i] = 'c';
		}
		else if(random_event < p_a+p_c+p_g){
			result[i] = 'g';
		}
		else{
			result[i] = 't';
		}
	}
	return result;
}

char LocusSequenceSimulator::apply_transition(const char nuc) const{
	switch(nuc){
		case 'a': 
			return 'g';
		case 'g': 
			return 'a';
		case 'c': 
			return 't';
		case 't': 
			return 'c';
		default:
			return '\0';
	}
}

char LocusSequenceSimulator::apply_transversion(const char nuc) const{
	double random_event = rand()*1.0/RAND_MAX;
	switch(nuc){
		case 'a':
		case 'g': 
			if(random_event < 0.5){
				return 'c';
			}
			else{
				return 't';
			}
		case 'c':
		case 't': 
			if(random_event < 0.5){
				return 'a';
			}
			else{
				return 'g';
			}
		default:
			return '\0';
	}
}



void LocusSequenceSimulator::check_data_validity3(void){
	for(map<size_t, string>::iterator iterator = _sequences.begin(); iterator != _sequences.end(); iterator++) {
		if(iterator->second.size() != _locus_len){
			cout << endl << "error!!!" << endl;
			cout << "iterator->second is " << iterator->second << endl;
			cout << "iterator->second.size() is " << iterator->second.size() << endl;
			cout << "_locus_len is " << _locus_len << endl;
		}
	} 
}
