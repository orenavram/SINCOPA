#include "LocusTreeStructure.h"

LocusTreeStructure::LocusTreeStructure(const vector<int> & mutations_history, vector<int> current_alleles, size_t total_num_of_bacteria, size_t pop_size){
	
	//current_alleles will maintain the alleles state at every step when going back in time

	int max_allele, max_allele_index, parent_allele;
	
	//for debugging:
	//int last_last_parent_allele=-1, last_last_max_allele=-1, last_parent_allele=-1, last_max_allele=-1;
	
	size_t prev_size;

	//remove duplicates
	sort(current_alleles.begin(), current_alleles.end());
	current_alleles.erase(unique(current_alleles.begin(), current_alleles.end() ), current_alleles.end() );

	vector<int> last_coalecsent = current_alleles; // maintain the alleles state that was in last coalecsent when going back in time

	vector<int> mutation_counters(current_alleles.size(), 0); // maintain how many mutations have happend between the last coalecsent until now when going back in time

	while(there_is_more_than_one_allele(last_coalecsent)){

		if(PRINT_MODE){cout << "last_coalecsent.size() = " << last_coalecsent.size() << endl;}
		max_allele = *max_element(current_alleles.begin(), current_alleles.end());
		max_allele_index = find(current_alleles.begin(), current_alleles.end(), max_allele) - current_alleles.begin();
		//current_allele = current_alleles[max_allele_index];
		parent_allele = mutations_history[max_allele];

		if(parent_allele == max_allele){
				cout << "max_allele == current_allele  !!!!"<< endl;
		}

		current_alleles[max_allele_index] = parent_allele; // change the allele to its parent
		mutation_counters[max_allele_index]++; // add one mutation to the counting

		//check if there is a coalecsent event and if so, update counters and tree structure
		int num_of_occurences = count(current_alleles.begin(), current_alleles.end(), parent_allele);
		if(num_of_occurences > 1){
			prev_size = last_coalecsent.size();
			current_alleles = coalesce(current_alleles, last_coalecsent, mutation_counters, parent_allele, num_of_occurences);
			if(last_coalecsent.size() == prev_size){
				cout << "size didnt change!!!!"<< endl;
			}
		}
		//if(DEBUG_MODE){
		//	last_last_parent_allele = last_parent_allele;
		//	last_last_max_allele = last_max_allele;
		//	last_parent_allele = parent_allele;
		//	last_max_allele = max_allele;
		//}
	}
}

vector<int> LocusTreeStructure::coalesce(vector<int> current_alleles, vector<int> & last_coalecsent, vector<int> & mutation_counters, int lca_allele, int num_of_occurences){
	int child, i = current_alleles.size()-1;
	int flag=0;
	while(i>=0){
		if(current_alleles[i] == lca_allele){

			child = last_coalecsent[i]; //get child of the coalescencing allele
			
			if (mutation_counters[i] != 0){ // avoid overrding with 0 in case of multiforcation
				if(_allelic_tree.find(child) != _allelic_tree.end()){
					cout << "Overriding child " << child << " !! already exists!!" << endl;
					cout << "_allelic_tree["<< child <<"] = " << _allelic_tree[child] << endl;
					cout << "_num_of_mutations_tree["<< child <<"] = " << _num_of_mutations_tree[child] << endl;
					flag=1;
				}
				_allelic_tree[child] = lca_allele; //map the child to its least common ansector
				_num_of_mutations_tree[child] = mutation_counters[i]; //track number of mutation between the child and its lca

				if(flag){
					cout << "Changed to:" << endl;
					cout << "_allelic_tree[child] = " << _allelic_tree[child] << endl;
					cout << "_num_of_mutations_tree[child] = " << _num_of_mutations_tree[child] << endl;
					flag=0;
				}

				//next 2 lines are relevant when handling the last occurence of parrent allele.
				//in all other cases they erased anyway.
				last_coalecsent[i] = lca_allele; // update the coalescent vector
				mutation_counters[i] = 0; // initialize counter
			}

			if(num_of_occurences == 1){ // don't remove last individuals data!!
				break;
			}

			// remove coalescenced individual data from lists
			current_alleles.erase(current_alleles.begin() + i);
			last_coalecsent.erase(last_coalecsent.begin() + i);
			mutation_counters.erase(mutation_counters.begin() + i);
			
			num_of_occurences--;
		}
		i--;
	}

	return current_alleles;
}

bool LocusTreeStructure::there_is_more_than_one_allele(vector<int> & last_coalecsent){
	return last_coalecsent.size() > 1;
}

int LocusTreeStructure::get_oldest_allele(){
 	int min_val = _allelic_tree.begin()->second;
	for(map<int, int>::iterator iterator = _allelic_tree.begin(); iterator != _allelic_tree.end(); iterator++) {
		if(min_val > iterator->second){
			min_val = iterator->second;
		}
	}
	return min_val;
}
