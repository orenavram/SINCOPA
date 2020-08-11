#include "WrightFisherProcess.h"



WrightFisherProcess::WrightFisherProcess(SuperPopulation & old_pops,
	long super_step, // how many steps between a super step event
	long max_num_of_steps,
	double mutation_rate_per_bacterium_per_step,
	double recombination_rate_per_step,
	double between_pops_recombination_rate_per_bacterium_per_step): _old_pops(NULL),_new_pops(NULL),_SuperPopulation1(old_pops), _SuperPopulation2(old_pops),
	_super_step(super_step),
	_max_num_of_steps(max_num_of_steps), 
	_mutation_rate_per_bacterium_per_step(mutation_rate_per_bacterium_per_step),
	_recombination_rate_per_step(recombination_rate_per_step),
	_between_pops_recombination_rate_per_bacterium_per_step(between_pops_recombination_rate_per_bacterium_per_step)
{
	const int startingMutationHistoryIndex = -1;
	_mutations_history.push_back(startingMutationHistoryIndex);
	
}

void WrightFisherProcess::runWF() {
	for (long step = 0; step < _max_num_of_steps; ++step ) {
		if (step % 2 == 0) {
			_old_pops = &_SuperPopulation1;
			_new_pops = &_SuperPopulation2;
		}
		else {
			_old_pops = &_SuperPopulation2;
			_new_pops = &_SuperPopulation1;
		}
		if (step % 1000 == 0) 
			cout << "step = " << step << endl;
		if (step%_super_step == 0) doSuperStep();
		
		
		wright_fisher_step(step);
		

	}
}

void WrightFisherProcess::wright_fisher_step(size_t step) {

	

	//cout << "printing old pops before WF step:" << endl;
	//_old_pops->print();
	//cout << "printing new pops before WF step:" << endl;
	//_new_pops->print();
	_new_pops->_num_of_different_alleles = _old_pops->_num_of_different_alleles;

	for (size_t popIndex = 0; popIndex < _old_pops->get_num_of_pops(); popIndex++) {
		double mutate, recombine, recombine_between_pops;

		(*_new_pops)[popIndex]._selective_pop_size = 0;
		(*_new_pops)[popIndex]._neutral_pop_size = 0;
		if (PRINT_MODE) _new_pops->check_superpop_validity();
		
		for (size_t bactIndex = 0; bactIndex < (*_old_pops)[popIndex].get_total_size(); bactIndex++) {
			int randomParentIndex = (*_old_pops)[popIndex].draw_bacterium_index_selectively();
			(*_new_pops)[popIndex].add_bacterium((*_old_pops)[popIndex][randomParentIndex]);
				
			if (PRINT_MODE) {
				cout << "bacteria " << bactIndex << " is the son of bacteria " << randomParentIndex << endl;
				if (static_cast <int> (_new_pops->get_num_of_different_alleles_at_locus()) < (*_new_pops)[popIndex][bactIndex].get_allele()) {
					cout << "the number of different alleles observed in the population = ";
					cout << _new_pops->get_num_of_different_alleles_at_locus() << endl;
					cout << "parental allele is ";
					cout << (*_new_pops)[popIndex][bactIndex].get_allele() << endl;
					exit(-1);
				}
			}
			
			//mutation
			mutate = rand()*1.0 / RAND_MAX;
			if (mutate < _mutation_rate_per_bacterium_per_step) {
				if (PRINT_MODE) {
					cout << "locus has BEFORE " << _new_pops->get_num_of_different_alleles_at_locus() << " different alleles" << endl;
				}
				
				mutate_bacterium(popIndex,bactIndex);
				
				_new_pops->add_new_allele_to_locus();

				if (PRINT_MODE) {
					cout << "locus has AFTER " << _new_pops->get_num_of_different_alleles_at_locus() << " different alleles" << endl;
				}
			}

			//recombination
			recombine = rand()*1.0 / RAND_MAX;
			if (recombine < _mutation_rate_per_bacterium_per_step) {
				recombine_bacterium(bactIndex, popIndex);
			}


			//between population recombination
			recombine_between_pops = rand()*1.0 / RAND_MAX;
			if (recombine_between_pops < _between_pops_recombination_rate_per_bacterium_per_step) {
				//if (step > 10000) {
				//	cout << (*_new_pops)[popIndex][bactIndex].get_allele()<<endl;
				//}
				recombination_between_populations(popIndex, bactIndex);
				//if (step > 10000) {
				//	cout << (*_new_pops)[popIndex][bactIndex].get_allele()<<endl;
				//}
			}

		


			if (PRINT_MODE) {
				_new_pops->check_superpop_validity();
				(*_new_pops)[popIndex][bactIndex].print(bactIndex);
			}

		}




		if (PRINT_MODE) {
			cout << "step = " << step << endl;
			(*_new_pops)[popIndex].print();
			(*_old_pops)[popIndex].print();
		}
	}
	//cout << "printing old pops:" << endl;
	//_old_pops->print();
	//cout << "printing new pops:" << endl;
	//_new_pops->print();
	//cout << endl;
}

void WrightFisherProcess::doSuperStep() {
/*	for (size_t i = 0; i < old_pops.get_num_of_pops(); i++) { // old_pops is more updated now. thus the super_step is applied on old_pops

		size_t randomParentIndex = new_pops.draw_pop_index_uniformly();

		old_pops.set_pop(new_pops.get_pop(randomParentIndex), i);

	}
	if (PRINT_MODE) { new_pops.check_superpop_validity(); }
	*/
}

void WrightFisherProcess::mutate_bacterium(size_t popIndex, size_t bactIndex) {
	size_t num_of_different_alleles = _new_pops->get_num_of_different_alleles_at_locus();
	size_t old_allele = (*_new_pops)[popIndex][bactIndex].get_allele();

	if (PRINT_MODE) {
		cout << "MUTATION: num_of_different_alleles is:" << num_of_different_alleles << endl;
	}

	if (num_of_different_alleles < old_allele) {
		cout << "new < old!!!!" << endl;
		exit(1);
	}

	_mutations_history.push_back(old_allele);

	(*_new_pops)[popIndex][bactIndex].set_allele(num_of_different_alleles);//set mutation
}

void WrightFisherProcess::recombine_bacterium(	const size_t current_bacterium_index,
												const size_t popIndex) {
	size_t random_other_bacterium_index = (*_old_pops)[0].draw_other_bacterium_index_uniformly(current_bacterium_index);
	int new_allele = (*_old_pops)[popIndex][random_other_bacterium_index].get_allele();

	if (PRINT_MODE) {
		int old_allele = (*_new_pops)[popIndex][current_bacterium_index].get_allele();
		cout << "RECOMBINATION: Bacterium " << current_bacterium_index << " from Bacterium " << random_other_bacterium_index << " received allele number " << new_allele << " (had " << old_allele << ")" << endl;
	}
	(*_new_pops)[popIndex][current_bacterium_index].set_allele(new_allele);//set recombination
}

void WrightFisherProcess::recombination_between_populations(const size_t target_pop_index,
															const size_t target_bacterium_index) {
	if (_old_pops->get_num_of_pops()<2) { //superpop is degenerate. No "other" pop to get allele from..
		cout << "less than two popiulations!" << endl;
		return;
	}

	//this bacterium will *recieve* the allele
	//size_t random_target_population_index = rand() % old_pops.get_num_of_pops();
	//size_t random_target_bacterium_index = rand() % old_pops.get_pop(random_target_population_index).get_total_size();
	//Bacterium target_bacterium = Bacterium(old_pops.get_pop(random_target_population_index).get_bacterium(random_target_bacterium_index)); // MUST be a copy of the target bacterium

	//this bacterium will *donate* the allele
	size_t random_source_population_index = rand() % (_old_pops->get_num_of_pops() - 1); //(excluding the source_pop)
	if (random_source_population_index >= target_pop_index) {
		random_source_population_index++;
	}
	//cout << "target population = " << target_pop_index << " source rand population= " << random_source_population_index << endl;
	size_t random_source_bacterium_index = rand() % ((*_old_pops)[random_source_population_index].get_total_size()); // index can be the same as the source bacterium since they are from different pops
	//Bacterium source_bacterium(_old_pops.get_pop(random_source_population_index).get_bacterium(random_source_bacterium_index));

	int old_allele = (*_new_pops)[target_pop_index][target_bacterium_index].get_allele();
	int new_allele = (*_old_pops)[random_source_population_index][random_source_bacterium_index].get_allele();

	if (PRINT_MODE) {
		cout << "BETWEEN POPULATION RECOMBINATION: Pop " << target_pop_index;
		cout << " Bacterium " << target_bacterium_index;
		cout << " from Pop " << random_source_population_index;
		cout << " Bacterium " << random_source_bacterium_index;
		cout << " received allele number " << new_allele;
		cout<< " (had " << old_allele << ")" << endl;
	}
	(*_new_pops)[target_pop_index][target_bacterium_index].set_allele(new_allele);//set recombination
}