#pragma once
#include "SuperPopulation.h"


class WrightFisherProcess {
public:

	explicit WrightFisherProcess(SuperPopulation & old_pops,
		long super_step, // how many steps between a super step event
		long max_num_of_steps,
		double mutation_rate_per_bacterium_per_step,
		double recombination_rate_per_step,
		double between_pops_recombination_rate_per_bacterium_per_step);
	void runWF();
	void doSuperStep();
	void wright_fisher_step(size_t step);
	void mutate_bacterium(size_t popIndex, size_t bactIndex);
		
	void recombine_bacterium(const size_t current_bacterium_index,const size_t popIndex);
	void recombination_between_populations(const size_t target_pop_index,
		const size_t target_bacterium_index);

	long _super_step; // how many steps between a super step event
	long _max_num_of_steps;
	double _mutation_rate_per_bacterium_per_step;
	double _recombination_rate_per_step;
	double _between_pops_recombination_rate_per_bacterium_per_step;
	vector<int> _mutations_history; // to be populated with the mutation history.
	SuperPopulation *_old_pops;
	SuperPopulation *_new_pops;
	SuperPopulation _SuperPopulation1; // these will hold the actual data
	SuperPopulation _SuperPopulation2; // // these will hold the actual data
	
	
	

};