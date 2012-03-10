/* 

  GA.h : header file for the GA library.
  Copyright (C) 2010 Nuno Alves

*/

#ifndef GA_h
#define GA_h

//#define        GPOPSIZE 100


class GA
{
	public:
		GA(int _popsize, int _ngenerations, int _bitmutation, int _numbitsgene);
		void print_all_tx_population(bool populationarray);
		void process_generation();
		unsigned int return_individual(int ind, int geneb);
		bool bitRead(unsigned int groupOfBits, int bitLocation);
		void write_t0_fitness(int idx, float value) { assert(value<65000); t0_fitness[idx]=value; }
		void reportFitnessStats(int verbosityLevel);
		void returnMemory() { delete []t0_population; delete []t1_population ; delete []t0_fitness; }
		void prepare_next_generation();
		int return_id_best_individual();
		
	private:
		void mate_t1();
		void randomize_t0_population();
		unsigned int return_tx_population(bool populationarray, int individual, int geneblock);
		unsigned int randomRange(int startVal, int endVal);
		unsigned int bitWrite(unsigned int groupOfBits, int bitLocation,bool bitValue);
		void evaluate_t0_fitness();
		void expected_count_t1();
		void populate_t1();

		int bitmutation; 
		int popsize;
		int ngenerations;
		int numbitsgene;

		unsigned int *t0_population;
		unsigned int *t1_population;
		float *t0_fitness;
		float *t0_cdf;

		float top_fitness_val;
		float sum_fitness_val;
		int most_fit_individual;
};

/*
class GA
{
	public:
		GA(int _popsize, int _ngenerations, int _bitmutation);
		void process_generation(int generationid);
		void reportStatistics(int generationid, int verbosityLevel);
		void prepare_next_generation();

		void write_t0_fitness(int idx, unsigned int value) { t0_fitness[idx]=value; }
		unsigned int read_t0_a_population(int idx) { return(t0_a_population[idx]); }
		unsigned int read_t0_b_population(int idx) { return(t0_b_population[idx]); }

		bool bitRead(unsigned int groupOfBits, int bitLocation);

	private:
		unsigned int t0_fitness[GPOPSIZE];
		unsigned int top_fitness_val;
		unsigned int sum_fitness_val;
		unsigned int t1_a_population[GPOPSIZE];
		unsigned int t1_b_population[GPOPSIZE];
		unsigned int bestcandidate_a;
		unsigned int bestcandidate_b;
		unsigned int t0_a_population[GPOPSIZE];
		unsigned int t0_b_population[GPOPSIZE];
		unsigned int next_gen_expected_count[GPOPSIZE];
		
		unsigned int randomRange(int startVal, int endVal);
		unsigned int bitWrite(unsigned int groupOfBits, int bitLocation,  bool bitValue);
		void randomize_t0_population();
		void evaluate_t0_fitness();
		void expected_count_t1();
		void populate_t1();
		void mate_t1();

		int bitmutation; //should be within the range [0..1000]
		int popsize;
		int ngenerations;
		


};
*/
#endif
