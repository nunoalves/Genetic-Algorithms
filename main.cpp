//rm -f GAexec GA.o  ; g++ -c GA.cpp ; g++ main.cpp GA.o -o GAexec ; ./GAexec

#include <cstdlib> 
#include <ctime> 
#include <iostream>
#include <stdio.h>
#include <assert.h>


#include "GA.h"

using namespace std;

int main()
{

	//population must be less than 65k
	int populationsize=500;   
	
	//number of generations cannot exceed 65k
	int numgenerations=50;     

	//the unit of bit mutation is per thousand. For example, if we set the
	//bitmutation to 10, approximately every 10 out of 10k bits will be mutated
	//bitmutation cannot exceed 10k
	int bitmutation=1;       

	//number of bits in the gene (in multiples of 32 please)
	//num_bits_gene cannot exceed 32k
	int num_bits_gene=64;
	
	
	//instanitating the class. this will initialize a set of random genes.
	GA ga(populationsize,numgenerations,bitmutation,num_bits_gene);

	//we start in the generation 0
	int generation=0;

	//printing all individuals
	//ga.print_all_tx_population(0);

	while (generation<numgenerations)
	{ 
       		cout<<"generation #"<<generation<<endl;
	       	//======================================================================
	       	//Fitness Evaluation
	       	//======================================================================
	       	//This is the part where we "how good" a particular combination of
	       	//genes is. In this simple case, we just want to maximize the number of 1's 
	       	//in the genes. That is to say, if a particular individual has more
		//1's than another, he will have a larger fitness.
	       	//======================================================================
	      	//begin: Fitness Evaluation

		//cycle through all individuals		
	       for (int individual=0 ; individual < populationsize ; individual++)
	       {	
			float individual_fitness=0;
		
			//cycle through all geneblocks		
			for (int j=0 ; j < num_bits_gene/32 ; j++)
			{	
				unsigned int geneblock=ga.return_individual(individual,j);
			
				//count the number of 1's
		   		for (int k=0 ; k < 32 ; k++)
				{
				  if (ga.bitRead(geneblock, k)) individual_fitness++;
				}

				//cout	<<"("<<individual<<","<<j<<") "
				//	<<ga.return_individual(individual,j)
				//	<<" = "<<individual_fitness
				//	<<endl;
			}

			individual_fitness=individual_fitness/(num_bits_gene/32);
	 		//cout<<"individual "<<individual<<" = "<<individual_fitness<<endl;
			ga.write_t0_fitness(individual,individual_fitness);
		  }	//end: Fitness Evaluation

		//the next steps perform the genetic optimization process
		//behind the scenes. It mates the best individuals and create the 
		//next generation with a new set of genes.   

		ga.process_generation();
		ga.reportFitnessStats(0);
	
		//cout<<endl;
		//ga.print_all_tx_population(0);
		//cout<<endl;
	
		ga.prepare_next_generation();
		
		//cout<<endl;
		//ga.print_all_tx_population(0);
		//cout<<endl;

		//printing the genes of the best individual
		cout<<"best genes = ";
		for (int j=0 ; j < num_bits_gene/32 ; j++)
		{	
			unsigned int geneblock;
			geneblock=ga.return_individual(ga.return_id_best_individual(),j);
		
			cout<<geneblock<<" ";
		}
		cout<<endl;

		generation=generation+1;
	}

	ga.returnMemory();
}//end of main()


