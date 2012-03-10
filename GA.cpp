#include <cstdlib> 
#include <ctime> 
#include <iostream>
#include <stdio.h>
#include <assert.h>

#include "GA.h"

using namespace std;


GA::GA(int _popsize, int _ngenerations, int _bitmutation, int _numbitsgene)
{
	srand((unsigned)time(0)); 

  	popsize 	= _popsize;
	ngenerations  	= _ngenerations;
	bitmutation	= _bitmutation;
	numbitsgene	= _numbitsgene;
	
	assert (popsize 	< 65535);
	assert (ngenerations 	< 65535);
	assert (bitmutation 	< 10000);
	assert (numbitsgene 	< 32000);
	assert (numbitsgene%32 	== 0);


	t0_population  = new unsigned int [ int(popsize*(numbitsgene / 32)) ];
	t1_population  = new unsigned int [ int(popsize*(numbitsgene / 32)) ];
	t0_fitness     = new float [ int(popsize) ];
	t0_cdf		   = new float [ int(popsize) ];

	//once i initialize things i can also initialize the random population
	randomize_t0_population();

	top_fitness_val		= 0;
	sum_fitness_val		= 0;
	most_fit_individual	= 0;
}

//randomize all genes of each individual in the t0_population
void GA::randomize_t0_population()
{
	for (int i=0 ; i < popsize ; i++)
  	{

  		for (int j=0 ; j < (numbitsgene / 32) ; j++)
 	 	{
		   	unsigned int randNumberA = random(); 
			t0_population[j+(i*(numbitsgene / 32))]=randNumberA;
		}
	
	}
}

int GA::return_id_best_individual()
{
	return(most_fit_individual);
}


unsigned int GA::return_individual(int ind, int geneb)
{
	return_tx_population(0, ind, geneb);
}

//polulationarray = 0 is the array of the current generation
//polulationarray = 1 is the array of the next generation
unsigned int GA::return_tx_population(bool polulationarray, int individual, int geneblock)
{
	if (polulationarray==0)
	return (t0_population[geneblock+(individual*(numbitsgene / 32))]);

	if (polulationarray==1)
	return (t1_population[geneblock+(individual*(numbitsgene / 32))]);

}

unsigned int GA::bitWrite(unsigned int groupOfBits, int bitLocation,bool bitValue)
{
	if (bitValue==1) groupOfBits |= 1 << bitLocation;
	else groupOfBits &= ~(1 << bitLocation); 	
	return(groupOfBits);
	
}

bool GA::bitRead(unsigned int groupOfBits, int bitLocation)
{
	bool bit = groupOfBits & (1 << bitLocation);
	return(bit);
}

//if startVal==1 and endVal==10 randomNumber[1..10]
unsigned int GA::randomRange(int startVal, int endVal)
{
	unsigned int random_integer; 
	
	random_integer = (rand()%endVal)+startVal; 
}


void GA::print_all_tx_population(bool polulationarray)
{
	for (int i=0 ; i < popsize ; i++)
	{
		cout<<"individual "<<i<<": ";
		for (int j=0 ; j < (numbitsgene / 32) ; j++)
		{	
			if (polulationarray==0)	
			cout<<return_tx_population(0,i,j)<<" ";
			if (polulationarray==1)	
			cout<<return_tx_population(1,i,j)<<" ";
		}
		cout<<endl;
	}

}

//store just the current generation fitness to save memory
void GA::evaluate_t0_fitness()
{
	top_fitness_val		= 0;
	sum_fitness_val		= 0;
	most_fit_individual	= 0;

	for (int ii=0 ; ii < popsize ; ii++)
	{
		sum_fitness_val = t0_fitness[ii]+sum_fitness_val;
		assert(sum_fitness_val<65000);
  
		if (t0_fitness[ii]>top_fitness_val) 
		{
			top_fitness_val=t0_fitness[ii];
			most_fit_individual=ii;	
		}
	}
}

//verbosityLevel=0,1,2
void GA::reportFitnessStats(int verbosityLevel)
{
	
	if (verbosityLevel>=0)
	{
    	cout	<<"top fitness="
    		<<top_fitness_val
		<<", top candidate="
		<<most_fit_individual
    		<<", sum all fitnesses="
    		<<sum_fitness_val<<endl;
	}
	
	if (verbosityLevel>=1)
	{
		//cycle through all individuals		
	       for (int ind=0 ; ind < popsize ; ind++)
	       {	
			cout<<"indivual #"<<ind<<" = ";		
			//cycle through all geneblocks		
			for (int j=0 ; j < numbitsgene/32 ; j++)
			{	
				unsigned int geneblock=return_individual(ind,j);
				cout<<geneblock<<" ";			
			}
			cout<<", " <<t0_fitness[ind];				
			cout<<", " <<t0_cdf[ind]<<endl;				
		  }
		
	}	
}

void GA::expected_count_t1()
{   
	//normalize the fitness
	for (int x=0; x < popsize; x++)
	{
		t0_cdf[x]=(t0_fitness[x]*100)/sum_fitness_val;
		if (x!=0) t0_cdf[x]=t0_cdf[x]+t0_cdf[x-1];
	}	

}

void GA::populate_t1()
{
	for (int i=0 ; i < popsize ; i++)
	{
		int randomValue=randomRange(1,100); 
		//cout<<randomValue<<endl;
		int select_index=0;
	
		for (int j=1 ; j < popsize ; j++)
		{
			if ((randomValue>=t0_cdf[j-1]) && (randomValue<t0_cdf[j]))
			{
				select_index=j; 
				j=popsize; //exit the loop
			}  
		}
	
		for (int k=0; k < numbitsgene/32 ; k++)
		{		
			t1_population[k+(i*(numbitsgene / 32))]= t0_population[k+(select_index*(numbitsgene / 32))];
		}
	
	}

}

//switch bits between two individuals at a particular cross point
void GA::mate_t1()
{
  for (int popindex=0 ; popindex < popsize ; popindex=popindex+2)
  {
		
		
		int crossoverSite=randomRange(0,numbitsgene-1);    
		int crossoverBlock=(crossoverSite+1)/32;
		int crossoverBit=(crossoverSite+1)%32;
		/*
		cout <<"======================================="<<endl;
		cout <<"crossoverSite="<<crossoverSite<<endl;
		cout<<"block="<<crossoverBlock<<endl;
		cout<<"crossoverBit="<<crossoverBit<<endl;
	

		cout<<"before crossover"<<endl; 
		cout<<"popindex="<<popindex<<" ";
		for (int k=0; k < numbitsgene/32 ; k++)	
			cout<<t1_population[k+(popindex*(numbitsgene / 32))]<<" ";				
			cout<<endl;
		
		cout<<"popindex="<<popindex+1<<" ";
		for (int k=0; k < numbitsgene/32 ; k++)
			cout<<t1_population[k+((popindex+1)*(numbitsgene / 32))]<<" ";
			cout<<endl;
		*/

		//start exchanging the bits on the crossover block 
		unsigned int t1_i=t1_population[crossoverBlock+(popindex*(numbitsgene / 32))];
		unsigned int t1_i_plus_1=t1_population[crossoverBlock+((popindex+1)*(numbitsgene / 32))];
		
		for (int bit=0; bit<crossoverBit;bit++)
		{
			//cout<<bit<<endl;
			bool siteA=bitRead(t1_i,bit);
			bool siteB=bitRead(t1_i_plus_1,bit);
	
			t1_population[crossoverBlock+(popindex*(numbitsgene / 32))]=bitWrite(t1_population[crossoverBlock+(popindex*(numbitsgene / 32))],bit,siteB);	
			t1_population[crossoverBlock+((popindex+1)*(numbitsgene / 32))]=bitWrite(t1_population[crossoverBlock+((popindex+1)*(numbitsgene / 32))],bit,siteA);	
	
		}
		//exchange the remaining blocks
		
		for (int block=crossoverBlock+1 ; block < numbitsgene / 32 ; block++)
		{
			unsigned int t1_i=t1_population[block+(popindex*(numbitsgene / 32))];
			unsigned int t1_i_plus_1=t1_population[block+((popindex+1)*(numbitsgene / 32))];
			
			t1_population[block+(popindex*(numbitsgene / 32))]=t1_i_plus_1;
			t1_population[block+((popindex+1)*(numbitsgene / 32))]=t1_i;
		}

		/*
		cout<<"after crossover"<<endl; 
		cout<<"popindex="<<popindex<<" ";
		for (int k=0; k < numbitsgene/32 ; k++)	
			cout<<t1_population[k+(popindex*(numbitsgene / 32))]<<" ";				
			cout<<endl;
		
		cout<<"popindex="<<popindex+1<<" ";
		for (int k=0; k < numbitsgene/32 ; k++)
			cout<<t1_population[k+((popindex+1)*(numbitsgene / 32))]<<" ";
			cout<<endl;
		cout<<endl;
		*/
		
  } //for (int i=0 ; i < popsize ; i=i+2)

}


//processes all the steps for current generation
void GA::process_generation()
{
	//calculates the average and find the best individual
	evaluate_t0_fitness();

	//update t0_probability with the probability of a particular individual
	//show up on t1?
	expected_count_t1();

	//populate the next generation array with the best performing individuals
	populate_t1();
	
	//mate every two elements at a random crossover point
	mate_t1();
}

//in essence, just copies each element of t1 into t0
void GA::prepare_next_generation()
{

	for (int i=0 ; i < popsize ; i++)
	{	
		for (int k=0; k < numbitsgene/32 ; k++)
		{		
			t0_population[k+(i*(numbitsgene / 32))]= t1_population[k+(i*(numbitsgene / 32))];
		}
	
	}
  
}

