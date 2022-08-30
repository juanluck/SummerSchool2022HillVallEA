/*

HillVallEA 

Real-valued Multi-Modal Evolutionary Optimization

By S.C. Maree
s.c.maree[at]amc.uva.nl
github.com/SCMaree/HillVallEA

Example script to demonstrate the usage of HillVallEA
 on the well-known 2D Six Hump Camel Back function

*/

#include "HillVallEA/hillvallea.hpp"
#include "HillVallEA/fitness.h"

namespace hillvallea
{
class example_problem_t : public fitness_t
{
	public:

	example_problem_t()
	{
	  number_of_parameters = 2;
	  maximum_number_of_evaluations = 10000;
	}
	~example_problem_t() {}

	void get_param_bounds(vec_t & lower, vec_t & upper) const
	{
	  lower.resize(number_of_parameters, 0);
	  upper.resize(number_of_parameters, 0);
	  
	  lower[0] = -3.0;
	  lower[1] = -2.0;
	  upper[0] = 3.0;
	  upper[1] = 2.0;
	}

	void define_problem_evaluation(solution_t & sol)
	{
	  double p0s = sol.param[0]*sol.param[0]; // param 0 squared
	  double p1s = sol.param[1]*sol.param[1]; // param 1 squared
	  
	  sol.f = (4.0-2.1*p0s + p0s*p0s/3.0) * p0s + sol.param[0]*sol.param[1] + (-4.0 + 4.0*p1s)*p1s;
	  sol.penalty = 0.0;
	}


	std::string name() const { return "SixHumpCamelBack"; }
};
}


double fitness_aggregate(std::vector<hillvallea::solution_pt> elitist_archive){
	
  double sum             = 1.0;
  double square_of_peaks =  elitist_archive.size() * elitist_archive.size();
  
  if ( square_of_peaks == 0) // we cannot divide by 0
	  return 9999999;

  for(size_t i = 0; i < elitist_archive.size(); ++i)
  {
     sum += -1 * elitist_archive[i]->f; // -1 is to change from a maximization to a minimization problem
  }
  
  return sum / square_of_peaks;
	
}

	
// Main: Run the CEC2013 niching benchmark
//--------------------------------------------------------
int main(int argc, char **argv)
{
  
	
	//std::cout << "You have entered " << argc
    //     << " arguments:" << "\n";
	double popinit   = std::stod(argv[2]);
	double popinc    = std::stod(argv[4]);
    double clustinit = std::stod(argv[6]);
    double clustinc  = std::stod(argv[8]);
    int    maxelites = std::stoi(argv[10]);
        
	
	
  // Problem definition
  // Note: define as minimization problem!
  //-----------------------------------------
  hillvallea::fitness_pt fitness_function = std::make_shared<hillvallea::example_problem_t>();
  hillvallea::vec_t lower_range_bounds, upper_range_bounds;
  fitness_function->get_param_bounds(lower_range_bounds, upper_range_bounds);
  
  // HillVallEA Settings
  //-----------------------------------------
  // Type of local optimizer to be used.
  // 0 = AMaLGaM, 1 = AMaLGaM-Univariate, 20 = iAMaLGaM, 21 = iAMaLGaM-Univariate
  size_t local_optimizer_index = 1; // AMaLGaM-Univariate (1) is suggested
  
  int maximum_number_of_evaluations = 10000; // maximum number of evaluations
  int maximum_number_of_seconds = 3600; // maximum runtime in seconds
  
  // if the optimum is known, you can terminate HillVallEA if it found a solution
  // with fitness below the value_to_reach (vtr)
  double value_to_reach = 0;
  bool use_vtr = false;
  
  // random seed initialization for reproducibility
  srand(time(NULL));
  int random_seed = rand();
  
  // Output to test files
  bool write_generational_solutions = false;
  bool write_generational_statistics = true;
  std::string write_directory = "./";
  std::string file_appendix = ""; // can be used when multiple runs are outputted in the same directory
  
  // Initialization of HillVallEA
  //-----------------------------------------
  hillvallea::hillvallea_t opt(
     fitness_function,
     (int) local_optimizer_index,
     (int) fitness_function->number_of_parameters,
     lower_range_bounds,
     upper_range_bounds,
     lower_range_bounds,
     upper_range_bounds,
     maximum_number_of_evaluations,
     maximum_number_of_seconds,
     value_to_reach,
     use_vtr,
     random_seed,
     write_generational_solutions,
     write_generational_statistics,
     write_directory,
     file_appendix
  );
  

  opt.set_params(popinit,popinc,clustinit,clustinc,maxelites);

  // Running HillVallEA
  std::cout << "Running HillVallEA on the Six Hump Camel back function" << std::endl;
  
  
  opt.run();
  
  std::cout << "HillVallEA finished" << std::endl;
  std::cout << "Generation statistics written to " << write_directory << "statistics" << file_appendix << ".dat" << std::endl;
  std::cout << "Elitist archive written to       " << write_directory << "elites" << file_appendix << ".dat" << std::endl;
  
  std::cout << "HillVallEA Obtained " << opt.elitist_archive.size() << " elites: " << std::endl;
  
  std::cout << "    Fitness      Penalty   Params" << std::endl;
  for(size_t i = 0; i < opt.elitist_archive.size(); ++i)
  {
    std::cout << std::setw(11) << std::scientific << std::setprecision(3) << opt.elitist_archive[i]->f << "  ";
    std::cout << std::setw(11) << std::scientific << std::setprecision(3) << opt.elitist_archive[i]->penalty << "  ";
    std::cout << std::setw(11) << std::scientific << std::setprecision(3) << opt.elitist_archive[i]->param << std::endl;
  }
  
  
  std::cout << std::endl << "Aggregated fitness -> ";
  
  std::cout << std::setw(11) << std::scientific << std::setprecision(3) << fitness_aggregate(opt.elitist_archive) << std::endl << std::endl;

  return(0);
}
