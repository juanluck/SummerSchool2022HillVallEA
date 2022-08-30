/*

AMaLGaM-Univariate as part of HillVallEA

Implementation by S.C. Maree
s.c.maree[at]amc.uva.nl
github.com/SCMaree/HillVallEA


*/
#include "mathfunctions.hpp"
#include "amalgam_univariate.hpp"


// init amalgam default parameters
hillvallea::amalgam_univariate_t::amalgam_univariate_t(const size_t number_of_parameters, const vec_t & lower_param_bounds, const vec_t & upper_param_bounds, double init_univariate_bandwidth, fitness_pt fitness_function, rng_pt rng) : optimizer_t(number_of_parameters, lower_param_bounds, upper_param_bounds, init_univariate_bandwidth, fitness_function, rng)
{

  // default values for AMaLGaM
  maximum_no_improvement_stretch = (int)(number_of_parameters + 25);
  selection_fraction = 0.35;
  st_dev_ratio_threshold = 1.0; 
  distribution_multiplier_decrease = 0.9;
  sample_succes_ratio_threshold = 0.1;
  param_std_tolerance = 1e-15; //  1e-15 needed to solve weierstrass! 
  fitness_std_tolerance = 1e-12;

  apply_ams = true;
  delta_ams = 2.0;

}

hillvallea::amalgam_univariate_t::~amalgam_univariate_t(){};

hillvallea::optimizer_pt hillvallea::amalgam_univariate_t::clone() const
{

  amalgam_univariate_pt opt = std::make_shared<amalgam_univariate_t>(number_of_parameters, lower_param_bounds, upper_param_bounds, init_univariate_bandwidth, fitness_function, rng);

  // Optimizer data members
  //-------------------------------------------
  opt->active = active;
  opt->number_of_parameters = number_of_parameters;
  opt->lower_param_bounds = lower_param_bounds;
  opt->upper_param_bounds = upper_param_bounds;
  opt->fitness_function = fitness_function;
  opt->number_of_generations = number_of_generations;
  opt->rng = rng;
  opt->pop = std::make_shared<population_t>(); //!! A copy of the contents, not the pointer
  opt->pop->addSolutions(*pop);
  opt->best = best;
  opt->average_fitness_history = average_fitness_history;
  opt->selection_fraction = selection_fraction;
  opt->init_univariate_bandwidth = init_univariate_bandwidth;

  // Stopping criteria
  //----------------------------------------------------------------------------------
  opt->maximum_no_improvement_stretch = maximum_no_improvement_stretch; ;
  opt->param_std_tolerance = param_std_tolerance;
  opt->fitness_std_tolerance = fitness_std_tolerance;

  // AMaLGaM data members
  //-------------------------------------------
  opt->mean = mean;
  opt->covariance = covariance;
  opt->cholesky = cholesky;
  opt->inverse_cholesky = inverse_cholesky;
  opt->no_improvement_stretch = no_improvement_stretch;
  opt->multiplier = multiplier;
  opt->old_mean = old_mean;
  opt->st_dev_ratio_threshold = st_dev_ratio_threshold;
  opt->distribution_multiplier_decrease = distribution_multiplier_decrease;
  opt->maximum_no_improvement_stretch = maximum_no_improvement_stretch;
  opt->sample_succes_ratio_threshold = sample_succes_ratio_threshold;
  opt->delta_ams = delta_ams;
  opt->apply_ams = apply_ams;

  return opt;
}

// Algorithm Name
std::string hillvallea::amalgam_univariate_t::name() const { return "AMaLGaM-Univariate"; }

// Initial initialization of the algorithm
// Population should be sorted on fitness (fittest first)
void hillvallea::amalgam_univariate_t::initialize_from_population(population_pt pop)
{
  this->pop = pop;
  multiplier = 1.0;
  no_improvement_stretch = 0;
  pop->mean(old_mean);
  mean = old_mean;
  pop->sort_on_fitness();
  best = *pop->sols[0]; 
}

size_t hillvallea::amalgam_univariate_t::recommended_popsize(const size_t problem_dimension) const
{
  return (size_t)std::max((double)((size_t)((2.0 / selection_fraction) + 1)), 10.0*pow((double)problem_dimension, 0.5));
}

void hillvallea::amalgam_univariate_t::update_distribution_multiplier(double & multiplier, const bool improvement, double & no_improvement_stretch, const double & sample_success_ratio, const double sdr) const
{

  // default variables;
  double sample_succes_ratio_threshold = 0.10;

  // if >90% of the samples is out of bounds, multiplier *0.5
  if (sample_success_ratio < sample_succes_ratio_threshold)
    multiplier *= 0.5;

  if (improvement)
  {
    no_improvement_stretch = 0.0;

    if (multiplier < 1.0)
      multiplier = 1.0;

    if (sdr > st_dev_ratio_threshold)
      multiplier /= distribution_multiplier_decrease;

  }
  else
  {

    if (multiplier <= 1.0)
      no_improvement_stretch++;

    if (multiplier > 1.0 || no_improvement_stretch >= maximum_no_improvement_stretch)
      multiplier *= distribution_multiplier_decrease;

    if (multiplier < 1.0 && no_improvement_stretch < maximum_no_improvement_stretch)
      multiplier = 1.0;

  }
}

// returns true if any of the termination criteria is satisfied
bool hillvallea::amalgam_univariate_t::checkTerminationCondition()
{

  if (number_of_generations == 0) {
    active = true;
    return !active;
  }

  // 1. if the cluster is empty, deactivate it.
  if (pop->size() == 0) {
    active = false;
    return !active;
  }

  // 2. check the maximum parameter variance
  // we scale the param std by it.
  double max_param_variance = 0.0;
  for (size_t i = 0; i < covariance.rows(); ++i) {
    if (covariance[i][i] > max_param_variance) {
      max_param_variance = covariance[i][i];
    }
  }

  vec_t mean;
  pop->mean(mean);

  // if the mean equals zero, we can't didivide by it, so terminate it when it is kinda small
  bool terminate_for_param_std_mean_zero = (mean.infinitynorm() <= 0 && sqrt(max_param_variance) < param_std_tolerance);
  bool terminate_on_parameter_std = sqrt(max_param_variance) / mean.infinitynorm() < param_std_tolerance;
  bool terminate_on_fitness_std = (pop->size() > 1) && (pop->relative_fitness_std() < fitness_std_tolerance);
  bool terminate_on_distribution_multiplier = (multiplier < 1e-10);

  if (terminate_for_param_std_mean_zero || terminate_on_parameter_std || terminate_on_fitness_std || terminate_on_distribution_multiplier)
  {
     active = false;
     return !active;
  }

  // if we have not terminated so far, the cluster is active.
  // if, due to selection, the cluster is shrunken, it is set to active again.
  active = true;
  return !active;

}

void hillvallea::amalgam_univariate_t::estimate_sample_parameters()
{

  // Compute sample mean and sample covariance
  old_mean = mean;

  // Change the focus of the search to the best solution 
  if (multiplier < 1.0)
    mean = pop->sols[0]->param;
  else
    pop->mean(mean);


  // if the population size is too small,
  // estimate a univariate covariance matrix
  if (pop->size() == 1)
  {
    covariance.setIdentity(mean.size(), mean.size());
    covariance.multiply(init_univariate_bandwidth*0.01);
  }
  else 
  {
    pop->covariance_univariate(mean, covariance);
  }
  // Cholesky decomposition
  choleskyDecomposition_univariate(covariance, cholesky);

  // apply the multiplier
  cholesky.multiply(sqrt(multiplier));

  // invert the cholesky decomposition
  int n = (int)covariance.rows();
  inverse_cholesky.setRaw(matrixLowerTriangularInverse(cholesky.toArray(), n), n, n);
  
}

// sample a new population
size_t hillvallea::amalgam_univariate_t::sample_new_population(const size_t sample_size)
{

  // Sample new population
  //----------------------------------------------------------------------------------------
  int number_of_samples = pop->fill_normal_univariate(sample_size, number_of_parameters, mean, cholesky, lower_param_bounds, upper_param_bounds, 1, rng);

  // apply the AMS
  if (apply_ams)
  {
    vec_t ams_direction = mean - old_mean;

    size_t number_of_ams_solutions = (size_t)(0.5*selection_fraction*sample_size); // alpha ams
    apply_ams_to_population(number_of_ams_solutions, delta_ams*multiplier, ams_direction); // ams, not to elite
  }

  // evaluate the population
  //---------------------------------------------------------------------------------------
  size_t number_of_evaluations = pop->evaluate(fitness_function, 1);
  pop->sort_on_fitness();

  // Update Params
  //---------------------------------------------------------------------------------------
  bool improvement = pop->improvement_over(best.f);
  double sdr = getSDR(best, mean, inverse_cholesky);
  double sample_success_ratio = (double)(sample_size - 1) / number_of_samples; // we do not sample the best.
  update_distribution_multiplier(multiplier, improvement, no_improvement_stretch, sample_success_ratio, sdr);
  best = *pop->first();

  number_of_generations++;

  return number_of_evaluations;
}





// Compute the SDR
//--------------------------------------------------------------------------------------
double hillvallea::amalgam_univariate_t::getSDR(const solution_t & best, const vec_t & mean, matrix_t & inverse_chol) const
{
  size_t i;

  // find improvements over the best.
  vec_t average_params(number_of_parameters, 0.0);
  for (i = 0; (i < pop->size()) && (pop->sols[i]->f < best.f); ++i) {
    average_params += pop->sols[i]->param;
  }

  if (i == 0)
    return 0.0;

  average_params /= (double)i;

  vec_t diff = average_params - mean;

  return inverse_chol.lowerProduct(diff).infinitynorm();

}


// Apply the Anticipated Mean Shift (AMS) to the first solutions (not the elite)
//------------------------------------------------------------------------
void hillvallea::amalgam_univariate_t::apply_ams_to_population(const size_t number_of_ams_solutions, const double ams_factor, const vec_t & ams_direction)
{

  // we have a shrink factor in case the shifts are
  // outside of the parameter bounds
  double shrink_factor;
  int attempts = 0;

  // loop over the first solutions to shift them,
  // but we save the elite.
  for (size_t i = 1; i < std::min(number_of_ams_solutions + 1, pop->sols.size()); ++i)
  {

    // try to sample within bounds
    shrink_factor = 2;

    // shift x.
    vec_t ams_params = pop->sols[i]->param;
    ams_params += shrink_factor * ams_factor * ams_direction;

    boundary_repair(pop->sols[i]->param, lower_param_bounds, upper_param_bounds);


    // try smaller shifts until the sol is within range
    while (attempts < 100 && !in_range(ams_params, lower_param_bounds, upper_param_bounds))
    {

      // if not, decrease the shrink_factor
      attempts++;
      shrink_factor *= 0.5;
      ams_params -= shrink_factor * ams_factor * ams_direction;

    }

    if (attempts < 100) {
      pop->sols[i]->param = ams_params;
    }

  }
}
