/*

HillVallEA

By S.C. Maree
s.c.maree[at]amc.uva.nl
github.com/SCMaree/HillVallEA

*/

#include "solution.hpp"

namespace hillvallea
{
  
  // initialize solution
  //----------------------------------------------
  solution_t::solution_t() {
    penalty = 0.0;
    elite = false;
    time_obtained = 0;
    feval_obtained = 0;
    generation_obtained = 0;
    cluster_number = -1;
    multiplier = 1.0;
    NormTabDis = 0.0;
  }

  solution_t::solution_t(size_t problem_size)
  {
    param.resize(problem_size,0.0);
    this->param_transformed.resize(problem_size, 0.0);
    penalty = 0.0;
    elite = false;
    time_obtained = 0;
    feval_obtained = 0;
    generation_obtained = 0;
    cluster_number = -1;
    multiplier = 1.0;
    NormTabDis = 0.0;
  }
  
  solution_t::solution_t(vec_t param)
  {
    penalty = 0.0;
    this->param = param;
    this->param_transformed.resize(this->param.size(), 0.0);
    elite = false;
    time_obtained = 0;
    feval_obtained = 0;
    generation_obtained = 0;
    cluster_number = -1;
    multiplier = 1.0;
    NormTabDis = 0.0;
  }
  
  solution_t::solution_t(const solution_t & other)
  {

    this->param = other.param;
    this->f = other.f;
    this->penalty = other.penalty;
    this->elite = other.elite;

    this->time_obtained = other.time_obtained;
    this->feval_obtained = other.feval_obtained;
    this->generation_obtained = other.generation_obtained;
    this->cluster_number = other.cluster_number;
    this->multiplier = other.multiplier;
    this->param_transformed = other.param_transformed;
    this->NormTabDis = other.NormTabDis;
  
  }
  
  // delete solution
  //----------------------------------------------
  solution_t::~solution_t() {}
  
  // comparison for solution_t pointers
  // is sol1 better than sol2?
  //-----------------------------------------------
  bool solution_t::better_solution_via_pointers(const solution_pt sol1, const solution_pt sol2) {
    return better_solution(*sol1, *sol2);
  }
  
  // defined as static!
  // returns true of the first solution is better than the second
  bool solution_t::better_solution(const solution_t & sol1, const solution_t & sol2)
  {
    
    if (sol1.penalty > 0)                     // sol1 is infeasible
    {
      if (sol2.penalty > 0)                   // both are infeasible
        return (sol1.penalty < sol2.penalty); // return the "most feasible"
      else
        return false;                         // only sol2 is feasible
    }
    else                                      // sol1 is feasible
    {
      if (sol2.penalty > 0)                   // only sol1 is feasible
        return true;
      else
        return (sol1.f < sol2.f);             // both are feasible
    }
  }
  
  // has sol1 a higher sample probability than sol2?
  // defined as static!
  bool solution_t::higher_probability(const std::shared_ptr<solution_t> sol1, const std::shared_ptr<solution_t> sol2)
  {
    return (sol1->probability > sol2->probability);
  }
  
  
  // computes the distance to another solution
  //------------------------------------------------------------------------------------
  double solution_t::param_distance(const solution_t & sol2) const
  {
    return (this->param - sol2.param).norm();
  }

  double solution_t::param_distance(const vec_t & param2) const
  {
    return (this->param - param2).norm();
  }
  
  
}







