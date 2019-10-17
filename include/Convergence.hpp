/*
 * Convergence.hpp
 *
 *  Created on: Sep 21, 2019
 *      Author: kon
 */

#ifndef INCLUDE_CONVERGENCE_HPP_
#define INCLUDE_CONVERGENCE_HPP_

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/vector.h>
#include <vector>
#include <numeric>
#include <utility>

using namespace dealii;

template <int dim>
	class Convergence
	{
		public:
			bool scale_factors_are_set;

			Convergence();
			Convergence(const DoFHandler<dim> & Poisson_dof_handler,const DoFHandler<dim> & continuity_dof_handler,
						const Vector<double>  & Poisson_solution   , const Vector<double> & Continuity_solution_electron,
																	 const Vector<double> & Continuity_solution_hole);

			void set_old_solutions(const Vector<double>  & Poisson_solution,
								   const Vector<double>  & Continuity_solution_electron,
								   const Vector<double>  & Continuity_solution_hole);

			void print_indexes();
			void print_residuums(unsigned int time_step);

			/** \brief Calculate the difference between solution in this step and previous step */
			/** The main purpose of this function is to calculate for each variable (E,phi,jx,jy,rho)
			 * the sum of an absolute values of the difference between solution in this step and
			 * previous step.
			 * The function loop over each degrees of freedom calculating:
			 * x_residuum += | x(tn) - x(t(n-1) |
			 *
			 * Additionally the scale factor for each variable is set on the first run of this function
			 * and the new solutions and residuums became after this function "old" solutions and residuums
			 */
			void calculate_residuum(const Vector<double> & Poisson_solution,
					                const Vector<double> & Continuity_solution_electron,
									const Vector<double> & Continuity_solution_hole);

			double sum_absolute_values(const Vector<double> & new_solution,
									   const Vector<double> & old_solution,
									   const unsigned int & begin,
									   const unsigned int & end);

			/** \brief check if solutions start to oscillate */
			/** I assumed that the residuums will decrease over the time
			 *  so the difference between old one and new one will be positive
			 *  if this not be a case it means that the solutions start to oscillate
			 *  so we can announce that we achieve steady state (or equilibrium state)
			 */
			bool check_steady_state();
			bool sanity_check();

			bool get_steady_state_idicator() const {return steady_state;}
			std::pair<unsigned int, unsigned int> get_density_indexes() const {return density_indexes;}
			std::pair<unsigned int, unsigned int> get_E_field_indexes() const {return electric_field_indexes;}


		private:

			std::pair<unsigned int, unsigned int> electric_field_indexes;
			std::pair<unsigned int, unsigned int> potential_indexes;

			std::pair<unsigned int, unsigned int> x_current_indexes;
			std::pair<unsigned int, unsigned int> y_current_indexes;
			std::pair<unsigned int, unsigned int> density_indexes;

			Vector<double> old_Poisson_solution;
			Vector<double> old_Continuity_solution_electrons;
			Vector<double> old_Continuity_solution_holes;

			double electric_field_residuum;
			double potential_residuum;

			double jx_hole_residuum;
			double jy_hole_residuum;
			double hole_density_residuum;

			double jx_electron_residuum;
			double jy_electron_residuum;
			double electron_density_residuum;

			double electric_field_scale_factor;
			double potential_scale_factor;

			double jx_hole_scale_factor;
			double jy_hole_scale_factor;
			double hole_density_scale_factor;

			double jx_electron_scale_factor;
			double jy_electron_scale_factor;
			double electron_density_scale_factor;

			//for checking if relative residuum is changed less by 1%
			double old_potential_res;
			double old_elec_density_res;
			double old_hole_density_res;

			bool steady_state;

	};


#endif /* INCLUDE_CONVERGENCE_HPP_ */
