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
			Convergence();
			Convergence(const DoFHandler<dim> & Poisson_dof_handler,const DoFHandler<dim> & continuity_dof_handler,
						const Vector<double>  & Poisson_solution   , const Vector<double> & Continuity_solution_electron,
																	 const Vector<double> & Continuity_solution_hole);
			void print_indexes();
			void print_residuums(unsigned int time_step);
			void calculate_residuum(const Vector<double> & Poisson_solution,
					                const Vector<double> & Continuity_solution_electron,
									const Vector<double> & Continuity_solution_hole);

			double sum_absolute_values(const Vector<double> & new_solution,
									   const Vector<double> & old_solution,
									   const unsigned int & begin,
									   const unsigned int & end);

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

			bool scale_factors_are_set;

	};


#endif /* INCLUDE_CONVERGENCE_HPP_ */
