/*
 * Convergence.cpp
 *
 *  Created on: Sep 21, 2019
 *      Author: kon
 */

#include "../include/Convergence.hpp"


/*double inner_product_accumulator(double x, double y)
{
	return x + std::abs(y);
}

double inner_product_product(double x, double y)
{
	return x - y;
}*/


template<int dim>
Convergence<dim>::
Convergence() :
scale_factors_are_set(false),
electric_field_residuum(0.0),
potential_residuum(0.0),
jx_hole_residuum(0.0),
jy_hole_residuum(0.0),
hole_density_residuum(0.0),
jx_electron_residuum(0.0),
jy_electron_residuum(0.0),
electron_density_residuum(0.0),
electric_field_scale_factor(0.0),
potential_scale_factor(0.0),
jx_hole_scale_factor(0.0),
jy_hole_scale_factor(0.0),
hole_density_scale_factor(0.0),
jx_electron_scale_factor(0.0),
jy_electron_scale_factor(0.0),
electron_density_scale_factor(0.0),
old_potential_res(10.0),
old_elec_density_res(10.0),
old_hole_density_res(10.0),
old_electric_field_res(10.0),
steady_state(false)
{
}

template<int dim>
Convergence<dim>::
Convergence(const DoFHandler<dim> & Poisson_dof_handler,const DoFHandler<dim> & LDG_dof_handler,
			const Vector<double>  & Poisson_solution   , const Vector<double> & Continuity_solution_electron,
														 const Vector<double> & Continuity_solution_hole):

scale_factors_are_set(false),
old_Poisson_solution(Poisson_solution),
old_Continuity_solution_electrons(Continuity_solution_electron),
old_Continuity_solution_holes(Continuity_solution_hole),
electric_field_residuum(0.0),
potential_residuum(0.0),
jx_hole_residuum(0.0),
jy_hole_residuum(0.0),
hole_density_residuum(0.0),
jx_electron_residuum(0.0),
jy_electron_residuum(0.0),
electron_density_residuum(0.0),
electric_field_scale_factor(0.0),
potential_scale_factor(0.0),
jx_hole_scale_factor(0.0),
jy_hole_scale_factor(0.0),
hole_density_scale_factor(0.0),
jx_electron_scale_factor(0.0),
jy_electron_scale_factor(0.0),
electron_density_scale_factor(0.0),
old_potential_res(10.0),
old_elec_density_res(10.0),
old_hole_density_res(10.0),
old_electric_field_res(10.0),
steady_state(false)
{


	std::vector<types::global_dof_index> dofs_per_component_Poisson(dim+1), dofs_per_component_LDG(dim+1);
	DoFTools::count_dofs_per_component(Poisson_dof_handler, dofs_per_component_Poisson);
	DoFTools::count_dofs_per_component(LDG_dof_handler, dofs_per_component_LDG);

	const unsigned int n_electric_field = dofs_per_component_Poisson[0];
	const unsigned int n_potential		= dofs_per_component_Poisson[dim];

	const unsigned int n_current        = dofs_per_component_LDG[0];
	const unsigned int n_density		= dofs_per_component_LDG[dim];

	electric_field_indexes.first  = 0;
	electric_field_indexes.second =   n_electric_field-1;

	potential_indexes.first         = n_electric_field;
	potential_indexes.second		= n_electric_field+n_potential-1;


	x_current_indexes.first  = 0;
	x_current_indexes.second =   n_current-1;

	y_current_indexes.first  =   n_current;
	y_current_indexes.second = 2*n_current-1;

	density_indexes.first    = 2*n_current;
	density_indexes.second	 = 2*n_current+n_density-1;

}


template<int dim>
void
Convergence<dim>::
set_old_solutions(const Vector<double>  & Poisson_solution,
				  const Vector<double>  & Continuity_solution_electron,
				  const Vector<double>  & Continuity_solution_hole)
{
	old_Poisson_solution 			  = Poisson_solution;
	old_Continuity_solution_electrons = Continuity_solution_electron;
	old_Continuity_solution_holes     = Continuity_solution_hole;
	old_potential_res    = 10.0;
	old_electric_field_res=10.0;
	old_elec_density_res = 10.0;
	old_hole_density_res = 10.0;
}

template<int dim>
void
Convergence<dim>::
calculate_residuum(const Vector<double> & Poisson_solution,
				   const Vector<double> & Continuity_solution_electron,
				   const Vector<double> & Continuity_solution_hole)
{
	electric_field_residuum		= sum_absolute_values(Poisson_solution, old_Poisson_solution, electric_field_indexes.first, electric_field_indexes.second);
	potential_residuum   		= sum_absolute_values(Poisson_solution, old_Poisson_solution, potential_indexes.first, potential_indexes.second);

	jx_electron_residuum		= sum_absolute_values(Continuity_solution_electron, old_Continuity_solution_electrons, x_current_indexes.first, x_current_indexes.second);
	jy_electron_residuum		= sum_absolute_values(Continuity_solution_electron, old_Continuity_solution_electrons, y_current_indexes.first, y_current_indexes.second);
	electron_density_residuum	= sum_absolute_values(Continuity_solution_electron, old_Continuity_solution_electrons, density_indexes.first, density_indexes.second);

	jx_hole_residuum			= sum_absolute_values(Continuity_solution_hole, old_Continuity_solution_holes, x_current_indexes.first, x_current_indexes.second);
	jy_hole_residuum			= sum_absolute_values(Continuity_solution_hole, old_Continuity_solution_holes, y_current_indexes.first, y_current_indexes.second);
	hole_density_residuum		= sum_absolute_values(Continuity_solution_hole, old_Continuity_solution_holes, density_indexes.first, density_indexes.second);


	if(!scale_factors_are_set)
	{
		std::cout << "Set up scale factors!" <<std::endl;

		print_residuums(0);
		std::cout<<std::endl;

		electric_field_scale_factor = electric_field_residuum;
		potential_scale_factor      = potential_residuum;

		jx_electron_scale_factor 	   = jx_electron_residuum;
		jy_electron_scale_factor 	   = jy_electron_residuum;
		electron_density_scale_factor  = electron_density_residuum;

		jx_hole_scale_factor 	   = jx_hole_residuum;
		jy_hole_scale_factor 	   = jy_hole_residuum;
		hole_density_scale_factor  = hole_density_residuum;
		scale_factors_are_set 	   = true;
	}


	electric_field_residuum /= electric_field_scale_factor;
	potential_residuum 	    /= potential_scale_factor;

	jx_electron_residuum 		/= jx_electron_scale_factor;
	jy_electron_residuum 		/= jy_electron_scale_factor;
	electron_density_residuum 	/= electron_density_scale_factor;

	jx_hole_residuum 		/= jx_hole_scale_factor;
	jy_hole_residuum 		/= jy_hole_scale_factor;
	hole_density_residuum 	/= hole_density_scale_factor;

	steady_state = check_steady_state();
	//std::cout<< "steady state:   " << steady_state << std::endl;

	old_potential_res      = potential_residuum;
	old_electric_field_res = electric_field_residuum;
	old_elec_density_res   = electron_density_residuum;
	old_hole_density_res   = hole_density_residuum;

	old_Poisson_solution = Poisson_solution;
	old_Continuity_solution_electrons = Continuity_solution_electron;
	old_Continuity_solution_holes = Continuity_solution_hole;

}

template<int dim>
double
Convergence<dim>::
sum_absolute_values(const Vector<double> & new_solution,
				    const Vector<double> & old_solution,
				    const unsigned int & begin,
				    const unsigned int & end)
{
	double sum=0.0;
	for(unsigned int i = begin; i <= end; i++)
	{
		sum += std::abs(new_solution[i] - old_solution[i]);
	}

	return sum;
}

template<int dim>
bool
Convergence<dim>::check_steady_state()
{
	return (
			(-potential_residuum        + old_potential_res)    < 0.0 /*&&
			 -electric_field_residuum   + old_electric_field_res< 0.0*/&&
			(-hole_density_residuum     + old_hole_density_res) < 0.0 /*||
			(-electron_density_residuum + old_elec_density_res) < 0.0*/
			);
}

template<int dim>
bool
Convergence<dim>::sanity_check()
{
	if(std::isnan(potential_residuum) ||
	   std::isnan(electric_field_residuum) ||
	   std::isnan(jx_electron_residuum) ||
	   std::isnan(jy_electron_residuum) ||
	   std::isnan(electron_density_residuum) ||
	   std::isnan(jx_hole_residuum) ||
	   std::isnan(jy_hole_residuum) ||
	   std::isnan(hole_density_residuum) )
		return false;
	else
		return true;
}

template<int dim>
void
Convergence<dim>::
print_residuums(unsigned int time_step)
{
	std::cout << "Time step:     "
			  << time_step
			  <<std::endl

			  << "  E:            "
			  << electric_field_residuum
			  << "    potential:  "
			  << potential_residuum
			  <<std::endl

			  << "  jx electron:  "
			  << jx_electron_residuum
			  << "   jy electron: "
			  << jy_electron_residuum
			  << " electron den:  "
			  << electron_density_residuum
			  << std::endl


			  << "  jx_hole:      "
			  << jx_hole_residuum
			  << "   jy_hole:     "
			  << jy_hole_residuum
			  << " hole density: "
			  << hole_density_residuum
			  << std::endl
			  << std::endl;
}

template<int dim>
void
Convergence<dim>::
print_indexes()
{
	std::cout << "Indexes of particular variables:" << std::endl;
	std::cout << "E  from:  " << electric_field_indexes.first <<"      to:" << electric_field_indexes.second << std::endl;
	std::cout << "Phi from: " << potential_indexes.first        <<"   to:" << potential_indexes.second        << std::endl;

	std::cout << std::endl;
	std::cout << "jx from:  " << x_current_indexes.first <<"      to:" << x_current_indexes.second << std::endl;
	std::cout << "jy from:  " << y_current_indexes.first <<"   to:" << y_current_indexes.second << std::endl;
	std::cout << "rho from: " << density_indexes.first   <<"   to:" << density_indexes.second   << std::endl;
	std::cout << std::endl;
}
