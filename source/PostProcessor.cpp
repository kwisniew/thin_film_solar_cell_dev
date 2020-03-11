#include "../include/PostProcessor.hpp"

using namespace dealii;


template<int dim>
PostProcessor<dim>::
PostProcessor(   const ParameterSpace::Parameters 	& sim_params,
				 const bool							& print_carrier,
				 const std::string					& name)
{
	scale_potential  = sim_params.thermal_voltage;
	scale_elec_field = sim_params.thermal_voltage/
					  (sim_params.characteristic_length*(sim_params.scaled_debye_length*sim_params.semiconductor_permittivity));
	scale_current	 = PhysicalConstants::electron_charge * sim_params.characteristic_denisty  /
			   	   	   (sim_params.characteristic_time * sim_params.characteristic_length);

	printing_carrier = print_carrier;	
	if(print_carrier)
	{
		density_name  = name.c_str();
		density_name += "Density";
		current_name  = name.c_str();
		current_name += "Current";
	}
}

template<int dim>
std::vector<std::string> 
PostProcessor<dim>::
get_names() const
{
	std::vector<std::string> solution_names;

	if(printing_carrier)
	{
		for(unsigned int d=0; d<dim; d++)
			solution_names.push_back(current_name.c_str());

		solution_names.push_back(density_name.c_str());	
	}
	else
	{
		for(unsigned int d=0; d<dim; d++)
			solution_names.push_back("Field");

		solution_names.push_back("Potential");
	} 
	return solution_names;
}

template<int dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation>
PostProcessor<dim>::
get_data_component_interpretation() const
{

	std::vector<DataComponentInterpretation::DataComponentInterpretation>
	interpretation;

	for(unsigned int d=0; d<dim; d++)
	{
		interpretation.push_back(DataComponentInterpretation::
					component_is_part_of_vector);
	}
	interpretation.push_back(DataComponentInterpretation::component_is_scalar);

	return interpretation;
}

template<int dim>
UpdateFlags
PostProcessor<dim>::
get_needed_update_flags() const
{
	return update_values;// | update_gradients | update_quadrature_points;
}


template<int dim>
void 
PostProcessor<dim>::
evaluate_vector_field
(const DataPostprocessorInputs::Vector<dim> &inputs,
       std::vector<Vector<double> >         &computed_quantities) const
{

	
	const unsigned int n_quadrature_points = inputs.solution_values.size();
	Assert(computed_quantities.size() == n_quadrature_points, ExcInternalError());
	if(printing_carrier)
	{
		for(unsigned int q=0; q<n_quadrature_points; q++)
		{
			// copy over current with scaling.
			for(unsigned int d=0; d<dim; d++)
				computed_quantities[q](d) = (scale_current * inputs.solution_values[q](d));
	
			// copy over densities without scaling
			computed_quantities[q](dim) = inputs.solution_values[q](dim);
		} // for q
	}
	else
	{
		for(unsigned int q=0; q<n_quadrature_points; q++)
		{
			// copy over electric field...
			for(unsigned int d=0; d<dim; d++)
				computed_quantities[q](d) = (scale_elec_field * inputs.solution_values[q](d));
	
			// copy over potential
			const double potential = scale_potential * inputs.solution_values[q](dim);
			computed_quantities[q](dim) = potential;
		} // for q
	} // else

}



/*
 *
 *
 *  CALCULATE ALL CURRENT
 *
 *
 */


template<int dim>
PostProcessor_currents<dim>::
PostProcessor_currents(   const ParameterSpace::Parameters 	& sim_params,
				 	  	  const std::string					& e_name,
						  const std::string					& h_name)
{
	elec_name  = e_name.c_str();
	hole_name  = h_name.c_str();

	double scale_elec_field = sim_params.thermal_voltage/
						  (sim_params.characteristic_length*(sim_params.scaled_debye_length*sim_params.semiconductor_permittivity));
	scale_dryf_current     = PhysicalConstants::electron_charge*sim_params.mobility*sim_params.characteristic_denisty*scale_elec_field;


	double scale_gradient	  =  sim_params.characteristic_denisty  / sim_params.characteristic_length;
	double diffusion_constant =  sim_params.mobility * sim_params.thermal_voltage;
	scale_diff_current   =  PhysicalConstants::electron_charge*diffusion_constant*scale_gradient;
}

template<int dim>
std::vector<std::string>
PostProcessor_currents<dim>::
get_names() const
{
	std::vector<std::string> solution_names;

	solution_names.push_back(elec_name + "_Dryf_Current");
	solution_names.push_back(elec_name + "_Dryf_Current");

	solution_names.push_back(elec_name + "_Diff_Current");
	solution_names.push_back(elec_name + "_Diff_Current");

	solution_names.push_back(elec_name + "_totalCurrent");
	solution_names.push_back(elec_name + "_totalCurrent");


	solution_names.push_back(hole_name + "_Dryf_Current");
	solution_names.push_back(hole_name + "_Dryf_Current");

	solution_names.push_back(hole_name + "_Diff_Current");
	solution_names.push_back(hole_name + "_Diff_Current");

	solution_names.push_back(hole_name + "_totalCurrent");
	solution_names.push_back(hole_name + "_totalCurrent");

	solution_names.push_back("Total_Current");
	solution_names.push_back("Total_Current");

	return solution_names;
}

template<int dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation>
PostProcessor_currents<dim>::
get_data_component_interpretation() const
{

	std::vector<DataComponentInterpretation::DataComponentInterpretation>
	interpretation;

	for(unsigned int i=0;i<7;i++)
	{
		for(unsigned int d=0; d<dim; d++)
		{
			interpretation.push_back(DataComponentInterpretation::
						component_is_part_of_vector);
		}
	}


	return interpretation;
}

template<int dim>
UpdateFlags
PostProcessor_currents<dim>::
get_needed_update_flags() const
{
	return update_values | update_gradients;// | update_quadrature_points;
}


template<int dim>
void
PostProcessor_currents<dim>::
evaluate_vector_field
(const DataPostprocessorInputs::Vector<dim> &inputs,
       std::vector<Vector<double> >         &computed_quantities) const
{


	const unsigned int n_quadrature_points = inputs.solution_values.size();
	if(computed_quantities.size() != n_quadrature_points) std::cerr << "Wielkosci sie nie zgadzaja!\n";
	Assert(computed_quantities.size() == n_quadrature_points, ExcInternalError());

	for(unsigned int q=0; q<n_quadrature_points; q++)
	{
		for(unsigned int d=0; d<dim; d++){

			computed_quantities[q](d)    = (scale_dryf_current * inputs.solution_values[q](d)* inputs.solution_values[q](5));
			computed_quantities[q](2+d)  = (scale_diff_current * inputs.solution_gradients[q][5][d]);
			computed_quantities[q](4+d)  = computed_quantities[q](d) + computed_quantities[q](2+d);

			computed_quantities[q](6+d)  = (scale_dryf_current * inputs.solution_values[q](d)* inputs.solution_values[q](8));
			computed_quantities[q](8+d)  = (-scale_diff_current * inputs.solution_gradients[q][8][d]);
			computed_quantities[q](10+d) = computed_quantities[q](6+d) + computed_quantities[q](8+d);

			computed_quantities[q](12+d) = computed_quantities[q](4+d) + computed_quantities[q](10+d);
		}


	} // for q


}
