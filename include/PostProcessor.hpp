#ifndef _POSTPROCESSOR_H__
#define _POSTPROCESSOR_H__

//#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_postprocessor.h>
#include <deal.II/lac/vector.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/exceptions.h>
#include "Parameters.hpp"
#include "CarrierPair.hpp"
#include "Assembly.hpp"
#include <string>

using namespace dealii;

/** \brief This will add back in units and calculate the current.*/
/** \note This is based off of step-32 in deal.ii examples.*/
template <int dim>
class PostProcessor : public DataPostprocessor<dim>
{
	public:

		/** \brief Constructor for Poisson post processor.*/
		/** The constructor will assign the scaling values to the object.*/
		PostProcessor(const ParameterSpace::Parameters 	& sim_parms,
			      const bool			& print_carrier,
			      const std::string			& name);

		/** \brief post processing done in this function.*/
		/** \note: This is called internally from DatOut::build_patches*/
		virtual
		void 
		evaluate_vector_field
		(const DataPostprocessorInputs::Vector<dim> &inputs,
			   std::vector<Vector<double> >         &computed_quantities) const;


		/** \brief Returns a vector containing the names of the data componetnt.*/
		virtual std::vector<std::string> get_names() const;
	
		/** \brief returns a vector which tells whether data comp. is scalar or vector.*/
		virtual
		std::vector<DataComponentInterpretation::DataComponentInterpretation>
		get_data_component_interpretation() const;

		/** \brief Returns the values which must be update every print.*/
		virtual UpdateFlags get_needed_update_flags() const;


	private:
		//double scale_density;
		double scale_current;
		double scale_potential;
		double scale_elec_field;
		//double material_permittivity;
		std::string density_name;
		std::string current_name;
		bool printing_carrier;

};


template <int dim>
class PostProcessor_currents : public DataPostprocessor<dim>
{
	public:

		/** \brief Constructor for Currents post processor.*/
		/** The constructor will assign the scaling values to the object.*/
		PostProcessor_currents(   const ParameterSpace::Parameters 	& sim_params,
								  const std::string					& e_name,
								  const std::string					& h_name);

		/** \brief post processing done in this function.*/
		/** \note: This is called internally from DatOut::build_patches*/
		virtual
		void
		evaluate_vector_field
		(const DataPostprocessorInputs::Vector<dim> &inputs,
			   std::vector<Vector<double> >         &computed_quantities) const;


		/** \brief Returns a vector containing the names of the data componetnt.*/
		virtual std::vector<std::string> get_names() const;

		/** \brief returns a vector which tells whether data comp. is scalar or vector.*/
		virtual
		std::vector<DataComponentInterpretation::DataComponentInterpretation>
		get_data_component_interpretation() const;

		/** \brief Returns the values which must be update every print.*/
		virtual UpdateFlags get_needed_update_flags() const;


	private:

		double scale_dryf_current;
		double scale_diff_current;
		std::string elec_name;
		std::string hole_name;

};



#endif 
