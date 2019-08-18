#include "../include/InitialConditions.hpp"

	using namespace dealii;

	template <int dim>
	double
	Electrons_Equilibrium<dim>::
	value(const Point<dim> &p,
		  const unsigned int component) const
	{
			// dim+1 components
		if(component < dim)
		{	
			// set the components of the current initially to zero
			return ZeroFunction<dim>(dim+1).value(p, component);
		}
		else //(component == dim)
		{
			if(p[0] < scaled_p_type_width)
				return scaled_p_type_donor_doping;
			else
				return scaled_n_type_donor_doping;
		}
	}

	template <int dim>
	void
	Electrons_Equilibrium<dim>::
	set_values(const double & n_type_donor_doping, const double & p_type_donor_doping, const double & p_type_width)
	{
		scaled_n_type_donor_doping = n_type_donor_doping;
		scaled_p_type_donor_doping = p_type_donor_doping;
		scaled_p_type_width  = p_type_width;

	}



	template <int dim>
	double
	Holes_Equilibrium<dim>::
	value(const Point<dim> &p,
		  const unsigned int component) const
	{
		// dim+1 components
		if(component < dim)
		{	
			// set the components of the current initially to zero
			return ZeroFunction<dim>(dim+1).value(p, component);
		}
		else //(component == dim)
		{
			if(p[0] < scaled_p_type_width)
				return scaled_p_type_acceptor_doping;
			else
				return scaled_n_type_acceptor_doping;
		}
	}


	template <int dim>
	void
	Holes_Equilibrium<dim>::
	set_values(const double & n_type_acceptor_doping, const double & p_type_acceptor_doping, const double & p_type_width)
	{
		scaled_n_type_acceptor_doping = n_type_acceptor_doping;
		scaled_p_type_acceptor_doping = p_type_acceptor_doping;
		scaled_p_type_width  = p_type_width;

	}

	template <int dim>
	double
	Reductants_Equilibrium<dim>::
	value(const Point<dim> &p,
		  const unsigned int component) const
	{
			// dim+1 components
		if(component < dim)
		{	
			// set the components of the current initially to zero
			return ZeroFunction<dim>(dim+1).value(p, component);
		}
		else //(component == dim)
		{
				return 30; 
		}
	}

	template <int dim>
	double
	Oxidants_Equilibrium<dim>::
	value(const Point<dim> &p,
		  const unsigned int component) const
	{
			// dim+1 components
		if(component < dim)
		{	
			// set the components of the current initially to zero
			return ZeroFunction<dim>(dim+1).value(p, component);
		}
		else //(component == dim)
		{
				return 29; 
		}
	}

