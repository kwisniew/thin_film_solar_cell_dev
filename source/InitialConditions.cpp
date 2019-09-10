#include "../include/InitialConditions.hpp"

	using namespace dealii;

	template <int dim>
	double
	Doping_profile_donors<dim>::
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
	Doping_profile_donors<dim>::
	set_values(const double & n_type_donor_doping, const double & p_type_donor_doping, const double & p_type_width)
	{
		scaled_n_type_donor_doping = n_type_donor_doping;
		scaled_p_type_donor_doping = p_type_donor_doping;
		scaled_p_type_width  = p_type_width;

	}



	template <int dim>
	double
	Doping_profile_acceptors<dim>::
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
	Doping_profile_acceptors<dim>::
	set_values(const double & n_type_acceptor_doping, const double & p_type_acceptor_doping, const double & p_type_width)
	{
		scaled_n_type_acceptor_doping = n_type_acceptor_doping;
		scaled_p_type_acceptor_doping = p_type_acceptor_doping;
		scaled_p_type_width  = p_type_width;

	}

	template <int dim>
	double
	Initial_condition_electrons<dim>::
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
			if(p[0] < (scaled_p_type_width - scaled_p_type_depletion_width) )
				return scaled_p_type_donor_density;
			else
			{
				if(p[0] < (scaled_p_type_width + scaled_n_type_depletion_width) )
					return 0;
				else
					return scaled_n_type_donor_density;
			}
		}
	}

	template <int dim>
	void
	Initial_condition_electrons<dim>::
	set_values( const double & n_type_donor_density,
				const double & p_type_donor_density,
				const double & p_type_width,
				const double & n_type_depletion_width,
				const double & p_type_depletion_width)
	{
		scaled_n_type_donor_density   = n_type_donor_density;
		scaled_p_type_donor_density   = p_type_donor_density;
		scaled_p_type_width  		  = p_type_width;
		scaled_n_type_depletion_width = n_type_depletion_width;
		scaled_p_type_depletion_width = p_type_depletion_width;

	}


	template <int dim>
	double
	Initial_condition_holes<dim>::
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
			if(p[0] < (scaled_p_type_width - scaled_p_type_depletion_width) )
				return scaled_p_type_acceptor_density;
			else
			{
				if(p[0] < (scaled_p_type_width + scaled_n_type_depletion_width) )
					return 0;
				else
					return scaled_n_type_acceptor_density;
			}
		}
	}

	template <int dim>
	void
	Initial_condition_holes<dim>::
	set_values( const double & n_type_acceptor_density,
				const double & p_type_acceptor_density,
				const double & p_type_width,
				const double & n_type_depletion_width,
				const double & p_type_depletion_width)
	{
		scaled_n_type_acceptor_density = n_type_acceptor_density;
		scaled_p_type_acceptor_density = p_type_acceptor_density;
		scaled_p_type_width  		   = p_type_width;
		scaled_n_type_depletion_width  = n_type_depletion_width;
		scaled_p_type_depletion_width  = p_type_depletion_width;

	}


	template <int dim>
	double
	LDG_Dirichlet_electron_density_bc<dim>::
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
			if(p[0] < 1e-6 )
				return scaled_p_type_donor_density;
			else
			{
				if(p[0] > (scaled_p_type_width + scaled_n_type_width - 1e-6) )
					return scaled_n_type_donor_density;
				else
				{
					//std::cout<<"This situation is not implemented!" <<std::endl;
					return 0;
				}

			}
		}
	}

	template <int dim>
	void
	LDG_Dirichlet_electron_density_bc<dim>::
	set_values( const double & n_type_donor_density,
				const double & p_type_donor_density,
				const double & p_type_width,
				const double & n_type_width)
	{
		scaled_n_type_donor_density = n_type_donor_density;
		scaled_p_type_donor_density = p_type_donor_density;
		scaled_p_type_width  		= p_type_width;
		scaled_n_type_width         = n_type_width;

	}

	template <int dim>
	double
	LDG_Dirichlet_hole_density_bc<dim>::
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
			if(p[0] < 1e-6 )
				return scaled_p_type_acceptor_density;
			else
			{
				if(p[0] > (scaled_p_type_width + scaled_n_type_width - 1e-6) )
					return scaled_n_type_acceptor_density;
				else
				{
					//std::cout<<"This situation is not implemented!" <<std::endl;
					return 0;
				}

			}
		}
	}

	template <int dim>
	void
	LDG_Dirichlet_hole_density_bc<dim>::
	set_values( const double & n_type_acceptor_density,
				const double & p_type_acceptor_density,
				const double & p_type_width,
				const double & n_type_width)
	{
		scaled_n_type_acceptor_density = n_type_acceptor_density;
		scaled_p_type_acceptor_density = p_type_acceptor_density;
		scaled_p_type_width  		= p_type_width;
		scaled_n_type_width         = n_type_width;
	}
