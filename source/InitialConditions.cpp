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
		/*std::cout<< scaled_n_type_donor_doping << "   "
				 << scaled_p_type_donor_doping << "   "
				 << scaled_p_type_width <<std::endl;*/
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
		if(!schottky)
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
		else
		{
			if(p[0] <= scaled_schottky_depletion_width)
				return scaled_schottky_electron_density;
			if(p[0] > scaled_schottky_depletion_width
					&&
			   p[0] < (scaled_p_type_width - scaled_p_type_depletion_width) )
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
				const double & p_type_depletion_width,
				const bool   & schottky_status,
				const double & schottky_depletion_width,
				const double & schottky_electron_density)
	{
		scaled_n_type_donor_density   	= n_type_donor_density;
		scaled_p_type_donor_density   	= p_type_donor_density;
		scaled_p_type_width  		  	= p_type_width;
		scaled_n_type_depletion_width 	= n_type_depletion_width;
		scaled_p_type_depletion_width   = p_type_depletion_width;
		schottky						= schottky_status;
		scaled_schottky_depletion_width = schottky_depletion_width;
		scaled_schottky_electron_density= schottky_electron_density;

		if(schottky)
			std::cout << "schottky depletion width:"
					  << scaled_schottky_depletion_width
					  << "\nschottky electron density:"
					  << scaled_schottky_electron_density
					  << std::endl;

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
		if(!schottky)
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
		else // schottky
		{
			if(p[0] <= scaled_schottky_depletion_width )
				return scaled_schottky_hole_density;
			if(p[0] > scaled_schottky_depletion_width
					&&
			   p[0] <= (scaled_p_type_width - scaled_p_type_depletion_width) )
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
				const double & p_type_depletion_width,
				const bool   & schottky_status,
				const double & schottky_depletion_width,
				const double & schottky_hole_density)
	{
		scaled_n_type_acceptor_density  = n_type_acceptor_density;
		scaled_p_type_acceptor_density  = p_type_acceptor_density;
		scaled_p_type_width  		    = p_type_width;
		scaled_n_type_depletion_width   = n_type_depletion_width;
		scaled_p_type_depletion_width   = p_type_depletion_width;
		schottky						= schottky_status;
		scaled_schottky_depletion_width = schottky_depletion_width;
		scaled_schottky_hole_density    = schottky_hole_density;

		if(schottky)
			std::cout << "schottky hole density:"
					  << scaled_schottky_hole_density
					  << std::endl;

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
			if(p[0] < 1e-10 )
				return scaled_p_type_electron_bc;
			else
			{
				if(p[0] > (scaled_p_type_width + scaled_n_type_width - 1e-10) )
					return scaled_n_type_electron_bc;
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
	set_values( const double & n_type_acceptor_density,
				const double & p_type_acceptor_density,
				const double & n_type_donor_density,
				const double & p_type_donor_density,
				const double & n_type_width,
				const double & p_type_width,
				const double & intrinsic_density)
	{
		scaled_n_type_acceptor_density = n_type_acceptor_density;
		scaled_p_type_acceptor_density = p_type_acceptor_density;

		scaled_n_type_donor_density    = n_type_donor_density;
		scaled_p_type_donor_density    = p_type_donor_density;

		scaled_n_type_width  		   = n_type_width;
		scaled_p_type_width            = p_type_width;

		scaled_intrinsic_density       = intrinsic_density;

		const double n_type_resultant_doping = (scaled_n_type_donor_density - scaled_n_type_acceptor_density);
		const double p_type_resultant_doping = (scaled_p_type_donor_density - scaled_p_type_acceptor_density);

		scaled_n_type_electron_bc = 0.5*( n_type_resultant_doping + std::sqrt(n_type_resultant_doping*n_type_resultant_doping + 4*scaled_intrinsic_density*scaled_intrinsic_density));
		scaled_p_type_electron_bc = 0.5*( p_type_resultant_doping + std::sqrt(p_type_resultant_doping*p_type_resultant_doping + 4*scaled_intrinsic_density*scaled_intrinsic_density));

		std::cout << "electrons on n_type boundary:   "
				  << scaled_n_type_electron_bc
				  << std::endl
				  << "electrons on p_type boundary:   "
				  << scaled_p_type_electron_bc
				  << std::endl;


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
			if(p[0] < 1e-10 )
				return scaled_p_type_hole_bc;
			else
			{
				if(p[0] > (scaled_p_type_width + scaled_n_type_width - 1e-10) )
					return scaled_n_type_hole_bc;
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
				const double & n_type_donor_density,
				const double & p_type_donor_density,
				const double & n_type_width,
				const double & p_type_width,
				const double & intrinsic_density)
	{
		scaled_n_type_acceptor_density = n_type_acceptor_density;
		scaled_p_type_acceptor_density = p_type_acceptor_density;

		scaled_n_type_donor_density    = n_type_donor_density;
		scaled_p_type_donor_density    = p_type_donor_density;

		scaled_n_type_width  		   = n_type_width;
		scaled_p_type_width            = p_type_width;

		scaled_intrinsic_density       = intrinsic_density;

		const double n_type_resultant_doping = (scaled_n_type_donor_density - scaled_n_type_acceptor_density);
		const double p_type_resultant_doping = (scaled_p_type_donor_density - scaled_p_type_acceptor_density);

		scaled_n_type_hole_bc = 0.5*( -n_type_resultant_doping + std::sqrt(n_type_resultant_doping*n_type_resultant_doping + 4*scaled_intrinsic_density*scaled_intrinsic_density));
		scaled_p_type_hole_bc = 0.5*( -p_type_resultant_doping + std::sqrt(p_type_resultant_doping*p_type_resultant_doping + 4*scaled_intrinsic_density*scaled_intrinsic_density));

		std::cout << "holes on n_type boundary:   "
				  << scaled_n_type_hole_bc
				  << std::endl
				  << "holes on p_type boundary:   "
				  << scaled_p_type_hole_bc
				  << std::endl;


	}


	template <int dim>
	double
	LDG_Dirichlet_electron_density_bc_gb<dim>::
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
			if(p[0] > 1e-10 )
			{
				std::cerr << "Dirichlet condition on grain boundary can be only on the left contact! That is: x coordinate must equal to 0";
				return -1;
			}
			else
			{
				return scaled_electron_bc;

			}
		}
	}

	template <int dim>
	void
	LDG_Dirichlet_electron_density_bc_gb<dim>::
	set_values( const double & scaled_electron_bc_precalculated)
	{
		scaled_electron_bc = scaled_electron_bc_precalculated;
		std::cout << "electrons on a left contact with grain boundary:   "
				  << scaled_electron_bc
				  << std::endl;


	}

	template <int dim>
	double
	LDG_Dirichlet_hole_density_bc_gb<dim>::
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
			if(p[0] > 1e-10 )
			{
				std::cerr << "Dirichlet condition on grain boundary can be only on the left contact! That is: x coordinate must equal to 0";
				return -1;
			}

			else
				return scaled_hole_bc;
		}
	}


	template <int dim>
	void
	LDG_Dirichlet_hole_density_bc_gb<dim>::
	set_values( const double & scaled_hole_bc_precalculated)
	{
		scaled_hole_bc = scaled_hole_bc_precalculated;
		std::cout << "hole on a left contact with grain boundary:   "
				  << scaled_hole_bc
				  << std::endl;


	}

