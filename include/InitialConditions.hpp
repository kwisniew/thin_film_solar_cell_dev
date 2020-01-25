#ifndef _INITIALCONDITIONS_H__
#define _INITIALCONDITIONS_H__

#include<deal.II/base/function.h>
#include <iostream>

using namespace dealii;

/** \brief Set doping profile for donors (e.g. 0 on p-type side and Nd for n-type side.*/
template<int dim>
class Doping_profile_donors : public Function<dim>
{
	public:
		/** \brief Default constructor. */
		Doping_profile_donors() : Function<dim>(dim+1)
		{};
		virtual ~Doping_profile_donors(){}

		void set_values(const double & n_type_donor_doping, const double & p_type_donor_doping, const double & p_type_width);

		/** \brief Returns value of \f$\rho_{n}^{e}\f$ at point p.*/
		/**
		 * Right now the component work only for value of "dim" - look at the function .value_list() eg. in
		 * assemble_local_semiconductor_rhs()
		 */
		virtual double value(const Point<dim> &p,
					 const unsigned int component = 0) const;

	private:
		double scaled_n_type_donor_doping;
		double scaled_p_type_donor_doping;
		double scaled_p_type_width;
};

	/** \brief Set doping profile for acceptors (e.g. Na on p-type side and 0 for n-type side.*/
	template<int dim>
	class Doping_profile_acceptors : public Function<dim>
	{
		public:
			/** \brief Default constructor. */
			Doping_profile_acceptors() : Function<dim>(dim+1)
			{}

			void set_values(const double & n_type_acceptor_doping, const double & p_type_acceptor_doping, const double & p_type_width);

			/** \brief Returns value of \f$\rho_{p}^{e}\f$ at point p.*/
			virtual double value(const Point<dim> & p,
					     const unsigned int component=0) const;

		private:
			double scaled_n_type_acceptor_doping;
			double scaled_p_type_acceptor_doping;
			double scaled_p_type_width;
	};

	/** \brief Initial conditions for electrons in full depletion approximation.
	 * That is: zero carrier density near to pn interface.
	 *  \f$ \rho_{r}^{\infty} \f$.*/
	template<int dim>
	class Initial_condition_electrons : public Function<dim>
	{
		public:
			/** \brief Default constructor. */
			Initial_condition_electrons() :
			Function<dim>(dim+1),
			scaled_n_type_donor_density(0),
			scaled_p_type_donor_density(0),
			scaled_p_type_width(0),
			scaled_n_type_depletion_width(0),
			scaled_p_type_depletion_width(0),
			scaled_schottky_depletion_width(0),
			scaled_schottky_electron_density(0),
			schottky(false)
			{}
			virtual ~Initial_condition_electrons(){}

			void set_values(const double & n_type_donor_density,
							const double & p_type_donor_density,
							const double & p_type_width,
							const double & n_type_depletion_width,
							const double & p_type_depletion_width,
							const bool   & schottky_status = false,
							const double & schottky_depletion_width  = 0,
							const double & schottky_electron_density = 0);


			/** \brief Returns value of \f$\rho_{r}^{\infty}\f$ at point p.*/
			virtual double value(const Point<dim> & p,
					     const unsigned int component=0) const;

		private:
			double scaled_n_type_donor_density;
			double scaled_p_type_donor_density;
			double scaled_p_type_width;
			double scaled_n_type_depletion_width;
			double scaled_p_type_depletion_width;
			double scaled_schottky_depletion_width;
			double scaled_schottky_electron_density;

			bool schottky;
	};



	/** \brief Initial conditions for electrons in full depletion approximation.
	 * That is: zero carrier density near to pn interface.
	 *  \f$ \rho_{o}^{\infty} \f$.*/
	template<int dim>
	class Initial_condition_holes : public Function<dim>
	{
		public:
			/** \brief Default constructor. */
			Initial_condition_holes() :
			Function<dim>(dim+1),
			scaled_n_type_acceptor_density(0),
			scaled_p_type_acceptor_density(0),
			scaled_p_type_width(0),
			scaled_n_type_depletion_width(0),
			scaled_p_type_depletion_width(0),
			scaled_schottky_depletion_width(0),
			scaled_schottky_hole_density(0),
			schottky(false)
			{}
			virtual ~Initial_condition_holes(){}
			void set_values(const double & n_type_acceptor_density,
							const double & p_type_acceptor_density,
							const double & p_type_width,
							const double & n_type_depletion_width,
							const double & p_type_depletion_width,
							const bool   & schottky_status = false,
							const double & schottky_depletion_width  = 0,
							const double & schottky_hole_density     = 0);


			/** \brief Returns value of \f$\rho_{r}^{\infty}\f$ at point p.*/
			virtual double value(const Point<dim> & p,
					     const unsigned int component=0) const;

		private:
			double scaled_n_type_acceptor_density;
			double scaled_p_type_acceptor_density;
			double scaled_p_type_width;
			double scaled_n_type_depletion_width;
			double scaled_p_type_depletion_width;
			double scaled_schottky_depletion_width;
			double scaled_schottky_hole_density;

			bool schottky;
	};

	template<int dim>
	class LDG_Dirichlet_electron_density_bc : public Function<dim>
	{
		public:
			/** \brief Default constructor. */
			LDG_Dirichlet_electron_density_bc() : Function<dim>(dim+1)
			{}
			virtual ~LDG_Dirichlet_electron_density_bc(){}

			void set_values(const double & n_type_acceptor_density,
							const double & p_type_acceptor_density,
							const double & n_type_donor_density,
							const double & p_type_donor_density,
							const double & p_type_width,
							const double & n_type_width,
							const double & intrinsic_density);


			/** \brief Returns value of \f$\rho_{r}^{\infty}\f$ at point p.*/
			virtual double value(const Point<dim> & p,
					     const unsigned int component=0) const;

		private:
			double scaled_n_type_acceptor_density;
			double scaled_p_type_acceptor_density;
			double scaled_n_type_donor_density;
			double scaled_p_type_donor_density;
			double scaled_p_type_width;
			double scaled_n_type_width;
			double scaled_intrinsic_density;

			double scaled_n_type_electron_bc;
			double scaled_p_type_electron_bc;
	};


	template<int dim>
	class LDG_Dirichlet_hole_density_bc : public Function<dim>
	{
		public:
			/** \brief Default constructor. */
			LDG_Dirichlet_hole_density_bc() : Function<dim>(dim+1)
			{}
			virtual ~LDG_Dirichlet_hole_density_bc(){}

			void set_values(const double & n_type_acceptor_density,
							const double & p_type_acceptor_density,
							const double & n_type_donor_density,
							const double & p_type_donor_density,
							const double & p_type_width,
							const double & n_type_width,
							const double & intrinsic_density);


			/** \brief Returns value of \f$\rho_{r}^{\infty}\f$ at point p.*/
			virtual double value(const Point<dim> & p,
					     const unsigned int component=0) const;

		private:
			double scaled_n_type_acceptor_density;
			double scaled_p_type_acceptor_density;
			double scaled_n_type_donor_density;
			double scaled_p_type_donor_density;
			double scaled_p_type_width;
			double scaled_n_type_width;
			double scaled_intrinsic_density;

			double scaled_n_type_hole_bc;
			double scaled_p_type_hole_bc;
	};


	template<int dim>
	class LDG_Dirichlet_electron_density_bc_gb : public Function<dim>
	{
		public:
			/** \brief Default constructor. */
			LDG_Dirichlet_electron_density_bc_gb() : Function<dim>(dim+1)
			{}
			virtual ~LDG_Dirichlet_electron_density_bc_gb(){}

			void set_values(const double & electron_bc);


			/** \brief Returns value of \f$\rho_{r}^{\infty}\f$ at point p.*/
			virtual double value(const Point<dim> & p,
					     const unsigned int component=0) const;

		private:
			double scaled_electron_bc;
	};


	template<int dim>
	class LDG_Dirichlet_hole_density_bc_gb : public Function<dim>
	{
		public:
			/** \brief Default constructor. */
			LDG_Dirichlet_hole_density_bc_gb() : Function<dim>(dim+1)
			{}
			virtual ~LDG_Dirichlet_hole_density_bc_gb(){}

			void set_values(const double & scaled_hole_bc_precalculated);


			/** \brief Returns value of \f$\rho_{r}^{\infty}\f$ at point p.*/
			virtual double value(const Point<dim> & p,
					     const unsigned int component=0) const;

		private:
			double scaled_hole_bc;
	};

#endif
