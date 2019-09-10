#ifndef _INITIALCONDITIONS_H__
#define _INITIALCONDITIONS_H__

#include<deal.II/base/function.h>

using namespace dealii;

/** \brief Set doping profile for donors (e.g. 0 on p-type side and Nd for n-type side.*/
template<int dim>
class Doping_profile_donors : public Function<dim>
{
	public:
		/** \brief Default constructor. */
		Doping_profile_donors() : Function<dim>(dim+1)
		{};

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
			Initial_condition_electrons() : Function<dim>(dim+1)
			{}

			void set_values(const double & n_type_donor_density,
							const double & p_type_donor_density,
							const double & p_type_width,
							const double & n_type_depletion_width,
							const double & p_type_depletion_width);


			/** \brief Returns value of \f$\rho_{r}^{\infty}\f$ at point p.*/
			virtual double value(const Point<dim> & p,
					     const unsigned int component=0) const;

		private:
			double scaled_n_type_donor_density;
			double scaled_p_type_donor_density;
			double scaled_p_type_width;
			double scaled_n_type_depletion_width;
			double scaled_p_type_depletion_width;
	};



	/** \brief Initial conditions for electrons in full depletion approximation.
	 * That is: zero carrier density near to pn interface.
	 *  \f$ \rho_{o}^{\infty} \f$.*/
	template<int dim>
	class Initial_condition_holes : public Function<dim>
	{
		public:
			/** \brief Default constructor. */
			Initial_condition_holes() : Function<dim>(dim+1)
			{}

			void set_values(const double & n_type_acceptor_density,
							const double & p_type_acceptor_density,
							const double & p_type_width,
							const double & n_type_depletion_width,
							const double & p_type_depletion_width);


			/** \brief Returns value of \f$\rho_{r}^{\infty}\f$ at point p.*/
			virtual double value(const Point<dim> & p,
					     const unsigned int component=0) const;

		private:
			double scaled_n_type_acceptor_density;
			double scaled_p_type_acceptor_density;
			double scaled_p_type_width;
			double scaled_n_type_depletion_width;
			double scaled_p_type_depletion_width;
	};

	template<int dim>
	class LDG_Dirichlet_electron_density_bc : public Function<dim>
	{
		public:
			/** \brief Default constructor. */
			LDG_Dirichlet_electron_density_bc() : Function<dim>(dim+1)
			{}

			void set_values(const double & n_type_donor_density,
							const double & p_type_donor_density,
							const double & p_type_width,
							const double & n_type_width);


			/** \brief Returns value of \f$\rho_{r}^{\infty}\f$ at point p.*/
			virtual double value(const Point<dim> & p,
					     const unsigned int component=0) const;

		private:
			double scaled_n_type_donor_density;
			double scaled_p_type_donor_density;
			double scaled_p_type_width;
			double scaled_n_type_width;
	};


	template<int dim>
	class LDG_Dirichlet_hole_density_bc : public Function<dim>
	{
		public:
			/** \brief Default constructor. */
			LDG_Dirichlet_hole_density_bc() : Function<dim>(dim+1)
			{}

			void set_values(const double & n_type_acceptor_density,
							const double & p_type_acceptor_density,
							const double & p_type_width,
							const double & n_type_width);


			/** \brief Returns value of \f$\rho_{r}^{\infty}\f$ at point p.*/
			virtual double value(const Point<dim> & p,
					     const unsigned int component=0) const;

		private:
			double scaled_n_type_acceptor_density;
			double scaled_p_type_acceptor_density;
			double scaled_p_type_width;
			double scaled_n_type_width;
	};

#endif
