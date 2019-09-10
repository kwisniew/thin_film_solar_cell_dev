#ifndef _PARAMETERS_H__
#define _PARAMETERS_H__

#include<deal.II/base/parameter_handler.h>
#include<iostream>
#include <math.h>

namespace PhysicalConstants
{
	const double electron_charge  	 = 1.6021765e-19;// [C]
	const double boltzman_constant   = 1.3806488e-23;// [J/K]=[kg*m^2/(s^2*K)]
	//permittivity mast be in cm!
	const double vacuum_permittivity = 8.8541878e-14;// [C V^{-1} cm^{-1}] = [A^2 s^4 kg^-1 m^-3]
	const double planck_constant     = 6.63e-34;// [m^{2} kg/s]
	const double free_electron_mass  = 9.11e-31;// [kg]
}

namespace ParameterSpace
{
	using namespace dealii;

	/// \brief Struct which holds the parameters used for simulations.
	
	
	/** These are the parameters which will be used by
	*		DriftDiffusionPoisson class to do simulations. They will be read
	*		read in using the ParameterReader class to read them from 
	*		input_file.prm and then use the ParameterHandler to 
	* 	set the variables and scale them if necessary.
	*/							 
	struct Parameters
	{
			// computational
			unsigned int n_global_refine;				
			unsigned int n_local_refine;				
			unsigned int time_stamps;
			double 	 	 	 h_max;
			double			 h_min;
			double 			 t_end;
			double 			 t_end_2;
			double 			 delta_t;	
			double 			 penalty;

			//mesh
			double scaled_boundary_layer;
			//double scaled_domain_length;
			double scaled_domain_height;
			double scaled_radius_one;
			double scaled_radius_two;
			double scaled_n_type_width;
			double scaled_p_type_width;
			double scaled_n_type_donor_density;
			double scaled_n_type_acceptor_density;
			double scaled_p_type_donor_density;
			double scaled_p_type_acceptor_density;


			//electron
			double scaled_electron_mobility;
			double scaled_electron_recombo_t;
			double scaled_electron_recombo_v;
			double electron_effective_mass;
			double scaled_n_type_depletion_width;
			//double scaled_k_et;


			//hole
			double scaled_hole_mobility;
			double scaled_hole_recombo_t;
			double scaled_hole_recombo_v;
			double hole_effective_mass;
			double scaled_p_type_depletion_width;
			//double scaled_k_ht;


			//Physical
			double scaled_intrinsic_density;
			double semiconductor_permittivity;
			double scaled_absorption_coeff;
			double scaled_photon_flux;
			//note: any time you use scaled_debye_length in a program it need to be multiply by semiconductor_permittivity
			double scaled_debye_length;
			double Nc_effective_dos;
			double Nv_effective_dos;
			double band_gap;
			double defect_energy;
			double temperatue;
			double thermal_voltage;

			double characteristic_length;
			double characteristic_time;
			double characteristic_denisty;


			double scaled_schottky_hole_density;
			double scaled_schottky_electron_density;
			double scaled_srh_electron_density;
			double scaled_srh_hole_density;


			bool illum_or_dark;
			bool insulated;
			bool restart_status;
			bool schottky_status;
	

			double scaled_applied_bias;
			double scaled_built_in_bias;
			double scaled_schottky_bias;

			double rescale_current;

			/** This functions makes almost everythingn 1.0 that is relevant for testing.*/
			void set_params_for_testing(const unsigned int & n_refine)
			{
				n_global_refine 	 = n_refine;
				t_end 			 = 1.0;
				scaled_electron_mobility = 1.0;
				scaled_hole_mobility	 = 1.0;
				scaled_absorption_coeff  = 0.0;
	
				scaled_domain_height  = 1.0;
				scaled_radius_one     = 0.5;
				scaled_radius_two     = 0.5;

				scaled_debye_length   = 1.0;
				characteristic_length  = 1.0;
				characteristic_denisty = 1.0;
				characteristic_time    = 1.0;
			}

			/** This function opens <code>input_file.prm<code/>, 
			*		reads in the specified 
			* 	parameter values and scales appropriate values for singular
			* 	perturbation scaling.
			* 
			*   NOTE: This must be called from the constructor: 
			*		SOLARCELL::SolarCellProblem()
			*/
			void parse_and_scale_parameters(ParameterHandler & prm)
			{

				// read in the parameters
				prm.enter_subsection("computational");
				this->n_global_refine = prm.get_integer("global refinements");
				this->n_local_refine  = prm.get_integer("local refinements");
				this->delta_t  = prm.get_double("time step size");
				this->t_end = prm.get_double("end time");
				this->t_end_2 = prm.get_double("end time 2");
				this->time_stamps = prm.get_integer("time stamps");
				this->restart_status = prm.get_bool("restart status");
				prm.leave_subsection();

				prm.enter_subsection("mesh");
				this->scaled_domain_height = prm.get_double("mesh height");
				this->scaled_radius_one    = prm.get_double("radius one");
				this->scaled_radius_two    = prm.get_double("radius two");
				this->scaled_boundary_layer= prm.get_double("boundary layer");
				this->scaled_n_type_width  = prm.get_double("n_type width");
				this->scaled_p_type_width  = prm.get_double("p_type width");
				prm.leave_subsection();

				prm.enter_subsection("electrons");
				this->scaled_electron_mobility = prm.get_double("mobility");
				this->scaled_electron_recombo_t = prm.get_double("recombination time");
				this->scaled_electron_recombo_v = prm.get_double("recombination velocity");
				this->electron_effective_mass   = prm.get_double("electron effective mass");
				prm.leave_subsection();

				prm.enter_subsection("holes");
				this->scaled_hole_mobility = prm.get_double("mobility");
				this->scaled_hole_recombo_t = prm.get_double("recombination time");
				this->scaled_hole_recombo_v = prm.get_double("recombination velocity");
				this->hole_effective_mass   = prm.get_double("hole effective mass");
				prm.leave_subsection();

				prm.enter_subsection("physical");	
				this->illum_or_dark = prm.get_bool("illumination status");
				this->insulated = prm.get_bool("insulated");
				this->schottky_status = prm.get_bool("schottky status");
				this->scaled_applied_bias = prm.get_double("applied bias");
				this->band_gap = prm.get_double("band gap");
				this->defect_energy = prm.get_double("defect energy level");
				this->temperatue = prm.get_double("temperature");
				this->scaled_schottky_bias = prm.get_double("schottky bias");
				this->characteristic_length = prm.get_double("characteristic length");
				this->characteristic_time = prm.get_double("characteristic time");
				this->scaled_n_type_donor_density     = prm.get_double("n_type donor density");
				this->scaled_n_type_acceptor_density  = prm.get_double("n_type acceptor density");
				this->scaled_p_type_donor_density     = prm.get_double("p_type donor density");
				this->scaled_p_type_acceptor_density  = prm.get_double("p_type acceptor density");
				this->scaled_photon_flux = prm.get_double("photon flux");
				this->scaled_absorption_coeff = prm.get_double("absorption coefficient"); 
				this->semiconductor_permittivity = prm.get_double("semiconductor permittivity");
				prm.leave_subsection();



/*				prm.enter_subsection("reductants");
				this->scaled_reductant_mobility = prm.get_double("mobility");
				prm.leave_subsection();*/


/*				prm.enter_subsection("oxidants");
				this->scaled_oxidant_mobility = prm.get_double("mobility");
				prm.leave_subsection();*/

				//compute parameters:


				this->thermal_voltage = PhysicalConstants::boltzman_constant*temperatue/PhysicalConstants::electron_charge;
				std::cout << "thermal voltage:    " << this->thermal_voltage << std::endl;

				//calculate intrinsic density:
				this->Nc_effective_dos = 2*pow(2*M_PI*PhysicalConstants::free_electron_mass*this->electron_effective_mass
												*PhysicalConstants::boltzman_constant*this->temperatue
												/(PhysicalConstants::planck_constant*PhysicalConstants::planck_constant)
											  ,1.5);

				std::cout << "Nc [cm^-3]:    " << this->Nc_effective_dos*1.0e-6 << std::endl;

				this->Nv_effective_dos = 2*pow(2*M_PI*PhysicalConstants::free_electron_mass*this->hole_effective_mass
												*PhysicalConstants::boltzman_constant*this->temperatue
												/(PhysicalConstants::planck_constant*PhysicalConstants::planck_constant)
											  ,1.5);

				std::cout << "Nv [cm^-3]:    " << this->Nv_effective_dos*1.0e-6 << std::endl;

				//band gap is multiplied by electron charge to convert eV [eV] to Joules [J].
				this->scaled_intrinsic_density = sqrt(this->Nv_effective_dos*this->Nc_effective_dos)
												*std::exp(-this->band_gap*PhysicalConstants::electron_charge
														  /(2*PhysicalConstants::boltzman_constant*this->temperatue));



				//change intrinsic density to [cm^-3]
				this->scaled_intrinsic_density *= 1.0e-6;

				std::cout << "Przeskalowana ni do centymetrów:   "
						  << this->scaled_intrinsic_density
						  << "cm^-3"
						  << std::endl;

				//calculate build in bias (not scaled yet):
				//Note: all densities (n_type_doping, p_type_doping and intrinsic_density)
				//      need to be unscaled and in [cm^-3]
				this->scaled_built_in_bias =  this->thermal_voltage
											 *log(this->scaled_n_type_donor_density*this->scaled_p_type_acceptor_density
												  /(this->scaled_intrinsic_density*this->scaled_intrinsic_density));

				std::cout << "Potencjał wbudowany:   "
						  << this->scaled_built_in_bias
						  << "[V]"
						  << std::endl;

				//calculate carrier densities on schottky contact
				//Note: all biases need to be without scaling, that is in volts: [V]
				this->scaled_schottky_electron_density = this->Nc_effective_dos
														*std::exp( -(this->scaled_built_in_bias - this->scaled_schottky_bias)
																	/this->thermal_voltage);

				this->scaled_schottky_hole_density     = this->Nv_effective_dos
														 *std::exp( -(this->scaled_schottky_bias)
																 	 /this->thermal_voltage);
				//change Schottky densities to [cm^-3]
				this->scaled_schottky_electron_density *= 1e-6;
				this->scaled_schottky_hole_density     *= 1e-6;


				//Initial condition depends on depletion width,
				//in case when Vapp > Vbi other initial value need to be impose eg. very low free carrier density at the beginning
				// depletion width will be in [cm] since dopings and permittivity ~ [cm]
				this->scaled_n_type_depletion_width = sqrt(2*PhysicalConstants::vacuum_permittivity*this->semiconductor_permittivity
															*this->scaled_p_type_acceptor_density*(this->scaled_built_in_bias - this->scaled_applied_bias)
															/
															(PhysicalConstants::electron_charge*this->scaled_n_type_donor_density
															* (this->scaled_n_type_donor_density + this->scaled_p_type_acceptor_density)));

				this->scaled_p_type_depletion_width = sqrt(2*PhysicalConstants::vacuum_permittivity*this->semiconductor_permittivity
															*this->scaled_n_type_donor_density*(this->scaled_built_in_bias - this->scaled_applied_bias)
															/
															(PhysicalConstants::electron_charge*this->scaled_p_type_acceptor_density
															* (this->scaled_n_type_donor_density + this->scaled_p_type_acceptor_density)));

				std::cout << "scaled_n_type_depletion_width:   "
						  << this->scaled_n_type_depletion_width*10000
						  << "um"
						  << std::endl
						  << "scaled_p_type_depletion_width:   "
						  << this->scaled_p_type_depletion_width*10000
						  << "um"
						  << std::endl;

				//defect energy in [eV]:
				this->defect_energy *= this->band_gap;
				//defect energy is eV, so we need to change them to [J] by multiplying it by e
				// see eg. p.105 Selberherr  (srh_electron_density=n1)
				this->scaled_srh_electron_density = Nc_effective_dos
											*std::exp(-this->defect_energy*PhysicalConstants::electron_charge
													  /(PhysicalConstants::boltzman_constant*this->temperatue));

				this->scaled_srh_hole_density = Nv_effective_dos
										*std::exp( (this->defect_energy - this->band_gap)*PhysicalConstants::electron_charge
												  /(PhysicalConstants::boltzman_constant*this->temperatue));

				std::cout << "srh_electron_density:   "
						  << this->scaled_srh_electron_density
						  << "cm^-3"
						  << std::endl
						  << "srh_hole_density:       "
						  << this->scaled_srh_hole_density
						  << "cm^-3"
						  << std::endl;

				// scale the parameters
				this->characteristic_denisty =
				// MAX(a, b) == ((a > b) ? a : b)
				(this->scaled_p_type_acceptor_density < this->scaled_n_type_donor_density ? this->scaled_p_type_acceptor_density : this->scaled_n_type_donor_density);

				this->scaled_n_type_donor_density      /= this->characteristic_denisty;
				this->scaled_n_type_acceptor_density   /= this->characteristic_denisty;
				this->scaled_p_type_donor_density      /= this->characteristic_denisty;
				this->scaled_p_type_acceptor_density   /= this->characteristic_denisty;
				this->scaled_intrinsic_density         /= this->characteristic_denisty;
				this->scaled_schottky_electron_density /= this->characteristic_denisty;
				this->scaled_schottky_hole_density     /= this->characteristic_denisty;
				this->scaled_srh_electron_density      /= this->characteristic_denisty;
				this->scaled_srh_hole_density		   /= this->characteristic_denisty;

				//CHRACTERISTIC LENGTH MUST BE IN CM!
				this->scaled_n_type_depletion_width /= (this->characteristic_length);
				this->scaled_p_type_depletion_width /= (this->characteristic_length);


				this->scaled_electron_recombo_t /= this->characteristic_time;
				this->scaled_hole_recombo_t /= this->characteristic_time;

				this->scaled_electron_recombo_v *= (this->characteristic_time /
								this->characteristic_length);

				this->scaled_hole_recombo_v *= (this->characteristic_time /
								this->characteristic_length);

				this->scaled_photon_flux *= (this->characteristic_time /
							     this->characteristic_denisty);

				this->scaled_absorption_coeff *= this->characteristic_length;

				// DONT INCLUDE THE MATERIAL DIELECTRIC IN THE DEBYE LENGTH
				this->scaled_debye_length = 
						(this->thermal_voltage *
						PhysicalConstants::vacuum_permittivity) /
//						this->material_permittivity) /
						(PhysicalConstants::electron_charge * 
						this->characteristic_denisty * 
						this->characteristic_length *
						this->characteristic_length);

				std::cout << "Debye length (dimensonless):   "
						  << this->scaled_debye_length
						  << std::endl;

				double mobility_scale = (this->characteristic_time *
						this->thermal_voltage) /
						(this->characteristic_length *
						this->characteristic_length);

				this->scaled_electron_mobility *= mobility_scale;
				this->scaled_hole_mobility *= mobility_scale;

				this->scaled_applied_bias  /= this->thermal_voltage;
				this->scaled_built_in_bias /= this->thermal_voltage;
				this->scaled_schottky_bias /= this->thermal_voltage;

				this->rescale_current = (PhysicalConstants::electron_charge 
							* this->characteristic_denisty 
							* this->characteristic_length)
							/ this->characteristic_time;


/*
				std::cout << "debeye length = " 
					<< this->scaled_debeye_length 
					<< std::endl;
				std::cout << "semiconductor perm = "
					<< this->semiconductor_permittivity
					<< std::endl;
				std::cout << "scaled electron mobility = " 
					<< this->scaled_electron_mobility 
					<< std::endl;
				std::cout << "scaled hole mobility = " 
					<< this->scaled_hole_mobility 
					<< std::endl;
				std::cout << "built in bias = " 
					<< this->scaled_built_in_bias 
					<< std::endl;
				std::cout << "rescale current = "
					<< this->rescale_current
					<< std::endl;
				std::cout << "photon flux = " 
					<< this->scaled_photon_flux
					<< std::endl;
				std::cout << "absorption coeff = "
					<< this->scaled_absorption_coeff 
					<< std::endl;
*/

		} // parse_and_scale_parameters(prm)
		
	};

	/// @author Michael Harmon
}


#endif
