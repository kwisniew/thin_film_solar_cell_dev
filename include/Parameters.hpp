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
	const double richardson_constant = 4*PhysicalConstants::electron_charge*M_PI*PhysicalConstants::free_electron_mass
								  	    *PhysicalConstants::boltzman_constant*PhysicalConstants::boltzman_constant
										/(PhysicalConstants::planck_constant*PhysicalConstants::planck_constant*PhysicalConstants::planck_constant);
										//C s^-1 m^-2 K^-2 = A m^-2 K^-2
}

namespace ParameterSpace
{
	using namespace dealii;

	/// \brief Struct which holds the parameters used for simulations.
	
	
	/** These are the parameters which will be used by
	*		SolarCell class to do simulations. They will be read
	*		read in using the ParameterReader class to read them from 
	*		input_file.prm and then use the ParameterHandler to 
	* 	set the variables and scale them if necessary.
	*/							 
	struct Parameters
	{
			// computational
			bool 		 calculate_steady_state;
			bool 		 calculate_equilibrium_state;
			bool 		 calculate_DLTS;
			bool		 calculate_IV_curve;
			bool		 calculate_CV_curve;
			bool		 join_pn_gb;
			std::string  type_of_simulation;
			unsigned int n_global_refine;				
			unsigned int n_local_refine;				
			unsigned int time_stamps;
			double 	 	 h_max;
			double		 h_min;
			double		 t_end;
			double 		 t_end_steady_state;
			double 		 t_end_equilibrium_state;
			double 		 t_end_DLTS;
			double 		 t_end_2;
			double 		 delta_t;
			double 		 penalty;
			bool 		 restart_status;
			bool		 restart_from_steady_state;
			std::string  type_of_restart;

			//mesh
			double scaled_grain_boundary_width;
			//double scaled_domain_length;
			double scaled_domain_height;
			double scaled_top_point_x;
			double scaled_top_point_y;
			double scaled_bottom_point_x;
			double scaled_bottom_point_y;
			double scaled_n_type_width;
			double scaled_p_type_width;
			double scaled_estimated_gb_depletion_width;


			//electron
			double scaled_electron_mobility;
			double scaled_electron_recombo_t;
			double scaled_electron_recombo_t_gb;
			double scaled_electron_schottky_recombo_v;
			double electron_effective_mass;
			double scaled_n_type_depletion_width;
			//double scaled_k_et;


			//hole
			double scaled_hole_mobility;
			double scaled_hole_recombo_t;
			double scaled_hole_recombo_t_gb;
			double scaled_hole_schottky_recombo_v;
			double hole_effective_mass;
			double scaled_p_type_depletion_width;
			double scaled_p_type_schottky_dwidth;

			//for calculating outputs currents in postprocesor, should be deleted in near future
			double mobility;

			//Physical
			double real_domain_height;
			double device_thickness;
			double scaled_intrinsic_density;
			double semiconductor_permittivity;
			double scaled_absorption_coeff;
			double scaled_photon_flux;
			//NOTE: any time you use scaled_debye_length in a program it need to be multiply by semiconductor_permittivity
			double scaled_debye_length;
			double Nc_effective_dos;
			double Nv_effective_dos;
			double band_gap;
			double defect_energy;
			double defect_energy_gb;
			double temperature;
			double thermal_voltage;

			double scaled_n_type_donor_density;
			double scaled_n_type_acceptor_density;
			double scaled_p_type_donor_density;
			double scaled_p_type_acceptor_density;

			double characteristic_length;
			double characteristic_denisty;
			double characteristic_time;
			double characteristic_time_steady_state;
			double characteristic_time_equilibrium_state;
			double characteristic_time_DLTS;

			double scaled_schottky_hole_density;
			double scaled_schottky_electron_density;
			double scaled_srh_electron_density;
			double scaled_srh_hole_density;
			double scaled_gb_defect_density;
			double scaled_gb_electron_density_Dirichlet;
			double scaled_gb_hole_density_Dirichlet;
			double scaled_srh_electron_density_gb;
			double scaled_srh_hole_density_gb;


			bool illum_or_dark;
			bool insulated;
			bool schottky_status;
			bool horizontal_gb;
			bool vertical_gb;
			bool grain_boundary_status;
	

			double scaled_applied_bias;
			double scaled_built_in_bias;
			double scaled_schottky_bias;

			double rescale_current;

			//IV
			double scaled_IV_min_V;
			double scaled_IV_max_V;
			unsigned int IV_n_of_data_points;
			double scaled_delta_V;
			double charge_scaling_factor;


			/** This functions makes almost everythingn 1.0 that is relevant for testing.*/
			void set_params_for_testing(const unsigned int & n_refine)
			{
				n_global_refine 	 = n_refine;
				t_end 			 = 1.0;
				scaled_electron_mobility = 1.0;
				scaled_hole_mobility	 = 1.0;
				scaled_absorption_coeff  = 0.0;
	
				scaled_domain_height  = 1.0;
				scaled_top_point_x     = 0.5;
				scaled_bottom_point_x     = 0.5;

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
				calculate_steady_state          = prm.get_bool("steady state");
				calculate_equilibrium_state     = prm.get_bool("equilibrium state");
				calculate_DLTS 				  = prm.get_bool("DLTS");
				calculate_IV_curve			  = prm.get_bool("IV curve");
				calculate_CV_curve			  = prm.get_bool("CV curve");
				join_pn_gb							  = prm.get_bool("join grain boundary and junction");

				n_global_refine = prm.get_integer("global refinements");
				n_local_refine  = prm.get_integer("local refinements");
				delta_t         = prm.get_double("time step size");

				t_end_steady_state          = prm.get_double("end time steady state");
				t_end_equilibrium_state     = prm.get_double("end time equilibrium state");
				t_end_DLTS = prm.get_double("end time DLTS");

				t_end_2 = prm.get_double("end time 2");
				time_stamps = prm.get_integer("time stamps");
				restart_status            = prm.get_bool("restart status");
				restart_from_steady_state = prm.get_bool("restart from steady state");
				prm.leave_subsection();

				prm.enter_subsection("mesh");
				scaled_domain_height = prm.get_double("mesh height");
				scaled_top_point_x    = prm.get_double("top point x");
				scaled_top_point_y    = prm.get_double("top point y");
				scaled_bottom_point_x    = prm.get_double("bottom point x");
				scaled_bottom_point_y    = prm.get_double("bottom point y");
				scaled_grain_boundary_width= prm.get_double("grain boundary width");
				scaled_n_type_width  = prm.get_double("n_type width");
				scaled_p_type_width  = prm.get_double("p_type width");
				scaled_estimated_gb_depletion_width = prm.get_double("grain boundary depletion width");
				prm.leave_subsection();

				prm.enter_subsection("electrons");
				scaled_electron_mobility = prm.get_double("mobility");
				scaled_electron_recombo_t = prm.get_double("recombination time");
				scaled_electron_recombo_t_gb = prm.get_double("grain boundary recombination time");
				//scaled_electron_recombo_v = prm.get_double("recombination velocity");
				electron_effective_mass   = prm.get_double("electron effective mass");
				prm.leave_subsection();

				prm.enter_subsection("holes");
				scaled_hole_mobility = prm.get_double("mobility");
				scaled_hole_recombo_t = prm.get_double("recombination time");
				scaled_hole_recombo_t_gb = prm.get_double("grain boundary recombination time");
				//scaled_hole_recombo_v = prm.get_double("recombination velocity");
				hole_effective_mass   = prm.get_double("hole effective mass");
				prm.leave_subsection();

				prm.enter_subsection("physical");
				real_domain_height = prm.get_double("real domain height");
				device_thickness   = prm.get_double("device thickness");
				illum_or_dark = prm.get_bool("illumination status");
				insulated = prm.get_bool("insulated");
				schottky_status = prm.get_bool("schottky status");
				grain_boundary_status = prm.get_bool("grain boundary status");
				vertical_gb   = prm.get_bool("vertical grain boundary");
				horizontal_gb = prm.get_bool("horizontal grain boundary");
				scaled_applied_bias = prm.get_double("applied bias");
				band_gap = prm.get_double("band gap");
				defect_energy = prm.get_double("defect energy level");
				defect_energy_gb = prm.get_double("grain boundary defect energy level");
				temperature = prm.get_double("temperature");
				scaled_schottky_bias = prm.get_double("schottky bias");
				characteristic_length = prm.get_double("characteristic length");

				characteristic_time_steady_state          = prm.get_double("characteristic time steady state");
				characteristic_time_equilibrium_state     = prm.get_double("characteristic time equilibrium state");
				characteristic_time_DLTS = prm.get_double("characteristic time DLTS");

				scaled_n_type_donor_density     = prm.get_double("n_type donor density");
				scaled_n_type_acceptor_density  = prm.get_double("n_type acceptor density");
				scaled_p_type_donor_density     = prm.get_double("p_type donor density");
				scaled_p_type_acceptor_density  = prm.get_double("p_type acceptor density");
				scaled_gb_defect_density = prm.get_double("grain boundary defect density");
				scaled_gb_electron_density_Dirichlet  = prm.get_double("gb electron density on dirichlet");
				scaled_gb_hole_density_Dirichlet      = prm.get_double("gb hole density on dirichlet");
				scaled_photon_flux = prm.get_double("photon flux");
				scaled_absorption_coeff = prm.get_double("absorption coefficient");
				semiconductor_permittivity = prm.get_double("semiconductor permittivity");
				prm.leave_subsection();

				prm.enter_subsection("IV");
				scaled_IV_min_V = prm.get_double("IV minimal voltage");
				scaled_IV_max_V = prm.get_double("IV maximal voltage");
				IV_n_of_data_points = prm.get_double("IV number of data points");
				scaled_delta_V        = prm.get_double("delta V");
				prm.leave_subsection();


				std::cout << "SZEROKOSC GRANICY ZIAREN:   " << scaled_grain_boundary_width << "\n";
				//TODO: NEED TO CHANGE IT AS SOON AS POSSIBLE
				//it is for scaling of outputs values
				mobility = scaled_electron_mobility;

				/*----------------------------------------------------------------------*/
				// 	Compute parameters
				/*----------------------------------------------------------------------*/

				std::cout <<std::endl;
				if(calculate_steady_state  && !calculate_equilibrium_state && !calculate_DLTS)
				{
					std::cout << "STEADY STATE CALCULATION" <<std::endl;

					type_of_simulation = "_steady_state";
					characteristic_time = characteristic_time_steady_state;
					t_end               = t_end_steady_state;
					scaled_applied_bias = 0;
				}
				else if(!calculate_steady_state  && calculate_equilibrium_state && !calculate_DLTS)
				{
					std::cout << "EQUILIBRIUM STATE CALCULATION" <<std::endl;

					type_of_simulation = "";
					characteristic_time = characteristic_time_equilibrium_state;
					t_end               = t_end_equilibrium_state;
				}
				else if(!calculate_steady_state  && !calculate_equilibrium_state && calculate_DLTS)
				{
					std::cout << "DLTS SIGNAL CALCULATION" <<std::endl;

					type_of_simulation = "";
					characteristic_time = characteristic_time_DLTS;
					t_end               = t_end_DLTS;
				}
				else
				{
					std::cout << "There need to be one and ONLY ONE type of SIMULATION AT THE SAME TIME!!! \n The steady state time will be used!" <<std::endl;
					type_of_simulation = "steady_state";
					characteristic_time = characteristic_time_steady_state;
					t_end               = t_end_steady_state;
				}

				if(restart_from_steady_state)
				{
					type_of_restart = "_steady_state";
				}
				else   type_of_restart = "";

				if(scaled_gb_defect_density == 0)
				{
					scaled_electron_recombo_t_gb = scaled_electron_recombo_t;
					scaled_hole_recombo_t_gb     = scaled_hole_recombo_t;
				}
				if(!grain_boundary_status)
				{
					vertical_gb= false;
					horizontal_gb=false;
				}

				std::cout 	<< std::endl
							<< "Parameters:    "
							<< std::endl;


				thermal_voltage = PhysicalConstants::boltzman_constant*temperature/PhysicalConstants::electron_charge;
				std::cout << "thermal voltage:    " << thermal_voltage << std::endl;

				//calculate intrinsic density [m^3]:
				Nc_effective_dos = 2*pow(2*M_PI*PhysicalConstants::free_electron_mass*electron_effective_mass
												*PhysicalConstants::boltzman_constant*temperature
												/(PhysicalConstants::planck_constant*PhysicalConstants::planck_constant)
											  ,1.5);

				std::cout << "Nc [cm^-3]:    " << Nc_effective_dos*1.0e-6 << std::endl;
				//[m^3]
				Nv_effective_dos = 2*pow(2*M_PI*PhysicalConstants::free_electron_mass*hole_effective_mass
												*PhysicalConstants::boltzman_constant*temperature
												/(PhysicalConstants::planck_constant*PhysicalConstants::planck_constant)
											  ,1.5);

				std::cout << "Nv [cm^-3]:    " << Nv_effective_dos*1.0e-6 << std::endl;

				//band gap is multiplied by electron charge to convert eV [eV] to Joules [J].
				//[m^3]
				scaled_intrinsic_density = sqrt(Nv_effective_dos*Nc_effective_dos)
												*std::exp(-band_gap*PhysicalConstants::electron_charge
														  /(2*PhysicalConstants::boltzman_constant*temperature));

				//change intrinsic density to [cm^-3]
				scaled_intrinsic_density *= 1.0e-6;

				std::cout << "Przeskalowana ni do centymetrów:   "
						  << scaled_intrinsic_density
						  << "cm^-3"
						  << std::endl;

				//calculate build in bias (not scaled yet):
				//Note: all densities (n_type_doping, p_type_doping and intrinsic_density)
				//      need to be "unscaled" and in [cm^-3]
				if(scaled_n_type_donor_density > 0 && scaled_p_type_acceptor_density > 0)
				{
					scaled_built_in_bias =  thermal_voltage
												 *log(scaled_n_type_donor_density*scaled_p_type_acceptor_density
													  /(scaled_intrinsic_density*scaled_intrinsic_density));
				}
				/*else if(schottky_status==true)
					scaled_built_in_bias = -scaled_schottky_bias;*/
				else
					scaled_built_in_bias = 0;


				std::cout << "Potencjał wbudowany:   "
						  << scaled_built_in_bias
						  << "[V]"
						  << std::endl;


				//const double n_type_resultant_doping = (scaled_n_type_donor_density - scaled_n_type_acceptor_density);
				const double p_type_resultant_doping = (scaled_p_type_donor_density - scaled_p_type_acceptor_density);

				scaled_schottky_hole_density     = 0.5*(-p_type_resultant_doping + std::sqrt(p_type_resultant_doping*p_type_resultant_doping + 4*scaled_intrinsic_density*scaled_intrinsic_density));
				scaled_schottky_electron_density = 0.5*( p_type_resultant_doping + std::sqrt(p_type_resultant_doping*p_type_resultant_doping + 4*scaled_intrinsic_density*scaled_intrinsic_density));

				scaled_schottky_hole_density    *= std::exp( -(scaled_schottky_bias)
					 	 	 	 	 	 	 	 	 	 	  /thermal_voltage);

				scaled_schottky_electron_density*= std::exp( (scaled_schottky_bias)
									 	 	 	 	 	 	 /thermal_voltage);

				std::cout << "schottky_electron_density:   "
						  << scaled_schottky_electron_density
						  << "cm^-3"
						  << std::endl
						  << "schottky_hole_density:       "
						  << scaled_schottky_hole_density
						  << "cm^-3"
						  << std::endl;


				//Initial condition depends on depletion width,
				//in case when Vapp > Vbi other initial value need to be impose eg. very low free carrier density at the beginning
				// depletion width will be in [cm] since dopings and permittivity ~ [cm]
				if(scaled_n_type_donor_density > 0)
				{
					scaled_n_type_depletion_width = sqrt(2*PhysicalConstants::vacuum_permittivity*semiconductor_permittivity
																*scaled_p_type_acceptor_density*(scaled_built_in_bias - scaled_applied_bias)
																/
																(PhysicalConstants::electron_charge*scaled_n_type_donor_density
																* (scaled_n_type_donor_density + scaled_p_type_acceptor_density)));
				}
				else
					scaled_n_type_depletion_width = 0.0;

				if(scaled_p_type_acceptor_density > 0)
				{
					scaled_p_type_depletion_width = sqrt(2*PhysicalConstants::vacuum_permittivity*semiconductor_permittivity
																*scaled_n_type_donor_density*(scaled_built_in_bias - scaled_applied_bias)
																/
																(PhysicalConstants::electron_charge*scaled_p_type_acceptor_density
																* (scaled_n_type_donor_density + scaled_p_type_acceptor_density)));
				}
				else
					scaled_p_type_depletion_width = 0.0;

				scaled_p_type_schottky_dwidth = sqrt(2*PhysicalConstants::vacuum_permittivity*semiconductor_permittivity
															*(scaled_schottky_bias)
															/
															(PhysicalConstants::electron_charge*scaled_p_type_acceptor_density));

				std::cout << "scaled_n_type_depletion_width:   "
						  << scaled_n_type_depletion_width*10000
						  << "um"
						  << std::endl
						  << "scaled_p_type_depletion_width:   "
						  << scaled_p_type_depletion_width*10000
						  << "um"
						  << std::endl
						  << "scaled_p_type_schottky_dwidth:   "
						  << scaled_p_type_schottky_dwidth*10000
						  << "um"
						  << std::endl;

				//defect energy in [eV]:
				defect_energy *= band_gap;
				//defect energy on grain boundary
				defect_energy_gb *= band_gap;

				//defect energy is eV, so we need to change them to [J] by multiplying it by e
				// see eg. p.105 Selberherr  (srh_electron_density=n1)
				//[m^3]
				scaled_srh_electron_density = Nc_effective_dos
											*std::exp(-defect_energy*PhysicalConstants::electron_charge
													  /(PhysicalConstants::boltzman_constant*temperature));
				//[cm^3]
				scaled_srh_electron_density *= 1.0e-6;

				//[m^3]
				scaled_srh_hole_density = Nv_effective_dos
										*std::exp( (defect_energy - band_gap)*PhysicalConstants::electron_charge
												  /(PhysicalConstants::boltzman_constant*temperature));
				//[cm^3]
				scaled_srh_hole_density *= 1.0e-6;



				// electron and hole density as if fermi level was equal defects energy levels on grain boundary
				// the same as for srh_hole/electron_density
				//[m^3]
				scaled_srh_electron_density_gb = Nc_effective_dos
											*std::exp(-defect_energy_gb*PhysicalConstants::electron_charge
													  /(PhysicalConstants::boltzman_constant*temperature));
				//[cm^3]
				scaled_srh_electron_density_gb *= 1.0e-6;

				//[m^3]
				scaled_srh_hole_density_gb = Nv_effective_dos
										*std::exp( (defect_energy_gb - band_gap)*PhysicalConstants::electron_charge
												  /(PhysicalConstants::boltzman_constant*temperature));
				//[cm^3]
				scaled_srh_hole_density_gb *= 1.0e-6;

				if(scaled_gb_defect_density == 0)
				{
					scaled_srh_electron_density_gb = scaled_srh_electron_density;
					scaled_srh_hole_density_gb     = scaled_srh_hole_density;
				}


				std::cout << "srh_electron_density:   "
						  << scaled_srh_electron_density
						  << "cm^-3"
						  << std::endl
						  << "srh_hole_density:       "
						  << scaled_srh_hole_density
						  << "cm^-3"
						  << std::endl
						  << "srh_grain_boundary_elec_density:   "
						  << scaled_srh_electron_density_gb
						  << "cm^-3"
						  << std::endl
						  << "srh_grain_boundary_hole_density:   "
						  << scaled_srh_hole_density_gb
						  << "cm^-3"
						  << std::endl;

				std::cout << "Richardson constant:   "
						  << PhysicalConstants::richardson_constant/1e4
						  << std::endl;
				//[m/s]
				scaled_hole_schottky_recombo_v     = PhysicalConstants::richardson_constant*temperature*temperature
													/(PhysicalConstants::electron_charge*Nv_effective_dos);
				scaled_electron_schottky_recombo_v = PhysicalConstants::richardson_constant*temperature*temperature
													/(PhysicalConstants::electron_charge*Nv_effective_dos);

				//[cm/s]
				scaled_hole_schottky_recombo_v     *=1e2;
				scaled_electron_schottky_recombo_v *=1e2;

				std::cout << "Recombination velocity on schottky[cm/s]:   "
						  << scaled_electron_schottky_recombo_v
						  << std::endl;


				/*----------------------------------------------------------------------*/
				// 	Scale parameters
				/*----------------------------------------------------------------------*/


				characteristic_denisty =
				// MAX(a, b) == ((a > b) ? a : b)
				(scaled_p_type_acceptor_density < scaled_n_type_donor_density ? scaled_p_type_acceptor_density : scaled_n_type_donor_density);
				if(characteristic_denisty == 0)
					characteristic_denisty = scaled_p_type_acceptor_density;

				scaled_n_type_donor_density      /= characteristic_denisty;
				scaled_n_type_acceptor_density   /= characteristic_denisty;
				scaled_p_type_donor_density      /= characteristic_denisty;
				scaled_p_type_acceptor_density   /= characteristic_denisty;
				scaled_intrinsic_density         /= characteristic_denisty;
				scaled_gb_defect_density         /= characteristic_denisty;
				scaled_schottky_electron_density /= characteristic_denisty;
				scaled_schottky_hole_density     /= characteristic_denisty;
				scaled_srh_electron_density      /= characteristic_denisty;
				scaled_srh_electron_density_gb   /= characteristic_denisty;
				scaled_srh_hole_density		     /= characteristic_denisty;
				scaled_srh_hole_density_gb	   	 /= characteristic_denisty;
				scaled_gb_electron_density_Dirichlet /= characteristic_denisty;
				scaled_gb_hole_density_Dirichlet     /= characteristic_denisty;

				//CHRACTERISTIC LENGTH MUST BE IN CM!
				scaled_n_type_depletion_width /= (characteristic_length);
				scaled_p_type_depletion_width /= (characteristic_length);
				scaled_p_type_schottky_dwidth /= (characteristic_length);


				scaled_electron_recombo_t    /= characteristic_time;
				scaled_electron_recombo_t_gb /= characteristic_time;
				scaled_hole_recombo_t        /= characteristic_time;
				scaled_hole_recombo_t_gb     /= characteristic_time;

				scaled_electron_schottky_recombo_v *= (characteristic_time /
								characteristic_length);

				scaled_hole_schottky_recombo_v *= (characteristic_time /
								characteristic_length);

				scaled_photon_flux *= (characteristic_time /
							     characteristic_denisty);

				scaled_absorption_coeff *= characteristic_length;

				// DONT INCLUDE THE MATERIAL DIELECTRIC IN THE DEBYE LENGTH
				scaled_debye_length =
						(thermal_voltage *
						PhysicalConstants::vacuum_permittivity) /
//						material_permittivity) /
						(PhysicalConstants::electron_charge * 
						characteristic_denisty *
						characteristic_length *
						characteristic_length);

				std::cout << "Debye length (dimensonless):   "
						  << scaled_debye_length
						  << std::endl;

				double mobility_scale = (characteristic_time *
						thermal_voltage) /
						(characteristic_length *
						characteristic_length);

				scaled_electron_mobility *= mobility_scale;
				scaled_hole_mobility *= mobility_scale;

				scaled_applied_bias  /= thermal_voltage;
				scaled_built_in_bias /= thermal_voltage;
				scaled_schottky_bias /= thermal_voltage;
				scaled_IV_max_V      /= thermal_voltage;
				scaled_IV_min_V      /= thermal_voltage;
				scaled_delta_V              /= thermal_voltage;


				rescale_current = (PhysicalConstants::electron_charge
							* characteristic_denisty
							* characteristic_length)
							/ characteristic_time;

				charge_scaling_factor =   real_domain_height*device_thickness*characteristic_length*
												characteristic_denisty*PhysicalConstants::electron_charge
											   /scaled_domain_height;



				std::cout 	<< std::endl
							<< std::endl;
/*
				std::cout << "debeye length = " 
					<< scaled_debeye_length
					<< std::endl;
				std::cout << "semiconductor perm = "
					<< semiconductor_permittivity
					<< std::endl;
				std::cout << "scaled electron mobility = " 
					<< scaled_electron_mobility
					<< std::endl;
				std::cout << "scaled hole mobility = " 
					<< scaled_hole_mobility
					<< std::endl;
				std::cout << "built in bias = " 
					<< scaled_built_in_bias
					<< std::endl;
				std::cout << "rescale current = "
					<< rescale_current
					<< std::endl;
				std::cout << "photon flux = " 
					<< scaled_photon_flux
					<< std::endl;
				std::cout << "absorption coeff = "
					<< scaled_absorption_coeff
					<< std::endl;
*/

		} // parse_and_scale_parameters(prm)
		
	};

	/// @author Michael Harmon & Konrad Wisniewski
}


#endif
