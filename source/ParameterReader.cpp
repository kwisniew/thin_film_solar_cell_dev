#include "../include/ParameterReader.hpp"

namespace ParameterSpace
{
	using namespace dealii;
	
	ParameterReader::
	ParameterReader(ParameterHandler & param_handler) 
	:
	prm(param_handler)
	{}

	void
	ParameterReader:: 
	read_parameters(const std::string parameter_file)
	{
		declare_parameters();
		prm.parse_input(parameter_file);
	}

	void
	ParameterReader:: 
	read_test_parameters(const std::string parameter_file)
	{
		declare_test_parameters();
		prm.parse_input(parameter_file);
	}

	void 
	ParameterReader::
	declare_parameters()
	{
		prm.enter_subsection("computational");

		prm.declare_entry("steady state", "true",
				Patterns::Bool(),
				"whether to calculate a steady state situation");

		prm.declare_entry("equilibrium state", "false",
				Patterns::Bool(),
				"whether to calculate a equilibrium state (eg. with applied voltage)");

		prm.declare_entry("DLTS", "false",
				Patterns::Bool(),
				"whether to calculate a DLTS (transient)");

		prm.declare_entry("IV curve", "false",
				Patterns::Bool(),
				"whether to calculate an IV curve");

		prm.declare_entry("CV curve", "false",
				Patterns::Bool(),
				"whether to calculate an CV curve");

		prm.declare_entry("join grain boundary and junction", "false",
				Patterns::Bool(),
				"whether to join solutions from pn and grain boundary calculations");

		prm.declare_entry("global refinements", "4",
				 Patterns::Integer(0,10),
				 "number of global refinements");

		prm.declare_entry("local refinements", "0",
				 Patterns::Integer(0,10),
				"number of local refinements in predefined critical areas");

		prm.declare_entry("time step size", "0.05",
				Patterns::Double(0,1),
				"scaled time step size of both subdomains");

		prm.declare_entry("end time steady state", "1e3",
				Patterns::Double(0),
				"time to steady state (in terms of characteristic time)");

		prm.declare_entry("end time equilibrium state", "1e3",
				Patterns::Double(0),
				"time to equilibrium state (in terms of characteristic time)");

		prm.declare_entry("end time DLTS", "1e3",
				Patterns::Double(0),
				"time to DLTS equilibrium state (in terms of characteristic time)");

		prm.declare_entry("end time 2", "20",
				Patterns::Double(0),
				"2nd end time for restarting at end time");

		prm.declare_entry("time stamps", "50",
				Patterns::Integer(0,1000),
				"number of output files until end time");

		prm.declare_entry("restart status", "false",
				Patterns::Bool(),
				"whether to use end of last run as start for this one");

		prm.declare_entry("restart from steady state", "false",
				Patterns::Bool(),
				"whether to use end of last steady state simulation for this one");

		prm.leave_subsection();				




		prm.enter_subsection("mesh");

		prm.declare_entry("mesh height", "1",
				Patterns::Double(0),
				"scaled domain height in y coordinate");

		prm.declare_entry("top point x", "0.5",
				Patterns::Double(0),
				"top radius of wire");

		prm.declare_entry("top point y", "5.0",
				Patterns::Double(0),
				"top radius of wire");

		prm.declare_entry("bottom point x", "1.0",
				Patterns::Double(0),
				"bottom radius of wire");

		prm.declare_entry("bottom point y", "0.0",
				Patterns::Double(0),
				"bottom radius of wire");

		prm.declare_entry("grain boundary width", "0.1",
				Patterns::Double(0),
				"width of boundary layer");

		prm.declare_entry("n_type width", "0.6",
						Patterns::Double(0),
						"width of n type layer");

		prm.declare_entry("p_type width", "4.0",
						Patterns::Double(0),
						"width of p type layer");

		prm.declare_entry("grain boundary depletion width", "1.0",
						Patterns::Double(0),
						"estimated grain boundary depletion width (probably by try and error process) "
						"- for the purpose of mesh refinement");

		prm.leave_subsection();



		prm.enter_subsection("electrons");
		prm.declare_entry("mobility", "1350.0",
				Patterns::Double(0),
				"electron mobility [cm^{2}/(Vs)]");

		prm.declare_entry("recombination time", "5e-5",
				Patterns::Double(0),
				"Recombination rate/time of electrons [s]");

		prm.declare_entry("grain boundary recombination time", "5e-5",
				Patterns::Double(0),
				"Recombination rate/time of electrons [s]");

		prm.declare_entry("electron effective mass", "0.2",
				 Patterns::Double(0),
				 "electron effective mass in conduction band [cm/s]");

		prm.leave_subsection();




		prm.enter_subsection("holes");
		prm.declare_entry("mobility", "480.0",
				Patterns::Double(0),
				"hole mobility [cm^{2}/(Vs)]");

		prm.declare_entry("recombination time", "5e-5",
				Patterns::Double(0),
				"Recombination rate/time of holes [s]");

		prm.declare_entry("grain boundary recombination time", "5e-5",
				Patterns::Double(0),
				"Recombination time of holes [s]");

		prm.declare_entry("hole effective mass", "0.2",
				 Patterns::Double(0),
				 "hole effective mass in valence band [cm/s]");
		prm.leave_subsection();


		prm.enter_subsection("physical");
		prm.declare_entry("real domain height", "0.5",
				Patterns::Double(0.1,5.0),
				"real domain height (or more correctly width) [cm]");

		prm.declare_entry("device thickness", "0.5",
				Patterns::Double(0.1,5.0),
				"device thickness [cm]");

		prm.declare_entry("applied bias", "0.0",
				Patterns::Double(-10.0,10.0),
				"the applied bias [v]");

		prm.declare_entry("band gap", "1.1",
				Patterns::Double(0),
				"the band gap of semiconductor [eV]");

		prm.declare_entry("defect energy level", "0.3",
				Patterns::Double(0.0,1.0),
				"defect level [in faction of band gap]");

		prm.declare_entry("grain boundary defect energy level", "0.2",
				Patterns::Double(0.0,1.0),
				"grain boundary defect energy level [in faction of band gap]");

		prm.declare_entry("temperature", "300.0",
				Patterns::Double(0),
				"temperature of the device [K]");

		prm.declare_entry("schottky bias", "0.0",
				Patterns::Double(0),
				"the Schottky barrier/applied  bias [V]");

		prm.declare_entry("illumination status", "false",
				Patterns::Bool(),
				"true means that the cell is illuminated, false means its not.");

		prm.declare_entry("insulated", "false",
				Patterns::Bool(),
				"whether device is insulated. See Grid::make_Neumann_boundaries");

		prm.declare_entry("schottky status", "false",
				Patterns::Bool(),
				"whether to have schottky contact. See Grid::make_Schottky_boundaries");

		prm.declare_entry("grain boundary status", "false",
				Patterns::Bool(),
				"whether to have grain boundary.");

		prm.declare_entry("vertical grain boundary", "true",
				Patterns::Bool(),
				"whether grain boundary is vertical (perpendicular to the current flow)");

		prm.declare_entry("horizontal grain boundary", "false",
				Patterns::Bool(),
				"whether grain boundary is horizontal (parallel to current flow)");

		prm.declare_entry("characteristic length", "1.0e-4",
				Patterns::Double(0),
				"the characteristic length scale [cm]");

		prm.declare_entry("n_type donor density", "1.0e18",
				Patterns::Double(0),
				"doping in n type region [cm^{-3}]");

		prm.declare_entry("n_type acceptor density", "0.0",
				Patterns::Double(0),
				"doping in n type region [cm^{-3}]");

		prm.declare_entry("p_type donor density", "0.0",
				Patterns::Double(0),
				"doping in p type region [cm^{-3}]");

		prm.declare_entry("p_type acceptor density", "1.0e15",
				Patterns::Double(0),
				"doping in p type region [cm^{-3}]");

		prm.declare_entry("grain boundary defect density", "1.0e14",
				Patterns::Double(0),
				"defect density within grain boundary [cm^{-3}]");

		prm.declare_entry("gb electron density on dirichlet", "1.0e10",
				Patterns::Double(0),
				"electron density on ohmic contact between back contact (left side of the device in paraview) and grain boundary [cm^{-3}]");

		prm.declare_entry("gb hole density on dirichlet", "1.0e10",
				Patterns::Double(0),
				"hole density on ohmic contact between back electrode (left side of the device in paraview) and grain boundary [cm^{-3}]");


		prm.declare_entry("characteristic time steady state", "1.0e-12",
				Patterns::Double(0),
				"the characteristic time scale [s]");

		prm.declare_entry("characteristic time equilibrium state", "1.0e-12",
				Patterns::Double(0),
				"the characteristic time scale [s]");

		prm.declare_entry("characteristic time DLTS", "1.0e-12",
				Patterns::Double(0),
				"the characteristic time scale [s]");

		prm.declare_entry("intrinsic density", "2.564e9",
				Patterns::Double(0),
				"the intrinsic density scale [cm^{-3}]");

		prm.declare_entry("semiconductor permittivity", "11.9",
				Patterns::Double(0),
				"semiconductor permittivity const []");

		prm.declare_entry("electrolyte permittivity", "100",
				Patterns::Double(0),
				"electrolyte permittivity const []");

		prm.declare_entry("photon flux", "1.2e17",
				Patterns::Double(0),
				"intensity of light  [cm^{-2}s^{-1} ]");
	
		prm.declare_entry("absorption coefficient", "1.74974e5",
				Patterns::Double(0),
				"absoprtion coefficient averaged over all energies  [cm^{-1} ]");

		prm.leave_subsection();


		prm.enter_subsection("IV");
		prm.declare_entry("IV minimal voltage", "-1.0",
				Patterns::Double(-10,10),
				"The lowest voltage value on IV curve [V]");

		prm.declare_entry("IV maximal voltage", "1.0",
				Patterns::Double(-10,10),
				"Recombination rate/time of electrons [s]");

		prm.declare_entry("IV number of data points", "20",
				 Patterns::Integer(0,100),
				 "number of points on the curve (without 0V point)");

		prm.declare_entry("delta V", "0.001",
				 Patterns::Double(0.0001,0.1),
				 "small voltage increment for capacitance calculations");

		prm.leave_subsection();
	}
	
	
	void 
	ParameterReader::
	declare_test_parameters()
	{
		prm.enter_subsection("computational");

		prm.declare_entry("global refinements", "2",
				Patterns::Integer(1,10),
				"number of global refinements");

		prm.declare_entry("local refinements", "0",
				Patterns::Integer(0,10),
				"number of local refinements in predefined critical areas");

		prm.declare_entry("time step size", "0.01",
				Patterns::Double(0,1),
				"scaled time step size of both subdomains");

		prm.declare_entry("end time", "1",
				Patterns::Double(0),
				"time to steady state (in terms of characteristic time)");

		prm.declare_entry("end time 2", "20",
				Patterns::Double(0),
				"2nd end time for restarting at end time");

		prm.declare_entry("time stamps", "100",
				Patterns::Integer(0,1000),
				"number of output files until end time");


		prm.declare_entry("restart status", "false",
				Patterns::Bool(),
				"whether to use end of last run as start for this one");

		prm.leave_subsection();				

		prm.enter_subsection("physical");

		prm.declare_entry("applied bias", "0.0",
				Patterns::Double(0),
				"the applied bias [v]");

		prm.declare_entry("built in bias", "0.0",
				Patterns::Double(0),
				"the built in / SCR  bias [V]");

		prm.declare_entry("schottky bias", "0.0",
				Patterns::Double(0),
				"the Schottky barrier/applied  bias [V]");


		prm.declare_entry("illumination status", "false",
				Patterns::Anything(),
				"On means that the cell is illuminated, off means its not.");

		prm.declare_entry("insulated", "false",
				Patterns::Bool(),
				"whether to device is insualted.  See Grid::make_Neumann_boundaries.");

		prm.declare_entry("schottky status", "false",
				Patterns::Bool(),
				"whether to have schottky contact. See Grid::make_Schottky_boundaries");
	
		prm.declare_entry("characteristic length", "1.0",
				Patterns::Double(0),
				"the characteristic length scale [cm]");

		prm.declare_entry("characteristic density", "1.0",
				Patterns::Double(0),
				"the characteristic density scale [cm^{-3}]");

		prm.declare_entry("characteristic time", "1.0",
				Patterns::Double(0),
				"the characteristic time scale [s]");

		prm.declare_entry("end time", "1",
				Patterns::Double(0),
				"physical end time (in terms of characteristic time)");

		prm.declare_entry("intrinsic density", "1.0",
				Patterns::Double(0),
				"the intrinsic density scale [cm^{-3}]");

		prm.declare_entry("semiconductor permittivity", "1.0",
				Patterns::Double(0),
				"semiconductor permittivity const []");

		prm.declare_entry("electrolyte permittivity", "1.0",
				Patterns::Double(0),
				"electrolyte permittivity const []");

		prm.declare_entry("photon flux", "0.0",
				Patterns::Double(0),
				"intensity of light  [cm^{-2}s^{-1} ]");
	
		prm.declare_entry("absorption coefficient", "0.0",
				Patterns::Double(0),
				"absoprtion coefficient averaged over all energies  [cm^{-1} ]");

		prm.leave_subsection();

		prm.enter_subsection("mesh");
		prm.declare_entry("mesh length", "2",
				Patterns::Double(0),
				"scaled domain length in x coordinate");

		prm.declare_entry("mesh height", "1",
				Patterns::Double(0),
				"scaled domain height in y coordinate");
	
		prm.declare_entry("radius one", "1.0",
				Patterns::Double(0),
				"top radius of wire");

		prm.declare_entry("radius two", "1.0",
				Patterns::Double(0),
				"bottom radius of wire");

		prm.declare_entry("boundary layer", "0.1",
				Patterns::Double(0),
				"width of boundary layer");
		prm.leave_subsection();	

		prm.enter_subsection("electrons");

		prm.declare_entry("mobility", "1.0", 
				Patterns::Double(0),
				"electron mobility [v/cm^{2}]");	

		prm.declare_entry("transfer rate", "0",
				Patterns::Double(0),
				"transfer rate of electrons [cm^{4} s^{-1}]");


		prm.declare_entry("recombination time", "0",
				 Patterns::Double(0),
				"Recombination rate/time of electrons [s]");

		prm.declare_entry("recombination velocity", "3e5",
				 Patterns::Double(0),
				"Recombination velocity of electrons [cm/s]");

		prm.leave_subsection();	
	
		prm.enter_subsection("holes");

		prm.declare_entry("mobility", "1.0", 
				Patterns::Double(0),
				"hole mobility [v/cm^{2}]");	

		prm.declare_entry("transfer rate", "0",
				 Patterns::Double(0),
				"transfer rate of holes [cm^{4} s^{-1}]");

		prm.declare_entry("recombination time", "0",
				 Patterns::Double(0),
				"Recombination rate/time of holes [s]");

		prm.declare_entry("recombination velocity", "2.9e-2",
				 Patterns::Double(0),
				"Recombination velocity of holes [cm/s]");

		prm.leave_subsection();	

		prm.enter_subsection("reductants");

		prm.declare_entry("mobility", "1.0", 
				Patterns::Double(0),
				"reductant mobility [v/cm^{2}]");	
	
		prm.leave_subsection();
	
		prm.enter_subsection("oxidants");

		prm.declare_entry("mobility", "1.0", 
				Patterns::Double(0),
				"oxidant mobility [v/cm^{2}]");	
	
		prm.leave_subsection();



	


	}
										 
}


