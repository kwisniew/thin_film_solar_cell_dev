#include "../include/SolarCell.hpp"
#include "../include/Grid.hpp"

#include "Grid.cpp"
#include "Assembly.cpp"
#include "MixedFEM.cpp"
#include "LDG.cpp"
#include "Generation.cpp"
#include "InitialConditions.cpp"
#include "BiasValues.cpp"
#include "CarrierPair.cpp"
#include "Convergence.cpp"


namespace SOLARCELL
{
	using namespace dealii;

	template<int dim>
	SolarCellProblem<dim>::
	SolarCellProblem(const unsigned int degree,
			 ParameterHandler & param)
	:
	degree(degree),
	delta_t(0.1),
	full_system(false),
	prm(param),
	Poisson_dof_handler(Poisson_triangulation),
	Poisson_fe(FE_RaviartThomas<dim>(degree-1),  1,
			   FE_DGQ<dim>(degree-1)          ,  1),
	Poisson_object(),
	semiconductor_dof_handler(semiconductor_triangulation),
	electron_hole_pair(),
	carrier_fe(FESystem<dim>(FE_DGQ<dim>(degree), dim), 1,
		       	   	   	   	 FE_DGQ<dim>(degree),       1),
	Mixed_Assembler(),
	donor_doping_profile(),
	acceptor_doping_profile(),
	schottky_p_type_electrons_eq(),
	schottky_p_type_holes_eq(),
	electrons_initial_condition(),
	holes_initial_condition(),
	electron_density_bc(),
	hole_density_bc(),
	built_in_bias(),
	applied_bias(),
	bulk_bias(),
	schottky_bias(),
	generation()
	{
		// set the parameters
		sim_params.parse_and_scale_parameters(prm);
	
		// set the voltage bias functions
		applied_bias.set_value(sim_params.scaled_applied_bias);
		//std::cout<<"applied bias:  " << sim_params.scaled_applied_bias << std::endl;
		built_in_bias.set_value(sim_params.scaled_built_in_bias);
		schottky_bias.set_value(sim_params.scaled_domain_height);
		schottky_bias.set_value(sim_params.scaled_schottky_bias);

		donor_doping_profile.set_values(
								sim_params.scaled_n_type_donor_density,
								sim_params.scaled_n_type_acceptor_density,
								sim_params.scaled_p_type_width);

		acceptor_doping_profile.set_values(
								sim_params.scaled_p_type_donor_density,
								sim_params.scaled_p_type_acceptor_density,
								sim_params.scaled_p_type_width);

		electrons_initial_condition.set_values(
									sim_params.scaled_n_type_donor_density,
									sim_params.scaled_p_type_donor_density,
									sim_params.scaled_p_type_width,
									sim_params.scaled_n_type_depletion_width,
									sim_params.scaled_p_type_depletion_width);

		holes_initial_condition.set_values(
									sim_params.scaled_n_type_acceptor_density,
									sim_params.scaled_p_type_acceptor_density,
									sim_params.scaled_p_type_width,
									sim_params.scaled_n_type_depletion_width,
									sim_params.scaled_p_type_depletion_width);

		electron_density_bc.set_values(
									sim_params.scaled_n_type_acceptor_density,
									sim_params.scaled_p_type_acceptor_density,
									sim_params.scaled_n_type_donor_density,
									sim_params.scaled_n_type_acceptor_density,
									sim_params.scaled_p_type_width,
									sim_params.scaled_n_type_width,
									sim_params.scaled_intrinsic_density);

		hole_density_bc.set_values(
									sim_params.scaled_n_type_acceptor_density,
									sim_params.scaled_p_type_acceptor_density,
									sim_params.scaled_n_type_donor_density,
									sim_params.scaled_n_type_acceptor_density,
									sim_params.scaled_p_type_width,
									sim_params.scaled_n_type_width,
									sim_params.scaled_intrinsic_density);

		// set the charges name, charge sign, and mobility
		electron_hole_pair.carrier_1.set_name("Electrons");
		electron_hole_pair.carrier_1.charge_number = -1.0;
		electron_hole_pair.carrier_1.scaled_mobility =
						sim_params.scaled_electron_mobility;
	
		electron_hole_pair.carrier_2.set_name("Holes");
		electron_hole_pair.carrier_2.charge_number = 1.0;
		electron_hole_pair.carrier_2.scaled_mobility =
						sim_params.scaled_hole_mobility;

		// set material name
		electron_hole_pair.set_name("Semiconductor-");

		// set the material permittivities
		electron_hole_pair.material_permittivity = 
						sim_params.semiconductor_permittivity;

		// set the generation function to be on or off
		if(sim_params.illum_or_dark)
			generation.set_illuminated_params(sim_params);
		else
			generation.set_dark_params();
	}	//SolarCellProblem

	// destructor
	template<int dim>
	SolarCellProblem<dim>::
	~SolarCellProblem()
	{
		Poisson_dof_handler.clear();
		semiconductor_dof_handler.clear();

	} // ~SolarCellProblem 

	template<int dim>
	void 
	SolarCellProblem<dim>::
	setup_dofs()
	{
		Poisson_object.setup_dofs(Poisson_fe,
					  Poisson_dof_handler);

		electron_hole_pair.setup_dofs(carrier_fe,
					  semiconductor_dof_handler);

	}	// setup+dos

	/////////////////////////////////////////////////////////////////////////////
	//				PRINT INFO
	/////////////////////////////////////////////////////////////////////////////

	template<int dim>
	void
	SolarCellProblem<dim>::
	print_sim_info()
	{
		std::cout << "---------------------------------------------------------\n"
				<< "Triangulation Info\n"
				<< "---------------------------------------------------------\n"
				<< "Number of active cells : " 
				<< Poisson_triangulation.n_active_cells()
				<< std::endl
				<< "Total number of cells: " 
				<< Poisson_triangulation.n_cells() << std::endl
				<< "h_min: "
				<< GridTools::minimal_cell_diameter(Poisson_triangulation) << std::endl
				<< "h_max: "
				<< GridTools::maximal_cell_diameter(Poisson_triangulation) << std::endl
				<< std::endl;
	
		Poisson_object.print_info(Poisson_dof_handler);
		electron_hole_pair.print_info(semiconductor_dof_handler);
	}

	/////////////////////////////////////////////////////////////////////////////
	//				MAPPINGS
	/////////////////////////////////////////////////////////////////////////////

	template<int dim>
	void
	SolarCellProblem<dim>::
	setup_mappings()
	{
		/*---------------------------------------------------------------------*/
		/*                   Semiconductor-Poisson Mapping                 	 */
		/*---------------------------------------------------------------------*/

		std::vector<std::pair<unsigned int, unsigned int>> semiconductor_cells;
		std::vector<Point<dim>> 			   semiconductor_cell_centers;
		std::vector<std::pair<unsigned int, unsigned int>> Poisson_cells;
		std::vector<Point<dim>> 			   Poisson_cell_centers;
		std::vector<std::pair<unsigned int, unsigned int>> electrolyte_cells;
		std::vector<Point<dim>>				   electrolyte_cell_centers;	

		// make mapping of Poisson cells to global poisson index							 
		typename DoFHandler<dim>::active_cell_iterator
		Poisson_cell	=	Poisson_dof_handler.begin_active(),
		Poisson_endc	=	Poisson_dof_handler.end();		
		for(; Poisson_cell != Poisson_endc; Poisson_cell++)	
		{
			Poisson_cells.push_back(
				std::pair<unsigned int, unsigned int>(Poisson_cell->level(),
				Poisson_cell->index()));
			Poisson_cell_centers.push_back(Poisson_cell->center());
		
		} // end Poisson_cell
	
		// make list of semiconductor cells	
		typename DoFHandler<dim>::active_cell_iterator
		semi_cell	=	semiconductor_dof_handler.begin_active(),
		semi_endc	=	semiconductor_dof_handler.end();
		for(; semi_cell != semi_endc; semi_cell++)
		{
			semiconductor_cells.push_back(
					std::pair<unsigned int, unsigned int>(semi_cell->level(),
					semi_cell->index()));
			semiconductor_cell_centers.push_back(semi_cell->center());
		}

		std::pair<unsigned int, unsigned int>		semi_pair;
		std::pair<unsigned int, unsigned int>		Poisson_pair;
		
		std::pair<std::pair<unsigned int, unsigned int>, 
				  std::pair<unsigned int, unsigned int>> mapping_pair;


		// find cells which are in the semiconductor and poisson and make mapping between them
		for(unsigned int i=0; i < Poisson_cells.size(); i++)
		{
			for(unsigned int j=0; j < semiconductor_cells.size(); j++)
			{
				if(Poisson_cell_centers[i].distance(semiconductor_cell_centers[j]) < 1e-13)
				{
					semi_pair = std::pair<unsigned int, unsigned int>(
								 semiconductor_cells[j].first, // level
								 semiconductor_cells[j].second); // index
					
					Poisson_pair = std::pair<unsigned int, unsigned int>(
								Poisson_cells[i].first, // level
								Poisson_cells[i].second); // index

					// mapping from semiconductor cell to Poisson cell
					mapping_pair.first  = semi_pair;				
					mapping_pair.second = Poisson_pair;				
					s_2_p_map.insert(mapping_pair);
					
				} // if
			} // end j
		} // end i
	} // setup_mappings

	////////////////////////////////////////////////////////////////////////////////
	// 			MIXED METHOD ROUTINES
	////////////////////////////////////////////////////////////////////////////////
	template<int dim>
	void 
	SolarCellProblem<dim>::
	assemble_Poisson_matrix()
	{
		// this is shown in deal.II on shared memory parallelism 
		// it builds the above matrix by distributing over multi threads locally
		// and then building up the global matrix sequentially (kinda)
		WorkStream::run(Poisson_dof_handler.begin_active(),
				Poisson_dof_handler.end(),
				std_cxx11::bind(&MixedPoisson::MixedFEM<dim>::
						assemble_local_Poisson_matrix, 
						Mixed_Assembler, // this object object
						std_cxx11::_1,
						std_cxx11::_2,
						std_cxx11::_3,
						sim_params.semiconductor_permittivity,
						//sim_params.electrolyte_permittivity,
						sim_params.scaled_debye_length), 
						std_cxx11::bind(&SOLARCELL::SolarCellProblem<dim>::
						copy_local_to_global_Poisson_matrix,
						this, // this object
						std_cxx11::_1),
				Assembly::AssemblyScratch<dim>(Poisson_fe,
							carrier_fe,
							QGauss<dim>(degree+2),
							QGauss<dim-1>(degree+2)),
				Assembly::Poisson::CopyData<dim>(Poisson_fe)
				);
	} // end assemble_Poisson_matrix

	template <int dim>
	void 
	SolarCellProblem<dim>::
	copy_local_to_global_Poisson_matrix(const Assembly::Poisson::CopyData<dim> & data)
	{
		// distribute local matrix to global Poisson matrix
		Poisson_object.constraints.distribute_local_to_global(data.local_matrix,
								data.local_dof_indices,
								Poisson_object.system_matrix);
	}	// copy_local_to_global_poisson

	template <int dim>
	void 
	SolarCellProblem<dim>::
	copy_local_to_global_Poisson_rhs(const Assembly::Poisson::CopyData<dim> & data)
	{
		// copy the local RHS into the global RHS for Poisson
		Poisson_object.constraints.distribute_local_to_global(data.local_rhs,
						  data.local_dof_indices,
						  Poisson_object.system_rhs);
/*		if(std::accumulate(data.local_rhs.begin(),data.local_rhs.end(),0)>30)
		{
			for (auto i: Poisson_object.system_rhs)
			{
					std::cout << i << "   ";
			}
			for (auto i: data.local_rhs)
			{
					std::cout << i << "   ";
			}
			std::cout<<std::endl;
			for (auto i: data.local_dof_indices)
			{
					std::cout << i << "   ";
			}
			std::cout<<std::endl;
			std::cout<<"Wartość elementu lokalnie: " << data.local_rhs[0] << ", a globalnie:  " <<Poisson_object.system_rhs[data.local_dof_indices[0]]<<std::endl;
			std::cout<<std::endl;
			std::cout<<std::endl;
		}*/
	}	

	template <int dim>
	void
	SolarCellProblem<dim>::
	assemble_Poisson_rhs()
	{
		// reset rhs to zero
		Poisson_object.system_rhs = 0;

		/*-------------------------------------------------------------*/
		// This is the one coupled to the semiconductor
		/*-------------------------------------------------------------*/
		WorkStream::run(semiconductor_dof_handler.begin_active(),
				semiconductor_dof_handler.end(),
				std_cxx11::bind(&SOLARCELL::SolarCellProblem<dim>::
						assemble_local_Poisson_rhs_for_semiconductor,
						this, // this object
						std_cxx11::_1,
						std_cxx11::_2,
						std_cxx11::_3),
				std_cxx11::bind(&SOLARCELL::SolarCellProblem<dim>::
						copy_local_to_global_Poisson_rhs,
						this, // this object
						std_cxx11::_1),
				Assembly::AssemblyScratch<dim>(Poisson_fe,
								carrier_fe,
								QGauss<dim>(degree+2),
								QGauss<dim-1>(degree+2)),
				Assembly::Poisson::CopyData<dim>(Poisson_fe)
				);
		
/*
		if(full_system)
		{
			-------------------------------------------------------------
			// This is the one coupled to the electrolyte
			-------------------------------------------------------------
			WorkStream::run(electrolyte_dof_handler.begin_active(),
					electrolyte_dof_handler.end(),
					std_cxx11::bind(&SOLARCELL::SolarCellProblem<dim>::
							assemble_local_Poisson_rhs_for_electrolyte,
							this, // this object
							std_cxx11::_1,
							std_cxx11::_2,
							std_cxx11::_3),
					std_cxx11::bind(&SOLARCELL::SolarCellProblem<dim>::
							copy_local_to_global_Poisson_rhs,
							this, // this object
							std_cxx11::_1),
					Assembly::AssemblyScratch<dim>(Poisson_fe,
									carrier_fe,
									QGauss<dim>(degree+2),
									QGauss<dim-1>(degree+2)),
					Assembly::Poisson::CopyData<dim>(Poisson_fe)
					);

		} // end if full_system
*/
	}

	template <int dim>
	void 
	SolarCellProblem<dim>::
	assemble_local_Poisson_rhs_for_semiconductor(
				const typename DoFHandler<dim>::active_cell_iterator 	& cell,
				Assembly::AssemblyScratch<dim>				& scratch,
				Assembly::Poisson::CopyData<dim>			& data)
	{
		const unsigned int	dofs_per_cell = scratch.Poisson_fe_values.dofs_per_cell;
		//std::cout<<"stopnie swobody w komórce Poisson:   " << dofs_per_cell <<std::endl;

		const unsigned int 	n_q_points    = scratch.Poisson_fe_values.n_quadrature_points;

		const unsigned int	n_face_q_points =	
							scratch.Poisson_fe_face_values.n_quadrature_points;


		// Get the actual values for vector field and potential from FEValues
		// Use Extractors instead of having to deal with shapefunctions directly
		const FEValuesExtractors::Vector VectorField(0); // Vector as in Vector field
		const FEValuesExtractors::Scalar Potential(dim);
		const FEValuesExtractors::Scalar Density(dim);

		// reset the local_rhs vector to be zero
		data.local_rhs=0;

		// get the Poisson cell info
		std::pair<unsigned int, unsigned int> Poisson_cell_info = 
					s_2_p_map[std::pair<unsigned int, unsigned int>(cell->level(),
											cell->index())];

		// get the Poisson cell
		typename DoFHandler<dim>::active_cell_iterator
					Poisson_cell(&Poisson_triangulation,
						Poisson_cell_info.first,
						Poisson_cell_info.second,
						&Poisson_dof_handler);

		Poisson_cell->get_dof_indices(data.local_dof_indices);
		scratch.Poisson_fe_values.reinit(Poisson_cell);
		// get the test rhs for poisson
		// Assemble the right hand side on semiconductor side
		scratch.carrier_fe_values.reinit(cell);
			
		// get the electron density values at the previous time step
		scratch.carrier_fe_values[Density].get_function_values(
						electron_hole_pair.carrier_1.solution,
						scratch.old_carrier_1_density_values);

		// get the hole density values at the previous time step
		scratch.carrier_fe_values[Density].get_function_values(
						electron_hole_pair.carrier_2.solution,
						scratch.old_carrier_2_density_values);
	

		// get doping profiles values on this cell
		donor_doping_profile.value_list(scratch.carrier_fe_values.get_quadrature_points(),
					scratch.donor_doping_values,
					dim); // calls the density values of the donor profile

		acceptor_doping_profile.value_list(scratch.carrier_fe_values.get_quadrature_points(),
				scratch.acceptor_doping_values,
				dim); // calls the density values of the acceptor profile

		// Loop over all the quadrature points in this cell
		for(unsigned int q=0; q<n_q_points; q++)
		{
			// copy over the test functions
			for(unsigned int k = 0; k<dofs_per_cell; k++)
				scratch.psi_i_potential[k] = scratch.Poisson_fe_values[Potential].value(k,q);

			// loop over the test function dofs for this cell
			for(unsigned int i=0; i<dofs_per_cell; i++)
			{
				// i-th potential basis functions at the point q
		
				// get the local RHS values for this cell
				// = -int_{Omega_{e}} (1/lambda^{2}) v	(N_{D} - N_{A}) - (n - p) d
				data.local_rhs(i) += scratch.psi_i_potential[i] * //-psi_i_potential *
						(
						(scratch.donor_doping_values[q] 
						- 
						scratch.acceptor_doping_values[q])
						+
						(electron_hole_pair.carrier_1.charge_number *
						scratch.old_carrier_1_density_values[q]	
						+
						electron_hole_pair.carrier_2.charge_number *
						scratch.old_carrier_2_density_values[q])
						) * scratch.Poisson_fe_values.JxW(q);
				if(i==4)
				{
					//std::cout<<"Numer q-pointa:   " << q <<"   Numer wierzchołka: " << i <<"Wartości potensjału:   "<< scratch.psi_i_potential[i] <<std::endl;
				}
			} // for i
		} // for q

		// loop over all the faces of this cell to calculate the vector
		// from the dirichlet boundary conditions if the face is on the 
		// Dirichlet portion of the boundary
		for(unsigned int face_no=0;
				face_no<GeometryInfo<dim>::faces_per_cell;
				face_no++)
		{
			// obtain the face_no-th face of this cell
			typename DoFHandler<dim>::face_iterator 	face = Poisson_cell->face(face_no);

			// apply Dirichlet boundary conditions.. 
			// Since we are in a semicondcutor cell we know to apply the 
			// biases
			if(face->at_boundary())
			{
				if(face->boundary_id() == Dirichlet)
				{	
					//std::cout<< "Dirichlet!" <<std::endl;
					// get the values of the shape functions at this boundary face
					scratch.Poisson_fe_face_values.reinit(Poisson_cell,face_no);
			
					// get the values of the dirichlet boundary conditions evaluated
					// on the quadrature points of this face
					built_in_bias.value_list(
							scratch.Poisson_fe_face_values.get_quadrature_points(),
							scratch.Poisson_bi_values);


/*					for (auto i: scratch.Poisson_bi_values)
					{
						if(i > 0)
						{
							std::cout << i << ' ';
						}

					}*/

					applied_bias.value_list(
							scratch.Poisson_fe_face_values.get_quadrature_points(),
							scratch.Poisson_bc_values);
			
					// copy over the normal vectors
					for(unsigned int k=0; k < n_face_q_points; k++)
						scratch.normals[k] = scratch.Poisson_fe_face_values.normal_vector(k);

					// loop over all the quadrature points of this face
					for(unsigned int q=0; q<n_face_q_points; q++)
					{
						// copy over the test functions
						for(unsigned int k = 0; k < dofs_per_cell; k++)
						{
							scratch.psi_i_field[k]=	
								scratch.Poisson_fe_face_values[VectorField].value(k,q);
						}

						// loop over all the test function dofs of this face
						for(unsigned int i=0; i<dofs_per_cell; i++)
						{
							// - \int_{face} p * n * (phi_{Dirichlet}) dx
							data.local_rhs(i) +=	
									-(scratch.psi_i_field[i] *
									scratch.normals[q] *
									(scratch.Poisson_bi_values[q] 
									-
									scratch.Poisson_bc_values[q]) *
									scratch.Poisson_fe_face_values.JxW(q));
							//std::cout<<"Numer q-pointa:   " << q <<"   Numer wierzchołka: " << i << "  wartość W: "  << scratch.psi_i_field[i] << std::endl;
							/*std::cout<< -(scratch.psi_i_field[i] *
									scratch.normals[q] *
									(scratch.Poisson_bi_values[q]
									-
									scratch.Poisson_bc_values[q]) *
									scratch.Poisson_fe_face_values.JxW(q)) <<std::endl;*/
						} // for i
						//std::cout<<scratch.Poisson_bi_values[q] << "   " << scratch.normals[q] <<std::endl;
					} // for q
				} // end if dirichlet

				if(face->boundary_id() == Schottky)
				{	
					// get the values of the shape functions at this boundary face
					scratch.Poisson_fe_face_values.reinit(Poisson_cell,face_no);
			
					// get the values of the dirichlet boundary conditions evaluated
					// on the quadrature points of this face
					built_in_bias.value_list(
								scratch.Poisson_fe_face_values.get_quadrature_points(),
								scratch.Poisson_bi_values);
					
					schottky_bias.value_list(
								scratch.Poisson_fe_face_values.get_quadrature_points(),
								scratch.Poisson_bc_values);
			
					// copy over the normal vectors
					for(unsigned int k=0; k < n_face_q_points; k++)
						scratch.normals[k] = scratch.Poisson_fe_face_values.normal_vector(k);

					// loop over all the quadrature points of this face
					for(unsigned int q=0; q<n_face_q_points; q++)
					{
						// copy over the test functions
						for(unsigned int k = 0; k < dofs_per_cell; k++)
						{
							scratch.psi_i_field[k]=	
								scratch.Poisson_fe_face_values[VectorField].value(k,q);
						}

						// loop over all the test function dofs of this face
						for(unsigned int i=0; i<dofs_per_cell; i++)
						{
							// - \int_{face} p * n * (phi_{Dichlet}) dx
							data.local_rhs(i) +=	
								 -(scratch.psi_i_field[i] *
								 scratch.normals[q] *
								(scratch.Poisson_bi_values[q] 
								 -
								scratch.Poisson_bc_values[q]) *
								scratch.Poisson_fe_face_values.JxW(q));
						} // for i
					} // for q
				} // end if schottky
			} // if boundary
		} // end for face_no
		//std::cout<< data.local_rhs <<std::endl;
	} // assemble_local_Poisson_rhs_for_semiconductor

	/////////////////////////////////////////////////////////////////////////////////
	// 		LOCAL DISCONTINUOUS GALERKIN ROUTINES
	/////////////////////////////////////////////////////////////////////////////////

	template <int dim>
	void 
	SolarCellProblem<dim>::
	assemble_LDG_system(const double & transient_or_steady)
	{	
		// assemble mass matrix for the electron_hole pair
//		std::cout<<"MACIERZ PRZED"<<std::endl;
//
//		electron_hole_pair.carrier_1.system_matrix.print_formatted(std::cout,0);

//		std::cout << "Czy chodzi o bulk?" << std::endl;
		WorkStream::run(semiconductor_dof_handler.begin_active(),
				semiconductor_dof_handler.end(),
				std_cxx11::bind(&LDG_System::LDG<dim>::
						assemble_local_LDG_mass_matrix,
						LDG_Assembler, // the Assembler object
						std_cxx11::_1,
						std_cxx11::_2,
						std_cxx11::_3,
						delta_t),
				std_cxx11::bind(&SOLARCELL::SolarCellProblem<dim>::
						copy_local_to_global_semiconductor_mass_matrix,
						this, 			// this object
						std_cxx11::_1),
				Assembly::AssemblyScratch<dim>(Poisson_fe,
							  carrier_fe,
							  QGauss<dim>(degree+2),
							  QGauss<dim-1>(degree+2)),
				Assembly::DriftDiffusion::CopyData<dim>(carrier_fe)
				);
//		std::cout<<"MACIERZ Po BULKU"<<std::endl;
//
//		electron_hole_pair.carrier_1.system_matrix.print_formatted(std::cout,0);
		// system matrix for the electron_hole pair
//		std::cout << "Czy o boundary condition?" << std::endl;
		WorkStream::run(semiconductor_dof_handler.begin_active(),
				semiconductor_dof_handler.end(),
				std_cxx11::bind(&LDG_System::LDG<dim>::
						assemble_local_LDG_cell_and_bc_terms,
						LDG_Assembler, // the Assembler object
						std_cxx11::_1,
						std_cxx11::_2,
						std_cxx11::_3,
						electron_hole_pair.carrier_1.scaled_mobility,
						electron_hole_pair.carrier_2.scaled_mobility,
						delta_t,
						transient_or_steady,
						electron_hole_pair.penalty),
				std_cxx11::bind(&SOLARCELL::SolarCellProblem<dim>::
						copy_local_to_global_semiconductor_system_matrix,
						this, 			// tthis object
						std_cxx11::_1),
				Assembly::AssemblyScratch<dim>(Poisson_fe,
							 carrier_fe,
							 QGauss<dim>(degree+2),
							 QGauss<dim-1>(degree+2)),
				Assembly::DriftDiffusion::CopyData<dim>(carrier_fe)
				);
//		std::cout<<"MACIERZ PO BOUNDARY CONDITION"<<std::endl;
//
//		electron_hole_pair.carrier_1.system_matrix.print_formatted(std::cout,0);

//		std::cout << "A może to dopiero fluxy?" << std::endl;
		// LDG FLLUXES
		LDG_Assembler.assemble_flux_terms(semiconductor_dof_handler,
						electron_hole_pair,
						Poisson_fe,
						carrier_fe);	
//		std::cout << "A może wszystko działa" << std::endl;
//		std::cout<<"MACIERZ PO FLUXACH"<<std::endl;
//		electron_hole_pair.carrier_1.system_matrix.print_formatted(std::cout,0);

	}

	template<int dim>
	void
	SolarCellProblem<dim>::
	copy_local_to_global_semiconductor_mass_matrix(
				const Assembly::DriftDiffusion::CopyData<dim>	& data)
	{
		// copy local mass matrix into the global mass matrix.. 
		electron_hole_pair.constraints.distribute_local_to_global(
								data.local_mass_matrix,
								data.local_dof_indices,
								electron_hole_pair.mass_matrix);
	}

	template<int dim>
	void
	SolarCellProblem<dim>::
	copy_local_to_global_semiconductor_system_matrix(
				const Assembly::DriftDiffusion::CopyData<dim>	& data)
	{
		// copy local system matrix into the global system matrix.. 
		electron_hole_pair.constraints.distribute_local_to_global(
								data.local_matrix_1,
								data.local_dof_indices,
								electron_hole_pair.carrier_1.system_matrix);

		// copy local system matrix into the global system matrix.. 
		electron_hole_pair.constraints.distribute_local_to_global(
								data.local_matrix_2,
								data.local_dof_indices,
								electron_hole_pair.carrier_2.system_matrix);
	}

	template<int dim>
	void
	SolarCellProblem<dim>::
	copy_local_to_global_semiconductor_system_rhs(						
						const Assembly::DriftDiffusion::CopyData<dim>	& data)
	{
		// copy local system rhs into the global system rhs.. 
		electron_hole_pair.constraints.distribute_local_to_global(
								data.local_carrier_1_rhs,
								data.local_dof_indices,
								electron_hole_pair.carrier_1.system_rhs);

		// copy local syste, matrix into the global system matrix.. 
		electron_hole_pair.constraints.distribute_local_to_global(
								data.local_carrier_2_rhs,
								data.local_dof_indices,
								electron_hole_pair.carrier_2.system_rhs);
	}

	template<int dim>
	void
	SolarCellProblem<dim>::
	assemble_semiconductor_rhs()
	{
		// set carrier_system_rhs = M * u^{n-1}
		electron_hole_pair.mass_matrix.vmult(electron_hole_pair.carrier_1.system_rhs,
						 electron_hole_pair.carrier_1.solution);
		//std::cout<<"1.1"<<std::endl;

		electron_hole_pair.mass_matrix.vmult(electron_hole_pair.carrier_2.system_rhs,
						 electron_hole_pair.carrier_2.solution);
		//std::cout<<"1.2"<<std::endl;
		// now run through the semiconductor cells and assemble rhs for the 
		// LDG equations.
		WorkStream::run(semiconductor_dof_handler.begin_active(),
				semiconductor_dof_handler.end(),
				std_cxx11::bind(&SOLARCELL::SolarCellProblem<dim>::
						assemble_local_semiconductor_rhs,
						this, // this object
						std_cxx11::_1,
						std_cxx11::_2,
						std_cxx11::_3,
						electron_hole_pair.penalty),
				std_cxx11::bind(&SOLARCELL::SolarCellProblem<dim>::
						copy_local_to_global_semiconductor_system_rhs,
						this, // this object
						std_cxx11::_1),
				Assembly::AssemblyScratch<dim>(Poisson_fe,
								carrier_fe,
								QGauss<dim>(degree+2),
								QGauss<dim-1>(degree+2)),
				Assembly::DriftDiffusion::CopyData<dim>(carrier_fe)
				);

	}		

	template<int dim>
	void 
	SolarCellProblem<dim>::
	assemble_local_semiconductor_rhs(
			const typename DoFHandler<dim>::active_cell_iterator & cell,
			Assembly::AssemblyScratch<dim>			     & scratch,
			Assembly::DriftDiffusion::CopyData<dim>	  	     & data,
			const double  					     & penalty)									
	{
		/*-------------------------------------------------------------------*/
		// This is for the case of a coupled Poisson's equation and
		// reactive interface
		/*-------------------------------------------------------------------*/

		// this assembles the drift term in the ldg formulation.  it uses the 
		// electric field at the current iteration and the density of the 
		// carrier at the previous time step
		const unsigned int dofs_per_cell	= 
							scratch.carrier_fe_values.dofs_per_cell;

		//std::cout<< "stopnie swobody LDG!:  " << dofs_per_cell << std::endl;
		const unsigned int n_q_points		=  
							scratch.carrier_fe_values.n_quadrature_points;
		const unsigned int n_face_q_points	=	
							scratch.carrier_fe_face_values.n_quadrature_points;

		cell->get_dof_indices(data.local_dof_indices);

		double h = cell->diameter();

		// get the Poisson cell info
		std::pair<unsigned int, unsigned int> Poisson_cell_info = 
				s_2_p_map[std::pair<unsigned int, unsigned int>(cell->level(),
										cell->index())];
		// get the Poisson cell
		typename DoFHandler<dim>::active_cell_iterator
					Poisson_cell(&Poisson_triangulation,
						Poisson_cell_info.first,
						Poisson_cell_info.second,
						&Poisson_dof_handler);

		// reinitialize the fe_values 
		scratch.carrier_fe_values.reinit(cell);
		scratch.Poisson_fe_values.reinit(Poisson_cell);
		
		// reset the local_rhs to be zero
		data.local_carrier_1_rhs=0;
		data.local_carrier_2_rhs=0;
	
		const FEValuesExtractors::Vector Current(0);
		const FEValuesExtractors::Scalar Density(dim);

		const FEValuesExtractors::Vector ElectricField(0);

		// get the values of carrier_1 and carrier_2 densities at 
		// the previous time step
		scratch.carrier_fe_values[Density].get_function_values(
							electron_hole_pair.carrier_1.solution,
							scratch.old_carrier_1_density_values);

		scratch.carrier_fe_values[Density].get_function_values(
							electron_hole_pair.carrier_2.solution,
							scratch.old_carrier_2_density_values);

		generation.value_list(scratch.carrier_fe_values.get_quadrature_points(),
					 scratch.generation_values);

		// get the electric field values at the previous time step
		scratch.Poisson_fe_values[ElectricField].get_function_values(	
									Poisson_object.solution,
									scratch.electric_field_values);

		const double inverse_perm		  =	1.0/sim_params.semiconductor_permittivity;
		const double inverse_debye_length = 1.0/sim_params.scaled_debye_length;

		// loop over all the quadrature points in this cell and compute body integrals
		for(unsigned int q=0; q<n_q_points; q++)
		{
			// copy over the test functions
			for(unsigned int k=0; k<dofs_per_cell; k++)
				scratch.psi_i_density[k] = scratch.carrier_fe_values[Density].value(k,q);

			for(unsigned int k=0; k<dofs_per_cell; k++)
				scratch.psi_i_current[k] = scratch.carrier_fe_values[Current].value(k,q);

			// loop over all the test function dofs and get the test functions
			for(unsigned int i=0; i<dofs_per_cell; i++)
			{
				// contribution from RHS function + Drift
				// int_{Omega} v * R dx
				data.local_carrier_1_rhs(i) += ( 
						(scratch.psi_i_density[i] * scratch.generation_values[q])
						+
						(scratch.psi_i_density[i] *
						SRH_Recombination(scratch.old_carrier_1_density_values[q],
										  scratch.old_carrier_2_density_values[q],
										  sim_params))
						+
						electron_hole_pair.carrier_1.charge_number *
				   	  	scratch.psi_i_current[i] *
						inverse_perm *
						inverse_debye_length *
						scratch.electric_field_values[q] *
						scratch.old_carrier_1_density_values[q]
						) * scratch.carrier_fe_values.JxW(q);

				data.local_carrier_2_rhs(i) += ( 
						(scratch.psi_i_density[i] * scratch.generation_values[q])
						+
						(scratch.psi_i_density[i] *
						SRH_Recombination(scratch.old_carrier_1_density_values[q],
									      scratch.old_carrier_2_density_values[q],
										  sim_params))
		 			  	+
					  	electron_hole_pair.carrier_2.charge_number *
					  	scratch.psi_i_current[i] *
					    inverse_perm *
						inverse_debye_length *
						scratch.electric_field_values[q] *
						scratch.old_carrier_2_density_values[q]
						) * scratch.carrier_fe_values.JxW(q);

			} // for i
		}	// for q
		
		// loop over all the faces of this cell and compute the contribution 
		// from the boundary conditions
		for(unsigned int face_no=0; 
				face_no< GeometryInfo<dim>::faces_per_cell;
				face_no++)
		{
			// get the face_no-th face of this cell
			typename DoFHandler<dim>::face_iterator 	face = cell->face(face_no);
			
			// if on boundary apply boundary conditions
			if(face->at_boundary() )
			{
				// reinitialize the fe_face_values for this cell ONLY if it is as the
				//boundary otherwise its a waste.  then assemble the appropriate
				//boundary conditions
				scratch.carrier_fe_face_values.reinit(cell, face_no);
			
				if(face->boundary_id() == Dirichlet)
				{
					// Get the doping profile values for the boundary conditions
					electron_density_bc.value_list(
								scratch.carrier_fe_face_values.get_quadrature_points(),
								scratch.carrier_1_bc_values,
								dim); // calls the density values of the donor profile
								     // not the current ones
					hole_density_bc.value_list(
								scratch.carrier_fe_face_values.get_quadrature_points(),
								scratch.carrier_2_bc_values,
								dim); // calls the density values of the donor profile
								     // not the current ones
					// copy over normal vectors
					for(unsigned int k=0; k<n_face_q_points; k++)
						scratch.normals[k] = scratch.carrier_fe_face_values.normal_vector(k);
	
					// loop over all the quadrature points on this face
					for(unsigned int q=0; q<n_face_q_points; q++)
					{
						// copy over the test functions
						for(unsigned int k=0; k<dofs_per_cell; k++)
						{
							scratch.psi_i_density[k] = 
								scratch.carrier_fe_face_values[Density].value(k,q);
						}
						for(unsigned int k=0; k<dofs_per_cell; k++)
						{
							scratch.psi_i_current[k] =
								 scratch.carrier_fe_face_values[Current].value(k,q);
						}
						// loop over all the test function dofs on this face
						for(unsigned int i=0; i<dofs_per_cell; i++)
						{
						// int_{\Gamma_{D}} ( -p^{-} n^{-} + penalty/h * v^{-}) * u_{D} ds
							data.local_carrier_1_rhs(i) += 
								    (-1.0 * scratch.psi_i_current[i] *
								    scratch.normals[q]
								    + 
								    (penalty/h) * scratch.psi_i_density[i]) *
								    scratch.carrier_1_bc_values[q] *
								    scratch.carrier_fe_face_values.JxW(q);				
						
							data.local_carrier_2_rhs(i) +=  
								   (-1.0 * scratch.psi_i_current[i] *
								   scratch.normals[q] 
								   + 
								   (penalty/h) * scratch.psi_i_density[i]) *
								   scratch.carrier_2_bc_values[q] *
								   scratch.carrier_fe_face_values.JxW(q);				
						} // for i
					}	// for q
				} // end Dirichlet
				else if(face->boundary_id() == Schottky)
				{
					// Get the doping profile values for the boundary conditions
					electron_density_bc.value_list(
								scratch.carrier_fe_face_values.get_quadrature_points(),
								scratch.carrier_1_bc_values,
								dim); // calls the density values of the donor profile
									  // not the current ones
					hole_density_bc.value_list(
								scratch.carrier_fe_face_values.get_quadrature_points(),
								scratch.carrier_2_bc_values,
								dim); // calls the density values of the donor profile
									  // not the current ones
					
	
					// get the electrona and hole densities at the pevious time step
					scratch.carrier_fe_face_values[Density].get_function_values(
									electron_hole_pair.carrier_1.solution,
									scratch.electron_interface_values);

					scratch.carrier_fe_face_values[Density].get_function_values(
									electron_hole_pair.carrier_2.solution,
									scratch.hole_interface_values);

	
					// loop over all the quadrature points on this face
					for(unsigned int q=0; q<n_face_q_points; q++)
					{
						// copy over the test functions
						for(unsigned int k=0; k<dofs_per_cell; k++)
						{
							scratch.psi_i_density[k] = 
								scratch.carrier_fe_face_values[Density].value(k,q);
						}
						// loop over all the test function dofs on this face
						for(unsigned int i=0; i<dofs_per_cell; i++)
						{
							// int_{\Sigma} -v^{-} ket (rho_n - rho_n^e) rho_o ds
							data.local_carrier_1_rhs(i) += 
									-1.0 * scratch.psi_i_density[i] *
									sim_params.scaled_electron_recombo_v*
									(scratch.electron_interface_values[q]
									-
									scratch.carrier_1_bc_values[q] ) * 
									scratch.carrier_fe_face_values.JxW(q);			
					
							// int_{\Sigma}  v^{-} kht (rho_p - rho_p^e) rho_r ds
							data.local_carrier_2_rhs(i) +=  
									+1.0 * scratch.psi_i_density[i] *
									sim_params.scaled_hole_recombo_v *
									(scratch.hole_interface_values[q]
									-
									scratch.carrier_2_bc_values[q] ) *
									scratch.carrier_fe_face_values.JxW(q);						
						} // for i
					} // for q
				}
				else if(face->boundary_id() == Neumann)
				{
					// NOTHIN TO DO IF INSULATING
				}
				else
					Assert(false, ExcNotImplemented() );
			} // end at boundary
		} // end for face_no
	} // end assemble_local_semiconductor_rhs

	////////////////////////////////////////////////////////////////////////
	//			SOLVERS
	////////////////////////////////////////////////////////////////////////

	template <int dim>
	void 
	SolarCellProblem<dim>::
	set_solvers()
	{
		//std::cout << "Ustawiam Solver: Poisson" << std::endl;
//		freopen("output.txt","w",stdout);
//		Poisson_object.system_matrix.print_formatted(std::cout,0);
		//fclose(stdout);

		Poisson_object.set_solver();
		//std::cout << "Ustawim Solver: Elektrony" << std::endl;
		//electron_hole_pair.carrier_1.system_matrix.print_formatted(std::cout,0);
		electron_hole_pair.carrier_1.set_solver();
		//std::cout << "Usatwiam Solver: Dziury" << std::endl;
		electron_hole_pair.carrier_2.set_solver();
	}
	
	template <int dim>
	void 
	SolarCellProblem<dim>::
	solve_Poisson()
	{
		Poisson_object.solve();	
	}		
	
	template<int dim>
	void
	SolarCellProblem<dim>::
	solve_full_system()
	{
		Threads::TaskGroup<void> task_group;
		
		task_group += Threads::new_task(&ChargeCarrierSpace::
						 Carrier<dim>::solve,
						 electron_hole_pair.carrier_1);

		task_group += Threads::new_task(&ChargeCarrierSpace::
						 Carrier<dim>::solve,
						 electron_hole_pair.carrier_2);

/*		task_group += Threads::new_task(&ChargeCarrierSpace::
						 Carrier<dim>::solve,
						 redox_pair.carrier_1);

		task_group += Threads::new_task(&ChargeCarrierSpace::
						 Carrier<dim>::solve,
						 redox_pair.carrier_2);*/

		task_group.join_all();
	} // solve_full_system

	template<int dim>
	void
	SolarCellProblem<dim>::
	solve_semiconductor_system()
	{
		Threads::TaskGroup<void> task_group;
		
		task_group += Threads::new_task(&ChargeCarrierSpace::
						Carrier<dim>::solve,
						electron_hole_pair.carrier_1);

		task_group += Threads::new_task(&ChargeCarrierSpace::
						 Carrier<dim>::solve,
						 electron_hole_pair.carrier_2);
		task_group.join_all();
	} // solve_semiconductor_system

	////////////////////////////////////////////////////////////////////////
	//			PRINTING & POSTPROCESSING
	////////////////////////////////////////////////////////////////////////
	template<int dim>
	double
	SolarCellProblem<dim>::
	calculate_uncompensated_charge()
	{
		QGauss<dim> quad_rule(degree+2);
		FEValues<dim> fe_values(carrier_fe,
		                      	quad_rule,
		                        update_values | update_JxW_values | update_quadrature_points);

		const unsigned int dofs_per_cell = carrier_fe.dofs_per_cell;
		const unsigned int n_q_points    = quad_rule.size();

		double uncompendated_charge = 0;

		std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
		std::vector<double>			donor_doping_values(n_q_points);
		std::vector<double>			acceptor_doping_values(n_q_points);
		std::vector<double>			carrier_1_density_values(n_q_points);
		std::vector<double>			carrier_2_density_values(n_q_points);

		const FEValuesExtractors::Scalar Density(dim);

		for (const auto &cell : semiconductor_dof_handler.active_cell_iterators())
	    {
			fe_values.reinit(cell);
			cell->get_dof_indices(local_dof_indices);

			// get the electron density values at the last time step
			fe_values[Density].get_function_values(
							electron_hole_pair.carrier_1.solution,
							carrier_1_density_values);

			// get the hole density values at the last time step
			fe_values[Density].get_function_values(
							electron_hole_pair.carrier_2.solution,
							carrier_2_density_values);

			donor_doping_profile.value_list(fe_values.get_quadrature_points(),
						donor_doping_values,
						dim); // calls the density values of the donor profile

			acceptor_doping_profile.value_list(fe_values.get_quadrature_points(),
						acceptor_doping_values,
						dim); // calls the density values of the acceptor profile

			for (unsigned int q_index = 0; q_index < n_q_points; ++q_index)
			{
				uncompendated_charge += std::abs(carrier_2_density_values[q_index]-
										 carrier_1_density_values[q_index]+
										 donor_doping_values[q_index]-
										 acceptor_doping_values[q_index])*
										fe_values.JxW(q_index);
			} //quad points
	    }// cells
		std::cout<< "\nChecking uncompensated charge!\n"
				 << "Full depletion approximation:       "
				 << (sim_params.scaled_n_type_depletion_width + sim_params.scaled_p_type_depletion_width)*
				    sim_params.scaled_domain_height * 1
				 << "\nCalculated uncompensated charge:    "
				 << uncompendated_charge;

		return uncompendated_charge/2;
	}

	template<int dim>
	void
	SolarCellProblem<dim>::
	print_results(unsigned int time_step_number)
	{
		// TODO: Make so that rescaled version will work in debug mode
		Threads::TaskGroup<void> task_group;
		task_group += Threads::new_task(&MixedPoisson::
						MixedFEM<dim>::output_rescaled_results,
						Mixed_Assembler,
						Poisson_dof_handler,
						Poisson_object.solution,
						sim_params,
						time_step_number);
			
		task_group += Threads::new_task(&LDG_System::
						LDG<dim>::output_rescaled_results,
						LDG_Assembler,
						semiconductor_dof_handler,
						electron_hole_pair,
						sim_params,
						time_step_number);

		task_group.join_all();
	} // print_results


	template<int dim>
	void
	SolarCellProblem<dim>::
	print_results_on_boundary(unsigned int time_step_number)
	{
		// TODO: Make so that rescaled version will work in debug mode
		Threads::TaskGroup<void> task_group;
		task_group += Threads::new_task(&LDG_System::
						LDG<dim>::output_unscaled_results_on_boundary,
						LDG_Assembler,
						semiconductor_dof_handler,
						electron_hole_pair,
		//				sim_params,
						time_step_number);


		task_group.join_all();
		}


	template<int dim>
	void
	SolarCellProblem<dim>::
	print_currents(unsigned int time_step_number)
	{

		const FESystem<dim> joint_fe(Poisson_fe, 1, carrier_fe, 1, carrier_fe, 1);
		DoFHandler<dim> joint_dof_handler(Poisson_triangulation);
		joint_dof_handler.distribute_dofs(joint_fe);
		std::cout << "Sanity check!\n"
				  << "Number of dofs in joint dof_handler:              "
				  << joint_dof_handler.n_dofs()
				  << "\n"
				  << "Number of sum of dofs in poisson and continuity:  "
				  << Poisson_dof_handler.n_dofs()+2*semiconductor_dof_handler.n_dofs()
				  << std::endl;
		Vector<double> joint_solution(joint_dof_handler.n_dofs());


		  {
			std::vector<types::global_dof_index> local_joint_dof_indices(
			  joint_fe.dofs_per_cell);
			std::vector<types::global_dof_index> local_poisson_dof_indices(
			  Poisson_fe.dofs_per_cell);
			std::vector<types::global_dof_index> local_carrier_dof_indices(
			  carrier_fe.dofs_per_cell);
			typename DoFHandler<dim>::active_cell_iterator
			  joint_cell       = joint_dof_handler.begin_active(),
			  joint_endc       = joint_dof_handler.end(),
			  poisson_cell     = Poisson_dof_handler.begin_active(),
			  semicon_cell     = semiconductor_dof_handler.begin_active();
			for (; joint_cell != joint_endc;
				 ++joint_cell, ++poisson_cell, ++semicon_cell)
			{
				  joint_cell->get_dof_indices(local_joint_dof_indices);
				  poisson_cell->get_dof_indices(local_poisson_dof_indices);
				  semicon_cell->get_dof_indices(local_carrier_dof_indices);
				  for (unsigned int i = 0; i < joint_fe.dofs_per_cell; ++i)
					if (joint_fe.system_to_base_index(i).first.first == 0)
					  {
						Assert(joint_fe.system_to_base_index(i).second < local_poisson_dof_indices.size(),ExcInternalError());

						joint_solution(local_joint_dof_indices[i]) = Poisson_object.solution(
									local_poisson_dof_indices[joint_fe.system_to_base_index(i).second]);
					  }
					else if(joint_fe.system_to_base_index(i).first.first == 1)
					  {
						Assert(joint_fe.system_to_base_index(i).first.first == 1,ExcInternalError());
						Assert(joint_fe.system_to_base_index(i).second <local_carrier_dof_indices.size(),ExcInternalError());

						joint_solution(local_joint_dof_indices[i]) = electron_hole_pair.carrier_1.solution(
							local_carrier_dof_indices[joint_fe.system_to_base_index(i).second]);
					  }
					else
					{
						Assert(joint_fe.system_to_base_index(i).first.first == 2,ExcInternalError());
						Assert(joint_fe.system_to_base_index(i).second <local_carrier_dof_indices.size(),ExcInternalError());

						joint_solution(local_joint_dof_indices[i]) = electron_hole_pair.carrier_2.solution(
							local_carrier_dof_indices[joint_fe.system_to_base_index(i).second]);
					}
			}
		  }


		DataOut<dim>  currents_out;
		PostProcessor_currents<dim> postprocessor_currents(sim_params,
														  electron_hole_pair.carrier_1.name.c_str(),
														  electron_hole_pair.carrier_2.name.c_str());

		currents_out.attach_dof_handler(joint_dof_handler);
		currents_out.add_data_vector(joint_solution,
									 postprocessor_currents);

		currents_out.build_patches();

		std::string currents_file = "Currents";
		currents_file += Utilities::int_to_string(time_step_number,3);
		//continuity_file += "grad.gp";
		currents_file += ".vtu";

		std::ofstream output_current(currents_file.c_str());
		currents_out.write_vtu(output_current);
//		data_out.write_gnuplot(output);
		output_current.close();
	}
	/*--------------------------------------------------------------------*/

	/*---------------------------------------------------------------------*/
	/*			RUN		  			       */
	/*---------------------------------------------------------------------*/

	template<int dim>
	void
	SolarCellProblem<dim>::
	run_full_system()
	{

		TimerOutput	timer(std::cout,
					TimerOutput::summary,
					TimerOutput::wall_times);

		timer.enter_subsection("Make Grids");
	
		// make the triangulations
		Grid_Maker::Grid<dim> grid_maker(sim_params);

		// make the grids and set the boundary conditions
		grid_maker.make_grids(semiconductor_triangulation,
					Poisson_triangulation
					);

		grid_maker.print_grid(Poisson_triangulation,"Grid.eps");
		grid_maker.print_grid(semiconductor_triangulation,"Semi.eps");

		timer.leave_subsection("Make Grids");

		// allocate the memory
		timer.enter_subsection("Allocate Memory");
		setup_dofs();
		timer.leave_subsection("Allocate Memory");

		
		timer.enter_subsection("Build Mappings");
		setup_mappings();
		timer.leave_subsection("Build Mappings");

		print_sim_info();

		// dont remove
		electron_hole_pair.penalty = 1.0e4;
	
		
		// assemble the global matrices
		timer.enter_subsection("Assemble Poisson Matrix");
		assemble_Poisson_matrix();
		timer.leave_subsection("Assemble Poisson Matrix");
	
		//assign delta_t
		delta_t = sim_params.delta_t;
		std::cout << "delta_t = " << delta_t << std::endl;

		// the 1.0 means we are calculating transients -> 1.0/delta_t
		timer.enter_subsection("Assemble LDG Matrices");
		assemble_LDG_system(1.0);
		timer.leave_subsection("Assemble LDG Matrices");

		// factor the matrices
		timer.enter_subsection("Factor Matrices");
		set_solvers();
		timer.leave_subsection("Factor Matrices");
	

		// make time stepping stuff... 
		unsigned int counter 	= 0;
		unsigned int time_step_number;
		double time;

		unsigned int number_outputs = sim_params.time_stamps;
		std::vector<double>	timeStamps(number_outputs);


		// if we are restarting the simulation, read in the dofs
		// from the end of the last simulation as the initial conditions of this one
		if(sim_params.restart_status)
		{
			unsigned int last_run_number;
			std::ifstream last_run;
			if(sim_params.restart_from_steady_state)
			{
				last_run.open("last_run_ss.txt");
			}
			else last_run.open("last_run.txt");

			if (last_run.is_open())
			{
			last_run >> last_run_number;
			std::cout<< "Next output file index will be: " <<last_run_number + 1 << std::endl;
			last_run.close();
			}
			else
			{
				std::cout << "UNABLE TO OPEN THE FILE. \n New outputs fill start from index 1!\n";
				last_run_number=1;
			}

			electron_hole_pair.read_dofs(sim_params.type_of_restart);
			for(unsigned int i=0; i<number_outputs; i++)
				timeStamps[i] = (i+1) * sim_params.t_end / number_outputs;
			// make the time stamps
			time_step_number	= last_run_number/*number_outputs*/;
			time 			= 0.0/*sim_params.t_end*/;
			std::cout << "running until t = " << sim_params.t_end << std::endl;
		}
		else // use defined initial condition functions
		{
			// Get the initial conditions
			VectorTools::project(semiconductor_dof_handler,
						electron_hole_pair.constraints,
						QGauss<dim>(degree+1),
						electrons_initial_condition,
						electron_hole_pair.carrier_1.solution);

/*			for (auto i: electron_hole_pair.carrier_1.solution)
			{
					std::cout << i << ' ' << std::endl;
			}*/

			VectorTools::project(semiconductor_dof_handler,
						electron_hole_pair.constraints,
						QGauss<dim>(degree+1),
						holes_initial_condition,
						electron_hole_pair.carrier_2.solution);

			// make the time stamps
			for(unsigned int i=0; i<number_outputs; i++)
				timeStamps[i] = (i+1) * sim_params.t_end / number_outputs;	

			time_step_number	= 0;
			time			= 0.0;
			std::cout << "running until t = " << sim_params.t_end << std::endl;
		}

		// get the intitial potential and electric field
		assemble_Poisson_rhs();
		solve_Poisson();
	
		// print the initial values
		print_results(time_step_number);

		// for testing convergence to steady state
		Convergence<dim> ConverganceCheck(Poisson_dof_handler,
									semiconductor_dof_handler,
										Poisson_object.solution,
										electron_hole_pair.carrier_1.solution,
										electron_hole_pair.carrier_2.solution);
		ConverganceCheck.print_indexes();

		// time stepping until the semiconductor converges at time 
		// t = t_end_1
		for(unsigned int k = 0;
				k < number_outputs;
				k++)
		{
			// while time < next print stamp time
			while(time < timeStamps[k])
			{
				timer.enter_subsection("Assemble semiconductor rhs");
				//std::cout<<"1"<<std::endl;
				assemble_semiconductor_rhs();
				timer.leave_subsection("Assemble semiconductor rhs");

				//std::cout<<"2"<<std::endl;
				timer.enter_subsection("Solve LDG Systems");
				solve_full_system();
				timer.leave_subsection("Solve LDG Systems");
	
				//std::cout<<"Jedziemy Poissona!"<<std::endl;
				timer.enter_subsection("Assemble Poisson rhs");
				assemble_Poisson_rhs();
				timer.leave_subsection("Assemble Poisson rhs");

				timer.enter_subsection("Solve Poisson system");
				solve_Poisson();			
				timer.leave_subsection("Solve Poisson system");

				time += delta_t;
				counter++;

			} // while

			// print the results
			timer.enter_subsection("Printing");
			time_step_number++;
			print_results(time_step_number);
			ConverganceCheck.calculate_residuum(Poisson_object.solution,
												electron_hole_pair.carrier_1.solution,
												electron_hole_pair.carrier_2.solution);
			ConverganceCheck.print_residuums(time_step_number);
			if(ConverganceCheck.get_steady_state_idicator())
			{
				std::cout << "STEADY STATE WAS ACHIEVED in "
						  << time_step_number
						  << " time step."
						  <<std::endl;
				break;
			}
			timer.leave_subsection("Printing");
		} // end for

		print_currents(time_step_number);
		calculate_uncompensated_charge();
		// print the dofs of the end state for restart situation.
		std::ofstream last_run;
		if(sim_params.calculate_steady_state)
		{
			last_run.open("last_run_ss.txt", std::ios::out | std::ios::trunc);
		}
		else last_run.open("last_run.txt", std::ios::out | std::ios::trunc);
		last_run << time_step_number;
		last_run.close();
		electron_hole_pair.print_dofs(sim_params.type_of_simulation);

	} // run


///////////////////////////////////////////////////////////////////////////////
//			TESTING ROUTINES
///////////////////////////////////////////////////////////////////////////////	

	/*----------------------------------------------------------------------*/
	// 	COUPLED MIXED-LDG TEST
	/*----------------------------------------------------------------------*/
	template<int dim>
	void
	SolarCellProblem<dim>::
	assemble_local_coupled_Poisson_test_rhs(
		const typename DoFHandler<dim>::active_cell_iterator & cell,
		Assembly::AssemblyScratch<dim>			 & scratch,
		Assembly::Poisson::CopyData<dim>		 & data,
		const double					 & time)	
	{
		const unsigned int dofs_per_cell = scratch.Poisson_fe_values.dofs_per_cell;

		const unsigned int n_q_points = scratch.Poisson_fe_values.n_quadrature_points;

		const unsigned int n_face_q_points =	
						scratch.Poisson_fe_face_values.n_quadrature_points;

		// SET RHS TIME
		Mixed_Assembler.test_coupling_Poisson_rhs.set_time(time);

		// Get the actual values for vector field and potential from FEValues
		// Use Extractors instead of having to deal with shapefunctions directly
		const FEValuesExtractors::Vector VectorField(0); // Vector as in Vector field
		const FEValuesExtractors::Scalar Potential(dim);
		const FEValuesExtractors::Scalar Density(dim);

		// reset the local_rhs vector to be zero
		data.local_rhs=0;

		// get the Poisson cell info
		std::pair<unsigned int, unsigned int> Poisson_cell_info = 
			s_2_p_map[std::pair<unsigned int, unsigned int>(cell->level(),																			cell->index())];

		// get the Poisson cell
		typename DoFHandler<dim>::active_cell_iterator
					Poisson_cell(&Poisson_triangulation,
						Poisson_cell_info.first,
						Poisson_cell_info.second,
						&Poisson_dof_handler);

		Poisson_cell->get_dof_indices(data.local_dof_indices);
		scratch.Poisson_fe_values.reinit(Poisson_cell);
		// get the test rhs for poisson
		// Assemble the right hand side on semiconductor side
		scratch.carrier_fe_values.reinit(cell);
			
		// get the electron density values at the previous time step
		scratch.carrier_fe_values[Density].get_function_values(
							electron_hole_pair.carrier_1.solution,
							scratch.old_carrier_1_density_values);
		// get the test rhs for poisson
		Mixed_Assembler.test_coupling_Poisson_rhs.value_list(
						scratch.Poisson_fe_values.get_quadrature_points(),
						scratch.Poisson_rhs_values);

		// Loop over all the quadrature points in this cell
		for(unsigned int q=0; q<n_q_points; q++)
		{
			// loop over the test function dofs for this cell
			for(unsigned int i=0; i<dofs_per_cell; i++)
			{
				// i-th potential basis functions at the point q
				const double psi_i_potential = 
						scratch.Poisson_fe_values[Potential].value(i,q);
		
				// get the local RHS values for this cell
				// = -int_{Omega_{e}} (C - u)
				data.local_rhs(i) += -psi_i_potential
					* ( scratch.Poisson_rhs_values[q]	
					- 	
					scratch.old_carrier_1_density_values[q]
					)	
					* scratch.Poisson_fe_values.JxW(q);
		
			} // for i
		} // for q

		// loop over all the faces of this cell to calculate the vector
		// from the dirichlet boundary conditions if the face is on the 
		// Dirichlet portion of the boundary
		for(unsigned int face_no=0;
				face_no<GeometryInfo<dim>::faces_per_cell;
				face_no++)
		{
			// obtain the face_no-th face of this cell
			typename DoFHandler<dim>::face_iterator face = Poisson_cell->face(face_no);

			// apply Dirichlet boundary conditions.. 
			// Since we are in a semicondcutor cell we know to apply the 
			// biases
			if((face->at_boundary()) && (face->boundary_id() == Dirichlet))
			{	
				// get the values of the shape functions at this boundary face
				scratch.Poisson_fe_face_values.reinit(Poisson_cell,face_no);
			
				// get the values of the dirichlet boundary conditions evaluated
				// on the quadrature points of this face	
				Mixed_Assembler.test_Poisson_bc.value_list(
						scratch.Poisson_fe_face_values.get_quadrature_points(),
						scratch.Poisson_bc_values);
		

				// loop over all the quadrature points of this face
				for(unsigned int q=0; q<n_face_q_points; q++)
				{
					// loop over all the test function dofs of this face
					for(unsigned int i=0; i<dofs_per_cell; i++)
					{
						// - \int_{face} p * n * (phi_{Dichlet}) dx
						data.local_rhs(i) +=	
								-(scratch.Poisson_fe_face_values[VectorField].value(i,q) *
								scratch.Poisson_fe_face_values.normal_vector(q) *
								scratch.Poisson_bc_values[q] *
								scratch.Poisson_fe_face_values.JxW(q));
					} // for i
				} // for q
			} // end if
		} // end for face_no

	} // assemble_local_coupled_Poisson_test_rhs

	template<int dim>
	void
	SolarCellProblem<dim>::
	assemble_local_coupled_DD_test_rhs(
				const typename DoFHandler<dim>::active_cell_iterator & cell,
				Assembly::AssemblyScratch<dim>			 & scratch,
				Assembly::DriftDiffusion::CopyData<dim>		 & data,
				const double					 & time,
				const double					 & penalty)	
	{

		const unsigned int dofs_per_cell = scratch.carrier_fe_values.dofs_per_cell;
		const unsigned int n_q_points = scratch.carrier_fe_values.n_quadrature_points;
		const unsigned int n_face_q_points = 
				scratch.carrier_fe_face_values.n_quadrature_points;


		cell->get_dof_indices(data.local_dof_indices);

		// get the Poisson cell info
		std::pair<unsigned int, unsigned int> Poisson_cell_info = 
					s_2_p_map[std::pair<unsigned int, unsigned int>(cell->level(),
											cell->index())];

		// get the Poisson cell
		typename DoFHandler<dim>::active_cell_iterator
					Poisson_cell(&Poisson_triangulation,
						Poisson_cell_info.first,
						Poisson_cell_info.second,
						&Poisson_dof_handler);

		scratch.Poisson_fe_values.reinit(Poisson_cell);

		scratch.carrier_fe_values.reinit(cell);
		
		// reset the local_rhs to be zero
		data.local_carrier_1_rhs=0;
		data.local_carrier_2_rhs=0;
	
		const FEValuesExtractors::Vector Current(0);
		const FEValuesExtractors::Scalar Density(dim);
		const FEValuesExtractors::Vector ElectricField(0);

		LDG_Assembler.test_DD_rhs.set_time(time);
		LDG_Assembler.test_DD_bc.set_time(time);

		// rhs values
		LDG_Assembler.test_DD_rhs.value_list(
				scratch.carrier_fe_values.get_quadrature_points(),
				scratch.generation_values);
	
		// get the values of carrier_1 density
		scratch.carrier_fe_values[Density].get_function_values(
				electron_hole_pair.carrier_1.solution,															scratch.old_carrier_1_density_values);


		// get the electric field values at the previous time step
		scratch.Poisson_fe_values[ElectricField].get_function_values(
					Poisson_object.solution,																	scratch.electric_field_values);


		double h = cell->diameter();
		// loop over all the quadrature points in this cell and compute body integrals
		for(unsigned int q=0; q<n_q_points; q++)
		{
			// loop over all the test function dofs and get the test functions
			for(unsigned int i=0; i<dofs_per_cell; i++)
			{
				const double psi_i_density = 
							scratch.carrier_fe_values[Density].value(i,q);
		
				const Tensor<1,dim> psi_i_field	=
							scratch.carrier_fe_values[Current].value(i,q);
					
				// contribution from RHS function + Drift
				// int_{Omega} v * R +  p * E * u dx
				data.local_carrier_1_rhs(i) += ( 
					(psi_i_density * scratch.generation_values[q])
					-
					psi_i_field *
					scratch.electric_field_values[q] *
					scratch.old_carrier_1_density_values[q]
					) *
					scratch.carrier_fe_values.JxW(q);

			} // for i
		}	// for q

		// loop over all the faces of this cell and compute the contribution 
		// from the boundary conditions
		for(unsigned int face_no=0; 
		    face_no< GeometryInfo<dim>::faces_per_cell;
		    face_no++)
		{
			// get the face_no-th face of this cell
			typename DoFHandler<dim>::face_iterator face = cell->face(face_no);
			
			// if on boundary apply boundayr conditions	
			if(face->at_boundary() )
			{
				// reinitialize the fe_face_values for this cell ONLY if it is as the
				//boundary otherwise its a waste.  then assemble the appropriate
				//boundary conditions
				scratch.carrier_fe_face_values.reinit(cell, face_no);
			
				if(face->boundary_id() == Dirichlet)
				{
					// get the density
					LDG_Assembler.test_DD_bc.value_list(
								scratch.carrier_fe_face_values.get_quadrature_points(),
								scratch.carrier_1_bc_values);

					// loop over all the quadrature points on this face
					for(unsigned int q=0; q<n_face_q_points; q++)
					{
						// loop over all the test function dofs on this face
						for(unsigned int i=0; i<dofs_per_cell; i++)
						{
							// get the test function
							const Tensor<1, dim>  psi_i_field = 
												scratch.carrier_fe_face_values[Current].value(i,q);
							const double 	psi_i_density = 
												scratch.carrier_fe_face_values[Density].value(i,q);
 
							// int_{\Gamma_{D}} -p^{-} n^{-} u_{D} ds
							data.local_carrier_1_rhs(i) += 
									(-1.0 * psi_i_field *
									 scratch.carrier_fe_face_values.normal_vector(q) 
									+ 
									(penalty/h) * psi_i_density) * 
									scratch.carrier_1_bc_values[q] *
									scratch.carrier_fe_face_values.JxW(q);	
						} // for i
					}	// for q
				} // end Dirichlet
				else if(face->boundary_id() == 666/*Interface*/)
				{
				}
				else if(face->boundary_id() == Neumann)
				{
					// NOTHIN TO DO IF INSULATING
				}
				else
					Assert(false, ExcNotImplemented() );
			} // end at boundary
		} // for face_no
	
	} // assemble_local_coupling_DD_rhs

	/*-----------------------------------------------------------------------*/
	/*			INTERFACE ASSMBLY ROUTINES			 */
	/*-----------------------------------------------------------------------*/

	template<int dim>
	void 
	SolarCellProblem<dim>::
	assemble_local_test_semiconductor_rhs(
			const typename DoFHandler<dim>::active_cell_iterator & cell,
			Assembly::AssemblyScratch<dim>			     & scratch,
			Assembly::DriftDiffusion::CopyData<dim>		     & data,
			const double 					     & time,
			const double 					     & penalty)									
	{
		const unsigned int dofs_per_cell  = scratch.carrier_fe_values.dofs_per_cell;
		const unsigned int n_q_points = scratch.carrier_fe_values.n_quadrature_points;
		const unsigned int n_face_q_points = 
					scratch.carrier_fe_face_values.n_quadrature_points;

		cell->get_dof_indices(data.local_dof_indices);

		// reinitialize the fe_values 
		scratch.carrier_fe_values.reinit(cell);
		
		// reset the local_rhs to be zero
		data.local_carrier_1_rhs=0;
		data.local_carrier_2_rhs=0;
	
		const FEValuesExtractors::Vector Current(0);
		const FEValuesExtractors::Scalar Density(dim);
		Tensor<1,dim>	field_values;
		field_values[0]	= 0.0; 
		field_values[1] = 0.0; 

		LDG_Assembler.test_interface_rhs.set_time(time);
		LDG_Assembler.test_interface_bc.set_time(time);
		LDG_Assembler.test_interface_function.set_time(time);

		// the test version	
		LDG_Assembler.test_interface_rhs.value_list(
						scratch.carrier_fe_values.get_quadrature_points(),
						scratch.generation_values);
	
		// get the values of carrier_1 and carrier_2 densities at the pevious time step
		scratch.carrier_fe_values[Density].get_function_values(
						electron_hole_pair.carrier_1.solution,
						scratch.old_carrier_1_density_values);


		double h = cell->diameter();
		// loop over all the quadrature points in this cell and compute body integrals
		for(unsigned int q=0; q<n_q_points; q++)
		{
			// loop over all the test function dofs and get the test functions
			for(unsigned int i=0; i<dofs_per_cell; i++)
			{
				const double psi_i_density = scratch.carrier_fe_values[Density].value(i,q);
		
				const Tensor<1,dim> psi_i_field	 = scratch.carrier_fe_values[Current].value(i,q);
					
				// contribution from RHS function + Drift
				// int_{Omega} v * R +  p * E * u dx
				data.local_carrier_1_rhs(i) += ( 
						(psi_i_density * scratch.generation_values[q])
						-
						psi_i_field *
						field_values *
						scratch.old_carrier_1_density_values[q]
						) *
						scratch.carrier_fe_values.JxW(q);

			} // for i
		}	// for q

		// loop over all the faces of this cell and compute the contribution 
		// from the boundary conditions
		for(unsigned int face_no=0; 
				face_no< GeometryInfo<dim>::faces_per_cell;
				face_no++)
		{
			// get the face_no-th face of this cell
			typename DoFHandler<dim>::face_iterator 	face = cell->face(face_no);
			
			// if on boundary apply boundayr conditions	
			if(face->at_boundary() )
			{
				// reinitialize the fe_face_values for this cell ONLY if it is as the
				//boundary otherwise its a waste.  then assemble the appropriate
				//boundary conditions
				scratch.carrier_fe_face_values.reinit(cell, face_no);
			
				if(face->boundary_id() == Dirichlet)
				{
					// get the density
					LDG_Assembler.test_interface_bc.value_list(
							scratch.carrier_fe_face_values.get_quadrature_points(),
							scratch.carrier_1_bc_values);

					// loop over all the quadrature points on this face
					for(unsigned int q=0; q<n_face_q_points; q++)
					{
						// loop over all the test function dofs on this face
						for(unsigned int i=0; i<dofs_per_cell; i++)
						{
							// get the test function
							const Tensor<1, dim>  psi_i_field = 
									scratch.carrier_fe_face_values[Current].value(i,q);
							const double psi_i_density = 
									scratch.carrier_fe_face_values[Density].value(i,q);
 
							// int_{\Gamma_{D}} -p^{-} n^{-} u_{D} ds
							data.local_carrier_1_rhs(i) += 
								(-1.0 * psi_i_field *
								scratch.carrier_fe_face_values.normal_vector(q) 
								+ 
								(penalty/h) * psi_i_density) * 
								scratch.carrier_1_bc_values[q] *
							 	scratch.carrier_fe_face_values.JxW(q);				
						} // for i
					}	// for q
				} // end Dirichlet
				else if(face->boundary_id() == Neumann)
				{
					// NOTHIN TO DO IF INSULATING
				}
				else
					Assert(false, ExcNotImplemented() );
			} // end at boundary
		} // for face_no
	}


	/*---------------------------------------------------------------------*/
	/*			STEADY STATE TEST					 */
	/*---------------------------------------------------------------------*/

	template<int dim>
	void
	SolarCellProblem<dim>::
	test_steady_state(const unsigned int & n_refine,
			ConvergenceTable & Mixed_table,
			ConvergenceTable & LDG_table)
	{
		full_system = false;

		double 	primary_error, flux_error;
		sim_params.set_params_for_testing(n_refine);

		Grid_Maker::Grid<dim> grid_maker(sim_params);
		grid_maker.make_test_grid(Poisson_triangulation, n_refine);
		grid_maker.make_test_grid(semiconductor_triangulation, n_refine);
//		grid_maker.refine_test_grid(Poisson_triangulation,2);
//		grid_maker.refine_test_grid(semiconductor_triangulation, 2);

		setup_dofs();
		electron_hole_pair.set_semiconductor_for_testing(sim_params);

		assemble_Poisson_matrix();
		assemble_LDG_system(0.0);

		set_solvers();

		unsigned int n_active_cells = Poisson_triangulation.n_active_cells();
		unsigned int n_dofs		= Poisson_dof_handler.n_dofs();
		double h =  GridTools::maximal_cell_diameter(Poisson_triangulation);
		Mixed_table.add_value("h", h);		
		Mixed_table.add_value("cells", n_active_cells);
		Mixed_table.add_value("dofs", n_dofs);

		n_active_cells = semiconductor_triangulation.n_active_cells();
		n_dofs		 = semiconductor_dof_handler.n_dofs();

		LDG_table.add_value("h", h);		
		LDG_table.add_value("cells", n_active_cells);
		LDG_table.add_value("dofs", n_dofs);
		
		// assemble poisson rhs
		WorkStream::run(Poisson_dof_handler.begin_active(),
				Poisson_dof_handler.end(),
				std_cxx11::bind(&MixedPoisson::MixedFEM<dim>::
						assemble_local_test_rhs,
						Mixed_Assembler, // the object
						std_cxx11::_1,
						std_cxx11::_2,
						std_cxx11::_3),
				std_cxx11::bind(&SOLARCELL::SolarCellProblem<dim>::
						copy_local_to_global_Poisson_rhs,
						this, // this object
						std_cxx11::_1),
				Assembly::AssemblyScratch<dim>(Poisson_fe,
								 carrier_fe,
								 QGauss<dim>(degree+2),
								 QGauss<dim-1>(degree+2)),
				Assembly::Poisson::CopyData<dim>(Poisson_fe)
				);

		//assemble ldg rhs
		WorkStream::run(semiconductor_dof_handler.begin_active(),
				semiconductor_dof_handler.end(),
				std_cxx11::bind(&LDG_System::LDG<dim>::
						assemble_local_test_rhs,
						LDG_Assembler, // Assembler object
						std_cxx11::_1,
						std_cxx11::_2,
						std_cxx11::_3,
						electron_hole_pair.penalty),
				std_cxx11::bind(&SOLARCELL::SolarCellProblem<dim>::
						copy_local_to_global_semiconductor_system_rhs,
						this, // this object
						std_cxx11::_1),
				Assembly::AssemblyScratch<dim>(Poisson_fe,																		 carrier_fe,
								 QGauss<dim>(degree+2),
								 QGauss<dim-1>(degree+2)),
								 Assembly::DriftDiffusion::CopyData<dim>(carrier_fe)
								 );

		solve_Poisson();
		solve_semiconductor_system();

		Mixed_Assembler.compute_errors(Poisson_triangulation, 
						Poisson_dof_handler,
						Poisson_object.solution,
						primary_error,
						flux_error);
	
		Mixed_table.add_value("Phi", primary_error);
		Mixed_table.add_value("D", flux_error);

	
		LDG_Assembler.compute_errors(semiconductor_triangulation,
						semiconductor_dof_handler,
						electron_hole_pair.carrier_1.solution,
						primary_error,
						flux_error,
						true,
						0.0); // no time

		LDG_table.add_value("u", primary_error);
		LDG_table.add_value("J", flux_error);
		
		Mixed_Assembler.output_unscaled_results(Poisson_dof_handler,
							Poisson_object.solution,
							n_refine);

		LDG_Assembler.output_unscaled_results(semiconductor_dof_handler,
							electron_hole_pair,
							n_refine);
	}

	/*---------------------------------------------------------------------*/
	/*			LDG-IMEX TEST						*/
	/*---------------------------------------------------------------------*/
	template<int dim>
	void
	SolarCellProblem<dim>::
	test_transient(const unsigned int & n_refine,
		  	 ConvergenceTable   & LDG_table)
	{
		full_system = false;
		double primary_error, flux_error;

		sim_params.set_params_for_testing(n_refine);

		Grid_Maker::Grid<dim> grid_maker(sim_params);
		grid_maker.make_test_grid(Poisson_triangulation, n_refine);
		grid_maker.make_test_tran_grid(semiconductor_triangulation, n_refine);
	//	grid_maker.refine_test_grid(semiconductor_triangulation, 2);

		setup_dofs();

		double h = GridTools::maximal_cell_diameter(semiconductor_triangulation);

		electron_hole_pair.set_semiconductor_for_testing(sim_params);
 		
		delta_t = 1.0;
		for(unsigned int i=0; i < carrier_fe.degree+1; i++)
			delta_t *= h;

		std::cout << "dt = " << delta_t << std::endl;

		assemble_LDG_system(1.0);

		electron_hole_pair.carrier_1.set_solver();

		// Get the initial conditions
		VectorTools::project(semiconductor_dof_handler,
					electron_hole_pair.constraints,
					QGauss<dim>(degree+1),
					LDG_Assembler.test_interface_initial,
					electron_hole_pair.carrier_1.solution);

		double t_end = 1.0;
		double time  = 0.0;
		
		while(time < t_end)
		{
			// set carrier_system_rhs = M * u^{n-1}
			electron_hole_pair.mass_matrix.vmult(electron_hole_pair.carrier_1.system_rhs,
												electron_hole_pair.carrier_1.solution);
			// set rhs for time step n-1
			WorkStream::run(semiconductor_dof_handler.begin_active(),
					  semiconductor_dof_handler.end(),
					  std_cxx11::bind(&LDG_System::LDG<dim>::
							assemble_local_test_transient_rhs,
							LDG_Assembler, // Assembler object
							std_cxx11::_1,
							std_cxx11::_2,
							std_cxx11::_3,
							electron_hole_pair.carrier_1.solution,
							time,
							electron_hole_pair.penalty),
					  std_cxx11::bind(&SOLARCELL::SolarCellProblem<dim>::
							copy_local_to_global_semiconductor_system_rhs,
							this, // this object
							std_cxx11::_1),
					  Assembly::AssemblyScratch<dim>(Poisson_fe,
									carrier_fe,
									QGauss<dim>(degree+2),
									QGauss<dim-1>(degree+2)),
					  Assembly::DriftDiffusion::CopyData<dim>(carrier_fe)
					 );
			electron_hole_pair.carrier_1.solve();
	
			// UPDATE TIME
			time += delta_t;
		}	// end while																		


		// LDG ERRORS
		LDG_Assembler.compute_errors(semiconductor_triangulation,
						 semiconductor_dof_handler,
						 electron_hole_pair.carrier_1.solution,
						 primary_error,
						 flux_error,
						 false,
						 time);

		unsigned int n_active_cells = semiconductor_triangulation.n_active_cells();
		unsigned int n_dofs		= semiconductor_dof_handler.n_dofs();

	
		LDG_table.add_value("h",h);
		LDG_table.add_value("cells", n_active_cells);
		LDG_table.add_value("dofs", n_dofs);

		LDG_table.add_value("u", primary_error);
		LDG_table.add_value("J", flux_error);
		
		LDG_Assembler.output_unscaled_results(semiconductor_dof_handler,
											electron_hole_pair,
											n_refine);

	} // end test_transient


	/*---------------------------------------------------------------------*/
	/*		LDG-MIXED TEST						 */
	/*---------------------------------------------------------------------*/

	template<int dim>
	void
	SolarCellProblem<dim>::
	test_DD_Poisson(const unsigned int & n_refine,
			ConvergenceTable	& Mixed_table,
			ConvergenceTable	& LDG_table)
	{
		full_system = false;
		double 	primary_error, flux_error;

		sim_params.set_params_for_testing(n_refine);

		Grid_Maker::Grid<dim> grid_maker(sim_params);
		grid_maker.make_DD_Poisson_grid(Poisson_triangulation, n_refine);
		grid_maker.make_DD_Poisson_grid(semiconductor_triangulation, n_refine);
	//	grid_maker.refine_test_grid(semiconductor_triangulation, 2);

		setup_dofs();
		setup_mappings();

		double h = GridTools::maximal_cell_diameter(semiconductor_triangulation);

		electron_hole_pair.set_semiconductor_for_testing(sim_params);
 		
		delta_t = 1.0;
		for(unsigned int i=0; i < carrier_fe.degree+1; i++)
			delta_t *= h;

		std::cout << "dt = " << delta_t << std::endl;

		assemble_Poisson_matrix();
		assemble_LDG_system(1.0);

		electron_hole_pair.carrier_1.set_solver();

		// Get the initial conditions
		VectorTools::project(semiconductor_dof_handler,
					electron_hole_pair.constraints,
					QGauss<dim>(degree+1),
					LDG_Assembler.test_interface_initial,
					electron_hole_pair.carrier_1.solution);
		
		set_solvers();

		double t_end = 1.0;
		double time  = 0.0;
	
		while(time < t_end)
		{
			///////////////////////////////////////
			// UPDATE POISSON RHS AND SOLVE
			/////////////////////////////////////////
			Poisson_object.system_rhs = 0;

			WorkStream::run(semiconductor_dof_handler.begin_active(),
					semiconductor_dof_handler.end(),
					std_cxx11::bind(&SOLARCELL::SolarCellProblem<dim>::
							assemble_local_coupled_Poisson_test_rhs,
							this, // this object
							std_cxx11::_1,
							std_cxx11::_2,
							std_cxx11::_3,
							time),
					std_cxx11::bind(&SOLARCELL::SolarCellProblem<dim>::
							copy_local_to_global_Poisson_rhs,
							this, // this object
							std_cxx11::_1),
					Assembly::AssemblyScratch<dim>(Poisson_fe,
									 carrier_fe,
									 QGauss<dim>(degree+2),
									 QGauss<dim-1>(degree+2)),
					Assembly::Poisson::CopyData<dim>(Poisson_fe)
					);
	
			solve_Poisson();

			///////////////////////////////////////
			// UPDATE DD RHS AND SOLVE
			/////////////////////////////////////////
			// set carrier_system_rhs = M * u^{n-1}
			electron_hole_pair.mass_matrix.vmult(electron_hole_pair.carrier_1.system_rhs,
								  electron_hole_pair.carrier_1.solution);
		
			// set rhs for time step n-1
			WorkStream::run(semiconductor_dof_handler.begin_active(),
					  semiconductor_dof_handler.end(),
					  std_cxx11::bind(&SOLARCELL::SolarCellProblem<dim>::
							   assemble_local_coupled_DD_test_rhs,		
							   this, // this object
							   std_cxx11::_1,
							   std_cxx11::_2,
							   std_cxx11::_3,
							   time,
							   electron_hole_pair.penalty),
					  std_cxx11::bind(&SOLARCELL::SolarCellProblem<dim>::
							    copy_local_to_global_semiconductor_system_rhs,
							    this, // this object
							    std_cxx11::_1),
					  Assembly::AssemblyScratch<dim>(Poisson_fe,
									    carrier_fe,
									    QGauss<dim>(degree+2),
									    QGauss<dim-1>(degree+2)),
					  Assembly::DriftDiffusion::CopyData<dim>(carrier_fe)
					);

			electron_hole_pair.carrier_1.solve();

			// UPDATE TIME
			time += delta_t;
		}	// end while																		

		// LDG ERRORS
		LDG_Assembler.compute_coupled_errors(semiconductor_triangulation,
							semiconductor_dof_handler,
							electron_hole_pair.carrier_1.solution,
							primary_error,
							flux_error,
							time);

		unsigned int n_active_cells = semiconductor_triangulation.n_active_cells();
		unsigned int n_dofs	   	= semiconductor_dof_handler.n_dofs();

		LDG_table.add_value("h",h);
		LDG_table.add_value("cells", n_active_cells);
		LDG_table.add_value("dofs", n_dofs);

		LDG_table.add_value("u", primary_error);
		LDG_table.add_value("J", flux_error);
		
		// MIXED ERRORS
		Mixed_Assembler.compute_errors(Poisson_triangulation, 
					Poisson_dof_handler,
					Poisson_object.solution,
					primary_error,
					flux_error);
		
		n_active_cells = Poisson_triangulation.n_active_cells();
		n_dofs  = Poisson_dof_handler.n_dofs();

		Mixed_table.add_value("h",h);
		Mixed_table.add_value("cells", n_active_cells);
		Mixed_table.add_value("dofs", n_dofs);

		Mixed_table.add_value("Phi", primary_error);
		Mixed_table.add_value("D", flux_error);
	}

} // namespace SOLARCELL
 // namespace SOLARCELL
