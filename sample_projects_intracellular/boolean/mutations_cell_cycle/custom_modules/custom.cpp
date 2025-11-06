/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "custom.h"
#include "../BioFVM/BioFVM.h"  
#include "../addons/PhysiBoSS/src/maboss_intracellular.h"
using namespace BioFVM;

// declare cell definitions here 

std::vector<bool> nodes;

void create_cell_types( void )
{
	// set the random seed 
	if (parameters.ints.find_index("random_seed") != -1)
	{
		SeedRandom(parameters.ints("random_seed"));
	}
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 

	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 

	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = NULL;
	cell_defaults.functions.update_phenotype = NULL; 
	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.pre_update_intracellular = pre_update_intracellular; 
	cell_defaults.functions.post_update_intracellular = post_update_intracellular; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	cell_defaults.custom_data.add_variable(parameters.strings("node_to_visualize"), "dimensionless", 0.0 ); //for paraview visualization
	
	// Add apoptosis timer variables
	cell_defaults.custom_data.add_variable("apoptosis_commitment", "dimensionless", 0.0);
	cell_defaults.custom_data.add_variable("apoptosis_threshold", "dimensionless", 10.0);

	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 
	
	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 

	get_cell_definition("default").functions.update_phenotype = wt_phenotype; 

	// Apply custom functions to specific cell types
    std::vector<std::string> cell_types_with_custom_functions = {"default"};
    
    for (const auto& cell_type : cell_types_with_custom_functions)
    {
        Cell_Definition* pCD = find_cell_definition(cell_type);
        if(pCD)
        {
            Cycle_Model* cycle_model = &(pCD->phenotype.cycle.model());
            
            // Apply custom phase entry functions for Flow Cytometry model
            if (cycle_model->code == PhysiCell_constants::flow_cytometry_separated_cycle_model)
            {
                // Find phase indices
                int G0G1_index = cycle_model->find_phase_index(PhysiCell_constants::G0G1_phase);
                int S_index = cycle_model->find_phase_index(PhysiCell_constants::S_phase);
                int G2_index = cycle_model->find_phase_index(PhysiCell_constants::G2_phase);
                int M_index = cycle_model->find_phase_index(PhysiCell_constants::M_phase);
                
                // Assign custom entry functions
                if (G0G1_index >= 0)
                    cycle_model->phases[G0G1_index].entry_function = custom_G0G1_phase_entry_function;
                if (S_index >= 0)
                    cycle_model->phases[S_index].entry_function = custom_S_phase_entry_function;
                if (G2_index >= 0)
                    cycle_model->phases[G2_index].entry_function = custom_G2_phase_entry_function;
                if (M_index >= 0)
                    cycle_model->phases[M_index].entry_function = custom_M_phase_entry_function;
                
                std::cout << "Applied custom phase entry functions to " << cell_type << " cells (Flow Cytometry model)." << std::endl;
            }
            else
            {
                // For other models, apply the mutation exit function to the live phase
                int live_phase_index = cycle_model->find_phase_index(PhysiCell_constants::live);
                if (live_phase_index >= 0)
                {
                    cycle_model->phase_link(live_phase_index, live_phase_index).exit_function = phase_exit_mutation_function;
                    std::cout << "Applied mutation exit function to " << cell_type << " cells." << std::endl;
                }
            }
        }
        else
        {
            std::cout << "Warning: Cell type " << cell_type << " not found." << std::endl;
        }
    }
	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/

	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries();

	/*
	   This summarizes the setup. 
	*/
	
	display_cell_definitions( std::cout ); 


	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
	// load cells from your CSV file
	load_cells_from_pugixml(); 	
}

void pre_update_intracellular( Cell* pCell, Phenotype& phenotype, double dt )
{
	// Update GF Boolean network nodes based on internalized GF concentration
	// This must happen before the Boolean network runs
	update_gf_boolean_nodes(pCell);
	
	// Commented out to avoid crash - parameter "$time_scale" doesn't exist in the network -- SHOULD BE FIXED
	// if (PhysiCell::PhysiCell_globals.current_time >= 100.0 
	// 	&& pCell->phenotype.intracellular->get_parameter_value("$time_scale") == 0.0
	// ){
	// 	pCell->phenotype.intracellular->set_parameter_value("$time_scale", 0.1);
	// }

}



void post_update_intracellular( Cell* pCell, Phenotype& phenotype, double dt )
{
	color_node(pCell);

	std::cout << "Cell " << pCell->ID << " generation " << pCell->generation << " parent " << pCell->parent_ID << std::endl;

	// Decay apoptosis commitment if no caspase signals are active
	if (pCell->phenotype.intracellular &&
		pCell->phenotype.intracellular->intracellular_type == "maboss")
	{
		MaBoSSIntracellular* maboss_model = static_cast<MaBoSSIntracellular*>(pCell->phenotype.intracellular);
		bool Casp3 = maboss_model->maboss.get_node_value("Casp3");
		bool Casp8 = maboss_model->maboss.get_node_value("Casp8");
		bool Casp9 = maboss_model->maboss.get_node_value("Casp9");
		
		// Only decay if no caspase signals are active
		if (!Casp3 && !Casp8 && !Casp9)
		{
			pCell->custom_data["apoptosis_commitment"] *= 0.95; // 5% decay per time step
		}
	}

	// ADD NODE C RESISTANCE
	// update_cell_from_boolean_model(pCell, phenotype, dt);
	
}

void update_cell_from_boolean_model(Cell* pCell, Phenotype& phenotype, double dt)
{
	static int death_decay_idx = pCell->custom_data.find_variable_index( "death_commitment_decay" );
	static int apoptosis_index = phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model );
	static float apoptosis_rate = pCell->custom_data["apoptosis_rate"];
	static float death_commitment_decay = pCell->custom_data["death_decay_idx"];
	bool S_entry_node = pCell->phenotype.intracellular->get_boolean_variable_value( "S_entry" );
	bool G2M_entry_node = pCell->phenotype.intracellular->get_boolean_variable_value( "G2M_entry" );
	bool G0G1_entry_node = pCell->phenotype.intracellular->get_boolean_variable_value( "G0G1_entry" );

	static int density_idx = microenvironment.find_density_index("drug");
	double drug_density_ext = pCell->nearest_density_vector()[density_idx]; // A density (mM)
	double density_int = pCell->phenotype.molecular.internalized_total_substrates[density_idx];
	// density_int /= cell_volume; // divide int tot substrate to get concentration


	float basal_apoptosis = 1e-06;
	float final_apoptosis_rate = basal_apoptosis + apoptosis_rate*drug_density_ext;



	if ( S_entry_node || G2M_entry_node || G0G1_entry_node )
	{
		pCell-> phenotype.death.rates[apoptosis_index] = basal_apoptosis;
	}
	else
	{
		pCell-> phenotype.death.rates[apoptosis_index] = final_apoptosis_rate;
	}

	return;
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
	std::vector< std::string > output( 4 , "rgb(0,0,0)" );
	
	if ( !pCell->phenotype.intracellular->get_boolean_variable_value( parameters.strings("node_to_visualize") ) )
	{
		output[0] = "rgb(255,0,0)";
		output[2] = "rgb(125,0,0)";
		
	}
	else{
		output[0] = "rgb(0, 255,0)";
		output[2] = "rgb(0, 125,0)";
	}
	
	return output;
}

void color_node(Cell* pCell){
	std::string node_name = parameters.strings("node_to_visualize");
	pCell->custom_data[node_name] = pCell->phenotype.intracellular->get_boolean_variable_value(node_name);
}

void wt_phenotype( Cell* pCell, Phenotype& phenotype, double dt )
{ 
	// Default phenotype, doesn't do anything special

	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL; 
		return; 
	}

	// Add info on parent cell ID and generation
	// change_custom_data_var(pCell, "parent_ID", pCell->parent_ID);
	// change_custom_data_var(pCell, "generation", pCell->generation);

	// First check O2 availability
	// update_cell_and_death_parameters_O2_based(pCell, phenotype, dt);

	return; 
}

void add_custom_cycle_function()
{

	for( int i=0; i < all_cells->size() ; i++ )
	{
		// to call each cell, the pointer is *all_cells)[i] instead of pCell
		(*all_cells)[i]->phenotype.cycle.model().phases[0].entry_function = phase_exit_mutation_function;
	}
	
	return;
}


void phase_exit_mutation_function( Cell* pCell, Phenotype& phenotype, double dt )
{

	// Choose a random node
	static std::default_random_engine generator;
    
	
	// Only proceed if the cell has a MaBoSS model
    if (pCell->phenotype.intracellular &&
        pCell->phenotype.intracellular->intracellular_type == "maboss")
    {
        MaBoSSIntracellular* maboss_model = static_cast<MaBoSSIntracellular*>(pCell->phenotype.intracellular);
        std::vector<std::string> node_names = maboss_model->maboss.get_all_node_names();

        if (!node_names.empty())
        {			
			if (uniform_random() < parameters.doubles("mutation_rate_threshold")){

				std::uniform_int_distribution<size_t> node_dist(0, node_names.size() - 1);
				size_t node_idx = node_dist(generator);
				std::string node_name = node_names[node_idx];

				// This is where we list the nodes that we want to mutate

				if (node_name == "S_entry" || node_name == "G2M_entry" || node_name == "G0G1_entry"){

					// Flip its value
					bool current_value = maboss_model->maboss.get_node_value(node_name);
					bool new_value = !current_value;
					maboss_model->maboss.set_node_value(node_name, new_value);

					// Record the mutation in the cell's mutations vector
					std::string mutation_record = node_name + "_" + std::to_string(new_value);
					pCell->custom_data.mutations.push_back(mutation_record);

					// std::cout << "Cell " << pCell->ID << ": flipped node " << node_name
					//           << " from " << current_value << " to " << new_value << std::endl;

				}
			}

        }
        else
        {
            std::cout << "Cell " << pCell->ID << " has no nodes in its MaBoSS model!" << std::endl;
        }
    }
    else
    {
        std::cout << "Cell " << pCell->ID << " does not have a MaBoSS model!" << std::endl;
    }
}


double get_custom_data_variable(Cell* pCell, std::string variable_name){
	int tmp_variable_idx = pCell->custom_data.find_variable_index(variable_name);
	double tmp_variable_value = pCell->custom_data[tmp_variable_idx];
	return tmp_variable_value;
}

void change_custom_data_var(Cell* pCell, std::string variable_name, double variable_new_value){
	int tmp_variable_idx = pCell->custom_data.find_variable_index(variable_name);
	pCell->custom_data[tmp_variable_idx] = variable_new_value;
	return;
}

// Helper function to check Boolean network for quiescence signals
bool check_boolean_network_quiescence( Cell* pCell )
{
	// Only proceed if the cell has a MaBoSS model
	if (pCell->phenotype.intracellular &&
		pCell->phenotype.intracellular->intracellular_type == "maboss")
	{
		MaBoSSIntracellular* maboss_model = static_cast<MaBoSSIntracellular*>(pCell->phenotype.intracellular);
		
		// Check cyclin levels for quiescence - based on Sizek Boolean network logic
		bool CyclinD1 = maboss_model->maboss.get_node_value("CyclinD1");
		bool CyclinA = maboss_model->maboss.get_node_value("CyclinA");
		bool CyclinE = maboss_model->maboss.get_node_value("CyclinE");
		bool CyclinB = maboss_model->maboss.get_node_value("CyclinB");
		
		// Check apoptosis signals - if any caspase is active, cell should not cycle
		bool Casp3 = maboss_model->maboss.get_node_value("Casp3");
		bool Casp8 = maboss_model->maboss.get_node_value("Casp8");
		bool Casp9 = maboss_model->maboss.get_node_value("Casp9");
		
		// If apoptosis is active, cell should not cycle
		if (Casp3 || Casp8 || Casp9)
		{
			return true; // Cell should remain quiescent due to apoptosis
		}
		
		// Quiescence logic based on cyclin levels
		// Cell is quiescent if CyclinD1 is low and no other cyclins are active
		if (!CyclinD1 && !CyclinA && !CyclinE && !CyclinB)
		{
			return true; // Cell should remain quiescent
		}
	}
	
	return false; // Cell can proceed with normal cell cycle
}

// Helper function to check if cell should enter apoptosis
bool check_boolean_network_apoptosis( Cell* pCell )
{
	// Only proceed if the cell has a MaBoSS model
	if (pCell->phenotype.intracellular &&
		pCell->phenotype.intracellular->intracellular_type == "maboss")
	{
		MaBoSSIntracellular* maboss_model = static_cast<MaBoSSIntracellular*>(pCell->phenotype.intracellular);

		// Check apoptosis signals
		bool Casp3 = maboss_model->maboss.get_node_value("Casp3");
		bool Casp8 = maboss_model->maboss.get_node_value("Casp8");
		bool Casp9 = maboss_model->maboss.get_node_value("Casp9");

		// Get apoptosis timer variables
		double& apoptosis_commitment = pCell->custom_data["apoptosis_commitment"];
		double apoptosis_threshold = pCell->custom_data["apoptosis_threshold"];

		// If any caspase is active, increase apoptosis commitment
		if (Casp3 || Casp8 || Casp9)
		{
			apoptosis_commitment += 1.0; // Increase commitment by 1 each time step
			std::cout << "Cell " << pCell->ID << " apoptosis commitment: " << apoptosis_commitment << "/" << apoptosis_threshold << std::endl;
		}

		// Only trigger apoptosis if commitment exceeds threshold
		if (apoptosis_commitment >= apoptosis_threshold)
		{
			return true;
		}
	}

	return false;
}

// Helper function to get cyclin-based transition readiness
bool check_cyclin_transition_readiness( Cell* pCell, const std::string& target_phase )
{
	// Only proceed if the cell has a MaBoSS model
	if (pCell->phenotype.intracellular &&
		pCell->phenotype.intracellular->intracellular_type == "maboss")
	{
		MaBoSSIntracellular* maboss_model = static_cast<MaBoSSIntracellular*>(pCell->phenotype.intracellular);
		
		// Get cyclin levels
		bool CyclinD1 = maboss_model->maboss.get_node_value("CyclinD1");
		bool CyclinA = maboss_model->maboss.get_node_value("CyclinA");
		bool CyclinE = maboss_model->maboss.get_node_value("CyclinE");
		bool CyclinB = maboss_model->maboss.get_node_value("CyclinB");
		
		// Check transition readiness based on target phase
		if (target_phase == "S")
		{
			// G0/G1 → S: Need CyclinD1 or CyclinE or CyclinA
			return (CyclinD1 || CyclinE || CyclinA);
		}
		else if (target_phase == "G2")
		{
			// S → G2: Need CyclinA (S phase completion)
			return CyclinA;
		}
		else if (target_phase == "M")
		{
			// G2 → M: Need CyclinB
			return CyclinB;
		}
		else if (target_phase == "G0G1")
		{
			// M → G0/G1: CyclinB should be low (mitosis completion)
			return !CyclinB;
		}
	}
	
	return false;
}

// G0/G1 phase entry function - handles quiescence based on Boolean network
void custom_G0G1_phase_entry_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	// Check for apoptosis first - if apoptosis is active, cell should not cycle
	if (check_boolean_network_apoptosis(pCell))
	{
		std::cout << "Cell " << pCell->ID << " entering G0/G1 phase but apoptosis is active - cell will remain quiescent" << std::endl;
		
		// Set very low transition rates to keep cell in G0/G1
		Cycle_Model* cycle_model = &(phenotype.cycle.model());
		int G0G1_index = cycle_model->find_phase_index(PhysiCell_constants::G0G1_phase);
		int S_index = cycle_model->find_phase_index(PhysiCell_constants::S_phase);
		cycle_model->transition_rate(G0G1_index, S_index) = 0.00001; // Very low rate
		return;
	}
	
	// Check if the Boolean network indicates quiescence based on cyclin levels
	if (check_boolean_network_quiescence(pCell))
	{
		std::cout << "Cell " << pCell->ID << " entering quiescent G0/G1 phase - low cyclin levels" << std::endl;
		
		// Reduce transition rate from G0/G1 to S phase when quiescent
		Cycle_Model* cycle_model = &(phenotype.cycle.model());
		int G0G1_index = cycle_model->find_phase_index(PhysiCell_constants::G0G1_phase);
		int S_index = cycle_model->find_phase_index(PhysiCell_constants::S_phase);
		cycle_model->transition_rate(G0G1_index, S_index) = 0.0001; // Very low rate
	}
	else
	{
		// Check if cell is ready to transition to S phase based on cyclins
		if (check_cyclin_transition_readiness(pCell, "S"))
		{
			std::cout << "Cell " << pCell->ID << " entering active G0/G1 phase - cyclins ready for S phase" << std::endl;
			
			// Normal transition rate from G0/G1 to S phase
			Cycle_Model* cycle_model = &(phenotype.cycle.model());
			int G0G1_index = cycle_model->find_phase_index(PhysiCell_constants::G0G1_phase);
			int S_index = cycle_model->find_phase_index(PhysiCell_constants::S_phase);
			cycle_model->transition_rate(G0G1_index, S_index) = 0.00324; // Normal rate
		}
		else
		{
			std::cout << "Cell " << pCell->ID << " entering G0/G1 phase - waiting for cyclin activation" << std::endl;
			
			// Reduced transition rate - waiting for cyclins
			Cycle_Model* cycle_model = &(phenotype.cycle.model());
			int G0G1_index = cycle_model->find_phase_index(PhysiCell_constants::G0G1_phase);
			int S_index = cycle_model->find_phase_index(PhysiCell_constants::S_phase);
			cycle_model->transition_rate(G0G1_index, S_index) = 0.0005; // Reduced rate
		}
	}
}

// S phase entry function - DNA synthesis phase
void custom_S_phase_entry_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	// Check for apoptosis first
	if (check_boolean_network_apoptosis(pCell))
	{
		std::cout << "Cell " << pCell->ID << " entering S phase but apoptosis is active - cell will remain in S phase" << std::endl;
		
		// Set very low transition rate to G2
		Cycle_Model* cycle_model = &(phenotype.cycle.model());
		int S_index = cycle_model->find_phase_index(PhysiCell_constants::S_phase);
		int G2_index = cycle_model->find_phase_index(PhysiCell_constants::G2_phase);
		cycle_model->transition_rate(S_index, G2_index) = 0.00001; // Very low rate
		return;
	}
	
	std::cout << "Cell " << pCell->ID << " entering S phase (DNA synthesis)" << std::endl;
	
	// Check cyclin levels for S phase progression
	if (pCell->phenotype.intracellular &&
		pCell->phenotype.intracellular->intracellular_type == "maboss")
	{
		MaBoSSIntracellular* maboss_model = static_cast<MaBoSSIntracellular*>(pCell->phenotype.intracellular);
		bool CyclinA = maboss_model->maboss.get_node_value("CyclinA");
		bool CyclinE = maboss_model->maboss.get_node_value("CyclinE");
		
		if (CyclinA)
		{
			std::cout << "Cell " << pCell->ID << " CyclinA is active during S phase - ready for G2 transition" << std::endl;
			
			// Normal transition rate to G2
			Cycle_Model* cycle_model = &(phenotype.cycle.model());
			int S_index = cycle_model->find_phase_index(PhysiCell_constants::S_phase);
			int G2_index = cycle_model->find_phase_index(PhysiCell_constants::G2_phase);
			cycle_model->transition_rate(S_index, G2_index) = 0.00208; // Normal rate
		}
		else if (CyclinE)
		{
			std::cout << "Cell " << pCell->ID << " CyclinE is active during S phase - early S phase" << std::endl;
			
			// Reduced transition rate - still in early S phase
			Cycle_Model* cycle_model = &(phenotype.cycle.model());
			int S_index = cycle_model->find_phase_index(PhysiCell_constants::S_phase);
			int G2_index = cycle_model->find_phase_index(PhysiCell_constants::G2_phase);
			cycle_model->transition_rate(S_index, G2_index) = 0.001; // Reduced rate
		}
		else
		{
			std::cout << "Cell " << pCell->ID << " No cyclins active during S phase - waiting for CyclinA" << std::endl;
			
			// Very low transition rate - waiting for CyclinA
			Cycle_Model* cycle_model = &(phenotype.cycle.model());
			int S_index = cycle_model->find_phase_index(PhysiCell_constants::S_phase);
			int G2_index = cycle_model->find_phase_index(PhysiCell_constants::G2_phase);
			cycle_model->transition_rate(S_index, G2_index) = 0.0002; // Very low rate
		}
	}
	
	// Standard S phase behavior - double nuclear volume
	// This is already handled by the standard S_phase_entry_function in PhysiCell
}

// G2 phase entry function - Gap 2 phase
void custom_G2_phase_entry_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	// Check for apoptosis first
	if (check_boolean_network_apoptosis(pCell))
	{
		std::cout << "Cell " << pCell->ID << " entering G2 phase but apoptosis is active - cell will remain in G2 phase" << std::endl;
		
		// Set very low transition rate to M
		Cycle_Model* cycle_model = &(phenotype.cycle.model());
		int G2_index = cycle_model->find_phase_index(PhysiCell_constants::G2_phase);
		int M_index = cycle_model->find_phase_index(PhysiCell_constants::M_phase);
		cycle_model->transition_rate(G2_index, M_index) = 0.00001; // Very low rate
		return;
	}
	
	std::cout << "Cell " << pCell->ID << " entering G2 phase" << std::endl;
	
	// Check cyclin levels for G2 phase progression
	if (pCell->phenotype.intracellular &&
		pCell->phenotype.intracellular->intracellular_type == "maboss")
	{
		MaBoSSIntracellular* maboss_model = static_cast<MaBoSSIntracellular*>(pCell->phenotype.intracellular);
		bool CyclinB = maboss_model->maboss.get_node_value("CyclinB");
		bool CyclinA = maboss_model->maboss.get_node_value("CyclinA");
		
		if (CyclinB)
		{
			std::cout << "Cell " << pCell->ID << " CyclinB is active during G2 phase - ready for mitosis" << std::endl;
			
			// Normal transition rate to M phase
			Cycle_Model* cycle_model = &(phenotype.cycle.model());
			int G2_index = cycle_model->find_phase_index(PhysiCell_constants::G2_phase);
			int M_index = cycle_model->find_phase_index(PhysiCell_constants::M_phase);
			cycle_model->transition_rate(G2_index, M_index) = 0.00417; // Normal rate
		}
		else if (CyclinA)
		{
			std::cout << "Cell " << pCell->ID << " CyclinA is active during G2 phase - early G2 phase" << std::endl;
			
			// Reduced transition rate - still in early G2 phase
			Cycle_Model* cycle_model = &(phenotype.cycle.model());
			int G2_index = cycle_model->find_phase_index(PhysiCell_constants::G2_phase);
			int M_index = cycle_model->find_phase_index(PhysiCell_constants::M_phase);
			cycle_model->transition_rate(G2_index, M_index) = 0.002; // Reduced rate
		}
		else
		{
			std::cout << "Cell " << pCell->ID << " No cyclins active during G2 phase - waiting for CyclinB" << std::endl;
			
			// Very low transition rate - waiting for CyclinB
			Cycle_Model* cycle_model = &(phenotype.cycle.model());
			int G2_index = cycle_model->find_phase_index(PhysiCell_constants::G2_phase);
			int M_index = cycle_model->find_phase_index(PhysiCell_constants::M_phase);
			cycle_model->transition_rate(G2_index, M_index) = 0.0003; // Very low rate
		}
	}
}

// M phase entry function - Mitosis phase
void custom_M_phase_entry_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	// Check for apoptosis first
	if (check_boolean_network_apoptosis(pCell))
	{
		std::cout << "Cell " << pCell->ID << " entering M phase but apoptosis is active - cell will remain in M phase" << std::endl;
		
		// Set very low transition rate back to G0/G1
		Cycle_Model* cycle_model = &(phenotype.cycle.model());
		int M_index = cycle_model->find_phase_index(PhysiCell_constants::M_phase);
		int G0G1_index = cycle_model->find_phase_index(PhysiCell_constants::G0G1_phase);
		cycle_model->transition_rate(M_index, G0G1_index) = 0.00001; // Very low rate
		return;
	}
	
	std::cout << "Cell " << pCell->ID << " entering M phase (mitosis)" << std::endl;
	
	// Check cyclin levels for mitosis progression
	if (pCell->phenotype.intracellular &&
		pCell->phenotype.intracellular->intracellular_type == "maboss")
	{
		MaBoSSIntracellular* maboss_model = static_cast<MaBoSSIntracellular*>(pCell->phenotype.intracellular);
		bool CyclinB = maboss_model->maboss.get_node_value("CyclinB");
		
		if (CyclinB)
		{
			std::cout << "Cell " << pCell->ID << " CyclinB is active during M phase - mitosis in progress" << std::endl;
			
			// Normal transition rate back to G0/G1 (mitosis completion)
			Cycle_Model* cycle_model = &(phenotype.cycle.model());
			int M_index = cycle_model->find_phase_index(PhysiCell_constants::M_phase);
			int G0G1_index = cycle_model->find_phase_index(PhysiCell_constants::G0G1_phase);
			cycle_model->transition_rate(M_index, G0G1_index) = 0.01667; // Normal rate
		}
		else
		{
			std::cout << "Cell " << pCell->ID << " CyclinB is low during M phase - mitosis completion" << std::endl;
			
			// Higher transition rate - mitosis is completing
			Cycle_Model* cycle_model = &(phenotype.cycle.model());
			int M_index = cycle_model->find_phase_index(PhysiCell_constants::M_phase);
			int G0G1_index = cycle_model->find_phase_index(PhysiCell_constants::G0G1_phase);
			cycle_model->transition_rate(M_index, G0G1_index) = 0.05; // Higher rate
		}
	}
	
	// Division will occur at the exit of this phase
	// This is handled automatically by PhysiCell
}

void update_gf_boolean_nodes( Cell* pCell )
{
	// Only proceed if the cell has a MaBoSS model
	if (pCell->phenotype.intracellular &&
		pCell->phenotype.intracellular->intracellular_type == "maboss")
	{
		MaBoSSIntracellular* maboss_model = static_cast<MaBoSSIntracellular*>(pCell->phenotype.intracellular);
		
		// Get the internalized GF concentration
		// GF is substrate ID 2 (oxygen=0, drug=1, GF=2)
		double gf_concentration = pCell->phenotype.molecular.internalized_total_substrates[2];
		
		// Thresholds for GF activation
		double gf_threshold = 0.25;      // 0.25 mM to activate GF
		double gf_high_threshold = 0.75; // 0.75 mM to activate GF_High
		
		// Update GF node based on concentration
		if (gf_concentration >= gf_threshold)
		{
			if (!maboss_model->maboss.get_node_value("GF"))
			{
				maboss_model->maboss.set_node_value("GF", true);
				std::cout << "Cell " << pCell->ID << " activated GF node (concentration: " << gf_concentration << " mM)" << std::endl;
			}
		}
		else
		{
			if (maboss_model->maboss.get_node_value("GF"))
			{
				maboss_model->maboss.set_node_value("GF", false);
				std::cout << "Cell " << pCell->ID << " deactivated GF node (concentration: " << gf_concentration << " mM)" << std::endl;
			}
		}
		
		// Update GF_High node based on concentration
		if (gf_concentration >= gf_high_threshold)
		{
			if (!maboss_model->maboss.get_node_value("GF_High"))
			{
				maboss_model->maboss.set_node_value("GF_High", true);
				std::cout << "Cell " << pCell->ID << " activated GF_High node (concentration: " << gf_concentration << " mM)" << std::endl;
			}
		}
		else
		{
			if (maboss_model->maboss.get_node_value("GF_High"))
			{
				maboss_model->maboss.set_node_value("GF_High", false);
				std::cout << "Cell " << pCell->ID << " deactivated GF_High node (concentration: " << gf_concentration << " mM)" << std::endl;
			}
		}
	}
}