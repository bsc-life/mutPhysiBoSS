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
#include <algorithm>
#include <cctype>
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
	cell_defaults.functions.cell_division_function = custom_cell_division_function; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	cell_defaults.custom_data.add_variable(parameters.strings("node_to_visualize"), "dimensionless", 0.0 ); //for paraview visualization
	
	// Add PLK1 node state tracking for debugging and analysis
	cell_defaults.custom_data.add_variable("Plk1", "dimensionless", 0.0);
	
	// Add entry node state tracking for debugging
	cell_defaults.custom_data.add_variable("S_entry", "dimensionless", 0.0);
	cell_defaults.custom_data.add_variable("G2M_entry", "dimensionless", 0.0);
	cell_defaults.custom_data.add_variable("G0G1_entry", "dimensionless", 0.0);
	
	// Add apoptosis timer variables
	cell_defaults.custom_data.add_variable("apoptosis_commitment", "dimensionless", 0.0);
	cell_defaults.custom_data.add_variable("apoptosis_threshold", "dimensionless", 100.0);

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

                // Assign arrest functions to prevent transitions when conditions aren't met
                if (G0G1_index >= 0 && S_index >= 0)
                    cycle_model->phase_link(G0G1_index, S_index).arrest_function = arrest_G0G1_to_S;
                if (S_index >= 0 && G2_index >= 0)
                    cycle_model->phase_link(S_index, G2_index).arrest_function = arrest_S_to_G2;
                if (G2_index >= 0 && M_index >= 0)
                    cycle_model->phase_link(G2_index, M_index).arrest_function = arrest_G2_to_M;
                if (M_index >= 0 && G0G1_index >= 0)
                {
                    cycle_model->phase_link(M_index, G0G1_index).arrest_function = arrest_M_to_G0G1;
                    cycle_model->phase_link(M_index, G0G1_index).exit_function = phase_exit_mutation_function;
                }
                
                std::cout << "Configured flow-cytometry separated cycle customizations for " << cell_type << " cells." << std::endl;
            }
            else if (cycle_model->code == PhysiCell_constants::flow_cytometry_cycle_model)
            {
                int G2M_index = cycle_model->find_phase_index(PhysiCell_constants::G2M_phase);
                int G0G1_index = cycle_model->find_phase_index(PhysiCell_constants::G0G1_phase);
                if (G2M_index >= 0 && G0G1_index >= 0)
                {
                    cycle_model->phase_link(G2M_index, G0G1_index).exit_function = phase_exit_mutation_function;
                    std::cout << "Applied mutation exit function to " << cell_type << " cells (Flow Cytometry basic)." << std::endl;
                }
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
	
	// Track PLK1 node state for debugging and analysis
	if (pCell->phenotype.intracellular &&
		pCell->phenotype.intracellular->intracellular_type == "maboss")
	{
		MaBoSSIntracellular* maboss_model = static_cast<MaBoSSIntracellular*>(pCell->phenotype.intracellular);
		if (maboss_model->maboss.has_node("Plk1"))
		{
			bool plk1_state = maboss_model->maboss.get_node_value("Plk1");
			pCell->custom_data["Plk1"] = plk1_state ? 1.0 : 0.0;
		}
	}

	// std::cout << "Cell " << pCell->ID << " generation " << pCell->generation << " parent " << pCell->parent_ID << std::endl;

	// Update apoptosis rate based on caspase activity and commitment
	if (pCell->phenotype.intracellular &&
		pCell->phenotype.intracellular->intracellular_type == "maboss")
	{
		MaBoSSIntracellular* maboss_model = static_cast<MaBoSSIntracellular*>(pCell->phenotype.intracellular);
		bool Casp3 = maboss_model->maboss.get_node_value("Casp3");
		bool Casp8 = maboss_model->maboss.get_node_value("Casp8");
		bool Casp9 = maboss_model->maboss.get_node_value("Casp9");
		
		// Get apoptosis rate index
		int apoptosis_index = phenotype.death.find_death_model_index(PhysiCell_constants::apoptosis_death_model);
		
		// Get apoptosis commitment
		double& apoptosis_commitment = pCell->custom_data["apoptosis_commitment"];
		double apoptosis_threshold = pCell->custom_data["apoptosis_threshold"];
		
		// If any caspase is active, increase apoptosis commitment
		if (Casp3 || Casp8 || Casp9)
		{
			// Increase commitment at a rate of 0.1 per minute (scaled by dt)
			double commitment_rate = 0.1; // per minute
			apoptosis_commitment += commitment_rate * dt;
			
			// Calculate apoptosis rate based on commitment level
			// Base rate when caspases are active
			double base_apoptosis_rate = 1e-6; // Very low base rate
			
			// Scale apoptosis rate with commitment level
			// The longer caspases are active, the higher the rate
			// Rate increases linearly with commitment up to a maximum
			double max_apoptosis_rate = 0.01; // Maximum rate (1/min)
			double commitment_factor = apoptosis_commitment / apoptosis_threshold; // 0 to 1+
			commitment_factor = std::min(commitment_factor, 1.0); // Cap at 1.0
			
			// Linear scaling: rate = base + (max - base) * commitment_factor
			double apoptosis_rate = base_apoptosis_rate + (max_apoptosis_rate - base_apoptosis_rate) * commitment_factor;
			
			// Set the apoptosis rate
			phenotype.death.rates[apoptosis_index] = apoptosis_rate;
			
			// Only print occasionally to avoid spam
			static int print_counter = 0;
			if (print_counter % 100 == 0)
			{
				std::cout << "Cell " << pCell->ID << " apoptosis commitment: " << apoptosis_commitment 
						  << "/" << apoptosis_threshold << ", apoptosis rate: " << apoptosis_rate << " 1/min" << std::endl;
			}
			print_counter++;
		}
		else
		{
			// No caspases active - decay commitment and reduce apoptosis rate
			pCell->custom_data["apoptosis_commitment"] *= 0.95; // 5% decay per time step
			
			// Set very low apoptosis rate when caspases are not active
			double base_apoptosis_rate = 1e-6;
			phenotype.death.rates[apoptosis_index] = base_apoptosis_rate;
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
	
	// Track entry node states for debugging
	if (pCell->phenotype.intracellular &&
		pCell->phenotype.intracellular->intracellular_type == "maboss")
	{
		MaBoSSIntracellular* maboss_model = static_cast<MaBoSSIntracellular*>(pCell->phenotype.intracellular);
		if (maboss_model->maboss.has_node("S_entry"))
			pCell->custom_data["S_entry"] = maboss_model->maboss.get_node_value("S_entry") ? 1.0 : 0.0;
		if (maboss_model->maboss.has_node("G2M_entry"))
			pCell->custom_data["G2M_entry"] = maboss_model->maboss.get_node_value("G2M_entry") ? 1.0 : 0.0;
		if (maboss_model->maboss.has_node("G0G1_entry"))
			pCell->custom_data["G0G1_entry"] = maboss_model->maboss.get_node_value("G0G1_entry") ? 1.0 : 0.0;
	}
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
	static bool plk1_node_missing_reported = false;
	static bool plk1_target_warning_reported = false;
	static bool plk1_effect_warning_reported = false;
    
	
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
			
			double plk1_rate_mut_prob = parameters.doubles("plk1_rate_mutation_probability");
			if (plk1_rate_mut_prob > 0.0 && uniform_random() < plk1_rate_mut_prob)
			{
				const std::string plk1_node_name = "Plk1";
				if (!maboss_model->maboss.has_node(plk1_node_name))
				{
					if (!plk1_node_missing_reported)
					{
						std::cout << "Warning: Plk1 node not found in MaBoSS model; skipping PLK1 rate mutation." << std::endl;
						plk1_node_missing_reported = true;
					}
				}
				else
				{
					std::string target = parameters.strings("plk1_rate_mutation_target");
					std::string target_lower = target;
					std::transform(target_lower.begin(), target_lower.end(), target_lower.begin(),
						[](unsigned char c){ return static_cast<char>(std::tolower(c)); });
					
					const char* rate_symbol_cstr = nullptr;
					if (target_lower == "up" || target_lower == "activation" || target_lower == "rate_up")
					{
						rate_symbol_cstr = "$u_Plk1";
					}
					else if (target_lower == "down" || target_lower == "inhibition" || target_lower == "rate_down")
					{
						rate_symbol_cstr = "$d_Plk1";
					}
					else
					{
						if (!plk1_target_warning_reported)
						{
							std::cout << "Warning: Unsupported plk1_rate_mutation_target '" << target
									  << "'. Expected 'up' or 'down'. Defaulting to 'down'." << std::endl;
							plk1_target_warning_reported = true;
						}
						rate_symbol_cstr = "$d_Plk1";
					}
					
					double effect_size = parameters.doubles("plk1_rate_effect_size");
					if (effect_size < 0.0)
					{
						if (!plk1_effect_warning_reported)
						{
							std::cout << "Warning: plk1_rate_effect_size is negative; clamping to 0." << std::endl;
							plk1_effect_warning_reported = true;
						}
						effect_size = 0.0;
					}
					
					std::string rate_symbol(rate_symbol_cstr);
					double current_rate = maboss_model->maboss.get_parameter_value(rate_symbol);
					// Use ADDITIVE effect instead of multiplicative to avoid exponential explosion
					// new_rate = current_rate + effect_size (additive)
					// This provides controlled, linear growth instead of exponential
					double new_rate = current_rate + effect_size;
					maboss_model->maboss.set_parameter_value(rate_symbol, new_rate);
					
					// Debug output to verify rate changes
					std::cout << "Cell " << pCell->ID << " (gen " << pCell->generation << "): "
							  << "PLK1 rate mutation - " << rate_symbol 
							  << " changed from " << current_rate 
							  << " to " << new_rate 
							  << " (additive effect=" << effect_size << ")" << std::endl;
					
					std::string mutation_record = "Plk1_rate_" + std::string(rate_symbol == "$u_Plk1" ? "up" : "down")
						+ "_" + std::to_string(effect_size);
					pCell->custom_data.mutations.push_back(mutation_record);
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
bool check_boolean_network_apoptosis( Cell* pCell, double dt )
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
		// Use a rate-based approach: commitment increases at a rate per minute
		// Typical apoptosis commitment takes 30-60 minutes, so use a rate of ~0.1-0.2 per minute
		if (Casp3 || Casp8 || Casp9)
		{
			// Increase commitment at a rate of 0.1 per minute (scaled by dt)
			// This means it takes ~100 minutes to reach threshold of 10.0
			// Adjust the rate (0.1) to make it faster or slower as needed
			double commitment_rate = 0.1; // per minute
			apoptosis_commitment += commitment_rate * dt;
			
			// Only print occasionally to avoid spam
			static int print_counter = 0;
			if (print_counter % 100 == 0)
			{
				std::cout << "Cell " << pCell->ID << " apoptosis commitment: " << apoptosis_commitment 
						  << "/" << apoptosis_threshold << " (rate: " << commitment_rate << "/min)" << std::endl;
			}
			print_counter++;
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

// Arrest function: Prevent G0/G1 -> S transition if conditions aren't met
bool arrest_G0G1_to_S( Cell* pCell, Phenotype& phenotype, double dt )
{
	// Check quiescence
	if (check_boolean_network_quiescence(pCell))
	{
		static int debug_counter = 0;
		if (debug_counter % 100 == 0)
		{
			std::cout << "Cell " << pCell->ID << " arrested G0G1->S: quiescent" << std::endl;
		}
		debug_counter++;
		return true; // Arrest transition
	}
	
	// Get Boolean network entry nodes and cyclin levels
	bool S_entry_node = false;
	bool cyclins_ready = false;
	
	if (pCell->phenotype.intracellular &&
		pCell->phenotype.intracellular->intracellular_type == "maboss")
	{
		MaBoSSIntracellular* maboss_model = static_cast<MaBoSSIntracellular*>(pCell->phenotype.intracellular);
		S_entry_node = maboss_model->maboss.get_node_value("S_entry");
		cyclins_ready = check_cyclin_transition_readiness(pCell, "S");
		
		// Debug output occasionally
		static int debug_counter2 = 0;
		if (debug_counter2 % 100 == 0 && !S_entry_node)
		{
			std::cout << "Cell " << pCell->ID << " arrested G0G1->S: S_entry=" << (S_entry_node ? "ON" : "OFF") 
					  << ", cyclins_ready=" << (cyclins_ready ? "YES" : "NO") << std::endl;
		}
		debug_counter2++;
	}
	
	// Only allow transition if BOTH S_entry node is ON AND cyclins are ready
	if (S_entry_node && cyclins_ready)
	{
		return false; // Don't arrest - allow transition
	}
	
	return true; // Arrest transition - conditions not met
}

// G0/G1 phase entry function - simplified to just log entry
void custom_G0G1_phase_entry_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	std::cout << "Cell " << pCell->ID << " entering G0/G1 phase" << std::endl;
}

// Arrest function: Prevent S -> G2 transition if conditions aren't met
bool arrest_S_to_G2( Cell* pCell, Phenotype& phenotype, double dt )
{
	// Get Boolean network entry nodes and cyclin levels
	bool G2M_entry_node = false;
	bool CyclinA = false;
	
	if (pCell->phenotype.intracellular &&
		pCell->phenotype.intracellular->intracellular_type == "maboss")
	{
		MaBoSSIntracellular* maboss_model = static_cast<MaBoSSIntracellular*>(pCell->phenotype.intracellular);
		G2M_entry_node = maboss_model->maboss.get_node_value("G2M_entry");
		CyclinA = maboss_model->maboss.get_node_value("CyclinA");
	}
	
	// Only allow transition if BOTH G2M_entry node is ON AND CyclinA is active
	if (G2M_entry_node && CyclinA)
	{
		return false; // Don't arrest - allow transition
	}
	
	return true; // Arrest transition - conditions not met
}

// S phase entry function - simplified to just log entry
void custom_S_phase_entry_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	std::cout << "Cell " << pCell->ID << " entering S phase (DNA synthesis)" << std::endl;
}

// Arrest function: Prevent G2 -> M transition if conditions aren't met
bool arrest_G2_to_M( Cell* pCell, Phenotype& phenotype, double dt )
{
	// Get Boolean network entry nodes and cyclin levels
	bool G2M_entry_node = false;
	bool CyclinB = false;
	bool cyclins_ready_M = false;
	
	if (pCell->phenotype.intracellular &&
		pCell->phenotype.intracellular->intracellular_type == "maboss")
	{
		MaBoSSIntracellular* maboss_model = static_cast<MaBoSSIntracellular*>(pCell->phenotype.intracellular);
		G2M_entry_node = maboss_model->maboss.get_node_value("G2M_entry");
		CyclinB = maboss_model->maboss.get_node_value("CyclinB");
		cyclins_ready_M = check_cyclin_transition_readiness(pCell, "M");
	}
	
	// Only allow transition if G2M_entry node is ON AND CyclinB is active AND cyclins are ready
	if (G2M_entry_node && CyclinB && cyclins_ready_M)
	{
		return false; // Don't arrest - allow transition
	}
	
	return true; // Arrest transition - conditions not met
}

// G2 phase entry function - simplified to just log entry
void custom_G2_phase_entry_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	std::cout << "Cell " << pCell->ID << " entering G2 phase" << std::endl;
}

// Arrest function: Prevent M -> G0/G1 transition if conditions aren't met
bool arrest_M_to_G0G1( Cell* pCell, Phenotype& phenotype, double dt )
{
	// Get Boolean network entry nodes and cyclin levels
	bool G0G1_entry_node = false;
	bool CyclinB = false;
	bool cyclins_ready_G0G1 = false;
	
	if (pCell->phenotype.intracellular &&
		pCell->phenotype.intracellular->intracellular_type == "maboss")
	{
		MaBoSSIntracellular* maboss_model = static_cast<MaBoSSIntracellular*>(pCell->phenotype.intracellular);
		G0G1_entry_node = maboss_model->maboss.get_node_value("G0G1_entry");
		CyclinB = maboss_model->maboss.get_node_value("CyclinB");
		cyclins_ready_G0G1 = check_cyclin_transition_readiness(pCell, "G0G1");
	}
	
	// Allow transition if CyclinB is low (mitosis completing) AND (G0G1_entry ON OR cyclins ready)
	if (!CyclinB && (G0G1_entry_node || cyclins_ready_G0G1))
	{
		return false; // Don't arrest - allow transition
	}
	
	return true; // Arrest transition - conditions not met
}

// M phase entry function - simplified to just log entry
void custom_M_phase_entry_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	std::cout << "Cell " << pCell->ID << " entering M phase (mitosis)" << std::endl;
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
		double gf_concentration = pCell->nearest_density_vector()[microenvironment.find_density_index("GF")];
		// std::cout << "Cell " << pCell->ID << " GF concentration: " << gf_concentration << " mM" << std::endl;
		
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

// Custom cell division function to reset Boolean model to initial state (t=0)
void custom_cell_division_function( Cell* parent, Cell* child )
{
	// Reset the child cell's Boolean model to initial state (t=0)
	// This ensures each daughter cell starts fresh, regardless of inheritance settings
	if (child->phenotype.intracellular &&
		child->phenotype.intracellular->intracellular_type == "maboss")
	{
		MaBoSSIntracellular* child_maboss = static_cast<MaBoSSIntracellular*>(child->phenotype.intracellular);
		
		// Reset the Boolean model to initial state
		// This calls restart_node_values() which resets all nodes to their initial values
		child_maboss->start();
		
		// After reset, restore GF and GF_High nodes based on current concentration
		// This is important because resetting the model may have turned them OFF
		update_gf_boolean_nodes(child);
		
		// Also reset apoptosis commitment for the new cell
		child->custom_data["apoptosis_commitment"] = 0.0;
		
		std::cout << "Cell " << child->ID << " (daughter of " << parent->ID << "): Boolean model reset to initial state (t=0), GF nodes restored" << std::endl;
	}
}