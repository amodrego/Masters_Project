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

#include "./custom.h"
#include "../core/PhysiCell_cell.h"
using namespace std;

// declare cell definitions here 

Cell_Definition motile_cell; 

void create_cell_types( void )
{
	// use the same random seed so that future experiments have the 
	// same initial histogram of oncoprotein, even if threading means 
	// that future division and other events are still not identical 
	// for all runs 
	
	SeedRandom( parameters.ints("random_seed") ); // or specify a seed here 
	
	// housekeeping 
	
	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	// Name the default cell type 
	
	cell_defaults.type = 0; 
	cell_defaults.name = "tumor cell"; 
	
	// set default cell cycle model 

	cell_defaults.functions.cycle_model = flow_cytometry_separated_cycle_model; 
	
	// set default_cell_functions; 
	
	cell_defaults.functions.update_phenotype = update_cell_and_death_parameters_O2_based; 
	
	// only needed for a 2-D simulation: 
	
	/*
	cell_defaults.functions.set_orientation = up_orientation; 
	cell_defaults.phenotype.geometry.polarity = 1.0;
	cell_defaults.phenotype.motility.restrict_to_2D = true; 
	*/
	
	// make sure the defaults are self-consistent. 
	
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	cell_defaults.phenotype.molecular.sync_to_microenvironment( &microenvironment );
	cell_defaults.phenotype.sync_to_functions( cell_defaults.functions );

	// set the rate terms in the default phenotype 

	// first find index for a few key variables. 
	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	int necrosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Necrosis" );
	int oxygen_substrate_index = microenvironment.find_density_index( "oxygen" ); 
	int glucose_substrate_index = microenvironment.find_density_index("glucose");
	int lactate_substrate_index = microenvironment.find_density_index("lactate");

	int G0G1_index = flow_cytometry_separated_cycle_model.find_phase_index( PhysiCell_constants::G0G1_phase );
	int S_index = flow_cytometry_separated_cycle_model.find_phase_index( PhysiCell_constants::S_phase );

	// initially no necrosis 
	cell_defaults.phenotype.death.rates[necrosis_model_index] = 0.0; 

	// set oxygen uptake / secretion parameters for the default cell type 
	cell_defaults.phenotype.secretion.uptake_rates[oxygen_substrate_index] = 10; // [1/min]
	cell_defaults.phenotype.secretion.secretion_rates[oxygen_substrate_index] = 0; // []
	cell_defaults.phenotype.secretion.saturation_densities[oxygen_substrate_index] = 38; // [mmHg]

	// set glucose uptake / secretion parameters for default cell type
	// This value is 10/6 due to the stechiometric relation glucose-oxygen to produce ATP
	cell_defaults.phenotype.secretion.uptake_rates[glucose_substrate_index] = 1.67; // [1/min]
	cell_defaults.phenotype.secretion.secretion_rates[glucose_substrate_index] = 0;
	cell_defaults.phenotype.secretion.saturation_densities[glucose_substrate_index] = 4.5; // [mg/ml]

	// Lactate values
	cell_defaults.phenotype.secretion.uptake_rates[lactate_substrate_index] = 0.0; // [1/min]
	cell_defaults.phenotype.secretion.secretion_rates[lactate_substrate_index] = 0.0;
	cell_defaults.phenotype.secretion.saturation_densities[lactate_substrate_index] = 20.0; // [mM]

	// cell velocity takes into account ECM drag forces
	cell_defaults.functions.update_velocity = drag_update_cell_velocity;
	
	// add custom data here, if any 	
	
	// Now, let's define another cell type. 
	// It's best to just copy the default and modify it. 
	
	// make this cell type randomly motile, less adhesive, greater survival, 
	// and less proliferative 
	
	motile_cell = cell_defaults; 
	motile_cell.type = 1; 
	motile_cell.name = "motile tumor cell"; 
	
	// make sure the new cell type has its own reference phenotyhpe

	motile_cell.parameters.pReference_live_phenotype = &( motile_cell.phenotype ); 
	
	// enable random motility 
	motile_cell.phenotype.motility.is_motile = true; 
	motile_cell.phenotype.motility.persistence_time = parameters.doubles( "motile_cell_persistence_time" ); // 15.0; // 15 minutes
	motile_cell.phenotype.motility.migration_speed = parameters.doubles( "motile_cell_migration_speed" ); // 0.25; // 0.25 micron/minute 
	motile_cell.phenotype.motility.migration_bias = 0.0;// completely random 
	
	// Set cell-cell adhesion to 5% of other cells 
	//motile_cell.phenotype.mechanics.cell_cell_adhesion_strength *= 
	//	parameters.doubles( "motile_cell_relative_adhesion" ); // 0.05;
	
	// Set apoptosis to zero 
	motile_cell.phenotype.death.rates[apoptosis_model_index] = 
		parameters.doubles( "motile_cell_apoptosis_rate" ); // 0.0; 
	
	// Set proliferation to 10% of other cells. 
	// Alter the transition rate from G0G1 state to S state
	//motile_cell.phenotype.cycle.data.transition_rate(G0G1_index,S_index) *= 
	//	parameters.doubles( "motile_cell_relative_cycle_entry_rate" ); // 0.1; 
	
	build_cell_definitions_maps(); 
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
/*	
	default_microenvironment_options.X_range = {-500, 500}; 
	default_microenvironment_options.Y_range = {-500, 500}; 
	default_microenvironment_options.Z_range = {-500, 500}; 
*/	
	// make sure to override and go back to 2D 
	if( default_microenvironment_options.simulate_2D == true )
	{
		std::cout << "Warning: overriding XML config option and setting to 3D!" << std::endl; 
		default_microenvironment_options.simulate_2D = false; 
	}	
	
	
/*
	// all this is in XML as of August 2019 (1.6.0)
	// no gradients need for this example 

	default_microenvironment_options.calculate_gradients = false; 
	
	// set Dirichlet conditions 

	default_microenvironment_options.outer_Dirichlet_conditions = true;
	
	// if there are more substrates, resize accordingly 
	std::vector<double> bc_vector( 1 , 38.0 ); // 5% o2
	default_microenvironment_options.Dirichlet_condition_vector = bc_vector;
	
	// set initial conditions 
	default_microenvironment_options.initial_condition_vector = { 38.0 }; 
*/
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
	// create some cells near the origin
	
	Cell* pC;
	//create 4 motile cells in order to have all of them proliferating and migrating normally
	//pC = create_cell(motile_cell); 
	//pC->assign_position( 300.0, 300.0, 300.0 );

	//pC = create_cell(motile_cell); 
	//pC->assign_position( -300.0, 300.0, 300.0 );
	//
	//pC = create_cell(motile_cell);
	//pC->assign_position( 300.0, -300.0, 300.0 );
	
	//pC = create_cell(motile_cell); 
	//pC->assign_position( -300.0, -300.0, 300.0 );

	pC = create_cell(motile_cell);
	pC->assign_position(0.0, 0.0, 0.0);

	pC = create_cell(motile_cell);
	pC->assign_position( 300.0, 300.0, 0.0);
	//
	pC = create_cell(motile_cell);
	pC->assign_position( -300.0, 300.0, 0.0);

	pC = create_cell(motile_cell);
	pC->assign_position( 300.0, -300.0, 0.0);

	pC = create_cell(motile_cell);
	pC->assign_position(-300.0, -300.0, 0.0);
	

	// loop to position a number of cells randomly
	/*int cont = 0;
	while (cont < 8) {
		pC = create_cell(motile_cell);
		pC->assign_position(rand() % 400 + -400, rand() % 400 + -400, rand() % 400 + -400);
		cont++;
	}*/

	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
	// start with flow cytometry coloring 
	
	std::vector<std::string> output = false_cell_coloring_cytometry(pCell); 
	
	// if the cell is motile and not dead, paint it black 
	
	if( pCell->phenotype.death.dead == false && 
		pCell->type == 1 )
	{
		 output[0] = "black"; 
		 output[2] = "black"; 	
	}
	
	return output; 
}

double locomotive_forces_generator( )
{
	// random number generator to define cell locomotive forces
	// based on anempirically obtained velocity distribution

	double random_value, force_value;

	random_value = UniformRandom();

  force_value =
		0.1157 * pow(random_value, 3) + 0.2422 * pow(random_value, 2) +
		0.0053 * random_value + 0.0048;

	force_value = 13.5*force_value;

	return force_value;
}

void drag_update_cell_velocity( Cell* pCell, Phenotype& phenotype, double dt ) {
    // sample ECM
    int ECM_density_index = microenvironment.find_density_index( "ECM" );
    double ECM_density = pCell->nearest_density_vector()[ECM_density_index];
    double dyn_viscosity;
    // get viscosity based on concentration
    if(ECM_density == 2.5) {
        dyn_viscosity = 7.96;
    }
    else if(ECM_density == 4.0){
        dyn_viscosity = 18.42;
    }
    else if(ECM_density == 6.0) {
        dyn_viscosity = 39.15;
    }
    // update the speed value
    pCell->phenotype.motility.migration_speed = 
        locomotive_forces_generator();
    // update velocity
    standard_update_cell_velocity(pCell, phenotype, dt);
    // include the 1/vu (1/ECM density) term to consider friction
    pCell->velocity /= dyn_viscosity;
    return;
}

void get_proliferation_necrosis_rates(Cell* pCell, Phenotype& phenotype, double dt, int start_phase_index, int end_phase_index) {

	int oxygen_substrate_index = pCell->get_microenvironment()->find_density_index("oxygen");
	int glucose_substrate_index = pCell->get_microenvironment()->find_density_index("glucose");
	int lactate_substrate_index = pCell->get_microenvironment()->find_density_index("lactate");

	//COMBUSTION// OXYGEN:GLUCOSE RATIO -> 1Gl:6O2

		// Regla de 3:
		//		20% O2 ----------- 0.28 mM
		//		5% O2  ----------- 0.07 mM (38 mmHg)

		// Regla de 3 para obtener % a partir de los mmHg que tenemos (38 mmHg = 5% O2)

	double O2_perc = (pCell->parameters.pO2 * 5) / 38;
	pCell->parameters.mMglucose = (pCell->parameters.pGlucose / 180) * 1000; // glucose molecular weight [mM] 

	pCell->parameters.mMo2 = (O2_perc * 0.28) / 20; // [mM]


	pCell->parameters.Vg_max = pCell->parameters.qg_max *
		(1 - (1 - pCell->parameters.qg_min / pCell->parameters.qg_max) * pCell->parameters.mMglucose / (pCell->parameters.mMglucose + pCell->parameters.k_g_o));
	pCell->parameters.Vo_max = pCell->parameters.qo_max *
		(1 - (1 - pCell->parameters.qo_min / pCell->parameters.qo_max) * pCell->parameters.mMo2 / (pCell->parameters.mMo2 + pCell->parameters.k_o_g));

	// SE OBTIENEN LOS RATIOS DE CONSUMO QUE DEPENDEN DE LA CONCENTRACION DE AMBAS SUSTANCIAS, SEGUN LA DISTRIBUCION MICHAELIS-MENTEN [5]
	pCell->parameters.qg = pCell->parameters.Vg_max * pCell->parameters.mMglucose / (pCell->parameters.mMglucose * pCell->parameters.k_g);
	pCell->parameters.qo2 = pCell->parameters.Vo_max * pCell->parameters.mMo2 / (pCell->parameters.mMo2 * pCell->parameters.k_o2);


	// Calculates production of ATP with the available amount of oxygen and glucose in a given moment
	pCell->parameters.P_atp = 2 * pCell->parameters.qg + 17.0 / 3.0 * pCell->parameters.qo2;

	// K PROLIFERATION
	pCell->parameters.k_prolif = Ki67_basic.transition_rate(1, 0) * (pCell->parameters.P_atp - pCell->parameters.P_atp_min) / 60.0;
	
	if (pCell->parameters.k_prolif < 0)
	{
		pCell->parameters.k_prolif = 0.0;
	}

	// K NECROSIS
	// For necrosis max value lysed cells ratio is used  
	pCell->parameters.k_necrosis = necrosis.transition_rate(0, 1) * (pCell->parameters.P_atp_min - pCell->parameters.P_atp) / 60.0;
	if (pCell->parameters.k_necrosis < 0)
	{
		pCell->parameters.k_necrosis = 0.0;
	}

	// if oxygen is under hypoxic threshold (15 mmHg) anaerobic cycle begins but cells still proliferate even lacking oxygen. 
	// If it decreases under necrosis threshold, it stops proliferating and only lives
		

	if ((pCell->nearest_density_vector())[oxygen_substrate_index] < pCell->parameters.o2_hypoxic_threshold)
	{
		//pCell->phenotype.secretion.uptake_rates[oxygen_substrate_index] = 0.0;
		pCell->phenotype.secretion.secretion_rates[lactate_substrate_index] = pCell->phenotype.secretion.uptake_rates[glucose_substrate_index] * 2;
		//cell_defaults.phenotype.secretion.secretion_rates[lactate_substrate_index] = phenotype.secretion.uptake_rates[glucose_substrate_index] * 2;
		
		(pCell->nearest_density_vector())[lactate_substrate_index] = (pCell->nearest_density_vector())[lactate_substrate_index] 
			+ pCell->phenotype.secretion.secretion_rates[lactate_substrate_index];


		pCell->parameters.k_necrosis = necrosis.transition_rate(0, 1) * pow((pCell->nearest_density_vector())[oxygen_substrate_index], 2) /
			(pow(pCell->parameters.L_max, 2) - pow((pCell->nearest_density_vector())[oxygen_substrate_index], 2));
		
		if ((pCell->nearest_density_vector())[lactate_substrate_index] > pCell->parameters.L_max)
		{
			pCell->parameters.k_necrosis = 9e99;
			//phenotype.death.trigger_death(necrosis_index);
			
		}

		if (pCell->parameters.pO2 < pCell->parameters.o2_necrosis_threshold)
		{
			pCell->parameters.k_prolif = 0.0;
		}

		
	}

	return;

}

