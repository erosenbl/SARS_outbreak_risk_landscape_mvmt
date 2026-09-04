# SARS_outbreak_risk_landscape_mvmt
R code for manuscript titled "Landscape structure and animal movement influence 
the fate of SARS-CoV-2 spillover in white-tailed deer"

**Elias Rosenblatt**

## Scripts folder
### Function script
**function_definitions.R** - A script to load packages and define custom 
functions for simulations. If interested in running contents of Example 
simulation folder, *run this script first*. DONE

### Contact rate calculation and selection subfolder
#### Contact rate calculation example
**example_contact_rate_calculation.R** - A script to show how raw pair-wise 
contact rates were calculated using collar data. DONE

#### Contact kernel selection script
**select_slow_and_fast_decay_contact_kernels.R** - A script to identify one slow
decay contact kernel and one fast decay contact kernel from the nine study areas
considered. DONE

## Example simulation folder
### Study area creation script
**study_area_grid_creation.R** - A script to demonstration the delineation of a 
study area grid, for use in **Example_simulation_run.R**. *This script is not functional*, but 
the resulting file (**Inputs/ECRP_grid.RDS**) has been included in this repo for use in 
**Example_simulation_run.R**.

### Simulation input file creation script
**simulation_input_creation.R** - A script to create an input file necessary 
for use in **Example_simulation_run.R**. *This script is not functional*, but 
the resulting file (**Inputs/ECRP_sim_input_fall.RDS**) has been included in this repo for use in 
**Example_simulation_run.R**.

### Example simulation script
**Example_simulation_run.R** - An example script running simulations for 
Elm Creek, MN with the density and contact kernel combinations used in 
this study. All inputs are linked to RDS files in the 
*Example_simulation/Inputs* folder. This script is functional, but will take a 
long time to run (see caution note for details).

### Inputs subfolder 
This folder includes RDS files necessary to run **Example_simulation_run.R**. 
Raw location and collar metadata are not included here, but may be available 
upon request.

**contact_rates_annual_CP.RDS** - glmTMB model describing how pair-wise contact 
rates decay with distance for the study area selected for the fast decay spatial
contact kernel condition in simulations. 

**contact_rates_annual_lone_oaks.RDS** - glmTMB model describing how pair-wise 
contact  rates decay with distance for the study area selected for the slow 
decay spatial contact kernel condition in simulations. 

**ECRP_grid.RDS** - spatVector (terra) object with example grid for the
Elm Creek, MN study area.

**ECRP_median_mvmt.RDS** - Example dataframe summarizing median daily 
displacement for female and male deer, for the Elm Creek, MN study area.

**ECRP_mvmt_shape_params.RDS** - Example dataframe containing gamma distribution
shape parameters for daily displacement, by sex, for Elm Creek, MN study area.

**ECRP_sim_input_fall.RDS** - Example list of simulation inputs for the Elm Creek, MN study area, for use in **Example_simulation_run.R**. List includes: population size (pop_size; empty to be populated in simulation), average group size (ave_group_size), sex ratio (sex_ratio), collared individuals for mimic habitat preferences (collar_roster), selection probabilities for daily movement (selection_probs; a data frame that lists grid cells and the 3rd order selection probability for each collared animal), selection probabilities for activity centroid selection (origin_probs; same format as selection_probs), and sex-specific centroid exponetial decay models (centroid_decay).

## Random forest analysis folder

