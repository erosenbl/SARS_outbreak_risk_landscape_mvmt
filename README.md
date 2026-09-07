# SARS_outbreak_risk_landscape_mvmt
R code for manuscript titled "Landscape structure and animal movement influence 
the fate of SARS-CoV-2 spillover in white-tailed deer"

Elias Rosenblatt<sup>1,*</sup>, Guillaume Bastille-Rousseau<sup>2</sup>, Michael
Egan<sup>2</sup>, James Forester<sup>3</sup>, Billy J. Gardner<sup>4</sup>, 
Tyler Garwood<sup>1</sup>, Daniel M. Grove<sup>4</sup>, 
Tadao Kishimoto<sup>2</sup>, Justin Kosiewska<sup>4</sup>, 
Sung-Joo Lee<sup>5</sup>, Kezia R. Manlove<sup>6</sup>, 
Cameron Mitchell<sup>4</sup>, Lisa I. Muller<sup>4</sup>, 
Laura D. Plimpton<sup>7</sup>, Meredith C.VanAcker<sup>7</sup>, 
Mark Q. Wilber<sup>4</sup>, W. David Walter<sup>8</sup>, 
Maria A. Diuk-Wasser<sup>5</sup>, Grete Wilson-Henjum<sup>9</sup>, 
Tiffany Wolf<sup>10</sup>, Jacob Wyrick<sup>4</sup>, Kim M. Pepin<sup>9</sup>, 
and Meggan E. Craft<sup>1</sup>.

**Affiliations**:
<sup>1</sup> Department of Ecology, Evolution, and Behavior, University of Minnesota, 1987 Upper Buford Circle, St. Paul, MN 55108, USA
<sup>2</sup> Center for Wildlife Sustainability Research, Southern Illinois University, 1125 Lincoln Drive, Carbondale, IL 62901, USA
<sup>3</sup> Department of Fisheries, Wildlife, and Conservation Biology, University of Minnesota, 2003 Upper Buford Circle, Suite 135, St. Paul, MN 55108, USA 
<sup>4</sup> School of Natural Resources, University of Tennessee Institute of Agriculture, 2431 Joe Johnson Dr, Knoxville, TN 37996, USA 
<sup>5</sup> Department of Ecology, Evolution, and Environmental Biology, Columbia University, 1200 Amsterdam Ave, New York, NY 10027, USA
<sup>6</sup> Department of Wildland Resources and Ecology Center, Utah State University, 5230 Old Main Hill, Logan, UT 84322, USA
<sup>7</sup> Department of Evolution, Ecology, and Organismal Biology, University of California - Riverside, 900 University Ave, Riverside, CA 92521, USA
<sup>8</sup> Department of Veterinary Microbiology and Pathology, College of Veterinary Medicine, Washington State University, P.O. Box 647040, Pullman, WA 99164-7040
<sup>9</sup> United States Department of Agriculture, Animal and Plant Health Inspection Service,  Wildlife Services, National Wildlife Research Center,4101 Laporte Ave, Fort Collins, CO 80521, USA
<sup>10</sup> Department of Veterinary Population Medicine, University of Minnesota, 1365 Gortner Ave, Saint Paul, MN 55108, USA

<sup>*</sup>Corresponding author: erosenbl@umn.edu

Mention of commercial products does not represent an endorsement by the US government. The findings and conclusions in this publication are those of the authors and should not be construed to represent any official USDA or US Government determination or policy.

## Scripts folder
### Function script
**function_definitions.R** - A script to load packages and define custom 
functions for simulations. If interested in running contents of Example 
simulation folder, *run this script first*.

### contact_rate_calculation_and_selection subfolder

**example_contact_rate_calculation.R** - A script to show how raw pair-wise 
contact rates were calculated using collar data.

**select_slow_and_fast_decay_contact_kernels.R** - A script to identify one slow
decay contact kernel and one fast decay contact kernel from the nine study areas
considered.

## example_simulation folder
### Study area creation script
**study_area_grid_creation.R** - A script to demonstration the delineation of a 
study area grid, for use in **example_simulation_run.R**. 
*This script is not functional*, but the resulting file 
(**Inputs/ECRP_grid.RDS**) has been included in this repo for use in **example_simulation_run.R**.

### Simulation input file creation script
**simulation_input_creation.R** - A script to create an input file necessary 
for use in **example_simulation_run.R**. *This script is not functional*, but 
the resulting file (**Inputs/ECRP_sim_input_fall.RDS**) has been included in 
this repo for use in 
**example_simulation_run.R**.

### Example simulation script
**example_simulation_run.R** - An example script running simulations for 
Elm Creek, MN with the density and contact kernel combinations used in 
this study. All inputs are linked to RDS files in the 
*example_simulation/Inputs* folder. This script is functional, but will take a 
long time to run (see caution note for details).

### Inputs subfolder 
This folder includes RDS files necessary to run **example_simulation_run.R**. 
Raw location and collar metadata are not included here, but may be available 
upon request.

**contact_rates_annual_CP.RDS** - glmTMB model describing how pair-wise contact 
rates decay with distance for the study area selected for the fast decay spatial
contact kernel condition in simulations. 

**contact_rates_annual_lone_oaks.RDS** - glmTMB model describing how pair-wise 
contact  rates decay with distance for the study area selected for the slow 
decay spatial contact kernel condition in simulations. 

**ECRP_grid_dist.RDS** - A data frame containing pairwise distances between grid 
cell centroids for the Elm Creek, MN study area 
(from the grid **ECRP_grid.RDS**).

**ECRP_grid.RDS** - spatVector (terra) object with example grid for the
Elm Creek, MN study area.

**ECRP_median_mvmt.RDS** - Example dataframe summarizing median daily 
displacement for female and male deer, for the Elm Creek, MN study area.

**ECRP_mvmt_shape_params.RDS** - Example dataframe containing gamma distribution
shape parameters for daily displacement, by sex, for Elm Creek, MN study area.

**ECRP_NLCD.RDS** - A SpatRaster with landcover categories from the 30m-pixel 
National Land Cover Database (NLCD) for for Elm Creek, MN study area. This 
raster file is used in the creation of several summary dataframes here, and used
in **Random forest analysis/example_landscape_metrics_calculation.R**.

**ECRP_sim_input_fall.RDS** - Example list of simulation inputs for the Elm 
Creek, MN study area, for use in **Example_simulation_run.R**. List includes: 
population size (pop_size; empty to be populated in simulation), 
average group size (ave_group_size), sex ratio (sex_ratio), 
collared individuals for mimic habitat preferences (collar_roster), 
selection probabilities for daily movement (selection_probs; a data frame 
that lists grid cells and the 3rd order selection probability for each 
collared animal), 
selection probabilities for activity centroid selection (origin_probs; 
same format as selection_probs), 
and sex-specific centroid exponetial decay models (centroid_decay).

## random_forest_analysis folder
This folder contains the necessary summary scripts, analytical scripts, and 
summary data to demonstrate the random forest analysis incuded in this study.

### Landscape metrics summary script
**example_landscape_metrics_calculation.R** - Example extraction of relevant 
landscape variables tested for influencing outbreak probability, epidemic growth
rate, and spatial spread rate. Script ends with function demonstration for 
summarizing landscape metrics for Elm Creek, MN.

### Random forest model fitting scripts
**random_forest_model_fits.R** - Script to take simulation output 
(**random_forest_ _analysis/Data/df_RF_final.RDS**) and fit random forest models
for outbreak probability, maximum incidence rate, and maximum spread rate, 
resulting in Figure 3 in the manuscript. Users could modify the output from 
these fitted models to recreate Figure 4 in the manusctipt.

###Data subfolder
**df_RF_final.RDS** - A data frame containing all simulations across all study 
areas, with relevant landscape and movement metrics used in 
**random_forest_model_fits.R**. Note that users can use these summary data to 
recreate Figure 2, Figure S3 and S4.

