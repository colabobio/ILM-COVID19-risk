/***
* Part of the GAMA CoVid19 Modeling Kit
* see http://gama-platform.org/covid19
* Author: Benoit Gaudou, Patrick Taillandier
* Tags: covid19,epidemiology
***/

@no_experiment

model CoVid19


global {
	//app based model
	float subset_probability <- 1.0; // for taking subset of all the contacts being saved for parameter estimation
	bool quaranatine_enabled <- true;
	int test_delay <- 0;
	int quarantine_time <- 14;

    // Par set 1
//	float A0<-0.2;
//	float A1<-2.0;
//	float B0<-0.2;
//	float B1<-2.0;
//	float eA0<-0.3;
//	float eA1<-2.3;
//	float eB0<-0.18;
//	float eB1<-2.0;
//	float quarantine_risk<-0.18;

    // Par set 2
//	float A0<-0.5; 
//	float A1<-1.5;
//	float B0<-0.5;
//	float B1<-0.8;
//	float eA0<-0.4; 
//	float eA1<-1.3;
//	float eB0<-0.65;
//	float eB1<-0.43;
//	float quarantine_risk<-0.12;

    // Par set 3
	float A0<-1.0; 
	float A1<-2.0;
	float B0<-1.0;
	float B1<-1.0;
	float eA0<-1.3; 
	float eA1<-3.1;
	float eB0<-0.7;
	float eB1<-1.1;
	float quarantine_risk<-0.32;	

	//Epidemiological status of the individual
	string susceptible <- "S";
	string latent <- "L";
	string asymptomatic <- "Ia";
	string presymptomatic <- "Ip";
	string symptomatic <- "Is";
	string recovered <- "R";
	string dead <- "D";
	
	//Diagnostic status of the individual
	string not_tested <- "Not tested";
	string tested_positive <- "Positive";
	string tested_negative <- "Negative";
	
	//Hospitalized status of the individual
	string no_need_hospitalisation <- "Not needed";
	string need_hospitalisation <- "Need hospitalisation";
	string need_ICU <- "Need ICU";

	
	string act_neighbor <- "visiting neighbor";
	string act_friend <- "visiting friend";
	string act_home <- "staying at home";
	string act_working <- "working";
	string act_studying <- "studying";
	string act_eating <- "eating";
	string act_shopping <- "shopping";
	string act_leisure <- "leisure";
	string act_outside <- "outside activity";
	string act_sport <- "sport";
	string act_other <- "other activity";
	
	//Type of model for building choice during activity
	string random <- "random";
	string gravity <- "gravity";
	string closest <- "closest";
	
	// OSM Constant (type of building)
	list<string> OSM_eat <- ["restaurant","bakery"];
	list<string> OSM_home <- ["yes","house", "manor","apartments",'chocolate','shoes',"caravan"];
	list<string> OSM_shop <- ['commercial','supermarket',"bakery","frozen_food","alcohol","retail","furniture","bicycle"];
	list<string> OSM_outside_activity <- [];
	list<string> OSM_leisure <- [];
	list<string> OSM_sport <- ['tennis','multi','basketball','soccer','rugby_league','swimming','cycling','pelota','boules','skateboard','beachvolleyball','athletics'];
	list<string> OSM_other_activity <- ['car_repair','garages','church','hairdresser',"chapel","memorial","ruins"];
	list<string> OSM_work_place <- ['office',"estate_agent","public","civic","government","manufacture","company"];
	list<string> OSM_school <- ["school"];
	
	//number of the column for the epidemiological parameters CSV file
	int epidemiological_csv_column_name <- 0; //Name of the parameter
	int epidemiological_csv_column_age <- 1; //Lower bound of the age category
	int epidemiological_csv_column_detail <- 2; //Detail of the parameter (i.e. Fixed, or following a distribution)
	int epidemiological_csv_column_parameter_one <- 3; //Value of the parameter (only this one is used for fixed, else it is the first parameter of the distribution)
	int epidemiological_csv_column_parameter_two <- 4; //Value of the parameter (only used as the second parameter for distribution)
	
	//Keys of the map of epidemiological parameters, must also be used in the CSV
	string epidemiological_transmission_human <- "Transmission_human";
	string epidemiological_transmission_building <- "Transmission_building";
	string epidemiological_basic_viral_decrease <- "Basic_viral_decrease";
	string epidemiological_fixed <- "Fixed";
	string epidemiological_lognormal <- "Lognormal";
	string epidemiological_normal <- "Normal";
	string epidemiological_weibull <- "Weibull";
	string epidemiological_gamma <- "Gamma";
	string epidemiological_uniform <- "Uniform";
	string epidemiological_successful_contact_rate_human <- "Successful_contact_rate_human";
	string epidemiological_successful_contact_rate_building <- "Successful_contact_rate_building";
	string epidemiological_factor_asymptomatic <-"Factor_asymptomatic";
	string epidemiological_proportion_asymptomatic <- "Proportion_asymptomatic";
	string epidemiological_basic_viral_release <- "Basic_viral_release";
	string epidemiological_probability_true_positive <- "Probability_true_positive";
	string epidemiological_probability_true_negative <- "Probability_true_negative";
	string epidemiological_proportion_wearing_mask <- "Proportion_wearing_mask";
	string epidemiological_factor_wearing_mask <- "Factor_wearing_mask";
	string epidemiological_incubation_period <-"Incubation_period";
	string epidemiological_serial_interval <- "Serial_interval";
	string epidemiological_proportion_hospitalisation <- "Proportion_hospitalisation";
	string epidemiological_onset_to_hospitalisation <- "Onset_to_hospitalisation";
	string epidemiological_proportion_icu <- "Proportion_icu";
	string epidemiological_hospitalisation_to_ICU <- "Hospitalisation_to_ICU";
	string epidemiological_stay_ICU <- "Stay_ICU";
	string epidemiological_proportion_death_symptomatic <- "Proportion_death_symptomatic";
	string epidemiological_onset_to_recovery <- "Onset_to_recovery";
	
}