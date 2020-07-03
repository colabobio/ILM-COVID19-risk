/***
* Part of the GAMA CoVid19 Modeling Kit
* see http://gama-platform.org/covid19
* Author: Huynh Quang Nghi, Patrick Taillandier, Damien Philippon
* Tags: covid19,epidemiology
***/


@no_experiment

model CoVid19

import "../Functions.gaml"
import "Activity.gaml"
import "Building.gaml"
import "../Constants.gaml"

global
{
	int total_number_of_infected <- 0;
	int total_number_reported <- 0;
	int total_number_individual <- 0;
	
	map<string, int> building_infections;
	map<int,int> total_incidence_age;
}

species Individual schedules: shuffle(Individual where (each.status != dead)){
	// app based model. only considering 2 covariates as of now. Can be expanded to more easily.
	float susceptibility <- get_susceptibility();
	float estimated_susceptibility <- get_estimated_susceptibility();
	
	float infectivity;
	float estimated_infectivity;
	
	bool has_preexisting_condition;
	float X <- ((65<age) or has_preexisting_condition) ? 1.0 :0.0;
	float Y <- 0.0;	
	int plen <- 14 * 4;
	list<float> noninf_prob <- list_with(plen, 1.0);
	float inf_risk <- 0.0;
	list<float> inf_risk_list <- list_with(plen, 0.0);
	
	float quarantined <- 0.0;
	int quarantine_count <- 0;
	
	// Symptoms predictive of COVID-19:
	// https://www.nature.com/articles/s41591-020-0916-2
    float smell_and_taste_loss <- 0.0;
    float cough <- 0.0;
    float fatigue <- 0.0;
    float skipped_meals <- 0.0;
	
	//#############################################################
	//App based model functions
	//#############################################################
	//Returns the susceptibility of the Individual
	float get_susceptibility{
		float individual_susceptibility <- A0 + A1*X;
		return individual_susceptibility;
	}
	float get_estimated_susceptibility {
		float individual_susceptibility <- eA0 + eA1*X;
		return individual_susceptibility;
	}
	//Age of the individual
	int age;
	//Sex of the individual
	int sex; //0 M 1 F
	
	//#############################################################
	//Location related attributes
	//#############################################################
	//ID of the household of the individual
	string household_id;
	//Bool to consider if the individual is at home
	bool is_at_home <- true;
	//Home building of the individual
	Building home;
	//School building of the individual (if student)
	Building school;
	//Working place of the individual (if working)
	Building working_place;
	//Relatives (i.e. same household) of the individual
	list<Individual> relatives;
	//Friends (i.e. possibility of leisure activities together) of the individual
	list<Individual> friends;
	//Colleagues (i.e. same working place) of the individual
	list<Individual> colleagues;	
	//List of contacts of the individual in the last cycle
	list<Individual> contacts;	
	//Current building of the individual
	Building current_place;
	//Bool to consider if the individual is outside of the commune
	bool is_outside <- false;
	
	
	//#############################################################
	//Epidemiologial attributes
	//#############################################################
	//Biological status of the individual (susceptible, latent, presymptomatic, symptomatic, asymptomatic, recovered, dead)
	string status <- susceptible; //S,L,Ip,Is,Ia,R,D
	//Report status of the individual (not-tested, negative, positive)
	string report_status <- not_tested; //Not-tested, Negative, Positive
	//Bool to represent the infectious potential of an individual (i.e. presymptomatic, symptomatic and asymptomatic)
	bool is_infectious <- define_is_infectious();
	//Bool to represent the infection status of an individual (i.e. latent, presymptomatic, symptomatic and asymptomatic)
	bool is_infected <- define_is_infected();
	//Bool to represent the asymptomatic status of an individual (i.e. presymptomatic and asymptomatic)
	bool is_asymptomatic <- define_is_asymptomatic();
	//Time in the latent period (corresponds to the incubation period for non-presymptomatic individuals)
	float latent_time; 
	//Time being infectious
	float infectious_time;
	//Serial interval (used for computation of the latent period)
	float serial_interval;
	//Contact rate with human individuals
	float contact_rate_human;
	//Viral release of the individual in the environment
	float basic_viral_release;
	//Variable for time related events
	int tick <- 0;
	
	//#############################################################
	//Agenda and activities attributes
	//#############################################################
	list<map<int, pair<Activity,list<Individual>>>> agenda_week;
	list<Individual> activity_fellows;
	Activity last_activity;
	map<Activity, list<Building>> building_targets;
	
	//#############################################################
	//Intervention related attributes
	//#############################################################
	//Reduction in the transmission when being asymptomatic (not coughing, etc)
	float factor_contact_rate_asymptomatic;
	//Reduction in the transmission when wearing a mask (coughing prevented)
	float factor_contact_rate_wearing_mask;
	//Bool to consider not following interventions
	bool free_rider <- false;
	//Probability of wearing a mask per time step
	float proba_wearing_mask;
	//Bool to represent wearing a mask
	bool is_wearing_mask;
	
	//#############################################################
	//Hospitalisation related attributes
	//#############################################################
	//Hospitalisation status of the individual (healthy, needing_hospitalisation, needing_ICU)
	string hospitalisation_status <- no_need_hospitalisation;
	//Bool to represent the fact of being hospitalised (not ICU)
	bool is_hospitalised <- false;
	//Bool to represent the fact of being in ICU (not hospitalised)
	bool is_ICU <- false;
	//Time before needing hospitalisation
	float time_before_hospitalisation;
	//Time before needing ICU
	float time_before_ICU;
	//Time to stay in ICU
	float time_stay_ICU;	
	
	//#############################################################
	//Tool functions
	//#############################################################
	//Return the fact of the Individual being infectious (i.e. asymptomatic, presymptomatic, symptomatic)
	bool define_is_infectious {
		return [asymptomatic,presymptomatic, symptomatic] contains status;
	}
	
	//Return the fact of the Individual not being infectious yet but infected (i.e. latent)
	bool is_latent {
		return status = latent;
	}
	
	//Return the fact of the Individual being infected (i.e. latent, asymptomatic, presymptomatic or symptomatic)
	bool define_is_infected {
		return is_infectious or self.is_latent();
	}
	
	//Return the fact of the Individual not showing any symptoms (i.e. asymptomatic, presymptomatic)
	bool define_is_asymptomatic {
		return [asymptomatic,presymptomatic] contains status;
	}
	
	//Action to call to change the status of the individual
	action set_status(string new_status) {
		status <- new_status;
		is_infectious <- define_is_infectious();
		is_infected <- define_is_infected();
		is_asymptomatic <- define_is_asymptomatic();
		if(new_status = symptomatic){ // symptomatic people are supposed to transmit at a higher rate.
			//infectivity = b0 + b1 * y;
			Y <- 1.0;
			infectivity <- B0 + B1;
			estimated_infectivity <- eB0 + eB1; 
			
            // Symptoms for positive cases (US)
            smell_and_taste_loss <- flip(0.67) ? 1.0 : 0.0;
            cough <- flip(0.45) ? 1.0 : 0.0;
            fatigue <- flip(0.23) ? 1.0 : 0.0;
            skipped_meals <- flip(0.47) ? 1.0 : 0.0;
		}
		else if ([asymptomatic,presymptomatic,latent] contains status){ // asymptomatic,latent and presymptomatic people can still transmit.
			Y <- 0.0;
			infectivity <- B0;
			estimated_infectivity <- eB0;
		}
		else{
			Y <- 0.0;
			infectivity <- 0.0;
			estimated_infectivity <- 0.0;
		}
	}
	
	
	//#############################################################
	//Actions
	//#############################################################
	//Initialise epidemiological parameters according to the age of the Individual
	action initialise_epidemio {
		factor_contact_rate_asymptomatic <- world.get_factor_contact_rate_asymptomatic(age);
		factor_contact_rate_wearing_mask <- world.get_factor_contact_rate_wearing_mask(age);
		basic_viral_release <- world.get_basic_viral_release(age);
		contact_rate_human <- world.get_contact_rate_human(age);
		proba_wearing_mask <- world.get_proba_wearing_mask(age);
	}
	
	//Action to call when performing a test on a individual
	action test_individual
	{
		//If the Individual is infected, we check for true positive
		if(self.is_infected)
		{
			if(world.is_true_positive(self.age))
			{
				report_status <- tested_positive;
				total_number_reported <- total_number_reported+1;
			}
			else
			{
				report_status <- tested_negative;
			}
		}
		else
		{
			//If the Individual is not infected, we check for true negative
			if(world.is_true_negative(self.age))
			{
				report_status <- tested_negative;
			}
			else
			{
				report_status <- tested_positive;
				total_number_reported <- total_number_reported+1;
			}
		}
	}

	//Action to call to update wearing a mask for a time step
	action update_wear_mask
	{
		//If the Individual is a free rider, it will not care for masks
		if(free_rider)
		{
			is_wearing_mask <- false;
		}
		else
		{
			if(flip(proba_wearing_mask))
			{
				is_wearing_mask <- true;
			}
			else
			{
				is_wearing_mask <- false;
			}
		}
	}
	
	action update_risk_prob
	{
	  // https://www.nature.com/articles/s41591-020-0916-2
	  float x <- - 1.32 - 0.01 * age + 0.44 * sex + 1.75 * smell_and_taste_loss + 0.31 * cough + 0.49 * fatigue + 0.39 * skipped_meals;	
	  float inf_pred <- 1/(1+exp(-x));
	  float inf_prev <- max(0.2, inf_prevalence);
	  float rs <- 1.0;
	  if (0 < inf_prev) {
	    rs <- inf_pred / inf_prev;
	  }
		
      inf_risk <- 0.0;
      
//	  loop n from: 0 to: length(noninf_prob) - 1 {
      
        int n <- length(noninf_prob) - 1;
	  	float qprod <- 1.0;
	    loop m from: 0 to: n - 1 {
	  	    qprod <- qprod * noninf_prob[m];
	    }
	    qprod <- qprod * (1 - noninf_prob[n]);
	    inf_risk <- inf_risk + qprod;
	     
//	  }

      inf_risk <- min(1.0, rs * inf_risk);
      
      // Shift risks one entry to the left and add new one
      inf_risk_list <- copy_between(inf_risk_list, 1, length(inf_risk_list));
      inf_risk_list <<+ [inf_risk];
		
	  if (0 < inf_risk) {
        string day <- string(int((current_date - starting_date) / #day));
	    save (day + "," + string(cycle) + "," + inf_prevalence + ","  + name + "," + 
	    	  smell_and_taste_loss + "," + cough + "," + fatigue + "," + skipped_meals + "," +
	    	  status + "," + quarantined + "," + inf_risk
	    ) to: "infection_risk.txt" type: "text" rewrite: false;	  	
	  }	
		
	  // Shift probabilities one entry to the left
	  noninf_prob <- copy_between(noninf_prob, 1, length(noninf_prob));
	  noninf_prob <<+ [1.0];
	}
	
	action update_quarantine
	{
		if (quaranatine_enabled) {
			if (quarantined < 1)  {		
	            int obs_idx <- plen - test_delay * 4 - 1;  
	            float obs_risk <- inf_risk_list[obs_idx];		
	            if (quarantine_risk <= obs_risk) {
				    quarantined <- 1.0;
				    quarantine_count <- 0;
			    }
			} else {				
				quarantine_count <- quarantine_count + 1;
				if (quarantine_time * 4 <= quarantine_count) {
					quarantined <- 0.0;
					quarantine_count <- 0;
				} 
			}			
		}
	}
	
	action update_symptoms
	{
		if [susceptible, latent, asymptomatic, presymptomatic] contains status {
	      // 	
          // General symptom prevalence in the population (US)
          smell_and_taste_loss <- flip(0.12) ? 1.0 : 0.0;
          cough <- flip(0.22) ? 1.0 : 0.0;
          fatigue <- flip(0.082) ? 1.0 : 0.0;
          skipped_meals <- flip(0.21) ? 1.0 : 0.0;   
		}
	}
	
	//Action to call to define a new case, obtaining different time to key events
	action define_new_case
	{
		//Add the new case to the total number of infected (not mandatorily known)
		total_number_of_infected <- total_number_of_infected +1;
		//Add the infection to the infections having been caused in the building
		if(building_infections.keys contains(current_place.type))
		{
			building_infections[current_place.type] <- building_infections[current_place.type] +1;
		}
		else
		{
			add 1 to: building_infections at: current_place.type;
		}
		//Add the infection to the infections of the same age
		if(total_incidence_age.keys contains(self.age))
		{
			total_incidence_age[self.age] <- total_incidence_age[self.age] +1;
		}
		else
		{
			add 1 to: total_incidence_age at: self.age;
		}
		
		
		//Set the status of the Individual to latent (i.e. not infectious)
		do set_status(latent);
		//Computation of the latent period as the incubation period for a given age
		self.latent_time <- world.get_incubation_time(self.age);
		//Computation of the serial interval for a given age, used to determine the latent period
		self.serial_interval <- world.get_serial_interval(self.age);
		//Computation of the infectious period for a given age, does not take into account presymptomatic period yet
		self.infectious_time <- world.get_infectious_time(self.age);
		
		
		//If the serial is negative, the Individual is considered being presymptomatic and extends its infectious period, reducing
		// the latent period to something lower than the incubation period
		if(serial_interval<0)
		{
			if(abs(serial_interval)>latent_time)
			{
				serial_interval <- -latent_time;
			}
			self.infectious_time <- max(0,self.infectious_time - self.serial_interval);
			self.latent_time <- max(0,self.latent_time + self.serial_interval);
		}
		
		//Reinitialise the tick to 0
		self.tick <- 0;
	}
	
	//Initialiase social network of the agents (colleagues, friends)
	action initialise_social_network(map<Building,list<Individual>> working_places, map<Building,list<Individual>> schools, map<int,list<Individual>> ind_per_age_cat) {
		
		int nb_friends <- max(0,round(gauss(nb_friends_mean,nb_friends_std)));
		loop i over: ind_per_age_cat.keys {
			if age < i {
				friends <- nb_friends among ind_per_age_cat[i];
				friends <- friends - self;
				break;
			}
		}
		
		if (working_place != nil) {
			int nb_colleagues <- max(0,int(gauss(nb_work_colleagues_mean,nb_work_colleagues_std)));
			if nb_colleagues > 0 {
				colleagues <- nb_colleagues among (working_places[working_place] - self);
			}
		} 
		if (school != nil) {
			int nb_classmates <- max(0,int(gauss(nb_classmates_mean,nb_classmates_std)));
			if nb_classmates > 0 {
				colleagues <- nb_classmates among (schools[school] - self);
			}
		}
 	}
	
	
	//Action to call when entering a new building to update the list of individuals of the buildings
	action enter_building(Building b) {
		if (current_place != nil ){
			current_place.individuals >> self;
		}	
		current_place <- b;
		is_at_home <- current_place = home;
		current_place.individuals << self;
		location <- any_location_in(current_place);
	}
	
	
	//Action to compute the hospitalisation time for symptomatic Individuals
	action set_hospitalisation_time{
		if(world.is_hospitalised(self.age))
		{
			//Compute the time before hospitalisation knowing the current biological status of the agent
			time_before_hospitalisation <- status=presymptomatic? abs(self.serial_interval)+world.get_time_onset_to_hospitalisation(self.age,self.infectious_time):world.get_time_onset_to_hospitalisation(self.age,self.infectious_time);
			if(time_before_hospitalisation>infectious_time)
			{
				time_before_hospitalisation <- infectious_time;
			}
			//Check if the Individual will need to go to ICU
			if(world.is_ICU(self.age))
			{
				//Compute the time before going to ICU once hospitalised
				time_before_ICU <- world.get_time_hospitalisation_to_ICU(self.age, self.time_before_hospitalisation);
				time_stay_ICU <- world.get_time_ICU(self.age);
				if(time_before_hospitalisation+time_before_ICU>=infectious_time)
				{
					time_before_hospitalisation <- infectious_time-time_before_ICU;
				}
			}
		}
	}
	
	
	//#############################################################
	//Reflexes
	//#############################################################
	//Reflex to trigger infection when outside of the commune
	reflex become_infected_outside when: is_outside {
		
//		if flip(proba_outside_contamination_per_hour) {
//			string day <- string(int((current_date - starting_date) /  #day));
//			save ("Day" + day + "-Cycle" + string(cycle) + "="  + self.name + " became infected outside ") to: "contact_data.txt" type: "text" rewrite: false;
//			do define_new_case;
//		}
		
	}
	//Reflex to trigger transmission to other individuals and environmental contamination
	reflex infect_others when: not is_outside and is_infectious
	{
		//Computation of the reduction of the transmission when being asymptomatic/presymptomatic and/or wearing mask
		float reduction_factor <- 1.0;
		if(is_asymptomatic)
		{
			reduction_factor <- reduction_factor * factor_contact_rate_asymptomatic;
		}
		if(is_wearing_mask)
		{
			reduction_factor <- reduction_factor * factor_contact_rate_wearing_mask;
		}
		
		//Performing environmental contamination
		if(current_place!=nil)and(allow_transmission_building)
		{
			ask current_place
			{
				do add_viral_load(reduction_factor*myself.basic_viral_release);
			}
		}
		
		//Perform human to human transmission
		contacts <- [];
		if allow_transmission_human {
			//If the Individual is at home, perform transmission on the household level with a higher factor
			if (is_at_home) {
				// as of now, we are using contact probability but we should use distance function in future to define the kernal K(i,j,t).
				float contact_proba <- contact_rate_human*reduction_factor;				
				list<Individual> new_contacts <- relatives where ((each.is_at_home and flip(contact_proba))and (each.status = susceptible));
				contacts <<+ new_contacts;
				ask new_contacts{
					save("asking relative - "+self+" transmitter - "+myself) to: "contact_data.txt" type: "text" rewrite: false;
					float q <- exp(-self.susceptibility * (1-self.quarantined) * myself.infectivity);
					float eq <- exp(-self.estimated_susceptibility * myself.estimated_infectivity);					
					float transmission_proba <- 1 - q;					
					self.noninf_prob[plen-1] <- self.noninf_prob[plen-1] * eq;
					save("transmission_proba for relative - "+transmission_proba) to: "contact_data.txt" type: "text" rewrite: false;
					if flip(transmission_proba){
						do define_new_case;
						string day <- string(int((current_date - starting_date) /  #day));
						save("transmitting to relatives in the same house") to: "contact_data.txt" type: "text" rewrite: false;
						save ("Day" + day + "-Cycle" + string(cycle) + "="  + myself.name + ":" + myself.location + "Y -"+ myself.Y + "->" + self.name + ":" + self.location + "X -"+ self.X) to: "contact_data.txt" type: "text" rewrite: false;
					}	 			
				}
				if (current_place.nb_households > 1) {
				    contact_proba <- contact_proba * reduction_coeff_all_buildings_inhabitants;
					list<Individual> new_contacts <- current_place.individuals where (each.status = susceptible and flip(contact_proba));
				    contacts <<+ new_contacts;
					ask new_contacts
			 		{
			 			save("asking individual in neighbourhood households - "+self+" transmitter - "+myself) to: "contact_data.txt" type: "text" rewrite: false;
			 			float q <- exp(-self.susceptibility * (1-self.quarantined) * myself.infectivity);
			 			float eq <- exp(-self.estimated_susceptibility * myself.estimated_infectivity);
			 			float transmission_proba <- 1 - q;
					    self.noninf_prob[plen-1] <- self.noninf_prob[plen-1] * eq;
						save("transmission_proba for neighbourhood households- "+transmission_proba) to: "contact_data.txt" type: "text" rewrite: false;
						if flip(transmission_proba){
							do define_new_case;
							string day <- string(int((current_date - starting_date) /  #day));
							save("transmitting to individual in neighbourhood households") to: "contact_data.txt" type: "text" rewrite: false;
							save ("Day" + day + "-Cycle" + string(cycle) + "="  + myself.name + ":" + myself.location + "Y -"+ myself.Y + "->" + self.name + ":" + self.location + "X -"+ self.X) to: "contact_data.txt" type: "text" rewrite: false;
							}
			 		}
				}
				
			}
			else {
				//Perform transmission with people doing the activity explicitly with the Individual
				float contact_proba <- contact_rate_human*reduction_factor;
				list<Individual> fellows <- activity_fellows where (flip(contact_proba) and (each.status = susceptible));
				if (species(last_activity) != Activity) {
					fellows <- fellows where (each.current_place = current_place); 
				}
				contacts <<+ fellows;	
				ask fellows {
					save("asking fellow from the same activity - "+self+" transmitter - "+myself) to: "contact_data.txt" type: "text" rewrite: false;
					float q <- exp(-self.susceptibility * (1-self.quarantined) * myself.infectivity);
					float eq <- exp(-self.estimated_susceptibility * myself.estimated_infectivity);
					float transmission_proba <- 1 - q;
					self.noninf_prob[plen-1] <- self.noninf_prob[plen-1] * eq;
					save("transmission_proba for fellows- "+transmission_proba) to: "contact_data.txt" type: "text" rewrite: false;
					if flip(transmission_proba){
						do define_new_case;
						string day <- string(int((current_date - starting_date) /  #day));
						save("transmitting to fellow from the same activity") to: "contact_data.txt" type: "text" rewrite: false;
						save ("Day" + day + "-Cycle" + string(cycle) + "="  + myself.name + ":" + myself.location + "Y -"+ myself.Y + "->" + self.name + ":" + self.location + "X -"+ self.X) to: "contact_data.txt" type: "text" rewrite: false;
						}
				}
				
				//Perform slightly reduced transmission with people not being involved in the activity but still being present
				contact_proba <- contact_proba * reduction_coeff_all_buildings_individuals;				
				list<Individual> new_contacts <- current_place.individuals where (each.status = susceptible and flip(contact_proba));				
				contacts <<+ new_contacts;
				ask new_contacts
		 		{
		 			save("asking indirect fellows from the same activity place- "+self+" transmitter - "+myself) to: "contact_data.txt" type: "text" rewrite: false;
		 			float q <- exp(-self.susceptibility * (1-self.quarantined) * myself.infectivity);
		 			float eq <- exp(-self.estimated_susceptibility * myself.estimated_infectivity);
					float transmission_proba <- 1 - q;
					self.noninf_prob[plen-1] <- self.noninf_prob[plen-1] * eq;
					save("transmission_proba for indirect fellows- "+transmission_proba) to: "contact_data.txt" type: "text" rewrite: false;
					if flip(transmission_proba){
						do define_new_case;
						string day <- string(int((current_date - starting_date) /  #day));
						save("transmitting to indirect fellows from the same activity place") to: "contact_data.txt" type: "text" rewrite: false;
						save ("Day" + day + "-Cycle" + string(cycle) + "="  + myself.name + ":" + myself.location + "Y -"+ myself.Y + "->" + self.name + ":" + self.location + "X -"+ self.X) to: "contact_data.txt" type: "text" rewrite: false;
						}
		 		}
		 	}
		}
		if (length(contacts) >0){
			loop i from: 0 to: length(contacts) -1 {
//				if (flip(subset_probability)){
					string day <- string(int((current_date - starting_date) /  #day));
					save (day + "," + string(cycle) + ","  + contacts[i].name + "," + contacts[i].X + "," + self.name + "," + self.Y +","+
						self.latent_time/nb_step_for_one_day+","+self.infectious_time/nb_step_for_one_day+"," +
						self.serial_interval/nb_step_for_one_day+","+self.time_before_hospitalisation/nb_step_for_one_day
					) to: "contacts_data_compact.txt" type: "text" rewrite: false;
					
//				}
			}			
		}
	}
	
	//Reflex to trigger to change the status and become infectious once the latent period is expired
	reflex become_infectious when: self.is_latent() and(tick >= latent_time)
	{
		//Set the individual as asymptomatic
		if(world.is_asymptomatic(self.age))
		{
			do set_status(asymptomatic);
			tick <- 0;
		}
		else
		{
			//If the individual was supposed to be infectious before symptom onset, become infectious
			if(serial_interval<0)
			{
				do set_status(presymptomatic);
				tick <- 0;
				do set_hospitalisation_time;
			}
			else
			{
				do set_status(symptomatic);
				tick <- 0;
				do set_hospitalisation_time;
			}
		}
	}
	
	//Reflex to become symptomatic when presymptomatic and presymptomatic period is expired
	reflex become_symptomatic when: status=presymptomatic and self.tick>=abs(self.serial_interval) {
		do set_status(symptomatic);
	}
	
	//Reflex to trigger the need of hospitalisation
	reflex update_before_hospitalisation when: (time_before_hospitalisation>0)
	{
		if(time_before_hospitalisation>0)
		{
			time_before_hospitalisation <-time_before_hospitalisation -1;
			if(time_before_hospitalisation<=0)
			{
				hospitalisation_status <- need_hospitalisation;
			}
		}
	}
	
	//Reflex to trigger the need of ICU
	reflex update_before_ICU when: (hospitalisation_status = need_hospitalisation) and (time_before_ICU>0)
	{
		if(time_before_ICU>0)
		{
			time_before_ICU <-time_before_ICU -1;
			if(time_before_ICU<=0)
			{
				self.hospitalisation_status <- need_ICU;
			}
		}
	}
	
	//Reflex to update the ICU status
	reflex update_stay_ICU when: (time_before_ICU<=0)and(time_stay_ICU>0)and(self.is_ICU=true)
	{
		if(time_stay_ICU>0)
		{
			time_stay_ICU <-time_stay_ICU -1;
			if(time_stay_ICU<=0)
			{
				if(world.is_fatal(self.age))
				{
					hospitalisation_status <- no_need_hospitalisation;
					do set_status(dead);
				}
				else
				{
					hospitalisation_status <- need_hospitalisation;
				}
			}
		}
	}
	
	//Reflex to become recovered or dead once the infectious period is expired
	reflex becomes_not_infectious when: ((status=symptomatic) or (status=asymptomatic))and(tick>=infectious_time)
	{
		//If the Individual is symptomatic, check for its hospitalisation needs
		if(self.status = symptomatic)
		{
			//If the Individual needed to be in ICU and has not been in ICU, dies
			if(self.hospitalisation_status=need_ICU)and(self.is_ICU=false)
			{
				do set_status(dead);
				self.hospitalisation_status <- no_need_hospitalisation;
			}
			else
			{
				//If the Individual needed to be in ICU, has been in ICU, might die according to a given age-related probability
				if(self.hospitalisation_status=need_ICU)and(self.is_ICU=true)and(world.is_fatal(self.age))
				{
					do set_status(dead);
					self.hospitalisation_status <- no_need_hospitalisation;
				}
				else
				{
					do set_status(recovered);
					self.hospitalisation_status <- no_need_hospitalisation;
				}
			}
		}
		else
		{
			do set_status(recovered);
			self.hospitalisation_status <- no_need_hospitalisation;
		}
	}

	//Reflex to execute the agenda	
	reflex execute_agenda {
		pair<Activity,list<Individual>> act <- agenda_week[current_date.day_of_week - 1][current_date.hour];
		if (act.key != nil) {
			if (Authority[0].allows(self, act.key)) {
				int nb_fellows <- Authority[0].limitGroupActivity(self, act.key) - 1;
					if (nb_fellows > 0) {
					activity_fellows <-nb_fellows among act.value;
				} else {
					activity_fellows <- [];
				}
				
				map<Building,list<Individual>> bds_ind <-  act.key.find_target(self);
				if not empty(bds_ind) {
					Building bd <- any(bds_ind.keys);
					list<Individual> inds <- bds_ind[bd];
					activity_fellows <- activity_fellows + inds;
					last_activity <- act.key;
					do enter_building(bd);
					is_outside <- current_place = the_outside;
				} else {
					activity_fellows <- [];
				}
			}
		}
	}
	
	//Reflex to update disease cycle
	reflex update_disease_cycle when:(status!=recovered)and(status!=dead) {
		tick <- tick + 1;
		
//		if(allow_transmission_building and (not is_infected)and(self.current_place!=nil))
//		{
//			if(flip(current_place.viral_load*successful_contact_rate_building))
//			{
//				do define_new_case();
//			}
//		}
		do update_wear_mask();
		do update_symptoms();
		do update_risk_prob();
		do update_quarantine();		
	}

	aspect default {
		if not is_outside {
			draw shape color: status = latent ? #pink : ((status = symptomatic)or(status=asymptomatic)or(status=presymptomatic)? #red : #green);
		}
		else{
			draw circle(10) color: status = latent ? #pink : ((status = symptomatic)or(status=asymptomatic)or(status=presymptomatic)? #red : #green);
		}
	}
}