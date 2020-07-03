/***
* Part of the GAMA CoVid19 Modeling Kit
* see http://gama-platform.org/covid19
* Author: Huynh Quang Nghi
* Tags: covid19,epidemiology
***/

model CoVid19

import "../../Model/Global.gaml"
import "../Abstract Experiment.gaml"

global {
    
	action define_policy{   
		/*ask Authority {
			policy <- create_lockdown_policy();
		}*/
	}
}

experiment "No Containment" parent: "Abstract Experiment" {
	output {
		display "Main" parent: default_display {
		}
		display "Chart" parent: default_white_chart {
		}
		display "Cumulative incidence" parent: cumulative_incidence {
		}
	}
	reflex save_results {		
		int inf <- length(Individual where (each.is_latent())) + length(Individual where (each.is_infectious));
		int tot <- length(Individual where (each.status=susceptible)) + 
		             length(Individual where (each.is_latent())) + 
		             length(Individual where (each.is_infectious));
		inf_prevalence <- float(inf)/float(tot);
		
		string day <- string(int((current_date - starting_date) /  #day));
		save(day + "," + cycle + "," + length(Individual where (each.status=susceptible)) + "," + length(Individual where (each.is_latent())) + "," + length(Individual where (each.is_infectious)) + "," + length(Individual where (each.status = recovered)) + "," + length(Individual where (each.status = dead)) ) to: "infected_number.txt" type: "text" rewrite: false;
	}
}