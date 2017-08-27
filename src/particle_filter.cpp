/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random> // Need this for sampling from distributions
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	
	num_particles = 1; // No. of particles. Keep size such that the filter deos not become slow nor less accurate
	
	default_random_engine gen;
	
	// This line creates a normal (Gaussian) distribution for x.
	normal_distribution<double> dist_x(x, std[0]);	
	// This line creates a normal (Gaussian) distribution for y.
	normal_distribution<double> dist_y(y, std[1]);
	// This line creates a normal (Gaussian) distribution for theta.
	normal_distribution<double> dist_theta(theta, std[2]);

	for(int i=0; i < num_particles; i++){
		
		Particle particle;
		particle.id = i;
		particle.x = dist_x(gen);
		particle.y = dist_y(gen);
		particle.theta = dist_theta(gen);
		particle.weight = 1.0;

		particles.push_back(particle);
		weights.push_back(1);
	}

	is_initialized = 1;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine gen;

	for(int i=0; i < num_particles; i++){

		double new_x;
		double new_y;
		double new_theta;

		if(fabs(yaw_rate) < 0.00001){
			new_x = particles[i].x + velocity*delta_t*cos(particles[i].theta);
			new_y = particles[i].y + velocity*delta_t*sin(particles[i].theta);
			new_theta = particles[i].theta;
		}
		else{
			new_x = particles[i].x + velocity/yaw_rate*(sin(particles[i].theta+yaw_rate*delta_t) - sin(particles[i].theta));
			new_y = particles[i].y + velocity/yaw_rate*(cos(particles[i].theta) - cos(particles[i].theta+yaw_rate*delta_t));
			new_theta = particles[i].theta + yaw_rate*delta_t;
		}

		// This line creates a normal (Gaussian) distribution for x.
		normal_distribution<double> dist_x(new_x, std_pos[0]);	
		// This line creates a normal (Gaussian) distribution for y.
		normal_distribution<double> dist_y(new_y, std_pos[1]);
		// This line creates a normal (Gaussian) distribution for theta.
		normal_distribution<double> dist_theta(new_theta, std_pos[2]);

		particles[i].x = dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_theta(gen);

	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html	

	double obs_x, obs_y;

	for(int i=0; i < num_particles; i++){
		
		cout << "" << endl;
		cout << "NEW PARTICLE" << endl;
		cout << "i = " << i << endl;
		
		particles[i].weight = 1.0;

		for(int j = 0; j < observations.size(); j++){
			
			//cout << "j = " << j << endl;

			// Transforming observations from VEHICLE cordinates to MAP cordinates
			if(particles[i].theta > 0){ //Counterclockwise
				
				observations[j].x = particles[i].x + cos(particles[i].theta)*observations[j].x - sin(particles[i].theta)*observations[j].y;
				observations[j].y = particles[i].y + sin(particles[i].theta)*observations[j].x + cos(particles[i].theta)*observations[j].y;
			}
			else{						//Clockwise

				observations[j].x = particles[i].x + cos(particles[i].theta)*observations[j].x + sin(particles[i].theta)*observations[j].y;
				observations[j].y = particles[i].y - sin(particles[i].theta)*observations[j].x + cos(particles[i].theta)*observations[j].y;	
			}

			obs_x = observations[j].x;
			obs_y = observations[j].y;

			// Association of an observation with its landmark
			double new_dist = 0.0;			
			double min_dist = 0.0;

			double mu_x = 0.0;
			double mu_y = 0.0;

			int flag = 1;

			for(int k=0; k < map_landmarks.landmark_list.size(); k++){	

				//cout << "k = " << k << endl;
				
				double landmk_x = map_landmarks.landmark_list[k].x_f;
				double landmk_y = map_landmarks.landmark_list[k].y_f;
						
				if(fabs(landmk_x - particles[i].x) <= sensor_range){

					if(fabs(landmk_y - particles[i].y) <= sensor_range){
			
						new_dist = dist(obs_x, obs_y, landmk_x, landmk_y);						
						
						if(flag == 1){
						//	cout << "min set" << endl;
							min_dist = new_dist;
							mu_x = landmk_x;
							mu_y = landmk_y;
							flag = 0;
						}

						if(new_dist < min_dist){
							min_dist = new_dist;							
							mu_x = landmk_x;
							mu_y = landmk_y;
							//cout << "updated min = " << min_dist << endl;							
						}
					}					
				}
			}
			// Particle weight calc.
			// The weight of the particle is calculated as the product of each measurement's Multi-variate Gaussian  probability. 
			// Hence for every particle, calc Multi-variate Gaussian probability for each transformed observation (formula below) 
			// using its associated landmark. Then after that is calc for each transformed observation, multiply all 
			// these probababilites to get the final weight of the particle.
			// Final particle weight = prob(obs1) * prob(obs2)* ..... *prob(obsn)
			//
			// FORMULA:
			//# calculate normalization term
			// gauss_norm= (1/(2 * pi * sig_x * sig_y))
			//
			// # calculate exponent
			// exponent= ((x_obs - mu_x)**2)/(2 * sig_x**2) + 
			//					((y_obs - mu_y)**2)/(2 * sig_y**2)
			//
			// # calculate weight using normalization terms and exponent
			// weight= gauss_norm * exp(-exponent)			

			double gauss_norm = (1/(2 * M_PI * std_landmark[0] * std_landmark[1]));
			double x_exponent = (pow((obs_x - mu_x),2)) / (2*pow(std_landmark[0],2)); 
			double y_exponent = (pow((obs_y - mu_y),2)) / (2*pow(std_landmark[1],2));
			double exponent = x_exponent + y_exponent;
			particles[i].weight = particles[i].weight * gauss_norm * exp(-exponent);
			
			cout << " obs x, obs_y " << obs_x << "," << obs_y << endl;
			cout << " mu x, mu_y " << mu_x << "," << mu_y << endl;
			cout << "x_exponent, y_exponent " << x_exponent << "," << y_exponent << endl;	
			cout << "weight" << gauss_norm * exp(-exponent) << endl;
			cout << "" << endl;
		}		
	}	
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	default_random_engine gen;
	discrete_distribution<int> distribution(weights.begin(), weights.end());

	vector<Particle> resample_particles;

	for(int i=0; i < num_particles; i++){
		resample_particles.push_back(particles[distribution(gen)]);
	}

	particles = resample_particles;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
