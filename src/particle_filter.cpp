/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
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

  	if (is_initialized) {
  		return;
  	}

  	std::default_random_engine gen;

  	// create normal distributions for x, y, and theta
  	std::normal_distribution<double> N_x(x, std[0]);
  	std::normal_distribution<double> N_y(y, std[1]);
  	std::normal_distribution<double> N_theta(theta, std[2]);

  	// resize the vectors of particles and weights
  	num_particles = 100;
  	particles.resize(num_particles);
  	//weights.resize(num_particles);

  	// generate the particles
  	for(Particle& p: particles) {
    	p.x = N_x(gen);
    	p.y = N_y(gen);
    	p.theta = N_theta(gen);
    	p.weight = 1;
  	}

  	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

  	std::default_random_engine gen;

  	for(Particle& p: particles){

	    //compute new position and theta
    	double new_theta = 0;
    	double new_x = 0;
    	double new_y = 0;
    	
    	// add measurements to each particle
    	if( fabs(yaw_rate) == 0) {  // constant velocity
    	//if( fabs(yaw_rate) < 0.0001) {  // constant velocity
      		new_x = p.x + velocity * delta_t * cos(p.theta);
      		new_y = p.y + velocity * delta_t * sin(p.theta);
      		new_theta = p.theta;
    	} else {
		    new_theta = p.theta + yaw_rate*delta_t;
     		new_x = p.x + velocity / yaw_rate * ( sin(new_theta) - sin(p.theta) );
      		new_y = p.y + velocity / yaw_rate * ( cos(p.theta) - cos(new_theta) );
    	}

	  	// generate random Gaussian noise
  		std::normal_distribution<double> N_x(new_x, std_pos[0]);
  		std::normal_distribution<double> N_y(new_y, std_pos[1]);
  		std::normal_distribution<double> N_theta(new_theta, std_pos[2]);

    	// predicted particles with added sensor noise
    	p.x = N_x(gen);
    	p.y = N_y(gen);
    	p.theta = N_theta(gen);
  	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

  	for(LandmarkObs& obs: observations) {
    	double minD = std::numeric_limits<float>::max();

    	for(LandmarkObs& pred: predicted) {
      		double distance = dist(obs.x, obs.y, pred.x, pred.y);
      		if( minD > distance){
        		minD = distance;
        		obs.id = pred.id;
      		}
    	}
  	}

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
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

  	for(Particle& p: particles)	{
    	p.weight = 1.0;

    	// step 1: collect valid landmarks
    	vector<LandmarkObs> predictions;
    	for(const Map::single_landmark_s& lm: map_landmarks.landmark_list) {
      		double distance = dist(p.x, p.y, lm.x_f, lm.y_f);
      		if ( distance < sensor_range) { // if the landmark is within the sensor range, save it to predictions
        		predictions.push_back(LandmarkObs{lm.id_i, lm.x_f, lm.y_f});
      		}
    	}

    	// step 2: convert observations coordinates from vehicle to map
    	vector<LandmarkObs> observations_map;
    	double cos_theta = cos(p.theta);
    	double sin_theta = sin(p.theta);

    	for(const LandmarkObs& obs: observations){
      		LandmarkObs tmp;
      		tmp.x = obs.x * cos_theta - obs.y * sin_theta + p.x;
      		tmp.y = obs.x * sin_theta + obs.y * cos_theta + p.y;
      		//tmp.id = obs.id; // maybe an unnecessary step, since the each obersation will get the id from dataAssociation step.
      		observations_map.push_back(tmp);
    	}

    	// step 3: find landmark index for each observation
    	dataAssociation(predictions, observations_map);

	    //record associations
    	vector<int> associations;
    	vector<double> sense_x;
    	vector<double> sense_y;

    	// step 4: compute the particle's weight:
    	// see equation this link:
    	for(LandmarkObs& obs_m: observations_map){

      		Map::single_landmark_s landmark = map_landmarks.landmark_list.at(obs_m.id-1);
      		double x_term = pow(obs_m.x - landmark.x_f, 2) / (2 * pow(std_landmark[0], 2));
      		double y_term = pow(obs_m.y - landmark.y_f, 2) / (2 * pow(std_landmark[1], 2));
      		double w = exp(-(x_term + y_term)) / (2 * M_PI * std_landmark[0] * std_landmark[1]);
      		p.weight *=  w;
	      	
	      	//record association
    	  	associations.push_back(obs_m.id);
      		sense_x.push_back(obs_m.x);
      		sense_y.push_back(obs_m.y);      
    	}

    	weights.push_back(p.weight);
    	SetAssociations(p, associations, sense_x, sense_y);    

  	}

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
  
  	// generate distribution according to weights
  	std::random_device rd;
  	std::mt19937 gen(rd());
  	std::discrete_distribution<> dist(weights.begin(), weights.end());

  	// create resampled particles
  	vector<Particle> resampled_particles;
  	resampled_particles.resize(num_particles);

  	// resample the particles according to weights
  	for(int i=0; i<num_particles; i++) {
    	int idx = dist(gen);
    	resampled_particles[i] = particles[idx];
  	}

  	// assign the resampled_particles to the previous particles
  	particles = resampled_particles;

  	// clear the weight vector for the next round
  	weights.clear();
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
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
