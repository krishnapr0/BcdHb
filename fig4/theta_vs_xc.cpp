/*
 * bcdHbTrials1D.cpp
 *
 *  Created on: May 26, 2021
 *      Author: michael
 */

// In heatmap with A = xc, theta* = -0.05, gamma*/nu = 10^(0.75). 
// Fix gamma/nu* and find how theta* varies with xc (by varying HHalf).

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <random>
#include <chrono>

#include "xoshiro256.hpp"

int main(){
	// Degradation ultimately sets the timescale
	// double nu = std::log(2.0)/(8*60); // Jaeger
	double nu = std::log(2.0)/(8*60);

	// Read in the relevant parameters
	double theta;
	double gammaexp;
	double gamma;
	int numSites = 40;
	bool bursts = false;
	double meanBurst = 1.0;
	double hiavg, hjavg, h2iavg, h2javg;
	double stdi, stdj;
	double sens;

	// Set the constants
	double pGeo = 1/(1+meanBurst);
	double BHalf = 690; // Erdmann
	double hillB = 3.0;
	double HHalf; // Holloway states 820-1300 is biased low (other studies)
	double hillH = 2.0;
	double BLenScale = 120.0; // Erdmann
	double spacing = 8.5; // Erdmann
	int bdyPosition = (numSites-1)/2; // Always put the boundary in the middle, but site goes from 0 to numSites-1
	double BExpAmp = BHalf*std::exp(bdyPosition*spacing/BLenScale);
	int A; // Center and amplitude of the tanh

	std::vector<int> HHalfs = {164, 328, 492, 656, 820, 984, 1148, 1312};

	// Set the simulation variables
	double t, dt; // time and increment
	double r; // random number for determining reaction
	int source, target, rxnIndex;
	double propSum, partSum;
	double propVec[3*numSites];

	// How long are the simulations? How many?
	int numTrials = 1000;
	double simTime = 10*25*60*4; // Cell cycle 13 is 25 minutes... the first part of cell cycle 14 is 75 minutes.
	// Running for 1000 minutes

	// Create the generator
	xoshiro256p gen;

	// Create the file
	std::ofstream strm;
	strm.open("hmap_theta_vs_xc.csv", std::ofstream::out|std::ofstream::trunc);

	for(theta = 0.3; theta > -0.31; theta -= 0.05) {
		gammaexp = -0.75;
		gamma = std::pow(10, gammaexp)*nu;

		for (int i  = 0;i < (int)HHalfs.size(); i++) {
			HHalf = HHalfs[i];
			
			// Use formulas to set the other parameters
			double xc = HHalf*std::pow((hillH-1)/(hillH+1),1/hillH);  // Modified, no xc-1
			double kH = 16*hillH*nu*xc/((hillH*hillH - 1)*((hillH*hillH - 1)*theta + 4)); // Modified from Mike's original
			double kB = (nu*xc -(hillH-1)*kH/(2*hillH))*2; // Doesn't change

			// Construct the molecule number profiles, set the Bicoid profile and initial Hunchback profile
			int BcdTemp;
			double BProdProp[numSites];
			int hNum[numSites];
			int hInitial[numSites];
			int hbMidInitial = (int)std::round(xc); // Tanh centered at xc
			for(int site=0;site<numSites;site++){
				BcdTemp = (int)std::round(BExpAmp*std::exp(-site*spacing/BLenScale));
				BProdProp[site] = kB*std::pow(BcdTemp/BHalf,hillB)/(std::pow(BcdTemp/BHalf,hillB)+1);
				hInitial[site] = (int)std::round(hbMidInitial*(1-std::tanh(spacing*(site-bdyPosition)/BLenScale))); // Tanh initialization
			}

			hiavg = 0.0;
			hjavg = 0.0;
			h2iavg = 0.0;
			h2javg = 0.0;
			sens = 0.0;

			// Loop over the trials
			for(int trial=0;trial<numTrials;trial++){
				// Initialize the system
				for(int site=0;site<numSites;site++){
					hNum[site] = hInitial[site];
				}
				t=0;
				dt=0;
				// Set the vector of propensities and the sum of all propensities
				propSum = 0;
				for(int site=0;site<numSites;site++){
					// Production of Hb
					propVec[site] = BProdProp[site] + kH*std::pow(hNum[site]/HHalf,hillH)/(std::pow(hNum[site]/HHalf,hillH)+1);
					// Degradation of Hb
					propVec[site+numSites] = nu*hNum[site];
					// Hb diffusion
					if((site==0)||(site==(numSites-1))){
						propVec[site+2*numSites] = gamma*hNum[site];
					} else {
						propVec[site+2*numSites] = 2*gamma*hNum[site];
					}
					propSum += propVec[site] + propVec[site+numSites] + propVec[site+2*numSites];
				}
				// Perform the trial
				while(t < simTime) {                                                 // Long run to steady state
					// draw a time increment and update
					dt = gen.exponential(1.0/propSum);
					t += dt;
					// determine which reaction occurs
					r = gen.uniform(0.0, 1.0);
					rxnIndex = 0;
					partSum = propVec[0];
					while(r > partSum/propSum){
						rxnIndex++;
						partSum += propVec[rxnIndex];
					}
					// Indentify the source site
					source = rxnIndex % numSites;
					// Determine whether it was birth, death, or diffusion and change things appropriately
					if(rxnIndex/numSites ==0){
						// Hb is produced. Subtract the relevant terms from propSum first
						propSum -= propVec[source]+propVec[source+numSites]+propVec[source+2*numSites];
						// Update the molecule numbers
						if(bursts){
							hNum[source] += gen.geometric(pGeo);
						} else {
							hNum[source]++;
						}
						// Update propVec, then the sum
						propVec[source] = BProdProp[source] + kH*std::pow(hNum[source]/HHalf,hillH)/(std::pow(hNum[source]/HHalf,hillH)+1);
						propVec[source+numSites] = nu*hNum[source];
						if((source==0)||(source==(numSites-1))){
							propVec[source+2*numSites] = gamma*hNum[source];
						} else {
							propVec[source+2*numSites] = 2*gamma*hNum[source];
						}
						propSum += propVec[source]+propVec[source+numSites]+propVec[source+2*numSites];
					} else if (rxnIndex/numSites ==1){
						// Hb is produced. Subtract the relevant terms from propSum first
						propSum -= propVec[source]+propVec[source+numSites]+propVec[source+2*numSites];
						hNum[source]--;
						// Update propVec, then the sum
						propVec[source] = BProdProp[source] + kH*std::pow(hNum[source]/HHalf,hillH)/(std::pow(hNum[source]/HHalf,hillH)+1);
						propVec[source+numSites] = nu*hNum[source];
						if((source==0)||(source==(numSites-1))){
							propVec[source+2*numSites] = gamma*hNum[source];
						} else {
							propVec[source+2*numSites] = 2*gamma*hNum[source];
						}
						propSum += propVec[source]+propVec[source+numSites]+propVec[source+2*numSites];
					} else {
						// Hb diffuses. Which direction does it go? Check for edge cases, flip a coin if in the middle
						if(source==0){
							target=1;
						} else if (source==(numSites-1)){
							target =numSites-2;
						} else {
							if(gen.uniform(0.0,1.0)<0.5){
								target=source-1;
							} else {
								target = source+1;
							}
						}
						// Subtract the relevant terms from propSum
						propSum -= propVec[source]+propVec[source+numSites]+propVec[source+2*numSites];
						propSum -= propVec[target]+propVec[target+numSites]+propVec[target+2*numSites];
						// Update the molecule numbers
						hNum[source]--;
						hNum[target]++;
						// Update propVec, then the sum. First the source
						propVec[source] = BProdProp[source] + kH*std::pow(hNum[source]/HHalf,hillH)/(std::pow(hNum[source]/HHalf,hillH)+1);
						propVec[source+numSites] = nu*hNum[source];
						if((source==0)||(source==(numSites-1))){
							propVec[source+2*numSites] = gamma*hNum[source];
						} else {
							propVec[source+2*numSites] = 2*gamma*hNum[source];
						}
						// Next the target
						propVec[target] = BProdProp[target] + kH*std::pow(hNum[target]/HHalf,hillH)/(std::pow(hNum[target]/HHalf,hillH)+1);
						propVec[target+numSites] = nu*hNum[target];
						if((target==0)||(target==(numSites-1))){
							propVec[target+2*numSites] = gamma*hNum[target];
						} else {
							propVec[target+2*numSites] = 2*gamma*hNum[target];
						}
						// Finally, update the sum
						propSum += propVec[source]+propVec[source+numSites]+propVec[source+2*numSites];
						propSum += propVec[target]+propVec[target+numSites]+propVec[target+2*numSites];
					}
				}
				
				hiavg += hNum[bdyPosition];
				hjavg += hNum[bdyPosition - 1];
				h2iavg += std::pow(hNum[bdyPosition],2);
				h2javg += std::pow(hNum[bdyPosition - 1],2);
			}

			hiavg /= numTrials;
			hjavg /= numTrials;
			h2iavg /= numTrials;
			h2javg /= numTrials;
			stdi = std::pow(h2iavg - std::pow(hiavg, 2), 0.5);
			stdj = std::pow(h2javg - std::pow(hjavg, 2), 0.5);

			sens = 2*(hjavg - hiavg)/(stdi + stdj);

			strm << sens << ",";
		}

		strm << std::endl;
	}
	
	strm.close();

	// terminate the program
	return 0;
}
