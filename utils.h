/*
Copyright ESPCI Paris (2016)
Contributors: Alexis Poncet, Vincent Démery

alexis.poncet@ens.fr

This software is a computer program whose purpose is to simulate the
motion of interacting Brownian particle, some of them being driven by
an external force.

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/
/*
 * BrownianTracers
 * utils.h
 * 
 * Author: Alexis Poncet
 * Email: alexis.poncet@ens.fr
 * Version: 1.0
 *
 * Define helper functions, parameter data structure, and dimension if needed.
 */

#ifndef BROWNIANTRACERS_UTILS_H_
#define BROWNIANTRACERS_UTILS_H_

#include <cmath>
#include <string>

struct Parameters {
	long nbParticles;  // Number of particles
	long nbTracers;  // Number of tracers
	double density;  // Total density
	double temperature;  // Temperature
	double force;  // Force
	double timestep;  // Time step
	long nbIters;  // Simulation duration
	long nbItersTh;  // Thermalization duration

	double divLength;  // Resolution for correlations
	long timeCalcCorrel;  // Iterations between computations of correlations

	bool outputMobility, outputCorrel;  // Output mobility / correlations
	std::string path;  // Directory to store the results
	std::string fnameMobility, fnameCorrel;  // Filenames for mobility / correl

	bool gzip;  // Compress the output using gzip
	bool verbose;  // Verbose mode
};

// Translate x into interval [-L/2, L/2[ 
inline double periodicBC(const double x, const double L) {
	return x - L * std::round(x / L);
}

// Translate x into interval [0, L[ 
inline double periodicBCpos(const double x, const double L) {
	return x - L * std::round(x / L) + L / 2;
}

// Compute the number of translations of length L (to the left) 
// to put x in interval [-L/2, L/2[
inline long calcWinding(const double x, const double L) {
	return std::lround(x/L);
}

// Compute pow(a, b) with b a positive integer
template<typename T, typename U>
int mypow(const T a, const U b) {
	if (b == 2)  // Most common case
		return a*a;
	if (b == 1)
		return a;
	if (b <= 0)
		return 1;
	if (b % 2 == 1)
		return a*mypow(a,b-1);
	T c = mypow(a,b/2);
	return c*c;
}

#endif // BROWNIANTRACERS_UTILS_H_
