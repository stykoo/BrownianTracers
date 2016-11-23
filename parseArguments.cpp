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
 * parseArguments.cpp
 *
 * Author: Alexis Poncet
 * Email: alexis.poncet@ens.fr
 * Version: 1.0
 *
 * Contain the function for parsing and checking command-line arguments.
 */

#include <iostream>
#include <fstream>
#include <string>
#include <boost/program_options.hpp>
#include "parseArguments.h"

using namespace std;
namespace po = boost::program_options;

// Parse command line arguments, store the values into structure Parameters
// and check if they are consistent.
// Return 1 if displaying help, 2 if an inconsistency is found, 0 otherwise.
int parseAndCheckArguments(int argc, char **argv, Parameters &p) {
    // Command-line only
    po::options_description genericOpts("Command-line only");
    genericOpts.add_options()
        ("config,c", po::value<string>(), "Use configuration file")
        ("verbose,v", po::bool_switch(&p.verbose), "Verbose mode")
        ("help,h", "Print help message and exit")
        ;

    // Both command line (prioritary) and configuration file
    po::options_description configOpts(
        "Configuration (configuration file can be used)");
    configOpts.add_options()
        ("particles,n", po::value<long>(&p.nbParticles)->required(),
         "Number of particles")
        ("tracers,m", po::value<long>(&p.nbTracers)->required(),
         "Number of tracers")
        ("density,r", po::value<double>(&p.density)->required(),
         "Density of the system")
        ("temperature,T", po::value<double>(&p.temperature)->required(),
         "Temperature")
        ("force,f", po::value<double>(&p.force)->required(), "Force")
        ("timestep,t", po::value<double>(&p.timestep)->required(), "Time step")
        ("iters,i", po::value<long>(&p.nbIters)->required(),
         "Number of time iterations")
        ("thermalization,j", po::value<long>(&p.nbItersTh)->default_value(
            DEFAULT_ITERS_TH),
         "Number of time iterations for thermalization")
        ("resolution,R",
         po::value<double>(&p.divLength)->default_value(DEFAULT_DIV_LENGTH),
         "Resolution for correlations")
        ("skip,s", po::value<long>(&p.timeCalcCorrel)->default_value(
             DEFAULT_TIME_CALC_CORREL),
         "Iterations between computations of correlations")
        ("mobility", po::bool_switch(&p.outputMobility), "Output mobility")
        ("correl", po::bool_switch(&p.outputCorrel), "Output correlations")
        ("livepos", po::value<long>(&p.outputPosIters)->default_value(-1),
         "Output positions of particles every given number of iterations")
        ("path,p", po::value<string>(&p.path)->default_value(DEFAULT_PATH),
         "Directory to export results")
        ("fileMobility",
          po::value<string>(&p.fnameMobility)->default_value(""),
         "File name (without path) for mobility (default format if none)")
        ("fileCorrel", po::value<string>(&p.fnameCorrel)->default_value(""),
         "File name (without path) for correlations (default format if none)")
        ("filePos", po::value<string>(&p.fnamePos)->default_value(""),
         "File name (without path and extension) for positions "
         "(default format if none)")
        ("gzip,z", po::bool_switch(&p.gzip), "Compress output using gzip")
        ;

    po::options_description cmdlineOpts;
    cmdlineOpts.add(configOpts).add(genericOpts);

    try {
        po::variables_map vars;
        po::store(po::parse_command_line(argc, argv, cmdlineOpts), vars);

        // Display help and exit
        if (vars.count("help")) {
            cout << "Usage: " << argv[0] << " options\n";
            cout << cmdlineOpts << endl;
            return 1;
        }

        // Use configuration file if one is defined
        if (vars.count("config")) {
            ifstream confFile(vars["config"].as<string>());
            if (confFile.is_open()) {
                po::store(po::parse_config_file(confFile, configOpts), vars);
            } else {
                cerr << "Warning: could not open configuration file "
                    << vars["config"].as<string>() << endl;
            }
        }

        po::notify(vars);
    } catch (std::exception& e) {
        std::cerr << "Error: " << e.what() << endl;
        return 2;
    }

    // Look for inconstencies
    if (checkParameters(p)) {
        return 2;
    }

    return 0;
}

// Check if parameters are consistent. Return 0 if they are, 1 else.
int checkParameters(const Parameters &p) {
    if (checkPositive(p.nbParticles, "the number of particles") ||
        checkPositive(p.density, "the density") ||
        checkPositive(p.temperature, "the temperature") ||
        checkPositive(p.force, "the force") ||
        checkPositive(p.timestep, "the time step") ||
        checkPositive(p.nbIters, "the number of iterations") ||
        checkPositive(p.divLength, "the resolution") ||
        checkPositive(p.timeCalcCorrel, "skip")) {
        return 1;
    }

    if (p.nbTracers < 0 || p.nbTracers > p.nbParticles) {
        cerr << "Error: the number of tracers should be between 0 and"
            << " the number of particles!" << endl;
        return 1;
    }

#ifndef VISU2D
    if (!p.outputMobility && !p.outputCorrel && p.outputPosIters < 1) {
        cerr << "Warning: you should use either '--mobility', '--correl "
            "or '--livepos'"
            << endl;
        return 1;
    }
#endif

    return 0;
}

// Return 1 and print error if a is negative. Return 0 otherwise.
int checkPositive(const long a, const string label) {
    if (a <= 0) {
        cerr << "Error: " << label << " should be strictly positive."
            << std::endl;
        return 1;
    }
    return 0;
}
int checkPositive(const double a, const string label) {
    if (a <= 0.) {
        cerr << "Error: " << label << " should be strictly positive."
            << std::endl;
        return 1;
    }
    return 0;
}

// Print parameters to standard output.
void printParameters(const Parameters &p) {
    cout << p.nbParticles << " particles, " << p.nbTracers
        << " tracers" << endl;
    cout << "Density: " << p.density << ", temperature: " << p.temperature
        << ", force: " << p.force << endl;
    cout << "Timestep: " << p.timestep << ", iterations: " << p.nbIters
        << ", thermalization: " << p.nbItersTh << endl;
    if (p.path != DEFAULT_PATH) {
        cout << "path: " << p.path << endl;
    }
    if (p.gzip) {
        cout << "gzip compression enabled" << endl;
    }

    if (p.outputMobility) {
        cout << "- Output mobility";
        if(!p.fnameMobility.empty()){
            cout << ", file: " << p.fnameMobility;
        }
        cout << endl;
    }
    if (p.outputCorrel) {
        cout << "- Output correlations, resolution: " << p.divLength
           << ", skip: " << p.timeCalcCorrel;
        if(!p.fnameCorrel.empty()){
            cout << ", file: " << p.fnameCorrel;
        }
        cout << endl;
    }
#ifdef VISU2D
    cout << "- Visualization in 2d" << endl;
#endif
    cout << endl;
}
