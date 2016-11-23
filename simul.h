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
 * simul.h
 *
 * Author: Alexis Poncet
 * Email: alexis.poncet@ens.fr
 * Version: 1.0
 *
 * Define the functions of simul.cpp.
 */

#ifndef BROWNIANTRACERS_SIMUL_H_
#define BROWNIANTRACERS_SIMUL_H_

#include <vector>
#include <forward_list>
#include <array>
#include <random>
#include <iostream>
#include <fstream>

#define DEFAULT_NB_RND_CHAR 6
#define DEFAULT_OUT_PRECISION 15

int runSimulation(const Parameters &p);

// Core functions for simulation
void calcForcesBetweenParticlesNaive(
    const std::vector< std::array<double, DIM> > &positions,
    std::vector< std::array<double, DIM> > &forces, const long nbParticles,
    const double boxLength);

void calcForcesBetweenParticles(
    const std::vector< std::array<double, DIM> > &positions,
    std::vector< std::array<double, DIM> > &forces,
    const long nbParticles, const double boxLength,
    const double subBoxesLength, const long nbSubBoxesLine,
    const long nbSubBoxes,
    const std::vector< std::forward_list<long> > &nbrLists);

void makeSubBoxes(const std::vector< std::array<double, DIM> > &particles,
                  std::vector<long> &boxOfParticle,
                  std::vector< std::forward_list<long> > &particlesOfBox,
                  const long nbParticles, const double subBoxLength,
                  const long nbSubBoxesLine);

void initCorrel(std::vector<long long> &countInDivTNT,
                std::vector<long long> &countInDivTT,
                long &nbDivLine, double &divLength, const double boxSize);

void calcCorrel(const std::vector< std::array<double, DIM> > &pos,
                std::vector<long long> &countInDivTNT,
                std::vector<long long> &countInDivTT,
                const long nbParticles, const long nbTracers,
                const double boxSize, const long nbDivLine);

std::forward_list<long> makeNeighborsList(const long L, const long i);

std::vector< std::forward_list<long> > makeNeighborsLists(const long L);


// Functions concerning export
void printPos(const std::vector< std::array<double, DIM> > &pos);

void exportLineMobility(const std::vector< std::array<double, DIM> > &pos,
                        const std::vector< std::array<double, DIM> > &forces,
                        const std::vector< std::array<double, DIM> > &windings,
                        const std::vector<double> initPosX,
                        const long nbParticles, const long nbTracers,
                        const double boxLength, const double time0,
                        std::ostream &stream);

void exportCorrel(const std::vector<long long> &countInDivTNT,
                  const std::vector<long long> &countInDivTT,
                  const long nbDivLine, const double divLength,
                  const double boxLength, std::ostream &stream);


// Functions taking care of files
int openFiles(const Parameters &p, std::ostream *&fileMobility,
              std::ostream *&fileCorrel, std::mt19937 &rndGen);

int openFile(const std::string &fname, std::ostream *&file, bool gzip,
             bool verbose);

void closeFiles(const Parameters &p, std::ostream *&fileMobility,
                std::ostream *&fileCorrel);

void genFileNames(const Parameters &p, std::string &fnameMobility,
                  std::string &fnameCorrel, std::mt19937 &rndGen);

void writeHeaders(const Parameters &p, std::ostream &fileMobility,
                  std::ostream &fileCorrel);

#endif  // BROWNIANTRACERS_SIMUL_H_
