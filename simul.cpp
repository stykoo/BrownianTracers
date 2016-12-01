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
 * simul.cpp
 *
 * Author: Alexis Poncet
 * Email: alexis.poncet@ens.fr
 * Version: 1.0
 *
 * Contain the functions for simulating the motion of the particles.
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <iomanip>
#include <ctime>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include "utils.h"
#include "simul.h"
#ifdef VISU2D
#include <thread>
#include "visu2d.h"
#endif

using namespace std;
namespace io = boost::iostreams;

// Run the simulation.
int runSimulation(const Parameters &p) {
// Beginning of initializations
    const double boxLength = pow(p.nbParticles / p.density,
                                 1.0 / ((double) DIM));  // Length of the box
    double divLength = p.divLength;  // Can be modify to divide evenly the box
    int time0 = 0;

    // Length of a sub-box (approximately 1.)
    const double subBoxLength = boxLength / floor(boxLength);
    // Number of sub-boxes in a given direction
    const long nbSubBoxesLine = (long) floor(boxLength);
    // Total number of sub-boxes
    const long nbSubBoxes = mypow(nbSubBoxesLine, DIM);
    // Lists of neighboring boxes of each box
    const vector< forward_list<long> > nbrLists = \
        makeNeighborsLists(nbSubBoxesLine);

    // Random generators
    random_device rd;
    mt19937 rndGen(rd());
    uniform_real_distribution<double> rndForPos(-boxLength/2., boxLength/2.);
    normal_distribution<double> rndForNoise(0.,sqrt(2. * p.temperature
                                                    * p.timestep));
    // Files for export
    ostream *fileMobility = nullptr, *fileCorrel = nullptr;
    string fullnamePos;
    if(openFiles(p, fileMobility, fileCorrel, fullnamePos, rndGen)){
        return 1;
    }
    writeHeaders(p, *fileMobility, *fileCorrel);

    // Random initial positions
    vector< array<double,DIM> > positions(p.nbParticles);
    for (long i = 0 ; i < p.nbParticles ; i++) {
        for (int a = 0 ; a < DIM ; a++) {
            positions[i][a] = rndForPos(rndGen);
        }
    }

    // Thread for visualization in 2d
#ifdef VISU2D
    DataThreadVisu visuData;
    visuData.pos = &positions;
    visuData.t = &time0;
    thread visuThread(threadVisu, visuData, p.nbParticles, p.nbTracers,
                      boxLength, p.timestep);
#endif

    vector< array<double,DIM> > forces(p.nbParticles);
// End of initializations

// Thermalization loop
    for (int t = 0 ; t < p.nbItersTh ; t++) {
        calcForcesBetweenParticles(positions, forces, p.nbParticles,
                                   boxLength, subBoxLength, nbSubBoxesLine,
                                   nbSubBoxes, nbrLists);

        // Force on the tracers (along the first coordinate)
        for (long i = 0 ; i < p.nbTracers ; i++)
            positions[i][0] += p.timestep * p.force;

        // Add forces, noise and enforce periodic boundary conditions
        for (long i = 0 ; i < p.nbParticles ; i++) {
            for (int a = 0 ; a < DIM ; a++) {
                positions[i][a] += p.timestep * forces[i][a];
                positions[i][a] += rndForNoise(rndGen);
                positions[i][a] = periodicBC(positions[i][a], boxLength);
            }
        }
    }
// End of thermalization loop

// Beginning of other intializations
    vector<long long> countInDivTNT, countInDivTT;
    long nbDivLine;
    vector<double> initPosX(p.nbTracers);
    vector< array<double,DIM> > windings(p.nbParticles);

    if (p.outputMobility) {
        // Initial position of the tracers along the first coordinate
        for (long i = 0 ; i < p.nbTracers ; i++)
            initPosX[i] = positions[i][0];

        // Number of turns for each particles
        for (long i = 0 ; i < p.nbParticles ; i++)
            for (int a = 0 ; a < DIM ; a++)
                windings[i][a] = 0;
    }

    if (p.outputCorrel) {
        initCorrel(countInDivTNT, countInDivTT, nbDivLine, divLength,
                   boxLength);
    }
// End of other initializations

// Main loop
    for (time0 = 0 ; time0 < p.nbIters ; time0++) {
        calcForcesBetweenParticles(positions, forces, p.nbParticles,
                                   boxLength, subBoxLength, nbSubBoxesLine,
                                   nbSubBoxes, nbrLists);

        // Force on the tracers (along the first coordinate)
        for (long i = 0 ; i < p.nbTracers ; i++)
            positions[i][0] += p.timestep * p.force;

        // Add forces, noise and enforce periodic boundary conditions
        for (long i = 0 ; i < p.nbParticles ; i++) {
            for (int a = 0 ; a < DIM ; a++) {
                positions[i][a] += p.timestep * forces[i][a];
                positions[i][a] += rndForNoise(rndGen);
                windings[i][a] += calcWinding(positions[i][a], boxLength);
                positions[i][a] = periodicBC(positions[i][a], boxLength);
            }
        }

       if (p.outputMobility) {
           exportLineMobility(positions, forces, windings, initPosX,
                              p.nbParticles, p.nbTracers, boxLength,
                              (time0 + 1) * p.timestep, *fileMobility);
        }
        if (p.outputCorrel && time0 % p.timeCalcCorrel == 0) {
            calcCorrel(positions, countInDivTNT, countInDivTT, p.nbParticles,
                       p.nbTracers, boxLength, nbDivLine);
        }
        if (p.outputPosIters > 0 && time0 % p.outputPosIters == 0) {
            exportPos(p, positions, fullnamePos, time0 / p.outputPosIters);
        }

    }
// End of main loop

    if (p.outputCorrel) {
        exportCorrel(countInDivTNT, countInDivTT, nbDivLine, divLength,
                     boxLength, *fileCorrel);
    }

    closeFiles(p, fileMobility, fileCorrel);

#ifdef VISU2D
    visuThread.join();
#endif

    return 0;
}


// Compute inter-particle forces (naive algorithm).
void calcForcesBetweenParticlesNaive(
        const vector< array<double,DIM> > &positions,
        vector< array<double,DIM> > &forces, const long nbParticles,
        const double boxLength) {
    for (long i = 0 ; i < nbParticles ; i++) {
        for (int a = 0 ; a < DIM ; a++) {
            forces[i][a] = 0;
        }
    }

    for (long i = 0 ; i < nbParticles ; i++) {
        for (long j = i+1 ; j<nbParticles ; j++) {
            double dr[DIM];
            double distsq = 0.;
            for (int a = 0 ; a < DIM ; a++) {
                dr[a] = periodicBC(positions[i][a] - positions[j][a],
                                   boxLength);
                distsq += dr[a] * dr[a];
            }
            if(distsq < 1. && distsq > 0.){
                for (int a = 0 ; a < DIM ; a++) {
                    double t = (1. / sqrt(distsq) - 1.) * dr[a];
                    forces[i][a] += t;
                    forces[j][a] -= t;  // f_ij = -f_ji
                }
            }
        }
    }
}

// Compute inter-particle forces (with decomposition in sub-boxes).
void calcForcesBetweenParticles(
        const vector< array<double,DIM> > &positions,
        vector< array<double,DIM> > &forces, const long nbParticles,
        const double boxLength, const double subBoxLength,
        const long nbSubBoxesLine, const long nbSubBoxes,
        const vector< forward_list<long> > &nbrLists) {
    for (long i = 0 ; i < nbParticles ; i++)
        for (int a = 0 ; a < DIM ; a++)
            forces[i][a] = 0;

    vector<long> boxOfParticle(nbParticles);
    vector< forward_list<long> > particlesOfBox(nbSubBoxes);

    // Classify the particles in sub-boxes
    makeSubBoxes(positions, boxOfParticle, particlesOfBox, nbParticles,
                 subBoxLength, nbSubBoxesLine);

    // For every particle, for every neighboring boxes of its box,
    // for every particle in this later box,
    // compute the force between the two particles.
    for (long i = 0 ; i < nbParticles ; i++) {
        for (long box : nbrLists[boxOfParticle[i]]) {
            for (auto it_j = particlesOfBox[box].begin() ;
                 it_j != particlesOfBox[box].end() && (*it_j) > i ;
                 it_j++) {
                double dr[DIM];
                double distsq = 0.;
                for (int a = 0 ; a < DIM ; a++) {
                    dr[a] = periodicBC(positions[i][a] - positions[*it_j][a],
                                       boxLength);
                    distsq += dr[a] * dr[a];
                }
                if (distsq < 1. && distsq > 0.) {
                    for (int a = 0 ; a < DIM ; a++) {
                        double t = (1. / sqrt(distsq) - 1.) * dr[a];
                        forces[i][a] += t;
                        forces[*it_j][a] -= t;  // f_ij = -f_ji
                    }
                }
            }
        }
    }
}

// Classify the particles in boxes of size approximately 1.
void makeSubBoxes(const vector< array<double,DIM> > &positions,
                  vector<long> &boxOfParticle,
                  vector< forward_list<long> > &particlesOfBox,
                  const long nbParticles, const double subBoxLength,
                  const long nbSubBoxesLine) {
    for (long i=0 ; i < nbParticles ; i++) {
        long box = 0;
        for (int a = 0 ; a < DIM ; a++) {
            box += mypow(nbSubBoxesLine,a) * \
                   ((long) floor(positions[i][a]/subBoxLength +
                                 (double) nbSubBoxesLine/2.));
        }
        boxOfParticle[i] = box;
        particlesOfBox[box].push_front(i);
    }
}

// Initialize parameters needed for the computation of correlations.
void initCorrel(vector<long long> &countInDivTNT,
                vector<long long> &countInDivTT,
                long &nbDivLine, double &divLength, const double boxSize) {
    nbDivLine = (long) floor(boxSize / divLength);
    divLength = boxSize / nbDivLine;
    const long nbDiv = mypow(nbDivLine, DIM);
    countInDivTNT.assign(nbDiv, 0);
    countInDivTT.assign(nbDiv, 0);
}

// Compute correlations between particules.
// Each tracer-non tracer (resp. tracer-tracer) occurence is added to the
// correponding division of countInDivTNT (resp. countInDivTT)
void calcCorrel(const vector< array<double,DIM> > &pos,
                vector<long long> &countInDivTNT,
                vector<long long> &countInDivTT,
                const long nbParticles, const long nbTracers,
                const double boxSize, const long nbDivLine) {
    for (long i = 0 ; i < nbTracers ; i++) {
        // Loop on tracer-non tracer pairs
        for (long j = nbTracers ; j < nbParticles ; j++) {
            long div = 0;
            for (int a=0 ; a<DIM ; a++) {
                double dx = periodicBC(pos[j][a] - pos[i][a], boxSize);
                long u = (long) floor((dx + boxSize/2.) * nbDivLine / boxSize);
                div += mypow(nbDivLine,a) * u;
            }
            countInDivTNT[div]++;
        }
        // Loop on tracer-tracer pairs
        for (long j = i+1 ; j < nbTracers ; j++) {
            long div = 0;
            for (int a = 0 ; a < DIM ; a++) {
                double dx = periodicBC(pos[j][a] - pos[i][a], boxSize);
                long u = (long) floor((dx + boxSize/2.) * nbDivLine / boxSize);
                div += mypow(nbDivLine,a) * u;
            }
            countInDivTT[div]++;
        }
    }
}

// Create the list of neighbors of a given box.
forward_list<long> makeNeighborsList(const long L, long i) {
    forward_list<long> nbrList = {i};
    long ia, ia1, ia2;

    for (int a = DIM-1 ; a >= 0 ; a--) {
        forward_list<long> tmp = nbrList;
        const long p = mypow(L, a);
        ia = i / p;
        i = i - p * ia;
        ia1 = (ia - 1 + L) % L - ia;
        ia2 = (ia + 1) % L - ia;

        for (long x : tmp) {
            nbrList.push_front(x + p*ia1);
            nbrList.push_front(x + p*ia2);
        }
    }
    return nbrList;
}

// Create the vector of the list of neighbors of each box.
vector< forward_list<long> > makeNeighborsLists(const long L) {
    const long N = mypow(L, DIM);
    vector< forward_list<long> > nbrLists(N);
    for (long i = 0 ; i < N ; i++)
        nbrLists[i] = makeNeighborsList(L, i);
    return nbrLists;
}


// Write one line in the output file for mobility
void exportLineMobility(const std::vector< std::array<double, DIM> > &pos,
                        const std::vector< std::array<double, DIM> > &forces,
                        const std::vector< std::array<double, DIM> > &windings,
                        const std::vector<double> initPosX,
                        const long nbParticles, const long nbTracers,
                        const double boxLength, const double time0,
                        ostream &stream) {
    // Force on the tracers from the other particles
    // and average movement of the tracers
    double fIntPar1 = 0., fIntPar2 = 0.,  displ = 0.;
    for (long i = 0 ; i < nbTracers ; i++) {
        fIntPar1 += forces[i][0];
        displ += windings[i][0] * boxLength + pos[i][0] - initPosX[i];
    }
    for (long i = nbTracers ; i < nbParticles ; i++) {
        fIntPar2 += forces[i][0];
    }

    stream << time0 << "\t" << fIntPar1 / nbTracers << "\t" << displ / nbTracers
        << "\t" << fIntPar2 / (nbParticles - nbTracers) << "\t"
        << (fIntPar1 + fIntPar2) / nbParticles << "\n";
}

// Export correlations: "dx1 dx2 ... countTNT countTT".
void exportCorrel(const vector<long long> &countInDivTNT,
                 const vector<long long> &countInDivTT, const long nbDivLine,
                 const double divLength, const double boxLength,
                 ostream &stream) {
    const long nbDiv = mypow(nbDivLine, DIM);
    long long sumTNT = 0;

    for (long long u : countInDivTNT) {
        sumTNT += u;
    }

    for (long k=0 ; k<nbDiv ; k++) {
        for (int a=0 ; a<DIM ; a++) {
            long u = (k / mypow(nbDivLine, a)) % nbDivLine;
            double dx = u * boxLength / nbDivLine - boxLength/2. + \
                        divLength/2.;
            stream << dx << " ";
        }
        stream << countInDivTNT[k]  << " "  << countInDivTT[k] << "\n";
    }
}


// Open files for mobility and correlations
// Use boost::iostreams for gzip compression if needed
int openFiles(const Parameters &p, ostream *&fileMobility,
              ostream *&fileCorrel, string &fullnamePos, mt19937 &rndGen) {
    std::string fullnameMobility, fullnameCorrel;
    genFileNames(p, fullnameMobility, fullnameCorrel, fullnamePos, rndGen);

    if (p.outputMobility) {
        if (openFile(fullnameMobility, fileMobility, p.gzip, p.verbose))
            return 1;
    }
    if (p.outputCorrel) {
        if (openFile(fullnameCorrel, fileCorrel, p.gzip, p.verbose))
            return 1;
    }

    return 0;
}

// Open a file with or without gzip compression
// file becomes either a io::fitering_ostream or a std::ofstream
int openFile(const string &fname, ostream *&file, bool gzip, bool verbose) {
    if (gzip) {
        file = new io::filtering_ostream();
        io::filtering_ostream *f = static_cast<io::filtering_ostream*>(file);
        f->push(io::gzip_compressor());

        io::file_sink fs = io::file_sink(fname);
        if(!fs.is_open()) {
            cerr << "Error: could not open " << fname << endl;
            return 1;
        }
        if (verbose) {
            cout << fname << " successfully opened" << endl;
        }
        f->push(fs);
    } else {
        file = new ofstream(fname);
        ofstream *f = static_cast<ofstream*>(file);
        if(!f->is_open()) {
            cerr << "Error: could not open " << fname << endl;
            return 1;
        }
        if (verbose) {
            cout << fname << " successfully opened" << endl;
        }
    }
    file->precision(DEFAULT_OUT_PRECISION);
    return 0;
}

// Close the files that were previously open
void closeFiles(const Parameters &p, ostream *&fileMobility,
                ostream *&fileCorrel) {
    if (p.outputMobility)
        delete fileMobility;
    if (p.outputCorrel)
        delete fileCorrel;
}

// Generate file names
void genFileNames(const Parameters &p, string &fnameMobility,
                  string &fnameCorrel, string &fnamePos, mt19937 &rndGen) {
    // Defaut file names
    // eg BT2d-20150320-102714-sUfrH-mobility.dat
    ostringstream oss;
    oss << "BT" << DIM << "d-";

    std::time_t t = std::time(nullptr);
    oss << std::put_time(std::localtime(&t), "%Y%m%d-%H%M%S") << "-";

    const char charset[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz";
    uniform_int_distribution<size_t> rndInt(0, sizeof(charset) - 1);
    for (int i = 0 ; i < DEFAULT_NB_RND_CHAR ; i++){
        oss << charset[rndInt(rndGen)];
    }

    fnameMobility = p.path + "/";
    fnameCorrel = p.path + "/";
    fnamePos = p.path + "/";
    if (p.fnameMobility.empty()) {
        fnameMobility += oss.str() + "-mobility.dat";
    } else {
        fnameMobility += p.fnameMobility;
    }
    if (p.fnameCorrel.empty()) {
        fnameCorrel += oss.str() + "-correl.dat";
    } else {
        fnameCorrel += p.fnameCorrel;
    }
    if (p.fnamePos.empty()) {
        fnamePos += oss.str() + "-pos";
    } else {
        fnameCorrel += p.fnameCorrel;
    }
    if (p.gzip) {
        fnameMobility += ".gz";
        fnameCorrel += ".gz";
    }
}

// Write file headers
void writeHeaders(const Parameters &p, ostream &fileMobility,
                  ostream &fileCorrel) {
    ostringstream oss;
    oss << "# BrownianTracers" << DIM << "d: "
        << "nbParticles=" << p.nbParticles << ", nbTracers=" << p.nbTracers
        << ", density=" << p.density << ", temperature=" << p.temperature
        << ", force=" << p.force << ", timestep=" << p.timestep
        << ", nbIters=" << p.nbIters << ", nbItersTh=" << p.nbItersTh;

    if (p.outputMobility) {
        fileMobility << oss.str() << "\n";
        fileMobility << "# time forceIntOnTracers displacementTracers"
            << " forceIntOnNonTracers forceIntAll" << endl;
    }
    if (p.outputCorrel) {
        fileCorrel << oss.str() << ", resolution=" << p.divLength
            << ", skip=" << p.timeCalcCorrel << "\n";
        fileCorrel << "# ";
        for(int a = 0 ; a < DIM ; a++){
            fileCorrel << "x" << a + 1 << " ";
        }
        fileCorrel << "correlTNT correlTT" << endl;
    }
}

// Write header for position file
void writeHeaderPos(const Parameters &p, ostream &filePos, const long k) {
    filePos << "# BrownianTracers" << DIM << "d: "
        << "nbParticles=" << p.nbParticles << ", nbTracers=" << p.nbTracers
        << ", density=" << p.density << ", temperature=" << p.temperature
        << ", force=" << p.force << ", timestep=" << p.timestep
        << ", nbIters=" << p.nbIters << ", nbItersTh=" << p.nbItersTh
        << ", time=" << k * p.outputPosIters << "\n";
    filePos << "#";
    for(int a = 0 ; a < DIM ; a++){
        filePos << " x" << a + 1;
    }
    filePos << endl;
}

// Export the positions of the particles.
void exportPos(const Parameters &p, const vector< array<double,DIM> > &pos,
               const string &fname, const long k) {
    ostringstream oss;
    oss << fname << "-" << setfill('0') << setw(DEFAULT_NB_DIGITS) << k
        << ".dat";
    if (p.gzip) {
        oss << ".gz";
    }

    ostream *filePos;
    if (openFile(oss.str(), filePos, p.gzip, false))
        return;

    writeHeaderPos(p, *filePos, k);
    for (long i = 0 ; i < (long) pos.size() ; i++) {
        for (int a = 0 ; a < DIM ; a++) {
            *filePos << pos[i][a] << "\t";
        }
        *filePos << endl;
    }

    delete filePos;
}
