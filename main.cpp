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
 * main.cpp
 *
 * Author: Alexis Poncet
 * Email: alexis.poncet@ens.fr
 * Version: 1.0
 *
 * Simulation of interacting Brownian particles,
 * some of them are driven by a force.
 */

#include <iostream>
#include "utils.h"
#include "parseArguments.h"
#include "simul.h"

using namespace std;

int main(int argc, char**argv) {
    Parameters params;
    int status = parseAndCheckArguments(argc, argv, params);
    if (status == 1)
        return 0;
    else if (status == 2)
        return 1;

    if (params.verbose) {
        cout << "--- This is " << argv[0] << " (DIM = " << DIM
            << "), compiled on " << __DATE__ << " at " << __TIME__ << " ---"
            << endl;
        printParameters(params);
    }

    status = runSimulation(params);

    return status;
}
