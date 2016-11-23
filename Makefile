# Copyright ESPCI Paris (2016)
# Contributors: Alexis Poncet, Vincent Démery

# alexis.poncet@ens.fr

# This software is a computer program whose purpose is to simulate the
# motion of interacting Brownian particle, some of them being driven by
# an external force.

# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use, 
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info". 

# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability. 

# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or 
# data to be ensured and,  more generally, to use and operate it in the 
# same conditions as regards security. 

# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.

CC=g++
CFLAGS=-W -Wall -ansi -pedantic -std=c++11 -O3
LDFLAGS=-lboost_program_options -lboost_iostreams
LDFLAGS_VISU2D=-lsfml-graphics -lsfml-window -lsfml-system -pthread
EXEC=brownianTracers-2d brownianTracers-3d brownianTracers-2d-visu

all: $(EXEC)

brownianTracers-2d: main-2d.o simul-2d.o parseArguments.o
		$(CC) -o $@ $^ $(LDFLAGS)

brownianTracers-3d: main-3d.o simul-3d.o parseArguments.o
		$(CC) -o $@ $^ $(LDFLAGS)

brownianTracers-2d-visu: main-2d-visu.o simul-2d-visu.o visu2d.o parseArguments-visu.o
		$(CC) -o $@ $^ $(LDFLAGS) $(LDFLAGS_VISU2D)


main-2d.o: main.cpp simul.h utils.h parseArguments.h
		$(CC) -o $@ -c $< -DDIM=2 $(CFLAGS)

main-3d.o: main.cpp simul.h utils.h parseArguments.h
		$(CC) -o $@ -c $< -DDIM=3 $(CFLAGS)

main-2d-visu.o: main.cpp simul.h utils.h parseArguments.h
		$(CC) -o $@ -c $< -DDIM=2 -DVISU2D $(CFLAGS)


parseArguments.o: parseArguments.cpp utils.h parseArguments.h
		$(CC) -o $@ -c $< $(CFLAGS)

parseArguments-visu.o: parseArguments.cpp utils.h parseArguments.h
		$(CC) -o $@ -c $< -DVISU2D $(CFLAGS)


simul-2d.o: simul.cpp utils.h simul.h
		$(CC) -o $@ -c $< -DDIM=2 $(CFLAGS)

simul-3d.o: simul.cpp utils.h simul.h
		$(CC) -o $@ -c $< -DDIM=3 $(CFLAGS)

simul-2d-visu.o: simul.cpp utils.h simul.h
		$(CC) -o $@ -c $< -DDIM=2 -DVISU2D $(CFLAGS)


visu2d.o: visu2d.cpp utils.h
		$(CC) -o $@ -c $< -DDIM=2 $(CFLAGS)


.PHONY: clean mrproper

clean:
		rm -f *.o

mrproper: clean
		rm -rf $(EXEC)
