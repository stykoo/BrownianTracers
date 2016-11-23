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
 * Define constants and data structure associated to visualization.
 * Define functions of visu2d.cpp.
 */

#ifndef BROWNIANTRACERS_VISU2D_H_
#define BROWNIANTRACERS_VISU2D_H_

#include <vector>
#include <array>
#include <SFML/Graphics.hpp>

#define WINDOW_SIZE 800
#define CIRCLE_RADIUS 0.5
#define FPS 24
#define FONT_FILE "font.ttf"

// Data structure to be passed to the visualization thread.
struct DataThreadVisu {
    std::vector< std::array<double, DIM> > *pos;  // Positions
    int *t;  // Current iteration of the simulation
};

void drawStateStatic(const std::vector< std::array<double, DIM> > &positions,
                     const long nbParticles, const long nbTracers,
                     const double boxLength);

void threadVisu(const DataThreadVisu &data, const long nbParticles,
                const long nbTracers, const double boxLength,
                const double dt);

void initVisu(sf::RenderWindow &window, sf::CircleShape &circle1,
              sf::CircleShape &circle2, sf::Font &font, sf::Text &text,
              const double boxLength, const bool fill = false);

void updateVisu(sf::RenderWindow &window, sf::CircleShape &circle1,
                sf::CircleShape &circle2, sf::Text &text,
                const std::vector< std::array<double, DIM> > &positions,
                const long nbParticles, const long nbTracers,
                const double boxLength, const double time,
                const int screenshotNum = -1);

void drawOneParticle(const std::array<double,DIM> &pos, const double boxLength,
                     sf::RenderTarget &target, sf::CircleShape &circle);

#endif  // BROWNIANTRACERS_VISU2D_H_
