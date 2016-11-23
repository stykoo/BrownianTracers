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
 * visu2d.cpp
 *
 * Author: Alexis Poncet
 * Email: alexis.poncet@ens.fr
 * Version: 1.0
 *
 * Contain the functions to visualize the particles in dimension 2.
 * Use SFML library.
 */

#include <iostream>
#include <cstdlib>
#include "utils.h"
#include "visu2d.h"

using namespace std;

// Draw a given static state of the system.
void drawStateStatic(const vector< array<double, DIM> > &positions,
                     const long nbParticles, const long nbTracers,
                     const double boxLength){
    if(DIM != 2){
        cerr << "Error: visualization in dimension " << DIM
            << " is not available." << endl;
        return;
    }

    sf::RenderWindow window(sf::VideoMode(WINDOW_SIZE, WINDOW_SIZE),
                            "Tracers");
    sf::CircleShape circle1((int) (WINDOW_SIZE/boxLength));
    sf::CircleShape circle2((int) (WINDOW_SIZE/boxLength));
    circle1.setFillColor(sf::Color::Transparent);
    circle2.setFillColor(sf::Color::Transparent);
    circle1.setOutlineThickness(-3);
    circle2.setOutlineThickness(-3);
    circle1.setOutlineColor(sf::Color::Red);
    circle2.setOutlineColor(sf::Color::Blue);

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        window.clear(sf::Color::White);

        for (long i = nbTracers ; i < nbParticles ; i++)
            drawOneParticle(positions[i], boxLength, window, circle2);
        for (long i = 0 ; i < nbTracers ; i++)
            drawOneParticle(positions[i], boxLength, window, circle1);

        window.display();
    }
}

// Draw the particules on the screen and update at a given FPS rate.
void threadVisu(const DataThreadVisu &data, const long nbParticles,
                const long nbTracers, const double boxLength,
                const double dt) {
    sf::RenderWindow window;
    sf::CircleShape circle1, circle2;
    sf::Font font;
    sf::Text text;

    initVisu(window, circle1, circle2, font, text, boxLength);

    window.setFramerateLimit(FPS);

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        updateVisu(window, circle1, circle2, text, *data.pos, nbParticles,
                   nbTracers, boxLength, (*data.t) * dt);
    }
}

// Initialize data for visualization.
void initVisu(sf::RenderWindow &window, sf::CircleShape &circle1,
              sf::CircleShape &circle2, sf::Font &font, sf::Text &text,
              const double boxLength, const bool fill) {
    window.create(sf::VideoMode(WINDOW_SIZE, WINDOW_SIZE), "Tracers");
    circle1.setRadius((int) (CIRCLE_RADIUS * WINDOW_SIZE / boxLength));
    circle2.setRadius((int) (CIRCLE_RADIUS * WINDOW_SIZE / boxLength));

    if (fill) {
        circle1.setFillColor(sf::Color::Red);
        circle2.setFillColor(sf::Color::Blue);
    } else {
        circle1.setFillColor(sf::Color::Transparent);
        circle2.setFillColor(sf::Color::Transparent);
        circle1.setOutlineThickness(-3);
        circle2.setOutlineThickness(-3);
        circle1.setOutlineColor(sf::Color::Red);
        circle2.setOutlineColor(sf::Color::Blue);
    }

    if (!font.loadFromFile(FONT_FILE))
        cerr << "Error: could not load font from " << FONT_FILE << endl;

    text.setFont(font);
    text.setString("");
    text.setCharacterSize(30);
    text.setFillColor(sf::Color::Black);
    text.setStyle(sf::Text::Bold);
}

// Update the screen with the new positions of the particles.
// If screenshotNum >= 0, save a screenshot.
void updateVisu(sf::RenderWindow &window, sf::CircleShape &circle1,
                sf::CircleShape &circle2, sf::Text &text,
                const vector< array<double, DIM> > &positions,
                const long nbParticles, const long nbTracers,
                const double boxLength, const double time,
                const int screenshotNum) {
    // With or without screenshot
    if(screenshotNum < 0) {
        window.clear(sf::Color::White);

        for (long i = nbTracers ; i < nbParticles ; i++)
            drawOneParticle(positions[i], boxLength, window, circle2);
        for (long i = 0 ; i < nbTracers ; i++)
            drawOneParticle(positions[i], boxLength, window, circle1);

        text.setString(to_string(time));
        window.draw(text);
        window.display();
    } else {
        sf::RenderTexture texture;
        texture.create(WINDOW_SIZE, WINDOW_SIZE);
        texture.clear(sf::Color::White);

        for (long i = nbTracers ; i < nbParticles ; i++)
            drawOneParticle(positions[i], boxLength, texture, circle2);
        for (long i = 0 ; i < nbTracers ; i++)
            drawOneParticle(positions[i], boxLength, texture, circle1);

        text.setString(to_string(time));
        texture.draw(text);
        texture.display();

        window.clear();
        window.draw(sf::Sprite(texture.getTexture()));
        window.display();

        sf::Image im = texture.getTexture().copyToImage();
        char buf[100];
        sprintf(buf, "screenshots/screenshot%07d.png", screenshotNum);
        im.saveToFile(buf);
    }
}

// Draw a single particule at a given position.
// Take periodic boundary conditions into account.
void drawOneParticle(const array<double,DIM> &pos, const double boxLength,
                     sf::RenderTarget &target, sf::CircleShape &circle) {
    const float x = WINDOW_SIZE * periodicBCpos(pos[0], boxLength) / boxLength;
    const float y = WINDOW_SIZE * periodicBCpos(pos[1], boxLength) / boxLength;
    const float diam = 2. * circle.getRadius();
    int per1 = 0, per2 = 0;

    if (x < diam)
        per1 = 1;
    else if(x > WINDOW_SIZE - diam)
        per1 = -1;
    if (y < diam)
        per2 = 1;
    else if(y > WINDOW_SIZE - diam)
        per2 = -1;

    for(int i = 0 ; i <= per1*per1 ; i++){
        for(int j = 0 ; j <= per2*per2 ; j++){
            circle.setPosition(x + per1*i*WINDOW_SIZE, y + per2*j*WINDOW_SIZE);
            target.draw(circle);
        }
    }
}
