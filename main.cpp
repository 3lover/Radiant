#include <SFML/Graphics.hpp>
#include <iostream>
#include <vector>
#include <random>
#include <chrono>

using namespace std;
using namespace sf;

const int BOUNDX = 1600;
const int BOUNDY = 900;
const double PI = 3.14159265358979323846;

// generates a random number in [0, 1)
double randDouble() {
    return (((double)rand() - 1) / RAND_MAX);
}

// draws a line between 2 points
void drawLine(RenderWindow& window, double ax, double ay, double bx, double by, double thickness, vector<int> color) {
    Vector2f point1(ax, ay);
    Vector2f point2(bx, by);

    // Calculate length and angle
    double dx = point2.x - point1.x;
    double dy = point2.y - point1.y;
    double length = sqrt(dx*dx + dy*dy);
    double angle = atan2(dy, dx) * 180 / PI;

    RectangleShape line(Vector2f(length, thickness));
    line.setOrigin(0, thickness / 2);
    line.setPosition(point1);
    line.setRotation(angle);
    line.setFillColor(Color(color[0], color[1], color[2]));

    window.draw(line);
}

class Projector {
public:
    CircleShape shape = CircleShape(50);
    vector<double> pos = { 0, 0 };
    vector<int> borderColor = { 255, 255, 255 };
    vector<int> laserColor = { 255, 255, 255 };
    double rotation = 0;
    double radius = 50;
    bool requireCharge = false;
    bool locked = false;

    vector<vector<double>> bounces = {};
    // full constructor in case
    Projector(vector<double> pos, vector<int> borderColor, vector<int> laserColor, double rotation, double radius, bool requireCharge, bool locked) {
        shape.setRadius(radius);
        this->pos = { pos[0], pos[1] };
        this->borderColor = borderColor;
        this->laserColor = laserColor;
        this->rotation = rotation;
        this->radius = radius;
        this->requireCharge = requireCharge;
        this->locked = locked;
    }
    // practically, we only need these from our constructor
    Projector(vector<double> pos, vector<int> laserColor, double rotation, bool requireCharge, bool locked) {
        shape.setRadius(radius);
        this->pos = { pos[0], pos[1] };
        this->borderColor = { 255, 255, 255 };
        this->laserColor = laserColor;
        this->rotation = rotation;
        this->radius = 50;
        this->requireCharge = requireCharge;
        this->locked = locked;
    }
    // renders the laser and all of its bounces, called in the update loop at 60Hz
    void render(RenderWindow& window) {
        // draw the shape first
        window.draw(shape);

        // draw the laser based on the segments calculated between bounces
        for (int i = 0; i < bounces.size() - 1; i++) {
            drawLine(window, bounces[i][0], bounces[i][1], bounces[i + 1][0], bounces[i + 1][1], 5, laserColor);
        }
    }
    // calculates the bounces of the laser, and alerts hit objects
    void update() {
        // calculate the bounces of the laser from various surfaces
        bounces = { {150, 150}, { 200, 200 }, {400, 100}, {400, 700 }, {500, 900 } };

        // update color and position of the render object
        shape.setFillColor(Color(laserColor[0], laserColor[1], laserColor[2]));
        shape.setPosition(pos[0], pos[1]);
        shape.setOutlineThickness(radius/5);
        shape.setOutlineColor(Color(borderColor[0], borderColor[1], borderColor[2]));
    }
    // if hit, check the history to see if we should change our color and propogate the hit
    void hit(vector<Projector*> hits) {
        // if we find ourselves anywhere on the list, stop the propogation and assign all lasers on the list the same color, since it is a loop
        bool infiniteLoop = false;
        for (Projector* hit : hits) if (hit == this) infiniteLoop = true;
        for (Projector* hit : hits) {
            Projector &proj = *hit;
            proj.laserColor = { 255, 0, 0 }; // set to red for testing
        }
    }
};

int main()
{
    // create the window
    sf::RenderWindow window(sf::VideoMode(BOUNDX, BOUNDY), "Radiant");

    vector<double> mouseLoc = {0, 0};
    auto lastTick = chrono::steady_clock::now();

    vector<Projector> projectors = {};
    projectors.push_back(Projector({100, 100}, {0,255,0}, 0, false, false));

    // run the program as long as the window is open
    while (window.isOpen())
    {
        // check all the window's events that were triggered since the last iteration of the loop
        sf::Event event;
        while (window.pollEvent(event)) {
            // check the type of the event
            switch (event.type) {
                // window closed
                case sf::Event::Closed: {
                    window.close();
                    break;
                }

                // key pressed
                case sf::Event::KeyPressed: {

                    break;
                }

                // mouse clicked
                case sf::Event::MouseButtonPressed: {
                    if (event.mouseButton.button == sf::Mouse::Left)
                    {
                        
                    }
                    break;
                }

                case sf::Event::MouseMoved: {
                    mouseLoc = { (double)event.mouseMove.x, (double)event.mouseMove.y };
                    break;
                }

                default: {
                    break;
                }
            }
        }

        // run the update loop every 17ms, for an update loop at about 58.8Hz
        double dt = chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - lastTick).count();
        if (dt < 17) continue;
        lastTick = chrono::steady_clock::now();

        // update balls positions and velocities
        for (int i = 0; i < projectors.size(); i++) {
            projectors[i].update();
        }

        // clear the window with black color
        window.clear(sf::Color::Black);

        // render every ball
        for (int i = 0; i < projectors.size(); i++) {
            projectors[i].render(window);
        }

        // end the current frame
        window.display();
    }

    return 0;
}