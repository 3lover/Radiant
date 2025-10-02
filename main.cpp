#include <SFML/Graphics.hpp>
#include <iostream>
#include <vector>
#include <random>
#include <chrono>

using namespace std;
using namespace sf;

const int BOUNDX = 1600;
const int BOUNDY = 900;

double randDouble() {
    return ((double)rand() / RAND_MAX);
}

class Ball {
public:
    CircleShape shape = CircleShape(50);
    vector<double> pos = { 0, 0 };
    vector<int> color = { 255, 255, 255 };
    vector<double> velocity = { 0, 0 };
    double radius = 1;
    Ball(double radius, vector<double> pos, vector<int> color, vector<double> velocity) {
        shape.setRadius(radius);
        this->pos = {pos[0], pos[1]};
        this->color = color;
        this->velocity = velocity;
        this->radius = radius;
    }
    // renders the ball, called in the update loop at 60Hz
    void render(RenderWindow& window) {
        shape.setOutlineThickness(10.f);
        shape.setOutlineColor(sf::Color(255, 255, 255));
        window.draw(shape);
    }
    // moves/colors the balls, reversing their velocity if they hit a wall
    void update() {
        // move the shape in the direction of its velocity
        pos[0] += velocity[0];
        pos[1] += velocity[1];

        // if we are out of bounds/hitting a wall, bounce us
        if ((pos[0] < 0 && velocity[0] < 0) || (pos[0] + radius * 2 > BOUNDX && velocity[0] > 0)) velocity[0] *= -1;
        if ((pos[1] < 0 && velocity[1] < 0) || (pos[1] + radius * 2 > BOUNDY && velocity[1] > 0)) velocity[1] *= -1;

        // update color and position of the render object
        shape.setFillColor(sf::Color(color[0], color[1], color[2]));
        shape.setPosition(pos[0], pos[1]);
    }
};

int main()
{
    // create the window
    sf::RenderWindow window(sf::VideoMode(BOUNDX, BOUNDY), "Radiant");

    vector<Ball> balls = {};
    for (int i = 0; i < 100; i++) {
        balls.push_back(Ball(
            randDouble() * 90 + 10,
            {randDouble() * (BOUNDX - 100), randDouble() * (BOUNDY - 100)},
            {(int)(randDouble() * 255), (int)(randDouble() * 255), (int)(randDouble() * 255)},
            {randDouble() * 30 - 15, randDouble() * 30 - 15}
        ));
    }

    auto lastTick = chrono::steady_clock::now();

    // run the program as long as the window is open
    while (window.isOpen())
    {
        // check all the window's events that were triggered since the last iteration of the loop
        sf::Event event;
        while (window.pollEvent(event))
        {
            // "close requested" event: we close the window
            if (event.type == sf::Event::Closed)
                window.close();
        }

        // run the update loop every 17ms, for an update loop at about 58.8Hz
        double dt = chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - lastTick).count();
        if (dt < 17) continue;
        lastTick = chrono::steady_clock::now();

        // update balls positions and velocities
        for (int i = 0; i < balls.size(); i++) {
            balls[i].update();
        }

        // clear the window with black color
        window.clear(sf::Color::Black);

        // render every ball
        for (int i = 0; i < balls.size(); i++) {
            balls[i].render(window);
        }

        // end the current frame
        window.display();
    }

    return 0;
}