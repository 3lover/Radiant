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

// the Ray struct, used for casting
struct Ray {
    double x = 0;
    double y = 0;
    double vx = 0;
    double vy = 0;
};

// line segments, which reflect and/or absorb light
class LineSegment {
public:
    vector<double> a = { 0, 0 };
    vector<double> b = { 0, 0 };
    vector<int> color = { 255, 255, 255 };
    double reflectivity = 1;
    LineSegment(vector<double> a, vector<double> b, vector<int> color, double reflectivity) {
        this->a = a;
        this->b = b;
        this->color = color;
        this->reflectivity = reflectivity;
    }
    // renders the line segment
    void render(RenderWindow& window) {
        drawLine(window, a[0], a[1], b[0], b[1], 5, color);
    }
    // returns true if a given ray intersects this, or false otherwise
    vector<double> intersects(Ray r) {
        // we don't care about the case of infinite intersection, we treat these lasers/segments as infinitely thin
        if (atan2(b[1] - a[1], b[0] - a[0]) == atan2(r.vy, r.vx)) return {};
        // the difference between the slopes
        double slopeDiff = (r.vy / r.vx) - ((b[1] - a[1]) / (b[0] - a[0]));
        // align them at the same x point as this segment's A point, so we can work out y distance from it
        double rayYAtPoint = r.y + (r.vy / r.vx) * (a[0] - r.x);
        // calculate the distance in y values between the two points, and work backwards to get the intersection point
        double yDiff = a[1] - rayYAtPoint;
        double Xi = (yDiff / slopeDiff) + a[0];
        double Yi = r.y + (r.vy / r.vx) * (Xi - r.x);

        // if the intersection x or y are outside the bounds of our line segment, we don't intersect
        vector<double> bounds = { min(a[0], b[0]) - 1, min(a[1], b[1]) - 1, max(a[0], b[0]) + 1, max(a[1], b[1]) + 1 };
        if (bounds[0] > Xi || bounds[1] > Yi || bounds[2] < Xi || bounds[3] < Yi) {
            return {};
        }

        // if the intersection is in the wrong direction from our ray, we don't intersect
        if ((Xi - r.x <= 0 && r.vx > 0) || (Xi - r.x >= 0 && r.vx < 0) || (Yi - r.y <= 0 && r.vy > 0) || (Yi - r.y >= 0 && r.vy < 0)) {
            return {};
        }

        return {Xi, Yi};
    }

    // gets the direction a ray should bounce off this
    double getBounce(Ray r) {
        // get the normal of this line segment
        vector<double> normal = { -(b[1] - a[1]), (b[0] - a[0]) };
        // split the ray's directional vector into perpendicular components to this normal
        double uScalar = (r.vy * normal[1] + r.vx * normal[0]) / (normal[1] * normal[1] + normal[0] * normal[0]);
        vector<double> u = {uScalar * normal[0], uScalar * normal[1]};
        vector<double> w = {r.vx - u[0], r.vy - u[1]};
        vector<double> vPrime = {w[0] - u[0], w[1] - u[1]};
        double newAngle = atan2(vPrime[1], vPrime[0]);

        return newAngle;
    }
};

// a projector shoots a laser and bounces it
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

        drawLine(
            window,
            pos[0] + radius,
            pos[1] + radius,
            pos[0] + cos(rotation) * radius + radius,
            pos[1] + sin(rotation) * radius + radius,
            3,
            {255, 0, 0}
        );
    }
    // calculates the bounces of the laser, and alerts hit objects
    void update(vector<LineSegment> segments) {
        // calculate the bounces of the laser from various surfaces
        //bounces = { {150, 150}, { 200, 200 }, {400, 100}, {400, 700 }, {500, 900 } };
        vector<double> lastPos = { pos[0] + radius, pos[1] + radius };
        double lastDir = rotation;
        bounces = {lastPos};
        int lastBounce = -1;
        for (int m = 0; m < 100; m++) {
            int closestIndex = -1;
            double closestDist = -1;
            vector<double> closestIntersection = {};
            bool hit = false;
            for (int i = 0; i < segments.size(); i++) {
                if (lastBounce == i) continue;
                vector<double> intersection = segments[i].intersects({ lastPos[0], lastPos[1], cos(lastDir), sin(lastDir) });
                if (intersection.size() == 0) continue;
                // if we do intersect, check if we are closest
                double distX = intersection[0] - lastPos[0];
                double distY = intersection[1] - lastPos[1];
                double intersectDist = sqrt(pow(distX, 2) + pow(distY, 2));
                if (closestDist < 0 || intersectDist < closestDist) {
                    closestDist = intersectDist;
                    closestIndex = i;
                    closestIntersection = intersection;
                    hit = true;
                }
            }
            if (!hit) {
                // find where we land among the boundary
                break;
            }
            lastPos = closestIntersection;
            lastDir = segments[closestIndex].getBounce({lastPos[0], lastPos[1], cos(lastDir), sin(lastDir)});
            bounces.push_back(lastPos);
            lastBounce = closestIndex;
            //cout << "Bounced " << closestIndex << " with strand with colors " << segments[closestIndex].color[0] << ", " << segments[closestIndex].color[1] << ", " << segments[closestIndex].color[2] << endl;
            if (segments[closestIndex].reflectivity == 0) break;
        }

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
    projectors.push_back(Projector({BOUNDX/2 - 25, BOUNDY / 2 - 25 }, {0,255,0}, 0, false, false));

    vector<LineSegment> segments = {};
    for (int i = 0; i < 8; i++) {
        segments.push_back(LineSegment({BOUNDX * randDouble(), BOUNDY * randDouble()}, { BOUNDX * randDouble(), BOUNDY * randDouble() }, {rand() % 255,rand() % 255,rand() % 255 }, 1));
    }

    segments.push_back(LineSegment({ 0, 0 }, { 0, BOUNDY }, { 255,255,255 }, 0));
    segments.push_back(LineSegment({ 0, BOUNDY }, { BOUNDX, BOUNDY }, { 255,255,255 }, 0));
    segments.push_back(LineSegment({ BOUNDX, BOUNDY }, { BOUNDX, 0 }, { 255,255,255 }, 0));
    segments.push_back(LineSegment({ BOUNDX, 0 }, { 0, 0 }, { 255,255,255 }, 0));

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
                    for (int i = 0; i < projectors.size(); i++) {
                        projectors[i].rotation = atan2(mouseLoc[1] - projectors[i].pos[1] - projectors[i].radius, mouseLoc[0] - projectors[i].pos[0] - projectors[i].radius);
                    }
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

        // update projector's bounces and stuff
        for (int i = 0; i < projectors.size(); i++) {
            projectors[i].update(segments);
        }

        // clear the window with black color
        window.clear(sf::Color::Black);

        // render every ball
        for (int i = 0; i < projectors.size(); i++) {
            projectors[i].render(window);
        }

        // render every segment
        for (int i = 0; i < segments.size(); i++) {
            segments[i].render(window);
        }

        // end the current frame
        window.display();
    }

    return 0;
}