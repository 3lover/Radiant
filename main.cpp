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

// takes in 2 colors and merges them, then normalizes them
vector<int> mergeColors(vector<int> c1, vector<int> c2) {
    vector<double> merged = { (double)(c1[0] + c2[0]), (double)(c1[1] + c2[1]), (double)(c1[2] + c2[2]) };
    double highest = max(max(merged[0], merged[1]), merged[2]);
    if (highest == 0) return { 0, 0, 0 };
    return { (int)(merged[0] / highest * 255), (int)(merged[1] / highest * 255) , (int)(merged[2] / highest * 255) };
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
        if (atan2(b[1] - a[1], b[0] - a[0]) == atan2(r.vy, r.vx)) {
            return {};
        }

        double Xi, Yi;

        // vertical slope catching
        if (abs(r.vx) < 0.001) {
            Xi = r.x;
            // if this is a vertical line, return false, we don't hit that with a vertical laser
            if (b[0] - a[0] == 0) return {};
            double segSlope = (b[1] - a[1]) / (b[0] - a[0]);
            double intercept = a[1] - segSlope * a[0];
            Yi = segSlope * Xi + intercept;
        }
        // horizontal slope catching
        else if (abs(r.vy) < 0.001) {
            Yi = r.y;
            // if this is a horizontal line, return false, we don't hit that with a vertical laser
            if (b[1] - a[1] == 0) return {};
            double segSlope = (b[0] - a[0]) / (b[1] - a[1]);
            double intercept = a[0] - segSlope * a[1];
            Xi = segSlope * Yi + intercept;
        }
        else {
            // the difference between the slopes
            double slopeDiff = (r.vy / r.vx) - ((b[1] - a[1]) / (b[0] - a[0]));
            // align them at the same x point as this segment's A point, so we can work out y distance from it
            double rayYAtPoint = r.y + (r.vy / r.vx) * (a[0] - r.x);
            // calculate the distance in y values between the two points, and work backwards to get the intersection point
            double yDiff = a[1] - rayYAtPoint;
            Xi = (yDiff / slopeDiff) + a[0];
            Yi = r.y + (r.vy / r.vx) * (Xi - r.x);
        }

        // if the intersection x or y are outside the bounds of our line segment, we don't intersect
        vector<double> bounds = { min(a[0], b[0]) - 1, min(a[1], b[1]) - 1, max(a[0], b[0]) + 1, max(a[1], b[1]) + 1 };
        if (bounds[0] > Xi || bounds[1] > Yi || bounds[2] < Xi || bounds[3] < Yi) {
            return {};
        }

        // if the intersection is in the wrong direction from our ray, we don't intersect
        if ((Xi - r.x < 0 && r.vx > 0) || (Xi - r.x > 0 && r.vx < 0) || (Yi - r.y < 0 && r.vy > 0) || (Yi - r.y > 0 && r.vy < 0)) {
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

// as we go through and calculate loops, this helps us track input and output colors
struct InputColor {
    int colorId;
    bool resolved;
    vector<int> value;
};
vector<vector<int>> connectedLoops = {};

// a projector shoots a laser and bounces it
class Projector {
public:
    CircleShape shape = CircleShape(50);
    vector<double> pos = { 0, 0 };
    vector<int> borderColor = { 255, 255, 255 };
    vector<int> laserColor = { 0, 0, 0 };
    vector<int> producedColor = { 255, 255, 255 };
    vector<int> loopedColor = { -1, -1, -1 };
    vector<vector<Projector*>> infiniteLoops = {};
    double rotation = 0;
    double radius = 25;
    double borderSize = 0.2;
    bool requireCharge = false;
    bool locked = false;
    vector<Projector*> hitting = {};
    vector<InputColor*> inputs = {};
    InputColor* output = new InputColor{ -1, false, {} };
    int colorId = 0;
    bool selected = false;
    double rotating = 0;
    vector<double> beams = {0};
    bool loopChecked = false;

    vector<vector<vector<double>>> bounces = {};
    // full constructor in case
    Projector(vector<double> pos, vector<int> borderColor, vector<int> laserColor, double rotation, double radius, bool requireCharge, bool locked, vector<int> producedColor, vector<double> beams) {
        shape.setRadius(radius);
        this->pos = { pos[0], pos[1] };
        this->borderColor = borderColor;
        this->laserColor = laserColor;
        this->producedColor = producedColor;
        this->rotation = rotation;
        this->radius = radius;
        this->requireCharge = requireCharge;
        this->locked = locked;
        this->beams = beams;
    }
    // practically, we only need these from our constructor
    Projector(vector<double> pos, vector<int> producedColor, double rotation, bool requireCharge, bool locked, vector<double> beams) {
        shape.setRadius(radius);
        this->pos = { pos[0], pos[1] };
        this->borderColor = { 255, 255, 255 };
        this->producedColor = producedColor;
        if (!requireCharge) this->laserColor = producedColor;
        this->rotation = rotation;
        this->radius = 25;
        this->requireCharge = requireCharge;
        this->locked = locked;
        this->beams = beams;
    }
    // renders the laser and all of its bounces, called in the update loop at 60Hz
    void render(RenderWindow& window) {
        // draw the shape first
        window.draw(shape);

        // draw the laser based on the segments calculated between bounces
        for (int i = 0; i < bounces.size(); i++) for (int j = 0; j < bounces[i].size() - 1; j++) {
            drawLine(window, bounces[i][j][0], bounces[i][j][1], bounces[i][j + 1][0], bounces[i][j + 1][1], 5, laserColor);
        }

        drawLine(
            window,
            pos[0] + radius,
            pos[1] + radius,
            pos[0] + cos(rotation) * radius + radius,
            pos[1] + sin(rotation) * radius + radius,
            3,
            { 255 - laserColor[0], 255 - laserColor[0], 255 - laserColor[0] }
        );
    }

    // returns true if a point falls inside this circle
    bool inside(double x, double y) {
        return pow(x - pos[0] - radius, 2) + pow(y - pos[1] - radius, 2) <= pow(radius, 2);
    }

    // checks if a line segment intersects this shape, and sends each set of intersections if so
    vector<double> intersects(vector<double> A, vector<double> B) {
        // compute the euclidean distance between A and B
        double LAB = sqrt(pow((B[0] - A[0]), 2) + pow((B[1] - A[1]), 2));

        // compute the direction vector D from A to B
        double Dx = (B[0] - A[0]) / LAB;
        double Dy = (B[1] - A[1]) / LAB;

        // the equation of the line AB is x = Dx*t + Ax, y = Dy*t + Ay with 0 <= t <= LAB.

        // compute the distance between the points A and E, where
        // E is the point of AB closest the circle center (Cx, Cy)
        double t = Dx * (pos[0] + radius - A[0]) + Dy * (pos[1] + radius - A[1]);

        // compute the coordinates of the point E
        double Ex = t * Dx + A[0];
        double Ey = t * Dy + A[1];

        // compute the euclidean distance between E and C
        double LEC = sqrt(pow((Ex - pos[0] - radius), 2) + pow((Ey - pos[1] - radius), 2));
        vector<double> bounds = { min(A[0], B[0]) - 1, min(A[1], B[1]) - 1, max(A[0], B[0]) + 1, max(A[1], B[1]) + 1 };

        // test if the line intersects the circle
        if (LEC < radius * (1+ borderSize)) {
            // compute distance from t to circle intersection point
            double dt = sqrt(pow(radius * (1 + borderSize), 2) - pow(LEC, 2));

            // compute first intersection point
            double Fx = (t - dt) * Dx + A[0];
            double Fy = (t - dt) * Dy + A[1];

            // compute second intersection point
            double Gx = (t + dt) * Dx + A[0];
            double Gy = (t + dt) * Dy + A[1];

            vector<double> intersections = {};

            // check if both fall within our line segment
            if (bounds[0] <= Fx && bounds[1] <= Fy && bounds[2] >= Fx && bounds[3] >= Fy) {
                intersections.push_back(Fx);
                intersections.push_back(Fy);
            }
            if (bounds[0] <= Gx && bounds[1] <= Gy && bounds[2] >= Gx && bounds[3] >= Gy) {
                intersections.push_back(Gx);
                intersections.push_back(Gy);
            }

            return intersections;
        }
        if (LEC == radius * (1 + borderSize)) {
            // tangent point to circle is E
            if (bounds[0] <= Ex && bounds[1] <= Ey && bounds[2] >= Ex && bounds[3] >= Ey) {
                return { Ex, Ey };
            }
        }
        return {};
    }

    vector<vector<double>> calculateBounces(vector<LineSegment> &segments, vector<Projector> &projectors, double dir) {
        vector<double> lastPos = { pos[0] + radius, pos[1] + radius };
        double lastDir = rotation + dir;
        vector<vector<double>> beamBounces = { lastPos };
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
                hit = true;
                // if we do intersect, check if we are closest
                double distX = intersection[0] - lastPos[0];
                double distY = intersection[1] - lastPos[1];
                double intersectDist = sqrt(pow(distX, 2) + pow(distY, 2));
                if (closestDist < 0 || intersectDist < closestDist) {
                    closestDist = intersectDist;
                    closestIndex = i;
                    closestIntersection = intersection;
                }
            }
            if (!hit) {
                // find where we land among the boundary
                break;
            }
            // we have our next line segment, before we push it check if it intersects with a projector instead
            int projectorId = -1;
            for (int i = 0; i < projectors.size(); i++) {
                if (pos == projectors[i].pos && m == 0) continue;
                vector<double> intersection = projectors[i].intersects(lastPos, closestIntersection);
                if (intersection.size() == 0) continue;
                // if we do intersect another projector, check if it is the closest one and if so send a hit
                projectorId = i;
                if (intersection.size() == 2) {
                    closestIntersection = { intersection[0], intersection[1] };
                    continue;
                }
                double dist1 = pow(intersection[0] - lastPos[0], 2) + pow(intersection[1] - lastPos[1], 2);
                double dist2 = pow(intersection[2] - lastPos[0], 2) + pow(intersection[3] - lastPos[1], 2);
                if (dist1 < dist2) closestIntersection = { intersection[0], intersection[1] };
                else closestIntersection = { intersection[2], intersection[3] };
            }

            // if we intersected a projector, hit them then end our bouncing
            if (projectorId != -1) {
                beamBounces.push_back(closestIntersection);
                projectors[projectorId].hitting.push_back(this);
                break;
            }

            lastPos = closestIntersection;
            lastDir = segments[closestIndex].getBounce({ lastPos[0], lastPos[1], cos(lastDir), sin(lastDir) });
            beamBounces.push_back(lastPos);
            lastBounce = closestIndex;
            if (segments[closestIndex].reflectivity == 0) break;
        }

        return beamBounces;
    }

    // calculates the bounces of the laser, and alerts hit objects
    void update(vector<LineSegment> &segments, vector<Projector> &projectors) {
        // if rotating, rotate us in the right direction
        rotation += rotating * 3.0 / 60.0;

        // calculate the bounces of the laser from various surfaces
        bounces = {};
        for (int i = 0; i < beams.size(); i++) {
            bounces.push_back(calculateBounces(segments, projectors, beams[i]));
        }

        // remove our laser color, we'll resolve for it later anyways
        if (requireCharge) laserColor = { 0, 0, 0 };
        else laserColor = producedColor;
    }

    
    
    // initializes the inputs and outputs related to this laser from the hitting list
    void generateThroughput() {
        // add the output beam to our inputs for every laser hitting us
        for (int i = 0; i < hitting.size(); i++) {
            inputs.push_back(hitting[i]->output);
        }
        // once finished running for every laser, each laser will have the input and output color channels affecting them
    }

    // if the laser has every output defined, go through and collapse the input to a resolved value and propogate it
    bool collapseColorInputs() {
        // if we already have our output, or don't have our inputs settled, return false representing a failed collapse
        for (InputColor* c : inputs) if (!c->resolved) { 
            return false;
        }
        if (output->resolved) {
            return false;
        }
        // if we reach here, every color output we have is resolved so crunch those numbers
        vector<int> avgColor = { producedColor[0], producedColor[1], producedColor[2] };
        for (InputColor* c : inputs) {
            avgColor = { avgColor[0] + c->value[0], avgColor[1] + c->value[1], avgColor[2] + c->value[2] };
            laserColor = mergeColors(avgColor, { 0, 0, 0 });
        }
        output->resolved = true;
        output->value = laserColor;
        return true;
        // once this runs for every projector until they all return false, all non-looped projectors should be solved
    }

    // finds the highest order loop
    void findHighestLoop(vector<int> previous, int safety = 0) {
        if (safety > 99) {
            cout << "Max loop limit reached!";
            return;
        }

        for (int i = previous.size() - 1; i >= 0; i--) {
            // if we have found a completed loop, push our highest-order list of looped projectors for merging later
            if (previous[i] == colorId) {
                vector<int> connected = {};
                for (int j = i; j < previous.size(); j++) connected.push_back(previous[j]);
                for (int j = 0; j < connectedLoops.size(); j++) {
                    if (connectedLoops[j].size() != connected.size()) continue;
                    bool found = false;
                    for (int k = 0; k < connectedLoops[j].size(); k++) {
                        for (int l = 0; l < connected.size(); l++) {
                            if (connected[l] == connectedLoops[j][k]) found = true;
                        }
                        if (!found) break;
                    }
                    if (!found) continue;
                    return;
                }
                connectedLoops.push_back(connected);
                return;
            }
        }

        // otherwise, add ourselves to the loop and move to our next source of input
        previous.push_back(colorId);
        for (int i = 0; i < hitting.size(); i++) {
            if (hitting[i]->output->resolved) continue;
            hitting[i]->findHighestLoop(previous, safety + 1);
        }
    }

    // update the colors of our lasers
    void updateColors() {
        // update color and position of the render object
        shape.setFillColor(Color(laserColor[0], laserColor[1], laserColor[2]));
        shape.setPosition(pos[0], pos[1]);
        shape.setOutlineThickness(radius / 5);
        vector<int> usedBorder = selected ? mergeColors(borderColor, { 0, 255, 0 }) : borderColor;
        shape.setOutlineColor(Color(usedBorder[0], usedBorder[1], usedBorder[2]));
    }
};

// solves a looped color by finding every input to every projector, or using the stored color of the looped items if no inputs exist
vector<int> solveLoop(vector<int> loop, vector<Projector> projectors) {
    vector<int> totalInput = { 0, 0, 0 };
    vector<int> storedColor = { 0, 0, 0 };
    vector<int> defaultColor = { 0, 0, 0 };
    for (int i = 0; i < loop.size(); i++) {
        int node = -1;
        for (int j = 0; j < projectors.size(); j++) {
            if (projectors[j].colorId == loop[i]) {
                node = j;
                break;
            }
        }
        if (node == -1) {
            cout << "Failed to find projector in loop";
            return {};
        }
        // search for the resolved inputs of this node
        for (int j = 0; j < projectors[node].hitting.size(); j++) {
            if (projectors[node].hitting[j]->output->resolved) {
                totalInput = {
                    totalInput[0] + projectors[node].hitting[j]->output->value[0],
                    totalInput[1] + projectors[node].hitting[j]->output->value[1],
                    totalInput[2] + projectors[node].hitting[j]->output->value[2]
                };
            }
            if (projectors[node].hitting[j]->loopedColor[0] != -1) {
                storedColor = {
                    storedColor[0] + projectors[node].hitting[j]->loopedColor[0],
                    storedColor[1] + projectors[node].hitting[j]->loopedColor[1],
                    storedColor[2] + projectors[node].hitting[j]->loopedColor[2]
                };
            }
            if (!projectors[node].hitting[j]->requireCharge) {
                defaultColor = {
                    defaultColor[0] + projectors[node].hitting[j]->producedColor[0],
                    defaultColor[1] + projectors[node].hitting[j]->producedColor[1],
                    defaultColor[2] + projectors[node].hitting[j]->producedColor[2]
                };
            }
        }
    }
    // we now have our total input from all sources of the loop, now either use the input, or our stored value if input is 0
    if (totalInput[0] != 0 || totalInput[1] != 0 || totalInput[2] != 0) {
        return mergeColors({ 0, 0, 0 }, totalInput);
    }

    if (storedColor[0] != 0 || storedColor[1] != 0 || storedColor[2] != 0) {
        return mergeColors({ 0, 0, 0 }, storedColor);
    }
    // if we don't have input and also have no stored value yet, find the stored value from the sum of the loop's inner contents
    return mergeColors({ 0, 0, 0 }, defaultColor);
}

int main()
{
    // create the window
    sf::RenderWindow window(sf::VideoMode(BOUNDX, BOUNDY), "Radiant");

    vector<double> mouseLoc = {0, 0};
    auto lastTick = chrono::steady_clock::now();

    vector<Projector> projectors = {};

    // the top left loop
    projectors.push_back(Projector({ 250, 250 }, { 0,0,0 }, 0, true, false, { 0 }));
    projectors.push_back(Projector({ 400, 250 }, { 0,0,255 }, 0, false, false, { 0 }));
    projectors.push_back(Projector({ 400, 400 }, { 0,0,0 }, 0, true, false, { PI / 2, 3 * PI / 2 }));

    // the bottom right loop
    projectors.push_back(Projector({ BOUNDX - 250, BOUNDY - 250 }, { 0,0,0 }, 0, true, false, { 0 }));
    projectors.push_back(Projector({ BOUNDX - 400, BOUNDY - 250 }, { 0,255,0 }, 0, false, false, { 0 }));
    projectors.push_back(Projector({ BOUNDX - 400, BOUNDY - 400 }, { 0,0,0 }, 0, true, false, { PI / 2, 3 * PI / 2 }));

    // the input
    projectors.push_back(Projector({ 250, BOUNDY - 250 }, { 255,0,0 }, 0, false, false, { 0 }));

    vector<LineSegment> segments = {};
    for (int i = 0; i < 0; i++) {
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
        Event event;
        while (window.pollEvent(event)) {
            // check the type of the event
            switch (event.type) {
                // window closed
                case Event::Closed: {
                    window.close();
                    break;
                }

                // key pressed
                case Event::KeyPressed: {
                    bool left = (event.key.code == Keyboard::Left) || (event.key.code == Keyboard::A);
                    bool right = (event.key.code == Keyboard::Right) || (event.key.code == Keyboard::D);
                    if (left || right) {
                        // rotate selected lasers when buttons are pressed
                        for (int i = 0; i < projectors.size(); i++) {
                            if (!projectors[i].selected) continue;
                            projectors[i].rotating = left ? -1 : 1;
                        }
                    }
                    break;
                }

                // key released
                case Event::KeyReleased: {
                    bool left = (event.key.code == Keyboard::Left) || (event.key.code == Keyboard::A);
                    bool right = (event.key.code == Keyboard::Right) || (event.key.code == Keyboard::D);
                    if (left || right) {
                        for (int i = 0; i < projectors.size(); i++) {
                            if (!projectors[i].selected) continue;
                            projectors[i].rotating = 0;
                        }
                    }
                    break;
                }

                // mouse clicked
                case Event::MouseButtonPressed: {
                    if (event.mouseButton.button == Mouse::Left)
                    {
                        // check if we fall within any projectors, and select them if we do
                        for (int i = 0; i < projectors.size(); i++) {
                            if (projectors[i].inside(mouseLoc[0], mouseLoc[1])) {
                                projectors[i].selected = true;
                            }
                            else projectors[i].selected = false;
                        }
                    }
                    break;
                }

                case Event::MouseMoved: {
                    mouseLoc = { (double)event.mouseMove.x, (double)event.mouseMove.y };
                    for (int i = 0; i < projectors.size(); i++) {
                        //projectors[i].rotation = atan2(mouseLoc[1] - projectors[i].pos[1] - projectors[i].radius, mouseLoc[0] - projectors[i].pos[0] - projectors[i].radius);
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

        // clear data from the previous loop, and set up projectors for propogation and bounce detection down the line
        for (int i = 0; i < projectors.size(); i++) {
            projectors[i].inputs.clear();
            projectors[i].hitting.clear();
            projectors[i].colorId = i;
            projectors[i].loopChecked = false;
            delete projectors[i].output;
            projectors[i].output = new InputColor{ i, false, {} };
        }
        // update the bounces and hits for each projector
        for (int i = 0; i < projectors.size(); i++) {
            projectors[i].update(segments, projectors);
        }
        // create the input and output vectors for every projector
        for (int i = 0; i < projectors.size(); i++) {
            projectors[i].generateThroughput();
        }
        // collapse all color inputs until we either finish, or are only left with infinite loops
        bool done = true;
        for (int j = 0; j < 999; j++) {
            done = true;
            for (int i = 0; i < projectors.size(); i++) {
                if (projectors[i].collapseColorInputs()) {
                    projectors[i].loopedColor = { -1, -1, -1 };
                    done = false;
                }
            }
            if (done) {
                break;
            }
            if (j == 998) cout << "Failed to collapse color inputs" << endl;
        }

        for (int a = 0; a < 99; a++) {
            connectedLoops = {};

            // we now only have loops, so for each laser attempt to find the highest order loops
            for (int i = 0; i < projectors.size(); i++) {
                projectors[i].findHighestLoop({}, 0);
            }

            // we have the highest order loops pushed to our connected loops vector, solve these then repeat
            for (int i = 0; i < connectedLoops.size(); i++) {
                vector<int> loopColor = solveLoop(connectedLoops[i], projectors);
                for (int j = 0; j < connectedLoops[i].size(); j++) {
                    int node = -1;
                    for (int k = 0; k < projectors.size(); k++) {
                        if (projectors[k].colorId == connectedLoops[i][j]) {
                            node = k;
                            break;
                        }
                    }
                    if (node == -1) {
                        cout << "Failed to find projector in loop";
                        return {};
                    }
                    // assign its output and resolve it
                    projectors[node].laserColor = loopColor;
                    projectors[node].loopedColor = loopColor;
                    projectors[node].output->resolved = true;
                    projectors[node].output->value = loopColor;
                }
            }

            // check if every projector has an output
            bool done = true;
            for (int j = 0; j < projectors.size(); j++) {
                if (!projectors[j].output->resolved) done = false;
            }
            if (done) break;

            if (a == 98) cout << "Failed to calculate loops" << endl;
        }

        // update our render colors to match what we calculated
        for (int i = 0; i < projectors.size(); i++) {
            projectors[i].updateColors();
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