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
    vector<int> laserColor = { 0, 0, 0 };
    vector<int> producedColor = { 255, 255, 255 };
    vector<int> loopedColor = { -1, -1, -1 };
    vector<vector<Projector*>> infiniteLoops = {};
    double rotation = 0;
    double radius = 25;
    double borderSize = 0.2;
    bool requireCharge = false;
    bool locked = false;
    vector<Projector*> hittingProjectors = {};
    bool selected = false;
    double rotating = 0;

    vector<vector<double>> bounces = {};
    // full constructor in case
    Projector(vector<double> pos, vector<int> borderColor, vector<int> laserColor, double rotation, double radius, bool requireCharge, bool locked, vector<int> producedColor) {
        shape.setRadius(radius);
        this->pos = { pos[0], pos[1] };
        this->borderColor = borderColor;
        this->laserColor = laserColor;
        this->producedColor = producedColor;
        this->rotation = rotation;
        this->radius = radius;
        this->requireCharge = requireCharge;
        this->locked = locked;
    }
    // practically, we only need these from our constructor
    Projector(vector<double> pos, vector<int> producedColor, double rotation, bool requireCharge, bool locked) {
        shape.setRadius(radius);
        this->pos = { pos[0], pos[1] };
        this->borderColor = { 255, 255, 255 };
        this->producedColor = producedColor;
        if (!requireCharge) this->laserColor = producedColor;
        this->rotation = rotation;
        this->radius = 25;
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

    void calculateBounces(vector<LineSegment> &segments, vector<Projector> &projectors) {
        vector<double> lastPos = { pos[0] + radius, pos[1] + radius };
        double lastDir = rotation;
        bounces = { lastPos };
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
                bounces.push_back(closestIntersection);
                projectors[projectorId].hittingProjectors.push_back(this);
                break;
            }

            lastPos = closestIntersection;
            lastDir = segments[closestIndex].getBounce({ lastPos[0], lastPos[1], cos(lastDir), sin(lastDir) });
            bounces.push_back(lastPos);
            lastBounce = closestIndex;
            //cout << "Bounced " << closestIndex << " with strand with colors " << segments[closestIndex].color[0] << ", " << segments[closestIndex].color[1] << ", " << segments[closestIndex].color[2] << endl;
            if (segments[closestIndex].reflectivity == 0) break;
        }
    }

    // calculates the bounces of the laser, and alerts hit objects
    void update(vector<LineSegment> &segments, vector<Projector> &projectors) {
        // if rotating, rotate us in the right direction
        rotation += rotating * 3.0 / 60.0;

        // calculate the bounces of the laser from various surfaces
        calculateBounces(segments, projectors);

        // remove our laser color, we'll resolve for it later anyways
        if (requireCharge) laserColor = { 0, 0, 0 };
        else laserColor = producedColor;
    }

    // propogates the colors from all lasers hitting this by back-calculating
    /*vector<int> propogateColor(Projector* parentCall, bool embraceInfinity, bool firstLoop, vector<int> trackedIndex = {}) {
        if (firstLoop) parentCall = this;

        // if we aren't hit, return our color (or lack thereof), unless we are outside an infinite loop being calculated for
        if (hittingProjectors.size() == 0) {
            if (requireCharge) return { 0, 0, 0 };
            else return producedColor;
        }

        // otherwise, sum all colors hitting this, and if in an infinite loop only care about those within the loop
        // first get our own color as a basis for addition
        vector<int> totalColor;
        if (requireCharge) totalColor = { 0, 0, 0 };
        else totalColor = embraceInfinity ? laserColor : producedColor;

        // loop through, until we reach a basis or ourselves we sum everything in between
        for (int i = 0; i < hittingProjectors.size(); i++) {
            Projector proj = *hittingProjectors[i];
            vector<int> foundColor;
            // if we stare into the infinite and it stares back, look away quickly to preserve what little sanity we have left
            if (parentCall == hittingProjectors[i]) {
                if (embraceInfinity) foundColor = { 0, 0, 0 };
                else foundColor = {-1};
            }
            else foundColor = proj.propogateColor(parentCall, embraceInfinity, false, trackedIndex);

            // if an infinite loop propogates back up, rerun our propogation to only care about produced colors (since they originate in the loop)
            if (foundColor[0] == -1) {
                foundColor.push_back(i); // keep track of the indeces this occurs at
                if (firstLoop) return propogateColor(parentCall, true, false, foundColor);
                return foundColor;
            }

            // otherwise, sum them like normal
            totalColor = { totalColor[0] + foundColor[0], totalColor[1] + foundColor[1] , totalColor[2] + foundColor[2] };
        }
        
        // finally, once everyone has responded to us, return the total color we found, which will be clamped later
        return totalColor;
    }*/

    // checks down every hitting vector to see if a specific one is found, returning the path to it if so (with the target last in the path)
    // for example, if L1 hits L2, which hits L3, which then hits L1 again, and L1 is passed as a target, this would return loop of { L1, L3, L2 }
    vector<Projector*> findLoop(Projector* target) {
        // if we don't have it, send an empty response
        if (hittingProjectors.size() == 0) return {};
        for (int i = 0; i < hittingProjectors.size(); i++) {
            // if we have it, start propogating the infinite loop backwards
            if (target == hittingProjectors[i]) return { hittingProjectors[i] };
            // if we don't find it here, look deeper
            Projector proj = *hittingProjectors[i];
            vector<Projector*> propogation = proj.findLoop(target);
            if (propogation.size() == 0) continue;
            // if looking deeper found it, add ourselves to the list and keep passing it back up
            propogation.push_back(hittingProjectors[i]);
            return propogation;
        }
        // if none of our branches found anything, send an empty response
        return {};
    }

    // returns an empty vector if not, or a vector of vectors, with each separate loop of projectors
    vector<vector<Projector*>> inInfiniteLoop() {
        vector<vector<Projector*>> loops = {};
        for (int i = 0; i < hittingProjectors.size(); i++) {
            // for every projector hitting us, check if we appear in their propogated hit list, and if so add it as an infinite loop
            vector<Projector*> thisLoop = findLoop(this);
            if (thisLoop.size() > 0) loops.push_back(thisLoop);
        }
        return loops;
    }

    // propogates color all the way down, reassigning over loop color as needed
    vector<int> propogateColor(bool onlyStatics = false) {
        // looped objects already know what color they want to be from prior calculations, so just use that
        if (infiniteLoops.size() > 0) {
            if (onlyStatics) return {0, 0, 0};
            return loopedColor;
        }

        vector<int> totalColor = { 0, 0, 0 };
        if (!requireCharge) totalColor = producedColor;

        // go through every projector hitting this and request their color too
        for (int i = 0; i < hittingProjectors.size(); i++) {
            // combine the colors of all things hitting us to get our color
            Projector proj = *hittingProjectors[i];
            vector<int> mergedColor = proj.propogateColor();
            totalColor = { totalColor[0] + mergedColor[0], totalColor[1] + mergedColor[1], totalColor[2] + mergedColor[2] };
        }

        return totalColor;
    }

    // propogates color through loops
    void propogateLoopColor(bool firstLoop) {
        // if we aren't in an infinite loop, remove any loop color we have
        if (infiniteLoops.size() == 0) {
            loopedColor = { -1, -1, -1 };
            return;
        }

        // for each loop, if it doesn't already have a color, calculate the base color for the loop
        vector<int> totalColor = { 0, 0, 0 };
        vector<Projector*> uniqueProjectors = {};
        for (vector<Projector*> loop : infiniteLoops) for (Projector* projPointer : loop) {
            // don't include the same projector twice
            bool found = false;
            for (Projector* projPointer2 : uniqueProjectors) if (projPointer == projPointer2) found = true;
            if (found) continue;
            uniqueProjectors.push_back(projPointer);

            // if our projector produces light, average it into the loop's average
            Projector proj = *projPointer;
            if (!proj.requireCharge) {
                totalColor[0] += proj.producedColor[0];
                totalColor[1] += proj.producedColor[1];
                totalColor[2] += proj.producedColor[2];
            }
        }

        // we have found the loop average, now we need to check if our loop has inputs we will use instead
        vector<int> loopInputs = { 0, 0, 0 };
        for (int i = 0; i < uniqueProjectors.size(); i++) {
            Projector proj = *uniqueProjectors[i];
            // go through every projector in our loop and find all inputs of light to our loop, if any exist
            for (int j = 0; j < proj.hittingProjectors.size(); j++) {
                Projector child = *proj.hittingProjectors[j];
                vector<int> childStatics = child.propogateColor(true);
                loopInputs = { loopInputs[0] + childStatics[0], loopInputs[1] + childStatics[1], loopInputs[2] + childStatics[2] };
            }
        }
        if (loopInputs[0] == 0 && loopInputs[1] == 0 && loopInputs[2] == 0) {
            // if we don't have any inputs, clamp our averaged loop color, we are done, yay!
            if (loopedColor[0] == -1 && loopedColor[1] == -1 && loopedColor[2] == -1) loopedColor = mergeColors({ 0, 0, 0 }, totalColor);
        }
        else {
            // if we have input to our loop, our looped color becomes that input
            loopedColor = mergeColors({ 0, 0, 0 }, loopInputs);
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

int main()
{
    // create the window
    sf::RenderWindow window(sf::VideoMode(BOUNDX, BOUNDY), "Radiant");

    vector<double> mouseLoc = {0, 0};
    auto lastTick = chrono::steady_clock::now();

    vector<Projector> projectors = {};
    projectors.push_back(Projector({ BOUNDX/2 - 125, BOUNDY / 2 - 25 }, {255,0,0}, 0, false, false));
    projectors.push_back(Projector({ BOUNDX / 2 + 75, BOUNDY / 2 + 25 }, { 0,255,0 }, 0, false, false));

    projectors.push_back(Projector({ BOUNDX / 2 - 125, BOUNDY / 2 - 225 }, { 0,0,255 }, 0, false, false));
    projectors.push_back(Projector({ BOUNDX / 2 + 75, BOUNDY / 2 + 225 }, { 0,0,0 }, 0, true, false));

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

        // update projector's bounces and stuff
        for (int i = 0; i < projectors.size(); i++) {
            projectors[i].hittingProjectors.clear();
        }
        // update the bounces and hits for each projector
        for (int i = 0; i < projectors.size(); i++) {
            projectors[i].update(segments, projectors);
        }
        // calculate the loops every projector is in
        for (int i = 0; i < projectors.size(); i++) {
            projectors[i].infiniteLoops = projectors[i].inInfiniteLoop();
        }
        for (int i = 0; i < projectors.size(); i++) {
            projectors[i].propogateLoopColor(true);
        }
        // propogate colors from hit projectors
        vector<vector<int>> savedLaserColors = {};
        for (int i = 0; i < projectors.size(); i++) {
            savedLaserColors.push_back(projectors[i].propogateColor());
        }
        // go through and update every laser at once, that way order does not affect how colors mix
        for (int i = 0; i < projectors.size(); i++) {
            projectors[i].laserColor = mergeColors({ 0, 0, 0 }, savedLaserColors[i]);
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