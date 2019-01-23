#include <stdio.h>
#include <cmath>
#include<iostream>
#define none 0,0,dt
#define cw (-1)
#define ccw (+1)
using namespace std;

double fric_s = 0.6, fric_k = 0.5;

class Ball {
private:
	double rx, ry, vx, vy, ax, ay; // Kinematic quantities
	double ar, av, aa; int dir; // Rotational quantities
	double m, k, radius, I; // Ball constants
public:
	Ball(double rx, double ry, double vx, double vy, double ax, double ay, double ar, double av, double aa, double m, double k, double r, int dir) :
		rx(rx), ry(ry), vx(vx), vy(vy), ax(ax), ay(ay), ar(ar), av(av), aa(aa), m(m), k(k), radius(r), dir(dir)
	{
		I = 2 * m*radius*radius / 5; // Moment of inertia for a solid sphere
	}
                    //  force   theta force   delta time
	void updateBall(double f, double t_f, double dt, int ball_num);

	// Accessors
	double getRadius() {
		return radius;
	}
	double getX() {
		return rx;
	}
	double getY() {
		return ry;
	}
	double getVx() {
	    return vx;
	}
	double getVy() {
        return vy;
	}
	double getK() {
		return k;
	}
	double getM() {
		return m;
	}
	double getV() {
        return hypot(vx, vy);
	}
	double getL() {
		// Angular momentum = L = I * av
		return I * av;
	}
	double getDir() {
	    return dir;
	}
};
void Ball::updateBall(double f, double t_f, double dt, int dir) {
	// F = ma -> a = f/m, calculate acceleration for this frame
	ax = f*cos(t_f) / m;
	ay = f*sin(t_f) / m;

	//s = s0 + v0t + 1/2at^2, calculate position for this frame
	rx = rx + vx*dt + ax*dt*dt / 2;
	ry = ry + vy*dt + ay*dt*dt / 2;

	// f, the sum of the spring forces, is also the normal force
	// F_f = fric_k * F_N -> T = r * F_f -> aa = T / I = r * frick_k * F_N / I
	// Note: r is ball radius minus compression
	aa = dir * abs (radius * fric_k * f / I);

	//ar = ar + av0*t + 1/2*aa*t^2
	ar = ar + av*dt + aa*dt*dt / 2;

	// v = v0 + at, set v0 for next frame
	vx = vx + ax*dt;
	vy = vy + ay*dt;

	// av = av0 + aa*t
	av = av + aa*dt;




}

// smh using global variables
Ball b1(0, 0, 0, 0, 0, 0, 0, 0, 0, 1.2, 1000, .5, cw); // x, y, v_x, v_y, a_x, a_y, ar, av, aa, m (kg), k (Nm-1), radius (m), spin dir
Ball b2(7, 0, -7.2, 0, 0, 0, 0, 0, 0, 3.6, 1000, .5, cw);

double c1 = 0, c2 = 0, t_c1 = 0, t_c2 = 0; // magnitude and direction of compression
bool col = false; // assume not initally colliding
double ttt=0;
int dir=0;
void collide(double dt) {
	cout << b1.getL() << "," << b2.getL() << ",";
	cout << b1.getX() << "," << b1.getY() << "," << b2.getX() << "," << b2.getY() <<","<< ttt << "\n";
	if (b1.getRadius() + b2.getRadius() <= hypot(b1.getX() - b2.getX(), b1.getY() - b2.getY())) { // If the radii of two balls added together (assuming perfect circles) is less than distance between the midpoints of the balls
		b1.updateBall(none,0); // No force applied to ball
		b2.updateBall(none,0);
		col = false; // The ball is not currently colliding
	}
	else {
		if (!col) { // If ball just started colliding
			b1.updateBall(none,1); // Allow simulation to continue for one more dt
			b2.updateBall(none,2);
			double d = b1.getRadius() + b2.getRadius() - hypot(b1.getX() - b2.getX(), b1.getY() - b2.getY()); // Find the total compression

			c1 = d * b2.getK() / (b1.getK() + b2.getK()); // After initally colliding, assume that each ball contributes a displacement proportional to the k to the collision
			c2 = d * b1.getK() / (b1.getK() + b2.getK());
			// The angle of compression is + to radical axis i.e. || to line segment connecting midpoints of balls
			t_c1 = atan2((b2.getY() - b1.getY()), (b2.getX() - b1.getX())); // Angle between midpoint of b2 to midpoint of b1
			t_c2 = atan2((b1.getY() - b2.getY()), (b1.getX() - b2.getX())); // Angle between midpoint of b1 to midpoint of b2

            // theta of movement relative to radical axis
            double t_d = atan2(b1.getVy(), b1.getVx()) - abs(t_c1 - asin(1));

            // determine direction
           // if (t_d < 0.001) {
             //   dir=0; // head on collision
            if (t_d > asin(1)) {
                dir=-1; // cw
            } else dir=1; // ccw


			col = true;
		}
		else { // If it is currently colliding
			double x1i = b1.getX(), y1i = b1.getY(), x2i = b2.getX(), y2i = b2.getY(); // Store initial conditions of ball

																					   // N3 states equal and opposite reaction, sum the spring forces of the two balls
																					   // F = Fs1 + Fs2 = (-k1x1) + (-k2x2) = -(k1x1+k2x2)

			double f = -(b1.getK()*c1 + b2.getK()*c2); // Calculate force

													   // Update conditions of balls
			b1.updateBall(f, t_c1, dt, dir);
			b2.updateBall(f, t_c2, dt, -dir);
			double t_s1 = atan2(b1.getY() - y1i, b1.getX() - x1i); // Calculate angles of displacement for each ball
			double t_s2 = atan2(b2.getY() - y2i, b2.getX() - x2i);

			// Calculate new compression, find the component of displacement that lies on compression angle
			c1 += hypot(b1.getX() - x1i, b1.getY() - y1i) * cos(t_c1 - t_s1);
			c2 += hypot(b2.getX() - x2i, b2.getY() - y2i) * cos(t_c2 - t_s2);
		}
	}

}

int main()
{
	double t = 2*hypot(b1.getY()-b2.getY(), b1.getX()-b2.getX())/((b1.getV()+b2.getV())/2);
	double dt = t/250;
	freopen("output.csv", "w", stdout);
	// NEED A CASE FOR BALLS ALREADY TOUCHING
	cout << "angular momentum b1, angular momentum b2, b1x, b1y, b2x, b2y, t\n";
	for (int i = 0; i <= t / dt; i++) {
        ttt+=dt;
		//cout << (i*dt) << " ";
		collide(dt);
	}
}

