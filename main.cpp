#include <iostream>

#include <vector>
#include <armadillo>

#include "zignor.h"
#include "zigrandom.h"

using namespace std;
using namespace arma;

struct Particle {
    double x;
    double y;
    double z;
};

struct System {
    int L;
    int W;
    int H;
    double prefac;
};

bool diffuse(Particle *particle, const System *system, const int ndiff = 100) {

    const int &L = system->L;
    const int &W = system->W;
    const int &H = system->H;

    const double &x0 = particle->x;
    const double &y0 = particle->y;
    const double &z0 = particle->z;

    double x1 = x0;
    double y1 = y0;
    double z1 = z0;

    for (int n = 0; n < ndiff; ++n) {
        const double dx = system->prefac*DRanNormalZig32();
        const double dy = system->prefac*DRanNormalZig32();

        x1 = x0 + dx;
        y1 = y0 + dy;

        if (x1 > L - 1) {
            const double delta = x1 - (L-1);
            x1 = (L-1) - delta;
        }
        else if (x1 < 0)
        {
            x1 = -x1;
        }

        if (y1 > W - 1) {
            return false;
        }
        else if (y1 < 0) {
            return false;
        }

        bool z_legal = false;

        while (!z_legal) {
            const double dz = system->prefac*DRanNormalZig32();
            z1 = z0 + dz;
            z_legal = (z1 < H - 1) && (z1 > 0);
        }
    }

    particle->x = x1;
    particle->y = y1;
    particle->z = z1;

    return true;
}

Particle *getNewParticle(const System *system, const bool yZero = false) {
    Particle *particle = new Particle;
    particle->x = DRan_MWC8222()*(system->L-1);
    particle->z = DRan_MWC8222()*(system->H-1);

    if (yZero) {
        particle->y = 0.5;
    }
    else {
        particle->y = DRan_MWC8222()*(system->W-1);
    }

    return particle;
}

int main() {

    int inseed = time(nullptr);
    int cseed = 100;
    int seed2 = inseed * 3;
    RanSetSeed_MWC8222(&seed2, cseed);
    RanNormalSetSeedZig32(&inseed, 5);

    System system;
    system.L = 100;
    system.W = 100;
    system.H = 15;

    const int nb = 30;
    mat hist(nb, nb, fill::zeros);
    const double deltax = (system.L-1)/double(nb);
    const double deltay = (system.W-1)/double(nb);

    const double dt = 0.01;
    const double D = 1.0;
    system.prefac = sqrt(2*D*dt);

    const int N0 = 0;
    vector<Particle *> particles;

    for (int i = 0; i < N0; ++i) {
        particles.push_back(getNewParticle(&system));
    }

    const int nSteps = 10000000;
    const int interval = 1000;
    const int fluxInterval = 10;

    vector<Particle*> removalQueue;
    for (int step = 0; step < nSteps; ++step) {
        if (step % fluxInterval == 0) {
            particles.push_back(getNewParticle(&system, true));
        }

        if (step % interval == 0) {
            cout.flush();
            cout << "\r" << step+1 << "/" << nSteps << " " << particles.size();

            for (const Particle* particle : particles) {
                hist(int(particle->x/deltax), int(particle->y/deltay))++;
            }

            stringstream s;
            s << "/tmp/fcav_evo_" << step/interval << ".arma";
            hist.save(s.str());
        }

        for (int i = particles.size() - 1; i >= 0; i--) {
            Particle *particle = particles.at(i);
            bool success = diffuse(particle, &system, 1);

            if (!success) {
                removalQueue.push_back(particle);
            }
        }

        for (int i = removalQueue.size() - 1; i >= 0; --i) {
            particles.erase(remove(particles.begin(),
                                   particles.end(),
                                   removalQueue.at(i)),
                            particles.end());
        }

        removalQueue.clear();

    }

    cout << endl;
    return 0;
}

























