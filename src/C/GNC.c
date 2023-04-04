#include <stdbool.h>
#include <math.h>

struct input;
struct input
{
    double r[3];
    double v[3];
    double m_prop;
};

struct commands;
struct commands
{
    bool NOOP;
};

struct parameters;
struct parameters
{
    double mu;
};

struct output;
struct output
{
    bool thr_on;
    double thr_dir[3];
    double specific_energy;
};

struct output step(const struct input *i, struct commands *c, const struct parameters *p)
{
    bool thr_on;
    double thr_dir[3];
    double specific_energy;

    double R2 = 0.0;
    double V2 = 0.0;
    for (int j = 0; j < 3; j++)
    {
        R2 += i->r[j] * i->r[j];
        V2 += i->v[j] * i->v[j];
    }

    specific_energy = V2 / 2 - p->mu / sqrt(R2);

    thr_on = (i->m_prop > 0) && (specific_energy < 0);

    if (thr_on && V2 > 0)
    {
        double V = sqrt(V2);

        for (int j = 0; j < 3; j++)
        {
            thr_dir[j] = i->v[j] / V;
        }
    }
    else
    {
        thr_dir[0] = 0;
        thr_dir[1] = 0;
        thr_dir[2] = 0;
    }

    struct output o = {
        .thr_on = thr_on,
        .thr_dir = *thr_dir,
        .specific_energy = specific_energy,
    };

    return o;
}
