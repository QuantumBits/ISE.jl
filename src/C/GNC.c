#include <stdbool.h>
#include <math.h>

struct gnc_i;
struct gnc_i
{
    double r[3];
    double v[3];
    double m_prop;
};

struct gnc_c;
struct gnc_c
{
    bool NOOP;
};

struct gnc_p;
struct gnc_p
{
    double mu;
};

struct gnc_o;
struct gnc_o
{
    bool thr_on;
    double thr_dir[3];
    double specific_energy;
};

double specific_energy(const double mu, const double r[3], const double v[3]);
double specific_energy(const double mu, const double r[3], const double v[3])
{
    double R = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
    double V = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);

    return V * V / 2 - mu / R;
}

void gnc_step(const struct gnc_i *i, struct gnc_c *c, const struct gnc_p *p, struct gnc_o *o)
{
    o->specific_energy = specific_energy(p->mu, i->r, i->v);

    o->thr_on = (i->m_prop > 0) && (o->specific_energy < 0);

    if (o->thr_on)
    {
        double V = sqrt(i->v[0] * i->v[0] + i->v[1] * i->v[1] + i->v[2] * i->v[2]);

        o->thr_dir[0] = i->v[0] / V;
        o->thr_dir[1] = i->v[1] / V;
        o->thr_dir[2] = i->v[2] / V;
    }
    else
    {
        o->thr_dir[0] = 0;
        o->thr_dir[1] = 0;
        o->thr_dir[2] = 0;
    }
}

int main()
{
    struct gnc_i i;
    struct gnc_c c;
    struct gnc_p p;
    struct gnc_o o;

    i.r[0] = 6371009.0; // m
    i.r[1] = 0;
    i.r[2] = 0;

    i.v[0] = 0;
    i.v[1] = 7672.59355; // m/s
    i.v[2] = 0;

    i.m_prop = 100.0;

    p.mu = 3.986004418e14;

    gnc_step(&i, &c, &p, &o);

    return 0;
}