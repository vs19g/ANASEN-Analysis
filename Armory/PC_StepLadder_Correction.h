#include <TF1.h>

double model_invert(double *y, double *p)
{
    double result = y[0];
    double slope = 0.52;
    double z_grid[8] = {147.998, 101.946, 59.7634, 19.6965, -19.6965, -59.7634, -101.946, -147.998};
    for (int i = 0; i < 7; i++)
    {
        if (y[0] <= z_grid[i] && y[0] > z_grid[i + 1])
        {
            double zavg = (z_grid[i] + z_grid[i + 1]) * 0.5; // midpoint about which we pivot
            result = (y[0] - zavg) / slope + zavg;
            break;
        }
    }
    return result;
}
// ---- Parabola variant (now takes optional slope; out-of-band -> Si, no clamp-snap) ----
inline double model_invert_a1c1(double cfrac, double z_si, double z_pivot,
                                double C, const double* c0_table, bool& ok,
                                double slope = 1.0)
{
    ok = false;
    static const double z_grid[8] =
        {147.998, 101.946, 59.7634, 19.6965, -19.6965, -59.7634, -101.946, -147.998};
    const double cfrac_hi = 0.90, cfrac_lo = 0.12;
    if (cfrac > cfrac_hi || cfrac < cfrac_lo) return z_si;

    const double sgn = (z_si >= z_pivot) ? +1.0 : -1.0;
    int gnb = -1; double z_nb = (sgn > 0) ? +1.0e30 : -1.0e30;
    for (int i = 0; i < 8; ++i) {
        if (sgn > 0 && z_grid[i] > z_pivot + 1e-6 && z_grid[i] < z_nb) { z_nb = z_grid[i]; gnb = i; }
        if (sgn < 0 && z_grid[i] < z_pivot - 1e-6 && z_grid[i] > z_nb) { z_nb = z_grid[i]; gnb = i; }
    }
    if (gnb < 0) return z_si;
    const double W = TMath::Abs(z_nb - z_pivot);
    if (W <= 0.0) return z_si;
    int cell = (sgn > 0) ? gnb : gnb - 1;
    if (cell < 0) cell = 0; if (cell > 6) cell = 6;

    double s2 = (c0_table[cell] - cfrac) / C;
    if (s2 < 0.0 || s2 > 1.0) return z_si;          // out-of-band -> Si (no boundary snap)
    double s = TMath::Sqrt(s2);
    double z_rec = z_pivot + sgn * s * W;
    if (slope > 0.0 && slope != 1.0) {
        const double zc = 0.5 * (z_grid[cell] + z_grid[cell + 1]);
        z_rec = (z_rec - zc) / slope + zc;
    }
    if (TMath::Abs(z_rec - z_si) > W) return z_si;
    ok = true;
    return z_rec;
}

// ---- Linear centre-fold variant ----
inline double model_invert_a1c1_linear(double cfrac, double z_si, double z_pivot,
                                        double k, const double* cfmin_table, bool& ok)
{
    ok = false;
    static const double z_grid[8] =
        {147.998, 101.946, 59.7634, 19.6965, -19.6965, -59.7634, -101.946, -147.998};
    const double cfrac_hi = 0.90, cfrac_lo = 0.12;
    if (cfrac > cfrac_hi || cfrac < cfrac_lo) return z_si;

    int cell = -1;
    for (int i = 0; i < 7; ++i)
        if (z_si <= z_grid[i] && z_si > z_grid[i + 1]) { cell = i; break; }
    if (cell < 0) return z_si;
    const double zc   = 0.5 * (z_grid[cell] + z_grid[cell + 1]);
    const double half = 0.5 * (z_grid[cell] - z_grid[cell + 1]);
    if (half <= 0.0 || k <= 0.0) return z_si;

    double f = (cfrac - cfmin_table[cell]) / k;
    if (f < 0.0) f = 0.0;
    if (f > 1.0) return z_si;
    const double d   = f * half;
    const double sgn = (z_pivot >= zc) ? +1.0 : -1.0;
    const double z_rec = zc + sgn * d;
    if (TMath::Abs(z_rec - z_si) > half) return z_si;
    ok = true;
    return z_rec;
}