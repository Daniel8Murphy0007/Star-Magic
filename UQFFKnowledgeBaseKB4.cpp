#include "UQFFKnowledgeBaseKB4.h"
#include <iostream>
#include <string>

/*
PAPER_719: UQFF KB4 -- Red Dwarf Compression_B (43.b)
Source: grok_share_ba508f76c8e.txt entry #68

=== Drawing 32: Nebular Cloud Photo ===
  U_g4 = k4 * rho_SCm * M_BH / d_g * exp(-alpha*t) * cos(pi*t_n) * (1+f_feedback)
  U_g4 ~= 1.69e-2 * exp(-0.001t) * cos(pi*t) J/m^3
  Star positions: (100,900), (500,900), (900,900), (500,100), (200,100)
  Distances: d_12=400, d_23=400, d_13=800, d_24=800, d_45=300
  Angles: theta_-4=180 deg, theta_-5=90 deg, theta_-7=90 deg
  rho_SCm = 2.39e-22 J/m^3 (nebula scale level 13)

=== Drawing 33: Shock-Induced Star Formation ===
  U_g4 = k4 * rho_SCm * M_star / d_g * exp(-alpha*t) * cos(pi*t_n) * (1+f_shock)
  U_g4 ~= 3.49e-6 * exp(-0.001t) * cos(pi*t) J/m^3
  SiO / formamide astrochemical tracers, protostellar jets

=== LENR/Collider (referenced from 43.b) ===
  Same LENR + Higgs equations as 43.c (W+e-+p->n+nu_e, Higgs ~125 GeV)

self.version = "Session176"
*/

double UQFFKnowledgeBaseKB4::U_g4_nebular_BH(double t) const {
    // U_g4 = k4 * rho_SCm * M_BH / d_g * exp(-alpha*t) * cos(pi*t_n) * (1+f_feedback)
    double t_n = 0.0;
    return 1.0 * rho_SCm_nebula * M_BH_nebula / d_g_nebula
           * std::exp(-alpha_decay * t) * std::cos(M_PI * t_n) * (1.0 + f_feedback);
}

double UQFFKnowledgeBaseKB4::U_g4_shock_SF(double t) const {
    // U_g4 for shock-induced star formation (Drawing 33)
    double t_n = 0.0;
    return 1.0 * rho_SCm_nebula * M_star_SF / d_g_SF
           * std::exp(-alpha_decay * t) * std::cos(M_PI * t_n) * (1.0 + f_shock);
}

double UQFFKnowledgeBaseKB4::geometric_distance(int i, int j) const {
    if (i < 0 || i >= (int)star_positions.size()) return 0.0;
    if (j < 0 || j >= (int)star_positions.size()) return 0.0;
    double dx = star_positions[j].x - star_positions[i].x;
    double dy = star_positions[j].y - star_positions[i].y;
    return std::sqrt(dx*dx + dy*dy);
}

double UQFFKnowledgeBaseKB4::geometric_angle(int apex, int from, int to) const {
    // Angle at vertex 'apex' in the triangle from-apex-to
    if (apex < 0 || from < 0 || to < 0) return 0.0;
    if (apex >= (int)star_positions.size()) return 0.0;
    double ax = star_positions[from].x - star_positions[apex].x;
    double ay = star_positions[from].y - star_positions[apex].y;
    double bx = star_positions[to].x   - star_positions[apex].x;
    double by = star_positions[to].y   - star_positions[apex].y;
    double dot = ax*bx + ay*by;
    double mag = std::sqrt(ax*ax+ay*ay) * std::sqrt(bx*bx+by*by);
    if (mag < 1e-30) return 0.0;
    double cosA = dot / mag;
    if (cosA >  1.0) cosA =  1.0;
    if (cosA < -1.0) cosA = -1.0;
    return std::acos(cosA) * 180.0 / M_PI;
}

double UQFFKnowledgeBaseKB4::primary_equation() const {
    double Ug4_bh = U_g4_nebular_BH(curr_t);
    double Ug4_sf = U_g4_shock_SF(curr_t);
    return Ug4_bh + Ug4_sf;
}

std::string UQFFKnowledgeBaseKB4::description() const {
    return "PAPER_719: UQFF KB4 -- Drawing 32 nebular BH + Drawing 33 shock SF | Session176";
}

void UQFFKnowledgeBaseKB4::self_update() {
    curr_t += time_step;
}

void UQFFKnowledgeBaseKB4::self_expand() {
    f_feedback *= 1.001;
    f_shock    *= 1.001;
}

void UQFFKnowledgeBaseKB4::simulate(int num_steps) {
    for (int step = 0; step < num_steps; ++step) {
        double t = curr_t + step * time_step;
        double P = primary_equation();
        std::cout << "Step " << step
                  << "  U_g4_BH=" << U_g4_nebular_BH(t)
                  << "  U_g4_SF=" << U_g4_shock_SF(t)
                  << "  total=" << P << "\n";
        self_update();
    }
}

#ifdef STANDALONE_UQFFKNOWLEDGEBASEKB4
int main() {
    UQFFKnowledgeBaseKB4 kb4;
    std::cout << "UQFF KB4 -- Red Dwarf Compression_B Analysis\n";
    std::cout << kb4.description() << "\n";
    std::cout << "U_g4_BH(t=0) = " << kb4.U_g4_nebular_BH(0.0) << " J/m^3\n";
    std::cout << "U_g4_SF(t=0) = " << kb4.U_g4_shock_SF(0.0)   << " J/m^3\n";
    std::cout << "d_12 = " << kb4.geometric_distance(0,1) << " (normalized)\n";
    std::cout << "angle_1-2-3 = " << kb4.geometric_angle(1,0,2) << " deg\n";
    kb4.simulate(3);
    return 0;
}
#endif
