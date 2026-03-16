/**
 * ================================================================================================
 * Header: UQFFLearningAssessment.h
 * 
 * Description: C++ Module for UQFF Learning Assessment Evolution_B Class
 *              This is the ninth module in a series of 500+ code files for the Universal Quantum
 *              Field Framework (UQFF) simulations, focusing on framework assessment and advancement
 *              metrics derived from recent examples (Westerlund 2, Pillars of Creation, Rings of Relativity),
 *              using insights from Hubble datasets, high-energy lab simulations, and UQFF refinements
 *              (dated May 09, 2025, updated for full integration on October 08, 2025).
 * 
 * Purpose: Encapsulates assessment of UQFF learning and advancement from the last three examples.
 *          Computes metrics for diverse regimes, dynamic processes, scalability, and overall advancement.
 *          Includes parameters from each example for comparative analysis. Supports dynamic variable updates.
 * 
 * Integration: Designed for inclusion in base program 'ziqn233h.cpp' (not present here).
 *              Instantiate class in main: UQFFLearningAssessment assess;
 *              Compute: double adv = assess.compute_advancement();
 * 
 * Key Features:
 *   - Default values aggregated from previous UQFF examples (e.g., M_wd2 = 30,000 Msun, etc.).
 *   - Metrics: diversity_score (number of regimes), dynamic_score (new terms introduced), scalability_score.
 *   - Advancement = (diversity + dynamic + scalability) / 3.0 * 100.0 (%).
 *   - Setter methods for updates: setVar(double new_val) or addToVar(double delta)/subtractFromVar(double delta).
 * 
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript.
 * Date: October 08, 2025
 * Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
 * ================================================================================================
 */

#ifndef UQFF_LEARNING_ASSESSMENT_H
#define UQFF_LEARNING_ASSESSMENT_H

#include <iostream>
#include <cmath>
#include <iomanip>

class UQFFLearningAssessment {
private:
    // Core assessment parameters (mutable for updates)
    double diversity_score;     // Number of physical regimes (default 3)
    double dynamic_score;       // Number of new dynamic terms (default 3: wind, erosion, lensing)
    double scalability_score;   // Adaptability across scales (default 8.0/10)
    double f_TRZ;               // Time-reversal factor (shared across examples)
    double rho_vac_UA;          // UA vacuum density (J/m^3)
    double rho_vac_SCm;         // SCm vacuum density (J/m^3)

    // Parameters from Westerlund 2 example
    double M_wd2;               // Initial mass (kg)
    double r_wd2;               // Radius (m)
    double tau_SF_wd2;          // Star formation timescale (s)
    double rho_wind_wd2;        // Wind density (kg/m^3)
    double v_wind_wd2;          // Wind velocity (m/s)

    // Parameters from Pillars of Creation example
    double M_pillars;           // Initial mass (kg)
    double r_pillars;           // Radius (m)
    double tau_SF_pillars;      // Star formation timescale (s)
    double E_0_pillars;         // Initial erosion factor
    double tau_erosion_pillars; // Erosion timescale (s)
    double rho_wind_pillars;    // Wind density (kg/m^3)
    double v_wind_pillars;      // Wind velocity (m/s)

    // Parameters from Rings of Relativity example
    double M_rings;             // Lensing mass (kg)
    double r_rings;             // Einstein radius (m)
    double Hz_rings;            // Hubble parameter at z (s^-1)
    double L_factor_rings;      // Lensing factor
    double rho_wind_rings;      // Wind density (kg/m^3)
    double v_wind_rings;        // Wind velocity (m/s)

    // Shared constants
    double G;                   // Gravitational constant
    double Lambda;              // Cosmological constant
    double c_light;             // Speed of light
    double hbar;                // Reduced Planck's constant
    double t_Hubble;            // Hubble time (s)
    double t_Hubble_gyr;        // Hubble time in Gyr
    double M_DM_factor;         // Dark matter mass fraction
    double delta_rho_over_rho;  // Density perturbation fraction
    double rho_fluid_shared;