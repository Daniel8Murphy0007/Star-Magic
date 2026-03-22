"""
gen_muge_assessment.py — Generator for UQFFLearningAssessment.h
Assessment metrics module for UQFF framework advancement tracking.
Computes diversity_score, dynamic_score, scalability_score, and
overall advancement percentage.

Run:  python gen_muge_assessment.py
Output: UQFFLearningAssessment.h
"""
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_FILE = os.path.join(SCRIPT_DIR, "UQFFLearningAssessment.h")


def get_content():
    return r"""/**
 * ================================================================================================
 * Header: UQFFLearningAssessment.h
 *
 * Description: UQFF Framework Advancement Metrics Module
 *              Aggregates parameters from MUGE modules to assess framework advancement:
 *              diversity, dynamic range, and scalability scores.
 *
 * Purpose: Compute advancement = (diversity + dynamic + scalability) / 3 * 100 (%)
 *          Inputs drawn from modules: Westerlund2, PillarsOfCreation, RingsOfRelativity.
 *
 * Metrics:
 *   diversity_score  = count of unique physics term types / total possible
 *   dynamic_score    = count of time-varying terms / total terms
 *   scalability_score = log10(M_max / M_min) / log10(1e20)  (mass range coverage)
 *
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript.
 * Date: October 08, 2025
 * Copyright: Daniel T. Murphy
 * ================================================================================================
 */

#ifndef UQFF_LEARNING_ASSESSMENT_H
#define UQFF_LEARNING_ASSESSMENT_H

#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <string>
#include <map>

namespace UQFF {

class UQFFLearningAssessment {
private:
    // Unique physics term type counts
    int n_unique_term_types;    // Number of distinct physics term categories present
    int n_total_possible;       // Total possible term categories in framework
    int n_dynamic_terms;        // Time-varying terms count
    int n_total_terms;          // Total terms count

    // Mass range for scalability (kg)
    double M_min;               // Smallest mass in assessed module suite
    double M_max;               // Largest mass in assessed module suite
    double M_scale_reference;   // Reference mass range scale factor (1e20)

    // Computed scores
    double diversity_score;     // [0,1]
    double dynamic_score;       // [0,1]
    double scalability_score;   // [0,1]
    double advancement_pct;     // Percent (0-100)

    // Per-module assessment records
    struct ModuleRecord {
        std::string name;
        int term_count;
        int dynamic_count;
        double mass_kg;
        std::vector<std::string> unique_terms;
    };
    std::vector<ModuleRecord> module_records;

public:
    UQFFLearningAssessment() { initializeDefaults(); }
    ~UQFFLearningAssessment() {}

    void initializeDefaults() {
        // Based on first 9 assessed modules (Westerlund2, Pillars, Rings + base suite)
        n_unique_term_types = 9;   // base, Hz, B-decay, erosion, expansion, lensing, wind, merge, SC
        n_total_possible = 15;     // Full UQFF term taxonomy
        n_dynamic_terms = 6;       // E(t), B(t), M(t), I(t), F(t), cum_D
        n_total_terms = 12;        // Standard MUGE term set

        double M_sun = 1.989e30;
        M_min = 1.4 * M_sun;       // Magnetar
        M_max = 1.0e14 * M_sun;    // Rings cluster
        M_scale_reference = 1.0e20;

        computeScores();
    }

    void computeScores() {
        diversity_score = static_cast<double>(n_unique_term_types) / n_total_possible;
        dynamic_score = static_cast<double>(n_dynamic_terms) / n_total_terms;
        double mass_range = std::log10(M_max / M_min);
        double mass_ref = std::log10(M_scale_reference);
        scalability_score = std::min(mass_range / mass_ref, 1.0);
        advancement_pct = (diversity_score + dynamic_score + scalability_score) / 3.0 * 100.0;
    }

    // Add a module record for aggregated assessment
    void addModule(const std::string& name, int terms, int dynamic_terms,
                   double mass_kg, const std::vector<std::string>& unique_terms) {
        ModuleRecord rec;
        rec.name = name;
        rec.term_count = terms;
        rec.dynamic_count = dynamic_terms;
        rec.mass_kg = mass_kg;
        rec.unique_terms = unique_terms;
        module_records.push_back(rec);

        // Update mass range
        if (mass_kg < M_min) M_min = mass_kg;
        if (mass_kg > M_max) M_max = mass_kg;

        // Update dynamic count
        n_dynamic_terms += dynamic_terms;
        n_total_terms += terms;

        computeScores();
    }

    // Setters
    bool setVariable(const std::string& varName, double newValue) {
        if (varName == "n_unique_term_types")
            { n_unique_term_types = static_cast<int>(newValue); }
        else if (varName == "n_total_possible")
            { n_total_possible = static_cast<int>(newValue); }
        else if (varName == "n_dynamic_terms")
            { n_dynamic_terms = static_cast<int>(newValue); }
        else if (varName == "n_total_terms")
            { n_total_terms = static_cast<int>(newValue); }
        else if (varName == "M_min") { M_min = newValue; }
        else if (varName == "M_max") { M_max = newValue; }
        else {
            std::cerr << "Error: Unknown variable '" << varName << "'.\n";
            return false;
        }
        computeScores();
        return true;
    }

    double getVariable(const std::string& varName) const {
        if (varName == "diversity_score") return diversity_score;
        if (varName == "dynamic_score") return dynamic_score;
        if (varName == "scalability_score") return scalability_score;
        if (varName == "advancement_pct") return advancement_pct;
        if (varName == "M_min") return M_min;
        if (varName == "M_max") return M_max;
        std::cerr << "Unknown variable '" << varName << "'.\n";
        return 0.0;
    }

    // Print full assessment report
    void printAssessment(std::ostream& os = std::cout) const {
        os << std::fixed << std::setprecision(4);
        os << "===== UQFF Framework Advancement Assessment =====\n";
        os << "  Unique term types: " << n_unique_term_types << " / " << n_total_possible << "\n";
        os << "  Dynamic terms:     " << n_dynamic_terms << " / " << n_total_terms << "\n";
        os << "  Mass range:        " << M_min << " kg to " << M_max << " kg\n";
        os << "\n  Scores:\n";
        os << "    Diversity:    " << diversity_score * 100.0 << " %\n";
        os << "    Dynamic:      " << dynamic_score * 100.0 << " %\n";
        os << "    Scalability:  " << scalability_score * 100.0 << " %\n";
        os << "\n  ADVANCEMENT:  " << advancement_pct << " %\n";
        os << "================================================\n";

        if (!module_records.empty()) {
            os << "\n  Modules assessed: " << module_records.size() << "\n";
            for (const auto& r : module_records) {
                os << "    " << r.name << ": " << r.term_count << " terms, "
                   << r.dynamic_count << " dynamic, M=" << r.mass_kg << " kg\n";
            }
        }
    }

    double getAdvancementPct() const { return advancement_pct; }
    double getDiversityScore() const { return diversity_score; }
    double getDynamicScore() const { return dynamic_score; }
    double getScalabilityScore() const { return scalability_score; }
};

}  // namespace UQFF

#endif  // UQFF_LEARNING_ASSESSMENT_H
"""


def main():
    print("gen_muge_assessment.py — Generating UQFFLearningAssessment.h ...")
    content = get_content()
    with open(OUTPUT_FILE, "w", encoding="utf-8") as f:
        f.write(content)
    print(f"  Written: {OUTPUT_FILE}  ({len(content.splitlines())} lines)")


if __name__ == "__main__":
    main()
