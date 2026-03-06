#!/usr/bin/env python3
"""
Production Pipeline - Unified UQFF Query to Analysis to Export

Integrates APIFetch.py, QCalc.py, and QCalc_stat.py into single production pipeline.
Supports batch processing, LaTeX export, and JSON persistence.

PHASE 6: QCalc automatically includes Phase 6 galaxy physics (M51, NGC1316, SMBH
binaries) when appropriate parameters are detected.

Python Ecosystem Status (March 6, 2026):
- CondensedPhysics2.py: 529 classes (~42,000 lines)
  - Thread 10220801: 10 solar UQFF calibration calculators (Ug1-4 full forms,
    11-year B_s(t) cycle, 2025 EHT M_bh=8.55e36 kg Sgr A*)
  - Thread 9c366646: GrokThreadUQFFExtensions imported (14 classes)
  - Thread 3a469fcc: 8 canonical UQFF calculators (Star Magic 14Apr2025 derivation)
    ReactorEfficiencyUQFFCanonical, FUPiNegativeTimeCanonical, QuasarJetNavierStokes,
    PlanetaryCoreHamiltonian, StellarAgeHelioCorrelation, DifferentialRotationDisk,
    SCmDipoleAmplified, YangMillsMassGap
  - Thread ff01cb3a: 5 full-reconstruction calculators (Star Magic 14Apr2025)
    SCmDerivativeHierarchy, Ug2SolarWindTransmutation, Ug4GalacticNonInteractive,
    SolarCycleCoupledFU, FrozenPlanetSolarWind
  - Thread f3c55f52: 5 vacuum-mediated calculators (Superconductivity Unifies Quantum & Gravity 09Sept2025)
    Ug4VacuumMediated, AGNFeedbackFactor, InflationEpochStructure,
    DiPseudoMonopoleOrigin, VacuumEnergyComponentDensity
- GrokThreadUQFFExtensions.py: 2,229 lines, 14 classes, GROK_THREAD_UQFF_CALCULATORS
  registry (13-term g_res, AsymCap, FractalTime, Monte Carlo, 17 buoyancy proofs)
- CondensedPhysicsAggregator.py: v1.2.0, 9 modules, ALL_CALCULATORS unified dict
- Whitepapers: 19 papers (PAPER_001-018 + PAPER_UQFF_VacuumEnergy; GW physics,
  quantum entanglement, LISA noise spectrum, vacuum/dark energy connection)

Routing note: Queries matching Solar Cycle / 11-year / GW BNS / SGWB / Magnetar /
ReactorEfficiency / QuasarJet / YangMills keywords are routed to CP2 (trigger-based, qcalc_cp2_hybrid.py).

Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
Copyright: © 2025-2026 Daniel T. Murphy - All Rights Reserved
Created: February 13, 2026
Updated: March 5, 2026
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Union
from datetime import datetime
import warnings

# Import core modules
from APIFetch import fetch, fetch_and_save
from QCalc import UnifiedFieldSolver, ComputeParams, CONSTANTS
from QCalc_stat import (
    compute_triple_point_analysis,
    print_analysis_summary,
    TriplePointAnalysis
)
from OPData import store as save_query_result


class ProductionPipeline:
    """
    Production-ready pipeline for UQFF analysis workflow.
    
    Workflow:
        1. Query → APIFetch.py (SIMBAD/NASA/Grok)
        2. Parameters → QCalc.py (Physics computation)
        3. Results → QCalc_stat.py (Statistical analysis)
        4. Export → JSON/LaTeX/CSV
    
    Example:
        >>> pipeline = ProductionPipeline()
        >>> result = pipeline.run_query("Sagittarius A*")
        >>> pipeline.export_latex(result, "output.tex")
    """
    
    def __init__(self, output_dir: str = "production_output"):
        """
        Initialize production pipeline.
        
        Args:
            output_dir: Directory for output files (JSON, LaTeX, CSV)
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.solver = UnifiedFieldSolver()
        self.results_cache: Dict[str, Dict] = {}
        
    def run_query(
        self,
        query_name: str,
        fetch_api: bool = True,
        unit_filter: Optional[str] = None,
        save_to_opdata: bool = True
    ) -> Dict:
        """
        Run complete pipeline for single astrophysical query.
        
        Args:
            query_name: Object name (e.g., "Sagittarius A*", "M87", "Betelgeuse")
            fetch_api: If True, fetch from API; if False, use existing bodies.csv
            unit_filter: Filter statistical analysis by unit (e.g., 'm/s²')
            save_to_opdata: Save result to OPData.py for recall
        
        Returns:
            Dictionary with:
                - query_name: Input query
                - api_data: Fetched parameters (M, r, B, z, etc.)
                - qcalc_result: Physics equations and results
                - statistical_analysis: Triple point analysis
                - timestamp: ISO timestamp
        
        Example:
            >>> result = pipeline.run_query("SGR 0501+4516", unit_filter='m/s²')
            >>> print(f"Analyzed {result['statistical_analysis'].num_equations} equations")
        """
        print(f"\n{'='*80}")
        print(f"PRODUCTION PIPELINE: {query_name}")
        print(f"{'='*80}\n")
        
        timestamp = datetime.now().isoformat()
        
        # Step 1: Fetch astrophysical parameters
        print("[1/4] Fetching astrophysical parameters...")
        if fetch_api:
            api_data = fetch(query_name)
            if not api_data or 'error' in api_data:
                raise ValueError(f"API fetch failed for '{query_name}': {api_data.get('error', 'Unknown error')}")
        else:
            # Use existing bodies.csv
            api_data = self._load_from_csv(query_name)
        
        print(f"    ✓ Fetched: M={api_data.get('M', 'N/A')}, r={api_data.get('r', 'N/A')}, B={api_data.get('B', 'N/A')}")
        
        # Step 2: Convert to ComputeParams and run QCalc solver
        print("[2/4] Running UQFF physics calculations...")
        params = self._api_to_compute_params(query_name, api_data)
        qcalc_result = self.solver.solve(params)
        
        num_equations = len(qcalc_result.get('long_form_equations', []))
        print(f"    ✓ Computed {num_equations} physics equations")
        
        # Step 3: Statistical analysis
        print("[3/4] Computing triple point statistical analysis...")
        try:
            # Remove unit filter to analyze all equations
            analysis = compute_triple_point_analysis(
                qcalc_result,
                unit_filter=None  # Analyze all equations regardless of unit
            )
            print(f"    ✓ Range: {analysis.range.min_value:.3e} to {analysis.range.max_value:.3e}")
            print(f"    ✓ Scale: 10^{analysis.scale.order_of_magnitude}")
            print(f"    ✓ Distribution: {analysis.probability.distribution_type} (R²={analysis.probability.fit_quality:.3f})")
        except Exception as e:
            warnings.warn(f"Statistical analysis failed: {e}")
            analysis = None
        
        # Step 4: Package results
        print("[4/4] Packaging results...")
        result = {
            'query_name': query_name,
            'timestamp': timestamp,
            'api_data': api_data,
            'qcalc_result': qcalc_result,
            'statistical_analysis': analysis.to_dict() if analysis else None,
            'metadata': {
                'num_equations': num_equations,
                'unit_filter': unit_filter,
                'solver_version': 'QCalc v1.0',
                'pipeline_version': 'Production v1.1',
                'cp2_classes': 548,
                'grok_thread_classes': 14,
                'aggregator_version': '1.2.0',
                'git_head': 'TBD',
                'thread_1a2726a4': 'UQFF Full Document Assimilation: Q_wave_47-81 Stats, H2O-H2 Rotor CS, DPM-THz MUGE, BEC Alpha-Clustering, Superconductive Complex Ui (5 classes; IPC 0x0A00-0x0A04)'
            }
        }
        
        # Cache result
        self.results_cache[query_name] = result
        
        # Save to OPData
        if save_to_opdata:
            # OPData.store() expects a single dict with query_name inside
            save_query_result({
                'query_name': query_name,
                **qcalc_result
            })
            print(f"    ✓ Saved to OPData.py")
        
        print(f"\n{'='*80}")
        print(f"PIPELINE COMPLETE: {query_name}")
        print(f"{'='*80}\n")
        
        return result
    
    def run_batch(
        self,
        query_list: List[str],
        unit_filter: Optional[str] = None,
        continue_on_error: bool = True
    ) -> Dict[str, Dict]:
        """
        Run pipeline for multiple queries in batch.
        
        Args:
            query_list: List of object names
            unit_filter: Optional unit filter for all queries
            continue_on_error: If True, continue batch even if one query fails
        
        Returns:
            Dictionary mapping query_name → result
        
        Example:
            >>> queries = ["Sgr A*", "M87", "NGC 1365", "Vela Pulsar"]
            >>> batch_results = pipeline.run_batch(queries, unit_filter='m/s²')
            >>> print(f"Processed {len(batch_results)} objects")
        """
        print(f"\n{'='*80}")
        print(f"BATCH PROCESSING: {len(query_list)} queries")
        print(f"{'='*80}\n")
        
        results = {}
        for i, query_name in enumerate(query_list, 1):
            print(f"\n[{i}/{len(query_list)}] Processing: {query_name}")
            try:
                result = self.run_query(query_name, unit_filter=unit_filter)
                results[query_name] = result
                print(f"    ✓ SUCCESS")
            except Exception as e:
                print(f"    ✗ FAILED: {e}")
                if not continue_on_error:
                    raise
                results[query_name] = {'error': str(e), 'query_name': query_name}
        
        print(f"\n{'='*80}")
        print(f"BATCH COMPLETE: {len(results)}/{len(query_list)} successful")
        print(f"{'='*80}\n")
        
        return results
    
    def export_json(self, result: Dict, filename: Optional[str] = None) -> Path:
        """
        Export result to JSON file.
        
        Args:
            result: Pipeline result dictionary
            filename: Output filename (default: <query_name>_<timestamp>.json)
        
        Returns:
            Path to saved JSON file
        """
        if filename is None:
            query_name = result['query_name'].replace(' ', '_').replace('*', 'star')
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            filename = f"{query_name}_{timestamp}.json"
        
        output_path = self.output_dir / filename
        
        # Custom JSON encoder for enum and other non-serializable types
        def json_serializer(obj):
            if hasattr(obj, '__class__') and hasattr(obj.__class__, '__name__'):
                if 'UQFFScale' in obj.__class__.__name__:
                    return str(obj.name) if hasattr(obj, 'name') else str(obj)
                if hasattr(obj, 'value'):
                    return obj.value
            return str(obj)
        
        with open(output_path, 'w') as f:
            json.dump(result, f, indent=2, default=json_serializer)
        
        print(f"✓ Exported JSON: {output_path}")
        return output_path
    
    def export_latex_table(
        self,
        result: Dict,
        filename: Optional[str] = None,
        include_stats: bool = True
    ) -> Path:
        """
        Export result as LaTeX table for manuscript.
        
        Args:
            result: Pipeline result dictionary
            filename: Output filename (default: <query_name>_table.tex)
            include_stats: Include statistical analysis in table
        
        Returns:
            Path to saved LaTeX file
        
        Example LaTeX output:
            \\begin{table}[ht]
            \\centering
            \\caption{UQFF Analysis: Sagittarius A*}
            \\begin{tabular}{lll}
            \\toprule
            Equation & Result & Unit \\\\
            \\midrule
            $F_U$ & $1.234 \\times 10^{10}$ & m/s² \\\\
            ...
        """
        if filename is None:
            query_name = result['query_name'].replace(' ', '_').replace('*', 'star')
            filename = f"{query_name}_table.tex"
        
        output_path = self.output_dir / filename
        
        with open(output_path, 'w') as f:
            f.write("% UQFF Analysis Table\n")
            f.write(f"% Generated: {datetime.now().isoformat()}\n")
            f.write(f"% Query: {result['query_name']}\n\n")
            
            f.write("\\begin{table}[ht]\n")
            f.write("\\centering\n")
            f.write(f"\\caption{{UQFF Analysis: {result['query_name']}}}\n")
            f.write("\\label{tab:" + result['query_name'].replace(' ', '_').lower() + "}\n")
            f.write("\\begin{tabular}{lll}\n")
            f.write("\\toprule\n")
            f.write("Equation & Result & Unit \\\\\n")
            f.write("\\midrule\n")
            
            # Add equations (limit to first 20 for readability)
            equations = result['qcalc_result'].get('long_form_equations', [])[:20]
            for eq in equations:
                name = eq['name'].replace('_', '\\_')
                result_val = f"{eq['result']:.3e}"
                unit = eq.get('unit', '')
                f.write(f"${name}$ & ${result_val}$ & {unit} \\\\\n")
            
            if len(result['qcalc_result'].get('long_form_equations', [])) > 20:
                f.write("\\multicolumn{3}{c}{...} \\\\\n")
            
            f.write("\\bottomrule\n")
            f.write("\\end{tabular}\n")
            f.write("\\end{table}\n")
            
            # Add statistical summary
            if include_stats and result.get('statistical_analysis'):
                f.write("\n% Statistical Summary\n")
                f.write("% Range: {:.3e} to {:.3e}\n".format(
                    result['statistical_analysis']['range']['min_value'],
                    result['statistical_analysis']['range']['max_value']
                ))
                f.write("% Scale: 10^{}\n".format(
                    result['statistical_analysis']['scale']['order_of_magnitude']
                ))
                f.write("% Distribution: {}\n".format(
                    result['statistical_analysis']['probability']['distribution_type']
                ))
        
        print(f"✓ Exported LaTeX table: {output_path}")
        return output_path
    
    def export_csv(self, result: Dict, filename: Optional[str] = None) -> Path:
        """
        Export equations to CSV file.
        
        Args:
            result: Pipeline result dictionary
            filename: Output filename (default: <query_name>_equations.csv)
        
        Returns:
            Path to saved CSV file
        """
        if filename is None:
            query_name = result['query_name'].replace(' ', '_').replace('*', 'star')
            filename = f"{query_name}_equations.csv"
        
        output_path = self.output_dir / filename
        
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write("Equation,Result,Unit,Description\n")
            equations = result['qcalc_result'].get('long_form_equations', [])
            for eq in equations:
                name = eq['name']
                result_val = eq['result']
                unit = eq.get('unit', '')
                desc = eq.get('description', '')
                f.write(f"{name},{result_val},{unit},{desc}\n")
        
        print(f"✓ Exported CSV: {output_path}")
        return output_path
    
    def _api_to_compute_params(self, query_name: str, api_data: Dict) -> ComputeParams:
        """Convert API data to ComputeParams object."""
        return ComputeParams(
            query_name=query_name,
            M=api_data.get('M', 1e30),
            r=api_data.get('r', 1e6),
            B=api_data.get('B', 1e-5),
            t=api_data.get('t', 1e8),
            z=api_data.get('z', 0.0),
            T=api_data.get('T', None),
            L=api_data.get('L', None),
            R=api_data.get('R', None),
            d=api_data.get('d', None),
            v=api_data.get('v', None),
            omega=api_data.get('omega', None),
            P=api_data.get('P', None),
            t_n=api_data.get('tn_years', 0.0)  # Map tn_years → t_n
        )
    
    def _load_from_csv(self, query_name: str) -> Dict:
        """Load parameters from bodies.csv."""
        csv_path = Path("bodies.csv")
        if not csv_path.exists():
            raise FileNotFoundError("bodies.csv not found")
        
        # Simple CSV parsing (assumes first line is header)
        with open(csv_path) as f:
            lines = f.readlines()
        
        header = lines[0].strip().split(',')
        for line in lines[1:]:
            values = line.strip().split(',')
            if values[0] == query_name:
                return dict(zip(header, values))
        
        raise ValueError(f"Query '{query_name}' not found in bodies.csv")


def main():
    """
    CLI interface for production pipeline.
    
    Usage:
        python production_pipeline.py "Sagittarius A*"
        python production_pipeline.py "M87" "NGC 1365" "Vela Pulsar"
    """
    if len(sys.argv) < 2:
        print("Usage: python production_pipeline.py <query1> [query2] [query3] ...")
        print("\nExample:")
        print('  python production_pipeline.py "Sagittarius A*"')
        print('  python production_pipeline.py "M87" "NGC 1365" "Vela Pulsar"')
        sys.exit(1)
    
    pipeline = ProductionPipeline()
    queries = sys.argv[1:]
    
    if len(queries) == 1:
        # Single query
        query_name = queries[0]
        result = pipeline.run_query(query_name, unit_filter='m/s²')
        
        # Export all formats
        pipeline.export_json(result)
        pipeline.export_latex_table(result)
        pipeline.export_csv(result)
        
        # Print summary
        if result.get('statistical_analysis'):
            print("\nStatistical Summary:")
            print_analysis_summary(
                TriplePointAnalysis(**result['statistical_analysis'])
            )
    else:
        # Batch processing
        batch_results = pipeline.run_batch(queries, unit_filter='m/s²')
        
        # Export each
        for query_name, result in batch_results.items():
            if 'error' not in result:
                pipeline.export_json(result)
                pipeline.export_latex_table(result)
                pipeline.export_csv(result)
        
        print(f"\nBatch processing complete: {len(batch_results)} queries")


if __name__ == '__main__':
    main()
