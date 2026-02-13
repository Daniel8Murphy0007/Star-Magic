#!/usr/bin/env python3
"""
Batch Runner - Process multiple queries from JSON configuration

Reads manuscript_batch_config.json and runs production pipeline for all queries.
Generates organized output for manuscript figures and tables.

Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
Created: February 13, 2026
"""

import json
import sys
from pathlib import Path
from datetime import datetime
from typing import Dict, List
import traceback

from production_pipeline import ProductionPipeline


class BatchRunner:
    """
    Batch processor for manuscript queries.
    
    Example:
        >>> runner = BatchRunner("manuscript_batch_config.json")
        >>> results = runner.run_all()
        >>> runner.generate_summary_report()
    """
    
    def __init__(self, config_file: str = "manuscript_batch_config.json"):
        """
        Initialize batch runner.
        
        Args:
            config_file: Path to JSON configuration file
        """
        self.config_file = Path(config_file)
        if not self.config_file.exists():
            raise FileNotFoundError(f"Config file not found: {config_file}")
        
        with open(self.config_file) as f:
            self.config = json.load(f)
        
        # Create timestamped output directory
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        batch_name = self.config.get('batch_name', 'batch').replace(' ', '_')
        output_dir = f"production_output/{batch_name}_{timestamp}"
        
        self.pipeline = ProductionPipeline(output_dir=output_dir)
        self.results: Dict[str, Dict] = {}
        self.errors: Dict[str, str] = {}
        
        print(f"Batch Runner initialized")
        print(f"  Config: {config_file}")
        print(f"  Output: {output_dir}")
        print(f"  Queries: {len(self.config['queries'])}")
    
    def run_all(self, stop_on_error: bool = False) -> Dict[str, Dict]:
        """
        Run pipeline for all queries in configuration.
        
        Args:
            stop_on_error: If True, stop batch on first error
        
        Returns:
            Dictionary mapping query_name → result
        """
        queries = self.config['queries']
        total = len(queries)
        
        print(f"\n{'='*80}")
        print(f"BATCH RUN: {self.config.get('batch_name', 'Unnamed Batch')}")
        print(f"Processing {total} queries...")
        print(f"{'='*80}\n")
        
        for i, query_config in enumerate(queries, 1):
            query_name = query_config['name']
            category = query_config.get('category', 'Unknown')
            unit_filter = query_config.get('unit_filter') or \
                          self.config.get('analysis_options', {}).get('unit_filter')
            
            print(f"\n[{i}/{total}] {query_name} ({category})")
            print("-" * 80)
            
            try:
                # Run pipeline
                result = self.pipeline.run_query(
                    query_name=query_name,
                    fetch_api=True,
                    unit_filter=unit_filter,
                    save_to_opdata=self.config.get('analysis_options', {}).get('save_to_opdata', True)
                )
                
                self.results[query_name] = result
                
                # Export all configured formats
                for export_format in self.config.get('export_formats', ['json']):
                    if export_format == 'json':
                        self.pipeline.export_json(result)
                    elif export_format == 'latex':
                        self.pipeline.export_latex_table(result)
                    elif export_format == 'csv':
                        self.pipeline.export_csv(result)
                
                print(f"✓ SUCCESS: {query_name}")
                
            except Exception as e:
                error_msg = f"{type(e).__name__}: {str(e)}"
                self.errors[query_name] = error_msg
                print(f"✗ FAILED: {query_name}")
                print(f"  Error: {error_msg}")
                
                if stop_on_error:
                    print("\nStopping batch due to error (stop_on_error=True)")
                    break
                else:
                    # Log full traceback for debugging
                    print(f"\n  Full traceback:")
                    traceback.print_exc()
                    print()
        
        print(f"\n{'='*80}")
        print(f"BATCH COMPLETE")
        print(f"  Successful: {len(self.results)}/{total}")
        print(f"  Failed: {len(self.errors)}/{total}")
        print(f"{'='*80}\n")
        
        return self.results
    
    def run_priority(self, priority: int) -> Dict[str, Dict]:
        """
        Run only queries matching specified priority.
        
        Args:
            priority: Priority level (1=highest, 2=medium, 3=low)
        
        Returns:
            Dictionary of results for matching queries
        """
        matching_queries = [
            q for q in self.config['queries']
            if q.get('priority', 999) == priority
        ]
        
        if not matching_queries:
            print(f"No queries found with priority {priority}")
            return {}
        
        print(f"Running {len(matching_queries)} priority-{priority} queries...")
        
        # Temporarily modify config
        original_queries = self.config['queries']
        self.config['queries'] = matching_queries
        
        results = self.run_all(stop_on_error=False)
        
        # Restore original config
        self.config['queries'] = original_queries
        
        return results
    
    def generate_summary_report(self, output_file: str = "batch_summary_report.md"):
        """
        Generate markdown summary report of batch results.
        
        Args:
            output_file: Output markdown filename
        """
        output_path = self.pipeline.output_dir / output_file
        
        with open(output_path, 'w') as f:
            f.write(f"# Batch Processing Summary Report\n\n")
            f.write(f"**Batch:** {self.config.get('batch_name', 'Unnamed')}\n")
            f.write(f"**Date:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"**Total Queries:** {len(self.config['queries'])}\n")
            f.write(f"**Successful:** {len(self.results)}\n")
            f.write(f"**Failed:** {len(self.errors)}\n\n")
            
            f.write("---\n\n")
            
            # Success summary
            if self.results:
                f.write("## Successful Queries\n\n")
                f.write("| Query | Category | Equations | Distribution | Scale |\n")
                f.write("|-------|----------|-----------|--------------|-------|\n")
                
                for query_name, result in self.results.items():
                    category = next(
                        (q['category'] for q in self.config['queries'] if q['name'] == query_name),
                        'Unknown'
                    )
                    num_eq = result['metadata']['num_equations']
                    
                    if result.get('statistical_analysis'):
                        dist = result['statistical_analysis']['probability']['distribution_type']
                        scale = f"10^{result['statistical_analysis']['scale']['order_of_magnitude']}"
                    else:
                        dist = "N/A"
                        scale = "N/A"
                    
                    f.write(f"| {query_name} | {category} | {num_eq} | {dist} | {scale} |\n")
                
                f.write("\n")
            
            # Error summary
            if self.errors:
                f.write("## Failed Queries\n\n")
                f.write("| Query | Error |\n")
                f.write("|-------|-------|\n")
                
                for query_name, error in self.errors.items():
                    f.write(f"| {query_name} | {error} |\n")
                
                f.write("\n")
            
            # Detailed results
            f.write("## Detailed Results\n\n")
            
            for query_name, result in self.results.items():
                f.write(f"### {query_name}\n\n")
                
                # Metadata
                f.write(f"- **Timestamp:** {result['timestamp']}\n")
                f.write(f"- **Equations:** {result['metadata']['num_equations']}\n")
                
                # API data
                api_data = result['api_data']
                f.write(f"- **Mass:** {api_data.get('M', 'N/A')}\n")
                f.write(f"- **Radius:** {api_data.get('r', 'N/A')}\n")
                f.write(f"- **Magnetic Field:** {api_data.get('B', 'N/A')}\n")
                
                # Statistical analysis
                if result.get('statistical_analysis'):
                    stats = result['statistical_analysis']
                    f.write(f"\n**Statistical Analysis:**\n")
                    f.write(f"- Range: {stats['range']['min_value']:.3e} to {stats['range']['max_value']:.3e}\n")
                    f.write(f"- Scale: 10^{stats['scale']['order_of_magnitude']}\n")
                    f.write(f"- Distribution: {stats['probability']['distribution_type']}\n")
                    f.write(f"- Fit Quality: R²={stats['probability']['fit_quality']:.3f}\n")
                    
                    if stats['probability']['outliers']:
                        f.write(f"- Outliers: {', '.join(stats['probability']['outliers'])}\n")
                
                f.write("\n")
        
        print(f"✓ Summary report saved: {output_path}")
    
    def generate_latex_summary_table(self, output_file: str = "batch_summary_table.tex"):
        """
        Generate LaTeX summary table for manuscript.
        
        Args:
            output_file: Output LaTeX filename
        """
        output_path = self.pipeline.output_dir / output_file
        
        with open(output_path, 'w') as f:
            f.write("% Batch Processing Summary Table\n")
            f.write(f"% Generated: {datetime.now().isoformat()}\n\n")
            
            f.write("\\begin{table*}[ht]\n")
            f.write("\\centering\n")
            f.write("\\caption{UQFF Analysis Summary for " + str(len(self.results)) + " Astrophysical Systems}\n")
            f.write("\\label{tab:uqff_batch_summary}\n")
            f.write("\\begin{tabular}{llcccc}\n")
            f.write("\\toprule\n")
            f.write("Object & Category & $N_{\\mathrm{eq}}$ & Distribution & Scale & $R^2$ \\\\\n")
            f.write("\\midrule\n")
            
            for query_name, result in self.results.items():
                category = next(
                    (q['category'] for q in self.config['queries'] if q['name'] == query_name),
                    'Unknown'
                )
                
                # Escape LaTeX special characters
                safe_name = query_name.replace('*', '$^*$').replace('_', '\\_')
                
                num_eq = result['metadata']['num_equations']
                
                if result.get('statistical_analysis'):
                    stats = result['statistical_analysis']
                    dist = stats['probability']['distribution_type']
                    scale = f"$10^{{{stats['scale']['order_of_magnitude']}}}$"
                    r2 = f"{stats['probability']['fit_quality']:.3f}"
                else:
                    dist = "---"
                    scale = "---"
                    r2 = "---"
                
                f.write(f"{safe_name} & {category} & {num_eq} & {dist} & {scale} & {r2} \\\\\n")
            
            f.write("\\bottomrule\n")
            f.write("\\end{tabular}\n")
            f.write("\\end{table*}\n")
        
        print(f"✓ LaTeX summary table saved: {output_path}")


def main():
    """
    CLI interface for batch runner.
    
    Usage:
        python batch_runner.py                    # Run all queries
        python batch_runner.py --priority 1       # Run priority 1 only
        python batch_runner.py --config custom.json
    """
    import argparse
    
    parser = argparse.ArgumentParser(description='Run batch processing from JSON config')
    parser.add_argument('--config', default='manuscript_batch_config.json',
                        help='Configuration file (default: manuscript_batch_config.json)')
    parser.add_argument('--priority', type=int, help='Run only specified priority level')
    parser.add_argument('--stop-on-error', action='store_true',
                        help='Stop batch on first error')
    
    args = parser.parse_args()
    
    try:
        runner = BatchRunner(config_file=args.config)
        
        if args.priority:
            results = runner.run_priority(args.priority)
        else:
            results = runner.run_all(stop_on_error=args.stop_on_error)
        
        # Generate reports
        runner.generate_summary_report()
        runner.generate_latex_summary_table()
        
        print(f"\n✓ Batch processing complete!")
        print(f"  Results: {len(results)}")
        print(f"  Output directory: {runner.pipeline.output_dir}")
        
    except Exception as e:
        print(f"\n✗ Batch runner failed: {e}")
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
