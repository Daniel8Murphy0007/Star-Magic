#!/usr/bin/env python3
"""
QCalc_stat.py - Statistical Analysis Module for UQFF Query Results
====================================================================

Comprehensive statistical analysis of QCalc.py equation results with triple point
analysis: Range, Scale, and Probability distribution fitting.

Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
Framework: UQFF 99.9% Solvability (Star-Magic)
Copyright: © 2025-2026 Daniel T. Murphy - All Rights Reserved
Created: February 13, 2026
"""

import numpy as np
import json
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass, field
from scipy import stats
from scipy.optimize import curve_fit
from pathlib import Path
import warnings

# Import data layer for accessing query results
from OPData import OutputDataStore, QUERY_RESULTS


# ═══════════════════════════════════════════════════════════════════════════════
# DATA STRUCTURES
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class RangeAnalysis:
    """Range analysis (min/max) for equation results."""
    min_value: float
    max_value: float
    range_span: float  # max - min
    dynamic_range_db: float  # 20*log10(max/min) if min > 0
    min_equation: str  # Name of equation with minimum value
    max_equation: str  # Name of equation with maximum value
    unit: str
    
    def to_dict(self) -> Dict:
        return {
            'min_value': self.min_value,
            'max_value': self.max_value,
            'range_span': self.range_span,
            'dynamic_range_db': self.dynamic_range_db,
            'min_equation': self.min_equation,
            'max_equation': self.max_equation,
            'unit': self.unit
        }


@dataclass
class ScaleAnalysis:
    """Scale analysis (order of magnitude) for equation results."""
    characteristic_scale: float  # Geometric mean
    order_of_magnitude: int  # Floor(log10(characteristic_scale))
    scale_factor: float  # log10(max/min) if min > 0
    median_value: float
    geometric_mean: float
    harmonic_mean: float  # If all values positive
    
    def to_dict(self) -> Dict:
        return {
            'characteristic_scale': self.characteristic_scale,
            'order_of_magnitude': self.order_of_magnitude,
            'scale_factor': self.scale_factor,
            'median_value': self.median_value,
            'geometric_mean': self.geometric_mean,
            'harmonic_mean': self.harmonic_mean
        }


@dataclass
class ProbabilityAnalysis:
    """Probability distribution analysis for equation results."""
    distribution_type: str  # 'gaussian', 'lognormal', 'uniform', 'exponential', 'power_law'
    mean: float
    std_dev: float
    skewness: float
    kurtosis: float
    fit_quality: float  # R² or Kolmogorov-Smirnov statistic
    fit_parameters: Dict[str, float]  # Distribution-specific parameters
    confidence_interval_95: Tuple[float, float]
    outliers: List[str]  # Names of equations with outlier values (>3σ)
    
    def to_dict(self) -> Dict:
        return {
            'distribution_type': self.distribution_type,
            'mean': self.mean,
            'std_dev': self.std_dev,
            'skewness': self.skewness,
            'kurtosis': self.kurtosis,
            'fit_quality': self.fit_quality,
            'fit_parameters': self.fit_parameters,
            'confidence_interval_95': self.confidence_interval_95,
            'outliers': self.outliers
        }


@dataclass
class FineFitRatio:
    """Fine fit ratio - precision metric for observed vs predicted values."""
    ratio_mean: float  # Mean of (observed / predicted)
    ratio_std: float  # Std dev of ratio
    percent_error_mean: float  # Mean abs percent error
    percent_error_std: float  # Std dev of percent error
    max_deviation: float  # Max abs(observed - predicted)
    max_deviation_equation: str  # Name of equation with max deviation
    chi_squared: float  # χ² goodness of fit
    reduced_chi_squared: float  # χ²/dof
    
    def to_dict(self) -> Dict:
        return {
            'ratio_mean': self.ratio_mean,
            'ratio_std': self.ratio_std,
            'percent_error_mean': self.percent_error_mean,
            'percent_error_std': self.percent_error_std,
            'max_deviation': self.max_deviation,
            'max_deviation_equation': self.max_deviation_equation,
            'chi_squared': self.chi_squared,
            'reduced_chi_squared': self.reduced_chi_squared
        }


@dataclass
class TriplePointAnalysis:
    """Complete triple point (range, scale, probability) statistical analysis."""
    query_name: str
    num_equations: int
    range: RangeAnalysis
    scale: ScaleAnalysis
    probability: ProbabilityAnalysis
    fine_fit: Optional[FineFitRatio] = None  # Requires reference/predicted values
    correlation_matrix: Optional[Dict[str, Dict[str, float]]] = None  # Cross-equation correlations
    
    def to_dict(self) -> Dict:
        result = {
            'query_name': self.query_name,
            'num_equations': self.num_equations,
            'range': self.range.to_dict(),
            'scale': self.scale.to_dict(),
            'probability': self.probability.to_dict()
        }
        if self.fine_fit:
            result['fine_fit'] = self.fine_fit.to_dict()
        if self.correlation_matrix:
            result['correlation_matrix'] = self.correlation_matrix
        return result
    
    @classmethod
    def from_dict(cls, d: Dict) -> 'TriplePointAnalysis':
        """Reconstruct a TriplePointAnalysis from a to_dict() result."""
        obj = cls.__new__(cls)
        obj.query_name = d['query_name']
        obj.num_equations = d['num_equations']
        obj.range = RangeAnalysis(**d['range'])
        obj.scale = ScaleAnalysis(**d['scale'])
        obj.probability = ProbabilityAnalysis(**d['probability'])
        obj.fine_fit = FineFitRatio(**d['fine_fit']) if d.get('fine_fit') else None
        obj.correlation_matrix = d.get('correlation_matrix')
        return obj

    def to_json(self, indent: int = 2) -> str:
        """Convert to JSON string."""
        return json.dumps(self.to_dict(), indent=indent)
    
    def save(self, filepath: str):
        """Save analysis to JSON file."""
        with open(filepath, 'w') as f:
            f.write(self.to_json())


# ═══════════════════════════════════════════════════════════════════════════════
# STATISTICAL ANALYSIS FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════

def compute_range_analysis(
    equations: List[Dict],
    unit_filter: Optional[str] = None
) -> RangeAnalysis:
    """
    Compute range analysis (min/max) for equation results.
    
    Args:
        equations: List of equation dicts with 'name', 'result', 'unit'
        unit_filter: Optional unit filter (e.g., 'm/s²') to analyze only matching equations
        
    Returns:
        RangeAnalysis object with min/max statistics
    """
    if unit_filter:
        equations = [eq for eq in equations if eq.get('unit') == unit_filter]
    
    if not equations:
        raise ValueError("No equations match the filter criteria")
    
    # Extract results (skip non-numeric or None values)
    results = []
    names = []
    for eq in equations:
        result = eq.get('result')
        if isinstance(result, (int, float, np.number)) and not np.isnan(result) and not np.isinf(result):
            results.append(float(result))
            names.append(eq['name'])
    
    if not results:
        raise ValueError("No valid numeric results found")
    
    results = np.array(results)
    
    min_val = np.min(results)
    max_val = np.max(results)
    min_idx = np.argmin(results)
    max_idx = np.argmax(results)
    
    range_span = max_val - min_val
    
    # Dynamic range in dB (20*log10(max/min)) if min > 0
    if min_val > 0:
        dynamic_range_db = 20 * np.log10(max_val / min_val)
    else:
        dynamic_range_db = float('inf')
    
    unit = equations[0].get('unit', '(dimensionless)')
    
    return RangeAnalysis(
        min_value=min_val,
        max_value=max_val,
        range_span=range_span,
        dynamic_range_db=dynamic_range_db,
        min_equation=names[min_idx],
        max_equation=names[max_idx],
        unit=unit
    )


def compute_scale_analysis(
    equations: List[Dict],
    unit_filter: Optional[str] = None
) -> ScaleAnalysis:
    """
    Compute scale analysis (order of magnitude) for equation results.
    
    Args:
        equations: List of equation dicts with 'name', 'result', 'unit'
        unit_filter: Optional unit filter
        
    Returns:
        ScaleAnalysis object with characteristic scales
    """
    if unit_filter:
        equations = [eq for eq in equations if eq.get('unit') == unit_filter]
    
    if not equations:
        raise ValueError("No equations match the filter criteria")
    
    # Extract positive results only (for geometric/harmonic means)
    results = []
    for eq in equations:
        result = eq.get('result')
        if isinstance(result, (int, float, np.number)) and not np.isnan(result) and not np.isinf(result):
            results.append(float(result))
    
    if not results:
        raise ValueError("No valid numeric results found")
    
    results = np.array(results)
    results_positive = np.abs(results[results != 0])  # Abs for log operations
    
    # Characteristic scale = geometric mean of absolute values
    if len(results_positive) > 0:
        geometric_mean = np.exp(np.mean(np.log(results_positive)))
        characteristic_scale = geometric_mean
        order_of_magnitude = int(np.floor(np.log10(characteristic_scale)))
        
        # Scale factor = log10(max/min) for positive values
        if len(results_positive) > 1:
            scale_factor = np.log10(np.max(results_positive) / np.min(results_positive))
        else:
            scale_factor = 0.0
        
        # Harmonic mean (only for all-positive values)
        if np.all(results > 0):
            harmonic_mean = len(results) / np.sum(1.0 / results)
        else:
            harmonic_mean = float('nan')
    else:
        geometric_mean = 0.0
        characteristic_scale = 0.0
        order_of_magnitude = 0
        scale_factor = 0.0
        harmonic_mean = 0.0
    
    median_value = np.median(results)
    
    return ScaleAnalysis(
        characteristic_scale=characteristic_scale,
        order_of_magnitude=order_of_magnitude,
        scale_factor=scale_factor,
        median_value=median_value,
        geometric_mean=geometric_mean,
        harmonic_mean=harmonic_mean
    )


def compute_probability_analysis(
    equations: List[Dict],
    unit_filter: Optional[str] = None,
    outlier_threshold: float = 3.0
) -> ProbabilityAnalysis:
    """
    Compute probability distribution analysis for equation results.
    
    Fits multiple distributions (Gaussian, log-normal, exponential, power-law) and
    selects the best fit based on Kolmogorov-Smirnov test.
    
    Args:
        equations: List of equation dicts with 'name', 'result', 'unit'
        unit_filter: Optional unit filter
        outlier_threshold: Number of std deviations for outlier detection (default: 3.0)
        
    Returns:
        ProbabilityAnalysis object with distribution fit and statistics
    """
    if unit_filter:
        equations = [eq for eq in equations if eq.get('unit') == unit_filter]
    
    if not equations:
        raise ValueError("No equations match the filter criteria")
    
    # Extract results and names
    results = []
    names = []
    for eq in equations:
        result = eq.get('result')
        if isinstance(result, (int, float, np.number)) and not np.isnan(result) and not np.isinf(result):
            results.append(float(result))
            names.append(eq['name'])
    
    if len(results) < 3:
        raise ValueError("Need at least 3 valid results for distribution fitting")
    
    results = np.array(results)
    
    # Basic statistics
    mean = np.mean(results)
    std_dev = np.std(results, ddof=1)
    skewness = stats.skew(results)
    kurtosis = stats.kurtosis(results)
    
    # 95% confidence interval (normal approximation)
    ci_95 = stats.norm.interval(0.95, loc=mean, scale=std_dev / np.sqrt(len(results)))
    
    # Outlier detection (>3σ from mean)
    z_scores = np.abs((results - mean) / std_dev) if std_dev > 0 else np.zeros_like(results)
    outlier_mask = z_scores > outlier_threshold
    outliers = [names[i] for i in np.where(outlier_mask)[0]]
    
    # Try fitting multiple distributions
    best_dist = 'gaussian'
    best_ks_stat = float('inf')
    best_params = {}
    
    # 1. Gaussian (normal) distribution
    try:
        ks_stat, _ = stats.kstest(results, 'norm', args=(mean, std_dev))
        if ks_stat < best_ks_stat:
            best_ks_stat = ks_stat
            best_dist = 'gaussian'
            best_params = {'mean': mean, 'std_dev': std_dev}
    except Exception:
        pass
    
    # 2. Log-normal distribution (if all positive)
    if np.all(results > 0):
        try:
            log_results = np.log(results)
            log_mean = np.mean(log_results)
            log_std = np.std(log_results, ddof=1)
            ks_stat, _ = stats.kstest(results, 'lognorm', args=(log_std, 0, np.exp(log_mean)))
            if ks_stat < best_ks_stat:
                best_ks_stat = ks_stat
                best_dist = 'lognormal'
                best_params = {'log_mean': log_mean, 'log_std': log_std}
        except Exception:
            pass
    
    # 3. Exponential distribution (if all positive)
    if np.all(results > 0):
        try:
            scale = np.mean(results)
            ks_stat, _ = stats.kstest(results, 'expon', args=(0, scale))
            if ks_stat < best_ks_stat:
                best_ks_stat = ks_stat
                best_dist = 'exponential'
                best_params = {'scale': scale}
        except Exception:
            pass
    
    # 4. Power-law distribution (if wide range of positive values)
    if np.all(results > 0) and np.max(results) / np.min(results) > 10:
        try:
            # Fit power law: P(x) ∝ x^(-alpha)
            log_results = np.log(results)
            alpha, _ = np.polyfit(np.log(results), np.log(np.arange(1, len(results) + 1)[::-1]), 1)
            alpha = -alpha  # Convention: P(x) ∝ x^(-alpha)
            # Use KS test approximation
            sorted_results = np.sort(results)
            empirical_cdf = np.arange(1, len(results) + 1) / len(results)
            theoretical_cdf = 1 - (sorted_results[0] / sorted_results) ** (alpha - 1)
            ks_stat = np.max(np.abs(empirical_cdf - theoretical_cdf))
            if ks_stat < best_ks_stat:
                best_ks_stat = ks_stat
                best_dist = 'power_law'
                best_params = {'alpha': alpha}
        except Exception:
            pass
    
    # R² quality (1 - KS statistic as approximation)
    fit_quality = max(0.0, 1.0 - best_ks_stat)
    
    return ProbabilityAnalysis(
        distribution_type=best_dist,
        mean=mean,
        std_dev=std_dev,
        skewness=skewness,
        kurtosis=kurtosis,
        fit_quality=fit_quality,
        fit_parameters=best_params,
        confidence_interval_95=ci_95,
        outliers=outliers
    )


def compute_fine_fit_ratio(
    observed_equations: List[Dict],
    predicted_equations: List[Dict]
) -> FineFitRatio:
    """
    Compute fine fit ratio - precision metric comparing observed vs predicted values.
    
    Args:
        observed_equations: List of observed equation results
        predicted_equations: List of predicted/reference equation results
        
    Returns:
        FineFitRatio object with precision metrics
    """
    # Match equations by name
    observed_dict = {eq['name']: eq['result'] for eq in observed_equations 
                     if isinstance(eq.get('result'), (int, float, np.number))}
    predicted_dict = {eq['name']: eq['result'] for eq in predicted_equations 
                      if isinstance(eq.get('result'), (int, float, np.number))}
    
    # Find common equations
    common_names = set(observed_dict.keys()) & set(predicted_dict.keys())
    
    if not common_names:
        raise ValueError("No common equations found between observed and predicted")
    
    ratios = []
    percent_errors = []
    deviations = []
    names = []
    
    for name in common_names:
        obs = observed_dict[name]
        pred = predicted_dict[name]
        
        if pred != 0 and not np.isnan(obs) and not np.isnan(pred):
            ratio = obs / pred
            percent_error = 100 * abs(obs - pred) / abs(pred)
            deviation = abs(obs - pred)
            
            ratios.append(ratio)
            percent_errors.append(percent_error)
            deviations.append(deviation)
            names.append(name)
    
    if not ratios:
        raise ValueError("No valid comparison pairs found")
    
    ratios = np.array(ratios)
    percent_errors = np.array(percent_errors)
    deviations = np.array(deviations)
    
    ratio_mean = np.mean(ratios)
    ratio_std = np.std(ratios, ddof=1)
    percent_error_mean = np.mean(percent_errors)
    percent_error_std = np.std(percent_errors, ddof=1)
    
    max_deviation = np.max(deviations)
    max_deviation_idx = np.argmax(deviations)
    max_deviation_equation = names[max_deviation_idx]
    
    # Chi-squared goodness of fit
    # χ² = Σ((observed - predicted)² / predicted²)
    chi_squared = np.sum((deviations / np.array([abs(predicted_dict[n]) for n in names])) ** 2)
    dof = len(ratios) - 1  # Degrees of freedom
    reduced_chi_squared = chi_squared / dof if dof > 0 else float('inf')
    
    return FineFitRatio(
        ratio_mean=ratio_mean,
        ratio_std=ratio_std,
        percent_error_mean=percent_error_mean,
        percent_error_std=percent_error_std,
        max_deviation=max_deviation,
        max_deviation_equation=max_deviation_equation,
        chi_squared=chi_squared,
        reduced_chi_squared=reduced_chi_squared
    )


def compute_correlation_matrix(
    equations: List[Dict],
    unit_filter: Optional[str] = None
) -> Dict[str, Dict[str, float]]:
    """
    Compute pairwise correlation matrix between equations.
    
    Requires multiple queries (time-series or parameter sweep) to compute correlations.
    For single query, returns None.
    
    Args:
        equations: List of equation dicts (must have multiple samples per equation)
        unit_filter: Optional unit filter
        
    Returns:
        Nested dict: {eq1_name: {eq2_name: correlation_coefficient}}
    """
    # This would require multiple query results to compute correlations
    # For now, return empty dict (placeholder for future multi-query analysis)
    warnings.warn("Correlation matrix requires multiple query results (time-series or parameter sweep)")
    return {}


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN ANALYSIS FUNCTION
# ═══════════════════════════════════════════════════════════════════════════════

def compute_triple_point_analysis(
    query_results: Dict,
    unit_filter: Optional[str] = None,
    predicted_results: Optional[Dict] = None,
    outlier_threshold: float = 3.0
) -> TriplePointAnalysis:
    """
    Compute complete triple point (range, scale, probability) statistical analysis.
    
    Args:
        query_results: Dict from QCalc.solve() or OPData.QUERY_RESULTS
        unit_filter: Optional unit filter (e.g., 'm/s²') to analyze only matching equations
        predicted_results: Optional predicted/reference results for fine fit ratio
        outlier_threshold: Number of std deviations for outlier detection
        
    Returns:
        TriplePointAnalysis object with complete statistics
        
    Example:
        >>> from QCalc import UnifiedFieldSolver, ComputeParams, CONSTANTS
        >>> from QCalc_stat import compute_triple_point_analysis
        >>> 
        >>> solver = UnifiedFieldSolver()
        >>> params = ComputeParams(query_name="Test", M=1e30, r=1e6)
        >>> result = solver.solve(params)
        >>> 
        >>> analysis = compute_triple_point_analysis(result, unit_filter='m/s²')
        >>> print(f"Range: {analysis.range.min_value:.3e} to {analysis.range.max_value:.3e}")
        >>> print(f"Scale: 10^{analysis.scale.order_of_magnitude}")
        >>> print(f"Distribution: {analysis.probability.distribution_type}")
    """
    query_name = query_results.get('query_name', 'Unknown')
    equations = query_results.get('long_form_equations', [])
    
    if not equations:
        raise ValueError("No equations found in query results")
    
    # Filter by unit if specified
    if unit_filter:
        filtered_equations = [eq for eq in equations if eq.get('unit') == unit_filter]
        if not filtered_equations:
            raise ValueError(f"No equations found with unit '{unit_filter}'")
    else:
        filtered_equations = equations
    
    # Compute range analysis
    range_analysis = compute_range_analysis(filtered_equations, unit_filter=None)
    
    # Compute scale analysis
    scale_analysis = compute_scale_analysis(filtered_equations, unit_filter=None)
    
    # Compute probability analysis
    probability_analysis = compute_probability_analysis(
        filtered_equations, 
        unit_filter=None,
        outlier_threshold=outlier_threshold
    )
    
    # Compute fine fit ratio (if predicted results provided)
    fine_fit = None
    if predicted_results:
        try:
            predicted_equations = predicted_results.get('long_form_equations', [])
            fine_fit = compute_fine_fit_ratio(filtered_equations, predicted_equations)
        except Exception as e:
            warnings.warn(f"Could not compute fine fit ratio: {str(e)}")
    
    # Compute correlation matrix (placeholder for multi-query analysis)
    correlation_matrix = None  # Requires time-series data
    
    return TriplePointAnalysis(
        query_name=query_name,
        num_equations=len(filtered_equations),
        range=range_analysis,
        scale=scale_analysis,
        probability=probability_analysis,
        fine_fit=fine_fit,
        correlation_matrix=correlation_matrix
    )


# ═══════════════════════════════════════════════════════════════════════════════
# CONVENIENCE FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════

def analyze_query_by_name(query_name: str, **kwargs) -> TriplePointAnalysis:
    """
    Analyze a query from OPData.QUERY_RESULTS by name.
    
    Args:
        query_name: Name of query to analyze
        **kwargs: Additional arguments passed to compute_triple_point_analysis
        
    Returns:
        TriplePointAnalysis object
    """
    if query_name not in QUERY_RESULTS:
        raise ValueError(f"Query '{query_name}' not found in QUERY_RESULTS")
    
    query_results = QUERY_RESULTS[query_name]
    return compute_triple_point_analysis(query_results, **kwargs)


def analyze_latest_query(**kwargs) -> TriplePointAnalysis:
    """
    Analyze the most recent query from OPData.
    
    Args:
        **kwargs: Additional arguments passed to compute_triple_point_analysis
        
    Returns:
        TriplePointAnalysis object
    """
    store = OutputDataStore()
    latest = store.get_latest_query()
    if not latest:
        raise ValueError("No queries found in OutputDataStore")
    
    return compute_triple_point_analysis(latest, **kwargs)


def print_analysis_summary(analysis: TriplePointAnalysis):
    """
    Print human-readable summary of triple point analysis.
    
    Args:
        analysis: TriplePointAnalysis object
    """
    print("=" * 80)
    print(f"TRIPLE POINT ANALYSIS: {analysis.query_name}")
    print("=" * 80)
    print(f"Number of equations: {analysis.num_equations}")
    print()
    
    print("RANGE ANALYSIS")
    print("-" * 80)
    print(f"  Min: {analysis.range.min_value:.6e} {analysis.range.unit} ({analysis.range.min_equation})")
    print(f"  Max: {analysis.range.max_value:.6e} {analysis.range.unit} ({analysis.range.max_equation})")
    print(f"  Span: {analysis.range.range_span:.6e} {analysis.range.unit}")
    if np.isfinite(analysis.range.dynamic_range_db):
        print(f"  Dynamic Range: {analysis.range.dynamic_range_db:.2f} dB")
    print()
    
    print("SCALE ANALYSIS")
    print("-" * 80)
    print(f"  Characteristic Scale: {analysis.scale.characteristic_scale:.6e}")
    print(f"  Order of Magnitude: 10^{analysis.scale.order_of_magnitude}")
    print(f"  Scale Factor: {analysis.scale.scale_factor:.2f} (log10 range)")
    print(f"  Median: {analysis.scale.median_value:.6e}")
    print(f"  Geometric Mean: {analysis.scale.geometric_mean:.6e}")
    if np.isfinite(analysis.scale.harmonic_mean):
        print(f"  Harmonic Mean: {analysis.scale.harmonic_mean:.6e}")
    print()
    
    print("PROBABILITY ANALYSIS")
    print("-" * 80)
    print(f"  Distribution Type: {analysis.probability.distribution_type}")
    print(f"  Mean: {analysis.probability.mean:.6e}")
    print(f"  Std Dev: {analysis.probability.std_dev:.6e}")
    print(f"  Skewness: {analysis.probability.skewness:.3f}")
    print(f"  Kurtosis: {analysis.probability.kurtosis:.3f}")
    print(f"  Fit Quality (R²): {analysis.probability.fit_quality:.4f}")
    print(f"  95% CI: [{analysis.probability.confidence_interval_95[0]:.6e}, {analysis.probability.confidence_interval_95[1]:.6e}]")
    if analysis.probability.outliers:
        print(f"  Outliers (>3σ): {', '.join(analysis.probability.outliers)}")
    print()
    
    if analysis.fine_fit:
        print("FINE FIT RATIO")
        print("-" * 80)
        print(f"  Ratio Mean: {analysis.fine_fit.ratio_mean:.4f}")
        print(f"  Ratio Std Dev: {analysis.fine_fit.ratio_std:.4f}")
        print(f"  Percent Error Mean: {analysis.fine_fit.percent_error_mean:.2f}%")
        print(f"  Percent Error Std Dev: {analysis.fine_fit.percent_error_std:.2f}%")
        print(f"  Max Deviation: {analysis.fine_fit.max_deviation:.6e} ({analysis.fine_fit.max_deviation_equation})")
        print(f"  Chi-Squared: {analysis.fine_fit.chi_squared:.4f}")
        print(f"  Reduced Chi-Squared: {analysis.fine_fit.reduced_chi_squared:.4f}")
        print()
    
    print("=" * 80)


# ═══════════════════════════════════════════════════════════════════════════════
# COMMAND-LINE INTERFACE
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    import sys
    
    print("QCalc_stat.py - Statistical Analysis Module")
    print("=" * 80)
    print()
    
    # Example usage with test data
    print("Running example analysis...")
    print()
    
    try:
        # Create test data (simulate QCalc results)
        from QCalc import UnifiedFieldSolver, ComputeParams, CONSTANTS
        
        solver = UnifiedFieldSolver()
        params = ComputeParams(
            query_name="Statistical Analysis Test",
            M=1.4 * CONSTANTS['M_sun'],
            r=20e3,
            B=1e10,
            t=1e8
        )
        
        result = solver.solve(params)
        
        # Analyze results (filter by m/s² unit for consistency)
        analysis = compute_triple_point_analysis(result, unit_filter='m/s²')
        
        # Print summary
        print_analysis_summary(analysis)
        
        # Save to JSON
        output_file = "triple_point_analysis.json"
        analysis.save(output_file)
        print(f"Analysis saved to {output_file}")
        
    except Exception as e:
        print(f"Error: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
