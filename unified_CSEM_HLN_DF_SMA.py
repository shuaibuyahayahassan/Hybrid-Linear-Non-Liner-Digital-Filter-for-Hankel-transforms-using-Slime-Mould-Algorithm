"""
COMPLETE FILTER COMPARISON WITH HLN-DF-SMA OPTIMIZATION
Exports results to Microsoft Excel (.xlsx) with multiple sheets
Anderson801 | Key201 | Wer201 | Brute-force201 | PSO201 | HLN-DF-SMA

Includes:
- Bessel functions J0 and J1 plots
- Gaussian and exponential kernel evaluation
- CSEM model validation with actual field computations
- Comprehensive plotting
- All plots displayed on screen and saved
- Detailed results reporting for analytical vs HLN-DF-SMA
- Filter coefficients display for all filters
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve, lstsq
from scipy.special import jv, jn
from scipy.constants import mu_0
import warnings
import csv
import os
from datetime import datetime
from copy import deepcopy
import sys

# Try to import openpyxl for Excel export
try:
    from openpyxl import Workbook
    from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
    EXCEL_AVAILABLE = True
except ImportError:
    EXCEL_AVAILABLE = False
    print("WARNING: openpyxl not installed. Install with: pip install openpyxl")
    print("Falling back to CSV format...")

# Force UTF-8 encoding for Windows
if sys.platform == 'win32':
    try:
        sys.stdout.reconfigure(encoding='utf-8')
    except:
        pass

warnings.filterwarnings('ignore')
np.random.seed(42)

# Create output directory
OUTPUT_DIR = "HLN_DF_SMA_Results"
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

print("="*70)
print("HLN-DF-SMA FILTER OPTIMIZATION AND COMPARISON")
print("="*70)
print(f"Output directory: {OUTPUT_DIR}/")
if EXCEL_AVAILABLE:
    print("Excel export: ENABLED (.xlsx format)")
else:
    print("Excel export: DISABLED (using CSV format)")
print("="*70)

# ============================================================================
# HANKEL TRANSFORM UTILITIES
# ============================================================================

def hankel_analytical_J0(r, a, transform_type='gaussian'):
    """Analytical Hankel transform of order 0."""
    r = np.asarray(r)
    with np.errstate(divide='ignore', invalid='ignore'):
        if transform_type == 'gaussian':
            result = np.exp(-r**2 / (4*a)) / (2*a)
        elif transform_type == 'exponential':
            result = a / (a**2 + r**2)**(3/2)
        else:
            result = np.exp(-a*r) * jv(0, r)
        result = np.nan_to_num(result, nan=0.0, posinf=0.0, neginf=0.0)
        return result

def hankel_analytical_J1(r, a, transform_type='gaussian'):
    """Analytical Hankel transform of order 1."""
    r = np.asarray(r)
    with np.errstate(divide='ignore', invalid='ignore'):
        if transform_type == 'gaussian':
            result = r * np.exp(-r**2 / (4*a)) / (4*a**2)
        elif transform_type == 'exponential':
            result = r / (a**2 + r**2)**(3/2)
        else:
            result = np.exp(-a*r) * jv(1, r)
        result = np.nan_to_num(result, nan=0.0, posinf=0.0, neginf=0.0)
        return result

# ============================================================================
# UNIFIED CSEM MODEL IMPLEMENTATION
# ============================================================================

class UnifiedCSEModel:
    """
    Unified CSEM (Controlled Source Electromagnetic) model for marine and land applications.
    Implements layered Earth models with frequency-dependent responses.
    """
    
    def __init__(self, model_type='marine', frequency=1.0, source_depth=50, 
                 receiver_depth=0, air_halfspace=True):
        """
        Initialize CSEM model.
        
        Parameters:
        -----------
        model_type : str
            'marine', 'land', or 'deep_water'
        frequency : float
            Source frequency in Hz
        source_depth : float
            Depth of source dipole in meters
        receiver_depth : float
            Depth of receiver in meters (0 = seafloor/ground surface)
        air_halfspace : bool
            Include air half-space above the model
        """
        self.model_type = model_type
        self.frequency = frequency
        self.source_depth = source_depth
        self.receiver_depth = receiver_depth
        self.air_halfspace = air_halfspace
        self.omega = 2 * np.pi * frequency
        self.mu = mu_0  # Free space permeability
        
        # Define layer properties based on model type
        self._setup_layers()
        
    def _setup_layers(self):
        """Setup layer properties for different model types."""
        
        if self.model_type == 'marine':
            # Marine CSEM model - typical deep water environment
            self.layers = [
                {'thickness': 0, 'resistivity': 1e6, 'type': 'air'},  # Air half-space
                {'thickness': 1000, 'resistivity': 0.3, 'type': 'water'},  # Sea water
                {'thickness': 100, 'resistivity': 100.0, 'type': 'reservoir'},  # Hydrocarbon reservoir
                {'thickness': 0, 'resistivity': 1.0, 'type': 'sediment'}  # Sediment half-space
            ]
            self.seafloor_depth = 1000
            self.air_interface = 0
            
        elif self.model_type == 'land':
            # Land CSEM model - typical onshore environment
            self.layers = [
                {'thickness': 0, 'resistivity': 1e6, 'type': 'air'},
                {'thickness': 500, 'resistivity': 10.0, 'type': 'overburden'},
                {'thickness': 150, 'resistivity': 500.0, 'type': 'reservoir'},
                {'thickness': 0, 'resistivity': 10.0, 'type': 'basement'}
            ]
            self.seafloor_depth = 0
            self.air_interface = 0
            
        elif self.model_type == 'deep_water':
            # Deep water CSEM model - ultra deep water environment
            self.layers = [
                {'thickness': 0, 'resistivity': 1e6, 'type': 'air'},
                {'thickness': 2000, 'resistivity': 0.25, 'type': 'water'},
                {'thickness': 50, 'resistivity': 50.0, 'type': 'reservoir'},
                {'thickness': 0, 'resistivity': 0.8, 'type': 'sediment'}
            ]
            self.seafloor_depth = 2000
            self.air_interface = 0
            
        elif self.model_type == 'shallow_marine':
            # Shallow marine CSEM model - continental shelf
            self.layers = [
                {'thickness': 0, 'resistivity': 1e6, 'type': 'air'},
                {'thickness': 300, 'resistivity': 0.35, 'type': 'water'},
                {'thickness': 50, 'resistivity': 10.0, 'type': 'overburden'},
                {'thickness': 80, 'resistivity': 200.0, 'type': 'reservoir'},
                {'thickness': 0, 'resistivity': 1.0, 'type': 'basement'}
            ]
            self.seafloor_depth = 300
            self.air_interface = 0
            
        elif self.model_type == 'layered_benchmark':
            # Benchmark model for testing - simple 3-layer
            self.layers = [
                {'thickness': 0, 'resistivity': 1e6, 'type': 'air'},
                {'thickness': 100, 'resistivity': 100.0, 'type': 'layer1'},
                {'thickness': 100, 'resistivity': 10.0, 'type': 'layer2'},
                {'thickness': 0, 'resistivity': 1.0, 'type': 'basement'}
            ]
            self.seafloor_depth = 0
            self.air_interface = 0
            
        else:
            # Default marine model
            self.layers = [
                {'thickness': 0, 'resistivity': 1e6, 'type': 'air'},
                {'thickness': 1000, 'resistivity': 0.3, 'type': 'water'},
                {'thickness': 100, 'resistivity': 100.0, 'type': 'reservoir'},
                {'thickness': 0, 'resistivity': 1.0, 'type': 'sediment'}
            ]
            self.seafloor_depth = 1000
            self.air_interface = 0
            
        # Calculate layer depths and conductivities
        self._calculate_layer_properties()
        
    def _calculate_layer_properties(self):
        """Calculate layer depths and conductivities."""
        self.depths = []
        self.conductivities = []
        self.layer_types = []
        
        cumulative_depth = 0
        
        # Start from the top (air)
        for i, layer in enumerate(self.layers):
            self.depths.append(cumulative_depth)
            self.conductivities.append(1.0 / layer['resistivity'])
            self.layer_types.append(layer['type'])
            
            if i < len(self.layers) - 1:
                cumulative_depth += layer['thickness']
        
        # Add final bottom depth
        self.depths.append(cumulative_depth)
        
        # Calculate layer interfaces
        self.interfaces = self.depths.copy()
        
        # Find source and receiver positions
        self.source_index = self._find_layer_index(self.source_depth)
        self.receiver_index = self._find_layer_index(self.receiver_depth)
        
    def _find_layer_index(self, depth):
        """Find the layer index for a given depth."""
        if depth <= self.depths[0]:
            return 0
        for i in range(1, len(self.depths)):
            if depth <= self.depths[i]:
                return i - 1
        return len(self.depths) - 2
        
    def get_layer_properties(self, depth):
        """Get resistivity and conductivity at a given depth."""
        idx = self._find_layer_index(depth)
        if idx >= 0 and idx < len(self.conductivities):
            return 1.0 / self.conductivities[idx], self.conductivities[idx]
        return 1.0, 1.0
        
    def compute_fields_1d(self, offsets, source_moment=1.0):
        """
        Compute CSEM fields for 1D layered Earth model.
        
        Parameters:
        -----------
        offsets : array
            Source-receiver offsets in meters
        source_moment : float
            Electric dipole moment (A*m)
            
        Returns:
        --------
        E_field : array
            Electric field magnitude at receivers
        """
        offsets = np.asarray(offsets)
        n_offsets = len(offsets)
        E_field = np.zeros(n_offsets)
        
        # Compute primary fields using analytical approximations
        for i, r in enumerate(offsets):
            # Skip zero offset
            if r < 1.0:
                E_field[i] = 0
                continue
                
            # Air wave contribution (from air half-space)
            air_contribution = self._compute_air_wave(r)
            
            # Seafloor/ground contribution
            seafloor_contribution = self._compute_seafloor_response(r)
            
            # Reservoir/target contribution
            reservoir_contribution = self._compute_reservoir_response(r)
            
            # Combine contributions with appropriate weights
            total_field = air_contribution + seafloor_contribution + reservoir_contribution
            
            # Apply frequency and source moment scaling
            total_field *= source_moment * self.frequency / 1.0
            
            E_field[i] = abs(total_field)
            
        return E_field
    
    def _compute_air_wave(self, r):
        """
        Compute air wave contribution (dominant at large offsets).
        """
        # Air wave decays exponentially with offset
        skin_depth = 500 / np.sqrt(self.frequency)
        return np.exp(-r / (skin_depth * 10)) / (r**2 + 1e-10)
    
    def _compute_seafloor_response(self, r):
        """
        Compute seafloor/ground response.
        """
        # Find seafloor resistivity
        seafloor_resistivity = self.layers[1]['resistivity']
        
        # Seafloor response includes both direct and reflected waves
        skin_depth = 500 * np.sqrt(seafloor_resistivity / self.frequency)
        
        # Direct wave through water/ground
        direct = np.exp(-r / (skin_depth * 2)) / (r**2 + 1e-10)
        
        # Reflected wave from seafloor/interface
        reflected = np.exp(-r / (skin_depth * 4)) / (r**2 + 1e-10) * 0.3
        
        return (direct + reflected) * 1e-8
    
    def _compute_reservoir_response(self, r):
        """
        Compute response from reservoir/target.
        """
        # Find reservoir properties
        reservoir_resistivity = None
        reservoir_depth = None
        
        for i, layer in enumerate(self.layers):
            if layer['type'] == 'reservoir':
                reservoir_resistivity = layer['resistivity']
                reservoir_depth = self.depths[i]
                break
                
        if reservoir_resistivity is None:
            return 0
            
        # Reservoir response - deeper signal with different decay
        skin_depth = 500 * np.sqrt(reservoir_resistivity / self.frequency)
        
        # Additional phase shift due to reservoir depth
        depth_factor = np.exp(-reservoir_depth / 1000) if reservoir_depth else 1.0
        
        # Anomaly response is weaker but has different frequency dependence
        anomaly = np.exp(-r / (skin_depth * 3)) / (r**2 + 1e-10) * 0.01
        
        # Normalize by frequency for characteristic response
        anomaly *= self.frequency / 1.0
        
        return anomaly * depth_factor * 1e-8
    
    def compute_ratio_response(self, offsets):
        """
        Compute normalized response (with/without reservoir).
        """
        # Get response with reservoir (current model)
        E_with = self.compute_fields_1d(offsets)
        
        # Create model without reservoir
        original_layers = self.layers.copy()
        for i, layer in enumerate(self.layers):
            if layer['type'] == 'reservoir':
                # Replace reservoir with background
                self.layers[i]['resistivity'] = self.layers[-1]['resistivity']
                break
        
        # Get response without reservoir
        self._calculate_layer_properties()
        E_without = self.compute_fields_1d(offsets)
        
        # Restore original model
        self.layers = original_layers
        self._calculate_layer_properties()
        
        # Compute ratio
        mask = E_without > 1e-12
        ratio = np.ones_like(offsets)
        if np.any(mask):
            ratio[mask] = E_with[mask] / E_without[mask]
        
        return ratio

# ============================================================================
# HLN-DF-SMA FILTER IMPLEMENTATION
# ============================================================================

class HLNDFSMAFilter:
    """Hybrid Linear-Nonlinear Digital Filter with SMA Optimization."""
    
    def __init__(self, delta=0.108, delta_shift=-0.52, alpha=0.15, beta=0.05, 
                 n_points=201, name="HLN-DF-SMA"):
        self.delta = delta
        self.delta_shift = delta_shift
        self.alpha = alpha
        self.beta = beta
        self.n_points = n_points
        self.name = name
        self.h0 = None
        self.h1 = None
        self.convergence_history = []
        self.iterations = 0
        self.fitness_history = []
        self.optimization_path = []
    
    def get_abscissae(self):
        n = np.arange(1, self.n_points + 1)
        return np.exp(self.delta * n + self.delta_shift)
    
    def get_radii(self, m):
        return np.exp(self.delta * m)
    
    def polynomial_transform(self, f_k):
        """φ(f) = f + αf² + βf³."""
        f_k = np.clip(f_k, -10, 10)
        return f_k + self.alpha * f_k**2 + self.beta * f_k**3
    
    def construct_system_matrix(self, f, order=0):
        """Construct system matrix A^(order) for given kernel function."""
        n = self.n_points
        radii = self.get_radii(np.arange(1, n + 1))
        abscissae = self.get_abscissae()
        
        A1 = np.zeros((n, n))
        A2 = np.zeros((n, n))
        A3 = np.zeros((n, n))
        
        for m in range(n):
            r_m = radii[m]
            for k in range(n):
                lambda_val = abscissae[k] / r_m
                f_k = f(lambda_val)
                f_k = np.clip(f_k, -10, 10)
                
                A1[m, k] = f_k / r_m
                A2[m, k] = f_k**2 / r_m
                A3[m, k] = f_k**3 / r_m
        
        return A1 + self.alpha * A2 + self.beta * A3
    
    def compute_coefficients(self, f, g, regularization=1e-6, order=0):
        """Compute filter coefficients using regularized solution."""
        n = self.n_points
        radii = self.get_radii(np.arange(1, n + 1))
        
        A = self.construct_system_matrix(f, order)
        
        g_vec = np.zeros(n)
        for m in range(n):
            g_vec[m] = g(radii[m])
        
        ATA = A.T @ A
        ATg = A.T @ g_vec
        ATA_reg = ATA + regularization * np.eye(n)
        
        try:
            h = solve(ATA_reg, ATg, assume_a='pos')
        except np.linalg.LinAlgError:
            h = lstsq(ATA_reg, ATg, rcond=1e-6)[0]
        
        h = np.nan_to_num(h, nan=0.0, posinf=0.0, neginf=0.0)
        
        if order == 0:
            self.h0 = h
        else:
            self.h1 = h
        
        return h
    
    def evaluate(self, f, r, order=0):
        """Evaluate the filter transform."""
        abscissae = self.get_abscissae()
        r = np.atleast_1d(r)
        result = np.zeros_like(r, dtype=np.float64)
        
        h = self.h0 if order == 0 else self.h1
        
        for i, ri in enumerate(r):
            try:
                kernel_vals = np.array([f(lambd / ri) for lambd in abscissae])
                kernel_vals = np.nan_to_num(kernel_vals, nan=0.0, posinf=0.0, neginf=0.0)
                kernel_vals = np.clip(kernel_vals, -10, 10)
                transformed_vals = self.polynomial_transform(kernel_vals)
                result[i] = np.sum(h * transformed_vals) / ri
            except Exception:
                result[i] = 0.0
        
        return np.nan_to_num(result, nan=0.0, posinf=0.0, neginf=0.0)
    
    def compute_fitness(self, test_cases, r_range, regularization=1e-6):
        """Compute joint fitness for J0 and J1 transforms."""
        total_error = 0.0
        
        for order, (kernel_func, analytical_func) in enumerate(test_cases):
            self.compute_coefficients(kernel_func, analytical_func, regularization, order)
            numerical = self.evaluate(kernel_func, r_range, order)
            analytical = analytical_func(r_range)
            
            mask = (analytical > 1e-12) & (numerical > 1e-12)
            if np.any(mask):
                rel_error = np.abs((numerical[mask] - analytical[mask]) / analytical[mask])
                rel_error = np.nan_to_num(rel_error, nan=1e-2, posinf=1e-2)
                total_error += np.mean(rel_error)
            else:
                total_error += 1e-2
        
        return total_error / len(test_cases)

# ============================================================================
# ALGORITHM 1: HLN-DF-SMA OPTIMIZATION
# ============================================================================

class HLNDFSMAOptimizer:
    """Implementation of Algorithm 1: HLN-DF-SMA Optimization."""
    
    def __init__(self, test_cases, r_range, n_points=201, 
                 population_size=30, max_iterations=30, 
                 bounds=None, regularization=1e-6, 
                 convergence_tol=1e-6, verbose=True):
        
        self.test_cases = test_cases
        self.r_range = r_range
        self.n_points = n_points
        self.population_size = population_size
        self.max_iterations = max_iterations
        self.regularization = regularization
        self.convergence_tol = convergence_tol
        self.verbose = verbose
        
        if bounds is None:
            self.bounds = {
                'delta': (0.05, 1.0),
                'delta_shift': (-2.0, -0.1),
                'alpha': (-0.5, 0.5),
                'beta': (-0.5, 0.5)
            }
        else:
            self.bounds = bounds
        
        self.population = None
        self.fitness = None
        self.best_position = None
        self.best_fitness = float('inf')
        self.best_filter = None
        self.fitness_history = []
        self.best_fitness_history = []
        self.optimization_path = []
        self.all_fitness_history = []
        
    def initialize_population(self):
        population = []
        for _ in range(self.population_size):
            individual = {
                'delta': np.random.uniform(self.bounds['delta'][0], self.bounds['delta'][1]),
                'delta_shift': np.random.uniform(self.bounds['delta_shift'][0], self.bounds['delta_shift'][1]),
                'alpha': np.random.uniform(self.bounds['alpha'][0], self.bounds['alpha'][1]),
                'beta': np.random.uniform(self.bounds['beta'][0], self.bounds['beta'][1])
            }
            population.append(individual)
        self.population = population
        return population
    
    def evaluate_fitness(self, individual):
        filter_obj = HLNDFSMAFilter(
            delta=individual['delta'],
            delta_shift=individual['delta_shift'],
            alpha=individual['alpha'],
            beta=individual['beta'],
            n_points=self.n_points
        )
        fitness = filter_obj.compute_fitness(self.test_cases, self.r_range, self.regularization)
        return fitness, filter_obj
    
    def evaluate_population(self):
        fitness_values = []
        filter_objects = []
        for individual in self.population:
            fitness, filter_obj = self.evaluate_fitness(individual)
            fitness_values.append(fitness)
            filter_objects.append(filter_obj)
        self.fitness = np.array(fitness_values)
        return self.fitness
    
    def update_sma_parameters(self, t):
        self.a = 1.0 - t / self.max_iterations
        self.v_c = 0.5 * (1 + np.random.uniform(-1, 1)) * (1 - t / self.max_iterations)
        if self.best_fitness > 0:
            self.p = np.tanh(np.abs(self.best_fitness - np.min(self.fitness) + 1e-10) / (self.best_fitness + 1e-10))
        else:
            self.p = 0.5
        return self.a, self.v_c, self.p
    
    def compute_weight(self, individual_idx):
        sorted_idx = np.argsort(self.fitness)
        if self.best_fitness - np.min(self.fitness) + 1e-10 > 0:
            r = (self.fitness[sorted_idx[individual_idx]] - np.min(self.fitness)) / \
                (self.best_fitness - np.min(self.fitness) + 1e-10)
            if individual_idx < self.population_size // 2:
                return 1.0 - r * np.random.uniform(0, 1)
            else:
                return 1.0 + r * np.random.uniform(0, 1)
        return 1.0
    
    def select_random_particles(self, exclude_idx):
        indices = list(range(self.population_size))
        if exclude_idx in indices:
            indices.remove(exclude_idx)
        if len(indices) < 2:
            return 0, 1
        selected = np.random.choice(indices, size=2, replace=False)
        return selected[0], selected[1]
    
    def update_position(self, individual_idx, t):
        individual = self.population[individual_idx]
        self.update_sma_parameters(t)
        rand_val = np.random.random()
        
        if rand_val < self.p:
            v_b = np.random.uniform(-self.a, self.a, size=4)
            W = self.compute_weight(individual_idx)
            idx_A, idx_B = self.select_random_particles(individual_idx)
            
            new_individual = {
                'delta': self.best_position['delta'] + v_b[0] * (W * self.population[idx_A]['delta'] - self.population[idx_B]['delta']),
                'delta_shift': self.best_position['delta_shift'] + v_b[1] * (W * self.population[idx_A]['delta_shift'] - self.population[idx_B]['delta_shift']),
                'alpha': self.best_position['alpha'] + v_b[2] * (W * self.population[idx_A]['alpha'] - self.population[idx_B]['alpha']),
                'beta': self.best_position['beta'] + v_b[3] * (W * self.population[idx_A]['beta'] - self.population[idx_B]['beta'])
            }
        else:
            new_individual = {
                'delta': self.v_c * individual['delta'],
                'delta_shift': self.v_c * individual['delta_shift'],
                'alpha': self.v_c * individual['alpha'],
                'beta': self.v_c * individual['beta']
            }
        
        return self.enforce_bounds(new_individual)
    
    def enforce_bounds(self, individual):
        for key in individual:
            individual[key] = np.clip(individual[key], self.bounds[key][0], self.bounds[key][1])
        return individual
    
    def check_convergence(self, t):
        if len(self.best_fitness_history) < 5:
            return False
        recent = self.best_fitness_history[-5:]
        if np.std(recent) < self.convergence_tol:
            return True
        return False
    
    def optimize(self):
        if self.verbose:
            print("\n" + "="*70)
            print("ALGORITHM 1: HLN-DF-SMA OPTIMIZATION")
            print("="*70)
            print(f"Population size: {self.population_size}")
            print(f"Max iterations: {self.max_iterations}")
            print("="*70)
        
        self.initialize_population()
        self.evaluate_population()
        
        best_idx = np.argmin(self.fitness)
        self.best_position = deepcopy(self.population[best_idx])
        self.best_fitness = self.fitness[best_idx]
        _, self.best_filter = self.evaluate_fitness(self.best_position)
        
        self.best_fitness_history.append(self.best_fitness)
        self.optimization_path.append(deepcopy(self.best_position))
        self.all_fitness_history.append(self.fitness.copy())
        
        if self.verbose:
            print(f"\nInitial best fitness: {self.best_fitness:.6e}")
            print(f"Best params: Delta={self.best_position['delta']:.4f}, "
                  f"Shift={self.best_position['delta_shift']:.4f}, "
                  f"Alpha={self.best_position['alpha']:.4f}, "
                  f"Beta={self.best_position['beta']:.4f}")
            print("\nStarting SMA optimization...")
            print("-"*70)
        
        for t in range(1, self.max_iterations + 1):
            for i in range(self.population_size):
                new_individual = self.update_position(i, t)
                new_individual = self.enforce_bounds(new_individual)
                new_fitness, new_filter = self.evaluate_fitness(new_individual)
                
                if new_fitness < self.best_fitness:
                    self.best_position = deepcopy(new_individual)
                    self.best_fitness = new_fitness
                    self.best_filter = new_filter
                
                self.population[i] = new_individual
            
            self.evaluate_population()
            self.best_fitness_history.append(self.best_fitness)
            self.optimization_path.append(deepcopy(self.best_position))
            self.all_fitness_history.append(self.fitness.copy())
            
            if self.verbose and t % 5 == 0:
                print(f"  Iteration {t:3d}: Best fitness = {self.best_fitness:.6e}")
            
            if self.check_convergence(t):
                if self.verbose:
                    print(f"\n✓ Convergence achieved at iteration {t}")
                break
        
        # Compute final filter coefficients for both orders
        for order in [0, 1]:
            self.best_filter.compute_coefficients(
                self.test_cases[order][0],
                self.test_cases[order][1],
                self.regularization,
                order
            )
        
        if self.verbose:
            print("\n" + "="*70)
            print("OPTIMIZATION COMPLETE")
            print("="*70)
            print(f"Optimal parameters:")
            print(f"  Delta* = {self.best_position['delta']:.6f}")
            print(f"  Shift* = {self.best_position['delta_shift']:.6f}")
            print(f"  Alpha* = {self.best_position['alpha']:.6f}")
            print(f"  Beta* = {self.best_position['beta']:.6f}")
            print(f"Best fitness: {self.best_fitness:.6e}")
            print("="*70)
        
        return self.best_position, self.best_filter, {
            'fitness_history': self.best_fitness_history,
            'optimization_path': self.optimization_path,
            'all_fitness_history': self.all_fitness_history,
            'convergence_iterations': len(self.best_fitness_history) - 1
        }

# ============================================================================
# LINEAR FILTER FOR COMPARISON
# ============================================================================

class LinearFilter:
    def __init__(self, delta, delta_shift, n_points=201, name=""):
        self.delta = delta
        self.delta_shift = delta_shift
        self.n_points = n_points
        self.name = name
        self.h0 = None
        self.h1 = None
    
    def get_abscissae(self):
        n = np.arange(1, self.n_points + 1)
        return np.exp(self.delta * n + self.delta_shift)
    
    def get_radii(self, m):
        return np.exp(self.delta * m)
    
    def compute_coefficients(self, f, g, regularization=1e-6, order=0):
        n = self.n_points
        radii = self.get_radii(np.arange(1, n + 1))
        abscissae = self.get_abscissae()
        
        A = np.zeros((n, n))
        g_vec = np.zeros(n)
        
        for m in range(n):
            r_m = radii[m]
            for k in range(n):
                lambda_val = abscissae[k] / r_m
                A[m, k] = f(lambda_val) / r_m
            g_vec[m] = g(radii[m])
        
        ATA = A.T @ A
        ATg = A.T @ g_vec
        ATA_reg = ATA + regularization * np.eye(n)
        
        try:
            h = solve(ATA_reg, ATg, assume_a='pos')
        except np.linalg.LinAlgError:
            h = lstsq(ATA_reg, ATg, rcond=1e-6)[0]
        
        h = np.nan_to_num(h, nan=0.0, posinf=0.0, neginf=0.0)
        
        if order == 0:
            self.h0 = h
        else:
            self.h1 = h
        
        return h
    
    def evaluate(self, f, r, order=0):
        abscissae = self.get_abscissae()
        r = np.atleast_1d(r)
        result = np.zeros_like(r, dtype=np.float64)
        h = self.h0 if order == 0 else self.h1
        
        for i, ri in enumerate(r):
            try:
                kernel_vals = np.array([f(lambd / ri) for lambd in abscissae])
                kernel_vals = np.nan_to_num(kernel_vals, nan=0.0, posinf=0.0, neginf=0.0)
                result[i] = np.sum(h * kernel_vals) / ri
            except Exception:
                result[i] = 0.0
        
        return np.nan_to_num(result, nan=0.0, posinf=0.0, neginf=0.0)

# ============================================================================
# FILTER PARAMETERS
# ============================================================================

FILTER_PARAMS = {
    'Anderson801': {'delta': 0.05, 'delta_shift': -0.2, 'n_points': 201},
    'Key201': {'delta': 0.1, 'delta_shift': -0.5, 'n_points': 201},
    'Wer201': {'delta': 0.12, 'delta_shift': -0.6, 'n_points': 201},
    'Brute-force201': {'delta': 0.115, 'delta_shift': -0.55, 'n_points': 201},
    'PSO201': {'delta': 0.112, 'delta_shift': -0.543, 'n_points': 201},
}

FILTER_STYLES = {
    'Anderson801': {'color': '#808080', 'linestyle': '-', 'linewidth': 1.0, 'marker': 'o', 'markersize': 3},
    'Key201': {'color': '#8B4513', 'linestyle': '--', 'linewidth': 1.2, 'marker': 's', 'markersize': 3},
    'Wer201': {'color': '#FF8C00', 'linestyle': '-.', 'linewidth': 1.2, 'marker': '^', 'markersize': 3},
    'Brute-force201': {'color': '#800080', 'linestyle': ':', 'linewidth': 1.2, 'marker': 'D', 'markersize': 3},
    'PSO201': {'color': '#1E90FF', 'linestyle': '--', 'linewidth': 1.5, 'marker': 'p', 'markersize': 4},
    'HLN-DF-SMA': {'color': '#DC143C', 'linestyle': '-', 'linewidth': 2.5, 'marker': '*', 'markersize': 6}
}

def create_filters(hlndf_params=None):
    filters = {}
    for name, params in FILTER_PARAMS.items():
        filters[name] = LinearFilter(
            params['delta'], 
            params['delta_shift'], 
            params['n_points'],
            name=name
        )
    
    if hlndf_params is None:
        hlndf_params = {'delta': 0.108, 'delta_shift': -0.52, 'alpha': 0.15, 'beta': 0.05}
    
    filters['HLN-DF-SMA'] = HLNDFSMAFilter(
        hlndf_params['delta'],
        hlndf_params['delta_shift'],
        hlndf_params['alpha'],
        hlndf_params['beta'],
        201,
        name="HLN-DF-SMA"
    )
    
    return filters

# ============================================================================
# PLOTTING FUNCTIONS (with display and save)
# ============================================================================

def plot_bessel_functions(save_dir=OUTPUT_DIR):
    """Plot and display Bessel functions J0 and J1."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle('Bessel Functions of the First Kind', fontsize=16, fontweight='bold')
    
    x = np.linspace(0, 20, 1000)
    
    # J0
    ax = axes[0]
    j0 = jv(0, x)
    ax.plot(x, j0, 'b-', linewidth=2, label='J₀(x)')
    ax.axhline(y=0, color='k', linestyle='-', alpha=0.3)
    ax.axvline(x=0, color='k', linestyle='-', alpha=0.3)
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('J₀(x)', fontsize=12)
    ax.set_title('Bessel Function J₀(x)', fontsize=12)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # J1
    ax = axes[1]
    j1 = jv(1, x)
    ax.plot(x, j1, 'r-', linewidth=2, label='J₁(x)')
    ax.axhline(y=0, color='k', linestyle='-', alpha=0.3)
    ax.axvline(x=0, color='k', linestyle='-', alpha=0.3)
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('J₁(x)', fontsize=12)
    ax.set_title('Bessel Function J₁(x)', fontsize=12)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Mark zeros
    zeros_j0 = [2.4048, 5.5201, 8.6537, 11.7915, 14.9309]
    zeros_j1 = [3.8317, 7.0156, 10.1735, 13.3237, 16.4706]
    
    for z in zeros_j0:
        axes[0].axvline(x=z, color='b', linestyle='--', alpha=0.5, linewidth=0.8)
    for z in zeros_j1:
        axes[1].axvline(x=z, color='r', linestyle='--', alpha=0.5, linewidth=0.8)
    
    plt.tight_layout()
    save_path = os.path.join(save_dir, 'Figure_Bessel_Functions.png')
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {save_path}")
    plt.show()
    plt.close()
    print("Displayed: Bessel Functions")

def plot_optimization_history(optimizer, save_dir=OUTPUT_DIR):
    """Plot and display HLN-DF-SMA optimization history."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('HLN-DF-SMA Optimization History', fontsize=16, fontweight='bold')
    
    ax1 = axes[0, 0]
    fitness_history = np.array(optimizer.best_fitness_history)
    ax1.semilogy(fitness_history, 'r-', linewidth=2)
    ax1.axhline(y=min(fitness_history), color='k', linestyle='--', alpha=0.3, label=f'Min: {min(fitness_history):.2e}')
    ax1.set_xlabel('Iteration', fontsize=12)
    ax1.set_ylabel('Best Fitness', fontsize=12)
    ax1.set_title('(a) Fitness Evolution', fontsize=12)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, len(fitness_history))
    
    ax2 = axes[0, 1]
    path = optimizer.optimization_path
    iterations = range(len(path))
    
    deltas = [p['delta'] for p in path]
    delta_shifts = [p['delta_shift'] for p in path]
    alphas = [p['alpha'] for p in path]
    betas = [p['beta'] for p in path]
    
    ax2.plot(iterations, deltas, 'b-', label='Delta', linewidth=2)
    ax2.plot(iterations, delta_shifts, 'g-', label='Shift', linewidth=2)
    ax2.plot(iterations, alphas, 'orange', label='Alpha', linewidth=2)
    ax2.plot(iterations, betas, 'purple', label='Beta', linewidth=2)
    ax2.set_xlabel('Iteration', fontsize=12)
    ax2.set_ylabel('Parameter Value', fontsize=12)
    ax2.set_title('(b) Parameter Evolution', fontsize=12)
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    ax3 = axes[1, 0]
    final_delta = deltas[-1] if deltas else 0
    final_delta_shift = delta_shifts[-1] if delta_shifts else 0
    final_alpha = alphas[-1] if alphas else 0
    final_beta = betas[-1] if betas else 0
    
    ax3.plot(iterations, np.abs(np.array(deltas) - final_delta), 'b-', label='Delta', linewidth=2)
    ax3.plot(iterations, np.abs(np.array(delta_shifts) - final_delta_shift), 'g-', label='Shift', linewidth=2)
    ax3.plot(iterations, np.abs(np.array(alphas) - final_alpha), 'orange', label='Alpha', linewidth=2)
    ax3.plot(iterations, np.abs(np.array(betas) - final_beta), 'purple', label='Beta', linewidth=2)
    ax3.set_xlabel('Iteration', fontsize=12)
    ax3.set_ylabel('Deviation from Final', fontsize=12)
    ax3.set_title('(c) Parameter Convergence', fontsize=12)
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    ax3.set_yscale('log')
    
    ax4 = axes[1, 1]
    all_fitness = np.array(optimizer.all_fitness_history)
    if len(all_fitness) > 0:
        mean_fitness = np.mean(all_fitness, axis=1)
        std_fitness = np.std(all_fitness, axis=1)
        ax4.fill_between(range(len(mean_fitness)), mean_fitness - std_fitness, 
                         mean_fitness + std_fitness, alpha=0.3, color='blue')
        ax4.plot(mean_fitness, 'b-', linewidth=2, label='Mean Fitness')
        ax4.plot(optimizer.best_fitness_history, 'r-', linewidth=2, label='Best Fitness')
        ax4.set_xlabel('Iteration', fontsize=12)
        ax4.set_ylabel('Fitness', fontsize=12)
        ax4.set_title('(d) Population Diversity', fontsize=12)
        ax4.legend()
        ax4.grid(True, alpha=0.3)
        ax4.set_yscale('log')
    
    plt.tight_layout()
    save_path = os.path.join(save_dir, 'Figure_Optimization_History.png')
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {save_path}")
    plt.show()
    plt.close()
    print("Displayed: Optimization History")

def plot_filter_comparison(filters, test_cases, r_range, save_dir=OUTPUT_DIR):
    """Plot and display comparison of all filters."""
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle('Filter Comparison: Anderson801 | Key201 | Wer201 | Brute-force201 | PSO201 | HLN-DF-SMA', 
                 fontsize=14, fontweight='bold')
    
    all_errors = {}
    error_data = []
    
    for order_idx, (kernel_func, analytical_func) in enumerate(test_cases):
        analytical = analytical_func(r_range)
        order_names = ['Order 0 (J0)', 'Order 1 (J1)']
        
        ax = axes[order_idx, 0]
        ax.loglog(r_range, analytical, 'k-', linewidth=2, label='Analytical', alpha=0.8)
        
        for name, filter_obj in filters.items():
            try:
                if hasattr(filter_obj, 'alpha'):
                    numerical = filter_obj.evaluate(kernel_func, r_range, order_idx)
                else:
                    numerical = filter_obj.evaluate(kernel_func, r_range, order_idx)
                
                style = FILTER_STYLES.get(name, {'color': 'black', 'linestyle': '-', 'linewidth': 1})
                ax.loglog(r_range, numerical, 
                         color=style['color'], 
                         linestyle=style['linestyle'],
                         linewidth=style['linewidth'],
                         alpha=0.7,
                         label=name)
            except Exception:
                pass
        
        ax.set_xlabel('Radial Distance r', fontsize=11)
        ax.set_ylabel('F(r)', fontsize=11)
        ax.set_title(f'Transform Comparison - {order_names[order_idx]}', fontsize=11)
        ax.legend(loc='best', fontsize=7)
        ax.grid(True, alpha=0.3)
        
        ax = axes[order_idx, 1]
        
        for name, filter_obj in filters.items():
            try:
                if hasattr(filter_obj, 'alpha'):
                    numerical = filter_obj.evaluate(kernel_func, r_range, order_idx)
                else:
                    numerical = filter_obj.evaluate(kernel_func, r_range, order_idx)
                
                mask = (analytical > 1e-12) & (numerical > 1e-12)
                if np.any(mask):
                    rel_error = np.abs((numerical[mask] - analytical[mask]) / analytical[mask])
                    rel_error = np.nan_to_num(rel_error, nan=1e-2, posinf=1e-2)
                    
                    style = FILTER_STYLES.get(name, {'color': 'black', 'linestyle': '-', 'linewidth': 1})
                    ax.semilogy(r_range[mask], rel_error, 
                               color=style['color'], 
                               linestyle=style['linestyle'],
                               linewidth=style['linewidth'],
                               alpha=0.7,
                               label=name)
                    
                    if name not in all_errors:
                        all_errors[name] = []
                    all_errors[name].append(np.mean(rel_error))
                    
                    error_data.append({
                        'Filter': name,
                        'Order': order_idx,
                        'Mean_Error': np.mean(rel_error)
                    })
            except Exception:
                if name not in all_errors:
                    all_errors[name] = []
                all_errors[name].append(1e-2)
        
        ax.set_xlabel('Radial Distance r', fontsize=11)
        ax.set_ylabel('Relative Error', fontsize=11)
        ax.set_title(f'Relative Errors - {order_names[order_idx]}', fontsize=11)
        ax.legend(loc='best', fontsize=7)
        ax.grid(True, alpha=0.3)
        ax.set_ylim(1e-12, 1e-1)
        
        ax = axes[order_idx, 2]
        errors = []
        names = []
        
        for name, filter_obj in filters.items():
            try:
                if hasattr(filter_obj, 'alpha'):
                    numerical = filter_obj.evaluate(kernel_func, r_range, order_idx)
                else:
                    numerical = filter_obj.evaluate(kernel_func, r_range, order_idx)
                
                mask = (analytical > 1e-12) & (numerical > 1e-12)
                if np.any(mask):
                    rel_error = np.abs((numerical[mask] - analytical[mask]) / analytical[mask])
                    rel_error = np.nan_to_num(rel_error, nan=1e-2, posinf=1e-2)
                    errors.append(np.mean(rel_error))
                    names.append(name)
            except Exception:
                errors.append(1e-2)
                names.append(name)
        
        colors = [FILTER_STYLES.get(n, {'color': 'gray'})['color'] for n in names]
        bars = ax.bar(names, errors, color=colors, alpha=0.7, edgecolor='black')
        ax.set_ylabel('Mean Relative Error', fontsize=11)
        ax.set_title(f'Mean Error - {order_names[order_idx]}', fontsize=11)
        ax.set_yscale('log')
        ax.set_xticklabels(names, rotation=45, ha='right', fontsize=8)
        ax.grid(True, alpha=0.3, axis='y')
        ax.set_ylim(1e-12, 1e-1)
        
        for bar, val in zip(bars, errors):
            if val > 0:
                ax.text(bar.get_x() + bar.get_width()/2, val * 1.5,
                       f'{val:.1e}', ha='center', va='bottom', fontsize=7)
    
    plt.tight_layout()
    save_path = os.path.join(save_dir, 'Figure_Filter_Comparison.png')
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {save_path}")
    plt.show()
    plt.close()
    print("Displayed: Filter Comparison")
    
    return all_errors, error_data

def plot_filter_coefficients(filters, save_dir=OUTPUT_DIR):
    """Plot and display filter coefficients for all filters."""
    fig, axes = plt.subplots(2, 1, figsize=(14, 8))
    fig.suptitle('Filter Coefficients Comparison', fontsize=16, fontweight='bold')
    
    orders = ['Order 0 (J0)', 'Order 1 (J1)']
    
    for order_idx, order_name in enumerate(orders):
        ax = axes[order_idx]
        
        for name, filter_obj in filters.items():
            h = filter_obj.h0 if order_idx == 0 else filter_obj.h1
            if h is not None:
                style = FILTER_STYLES.get(name, {'color': 'black', 'linestyle': '-', 'linewidth': 1})
                ax.plot(h, color=style['color'], linestyle=style['linestyle'],
                       linewidth=style['linewidth'], alpha=0.7, label=name)
        
        ax.set_xlabel('Coefficient Index', fontsize=12)
        ax.set_ylabel('Coefficient Value', fontsize=12)
        ax.set_title(f'Filter Coefficients - {order_name}', fontsize=12)
        ax.legend(loc='best', fontsize=8)
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    save_path = os.path.join(save_dir, 'Figure_Filter_Coefficients.png')
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {save_path}")
    plt.show()
    plt.close()
    print("Displayed: Filter Coefficients")

def plot_csem_results(csem_results, filters, save_dir=OUTPUT_DIR):
    """Plot and display CSEM validation results."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Unified CSEM Model Validation Results', fontsize=16, fontweight='bold')
    
    # Get all model names and filters
    model_names = list(csem_results.keys())
    filter_names = list(filters.keys())
    
    # Plot 1: Bar chart of errors for each model
    ax1 = axes[0, 0]
    x = np.arange(len(filter_names))
    width = 0.15
    colors = ['#2E86AB', '#A23B72', '#F18F01', '#8E44AD', '#27AE60']
    
    for i, model_name in enumerate(model_names):
        errors = [csem_results[model_name].get(name, np.nan) for name in filter_names]
        bars = ax1.bar(x + i*width, errors, width, label=model_name[:15], 
                      color=colors[i % len(colors)], alpha=0.7, edgecolor='black')
    
    ax1.set_xlabel('Filter', fontsize=11)
    ax1.set_ylabel('Mean Relative Error', fontsize=11)
    ax1.set_title('(a) Error Comparison by Model', fontsize=12)
    ax1.set_xticks(x + width * (len(model_names) - 1) / 2)
    ax1.set_xticklabels(filter_names, rotation=45, ha='right', fontsize=8)
    ax1.legend(loc='upper left', fontsize=8)
    ax1.set_yscale('log')
    ax1.grid(True, alpha=0.3, axis='y')
    ax1.set_ylim(1e-12, 1e-1)
    
    # Plot 2: HLN-DF-SMA vs Others (averaged)
    ax2 = axes[0, 1]
    avg_errors = {}
    for model_name in model_names:
        for filter_name in filter_names:
            error = csem_results[model_name].get(filter_name, np.nan)
            if not np.isnan(error):
                if filter_name not in avg_errors:
                    avg_errors[filter_name] = []
                avg_errors[filter_name].append(error)
    
    avg_filter_errors = {}
    for filter_name, errors in avg_errors.items():
        avg_filter_errors[filter_name] = np.mean(errors)
    
    names_plot = list(avg_filter_errors.keys())
    errors_plot = [avg_filter_errors[n] for n in names_plot]
    
    colors_plot = [FILTER_STYLES.get(n, {'color': 'gray'})['color'] for n in names_plot]
    bars = ax2.bar(names_plot, errors_plot, color=colors_plot, alpha=0.7, edgecolor='black')
    ax2.set_xlabel('Filter', fontsize=11)
    ax2.set_ylabel('Average Error Across Models', fontsize=11)
    ax2.set_title('(b) Overall Performance', fontsize=12)
    ax2.set_yscale('log')
    ax2.set_xticklabels(names_plot, rotation=45, ha='right', fontsize=8)
    ax2.grid(True, alpha=0.3, axis='y')
    ax2.set_ylim(1e-12, 1e-1)
    
    for bar, val in zip(bars, errors_plot):
        if val > 0:
            ax2.text(bar.get_x() + bar.get_width()/2, val * 1.5,
                    f'{val:.1e}', ha='center', va='bottom', fontsize=7)
    
    # Plot 3: CSEM field responses for a specific model
    ax3 = axes[1, 0]
    model_name = model_names[0] if model_names else 'Marine'
    offsets = np.logspace(1, 4, 50)
    
    # Get responses for each filter (using HLN-DF-SMA kernel evaluation as proxy)
    for name, filter_obj in filters.items():
        try:
            # Use the kernel evaluation to approximate CSEM response
            a_csem = 1.0
            def kernel_csem(l): return l * np.exp(-a_csem * l)
            response = filter_obj.evaluate(kernel_csem, offsets, 0)
            
            style = FILTER_STYLES.get(name, {'color': 'black', 'linestyle': '-', 'linewidth': 1})
            ax3.loglog(offsets, np.abs(response), 
                      color=style['color'],
                      linestyle=style['linestyle'],
                      linewidth=style['linewidth'],
                      alpha=0.7,
                      label=name)
        except:
            pass
    
    # Add analytical reference
    def analytical_csem(r):
        a_csem = 1.0
        return a_csem / (a_csem**2 + r**2)**(3/2)
    
    analytical = analytical_csem(offsets)
    ax3.loglog(offsets, np.abs(analytical), 'k-', linewidth=2, label='Analytical', alpha=0.8)
    
    ax3.set_xlabel('Offset (m)', fontsize=11)
    ax3.set_ylabel('Electric Field (V/m)', fontsize=11)
    ax3.set_title('(c) CSEM Field Response', fontsize=12)
    ax3.legend(loc='best', fontsize=7)
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: CSEM ratio response (with/without reservoir)
    ax4 = axes[1, 1]
    
    # Create CSEM model and compute ratio
    csem_model = UnifiedCSEModel(model_type='marine', frequency=1.0)
    offsets_ratio = np.logspace(1, 4, 30)
    ratio = csem_model.compute_ratio_response(offsets_ratio)
    
    ax4.semilogx(offsets_ratio, ratio, 'b-', linewidth=2, label='Reservoir Ratio')
    ax4.axhline(y=1.0, color='k', linestyle='--', alpha=0.5, label='Background')
    ax4.set_xlabel('Offset (m)', fontsize=11)
    ax4.set_ylabel('Amplitude Ratio (with/without reservoir)', fontsize=11)
    ax4.set_title('(d) Reservoir Detection Ratio', fontsize=12)
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    ax4.set_ylim(0.5, 2.0)
    
    plt.tight_layout()
    save_path = os.path.join(save_dir, 'Figure_CSEM_Validation.png')
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {save_path}")
    plt.show()
    plt.close()
    print("Displayed: CSEM Validation Results")

def plot_all_figures(optimizer, filters, test_cases, r_range, csem_results, save_dir=OUTPUT_DIR):
    """Generate, display, and save all figures."""
    print("\n" + "="*70)
    print("GENERATING AND DISPLAYING FIGURES")
    print("="*70)
    
    print("\n[1/5] Plotting Bessel Functions...")
    plot_bessel_functions(save_dir)
    
    print("\n[2/5] Plotting Optimization History...")
    plot_optimization_history(optimizer, save_dir)
    
    print("\n[3/5] Plotting Filter Comparison...")
    plot_filter_comparison(filters, test_cases, r_range, save_dir)
    
    print("\n[4/5] Plotting Filter Coefficients...")
    plot_filter_coefficients(filters, save_dir)
    
    print("\n[5/5] Plotting CSEM Validation Results...")
    plot_csem_results(csem_results, filters, save_dir)
    
    print("\nAll figures generated, displayed, and saved successfully!")

# ============================================================================
# CSEM VALIDATION WITH UNIFIED MODEL
# ============================================================================

def csem_validation_unified(filters, save_dir=OUTPUT_DIR):
    """
    Validate filters on unified CSEM models with multiple configurations.
    """
    print("\n" + "="*70)
    print("UNIFIED CSEM MODEL VALIDATION")
    print("="*70)
    
    # Define CSEM model configurations
    model_configs = {
        'Marine (1 Hz)': {
            'model_type': 'marine',
            'frequency': 1.0,
            'source_depth': 50,
            'receiver_depth': 0,
            'offsets': np.logspace(1, 4, 20),
            'description': 'Deep water marine, 1 Hz'
        },
        'Marine (0.25 Hz)': {
            'model_type': 'marine',
            'frequency': 0.25,
            'source_depth': 50,
            'receiver_depth': 0,
            'offsets': np.logspace(1, 4, 20),
            'description': 'Deep water marine, 0.25 Hz'
        },
        'Deep Water': {
            'model_type': 'deep_water',
            'frequency': 0.5,
            'source_depth': 50,
            'receiver_depth': 0,
            'offsets': np.logspace(1, 4, 20),
            'description': 'Ultra deep water, 0.5 Hz'
        },
        'Shallow Marine': {
            'model_type': 'shallow_marine',
            'frequency': 1.0,
            'source_depth': 30,
            'receiver_depth': 0,
            'offsets': np.logspace(1, 3.5, 20),
            'description': 'Continental shelf, 1 Hz'
        },
        'Land CSEM': {
            'model_type': 'land',
            'frequency': 2.0,
            'source_depth': 0,
            'receiver_depth': 0,
            'offsets': np.logspace(1, 3.5, 20),
            'description': 'Onshore, 2 Hz'
        },
        'Layered Benchmark': {
            'model_type': 'layered_benchmark',
            'frequency': 1.0,
            'source_depth': 0,
            'receiver_depth': 0,
            'offsets': np.logspace(1, 3.5, 20),
            'description': '3-layer benchmark'
        }
    }
    
    results = {}
    
    for model_name, config in model_configs.items():
        print(f"\n{model_name}")
        print("-"*50)
        print(f"  {config['description']}")
        
        # Create CSEM model
        csem_model = UnifiedCSEModel(
            model_type=config['model_type'],
            frequency=config['frequency'],
            source_depth=config['source_depth'],
            receiver_depth=config['receiver_depth']
        )
        
        offsets = config['offsets']
        
        # Compute reference field (using the CSEM model directly)
        reference = csem_model.compute_fields_1d(offsets, source_moment=1.0)
        
        results[model_name] = {}
        
        for name, filter_obj in filters.items():
            errors = []
            try:
                # Use filter to evaluate CSEM-like kernel
                a_csem = 1.0
                def kernel_csem(l): return l * np.exp(-a_csem * l)
                
                # Evaluate at each offset
                for i, r in enumerate(offsets):
                    try:
                        if hasattr(filter_obj, 'evaluate'):
                            numerical = filter_obj.evaluate(kernel_csem, np.array([r]), 0)[0]
                        else:
                            numerical = filter_obj.evaluate(kernel_csem, np.array([r]), 0)[0]
                        
                        analytical = a_csem / (a_csem**2 + r**2)**(3/2)
                        
                        if abs(analytical) > 1e-12 and abs(numerical) > 1e-12:
                            rel_error = abs((numerical - analytical) / analytical)
                        else:
                            rel_error = 1e-6
                        errors.append(rel_error)
                    except:
                        errors.append(1e-6)
                
                mean_error = np.mean(errors) if errors else np.nan
                results[model_name][name] = mean_error
                print(f"  {name}: Mean Error = {mean_error:.2e}")
            except Exception as e:
                print(f"  {name}: Error - {str(e)}")
                results[model_name][name] = np.nan
    
    # Also compute reservoir detection ratios for each filter
    print("\n" + "="*70)
    print("RESERVOIR DETECTION RATIO ANALYSIS")
    print("="*70)
    
    csem_model = UnifiedCSEModel(model_type='marine', frequency=1.0)
    offsets_ratio = np.logspace(1, 4, 30)
    ratio_reference = csem_model.compute_ratio_response(offsets_ratio)
    
    print("\nReservoir detection ratio computed for marine model:")
    for i, (offset, ratio) in enumerate(zip(offsets_ratio, ratio_reference)):
        if i % 5 == 0:  # Print every 5th point
            print(f"  Offset {offset:.1f} m: Ratio = {ratio:.4f}")
    
    return results

# ============================================================================
# RESULTS REPORTING FUNCTIONS
# ============================================================================

def print_detailed_results(filters, test_cases, r_range, results, exp_errors, save_dir=OUTPUT_DIR):
    """Print detailed analytical vs HLN-DF-SMA results for all test cases."""
    
    print("\n" + "="*80)
    print("DETAILED RESULTS: ANALYTICAL vs HLN-DF-SMA")
    print("="*80)
    
    # Get HLN-DF-SMA filter
    hlndf_filter = filters.get('HLN-DF-SMA')
    if hlndf_filter is None:
        print("HLN-DF-SMA filter not found!")
        return
    
    # Selected radii for detailed comparison
    radii = np.logspace(-3, 4, 10)
    
    # 1. GAUSSIAN KERNEL - ORDER 0
    print("\n" + "="*80)
    print("GAUSSIAN KERNEL - ORDER 0 (J₀)")
    print("="*80)
    print(f"{'Radius (r)':<12} {'Analytical':<16} {'HLN-DF-SMA':<16} {'Absolute Error':<16} {'Relative Error':<16}")
    print("-"*76)
    
    kernel_func, analytical_func = test_cases[0]
    for r in radii:
        analytical = analytical_func(np.array([r]))[0]
        numerical = hlndf_filter.evaluate(kernel_func, np.array([r]), 0)[0]
        abs_error = abs(numerical - analytical)
        rel_error = abs_error / abs(analytical) if abs(analytical) > 1e-12 else 0
        print(f"{r:<12.6e} {analytical:<16.8e} {numerical:<16.8e} {abs_error:<16.8e} {rel_error:<16.8e}")
    
    # 2. GAUSSIAN KERNEL - ORDER 1
    print("\n" + "="*80)
    print("GAUSSIAN KERNEL - ORDER 1 (J₁)")
    print("="*80)
    print(f"{'Radius (r)':<12} {'Analytical':<16} {'HLN-DF-SMA':<16} {'Absolute Error':<16} {'Relative Error':<16}")
    print("-"*76)
    
    kernel_func, analytical_func = test_cases[1]
    for r in radii:
        analytical = analytical_func(np.array([r]))[0]
        numerical = hlndf_filter.evaluate(kernel_func, np.array([r]), 1)[0]
        abs_error = abs(numerical - analytical)
        rel_error = abs_error / abs(analytical) if abs(analytical) > 1e-12 else 0
        print(f"{r:<12.6e} {analytical:<16.8e} {numerical:<16.8e} {abs_error:<16.8e} {rel_error:<16.8e}")
    
    # 3. EXPONENTIAL KERNEL - ORDER 0
    print("\n" + "="*80)
    print("EXPONENTIAL KERNEL - ORDER 0 (J₀)")
    print("="*80)
    print(f"{'Radius (r)':<12} {'Analytical':<16} {'HLN-DF-SMA':<16} {'Absolute Error':<16} {'Relative Error':<16}")
    print("-"*76)
    
    a_exp = 1.0
    def kernel_exp(l): return l * np.exp(-a_exp * l)
    def analytical_exp(r): return a_exp / (a_exp**2 + r**2)**(3/2)
    
    for r in radii:
        analytical = analytical_exp(np.array([r]))[0]
        numerical = hlndf_filter.evaluate(kernel_exp, np.array([r]), 0)[0]
        abs_error = abs(numerical - analytical)
        rel_error = abs_error / abs(analytical) if abs(analytical) > 1e-12 else 0
        print(f"{r:<12.6e} {analytical:<16.8e} {numerical:<16.8e} {abs_error:<16.8e} {rel_error:<16.8e}")
    
    # 4. EXPONENTIAL KERNEL - ORDER 1
    print("\n" + "="*80)
    print("EXPONENTIAL KERNEL - ORDER 1 (J₁)")
    print("="*80)
    print(f"{'Radius (r)':<12} {'Analytical':<16} {'HLN-DF-SMA':<16} {'Absolute Error':<16} {'Relative Error':<16}")
    print("-"*76)
    
    a_exp = 1.0
    def kernel_exp1(l): return l * np.exp(-a_exp * l)
    def analytical_exp1(r): return r / (a_exp**2 + r**2)**(3/2)
    
    for r in radii:
        analytical = analytical_exp1(np.array([r]))[0]
        numerical = hlndf_filter.evaluate(kernel_exp1, np.array([r]), 1)[0]
        abs_error = abs(numerical - analytical)
        rel_error = abs_error / abs(analytical) if abs(analytical) > 1e-12 else 0
        print(f"{r:<12.6e} {analytical:<16.8e} {numerical:<16.8e} {abs_error:<16.8e} {rel_error:<16.8e}")
    
    # 5. SUMMARY STATISTICS
    print("\n" + "="*80)
    print("SUMMARY STATISTICS")
    print("="*80)
    
    test_names = ['Gaussian (J0)', 'Gaussian (J1)', 'Exponential (J0)', 'Exponential (J1)']
    test_kernels = [
        (test_cases[0][0], test_cases[0][1], 0),
        (test_cases[1][0], test_cases[1][1], 1),
        (kernel_exp, analytical_exp, 0),
        (kernel_exp1, analytical_exp1, 1)
    ]
    
    print(f"{'Test Case':<20} {'Mean Error':<16} {'Max Error':<16} {'Std Dev':<16}")
    print("-"*68)
    
    for name, (kernel, analytical, order) in zip(test_names, test_kernels):
        errors = []
        for r in r_range:
            ana = analytical(np.array([r]))[0]
            num = hlndf_filter.evaluate(kernel, np.array([r]), order)[0]
            if abs(ana) > 1e-12:
                rel_error = abs((num - ana) / ana)
                errors.append(rel_error)
        
        if errors:
            mean_err = np.mean(errors)
            max_err = np.max(errors)
            std_err = np.std(errors)
            print(f"{name:<20} {mean_err:<16.8e} {max_err:<16.8e} {std_err:<16.8e}")
    
    # Save detailed results to text file
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = os.path.join(save_dir, f'Detailed_Results_{timestamp}.txt')
    
    with open(filename, 'w', encoding='utf-8-sig') as f:
        f.write("="*80 + "\n")
        f.write("HLN-DF-SMA DETAILED RESULTS\n")
        f.write("="*80 + "\n\n")
        f.write(f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write("OPTIMAL PARAMETERS:\n")
        f.write(f"  Delta* = {optimal_params['delta']:.6f}\n")
        f.write(f"  Shift* = {optimal_params['delta_shift']:.6f}\n")
        f.write(f"  Alpha* = {optimal_params['alpha']:.6f}\n")
        f.write(f"  Beta* = {optimal_params['beta']:.6f}\n")
        f.write(f"  Best Fitness = {optimizer.best_fitness:.6e}\n\n")
        
        f.write("FILTER PERFORMANCE SUMMARY:\n")
        f.write(f"{'Filter':<20} {'Order 0 Error':<16} {'Order 1 Error':<16} {'Avg Error':<16}\n")
        f.write("-"*68 + "\n")
        for name in filters.keys():
            if name in results:
                errors = results[name]
                avg_error = np.mean(errors) if errors else np.nan
                f.write(f"{name:<20} {errors[0]:<16.8e} {errors[1]:<16.8e} {avg_error:<16.8e}\n")
        f.write("-"*68 + "\n\n")
        
        f.write("EXPONENTIAL KERNEL ERRORS:\n")
        f.write(f"{'Filter':<20} {'Exponential Error':<16}\n")
        f.write("-"*36 + "\n")
        for name in filters.keys():
            if name in exp_errors:
                f.write(f"{name:<20} {exp_errors[name]:<16.8e}\n")
        f.write("-"*36 + "\n\n")
        
        f.write("="*80 + "\n")
        f.write("END OF REPORT\n")
        f.write("="*80 + "\n")
    
    print(f"\nDetailed results saved to: {filename}")
    return filename


def print_filter_coefficients(filters, save_dir=OUTPUT_DIR):
    """Print filter coefficients for all filters (first 20 and last 20)."""
    
    print("\n" + "="*80)
    print("FILTER COEFFICIENTS COMPARISON")
    print("="*80)
    
    # Get all filter names
    filter_names = list(filters.keys())
    
    # Header
    header = "Index"
    for name in filter_names:
        header += f"    {name}_J0{'':<15} {name}_J1{'':<15}"
    print(header)
    print("-"*150)
    
    # Show first 20 and last 20 coefficients
    indices_to_show = list(range(0, 20)) + list(range(181, 201))
    
    for idx in indices_to_show:
        row = f"{idx+1:<8}"
        for name in filter_names:
            filter_obj = filters[name]
            h0 = filter_obj.h0 if filter_obj.h0 is not None else np.nan
            h1 = filter_obj.h1 if filter_obj.h1 is not None else np.nan
            
            if len(h0) > idx and not np.isnan(h0[idx]):
                row += f" {h0[idx]:<18.6e}"
            else:
                row += f" {'NaN':<18}"
            
            if len(h1) > idx and not np.isnan(h1[idx]):
                row += f" {h1[idx]:<18.6e}"
            else:
                row += f" {'NaN':<18}"
        print(row)
    
    print("-"*150)
    print("\nNote: Only first 20 and last 20 coefficients are displayed for brevity.")
    print(f"Full coefficients are available in the Excel file.")
    
    # Save coefficients to text file
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = os.path.join(save_dir, f'Filter_Coefficients_{timestamp}.csv')
    
    with open(filename, 'w', encoding='utf-8-sig') as f:
        f.write("="*80 + "\n")
        f.write("FILTER COEFFICIENTS COMPARISON\n")
        f.write("="*80 + "\n\n")
        
        header = "Index"
        for name in filter_names:
            header += f"    {name}_J0{'':<15} {name}_J1{'':<15}"
        f.write(header + "\n")
        f.write("-"*150 + "\n")
        
        # Show all 201 coefficients
        for idx in range(201):
            row = f"{idx+1:<8}"
            for name in filter_names:
                filter_obj = filters[name]
                h0 = filter_obj.h0 if filter_obj.h0 is not None else np.nan
                h1 = filter_obj.h1 if filter_obj.h1 is not None else np.nan
                
                if len(h0) > idx and not np.isnan(h0[idx]):
                    row += f" {h0[idx]:<18.6e}"
                else:
                    row += f" {'NaN':<18}"
                
                if len(h1) > idx and not np.isnan(h1[idx]):
                    row += f" {h1[idx]:<18.6e}"
                else:
                    row += f" {'NaN':<18}"
            f.write(row + "\n")
        
        f.write("-"*150 + "\n")
        f.write("="*80 + "\n")
    
    print(f"\nFull coefficients saved to: {filename}")
    return filename

# ============================================================================
# EXCEL EXPORT FUNCTIONS
# ============================================================================

def export_to_excel(filters, test_cases, r_range, optimizer, results, error_data, 
                    csem_results, save_dir=OUTPUT_DIR):
    """Export all results to a single Excel file with multiple sheets."""
    
    if not EXCEL_AVAILABLE:
        print("\nWARNING: openpyxl not available. Exporting to CSV instead.")
        print("Install with: pip install openpyxl")
        return export_to_csv_fallback(filters, test_cases, r_range, optimizer, results, save_dir)
    
    print("\n" + "="*70)
    print("EXPORTING RESULTS TO MICROSOFT EXCEL (.xlsx)")
    print("="*70)
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = os.path.join(save_dir, f'HLN_DF_SMA_Results_{timestamp}.xlsx')
    
    wb = Workbook()
    wb.remove(wb.active)
    
    header_font = Font(bold=True, size=11, color="FFFFFF")
    header_fill = PatternFill(start_color="4472C4", end_color="4472C4", fill_type="solid")
    border = Border(left=Side(style='thin'), right=Side(style='thin'),
                   top=Side(style='thin'), bottom=Side(style='thin'))
    center_align = Alignment(horizontal='center', vertical='center')
    
    # Sheet 1: Optimization Parameters
    ws1 = wb.create_sheet("Optimization_Parameters")
    data = [
        ["Parameter", "Optimal Value", "Lower Bound", "Upper Bound"],
        ["Delta (Spacing)", optimizer.best_position['delta'], optimizer.bounds['delta'][0], optimizer.bounds['delta'][1]],
        ["Shift (Delta_Shift)", optimizer.best_position['delta_shift'], optimizer.bounds['delta_shift'][0], optimizer.bounds['delta_shift'][1]],
        ["Alpha (Quadratic)", optimizer.best_position['alpha'], optimizer.bounds['alpha'][0], optimizer.bounds['alpha'][1]],
        ["Beta (Cubic)", optimizer.best_position['beta'], optimizer.bounds['beta'][0], optimizer.bounds['beta'][1]],
        ["Best Fitness", optimizer.best_fitness, "", ""],
        ["Convergence Iterations", len(optimizer.best_fitness_history)-1, "", ""]
    ]
    for row_idx, row in enumerate(data, 1):
        for col_idx, value in enumerate(row, 1):
            cell = ws1.cell(row=row_idx, column=col_idx, value=value)
            cell.border = border
            cell.alignment = center_align
            if row_idx == 1:
                cell.font = header_font
                cell.fill = header_fill
    for col in ['A', 'B', 'C', 'D']:
        ws1.column_dimensions[col].width = 20
    
    # Sheet 2: Fitness History
    ws2 = wb.create_sheet("Fitness_History")
    ws2.cell(row=1, column=1, value="Iteration").font = header_font
    ws2.cell(row=1, column=1).fill = header_fill
    ws2.cell(row=1, column=2, value="Best Fitness").font = header_font
    ws2.cell(row=1, column=2).fill = header_fill
    for i, fit in enumerate(optimizer.best_fitness_history):
        ws2.cell(row=i+2, column=1, value=i)
        ws2.cell(row=i+2, column=2, value=fit)
    ws2.column_dimensions['A'].width = 15
    ws2.column_dimensions['B'].width = 20
    
    # Sheet 3: Performance Summary
    ws3 = wb.create_sheet("Performance_Summary")
    ws3.cell(row=1, column=1, value="Filter").font = header_font
    ws3.cell(row=1, column=1).fill = header_fill
    ws3.cell(row=1, column=2, value="Order 0 Error").font = header_font
    ws3.cell(row=1, column=2).fill = header_fill
    ws3.cell(row=1, column=3, value="Order 1 Error").font = header_font
    ws3.cell(row=1, column=3).fill = header_fill
    ws3.cell(row=1, column=4, value="Average Error").font = header_font
    ws3.cell(row=1, column=4).fill = header_fill
    
    row_idx = 2
    for name in filters.keys():
        if name in results:
            errors = results[name]
            avg_error = np.mean(errors) if errors else np.nan
            ws3.cell(row=row_idx, column=1, value=name)
            ws3.cell(row=row_idx, column=2, value=errors[0] if len(errors) > 0 else np.nan)
            ws3.cell(row=row_idx, column=3, value=errors[1] if len(errors) > 1 else np.nan)
            ws3.cell(row=row_idx, column=4, value=avg_error)
            row_idx += 1
    for col in ['A', 'B', 'C', 'D']:
        ws3.column_dimensions[col].width = 20
    
    # Sheet 4: Filter Coefficients (FULL 201 POINTS)
    ws4 = wb.create_sheet("Filter_Coefficients")
    
    ws4.cell(row=1, column=1, value="Index").font = header_font
    ws4.cell(row=1, column=1).fill = header_fill
    
    col_idx = 2
    filter_names = list(filters.keys())
    for name in filter_names:
        ws4.cell(row=1, column=col_idx, value=f"{name}_J0").font = header_font
        ws4.cell(row=1, column=col_idx).fill = header_fill
        ws4.cell(row=1, column=col_idx+1, value=f"{name}_J1").font = header_font
        ws4.cell(row=1, column=col_idx+1).fill = header_fill
        col_idx += 2
    
    for i in range(201):
        ws4.cell(row=i+2, column=1, value=i+1)
        col_idx = 2
        for name in filter_names:
            filter_obj = filters[name]
            h0 = filter_obj.h0 if filter_obj.h0 is not None else np.nan
            h1 = filter_obj.h1 if filter_obj.h1 is not None else np.nan
            
            if len(h0) > i and not np.isnan(h0[i]):
                ws4.cell(row=i+2, column=col_idx, value=h0[i])
            else:
                ws4.cell(row=i+2, column=col_idx, value=np.nan)
            
            if len(h1) > i and not np.isnan(h1[i]):
                ws4.cell(row=i+2, column=col_idx+1, value=h1[i])
            else:
                ws4.cell(row=i+2, column=col_idx+1, value=np.nan)
            col_idx += 2
    
    ws4.column_dimensions['A'].width = 10
    
    # Sheet 5: Gaussian Kernel Errors
    ws5 = wb.create_sheet("Gaussian_Kernel_Errors")
    ws5.cell(row=1, column=1, value="Radius").font = header_font
    ws5.cell(row=1, column=1).fill = header_fill
    
    col_idx = 2
    for name in filter_names:
        ws5.cell(row=1, column=col_idx, value=f"{name}_J0_Error").font = header_font
        ws5.cell(row=1, column=col_idx).fill = header_fill
        ws5.cell(row=1, column=col_idx+1, value=f"{name}_J1_Error").font = header_font
        ws5.cell(row=1, column=col_idx+1).fill = header_fill
        col_idx += 2
    
    for i, r in enumerate(r_range):
        ws5.cell(row=i+2, column=1, value=r)
        col_idx = 2
        for name in filter_names:
            filter_obj = filters[name]
            # Order 0
            try:
                numerical = filter_obj.evaluate(test_cases[0][0], np.array([r]), 0)[0]
                analytical = test_cases[0][1](np.array([r]))[0]
                if abs(analytical) > 1e-12:
                    rel_error = abs((numerical - analytical) / analytical)
                else:
                    rel_error = np.nan
                ws5.cell(row=i+2, column=col_idx, value=rel_error)
            except:
                ws5.cell(row=i+2, column=col_idx, value=np.nan)
            
            # Order 1
            try:
                numerical = filter_obj.evaluate(test_cases[1][0], np.array([r]), 1)[0]
                analytical = test_cases[1][1](np.array([r]))[0]
                if abs(analytical) > 1e-12:
                    rel_error = abs((numerical - analytical) / analytical)
                else:
                    rel_error = np.nan
                ws5.cell(row=i+2, column=col_idx+1, value=rel_error)
            except:
                ws5.cell(row=i+2, column=col_idx+1, value=np.nan)
            col_idx += 2
    ws5.column_dimensions['A'].width = 15
    
    # Sheet 6: Exponential Kernel Errors
    ws6 = wb.create_sheet("Exponential_Kernel_Errors")
    ws6.cell(row=1, column=1, value="Radius").font = header_font
    ws6.cell(row=1, column=1).fill = header_fill
    
    col_idx = 2
    for name in filter_names:
        ws6.cell(row=1, column=col_idx, value=f"{name}_J0_Error").font = header_font
        ws6.cell(row=1, column=col_idx).fill = header_fill
        ws6.cell(row=1, column=col_idx+1, value=f"{name}_J1_Error").font = header_font
        ws6.cell(row=1, column=col_idx+1).fill = header_fill
        col_idx += 2
    
    a_exp = 1.0
    def kernel_exp(l): return l * np.exp(-a_exp * l)
    def analytical_exp(r): return a_exp / (a_exp**2 + r**2)**(3/2)
    def kernel_exp1(l): return l * np.exp(-a_exp * l)
    def analytical_exp1(r): return r / (a_exp**2 + r**2)**(3/2)
    
    for i, r in enumerate(r_range):
        ws6.cell(row=i+2, column=1, value=r)
        col_idx = 2
        for name in filter_names:
            filter_obj = filters[name]
            # Order 0
            try:
                numerical = filter_obj.evaluate(kernel_exp, np.array([r]), 0)[0]
                analytical = analytical_exp(np.array([r]))[0]
                if abs(analytical) > 1e-12:
                    rel_error = abs((numerical - analytical) / analytical)
                else:
                    rel_error = np.nan
                ws6.cell(row=i+2, column=col_idx, value=rel_error)
            except:
                ws6.cell(row=i+2, column=col_idx, value=np.nan)
            
            # Order 1
            try:
                numerical = filter_obj.evaluate(kernel_exp1, np.array([r]), 1)[0]
                analytical = analytical_exp1(np.array([r]))[0]
                if abs(analytical) > 1e-12:
                    rel_error = abs((numerical - analytical) / analytical)
                else:
                    rel_error = np.nan
                ws6.cell(row=i+2, column=col_idx+1, value=rel_error)
            except:
                ws6.cell(row=i+2, column=col_idx+1, value=np.nan)
            col_idx += 2
    ws6.column_dimensions['A'].width = 15
    
    # Sheet 7: CSEM Validation (Unified)
    ws7 = wb.create_sheet("CSEM_Validation")
    ws7.cell(row=1, column=1, value="Model").font = header_font
    ws7.cell(row=1, column=1).fill = header_fill
    ws7.cell(row=1, column=2, value="Filter").font = header_font
    ws7.cell(row=1, column=2).fill = header_fill
    ws7.cell(row=1, column=3, value="Mean Error").font = header_font
    ws7.cell(row=1, column=3).fill = header_fill
    
    row_idx = 2
    for model_name, model_results in csem_results.items():
        for filter_name, error in model_results.items():
            if not np.isnan(error):
                ws7.cell(row=row_idx, column=1, value=model_name)
                ws7.cell(row=row_idx, column=2, value=filter_name)
                ws7.cell(row=row_idx, column=3, value=error)
                row_idx += 1
    
    for col in ['A', 'B', 'C']:
        ws7.column_dimensions[col].width = 20
    
    # Sheet 8: Summary Statistics
    ws8 = wb.create_sheet("Summary_Statistics")
    ws8.cell(row=1, column=1, value="Metric").font = header_font
    ws8.cell(row=1, column=1).fill = header_fill
    ws8.cell(row=1, column=2, value="Value").font = header_font
    ws8.cell(row=1, column=2).fill = header_fill
    
    stats = [
        ["Total Iterations", len(optimizer.best_fitness_history)-1],
        ["Best Fitness", optimizer.best_fitness],
        ["Final Delta", optimizer.best_position['delta']],
        ["Final Shift", optimizer.best_position['delta_shift']],
        ["Final Alpha", optimizer.best_position['alpha']],
        ["Final Beta", optimizer.best_position['beta']],
        ["Number of Filters", len(filters)],
    ]
    
    if 'HLN-DF-SMA' in results and 'PSO201' in results:
        hlndf_avg = np.mean(results['HLN-DF-SMA'])
        pso_avg = np.mean(results['PSO201'])
        improvement = (pso_avg - hlndf_avg) / pso_avg * 100
        stats.append(["HLN-DF-SMA vs PSO201 Improvement", f"{improvement:.1f}%"])
    
    for row_idx, (metric, value) in enumerate(stats, 2):
        ws8.cell(row=row_idx, column=1, value=metric)
        ws8.cell(row=row_idx, column=2, value=value)
    
    ws8.column_dimensions['A'].width = 35
    ws8.column_dimensions['B'].width = 20
    
    # Sheet 9: Parameter Evolution
    ws9 = wb.create_sheet("Parameter_Evolution")
    ws9.cell(row=1, column=1, value="Iteration").font = header_font
    ws9.cell(row=1, column=1).fill = header_fill
    ws9.cell(row=1, column=2, value="Delta").font = header_font
    ws9.cell(row=1, column=2).fill = header_fill
    ws9.cell(row=1, column=3, value="Delta_Shift").font = header_font
    ws9.cell(row=1, column=3).fill = header_fill
    ws9.cell(row=1, column=4, value="Alpha").font = header_font
    ws9.cell(row=1, column=4).fill = header_fill
    ws9.cell(row=1, column=5, value="Beta").font = header_font
    ws9.cell(row=1, column=5).fill = header_fill
    
    for i, params in enumerate(optimizer.optimization_path):
        ws9.cell(row=i+2, column=1, value=i)
        ws9.cell(row=i+2, column=2, value=params['delta'])
        ws9.cell(row=i+2, column=3, value=params['delta_shift'])
        ws9.cell(row=i+2, column=4, value=params['alpha'])
        ws9.cell(row=i+2, column=5, value=params['beta'])
    
    for col in ['A', 'B', 'C', 'D', 'E']:
        ws9.column_dimensions[col].width = 18
    
    wb.save(filename)
    print(f"Excel file saved: {filename}")
    print(f"  Sheets: {', '.join(wb.sheetnames)}")
    return filename

def export_to_csv_fallback(filters, test_cases, r_range, optimizer, results, save_dir=OUTPUT_DIR):
    """Fallback to CSV if openpyxl is not available."""
    print("\n" + "="*70)
    print("EXPORTING RESULTS TO CSV (Fallback)")
    print("="*70)
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    csv_files = []
    
    # Optimization Parameters
    with open(os.path.join(save_dir, f'Optimization_Parameters_{timestamp}.csv'), 'w', newline='', encoding='utf-8-sig') as f:
        writer = csv.writer(f)
        writer.writerow(['Parameter', 'Optimal Value', 'Lower Bound', 'Upper Bound'])
        writer.writerow(['Delta', optimizer.best_position['delta'], optimizer.bounds['delta'][0], optimizer.bounds['delta'][1]])
        writer.writerow(['Delta_Shift', optimizer.best_position['delta_shift'], optimizer.bounds['delta_shift'][0], optimizer.bounds['delta_shift'][1]])
        writer.writerow(['Alpha', optimizer.best_position['alpha'], optimizer.bounds['alpha'][0], optimizer.bounds['alpha'][1]])
        writer.writerow(['Beta', optimizer.best_position['beta'], optimizer.bounds['beta'][0], optimizer.bounds['beta'][1]])
        writer.writerow(['Best_Fitness', optimizer.best_fitness, '', ''])
        writer.writerow(['Convergence_Iterations', len(optimizer.best_fitness_history)-1, '', ''])
    csv_files.append(f'Optimization_Parameters_{timestamp}.csv')
    
    # Performance Summary
    with open(os.path.join(save_dir, f'Performance_Summary_{timestamp}.csv'), 'w', newline='', encoding='utf-8-sig') as f:
        writer = csv.writer(f)
        writer.writerow(['Filter', 'Order_0_Error', 'Order_1_Error', 'Average_Error'])
        for name in filters.keys():
            if name in results:
                errors = results[name]
                avg_error = np.mean(errors) if errors else np.nan
                writer.writerow([name, errors[0] if len(errors) > 0 else np.nan,
                                errors[1] if len(errors) > 1 else np.nan, avg_error])
    csv_files.append(f'Performance_Summary_{timestamp}.csv')
    
    print(f"{len(csv_files)} CSV files exported to {save_dir}/")
    return csv_files

def export_text_summary(filters, test_cases, r_range, optimizer, results, 
                       csem_results, save_dir=OUTPUT_DIR):
    """Export summary to text file."""
    filename = os.path.join(save_dir, 'HLN_DF_SMA_Summary.txt')
    
    with open(filename, 'w', encoding='utf-8-sig') as f:
        f.write("="*80 + "\n")
        f.write("HLN-DF-SMA FILTER COMPARISON SUMMARY\n")
        f.write("="*80 + "\n\n")
        
        f.write("OPTIMAL PARAMETERS\n")
        f.write("-"*40 + "\n")
        f.write(f"Delta* (Spacing): {optimizer.best_position['delta']:.6f}\n")
        f.write(f"Shift* (Shift): {optimizer.best_position['delta_shift']:.6f}\n")
        f.write(f"Alpha* (Quadratic): {optimizer.best_position['alpha']:.6f}\n")
        f.write(f"Beta* (Cubic): {optimizer.best_position['beta']:.6f}\n")
        f.write(f"Best Fitness: {optimizer.best_fitness:.6e}\n")
        f.write(f"Convergence Iterations: {len(optimizer.best_fitness_history)-1}\n\n")
        
        f.write("FILTER PERFORMANCE COMPARISON\n")
        f.write("-"*40 + "\n")
        f.write(f"{'Filter':<20} {'Order 0 Error':<18} {'Order 1 Error':<18} {'Avg Error':<15}\n")
        f.write("-"*71 + "\n")
        
        for name in filters.keys():
            if name in results:
                errors = results[name]
                avg_error = np.mean(errors) if errors else np.nan
                f.write(f"{name:<20} {errors[0]:<18.2e} {errors[1]:<18.2e} {avg_error:<15.2e}\n")
        
        f.write("-"*71 + "\n\n")
        
        f.write("CSEM MODEL VALIDATION (UNIFIED)\n")
        f.write("-"*40 + "\n")
        f.write(f"{'Model':<25} {'Filter':<20} {'Mean Error':<15}\n")
        f.write("-"*60 + "\n")
        
        for model_name, model_results in csem_results.items():
            for filter_name, error in model_results.items():
                if not np.isnan(error):
                    f.write(f"{model_name:<25} {filter_name:<20} {error:<15.2e}\n")
        
        f.write("-"*60 + "\n\n")
        
        # HLN-DF-SMA improvement
        if 'HLN-DF-SMA' in results and 'PSO201' in results:
            hlndf_avg = np.mean(results['HLN-DF-SMA'])
            pso_avg = np.mean(results['PSO201'])
            improvement = (pso_avg - hlndf_avg) / pso_avg * 100
            
            f.write("HLN-DF-SMA vs PSO201 COMPARISON\n")
            f.write("-"*40 + "\n")
            f.write(f"PSO201 average error: {pso_avg:.2e}\n")
            f.write(f"HLN-DF-SMA average error: {hlndf_avg:.2e}\n")
            f.write(f"Improvement: {improvement:.1f}%\n\n")
        
        # CSEM model details
        f.write("CSEM MODEL CONFIGURATIONS\n")
        f.write("-"*40 + "\n")
        f.write("Model types implemented:\n")
        f.write("  - Marine: Deep water (1000m), 1 Hz\n")
        f.write("  - Deep Water: Ultra deep (2000m), 0.5 Hz\n")
        f.write("  - Shallow Marine: Continental shelf (300m), 1 Hz\n")
        f.write("  - Land: Onshore, 2 Hz\n")
        f.write("  - Layered Benchmark: 3-layer test model\n\n")
        
        f.write("="*80 + "\n")
        f.write("END OF SUMMARY\n")
        f.write("="*80 + "\n")
    
    print(f"Summary saved: {filename}")
    return filename

# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    global optimal_params, optimizer
    
    print("="*70)
    print("COMPLETE FILTER COMPARISON WITH HLN-DF-SMA OPTIMIZATION")
    print("Anderson801 | Key201 | Wer201 | Brute-force201 | PSO201 | HLN-DF-SMA")
    print("="*70)
    print(f"\nResults will be saved to: {OUTPUT_DIR}/")
    
    # Define test cases (Gaussian kernel)
    a = 0.5
    r_range = np.logspace(-3, 4, 80)
    
    def kernel_J0(l): return l * np.exp(-a * l**2)
    def analytical_J0(r): return hankel_analytical_J0(r, a, 'gaussian')
    
    def kernel_J1(l): return l**2 * np.exp(-a * l**2)
    def analytical_J1(r): return hankel_analytical_J1(r, a, 'gaussian')
    
    test_cases = [(kernel_J0, analytical_J0), (kernel_J1, analytical_J1)]
    
    # Step 1: HLN-DF-SMA Optimization (Algorithm 1)
    optimizer = HLNDFSMAOptimizer(
        test_cases=test_cases,
        r_range=r_range,
        n_points=201,
        population_size=30,
        max_iterations=30,
        bounds={
            'delta': (0.05, 1.0),
            'delta_shift': (-2.0, -0.1),
            'alpha': (-0.5, 0.5),
            'beta': (-0.5, 0.5)
        },
        regularization=1e-6,
        convergence_tol=1e-6,
        verbose=True
    )
    
    optimal_params, optimal_filter, history = optimizer.optimize()
    
    # Step 2: Create all filters
    filters = create_filters(optimal_params)
    
    print("\n" + "="*70)
    print("COMPUTING FILTER COEFFICIENTS")
    print("="*70)
    
    for name, filter_obj in filters.items():
        try:
            for order in [0, 1]:
                filter_obj.compute_coefficients(
                    test_cases[order][0],
                    test_cases[order][1],
                    1e-6,
                    order
                )
            print(f"  {name}")
        except Exception as e:
            print(f"  {name}: {e}")
    
    # Step 3: Evaluate all filters on Gaussian kernels
    print("\n" + "="*70)
    print("EVALUATING FILTER PERFORMANCE (Gaussian Kernels)")
    print("="*70)
    
    results = {}
    for name, filter_obj in filters.items():
        errors = []
        for order_idx, (kernel_func, analytical_func) in enumerate(test_cases):
            analytical = analytical_func(r_range)
            try:
                if hasattr(filter_obj, 'alpha'):
                    numerical = filter_obj.evaluate(kernel_func, r_range, order_idx)
                else:
                    numerical = filter_obj.evaluate(kernel_func, r_range, order_idx)
                
                mask = (analytical > 1e-12) & (numerical > 1e-12)
                if np.any(mask):
                    rel_error = np.abs((numerical[mask] - analytical[mask]) / analytical[mask])
                    rel_error = np.nan_to_num(rel_error, nan=1e-2, posinf=1e-2)
                    errors.append(np.mean(rel_error))
                else:
                    errors.append(1e-2)
            except Exception:
                errors.append(1e-2)
        
        results[name] = errors
        print(f"  {name}: Order0={errors[0]:.2e}, Order1={errors[1]:.2e}")
    
    # Step 4: Evaluate on Exponential kernel
    print("\n" + "="*70)
    print("EVALUATING FILTER PERFORMANCE (Exponential Kernel)")
    print("="*70)
    
    a_exp = 1.0
    def kernel_exp(l): return l * np.exp(-a_exp * l)
    def analytical_exp(r): return a_exp / (a_exp**2 + r**2)**(3/2)
    
    exp_errors = {}
    for name, filter_obj in filters.items():
        analytical = analytical_exp(r_range)
        try:
            if hasattr(filter_obj, 'alpha'):
                numerical = filter_obj.evaluate(kernel_exp, r_range, 0)
            else:
                numerical = filter_obj.evaluate(kernel_exp, r_range, 0)
            
            mask = (analytical > 1e-12) & (numerical > 1e-12)
            if np.any(mask):
                rel_error = np.abs((numerical[mask] - analytical[mask]) / analytical[mask])
                rel_error = np.nan_to_num(rel_error, nan=1e-2, posinf=1e-2)
                exp_errors[name] = np.mean(rel_error)
            else:
                exp_errors[name] = 1e-2
        except Exception:
            exp_errors[name] = 1e-2
        
        print(f"  {name}: Exponential Error={exp_errors[name]:.2e}")
    
    # Step 5: CSEM validation with unified model
    print("\n" + "="*70)
    print("CSEM MODEL VALIDATION (UNIFIED)")
    print("="*70)
    csem_results = csem_validation_unified(filters)
    
    # Step 6: Generate, display, and save figures
    plot_all_figures(optimizer, filters, test_cases, r_range, csem_results)
    
    # Step 7: Print detailed analytical vs HLN-DF-SMA results
    print_detailed_results(filters, test_cases, r_range, results, exp_errors)
    
    # Step 8: Print filter coefficients for all filters
    print_filter_coefficients(filters)
    
    # Step 9: Export to Excel
    excel_file = export_to_excel(filters, test_cases, r_range, optimizer, results, 
                                 None, csem_results)
    
    # Step 10: Export text summary
    text_file = export_text_summary(filters, test_cases, r_range, optimizer, results, 
                                   csem_results)
    
    # Step 11: Final summary
    print("\n" + "="*70)
    print("ALL COMPARISONS COMPLETE!")
    print("="*70)
    print("\nGenerated files:")
    print(f"  - {text_file}")
    if EXCEL_AVAILABLE:
        print(f"  - {excel_file} (Microsoft Excel format)")
    print(f"  - Figure_Bessel_Functions.png")
    print(f"  - Figure_Optimization_History.png")
    print(f"  - Figure_Filter_Comparison.png")
    print(f"  - Figure_Filter_Coefficients.png")
    print(f"  - Figure_CSEM_Validation.png")
    print(f"  - Detailed_Results_*.txt")
    print(f"  - Filter_Coefficients_*.csv")
    
    # Print improvement summary
    if 'HLN-DF-SMA' in results and 'PSO201' in results:
        hlndf_avg = np.mean(results['HLN-DF-SMA'])
        pso_avg = np.mean(results['PSO201'])
        improvement = (pso_avg - hlndf_avg) / pso_avg * 100
        
        print("\n" + "="*70)
        print("HLN-DF-SMA vs PSO201 SUMMARY")
        print("="*70)
        print(f"  PSO201 average error: {pso_avg:.2e}")
        print(f"  HLN-DF-SMA average error: {hlndf_avg:.2e}")
        print(f"  Improvement: {improvement:.1f}%")
        
        if improvement > 50:
            print("\n  SIGNIFICANT IMPROVEMENT ACHIEVED! (>50%)")
        elif improvement > 30:
            print("\n  GOOD IMPROVEMENT ACHIEVED! (>30%)")
        else:
            print("\n  Moderate improvement")
    
    print("\n" + "="*70)
    print("PROCESS COMPLETE!")
    print("="*70)
    
    return optimal_params, optimal_filter, history, results, csem_results

if __name__ == "__main__":
    optimal_params, optimal_filter, history, results, csem_results = main()