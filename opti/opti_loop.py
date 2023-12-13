import numpy as np
import subprocess
from opti_analysis import analyze_data
import pandas as pd
from scipy.optimize import minimize
import os
import csv


def generate_random_parameters():
    # Generate 12 random parameters
    # Adjust the range according to your specific requirements
    return np.random.uniform(low=0.0001, high=1000.0, size=12)

def call_matlab_function(params):
    # Convert parameters to a string format accepted by MATLAB
    param_str = ', '.join(map(str, params))
    matlab_command = f"matlab -batch \"generateControlFiles({param_str})\""
    subprocess.run(matlab_command, shell=True)
    
def call_rust_function():
    # Navigate to the Rust project directory
    rust_project_dir = "/home/bholder/gigabit_code"

    # Rust command to run
    rust_command = "cargo run --release -p sim"

    # Execute the command in the Rust project directory
    subprocess.run(rust_command, shell=True, cwd=rust_project_dir)


def analyze_output():
    file_path = '/home/bholder/data/out.csv'
    results = analyze_data(file_path)
    if os.path.exists('/home/bholder/data/out.csv'):
        os.remove('/home/bholder/data/out.csv')
    return results

# def optimize(max_iterations):
#     for iteration in range(max_iterations):
#         params = generate_random_parameters()
#         matlab_result = call_matlab_function(params)
#         rust_result = call_rust_function()
#         analysis = analyze_output(rust_result)
#         # Use the analysis to inform the next set of parameters
#         # This is where your optimization algorithm will come into play
        


# def optimize(max_iterations, output_csv, save_every):
#     # Initialize a DataFrame to store the results
#     results_df = pd.DataFrame()

#     for iteration in range(max_iterations):
#         params = generate_random_parameters()
#         call_matlab_function(params)
#         call_rust_function()
#         analysis_results = analyze_output()

#         # Store parameters and results in the DataFrame
#         iteration_data = {
#             'iteration': iteration,
#             **{f'param_{i+1}': param for i, param in enumerate(params)},
#             **{f'{col}_peak_to_peak': analysis_results[col]['peak_to_peak'] for col in analysis_results},
#             **{f'{col}_sigma': analysis_results[col]['sigma'] for col in analysis_results}
#         }
#         results_df = results_df.append(iteration_data, ignore_index=True)

#         # Check if it's time to save the results to CSV
#         if iteration % save_every == 0:
#             if iteration == save_every or not os.path.exists(output_csv):
#                 # If this is the first save or the file doesn't exist, include column names
#                 results_df.to_csv(output_csv, index=False, header=True)
#             else:
#                 # Append to the existing CSV file without writing column names again
#                 results_df.to_csv(output_csv, index=False, header=False, mode='a')

#             print(f'Saved results up to iteration {iteration} to {output_csv}')

#     # Write the remaining results to the CSV file
#     results_df.to_csv(output_csv, index=False)
#     print(f'Saved final results to {output_csv}')
    
#     if os.path.exists('/home/bholder/data/out.csv'):
#         os.remove('/home/bholder/data/out.csv')

# Define your custom objective function
def objective_function(params):
    # Call MATLAB, Rust, and analyze_output functions here with 'params'
    # Calculate and return the value of 'sigma' to minimize
    # Example: (replace with actual calls)
    call_matlab_function(params)
    call_rust_function()
    analysis_results = analyze_output()
    
    # Extract sigma values for 'ra', 'dec', and 'fr'
    sig_ra = analysis_results['ra']['sigma']
    sig_dec = analysis_results['dec']['sigma']
    sig_fr = analysis_results['fr']['sigma']
    
    # You can perform calculations with sig_ra, sig_dec, and sig_fr here if needed
    # Example calculation:
    sigma = 100000*(sig_ra + sig_dec + sig_fr)
    
     # Save parameters and sigmas to a results file
    results_file = 'results.csv'

    if not os.path.exists(results_file):
        # Create a new results file with headers
        with open(results_file, 'w', newline='') as file:
            writer = csv.writer(file)
            # Write headers
            headers = ['P1', 'P2','P3','P4','P5','P6','P7','P8','P9','P10','P11', 'Param_12', 'Sig_RA', 'Sig_Dec', 'Sig_FR']
            writer.writerow(headers)
    
    # Append parameters and sigmas as a new row in the results file
    with open(results_file, 'a', newline='') as file:
        writer = csv.writer(file)
        # Concatenate the parameters and sigmas individually
        
        params_list = params.tolist()
        row = params_list + [sig_ra, sig_dec, sig_fr]
        writer.writerow(row)
    
    
    
    return sigma

# Initialize starting parameters
initial_params = [450.0, 800.0, 1000.0, 1500.0, 7500.0, 1000.0, 0.005, 0.4, 0.1, 0.01, 0.001, 0.1]  # Replace with your initial parameter values

param_bounds = [
    (0.0001, 10000),  # Parameter 1
    (0.0001, 10000),  # Parameter 2
    (0.0001, 10000),  # Parameter 3
    (0.0001, 10000),  # Parameter 4
    (0.0001, 10000),  # Parameter 5
    (0.0001, 10000),  # Parameter 6
    (0.0001, 10000),  # Parameter 7
    (0.0001, 10000),  # Parameter 8
    (0.0001, 10000),  # Parameter 9
    (0.0001, 10000),  # Parameter 10
    (0.0001, 10000),  # Parameter 11
    (0.0001, 10000)   # Parameter 12
]

tolerance = 1e-4

# Use SciPy's minimize function to optimize
result = minimize(objective_function, initial_params, method='nelder-mead', bounds=param_bounds, tol=tolerance)

# Extract optimized parameters
optimized_params = result.x

# Store results in a DataFrame
results_df = pd.DataFrame({'Parameter': optimized_params, 'Sigma': result.fun})

# Write results to a CSV file
results_df.to_csv('optimization_results.csv', index=False)