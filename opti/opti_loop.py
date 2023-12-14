import numpy as np
import subprocess
from opti_analysis import analyze_data
import pandas as pd
from scipy.optimize import minimize
import os
import csv
import time
import json

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

    with open(os.devnull, 'w') as nullfile:
        subprocess.run(rust_command, shell=True, cwd=rust_project_dir, stdout=nullfile, stderr=nullfile)

# Call the function


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
def update_control_parameters(gains):
    kp_roll = gains[0]
    ki_roll = gains[1]
    kp_pitch = gains[2]
    ki_pitch = gains[3]
    kp_yaw = gains[4]
    ki_yaw = gains[5]
    try:
        file_path = '/home/bholder/gigabit_code/sim/src/gains.json'
        # Read the existing data from the file
        with open(file_path, 'r') as file:
            data = json.load(file)

        # Update the control parameters
        if "fine" in data:
            data["fine"]["kp_roll"] = float(kp_roll)
            data["fine"]["ki_roll"] = float(ki_roll)
            data["fine"]["kp_pitch"] = float(kp_pitch)
            data["fine"]["ki_pitch"] = float(ki_pitch)
            data["fine"]["kp_yaw"] = float(kp_yaw)
            data["fine"]["ki_yaw"] = float(ki_yaw)
        else:
            raise ValueError("JSON structure does not contain 'fine' key")

        # Write the updated data back to the file
        with open(file_path, 'w') as file:
            json.dump(data, file, indent=4)

    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
    except json.JSONDecodeError:
        print("Error: The file does not contain valid JSON.")
    except Exception as e:
        print(f"An error occurred: {e}")

# Example usage
# update_control_parameters('path_to_your_file.json', 400, 700, 300, 8500, 300, 11000)

def objective_function(params_lqr, params_pid):
    # Call MATLAB, Rust, and analyze_output functions here with 'params'
    # Calculate and return the value of 'sigma' to minimize
    # Example: (replace with actual calls)
    start = time.time()


    # Update control parameters for PID
    update_control_parameters(params_pid)

    # Call MATLAB function with LQR parameters
    call_matlab_function(params_lqr)
    
    call_rust_function()
    analysis_results = analyze_output()
    
    # Extract sigma values for 'ra', 'dec', and 'fr'
    sig_ra = analysis_results['ra']['sigma']
    sig_dec = analysis_results['dec']['sigma']
    sig_fr = analysis_results['fr']['sigma']
    
    p2p_ra = analysis_results['ra']['peak_to_peak']
    p2p_dec = analysis_results['dec']['peak_to_peak']
    p2p_fr = analysis_results['fr']['peak_to_peak']
    
    sig_ra = 206265 * sig_ra
    sig_dec = 206265 * sig_dec
    sig_fr = 206265 * sig_fr
    
    p2p_ra = 206265 * p2p_ra
    p2p_dec = 206265 * p2p_dec
    p2p_fr = 206265 * p2p_fr
    
    # You can perform calculations with sig_ra, sig_dec, and sig_fr here if needed
    # Example calculation:
    sigma = (sig_ra**2 + sig_dec**2 + sig_fr**2)
    
     # Save parameters and sigmas to a results file
    results_file = 'results.csv'

    if not os.path.exists(results_file):
        # Create a new results file with headers
        with open(results_file, 'w', newline='') as file:
            writer = csv.writer(file)
            # Write headers
            headers = ['P1', 'P2','P3','P4','P5','P6','P7','P8','P9','P10','P11', 'P12','kp_roll', 'ki_roll', 'kp_pitch', 'ki_pitch', 'kp_yaw', 'ki_yaw', 'Sig_RA', 'Sig_Dec', 'Sig_FR', 'P2P_RA','P2P_DEC','P2P_FR']
            writer.writerow(headers)
    
    # Append parameters and sigmas as a new row in the results file
    with open(results_file, 'a', newline='') as file:
        writer = csv.writer(file)
        # Concatenate the parameters and sigmas individually
        
        # params_list = params_lqr + params_pid
        params_list = np.concatenate([np.array(params_lqr), np.array(params_pid)]).tolist()

        row = params_list + [sig_ra, sig_dec, sig_fr, p2p_ra, p2p_dec, p2p_fr]
        writer.writerow(row)
    
    end = time.time()
    print("Time Elapsed: {}", end - start)
    return sigma

def append_to_csv(file_name, params, sigma):
    # Check if the file exists and create it with headers if it doesn't
    if not os.path.exists(file_name):
        with open(file_name, 'w', newline='') as file:
            writer = csv.writer(file)
            headers = ['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10', 'P11', 'P12', 'PID1', 'PID2', 'PID3', 'PID4', 'PID5', 'PID6', 'Sigma']
            writer.writerow(headers)

    # Append the new row
    with open(file_name, 'a', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(params + [sigma])


def optimize_pid(objective_function, initial_pid_params, pid_bounds, fixed_lqr_params, tolerance):
    pid_result = minimize(lambda pid_params: objective_function(fixed_lqr_params, pid_params), initial_pid_params, method='nelder-mead', bounds=pid_bounds, tol=tolerance)

    return pid_result.x, pid_result.fun

def optimize_lqr(objective_function, initial_lqr_params, lqr_bounds, fixed_pid_params, tolerance):

    lqr_result = minimize(lambda lqr_params: objective_function(lqr_params, fixed_pid_params), initial_lqr_params, method='nelder-mead', bounds=lqr_bounds, tol=tolerance)

    return lqr_result.x, lqr_result.fun

def optimize_all(objective_function, initial_params, param_bounds, tolerance):

    all_result = minimize(lambda params: objective_function(params[:12], params[12:]), initial_params, method='nelder-mead', bounds=param_bounds, tol=tolerance)

    return all_result.x, all_result.fun

####################################################################
#
#
#
#
#
#
#
#
#
#
####################################################################
def main():
    initial_params = [450.0, 800.0, 1000.0, 1500.0, 7500.0, 1000.0, 0.005, 0.4, 0.1, 0.01, 0.001, 0.1,
                        350.0, 600.0, 250.0, 8050.0, 260.0, 10000.0]

    param_bounds = [(0.01, 10000.0)] * 6 + [(0.00001, 10)] * 6 + [(0.1, 100000.0)] * 6

    tolerance = 1e-4

    [x,y] = optimize_pid(objective_function, initial_params[12:], param_bounds[12:], initial_params[:12], tolerance)
    # [x,y] = optimize_lqr(objective_function, initial_params[:12], param_bounds[:12], initial_params[12:], tolerance)
    # [x,y] = optimize_all(objective_function, initial_params, param_bounds, tolerance)

if __name__ == "__main__":
    main()
