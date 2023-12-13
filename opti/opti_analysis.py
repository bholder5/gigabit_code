import pandas as pd


def analyze_data(file_path):
    # Read the CSV file
    data = pd.read_csv(file_path)

    # Find the maximum time and filter data for the last 300 seconds
    max_time = data['t'].max()
    filtered_data = data[data['t'] > (max_time - 300)]

    # Columns to analyze
    columns = ['ra', 'dec', 'fr']

    analysis_results = {}
    for col in columns:
        # Peak-to-peak delta
        peak_to_peak = filtered_data[col].max() - filtered_data[col].min()

        # 1-sigma variation and subtract the mean
        sigma = filtered_data[col].std()
        mean_adjusted = filtered_data[col] - filtered_data[col].mean()

        analysis_results[col] = {
            'peak_to_peak': peak_to_peak,
            'sigma': sigma
        }

    return analysis_results

# # Call the function with your CSV file path
# results = analyze_data('/home/bholder/data/out.csv')

