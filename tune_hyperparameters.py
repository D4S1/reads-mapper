import argparse
import itertools
import mapper
import utils
import pandas as pd
import traceback
from datetime import datetime

def grid_search(param_grid: dict, reads_filename: str, genome_filename: str, output_prefix: str) -> pd.DataFrame:
    """
    Perform a grid search over the parameter grid and execute the mapper.
    
    Args:
        param_grid (dict): A dictionary where keys are parameter names and values are lists of possible values.
        reads_filename (str): Path to the reads file.
        genome_filename (str): Path to the genome file.
        output_prefix (str): Prefix for the output files.
    
    Returns:
        pd.DataFrame: A DataFrame containing parameter combinations and their corresponding accuracy scores.
    """
    log_file = f"{output_prefix}_log.txt"
    results = []
    keys, values = zip(*param_grid.items())

    with open(log_file, 'w') as log:
        log.write(f"Grid Search Log - {datetime.now()}\n")
        log.write("="*50 + "\n")

        for i, combination in enumerate(itertools.product(*values)):
            params = dict(zip(keys, combination))
            out_file = f"{output_prefix}_run_{i}.txt"  # Unique output file for each parameter combination

            try:
                log.write(f"\nRunning with parameters: {params}\n")

                # Call the mapper function with unpacked parameters
                process_time, max_time, avg_time = mapper.main(
                    reads_filename,
                    genome_filename,
                    params['wind_size'],
                    params['k'],
                    params['hash'],
                    params['err_max'],
                    params['delta'],
                    out_file
                )
                if process_time is None:
                    continue
                # Calculate accuracy
                accuracy = utils.accuracy(out_file, reads_filename.replace('fasta', 'txt'))

                # Log the results
                log.write(f"Process time: {process_time} min\n")
                log.write(f"Max time: {max_time} min\n")
                log.write(f"Average time per read: {avg_time} min\n")
                log.write(f"Accuracy: {accuracy:.4f}\n")

                # Append results to DataFrame
                print({**params, 'accuracy': accuracy, 'process_time':process_time, 'avg_time': avg_time})
                results.append({**params, 'accuracy': accuracy, 'process_time':process_time, 'avg_time': avg_time})

            except AssertionError as e:
                if "m cannot be 0" in str(e):
                    log.write(f"Error: {e}\n")
                    log.write(f"Skipping parameters: {params}\n")
                else:
                    log.write(f"Unexpected AssertionError: {e}\n")
                    log.write(traceback.format_exc())
            except Exception as e:
                log.write(f"Unexpected error: {e}\n")
                log.write(traceback.format_exc())

    # Convert results to a DataFrame
    results_df = pd.DataFrame(results)
    return results_df

if __name__ == '__main__':
    param_grid = {
        'wind_size': [80, 90, 100],
        'k': [9, 10, 11],
        'hash': [12345],
        'err_max': [0.075, 0.1, 0.125],
        'delta': [0.1, 0.125, 0.15, 0.175]
    }

    parser = argparse.ArgumentParser()
    parser.add_argument('reference', type=str, help="Path to the genome reference file")
    parser.add_argument('reads', type=str, help="Path to the reads file")
    parser.add_argument('output', type=str, help="Prefix for output files")
    args = parser.parse_args()

    # Perform the grid search
    results_df = grid_search(param_grid, args.reads, args.reference, args.output)

    # Save results to CSV
    results_csv = f"{args.output}_results.csv"
    results_df.to_csv(results_csv, index=False)
    print(f"Grid search results saved to {results_csv}")
