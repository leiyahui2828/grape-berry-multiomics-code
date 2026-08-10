import cooler
import pandas as pd
import numpy as np
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed

def calculate_valid_bins_percentage(file_path, threshold=1000, num_threads=1, chunk_size=100000, resolutions=None):
    if file_path.endswith('.mcool'):
        return calculate_mcool_valid_bins_percentage(file_path, threshold, num_threads, chunk_size, resolutions)
    else:
        return calculate_cool_valid_bins_percentage(file_path, threshold, chunk_size)

def calculate_cool_valid_bins_percentage(file_path, threshold=1000, chunk_size=100000):
    clr = cooler.Cooler(file_path)
    bins = clr.bins()[:]
    total_valid_bins = 0
    total_bins = 0

    for start in range(0, len(bins), chunk_size):
        end = min(start + chunk_size, len(bins))
        matrix = clr.matrix(balance=False)[start:end, start:end]

        bin_sums = np.sum(matrix, axis=0) + np.sum(matrix, axis=1) - np.diag(matrix)
        valid_bins = np.sum(bin_sums > threshold)
        total_valid_bins += valid_bins
        total_bins += bin_sums.size

    valid_bin_percentage = (total_valid_bins / total_bins) * 100 if total_bins > 0 else 0

    return {
        'File': file_path,
        'Valid Bin Percentage': valid_bin_percentage,
        'Total Valid Bins': total_valid_bins,
        'Total Bins': total_bins
    }

def calculate_mcool_valid_bins_percentage(file_path, threshold=1000, num_threads=1, chunk_size=100000, resolutions=None):
    available_resolutions = cooler.fileops.list_coolers(file_path)

    print(f"Available resolutions in {file_path}: {available_resolutions}")
    if resolutions:
        resolutions = [f'/resolutions/{res}' for res in resolutions.split(',') if f'/resolutions/{res}' in available_resolutions]
        print(resolutions)
        if not resolutions:
            raise ValueError(f"No matching resolutions found in the .mcool file for specified resolutions.")
    else:
        resolutions = available_resolutions

    results = []

    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        future_to_resolution = {
            executor.submit(process_resolution, file_path, resolution_path, threshold, chunk_size): resolution_path
            for resolution_path in resolutions
        }

        for future in as_completed(future_to_resolution):
            resolution_path = future_to_resolution[future]
            try:
                result = future.result()
                if result is not None:
                    results.append(result)
            except Exception as e:
                print(f"Error processing resolution {resolution_path}: {e}")

    return {
        'File': file_path,
        'Results': results
    }

def process_resolution(file_path, resolution_path, threshold, chunk_size):
    clr_res = cooler.Cooler(f"{file_path}::{resolution_path}")
    bins = clr_res.bins()[:]

    print(f"Processing resolution: {resolution_path}")

    total_valid_bins = 0
    total_bins = 0

    for start in range(0, len(bins), chunk_size):
        end = min(start + chunk_size, len(bins))
        matrix = clr_res.matrix(balance=False)[start:end, start:end]

        bin_sums = np.sum(matrix, axis=0) + np.sum(matrix, axis=1) - np.diag(matrix)
        valid_bins = np.sum(bin_sums > threshold)
        total_valid_bins += valid_bins
        total_bins += bin_sums.size

    valid_bin_percentage = (total_valid_bins / total_bins) * 100 if total_bins > 0 else 0

    return {
        'Resolution': resolution_path.split('/')[-1],
        'Valid Bin Percentage': valid_bin_percentage,
        'Total Valid Bins': total_valid_bins,
        'Total Bins': total_bins
    }

def save_results_to_csv(results, output_path):
    if 'Results' in results:
        df = pd.DataFrame(results['Results'])
        df.insert(0, 'File', results['File'])
    else:
        df = pd.DataFrame([results])

    df.to_csv(output_path, index=False)
    print(f"Results saved to {output_path}")

def main():
    parser = argparse.ArgumentParser(description="Calculate valid bins percentage from .cool or .mcool files.")
    parser.add_argument("file_path", type=str, help="Path to the .cool or .mcool file")
    parser.add_argument("--output", type=str, required=True, help="Path to save the results as CSV file")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads to use for processing")
    parser.add_argument("--chunk_size", type=int, default=100000, help="Chunk size for processing large matrices")
    parser.add_argument("--resolutions", type=str, help="Comma-separated list of resolutions to process in the .mcool file")

    args = parser.parse_args()

    try:
        results = calculate_valid_bins_percentage(args.file_path, num_threads=args.threads, chunk_size=args.chunk_size, resolutions=args.resolutions)
        save_results_to_csv(results, args.output)
    except Exception as e:
        print(f"Error processing file {args.file_path}: {e}")

if __name__ == "__main__":
    main()

