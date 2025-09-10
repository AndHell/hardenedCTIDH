import time
import sys
import json
import csv
import os
from datetime import datetime
from pathlib import Path

def save_optimization_results(cost, N0, m0, dac_paths, initial_batch_sizes, B,
                             primes, minkeyspace, cpu_time, wall_time, n_ells, cofactors_conf, output_dir="results"):
    """
    Save CSIDH optimization results to JSON and CSV files organized by number of primes (ells) and batch size.

    Args:
        cost: Optimization cost result
        N0: Final batch sizes tuple
        m0: Final batch parameters tuple
        dac_paths: Dictionary from build_dac_paths function
        initial_batch_sizes: Tuple from find_initial_batch_sizes function
        B: Number of batches
        primes: List of prime numbers used
        minkeyspace: Minimum keyspace requirement
        cpu_time: CPU time consumed in seconds
        wall_time: Wall clock time consumed in seconds
        n_ells: Number of primes (ells) used
        output_dir: Directory to save results (default: "results")
    """
    # Create output directory if it doesn't exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Create cofactors_conf-specific subdirectory
    cofactors_dir = Path(output_dir) / f"{cofactors_conf}"
    cofactors_dir.mkdir(parents=True, exist_ok=True)

    # Create ells-specific subdirectory
    ells_dir = Path(cofactors_dir) / f"ells_{n_ells}"
    ells_dir.mkdir(parents=True, exist_ok=True)

    # Create batch-specific subdirectory within ells directory
    batch_dir = ells_dir / f"batch_{B}"
    batch_dir.mkdir(parents=True, exist_ok=True)

    # Generate timestamp for unique filenames
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    # Prepare data structure for JSON
    results_data = {
        "metadata": {
            "timestamp": timestamp,
            "n_ells": n_ells,
            "batch_size": B,
            "num_primes": len(primes),
            "minkeyspace": minkeyspace,
            "optimization_cost": cost,
            "cpu_time_seconds": cpu_time,
            "wall_time_seconds": wall_time,
            "cpu_time_formatted": format_time(cpu_time),
            "wall_time_formatted": format_time(wall_time)
        },
        "configuration": {
            "initial_batch_sizes": list(initial_batch_sizes) if initial_batch_sizes else None,
            "final_batch_sizes": list(N0),
            "batch_parameters": list(m0)
        },
        "primes": primes,
        "dac_paths": {
            str(batch_idx): {
                "lengths": list(info["lengths"]),
                "primes": info["primes"]
            }
            for batch_idx, info in dac_paths.items()
        }
    }

    # Save JSON file
    json_filename = batch_dir / f"optimization_results_{timestamp}.json"
    with open(json_filename, 'w') as f:
        json.dump(results_data, f, indent=2)

    # Save CSV file with batch details
    csv_filename = batch_dir / f"batch_details_{timestamp}.csv"
    with open(csv_filename, 'w', newline='') as f:
        writer = csv.writer(f)

        # Write header
        writer.writerow([
            "Batch_Index", "Initial_Size", "Final_Size", "Batch_Parameter",
            "Num_Primes_in_Batch", "DAC_Path_Lengths", "Primes_in_Batch"
        ])

        # Write batch data
        for batch_idx in range(B):
            initial_size = initial_batch_sizes[batch_idx] if initial_batch_sizes else "N/A"
            final_size = N0[batch_idx] if batch_idx < len(N0) else "N/A"
            batch_param = m0[batch_idx] if batch_idx < len(m0) else "N/A"

            if batch_idx in dac_paths:
                dac_info = dac_paths[batch_idx]
                num_primes = len(dac_info["primes"])
                path_lengths = ";".join(map(str, sorted(dac_info["lengths"])))
                batch_primes = ";".join(map(str, dac_info["primes"]))
            else:
                num_primes = 0
                path_lengths = ""
                batch_primes = ""

            writer.writerow([
                batch_idx, initial_size, final_size, batch_param,
                num_primes, path_lengths, batch_primes
            ])

    # Save summary CSV file for this ells value (append mode for tracking multiple runs)
    summary_csv = ells_dir / f"ells_{n_ells}_summary.csv"
    file_exists = summary_csv.exists()

    with open(summary_csv, 'a', newline='') as f:
        writer = csv.writer(f)

        # Write header only if file doesn't exist
        if not file_exists:
            writer.writerow([
                "Timestamp", "Batch_Size", "Optimization_Cost", "Minkeyspace", "Num_Primes",
                "Initial_Batch_Sizes", "Final_Batch_Sizes", "Batch_Parameters",
                "Total_DAC_Paths", "Min_Path_Length", "Max_Path_Length",
                "CPU_Time_Seconds", "Wall_Time_Seconds", "CPU_Time_Formatted", "Wall_Time_Formatted"
            ])

        # Calculate path statistics
        all_path_lengths = []
        total_dac_paths = 0
        for info in dac_paths.values():
            all_path_lengths.extend(info["lengths"])
            total_dac_paths += len(info["lengths"])

        min_path_length = min(all_path_lengths) if all_path_lengths else 0
        max_path_length = max(all_path_lengths) if all_path_lengths else 0

        writer.writerow([
            timestamp, B, cost, minkeyspace, len(primes),
            ";".join(map(str, initial_batch_sizes)) if initial_batch_sizes else "N/A",
            ";".join(map(str, N0)),
            ";".join(map(str, m0)),
            total_dac_paths, min_path_length, max_path_length,
            cpu_time, wall_time, format_time(cpu_time), format_time(wall_time)
        ])

    # Save batch-specific summary CSV file (append mode for tracking multiple runs)
    batch_summary_csv = batch_dir / f"batch_{B}_summary.csv"
    batch_file_exists = batch_summary_csv.exists()

    with open(batch_summary_csv, 'a', newline='') as f:
        writer = csv.writer(f)

        # Write header only if file doesn't exist
        if not batch_file_exists:
            writer.writerow([
                "Timestamp", "Optimization_Cost", "Minkeyspace", "Num_Primes",
                "Initial_Batch_Sizes", "Final_Batch_Sizes", "Batch_Parameters",
                "Total_DAC_Paths", "Min_Path_Length", "Max_Path_Length",
                "CPU_Time_Seconds", "Wall_Time_Seconds", "CPU_Time_Formatted", "Wall_Time_Formatted"
            ])

        writer.writerow([
            timestamp, cost, minkeyspace, len(primes),
            ";".join(map(str, initial_batch_sizes)) if initial_batch_sizes else "N/A",
            ";".join(map(str, N0)),
            ";".join(map(str, m0)),
            total_dac_paths, min_path_length, max_path_length,
            cpu_time, wall_time, format_time(cpu_time), format_time(wall_time)
        ])

    # Save a consolidated report
    report_filename = batch_dir / f"optimization_report_{timestamp}.txt"
    with open(report_filename, 'w') as f:
        f.write(f"CSIDH Optimization Report - {n_ells} Primes, Batch Size {B}\n")
        f.write(f"=" * 60 + "\n\n")
        f.write(f"Timestamp: {timestamp}\n")
        f.write(f"Number of Primes (ells): {n_ells}\n")
        f.write(f"Optimization Cost: {cost:.6f}\n")
        f.write(f"Minimum Keyspace: {minkeyspace}\n")
        f.write(f"Number of Primes Used: {len(primes)}\n")
        f.write(f"Number of Batches: {B}\n")
        f.write(f"CPU Time: {format_time(cpu_time)} ({cpu_time:.6f} seconds)\n")
        f.write(f"Wall Time: {format_time(wall_time)} ({wall_time:.6f} seconds)\n")
        f.write(f"CPU Efficiency: {(cpu_time/wall_time)*100:.2f}%\n\n")

        f.write("Configuration:\n")
        f.write(f"  Initial Batch Sizes: {initial_batch_sizes}\n")
        f.write(f"  Final Batch Sizes: {N0}\n")
        f.write(f"  Batch Parameters: {m0}\n\n")

        f.write("Batch Details:\n")
        for batch_idx in range(B):
            f.write(f"  Batch {batch_idx}:\n")
            if batch_idx in dac_paths:
                info = dac_paths[batch_idx]
                f.write(f"    Primes: {info['primes']}\n")
                f.write(f"    DAC Path Lengths: {sorted(info['lengths'])}\n")
            f.write(f"    Final Size: {N0[batch_idx] if batch_idx < len(N0) else 'N/A'}\n")
            f.write(f"    Parameter: {m0[batch_idx] if batch_idx < len(m0) else 'N/A'}\n\n")

        f.write("Files Generated:\n")
        f.write(f"  JSON: {json_filename.name}\n")
        f.write(f"  CSV Details: {csv_filename.name}\n")
        f.write(f"  Batch Summary: {batch_summary_csv.name}\n")
        f.write(f"  Ells Summary: {summary_csv.name}\n")
        f.write(f"  Report: {report_filename.name}\n")

    print(f"Results saved to directory: {batch_dir}")
    print(f"Generated files:")
    print(f"  - JSON: {json_filename}")
    print(f"  - CSV Details: {csv_filename}")
    print(f"  - Batch Summary: {batch_summary_csv}")
    print(f"  - Ells Summary: {summary_csv}")
    print(f"  - Report: {report_filename}")
    print(f"Timing: CPU={format_time(cpu_time)}, Wall={format_time(wall_time)}")

    return {
        "json_file": str(json_filename),
        "csv_file": str(csv_filename),
        "batch_summary_file": str(batch_summary_csv),
        "ells_summary_file": str(summary_csv),
        "report_file": str(report_filename)
    }


def format_time(seconds):
    """
    Format time in seconds to a human-readable string.

    Args:
        seconds: Time in seconds

    Returns:
        str: Formatted time string
    """
    if seconds < 60:
        return f"{seconds:.2f}s"
    elif seconds < 3600:
        minutes = int(seconds // 60)
        secs = seconds % 60
        return f"{minutes}m {secs:.2f}s"
    else:
        hours = int(seconds // 3600)
        minutes = int((seconds % 3600) // 60)
        secs = seconds % 60
        return f"{hours}h {minutes}m {secs:.2f}s"


def get_cpu_time():
    """
    Get CPU time for current process and all child processes.

    Returns:
        float: CPU time in seconds
    """
    try:
        # Get current process
        process = psutil.Process(os.getpid())

        # Get CPU times for current process
        cpu_times = process.cpu_times()
        cpu_time = cpu_times.user + cpu_times.system

        # Add CPU time from child processes (including terminated ones)
        try:
            # This gets cumulative CPU time of all children (including terminated)
            children_times = process.cpu_times(percpu=False)
            # Try to get children CPU times more reliably
            for child in process.children(recursive=True):
                try:
                    child_cpu_times = child.cpu_times()
                    cpu_time += child_cpu_times.user + child_cpu_times.system
                except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
                    # Child process might have terminated
                    pass
        except Exception:
            pass

        return cpu_time
    except Exception:
        # Fallback to time.process_time() if psutil fails
        return time.process_time()


def get_process_cpu_time():
    """
    Alternative CPU time measurement using resource module.

    Returns:
        float: CPU time in seconds
    """
    try:
        import resource
        usage = resource.getrusage(resource.RUSAGE_SELF)
        children_usage = resource.getrusage(resource.RUSAGE_CHILDREN)

        # Sum user and system time for self and children
        total_cpu = (usage.ru_utime + usage.ru_stime +
                    children_usage.ru_utime + children_usage.ru_stime)
        return total_cpu
    except ImportError:
        # Fallback if resource module not available
        return time.process_time()


def load_optimization_results(json_file_path):
    """
    Load optimization results from a JSON file.

    Args:
        json_file_path: Path to the JSON file containing results

    Returns:
        dict: Loaded optimization results
    """
    with open(json_file_path, 'r') as f:
        return json.load(f)


def compare_batch_results(batch_sizes, n_ells, results_dir="results"):
    """
    Compare optimization results across different batch sizes for a specific number of ells.

    Args:
        batch_sizes: List of batch sizes to compare
        n_ells: Number of primes (ells) to compare
        results_dir: Directory containing results

    Returns:
        dict: Comparison data
    """
    comparison_data = {}
    ells_dir = Path(results_dir) / f"ells_{n_ells}"

    for B in batch_sizes:
        batch_dir = ells_dir / f"batch_{B}"
        summary_file = batch_dir / f"batch_{B}_summary.csv"

        if summary_file.exists():
            with open(summary_file, 'r') as f:
                reader = csv.DictReader(f)
                # Get the most recent result (last row)
                rows = list(reader)
                if rows:
                    latest_result = rows[-1]
                    comparison_data[B] = {
                        "cost": float(latest_result["Optimization_Cost"]),
                        "num_primes": int(latest_result["Num_Primes"]),
                        "min_path_length": int(latest_result["Min_Path_Length"]),
                        "max_path_length": int(latest_result["Max_Path_Length"]),
                        "cpu_time": float(latest_result.get("CPU_Time_Seconds", 0)),
                        "wall_time": float(latest_result.get("Wall_Time_Seconds", 0)),
                        "timestamp": latest_result["Timestamp"]
                    }

    return comparison_data


def compare_ells_results(n_ells_list, results_dir="results"):
    """
    Compare optimization results across different numbers of ells.

    Args:
        n_ells_list: List of ells numbers to compare
        results_dir: Directory containing results

    Returns:
        dict: Comparison data organized by ells
    """
    comparison_data = {}

    for n_ells in n_ells_list:
        ells_dir = Path(results_dir) / f"ells_{n_ells}"
        summary_file = ells_dir / f"ells_{n_ells}_summary.csv"

        if summary_file.exists():
            with open(summary_file, 'r') as f:
                reader = csv.DictReader(f)
                rows = list(reader)
                if rows:
                    # Find the best result (minimum cost) for this ells value
                    best_result = min(rows, key=lambda x: float(x["Optimization_Cost"]))
                    comparison_data[n_ells] = {
                        "cost": float(best_result["Optimization_Cost"]),
                        "batch_size": int(best_result["Batch_Size"]),
                        "num_primes": int(best_result["Num_Primes"]),
                        "min_path_length": int(best_result["Min_Path_Length"]),
                        "max_path_length": int(best_result["Max_Path_Length"]),
                        "cpu_time": float(best_result.get("CPU_Time_Seconds", 0)),
                        "wall_time": float(best_result.get("Wall_Time_Seconds", 0)),
                        "timestamp": best_result["Timestamp"]
                    }

    return comparison_data
