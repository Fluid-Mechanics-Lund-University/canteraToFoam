import os
import re
import argparse
import csv
import pandas as pd

def extract_uniform_value(file_path):
    """
    Reads a file and extracts the numeric value from the line matching:
    internalField   uniform   <value>;
    If the file doesn’t exist or no match is found, returns '0'.
    """
    pattern = re.compile(r'internalField\s+uniform\s+([0-9.eE+-]+);')
    if not os.path.isfile(file_path):
        return '0'
    with open(file_path, 'r') as f:
        for line in f:
            match = pattern.search(line)
            if match:
                return match.group(1)
    return '0'

def main():
    parser = argparse.ArgumentParser(
        description="Extract 'internalField uniform' values from files in numbered folders."
    )
    parser.add_argument(
        'species',
        help='Species names in parentheses, e.g. "(NH3 O2)"'
    )
    args = parser.parse_args()

    # Parse species list from the single argument
    species_list = args.species.strip('()').split()

    # Discover numeric directories in the current working directory
    timesteps = sorted(
        [d for d in os.listdir('.') if os.path.isdir(d) and re.match(r'^[0-9]+(\.[0-9]+)?$', d)],
        key=float
    )

    # Import pandas for DataFrame handling
    
    # Create a dictionary to store data
    data = {'timestep': []}
    for sp in species_list:
        data[sp] = []
    
    # Collect data for each timestep
    for ts in timesteps:
        data['timestep'].append(ts)
        for sp in species_list:
            file_path = os.path.join(ts, sp)
            value = extract_uniform_value(file_path)
            data[sp].append(value)
    
    # Create DataFrame
    df = pd.DataFrame(data)
    
    # Write to CSV
    output_file = "./output.csv"
    df.to_csv(output_file, index=False)
    

    print(f"Extraction complete. Results saved to {output_file}")

if __name__ == '__main__':
    main()

