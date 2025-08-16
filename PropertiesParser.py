import csv

def parse_records(input_file, output_file):
    # Read the entire file
    with open(input_file, 'r') as f:
        content = f.read()
    
    # Split records by blank lines
    records = [block.strip() for block in content.strip().split('\n\n') if block.strip()]
    
    # Parse each record into a dictionary
    data = []
    for record in records:
        record_dict = {}
        for line in record.split('\n'):
            if ':' in line:
                key, value = line.split(':', 1)
                record_dict[key.strip()] = value.strip()
        data.append(record_dict)
    
    # Get all unique column names (in case some records have missing fields)
    columns = sorted(set(key for d in data for key in d.keys()))
    
    # Write to CSV
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=columns)
        writer.writeheader()
        for row in data:
            writer.writerow(row)

# Example usage:
parse_records('new_properties.txt', 'cleaned_data.csv')