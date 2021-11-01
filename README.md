# ratgtex-server-data
Process RatGTEx pipeline results into data files for the RatGTEx site

## Requirements

- pandas (Python)
- gtfparse (Python)
- tidyverse (R)

## Usage

Run `run.sh`, which calls R and Python scripts that prepare the server data. The parameters are:

1. Path to the RatGTEx pipeline base directory
2. Output directory path
3. The remaining parameters are the list of tissues to include
