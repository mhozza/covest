# covest
Tool that estimates coverage (and genome size) of dna sequence from reads.

## Installation

### For development:
`pip install -e .` from the project directory

## Usage
type `covest --help` for the usage.

## Other included tools

- `geset.py` tool for estimation genome size from reads and known coverage
- `prepare_experient.py` tool for experiment pipeline setup
- `experiment_table.py` tool which collects data from experiment and create a nice table (html, tex, and csv formats are supported)
- `sam_to_fasta.py` tool for converting sam file to fasta
- Read simulator:
    - `generate_sequence.py` random sequence generator
    - `read_simulator.py` tool for generating random reads form the sequence


