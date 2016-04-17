# covest
Tool that estimates coverage (and genome size) of dna sequence from reads.

```
usage: covest.py [-h] [--load LOAD] [-k KMER_SIZE] [-r READ_LENGTH]
                 [--plot [PLOT]] [-rp] [-ll] [-t TRIM] [-M MAX_COVERAGE]
                 [-at [AUTO_TRIM]] [-g GRID] [-e ERROR_RATE] [-es ERROR_SCALE]
                 [-c COVERAGE] [-q1 Q1] [-mq1 MIN_Q1] [-q2 Q2] [-q Q] [-so]
                 [-f] [-T THREAD_COUNT] [-s GENOME_SIZE]
                 input_histogram

Simulate reads form random genome with errors

positional arguments:
  input_histogram       Input histogram

optional arguments:
  -h, --help            show this help message and exit
  --load LOAD           Load json
  -k KMER_SIZE, --kmer-size KMER_SIZE
                        Kmer size
  -r READ_LENGTH, --read-length READ_LENGTH
                        Read length
  --plot [PLOT]         Plot probabilities (use --plot 1 to plot probs * j)
  -rp, --repeats        Estimate with repeats
  -ll, --ll-only        Only compute log likelihood
  -t TRIM, --trim TRIM  Trim histogram at this value
  -M MAX_COVERAGE, --max-coverage MAX_COVERAGE
                        Upper coverage limit
  -at [AUTO_TRIM], --auto-trim [AUTO_TRIM]
                        Trim histogram automatically with this threshold
  -g GRID, --grid GRID  Grid search type: 0 - None, 1 - Pre-grid, 2 - Post-
                        grid
  -e ERROR_RATE, --error-rate ERROR_RATE
                        Error rate
  -es ERROR_SCALE, --error-scale ERROR_SCALE
                        Error scale
  -c COVERAGE, --coverage COVERAGE
                        Coverage
  -q1 Q1                q1
  -mq1 MIN_Q1, --min-q1 MIN_Q1
                        minimum single copy ratio
  -q2 Q2                q2
  -q Q                  q
  -so, --start-original
                        Start form given values
  -f, --fix             Fix some vars, optimize others
  -T THREAD_COUNT, --thread-count THREAD_COUNT
                        Thread count
  -s GENOME_SIZE, --genome-size GENOME_SIZE
                        Calculate genome size from reads
```
