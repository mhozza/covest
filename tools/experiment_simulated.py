#! /usr/bin/env python3
import os
import subprocess
from copy import deepcopy

generator = './generate_sequence.py {seq_name}.fa'
simulator = './read_simulator.py {seq_name}.fa {infile_base}.fa'\
            ' -f {infile_base_ef}.fa -e {error} -c {cov}'
jellyfish_count = 'jellyfish count -m {k} -s 500M -t 16 -C {infile_base}.fa -o {infile_base}.jf'
jellyfish_hist = 'jellyfish histo {infile_base}.jf_0 -o {infile_base}_k{k}.dist'
khmer_count = './khmer/scripts/load-into-counting.py -x 1e9'\
              ' -k {k} hash_table.kh {infile_base}.fa'
khmer_hist = './khmer/scripts/abundance-dist.py'\
             ' hash_table.kh {infile_base}.fa {infile_base}_k{k}.dist'
khmer_cov = './khmer-recipes/005-estimate-total-genome-size/estimate-total-genome-size.py'\
            ' {infile_base}.fa {infile_base}_k{k}.dist {khmer_cov}'
# estimator = './coverage_estimator2.py {infile_base}_k{k}.dist -e {error} -k {k}'
estimator = './covest.py {infile_base}_k{k}.dist -g -s {infile_base}.fa'\
            ' -t 100 -T 16 -e {error} -k {k} -c {cov} -rp'

path = 'experiment3_3'

seq_cnt = 5
error_rates = [0.01, 0.03, 0.05, 0.1]
#error_rates = [0.03]
coverages = [0.5, 1, 2, 4, 10, 50]
#coverages = [0.1]
ks = [21]

generate = False
compute_hist = False
run_khmer = False
use_jellyfish = True
VERBOSE = True


def run(command, output=None):
    if VERBOSE:
        print('executing command:', command)
    f = None
    if output:
        f = open(output, 'w')
    return subprocess.call(command.split(), stdout=f)


if __name__ == '__main__':
    for s in range(seq_cnt):
        seq_name = os.path.join(path, 'simulated{}'.format(s))
        if generate:
            run(generator.format(seq_name=seq_name))
        for c in coverages:
            for e in error_rates:
                params = {
                    'seq_name': seq_name,
                    'error': e,
                    'cov': c,
                    'khmer_cov': max(1, int(c)),
                }
                infile_base = '{seq_name}_c{cov}_e{error}'.format(**params)
                infile_base_ef = '{seq_name}_c{cov}_f{error}'.format(**params)
                params['infile_base'] = infile_base
                params['infile_base_ef'] = infile_base_ef
                if generate:
                    run(simulator.format(**params))
                    run(simulator.format(**params))
                for k in ks:
                    params['k'] = k
                    params2 = deepcopy(params)
                    params2['infile_base'] = params2['infile_base_ef']
                    for p in [params, params2]:
                        if compute_hist:
                            if use_jellyfish:
                                run(jellyfish_count.format(**p))
                                run(jellyfish_hist.format(**p))
                            else:
                                run(khmer_count.format(**p))
                                run(khmer_hist.format(**p))
                        run(estimator.format(**p),
                            '{infile_base}_k{k}.est.out'.format(**p))
                        if run_khmer:
                            run(khmer_cov.format(**p),
                                '{infile_base}_k{k}.khmer.out'.format(**p))
