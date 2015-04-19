#! /usr/bin/env python
import os
import subprocess
from copy import deepcopy

simulator = './simulator.py {infile_base}.fa -e {error} -c {cov} -f {infile_base_ef}.fa'
khmer_hash = './khmer/scripts/load-into-counting.py -x 1e8'\
             ' -k {k} hash_table.kh {infile_base}.fa'
khmer_hist = './khmer/scripts/abundance-dist.py'\
             ' hash_table.kh {infile_base}.fa {infile_base}_k{k}.dist'
khmer_cov = './khmer-recipes/005-estimate-total-genome-size/estimate-total-genome-size.py'\
            ' {infile_base}.fa {infile_base}_k{k}.dist {khmer_cov}'
# estimator = './coverage_estimator2.py {infile_base}_k{k}.dist -e {error} -k {k}'
estimator = './coverage_estimator3.py {infile_base}_k{k}.dist -e {error} -k {k} -c {cov}'

path = 'experiment2'

error_rates = [0.01, 0.03, 0.05, 0.1]
# error_rates = [0.05, 0.1]
coverages = [0.1, 0.5, 1, 2, 4, 10, 30, 50]
# coverages = [0.1, .5]
ks = [20, 30]

generate = False
compute_hist = generate
run_khmer = True


def run(command, output=None):
    f = None
    if output:
        f = open(output, 'w')
    return subprocess.call(command.split(), stdout=f)

if __name__ == '__main__':
    for c in coverages:
        for e in error_rates:
            params = {
                'error': e,
                'cov': c,
                'khmer_cov': max(1, int(c)),
            }
            infile_base = os.path.join(
                path, 'experiment1_c{cov}_e{error}'.format(**params)
            )
            infile_base_ef = os.path.join(
                path, 'experiment1_c{cov}_f{error}'.format(**params)
            )
            params['infile_base'] = infile_base
            params['infile_base_ef'] = infile_base_ef
            if generate:
                run(simulator.format(**params))
            for k in ks:
                params['k'] = k
                params2 = deepcopy(params)
                params2['infile_base'] = params2['infile_base_ef']
                for p in [params, params2]:
                    if compute_hist:
                        run(khmer_hash.format(**p))
                        run(khmer_hist.format(**p))
                    run(estimator.format(**p),
                        '{infile_base}_k{k}.est.out'.format(**p))
                    if run_khmer:
                        run(khmer_cov.format(**p),
                            '{infile_base}_k{k}.khmer.out'.format(**p))
