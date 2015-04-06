#! /usr/bin/env python
import os

khmer_hash = './khmer/scripts/load-into-counting.py -x 1e8'\
             ' -k {k} hash_table.kh {infile_base}.fa'
khmer_hist = './khmer/scripts/abundance-dist.py'\
             ' hash_table.kh {infile_base}.fa {infile_base}.k{k}.dist'
khmer_cov = './khmer-recipes/005-estimate-total-genome-size/estimate-total-genome-size.py'\
            ' {infile_base}.fa {infile_base}.k{k}.dist 10'\
            '>{infile_base}.k{k}.khmer.out'
simulator = './simulator.py {infile_base}.fa -e {error} -c {cov} -f {infile_base_ef}.fa'
estimator = './coverage_estimator2.py {infile_base}.k{k}.dist -e {error} -k {k}'\
            '>{infile_base}.k{k}.out'

error_rates = [0.01, 0.03, 0.05, 0.1, 0.5]
coverages = [1, 2, 4, 10, 16]
ks = [15, 20, 25, 30]

for c in coverages:
    for e in error_rates:
        path = 'experiment1'
        infile_base = os.path.join(path, 'experiment1.c{c}.e{e}')
        infile_base_ef = os.path.join(path, 'experiment1.c{c}.ef')
        params = {
            'error': e,
            'cov': c,
            'infile_base': infile_base,
            'infile_base_ef': infile_base_ef,
        }
        os.system(simulator.format(**params))
        for k in ks:
            params2 = params.copy()
            params2['infile_base'] = params2['infile_base_ef']
            for params in [params, params2]:
                params['k'] = k
                os.system(khmer_hash.format(**params))
                os.system(khmer_hist.format(**params))
                os.system(khmer_cov.format(**params))
                os.system(estimator.format(**params))
