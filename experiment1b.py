#! /usr/bin/env python
import os
import subprocess
from copy import deepcopy

source = 'data/chr14.fa'
# source = 'data/yeast/chrXV.fa'

simulator = 'art_bin_VanillaIceCream/art_illumina -i {src} -o {infile_base} -f {cov} -l 100 -ef'
simulator_ef = './sam_to_fasta.py {infile_base}_errFree.sam'
khmer_hash = './khmer/scripts/load-into-counting.py -x 1e9'\
             ' -k {k} hash_table.kh {infile}'
khmer_hist = './khmer/scripts/abundance-dist.py'\
             ' hash_table.kh {infile} {infile_base}_k{k}.dist'
khmer_cov = './khmer-recipes/005-estimate-total-genome-size/estimate-total-genome-size.py'\
            ' {infile} {infile_base}_k{k}.dist {khmer_cov}'
estimator = './coverage_estimator3.py {infile_base}_k{k}.dist -k {k} --repeats'

path = 'experiment2h'

coverages = [0.1, 0.5, 1, 2, 4, 10, 50]
ks = [20, 30]

generate = False


def run(command, output=None):
    f = None
    if output:
        f = open(output, 'w')
    return subprocess.call(command.split(), stdout=f)

if __name__ == '__main__':
    for c in coverages:
        params = {
            'src': source,
            'cov': c,
            'khmer_cov': max(1, int(c)),
        }
        infile_base = os.path.join(
            path, 'experiment1_c{cov}'.format(**params)
        )

        params['infile_base'] = infile_base
        params['infile_base_ef'] = infile_base + 'f'
        params['infile'] = '{}.fq'.format(params['infile_base'])
        params['infile_ef'] = '{}_errFree.fa'.format(params['infile_base'])

        if generate:
            run(simulator.format(**params))
            run(simulator_ef.format(**params))
            run('rm {path}/*.sam {path}/*.aln'.format(path=path))

        for k in ks:
            params['k'] = k
            params2 = deepcopy(params)
            params2['infile'] = params2['infile_ef']
            params2['infile_base'] = params2['infile_base_ef']
            for p in [params, params2]:
                run(khmer_hash.format(**p))
                run(khmer_hist.format(**p))
                run(khmer_cov.format(**p),
                    '{infile_base}_k{k}.khmer.out'.format(**p))
                run(estimator.format(**p),
                    '{infile_base}_k{k}.est.out'.format(**p))
