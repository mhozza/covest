#! /usr/bin/env python
import os
import subprocess
from copy import deepcopy

source = 'data/chr14.fa'
source = 'data/ecoli2.fa'
source = 'data/yeast/chrXV.fa'
source = 'data/prame-partial-hg38.fa'

simulator = 'art_bin_VanillaIceCream/art_illumina -i {src} -o {infile_base} -f {cov} -l 100 -ef'
simulator_ef = './sam_to_fasta.py {infile_base}_errFree.sam'
jellyfish_count = 'jellyfish count -m {k} -s 500M -t 16 -C {infile} -o table.jf'
jellyfish_hist = 'jellyfish histo table.jf -o {infile_base}_k{k}.dist'
khmer_count = './khmer/scripts/load-into-counting.py -x 2e9 -n 6'\
              ' -k {k} hash_table.kh {infile}'
khmer_hist = './khmer/scripts/abundance-dist.py'\
             ' hash_table.kh {infile} {infile_base}_k{k}.dist'
khmer_cov = './khmer-recipes/005-estimate-total-genome-size/estimate-total-genome-size.py'\
            ' {infile} {infile_base}_k{k}.dist {khmer_cov}'
estimator = './covest.py {infile_base}_k{k}.dist -g -s {infile}'\
            ' -t 100 -T 16 -k {k} -c {cov}'
estimator_r = './covest.py {infile_base}_k{k}.dist -g -s {infile}'\
              ' -t 100 -T 16 -k {k} -c {cov} -rp'

path = 'experiment3p'

coverages = [0.5, 1, 4, 10, 50]
# coverages = [4, 10]
ks = [21]

generate = True
generate_dist = generate
use_jellyfish = True
ef = False
run_khmer = False
run_rep = False
run_norep = False


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
            if ef:
                run(simulator_ef.format(**params))
            # run('rm {path}/*.sam {path}/*.aln'.format(path=path))

        for k in ks:
            params['k'] = k
            pp = [params]
            if ef:
                params2 = deepcopy(params)
                params2['infile'] = params2['infile_ef']
                params2['infile_base'] = params2['infile_base_ef']
                pp.append(params2)
            for p in pp:
                if generate_dist:
                    if use_jellyfish:
                        run(jellyfish_count.format(**p))
                        run(jellyfish_hist.format(**p))
                    else:
                        run(khmer_count.format(**p))
                        run(khmer_hist.format(**p))
                        if run_khmer:
                            run(khmer_cov.format(**p),
                                '{infile_base}_k{k}.khmer.out'.format(**p))
                if run_norep:
                    run(estimator.format(**p),
                        '{infile_base}_k{k}.est.out'.format(**p))
                if run_rep:
                    run(estimator_r.format(**p),
                        '{infile_base}_k{k}.est_r.out'.format(**p))
