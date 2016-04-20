#! /usr/bin/env python
import os
import subprocess
import argparse

source = 'data/chr14.fa'
source = 'data/ecoli2.fa'
source = 'data/yeast/chrXV.fa'
source = 'data/prame-partial-hg38.fa'

simulator = 'art_bin_VanillaIceCream/art_illumina -i {src} -o {infile_base} -f {cov} -l 100'
simulator_simple = './read_simulator.py {src} {infile_base}.fa'\
                   ' -e 0.03 -c {cov}'

jellyfish_count = 'jellyfish count -m {k} -s 500M -t 16 -C {infile} -o {infile}.jf'
jellyfish_hist = 'jellyfish histo {infile}.jf_0 -o {infile_base}_k{k}.dist'
khmer_count = './khmer/scripts/load-into-counting.py -x 2e9 -n 6'\
              ' -k {k} hash_table.kh {infile}'
khmer_hist = './khmer/scripts/abundance-dist.py'\
             ' hash_table.kh {infile} {infile_base}_k{k}.dist'
khmer_cov = './khmer-recipes/005-estimate-total-genome-size/estimate-total-genome-size.py'\
            ' {infile} {infile_base}_k{k}.dist {khmer_cov}'
estimator = './covest.py {infile_base}_k{k}.dist -g -s {infile}'\
            ' -t 100 -T 16 -k {k} -c {cov} -M 700'
estimator_r = './covest.py {infile_base}_k{k}.dist -g -s {infile}'\
              ' -t 100 -T 16 -k {k} -c {cov} -rp -M 700'

coverages = [0.5, 1, 2, 4, 10, 30, 50]
# coverages = [4, 10]
ks = [21]

USE_SIMPLE_SIMULATOR = False
generate = False
generate_dist = True
use_jellyfish = True
run_khmer = False
run_rep = False
run_norep = False


def run(command, output=None):
    f = None
    if output:
        f = open(output, 'w')
    return subprocess.call(command.split(), stdout=f)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('path')
    args = parser.parse_args()
    path = args.path

    for c in coverages:
        params = {
            'src': source,
            'cov': c,
            'khmer_cov': max(1, int(c)),
        }
        infile_base = os.path.join(
            path, 'ecoli-c{cov}'.format(**params)
        )

        if USE_SIMPLE_SIMULATOR:
            params['infile_base'] = infile_base
            params['infile'] = '{}.fa'.format(params['infile_base'])
        else:
            params['infile_base'] = infile_base
            params['infile'] = '{}.fastq'.format(params['infile_base'])

        if generate:
            if USE_SIMPLE_SIMULATOR:
                run(simulator_simple.format(**params))
            else:
                run(simulator.format(**params))
        for k in ks:
            params['k'] = k
            pp = [params]
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
