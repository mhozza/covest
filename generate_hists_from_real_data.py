#! /usr/bin/env python
import argparse
import os
import subprocess
from copy import deepcopy

# source = 'data/chr14.fa'
# source = 'data/yeast/chrXV.fa'
# source = 'data/prame-partial-hg38.fa'

simulator = 'art_bin_VanillaIceCream/art_illumina -i {src} -o {infile_base} -f {cov} -l 100 -ef'
simulator_ef = './sam_to_fasta.py {infile_base}_errFree.sam'
jellyfish_count = 'jellyfish count -m {k} -s 500M -t 16 -C {infile} -o table.jf'
jellyfish_hist = 'jellyfish histo table.jf -o {infile_base}_k{k}.dist'
khmer_count = './khmer/scripts/load-into-counting.py -x 2e9 -n 6'\
              ' -k {k} hash_table.kh {infile}'
khmer_hist = './khmer/scripts/abundance-dist.py'\
             ' hash_table.kh {infile} {infile_base}_k{k}.dist'


path = 'experiment2p'

coverages = [0.5, 1, 2, 4, 10, 30]
ks = [21]

generate = True
generate_dist = generate
use_jellyfish = True
ef = True


def run(command, output=None):
    f = None
    if output:
        f = open(output, 'w')
    return subprocess.call(command.split(), stdout=f)


def main(args):
    for c in coverages:
        path = args.path
        params = {
            'src': args.source,
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


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Simulate reads form random genome with errors')
    parser.add_argument('source', help='source sequence')
    parser.add_argument('path', help='target path')
    args = parser.parse_args()
    main(args)
