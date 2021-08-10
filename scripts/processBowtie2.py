import argparse
import sys
import os
import glob
import logging
import pandas as pd
import pysam

class App(object):

    def run(self, args):
        name = os.path.basename(args[0])
        parser = self.create_parser(name)
        opts = parser.parse_args(args[1:])
        return self.main(name, opts)

    def create_parser(self, name):
        p = argparse.ArgumentParser(
            prog = name,
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description = '''process bowtie-aligned bam files to get reads;
            specifically for the CASTx129 slide-seq data'''
        )
        p.add_argument(
            '-i', '--input',
            help = 'input bam', required=True
        )
        p.add_argument(
            '-o', '--output',
            help = 'name of the output csv', required=True
        )
        p.add_argument(
            '--transcript_mode',
            help = 'save transcripts instead of gene names',
            action='store_true'
        )
        p.add_argument(
            '--log_file',
            help = 'write log messages to file'
        )
        return p

    def main(self, name, opts):
        logging.basicConfig(filename = opts.log_file,
                            format = '%(levelname)s (%(asctime)s): %(message)s')
        log = logging.getLogger(name)
        log.setLevel(logging.INFO)

        infile = pysam.AlignmentFile(opts.input, 'rb')
        tx_gene = pd.read_csv('ensembl_tx_gene.csv')
        tx_gene = tx_gene.set_index('ensembl_transcript_id')
        tx_gene = tx_gene.to_dict('index')
        log.info('Parsing reads from pre-mapped bam')
        results = {} # dictionary of beads with dictionary of transcripts with dictionary of umis
        total_reads = 0
        aligned_reads = 0
        uniquely_aligned_reads = 0
        transcript = []
        for s in infile:
            if s.flag < 200: # then it's starting a new read
                total_reads += 1
                if total_reads % 2000000 == 0:
                    log.info('processed {} total reads'.format(total_reads))
                if len(transcript) > 0: # need to write a previous entry before continuing
                    contains_CAST = False
                    contains_129 = False
                    for t in transcript:
                        if 'CAST' in t:
                            contains_CAST = True
                        if '129' in t:
                            contains_129 = True
                    if not (contains_CAST and contains_129): # if it only has one
                        # loop through and figure out which gene(s) they correspond to
                        genes = []
                        txs = []
                        for t in transcript:
                            tx, allele = t.split('_')
                            if tx in tx_gene:
                                txs.append(tx)
                                genes.append(tx_gene[tx]['external_gene_name'])
                        if len(genes) > 0:
                            if opts.transcript_mode:
                                mycondition = all(x == genes[0] for x in genes) & all(x == txs[0] for x in txs)
                            else:
                                mycondition = all(x == genes[0] for x in genes)
                            if mycondition: # if the transcripts go to the same gene (or in transcript mode, if all go to the same transcript)
                                uniquely_aligned_reads += 1
                                if opts.transcript_mode: # save transcripts
                                    name = txs[0]
                                else:
                                    name = genes[0] # save genes
                                if not bead in results:
                                    results[bead] = {}
                                    results[bead][name] = {'CAST': [], '129': []}
                                    results[bead][name][allele] = [umi]
                                else:
                                    if not name in results[bead]:
                                        results[bead][name] = {'CAST': [], '129': []}
                                        results[bead][name][allele] = [umi]
                                    else:
                                        if not umi in results[bead][name][allele]: # unique umi
                                            results[bead][name][allele].append(umi)
                if s.reference_id != -1: # then it aligned
                    aligned_reads +=1
                    mismatches = [x[1] for x in s.tags if x[0] == 'NM']
                    if mismatches[0] <= 3: # only take those with 3 or less mismatches to reference
                        bead = s.query_name.split(':')[6]
                        umi = s.query_name.split(':')[8]
                        transcript = [infile.getrname(s.reference_id)]
                    else:
                        transcript = []
                else: # then it didn't align
                    transcript = []
            else: # then it's reporting additional alignments
                if ([x[1] for x in s.tags if x[0] == 'NM'][0] <= mismatches[0]) and (mismatches[0] <= 3):
                    # only save those that have less than or equal to number of mismatches as first
                    transcript.append(infile.getrname(s.reference_id))

        log.info('fraction aligned: {}, fraction uniquely aligned: {}'.format(round(aligned_reads/total_reads, 2), round(uniquely_aligned_reads/total_reads, 2)))
        log.info('writing output csv to {}'.format(opts.output))
        f = open(opts.output, 'w')
        f.write('bead,gene,CAST,129\n') # header
        for b in results:
            for g in results[b]:
                f.write('{},{},{},{}\n'.format(b,g,len(results[b][g]['CAST']),len(results[b][g]['129'])))

        f.close()
        log.info('Done')


if __name__ == '__main__':
    app = App()
    app.run(sys.argv)
