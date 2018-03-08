# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Classes to create, access, and/or populate contig links from long
    reads and other resources.
"""

import os
import itertools

import anvio
import anvio.tables as t
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.tables import contiglinks
from anvio.drivers.blast import BLAST

from anvio.errors import ConfigError


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


class LongReadLinks(object):
    def __init__(self, args=None, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.contigs_db_path = A('contigs_db')
        self.long_reads_fasta = A('long_reads_fasta')
        self.min_percent_identity = A('min_percent_identity')
        self.min_alignment_length = A('min_alignment_length')
        self.skip_alignment_position_check = A('skip_alignment_position_check')
        self.output_dir_path = A('output_dir')
        self.num_threads = A('num_threads')
        self.min_e_value = A('min_e_value')


    def check_params(self):
        utils.is_contigs_db(self.contigs_db_path)
        filesnpaths.is_file_fasta_formatted(self.long_reads_fasta)

        if self.output_dir_path:
            if not os.path.exists(self.output_dir_path):
                filesnpaths.gen_output_directory(self.output_dir_path)
            else:
                filesnpaths.is_output_dir_writable(self.output_dir_path)
        else:
            self.output_dir_path = filesnpaths.get_temp_directory_path()

        self.log_file_path = self.get_output_file_path('log.txt')

        if not isinstance(self.min_percent_identity, float):
            raise ConfigError("Minimum percent identity value must be of type float :(")

        if not isinstance(self.min_e_value, float):
            raise ConfigError("Minimum e-value must be of type float :(")

        if not isinstance(self.min_alignment_length, int):
            raise ConfigError("Minimum alignment length must be of type int :(")

        if not isinstance(self.num_threads, int):
            raise ConfigError("Number of threads to use must be of type int :(")

        if self.min_percent_identity < 0 or self.min_percent_identity > 100:
            raise ConfigError("Minimum percent identity must be between 0%% and 100%%. Although your %.2f%% is\
                               pretty cute, too." % self.min_percent_identity)

        if self.min_alignment_length < 20:
            raise ConfigError("Minimum alignment lenght less than 20 does not meet the minimum amount of sense\
                               a decision should make threshold for anvi'o :(")

        if self.min_alignment_length < 1:
            raise ConfigError("Number of threads to use can't be smaller than 1.")


    def check_programs(self):
        utils.is_program_exists('blastn')


    def get_output_file_path(self, file_name, delete_if_exists=False):
        output_file_path = os.path.join(self.output_dir_path, file_name)

        if delete_if_exists:
            if os.path.exists(output_file_path):
                os.remove(output_file_path)

        return output_file_path


    def run_blast(self, query, target):
        blast = BLAST(query, target_fasta=target, search_program='blastn', run=self.run, progress=self.progress)

        blast.log_file_path = self.log_file_path
        blast.num_threads = self.num_threads
        blast.evalue = self.min_e_value
        blast.search_output_path = self.get_output_file_path('blast-search-results.txt')

        return blast.get_blast_results()


    def gen_connections(self, blast_results, contigs_basic_info):
        self.progress.new('Processing search results')
        self.progress.update('...')

        # mapping for the fields in the blast output
        mapping = [str, str, float, int, int, int, int, int, int, int, float, float]

        long_read_commons = {}

        num_hits_passed = 0
        num_hits_considered = 0
        removed_due_to_percent_identity = 0
        removed_due_to_alignment_length = 0
        removed_due_to_position_of_alignment = 0
        for line in open(blast_results):
            fields = line.strip().split('\t')

            query_id, subject_id, perc_id, aln_length, mismatches, gaps, q_start, q_end, s_start, s_end, e_val, bit_score = \
                [mapping[i](fields[i]) for i in range(0, len(mapping))]
            s_length = contigs_basic_info[subject_id]['length']

            num_hits_considered += 1

            if num_hits_considered % 5000 == 0:
                self.progress.update('Lines processed %s ...' % pp(num_hits_considered))

            #
            # FILTERING BASED ON ALIGNMENT LENGTH
            #
            if aln_length < self.min_alignment_length:
                removed_due_to_alignment_length += 1
                continue

            #
            # FILTERING BASED ON PERCENT IDENTITY
            #
            if perc_id < self.min_percent_identity:
                removed_due_to_percent_identity += 1
                continue

            #
            # FILTERING BASED ON POSITION OF ALIGNMENT
            #
            if not self.skip_alignment_position_check:
                if (min(s_start, s_end) - self.min_alignment_length) < 0 or (max(s_start, s_end) + self.min_alignment_length) > s_length:
                    # this is OK. which means the alignment is either like this:
                    pass
                else:
                    removed_due_to_position_of_alignment += 1
                    continue

            if anvio.DEBUG:
                print("OKK %06d: %s vs %s: %.2f%% identity over %d nts. Contig length, aln start, aln end: %d:%d:%d" \
                                    % (num_hits_considered, query_id, subject_id, perc_id, aln_length, s_length, s_start, s_end))

            if query_id not in long_read_commons:
                long_read_commons[query_id] = set([])

            long_read_commons[query_id].add(subject_id)

            num_hits_passed += 1

        self.progress.end()

        self.run.warning(None, 'Hits', lc='cyan')
        self.run.info('Hits considered', '%s' % pp(num_hits_considered))
        self.run.info('Num hits removed due to percent identity', '%s' % pp(removed_due_to_percent_identity), mc='red')
        self.run.info('Num hits removed due to min alignment length', '%s' % pp(removed_due_to_alignment_length), mc='red')
        self.run.info('Num hits removed due to alignment position', '%s' % pp(removed_due_to_position_of_alignment), mc='red')
        self.run.info('Hits kept', '%s' % pp(num_hits_passed))

        self.run.warning(None, 'Connections', lc='cyan')
        self.run.info('Num long reads with contig hits', '%s' % pp(len(long_read_commons)))

        # remove long reads that has a hit only a single contig
        singletons = set([q_id for q_id in long_read_commons if len(long_read_commons[q_id]) == 1])
        for q_id in singletons:
            long_read_commons.pop(q_id)

        long_read_ids_sorted_by_num_connections = sorted([(len(long_read_commons[long_read_id]), long_read_id) for long_read_id in long_read_commons], reverse=True)
        num_connecting_contigs = [e[0] for e in long_read_ids_sorted_by_num_connections]

        self.run.info('Num singletons removed', '%s' % pp(len(singletons)), mc='red')
        self.run.info('Total num contig connections', sum(num_connecting_contigs))
        self.run.info('Min num connected contigs', min(num_connecting_contigs))
        self.run.info('Max num connected contigs', max(num_connecting_contigs))

        if anvio.DEBUG:
            self.run.warning(None, 'DEBUG: Top connecting long reads', lc='green')
            for i in range(0, 10):
                e = long_read_ids_sorted_by_num_connections[i]
                self.run.info(e[1], e[0], nl_after=(1 if i == 9 else 0))

        connections = set()

        for v in long_read_commons.values():
            # here we add a weight of 1 to all links
            for tpl in [tuple(sorted(x) + [1]) for x in itertools.combinations(v, r=2)]:
                connections.add(tpl)

        self.run.info('Final num connections', pp(len(connections)), nl_after=1)

        return connections


    def sanity_check(self):
        self.check_programs()

        self.check_params()

        self.run.log_file_path = self.log_file_path
        self.run.info('Args', (str(self.args)), quiet=True)


    def process(self):
        # check sanity
        self.sanity_check()

        self.run.info('Output directory', self.output_dir_path)
        self.run.info('Contigs database', self.contigs_db_path)
        self.run.info('Long reads', self.long_reads_fasta)
        self.run.info('Min percent identity', self.min_percent_identity)
        self.run.info('Min alignment length', self.min_alignment_length)
        self.run.info('Min e-value', self.min_e_value)

        # exporting contigs.
        contigs_fasta_path = self.get_output_file_path('contigs.fa')

        if os.path.exists(contigs_fasta_path):
            self.run.warning("Notice: A BLAST database is found in the output directory, and will be used!")
        else:
            self.progress.new('Exporting contigs')
            self.progress.update('...')
            utils.export_sequences_from_contigs_db(self.contigs_db_path, contigs_fasta_path, splits_mode=False)
            self.progress.end()

        ## run search
        BLAST_results = self.run_blast(self.long_reads_fasta, contigs_fasta_path)

        contigs_db = dbops.ContigsDatabase(self.contigs_db_path)
        contigs_basic_info = contigs_db.db.get_table_as_dict(t.contigs_info_table_name, string_the_key=True)

        # filter search results and generate a long_read_commons dict
        connections = self.gen_connections(BLAST_results, contigs_basic_info)

        # populate the table
        table_for_contig_links = contiglinks.TableForContigLinks(self.contigs_db_path)
        table_for_contig_links.populate_contig_links('test', connections)

        self.run.quit()


