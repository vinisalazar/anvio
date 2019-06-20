# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    A module to dealing with genome storages.

    Pangenomic workflow heavily uses this module.

    Ad hoc access to make sene of internal or external genome descriptions is also welcome.
"""

import time
import hashlib
import collections

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.fastalib as fastalib
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print



class GenomeStorage():
    def __init__(self, db_path, create_new=False, genome_names_to_focus=None, skip_init_functions=False, progress=progress, run=run):
        self.db_path = db_path
        self.db_type = 'genomestorage'
        self.version = anvio.__genomes_storage_version__
        self.genome_names_to_focus = genome_names_to_focus
        self.skip_init_functions = skip_init_functions

        self.progress = progress
        self.run = run

        # will be populated by init()
        self.genome_names = None

        if create_new:
            self.check_storage_path_for_create_new()

        self.db = db.DB(self.db_path, self.version, new_database=create_new)

        Attributes = collections.namedtuple('Attributes', 'primary_key has_numeric_key filter filter_source')

        # IMPORTANT: keys of this dictionary are not same as table names. as tables/__init__.py 
        # contains inconsistencies between variable names and table names, keys here follows variable names.
        # and will be referred as table_id inside this class.  
        self.tables = collections.OrderedDict({
                'contig_sequences': Attributes(primary_key='contig', has_numeric_key=False, 
                                                filter='contig', filter_source='contig'),

                'contigs_info':     Attributes(primary_key='contig', has_numeric_key=False,
                                                filter='contig', filter_source='contig'),

                'splits_info':      Attributes(primary_key='split', has_numeric_key=False, 
                                                filter='split', filter_source='split'),

                'genes_in_contigs': Attributes(primary_key='gene_callers_id', has_numeric_key=False, 
                                                filter='gene_callers_id', filter_source='gene_callers_id'),

                'genes_in_splits':  Attributes(primary_key='entry_id', has_numeric_key=True,
                                                filter='gene_callers_id', filter_source='gene_callers_id'),

                'gene_amino_acid_sequences':  Attributes(primary_key='gene_callers_id', has_numeric_key=False,
                                                        filter='gene_callers_id', filter_source='gene_callers_id'),

                'gene_function_calls':  Attributes(primary_key='entry_id', has_numeric_key=True,
                                                    filter='gene_callers_id', filter_source='gene_callers_id'),

                'hmm_hits_info':  Attributes(primary_key='source', has_numeric_key=False, filter=None, 
                                                filter_source=None),

                'hmm_hits':  Attributes(primary_key='entry_id', has_numeric_key=True, 
                                        filter='gene_callers_id', filter_source='gene_callers_id'),

                'hmm_hits_splits': Attributes(primary_key='entry_id', has_numeric_key=True, 
                                            filter='split', filter_source='split'),

                'nt_position_info': Attributes(primary_key='contig_name', has_numeric_key=False, 
                                            filter='contig_name', filter_source='contig'),

                'splits_taxonomy': Attributes(primary_key='split', has_numeric_key=False, filter='split', 
                                            filter_source='split'),

                'genes_taxonomy': Attributes(primary_key='gene_callers_id', has_numeric_key=False, 
                                            filter='gene_callers_id', filter_source='gene_callers_id'),

                'taxon_names': Attributes(primary_key='taxon_id', has_numeric_key=False, filter=None, 
                                            filter_source=None),
            })

        self.next_available_id = {}

        if create_new:
            self.create_tables()
            self.set_meta_values()
        else:
            self.init()

        self.populate_next_available_ids()


    def check_storage_path_for_create_new(self):
        if not self.db_path.endswith('GENOMES.db'):
            raise ConfigError("The genomes storage file must end with '-GENOMES.db'. Anvi'o developers do know how ridiculous\
                                this requirement sounds like, but if you have seen the things they did, you would totally\
                                understand why this is necessary.")

        filesnpaths.is_output_file_writable(self.db_path)


    def populate_next_available_ids(self):
        for table_id, attributes in self.tables.items():
            if attributes.has_numeric_key:
                name, structure, types = self.get_table_defs(table_id)

                self.next_available_id[table_id] = self.db.get_max_value_in_column(name, 
                    attributes.primary_key, value_if_empty=-1) + 1

    def init(self):
        genome_names_in_db = self.get_all_genome_names()

        if self.genome_names_to_focus:
            genome_names_to_focus_missing_from_db = [g for g in self.genome_names_to_focus if g not in genome_names_in_db]

            # make sure the user knows what they're doing
            if genome_names_to_focus_missing_from_db:
                raise ConfigError("%d of %d genome names you wanted to focus are missing from the genomes sotrage.\
                                 Although this may not be a show-stopper, anvi'o likes to be explicit, so here we\
                                 are. Not going anywhere until you fix this. For instance this is one of the missing\
                                 genome names: '%s', and this is one random genome name from the database: '%s'" % \
                                         (len(genome_names_to_focus_missing_from_db), len(self.genome_names_to_focus),\
                                         genome_names_to_focus_missing_from_db[0], ', '.join(genome_names_in_db)))

            self.genome_names = self.genome_names_to_focus
        else:
            self.genome_names = genome_names_in_db

        ## load the data
        self.progress.new('Recovering data from the db')
        self.progress.update('Loading genomes basic info...')
        where_clause = """genome_name IN (%s)""" % ",".join('"' + item + '"' for item in self.genome_names)
        self.genomes_info = self.db.get_some_rows_from_table_as_dict(t.genome_info_table_name, where_clause)

        self.gene_info = {}
        self.progress.update('Loading genes info for %s genomes...' % len(self.genomes_info))
        
        gene_aa_sequences_table_name, _, _ = self.get_table_defs('gene_amino_acid_sequences')
        gene_table_name, _, _ = self.get_table_defs('genes_in_contigs')
        gene_functions_table_name, _, _ = self.get_table_defs('gene_function_calls')

        gene_aa_sequences = self.db.get_some_rows_from_table_as_dict(gene_aa_sequences_table_name, where_clause)

        for gene_info_tuple in self.db.get_some_rows_from_table(gene_table_name, where_clause):
            gene_callers_id, contig, start, stop, direction, partial, version, source, genome_name = gene_info_tuple
            length = stop - start
            
            if genome_name not in self.gene_info:
                self.gene_info[genome_name] = {}

            self.gene_info[genome_name][gene_callers_id] = {
                'aa_sequence': gene_aa_sequences[gene_callers_id],
                'partial': partial,
                'length': length,
                'functions': {}
            }


        if not self.skip_init_functions:
            for gene_function_tuple in self.db.get_some_rows_from_table(gene_functions_table_name, where_clause):
                entry_id, gene_callers_id, source, accession, function, e_value, genome_name = gene_function_tuple
                
                if gene_callers_id not in self.gene_info[genome_name]:
                    continue

                function_dict = {
                    'accession': accession,
                    'function': function,
                    'e_value': e_value
                }

                self.gene_info[genome_name][gene_callers_id]['functions'][source] = function_dict

        self.progress.end()


    def get_table_defs(self, table_id):
        table_name      =   getattr(t, "%s_table_name"      % table_id)
        table_structure = [*getattr(t, "%s_table_structure" % table_id)]
        table_types     = [*getattr(t, "%s_table_types"     % table_id)]

        return (table_name, table_structure, table_types)


    def store_genomes(self, description):
        for name, genome_desc in description.genomes.items():
            self.add_genome(genome_desc)


    def add_genome(self, genome):
        if genome['name'] in self.get_all_genome_names():
            raise ConfigError("It seems the genome '%s' you are trying to add, already exists \
                               in the genome storage. If you see this message while creating \
                               a genome storage with multiple internal or external genomes make sure\
                               that you don't have duplicate entries." % genome['name'])
        
        source_db = db.DB(genome['contigs_db_path'], anvio.__contigs__version__)

        for table_id, attributes in self.tables.items():
            self.import_table(source_db, table_id, attributes, genome)

        # add genome information to genome_info table
        values = (genome['name'], )

        for column_name in t.genome_info_table_structure[1:]:
            if genome[column_name]:
                values += (genome[column_name], )
            else:
                # the following line will add a -1 for any `key` that has the value of `None`. the reason
                # we added this was to be able to work with contigs databases without any hmm hits for SCGs
                # which is covered in https://github.com/merenlab/anvio/issues/573
                values += (-1, )

        self.db.insert(t.genome_info_table_name, values=values)

        self.update_storage_hash()


    def import_table(self, source_db, table_id, attributes, genome):
        table_name, table_structure, table_types = self.get_table_defs(table_id)
        primary_key_index = table_structure.index(attributes.primary_key)

        filters = {
            'contig': genome['contigs_of_interest'],
            'split': genome['splits_of_interest'],
            'gene_callers_id': genome['gene_caller_ids'],
        }

        # it will be "WHERE 1" and match all, so no filter
        where_clause = '1' 

        if attributes.filter_source in filters:
            first_item = next(iter(filters[attributes.filter_source]))
            is_string = isinstance(first_item, str)

            if is_string:
                # add quotes
                keys = map(lambda x: "'%s'" % x, filters[attributes.filter_source])
            else:
                keys = map(str, filters[attributes.filter_source])

            where_clause = '%s IN (%s)' % (attributes.filter,
                                           ",".join(keys))

        table_content = source_db.get_some_rows_from_table(table_name, where_clause=where_clause)
        entries = []

        while table_content:
            # here we convert row entry tuple to array
            # because we may need to modify the primary key or add a genome_name column value.
            # we will later convert it back to tuple
            entry = [*table_content.pop()]

            if attributes.has_numeric_key:
                entry[primary_key_index] = self.next_available_id[table_id]

                self.next_available_id[table_id] += 1

            # add value for genome_name column
            entry.append(genome['name'])

            entries.append(tuple(entry))

        if len(entries):
            self.db._exec_many('''INSERT INTO %s VALUES(%s)''' % (table_name, ','.join(['?'] * len(entries[0]))), entries)
            self.db.commit()


    def create_tables(self):
        self.db.create_table(t.genome_info_table_name, t.genome_info_table_structure, t.genome_info_table_types)

        for table_id, attributes in self.tables.items():
            table_name, table_structure, table_types = self.get_table_defs(table_id)

            table_structure += ['genome_name']
            table_types     += [   'text'    ]

            self.db.create_table(table_name, table_structure, table_types)

            # create indexes
            if attributes.primary_key:
                self.db._exec("CREATE INDEX %s_primary_index ON %s (%s);" % 
                               (table_name, table_name, attributes.primary_key))

                self.db._exec("CREATE INDEX %s_primary_and_genome_index ON %s (%s, genome_name);" % 
                               (table_name, table_name, attributes.primary_key))


    def set_meta_values(self):
        self.db.set_meta_value('hash', 'hash_EMPTY_DATABASE')
        self.db.set_meta_value('db_type', self.db_type)


    def get_all_genome_names(self):
        return self.db.get_single_column_from_table(t.genome_info_table_name, 'genome_name')


    def get_genomes_dict(self):
        # we retrieve all table at once to avoid seperate sql queries
        all_genomes_dict = self.db.get_table_as_dict(t.genome_info_table_name)
        result = {}

        # copy genomes requested by user to result dictionary
        for genome_name in self.genome_names:
            result[genome_name] = all_genomes_dict[genome_name]

        return result




    def update_storage_hash(self):
        # here we create a signature for the storage itself by concatenating all hash values from all genomes. even if one
        # split is added or removed to any of these genomes will change this signature. since we will tie this information
        # to the profile database we will generate for the pangenome analysis, even if one split is added or removed from any
        # of the genomes will make sure that the profile databases from this storage and storage itself are not compatible:

        concatenated_genome_hashes = '_'.join(sorted(map(str, self.db.get_single_column_from_table(t.genome_info_table_name, 'genome_hash'))))
        new_hash = 'hash' + str(hashlib.sha224(concatenated_genome_hashes.encode('utf-8')).hexdigest()[0:8])

        self.db.set_meta_value('hash', new_hash)


    def close(self):
        self.db.disconnect()
