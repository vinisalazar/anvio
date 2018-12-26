# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Classes to define and work with anvi'o metapangenomics workflows.
"""


import anvio
import anvio.terminal as terminal

from anvio.workflows import WorkflowSuperClass
from anvio.workflows.contigs import ContigsDBWorkflow
from anvio.workflows.metagenomics import MetagenomicsWorkflow
from anvio.workflows.pangenomics import PangenomicsWorkflow


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Alon Shaiber"
__email__ = "alon.shaiber@gmail.com"


run = terminal.Run()
progress = terminal.Progress()


class MetaPangenomicsWorkflow(MetagenomicsWorkflow, PangenomicsWorkflow, ContigsDBWorkflow, WorkflowSuperClass):
    def __init__(self, args=None, run=terminal.Run(), progress=terminal.Progress()):
        self.init_workflow_super_class(args, workflow_name='metapangenomics')

        self.run = run
        self.progress = progress

        # Some items from inherited workflows are not relevant in metapangenomics
        self.config_params_to_remove = ['external_genomes', 'megahit', 'idba_ud', \
                                        'metaspades', 'merge_fastqs_for_co_assembly', \
                                        'merge_fastas_for_co_assembly', 'references_mode']

        # know thyself.
        self.name = 'metapangenomics'

        # initialize the base classes
        PangenomicsWorkflow.__init__(self)
        MetagenomicsWorkflow.__init__(self)

        # remove parameters that are not available for metapangenomics
        for param in self.config_params_to_remove:
            self.rules.clear(param)
            try:
                self.default_config.pop(param)
            except:
                pass

        self.rules.extend([])

        self.general_params.extend([''])

        rule_acceptable_params_dict = {}

        # defining the accesible params per rule
        rule_acceptable_params_dict[''] = []

        # updating the dict from upstream
        self.rule_acceptable_params_dict.update(rule_acceptable_params_dict)

        self.dirs_dict.update({"QC_DIR": "01_QC",
                               "FASTA_DIR": "02_FASTA",
                               "CONTIGS_DIR": "03_CONTIGS",
                               "MAPPING_DIR": "04_MAPPING",
                               "PROFILE_DIR": "05_ANVIO_PROFILE",
                               "MERGE_DIR": "06_MERGED"})

        self.default_config.update({'metapangenome_fastas_txt': 'metapangenome-fastas.txt'})


    def init(self):
        # Metapangenomics is only available for references mode of metagenomics workflow
        self.references_mode = True


    def get_fasta:
        # get_fasta is used when we want the fasta file that we will use
        # for things like anvi-gen-contigs-database and bowtie-build
        # hence we want the merged fasta file


        # NOTICE that we don't need to override the get_raw_fasta function
        # get_raw_fasta is only used for reformat fasta
        # in metapangenomics we want to reformat the individual fasta
        # separately and only then merge them. And this will be done naturally
