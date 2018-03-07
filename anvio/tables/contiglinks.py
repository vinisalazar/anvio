# -*- coding: utf-8
# pylint: disable=line-too-long

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.tables.tableops import Table


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Ozcan C. Esen"
__email__ = "ozcanesen@gmail.com"
__status__ = "Development"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


class TableForContigLinks(Table):
    def __init__(self, db_path, run=run, progress=progress):
        self.db_path = db_path

        utils.is_contigs_db(self.db_path)

        Table.__init__(self, self.db_path, anvio.__contigs__version__, run, progress)

        self.set_next_available_id(t.contig_links_table_name)
        self.set_next_available_id(t.contig_links_table_name)


    def remove_source(self, source):
        self.delete_entries_for_key('source', source, [t.contig_links_table_name])


    def populate_contig_links(self, source, connections_list):
        self.remove_source(source)

        db_entries = []
        for connection in connections_list:
            db_entries.append(tuple([self.next_id(t.contig_links_table_name)] + list(connection) + [source]))

        database = db.DB(self.db_path, utils.get_required_version_for_db(self.db_path))

        database._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?)''' % t.contig_links_table_name, db_entries)

        database.disconnect()
