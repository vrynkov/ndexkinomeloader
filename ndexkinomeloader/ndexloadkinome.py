#! /usr/bin/env python

import argparse
import sys
import logging
from logging import config
from ndexutil.config import NDExUtilConfig
import ndexkinomeloader

import requests
import os
import zipfile

import csv

SUCCESS = 0
ERROR = 2

logger = logging.getLogger(__name__)

TSV2NICECXMODULE = 'ndexutil.tsv.tsv2nicecx2'

LOG_FORMAT = "%(asctime)-15s %(levelname)s %(relativeCreated)dms " \
             "%(filename)s::%(funcName)s():%(lineno)d %(message)s"

def _parse_arguments(desc, args):
    """
    Parses command line arguments
    :param desc:
    :param args:
    :return:
    """
    help_fm = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=help_fm)

    parser.add_argument('datadir', help='Directory where BioGRID Kinome data downloaded to and processed from')

    parser.add_argument('--profile', help='Profile in configuration '
                                          'file to use to load '
                                          'NDEx credentials which means'
                                          'configuration under [XXX] will be'
                                          'used '
                                          '(default '
                                          'ndexkinomeloader)',
                        default='ndexkinomeloader')
    parser.add_argument('--logconf', default=None,
                        help='Path to python logging configuration file in '
                             'this format: https://docs.python.org/3/library/'
                             'logging.config.html#logging-config-fileformat '
                             'Setting this overrides -v parameter which uses '
                             ' default logger. (default None)')

    parser.add_argument('--conf', help='Configuration file to load '
                                       '(default ~/' +
                                       NDExUtilConfig.CONFIG_FILE)
    parser.add_argument('--verbose', '-v', action='count', default=0,
                        help='Increases verbosity of logger to standard '
                             'error for log messages in this module and'
                             'in ' + TSV2NICECXMODULE + '. Messages are '
                             'output at these python logging levels '
                             '-v = ERROR, -vv = WARNING, -vvv = INFO, '
                             '-vvvv = DEBUG, -vvvvv = NOTSET (default no '
                             'logging)')
    parser.add_argument('--version', action='version',
                        version=('%(prog)s ' + ndexkinomeloader.__version__))

    parser.add_argument('--biogridversion', help='Version of BioGRID Release', default='3.5.177')

    parser.add_argument('--skipdownload', action='store_true',
                        help='If set, skips download of  BioGRID Kinome and assumes data already reside in <datadir>'
                             'directory')

    return parser.parse_args(args)


def _setup_logging(args):
    """
    Sets up logging based on parsed command line arguments.
    If args.logconf is set use that configuration otherwise look
    at args.verbose and set logging for this module and the one
    in ndexutil specified by TSV2NICECXMODULE constant
    :param args: parsed command line arguments from argparse
    :raises AttributeError: If args is None or args.logconf is None
    :return: None
    """

    if args.logconf is None:
        level = (50 - (10 * args.verbose))
        logging.basicConfig(format=LOG_FORMAT,
                            level=level)
        logging.getLogger(TSV2NICECXMODULE).setLevel(level)
        logger.setLevel(level)
        return

    # logconf was set use that file
    logging.config.fileConfig(args.logconf,
                              disable_existing_loggers=False)


class NDExNdexkinomeloaderLoader(object):
    """
    Class to load content
    """
    def __init__(self, args):
        """

        :param args:
        """
        self._conf_file = args.conf
        self._profile = args.profile
        self._user = None
        self._pass = None
        self._server = None

        self._ndex = None

        self._biogrid_version = args.biogridversion
        self._datadir = os.path.abspath(args.datadir)
        self._skipdownload = args.skipdownload

        self._kinome_zip = os.path.join(self._datadir, self._get_kinome_zip_file_name())
        self._interactions = self._get_interactions_file_name()
        self._ptm = self._get_ptm_file_name()
        self._genes = self._get_genes_file_name()
        self._relations = self._get_relations_file_name()


        self._interaction_headers = ["#BIOGRID ID", "ENTREZ GENE ID", "INTERACTION COUNT", "PTM COUNT",
                   "CHEMICAL INTERACTION COUNT", "SOURCE", "CATEGORY VALUES", "SUBCATEGORY VALUES"]
        self._gene_lookup = {}



    #'BIOGRID-PROJECT-kinome_project_sc-INTERACTIONS-3.5.177.tab2.txt'
    #'BIOGRID-PROJECT-kinome_project_sc-PTM-3.5.177.ptmtab.txt'
    #'BIOGRID-PROJECT-kinome_project_sc-GENES-3.5.177.projectindex.txt'

    #'BIOGRID-PROJECT-kinome_project_sc-PTM-RELATIONSHIPS-3.5.177.ptmrel.txt'





    def _parse_config(self):
            """
            Parses config
            :return:
            """
            ncon = NDExUtilConfig(conf_file=self._conf_file)
            con = ncon.get_config()
            self._user = con.get(self._profile, NDExUtilConfig.USER)
            self._pass = con.get(self._profile, NDExUtilConfig.PASSWORD)
            self._server = con.get(self._profile, NDExUtilConfig.SERVER)


    def _get_kinome_prefix(self):
        return 'BIOGRID-PROJECT-kinome_project_sc-'

    def _get_kinome_zip_file_name(self):
        return self. _get_kinome_prefix() + self._biogrid_version + '.zip'

    def _get_interactions_file_name(self):
        return os.path.join(self._datadir, \
                            self._get_kinome_prefix() + 'INTERACTIONS-' + self._biogrid_version + '.tab2.txt')

    def _get_ptm_file_name(self):
        return os.path.join(self._datadir,\
                            self._get_kinome_prefix() + 'PTM-' + self._biogrid_version + '.ptmtab.txt')

    def _get_genes_file_name(self):
        return os.path.join(self._datadir, \
                            self._get_kinome_prefix() + 'GENES-' + self._biogrid_version + '.projectindex.txt')

    def _get_relations_file_name(self):
        return os.path.join(self._datadir, \
                            self._get_kinome_prefix() + 'PTM-RELATIONSHIPS-' + self._biogrid_version + '.ptmrel.txt')


    def _get_kinome_download_url(self):
        return 'https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-' + \
            self._biogrid_version + '/' + self._get_kinome_zip_file_name()


    def _download_file(self, url):

        #if not os.path.exists(self._datadir):
        #    os.makedirs(self._datadir)
        try:
            print(url)
            response = requests.get(url)

            if response.status_code // 100 == 2:
                with open(self._kinome_zip, "wb") as received_file:
                    received_file.write(response.content)
            else:
                return response.status_code

        except requests.exceptions.RequestException as e:
            logger.exception('Caught exception')
            print('\n\n\tException: {}\n'.format(e))
            return ERROR

        return SUCCESS


    def _download_kinome_files(self):
        url = self._get_kinome_download_url()
        download_status = self._download_file(url)
        return download_status


    def _check_if_data_dir_exists(self):
        data_dir_existed = True

        if not os.path.exists(self._datadir):
            data_dir_existed = False
            os.makedirs(self._datadir, mode=0o755)

        return data_dir_existed


    def _unzip_kinome(self):
        try:
            with zipfile.ZipFile(self._kinome_zip, "r") as zip_ref:
                zip_ref.extractall(self._datadir)
        except Exception as e:
            print('\n\n\tException: {}\n'.format(e))
            return ERROR

        return SUCCESS


    def _build_gene_lookup(self):

        try:
            _gene_lookup = pd.read_csv(self._genes, sep='\t')
        except:
            return ERROR

        for index, row in _gene_lookup.iterrows():
            self._gene_lookup[str(row['ENTREZ GENE ID'])] = \
                {
                    'INTERACTION COUNT': row['INTERACTION COUNT'],
                    'PTM COUNT': row['PTM COUNT'],
                    'CHEMICAL INTERACTION COUNT': row['CHEMICAL INTERACTION COUNT'],
                    'SOURCE': row['SOURCE'],
                    'CATEGORY VALUES': row['CATEGORY VALUES'],
                    'SUBCATEGORY VALUES': row['SUBCATEGORY VALUES']
                }

        return SUCCESS


    def _build_gene_tsv(self, entrez_gene_A_data, entrez_gene_B_data):
        ret_array = []

        ret_array.append(entrez_gene_A_data['INTERACTION COUNT'])
        ret_array.append(entrez_gene_B_data['INTERACTION COUNT'])
        ret_array.append(entrez_gene_A_data['PTM COUNT'])
        ret_array.append(entrez_gene_B_data['PTM COUNT'])
        ret_array.append(entrez_gene_A_data['CHEMICAL INTERACTION COUNT'])
        ret_array.append(entrez_gene_B_data['CHEMICAL INTERACTION COUNT'])
        ret_array.append(entrez_gene_A_data['SOURCE'])
        ret_array.append(entrez_gene_B_data['SOURCE'])
        ret_array.append(entrez_gene_A_data['CATEGORY VALUES'])
        ret_array.append(entrez_gene_B_data['CATEGORY VALUES'])
        ret_array.append(entrez_gene_A_data['SUBCATEGORY VALUES'])
        ret_array.append(entrez_gene_B_data['SUBCATEGORY VALUES'])

        ret_str = '\t'.join(str(element)  if element != '-' else '' for element in ret_array)
        return ret_str


    def _create_ppi_file(self):
        interactions_header = \
            ['#BioGRID Interaction ID', 'Entrez Gene Interactor A', 'Entrez Gene Interactor B',
             'BioGRID ID Interactor A','BioGRID ID Interactor B', 'Systematic Name Interactor A',
             'Systematic Name Interactor B', 'Official Symbol Interactor A', 'Official Symbol Interactor B',
             'Synonyms Interactor A', 'Synonyms Interactor B', 'Experimental System', 'Experimental System Type',
             'Author', 'Pubmed ID', 'Organism Interactor A', 'Organism Interactor B', 'Throughput', 'Score',
             'Modification', 'Phenotypes', 'Qualifications', 'Tags', 'Source Database'
             ]

        new_headers =  ['Interaction Count A', 'Interaction Count B',
                        'PTM Count A', 'PTM Count B',
                        'Chemical Interaction Count A', 'Chemical Interaction Count B',
                        'Source A', 'Source B', 'Category Values A', 'Category Values B',
                        'SubCategory Values A', 'SubCategory Values B']

        interactions_header_1 = interactions_header + new_headers

        ppi_network = os.path.join(self._datadir, 'ppi_network_1.txt')

        default_gene_data = {'INTERACTION COUNT':'',
                             'PTM COUNT':'',
                             'CHEMICAL INTERACTION COUNT': '',
                             'SOURCE': '',
                             'CATEGORY VALUES': '',
                             'SUBCATEGORY VALUES': ''}
        try:
            with open(self._interactions, 'r') as tsv:
                reader = csv.reader(tsv, delimiter='\t')


                with open(ppi_network, 'w') as o_f:

                    output_header = '\t'.join(h for h in interactions_header_1) + '\n'
                    o_f.write(output_header)

                    for row in reader:
                        # skip header since we already wrote it to output
                        break

                    for row in reader:
                        entrez_gene_interactor_a = row[1]
                        entrez_gene_interactor_b = row[2]

                        entrez_gene_A_data = self._gene_lookup.get(entrez_gene_interactor_a, default_gene_data)
                        entrez_gene_B_data = self._gene_lookup.get(entrez_gene_interactor_b, default_gene_data)

                        gene_tsv = self._build_gene_tsv(entrez_gene_A_data, entrez_gene_B_data)

                        output_tsv = '\t'.join(element if element != '-' else '' for element in row) + '\t' + gene_tsv + '\n'
                        o_f.write(output_tsv)

        except:
            return ERROR

        return SUCCESS


    def _create_ptm_file(self):
        ptm_header = \
            ['#PTM ID', 'Entrez Gene ID', 'BioGRID ID', 'Systematic Name', 'Official Symbol',
             'Synonymns', 'Sequence', 'Refseq ID', 'Position', 'Post Translational Modification',
             'Residue', 'Author', 'Pubmed ID', 'Organism ID', 'Organism Name', 'Has Relationships',
             'Notes', 'Source Database']

        new_header = ['Target Name', 'Target Represents']

        ptm_header_1 = ptm_header + new_header

        ppi_network = os.path.join(self._datadir, 'ptm_network_2.txt')

        try:
            with open(self._ptm, 'r') as tsv:
                reader = csv.reader(tsv, delimiter='\t')

                with open(ppi_network, 'w') as o_f:

                    for row in reader:
                        output_tsv = '\t'.join(e for e in row) + '\t' + '\t'.join(e for e in new_header) + '\n'
                        o_f.write(output_tsv)
                        # skip header since we already wrote it to output
                        break

                    for row in reader:
                        target_name = str(row[10]) + str(row[8])
                        target_represents = row[4] + '-' + str(row[10]) + '-' + str(row[8])

                        output_tsv = '\t'.join(e if e != '-' else '' for e in row) + '\t' + \
                                     target_name + '\t' + target_represents + '\n'

                        o_f.write(output_tsv)

        except:
            return ERROR

        return SUCCESS


    def run(self):
        """
        Runs content loading for NDEx KINOME Content Loader
        :param theargs:
        :return:
        """
        self._parse_config()

        data_dir_existed = self._check_if_data_dir_exists()

        if self._skipdownload is False or data_dir_existed is False:
            status_code = self._download_kinome_files()
            if status_code != 0:
                return ERROR

            status_code = self._unzip_kinome()
            if status_code != 0:
                return ERROR

        self._build_gene_lookup()

        # Step 1 - create PPI file from GENES and INTERACTIONS files
        self._create_ppi_file()

        # Step 2 - create PTM network file
        self._create_ptm_file()

        return SUCCESS


def main(args):
    """
    Main entry point for program
    :param args:
    :return:
    """
    desc = """
    Version {version}

    Loads NDEx KINOME Content Loader data into NDEx (http://ndexbio.org).
    
    To connect to NDEx server a configuration file must be passed
    into --conf parameter. If --conf is unset the configuration 
    the path ~/{confname} is examined. 
         
    The configuration file should be formatted as follows:
         
    [<value in --profile (default ncipid)>]
         
    {user} = <NDEx username>
    {password} = <NDEx password>
    {server} = <NDEx server(omit http) ie public.ndexbio.org>
    
    
    """.format(confname=NDExUtilConfig.CONFIG_FILE,
               user=NDExUtilConfig.USER,
               password=NDExUtilConfig.PASSWORD,
               server=NDExUtilConfig.SERVER,
               version=ndexkinomeloader.__version__)
    theargs = _parse_arguments(desc, args[1:])
    theargs.program = args[0]
    theargs.version = ndexkinomeloader.__version__

    try:
        _setup_logging(theargs)
        loader = NDExNdexkinomeloaderLoader(theargs)
        return loader.run()
    except Exception as e:
        logger.exception('Caught exception')
        return 2
    finally:
        logging.shutdown()


if __name__ == '__main__':  # pragma: no cover
    sys.exit(main(sys.argv))
