#! /usr/bin/env python

import argparse
import sys
import logging
from logging import config
from ndexutil.config import NDExUtilConfig
import ndexkinomeloader

import requests
import os
from ndexutil.tsv.streamtsvloader import StreamTSVLoader
import zipfile

import csv
import json
import pandas as pd
import ndexutil.tsv.tsv2nicecx2 as t2n

import ndex2
from ndex2.client import Ndex2

import re
import copy


SUCCESS = 0
ERROR = 2

STYLE = 'style.cx'
PTI_LOAD_PLAN = 'kinome_interactions-plan.json'
PTM_LOAD_PLAN = 'kinome_ptm-plan.json'

logger = logging.getLogger(__name__)

TSV2NICECXMODULE = 'ndexutil.tsv.tsv2nicecx2'

LOG_FORMAT = "%(asctime)-15s %(levelname)s %(relativeCreated)dms " \
             "%(filename)s::%(funcName)s():%(lineno)d %(message)s"

def get_package_dir():
    """
    Gets directory where package is installed
    :return:
    """
    return os.path.dirname(ndexkinomeloader.__file__)


def get_load_plan(load_plan):
    """
    Gets the load plan stored in this package
    :return: path to file
    :rtype: string
    """
    return os.path.join(get_package_dir(), load_plan)


def get_style():
    """
    Gets the style stored with this package

    :return: path to file
    :rtype: string
    """
    return os.path.join(get_package_dir(), STYLE)


def _parse_arguments(desc, args):
    """
    Parses command line arguments
    :param desc:
    :param args:
    :return:
    """
    help_fm = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=desc, formatter_class=help_fm)

    styling_group = parser.add_mutually_exclusive_group()

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

    parser.add_argument('--loadpti', help='Path to PTI load plan file in json format', default=get_load_plan(PTI_LOAD_PLAN))
    parser.add_argument('--loadptm', help='Path to PTM load plan file in json format', default=get_load_plan(PTM_LOAD_PLAN))

    styling_group.add_argument('--style', help='Path to NDEx CX file to use for styling networks', default=get_style())


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

    styling_group.add_argument('--template',
           help='UUID of network to use for styling networks (the same account where networks are located)')

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

        self._args = args

        self._ndex = None
        self._template_UUID = args.template

        self._biogrid_version = args.biogridversion
        self._datadir = os.path.abspath(args.datadir)
        self._skipdownload = args.skipdownload

        self._kinome_zip = os.path.join(self._datadir, self._get_kinome_zip_file_name())
        self._interactions = self._get_interactions_file_name()
        self._ptm = self._get_ptm_file_name()
        self._genes = self._get_genes_file_name()
        self._relations = self._get_relations_file_name()

        self._ppi_network_1 = os.path.join(self._datadir, 'ppi_network_1.txt')
        self._ptm_network_2 = os.path.join(self._datadir, 'ptm_network_2.txt')


        self._interaction_headers = ["#BIOGRID ID", "ENTREZ GENE ID", "INTERACTION COUNT", "PTM COUNT",
                   "CHEMICAL INTERACTION COUNT", "SOURCE", "CATEGORY VALUES", "SUBCATEGORY VALUES"]
        self._gene_lookup = {}

        self._pti_load_plan = args.loadpti
        self._ptm_load_plan = args.loadptm

        self._ppi_attributes = {}
        self._ptm_attributes = {}

        self._cx_pti = os.path.join(self._datadir, 'pti_1.cx')
        self._cx_ptm = os.path.join(self._datadir, 'ptm_2.cx')
        self._cx_merged = os.path.join(self._datadir, 'merged_3.cx')


    def _get_user_agent(self):
        """
        :return:
        """
        return 'kinome/' + self._biogrid_version


    def _create_ndex_connection(self):
        """
        creates connection to ndex
        :return:
        """
        if self._ndex is None:

            try:
                self._ndex = Ndex2(host=self._server, username=self._user,
                                   password=self._pass, user_agent=self._get_user_agent())
            except Exception as e:
                self._ndex = None

        return self._ndex

    def _load_style_template(self):
        """
        Loads the CX network specified by self._args.style into self._template
        :return:
        """
        self._template = ndex2.create_nice_cx_from_file(os.path.abspath(self._args.style))


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

        #ppi_network = os.path.join(self._datadir, 'ppi_network_1.txt')

        default_gene_data = {'INTERACTION COUNT':'',
                             'PTM COUNT':'',
                             'CHEMICAL INTERACTION COUNT': '',
                             'SOURCE': '',
                             'CATEGORY VALUES': '',
                             'SUBCATEGORY VALUES': ''}
        try:
            with open(self._interactions, 'r') as tsv:
                reader = csv.reader(tsv, delimiter='\t')


                with open(self._ppi_network_1, 'w') as o_f:

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

                        synonyms_interactor_a = '|ncbigene:' + row[1] + '|' + row[5]
                        synonyms_interactor_b = '|ncbigene:' + row[2] + '|' + row[6]

                        row[9] = row[9] + synonyms_interactor_a
                        row[10] = row[10] + synonyms_interactor_b

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

        try:
            with open(self._ptm, 'r') as tsv:
                reader = csv.reader(tsv, delimiter='\t')

                with open(self._ptm_network_2, 'w') as o_f:

                    for row in reader:
                        output_tsv = '\t'.join(e for e in row) + '\t' + '\t'.join(e for e in new_header) + '\n'
                        o_f.write(output_tsv)
                        # skip header since we already wrote it to output
                        break

                    for row in reader:
                        position_column_value = str(row[8]).rstrip()
                        if position_column_value == '-':
                            position_column_value = '?'
                            row[8] = 'undefined'

                        target_name = str(row[10]) + position_column_value
                        target_represents = row[4] + '-' + str(row[10]) + '-' + str(row[8])

                        synonyms = '|ncbigene:' + row[1] + '|' + row[3] + '|' + row[7]
                        row[5] = row[5] + synonyms

                        output_tsv = '\t'.join(e if e != '-' else '' for e in row) + '\t' + \
                                     target_name + '\t' + target_represents + '\n'

                        o_f.write(output_tsv)

        except:
            return ERROR

        return SUCCESS


    def _init_network_attributes(self, network, type='pti'):
        if type == 'pti':
            network.set_name('PTI - Step 1')
        elif type == 'ptm':
            network.set_name('PTM - Step 2')
        elif type == 'merged':
            network.set_name('FULLY MERGED - Step 3')

        network.set_network_attribute('prov:wasDerivedFrom', self._get_kinome_download_url())
        network.set_network_attribute('prov:wasGeneratedBy',
                '<a href="https://github.com/vrynkov/ndexkinomeloader" target="_blank">ndexkinomeloader ' \
                + str(ndexkinomeloader.__version__) + '</a>')

        network.set_network_attribute('__iconurl', 'https://home.ndexbio.org/img/biogrid_logo.jpg')

        network.apply_style_from_network(self._template)



    def _generate_CX_file(self, load_plan, network_tsv):

        with open(load_plan, 'r') as lp:
            plan = json.load(lp)

        dataframe = pd.read_csv(network_tsv,
                                dtype=str,
                                na_filter=False,
                                delimiter='\t',
                                engine='python')

        network = t2n.convert_pandas_to_nice_cx_with_load_plan(dataframe, plan)

        return network, SUCCESS


    def _write_nice_cx_to_file(self, network_in_cx, cx_file_path):

        with open(cx_file_path, 'w') as f:
            json.dump(network_in_cx.to_cx(), f, indent=4)



    def _upload_CX(self, path_to_network_in_CX, network_UUID):

        with open(path_to_network_in_CX, 'br') as network_out:
            try:
                if network_UUID is None:
                    self._ndex.save_cx_stream_as_new_network(network_out)
                else:
                    self._ndex.update_cx_network(network_out, network_UUID)

            except Exception as e:
                print(e)
                return ERROR

        return SUCCESS


    def _merge_attributes(self, attribute_list_1, attribute_list_2):

        for attribute1 in attribute_list_1:

            name1 = attribute1['n']

            found = False
            for attribute2 in attribute_list_2:
                if attribute2['n'] == name1:
                    found = True
                    break

            if not found:
                continue

            #if attribute1['v'] == attribute2['v']:
                # attribute with the same name and value; do not add
            #    continue

            if not 'd' in attribute1:
                attribute1['d'] = 'list_of_string'
            elif attribute1['d'] == 'boolean':
                attribute1['d'] = 'list_of_boolean'
            elif attribute1['d'] == 'double':
                attribute1['d'] = 'list_of_double'
            elif attribute1['d'] == 'integer':
                attribute1['d'] = 'list_of_integer'
            elif attribute1['d'] == 'long':
                attribute1['d'] = 'list_of_long'
            elif attribute1['d'] == 'string':
                attribute1['d'] = 'list_of_string'

            if not 'd' in attribute2:
                attribute2['d'] = 'list_of_string'
            elif attribute2['d'] == 'boolean':
                attribute2['d'] = 'list_of_boolean'

            new_list_of_values = []

            if isinstance(attribute1['v'], list):
                for value in attribute1['v']:
                    if (attribute2['d'] == 'list_of_boolean') or (value not in new_list_of_values):
                        new_list_of_values.append(value)
            else:
                if attribute1['v'] not in new_list_of_values and attribute1['v']:
                    new_list_of_values.append(attribute1['v'])

            if isinstance(attribute2['v'], list):
                for value in attribute2['v']:
                    if (attribute2['d'] == 'list_of_boolean') or (value not in new_list_of_values and value):
                        new_list_of_values.append(value)
            else:
                if attribute2['v'] not in new_list_of_values and attribute2['v']:
                    new_list_of_values.append(attribute2['v'])

            if attribute1['d'] == 'list_of_boolean':
                # if new_list_of_values contains a list of booleans and they all have the same value,
                # then replace all values with one
                set_of_booleans = set(attribute1['v'])
                if len(set_of_booleans) == 1:
                    new_list_of_values = list(set_of_booleans)

            attribute1['v'] = new_list_of_values


    def _collapse_edges(self, network_in_cx):

        unique_edges = {}

        # in the loop below, we build a map where key is a tuple (edge_source, interacts, edge_target)
        # and the value is a list of edge ids
        for edge_id, edge in network_in_cx.edges.items():

            edge_key = (edge['s'], edge['i'], edge['t'])
            edge_key_reverse = (edge['t'], edge['i'], edge['s'])


            if edge_key in unique_edges:
                if (edge_id not in unique_edges[edge_key]):
                    unique_edges[edge_key].append(edge_id)

            elif edge_key_reverse in unique_edges:
                if (edge_id not in unique_edges[edge_key_reverse]):
                    unique_edges[edge_key_reverse].append(edge_id)

            else:
                unique_edges[edge_key] = [edge_id]

        #print(len(unique_edges))

        # build collapsed edges and collapsed edges attributes
        # and then use them to replace network_in_cx.edges and network_in_cx.edgeAttributes
        collapsed_edges = {}
        collapsed_edgeAttributes = {}


        # create a new edges aspect in collapsed_edges
        for key, list_of_edge_attribute_ids in unique_edges.items():
            number_of_edges = len(list_of_edge_attribute_ids)
            edge_id = list_of_edge_attribute_ids.pop(0)
            collapsed_edges[edge_id] = network_in_cx.edges[edge_id]

            if not list_of_edge_attribute_ids:
                collapsed_edgeAttributes[edge_id] = network_in_cx.edgeAttributes[edge_id]
                del network_in_cx.edgeAttributes[edge_id]
                continue

            attribute_list = network_in_cx.edgeAttributes[edge_id]

            # here, the list of collapsed edges is not empty, we need to iterate over it
            # and add attributes of the edge(s) to already existing list of edge attributes
            for attribute_id in list_of_edge_attribute_ids:

                attribute_list_for_adding = network_in_cx.edgeAttributes[attribute_id]

                self._merge_attributes(attribute_list, attribute_list_for_adding)

                if number_of_edges > 1:
                    collapse_index = {
                        'po': edge_id,
                        'n': 'Collapse Index',
                        'v': number_of_edges,
                        'd': 'long'
                    }
                    attribute_list.append(collapse_index)

                collapsed_edgeAttributes[edge_id] = attribute_list

        del network_in_cx.edges
        network_in_cx.edges = collapsed_edges

        del network_in_cx.edgeAttributes
        network_in_cx.edgeAttributes = collapsed_edgeAttributes


    def _get_network_summaries_from_NDEx_server(self):

        try:
            network_summaries = self._ndex.get_network_summaries_for_user(self._user)
        except Exception as e:
            print("\n{}: {}".format(type(e).__name__, e))
            return None, ERROR

        return network_summaries, SUCCESS


    def _get_network_uuid(self, network_name, network_summaries):

        for summary in network_summaries:
            network_name_1 = summary.get('name')

            if network_name_1 is not None:
                if network_name_1 == network_name:
                    return summary.get('externalId'), True

        return None, False


    def _network_exists_on_server(self, ptm_CX_network, summaries):

        network_name = ptm_CX_network.get_name()

        for summary in summaries:
            network_name_1 = summary.get('name')

            if network_name_1 is not None:
                if network_name_1 == network_name:
                    return summary.get('externalId')

        return None


    def _rename_ptm_network_nodes(self, ptm_network_in_cx):

        pattern = re.compile("^([A-Za-z]+[0-9]*)-([A-Z]+)-([0-9]+|[A-Za-z]+)$")

        for index, node in ptm_network_in_cx.nodes.items():
            if node['n'] and pattern.match(node['n']):
                broken_name = node['n'].split('-')
                if len(broken_name) == 3:
                    if broken_name[2].strip().lower() == 'undefined':
                        node['n'] = broken_name[1] + '?'
                    else:
                        node['n'] = broken_name[1] + broken_name[2]


    def _get_ptm_ids_for_edge(self, edge_attributes):
        for edge_attribute in edge_attributes:
            if edge_attribute['n'] and edge_attribute['n'].strip().lower() == 'biogrid ptm id':
                if edge_attribute['v']:
                    return edge_attribute['v']

        return None


    def _add_ptm_ids_to_target_node(self, biogrid_ptm_ids, node_attributes):
        for node_attribute in node_attributes:
            if node_attribute['n'] and node_attribute['n'].strip().lower() == 'biogrid ptm id':
                node_attribute['v'] = biogrid_ptm_ids
                break

    def _add_BioGRID_PTM_IDs_to_ptm_nodes(self, ptm_network_in_cx):
        for index, edge in ptm_network_in_cx.edges.items():
            edge_attributes = ptm_network_in_cx.edgeAttributes[edge['@id']]
            biogrid_ptm_ids = self._get_ptm_ids_for_edge(edge_attributes)

            if biogrid_ptm_ids:
                target_node_attributes = ptm_network_in_cx.nodeAttributes[edge['t']]
                self._add_ptm_ids_to_target_node(biogrid_ptm_ids, target_node_attributes)



    def _build_pti_node_name_to_node_id_dictionary(self, pti_CX_network):
        pti_nodes = {}
        for node in pti_CX_network.get_nodes():
            node_obj = node[1]
            node_id = node_obj['@id']
            node_name = node_obj['n']
            if node_name not in pti_nodes:
                pti_nodes[node_name] = node_id
            else:
                raise Exception('Found duplicate node name in PTI network: ' + node_name + ' ids: ' +
                                pti_nodes[node_name] + ', ' + node_id)

        return pti_nodes




    def _build_ptm_node_name_to_node_id_dictionary(self, ptm_CX_network):
        ptm_nodes = {}
        for node in ptm_CX_network.get_nodes():
            node_obj = node[1]
            node_id = node_obj['@id']
            node_name = node_obj['n']

            node_type = ptm_CX_network.get_node_attribute(node_obj, 'type')

            if node_type and node_type['v'].strip().lower() == 'protein':

                # for ptm network we only want names of protein/gene nodes
                if node_name not in ptm_nodes:
                    ptm_nodes[node_name] = node_id
                else:
                    raise Exception('Found duplicate node name in PTM network: ' + node_name + ' ids: ' +
                                    ptm_nodes[node_name] + ', ' + node_id)
        return ptm_nodes


    def _get_all_edges_for_node(self, node_id, cx_network):
        edges = []

        for edge in cx_network.get_edges():
            if node_id == edge[1]['s'] or node_id == edge[1]['t']:
                edges.append(edge[1])

        return edges






    def _merge_ptm_onto_pti(self, pti_node_name_dict, ptm_node_name_dict, pti_CX_network,
                            ptm_CX_network, protein_id_to_ptm_ids_dict, src_target_edge_ptm_ids_dict):

        pti_CX_network.node_int_id_generator = max(pti_CX_network.nodes.keys()) + 1
        pti_CX_network.edge_int_id_generator = max(pti_CX_network.edges.keys()) + 1

        inv_ptm_node_name_dict = {v: k for k, v in ptm_node_name_dict.items()}

        # iterate over ptm protein nodes
        for protein_id, ptms in protein_id_to_ptm_ids_dict.items():
            # print(f'protein_id={protein_id}  ptms={ptms}')

            pti_protein_node_id = pti_node_name_dict[inv_ptm_node_name_dict[protein_id]]

            for ptm_id in ptms:
                # get ptm node and ptm nodes' properties
                ptm_node = ptm_CX_network.get_node(ptm_id)
                ptm_node_props = ptm_CX_network.get_node_attributes(ptm_id)

                # add this ptm node to pti network
                new_node_id = pti_CX_network.create_node(ptm_node['n'], ptm_node['r'])

                # modify nodes properties to correct node id
                for prop in ptm_node_props:
                    prop['po'] = new_node_id

                # set the node attributes to the node in pti network
                pti_CX_network.nodeAttributes[new_node_id] = copy.deepcopy(ptm_node_props)


                ptm_edge_id = src_target_edge_ptm_ids_dict.get((protein_id, ptm_id), None)
                if ptm_edge_id is None:
                    raise Exception('Unable to find edge with between nodes with Ids ' + protein_id + ' and ' + ptm_id)

                # from ptm network, get edge and edge attributes
                ptm_edge = ptm_CX_network.edges[ptm_edge_id]
                ptm_edge_props = ptm_CX_network.get_edge_attributes(ptm_edge_id)

                # add this edge to pti network between protein node and newly added ptm node
                new_edge_id = pti_CX_network.create_edge(pti_protein_node_id, new_node_id, ptm_edge['i'])

                # modify edge properties to correct edge id
                for prop in ptm_edge_props:
                    prop['po'] = new_edge_id

                # set the node attributes to the node in pti network
                pti_CX_network.edgeAttributes[new_edge_id] = copy.deepcopy(ptm_edge_props)

        return pti_CX_network


    def _build_protein_id_to_ptm_ids_dict(self, protein_name_dict, ptm_CX_network):

        protein_id_to_ptm_ids_dict = {}

        inv_protein_name_dict = {v: k for k, v in protein_name_dict.items()}


        for edge in ptm_CX_network.get_edges():
            edge_source_id = edge[1]['s']
            edge_target_id = edge[1]['t']

            if edge_source_id in inv_protein_name_dict:
                if edge_source_id not in protein_id_to_ptm_ids_dict:
                    protein_id_to_ptm_ids_dict[edge_source_id] = []

                protein_id_to_ptm_ids_dict[edge_source_id].append(edge_target_id)

            #elif edge_target_id in inv_protein_name_dict:
            #    if edge_target_id not in protein_id_to_ptm_ids_dict:
            #        protein_id_to_ptm_ids_dict[edge_target_id] = []

            #    protein_id_to_ptm_ids_dict[edge_source_id].append(edge_target_id)

        return protein_id_to_ptm_ids_dict


    def _build_src_target_edge_ptm_ids_dict(self, ptm_CX_network):
        src_target_edge_ptm_ids_dict = {}

        for edge_tuple in ptm_CX_network.get_edges():
            edge_id = edge_tuple[0]
            edge_source_node_id = edge_tuple[1]['s']
            edge_target_node_id = edge_tuple[1]['t']

            key = (edge_source_node_id, edge_target_node_id)
            if key not in src_target_edge_ptm_ids_dict:
                src_target_edge_ptm_ids_dict[key] = edge_id

        return src_target_edge_ptm_ids_dict


    def run(self):
        """
        Runs content loading for NDEx KINOME Content Loader
        :param theargs:
        :return:
        """
        self._parse_config()
        self._load_style_template()

        data_dir_existed = self._check_if_data_dir_exists()

        if self._skipdownload is False or data_dir_existed is False:
            status_code = self._download_kinome_files()
            if status_code != 0:
                return ERROR

            status_code = self._unzip_kinome()
            if status_code != 0:
                return ERROR

        self._build_gene_lookup()


        self._create_ndex_connection()


        summaries, ret_value = self._get_network_summaries_from_NDEx_server()
        if ret_value != SUCCESS:
            return ret_value


        # Step 1 - create PPI file from GENES and INTERACTIONS files
        self._create_ppi_file()

        pti_CX_network, ret_value = self._generate_CX_file(self._pti_load_plan, self._ppi_network_1)
        if ret_value != SUCCESS:
            return ret_value

        self._collapse_edges(pti_CX_network)
        self._init_network_attributes(pti_CX_network, 'pti')
        self._write_nice_cx_to_file(pti_CX_network, self._cx_pti)
        network_UUID = self._network_exists_on_server(pti_CX_network, summaries)
        self._upload_CX(self._cx_pti, network_UUID)


        # Step 2 - create PTM network file
        self._create_ptm_file()

        ptm_CX_network, ret_value = self._generate_CX_file(self._ptm_load_plan, self._ptm_network_2)
        self._rename_ptm_network_nodes(ptm_CX_network)
        if ret_value != SUCCESS:
            return ret_value

        self._collapse_edges(ptm_CX_network)
        self._add_BioGRID_PTM_IDs_to_ptm_nodes(ptm_CX_network)
        self._init_network_attributes(ptm_CX_network, 'ptm')
        self._write_nice_cx_to_file(ptm_CX_network, self._cx_ptm)
        network_UUID = self._network_exists_on_server(ptm_CX_network, summaries)
        self._upload_CX(self._cx_ptm, network_UUID)


        # Step 3 - merge PTM network with PTI network on protein/genes:
        # in essence, we add edges from PTM network to PTI based on node names
        pti_CX_network = ndex2.create_nice_cx_from_file(self._cx_pti)
        ptm_CX_network = ndex2.create_nice_cx_from_file(self._cx_ptm)

        # in this dictionary for pti network, key is protein node name, value is to node id:
        #   pti_node_name_dict: { 'CHD1': 0, 'CKA1': 1, 'CKA2': 2, ...}
        pti_node_name_dict = self._build_pti_node_name_to_node_id_dictionary(pti_CX_network)

        # in this dictionary for ptm network, key is protein node name, value is to node id:
        #   pti_node_name_dict: { 'ADK1': 0, 'ADR1': 2, 'AKL1': 40, ...}
        ptm_node_name_dict = self._build_ptm_node_name_to_node_id_dictionary(ptm_CX_network)

        #print(f'protein nodes in pti={len(pti_node_name_dict)}   protein nodes in ptm={len(ptm_node_name_dict)}')


        # in this dictionary for ptm network, key is protein node id, value is list of ptm ids:
        #   pti_node_name_dict: {0: [1, 3108, 3521, 3522, 3523], 2: [3, 4, 5, 6, 7, 8, 9], 40: [41, 42, 43, 44, 45], ...}
        protein_id_to_ptm_ids_dict = self._build_protein_id_to_ptm_ids_dict(ptm_node_name_dict, ptm_CX_network)

        # in this dictionary for ptm network, key is a tuple (source Id, target Id), and
        # value is edge id
        src_target_edge_ptm_ids_dict = self._build_src_target_edge_ptm_ids_dict(ptm_CX_network)


        merged_ptm_pti_network = self._merge_ptm_onto_pti(pti_node_name_dict, ptm_node_name_dict,
                  pti_CX_network, ptm_CX_network, protein_id_to_ptm_ids_dict, src_target_edge_ptm_ids_dict)

        self._init_network_attributes(merged_ptm_pti_network, 'merged')
        self._write_nice_cx_to_file(merged_ptm_pti_network, self._cx_merged)
        network_UUID = self._network_exists_on_server(merged_ptm_pti_network, summaries)
        self._upload_CX(self._cx_merged, network_UUID)


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
        return ERROR
    finally:
        logging.shutdown()


if __name__ == '__main__':  # pragma: no cover
    sys.exit(main(sys.argv))
