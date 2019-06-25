#!/usr/bin/env python

from abc import ABC, abstractmethod
import mygene
import json
from Bio import Entrez


class Gene(ABC):
    """
    The Subject interface declares common operations for both RealSubject and
    the Proxy. As long as the client works with RealSubject using this
    interface, you'll be able to pass it a proxy instead of a real subject.
    """

    @abstractmethod
    def map_ids(self, symbol):
        pass

    def _aggregate_ids_in_single_format(self, gene_ids_info):
        pass

    def _search_entrez(self, gene_schema):
        pass

class RealGene(Gene):
    """
    The RealSubject contains some core business logic. Usually, RealSubjects are
    capable of doing some useful work which may also be very slow or sensitive -
    e.g. correcting input data. A Proxy can solve these issues without any
    changes to the RealSubject's code.
    """
    def __init__(self):
        Entrez.email = 'ferrando.rosa@gmail.com'
        self._mg = mygene.MyGeneInfo()

    def map_ids(self, symbol):
        # print("RealGene: Handling request.")
        gene_info_schema = {}
        try:
            gene_info = self._mg.query(symbol, scopes='entrezgene', fields='symbol, entrezgene', species='human', verbose=False)

            if (gene_info['total'] > 0):
                gene_info_schema = self._aggregate_ids_in_single_format(gene_info)
            else:
                gene_info_schema['symbol'] = symbol

            gene_info_schema = self._search_entrez(gene_info_schema)
            # print("Gene Schema: ", gene_info_schema)
        except Exception as e:
            gene_info_schema['symbol'] = symbol
            print(symbol, e)

        return gene_info_schema

    def _aggregate_ids_in_single_format(self, gene_ids_info):
        gene_schema = {}

        json_with_entrezgene = next(((i,d)  for i,d in enumerate(gene_ids_info['hits']) if 'entrezgene' in d), None)
        if not(json_with_entrezgene == None):
            gene_schema['symbol'] = gene_ids_info['hits'][json_with_entrezgene[0]]['symbol']
            gene_schema['entrezgene'] = gene_ids_info['hits'][json_with_entrezgene[0]]['entrezgene']
            if (gene_ids_info['total'] == 2):
                json_index = abs(json_with_entrezgene[0]-1)
            else:
                json_index = json_with_entrezgene[0]

            gene_schema['ensemble'] = gene_ids_info['hits'][json_index]['_id']
        else:
            gene_schema['symbol'] = gene_ids_info['hits'][0]['symbol']
            gene_schema['ensemble'] = gene_ids_info['hits'][0]['_id']

        return gene_schema

    def _search_entrez(self, gene_schema):
        if not any('entrezgene' in d for d in gene_schema) and any('ensemble' in d for d in gene_schema):
            gene_name= gene_schema['ensemble']
            most_likely_entry = Entrez.esearch(db="gene",term="{gene_name} AND 9606 [Taxonomy ID]".format(gene_name=gene_name),retmode="json")
            most_likely_entry_json = json.loads(most_likely_entry.read())
            my_ids = most_likely_entry_json['esearchresult']['idlist']
            if (my_ids != []):
                gene_schema['entrezgene'] = my_ids[0]
        return gene_schema

class Proxy(Gene):
    """
    The Proxy has an interface identical to the RealSubject.
    """

    def __init__(self):
        self._real_gene = RealGene()
        self._gene_schema = []


    def map_ids(self,symbol):
        """
        The most common applications of the Proxy pattern are lazy loading,
        caching, controlling the access, logging, etc. A Proxy can perform one
        of these things and then, depending on the result, pass the execution to
        the same method in a linked RealSubject object.
        """
        info_cached = self.check_cached(symbol)
        if not(info_cached['cached']):
            gene_info_schema = self._real_gene.map_ids(symbol)
            self.write_cache(gene_info_schema)
        else:
            gene_info_schema = info_cached['gene_info']
#        print(gene_info_schema)
        return gene_info_schema


    def _get_cached_gene_by_symbol(self, symbol):
        return next((d for i,d in enumerate(self._gene_schema) if d['symbol'] == symbol), None)

    def check_cached(self, symbol):
        # print("Proxy: Checking gene cached prior to firing a real request: ")
        gene_info = self._get_cached_gene_by_symbol(symbol)
        if not(gene_info == None):
            return {'cached': True, 'gene_info': gene_info}
        else:
            return {'cached': False, 'gene_info': None}

    def write_cache(self, gene_info_schema):
        # print("Proxy: writing gene in cache: ", gene_info_schema)
        self._gene_schema.append(gene_info_schema)



def client_code(gene: Gene, symbol_vector):
    """
    The client code is supposed to work with all objects (both subjects and
    proxies) via the Subject interface in order to support both real subjects
    and proxies. In real life, however, clients mostly work with their real
    subjects directly. In this case, to implement the pattern more easily, you
    can extend your proxy from the real subject's class.
    """
    # ...
    all_symbol_info = list(map(gene.map_ids, symbol_vector))
    print(all_symbol_info)
    # ...


if __name__ == "__main__":
    print("Client: Executing the same client code with a proxy:")
    proxy = Proxy()

    my_genes = ['AC005165.1','HNRNPA3P1', 'LOC646588', 'LOC100289113', 'AC005165.1','HNRNPA3P1']
    client_code(proxy, my_genes)
