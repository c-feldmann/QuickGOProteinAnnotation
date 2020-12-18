import json
import warnings
from typing import *
import pandas as pd

import networkx as nx
import requests


def json_query(url):
    """Loads data from URL and returns them in json format."""
    r = requests.get(url, headers={"Accept": "application/json"})
    data = json.loads(r.text)
    if not r.ok:
        if r.status_code == 400:
            raise ValueError(r.status_code, f"Invalid URL: {url}")
        raise ConnectionError(r.status_code)
    return data


class GOTerm:
    """ This class represents a term used to annotate entries in EBI QuickGO. [1]

    QuickGO has hierarchically arranged terms for annotation. Each term has an ID, a name and usually also a short
    description.
    This class is designed for comfortable access to these data. Since the IDs may change, but outdated IDs map to the
    new assignment, the class aims to correct outdated IDs.

    Example Terms;
        "kinase activity": "GO:0016301"
        "cellular_component": "GO:0005575"

    Example Usage.
        kinase_activity = GOTerm("GO:0016301")
        print(kinase_activity.name)
            -> "kinase activity"

    Reference:
        [1] https://www.ebi.ac.uk/QuickGO/
    """

    def __init__(self, go_id: str, name: str, definition: str, aspect: str):
        """Initializes the object.

        Args:
            go_id (str): GO ID which is represented by this object.
        """
        self._go_id = go_id
        self._name = name
        self._definition = definition
        self._aspect = aspect

    def __str__(self):
        return self.name

    @classmethod
    def from_id(cls, go_id):
        """Queries QuickGO Webserver for Data"""
        if go_id[:3] != "GO:":
            raise ValueError('Every ID begins with: "GO:"')

        try:
            base_url = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms"
            incoming_data = json_query("{}/Go%3A{}".format(base_url, go_id.lstrip("GO:")))

            if len(incoming_data['results']) == 0:
                raise ValueError("Received no ID for {}".format(go_id))
            elif len(incoming_data['results']) != 1:
                raise ValueError("Received multiple possible IDs for {}".format(go_id))
            data = incoming_data['results'][0]

        except ValueError:
            warnings.warn(f"Not a valid ID: {go_id}")
            raise

        if data['isObsolete'] is True:
            raise AssertionError(f'{go_id} is obsolete!')

        data_go_id = data['id']
        if data_go_id != go_id:
            warnings.warn('{} is updated to {}!'.format(go_id, data_go_id))
        name = data['name']
        definition = data['definition']['text']
        aspect = data['aspect']
        return cls(data_go_id, name=name, definition=definition, aspect=aspect)

    @property
    def go_id(self) -> str:
        return self._go_id

    @property
    def name(self) -> str:
        return self._name

    @property
    def definition(self) -> str:
        return self._definition

    @property
    def aspect(self) -> str:
        return self._aspect


class GoMolecularFunction(GOTerm):
    """ "Molecular Functions" are a subset of GOTerms used to describe functions of a protein.

    Wrapper-class for GOTerm to ensure the aspect is "molecular_function"
    """

    def __init__(self, go_id: str, name: str, definition: str, aspect: str):
        """ Initializes the Object and ensures the correct aspect association.

        Args:
            go_id: ID of GOTerm represented by the object.

        Raises:
            TypeError: If the downloaded GOTerm has a aspect other than "molecular_function".
                Other aspects:  "Cellular Component" and "Biological Process".

        """
        super().__init__(go_id, name, definition, aspect)
        if self.aspect != 'molecular_function':
            raise TypeError(f'{go_id} does not refer to a molecular function!')


class GOMolecularFunctionHierarchy:
    """Class to manage hierarchically relations of GOMolecularFunctions.

     GOMolecularFunctions are arranged hierarchically and for a protein with a specific function, all supercategories of
     this function are present in the protein. E.g. A protein kinase is always a kinase, a transferase,
     and always an enzyme. To manage these dependencies a directed graphs is utilised.

    """

    def __init__(self):
        self._rerouted_function_dict: Dict[str, str] = dict()
        self._function_dict: Dict[str, GoMolecularFunction] = dict()
        self._graph = nx.DiGraph()

    @property
    def managed_function_ids(self) -> Set[str]:
        assert not set(self._rerouted_function_dict.keys()).intersection(self._function_dict.keys())
        return set(self._rerouted_function_dict.keys()).union(self._function_dict.keys())

    @property
    def graph(self):
        return self._graph

    def get_parent(self, go_id) -> List[GoMolecularFunction]:
        return [self._function_dict[node] for node in self.graph.neighbors(go_id)]

    def get_children(self, go_id) -> List[GoMolecularFunction]:
        return [self._function_dict[node] for node in self.graph.reverse().neighbors(go_id)]

    @staticmethod
    def _go_path(start_go_id, end_go_id):
        start_go_id = start_go_id.lstrip('GO:')
        end_go_id = end_go_id.lstrip('GO:')

        data = json_query(f"https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/GO%3A{start_go_id}/paths/"
                          f"GO%3A{end_go_id}?relations=is_a")
        if data['pageInfo'] is not None:
            raise NotImplementedError('This path has page info! There are probably multiple pages. Cannot be handled'
                                      'yet and action is aborted!')
        if not isinstance(data['results'], list):
            raise TypeError(f"Unexpected Type: {data['results']}, GO:{start_go_id}-GO:{end_go_id}")
        return data['results']

    def _add_path_to_top_level(self, go_id):
        top_level_id = "GO:0003674"  # GO ID of "Molecular Function"
        paths = self._go_path(go_id, top_level_id)
        if len(paths) == 0:
            warnings.warn(f'No path to "Molecular Function" {go_id} remains unconnected')
        for path in paths:
            for edge in path:

                # Child and parents are not added directly because the IDs might have changed. The Database is not
                # exactly up to date.
                if edge['child'] not in self.managed_function_ids:
                    self._add_go_function_without_path(edge['child'])
                child = self.get_go_function_from_id(edge['child']).go_id

                if edge['parent'] not in self.managed_function_ids:
                    self._add_go_function_without_path(edge['parent'])
                parent = self.get_go_function_from_id(edge['parent']).go_id

                if edge['relationship'] != 'is_a':
                    raise AssertionError("Unexpected relation: {}".format(edge['relationship']))

                self._graph.add_edge(child, parent)

    def _add_go_function_without_path(self, go_id):
        if go_id not in self.managed_function_ids:
            go_function = GoMolecularFunction.from_id(go_id)
            if go_id != go_function.go_id:
                self._rerouted_function_dict[go_id] = go_function.go_id
            if go_function.go_id not in self._function_dict:
                self._function_dict[go_function.go_id] = go_function

    def add_go_function(self, go_id: str) -> None:
        if go_id not in self.managed_function_ids:
            self._add_go_function_without_path(go_id)
            updated_id = self.get_go_function_from_id(go_id).go_id
            self._add_path_to_top_level(updated_id)

    def get_go_function_from_id(self, go_id) -> GoMolecularFunction:
        if go_id in self._rerouted_function_dict:
            query_id = self._rerouted_function_dict[go_id]
        else:
            query_id = go_id
        return self._function_dict[query_id]

    def get_nodes_upwards(self, go_id) -> Set[GoMolecularFunction]:
        r_nodes: Set[GoMolecularFunction] = set()
        processed_ids = set()
        iter_nodes: List[GoMolecularFunction] = [self.get_go_function_from_id(go_id)]
        while iter_nodes:
            node = iter_nodes.pop()
            r_nodes.add(node)
            processed_ids.add(node.go_id)
            for parent in self.get_parent(node.go_id):
                if parent.go_id in processed_ids:
                    continue
                iter_nodes.append(parent)
        return r_nodes

    def get_nodes_downwards(self, go_id) -> Set[GoMolecularFunction]:
        r_nodes: Set[GoMolecularFunction] = set()
        processed_ids = set()
        iter_nodes: List[GoMolecularFunction] = [self.get_go_function_from_id(go_id)]
        while iter_nodes:
            node = iter_nodes.pop()
            r_nodes.add(node)
            processed_ids.add(node.go_id)
            for child in self.get_children(node.go_id):
                if child.go_id in processed_ids:
                    continue
                iter_nodes.append(child)
        return r_nodes


class AllFunctionAnnotation:
    def __init__(self, alternative_name_dict=None, simplify_name=True):
        if alternative_name_dict is None:
            alternative_name_dict = dict()
        self._alternative_name_dict = alternative_name_dict
        self._simplify_name = simplify_name
        self._function_relations = GOMolecularFunctionHierarchy()

    def annotate_proteins(self, protein_list, as_dataframe=True) -> Union[List[Dict], pd.DataFrame]:
        result_df = []
        for uniprot_id in protein_list:
            result_df.extend(self.get_protein_functions(uniprot_id, as_dataframe=False))
        if as_dataframe:
            return pd.DataFrame(result_df)
        else:
            return result_df

    def get_protein_functions(self, uniprot_id, as_dataframe=True) -> Union[List[Dict], pd.DataFrame]:
        explicit_ids = self._get_protein_explicit_function_ids(uniprot_id)
        all_ids: Set[str] = set()
        for go_id in explicit_ids:
            all_ids.add(go_id)
            implicit_go_functions = self._function_relations.get_nodes_upwards(go_id)
            implicit_go_ids = {go_function.go_id for go_function in implicit_go_functions}
            all_ids.update(implicit_go_ids)
        # Remove "molecular_function" annotation
        all_ids.remove("GO:0003674")
        function_dict_list = []
        for go_id in all_ids:
            if go_id in self._alternative_name_dict:
                name = self._alternative_name_dict[go_id]
            elif self._simplify_name:
                name = self._function_relations.get_go_function_from_id(go_id).name
                name = name.rstrip(" activity")
            else:
                name = self._function_relations.get_go_function_from_id(go_id).name
            function_dict_list.append({"uniprot_id": uniprot_id,
                                       "go_id": go_id,
                                       "protein_function": name})
        if len(function_dict_list) == 0:
            function_dict_list.append({"uniprot_id": uniprot_id,
                                       "protein_function": "no_function"})
        if as_dataframe:
            return pd.DataFrame(function_dict_list)
        else:
            return function_dict_list

    def _get_protein_explicit_function_ids(self, uniprot_id: str):
        query_result = self._raw_function_annotations(uniprot_id)
        # Validate query
        for item in query_result:
            # Assert that query worked and returned proper go functions
            assert uniprot_id in item["geneProductId"], print(uniprot_id, item)
            assert item["goAspect"] == "molecular_function"
            assert item["qualifier"] == "enables"
        extracted_ids = [item["goId"] for item in query_result]

        # Add go ids and find potential ID updated
        for go_id in extracted_ids:
            self._function_relations.add_go_function(go_id)

        updated_ids = {self._function_relations.get_go_function_from_id(go_id).go_id for go_id in extracted_ids}
        return updated_ids

    def _raw_function_annotations(self, uniprot_id, page=1) -> list:
        """Currently all annotations are downloaded. Potentially it might be better to select a confidence level."""
        base_url = "https://www.ebi.ac.uk/QuickGO/services/annotation/"
        fixed_query = "search?selectedFields=geneProductId&"
        variable_query = f"geneProductId={uniprot_id}&aspect=molecular_function&qualifier=enables&page={page}"

        data = json_query(f"{base_url}{fixed_query}{variable_query}")

        # If no data: return Empty list
        if data["numberOfHits"] == 0:
            return []

        assert data["pageInfo"]["current"] == page
        results = data['results']
        # If not last page: get results of lower page(s) and add them to current results
        if page < data["pageInfo"]["total"]:
            lower_page_results = self._raw_function_annotations(uniprot_id, page + 1)
            results.extend(lower_page_results)
        elif page == 1:
            assert len(results) == data["numberOfHits"]
        return results


class SelectedFunctionAnnotation(AllFunctionAnnotation):
    def __init__(self, functions: set, alternative_name_dict: Optional[dict] = None, simplify_name=True):
        super(SelectedFunctionAnnotation, self).__init__(alternative_name_dict=alternative_name_dict,
                                                         simplify_name=simplify_name)
        self._functions = functions

    def get_protein_functions(self, uniprot_id, as_dataframe=True) -> Union[List[Dict], pd.DataFrame]:
        protein_function_list: List[Dict] = super().get_protein_functions(uniprot_id, as_dataframe=False)
        protein_function_list = [annotation for annotation in protein_function_list
                                 if annotation["go_id"] in self._functions]
        if len(protein_function_list) == 0:
            protein_function_list.append({"uniprot_id": uniprot_id,
                                          "protein_function": "no_function"})
        if as_dataframe:
            return pd.DataFrame(protein_function_list)
        else:
            return protein_function_list


class SpecialFunctionAnnotation(AllFunctionAnnotation):
    def __init__(self, function_specifications: List[Tuple[Set, Set, str]]):
        self._function_specifications = function_specifications
        super(SpecialFunctionAnnotation, self).__init__()

    def get_protein_functions(self, uniprot_id, as_dataframe=True) -> Union[List[Dict], pd.DataFrame]:
        protein_function_list: List[Dict] = super().get_protein_functions(uniprot_id, as_dataframe=False)
        protein_go_id_set = {annotation["go_id"] for annotation in protein_function_list}
        result_dict_list = []
        for specific_function in self._function_specifications:
            has_required = specific_function[0].issubset(protein_go_id_set)
            no_excluded = len(specific_function[1].intersection(protein_go_id_set)) == 0
            if has_required and no_excluded:
                result_dict_list.append({"uniprot_id": uniprot_id,
                                         "protein_function": specific_function[2]})
        if len(result_dict_list) == 0:
            result_dict_list.append({"uniprot_id": uniprot_id,
                                     "protein_function": "no_function"})
        if as_dataframe:
            return pd.DataFrame(result_dict_list)
        else:
            return result_dict_list
