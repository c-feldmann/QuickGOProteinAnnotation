import json
import warnings
from typing import *

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

    Example Useage.
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

    Wrapperclass for GOTerm to ensure the aspect is "molecular_function"
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

     GOMolecularFunctions are arranged hierarchically and for a protein with a specific function, all supercategorys of
     this function are present in the protein. E.g. A protein kinase is always a kinase, a transferase,
     and always an enzyme. To manage these dependencies a directed graphs is utilised.

    """

    def __init__(self):
        self._rerouted_function_dict: Dict[str, str] = dict()
        self._function_dict: Dict[str, GoMolecularFunction] = dict()
        self._graph = nx.DiGraph()

    @property
    def mangaged_funcition_ids(self) -> Set[str]:
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
                if edge['child'] not in self.mangaged_funcition_ids:
                    self._add_go_function_from_id(edge['child'])
                child = self.get_go_function_from_id(edge['child']).go_id

                if edge['parent'] not in self.mangaged_funcition_ids:
                    self._add_go_function_from_id(edge['parent'])
                parent = self.get_go_function_from_id(edge['parent']).go_id

                if edge['relationship'] != 'is_a':
                    raise AssertionError("Unexpected relation: {}".format(edge['relationship']))

                self._graph.add_edge(child, parent)

    def _add_go_function_from_id(self, go_id):
        if go_id not in self.mangaged_funcition_ids:
            go_function = GoMolecularFunction.from_id(go_id)
            if go_id != go_function.go_id:
                self._rerouted_function_dict[go_id] = go_function.go_id
            if go_function.go_id not in self._function_dict:
                self._function_dict[go_function.go_id] = go_function

    def add_go_function(self, go_id: str) -> None:
        if go_id not in self.mangaged_funcition_ids:
            self._add_go_function_from_id(go_id)
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


class Protein:
    def __init__(self, uniprot_id: str):
        self._uniprot_id = uniprot_id
        self._function_id_set = set()

    @property
    def uniprot_id(self):
        return self._uniprot_id

    @property
    def functions(self) -> Set[str]:
        return self._function_id_set

    @functions.setter
    def functions(self, function_dict):
        if not isinstance(function_dict, set):
            raise TypeError
        self._function_id_set = function_dict


class ProteinFactory:
    def __init__(self):
        self._hierarchy = GOMolecularFunctionHierarchy()
        self._created_proteins: Dict[str, Protein] = dict()
        self._managed_classes_dict: Dict[str, Dict[str, Set[str]]] = dict()
        self._protein_class_dict = dict()
        self._updated = False

    def _create_protein(self, uniprot_id) -> Protein:
        protein = Protein(uniprot_id)
        protein_functions = list(self.go_protein_function_annotations(uniprot_id))

        r_ids = set()
        while protein_functions:
            iter_id = protein_functions.pop()
            if iter_id in r_ids:
                continue
            r_ids.add(iter_id)
            upward_nodes = [node.go_id for node in self._hierarchy.get_nodes_upwards(iter_id)]
            r_ids.update(upward_nodes)

        protein.functions = r_ids
        return protein

    def _query_go_for_function(self, uniprot_id, page=1) -> list:
        """Currently all annotations are downloaded. Potentially it might be better to select a confidence level."""
        data = json_query(f"https://www.ebi.ac.uk/QuickGO/services/annotation/search?selectedFields=geneProductId&"
                          f"geneProductId={uniprot_id}&aspect=molecular_function&qualifier=enables&page={page}")

        if data["numberOfHits"] == 0:
            return []
        assert data["pageInfo"]["current"] == page
        if data["pageInfo"]["total"] == 1:
            return data['results']
        else:
            if page < data["pageInfo"]["total"]:
                # If not last page: get content of page and add content of next page, this works recursive.
                cur_data = data['results']
                next_data = self._query_go_for_function(uniprot_id, page + 1)
                cur_data.extend(next_data)
                if page == 1:
                    assert len(cur_data) == data["numberOfHits"]
                return cur_data
            else:
                # If last Page return content
                return data['results']

    def go_protein_function_annotations(self, uniprot_id):
        query_result = self._query_go_for_function(uniprot_id)
        for entry in query_result:
            assert uniprot_id in entry["geneProductId"], print(uniprot_id, entry)
            if entry["goAspect"] == "molecular_function":
                if entry["qualifier"] == "enables":
                    self._hierarchy.add_go_function(entry["goId"])
                    yield self._hierarchy.get_go_function_from_id(entry["goId"]).go_id

    def add_protein(self, uniprot_id) -> None:
        protein = self._create_protein(uniprot_id)
        self._created_proteins[uniprot_id] = protein

    def get_go_function(self, go_id) -> GoMolecularFunction:
        return self._hierarchy.get_go_function_from_id(go_id)

    def define_protein_class(self,
                             name: str,
                             present_functions: Set[str],
                             not_present_functions: Set[str],
                             check_definition_clash: bool = True):
        """Adds a protein_class for classification.

        TODO: Add explanation

        Args:
            check_definition_clash (bool):
            name:
            present_functions (set): All GO IDs which must be present
            not_present_functions (set): All GO IDs which must not be matched

        Returns:
            None
        """
        for go_id in present_functions.union(not_present_functions):
            self._hierarchy.add_go_function(go_id)
        self._updated = False

        if check_definition_clash:
            new_definion_covered_go_terms = set()
            for go_id in present_functions:
                new_definion_covered_go_terms.update([node.go_id for node in
                                                      self._hierarchy.get_nodes_downwards(go_id)])
            # Remove and ignore functions covered by "not_present_functions"
            for go_id in not_present_functions:
                new_definion_covered_go_terms.difference_update([node.go_id for node in
                                                                 self._hierarchy.get_nodes_downwards(go_id)])

            for protein_class, definition in self._managed_classes_dict.items():
                covered_go_terms = set()
                for go_id in definition["contains"]:
                    covered_go_terms.update([node.go_id for node in self._hierarchy.get_nodes_downwards(go_id)])
                for go_id in definition["contains_not"]:
                    covered_go_terms.difference_update([node.go_id for node in
                                                        self._hierarchy.get_nodes_downwards(go_id)])
                inter_class_go_intersecitons = covered_go_terms.intersection(new_definion_covered_go_terms)
                if len(inter_class_go_intersecitons) > 0:
                    raise ValueError(f"{inter_class_go_intersecitons} are already used in {protein_class}!")

        self._managed_classes_dict[name] = {"contains": present_functions,
                                            "contains_not": not_present_functions}

    def _functions_of_class(self, protein_class):
        class_functions = set()
        for go_id in self._managed_classes_dict[protein_class]["contains"]:
            class_functions.update([node.go_id for node in self._hierarchy.get_nodes_downwards(go_id)])
        for go_id in self._managed_classes_dict[protein_class]["contains_not"]:
            class_functions.difference_update([node.go_id for node in self._hierarchy.get_nodes_downwards(go_id)])
        return class_functions

    def check_classification_overlap(self):
        function_dict = dict()
        for p_class in self._managed_classes_dict.keys():
            function_dict[p_class] = self._functions_of_class(p_class)

        for i, p_class1 in enumerate(self._managed_classes_dict.keys()):
            p_class1_functions = function_dict[p_class1]

            for j, p_class2 in enumerate(self._managed_classes_dict.keys()):
                if i >= j:
                    continue
                p_class2_functions = function_dict[p_class2]
                overlap = p_class1_functions.intersection(p_class2_functions)
                if len(overlap) > 0:
                    warnings.warn(f"{overlap} occure in {p_class1} and {p_class2}")

    def _is_match(self, uniprot_id, class_name):
        must_occure = self._managed_classes_dict[class_name]["contains"]
        must_not_occure = self._managed_classes_dict[class_name]["contains_not"]
        protein_functions = set(self._created_proteins[uniprot_id].functions)
        if must_occure.issubset(protein_functions) and len(must_not_occure.intersection(protein_functions)) == 0:
            return True
        else:
            return False

    def _classify_protein(self, uniprot_id) -> Set[str]:
        matching_classes = set()
        for protein_class in self._managed_classes_dict.keys():
            if self._is_match(uniprot_id, protein_class):
                matching_classes.add(protein_class)
        return matching_classes

    def get_matching_classes(self, uniprot_id) -> Set[str]:
        if uniprot_id not in self._created_proteins:
            self.add_protein(uniprot_id)
        return self._classify_protein(uniprot_id)

    def get_protein_class(self, uniprot_id) -> str:
        if uniprot_id not in self._created_proteins:
            self.add_protein(uniprot_id)
        matching_classes = self._classify_protein(uniprot_id)
        if len(matching_classes) == 1:
            return list(matching_classes)[0]
        else:
            return "Multiclass_({})".format(", ".join(sorted(matching_classes)))
