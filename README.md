# QuckGOProteinAnnotation
The database [QuickGO](https://www.ebi.ac.uk/QuickGO/) provides protein function annotations for proteins, specified by
[UniProt](https://www.uniprot.org/) ID. Arranging proteins by function rather than family extends protein associations
beyond evolutionary relations. However, proteins may have multiple functions (e.g. receptor tyrosine kinases) and
are therefore not uniquely assigned.

Provided code can be used to extract (specified) annotations from QuickGO.

## Installation in Conda
If not already installed, install **pip** and **git**:
```
conda install git
conda install pip
```
Then install via pip:
```
pip install git+git://github.com/c-feldmann/QuickGOProteinAnnotation
```

## Quickstart
### From Terminal
The script `annotate_protein_list.py` takes an input-file (here:  *demo_data/demo_uniprot_ids.tsv*) where
proteins are specified in the column "uniprot_id". Results are saved to the file *demo_data/demo_output.tsv* as a
tab-separated file.
```
python annotate_protein_list.py -i demo_data/demo_uniprot_ids.tsv -o demo_data/demo_output.tsv -c "uniprot_id" -s tab
```
|   Argument    |   Explanation |
|:-------------:|:-------------:|
| -i | input file |
| -o | output file |
| -c | column name |
| -s | separator |

The default value for `-s` is "tab", whereas the default output-file is named *go_function_annotation.tsv*.
### In Python
A short example how this package could be used in a python code:


```python
from go_protein_annotation  import DefaultAnnotation
```


```python
test_proteins = ["Q16512", "P30085", "P25774"]
default_annotation = DefaultAnnotation()
protein_class_df = default_annotation.annotate_proteins(test_proteins)
```


```python
protein_class_df
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>uniprot_id</th>
      <th>protein_function</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>Q16512</td>
      <td>Transcription regulator</td>
    </tr>
    <tr>
      <th>1</th>
      <td>Q16512</td>
      <td>Kinase</td>
    </tr>
    <tr>
      <th>2</th>
      <td>P30085</td>
      <td>Kinase</td>
    </tr>
    <tr>
      <th>3</th>
      <td>P25774</td>
      <td>Peptidase</td>
    </tr>
  </tbody>
</table>
</div>



## Details
QuckGO functions are ordered hierarchically. E.g. an explicit annotation of
[peptidase activity](https://www.ebi.ac.uk/QuickGO/term/GO:0008233) implies also a
[hydrolase activity](https://www.ebi.ac.uk/QuickGO/term/GO:0016787). Provided code extracts
all explicit functional annotations and extends it with implicit annotations.
### All Protein Functions
To obtain *all* annotations for a protein the class `AllFunctionAnnotation` is used.


```python
from go_protein_annotation import  AllFunctionAnnotation
all_functions = AllFunctionAnnotation()

# For a single protein
all_functions_q16512 = all_functions.get_protein_functions("Q16512")

# For a list of proteins
protein_functions = all_functions.annotate_proteins(["Q16512", "P30085"])
```


```python
all_functions_q16512.head(10)
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>uniprot_id</th>
      <th>go_id</th>
      <th>protein_function</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>Q16512</td>
      <td>GO:0005515</td>
      <td>protein binding</td>
    </tr>
    <tr>
      <th>1</th>
      <td>Q16512</td>
      <td>GO:0035639</td>
      <td>purine ribonucleoside triphosphate binding</td>
    </tr>
    <tr>
      <th>2</th>
      <td>Q16512</td>
      <td>GO:0000166</td>
      <td>nucleotide binding</td>
    </tr>
    <tr>
      <th>3</th>
      <td>Q16512</td>
      <td>GO:1901363</td>
      <td>heterocyclic compound binding</td>
    </tr>
    <tr>
      <th>4</th>
      <td>Q16512</td>
      <td>GO:0050681</td>
      <td>androgen receptor binding</td>
    </tr>
    <tr>
      <th>5</th>
      <td>Q16512</td>
      <td>GO:0140110</td>
      <td>transcription regulator</td>
    </tr>
    <tr>
      <th>6</th>
      <td>Q16512</td>
      <td>GO:0017076</td>
      <td>purine nucleotide binding</td>
    </tr>
    <tr>
      <th>7</th>
      <td>Q16512</td>
      <td>GO:0019901</td>
      <td>protein kinase binding</td>
    </tr>
    <tr>
      <th>8</th>
      <td>Q16512</td>
      <td>GO:0042826</td>
      <td>histone deacetylase binding</td>
    </tr>
    <tr>
      <th>9</th>
      <td>Q16512</td>
      <td>GO:0035257</td>
      <td>nuclear hormone receptor binding</td>
    </tr>
  </tbody>
</table>
</div>




```python
protein_functions.groupby("uniprot_id").nunique()
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>go_id</th>
      <th>protein_function</th>
    </tr>
    <tr>
      <th>uniprot_id</th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>P30085</th>
      <td>30</td>
      <td>30</td>
    </tr>
    <tr>
      <th>Q16512</th>
      <td>55</td>
      <td>55</td>
    </tr>
  </tbody>
</table>
</div>



### A Subset of Protein Functions
Often it can be useful to extract only a subset of protein functions. This can be achieved using the class
`SelectedFunctionAnnotation`.


```python
from go_protein_annotation import SelectedFunctionAnnotation
```


```python
selected_functions = {"GO:0016301",  # Kinase activity
                      "GO:0140110",  # Transcription regulator activity
                      "GO:0008233",  # Peptidase activity
                      }
sel_function_extraction = SelectedFunctionAnnotation(selected_functions)
out = sel_function_extraction.get_protein_functions("Q16512")
```


```python
out
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>uniprot_id</th>
      <th>go_id</th>
      <th>protein_function</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>Q16512</td>
      <td>GO:0140110</td>
      <td>transcription regulator</td>
    </tr>
    <tr>
      <th>1</th>
      <td>Q16512</td>
      <td>GO:0016301</td>
      <td>kinase</td>
    </tr>
  </tbody>
</table>
</div>



### User defined Protein Annotations
Users can also specify groups based on personal preferences. Therefore three arguments need to be specified:
* Required functions: A set of functions which a protein must have to be assigned to this group.
* Permitted functions: A set of functions of which must not overlap with the protein functions.
* A name

This class is also used to define the class `DefaultAnnotation`. The individual definitions can be found in the file
*go_protein_annotation/default_use.py*.
A simple example to separate protein kinases from other kinases and non-kinases:


```python
from go_protein_annotation import SpecialFunctionAnnotation
# Must have 'GO:0004672' (protein kinase activity)
# No permitted functions
# Name: "Protein kinase"
protein_kinases = ({"GO:0004672"}, set(), "Protein kinase")

# Must have 'GO:0004672' (kinase activity)
# Must not have '"GO:0004672' (protein kinase activity)
# Name: "Other kinase"
other_kinases = ({"GO:0016301"}, {"GO:0004672"}, "Other kinase")

# No required functions (all proteins would match this)
# Must not have '"GO:0016301' (kinase activity)
# Name: "Non-kinase"
non_kinases = (set(), {"GO:0016301"}, "Non-kinase")

example_classification = SpecialFunctionAnnotation([protein_kinases, other_kinases, non_kinases])

test_protein_annotations = example_classification.annotate_proteins(test_proteins)
```


```python
test_protein_annotations
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>uniprot_id</th>
      <th>protein_function</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>Q16512</td>
      <td>Protein kinase</td>
    </tr>
    <tr>
      <th>1</th>
      <td>P30085</td>
      <td>Other kinase</td>
    </tr>
    <tr>
      <th>2</th>
      <td>P25774</td>
      <td>Non-kinase</td>
    </tr>
  </tbody>
</table>
</div>



### Default Function Definition
See *go_protein_annotation/default_use.py*. Explicit explanation will follow.
### Miscellaneous
* Only QuickGO protein functions are used. QuckGO also gives information about involvement in biological
processes. These annotations are not considered.
* The classes `AllFunctionAnnotation` and `SelectedFunctionAnnotation` accept the keyword `alternative_name_dict`
  * Keys: GO ID
  * Value: Alternative name
* The classes `AllFunctionAnnotation` and `SelectedFunctionAnnotation` accept the keyword `simplify_name`
  * True (default): " activity" is removed from each protein function name (e.g. "kinase activity" -> "kinase")
  * False: protein functions are named as given by QuickGO
