# QuckGOProteinAnnotation
Protein function annotations provided by [QuickGO](https://www.ebi.ac.uk/QuickGO/) can be used to
classify proteins specified by [UniProt](https://www.uniprot.org/) ID.
So far, currently class-assignment are non-unique, so a protein can be

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
```
python annotate_protein_list.py -i demo_data/demo_uniprot_ids.tsv -o demo_data/demo_output.tsv -c "uniprot_id" -s tab
```
### In Python


```python
from go_protein_annotation  import DefaultAnnotation
from go_protein_annotation import  AllFunctionAnnotation
```


```python
test_proteins = ["Q16512", "P30085", "P25774"]
```


```python
default_annotation = DefaultAnnotation()
protein_class_df = default_annotation.annotate_proteins(test_proteins)
```


```python
protein_class_df

```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>uniprot_id</th>
      <th>functions</th>
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
[peptidase activity](https://www.ebi.ac.uk/QuickGO/term/GO:0008233) implies a
[hydrolase activity](https://www.ebi.ac.uk/QuickGO/term/GO:0016787) as well. Provided code extracts
all explicit functional annotations and extends it with implicit annotations.


```python
all_functions = AllFunctionAnnotation()
all_functions.get_protein_functions("Q16512")

```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>uniprot_id</th>
      <th>go_id</th>
      <th>name</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>Q16512</td>
      <td>GO:0032559</td>
      <td>adenyl ribonucleotide binding</td>
    </tr>
    <tr>
      <th>1</th>
      <td>Q16512</td>
      <td>GO:0003682</td>
      <td>chromatin binding</td>
    </tr>
    <tr>
      <th>2</th>
      <td>Q16512</td>
      <td>GO:0009931</td>
      <td>calcium-dependent protein serine/threonine kinase</td>
    </tr>
    <tr>
      <th>3</th>
      <td>Q16512</td>
      <td>GO:0043168</td>
      <td>anion binding</td>
    </tr>
    <tr>
      <th>4</th>
      <td>Q16512</td>
      <td>GO:0004697</td>
      <td>protein kinase C</td>
    </tr>
    <tr>
      <th>5</th>
      <td>Q16512</td>
      <td>GO:0004672</td>
      <td>protein kinase</td>
    </tr>
    <tr>
      <th>6</th>
      <td>Q16512</td>
      <td>GO:0010857</td>
      <td>calcium-dependent protein kinase</td>
    </tr>
    <tr>
      <th>7</th>
      <td>Q16512</td>
      <td>GO:0016922</td>
      <td>nuclear receptor binding</td>
    </tr>
    <tr>
      <th>8</th>
      <td>Q16512</td>
      <td>GO:0005102</td>
      <td>signaling receptor binding</td>
    </tr>
    <tr>
      <th>9</th>
      <td>Q16512</td>
      <td>GO:0016740</td>
      <td>transferase</td>
    </tr>
    <tr>
      <th>10</th>
      <td>Q16512</td>
      <td>GO:0016301</td>
      <td>kinase</td>
    </tr>
    <tr>
      <th>11</th>
      <td>Q16512</td>
      <td>GO:0000166</td>
      <td>nucleotide binding</td>
    </tr>
    <tr>
      <th>12</th>
      <td>Q16512</td>
      <td>GO:0097159</td>
      <td>organic cyclic compound binding</td>
    </tr>
    <tr>
      <th>13</th>
      <td>Q16512</td>
      <td>GO:0032555</td>
      <td>purine ribonucleotide binding</td>
    </tr>
    <tr>
      <th>14</th>
      <td>Q16512</td>
      <td>GO:0035258</td>
      <td>steroid hormone receptor binding</td>
    </tr>
    <tr>
      <th>15</th>
      <td>Q16512</td>
      <td>GO:0098772</td>
      <td>molecular function regulator</td>
    </tr>
    <tr>
      <th>16</th>
      <td>Q16512</td>
      <td>GO:0140096</td>
      <td>catalytic activity, acting on a protein</td>
    </tr>
    <tr>
      <th>17</th>
      <td>Q16512</td>
      <td>GO:0050681</td>
      <td>androgen receptor binding</td>
    </tr>
    <tr>
      <th>18</th>
      <td>Q16512</td>
      <td>GO:0005515</td>
      <td>protein binding</td>
    </tr>
    <tr>
      <th>19</th>
      <td>Q16512</td>
      <td>GO:0008134</td>
      <td>transcription factor binding</td>
    </tr>
    <tr>
      <th>20</th>
      <td>Q16512</td>
      <td>GO:1901363</td>
      <td>heterocyclic compound binding</td>
    </tr>
    <tr>
      <th>21</th>
      <td>Q16512</td>
      <td>GO:0019901</td>
      <td>protein kinase binding</td>
    </tr>
    <tr>
      <th>22</th>
      <td>Q16512</td>
      <td>GO:0005488</td>
      <td>binding</td>
    </tr>
    <tr>
      <th>23</th>
      <td>Q16512</td>
      <td>GO:0042393</td>
      <td>histone binding</td>
    </tr>
    <tr>
      <th>24</th>
      <td>Q16512</td>
      <td>GO:0030374</td>
      <td>nuclear receptor coactivator</td>
    </tr>
    <tr>
      <th>25</th>
      <td>Q16512</td>
      <td>GO:0032553</td>
      <td>ribonucleotide binding</td>
    </tr>
    <tr>
      <th>26</th>
      <td>Q16512</td>
      <td>GO:0003824</td>
      <td>catal</td>
    </tr>
    <tr>
      <th>27</th>
      <td>Q16512</td>
      <td>GO:0035184</td>
      <td>histone threonine kinase</td>
    </tr>
    <tr>
      <th>28</th>
      <td>Q16512</td>
      <td>GO:0030554</td>
      <td>adenyl nucleotide binding</td>
    </tr>
    <tr>
      <th>29</th>
      <td>Q16512</td>
      <td>GO:0035639</td>
      <td>purine ribonucleoside triphosphate binding</td>
    </tr>
    <tr>
      <th>30</th>
      <td>Q16512</td>
      <td>GO:0005080</td>
      <td>protein kinase C binding</td>
    </tr>
    <tr>
      <th>31</th>
      <td>Q16512</td>
      <td>GO:0004674</td>
      <td>protein serine/threonine kinase</td>
    </tr>
    <tr>
      <th>32</th>
      <td>Q16512</td>
      <td>GO:0004698</td>
      <td>calcium-dependent protein kinase C</td>
    </tr>
    <tr>
      <th>33</th>
      <td>Q16512</td>
      <td>GO:0005524</td>
      <td>ATP binding</td>
    </tr>
    <tr>
      <th>34</th>
      <td>Q16512</td>
      <td>GO:0017076</td>
      <td>purine nucleotide binding</td>
    </tr>
    <tr>
      <th>35</th>
      <td>Q16512</td>
      <td>GO:0003713</td>
      <td>transcription coactivator</td>
    </tr>
    <tr>
      <th>36</th>
      <td>Q16512</td>
      <td>GO:0042826</td>
      <td>histone deacetylase binding</td>
    </tr>
    <tr>
      <th>37</th>
      <td>Q16512</td>
      <td>GO:0019899</td>
      <td>enzyme binding</td>
    </tr>
    <tr>
      <th>38</th>
      <td>Q16512</td>
      <td>GO:1901265</td>
      <td>nucleoside phosphate binding</td>
    </tr>
    <tr>
      <th>39</th>
      <td>Q16512</td>
      <td>GO:0019900</td>
      <td>kinase binding</td>
    </tr>
    <tr>
      <th>40</th>
      <td>Q16512</td>
      <td>GO:0035173</td>
      <td>histone kinase</td>
    </tr>
    <tr>
      <th>41</th>
      <td>Q16512</td>
      <td>GO:0140110</td>
      <td>transcription regulator</td>
    </tr>
    <tr>
      <th>42</th>
      <td>Q16512</td>
      <td>GO:0051020</td>
      <td>GTPase binding</td>
    </tr>
    <tr>
      <th>43</th>
      <td>Q16512</td>
      <td>GO:0043167</td>
      <td>ion binding</td>
    </tr>
    <tr>
      <th>44</th>
      <td>Q16512</td>
      <td>GO:0140297</td>
      <td>DNA-binding transcription factor binding</td>
    </tr>
    <tr>
      <th>45</th>
      <td>Q16512</td>
      <td>GO:0097367</td>
      <td>carbohydrate derivative binding</td>
    </tr>
    <tr>
      <th>46</th>
      <td>Q16512</td>
      <td>GO:0016772</td>
      <td>transferase activity, transferring phosphorus-...</td>
    </tr>
    <tr>
      <th>47</th>
      <td>Q16512</td>
      <td>GO:0061629</td>
      <td>RNA polymerase II-specific DNA-binding transcr...</td>
    </tr>
    <tr>
      <th>48</th>
      <td>Q16512</td>
      <td>GO:0035402</td>
      <td>histone kinase activity (H3-T11 specific)</td>
    </tr>
    <tr>
      <th>49</th>
      <td>Q16512</td>
      <td>GO:0016773</td>
      <td>phosphotransferase activity, alcohol group as ...</td>
    </tr>
    <tr>
      <th>50</th>
      <td>Q16512</td>
      <td>GO:0036094</td>
      <td>small molecule binding</td>
    </tr>
    <tr>
      <th>51</th>
      <td>Q16512</td>
      <td>GO:0031267</td>
      <td>small GTPase binding</td>
    </tr>
    <tr>
      <th>52</th>
      <td>Q16512</td>
      <td>GO:0003712</td>
      <td>transcription coregulator</td>
    </tr>
    <tr>
      <th>53</th>
      <td>Q16512</td>
      <td>GO:0051427</td>
      <td>hormone receptor binding</td>
    </tr>
    <tr>
      <th>54</th>
      <td>Q16512</td>
      <td>GO:0035257</td>
      <td>nuclear hormone receptor binding</td>
    </tr>
  </tbody>
</table>
</div>


