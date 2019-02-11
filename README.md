# COMP546_Project2018
An algorithm to achieve dataset reduction by Clustering 3D structures of protein-protein interaction interfaces

**Protein-Protein interactions (PPIs)** are crucial for many biological processes in living organisms. Their importance necessitates the identification and prediction of novel interactions. One of the most common computational approaches for PPI prediction is *template-based structural alignment* which uses a set of experimentally known PPI interfaces (i.e., templates) to predict new interactions. The bottleneck in template-based prediction is the computational cost due to the one-to-all comparison of the query protein against a database of all known interfaces which includes tens of thousands of templates. In this project, I implemented an algorithm which reduces the list of known interaction interfces templates without losing its structural variation; this results in a reduced set which then can be used for a faster template-based prediction.

A **_template interface_** represents an interaction interface between a pair of proteins. Thus, it is important to note that in the case of having similar 3D structures of two templates representing interactions between different protein pairs, both templates should be included in the set since it is important for researches and scientists to know about both possibilities. In order to ensure this, the full templates list was divided into PPI group; each group contains the template interfaces of only one protein pair. For example, group *A-B* is the set of all templates which represent interfaces of the interaction between protein *A* and protein *B*.

# Inputs:
The algorithm **CLUSTER** takes as an input a set of templates (one of the PPI groups) and their similarity matrix. These inputs are then used to construct an undirected graph *G(V, E)* where each vertex represents a template interface and an edge *(u, v)* exists iff the structural similarity between template v and template u is above a certain threshold SIMILARITY_THRESHOLD.

# Output:
**CLUSTER** outputs a set of vertices which provides a minimum structural representation for the whole graph G.
Thus, after the full template list is grouped by protein pairs, each subset of templates is input to CLUSTER and all representatives are collected to obtain the minimum representatives set for the full templates list.

# Example:

Assume a graph G has **5 nodes** (a list of 5 templates) with the following similarity scores matrix:

|     |  A  |  B  |  C  |  D  |  E  |
| --- |:---:|:---:|:---:|:---:|:---:|
|  A  |1.000|0.287|0.196|0.268|0.280|
|  B  |0.287|1.000|**0.885**|**0.879**|**0.851**|
|  C  |0.196|**0.885**|1.000|**0.826**|**0.802**|
|  D  |0.268|**0.879**|**0.825**|1.000|0.789|
|  E  |0.280|**0.851**|**0.802**|0.789|1.000|

If we apply a *SIMILARITY_THRESHOLD* of **_0.8_**, then we can obtain the input graph from the adjacency matrix as follows:

|     |  A  |  B  |  C  |  D  |  E  |             
| --- |:---:|:---:|:---:|:---:|:---:|
|  A  |     |0    |0    |    0|    0|
|  B  |0    |     |**1**|**1**|**1**|
|  C  |0    |**1**|     |**1**|**1**|
|  D  |0    |**1**|**1**|     |0    |
|  E  |0    |**1**|**1**|0    |     |

![alt text](https://github.com/ahakouz17/COMP546_Project2018/blob/master/GraphGExample.png "Examples Graph")

The algorithm **CLUSTER** will output a set of *representatives* **R**; R = {A, B}.
𝑠.𝑡: 𝑇ℎ𝑒𝑟𝑒 𝑖𝑠 𝑎𝑛 𝑒𝑑𝑔𝑒 (𝑢,𝑣) ∈𝐸; ∀ 𝑣 ∈(𝑉−𝑅), 𝑤ℎ𝑒𝑟𝑒 𝑢 ∈ 𝑅
i.e., every vertex in G which was not picked as a representative should be connected (i.e., structurally similar to) at least one of the representative vertices.
The main goal is to minimize the size of set **R** to remove redundancy and increase the speed of predictions.

![alt text](https://github.com/ahakouz17/COMP546_Project2018/blob/master/OutputExample.png "Example Output")


