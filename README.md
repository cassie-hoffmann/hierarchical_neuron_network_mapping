# Hierarchical mapping of neuron network morphology

This pipeline quantifies the network structure of in vitro neuron ensembles from optical microscopy images. It can be implemented to characterise neuron self-organisation behaviours across healthy or perturbed states. For example, the segmentations in Figure 1 show healthy neurons organising from random into modular network structure, which increases local and global network efficiency.
Within the pipeline, a network graph is generated from neuron morphology where somata are encoded as nodes (displayed in orange) and neurites as edges (displayed in grey). Graph theoretical analysis is then used to quantify key features of network topology.

Figure 1.
<img width="1035" height="419" alt="Screenshot 2026-04-09 at 3 37 13 pm" src="https://github.com/user-attachments/assets/9fecdcee-909c-456e-82a7-c48e3f649e87" />


Our pipeline has two modes that are compatible with different microscopy datasets. The first is a ‘hierarchical mode’ that generates multiresolution networks by combining micro- and meso-scale images of neuron architecture. The second is a ‘unimodal mode’ that generates single-resolution networks from mesoscale images alone, and is rapidly deployable in high-throughput datasets.  

Requirements:
-	Both modes require the Brain Connectivity Toolbox (https://sites.google.com/site/bctnet/), which should be downloaded to the same path as the pipeline.
  
Optional packages:
-	It is recommended to generate segmentations of somata and neurite structures with Ilastik (https://www.ilastik.org/) and use them as inputs into the pipeline. This machine learning-driven tool is highly adaptable to most fluorescence microscopy images, even with low signal to noise ratio. Our pipeline also accepts raw fluorescence images, but users should verify that thresholding procedures generate appropriate segmentations. 
-	In hierarchical mode, StarDist (https://github.com/stardist/stardist) may optionally be used to obtain segmentations where somata are individuated within clusters.
  
**Hierarchical mode**
In hierarchical mode, a multiscale network is constructed by embedding representative images of neurospheres within time-course images of neural populations. Although compatible with any optical microscopy modalities, the two input datasets used here are (i) immunostained neurospheres captured in 3D with super-resolution imaging and (ii) live fluorescent reporter populations captured in 2D with widefield imaging (Figure 1a). As clusters in population images have z-depths that are intractable with widefield techniques alone, these objects serve as sites for neurosphere insertion. On a network level, this involves embedding SR microscale graphs into widefield mesoscale frameworks, and generating a multiresolution graph as output (Figure 1b). Graph interconnectivity is represented as an adjacency matrix, which serves as a basis for network analysis (Figure 1c). 

Figure 2.
<img width="1607" height="864" alt="dskl" src="https://github.com/user-attachments/assets/05cdf32e-860b-483e-915e-30b7f0365799" />

For hierarchical mode, please run ‘network_reconstruction_multiscale.m’. The inputs are segmentations/raw images of 2D mesoscale neuron populations and 3D microscale neurospheres or clusters. The outputs are a .tif network visualisation and a .mat file of the following metrics: number of nodes, network edge density, node degree, community structure, community density, number of communities. More metrics may be optionally included from the Brain Connectivity Toolbox (https://sites.google.com/site/bctnet/).

**Unimodal mode**
For unimodal mode, please run ‘network_reconstruction_2d.m’. This mode is designed for rapid processing of single-modality 2D microscopy datasets, such as those containing widefield images alone. In the absence of 3D super-resolution imaging inputs, the pipeline interpolates cell configurations within clusters via morphological heuristics. 
The input is segmentations/raw images of 2D mesoscale neuron populations, and the outputs are the same as hierarchical mode above. To facilitate cluster interpolation, the user must specify (i) the expected size of a single neuron soma and (ii) the desired proportion of within-cluster connectivity. For example, Figure 2b shows a network ROI with 10% within-cluster connectivity, while Figure 2c shows 90% within-cluster connectivity. For proportions lower than 100%, edge placement follows a distance-decay law. 

Figure 3. 
<img width="530" height="488" alt="Screenshot 2026-04-09 at 5 15 31 pm" src="https://github.com/user-attachments/assets/53b3e8c8-7fe9-4c10-beb8-951f726c6f67" />





