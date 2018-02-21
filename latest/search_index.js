var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#BPJSpec.HierarchicalTransferMatrix",
    "page": "Home",
    "title": "BPJSpec.HierarchicalTransferMatrix",
    "category": "Type",
    "text": "struct HierarchicalTransferMatrix\n\nThis type represents the transfer matrix of an interferometer. This matrix effectively describes how an interferometer responds to the sky, including the antenna primary beam, bandpass, and baseline distribution.\n\n\"Hierarchical\" refers to the fact that we save on some computational and storage requirements by separating long baselines from short baselines.\n\nFields\n\npath points to the directory where the matrix is stored\nmetadata describes the properties of the interferometer\nhierarchy describes how the baselines are grouped\nfrequencies is an alias for metadata.frequencies\nbandwidth is an alias for metadata.bandwidth\nlmax is the maximum value of the total angular momentum quantum number l\nmmax is the maximum value of the azimuthal quantum number m\n\nImplementation\n\nAll of the data is stored on disk and only read into memory on-request. Generally, this approach is necessary because the entire transfer matrix is too large to entirely fit in memory, and because the matrix is block diagonal we can work with blocks individually.\n\n\n\n"
},

{
    "location": "#BPJSpec-Documentation-1",
    "page": "Home",
    "title": "BPJSpec Documentation",
    "category": "section",
    "text": "BPJSpec is a 21-cm power spectrum code developed for the OVRO-LWA based on the m-mode analysis formalism.CurrentModule = BPJSpec\nDocTestSetup = quote\n    using BPJSpec\nendHierarchicalTransferMatrix"
},

]}
