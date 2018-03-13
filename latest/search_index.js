var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#BPJSpec.TransferMatrix",
    "page": "Home",
    "title": "BPJSpec.TransferMatrix",
    "category": "type",
    "text": "struct TransferMatrix <: AbstractBlockMatrix{Matrix{Complex128}, 2}\n\nThis type represents the transfer matrix of an interferometer. This matrix effectively describes how an interferometer responds to the sky, including the antenna primary beam, bandpass, and baseline distribution.\n\nThis matrix is hierarchical in the sense that we save on some computational and storage requirements by separating long baselines from short baselines.\n\nFields\n\nstorage contains instructions on how to read the matrix from disk\ncache is used if we want to keep the matrix in memory\nmetadata stores the interferometer\'s metadata\nlmax is the largest value of the l quantum number used by the matrix\nmmax is the largest value of the m quantum number used by the matrix\n\n\n\n"
},

{
    "location": "#BPJSpec.MModes",
    "page": "Home",
    "title": "BPJSpec.MModes",
    "category": "type",
    "text": "struct MModes{S} <: AbstractBlockMatrix{Vector{Complex128}, 2}\n\nThis type represents the m-modes measured by the interferometer.\n\nFields\n\nstorage contains instructions on how to read the m-modes from disk\ncache is used if we want to keep the m-modes in memory\nmmax is the largest value of the m quantum number\nfrequencies is the list of frequencies\nbandwidth is the bandwidth associated with each frequency channel\n\n\n\n"
},

{
    "location": "#BPJSpec-Documentation-1",
    "page": "Home",
    "title": "BPJSpec Documentation",
    "category": "section",
    "text": "BPJSpec is a 21-cm power spectrum code developed for the OVRO-LWA based on the m-mode analysis formalism.CurrentModule = BPJSpec\nDocTestSetup = quote\n    using BPJSpec\nendTransferMatrix\nMModes"
},

{
    "location": "#BPJSpec.AbstractBlockMatrix",
    "page": "Home",
    "title": "BPJSpec.AbstractBlockMatrix",
    "category": "type",
    "text": "abstract type AbstractBlockMatrix{B, N}\n\nThis type represents a (potentially enormous) block-diagonal matrix. This type is designed to be general enough to handle large matrices that fit in memory as well as enormous matrices that do not fit in memory. In principle this type can also be used to store small matrices, but it would be relatively inefficient compared to the standard Array{T, N}.\n\nType Parameters\n\nB specifies the type of the blocks that compose the matrix\nN specifies the number of indices used to index into the matrix blocks\n\nRequired Fields\n\nstorage contains instructions on how to read matrix blocks\ncache is used if we want to read the matrix from disk and then keep it in memory for faster   access.\n\n\n\n"
},

{
    "location": "#Internals-1",
    "page": "Home",
    "title": "Internals",
    "category": "section",
    "text": "AbstractBlockMatrix"
},

]}
