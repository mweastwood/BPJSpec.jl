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

{
    "location": "wrappers/#",
    "page": "Wrappers",
    "title": "Wrappers",
    "category": "page",
    "text": ""
},

{
    "location": "wrappers/#Wrappers-1",
    "page": "Wrappers",
    "title": "Wrappers",
    "category": "section",
    "text": ""
},

{
    "location": "wrappers/#BPJSpec.CosmologyWrapper.comoving_distance",
    "page": "Wrappers",
    "title": "BPJSpec.CosmologyWrapper.comoving_distance",
    "category": "function",
    "text": "comoving_distance(z)\n\nCalculate the comoving distance (in units of Mpc) to the redshift z.\n\nUsage:\n\njulia> comoving_distance(1)\n3371.509961954628 Mpc\n\njulia> comoving_distance(10)\n9689.514711746533 Mpc\n\n\n\n"
},

{
    "location": "wrappers/#BPJSpec.CosmologyWrapper.age",
    "page": "Wrappers",
    "title": "BPJSpec.CosmologyWrapper.age",
    "category": "function",
    "text": "age(z)\n\nCalculate the age of the universe (in units of Gyr) to the redshift z.\n\nUsage:\n\njulia> age(1)\n5.918077173774152 Gyr\n\njulia> age(10)\n0.4785005773464139 Gyr\n\n\n\n"
},

{
    "location": "wrappers/#BPJSpec.CosmologyWrapper.frequency",
    "page": "Wrappers",
    "title": "BPJSpec.CosmologyWrapper.frequency",
    "category": "function",
    "text": "frequency(z)\n\nCalculate the frequency of the 21 cm line of Hydrogen at the redshift z.\n\nUsage:\n\njulia> frequency(1)\n710.202875885 MHz\n\njulia> frequency(10)\n129.12779561545454 MHz\n\n\n\n"
},

{
    "location": "wrappers/#BPJSpec.CosmologyWrapper.redshift",
    "page": "Wrappers",
    "title": "BPJSpec.CosmologyWrapper.redshift",
    "category": "function",
    "text": "redshift(ν)\n\nCalculate the redshift from which the emission originates if the 21 cm line is observed at the frequency ν.\n\nUsage:\n\njulia> redshift(100u\"MHz\")\n13.2040575177\n\njulia> redshift(200u\"MHz\")\n6.10202875885\n\n\n\n"
},

{
    "location": "wrappers/#CosmologyWrapper.jl-1",
    "page": "Wrappers",
    "title": "CosmologyWrapper.jl",
    "category": "section",
    "text": "CurrentModule = BPJSpec\nDocTestSetup = quote\n    using BPJSpec.CosmologyWrapper\n    using Unitful\nendcomoving_distance\nage\nfrequency\nredshift"
},

{
    "location": "wrappers/#FastTransformsWrapper.jl-1",
    "page": "Wrappers",
    "title": "FastTransformsWrapper.jl",
    "category": "section",
    "text": ""
},

]}
