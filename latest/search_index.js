var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#BPJSpec-Documentation-1",
    "page": "Home",
    "title": "BPJSpec Documentation",
    "category": "section",
    "text": "BPJSpec is a 21-cm power spectrum code developed for the OVRO-LWA based on the m-mode analysis formalism."
},

{
    "location": "noise-model/#",
    "page": "Noise Model",
    "title": "Noise Model",
    "category": "page",
    "text": ""
},

{
    "location": "noise-model/#Noise-Model-1",
    "page": "Noise Model",
    "title": "Noise Model",
    "category": "section",
    "text": ""
},

{
    "location": "noise-model/#BPJSpec.NoiseModel",
    "page": "Noise Model",
    "title": "BPJSpec.NoiseModel",
    "category": "type",
    "text": "struct NoiseModel\n\nThis type represents the thermal noise contributed to the measurement of a set of m-modes.\n\nA careful reading of Taylor, Carilli, Perley chapter 9 reveals that under the convention that Stokes-I is (rm xx + yy)2 we get the following expressions:\n\nUncertainty on single polarization visibilities (in flux density units):\n\n_rm xx = fracsqrt2 k T_rm sysA_e sqrt()\n\nUncertainty on Stokes I visibilities (in flux density units):\n\n_rm Stokes-I = frac_rm xxsqrt2 = frack T_rm sysA_e sqrt()\n\nWhere k is the Boltzmann constant, T_rm sys is the system temperature, A_e is the effective collecting area,  is the bandwidth, and  is the integration time.\n\nHowever, for a dipole antenna, the effective collecting area is not a very physically meaningful value. However, it turns out that we can relate the effective collecting are to the solid angle subtended by the primary beam :\n\nA_e = frac^2\n\nnote: Note\nThere seems to be some ambiguity in the literature in regards to notation. I believe we originally assumed that A_e refers to the maximum effective collecting area, and that we have normalized the beam to be unity in that direction.\n\nFinally we end  up with the following expression after including an additional contribution due to time smearing:\n\n_rm m-modes = frack T_rm sys ^2 sqrtN_rm int\n                  sincleft(fracmtautextsidereal dayright)\n\nFields:\n\nTsys specifies the system temperature\nτ specfies the length of a single integration\nNint specifies the total number of integrations used in the dataset\nΩ is the solid angle subtended by the primary beam\n\nUsage:\n\njulia> model = BPJSpec.NoiseModel(1000u\"K\", 13u\"s\", 6628, 2.41u\"sr\")\nNoiseModel(Tsys = 1000.0 K, τ = 13.0 s, Nint = 6628, Ω = 2.410 sr)\n\njulia> model(100, 74u\"MHz\", 24u\"kHz\")\n4.456470453155544 Jy\n\n\n\n"
},

{
    "location": "noise-model/#API-1",
    "page": "Noise Model",
    "title": "API",
    "category": "section",
    "text": "CurrentModule = BPJSpec\nDocTestSetup = quote\n    using BPJSpec\n    using Unitful, UnitfulAstro\nendNoiseModel"
},

{
    "location": "foreground-filtering/#",
    "page": "Foreground Filtering",
    "title": "Foreground Filtering",
    "category": "page",
    "text": ""
},

{
    "location": "foreground-filtering/#Foreground-Filtering-1",
    "page": "Foreground Filtering",
    "title": "Foreground Filtering",
    "category": "section",
    "text": ""
},

{
    "location": "foreground-filtering/#BPJSpec.ForegroundComponent",
    "page": "Foreground Filtering",
    "title": "BPJSpec.ForegroundComponent",
    "category": "type",
    "text": "struct ForegroundComponent <: SkyComponent\n\nThis type represents the contribution of a single foreground component to the multi-frequency angular power spectrum.\n\nC_l(_1 _2) = A left(fracl+11001right)^-\n                  left(frac_1 _2_0^2right)^-\n                  expleft(-fraclog^2(_1_2)2^2right)\n\nFields:\n\nν0 specifies the reference frequency\nA specifies the overall amplitude at the reference frequency\nα is the power-law index for the multipole moment l\nβ is the power-law index for frequency\nζ essentially determines how quickly the foreground component decorrelates in frequency\n\nUsage:\n\nSome quick foreground models can be constructed based on the work of Santos, Cooray & Knox 2005. These can be accessed with the functions extragalactic_point_sources, extragalactic_free_free, galactic_synchrotron, and galactic_free_free.\n\njulia> BPJSpec.extragalactic_point_sources()\nForegroundComponent(ν0 = 130.000 MHz, A = 57.000 mK², α = 1.100, β = 2.070, ζ = 1.000)\n\njulia> BPJSpec.extragalactic_free_free()\nForegroundComponent(ν0 = 130.000 MHz, A = 0.014 mK², α = 1.000, β = 2.100, ζ = 35.000)\n\njulia> BPJSpec.galactic_synchrotron()\nForegroundComponent(ν0 = 130.000 MHz, A = 700.000 mK², α = 2.400, β = 2.800, ζ = 4.000)\n\njulia> BPJSpec.galactic_free_free()\nForegroundComponent(ν0 = 130.000 MHz, A = 0.088 mK², α = 3.000, β = 2.150, ζ = 35.000)\n\njulia> component = BPJSpec.ForegroundComponent(100u\"MHz\", 1u\"K^2\", 1, 1, 100)\n       component(100, 74u\"MHz\", 76u\"MHz\")\n17.622494197510505 K^2\n\n\n\n"
},

{
    "location": "foreground-filtering/#BPJSpec.CylindricalPS",
    "page": "Foreground Filtering",
    "title": "BPJSpec.CylindricalPS",
    "category": "type",
    "text": "struct CylindricalPS <: PowerSpectrum <: SkyComponent\n\nThis type represents a cylindrically averaged power spectrum of the 21-cm brightness temperature over a given range in redshift. An instance of this type can be evaluated to compute an approximation of the multi-frequency angular power spectrum.\n\nC_l(_1 _2)  \n\nFields:\n\nzrange specifies the range over which this power spectrum is valid. It is used to define the range over which cosmological quantities can be approximated\nkpara is a list of wave numbers parallel to the line-of-sight\nkperp is a list of wave numbers perpendicular to the line-of-sight\npower is the amplitude of the power spectrum at the grid points specified by kpara and kperp\n\nUsage:\n\njulia> power_spectrum = BPJSpec.CylindricalPS((10, 30),\n                                              [0.0, 1.0].*u\"Mpc^-1\", # kpara\n                                              [0.0, 1.0].*u\"Mpc^-1\", # kperp\n                                              [1.0 1.0; 1.0 1.0].*u\"mK^2*Mpc^3\")\nCylindricalPS(10.0 < z < 30.0, k∥ = 0.00 Mpc⁻¹…1.00 Mpc⁻², k⟂ = 0.00 Mpc⁻¹…1.00 Mpc⁻², P ~ 1.0 mK²Mpc³)\n\njulia> power_spectrum(100, 74u\"MHz\", 74u\"MHz\")\n2.694677277820449e-15 K^2\n\njulia> power_spectrum(100, 74u\"MHz\", 75u\"MHz\")\n-4.315090053966302e-17 K^2\n\n\n\n"
},

{
    "location": "foreground-filtering/#BPJSpec.extragalactic_point_sources",
    "page": "Foreground Filtering",
    "title": "BPJSpec.extragalactic_point_sources",
    "category": "function",
    "text": "extragalactic_point_sources()\n\nConstructs a model of extragalactic point sources based on the work of Santos, Cooray & Knox 2005.\n\n\n\n"
},

{
    "location": "foreground-filtering/#BPJSpec.extragalactic_free_free",
    "page": "Foreground Filtering",
    "title": "BPJSpec.extragalactic_free_free",
    "category": "function",
    "text": "extragalactic_free_free()\n\nConstructs a model of extragalactic free-free emission based on the work of Santos, Cooray & Knox 2005.\n\n\n\n"
},

{
    "location": "foreground-filtering/#BPJSpec.galactic_synchrotron",
    "page": "Foreground Filtering",
    "title": "BPJSpec.galactic_synchrotron",
    "category": "function",
    "text": "galactic_synchrotron()\n\nConstructs a model of galactic synchrotron emission based on the work of Santos, Cooray & Knox 2005.\n\n\n\n"
},

{
    "location": "foreground-filtering/#BPJSpec.galactic_free_free",
    "page": "Foreground Filtering",
    "title": "BPJSpec.galactic_free_free",
    "category": "function",
    "text": "galactic_free_free()\n\nConstructs a model of galactic free-free emission based on the work of Santos, Cooray & Knox 2005.\n\n\n\n"
},

{
    "location": "foreground-filtering/#API-1",
    "page": "Foreground Filtering",
    "title": "API",
    "category": "section",
    "text": "CurrentModule = BPJSpec\nDocTestSetup = quote\n    using BPJSpec\n    using Unitful, UnitfulAstro\nendForegroundComponent\nCylindricalPS\nextragalactic_point_sources\nextragalactic_free_free\ngalactic_synchrotron\ngalactic_free_free"
},

{
    "location": "enormous-matrices/#",
    "page": "Enormous Matrices",
    "title": "Enormous Matrices",
    "category": "page",
    "text": ""
},

{
    "location": "enormous-matrices/#Enormous-Matrices-1",
    "page": "Enormous Matrices",
    "title": "Enormous Matrices",
    "category": "section",
    "text": ""
},

{
    "location": "enormous-matrices/#BPJSpec.AbstractBlockMatrix",
    "page": "Enormous Matrices",
    "title": "BPJSpec.AbstractBlockMatrix",
    "category": "type",
    "text": "abstract type AbstractBlockMatrix{B, N}\n\nThis type represents a (potentially enormous) block-diagonal matrix. This type is designed to be general enough to handle large matrices that fit in memory as well as enormous matrices that do not fit in memory. In principle this type can also be used to store small matrices, but it would be relatively inefficient compared to the standard Array{T, N}.\n\nType Parameters\n\nB specifies the type of the blocks that compose the matrix\nN specifies the number of indices used to index into the matrix blocks\n\nRequired Fields\n\nstorage contains instructions on how to read matrix blocks\ncache is used if we want to read the matrix from disk and then keep it in memory for faster   access.\n\n\n\n"
},

{
    "location": "enormous-matrices/#BPJSpec.TransferMatrix",
    "page": "Enormous Matrices",
    "title": "BPJSpec.TransferMatrix",
    "category": "type",
    "text": "struct TransferMatrix <: AbstractBlockMatrix{Matrix{Complex128}, 2}\n\nThis type represents the transfer matrix of an interferometer. This matrix effectively describes how an interferometer responds to the sky, including the antenna primary beam, bandpass, and baseline distribution.\n\nThis matrix is hierarchical in the sense that we save on some computational and storage requirements by separating long baselines from short baselines.\n\nFields\n\nstorage contains instructions on how to read the matrix from disk\ncache is used if we want to keep the matrix in memory\nmetadata stores the interferometer\'s metadata\nlmax is the largest value of the l quantum number used by the matrix\nmmax is the largest value of the m quantum number used by the matrix\n\n\n\n"
},

{
    "location": "enormous-matrices/#BPJSpec.MModes",
    "page": "Enormous Matrices",
    "title": "BPJSpec.MModes",
    "category": "type",
    "text": "struct MModes{S} <: AbstractBlockMatrix{Vector{Complex128}, 2}\n\nThis type represents the m-modes measured by the interferometer.\n\nFields\n\nstorage contains instructions on how to read the m-modes from disk\ncache is used if we want to keep the m-modes in memory\nmmax is the largest value of the m quantum number\nfrequencies is the list of frequencies\nbandwidth is the bandwidth associated with each frequency channel\n\n\n\n"
},

{
    "location": "enormous-matrices/#API-1",
    "page": "Enormous Matrices",
    "title": "API",
    "category": "section",
    "text": "CurrentModule = BPJSpec\nDocTestSetup = quote\n    using BPJSpec\nendAbstractBlockMatrix\nTransferMatrix\nMModes"
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
