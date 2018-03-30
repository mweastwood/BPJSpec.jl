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
    "text": "struct NoiseModel\n\nThis type represents the thermal noise contributed to the measurement of a set of m-modes.\n\nA careful reading of Taylor, Carilli, Perley chapter 9 reveals that under the convention that Stokes-I is (rm xx + yy)2 we get the following expressions:\n\nUncertainty on single polarization visibilities (in flux density units):\n\n_rm xx = fracsqrt2 k T_rm sysA_e sqrt\n\nUncertainty on Stokes I visibilities (in flux density units):\n\n_textStokes-I = frac_rm xxsqrt2 = frack T_rm sysA_e sqrt\n\nWhere k is the Boltzmann constant, T_rm sys is the system temperature, A_e is the effective collecting area,  is the bandwidth, and  is the integration time.\n\nHowever, for a dipole antenna, the effective collecting area is not a very physically meaningful value. However, it turns out that we can relate the effective collecting are to the solid angle subtended by the primary beam :\n\nA_e = frac^2\n\nnote: Note\nThere seems to be some ambiguity in the literature in regards to notation. I believe we originally assumed that A_e refers to the maximum effective collecting area, and that we have normalized the beam to be unity in that direction.\n\nFinally we end  up with the following expression after including an additional contribution due to time smearing:\n\n_textm-modes = frack T_rm sys ^2 sqrtN_rm int\n                   rm sincleft(fracmtautextsidereal dayright)\n\nFields:\n\nTsys specifies the system temperature\nτ specfies the length of a single integration\nNint specifies the total number of integrations used in the dataset\nΩ is the solid angle subtended by the primary beam\n\nUsage:\n\njulia> model = BPJSpec.NoiseModel(1000u\"K\", 13u\"s\", 6628, 2.41u\"sr\")\nNoiseModel(Tsys = 1000.0 K, τ = 13.0 s, Nint = 6628, Ω = 2.410 sr)\n\njulia> model(100, 74u\"MHz\", 24u\"kHz\")\n4.456470453155544 Jy\n\n\n\n"
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
    "location": "foreground-filtering/#BPJSpec.foreground_filter!",
    "page": "Foreground Filtering",
    "title": "BPJSpec.foreground_filter!",
    "category": "function",
    "text": "foreground_filter!(output_mmodes, output_transfermatrix, output_covariance,\n                   output_foreground_filter, output_noise_whitener,\n                   input_mmodes, input_transfermatrix, input_noisematrix,\n                   input_signalmatrix, input_foregroundmatrix;\n                   threshold=0.1, tempdir=\"/tmp\", cleanup=true)\n\nFilter foreground emission from the dataset and whiten the noise.\n\nArguments:\n\noutput_mmodes will be populated with the foreground filtered m-modes\noutput_transfermatrix will be populated with the foreground-filtered transfer-matrix\noutput_covariance will be populated with the covariance matrix of the output m-modes\noutput_foreground_filter is the matrix used to filter foreground emission\noutput_noise_whitener is the matrix used to whiten the noise covariance (the thermal noise and the foreground contribution)\ninput_mmodes is the measured input m-modes\ninput_transfermatrix is the transfer matrix describing the response of the interferometer\ninput_noisematrix is the thermal noise covariance matrix\ninput_signalmatrix is the covariance matrix due to the 21-cm power spectrum\ninput_foregroundmatrix is the covariance matrix due to the foreground radio emission\n\nwarn: Warn\nPlease pay attention to the argument order. Swapping arguments could lead to this function filtering all of the 21-cm signal while leaving the foreground emission in tact!\n\nKeyword Arguments:\n\nthreshold is the maximum allowed foreground-signal ratio\ntempdir is a path to a directory where temporary files can be placed\ncleanup determines whether or not files placed in the tempdir are automatically removed. This can be set to false if you wish to check some of the intermediate output.\n\n\n\n"
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
    "text": "CurrentModule = BPJSpec\nDocTestSetup = quote\n    using BPJSpec\n    using Unitful, UnitfulAstro\nendforeground_filter!\nForegroundComponent\nCylindricalPS\nextragalactic_point_sources\nextragalactic_free_free\ngalactic_synchrotron\ngalactic_free_free"
},

{
    "location": "imaging/#",
    "page": "Imaging",
    "title": "Imaging",
    "category": "page",
    "text": ""
},

{
    "location": "imaging/#Imaging-1",
    "page": "Imaging",
    "title": "Imaging",
    "category": "section",
    "text": ""
},

{
    "location": "imaging/#BPJSpec.tikhonov",
    "page": "Imaging",
    "title": "BPJSpec.tikhonov",
    "category": "function",
    "text": "tikhonov(transfermatrix, mmodes; regularization=1e-2, mfs=false, storage=NoFile())\n\nCreate a dirty image of the sky using Tikhonov regularization.\n\nArguments:\n\ntransfermatrix the interferometer\'s transfer matrix, describing its response to the sky\nmmodes the m-modes measured by the interferometer\n\nKeyword Arguments:\n\nregularization the amplitude of the Tikhonov regularization parameter\nmfs determines whether or not to perform Multi-Frequency Synthesis imaging. If this parameter is set to true, all frequency channels will be used to generate a single image of the sky. If this parameter is set to false, an image of the sky will be generated for each frequency channel.\nstorage determines how the computed spherical harmonic coefficients will be stored\n\n\n\n"
},

{
    "location": "imaging/#API-1",
    "page": "Imaging",
    "title": "API",
    "category": "section",
    "text": "CurrentModule = BPJSpec\nDocTestSetup = quote\n    using BPJSpec\nendtikhonov"
},

{
    "location": "power-spectrum-estimation/#",
    "page": "Power Spectrum Estimation",
    "title": "Power Spectrum Estimation",
    "category": "page",
    "text": ""
},

{
    "location": "power-spectrum-estimation/#Power-Spectrum-Estimation-1",
    "page": "Power Spectrum Estimation",
    "title": "Power Spectrum Estimation",
    "category": "section",
    "text": ""
},

{
    "location": "power-spectrum-estimation/#BPJSpec.q_estimator",
    "page": "Power Spectrum Estimation",
    "title": "BPJSpec.q_estimator",
    "category": "function",
    "text": "q_estimator(mmodes, transfermatrix, covariancematrix, basis)\n\nEvaluate the q estimator:\n\nq_a = v^* C^-1 B C_a B^* C^-1 v\n\nArguments:\n\nmmodes or v specifies the list of measured m-modes\ntransfermatrix or B specifies the interferometer\'s response to the sky\ncovariancematrix or C specifies the covariance of the measured m-modes\nbasis or C_a is a list of angular covariance matrices that represent the change in the covariance with respect to an increase in power of each 21-cm power spectrum bin\n\n\n\n"
},

{
    "location": "power-spectrum-estimation/#BPJSpec.fisher_information",
    "page": "Power Spectrum Estimation",
    "title": "BPJSpec.fisher_information",
    "category": "function",
    "text": "fisher_information(transfermatrix, covariancematrix, basis; iterations=10)\n\nCompute a Monte-Carlo approximation of the Fisher information matrix.\n\nF_ab = rm trleft( C^-1 C_a C^-1 C_b right)\n\nArguments:\n\ntransfermatrix or B specifies the interferometer\'s response to the sky\ncovariancematrix or C specifies the covariance of the measured m-modes\nbasis or C_a is a list of angular covariance matrices that represent the change in the covariance with respect to an increase in power of each 21-cm power spectrum bin\n\nKeyword Arguments:\n\niterations is the number of Monte Carlo simulations to perform\n\n\n\n"
},

{
    "location": "power-spectrum-estimation/#BPJSpec.noise_bias",
    "page": "Power Spectrum Estimation",
    "title": "BPJSpec.noise_bias",
    "category": "function",
    "text": "noise_bias(transfermatrix, covariancematrix, basis; iterations=10)\n\nCompute a Monte-Carlo approximation of the noise bias to the quadratic estimator.\n\n\n\n"
},

{
    "location": "power-spectrum-estimation/#BPJSpec.full_rank_compress!",
    "page": "Power Spectrum Estimation",
    "title": "BPJSpec.full_rank_compress!",
    "category": "function",
    "text": "full_rank_compress!(output_mmodes, output_transfermatrix, output_noisematrix,\n                    input_mmodes,  input_transfermatrix,  input_noisematrix;\n                    progress=false)\n\nIn the case where the interferometer has more baselines than there are spherical harmonic coefficients to measure, the transfer matrix is tall and skinny. This also indicates that we have made redundant measurements that can be averaged together with no information loss.\n\nIn this routine we use the singular value decomposition (SVD) of the transfer matrix to compress the measurements. However, the SVD is just as large as the transfer matrix itself, and will take a lot of disk space to store. Therefore we will compute the SVD, compress everything with it all at once so that there is no need to store the SVD as well.\n\nArguments:\n\noutput_mmodes the output compressed m-modes\noutput_transfermatrix the output compressed transfer matrix\noutput_noisematrix the output compressed noise covariance matrix\ninput_mmodes the input m-modes that will be compressed\ninput_transfermatrix the input transfer matrix that will be used to generate the compression\ninput_noisematrix the input noise covariance matrix\n\nKeyword Arguments:\n\nprogress if set to true, a progress bar will be displayed\n\n\n\n"
},

{
    "location": "power-spectrum-estimation/#API-1",
    "page": "Power Spectrum Estimation",
    "title": "API",
    "category": "section",
    "text": "CurrentModule = BPJSpec\nDocTestSetup = quote\n    using BPJSpec\nendq_estimator\nfisher_information\nnoise_bias\nfull_rank_compress!"
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
    "text": "Fundamentally, m-mode analysis involves operations on block-diagonal matrices. Although the block-diagonal structure of these matrices makes it possible to store and perform computations on these matrices, they can still be enormous. That is, while each block individually may fit in the system memory, multiple blocks may not. Here we introduce some of the abstractions developed for working with block diagonal matrices that may each be multiple terabytes in size."
},

{
    "location": "enormous-matrices/#BPJSpec.AbstractBlockMatrix",
    "page": "Enormous Matrices",
    "title": "BPJSpec.AbstractBlockMatrix",
    "category": "type",
    "text": "abstract type AbstractBlockMatrix{B, N}\n\nThis type represents a (potentially enormous) block-diagonal matrix. This type is designed to be general enough to handle large matrices that fit in memory as well as enormous matrices that do not fit in memory. In principle this type can also be used to store small matrices, but it would be relatively inefficient compared to the standard Array{T, N}.\n\nType Parameters\n\nB specifies the type of the blocks that compose the matrix\nN specifies the number of indices used to index into the matrix blocks\n\nRequired Fields\n\nstorage contains instructions on how to read matrix blocks\ncache is used if we want to read the matrix from disk and then keep it in memory for faster   access.\n\n\n\n"
},

{
    "location": "enormous-matrices/#BPJSpec.SimpleBlockVector",
    "page": "Enormous Matrices",
    "title": "BPJSpec.SimpleBlockVector",
    "category": "type",
    "text": "struct SimpleBlockVector <: AbstractBlockMatrix{Vector{Complex128}, 1}\n\nThis type represents a (potentially enormous) complex-valued vector that has been split into blocks. Each of these blocks is indexed by a number that varies from 1 to length.\n\nFields:\n\nstorage contains instructions on how to read the vector from disk\ncache is used if we want to keep the vector in memory\nlength determines the number of blocks the vector is divided into\n\nUsage:\n\njulia> x = create(SimpleBlockVector, 10)\nSimpleBlockVector(<no file>, cached=true, length=10)\n\njulia> x[5] = Complex128[1, 2, 3, 4, 5];\n\njulia> x[5]\n5-element Array{Complex{Float64},1}:\n 1.0+0.0im\n 2.0+0.0im\n 3.0+0.0im\n 4.0+0.0im\n 5.0+0.0im\n\nSee also: SimpleBlockMatrix, AbstractBlockMatrix\n\n\n\n"
},

{
    "location": "enormous-matrices/#BPJSpec.SimpleBlockMatrix",
    "page": "Enormous Matrices",
    "title": "BPJSpec.SimpleBlockMatrix",
    "category": "type",
    "text": "struct SimpleBlockMatrix <: AbstractBlockMatrix{Matrix{Complex128}, 1}\n\nThis type represents a (potentially enormous) complex-valued matrix that has been split into blocks. Each of these blocks is indexed by a number that varies from 1 to length.\n\nFields:\n\nstorage contains instructions on how to read the matrix from disk\ncache is used if we want to keep the matrix in memory\nlength determines the number of blocks the matrix is divided into\n\nUsage:\n\njulia> x = create(SimpleBlockMatrix, 10)\nSimpleBlockMatrix(<no file>, cached=true, length=10)\n\njulia> x[5] = Complex128[1 2; 3 4];\n\njulia> x[5]\n2×2 Array{Complex{Float64},2}:\n 1.0+0.0im  2.0+0.0im\n 3.0+0.0im  4.0+0.0im\n\nSee also: SimpleBlockVector, AbstractBlockMatrix\n\n\n\n"
},

{
    "location": "enormous-matrices/#BPJSpec.MBlockVector",
    "page": "Enormous Matrices",
    "title": "BPJSpec.MBlockVector",
    "category": "type",
    "text": "struct MBlockVector <: AbstractBlockMatrix{Vector{Complex128}, 1}\n\nThis type represents a (potentially enormous) complex-valued vector that has been split into blocks. Each of these blocks is indexed by its value of m that varies from 0 to mmax.\n\nFields:\n\nstorage contains instructions on how to read the matrix from disk\ncache is used if we want to keep the matrix in memory\nmmax determines the largest value of the m quantum number used by the matrix\n\nUsage:\n\njulia> x = create(MBlockVector, 10)\nMBlockVector(<no file>, cached=true, mmax=10)\n\njulia> x[0] = Complex128[1, 2, 3, 4, 5];\n\njulia> x[0]\n5-element Array{Complex{Float64},1}:\n 1.0+0.0im\n 2.0+0.0im\n 3.0+0.0im\n 4.0+0.0im\n 5.0+0.0im\n\nSee also: MBlockMatrix, AbstractBlockMatrix\n\n\n\n"
},

{
    "location": "enormous-matrices/#BPJSpec.MBlockMatrix",
    "page": "Enormous Matrices",
    "title": "BPJSpec.MBlockMatrix",
    "category": "type",
    "text": "struct MBlockMatrix <: AbstractBlockMatrix{Matrix{Complex128}, 1}\n\nThis type represents a (potentially enormous) complex-valued matrix that has been split into blocks. Each of these blocks is indexed by its value of m that varies from 0 to mmax.\n\nFields:\n\nstorage contains instructions on how to read the matrix from disk\ncache is used if we want to keep the matrix in memory\nmmax determines the largest value of the m quantum number used by the matrix\n\nUsage:\n\njulia> x = create(MBlockMatrix, 10)\nMBlockMatrix(<no file>, cached=true, mmax=10)\n\njulia> x[0] = Complex128[1 2; 3 4];\n\njulia> x[0]\n2×2 Array{Complex{Float64},2}:\n 1.0+0.0im  2.0+0.0im\n 3.0+0.0im  4.0+0.0im\n\nSee also: MBlockVector, AbstractBlockMatrix\n\n\n\n"
},

{
    "location": "enormous-matrices/#BPJSpec.FBlockVector",
    "page": "Enormous Matrices",
    "title": "BPJSpec.FBlockVector",
    "category": "type",
    "text": "struct FBlockVector <: AbstractBlockMatrix{Vector{Complex128}, 1}\n\nThis type represents a (potentially enormous) complex-valued vector that has been split into blocks. Each of these blocks is indexed by the index of the corresponding frequency channel, which varies from 1 to length(frequencies).\n\nFields:\n\nstorage contains instructions on how to read the matrix from disk\ncache is used if we want to keep the matrix in memory\nfrequencies is a list of the frequency channels represented by this matrix\nbandwidth is a list of the corresponding bandwidth of each frequency channel\n\nUsage:\n\njulia> x = create(FBlockVector, [74u\"MHz\", 100u\"MHz\"], [24u\"kHz\", 24u\"kHz\"])\nFBlockVector(<no file>, cached=true, frequencies=74.000 MHz…100.000 MHz, bandwidth~24 kHz)\n\njulia> x[1] = Complex128[1, 2, 3, 4, 5];\n\njulia> x[1]\n5-element Array{Complex{Float64},1}:\n 1.0+0.0im\n 2.0+0.0im\n 3.0+0.0im\n 4.0+0.0im\n 5.0+0.0im\n\nSee also: FBlockMatrix, AbstractBlockMatrix\n\n\n\n"
},

{
    "location": "enormous-matrices/#BPJSpec.FBlockMatrix",
    "page": "Enormous Matrices",
    "title": "BPJSpec.FBlockMatrix",
    "category": "type",
    "text": "struct FBlockMatrix <: AbstractBlockMatrix{Matrix{Complex128}, 1}\n\nThis type represents a (potentially enormous) complex-valued matrix that has been split into blocks. Each of these blocks is indexed by the index of the corresponding frequency channel, which varies from 1 to length(frequencies).\n\nFields:\n\nstorage contains instructions on how to read the matrix from disk\ncache is used if we want to keep the matrix in memory\nfrequencies is a list of the frequency channels represented by this matrix\nbandwidth is a list of the corresponding bandwidth of each frequency channel\n\nUsage:\n\njulia> x = create(FBlockMatrix, [74u\"MHz\", 100u\"MHz\"], [24u\"kHz\", 24u\"kHz\"])\nFBlockMatrix(<no file>, cached=true, frequencies=74.000 MHz…100.000 MHz, bandwidth~24 kHz)\n\njulia> x[1] = Complex128[1 2; 3 4];\n\njulia> x[1]\n2×2 Array{Complex{Float64},2}:\n 1.0+0.0im  2.0+0.0im\n 3.0+0.0im  4.0+0.0im\n\nSee also: FBlockVector, AbstractBlockMatrix\n\n\n\n"
},

{
    "location": "enormous-matrices/#BPJSpec.MFBlockVector",
    "page": "Enormous Matrices",
    "title": "BPJSpec.MFBlockVector",
    "category": "type",
    "text": "struct MFBlockVector <: AbstractBlockMatrix{Vector{Complex128}, 2}\n\nThis type represents a (potentially enormous) complex-valued vector that has been split into blocks. Each of these blocks is indexed by its value of m, which varies from 0 to mmax, and the index of the corresponding frequency channel, which varies from 1 to length(frequencies).\n\nFields:\n\nstorage contains instructions on how to read the matrix from disk\ncache is used if we want to keep the matrix in memory\nmmax determines the largest value of the m quantum number used by the matrix\nfrequencies is a list of the frequency channels represented by this matrix\nbandwidth is a list of the corresponding bandwidth of each frequency channel\n\nUsage:\n\njulia> x = create(MFBlockVector, 2, [74u\"MHz\", 100u\"MHz\"], [24u\"kHz\", 24u\"kHz\"])\nMFBlockVector(<no file>, cached=true, mmax=2, frequencies=74.000 MHz…100.000 MHz, bandwidth~24 kHz)\n\njulia> x[0, 1] = Complex128[1, 2, 3, 4, 5];\n\njulia> x[0, 1]\n5-element Array{Complex{Float64},1}:\n 1.0+0.0im\n 2.0+0.0im\n 3.0+0.0im\n 4.0+0.0im\n 5.0+0.0im\n\nSee also: MFBlockMatrix, AbstractBlockMatrix\n\n\n\n"
},

{
    "location": "enormous-matrices/#BPJSpec.MFBlockMatrix",
    "page": "Enormous Matrices",
    "title": "BPJSpec.MFBlockMatrix",
    "category": "type",
    "text": "struct MFBlockMatrix <: AbstractBlockMatrix{Matrix{Complex128}, 2}\n\nThis type represents a (potentially enormous) complex-valued matrix that has been split into blocks. Each of these blocks is indexed by its value of m, which varies from 0 to mmax, and the index of the corresponding frequency channel, which varies from 1 to length(frequencies).\n\nFields:\n\nstorage contains instructions on how to read the matrix from disk\ncache is used if we want to keep the matrix in memory\nmmax determines the largest value of the m quantum number used by the matrix\nfrequencies is a list of the frequency channels represented by this matrix\nbandwidth is a list of the corresponding bandwidth of each frequency channel\n\nUsage:\n\njulia> x = create(MFBlockMatrix, 2, [74u\"MHz\", 100u\"MHz\"], [24u\"kHz\", 24u\"kHz\"])\nMFBlockMatrix(<no file>, cached=true, mmax=2, frequencies=74.000 MHz…100.000 MHz, bandwidth~24 kHz)\n\njulia> x[0, 1] = Complex128[1 2; 3 4];\n\njulia> x[0, 1]\n2×2 Array{Complex{Float64},2}:\n 1.0+0.0im  2.0+0.0im\n 3.0+0.0im  4.0+0.0im\n\nSee also: MFBlockVector, AbstractBlockMatrix\n\n\n\n"
},

{
    "location": "enormous-matrices/#BPJSpec.LBlockMatrix",
    "page": "Enormous Matrices",
    "title": "BPJSpec.LBlockMatrix",
    "category": "type",
    "text": "struct LBlockMatrix <: AbstractBlockMatrix{Matrix{Float64}, 1}\n\nThis type represents a (potentially enormous) complex-valued matrix that has been split into blocks. Each of these blocks is indexed by its value of l, which varies from 0 to lmax.\n\nFields:\n\nstorage contains instructions on how to read the matrix from disk\ncache is used if we want to keep the matrix in memory\nlmax determines the largest value of the l quantum number used by the matrix\nfrequencies is a list of the frequency channels represented by this matrix\nbandwidth is a list of the corresponding bandwidth of each frequency channel\n\nUsage:\n\njulia> x = create(LBlockMatrix, 2, [74u\"MHz\", 100u\"MHz\"], [24u\"kHz\", 24u\"kHz\"])\nLBlockMatrix(<no file>, cached=true, lmax=2, frequencies=74.000 MHz…100.000 MHz, bandwidth~24 kHz)\n\njulia> l = BPJSpec.L(0);\n\njulia> x[l] = Float64[1 2; 3 4];\n\njulia> x[l]\n2×2 Array{Float64,2}:\n 1.0  2.0\n 3.0  4.0\n\nSee also: LMBlockVector, AbstractBlockMatrix\n\n\n\n"
},

{
    "location": "enormous-matrices/#BPJSpec.LMBlockVector",
    "page": "Enormous Matrices",
    "title": "BPJSpec.LMBlockVector",
    "category": "type",
    "text": "struct LMBlockVector <: AbstractBlockMatrix{Vector{Complex128}, 2}\n\nThis type represents a (potentially enormous) complex-valued vector that has been split into blocks. Each of these blocks is indexed by its value of l, which varies from 0 to lmax, and m, which varies from 0 to mmax with the restriction that m  l.\n\nFields:\n\nstorage contains instructions on how to read the matrix from disk\ncache is used if we want to keep the matrix in memory\nlmax determines the largest value of the l quantum number used by the matrix\nmmax determines the largest value of the m quantum number used by the matrix\nfrequencies is a list of the frequency channels represented by this matrix\nbandwidth is a list of the corresponding bandwidth of each frequency channel\n\nUsage:\n\njulia> x = create(LMBlockVector, 2, 2, [74u\"MHz\", 100u\"MHz\"], [24u\"kHz\", 24u\"kHz\"])\nLMBlockVector(<no file>, cached=true, lmax=2, mmax=2, frequencies=74.000 MHz…100.000 MHz, bandwidth~24 kHz)\n\njulia> x[0, 0] = Complex128[1, 2, 3, 4, 5];\n\njulia> x[0, 0]\n5-element Array{Complex{Float64},1}:\n 1.0+0.0im\n 2.0+0.0im\n 3.0+0.0im\n 4.0+0.0im\n 5.0+0.0im\n\nSee also: LBlockMatrix, AbstractBlockMatrix\n\n\n\n"
},

{
    "location": "enormous-matrices/#API-1",
    "page": "Enormous Matrices",
    "title": "API",
    "category": "section",
    "text": "CurrentModule = BPJSpec\nDocTestSetup = quote\n    using BPJSpec\n    using Unitful\nendAbstractBlockMatrix\nSimpleBlockVector\nSimpleBlockMatrix\nMBlockVector\nMBlockMatrix\nFBlockVector\nFBlockMatrix\nMFBlockVector\nMFBlockMatrix\nLBlockMatrix\nLMBlockVector"
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
    "location": "wrappers/#BPJSpec.GSLWrapper.Y",
    "page": "Wrappers",
    "title": "BPJSpec.GSLWrapper.Y",
    "category": "function",
    "text": "Y(l, m, θ, ϕ)\n\nThe spherical harmonic function (using the Condon-Shortley phase convention).\n\nUsage:\n\njulia> using BPJSpec\n\njulia> BPJSpec.Y(0, 0, 1, 2)\n0.28209479177387814 + 0.0im\n\njulia> 1/sqrt(4π)\n0.28209479177387814\n\n\n\n"
},

{
    "location": "wrappers/#GSLWrapper.jl-1",
    "page": "Wrappers",
    "title": "GSLWrapper.jl",
    "category": "section",
    "text": "CurrentModule = BPJSpec\nDocTestSetup = quote\n    using BPJSpec.GSLWrapper\nendY"
},

{
    "location": "wrappers/#FastTransformsWrapper.jl-1",
    "page": "Wrappers",
    "title": "FastTransformsWrapper.jl",
    "category": "section",
    "text": ""
},

]}
