using BPJSpec
using Base.Test

srand(123)

# Verify that the visibilities -> m-modes -> visibilities
# round-trip works.
let Nant = 3, mmax = 5
    Nbase = div(Nant*(Nant-1),2)
    data = zeros(Complex128,Nbase,2mmax+1)
    rand!(data)
    mmodes = MModes(data,mmax=mmax)
    data′ = visibilities(mmodes)
    @test_approx_eq data data′
end

# Check that the gridded data I/O is lossless.
let data = rand(Complex64,12,34), weights = rand(Float64,34)
    filename = tempname()*".h5"
    BPJSpec.write_data(filename,data,weights)
    newdata,newweights = BPJSpec.read_data(filename)
    @test data == newdata
    @test weights == newweights
end

# Check that the m-mode I/O is lossless.
let v = TransferMatrix(5,10)
    for m = 0:mmax(v)
        rand!(BPJSpec.block(v,m))
    end
    filename = tempname()*".h5"
    BPJSpec.write_mmodes(filename,v)
    newv = BPJSpec.read_mmodes(filename)
    @test v == newv
end

# Check that the transfer matrix I/O is lossless.
let B = TransferMatrix(5,10,10)
    for m = 0:mmax(B)
        rand!(BPJSpec.block(B,m))
    end
    filename = tempname()*".h5"
    BPJSpec.write_transfermatrix(filename,B)
    newB = BPJSpec.read_transfermatrix(filename)
    @test B == newB
end

