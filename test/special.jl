# Test the spherical harmonics
for i = 1:10
    θ =  π*rand()
    ϕ = 2π*rand()
    @test BPJSpec.Y(0,0,θ,ϕ) ≈ 1/2 * sqrt(1/π)
    @test BPJSpec.Y(1,0,θ,ϕ) ≈ 1/2 * sqrt(3/π) * cos(θ)
    @test BPJSpec.Y(1,+1,θ,ϕ) ≈ -1/2 * sqrt(3/(2π)) * sin(θ) * exp(+1im*ϕ)
    @test BPJSpec.Y(1,-1,θ,ϕ) ≈ +1/2 * sqrt(3/(2π)) * sin(θ) * exp(-1im*ϕ)
end

# Test the spherical Bessel function
for i = 1:10
    x = 100*randn() |> abs
    @test BPJSpec.j(0,x) ≈ sin(x)/x
    @test BPJSpec.j(1,x) ≈ sin(x)/x^2 - cos(x)/x
    @test BPJSpec.j(2,x) ≈ (3/x^2 - 1)*sin(x)/x - 3cos(x)/x^2
end

# Test the plane wave expansion
let u = 2.0, v = -3.0, w = 2.0
    alm_real, alm_imag = BPJSpec.planewave(u,v,w)
    map_real = alm2map(alm_real,nside=512)
    map_imag = alm2map(alm_imag,nside=512)
    map_test_real = HealpixMap(Float64,512)
    map_test_imag = HealpixMap(Float64,512)
    for i = 1:nside2npix(512)
        vec = LibHealpix.pix2vec_ring(512,i)
        phase = 2π*(u*vec[1]+v*vec[2]+w*vec[3])
        map_test_real[i] = cos(phase)
        map_test_imag[i] = sin(phase)
    end
    @test pixels(map_real) ≈ pixels(map_test_real)
    @test pixels(map_imag) ≈ pixels(map_test_imag)
end

# Test the matrix trace
for i = 1:10
    A = rand(Float64,50,50)
    B = rand(Float64,50,50)
    @test trace(A*B) ≈ BPJSpec.tr(A,B)
    A = rand(Complex128,50,50)
    B = rand(Complex128,50,50)
    @test trace(A*B) ≈ BPJSpec.tr(A,B)
end

# Test force_hermitian and force_posdef
let
    A = rand(10,10)
    @test ishermitian(BPJSpec.force_hermitian(A))
    B = diagm(linspace(1,-1e-11,10))
    @test isposdef(BPJSpec.force_posdef(B))
end

