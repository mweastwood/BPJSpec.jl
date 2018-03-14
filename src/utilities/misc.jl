"Convenience method for computing the conjugate transpose of a matrix."
T(A) = ctranspose(A)

"Convenience function for guaranteeing a matrix is Hermitian."
H(A) = 0.5*(A+A')

"Try to make sure a matrix that should be positive definite is in fact positive definite."
function fix(A)
    N = size(A, 1)
    N == 0 && return A
    B = H(A)
    λ = eigvals(B)
    λmin = minimum(λ)
    λmax = maximum(λ)
    if λmin ≤ 0
        factor = N * eps(Float64) * λmax
        return B + factor * I
    else
        return B
    end
end

"Useful little function that helps account for grouping of positive and negative m."
two(m) = ifelse(m != 0, 2, 1)

"Convert and strip units from the given quantity."
u(units, quantity) = ustrip(uconvert(units, quantity))

@enum Mode r a w
const mode = Dict(r => (false, false, false), a => (true, true, false), w => (true, true, true))

