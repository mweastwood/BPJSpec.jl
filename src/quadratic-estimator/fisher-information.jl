# Copyright (c) 2015-2017 Michael Eastwood
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

doc"""
    fisher_information(transfermatrix, covariancematrix, basis; iterations=10)

Compute a Monte-Carlo approximation of the Fisher information matrix.

```math
F_{ab} = {\rm tr}\left( C^{-1} C_a C^{-1} C_b \right)
```

**Arguments:**

* `transfermatrix` or $B$ specifies the interferometer's response to the sky
* `covariancematrix` or $C$ specifies the covariance of the measured $m$-modes
* `basis` or $C_a$ is a list of angular covariance matrices that represent the change in the
  covariance with respect to an increase in power of each 21-cm power spectrum bin

**Keyword Arguments:**

* `iterations` is the number of Monte Carlo simulations to perform
"""
function fisher_information(transfermatrix, covariancematrix, basis; iterations=10)
    μ, Σ = fisher_monte_carlo(:fisher, transfermatrix, covariancematrix, basis, iterations)
    Σ
end

"""
    noise_bias(transfermatrix, covariancematrix, basis; iterations=10)

Compute a Monte-Carlo approximation of the noise bias to the quadratic estimator.
"""
function noise_bias(transfermatrix, covariancematrix, basis; iterations=10)
    μ, Σ = fisher_monte_carlo(:bias, transfermatrix, covariancematrix, basis, iterations)
    μ
end

function fisher_monte_carlo(instruction, transfermatrix, covariancematrix, basis, iterations)
    queue = collect(1:iterations)
    lck = ReentrantLock()
    prg = Progress(length(queue))
    increment() = (lock(lck); next!(prg); unlock(lck))

    N = length(basis)
    Q = zeros(N, iterations)

    @sync for worker in workers()
        @async begin
            input  = RemoteChannel()
            output = RemoteChannel()
            remotecall(fisher_remote_processing_loop, worker, input, output,
                       transfermatrix, covariancematrix, basis)
            try
                while length(queue) > 0
                    iteration = shift!(queue)
                    put!(input, instruction)
                    Q[:, iteration] = take!(output)
                    increment()
                end
            finally
                put!(input, :quit)
            end
        end
    end

    μ = mean(Q, 2)
    Σ = zeros(N, N)
    for iteration = 1:iterations
        q = Q[:, iteration]
        Σ .+= (q.-μ)*(q.-μ)'
    end
    Σ ./= iterations - 1
    μ, Σ
end

function fisher_remote_processing_loop(input, output, transfermatrix, covariancematrix, basis)
    cache!(transfermatrix)
    cache!(covariancematrix)
    foreach(cache!, basis)
    while true
        instruction = take!(input)
        if instruction == :quit
            break
        elseif instruction == :fisher
            q = compute_fisher_q(transfermatrix, covariancematrix, basis)
            put!(output, q)
        elseif instruction == :bias
            q = compute_bias_q(transfermatrix, covariancematrix, basis)
            put!(output, q)
        else
            put!(output, :error)
        end
    end
end

function compute_fisher_q(transfermatrix, covariancematrix, basis)
    mmodes = RandomBlockVector(covariancematrix)
    q_estimator(mmodes, transfermatrix, covariancematrix, basis)
end

function compute_bias_q(transfermatrix, covariancematrix, basis)
    mmodes = WhiteNoiseBlockVector()
    q_estimator(mmodes, transfermatrix, covariancematrix, basis)
end

