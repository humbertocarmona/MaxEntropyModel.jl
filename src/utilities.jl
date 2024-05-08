function mean_1st_order_moments(M::Matrix{<:Number})
    return mean(M, dims=1)[:]
end

function mean_2nd_order_moments(M::Matrix{<:Number})
    MM = M' * M
    MM /= size(M, 1)

    n = size(M, 2)
    XY = Array{Float64}(undef, div(n * (n - 1), 2))
    t = 1
    @simd for i in 1:n-1
        @simd for j in i+1:n
            XY[t] = MM[i, j]
            t += 1
        end
    end

    return XY
end

function mean_3rd_order_moments(M::Matrix{<:Number})
    nrows, N = size(M) # samples, spins

    # Preallocate XYZ vector
    XYZ = zeros(Float64, div(N * (N - 1) * (N - 2), 6))

    # Iterate over each row of M
    @inbounds for row in eachrow(M)
        x = view(row, :) # Use view for efficiency
        t = 1
        @simd for i in 1:N-2
            @simd for j in i+1:N-1
                @simd for k in j+1:N
                    XYZ[t] += x[i] * x[j] * x[k]
                    t += 1
                end
            end
        end
    end

    # Normalize by the number of rows
    XYZ ./= nrows

    return XYZ
end

"""
	spin_permutations_iterator(n; spin_values=[+1, -1])

	Generate an iterator for all possible spin permutations of length `n` with specified spin values.

	# Arguments:
	- `n::Int`: Length of the spin permutation.
	- `spin_values::Vector{Int}`: Values representing the possible spins. Default is [+1, -1].

	# Returns:
	- `ProductIterator`: Iterator over all possible spin permutations.

	The function creates an iterator that generates all possible spin permutations of length `n` by combining the specified spin values.

	# Examples
	```julia
		iterator = spin_permutations_iterator(3)
		for spin_permutation in iterator
			println(spin_permutation)
		end
	```
"""
function spin_permutations_iterator(n; spin_values=[+1, -1])
    return Base.Iterators.product(fill(spin_values, n)...)
end

function ones_distribution(M::Matrix{Int64})
    nsamples, N = size(M)

    # Preallocate K array
    K = Vector{Int64}(undef, nsamples)

    # Iterate over each row of M
    @inbounds for i in 1:nsamples
        K[i] = count(isequal(1), M[i, :])
    end

    # Generate bin edges
    edges = collect(-0.5:1:N+0.5)

    # Compute histogram directly
    hist_counts = zeros(Int64, length(edges) - 1)
    for k in K
        if k >= 0 && k <= N
            hist_counts[k+1] += 1  # Adjusting index to align with bin edges
        end
    end

    # Normalize histogram
    p = hist_counts / nsamples

    # Generate bin centers
    centers = collect(0:N)

    return centers, p
end

function make_bonds(N::Int64)
    bond = zeros(Int64, N, N)
    t = 1
    for i in 1:N-1, j in i+1:N
        bond[i, j] = bond[j, i] = t
        t += 1
    end
    return bond
end
