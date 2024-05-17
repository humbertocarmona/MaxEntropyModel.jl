function mean_1st_order_moments(M::Matrix{<:Number})
    return mean(M, dims=1)[:]
end

function straighten(M::Matrix{<:Number})
    n = size(M, 2)
    XY = Array{Float64}(undef, div(n * (n - 1), 2))
    t = 1
    @simd for i in 1:n-1
        @simd for j in i+1:n
            XY[t] = M[i, j]
            t += 1
        end
    end
    return XY
end

function crunk(arr::Vector{<:Number})
    n = round(Int64, (1 + sqrt(1 + 8 * length(arr))) / 2)
    @assert length(arr) == n * (n - 1) ÷ 2 "wrong size array"
    M = zeros(n, n)
    t = 1
    for i in 1:n-1
        for j in i+1:n
            M[i, j] = arr[t]
            M[j, i] = M[i, j]
            t += 1
        end
    end
    return M
end

function mean_2nd_order_moments(M::Matrix{<:Number})
    MM = M' * M
    MM /= size(M, 1)

    return straighten(MM)
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

function spin_permutations_iterator2(n; spin_values=[+1, -1])
    return Base.Iterators.map(x -> [spin_values[bit+1] for bit in x], Base.Iterators.product(fill(0:1, n)...))
end

function gray_code_iterator(n; spin_values=[+1, -1])
    # Function to generate the nth Gray code number
    gray_code(k) = k ⊻ (k >>> 1)

    # Iterator to generate the Gray code sequence
    return Base.Iterators.map(
        k -> [spin_values[(gray_code(k)>>>i)&1+1] for i in (n-1):-1:0],
        0:(2^n-1)
    )
end
#==
n=1
m1.sj .= ones(Int,nspins)
en = energy(m1)
for s in MaxEntropyModel.gray_code_iterator(nspins)
	sj = collect(s)
	i = findfirst(x->x!=0, sj - m1.sj)
	if ~isnothing(i)
		en += deltaEnergy(m1,i)
		m1.sj .= sj
	end
	hj = energy(m1)
	println("$i, $(en), $(hj), $(hj -en)")
	println(sold, "->", m1.sj)
	n+=1
	(n>=5) && break;
end

==#

function ones_distribution(M::Matrix{<:Number})
    nsamples, N = size(M)

    # Preallocate K array
    K = Vector{Int64}(undef, nsamples)

    # Iterate over each row of M
    @inbounds for i in 1:nsamples
        K[i] = count(x -> x > 0, M[i, :])
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

function map_to_unit_interval(x, xmin, xmax)
    if xmin > xmax
        xmax, xmin = xmin, xmax
    end

    # slope
    m = 2 / (xmax - xmin)
    # y-intercept
    c = -(xmin + xmax) / (xmax - xmin)

    # linear transformation
    y = m * x .+ c

    # Ensure the output stays within the bounds of -1 to 1
    return clamp(y, -1, 1)
end

function centered_moments_obs(m::MaxEnt)

    covariance = zeros(size(m.xy_obs))
    pearson = zeros(size(m.xy_obs))
    xy_matrix = zeros(m.nspins, m.nspins)
    triplets = zeros(size(m.xyz_obs))

    stdev = sqrt.(1.0 .- m.x_obs .^ 2)
    t = 1
    for i in 1:m.nspins-1
        for j in i+1:m.nspins
            covariance[t] = m.xy_obs[t] - m.x_obs[i] * m.x_obs[j]
            pearson[t] = covariance[t] / (stdev[i] * stdev[j])
            xy_matrix[i, j] = m.xy_obs[t]
            xy_matrix[j, i] = xy_matrix[i, j]
            t += 1
        end
    end

    t = 1
    for i in 1:m.nspins-2
        for j in i+1:m.nspins-1
            for k in j+1:m.nspins
                triplets[t] = (m.xyz_obs[t] - xy_matrix[i, j] * m.x_obs[k] -
                               xy_matrix[i, k] * m.x_obs[j] -
                               xy_matrix[j, k] * m.x_obs[i] +
                               2 * m.x_obs[i] * m.x_obs[j] * m.x_obs[k])
                t += 1
            end
        end
    end

    return covariance, pearson, triplets

end

function centered_moments_mod(m::MaxEnt)

    covariance = zeros(size(m.xy_mod))
    pearson = zeros(size(m.xy_mod))
    xy_matrix = zeros(m.nspins, m.nspins)
    triplets = zeros(size(m.xyz_mod))

    stdev = sqrt.(1.0 .- m.x_mod .^ 2)
    t = 1
    for i in 1:m.nspins-1
        for j in i+1:m.nspins
            covariance[t] = m.xy_mod[t] - m.x_mod[i] * m.x_mod[j]
            pearson[t] = covariance[t] / (stdev[i] * stdev[j])
            xy_matrix[i, j] = m.xy_mod[t]
            xy_matrix[j, i] = xy_matrix[i, j]
            t += 1
        end
    end

    t = 1
    for i in 1:m.nspins-2
        for j in i+1:m.nspins-1
            for k in j+1:m.nspins
                triplets[t] = (m.xyz_mod[t] - xy_matrix[i, j] * m.x_mod[k] -
                               xy_matrix[i, k] * m.x_mod[j] -
                               xy_matrix[j, k] * m.x_mod[i] +
                               2 * m.x_mod[i] * m.x_mod[j] * m.x_mod[k])
                t += 1
            end
        end
    end

    return covariance, pearson, triplets

end

function pearson_mod!(m::MaxEnt)

    m.pearson_mod = zeros(size(m.xy_mod))
    s = sqrt.(1.0 .- m.x_mod .^ 2)

    t = 1
    for i in 1:m.nspins-1
        for j in i+1:m.nspins
            m.pearson_mod[t] = (m.xy_mod[t] - m.x_mod[i] * m.x_mod[j]) / (s[i] * s[j])
            t += 1
        end
    end

    return nothing
end


@inline function exp_q(x::Float64, q::Float64)
    (q == 1.0) && (return exp(x))
    y = (1 + (1 - q) * x)
    (y < 0.0) && (return 0.0)
    return y^(1.0 / (1.0 - q))
end

@inline function ln_q(x::Float64, q::Float64)
    (q == 1.0) && (return log(x))
    return (x^(1 - q) - 1) / (1 - q)
end

@inline function compute_energy_shift(m::MaxEnt, q)
    σi = sign.(m.h * (1.0 - q))
    σij = sign.(m.J * (1.0 - q))
    return -dot(m.h, σi) - dot(m.J, σij)
end
