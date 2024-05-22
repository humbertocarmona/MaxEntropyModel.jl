using MaxEntropyModel
using Test
using LinearAlgebra
using Random, Distributions
using StatsBase: rmsd

function testmodel()
    M = map(x -> x < 0.5 ? -1 : 1, rand(Int64, 1000, 20))
    model = MaxEntropyModel.MaxEnt(M)
    @test model.model == "MaxEnt"
    @test model.runid == "test"
    @test model.nspins == 20
end
# @testset "MaxEnt" testmodel()


function testio()
    M = map(x -> x < 0.5 ? -1 : 1, rand(Int64, 40, 20))
    model = MaxEntropyModel.MaxEnt(M)

    MaxEntropyModel.write_model(model, "test.json", force=true)
    @test isfile("test.json")

    model2 = MaxEntropyModel.read_model("test.json")

end
# @testset "io" testio()

function testenergy()
    M = map(x -> x < 0.5 ? -1 : 1, rand(Int64, 40, 20))
    model = MaxEntropyModel.MaxEnt(M)

    s = ones(Int64, 20)
    model.h .= 1.0
    model.J .= 0.0
    en = MaxEntropyModel.energy(model,s)
    @test en == -model.nspins

    model.h .= 0.0
    model.J .= 1.0
    en = MaxEntropyModel.energy(model,s)
    @test en == -model.nspins * (model.nspins - 1) / 2

    s = map(x -> x < 0.5 ? -1 : 1, rand(model.nspins))
    model.h = rand(Normal(0, 0.1), size(model.h))
    model.J = rand(Normal(0, 0.1), size(model.J))
    en0 = MaxEntropyModel.energy(model,s)
    for i in 1:model.nspins
        δe = MaxEntropyModel.deltaEnergy(model,s, i)
        s[i] = -s[i]
        en1 = MaxEntropyModel.energy(model,s)
        @test isapprox(δe, (en1 - en0), atol=1e-14)
        en0 = en1
    end

end
# @testset "energy" testenergy()

function testfull()
    rng = Xoshiro(4323)
    M = map(x -> x < 0.5 ? -1 : 1, rand(rng, Int64, 1000, 20))
    m0 = MaxEnt(M)
    m0.h = rand(rng, Normal(0.0, 0.10), size(m0.h))
    m0.J = rand(rng, Normal(0.0, 0.10), size(m0.J))

    samples = metropolis_iteration!(m0, true)

    m1 = MaxEnt(samples)
    m1.n_relax_steps = 100
    max_entropy_relax!(m1)

    @test rmsd(m0.h, m1.h, normalize=true) < 0.015
    @test rmsd(m0.J, m1.J, normalize=true) < 0.01

end
#@testset "full_iteration" testfull()
