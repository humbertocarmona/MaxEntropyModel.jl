using MaxEntropyModel
using Test
using LinearAlgebra
using Random, Distributions

function testmodel()
    M = map(x -> x < 0.5 ? -1 : 1, rand(Int64, 1000, 20))
    model = MaxEntropyModel.MaxEnt(M)
    @test model.model == "MaxEnt"
    @test model.runid == "test"
    @test model.nspins == 20
    @test model.S_obs == M
end
#@testset "MaxEnt" testmodel()


function testio()
    M = map(x -> x < 0.5 ? -1 : 1, rand(Int64, 40, 20))
    model = MaxEntropyModel.MaxEnt(M)

    MaxEntropyModel.write_model(model, "test.json", force=true)
    @test isfile("test.json")

    model2 = MaxEntropyModel.read_model("test.json")
    @test model.S_obs == model2.S_obs

end
@testset "io" testio()

function testenergy()
    M = map(x -> x < 0.5 ? -1 : 1, rand(Int64, 40, 20))
    model = MaxEntropyModel.MaxEnt(M)

    model.h .= 1.0
    model.J .= 0.0
    model.s = ones(Int64, 20)
    en = MaxEntropyModel.energy(model)
    @test en == -model.nspins

    model.h .= 0.0
    model.J .= 1.0
    en = MaxEntropyModel.energy(model)
    @test en == -model.nspins * (model.nspins - 1) / 2

    model.h = rand(Normal(0, 0.1), size(model.h))
    model.J = rand(Normal(0, 0.1), size(model.J))
    model.s = map(x -> x < 0.5 ? -1 : 1, rand(model.nspins))
    en0 = MaxEntropyModel.energy(model)
    for i in 1:model.nspins
        δe = MaxEntropyModel.deltaEnergy(model, i)
        model.s[i] = -model.s[i]
        en1 = MaxEntropyModel.energy(model)
        @test isapprox(δe, (en1 - en0), atol=1e-14)
        en0 = en1
    end

end
@testset "energy" testenergy()

