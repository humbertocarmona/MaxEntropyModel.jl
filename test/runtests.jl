using MaxEntropyModel
using Test
using LinearAlgebra

function testmodel()
    M = map(x -> x < 0.5 ? -1 : 1, rand(Int64, 1000, 20))
    model = MaxEntropyModel.MaxEnt(M)
    @test model.model == "MaxEnt"
    @test model.runid == "test"
    @test model.nspins == 20
    @test model.S_obs == M
end
@testset "MaxEnt" testmodel()


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

    δ = MaxEntropyModel.deltaEnergy(model, 1)
    model.s[1] = -model.s[1]
    fn = MaxEntropyModel.energy(model)
    δ1 = fn - en
    @test δ1 == δ

    δ2 = MaxEntropyModel.deltaEnergy(model, 20)
    model.s[20] = -model.s[20]
    en = MaxEntropyModel.energy(model)
    δ3 = en - fn
    @test δ2 == δ3
end
@testset "energy" testenergy()

