using InverseLaplace
using Aqua: Aqua

@testset "aqua deps compat" begin
    Aqua.test_deps_compat(InverseLaplace)
end

# This often gives false positive
@testset "aqua project toml formatting" begin
    Aqua.test_project_toml_formatting(InverseLaplace)
end

@testset "aqua unbound_args" begin
    Aqua.test_unbound_args(InverseLaplace)
end

@testset "aqua undefined exports" begin
    Aqua.test_undefined_exports(InverseLaplace)
end

# Depending on Optim causes many ambiguity errors outside our control
# @testset "aqua test ambiguities" begin
#     Aqua.test_ambiguities([InverseLaplace, Core, Base])
# end

@testset "aqua piracy" begin
    Aqua.test_piracy(InverseLaplace)
end

@testset "aqua project extras" begin
    Aqua.test_project_extras(InverseLaplace)
end

@testset "aqua state deps" begin
    Aqua.test_stale_deps(InverseLaplace)
end
