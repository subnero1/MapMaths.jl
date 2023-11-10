using MapMaths
using Documenter

DocMeta.setdocmeta!(
    MapMaths,
    :DocTestSetup,
    :(using MapMaths),
)

tmp = mktempdir()
open(joinpath(tmp, "README.md"); write = true) do test_readme
    open("../README.md"; read = true) do readme
        println(
            test_readme,
            """
            ```jldoctest MapMaths
            julia> using MapMaths, PrettyTables
            ```
            """
        )
        write(test_readme, readme)
    end
end
doctest(tmp, [MapMaths])