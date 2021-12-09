"""
    struct GlobalMatrices
Stores global matrices and some other things.
Everything that is very expensive to compute which stays
constant throughout the simulation should be stored here.
    GlobalMatrices(basis::Basis, filter::Filter, dimensions)
Initialize GlobalMatrices for `basis`, `filter` and `dimensions`.
"""
 struct GlobalMatrices
     normalsigns::Dict{Face, Int64}
     normalidxs::Dict{Face, Int64}
     oppositefaces::Dict{Face, Face}

     function GlobalMatrices(dimensions)
        normalsigns = Dict(left => -1, right => 1, top => 1, bottom => -1)
        normalidxs = Dict(left => 1, right => 1, top => 2, bottom => 2)
        oppositefaces = Dict(left => right, right => left, top => bottom, bottom => top)

         new(normalsigns,
            normalidxs,
            oppositefaces)
     end
 end