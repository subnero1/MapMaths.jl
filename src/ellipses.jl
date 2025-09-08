abstract type Ellipse end

abstract type EllipseParameter end
macro ellipse_parameter(Type::Symbol)
    return esc(quote
        struct $Type{T <: Number} <: EllipseParameter
            param::T
        end
        Base.iterate(p::$Type, state...) = iterate((p.param,), state...)
    end)
end

Base.getindex(e::Ellipse, P::Type{<:EllipseParameter}) = only(P(e))

@ellipse_parameter(SemiMajorAxis)
@ellipse_parameter(SemiMinorAxis)
@ellipse_parameter(Flattening)
@ellipse_parameter(InverseFlattening)
@ellipse_parameter(Eccentricity)
@ellipse_parameter(SquaredEccentricity)

struct EllipseAFi{T <: Number} <: Ellipse
    a::T
    fi::T
end
@symmetric Base.:|((a,)::SemiMajorAxis, (fi,)::InverseFlattening) = EllipseAFi(a, fi)
Base.iterate(e::EllipseAFi, state...) = iterate((e.a, e.fi), state...)
Base.show(io::IO, (a, fi)::EllipseAFi) = print(io, "Ellipse(", SemiMajorAxis(a), " | ", InverseFlattening(fi), ")")
(P::Type{<:SemiMajorAxis})((a, fi)::EllipseAFi) = P(a)
(P::Type{<:SemiMinorAxis})((a, fi)::EllipseAFi) = P(a * (1 - inv(fi)))
(P::Type{<:Flattening})((a, fi)::EllipseAFi) = P(inv(fi))
(P::Type{<:InverseFlattening})((a, fi)::EllipseAFi) = P(fi)
(P::Type{<:SquaredEccentricity})((a, fi)::EllipseAFi) = P((2 - inv(fi)) / fi)
(P::Type{<:Eccentricity})(e::EllipseAFi) = P(sqrt(e[SquaredEccentricity]))

"""
    normal(e::Ellipse, p) -> n

Compute a unit vector `n` such that the ray `p + [0, Inf) * n` intersects the
ellipse at a right angle. 

For any given `e` and `p`, there are between two to four such vectors.
`normal(e, p)` returns one of the vectors with the closest intersection.
"""
function normal(ellipse::Ellipse, (x, y))
    # Implements the SC Halley algorithm from
    # https://www.researchgate.net/profile/Toshio-Fukushima/publication/264118552_F06-gconv2/data/53ce803e0cf279d93530b524/F06-gconv2.pdf
    a = e[SemiMajorAxis]
    b = e[SemiMinorAxis]
    E = e[SquaredEccentricity]
    ec = b / a

    P = p / a
    Z = ec / a * abs(z)
    S = Z
    C = ec^2 * P # Paper says ec * P, but ec^2 * P is correct

    A = hypot(S, C)
    B = 1.5 * E^2 * P * (S * C)^2 * (A - ec)
    D = Z * A^3 + E * S^3
    F = P * A^3 - E * C^3
    S = D * F - B * S
    C = F^2 - B * C
    return ifelse(p > 0, normalize(SVector{2, Float64}(ec * C, copysign(S, z))), SVector{2, Float64}(0, copysign(1, z)))
end