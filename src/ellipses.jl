"""
    Ellipse

Supertype for ellipses. 

Notable operations: 

- New ellipse types can be defined using the [`@ellipse`](@ref) macro.
- Properties of an ellipse can be looked up using `ellipse[Parameter]` where
  Parameter <: EllipseParameter`. 
"""
abstract type Ellipse end

"""
    EllipseParameter

Supertype for ellipse parameters (e.g. `SemiMajorAxis`).
"""
abstract type EllipseParameter end

"""
    @ellipse_parameter(P)

Define a new ellipse parameter type. 
"""
macro ellipse_parameter(P)
    P = esc(P)
    return quote
        Core.@__doc__ struct $P{T <: Number} <: EllipseParameter
            value::T
            $P{T}(value::Number) where {T <: Number} = new{T}(value)
        end
        (::Type{<:$P})(value::Number) = $P{default_numtype(value)}(value)
        (C::Type{<:$P})((value,)::$P) = C(value)
        Base.convert(C::Type{<:$P}, c::$P) = C(c)

        numtype(::Type{<:$P}) = Number
        numtype(::Type{<:$P{T}}) where {T <: Number} = T

        Base.show(io::IO, p::$P) = print(io, $P, "(", p.value, ")")
    end
end

Base.getindex(ellipse::Ellipse, P::Type{<:EllipseParameter}) = P(ellipse).value

"""
    @ellipse(E, params...)

Define a new ellipse type for ellipses defined in terms of the given parameters. 

# Example
```julia
@ellipse(EllipseAB, a::SemiMajorAxis, b::SemiMinorAxis)
```
"""
macro ellipse(E, params...)
    @assert all(p -> isexpr(p, :(::)), params)
    E = esc(E)
    type_vars = ntuple(i -> Symbol(:P, i), length(params))
    vars = map(p -> esc(p.args[1]), params)
    types = map(p -> esc(p.args[2]), params)
    return quote
        Core.@__doc__ struct $E{T <: Number} <: Ellipse
            $(map(var -> :($var::T), vars)...)
        end
        @symmetric function MapMaths.Ellipse($(map((var, type) -> :($var::$type), vars, types)...))
            return Ellipse(default_numtype($(map(var -> :($var.value), vars)...)), $(vars...))
        end
        $(
            map(
                (permuted_vars, permuted_types) -> quote
                    function MapMaths.Ellipse(
                            T::Type{<:Number},
                            $(map((var, type) -> :($var::$type), permuted_vars, permuted_types)...),
                        )
                        return $E{T}($(map(var -> :($var.value), permuted_vars)...))
                    end
                end,
                permutations(vars),
                permutations(types),
            )...
        )
        MapMaths.numtype(::Type{<:$E}) = Number
        MapMaths.numtype(::Type{<:$E{T}}) where {T <: Number} = T
        function Base.show(io::IO, ellipse::$E)
            print(io, Ellipse, "(")
            join(
                io,
                ($(map((type, var) -> :($type(getfield(ellipse, $(QuoteNode(unescape(var)))))), types, vars)...),)
                , ", "
            )
            print(io, ")")
            return nothing
        end
    end
end

"""
    SemiMajorAxis

Semi-major axis of an ellipse.
"""
@ellipse_parameter(SemiMajorAxis)

"""
    SemiMinorAxis

Semi-minor axis of an ellipse.
"""
@ellipse_parameter(SemiMinorAxis)

"""
    AxisRatio

Semi-major axis divided by semi-minor axis of an ellipse.
"""
@ellipse_parameter(AxisRatio)

"""
    Flattening

Flattening of an ellipse: `(a - b) / a` where `a` and `b` are the semi-major and
semi-minor axes, respectively.
"""
@ellipse_parameter(Flattening)

"""
    InverseFlattening

Inverse flattening of an ellipse: `1 / f` where `f` is the flattening.
"""
@ellipse_parameter(InverseFlattening)

"""
    Eccentricity

Eccentricity of an ellipse: `sqrt(1 - (b / a)^2)` where `a` and `b` are the
semi-major and semi-minor axes, respectively.
"""
@ellipse_parameter(Eccentricity)

"""
    SquaredEccentricity

Squared eccentricity of an ellipse: `1 - (b / a)^2` where `a` and `b` are the
semi-major and semi-minor axes, respectively.
"""
@ellipse_parameter(SquaredEccentricity)

"""
    EllipseAFi

Ellipse defined by its semi-major axis and inverse flattening.
"""
@ellipse(EllipseAFi, a::SemiMajorAxis, fi::InverseFlattening)
(P::Type{<:SemiMajorAxis})((; a, fi)::EllipseAFi) = P(a)
(P::Type{<:SemiMinorAxis})((; a, fi)::EllipseAFi) = P(a * (1 - inv(fi)))
(P::Type{<:AxisRatio})((; a, fi)::EllipseAFi) = P(fi / (fi - 1))
(P::Type{<:Flattening})((; a, fi)::EllipseAFi) = P(inv(fi))
(P::Type{<:InverseFlattening})((; a, fi)::EllipseAFi) = P(fi)
(P::Type{<:SquaredEccentricity})((; a, fi)::EllipseAFi) = P((2 - inv(fi)) / fi)
(P::Type{<:Eccentricity})(ellipse::EllipseAFi) = P(sqrt(ellipse[SquaredEccentricity]))

"""
    EllipseAB

Ellipse defined by its semi-major and semi-minor axes.
"""
@ellipse(EllipseAB, a::SemiMajorAxis, b::SemiMinorAxis)
(P::Type{<:SemiMajorAxis})((; a, b)::EllipseAB) = P(a)
(P::Type{<:SemiMinorAxis})((; a, b)::EllipseAB) = P(b)
(P::Type{<:AxisRatio})((; a, b)::EllipseAB) = P(a / b)
(P::Type{<:Flattening})((; a, b)::EllipseAB) = P((a - b) / a)
(P::Type{<:InverseFlattening})((; a, b)::EllipseAB) = P(a / (a - b))
(P::Type{<:SquaredEccentricity})((; a, b)::EllipseAB) = P(1 - (b / a)^2)
(P::Type{<:Eccentricity})(ellipse::EllipseAB) = P(sqrt(ellipse[SquaredEccentricity]))

"""
    normal(ellipse, point) -> n

Compute a unit vector `n` such that the line `point + ℝ * n` intersects the
ellipse at a right angle. 
"""
normal(ellipse::Ellipse, point) = Tuple(normalize(SVector(unnormalized_normal(ellipse, point))))

"""
    unnormalized_normal(ellipse, point) -> n

Same as [`normal()`](@ref), but does not normalize the output vector for better
performance.
"""
function unnormalized_normal(ellipse::Ellipse, (p, z))
    # Implements the SC Halley algorithm from
    # https://www.researchgate.net/profile/Toshio-Fukushima/publication/264118552_F06-gconv2/data/53ce803e0cf279d93530b524/F06-gconv2.pdf
    a = ellipse[SemiMajorAxis]
    b = ellipse[SemiMinorAxis]
    E = ellipse[SquaredEccentricity]
    ec = b / a

    P = abs(p) / a
    Z = ec / a * abs(z)
    S = Z
    C = ec^2 * P # Paper says ec * P, but ec^2 * P is correct

    # TODO: Iterate to achieve machine precision?
    A = hypot(S, C)
    B = 1.5 * E^2 * P * (S * C)^2 * (A - ec)
    D = Z * A^3 + E * S^3
    F = P * A^3 - E * C^3
    S = D * F - B * S
    C = F^2 - B * C

    return ifelse(
        p != 0,
        (copysign(ec * C, p), copysign(S, z)),
        (zero(ec * C), copysign(one(S), z)),
    )
end
