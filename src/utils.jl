flatten() = ()
flatten(x, y...) = (x, flatten(y...)...)
flatten(x::Tuple, y...) = (flatten(x...)..., flatten(y...)...)

permutations(a) = ((a,))
permutations(a, b) = ((a, b), (b, a))
permutations(a, b, c) = ((a, b, c), (a, c, b), (b, a, c), (b, c, a), (c, a, b), (c, b, a))

numtype(x) = numtype(typeof(x))
numtype(T::Type) = throw(MethodError(numtype, (T,)))
numtype(T::Type{<:Number}) = T
numtype(T::Type{<:NTuple{<:Any, Number}}) = promote_numtype(fieldtypes(T)...)

promote_numtype(args...) = promote_type(numtype.(args)...)
default_numtype(args...) = default_numtype(promote_numtype(args...))
default_numtype(x::Number) = convert(default_numtype(typeof(x)), x)
default_numtype(T::Type{<:Number}) = T
default_numtype(T::Type{<:Integer}) = float(T)

maybe_only(x::Vararg{<:Any, 1}) = only(x)
maybe_only(x::Vararg) = x

abstract type Combineable end
struct Combo{P <: Tuple} <: Combineable
    parts::P
    Combo{P}(parts) where {P <: Tuple} = new{P}(parts)
end
struct ComboTV{T, V <: Combineable}
    value::V
    ComboTV{T, V}(value) where {T, V <: Combineable} = new{T, V}(value)
end
Combo(P::Type{<:Combineable}...) = Combo{<:Tuple{P...}}
Combo(T::Type{<:Type}, V::Type{<:Combineable}) = ComboTV{<:T.body.parameters[1].ub, <:V}
Base.:|(parts::Vararg{Combineable}) = Combo{typeof(parts)}(parts)
Base.:|(combo::Combo, part::Combineable) = |(combo.parts..., part)
Base.show(io::IO, combo::Combo) = join(io, combo.parts, " | ")
Base.show(io::IO, ::Type{Combo{P}}) where {P <: Tuple} = print(io, Combo, "(", join(P.parameters, ", "), ")")
Base.iterate(c::Combo, state...) = iterate(c.parts, state...)

Base.:|(T::Type, value::Combineable) = ComboTV{T, typeof(value)}(value)
Base.iterate(c::ComboTV{T}, state...) where {T} = iterate((T, c.value), state...)
Base.show(io::IO, (T, v)::ComboTV) = print(io, T, " | ", v)

macro symmetric(fun_def)
    @assert fun_def.head in (:function, :(=))
    @assert fun_def.args[1].head == :call
    name = fun_def.args[1].args[1]
    args = fun_def.args[1].args[2:end]
    return esc(
        quote
            $fun_def
            $(
                map(
                    permuted_args -> begin
                        if all(permuted_args .== args)
                            return nothing
                        end
                        return :($name($(permuted_args...)) = $name($(args...)))
                    end,
                    permutations(args...),
                )...
            )
        end,
    )
    return nothing
end