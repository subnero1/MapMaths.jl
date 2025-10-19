flatten() = ()
flatten(x, y...) = (x, flatten(y...)...)
flatten(x::Tuple, y...) = (flatten(x...)..., flatten(y...)...)

permutations(a) = ((a,))
permutations(a, b) = ((a, b), (b, a))
permutations(a, b, c) = ((a, b, c), (a, c, b), (b, a, c), (b, c, a), (c, a, b), (c, b, a))

"""
    numtype(x)

Concrete `Number` type associated with `x`. 

`numtype(x)` is similar and often identical to `eltype(x)` but differs in two
important ways:

- `numtype(x)` is a concrete, promoted type while `eltype(x)` is a potentially
  abstract, typejoined type. 
  ```jldoctest
  julia> numtype(Tuple{Int, Float64})
  Float64

  julia> eltype(Tuple{Int, Float64})
  Real
  ```

- `numtype(x)` may differ from `eltype(x)` for composite types. For example,
  given a type 
  ```julia
  struct Particles{T<:Number} 
      locations::Vector{SVector{3,T}} 
  end
  ```
  the appropriate definitions for `numtype` and `eltype` would be
  ```julia
  numtype(::Type{Particles{T}}) where {T<:Number} = T
  eltype(::Type{Particles{T}}) where {T<:Number} = SVector{3,T}
  ```

`numtype()` accepts both values and types. The type version should be defined
for new types. 
"""
function numtype end

numtype(x) = numtype(typeof(x))
numtype(T::Type) = throw(MethodError(numtype, (T,)))
numtype(T::Type{<:Number}) = T
numtype(T::Type{<:NTuple{<:Any, Number}}) = promote_numtype(fieldtypes(T)...)

"""
    promote_numtype(args...) = promote_type(numtype.(args)...)

Convenience function to promote the [`numtype`](@ref) of multiple arguments.
"""
promote_numtype(args...) = promote_type(numtype.(args)...)

"""
    default_numtype(args...)

The default [`numtype`](@ref) to be used in computations involving `args`.

Similar to `promote_numtype` but promotes `Integer` types to floating point types.
"""
default_numtype(args...) = default_numtype(promote_numtype(args...))
default_numtype(T::Type{<:Number}) = T
default_numtype(T::Type{<:Integer}) = float(T)

"""
    @number_tuple(Type [<: Supertype], N)

Define a new type `Type{T <: Number}` that (mostly) quacks like a `NTuple{N,
T}`. The number type `T` is derived using [`default_numtype`](@ref), i.e.
`Integer` types are promoted to floating point types by default.
"""
macro number_tuple(Type_maybe_sub_Supertype, N::Integer)
    if Type_maybe_sub_Supertype isa Symbol
        Type = esc(Type_maybe_sub_Supertype)
        type_signature = :($Type{T <: Number})
    elseif isexpr(Type_maybe_sub_Supertype, :(<:), 2)
        (Type, Supertype) = esc.(Type_maybe_sub_Supertype.args)
        type_signature = :($Type{T <: Number} <: $Supertype)
    else
        error("First argument to @number_tuple must be either `Type` or `Type <: Supertype`")
    end
    return quote
        struct $type_signature
            elements::NTuple{$N, T}
            $Type{T}(elements::Vararg{Number, $N}) where {T <: Number} = new{T}(elements)
        end

        $Type(numbers::Vararg{Number, $N}) = $Type{default_numtype(numbers)}(numbers)
        (C::Type{<:$Type})(numbers::NTuple{$N, Number}) = C(numbers...)
        (C::Type{<:$Type})(numbers::StaticVector{$N, <:Number}) = C(numbers...)
        (C::Type{<:$Type})(c::$Type) = C(c.elements)

        Base.length(::Union{$Type, Type{<:$Type}}) = $N
        Base.eltype(c::$Type{T}) where {T <: Number} = T
        Base.iterate(c::$Type, state...) = iterate(c.elements, state...)
        Base.Tuple(c::$Type) = c.elements
        MapMaths.numtype(::Type{$Type{T}}) where {T <: Number} = T
        function Base.show(io::IO, c::$Type)
            print(io, typeof(c), "(")
            join(io, Tuple(c), ", ")
            print(io, ")")
            return nothing
        end
    end
end

"""
    @symmetric function_definition

Symmetrize the given function definition. 

# Example
```jldoctest
julia> @symmetric foo(x::Int, y::Float64) = x + y;

julia> # Argument order doesn't matter:
       foo(1, 1.0) == foo(1.0, 1) == 2.0
true
```
"""
macro symmetric(fun_def)
    @assert isexpr(fun_def, :function) || isexpr(fun_def, :(=))
    signature_and_maybe_where_clause = fun_def.args[1]
    signature = if isexpr(signature_and_maybe_where_clause, :where)
        signature_and_maybe_where_clause.args[1]
    else
        signature_and_maybe_where_clause
    end
    @assert isexpr(signature, :call)
    name = signature.args[1]
    arg_defs = map(((i, arg_def),) -> begin
        @assert isexpr(arg_def, :(::))
        return :($(Symbol("x$i"))::$(arg_def.args[end]))
    end, enumerate(signature.args[2:end]))
    args = ntuple(i -> Symbol("x$i"), length(arg_defs))
    call = if isexpr(signature_and_maybe_where_clause, :where)
        where_clause = signature_and_maybe_where_clause.args[2:end]
        arg_defs -> :($name($(arg_defs...)) where {$(where_clause...)})
    else
        arg_defs -> :($name($(arg_defs...)))
    end
    return esc(
        quote
            $fun_def
            $(
                map(
                    permuted_arg_defs -> begin
                        if all(permuted_arg_defs .== arg_defs)
                            return nothing
                        end
                        return :($(call(permuted_arg_defs)) = $name($(args...)))
                    end,
                    permutations(arg_defs...),
                )...
            )
        end,
    )
    return nothing
end

strip_macros(expr) =
    if isexpr(expr, :macrocall, 3)
        return strip_macros(expr.args[3])
    else
        return expr
    end

"""
    @ampersand @combo struct Type [<: Supertype]
        (PartType1, PartType2, ...)
    end

Add an `&()` method which constructs a `Type` value from its parts.
"""
macro ampersand(full_struct_def::Expr)
    struct_def = strip_macros(full_struct_def)
    @assert isexpr(struct_def, :struct)
    @assert struct_def.args[1] == false # mutability
    Type_maybe_sub_Supertype = struct_def.args[2]
    body = struct_def.args[3]
    @assert isexpr(body, :block)
    @assert isexpr(body.args[2], :tuple)
    Parts = esc.(body.args[2].args)
    if Type_maybe_sub_Supertype isa Symbol
        Type = esc(Type_maybe_sub_Supertype)
    elseif isexpr(Type_maybe_sub_Supertype, :(<:), 2)
        Type = esc(Type_maybe_sub_Supertype.args[1])
    else
        error("Unexpected type signature: $Type_maybe_sub_Supertype")
    end
    return quote
        $(esc(full_struct_def))

        @symmetric Base.:&($(map(((i, Part),) -> :($(Symbol("c$i"))::$Part), enumerate(Parts))...)) =
            $Type($(ntuple(i -> Symbol("c$i"), length(Parts))...))
    end
end

"""
    @combo struct Type [<: Supertype]
        (PartType1, PartType2, ...)
    end

Define a new type for unordered combinations of `PartType` values. Such
combinations can be constructed using `Type(part1, part2...)`, and the parts can
be recovered using `PartTypeX(::Type)`.

# Example
```jldoctest
julia> struct Real x::Float64 end

julia> struct Imag x::Float64 end

julia> @combo struct Complex; (Real, Imag) end

julia> Complex(Real(1), Imag(2))
Complex(Real(1.0), Imag(2.0))

julia> Complex(Imag(2), Real(1)) # Order doesn't matter
Complex(Real(1.0), Imag(2.0))

julia> Complex(Real(1), Imag(2)) |> Imag
Imag(2.0)
```
"""
macro combo(struct_def)
    return combo(struct_def)
end

function parts end

function combo(struct_def)
    @assert isexpr(struct_def, :struct)
    @assert struct_def.args[1] == false # mutability
    Type_maybe_sub_Supertype = struct_def.args[2]
    body = struct_def.args[3]
    @assert isexpr(body, :block)
    @assert isexpr(body.args[2], :tuple)
    Parts = esc.(body.args[2].args)
    PartTypeVars = ntuple(i -> Symbol("P$i"), length(Parts))
    PartTypeVarDefs = map((V, T) -> :($V <: $T), PartTypeVars, Parts)
    part_vars = ntuple(i -> Symbol("part_$i"), length(Parts))
    part_var_defs = map((part, Part) -> :($part::$Part), part_vars, Parts)
    if Type_maybe_sub_Supertype isa Symbol
        Type = esc(Type_maybe_sub_Supertype)
        type_signature = :($Type{$(PartTypeVarDefs...)})
    elseif isexpr(Type_maybe_sub_Supertype, :(<:), 2)
        (Type, Supertype) = esc.(Type_maybe_sub_Supertype.args)
        type_signature = :($Type{$(PartTypeVarDefs...)} <: $Supertype)
    else
        error("Unexpected type signature: $Type_maybe_sub_Supertype")
    end
    return quote
        struct $type_signature
            parts::Tuple{$(PartTypeVars...)}
            @symmetric function $Type{$(PartTypeVars...)}(
                $(map((part_var, PartTypeVar) -> :($part_var::$PartTypeVar), part_vars, PartTypeVars)...),
            ) where {$(PartTypeVarDefs...)}
                return new{$(PartTypeVars...)}(tuple($(part_vars...)))
            end
        end

        @symmetric function $Type{$(PartTypeVars...)}($(part_var_defs...)) where {$(PartTypeVarDefs...)}
            return $Type{$(PartTypeVars...)}(
                $(map((part_var, PartTypeVar) -> :($PartTypeVar($part_var)), part_vars, PartTypeVars)...),
            )
        end

        @symmetric $Type($(part_var_defs...)) =
            $Type{$(map(part_var -> :(typeof($part_var)), part_vars)...)}($(part_vars...))
        # @symmetric $Type{$(PartTypeVars...)}($(part_var_defs...)) where {$(PartTypeVarDefs...)} =
        #     $Type{$(PartTypeVars...)}((
        #         $(map((part_var, PartTypeVar) -> :($PartTypeVar($part_var)), part_vars, PartTypeVars)...),
        #     ))
        # @symmetric (C::Type{<:$Type})($(part_var_defs...)) = C(tuple($(part_vars...)))
        (C::Type{<:$Type})(c::$Type) = C(c.parts...)

        $(map(((i, Part),) -> quote
            (P::Type{<:$Part})(c::$Type) = P(c.parts[$i])
        end, enumerate(Parts))...)

        MapMaths.parts(c::$Type) = c.parts

        function Base.show(io::IO, c::$Type)
            print(io, $Type, "(")
            join(io, c.parts, ", ")
            print(io, ")")
            return nothing
        end
    end
end

"""
    @conversion_constructible abstract type Type [<: Supertype] end

Define an abstract type `Type` along with an identity conversion constructor.

The conversion constructor added by this macro collides with the default outer
constructors of the subtypes. Therefore, subtypes should make sure to remove the
default outer constructor by defining an explicit inner constructor. 

See also [`@convertible`](@ref).
"""
macro conversion_constructible(abstract_type_def::Expr)
    @assert isexpr(abstract_type_def, :abstract)
    Type_maybe_sub_Supertype = abstract_type_def.args[1]
    if Type_maybe_sub_Supertype isa Symbol
        Type = esc(Type_maybe_sub_Supertype)
        type_signature = Type
    elseif isexpr(Type_maybe_sub_Supertype, :(<:), 2)
        (Type, Supertype) = esc.(Type_maybe_sub_Supertype.args)
        type_signature = :($Type <: $Supertype)
    else
        error("Unexpected type signature: $Type_maybe_sub_Supertype")
    end
    return quote
        abstract type $type_signature end
        (::Type{C})(c::C) where {C <: $Type} = c
        Base.convert(C::Type{<:$Type}, c::$Type) = C(c)
    end
end

"""
    @convertible abstract type Type [<: Supertype] end

Define an abstract type `Type` along with a conversion method and an identity
conversion constructor for itself and all its subtypes.

The conversion constructor added by this macro collides with the default outer
constructors of the subtypes. Therefore, subtypes should make sure to remove the
default outer constructor by defining an explicit inner constructor. 

See also [`@conversion_constructible`](@ref).

# Example
```jldoctest
julia> @convertible abstract type MyFloat end

julia> struct MyFloat32 <: MyFloat 
           x::Float32 
           MyFloat32(x::Number) = new(x) # Remove default outer constructor
       end

julia> # Now you can conversion-construct `MyFloat32` to itself and `MyFloat`:
       @assert MyFloat32(MyFloat32(1f0)) == MyFloat32(1f0)
       @assert MyFloat(MyFloat32(1f0)) == MyFloat32(1f0)

julia> # However, you can't convert `MyFloat32` to `MyFloat64` yet:
       struct MyFloat64 <: MyFloat 
           x::Float64 
           MyFloat64(x::Number) = new(x) # Remove default outer constructor
       end

julia> MyFloat64(MyFloat32(1f0)) 
ERROR: MethodError: no method matching MyFloat64(::MyFloat32)
The type `MyFloat64` exists, but no method is defined for this combination of
argument types when trying to construct it.
[...]

julia> # Let's fix that
       MyFloat64(float32::MyFloat32) = MyFloat64(float32.x);

julia> # Now we can conversion-construct and also convert
       @assert MyFloat64(MyFloat32(1f0)) == MyFloat64(1.0)
       @assert convert(MyFloat64, MyFloat32(1f0)) == MyFloat64(1.0)
```
"""
macro convertible(abstract_type_def::Expr)
    @assert isexpr(abstract_type_def, :abstract)
    Type_maybe_sub_Supertype = abstract_type_def.args[1]
    if Type_maybe_sub_Supertype isa Symbol
        Type = esc(Type_maybe_sub_Supertype)
        type_signature = Type
    elseif isexpr(Type_maybe_sub_Supertype, :(<:), 2)
        (Type, Supertype) = esc.(Type_maybe_sub_Supertype.args)
        type_signature = :($Type <: $Supertype)
    else
        error("Unexpected type signature: $Type_maybe_sub_Supertype")
    end
    return quote
        abstract type $type_signature end
        (::Type{C})(c::C) where {C <: $Type} = c
        Base.convert(C::Type{<:$Type}, c::$Type) = C(c)
    end
end

"""
    @static_length(Type, N)

Add a method to `Base.length` for `Type`.
"""
macro static_length(Type, N)
    Type = esc(Type)
    N = esc(N)
    return :(Base.length(::Union{$Type, Type{<:$Type}}) = $N)
end