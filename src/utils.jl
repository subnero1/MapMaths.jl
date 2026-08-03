"""
    numtype(x) -> T <: Number

Number type used by `x`. 

`numtype(x)` is similar and often identical to `eltype(x)` but differs in two
important ways:

- `numtype()` uses promotion while `eltype(x)` uses typejoins. 
  ```jldoctest
  julia> numtype(Tuple{Int, Float64})
  Float64

  julia> eltype(Tuple{Int, Float64})
  Real
  ```

- `numtype(x)` usually differs from `eltype(x)` for collections.
  ```julia
  julia> numtype([(1,2), (3,4)])
  Int64
  
  julia> eltype([(1,2), (3,4)])
  Tuple{Int64, Int64}
  ```

`numtype()` accepts both values and types. The type version should be defined
for new types. 
"""
function numtype end

numtype(x) = numtype(typeof(x))
numtype(T::Type) = throw(MethodError(numtype, (T,)))
numtype(T::Type{<:Number}) = T
numtype(T::Type{<:Tuple}) = promote_numtype(fieldtypes(T)...)
numtype(T::Type{<:AbstractArray}) = numtype(eltype(T))

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
    with_numtype(x, T::Type{<:Number})

Convert `x` to [`numtype`](@ref) `T`.
"""
function with_numtype end

"""
    is_set_equal(a::Tuple, b::Tuple)

Determine whether `a` and `b` have the same elements. Allocation-free
alternative to `issetequal(a, b)`.
"""
is_set_equal(a::Tuple, b::Tuple) = length(a) == length(b) && all(x -> x in b, a)


"""
    flatten(x::Tuple) -> Tuple

Strip away up to one level of nested tuples.

# Example
```jldoctest
julia> flatten(((1, 2), 3, (4, (5, 6))))
(1, 2, 3, 4, (5, 6))
```
"""
flatten(x::Tuple) = flatten_splat(x...)
flatten_splat() = ()
flatten_splat(x, y...) = (x, flatten_splat(y...)...)
flatten_splat(x::Tuple, y...) = (x..., flatten_splat(y...)...)


"""
    @symmetric function_definition

Symmetrize the given function definition. 

# Example
```jldoctest
julia> @symmetric foo(x::Int, y::Float64) = x + y;

julia> foo(1, 1.0) == foo(1.0, 1) == 2.0
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
    arg_defs = map(
        ((i, arg_def),) -> begin
            @assert isexpr(arg_def, :(::))
            return :($(Symbol("x$i"))::$(arg_def.args[end]))
        end, enumerate(signature.args[2:end])
    )
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
                    permutations(arg_defs),
                )...
            )
        end,
    )
    return nothing
end


"""
    @print_var(io, var)

Print `var` along with its module path if it is not imported into `Main`, and as
just `"\$var"` if it is imported. 

# Example
```julia
julia> module M
           using MapMaths: @print_var
           function foo() 
               @print_var(stdout, foo)
               println("()")
           end
       end;

julia> M.foo()
Main.M.foo()

julia> using .M: foo

julia> foo()
foo()
```
"""
macro print_var(io, var)
    io = esc(io)
    var = esc(var)
    symbol = Meta.quot(unescape(var))
    return quote
        if isdefined(Main, $symbol) && getfield(Main, $symbol) == $var
            print($io, $symbol)
        else
            print($io, @__MODULE__, ".", $symbol)
        end
    end
end
