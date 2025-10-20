"""
    @coordinate_combo struct Type [<: Supertype]
        (PartType1, PartType2, ...)
    end

Extension of [`@combo`](@ref) which defines additional methods specific to
coordinate types.
"""
macro coordinate_combo(struct_def::Expr)
    @assert isexpr(struct_def, :struct)
    @assert struct_def.args[1] == false # mutability
    Type_maybe_sub_Supertype = struct_def.args[2]
    body = struct_def.args[3]
    @assert isexpr(body, :block)
    @assert isexpr(body.args[2], :tuple)
    Parts = esc.(body.args[2].args)
    PartTypeVars = ntuple(i -> Symbol("P$i"), length(Parts))
    PartTypeVarDefs = map((V, T) -> :($V <: $T), PartTypeVars, Parts)
    if Type_maybe_sub_Supertype isa Symbol
        Type = esc(Type_maybe_sub_Supertype)
        Supertype = nothing
        type_signature = :($Type{$(PartTypeVarDefs...)})
    elseif isexpr(Type_maybe_sub_Supertype, :(<:), 2)
        (Type, Supertype) = esc.(Type_maybe_sub_Supertype.args)
        type_signature = :($Type{$(PartTypeVarDefs...)} <: $Supertype)
    else
        error("Unexpected type signature: $Type_maybe_sub_Supertype")
    end
    return quote
        $(
            if !isnothing(Supertype)
                quote
                    if hasmethod(length, Tuple{Type{$Supertype}})
                        @assert sum(length.(tuple($(Parts...)))) == length($Supertype)
                    end
                end
            end
        )
        $(combo(struct_def))
        @static_length($Type, sum(length.(tuple($(Parts...)))))
        MapMaths.coords(c::$Type) = flatten(coords.(parts(c)))
        function MapMaths.numtype(::Type{$Type{$(PartTypeVars...)}}) where {$(PartTypeVarDefs...)}
            return promote_numtype($(PartTypeVars...))
        end
    end
end

"""
    @coordinate_number_tuple Type [<: Supertype], N

Extension of [`@number_tuple`](@ref) with additional methods specific to
coordinate types.
"""
macro coordinate_number_tuple(Type_maybe_sub_Supertype, N::Integer)
    if Type_maybe_sub_Supertype isa Symbol
        Type = esc(Type_maybe_sub_Supertype)
    elseif isexpr(Type_maybe_sub_Supertype, :(<:), 2)
        Type = esc(Type_maybe_sub_Supertype.args[1])
    else
        error("First argument to @coordinate_number_tuple must be either `Type` or `Type <: Supertype`")
    end
    return quote
        $(number_tuple(Type_maybe_sub_Supertype, N))
        MapMaths.coords(c::$Type) = Tuple(c)
    end
end

"""
    @named_coordinate_combo Name = ComboType{PartType1, PartType2, ...}

Add a type alias for a combination of coordinate types, along with several
method definitions to make the alias behave like a [`@number_tuple`](@ref).
"""
macro named_coordinate_combo(expr::Expr)
    @assert isexpr(expr, :(=), 2)
    @assert expr.args[1] isa Symbol
    Name = esc(expr.args[1])
    ComboType_and_PartTypes = expr.args[2]
    @assert isexpr(ComboType_and_PartTypes, :curly)
    ComboType = esc(ComboType_and_PartTypes.args[1])
    PartTypes = esc.(ComboType_and_PartTypes.args[2:end])
    part_vars = ntuple(i -> Symbol("part_$i"), length(PartTypes))
    part_var_defs = map((part, Part) -> :($part::supertype($Part)), part_vars, PartTypes)
    return quote
        const $Name{T <: Number} = $ComboType{$(map(PartType -> :($PartType{T}), PartTypes)...)}

        function $Name{T}(parts::Vararg{Number, length($ComboType)}) where {T <: Number}
            tail_1 = parts
            $(
                map(
                    ((i, PartType),) -> quote
                        $(Symbol("part_$i")) = $PartType{T}($(Symbol("tail_$i"))[1:length($PartType)])
                        $(Symbol("tail_$(i+1)")) = $(Symbol("tail_$i"))[length($PartType)+1:end]
                    end,
                    enumerate(PartTypes),
                )...
            )
            return $ComboType($(ntuple(i -> Symbol("part_$i"), length(PartTypes))...))
        end
        $Name(numbers::Vararg{Number, length($ComboType)}) = $Name{default_numtype(numbers)}(numbers)
        (C::Type{<:$Name})(numbers::NTuple{length($ComboType), Number}) = C(numbers...)
        (C::Type{<:$Name})(numbers::StaticVector{length($ComboType), <:Number}) = C(numbers...)
        (C::Type{<:$Name})(c::$Name) = C(parts(c)...)
        $Name($(part_var_defs...)) = $Name{promote_numtype($(part_vars...))}($(part_vars...))

        function Base.show(io::IO, c::$Name)
            print(io, $Name, "{", numtype(c), "}(")
            join(io, coords(c), ", ")
            return print(io, ")")
        end
    end
end
