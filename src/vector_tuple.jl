function vector_tuple(Type::Symbol, SuperType::Symbol, N)
    return quote
        struct $Type{T <: Number} <: $SuperType
            numbers::NTuple{$N, T}
            $Type{T}(numbers::Vararg{Number, $N}) where {T <: Number} = new{T}(numbers)
        end

        $Type(numbers::Vararg{Number, $N}) = $Type{default_numtype(numbers)}(numbers)
        (::Type{V})(numbers::NTuple{$N, Number}) where {V <: $Type} = V(numbers...)
        (::Type{V})(numbers::StaticVector{$N, <:Number}) where {V <: $Type} = V(numbers...)
        (::Type{V})(v::$Type) where {V <: $Type} = V(v.numbers)

        numtype(::Type{$Type{T}}) where {T <: Number} = T
        # eltype(::Type{<:$Type{T}}) where {T <: Number} = T
        # length(::Type{<:$Type}) = $N

        function Base.show(io::IO, v::$Type)
            print(io, typeof(v), "(")
            join(io, v.numbers, ", ")
            return print(io, ")")
        end
    end
end

macro vector_tuple(Type_sub_SuperType::Expr, N)
    if Type_sub_SuperType.head == :(<:)
        (Type::Symbol, SuperType::Symbol) = Type_sub_SuperType.args
    elseif Type_sub_SuperType isa Symbol
        Type = Type_sub_SuperType
        SuperType = :Any
    else
        error(
            "First argument to @vector_tuple must be either a symbol or a subtype (<:) expression of two symbols",
        )
    end
    return esc(vector_tuple(Type, SuperType, N))
end