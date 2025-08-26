numtype(x) = numtype(typeof(x))
numtype(T::Type) = throw(MethodError(numtype, (T,)))
numtype(T::Type{<:Number}) = T
numtype(T::Type{<:NTuple{<:Any, Number}}) = promote_numtype(fieldtypes(T)...)

promote_numtype(args...) = promote_type(numtype.(args)...)
default_numtype(args...) = float(promote_numtype(args...))