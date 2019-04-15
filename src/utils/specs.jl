"""
This file is to test ideas on creating a structure with specs
"""

# import Base: getfield, getproperty

struct Spec
    # name::Char
    # shape::Tuple{Int, Int}
    # startID::Int
    # N::Int
    name
    shape
    startID
    N

    function Spec(name::String, shape::Tuple{Int, Int}, start::Int)
        new(name, shape, start, prod(shape))
    end
end

# abstract type FinancialModel end
#
# struct Model <: FinancialModel
#     x
#
#     function Model(x = [1 2 4 5 6 7])
#         new(x)
#     end
# end
#
# model = Model()
#
# function getproperty(obj::FinancialModel, spec::Spec)
#     return getfield(obj, spec)
# end
#
# function getfield(obj::FinancialModel, spec::Spec)
#     return view(obj, )
#     reshape(view(J, N+1:2*N, 2*N+1:3*N), N, N)
# end
#
# """
#     nelem(x)
# Number of elements that belong to Spec 'x'.
# """
# function nelem(x::Spec) prod(x.shape) end
#
# #---------------------------------------------------#
#
# struct MyType
#     x
# end
#
# function Base.getproperty(obj::MyType, sym::Symbol)
#     if sym === :special
#         return obj.x + 1
#     else # fallback to getfield
#         return getfield(obj, sym)
#     end
# end
