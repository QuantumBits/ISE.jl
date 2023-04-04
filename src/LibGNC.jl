module LibGNC

const size_t = Csize_t
const libGNC = joinpath(@__DIR__, "..", "lib", "libGNC.so")


struct input
    r::NTuple{3, Cdouble}
    v::NTuple{3, Cdouble}
    m_prop::Cdouble
end

struct commands
    NOOP::Bool
end

struct parameters
    mu::Cdouble
end

struct output
    thr_on::Bool
    thr_dir::NTuple{3, Cdouble}
    specific_energy::Cdouble
end

function step(i, c, p)
    ccall((:step, libGNC), output, (Ptr{input}, Ptr{commands}, Ptr{parameters}), i, c, p)
end

end # module
