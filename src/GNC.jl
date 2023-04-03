module GNC

using ..Types

using LinearAlgebra

using DynamicalSystems, Unitful

using BenchmarkTools


struct input
    r::Position     # Current position
    v::Velocity     # Current velocity
    m_prop::Mass    # Current propellant mass
end

struct commands
end

struct parameters
    μE::GravitationalParameter
end

struct output
    thr_on::Bool
    thr_dir::SVector{3,Float64}
    specific_energy::SpecificEnergy
end

function step(i::input, c::commands, p::parameters)::Tuple{output,commands}

    ϵ0 = ϵ(i.r, i.v, p.μE)

    thr_on = (i.m_prop > 0.0u"kg") && (ϵ0 < 0.0u"m^2/s^2")

    if thr_on
        thr_direction = normalize(ustrip.(i.v))
    else
        thr_direction = SVector(0.0, 0.0, 0.0)
    end

    return (
        output(thr_on, thr_direction, ϵ0),
        commands()
    )

end

"""
    ϵ(r::Position, v::Velocity, μ::GravitationalParameter)

Orbital specific energy (m^2/s^2)
"""
function ϵ(r::Position, v::Velocity, μ::GravitationalParameter)

    R = norm(r)
    V = norm(v)

    return upreferred(V^2 / 2 - μ / R)
end


function test_run()

    RE = upreferred(6371.009u"km")
    μE = upreferred(398600.4418u"km^3/s^2")
    R0 = upreferred(RE + 400u"km")
    V0 = upreferred(sqrt(μE / R0))

    r0 = Position(R0, 0.0u"m", 0.0u"m")
    v0 = Velocity(0.0u"m/s", V0, 0.0u"m/s")
    mp = 100.0u"kg"

    gnc_i = input(r0, v0, mp)
    gnc_c = commands()
    gnc_p = parameters(μE)

    gnc_o, gnc_c = step(gnc_i, gnc_c, gnc_p)

    return (gnc_i, gnc_c, gnc_p, gnc_o)

end

end