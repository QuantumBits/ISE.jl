module ORB

using ..Types

using LinearAlgebra

using DynamicalSystems, Unitful

using SatelliteToolbox

using BenchmarkTools

struct input
    r::Position
    v::Velocity
    a::Acceleration
    t::Time
    t_prev::Time
end

struct commands

end

struct parameters
    μE::GravitationalParameter
    # EGM::GravityModel_Coefs
end

struct output
    r::Position
    v::Velocity
    a::Acceleration
    t::Time
end

function step(i::input, c::commands, p::parameters)::Tuple{output,commands}

    a = -p.μE .* i.r ./ norm(i.r)^3 + i.a
    # a = compute_g(p.EGM, ustrip.(upreferred.(i.r)))u"m/s^2"
    v = i.v .+ a .* (i.t - i.t_prev)
    r = i.r .+ i.v .* (i.t - i.t_prev)
    t = i.t

    return (
        output(r, v, a, t),
        commands()
    )

end


function test_run()

    RE = upreferred(6371.009u"km")
    μE = upreferred(398600.4418u"km^3/s^2")
    R0 = upreferred(RE + 400u"km")
    V0 = upreferred(sqrt(μE / R0))

    r0 = Position(R0, 0.0u"m", 0.0u"m")
    v0 = Velocity(0.0u"m/s", V0, 0.0u"m/s")
    a0 = Acceleration(0.0u"m/s^2", 0.0u"m/s^2", 0.0u"m/s^2")
    t0 = 0.0u"s"

    orb_i = input(r0, v0, a0, t0 + 1.0u"s", t0)
    orb_c = commands()
    orb_p = parameters(μE)

    orb_o, orb_c = step(orb_i, orb_c, orb_p)

    return (orb_i, orb_c, orb_p, orb_o)

end

end