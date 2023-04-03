module GPS

using ..Types

using Random

using DynamicalSystems, Unitful

using BenchmarkTools

struct input
    r::Position
    v::Velocity
    t::Time
end

struct commands

end

struct parameters

    rng_r::MersenneTwister
    std_r::Length

    rng_v::MersenneTwister
    std_v::Speed

    rng_t::MersenneTwister
    std_t::Time

end

struct output
    r_est::Position
    v_est::Velocity
    t_est::Time
end

function step(i::input, c::commands, p::parameters)::Tuple{output,commands}

    r_est = i.r .+ p.std_r .* randn(p.rng_r, 3)
    v_est = i.v .+ p.std_v .* randn(p.rng_v, 3)
    t_est = i.t + p.std_t * randn(p.rng_t)

    return (
        output(r_est, v_est, t_est),
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
    t0 = 0.0u"s"

    gps_i = input(r0, v0, t0)
    gps_c = commands()
    gps_p = parameters(
        MersenneTwister(1), upreferred(10.0u"m"),
        MersenneTwister(2), upreferred(1.0u"mm/s"),
        MersenneTwister(3), upreferred(1.0u"μs"))

    gps_o, gps_c = step(gps_i, gps_c, gps_p)

    return (gps_i, gps_c, gps_p, gps_o)

end

end