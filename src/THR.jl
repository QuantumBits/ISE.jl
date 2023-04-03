module THR

using ..Types

using LinearAlgebra

using DynamicalSystems, Unitful

using BenchmarkTools


struct input
    m_prop::Mass    # Propellant mass
    t::Time         # Current time
    t_prev::Time    # Previous time
end

struct commands
    on::Bool
end

struct parameters
    m_dry::Mass                 # Dry mass
    Isp::Time                   # Specific impulse
    g0::AccelerationMagnitude   # Gravitational acceleration at sea level
    dmdt::typeof(1.0u"kg/s")    # Mass flow rate
end

struct output
    T::AccelerationMagnitude    # Imparted thrust
    m_prop::Mass                # Remaining propellant mass
    t::Time                     # Current time
end

function step(i::input, c::commands, p::parameters)::Tuple{output,commands}

    Δt = i.t - i.t_prev
    m_tot = p.m_dry + i.m_prop


    # Commanded to turn on AND
    # Total mass is greater than zero AND
    # Time step is greater than zero
    if c.on && (m_tot > 0.0u"kg") && Δt > 0.0u"s"

        # Propellant consumed
        dm = max(0.0u"kg", min(i.m_prop, p.dmdt * Δt))

        # Mean Acceleration 
        T = (p.g0 * p.Isp * p.dmdt) / (m_tot - 0.5 * dm)

        # Remaining propellant mass
        m_prop = i.m_prop - dm

    else
        T = 0.0u"m/s^2"
        m_prop = i.m_prop
    end

    return (
        output(T, m_prop, i.t),
        commands(c.on)
    )

end


function test_run()

    Isp = upreferred(300.0u"s")
    g0 = upreferred(9.80665u"m/s^2")
    dmdt = upreferred(1.0u"kg/s")
    m_dry = upreferred(200.0u"kg")
    m_prop = upreferred(100.0u"kg")

    t0 = 0.0u"s"

    thr_i = input(m_prop, t0 + 1.0u"s", t0)
    thr_c = commands(true)
    thr_p = parameters(m_dry, Isp, g0, dmdt)

    thr_o, thr_c = step(thr_i, thr_c, thr_p)

    return (thr_i, thr_c, thr_p, thr_o)

end


end