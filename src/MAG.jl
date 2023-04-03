module MAG

using SatelliteToolbox, Unitful

function step(u, p, t)

    t0 = p[1]
    r, λ, Ω = u

    return SVector(r, λ, Ω, uconvert.(u"T", igrf(t0 + t,
        ustrip(uconvert(u"m", r)),
        ustrip(uconvert(u"rad", λ)),
        ustrip(uconvert(u"rad", Ω)), Val(:geocentric))u"nT")...)

end

end