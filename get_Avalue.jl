using QuadGK

function α(v)
    4*π*v*k(v)
end

function A(ρ, v_min, v_max)
    (1/ρ) * quadgk(α, v_min, v_max, rtol=1e-3)
end
