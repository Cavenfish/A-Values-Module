using Plots
using LinearAlgebra

#Hertz–Knudsen equation
function sticking(α, ρ, M, T)
    R = 8.31446261815 #Gas Constant
    Avo = 6.02214076e23 #Avogadors Number
    (α * ρ * Avo) / √(2 * π * M * R * T)
end

#Snell's Law solved for θt
function snell_θt(n1, n2, θi)
    asin((n1 / n2) * sin(θi))
end

#(Heavens et al): interferance pattern reflection formula
function R(r1, r2, β)
    (r1 .+ r2 .* ℯ.^(-2im .* β)) ./ (1 .+ r1 .* r2 .* ℯ.^(-2im .* β))
end

#(Heavens et al): interferance pattern transmission formula
function T(t1, t2, r1, r2, β)
    (t1 .* t2 .*ℯ.^(-1im .* β)) ./ (1 .+ r1 .* r2 .* ℯ.^(-2im .* β))
end

#Phase Difference Relation
function δ(n, d, ϕ, λ)
    (2 .* π .* n .* d .* cos(ϕ)) ./ (λ)
end

#Generate the y-values for the pattern graph
function pattern(t_rate, k, λ, θi, n, l, I, front)
    df = t_rate .* k
    db = t_rate

    if front
        df = t_rate
        db = t_rate .* k
    end

    θ1 = snell_θt(n[1],n[2],θi)
    θ2 = snell_θt(n[2],n[3],θ1)
    θ3 = snell_θt(n[3],n[2],θ2)

    r1 = (n[1] - n[2]) / (n[1] + n[2])
    r2 = (n[2] - n[3]) / (n[2] + n[3])
    r3 = (n[3] - n[2]) / (n[3] + n[2])
    r4 = (n[2] - n[1]) / (n[2] + n[1])

    t1 = (2*n[1]) / (n[1] + n[2])
    t2 = (2*n[2]) / (n[2] + n[3])
    t3 = (2*n[3]) / (n[3] + n[2])
    t4 = (2*n[2]) / (n[2] + n[1])

    β_front = δ(n[2], df, θ1, λ)
    γ       = δ(n[3], l, θ2, λ)
    β_back  = δ(n[2], db, θ3, λ)
    T1      = (abs.(T(t1, t2, r1, r2, β_front))).^2
    T2      = (abs.(T(t3, t4, r3, r4, β_back))).^2
    R1      = (abs.(R(r1, r2, β_front))).^2
    R2      = (abs.(R(r3, r4, β_back))).^2

    [(R(R1,(T2 .* T1 .* R2),γ) ), (db .+ df), (R1), (df)]
end

Temp = 10
M    = .01801528
α    = 1
ρ    = 6.6661e-5
k   = (.1)
d    = 0.5e-3
λ1   = 405e-9
λ2   = 532e-9
θi   = 0.785398
n    = [1, 1.308, 1.527]
h2oV = 2.99e-29
time = 1:5500
I    = 20000

#Particles/Molecules per second times particle/molecule size = thickness rate
Φ = sticking(α, ρ, M, Temp) * h2oV


ice  = Φ .* time

p1 = pattern(ice, k, λ1, θi, n, d, I, true)
p2 = pattern(ice, k, λ1, θi, n, d, I, false)

y1 = (abs.(p1[1])).^2 .*I
y2 = (abs.(p2[1])).^2 .*I
y3 = (p1[3]) .*I

plot(layout=(3,1))
plot!(p1[2], y1, label="Front Laser", subplot=1)
plot!(p2[2], y2, label="Back Laser", color=:green,  subplot=2)
plot!(p1[4], y3, label="R1 Only", color=:red, subplot=3)
title!("Theoretical Laser Interference Pattern", subplot=1)
ylabel!("Laser Intensity (arbitrary units)", subplot=2)
xlabel!("Ice Thinkness (m)", subplot=3)
