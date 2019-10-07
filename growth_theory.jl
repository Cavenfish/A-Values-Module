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

#Fresnel equation for S-polarized light reflectance coeffiecient
function sR(n1, n2, θi, θt)
    abs((n1 * cos(θi) - n2 * cos(θt)) / (n1 * cos(θi) + n2 * cos(θt)))
end

#Fresnel equation for P-polarized light reflectance coefficient
function pR(n1, n2, θi, θt)
    abs((n1 * cos(θt) - n2 * cos(θi)) / (n1 * cos(θt) + n2 * cos(θi)))
end

#Fresnel equation for S-polarized light reflectance coeffiecient
function sT(n1, n2, θi, θt)
    abs((2*n1*cos(θi)) / (n1 * cos(θi) + n2 * cos(θt)))
end

#Fresnel equation for P-polarized light reflectance coefficient
function pT(n1, n2, θi, θt)
    abs((2*n1*cos(θi)) / (n1 * cos(θt) + n2 * cos(θi)))
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
function pattern(t_rate, k, λ, θi, n, l, I)
    df = t_rate
    db = t_rate .* k

    θ1 = snell_θt(n[1],n[2],θi)
    θ2 = snell_θt(n[2],n[3],θ1)
    θ3 = snell_θt(n[3],n[2],θ2)
    θ4 = snell_θt(n[2],n[1],θ3)

    r1 = (0.5)*(sR(n[1], n[2], θi, θ1) +  pR(n[1], n[2], θi, θ1))
    r2 = (0.5)*(sR(n[2], n[3], θ1, θ2) +  pR(n[2], n[3], θ1, θ2))
    r3 = (0.5)*(sR(n[3], n[2], θ2, θ3) +  pR(n[3], n[2], θ2, θ3))
    r4 = (0.5)*(sR(n[2], n[1], θ3, θ4) +  pR(n[2], n[1], θ3, θ4))

    t1 = (0.5)*(sT(n[1], n[2], θi, θ1) +  pT(n[1], n[2], θi, θ1))
    t2 = (0.5)*(sT(n[2], n[3], θ1, θ2) +  pT(n[2], n[3], θ1, θ2))
    t3 = (0.5)*(sT(n[3], n[2], θ2, θ3) +  pT(n[3], n[2], θ2, θ3))
    t4 = (0.5)*(sT(n[2], n[1], θ3, θ4) +  pT(n[2], n[1], θ3, θ4))

    β_front = δ(n[2], df, θ1, λ)
    γ       = δ(n[3], l, θ2, λ)
    β_back  = δ(n[2], db, θ3, λ)
    T1      = T(t1, t2, r1, r2, β_front)
    R1      = R(r1, r2, β_front)
    R2      = R(r3, r4, β_back)

    [(R(R1,(T1 .* R2),γ) .*I), (db .+ df), (R1 .* I)]
end

Temp = 10
M    = .01801528
α    = 1
ρ    = 6.6661e-5
k1   = (1/2.5)
k2   = (1/k1)
d    = 2e-3
λ1   = 405e-9
λ2   = 532e-9
θi   = 0.785398
n    = [1, 1.309, 1.527]
h2oV = 2.99e-29
time = 1:5500
I    = 20000

#Particles/Molecules per second times particle/molecule size = thickness rate
Φ = sticking(α, ρ, M, Temp) * h2oV


ice  = Φ .* time
ice1 = Φ .* time
ice2 = 0:0.1e-9:500e-9

p1 = pattern(ice, k1, λ1, θi, n, d, I)
p2 = pattern(ice1, k2, λ2, θi, n, d, I)


plot(layout=(3,1))
plot!(p1[2], real(p1[1]), label="Front Laser", subplot=1)
plot!(p2[2], real(p2[1]), label="Back Laser", color=:green,  subplot=2)
plot!(p1[2], real(p1[3]), label="R1 Only", color=:red, subplot=3)
title!("Theorectical Laser Interferance Pattern", subplot=1)
ylabel!("Laser Intensity (arbitrary units)", subplot=2)
xlabel!("Ice Thinkness (m)", subplot=3)
