using Plots

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

#Wave impedence equation
function Z(n)
    Z0 = 376.730313 #Wave Impedence of Free Space
    (Z0 / n)
end

#Fresnel equation for S-polarized light reflectance
function sR(n1, n2, θi, θt)
    abs((n1 * cos(θi) - n2 * cos(θt)) / (n1 * cos(θi) + n2 * cos(θt)))^2
end

#Fresnel equation for P-polarized light reflectance
function pR(n1, n2, θi, θt)
    abs((n1 * cos(θt) - n2 * cos(θi)) / (n1 * cos(θt) + n2 * cos(θi)))^2
end

#

#The sinesoidal equation that when plotted over time gives the interferance
#pattern graph
function pattern(t_rate, k, λ, θi, n, l, I)
    R = []
    ϕ1 = snell_θt(n[1], n[2], θi)
    ϕ2 = snell_θt(n[2], n[3], ϕ1)


    #find all reflection going forward until hit window
    for i = 1:2
        θt = snell_θt(n[i], n[i+1], θi)
        append!(R, 0.5*(sR(n[i], n[i+1], θi, θt) + pR(n[i], n[i+1], θi, θt)))
        θi = θt
        if i + 2 < 3
            θt = snell_θt(n[i+1], n[i+2], θi)
        end
    end

    #find all reflection going forward until back to vacuum
    for j = 1:2
        i = 3 - (j - 1)
        θt = snell_θt(n[i], n[i-1], θi)
        append!(R, 0.5*(sR(n[i], n[i-1], θi, θt) + pR(n[i], n[i-1], θi, θt)))
        θi = θt
        if i - 2 >= 1
            θt = snell_θt(n[i-1], n[i-2], θi)
        end
    end


    #Adjust coefficient transmission/reflection conservation
    a = R[1] * I
    b = (1-R[1])*R[2]*(1-R[4]) * I
    c = (1-R[1])*(1-R[2])*R[3]*(1-R[3])*(1-R[4]) * I
    d = (1-R[1])*(1-R[2])*(1-R[3])*R[4]*(1-R[2])*(1-R[3])*(1-R[4]) * I


    #Convert growth rate to angle change rate.
    θ1 =  (n[2] * 4 * π) / ( λ * cos(ϕ1)) .* t_rate
    θ2 = ((n[2] * 4 * π) / ( λ * cos(ϕ1)) .* t_rate) .* k
    γ  = cos((n[3] * 4 * π) / (λ * cos(ϕ2)) * l )

    println(a)
    println(c)
    println(γ)

    ( a .+ (b .* sin.(θ1)) ) .+ ( c .+ (d .* sin.(θ2)) ) .* γ
end

T    = 10
M    = .01801528
α    = 1
ρ    = 6.6661e-5
k1   = (1/6)
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
Φ = sticking(α, ρ, M, T) * h2oV


ice  = Φ .* time
ice1 = Φ .* time
ice2 = 0:0.1e-9:500e-9

y1 = pattern(ice, k1, λ1, θi, n, d, I)
y2 = pattern(ice1, k2, λ2, θi, n, d, I)


plot(layout=(2,1))
plot!(ice2, y1, label="Front Laser", subplot=1)
plot!(ice2, y2, label="Back Laser", color=:green,  subplot=2)
title!("Theorectical Laser Interferance Pattern", subplot=1)
ylabel!("Laser Intensity (arbitrary units)", subplot=1)
xlabel!("Ice Thinkness (m)", subplot=2)
