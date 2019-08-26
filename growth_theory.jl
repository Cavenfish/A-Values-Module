using Plots

#Hertz–Knudsen equation
function sticking(α, ρ, M, T)
    R   = 8.31446261815 #Gas Constant
    Avo = 6.02214076e23 #Avogadors Number
    (α*ρ*Avo) / √(2*π*M*R*T)
end

#Snell's Law solved for θ
function snell_θt(n1, n2, θi)
    asin( (n1/n2)*sin(θi) )
end

#Wave impedence equation
function Z(n)
    Z0 = 376.730313 #Wave Impedence of Free Space
    (Z0/n)
end

#Fresnel equation for S-polarized light reflectance
function sR(n1, n2, θi, θt)
    abs( (n2*cos(θi) - n1*cos(θt)) / (n2*cos(θi) + n1*cos(θt)) )^2
end

#Fresnel equation for P-polarized light reflectance
function pR(n1, n2, θi, θt)
    abs( (n2*cos(θt) - n1*cos(θi)) / (n2*cos(θt) + n1*cos(θi)) )^2
end

#The sinesoidal equation that when plotted over time gives the interferance
#pattern graph
function pattern(t_rate, c, λ, θi, n)
    R  = []
    ϕ  = snell_θt(n[1], n[2], θi)

    #find all reflection going forward until hit window
    for i = 1:2
        θt = snell_θt(n[i],n[i+1],θi)
        append!(R, (sR(n[i],n[i+1],θi,θt) + pR(n[i],n[i+1],θi,θt)))
        θi = θt
        if i+2 < 3
            θt = snell_θt(n[i+1],n[i+2],θi)
        end
    end

    #find all reflection going forward until back to vacuum
    for j = 1:2
        i = 3 - (j-1)
        θt = snell_θt(n[i],n[i-1],θi)
        append!(R, (sR(n[i],n[i-1],θi,θt) + pR(n[i],n[i-1],θi,θt)))
        θi = θt
        if i-2 >= 1
            θt = snell_θt(n[i-1],n[i-2],θi)
        end
    end

    #Adjust R3 and R4 to account for transmission/reflection conservation
    R[3] = (1-R[3])*R[3]  +  (1-R[4])*R[3]
    R[4] = (1-R[2])*R[4]  +  (1-R[3])*R[4]  +  (1-R[4])*R[4]

    #Convert growth rate to angle change rate
    θ1 =  (4*π)/(λ*cos(ϕ)) .* t_rate
    θ2 = ((4*π)/(λ*cos(ϕ)) .* t_rate) .* c

    ( R[1] .+ (R[2] .* cos.(θ1)) ) .+ ( R[3] .+ (R[4] .* cos.(θ2)) )
end

T    = 10
M    = .01801528
α    = 1
ρ    = 6.6661e-5
size = 2.5e-10
c    = 0.25
λ    = 405e-9
θi   = 0.785398
n    = [1, 1.309, 1.527]
time = 1:0.1:100

#Particles/Molecules per second times particle/molecule size = thickness rate
t_rate = (sticking(α, ρ, M, T) * size) .* time


y = pattern(t_rate, c, λ, θi, n)


plot(time,y, ylims=(-0.2, 0.3))
