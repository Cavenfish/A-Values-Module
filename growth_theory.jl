
#Hertz–Knudsen equation
function sticking(α, ρ, M, T)
    R   = 8.31446261815 #Gas Constant
    Avo = 6.02214076e23 #Avogadors Number
    (α*ρ*Avo) / √(2*π*M*R*T)
end

#Snell's Law solved for θ
function snell_θ2(n1, n2, θ1)
    asin( (n1/n2)*sin(θ1) )
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
function pattern(t_rate, R1, R2, R3, R4)
    θ1 =
    θ2 =

    ( R1 + R2*cos(θ1) ) - ( R3 + R4*cos(θ2) )
end

#Particles/Molecules per second times particle/molecule size = thickness rate
t_rate = sticking(α, ρ, M, T) * size
