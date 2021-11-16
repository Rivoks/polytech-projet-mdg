using Plots, Printf, LinearAlgebra, GeoData, NCDatasets, JLD, Interpolations

# On charge notre fichier GeoData...
# Puis on extrait les valeurs dans les variables...
var = load("BedMachineGreenland_96_176.jld"); # ultra low res data
Zbed = var["Zbed"]
xc = var["xc"]
Hice = var["Hice"]


S = zeros(96, 176) # Épaisseur de de la glace
B = Zbed.data # Niveau de la mer
H = Hice.data # Hauteur de la glace
S .= B .+ H


# display(B)
# display(B[1])

# Ici on plot les deuxs niveaux en même temps pour pouvoir
# visuellement définir un intervalle d'inteprolation
plot(xc, B[:, 80], color = :blue)
plot!(xc, S[:, 80], color = :red)


# Il faut interpoler entre -3*10^5 et 6*10^5
# Pour cela on on essaye de déterminer les indices dans la matrice...
# On utilise une boucle "for" pour cela (on prend les valeurs au dessus de S = 500)
for i = 30:96
    if (S[:, 80][i] < 500)
        display(i)
        break
    end
end

# On obtient l'intervalle [21;87]
# On lance la fonction d'interpolation (issue d'un package)
xc = 21:1:87
interp_linear = LinearInterpolation(xs, S)








#using Plots, JLD , GeoData ,LinearAlgebra, NCDatasets
#include(joinpath(@__DIR__, "helpers.jl"))




#d = load("BedMachineGreenland_96_176.jld")
#Zbed = d["Zbed"]

#B = Zbed.data
#display(B); 
# load the data
#print("Loading the data ... ")
#Zbed, Hice, Mask, dx, dy, xc, yc = load_data(; nx=96) # nx=96,160 are included in the repo
#                                                      # other numbers will trigger a 2GB download
#println("done.")

# Faire xc * 176 // yc * 96 