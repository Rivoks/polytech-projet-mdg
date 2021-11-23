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
# On utilise une boucle "for" pour cela (on prend les valeurs au dessus de S = 500
for i = 30:96
    if (S[:, 80][i] < 500)
        display(i)
        break
    end
end
# On obtient l'intervalle [21;87]
xs = 21:1:87

# On va préparer ici nos pramètres pour l'interpolation
# On définit une variable temporaire correspondant à la colonne 80 de S
tmp = S[:, 80]

# On va créer un masque pour tmp entre 21 et 87 
for i = 1:96
    if (i < 21 || i > 87)
        tmp[i] = 0 * tmp[i]
    end
end

# display(tmp)
# display(S[:, 80]);
# display(xs);
# A = [S[:, 80] for x in xs]
# display(A[1])

# On définit un intervalle sur lequel interpoler (ici tout l'intervalle)
xq = 1:1:96

# On lance la fonction d'interpolation (issue d'un package)
interp_linear = LinearInterpolation(xq, tmp)

# Obtient la fonction d'interpollation pour nos valeurs
# Ainsi il suffit d'entrer un paramètre un floatant pour obtenir la valeur interpollée
interp_linear(50)
# D'ailleurs le package nous permet d'avoir des approximations entre chaque valeurs
interp_linear(50.5)

# Notre but va désormais être de créer notre propre fonction d'interpolation
# À partir d'un vecteur et d'un intervalle on retournera un nouveau vecteur interpollé, 
# avec des approximations entre chaques valeurs

function interp(x, rg, new_rg)
    tmp = x
    res = zeros(length(new_rg))
    # Masque...
    for i = 1:1:length(rg)
        if (i < 21 || i > 87)
            tmp[i] = 0 * tmp[i]
        end
    end

    # Package d'interpolation...
    inter_lin = LinearInterpolation(rg, tmp)

    # On assigne la valeur interpollée à chaque coordonnées du vecteur...
    for j = 1:1:length(new_rg)
        res[j] = inter_lin(new_rg[j])
    end

    # Enfin on retourne notre nouveau vecteur interpollé...
    return res

end

# On définit l'intervalle issu de xc
display(xc);
display(xc[1]);
display(xc[end]);
display(xc[86] - xc[85]);

# On obtient notre nouvel intervalle
rg_xc = -627725:15600:854275


test = interp(S[:, 80], rg_xc, -627725:7800:854275)


