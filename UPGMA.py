#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import string
import sys
import ast
from FonctionsGenerales import myNoeud, creerListe_de_especes, creerMatrice, trouverMin, creerArbre, formatNewick

'''-----------------------------------------------------------------------------'''
'''                                    UPGMA                                    '''
'''-----------------------------------------------------------------------------'''

def ajouterNoeud(arbre, filsGauche, filsDroit, minDistance):
    '''Arguments:
            1. arbre, qui est l'arbre actuel auquel on veut ajouter un nouveau noeud
            2. filsGauche, le noeud qui va devenir le filsGauche du nouveau noeud
            3. filsDroit, le noeud qui va devenir le filsDroit du nouveau noeud
            4. minDistance, la distance entre filsGauche et filsDroit qui va être
                utilisée pour calculer la distance des edges entre le nouveau noeud
                et ses fils.
       Outputs:
            1. l'arbre avec le nouveau noeud et toute son information
       Purpose:
            Ajouter un nouveau noeud à l'arbre et y attacher toute l'information
                nécessaire.
    '''

    cardinaux = filsGauche.cardinaux + filsDroit.cardinaux

    #Il y a quelque chose assez bizarre ici, car il faut utiliser le .height de
        #fils droit pour calculer le lgFilsGauche et vice versa. Je ne sais pas
        #vraiment pourquoi, mais il marche donc...
    lgFilsGauche = (minDistance/2) - filsDroit.height
    lgFilsDroit = (minDistance/2) - filsGauche.height

    height = filsGauche.height + lgFilsDroit

    arbre["(" + filsGauche.sp + "," + filsDroit.sp + ")"] = myNoeud("(" + filsGauche.sp + "," + filsDroit.sp + ")", filsGauche, filsDroit, lgFilsGauche, lgFilsDroit, cardinaux, height)

    return [arbre, arbre["(" + filsGauche.sp + "," + filsDroit.sp + ")"]]

def MAJ_matrice(groupes, prevMatrice, espece1index, espece2index, nouveauNoeud):
    '''Arguments:
            1. groupes, la dictionnaire de toute les groupes inclus actuellement
                dans l'algorithme. 'groupes' existe seulement comme variable locale
                dans la fonctionne UPGMA
            2. prevMatrice, la matrice de la dernière étape qui va être utilisée
                pour construire la nouvelle matrice
            3. espece1index, l'indice de la première espèce qu'on veux supprimer de
                la matrice
            4. espece2index, l'indice de la deuxième espèce qu'on veux supprimer de
                la matrice
            5. nouveauNoeud, le noued qu'on a ajouté dans la fonctionne ajouterNoeud,
                avec toute l'information nécessaire du nouveau groupe qu'on veux
                ajouter à la nouvelle matrice.
       Output:
            1. la nouvelle matrice avec les distances de tous les groupes inclus
                actuellement dans l'algorithme.
       Purpose:
            Mettre à jour la matrice, ou plus correctement créer une nouvelle matrice
                qui guarde toutes les distances entre tous les groupes.
    '''

    liste_de_especes = list(groupes.keys())
    nb_especes = len(liste_de_especes)
    matrice = [[0 for x in range(nb_especes+1)] for y in range(nb_especes+1)]

    #creer la nouvelle matrice
    for i in range(1,len(liste_de_especes)+1):
        matrice[0][i] = liste_de_especes[i-1]
        matrice[i][0] = liste_de_especes[i-1]

    #établir les noms des espèces qu'on veut supprimer de la matrice:
    espece1 = nouveauNoeud.filsGauche.sp
    espece2 = nouveauNoeud.filsDroit.sp

    #remplir la matrice avec les valeurs de la matrice précédente,
        #mais sauter les espèces jointes dans cette itération.
    prevI = 2
    prevJ = 1
    for i in range(2, nb_especes):
        for j in range(1,i):
            while(prevI == espece1index or prevI == espece2index):
                prevI+=1
            while(prevJ == espece1index or prevJ == espece2index):
                prevJ+=1
                while(prevMatrice[prevI][prevJ] == 0 or prevI == espece1index or prevI == espece2index):
                    prevI+=1
            matrice[i][j] = prevMatrice[prevI][prevJ]
            prevJ += 1
        prevI += 1
        prevJ = 1


    #calculer les distances entre le nouveau noeud (A,B) et tous les autres:

    cardA = nouveauNoeud.filsGauche.cardinaux
    cardB = nouveauNoeud.filsDroit.cardinaux

    dernierRang = nb_especes

    #itérer au-dessus du dernier rang de la nouvelle matrice et ajouter les valeurs
        #selon la formule:
                # d[(A U B), X] = |A|*D(A,X) + |B|*D(B,X)
                #                 -----------------------
                #                         |A| + |B|
                #
                # soit "D" est prevMatrice et "d" est matrice

    prevCol = 1
    for col in range(1, nb_especes):
        while(prevCol == espece1index or prevCol == espece2index):
            prevCol += 1
        if(prevMatrice[0][prevCol] == nouveauNoeud.filsGauche.sp):
            cardA = nouveauNoeud.filsGauche.cardinaux
            cardB = nouveauNoeud.filsDroit.cardinaux
        else:
            cardB = nouveauNoeud.filsGauche.cardinaux
            cardA = nouveauNoeud.filsDroit.cardinaux
        prevValeur1 = max(prevMatrice[espece1index][prevCol], prevMatrice[prevCol][espece1index])
        prevValeur2 = max(prevMatrice[espece2index][prevCol], prevMatrice[prevCol][espece2index])
        matrice[dernierRang][col] = float(((cardA * prevValeur1) +
                                    (cardB * prevValeur2)) /
                                    (cardA + cardB))
        prevCol += 1

    return matrice

def UPGMA(matriceListe, liste_de_especes):
    '''Arguments:
            1. matriceListe, une matrice dans la forme d'une liste qui a toutes
                les différences entre les espèces en question. Cette matrice
                peut être plein, ou elle peut être que la moitié (étant donné
                qu'elle est symmetrique)
            2. liste_de_especes, un argument facultatif qui peut être une liste
                de toutes les espèces en question, ou peut être None. Dans le
                cas de None, l'algorithme construira une liste d'espèces de
                n'importe quel longeur par utiliser l'alphabet.
        Output:
            1. Un string dans le format Newick d'un arbre phylogénétique construit
                avec les données. L'arbre étiquettera les feuilles et les distances
                entre les noeuds.
        Purpose:
            Construire un arbre phylogénétique entre les espèces fournies, par
            utiliser l'algorithme de UPGMA.
    '''

    #CREER LA MATRICE PRELIMINAIRE
        #d'abord, on créer le variable nEspeces pour qu'on puisse changer le nombre
        #d'espèces
    nEspeces = len(matriceListe)

    #Il faut assurer que UPGMA peut soutenir n'importe quel nombre d'espèces,
        #donc ici on a une boucle qui utilise l'alphabet de nommer les especes.
        #Si le nombre d'espèces est >26, on commence à utiliser 'AA', 'BB', et
        #'so on'.
    if(liste_de_especes == None):
        liste_de_especes = creerListe_de_especes(nEspeces)
    elif(len(matriceListe) != len(liste_de_especes)):
        sys.exit("ERREUR UPGMA: matriceListe et liste_de_especes ne sont pas la même taille!")

    maMatrice = creerMatrice(matriceListe, liste_de_especes, nEspeces)
    monArbre = creerArbre(nEspeces, liste_de_especes)

    #creer une dictionnaire des groupes pour la mise a jour de la matrice
    groupes = {}
    for i in range(len(liste_de_especes)):
        groupes[liste_de_especes[i]] = liste_de_especes[i]

    while(len(maMatrice) >= 3):
        #"mins" devient une liste de trois choses:
            #0. le valeur minimum de la matrice
            #1. l'index de la première espèce trouvée d'avoir la distance minimale
            #2. l'index de la deuxième espèce trouvée d'avoir la distance minimale
        mins = trouverMin(maMatrice, nEspeces)
        minimum = mins[0]
        espece1 = maMatrice[0][mins[1]]
        espece2 = maMatrice[0][mins[2]]
        espece1index = mins[1]
        espece2index = mins[2]
        #supprimer les especes minimums de la dictionnaire "groupes" et ajouter une autre cle qui est
            #la combinaison des deux
        del groupes[espece1]
        del groupes[espece2]
        groupes["(" + espece1 + "," + espece2 + ")"] = "(" + espece1 + "," + espece2 + ")"
        nEspeces -= 1

        #ajouter un nouveau noeud à l'arbre
        monArbreAvecNouveauNoeud = ajouterNoeud(monArbre, monArbre[espece1], monArbre[espece2], minimum)
        monArbre = monArbreAvecNouveauNoeud[0]
        nouveauNoeud = monArbreAvecNouveauNoeud[1]

        #METTRE A JOUR LA MATRICE
        #D'abord, créer une nouveau liste de valeurs qui n'inclut pas les espèces combinés
        maMatrice = MAJ_matrice(groupes, maMatrice, espece1index, espece2index, nouveauNoeud)

    racine = maMatrice[0][1]

    print(formatNewick(monArbre[racine], None))

if __name__ == "__main__":
    # M1 = [[0,8,7,12], [8,0,9,14], [7,9,0,11], [12,14,11,0]]
    # M2 = [[0,2,3,8,14,18],[2,0,3,8,14,18],
    #       [3,3,0,8,14,18],[8,8,8,0,14,18],
    #       [14,14,14,14,0,18],[18,18,18,18,18,0]]
    # #UPGMA
    # M3 = [[0,19,27,8,33,18,13],[19,0,31,18,36,1,13],
    #           [27,31,0,26,41,32,29],[8,18,26,0,31,17,14],
    #           [33,36,41,31,0,35,28],[18,1,32,17,35,0,12],
    #           [13,13,29,14,28,12,0]]
    # #Neighbor Joining
    #M4 = [[0,2,4,6,6,8],[2,0,4,6,6,8],[4,4,0,6,6,8],[6,6,6,0,4,8],[6,6,6,4,0,8],[8,8,8,8,8,0]]
    # UPGMA(M3, None)
    # MWiki = [[0,5,9,9,8],
    #          [5,0,10,10,9],
    #          [9,10,0,8,7],
    #          [9,10,8,0,3],
    #          [8,9,7,3,0]]
    # MWiki2 = [[0,0,0,0],
    #           [8,0,0,0],
    #           [7,3,0,0],
    #           [7.0,7.0,6.0,0]]
    # liste_de_especes = ["C", "D", "E", "(B,A)"]
    # neighbor_joining(M4, None)
    inputMatrix = ast.literal_eval( sys.argv[1] )
    liste_de_especes = sys.argv[2]
    if sys.argv[2].strip()=='None':
        liste_de_especes = None 
    UPGMA(inputMatrix,liste_de_especes)
