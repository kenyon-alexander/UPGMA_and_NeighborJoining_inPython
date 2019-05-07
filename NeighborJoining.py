#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import string
import sys
import ast

'''-----------------------------------------------------------------------------'''
'''                              NEIGHBOR JOINING                               '''
'''-----------------------------------------------------------------------------'''

def creerMatriceQ(matriceM, nb_especes, liste_de_especes):
    '''Arguments:
            1. matriceM, la matrice actuelle qui a toutes les espèces et groupes
               matriceM, the current matrix which has all of the species and groups
            2. nb_especes, la nombre d'espèces actuellement dans la matrice
               nb_especies, the number of species currently in the matrix
            3. liste_de_especes, la liste des espèces inclues dans la matrice actuellement
               list_de_especes, the list of species included in the current matrix.
       Outputs:
            1. la matrice Q, qui est une des étapes de chaque itération de neighbor_joining
               the matrix Q, which is one of the steps of each iteration of neighbor_joining
       Purpose:
            Construire la matrice Q pour neighbor_joining
            Construct the matrix Q for neighbor_joining.
    '''

    Q = [[0 for x in range(nb_especes+1)] for y in range(nb_especes+1)]
    n = nb_especes
    d = matriceM

    for i in range(1,nb_especes+1):
        Q[0][i] = liste_de_especes[i-1]
        Q[i][0] = liste_de_especes[i-1]

    for i in range(1, nb_especes+1):
        for j in range(1, i):
            #sauter les diagonales de la matrice:
            #skip the diagonals of the matrix.
            if(i != j):

                #Ici, on calcule chaque valeur de Q selon la formule suivante:
                #Here, we calculate each value of Q according to the following
                    #formula:
                #                            n            n
                #                           ----         ----
                # Q(a,b) = (n-2) * d(a,b) - \   d(a,k) - \   d(b,k)
                #                           /            /
                #                           ----         ----
                #                           k=1          k=1

                #Ici on calcule les deux sommations dans la formule:
                #Here we calculate the two summations in the formula:
                sumDak = 0
                sumDbk = 0
                k = 1
                while d[i][k] != 0:
                    sumDak += d[i][k]
                    k += 1
                while k < nb_especes+1:
                    sumDak += d[k][i]
                    k += 1
                k = 1
                while d[j][k] != 0:
                    sumDbk += d[j][k]
                    k += 1
                while k < nb_especes+1:
                    sumDbk += d[k][j]
                    k+= 1

                #Maintenant on applique la formule:
                #Now we apply the formula.
                Q[i][j] = ((n-2) * d[i][j]) - sumDak - sumDbk

    return Q

def ajouterNoeud_NJ(arbre, filsGauche, filsDroit, indexFilsGauche, indexFilsDroit, matriceD):
    '''Arguments:
            1. arbre, l'arbre auquel on veut ajouter un nouveau noeud
               arbre, the tree to which we want to add a new node.
            2. filsGauche, un string qui a l'espèce nom d'un de noeuds qu'on veut connecter
               filsGuahce, a string which has the species name of one of the nodes that we
                    want to connect.
            3. filsGauche, un string qui a l'espèce nom de l'autre noeud qu'on veut connecter
               filsGauche, a string which has the species name of the other node that we
                    want to connect.
            4. indexFilsGauche, l'index du filsGauche dans la matrice D
               indexFilsGauche, the index of filsGauche in the matrix D
            5. indexFilsDroit, l'index du filsDroit dans la matrice D
               indexFilsDroit, the index of filsDroit in the matrix D
            6. matriceD, la matrice des distances entre tous les groupes.
               matriceD, the matrix of distances between all of the groups.
       Output:
            1. un arbre qui est l'arbre d'input, avec un nouveau noeud qui connecte les deux
                noeuds inputted.
               a tree which is the input tree, with a node which connects the two nodes inputted.
       Purpose:
            Connecter deux noeuds dans l'arbre actuel par créer un nouveau noeud comme parent
                pour les deux.
            Connect two nodes in the current tree by creating a new node as a parent of both.
    '''
    #Pour ajouter un noeud à l'arbre dans l'algorithme de neighbor-joining, on a besoin
        #de les formules suivantes pour calculer la distance entre le nouveau noeud et ses fils:
    #To add a new node to the tree in neighbor-joining, we need the following formula to calculate
        #the distance between the new node and its children:
        #                                 __  n             n       __
        #               d(a,b)      1    |   ----          ----       |
        # delta(a,u) = -------- + ------ |   \   d(a,k) -  \   d(b,k) |
        #                 2       2(n-2) |   /             /          |
        #                                |   ----          ----       |
        #                                 --  k=1           k=1     --
        #
        # delta(b,u) = d(a,b) - delta(a,u)

    n = len(matriceD) - 1

    cardinaux = filsGauche.cardinaux + filsDroit.cardinaux

    sumDak = 0
    sumDbk = 0
    k = 1
    while matriceD[indexFilsGauche][k] != 0:
        sumDak += matriceD[indexFilsGauche][k]
        k += 1
    while k < n+1:
        sumDak += matriceD[k][indexFilsGauche]
        k += 1
    k = 1
    while matriceD[indexFilsDroit][k] != 0:
        sumDbk += matriceD[indexFilsDroit][k]
        k += 1
    while k < n+1:
        sumDbk += matriceD[k][indexFilsDroit]
        k+= 1

    lgFilsGauche = (matriceD[indexFilsGauche][indexFilsDroit]/2) + ((1/(2*(n-2)))*(sumDbk - sumDak))
    lgFilsDroit = matriceD[indexFilsGauche][indexFilsDroit] - lgFilsGauche

    #On n'a pas besoin de cardinaux ni de height dans neighbor_joining, donc on utilise "None" ici.
        #Pour qu'on n'aie pas besoin d'une nouvelle classe de myNoeud pour neighbor-joining.
    #We don't need the cardinality nor the height in neighbor_joining, so we use "None" here.
        #So that we don't need a new class of myNoeud for neighbor_joining.
    arbre["(" + filsGauche.sp + "," + filsDroit.sp + ")"] = myNoeud("(" + filsGauche.sp + "," + filsDroit.sp + ")", filsGauche, filsDroit, lgFilsGauche, lgFilsDroit, cardinaux, None)

    return [arbre, arbre["(" + filsGauche.sp + "," + filsDroit.sp + ")"]]

def MAJ_matrice_NJ(groupes, prevMatrice, espece1index, espece2index, nouveauNoeud):
    '''Arguments:
            1. groupes, une dictionnaire des groupes dans l'algorithme actuel.
               groupes, a dictionary of the groups in the current iteration.
            2. prevMatrice, la matrice D qui va être "mise à jour"
               prevMatrix, the matrix D which will be "updated"
            3. espece1index, l'index d'une des espèces qu'on veut supprimer de la matrice
               espece1index, the index of one of the species that we want ot remove from
                    the matrix
            4. espece2index, l'index de l'autre espèce qu'on veut supprimer de la matrice
               espece2index, the index of the other species that we want to remove from
                    the matrix
            5. nouveauNoeud, le noeud qui était juste ajouté à l'arbre.
               nouveauNoeud, the node which was just added to the tree.
       Output:
            1. une nouvelle matriceD, mise à jour avec toute les distances entre tous les
                groupes actuels.
               a new matriceD, updated with all of the distances between all of the
                    current groups.
       Purpose:
            Mettre à jour la matrice D.
            Update matrix D.
    '''

    liste_de_especes = list(groupes.keys())
    nb_especes = len(liste_de_especes)
    matrice = [[0 for x in range(nb_especes+1)] for y in range(nb_especes+1)]

    #creer la nouvelle matrice
    #create the new matrix.
    for i in range(1,len(liste_de_especes)+1):
        matrice[0][i] = liste_de_especes[i-1]
        matrice[i][0] = liste_de_especes[i-1]

    #établir les noms des espèces qu'on veut supprimer de la matrice:
    #establish the names of the species that we want ot remove from the matrix:
    espece1 = nouveauNoeud.filsGauche.sp
    espece2 = nouveauNoeud.filsDroit.sp

    #remplir la matrice avec les valeurs de la matrice précédente,
        #mais sauter les espèces jointes dans cette itération.
    #fill the matrix with the values of the previous matrix, but skip the
        #joined species from this iteration.
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

    dernierRang = nb_especes

    #itérer au-dessus du dernier rang de la nouvelle matrice et ajouter les valeurs
        #selon la formule:
    #Iterate over the last row of the new matrix and add values according to the
        #following formula:
                #           1
                # d(u,k) = --- [d(f,k) + d(g,k) - d(f,g)]
                #           2

    prevCol = 1
    for col in range(1, nb_especes):
        while(prevCol == espece1index or prevCol == espece2index):
            prevCol += 1
        dfk = max(prevMatrice[espece1index][prevCol], prevMatrice[prevCol][espece1index])
        dgk = max(prevMatrice[espece2index][prevCol], prevMatrice[prevCol][espece2index])
        matrice[dernierRang][col] = .5*(dfk + dgk - prevMatrice[espece1index][espece2index])
        prevCol += 1

    return matrice

def connecter_noeuds(arbre, sp1, sp2, distance):
    '''Arguments:
            1. arbre, l'arbre dans lequel on veut connecter deux noeuds
               arbre, the tree in which we want to connect two nodes.
            2. sp1, un string avec le nom d'espèce d'un des noeuds qu'on veut connecter
               sp1, a string with the species name of one of the nodes we want to connect
            3. sp2, un string avec le nom d'espèce de l'autre noeud qu'on veut connecter
               sp2, a string wit the species name of the other node we want to connect
            4. distance, la distance entre les deux noeuds qu'on veut connecter.
               distance, the distance between the two nodes that we want to connect.
       Output:
            1. la même arbre d'input, mais avec les deux noeuds connectés.
               the same input tree, but with the two nodes connected.
       Purpose:
            Parfois, on ne veut pas ajouter un autre noeud, parce que les deux noeuds
                existent déjà. Pourtant, on veut quand même connecter deux noeuds,
                donc cette fonction existe pour cette raison.
           Sometimes, we don't want to add another node, because the two nodes already
                exist. However, we want to connect the two nodes, so this function
                exists for this reason.
    '''

    if arbre[sp1].filsGauche and arbre[sp1].filsDroit:
        if arbre[sp2].filsDroit:
            arbre[sp2].filsGauche = arbre[sp1]
            arbre[sp2].lgFilsGauche = distance
        else:
            arbre[sp2].filsDroit = arbre[sp1]
            arbre[sp2].lgFilsDroit = arbre[sp1]
    else:
        if arbre[sp1].filsDroit:
            arbre[sp1].filsGauche = arbre[sp2]
            arbre[sp1].lgFilsGauche = distance
        else:
            arbre[sp1].filsDroit = arbre[sp2]
            arbre[sp1].lgFilsDroit = distance

    return arbre

def trouver_racine(arbre):
    '''Arguments:
            1. arbre, l'arbre dans lequel on veut trouver la racine
               arbre, the tree in which we want to find the root.
       Output:
            1. racine, la racine de l'arbre d'input
               racine, the root of the input tree.
       Purpose:
            Dans neighbor_joining, il n'est pas évidant où est la racine de l'arbre
                à cause de la construction de l'arbre créé par l'algorithme. Donc
                cette fonction existe pour cette raison. La "racine" dans ce cas
                est le noeud dans l'arbre avec le cardinaux le plus grand.
            In neighbor_joining, it's not obvious where the root of the tree is
                because of the construction of the tree created by the algorithm.
                So this function exists for that reason. The "racine" in this case
                is the node in the tree with the greatest cardinality.
    '''
    maxCard = 0
    racine = None
    for noeud in arbre:
        if arbre[noeud].cardinaux > maxCard:
            maxCard = arbre[noeud].cardinaux
            racine = arbre[noeud]

    return racine

def neighbor_joining(matriceListe, liste_de_especes):
    '''Arguments:
            1. matriceListe, une matrice dans la forme d'une liste dans Python,
                qui a toutes les distances entre toutes les espèces qu'on veut utiliser
                pour construire l'arbre phylogénétique.
               matriceList, a matrix in the form of a list in Python, which as all of
                    the distances between all of the species that we want to use to
                    construct the phylogenetic tree.
            2. liste_de_especes, une liste facultatif qui à tous les noms de toutes les
                espèces inclues. Si liste_de_especes est None, l'algorithme construira
                une liste soi-même.
               list_de_especes, an optional list which has all of the names of all of the
                    included species. If the list_de_especes is None, the algorithm will
                    construct a list anyway.
       Output:
            1. Un string qui construit un arbre dans le format Newick.
               A string which constructs the phylogenetic tree in Newick format.
       Purpose:
            Construire un arbre phylogénétique donné une matrice de distances entre
                n espèces. Cet algorithme est une version de l'algorithme de Neighbor-Joining.
            Construct a phylogenetic tree given a matrix of distances between n species.
                This algorithm is a version of the Neighbor-Joining algorithm.
    '''

    #CREER LA MATRICE PRELIMINAIRE
        #d'abord, on créer le variable nEspeces pour qu'on puisse changer le nombre
        #d'espèces
    #CREATE THE PRELIMINARY MATRIX
        #First, we create the variable nEspeces so that we can change the number of
        #species.
    nEspeces = len(matriceListe)

    #Il faut assurer que neighbor_joining peut soutenir n'importe quel nombre d'espèces,
        #donc ici on a une boucle qui utilise l'alphabet de nommer les especes.
        #Si le nombre d'espèces est >26, on commence à utiliser 'AA', 'BB', et
        #'so on'.
    #We need to make sure that neighbor_joining can support any number of species,
        #so here we have a loop which uses the alphabet to name the species. If the
        #number of species is greater than 26, we start using 'AA', 'BB', and so on.
    if(liste_de_especes == None):
        liste_de_especes = creerListe_de_especes(nEspeces)
    elif(len(matriceListe) != len(liste_de_especes)):
        sys.exit("ERREUR neighbor_joining: matriceListe et liste_de_especes ne sont pas la même taille!")

    monArbre = creerArbre(nEspeces, liste_de_especes)
    maMatrice = creerMatrice(matriceListe, liste_de_especes, nEspeces)

    #creer une dictionnaire des groupes pour la mise a jour de la matrice
    #create a dictionary of the groups for the matrix update.
    groupes = {}
    for i in range(len(liste_de_especes)):
        groupes[liste_de_especes[i]] = liste_de_especes[i]

    while len(maMatrice) >= 3:

        if(len(maMatrice) == 3):
            monArbre = connecter_noeuds(monArbre, maMatrice[0][1], maMatrice[0][2], maMatrice[2][1])
            break
        else:
            matriceQ = creerMatriceQ(maMatrice, nEspeces, list(groupes.keys()))

            minList = trouverMin(matriceQ, nEspeces)
            minimum = minList[0]
            espece1 = maMatrice[0][minList[1]]
            espece2 = maMatrice[0][minList[2]]
            espece1index = minList[1]
            espece2index = minList[2]
            #supprimer les especes minimums de la dictionnaire "groupes" et ajouter une autre cle qui est
                #la combinaison des deux
            #remove the minimum species from the dictionary "groupes" and add a new key which is a
                #combination of the two species.
            del groupes[espece1]
            del groupes[espece2]
            groupes["(" + espece1 + "," + espece2 + ")"] = "(" + espece1 + "," + espece2 + ")"
            nEspeces -= 1

            monArbreAvecNouveauNoeud = ajouterNoeud_NJ(monArbre, monArbre[espece1], monArbre[espece2], espece1index, espece2index, maMatrice)
            monArbre = monArbreAvecNouveauNoeud[0]
            nouveauNoeud = monArbreAvecNouveauNoeud[1]

            maMatrice = MAJ_matrice_NJ(groupes, maMatrice, espece1index, espece2index, nouveauNoeud)

    racine = trouver_racine(monArbre)

    print(formatNewick(racine, None))

if __name__ == "__main__":
    inputMatrix = ast.literal_eval( sys.argv[1] )
    liste_de_especes = sys.argv[2]
    if sys.argv[2].strip()=='None':
        liste_de_especes = None
    neighbor_joining(inputMatrix,liste_de_especes)
