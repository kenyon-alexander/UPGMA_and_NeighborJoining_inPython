#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import string
import sys

'''-----------------------------------------------------------------------------'''
'''                          FONCTIONNES GENERALES                              '''
'''-----------------------------------------------------------------------------'''

class myNoeud:

    # Initializer / Instance Attributes
    #Initialize instance attributes
    def __init__(self, sp, filsDroit, filsGauche, lgFilsGauche, lgFilsDroit, cardinaux, height):
        self.sp = sp
        self.filsDroit = filsDroit
        self.filsGauche = filsGauche
        self.lgFilsDroit = lgFilsDroit
        self.lgFilsGauche = lgFilsGauche
        self.cardinaux = cardinaux
        self.height = height

def creerListe_de_especes(nb_especes):

    if(nb_especes == 0):
        sys.exit("ERREUR: creerListe_de_especes: nb_especes ne peut pas être 0!")

    liste_de_especes = []
    letterWrapper = 1
    compte = 0
    for i in range(nb_especes):
        if(compte == 26):
            letterWrapper += 1
            compte = 0
        newEspece = string.ascii_uppercase[compte] * letterWrapper
        liste_de_especes.append(newEspece)
        compte += 1

    return liste_de_especes

def creerMatrice(matriceListe, liste_de_especes, nb_especes):
    '''Arguments:
            1. matriceListe, une matrice de valeurs de différence entre especes
               matriceList, a matrix of difference values between species.
            2. liste_de_especes, un argument facultatif qui fournit les noms d'espèces
               liste_de_especes, an optional argument which provides species names.
            3. nb_especes, le nombre d'especes qu'on veut voir dans l'arbre phylogénétique.
               nb_especes, the number of species that we want to see in the phylogenetic
                    tree.
       Return:
            1. Une matrice avec tous les valeurs de différence entre espèces et les
                noms d'espèces
               A matrix with all of the difference values between species and the
                    species names.
            2. liste_de_especes. Si l'argument liste_de_especes est None, cette
                fonctionne construira une liste de caractères itératives de ASCII.
               list_de_especes. Si l'argument liste_de_especes est None, this function
                    will constrcut a list of iterative ASCII characters.
       Purpose:
            Créer la matrice originale pour commencer l'algorithme
            Create the original matrix to start the algorithm.
    '''

    #Il n'est pas nécessaire d'inclure le nombre d'espèces comme argument car
        #ce nombre est déjà implicite dans la matriceListe donnée.
    #It's not necessary to include the number of species as an argument because
        #this number is already implicit in the matriceList given.
    if(nb_especes < 2):
        sys.exit("ERREUR creerMatrice: UPGMA ne marche pas pour moins de 2 espèces!")

    #La matrice créée doit être une dimension plus grande que la liste donnée,
        #parce qu'il faut guarder les noms d'espèces dans le premier colomne et rang
    #The create matrix must be one dimension greater than the give list, because
        #we need to store the species names in the first row and column.
    matrice = [[0 for x in range(nb_especes+1)] for y in range(nb_especes+1)]

    #Remplir la matrice avec les noms d'espèces, quelque soit ou pas les noms
        #d'espèces sont fournis
    #Fill the matrix with the species names, whether or not the names are provided.
    for i in range(1,nb_especes+1):
        matrice[0][i] = liste_de_especes[i-1]
        matrice[i][0] = liste_de_especes[i-1]

    #Maintenant on remplit la matrice avec les valeurs dans la matriceListe.
    #Now we fill them matrix with the values in matriceList.
    for i in range(1, nb_especes+1):
        for j in range(1, nb_especes+1):
            matrice[i][j] = matriceListe[i-1][j-1]

    return matrice

def trouverMin(matrice, nb_especes):
    '''Arguments:
            1. matrice, la matrice d'espèces (ou groupes d'espèces) avec leurs
                différences.
               matrice, the matrix of species (or groups of species) with their
                differences.
            2. nb_especes, le nombre d'espèces actuellement inclues dans la matrice
               nb_especies, the number of species currently included in the matrix.
       Output:
            0. la valeur minimale de la matrice
               the minimum value in the matrix
            1. l'indice du rang de la valeur minimale
               the index of the row of the minimum value
            2. l'indice du colomne de la valeur minimale
               the index of the column of the minimum value
       Purpose:
            Une des étapes de UPGMA. Trouver la valeur minimale de la matrice
                pour choisir les deux groupes qui sont les plus proches.
            One of the steps of UPGMA and neighbor_joining. Find the minimum value
                of the matrix to chose the two groups which are the closest.
    '''

    #Si le nombre d'espèces est moins de 2, on sait qu'il y avait une erreur
    #If the number of species is less than 2, we know that there was an error.
    if nb_especes < 2:
        sys.exit("ERREUR trouverMin: nb_especes moins de 2!")

    minimum = matrice[2][1]
    iMin = 1
    jMin = 0

    #L'assomption ici est que la matrice soit symmétrique, donc on n'utilise
        #que la moitié.
    #The assumption here is that the matrix is symmetrical, so we use only half
        #of it.
    for i in range(2, nb_especes+1):
        for j in range(1,i):
            if float(matrice[i][j]) <= float(minimum):
                minimum = matrice[i][j]
                iMin = i
                jMin = j

    return [minimum, iMin, jMin]

def creerArbre(nb_especes, liste_de_especes):
    '''Arguments:
            1. nb_especes, le nombre d'espèces originalement dans l'arbre
               nb_especes, the number of species originally in the tree.
            2. liste_de_especes, une liste d'espèces fournies ou créés par
                la fonctionne creerMatrice.
               list_de_especes, the list of species provided or created by the
                the function creerMatrice.
       Output:
            1. une dictionnaire d'objets 'myNoeud', qui ont leurs nom d'espèce
                comme clé et un 'myNoeud' objet comme valeur.
               a dictionary of 'myNoeud' objects, which have their species names
                    as keys and a 'myNoeud' oject as a value.
       Purpose:
            Construire qu'une fois l'arbre, auquel on peut ajouter plus de noeuds
                dans la fonctionne UPGMA.
            Construct the tree one time, to which we can add more nodes in the function
                UPGMA or neighbor_joining
    '''

    dict_de_noeuds = {}

    for i in range(nb_especes):
        dict_de_noeuds[liste_de_especes[i]] = myNoeud(liste_de_especes[i], None, None, None, None, 1, 0)

    return dict_de_noeuds

def formatNewick(racine, distToParent):
    '''Arguments:
            1. racine, la racine de l'arbre qu'on veut formatter dans la
                format Newick
               racine, the root of the tree that we want to format in Newick format
            2. distToParent, la distance à la parent, seulement utilisé pour
                les feuilles. Si la racine n'est pas une feuille, distToParent
                est None.
               distToParent, the distance to the parent, only used for if 'racine'
                    is a leaf. If the 'racine' is not a leaf, distToParent is None.
       Output:
            1. Un string de l'arbre ou sous-arbre en format Newick.
               A string of the tree or sub-tree in Newick format.
       Purpose:
            Il s'agit d'un algorithme récursif qui traverse l'arbre et construit
                un string dans le format Newick pour afficher les feuilles et
                les distances entre noeuds.
            This is a recursive algorithm which traverses the tree and constructs a
                string in Newick format to display the leaves and the distances between
                nodes.
    '''


    #Base Case:
    if(racine.cardinaux == 1):
        return (racine.sp + distToParent)

    #confusingly, the below strings are only ever used if the child in question is
        #a leaf. In other words, the above base case uses these strings, but if
        #the base case is not triggered then these strings will never be used. Instead,
        #the strings blancSiFilsDroitEstFeuille and blancSiFilsGaucheEstFeuille are
        #used in their stead, and can be essentially erased if the current root is a
        #leaf. I had to install this confusing algorithm because the node object
        #that I created does not have a "distance to parent" variable, which prevents
        #acces to the edge length from the leaf itself. There are other work-arounds
        #for this problem, including simply adding a distanceToParent variable to the
        #myNoeud class, but this algorithm was an interesting challenge and it works well.
    distFilsGauche = (":" + str(racine.lgFilsGauche))
    distFilsDroit = (":" + str(racine.lgFilsDroit))

    blancSiFilsGaucheEstFeuille = distFilsGauche
    blancSiFilsDroitEstFeuille = distFilsDroit
    if(racine.filsGauche.cardinaux == 1):
        blancSiFilsGaucheEstFeuille = ""
    if(racine.filsDroit.cardinaux == 1):
        blancSiFilsDroitEstFeuille = ""

    return ("(" + formatNewick(racine.filsGauche, distFilsGauche) + blancSiFilsGaucheEstFeuille + "," + formatNewick(racine.filsDroit, distFilsDroit) + blancSiFilsDroitEstFeuille + ")")







    
