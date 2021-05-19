import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Slider
import xlrd
from datetime import date


def SIR_1( s0, i0, g, a, n, t_max): # Evolution Sains, Infectés et restants ( restants correspond au mort démographique, aux morts liésà la maladie et aux guéris )
    s = [s0]
    i = [i0]
    r = [0]
    t = [0]
    d = t_max/n

    for j in range(1,n+1) :
        t.append(j*d)

    for j in range(n):
        s.append( s[j] + d* ( -g * s[j] * i[j] ))
        i.append( i[j] + d* ( g * s[j] * i[j] - a*i[j]))
        r.append( r[j] + d* a*i[j])

    t = np.array(t)
    s = np.array(s)
    i = np.array(i)
    r = np.array(r)
    
    return { "Titre" : "Model SIR",
             "Sains" : s,
             "Infectés" : i,
             "Restants" : r,
             "Temps" : t }

def SIR_2( s0, i0, beta, b, gamma, nb_point, t_max) :

    # S' = b*N - b*S - beta/N*S*I
    # I' = beta/N * S * I - gamma * I - b * I

    s = [s0]
    i = [i0]
    t = [0]
    d = t_max / nb_point

    N = s0 + i0 # Population initial

    for j in range(1, nb_point+1) :
        t.append(j*d)

    for j in range(nb_point):
        s.append( s[j] + d * ( b*N - b*s[j] - beta/N * s[j] * i[j]))
        i.append( i[j] + d * ( beta/N * s[j] * i[j] - gamma * i[j] - b * i[j]))

    t = np.array(t)
    s = np.array(s)
    i = np.array(i)

    r = np.array( len(i)*[N] ) - s - i  # S + I + R = N constant

    return { "Titre" : "Model SIR",
             "Sains" : s,
             "Infectés" : i,
             "Restants" : r,
             "Temps" : t }

def SIR_cumule( s0, i0, beta, b, gamma, nb_point, t_max) :

    # S' = b*N - b*S - beta/N*S*I
    # I' = beta/N * S * I - gamma * I - b * I

    s = [s0]
    i = [i0]
    t = [0]
    i_cumule = [i0]
    d = t_max / nb_point

    N = s0 + i0 # Population initial

    for j in range(1, nb_point+1) :
        t.append(j*d)

    for j in range(nb_point):
        s.append( s[j] + d * ( b*N - b*s[j] - beta/N * s[j] * i[j]))
        i.append( i[j] + d * ( beta/N * s[j] * i[j] - gamma * i[j] - b * i[j]))

        i_cumule.append( i_cumule[j] + max( 0, i[j+1] - i[j]) ) #i_cumulé contient la somme des nouveaux cas depuis t0 à t donnée.
        

    t = np.array(t)
    s = np.array(s)
    i = np.array(i)
    i_cumule = np.array(i_cumule)
    
    r = np.array( len(i)*[N] ) - s - i  # S + I + R = N constant

    return { "Titre" : "Model SIR ( infecté cumulé) ",
             "Inféctés Cumulé" : i_cumule,
             "Temps" : t }

def SIGM ( s0, i0, beta, b, eta, mu,  nb_point, t_max): # modéle SIR contabilisant le nombre de mort

    # S' = b*N - b*S - beta/N*S*I
    # I' = beta/N * S * I - (eta + mu) * I - b * I
    # M' = mu * I

    
    s = [s0]
    i = [i0]
    m = [0]
    t = [0]
    d = t_max / nb_point

    N = s0 + i0 # Population initial

    for j in range(1, nb_point+1) :
        t.append(j*d)

    for j in range(nb_point):
        s.append( s[j] + d * ( b*N - b*s[j] - beta/N * s[j] * i[j]))
        i.append( i[j] + d * ( beta/N * s[j] * i[j] - (eta + mu) * i[j] - b * i[j]))
        m.append( m[j] + d * ( mu * i[j] ))

    t = np.array(t)
    s = np.array(s)
    i = np.array(i)
    m = np.array(m)
    
    g = np.array( len(i)*[N] ) - s - i - m # S + I + R + M + G = N constant

    return { "Titre" : "Model SIGM",
             "Sains" : s,
             "Infectés" : i,
             "Guéris" : g,
             "Morts" : m,
             "Temps" : t }
    
    
def SIRP (s0, i0, beta, b, gamma,rho_s, rho_i, phi_s, phi_i, nb_point, t_max): # Modèle SIR avec des protections ( masques, gants ) ralentissant l'épidémie

    # S' = b*N - b*S - alpha*S*I
    # I' = alpha * S * I - gamma * I - b * I

    s = [s0]
    i = [i0]
    t = [0]
    d = t_max / nb_point

    N = s0 + i0 # Population initial
    

    for j in range(1, nb_point+1) :
        t.append(j*d)

    # Coefficient de contagion avec protections
    alpha = beta/N * (rho_s*rho_i*phi_s*phi_i + (1 - rho_s)*rho_i*phi_i + rho_s*(1 - rho_i)*phi_s + (1 - rho_s)*(1 - rho_i))

    
    for j in range(nb_point):
        s.append( s[j] + d * ( b*N - b*s[j] - alpha * s[j] * i[j]))
        i.append( i[j] + d * ( alpha * s[j] * i[j] - gamma * i[j] - b * i[j]))

    t = np.array(t)
    s = np.array(s)
    i = np.array(i)

    r = np.array( len(i)*[N] ) - s - i  # S + I + R = N constant

    return { "Titre" : "Model SIR avec protection",
             "Sains" : s,
             "Infectés" : i,
             "Restants" : r,
             "Temps" : t }


def SIGMP (s0, i0, beta, b, eta, mu,rho_s, rho_i, phi_s, phi_i, nb_point, t_max): # Modèle SIR avec des protections contabilisant les morts
    
    # S' = b*N - b*S - alpha*S*I
    # I' = alpha * S * I - (eta + mu) * I - b * I
    # M' = mu * I
    
    s = [s0]
    i = [i0]
    m = [0] 
    t = [0]
    d = t_max / nb_point

    N = s0 + i0 # Population initial
    

    for j in range(1, nb_point+1) :
        t.append(j*d)

    # Coefficient de contagion avec protections
    alpha = beta/N * (rho_s*rho_i*phi_s*phi_i + (1 - rho_s)*rho_i*phi_i + rho_s*(1 - rho_i)*phi_s + (1 - rho_s)*(1 - rho_i))

    
    for j in range(nb_point):
        s.append( s[j] + d * ( b*N - b*s[j] - alpha * s[j] * i[j]))
        i.append( i[j] + d * ( alpha * s[j] * i[j] - (eta + mu) * i[j] - b * i[j]))
        m.append( m[j] + d * ( mu * i[j] ))


    t = np.array(t)
    s = np.array(s)
    i = np.array(i)
    m = np.array(m)

    g = np.array( len(i)*[N] ) - s - i - m     # S + I + R + M + G = N constant

    return { "Titre" : "Model SIGM avec protection",
             "Sains" : s,
             "Infectés" : i,
             "Guéris" : g,
             "Morts" : m,
             "Temps" : t }



def show_model( m) : #Affiche un modéle représenter par un dictionnaire de np.array

    fig ,ax = plt.subplots()
    ax.clear()
    
    for titre, array in m.items() : 
        if( titre != "Temps" and titre != "Titre" ):
            ax.plot( m["Temps"], array, label=titre)     # On affiche les courbes

    ax.legend()
    ax.set_title(m["Titre"])

#Slider non aboutis

##def init_slider_SIR_2():
##
##    dic_slider = {
##        "s0" : Slider( plt.axes([0.25, 0.03, 0.50, 0.02]) , "s0 :", 0, 500, valinit=CI["s0"]),
##        "i0" : Slider( plt.axes([0.25, 0.05, 0.50, 0.02]) , "i0 :", 0, 500, valinit=CI["i0"]),
##        "beta" : Slider( plt.axes([0.25, 0.07, 0.50, 0.02]) , "beta :", 0, 10, valinit=CI["beta"]),
##        "b" : Slider( plt.axes([0.25, 0.09, 0.50, 0.02]) , "b :", 0, 1, valinit=CI["b"]),
##        "gamma" : Slider( plt.axes([0.25, 0.11, 0.50, 0.02]) , "gamma : :", 0, 1, valinit=CI["gamma"]),
##        "nb_points" : Slider( plt.axes([0.25, 0.13, 0.50, 0.02]) , "nb_points :", 0, 500, valinit=CI["nb_points"]),
##        "t_max" : Slider( plt.axes([0.25, 0.15, 0.50, 0.02]) , "t_max :", 0, 40, valinit=CI["t_max"]) }
##
##    def uptade( val) :
##        for cle, slider in dic_slider :
##            print(slider)
##            #CI[nom] = slider.val
##
##        showModel(SIR_2 ( CI["s0"], CI["i0"], CI["beta"], CI["b"], CI["gamma"], CI["nb_points"], CI["t_max"]))
##    
##    for nom, slider in dic_slider.items():
##        slider.on_changed(uptade)
##
##    
##    plt.show()


def lire_donnees(pays): # Lis et représente les données d'une épidémie
    global date
    
    document = xlrd.open_workbook("ebola_chiffres.xls")
    feuille = document.sheet_by_index(0)

    n = feuille.nrows # nombre lignes de la feuille
    cas_cumules = [] # le document informe sur les nombre de cas cumulé à une date donnée

    for j in range(n):
        if ( feuille.cell_value(rowx=j, colx=0) == "Cumulative number of confirmed Ebola cases" and
             feuille.cell_value(rowx=j, colx=1) == pays ):
            
            cas_cumules.append(( feuille.cell_value(rowx=j, colx=2), feuille.cell_value(rowx=j, colx=3)))

    # on transforme les dates en nombre de jour depuis 01-01-2014
    date_init = date(2014,1,1)
    
    for j in range(len(cas_cumules)): 
        
        date_suivante = cas_cumules[j][0]
    
        annee = int(date_suivante[0:4])
        mois = 0
        jour = 0
        
        if( date_suivante[5] != '0' ):
            mois = int(date_suivante[5:7])
        else:
            mois = int(date_suivante[6])

        if( date_suivante[8] != '0' ):
            jour = int(date_suivante[8:10])
        else:
            jour = int(date_suivante[9])

        diff_jour = ( date(annee,mois,jour) - date_init ).days
        
        cas_cumules[j] = (diff_jour, cas_cumules[j][1])

    cas_cumules.sort(key= lambda x:x[0]) #on range ces nombres de jours dans l'ordre croissant

    for j in range(len(cas_cumules) - 1 , -1, -1) : # on fait demarrer les jours du J0 de l'épidémie
        cas_cumules[j] = ( cas_cumules[j][0] - cas_cumules[0][0], cas_cumules[j][1] )

    # on genere le nb de nouveau cas
    nouveau_cas = [(0,0)]

    for j in range(1,len(cas_cumules)):
        nouveau_cas.append((cas_cumules[j][0],
                            abs((cas_cumules[j][1] - cas_cumules[j-1][1] ))))

    # on genere le nb de nouveau cas par jour 
    nouveau_cas_par_jour = [(0,0)]

    for j in range(1, len(cas_cumules)):
        nouveau_cas_par_jour.append( ( nouveau_cas[j][0], nouveau_cas[j][1] / (nouveau_cas[j][0] - nouveau_cas[j-1][0])))


    #Affichage cas cumulé
    dates = [cas_cumules[j][0] for j in range(len(cas_cumules))]
    cas = [cas_cumules[j][1] for j in range(len(cas_cumules))]
    #ax.plot(date, cas)
    
    
    #Affichage nouveau cas par jour brute
    fig ,ax = plt.subplots()
    dates = [nouveau_cas_par_jour[j][0] for j in range(len(nouveau_cas_par_jour))]
    cas = [nouveau_cas_par_jour[j][1] for j in range(len(nouveau_cas_par_jour))]
    ax.plot(dates, cas, label = "Nouveau cas par jour sans lissage")

    #Affichage courbe lissé avec gaussienne
    from scipy.ndimage.filters import gaussian_filter1d # Lissage des courbes
    cas = gaussian_filter1d(cas, sigma=5)
    ax.plot(dates, cas, label = "Nouveau cas avec lissage gaussien")

    #Affichage lissage moyenne sur 7 jours
    nouveau_cas_lisse = lissage_sur_periode( 7, nouveau_cas)
    dates = [nouveau_cas_lisse[j][0] for j in range(len(nouveau_cas_lisse))]
    cas = [nouveau_cas_lisse[j][1] for j in range(len(nouveau_cas_lisse))]
    ax.plot(dates, cas, label = "Nouveau cas avec lissage sur 7 jours" , color = "red")
    
    ax.legend()
    ax.set_title("Cas par jour d'ébola en " + pays)
    plt.show()
    
    #return cas_cumules


def lissage_sur_periode ( periode , tab) : #on lisse en faisant la moyenne sur les *periode* derniers jours 
    
    decoupe = []
    for i in range((tab[-1][0]//periode + 1)): # initialisation
        decoupe.append([])

    for date, val in tab : # on remplit
        decoupe[date//periode].append( (date, val) )
    

    for i in range(1,len(decoupe)): # on met les cases vides égal à leurs cases précédente
        if( decoupe[i] == []):
            decoupe[i] = decoupe[i-1]

    moyennes = []
    for i in range(len(decoupe)) : # on fait la moyenne
        
        moyenne = 0
        for date, val in decoupe[i] :
            moyenne += val
        moyenne = moyenne / periode
        moyennes.append((periode*(i+1),moyenne))
    
    return moyennes

#Condition Initial
CI = {"s0" : 90,        #Nombre de personnes saines initials
      "i0" : 10,        #Nombr d'infectés initiaux
      "beta" : 2,       #Nombre de personnes recontrés pr unité de temps par un individu infecté
      "b" : 0,          #Taux de mortalité démographique ( on considère souvent une petite échelle de temps, ce qui revient à b = 0 )
      "gamma" : 0.5,    #taux de personnes échapant de la maladie par unité de temps ( morts ou guéris ) 

      "eta" : 0.1,      #Taux de personnes guerissants de la maladie
      "mu" : 0.4,       #Taux de mortalité
      
      "rho_s" : 1/3,    #Taux de population saine portant une protection ( masque, gants ... )
      "rho_i" : 2/3,    #Taux de la population infecté portant une protection
      "phi_i" : 0.2,    #Taux d'amortissement de la transmission du au port de protection parl les peronnes infectiés
      "phi_s" : 0.8,    #Taux d'amortissement de la transmission du au port de protection parl les peronnes sains

      "nb_points" : 100,#Nombre de point de calculs
      "t_max" : 10      #Unité de temps jusque laquelle on affiche la courbre
      }   

# MODEL INITIAL 
#m_SIR = SIR_1( 0.9, 0.1, 0.7, 0.2, 50, 25)

############################

#MODEL SANS MORT

m_SIR_2 = SIR_2 ( CI["s0"], CI["i0"], CI["beta"], CI["b"], CI["gamma"], CI["nb_points"], CI["t_max"])
show_model(m_SIR_2)


#m_SIR_cumule = SIR_cumule ( CI["s0"], CI["i0"], CI["beta"], CI["b"], CI["gamma"], CI["nb_points"], CI["t_max"])
#show_model(m_SIR_cumule)


m_SIRP = SIRP( CI["s0"], CI["i0"], CI["beta"], CI["b"], CI["gamma"], CI["rho_s"], CI["rho_i"], CI["phi_s"], CI["phi_i"], CI["nb_points"], CI["t_max"])
show_model(m_SIRP)

############################

# MODELS AVEC MORT 

m_SIGM = SIGM( CI["s0"], CI["i0"], CI["beta"], CI["b"], CI["eta"], CI["mu"], CI["nb_points"], CI["t_max"])
show_model(m_SIGM)

m_SIGMP = SIGMP( CI["s0"], CI["i0"], CI["beta"], CI["b"], CI["eta"],CI["mu"], CI["rho_s"], CI["rho_i"], CI["phi_s"], CI["phi_i"], CI["nb_points"], CI["t_max"])
show_model(m_SIGMP)

lire_donnees("Sierra Leone")
plt.show()
