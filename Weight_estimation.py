#Author:Liberti Marco
#Rotary Wing Aerodynamics course, prof. Renato Tognaccini
#University of Naples Federico II
#Academic Year 2023-2024
#Date:17/09/2024
#
#Version:1.0.0
# import numpy as np: Importa NumPy per operazioni matematiche avanzate e gestione degli array.

# Funzione di Conversione delle Masse (convmass):
# Converte le masse tra libbre (lbm) e chilogrammi (kg). La funzione usa un dizionario di fattori di conversione e restituisce le masse convertite.

# Funzione Principale (estimation_weight_prouty):
# Parametri: Prende un dizionario elicottero contenente i parametri necessari per il calcolo dei pesi.

# Calcolo dei Pesi:
# Ogni variabile rappresenta un peso specifico di una parte dell'elicottero, calcolato con una formula empirica basata sui parametri forniti.
# Le formule sono usate per stimare i pesi del rotore principale, del corpo, della coda, del carrello di atterraggio, del sistema di propulsione, ecc.

# Conversione delle Masse:
# Matrix_masses è un array contenente i pesi di tutte le componenti. È convertito da libbre a chilogrammi.
# M_empty_kg rappresenta il peso totale stimato dell'elicottero in chilogrammi.

# Output:
# Restituisce M_empty_kg (peso totale) e Matrix_masses_kg (array delle masse delle singole componenti).




import numpy as np

def convmass(masses, from_unit, to_unit):
    """
    Converte le masse tra unità di misura (lbm a kg e viceversa).
    Parametri:
        masses (array o float): Masse da convertire.
        from_unit (str): Unità di misura di origine (es. 'lbm').
        to_unit (str): Unità di misura di destinazione (es. 'kg').
    Ritorna:
        Array o float delle masse convertite.
    """
    # Fattori di conversione tra lbm e kg
    conversion_factors = {
        ('lbm', 'kg'): 0.453592,
        ('kg', 'lbm'): 2.20462
    }
    
    if from_unit == to_unit:
        return masses
    return masses * conversion_factors[(from_unit, to_unit)]

def estimation_weight_prouty(elicottero):
    """
    Stima i pesi di un elicottero basandosi su parametri forniti.
    Parametri:
        elicottero (dict): Dizionario contenente i parametri necessari per il calcolo.
    Ritorna:
        M_empty_kg (float): Peso totale stimato dell'elicottero in chilogrammi.
        Matrix_masses_kg (array): Array dei pesi delle componenti in chilogrammi.
    """

    # Parametri dell'elicottero (estratti dal dizionario 'elicottero')
    b = elicottero['b']  # Larghezza della pala del rotore principale (m)
    c = elicottero['c']  # Cord (larghezza) della pala del rotore principale (m)
    R = elicottero['R']  # Raggio del rotore principale (m)
    Omega = elicottero['Omega']  # Velocità angolare del rotore principale (rad/s)
    g = elicottero['g']  # Accelerazione gravitazionale (m/s²)
    A_H = elicottero['A_H']  # Area del rotore orizzontale (m²)
    AR_H = elicottero['AR_H']  # Rapporto d'aspetto del rotore orizzontale
    A_V = elicottero['A_V']  # Area del rotore verticale (m²)
    AR_V = elicottero['AR_V']  # Rapporto d'aspetto del rotore verticale
    n_tailgearboxes = elicottero['n_tailgearboxes']  # Numero di ingranaggi nella coda
    R_T = elicottero['R_T']  # Raggio del rotore della coda (m)
    transm_hp_rating = elicottero['transm_hp_rating']  # Potenza del sistema di trasmissione (hp)
    Omega_T = elicottero['Omega_T']  # Velocità angolare del rotore della coda (rad/s)
    GW = elicottero['GW']  # Peso lordo dell'elicottero (lbm)
    L_F = elicottero['L_F']  # Lunghezza del fusolotto (m)
    S_wetF = elicottero['S_wetF']  # Superficie bagnata del fusolotto (m²)
    n_wheellegs = elicottero['n_wheellegs']  # Numero di gambe del carrello (ruote)
    N_eng = elicottero['N_eng']  # Numero di motori
    Installed_wt_pereng = elicottero['Installed_wt_pereng']  # Peso installato per motore (lbm)
    Cap_In_Gal = elicottero['Cap_In_Gal']  # Capacità del serbatoio (galloni)
    N_tanks = elicottero['N_tanks']  # Numero di serbatoi
    RPM_eng = elicottero['RPM_eng']  # Velocità del motore (giri/minuto)
    tail_hp_rating = elicottero['tail_hp_rating']  # Potenza del motore di coda (hp)
    N_gearboxes = elicottero['N_gearboxes']  # Numero di ingranaggi
    # Note: 'b', 'c', 'Omega' e 'R' sono riutilizzati e non richiedono un'ulteriore definizione

    # Calcolo del peso delle pale del rotore principale
    W_bM = 0.026 * b**0.66 * c * R**1.3 * (Omega * R)**0.67
    
    # Calcolo del momento di inerzia delle pale del rotore principale
    J = 0.0311 * W_bM * R**2
    
    # Calcolo del peso del hub e degli snodi del rotore principale
    W_hM = 0.0037 * b**0.28 * R**1.5 * (Omega * R)**0.43 * (0.67 * W_bM + g * J / R**2)**0.55
    
    # Peso totale del rotore principale
    W_rotor = W_bM + W_hM
    
    # Calcolo del peso del rotore orizzontale
    W_H = 0.72 * A_H**1.2 * AR_H**0.32
    
    # Calcolo del peso del rotore verticale
    W_V = 1.05 * A_V**0.94 * AR_V**0.53 * n_tailgearboxes**0.71
    
    # Calcolo del peso del rotore di coda
    W_T = 1.4 * R_T**0.09 * (transm_hp_rating / Omega)**0.9
    
    # Peso totale della coda
    W_tail = W_H + W_V + W_T
    
    # Calcolo del peso del corpo dell'elicottero (fuselotto)
    W_body = 6.9 * (GW / 1000)**0.49 * L_F**0.61 * S_wetF**0.25
    
    # Calcolo del peso del carrello di atterraggio
    W_landinggear = 40 * (GW / 1000)**0.67 * n_wheellegs**0.54
    
    # Calcolo del peso dell'installazione del motore
    W_eng = N_eng * Installed_wt_pereng
    
    # Calcolo del peso dei sottosistemi di propulsione
    W_Pss = 2 * W_eng**0.59 * N_eng**0.79
    
    # Peso totale della propulsione
    W_prop = W_eng + W_Pss
    
    # Calcolo del peso del sistema di alimentazione
    W_FS = 0.43 * Cap_In_Gal**0.77 * N_tanks**0.59
    
    # Calcolo del peso del sistema di trasmissione
    W_DS = 13.6 * transm_hp_rating**0.82 * (RPM_eng / 1000)**0.037 * \
           ((tail_hp_rating / transm_hp_rating) * (Omega / Omega_T))**0.068 * \
           N_gearboxes**0.066 / Omega**0.64
    
    # Peso totale della trasmissione
    W_trasm = W_FS + W_DS
    
    # Calcolo del peso dei controlli della cabina
    W_CC = 11.5 * (GW / 1000)**0.40
    
    # Calcolo del peso del controllo dei sistemi
    W_SC = 36 * b * c**2.2 * (Omega * R / 1000)**3.2
    
    # Peso totale dei controlli
    W_controls = W_CC + W_SC
    
    # Peso dell'unità di potenza ausiliaria
    W_APP = 150
    W_aux = W_APP
    
    # Calcolo del peso degli strumenti
    W_I = 3.5 * (GW / 1000)**1.3
    W_instr = W_I
    
    # Calcolo del peso del sistema idraulico
    W_hyd = 37 * b**0.63 * c**1.3 * (Omega * R / 1000)**2.1
    
    # Calcolo del peso del sistema elettrico (escluso il peso idraulico)
    W_EL = 9.6 * transm_hp_rating**0.65 / (GW / 1000)**0.4 - W_hyd
    W_hpe = W_hyd + W_EL
    
    # Peso dell'avionica
    W_av = 50
    
    # Calcolo del peso degli arredi e dell'equipaggiamento
    W_FE = 6 * (GW / 1000)**1.3
    W_equip = W_FE
    
    # Peso dell'aria condizionata e del sistema anti-ghiaccio
    W_AC_AI = 8 * (GW / 1000)
    W_antiicing = W_AC_AI
    
    # Creazione di una matrice con i pesi di ciascuna componente
    Matrix_masses = np.array([W_rotor, W_tail, W_body, W_landinggear, W_prop,
                              W_trasm, W_controls, W_instr, W_hpe, W_av,
                              W_equip, W_antiicing])
    
    # Conversione delle masse dalla libbra (lbm) al chilogrammo (kg)
    Matrix_masses_kg = convmass(Matrix_masses, 'lbm', 'kg')
    
    # Calcolo del peso totale dell'elicottero e conversione in chilogrammi
    M_empty_kg = convmass(np.sum(Matrix_masses), 'lbm', 'kg')
    
    return M_empty_kg, Matrix_masses_kg
