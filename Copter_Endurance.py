import numpy as np

def Endurance(Pn, Wfuel, Voo):
    """
    Funzione che calcola il tempo di volo (endurance) conoscendo la quantità di carburante
    sull'elicottero in volo traslato a diverse velocità.
    
    Input:
    Pn    = Potenza in KWatt
    Wfuel = Peso del carburante in Kg
    Voo   = Lista o array delle velocità di volo traslato

    Restituisce:
    Un array con i tempi di volo (endurance) per ogni velocità fornita.
    """
    P_s = Pn            # Potenza in KW
    Wfuel_kg = Wfuel    # Peso del carburante in kg
    SFC = 0.458         # Consumo specifico di carburante (assunto costante)

    tempo = np.zeros(len(Voo))  # Inizializzazione dell'array per i tempi di volo

    for i, v in enumerate(Voo):
        t_i = Wfuel_kg / (SFC * P_s/60)  # Tempo di endurance in min
        tempo = t_i
        i +=1
    return tempo
