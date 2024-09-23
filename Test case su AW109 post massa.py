import unittest
import numpy as np
from ARR_Liberti_marco_ultimate import convmass, estimation_weight_prouty
from prettytable import PrettyTable


class TestHelicopterWeightEstimation(unittest.TestCase):

    def setUp(self):
        # Parametri reali del velivolo AW109
        self.elicottero = {
            'b': 0.25,  # Larghezza della pala del rotore principale (m)
            'c': 0.1,  # Corda del rotore principale (m)
            'R': 5.0,  # Raggio del rotore principale (m)
            'Omega': 5.0,  # Velocità angolare del rotore principale (rad/s)
            'g': 9.81,  # Accelerazione gravitazionale (m/s²)
            'A_H': 20.0,  # Area del rotore orizzontale (m²)
            'AR_H': 5.0,  # Rapporto d'aspetto del rotore orizzontale
            'A_V': 10.0,  # Area del rotore verticale (m²)
            'AR_V': 3.0,  # Rapporto d'aspetto del rotore verticale
            'n_tailgearboxes': 2,  # Numero di ingranaggi nella coda
            'R_T': 1.5,  # Raggio del rotore di coda (m)
            'transm_hp_rating': 800,  # Potenza del sistema di trasmissione (hp)
            'Omega_T': 5.0,  # Velocità angolare del rotore di coda (rad/s)
            'GW': 6000,  # Peso lordo dell'elicottero (lbm)
            'L_F': 12.0,  # Lunghezza del fusolotto (m)
            'S_wetF': 30.0,  # Superficie bagnata della fusoliera (m²)
            'n_wheellegs': 3,  # Numero di gambe del carrello
            'N_eng': 2,  # Numero di motori
            'Installed_wt_pereng': 400,  # Peso installato per motore (lbm)
            'Cap_In_Gal': 200,  # Capacità del serbatoio (galloni)
            'N_tanks': 2,  # Numero di serbatoi
            'RPM_eng': 3000,  # Velocità del motore (giri/minuto)
            'tail_hp_rating': 200,  # Potenza del motore di coda (hp)
            'N_gearboxes': 2,  # Numero di ingranaggi
        }

        # Parametri stimati tipici per questa classe di elicotteri (valori di riferimento)
        self.expected_M_empty_kg = 1500  # Peso a vuoto stimato in base alla classe di elicotteri
        self.expected_matrix_masses_kg = np.array([
            0.045,  # Rotore principale
            95.49,  # Coda (incluso il rotore di coda)
            80.24,  # Fusoliera
            109.08,  # Carrello di atterraggio
            443.85,  # Propulsione
            542.14,  # Trasmissione
            10.68,  # Controlli
            16.31,  # Strumenti
            163.94,  # Sistemi idraulici ed elettrici
            22.68,  # Avionica
            27.95,  # Equipaggiamento
            21.77  # Aria condizionata e sistema anti-ghiaccio
        ])

    def test_estimation_weight_prouty(self):
        # Eseguiamo la funzione
        M_empty_kg, Matrix_masses_kg = estimation_weight_prouty(self.elicottero)

        # Creiamo la tabella con i pesi stimati e i pesi attesi
        table = PrettyTable()
        table.field_names = ["Componente", "Peso Stimato (kg)", "Peso Atteso (kg)"]

        component_names = [
            "Rotore principale", "Coda", "Fusoliera", "Carrello di atterraggio",
            "Propulsione", "Trasmissione", "Controlli", "Strumenti",
            "Sistemi idraulici ed elettrici", "Avionica", "Equipaggiamento",
            "Aria condizionata e sistema anti-ghiaccio"
        ]

        for i, component_name in enumerate(component_names):
            table.add_row([component_name, round(Matrix_masses_kg[i], 2), round(self.expected_matrix_masses_kg[i], 2)])

        print(table)

        # Verifica del peso totale stimato
        self.assertAlmostEqual(M_empty_kg, self.expected_M_empty_kg, delta=50)

        # Verifica dei singoli componenti
        np.testing.assert_array_almost_equal(Matrix_masses_kg, self.expected_matrix_masses_kg, decimal=2)


if __name__ == '__main__':
    unittest.main()
