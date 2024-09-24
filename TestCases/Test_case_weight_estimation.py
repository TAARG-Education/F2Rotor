import numpy as np
from ARR_Liberti_marco_ultimate import estimation_weight_prouty
from prettytable import PrettyTable

# Real AW109 helicopter data
helicopter_data = {
    'b': 0.25,  # Main rotor blade width (m) - Estimated based on the helicopter class
    'c': 0.1,  # Main rotor chord (m) - Estimated based on the helicopter class
    'R': 5.0,  # Main rotor radius (m) - Estimated based on the helicopter class
    'Omega': 5.0,  # Main rotor angular speed (rad/s) - Estimated based on the helicopter class
    'g': 9.81,  # Gravitational acceleration (m/s²) - Actual
    'A_H': 20.0,  # Horizontal rotor area (m²) - Estimated based on the helicopter class
    'AR_H': 5.0,  # Horizontal rotor aspect ratio - Estimated based on the helicopter class
    'A_V': 10.0,  # Vertical rotor area (m²) - Estimated based on the helicopter class
    'AR_V': 3.0,  # Vertical rotor aspect ratio - Estimated based on the helicopter class
    'n_tailgearboxes': 2,  # Number of tail gearboxes - Estimated based on the helicopter class
    'R_T': 1.5,  # Tail rotor radius (m) - Estimated based on the helicopter class
    'transm_hp_rating': 800,  # Transmission horsepower rating (hp) - Estimated based on the helicopter class
    'Omega_T': 5.0,  # Tail rotor angular speed (rad/s) - Estimated based on the helicopter class
    'GW': 6000,  # Gross weight of the helicopter (lbm) - Actual
    'L_F': 12.0,  # Fuselage length (m) - Estimated based on the helicopter class
    'S_wetF': 30.0,  # Fuselage wetted area (m²) - Estimated based on the helicopter class
    'n_wheellegs': 3,  # Number of landing gear legs - Actual
    'N_eng': 2,  # Number of engines - Actual
    'Installed_wt_pereng': 400,  # Installed weight per engine (lbm) - Actual
    'Cap_In_Gal': 200,  # Fuel capacity (gallons) - Actual
    'N_tanks': 2,  # Number of tanks - Actual
    'RPM_eng': 3000,  # Engine speed (RPM) - Estimated based on the helicopter class
    'tail_hp_rating': 200,  # Tail engine horsepower rating (hp) - Estimated based on the helicopter class
    'N_gearboxes': 2,  # Number of gearboxes - Estimated based on the helicopter class
}

# Expected estimated data for this class of helicopters
expected_empty_weight_kg = 1500  # Estimated empty weight based on helicopter class
expected_component_weights_kg = np.array([
    0.045,  # Main rotor - Estimated
    95.49,  # Tail (including tail rotor) - Estimated
    80.24,  # Fuselage - Estimated
    109.08,  # Landing gear - Estimated
    443.85,  # Propulsion - Estimated
    542.14,  # Transmission - Estimated
    10.68,  # Controls - Estimated
    16.31,  # Instruments - Estimated
    163.94,  # Hydraulic and electric systems - Estimated
    22.68,  # Avionics - Estimated
    27.95,  # Equipment - Estimated
    21.77   # Air conditioning and anti-ice system - Estimated
])

def test_weight_estimation():
    # Execute the weight estimation function
    empty_weight_kg, component_weights_kg = estimation_weight_prouty(helicopter_data)

    # Create a table to display estimated vs expected weights
    table = PrettyTable()
    table.field_names = ["Component", "Estimated Weight (kg)", "Expected Weight (kg)"]

    component_names = [
        "Main rotor", "Tail", "Fuselage", "Landing gear",
        "Propulsion", "Transmission", "Controls", "Instruments",
        "Hydraulic and electric systems", "Avionics", "Equipment",
        "Air conditioning and anti-ice system"
    ]

    for i, component_name in enumerate(component_names):
        table.add_row([component_name, round(component_weights_kg[i], 2), round(expected_component_weights_kg[i], 2)])

    print(table)

    # Validate total estimated weight
    assert abs(empty_weight_kg - expected_empty_weight_kg) <= 50, "Total weight estimate is outside the acceptable range."

    # Validate individual component weights
    np.testing.assert_array_almost_equal(component_weights_kg, expected_component_weights_kg, decimal=2)


# Run the test
if __name__ == '__main__':
    test_weight_estimation()
