import math as ma
from parasite_rods_new import landing_rods

# Test class for the helicopter configuration
class HelicopterConfig:
    def __init__(self, config):
        self.config = config  # Initialize the configuration with the provided config dictionary

# Test data for the ParasiteArea configuration
test_config = {
    "ParasiteArea": {
        4: {  # Rear section (1 rod)
            "R": 0.14,        # Rod radius in meters
            "L": 1.07,        # Rod length in meters
            "alpha": ma.radians(29),  # Angle of attack in radians (converted from degrees)
            "S": 209.69        # Reference area in square meters
        },
        5: {  # Front section (2 rods)
            "R": 0.14,       # Rod radius in meters
            "L": 1.07,        # Rod length in meters
            "alpha": ma.radians(18),  # Angle of attack in radians (converted from degrees)
            "S": 209.69        # Reference area in square meters
        }
    }
}

# Execution of the test
if __name__ == "__main__":
    # Create an instance of HelicopterConfig with the test data
    helicopter = HelicopterConfig(test_config)

    # Call the landing_rods function to calculate the drag values
    f_A_lr, f_lr = landing_rods(helicopter)

    # Print the results of the calculations
    print(f"Combined parasite drag area (f_A_lr): {f_A_lr}")
    print(f"Total parasite drag (f_lr): {f_lr}")
