import math

def wrap_phi(angle):
    """Wraps an angle to be between -180 and +180 degrees"""
    return (angle + 180) % 360 - 180

# The name of the file you want to generate
filename = "detector_geometry.dat"

# Open the file in 'write' mode
with open(filename, "w") as f:
    
    # --- 1. SX3 Geometry ---
    f.write("=========================================================\n")
    f.write("          SX3 BARREL AZIMUTHAL ANGLES (Degrees)          \n")
    f.write("=========================================================\n")
    f.write(" Det ID | Strip 0 | Strip 1 | Strip 2 | Strip 3 | Det Center\n")
    f.write("---------------------------------------------------------\n")
    
    for det_id in range(12): # Assuming 12 barrel detectors
        strips = []
        for stripF in range(4):
            stripF_rev = 3 - stripF
            num = (2 * stripF_rev - 3) * 40.30
            den = 8.0 * 88.0 * math.cos(math.radians(15.0))
            beta_n = 15.0 + math.degrees(math.atan2(num, den))
            phi_n = ((-det_id + 0.5) * 30.0 + beta_n) + 45.0
            strips.append(wrap_phi(phi_n))
        
        det_center = wrap_phi(((-det_id + 0.5) * 30.0 + 15.0) + 45.0)
        f.write(f"   {det_id:2d}   |  {strips[0]:6.2f} |  {strips[1]:6.2f} |  {strips[2]:6.2f} |  {strips[3]:6.2f} |   {det_center:6.2f}\n")


    # --- 2. PC Wire Geometry ---
    f.write("\n\n=========================================================\n")
    f.write("     PROPORTIONAL COUNTER WIRE ANGLES (Degrees)          \n")
    f.write("=========================================================\n")
    f.write(" Wire ID | Anode (-z) | Anode (+z) | Cathode (-z) | Cathode (+z)\n")
    f.write("----------------------------------------------------------------\n")

    k = 360.0 / 24.0
    offset_a1 = -6*k - 3*k
    offset_c1 = -4*k - 2*k - (360.0/48.0)
    wireShift = 3
    offset_a2 = offset_a1 + wireShift*k
    offset_c2 = offset_c1 - wireShift*k

    for i in range(24):
        phi_a1 = wrap_phi(-k * i + offset_a1)
        phi_a2 = wrap_phi(-k * i + offset_a2)
        phi_c1 = wrap_phi(k * i + offset_c1)
        phi_c2 = wrap_phi(k * i + offset_c2)
        
        f.write(f"    {i:2d}   |  {phi_a1:7.2f}   |  {phi_a2:7.2f}   |   {phi_c1:7.2f}    |   {phi_c2:7.2f}\n")

print(f"Success! The geometry lookup table has been saved to '{filename}'.")