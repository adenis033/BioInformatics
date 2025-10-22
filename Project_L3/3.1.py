# The melting temp (tm) is the temp at which one half of a particular DNA will disociate
# and become a single stand of DNA primerland and seq are of critical importance
# in designing the parameters of a successful amplification, the tm of a nucleic acid duplex
# increseas both with its lenght and with increasing GC content

# tm = 4 ( G + C ) + 2 (A + T)
# tm = 81.5 + 16.6(log10([Nat])) + 0.41*(%GC) - 600/lenght

# Implement an app that computes the tm of a DNA seq by using one of those formulas or both of them

# Input = a string of DNA
# Output = temp in Celsius

import math

text = input("Enter DNA sequence: ").upper()

length = len(text)

A = text.count('A')
T = text.count('T')
G = text.count('G')
C = text.count('C')

tm_simple = 4 * (G + C) + 2 * (A + T)
gc = ((G + C) / length) * 100
Nat = 0.0001
tm_advanced = 81.5 + 16.6 * math.log10(Nat) + 0.41 * gc - (600 / length)

print(f"Tm (1st formula): {tm_simple:.2f} °C")
print(f"Tm (2nd formula): {tm_advanced:.2f} °C")
