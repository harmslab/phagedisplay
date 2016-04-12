__description__ = \
"""
Pre-defined weight-matrices for sequence distance calculations.
"""
__author__ = "Michael J. Harms"
__date__ = "2016-04-11"

blosum62 = { \
    '*' : {'*' :  1.0, 'A' : -4.0, 'B' : -4.0, 'C' : -4.0, 'D' : -4.0, 
           'E' : -4.0, 'F' : -4.0, 'G' : -4.0, 'H' : -4.0, 'I' : -4.0, 
           'K' : -4.0, 'L' : -4.0, 'M' : -4.0, 'N' : -4.0, 'P' : -4.0, 
           'Q' : -4.0, 'R' : -4.0, 'S' : -4.0, 'T' : -4.0, 'V' : -4.0, 
           'W' : -4.0, 'X' : -4.0, 'Y' : -4.0, 'Z' : -4.0, },

    'A' : {'*' : -4.0, 'A' :  4.0, 'B' : -2.0, 'C' :  0.0, 'D' : -2.0, 
           'E' : -1.0, 'F' : -2.0, 'G' :  0.0, 'H' : -2.0, 'I' : -1.0, 
           'K' : -1.0, 'L' : -1.0, 'M' : -1.0, 'N' : -2.0, 'P' : -1.0, 
           'Q' : -1.0, 'R' : -1.0, 'S' :  1.0, 'T' :  0.0, 'V' :  0.0, 
           'W' : -3.0, 'X' :  0.0, 'Y' : -2.0, 'Z' : -1.0, },

    'B' : {'*' : -4.0, 'A' : -2.0, 'B' :  4.0, 'C' : -3.0, 'D' :  4.0, 
           'E' :  1.0, 'F' : -3.0, 'G' : -1.0, 'H' :  0.0, 'I' : -3.0, 
           'K' :  0.0, 'L' : -4.0, 'M' : -3.0, 'N' :  3.0, 'P' : -2.0, 
           'Q' :  0.0, 'R' : -1.0, 'S' :  0.0, 'T' : -1.0, 'V' : -3.0, 
           'W' : -4.0, 'X' : -1.0, 'Y' : -3.0, 'Z' :  1.0, },

    'C' : {'*' : -4.0, 'A' :  0.0, 'B' : -3.0, 'C' :  9.0, 'D' : -3.0, 
           'E' : -4.0, 'F' : -2.0, 'G' : -3.0, 'H' : -3.0, 'I' : -1.0, 
           'K' : -3.0, 'L' : -1.0, 'M' : -1.0, 'N' : -3.0, 'P' : -3.0, 
           'Q' : -3.0, 'R' : -3.0, 'S' : -1.0, 'T' : -1.0, 'V' : -1.0, 
           'W' : -2.0, 'X' : -2.0, 'Y' : -2.0, 'Z' : -3.0, },

    'D' : {'*' : -4.0, 'A' : -2.0, 'B' :  4.0, 'C' : -3.0, 'D' :  6.0, 
           'E' :  2.0, 'F' : -3.0, 'G' : -1.0, 'H' : -1.0, 'I' : -3.0, 
           'K' : -1.0, 'L' : -4.0, 'M' : -3.0, 'N' :  1.0, 'P' : -1.0, 
           'Q' :  0.0, 'R' : -2.0, 'S' :  0.0, 'T' : -1.0, 'V' : -3.0, 
           'W' : -4.0, 'X' : -1.0, 'Y' : -3.0, 'Z' :  1.0, },

    'E' : {'*' : -4.0, 'A' : -1.0, 'B' :  1.0, 'C' : -4.0, 'D' :  2.0, 
           'E' :  5.0, 'F' : -3.0, 'G' : -2.0, 'H' :  0.0, 'I' : -3.0, 
           'K' :  1.0, 'L' : -3.0, 'M' : -2.0, 'N' :  0.0, 'P' : -1.0, 
           'Q' :  2.0, 'R' :  0.0, 'S' :  0.0, 'T' : -1.0, 'V' : -2.0, 
           'W' : -3.0, 'X' : -1.0, 'Y' : -2.0, 'Z' :  4.0, },

    'F' : {'*' : -4.0, 'A' : -2.0, 'B' : -3.0, 'C' : -2.0, 'D' : -3.0, 
           'E' : -3.0, 'F' :  6.0, 'G' : -3.0, 'H' : -1.0, 'I' :  0.0, 
           'K' : -3.0, 'L' :  0.0, 'M' :  0.0, 'N' : -3.0, 'P' : -4.0, 
           'Q' : -3.0, 'R' : -3.0, 'S' : -2.0, 'T' : -2.0, 'V' : -1.0, 
           'W' :  1.0, 'X' : -1.0, 'Y' :  3.0, 'Z' : -3.0, },

    'G' : {'*' : -4.0, 'A' :  0.0, 'B' : -1.0, 'C' : -3.0, 'D' : -1.0, 
           'E' : -2.0, 'F' : -3.0, 'G' :  6.0, 'H' : -2.0, 'I' : -4.0, 
           'K' : -2.0, 'L' : -4.0, 'M' : -3.0, 'N' :  0.0, 'P' : -2.0, 
           'Q' : -2.0, 'R' : -2.0, 'S' :  0.0, 'T' : -2.0, 'V' : -3.0, 
           'W' : -2.0, 'X' : -1.0, 'Y' : -3.0, 'Z' : -2.0, },

    'H' : {'*' : -4.0, 'A' : -2.0, 'B' :  0.0, 'C' : -3.0, 'D' : -1.0, 
           'E' :  0.0, 'F' : -1.0, 'G' : -2.0, 'H' :  8.0, 'I' : -3.0, 
           'K' : -1.0, 'L' : -3.0, 'M' : -2.0, 'N' :  1.0, 'P' : -2.0, 
           'Q' :  0.0, 'R' :  0.0, 'S' : -1.0, 'T' : -2.0, 'V' : -3.0, 
           'W' : -2.0, 'X' : -1.0, 'Y' :  2.0, 'Z' :  0.0, },

    'I' : {'*' : -4.0, 'A' : -1.0, 'B' : -3.0, 'C' : -1.0, 'D' : -3.0, 
           'E' : -3.0, 'F' :  0.0, 'G' : -4.0, 'H' : -3.0, 'I' :  4.0, 
           'K' : -3.0, 'L' :  2.0, 'M' :  1.0, 'N' : -3.0, 'P' : -3.0, 
           'Q' : -3.0, 'R' : -3.0, 'S' : -2.0, 'T' : -1.0, 'V' :  3.0, 
           'W' : -3.0, 'X' : -1.0, 'Y' : -1.0, 'Z' : -3.0, },

    'K' : {'*' : -4.0, 'A' : -1.0, 'B' :  0.0, 'C' : -3.0, 'D' : -1.0, 
           'E' :  1.0, 'F' : -3.0, 'G' : -2.0, 'H' : -1.0, 'I' : -3.0, 
           'K' :  5.0, 'L' : -2.0, 'M' : -1.0, 'N' :  0.0, 'P' : -1.0, 
           'Q' :  1.0, 'R' :  2.0, 'S' :  0.0, 'T' : -1.0, 'V' : -2.0, 
           'W' : -3.0, 'X' : -1.0, 'Y' : -2.0, 'Z' :  1.0, },

    'L' : {'*' : -4.0, 'A' : -1.0, 'B' : -4.0, 'C' : -1.0, 'D' : -4.0, 
           'E' : -3.0, 'F' :  0.0, 'G' : -4.0, 'H' : -3.0, 'I' :  2.0, 
           'K' : -2.0, 'L' :  4.0, 'M' :  2.0, 'N' : -3.0, 'P' : -3.0, 
           'Q' : -2.0, 'R' : -2.0, 'S' : -2.0, 'T' : -1.0, 'V' :  1.0, 
           'W' : -2.0, 'X' : -1.0, 'Y' : -1.0, 'Z' : -3.0, },

    'M' : {'*' : -4.0, 'A' : -1.0, 'B' : -3.0, 'C' : -1.0, 'D' : -3.0, 
           'E' : -2.0, 'F' :  0.0, 'G' : -3.0, 'H' : -2.0, 'I' :  1.0, 
           'K' : -1.0, 'L' :  2.0, 'M' :  5.0, 'N' : -2.0, 'P' : -2.0, 
           'Q' :  0.0, 'R' : -1.0, 'S' : -1.0, 'T' : -1.0, 'V' :  1.0, 
           'W' : -1.0, 'X' : -1.0, 'Y' : -1.0, 'Z' : -1.0, },

    'N' : {'*' : -4.0, 'A' : -2.0, 'B' :  3.0, 'C' : -3.0, 'D' :  1.0, 
           'E' :  0.0, 'F' : -3.0, 'G' :  0.0, 'H' :  1.0, 'I' : -3.0, 
           'K' :  0.0, 'L' : -3.0, 'M' : -2.0, 'N' :  6.0, 'P' : -2.0, 
           'Q' :  0.0, 'R' :  0.0, 'S' :  1.0, 'T' :  0.0, 'V' : -3.0, 
           'W' : -4.0, 'X' : -1.0, 'Y' : -2.0, 'Z' :  0.0, },

    'P' : {'*' : -4.0, 'A' : -1.0, 'B' : -2.0, 'C' : -3.0, 'D' : -1.0, 
           'E' : -1.0, 'F' : -4.0, 'G' : -2.0, 'H' : -2.0, 'I' : -3.0, 
           'K' : -1.0, 'L' : -3.0, 'M' : -2.0, 'N' : -2.0, 'P' :  7.0, 
           'Q' : -1.0, 'R' : -2.0, 'S' : -1.0, 'T' : -1.0, 'V' : -2.0, 
           'W' : -4.0, 'X' : -2.0, 'Y' : -3.0, 'Z' : -1.0, },

    'Q' : {'*' : -4.0, 'A' : -1.0, 'B' :  0.0, 'C' : -3.0, 'D' :  0.0, 
           'E' :  2.0, 'F' : -3.0, 'G' : -2.0, 'H' :  0.0, 'I' : -3.0, 
           'K' :  1.0, 'L' : -2.0, 'M' :  0.0, 'N' :  0.0, 'P' : -1.0, 
           'Q' :  5.0, 'R' :  1.0, 'S' :  0.0, 'T' : -1.0, 'V' : -2.0, 
           'W' : -2.0, 'X' : -1.0, 'Y' : -1.0, 'Z' :  3.0, },

    'R' : {'*' : -4.0, 'A' : -1.0, 'B' : -1.0, 'C' : -3.0, 'D' : -2.0, 
           'E' :  0.0, 'F' : -3.0, 'G' : -2.0, 'H' :  0.0, 'I' : -3.0, 
           'K' :  2.0, 'L' : -2.0, 'M' : -1.0, 'N' :  0.0, 'P' : -2.0, 
           'Q' :  1.0, 'R' :  5.0, 'S' : -1.0, 'T' : -1.0, 'V' : -3.0, 
           'W' : -3.0, 'X' : -1.0, 'Y' : -2.0, 'Z' :  0.0, },

    'S' : {'*' : -4.0, 'A' :  1.0, 'B' :  0.0, 'C' : -1.0, 'D' :  0.0, 
           'E' :  0.0, 'F' : -2.0, 'G' :  0.0, 'H' : -1.0, 'I' : -2.0, 
           'K' :  0.0, 'L' : -2.0, 'M' : -1.0, 'N' :  1.0, 'P' : -1.0, 
           'Q' :  0.0, 'R' : -1.0, 'S' :  4.0, 'T' :  1.0, 'V' : -2.0, 
           'W' : -3.0, 'X' :  0.0, 'Y' : -2.0, 'Z' :  0.0, },

    'T' : {'*' : -4.0, 'A' :  0.0, 'B' : -1.0, 'C' : -1.0, 'D' : -1.0, 
           'E' : -1.0, 'F' : -2.0, 'G' : -2.0, 'H' : -2.0, 'I' : -1.0, 
           'K' : -1.0, 'L' : -1.0, 'M' : -1.0, 'N' :  0.0, 'P' : -1.0, 
           'Q' : -1.0, 'R' : -1.0, 'S' :  1.0, 'T' :  5.0, 'V' :  0.0, 
           'W' : -2.0, 'X' :  0.0, 'Y' : -2.0, 'Z' : -1.0, },

    'V' : {'*' : -4.0, 'A' :  0.0, 'B' : -3.0, 'C' : -1.0, 'D' : -3.0, 
           'E' : -2.0, 'F' : -1.0, 'G' : -3.0, 'H' : -3.0, 'I' :  3.0, 
           'K' : -2.0, 'L' :  1.0, 'M' :  1.0, 'N' : -3.0, 'P' : -2.0, 
           'Q' : -2.0, 'R' : -3.0, 'S' : -2.0, 'T' :  0.0, 'V' :  4.0, 
           'W' : -3.0, 'X' : -1.0, 'Y' : -1.0, 'Z' : -2.0, },

    'W' : {'*' : -4.0, 'A' : -3.0, 'B' : -4.0, 'C' : -2.0, 'D' : -4.0, 
           'E' : -3.0, 'F' :  1.0, 'G' : -2.0, 'H' : -2.0, 'I' : -3.0, 
           'K' : -3.0, 'L' : -2.0, 'M' : -1.0, 'N' : -4.0, 'P' : -4.0, 
           'Q' : -2.0, 'R' : -3.0, 'S' : -3.0, 'T' : -2.0, 'V' : -3.0, 
           'W' : 11.0, 'X' : -2.0, 'Y' :  2.0, 'Z' : -3.0, },

    'X' : {'*' : -4.0, 'A' :  0.0, 'B' : -1.0, 'C' : -2.0, 'D' : -1.0, 
           'E' : -1.0, 'F' : -1.0, 'G' : -1.0, 'H' : -1.0, 'I' : -1.0, 
           'K' : -1.0, 'L' : -1.0, 'M' : -1.0, 'N' : -1.0, 'P' : -2.0, 
           'Q' : -1.0, 'R' : -1.0, 'S' :  0.0, 'T' :  0.0, 'V' : -1.0, 
           'W' : -2.0, 'X' : -1.0, 'Y' : -1.0, 'Z' : -1.0, },

    'Y' : {'*' : -4.0, 'A' : -2.0, 'B' : -3.0, 'C' : -2.0, 'D' : -3.0, 
           'E' : -2.0, 'F' :  3.0, 'G' : -3.0, 'H' :  2.0, 'I' : -1.0, 
           'K' : -2.0, 'L' : -1.0, 'M' : -1.0, 'N' : -2.0, 'P' : -3.0, 
           'Q' : -1.0, 'R' : -2.0, 'S' : -2.0, 'T' : -2.0, 'V' : -1.0, 
           'W' :  2.0, 'X' : -1.0, 'Y' :  7.0, 'Z' : -2.0, },

    'Z' : {'*' : -4.0, 'A' : -1.0, 'B' :  1.0, 'C' : -3.0, 'D' :  1.0, 
           'E' :  4.0, 'F' : -3.0, 'G' : -2.0, 'H' :  0.0, 'I' : -3.0, 
           'K' :  1.0, 'L' : -3.0, 'M' : -1.0, 'N' :  0.0, 'P' : -1.0, 
           'Q' :  3.0, 'R' :  0.0, 'S' :  0.0, 'T' : -1.0, 'V' : -2.0, 
           'W' : -3.0, 'X' : -1.0, 'Y' : -2.0, 'Z' :  4.0, },

}
