__description__ = \
"""
Simulation parameters specififed here.

The parameters are all the fraction of the (len(alphabet)**len(seq)) number of 
sequences that pass each round.
"""


all_param_sets = {
    "lwheeler_00" : {
        "company_tube_fx"       :  2.60e-03  ,  # number of sequences in tube from company
        "pipet_onto_plate_fx"   :  2.81e-05  ,  # number of sequences pipetted into well
        "protein_on_plate_fx"   :  1.50e-03  ,  # number of protein molecules in well
        "amplify_post_bind_fx"  :  6.40e-03  ,  # number of phage after amplifying
        "illumina_reads_fx"     :  2.44e-10  ,  # number of illumina reads per round
        "first_round_recov_fx"  :  9.77e-10  ,  # number of phage recovered after binding
    }
}

# -----------------------------------------------------------------------------
# NO USER EDITS REQUIRED BELOW THIS LINE 
# -----------------------------------------------------------------------------

class ParamHolder:
    """
    Mini class for holding onto parameter sets so they can be accessed via
    param_sets.lwheeler_00 (etc.)
    """

    def __init__(self,input_dict):
        """
        Map user-editable param_sets dictionary into self.__dict__
        """

        for k in input_dict.keys():
            setattr(self,k,input_dict[k])


param_sets = ParamHolder(all_param_sets)

