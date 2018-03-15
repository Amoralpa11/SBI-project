def str_comparison_superimposer(str1, str2):
    '''This function returns a list indicating which chains
    from str1 are the same in str2 based on % of identity.
    It requires 2 arguments, which are the 2 srtuctures we
    want to compare.
    '''
    ls = 0
    ppb = PPBuilder()
    ls_count = [[0, 0], [0, 0]]
    for pp1 in ppb.build_peptides(str1):
        seq1 = pp1.get_sequence()

        if ls == 0:
            ls += 1
        else:
            ls = 0
        ls2 = 0
        for pp2 in ppb.build_peptides(str2):
            seq2 = pp2.get_sequence()
            alignment = pairwise2.align.globalxx(seq1, seq2)
            score = alignment[0][2]
            length = max(len(seq1), len(seq2))
            ident_perc = score / length  # to look at, choose longest one?

            if ident_perc > 0.95:
                chain1 = list(str1[0].get_chains())[ls]
                chain2 = list(str2[0].get_chains())[ls2]
                sup = Superimposer()
                sup.set_atoms(chain1, chain2)
                sup.apply(str2)

            else:
                ls_count[ls][ls2] = 0

            ls2 += 1
    return ls_count