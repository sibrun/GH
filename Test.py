import HairyGraphBiComplex as HGBC

if __name__ == "__main__":

    gc = HGBC.HairyBiGC(range(6, 12), 0, True, False)
    gc.build_basis()
    # gc.build_matrix()

