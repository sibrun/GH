# list all expected ranks of triconnected 11-loop
import OrdinaryVariants


loops = 11
for v in range(3, 22):
    V = OrdinaryVariants.OrdinaryGVSTriconnected(v, loops, False)
    dim = V.get_dimension()
    print(f"v = {v}, dim = {dim}")
    

