"""Generate sample images of graphs to be used on the webpage."""

import SpecialGraphs
import Parameters
import OrdinaryGraphComplex
import HairyGraphComplex
import StoreLoad

imgdir = Parameters.web_dir + "/img/"
StoreLoad.makedirs(imgdir)

# Commutatie graph complex
G = SpecialGraphs.wheel_graph(7)
VS = OrdinaryGraphComplex.OrdinaryGVS(8, 7, False)
# p = VS.plot_graph(G)
p = G.plot(vertex_labels=False, transparent=True)
p.save(imgdir + "wheel7.png")


# hairy graph complex
G = SpecialGraphs.hedgehog_graph(5)
VS = HairyGraphComplex.HairyGraphVS(5, 1, 5, False, False)
p = G.plot(vertex_labels=False, transparent=True)
p.save(imgdir + "hedgehog5.png")
