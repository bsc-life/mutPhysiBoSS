import maboss
from maboss import Node, Network
from maboss import Simulation

sizek_bnd = "./config/sizek_BN/Sizek_model.bnd"
sizek_cfg = "./config/sizek_BN/Sizek_model.cfg"

sizek_model = maboss.load(sizek_bnd, sizek_cfg)


# running the model
run_model = sizek_model.run()