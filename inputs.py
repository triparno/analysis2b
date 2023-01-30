import numpy as np

model = "B-L_boson"

frame = 'lab'

coup_scal = 1  # 0.1

ch = 'el'

lum = 1E4

Nevents = 3


cl_st = 0.1

wd_dir = "../Data/Width_tot/"
cnt_file = ("../Data/cont/"
            + model + "_" + ch + f"_{lum}_{Nevents}")


mrange = np.logspace(-2, 0, 500)
grange = np.logspace(-8, -1, 10)

gscan_step = 20
mscan_step = 10
Eph_st = 10
Delz_st = 100
