# -*- coding: utf-8 -*-
import Gmode as Gm

gmode = Gm.Gmode()
gmode.load_data(filename="SDSSMOC/lists/MOC4_2quartile_refl.dat")
gmode.run(q1=1.3,mlim=0.15,ulim=1.0,grid=2,save='y', name="MOC2q")
gmode.evaluate(q2=1.3)
gmode.correspondence("SDSSMOC/bus_templates.pkl", "bus",q1=2.0, var = range(1,gmode.elems.shape[1])) 

