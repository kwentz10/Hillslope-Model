#Functions for Hillslope Model
#Created February 7 2016

import numpy as np
from matplotlib import pyplot as plt


class Hillslope:

    def initial_plot_conds(self, t_min, t_max, dt, x_min, x_max, dx, zu_max, zu_min, nz, zl_max, zl_min):
        self.timestep = np.arange(t_min, t_max+dt, dt).tolist() #0 to 100 days in time
        self.spacestep = np.arange(x_min, x_max, dx).tolist() #0 to 1500 m in space
        self.z_upper = np.linspace(zu_max, zu_min, nz).tolist() #initial upper terrace elevation (m)
        self.z_lower = np.linspace(zl_max, zl_min, nz).tolist() #initial lower terrace elevation (m)
        self.z = self.z_upper+self.z_lower #combine elevations arrays, list
        self.fig = plt.figure(figsize=(14,6)) #figure blueprint
        self.schematic = plt.subplot()
        self.schematic.plot(self.spacestep, self.z)
        plt.xlabel('Distance (m)')
        plt.ylabel('Marine Terrace Elevation (m)')
        plt.show()

    def dQdx_func (self, dx, k):
        self.s = np.diff(self.z)/dx #slope
        self.dQdx = (-k *(-(np.diff(self.s)/dx) )).tolist() #mass flux
        return self.dQdx

    def dzdt_func (self, dx, k, w, h, dens_rock, dens_soil):
        self.dzdt=[]
        self.dQdx=hillslope_model.dQdx_func(dx, k)
        for i in self.dQdx:   #dQdx is tiny because the slopes are identical for upper and lower terraces (so difference is tiny)--dz is identical so dzin-dzout = 0
            self.dzdt+= [(dens_rock/dens_soil * w) - (1/dens_soil*i)] #aggradation/degradation
        return self.dzdt

    def run (self, dx, k, dens_rock, dens_soil, dt):
        for t in range(len(self.timestep)):
            current_time=self.timestep[t]
            self.dQdx = hillslope_model.dQdx_func(dx, k)
            self.dzdt = hillslope_model.dzdt_func(dx, k, w, h, dens_rock, dens_soil)
            self.z=np.asarray(self.z) #to array
            self.dzdt=np.asarray(self.dzdt) #to array
            self.znew=[]
            self.znew=self.z[1: len(self.z)-1]+self.dzdt #beginning and ending points stay at same elevation (do not add dzdt)
            self.znew = self.znew.tolist()
            self.z=self.z.tolist() #to list
            plt.pause(0.2)
            self.znew_plt=[self.z[0]]+self.znew+[self.z[len(self.z)-1]]
            self.schematic.plot(self.spacestep, self.znew_plt)
            break

    def finalize(self):
        self.fig.savefig('kwentz_hs.png', bbox_inches='tight')
 

if __name__== "__main__":
        t_min=0 #yr
        t_max=350 # yr
        dt=1 #yr
        x_min=0.0 #m
        x_max=1500.0 #m
        dx=10 #m
        zu_max=100 #m
        zu_min=90 #m
        zl_max=60 #m
        zl_min=45 #m
        nz=(x_max/dx)/2 #total number of z data points; divide by two because I do upper and lower terrace (same dims as spacestep)
        kappa=0.001 #soil creep k ranges from 0.001 to 0.04 m2yr-1
        dens_rock=2650.0 #gcm-3
        dens_soil=2000.0 #gcm-3
        k=dens_soil*kappa
        w=10**(-6) #weathering rate myr-1(10-60)
        h=0.05 #m
    
        hillslope_model=Hillslope()

        #initilialize
        hillslope_model.initial_plot_conds(t_min, t_max, dt, x_min, x_max, dx, zu_max, zu_min, nz, zl_max, zl_min)
        #run
        hillslope_model.run(dx, k, dens_rock, dens_soil, dt)
        #finalize
        hillslope_model.finalize()

   

