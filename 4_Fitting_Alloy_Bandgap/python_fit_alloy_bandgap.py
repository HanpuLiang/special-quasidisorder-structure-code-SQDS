import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


class FitAlloyBandgap():
    def __init__(self, file_name):
        self.load_DOS_file(file_name)

    energy = 0
    DOS = 0
    bandgap = 0
    cb_st = 0
    cb_step = 0
    vb_st = 0
    vb_step = 0
    N_fit_point = 21
    cb_point = 0
    vb_point = 0
    fit_cb_point = 0
    fit_vb_point = 0

    def load_DOS_file(self, file_name):
        with open(file_name, 'r') as obj:
            ct = obj.readlines()
        DOS_data = np.array([[float(item) for item in line.split()] for line in ct[1:]])
        self.energy = DOS_data[:,0]
        self.DOS = DOS_data[:,1]

    def set_fit_range(self, cb_st, cb_step, vb_st, vb_step):
        self.cb_st = cb_st
        self.cb_step = cb_step
        self.vb_st = vb_st
        self.vb_step = vb_step
    
    def fit_bandedge(self):
        fit_cb_E, fit_cb_dos = self.get_fit_data('cbm')
        fit_vb_E, fit_vb_dos = self.get_fit_data('vbm')
        
        popt_cb, pcov_cb = curve_fit(self.band_equation, fit_cb_E, fit_cb_dos, bounds=([-1e5, -1e5, -1e5], [1e5, 1e5, 1e5]))
        popt_vb, pcov_vb = curve_fit(self.band_equation, fit_vb_E, fit_vb_dos, bounds=([-1e5, -1e5, -1e5], [1e5, 1e5, 1e5]))
        E_cb_i = np.linspace(min(fit_cb_E) - 0.5, max(fit_cb_E) + 0.5, 1000)
        E_vb_i = np.linspace(min(fit_vb_E) - 0.5, max(fit_vb_E) + 0.5, 1000)
        dos_cb_i = self.band_equation(E_cb_i, *popt_cb)
        dos_vb_i = self.band_equation(E_vb_i, *popt_vb)

        self.cb_point = np.concatenate((fit_cb_E.reshape(-1,1), fit_cb_dos.reshape(-1,1)), axis=1)
        self.vb_point = np.concatenate((fit_vb_E.reshape(-1,1), fit_vb_dos.reshape(-1,1)), axis=1)
        self.fit_cb_point = np.concatenate((E_cb_i.reshape(-1,1), dos_cb_i.reshape(-1,1)), axis=1)
        self.fit_vb_point = np.concatenate((E_vb_i.reshape(-1,1), dos_vb_i.reshape(-1,1)), axis=1)
    
    def band_equation(E, A, c1, C):
        return A * (E + 5*c1*E**2 + 8*c1**2*E**3 + 4*c1**3*E**4) + C    # second order energy correction
        
    def get_bandgap(self):
        dos_0_id_cb = np.argmin(np.abs(self.fit_cb_point[:,1]))
        dos_0_id_vb = np.argmin(np.abs(self.fit_vb_point[:,1]))
        E_cbm = self.fit_cb_point[dos_0_id_cb,0]
        E_vbm = self.fit_vb_point[dos_0_id_vb,0]
        self.bandgap = E_cbm - E_vbm
        return self.bandgap

    def get_fit_data(self, sign):
        if sign == 'vbm':
            cur_st, cur_ed = self.vb_st + self.vb_step, self.vb_st
        elif sign == 'cbm':
            cur_st, cur_ed = self.cb_st, self.cb_st + self.cb_step
        fit_E_point = np.linspace(cur_st, cur_ed, self.N_fit_point)
        fit_E = np.zeros((self.N_fit_point,))
        fit_dos = np.zeros((self.N_fit_point,))
        for i, e in enumerate(fit_E_point):
            cur_e_id = np.argmin(np.abs(self.energy - e))
            fit_E[i] = self.energy[cur_e_id]
            fit_dos[i] = self.DOS[cur_e_id]
        return fit_E, fit_dos

if __name__ == '__main__':
    cell221 = FitAlloyBandgap('dos-221-lro0-sro0-k444.dat')
    cell221.set_fit_range(1, -0.3, 0.5, -1.4)
    cell221.fit_bandedge()
    Eg221 = cell221.get_bandgap()

    cell442 = FitAlloyBandgap('dos-442-lro0-sro0-k111.dat')
    cell442.set_fit_range(0.8, -0.6, 0.8, -1.5)
    cell442.fit_bandedge()
    Eg442 = cell442.get_bandgap()
    print(Eg221)
    print(Eg442)