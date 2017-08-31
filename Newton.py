import numpy as np
import warnings
from semiconductor.material.intrinsic_carrier_density import IntrinsicCarrierDensity as ni
from semiconductor.electrical.ionisation import Ionisation as Ion
from semiconductor.material.thermal_velocity import ThermalVelocity as Vel_th
import scipy.constants as const
from semiconductor.general_functions import carrierfunctions as CF

kb = const.k / const.e


def generteCaldts(T, Ndop, doptypelist, ionauthor='Altermatt_2006_table1', vthauthor='Green_1990', niauthor='Couderc_2014', **kwarg):
    vth_e300, vth_h300 = Vel_th().update(temp=300, author=vthauthor)
    Caldts = []
    for i in range(len(T)):
        vth_e, vth_h = Vel_th().update(temp=T[i], author=vthauthor)
        Ni = ni().update(temp=T[i], author=niauthor)
        doptype = doptypelist[i]
        if doptype == 'n':
            Nidop = Ion(temp=T[i], ni_author=ni_author).update_dopant_ionisation(
                N_dop=Ndop[i], nxc=0, impurity='phosphorous', author=ionauthor)
            n0, p0 = CF.get_carriers(
                1, Nidop, 0, temp=T[i], ni_author=ni_author)
        elif doptype == 'p':
            Nidop = Ion(temp=T[i], ni_author=ni_author).update_dopant_ionisation(
                N_dop=Ndop[i], nxc=0, impurity='boron', author=ionauthor)
            n0, p0 = CF.get_carriers(
                Nidop, 0, 0, temp=T[i], ni_author=ni_author)
        dts = {
            'ni': Ni,
            'vth_e': vth_e,
            'vth_h': vth_h,
            'n0': n0,
            'p0': p0,
            'T': T[i],
            'doptype': doptype,
            'vth_e300': vth_e300,
            'vth_h300': vth_h300
        }
        Caldts.append(dts.copy())
    return Caldts


def Fm(x, m, Caldts, **kwarg):
    vth_e300 = Caldts['vth_e300']
    vth_h300 = Caldts['vth_h300']
    Et = x[0]
    tauminor = x[1]
    k = x[2]
    vth_e = Caldts['vth_e']
    vth_h = Caldts['vth_h']
    Ni = Caldts['ni']
    n0 = Caldts['n0']
    p0 = Caldts['p0']
    T = Caldts['T']
    doptype = Caldts['doptype']
    if doptype == 'n':
        res = -m + (tauminor * vth_h300 / vth_h) * (vth_h / vth_e / k *
                                                    (1 - Ni * np.exp(-Et / kb / T) / n0) - Ni * np.exp(Et / kb / T) / n0)
    elif doptype == 'p':
        res = -m + (tauminor * vth_e300 / vth_e) * (vth_e * k / vth_h) * (1 - Ni / p0 *
                                                                          np.exp(Et / kb / T)) - Ni * np.exp(-Et / kb / T) / p0 * (tauminor * vth_e300 / vth_e)
    return res


def Fb(x, b, Caldts, **kwarg):
    vth_e300 = Caldts['vth_e300']
    vth_h300 = Caldts['vth_h300']
    Et = x[0]
    tauminor = x[1]
    k = x[2]
    vth_e = Caldts['vth_e']
    vth_h = Caldts['vth_h']
    Ni = Caldts['ni']
    n0 = Caldts['n0']
    p0 = Caldts['p0']
    T = Caldts['T']
    doptype = Caldts['doptype']
    if doptype == 'n':
        res = -b + (tauminor * vth_h300 / vth_h) * (1 + Ni * np.exp(Et /
                                                                    kb / T) / n0 + Ni * np.exp(-Et / kb / T) / n0 / k * vth_h / vth_e)
    elif doptype == 'p':
        res = -b + (tauminor * vth_e300 / vth_e) * (1 + vth_e * k / vth_h / p0 * Ni *
                                                    np.exp(Et / kb / T)) + Ni * np.exp(-Et / kb / T) / p0 * (tauminor * vth_e300 / vth_e)
    return res


def dFbdtauminor(x, Caldts, **kwarg):
    vth_e300 = Caldts['vth_e300']
    vth_h300 = Caldts['vth_h300']
    Et = x[0]
    tauminor = x[1]
    k = x[2]
    vth_e = Caldts['vth_e']
    vth_h = Caldts['vth_h']
    Ni = Caldts['ni']
    n0 = Caldts['n0']
    p0 = Caldts['p0']
    T = Caldts['T']
    doptype = Caldts['doptype']
    if doptype == 'n':
        res = (1 + Ni * np.exp(Et / kb / T) / n0 + Ni * np.exp(-Et /
                                                               kb / T) / n0 / k * vth_h / vth_e) / vth_h * vth_h300
    elif doptype == 'p':
        res = (1 + vth_e * k / vth_h / p0 * Ni * np.exp(Et / kb / T) +
               Ni * np.exp(-Et / kb / T) / p0) / vth_e * vth_e300
    return res


def dFbdk(x, Caldts, **kwarg):
    vth_e300 = Caldts['vth_e300']
    vth_h300 = Caldts['vth_h300']
    Et = x[0]
    tauminor = x[1]
    k = x[2]
    vth_e = Caldts['vth_e']
    vth_h = Caldts['vth_h']
    Ni = Caldts['ni']
    n0 = Caldts['n0']
    p0 = Caldts['p0']
    T = Caldts['T']
    doptype = Caldts['doptype']
    if doptype == 'n':
        res = -1 * (tauminor * vth_h300 / vth_h) * Ni * \
            np.exp(-Et / kb / T) / n0 * vth_h / vth_e / (k**2)
    elif doptype == 'p':
        res = 1 * (tauminor * vth_e300 / vth_e) * vth_e / \
            vth_h / p0 * Ni * np.exp(Et / kb / T)
    return res


def dFbdEt(x, Caldts, **kwarg):
    vth_e300 = Caldts['vth_e300']
    vth_h300 = Caldts['vth_h300']
    Et = x[0]
    tauminor = x[1]
    k = x[2]
    vth_e = Caldts['vth_e']
    vth_h = Caldts['vth_h']
    Ni = Caldts['ni']
    n0 = Caldts['n0']
    p0 = Caldts['p0']
    T = Caldts['T']
    doptype = Caldts['doptype']
    if doptype == 'n':
        res = (tauminor * vth_h300 / vth_h) / n0 * (Ni * np.exp(Et / kb /
                                                                T) / kb / T - vth_h / vth_e / k * np.exp(-Et / kb / T) / kb / T)
    elif doptype == 'p':
        res = (tauminor * vth_e300 / vth_e) * k * vth_e / vth_h / p0 * Ni * np.exp(Et / kb / T) / \
            kb / T - Ni / p0 * np.exp(-Et / kb / T) / \
            kb / T * (tauminor * vth_e300 / vth_e)
    return res


def dFmdtauminor(x, Caldts, **kwarg):
    vth_e300 = Caldts['vth_e300']
    vth_h300 = Caldts['vth_h300']
    Et = x[0]
    tauminor = x[1]
    k = x[2]
    vth_e = Caldts['vth_e']
    vth_h = Caldts['vth_h']
    Ni = Caldts['ni']
    n0 = Caldts['n0']
    p0 = Caldts['p0']
    T = Caldts['T']
    doptype = Caldts['doptype']
    if doptype == 'n':
        res = (vth_h / vth_e / k * (1 - Ni * np.exp(-Et / kb / T) /
                                    n0) - Ni * np.exp(Et / kb / T) / n0) / vth_h * vth_h300
    elif doptype == 'p':
        res = (vth_e * k / vth_h * (1 - Ni / p0 * np.exp(Et / kb / T)) -
               Ni / p0 * np.exp(-Et / kb / T)) / vth_e * vth_e300
    return res


def dFmdk(x, Caldts, **kwarg):
    vth_e300 = Caldts['vth_e300']
    vth_h300 = Caldts['vth_h300']
    Et = x[0]
    tauminor = x[1]
    k = x[2]
    vth_e = Caldts['vth_e']
    vth_h = Caldts['vth_h']
    Ni = Caldts['ni']
    n0 = Caldts['n0']
    p0 = Caldts['p0']
    T = Caldts['T']
    doptype = Caldts['doptype']
    if doptype == 'n':
        res = -1 * (tauminor * vth_h300 / vth_h) / vth_e * vth_h * \
            (1 - Ni / n0 * np.exp(-Et / kb / T)) / (k**2)
        # print(res)
    elif doptype == 'p':
        res = (tauminor * vth_e300 / vth_e) * vth_e / \
            vth_h * (1 - Ni / p0 * np.exp(Et / kb / T))
    return res


def dFmdEt(x, Caldts, **kwarg):
    vth_e300 = Caldts['vth_e300']
    vth_h300 = Caldts['vth_h300']
    Et = x[0]
    tauminor = x[1]
    k = x[2]
    vth_e = Caldts['vth_e']
    vth_h = Caldts['vth_h']
    Ni = Caldts['ni']
    n0 = Caldts['n0']
    p0 = Caldts['p0']
    T = Caldts['T']
    doptype = Caldts['doptype']
    if doptype == 'n':
        res = -1 * (tauminor * vth_h300 / vth_h) * (-vth_h / vth_e / k / n0 * Ni *
                                                    np.exp(-Et / kb / T) / kb / T + Ni / n0 * np.exp(Et / kb / T) / kb / T)
    elif doptype == 'p':
        res = -1 * (tauminor * vth_e300 / vth_e) * k / vth_h * vth_e / p0 * Ni * np.exp(Et / kb / T) / \
            kb / T + Ni / p0 * np.exp(-Et / kb / T) / \
            kb / T * (tauminor * vth_e300 / vth_e)
    return res


def generatemb(m, b, **kwarg):
    res = []
    for i in range(len(m)):
        res.append(m[i])
        res.append(b[i])
    return res


def generatefunc(x, mb, CaldtsList, **kwarg):
    res = np.zeros((2 * len(CaldtsList), 1))
    for i in range(len(CaldtsList)):
        res[2 * i] = (Fm(x=x, m=mb[2 * i], Caldts=CaldtsList[i]))
        res[2 * i + 1] = (Fb(x=x, b=mb[2 * i + 1], Caldts=CaldtsList[i]))
    return res


def generatejacob(x, CaldtsList, **kwarg):
    res = np.zeros((2 * len(CaldtsList), 3))
    # print(res.shape)
    for i in range(len(CaldtsList)):
        res[2 * i, 0] = dFmdEt(x=x, Caldts=CaldtsList[i])
        res[2 * i, 1] = dFmdtauminor(x=x, Caldts=CaldtsList[i])
        res[2 * i, 2] = dFmdk(x=x, Caldts=CaldtsList[i])
        res[2 * i + 1, 0] = dFbdEt(x=x, Caldts=CaldtsList[i])
        res[2 * i + 1, 1] = dFbdtauminor(x=x, Caldts=CaldtsList[i])
        res[2 * i + 1, 2] = dFbdk(x=x, Caldts=CaldtsList[i])
    return res


def NewtonMethodUP(x0, mb, CaldtsList, tol=1e-12, MaxIntNum=100, ** kwarg):
    x0 = np.asarray(x0)
    x0 = x0.reshape(x0.shape[0], 1)
    Havewefoundsolution = False
    # print(x0)
    # print(mb)
    # print(CaldtsList)
    for i in range(MaxIntNum):
        if x0[0] < 0:
            x0[0] = np.random.random() * 0.6
            x0[1] = np.random.random()
            x0[2] = np.random.random()
        if x0[0] > 0.6:
            x0[0] = np.random.random() * 0.6
            x0[1] = np.random.random()
            x0[2] = np.random.random()
        if x0[1] < 0:
            x0[0] = np.random.random() * 0.6
            x0[1] = np.random.random()
            x0[2] = np.random.random()
        if x0[1] > 1:
            x0[0] = np.random.random() * 0.6
            x0[1] = np.random.random()
            x0[2] = np.random.random()
        if x0[2] > 1e3:
            x0[0] = np.random.random() * 0.6
            x0[1] = np.random.random()
            x0[2] = np.random.random()
        if x0[2] < 1e-3:
            x0[0] = np.random.random() * 0.6
            x0[1] = np.random.random()
            x0[2] = np.random.random()
        y = generatefunc(x=x0, mb=mb, CaldtsList=CaldtsList)
        # print(x0)
        yprime = generatejacob(x=x0, CaldtsList=CaldtsList)
        # print(np.linalg.matrix_rank(yprime))
        x1 = x0 - np.dot(np.linalg.pinv(yprime), y)
        if((abs(x1 - x0) <= tol * abs(x0)).all()):
            Havewefoundsolution = True
            break
        x0 = x1
    if not(Havewefoundsolution):
        warnings.warn('Result not converge')
        x0 = np.NaN
    return x0


def NewtonMethodDOWN(x0, mb, CaldtsList, tol=1e-10, MaxIntNum=100, ** kwarg):
    x0 = np.asarray(x0)
    x0 = x0.reshape(x0.shape[0], 1)
    Havewefoundsolution = False
    for i in range(MaxIntNum):
        if x0[0] < -0.6:
            x0[0] = np.random.random() * 0.6 - 0.6
            x0[1] = np.random.random()
            x0[2] = np.random.random()
        if x0[0] > 0:
            x0[0] = np.random.random() * 0.6 - 0.6
            x0[1] = np.random.random()
            x0[2] = np.random.random()
        if x0[1] < 0:
            x0[0] = np.random.random() * 0.6 - 0.6
            x0[1] = np.random.random()
            x0[2] = np.random.random()
        if x0[1] > 1:
            x0[0] = np.random.random() * 0.6 - 0.6
            x0[1] = np.random.random()
            x0[2] = np.random.random()
        if x0[2] > 1e3:
            x0[0] = np.random.random() * 0.6 - 0.6
            x0[1] = np.random.random()
            x0[2] = np.random.random()
        if x0[2] < 1e-3:
            x0[0] = np.random.random() * 0.6 - 0.6
            x0[1] = np.random.random()
            x0[2] = np.random.random()
        y = generatefunc(x=x0, mb=mb, CaldtsList=CaldtsList)
        yprime = generatejacob(x=x0, CaldtsList=CaldtsList)
        x1 = x0 - np.dot(np.linalg.pinv(yprime), y)
        # print(np.linalg.matrix_rank(yprime))
        # print(np.dot(np.linalg.inv(yprime), y))
        if((abs(x1 - x0) <= tol * abs(x0)).all()):
            Havewefoundsolution = True
            break
        x0 = x1
    if not(Havewefoundsolution):
        warnings.warn('Result not converge')
        x0 = np.NaN
    return x0
