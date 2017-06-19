import numpy as np
import sys
import os
from PyQt5.QtWidgets import QApplication, QMainWindow, QFileDialog, QListWidgetItem, QLineEdit, QLabel, QRadioButton, QGridLayout, QPushButton, QAction, QActionGroup, QMenu, QInputDialog, qApp, QVBoxLayout
from PyQt5 import QtGui, QtWidgets
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from uiLSanalysis import Ui_IDLS_analyzer_dpss
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt
from semiconductor.electrical.ionisation import Ionisation as Ion
from semiconductor.material.intrinsic_carrier_density import IntrinsicCarrierDensity as NI
from semiconductor.material.thermal_velocity import ThermalVelocity as the_vel
from semiconductor.material.densityofstates import DOS as dos
from semiconductor.general_functions import carrierfunctions as CF
from scipy.optimize import curve_fit
import scipy.constants as const


class fitres(object):
    def __init__(self, name=None, Ndop=None, temp=None, doptype=None, m=None, b=None, parentid=None, brotherid=None, uid=None, **kwarg):
        self.name = name
        self.Ndop = Ndop
        self.temp = temp
        self.doptype = doptype
        self.m = m
        self.b = b
        self.uid = uid
        self.parentid = parentid
        self.brotherid = brotherid


class LSana_dpss(QMainWindow, Ui_IDLS_analyzer_dpss):

    def __init__(self, rawdata, uidlist, mainwindow, parent=None, **kwarg):

        QMainWindow.__init__(self, parent)
        self.setupUi(self)

        self.Rawdat = rawdata
        self.uidlist = uidlist
        self.mainwindow = mainwindow
        self.currentfitres = None
        self.fitreslist = []
        self.currentvroid = 0

        self.listWidget_fit.setSelectionMode(
            QtWidgets.QAbstractItemView.ExtendedSelection)
        self.listWidget_fit.itemSelectionChanged.connect(self.enablewidgts)
        self.listWidget_fit.clicked.connect(
            self.listWidget_fitres.clearSelection)
        self.listWidget_fit.itemSelectionChanged.connect(self.updateplot)
        self.listWidget_fit.itemDoubleClicked.connect(self.opensetparam2)
        self.listWidget_fit.setContextMenuPolicy(Qt.CustomContextMenu)
        self.listWidget_fit.customContextMenuRequested.connect(
            self.itemrightclicked2)

        self.listWidget_fitres.setSelectionMode(
            QtWidgets.QAbstractItemView.ExtendedSelection)
        self.listWidget_fitres.itemSelectionChanged.connect(self.enablewidgts)
        self.listWidget_fitres.clicked.connect(
            self.listWidget_fit.clearSelection)
        self.listWidget_fitres.itemSelectionChanged.connect(self.updateplot)
        self.listWidget_fitres.setContextMenuPolicy(Qt.CustomContextMenu)
        self.listWidget_fitres.customContextMenuRequested.connect(
            self.itemrightclicked)

        self.comboBox_plotopt2.currentIndexChanged.connect(self.updateplot)

        self.lineEdit_fitini1.textChanged[str].connect(self.enablewidgts)
        self.lineEdit_fitini2.textChanged[str].connect(self.enablewidgts)

        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.verticalLayout.addWidget(self.canvas)
        self.verticalLayout.addWidget(self.toolbar)
        self.ax1 = self.figure.add_subplot(111)
        self.ax1.set_xlabel(r'Excess carrier density $[cm^{-3}]$')
        self.ax1.set_ylabel('Lifetime [s]')
        self.figure.tight_layout()

        self.pushButton_back.clicked.connect(self.close)

        self.pushButton_set2.clicked.connect(self.opensetparam2)
        self.pushButton_set2.setEnabled(False)

        self.pushButton_NT.clicked.connect(self.NT)
        self.pushButton_NT.setEnabled(False)

        self.pushButton_dpss.clicked.connect(self.dpss)
        self.pushButton_dpss.setEnabled(False)

        self.pushButton_twofit.clicked.connect(self.twofit)
        self.pushButton_twofit.setEnabled(False)

        self.pushButton_1dplot.clicked.connect(self.onedplot)
        self.pushButton_1dplot.setEnabled(False)

        self.pushButton_onefit.clicked.connect(self.onefit)
        self.pushButton_onefit.setEnabled(False)

        self.pushButton_2dplot.clicked.connect(self.twodplot)
        self.pushButton_2dplot.setEnabled(False)

        self.pushButton_simufit.clicked.connect(self.simufit)
        self.pushButton_simufit.setEnabled(False)

        self.pushButton_acceptfit.clicked.connect(self.acceptfit)
        self.pushButton_acceptfit.setEnabled(False)

        self.pushButton_delfitres.clicked.connect(self.delfit)
        self.pushButton_delfitres.setEnabled(False)

        self.radioButton_ind.toggled.connect(self.enablewidgts)

        for data in self.Rawdat:
            if data.uid in self.uidlist:
                dataitem = QListWidgetItem(parent=self.listWidget_fit)
                dataitem.setText(data.name)
                dataitem.setData(32, data.uid)
                self.listWidget_fit.addItem(dataitem)

        availablDOS = dos().available_models()
        availablNI = NI().available_models()
        availablIon = Ion().available_models()
        availablthervel = the_vel().available_models()

        self.menubar = self.menuBar()
        self.choosmodel = self.menubar.addMenu('Choose your models')
        self.dosmodel = self.choosmodel.addMenu('Density of states')
        self.nimodel = self.choosmodel.addMenu('ni models')
        self.ionmodel = self.choosmodel.addMenu('Ionisation')
        self.themodel = self.choosmodel.addMenu('thermal velocity')

        self.dosgroup = QActionGroup(self, exclusive=True)
        self.nigroup = QActionGroup(self, exclusive=True)
        self.iongroup = QActionGroup(self, exclusive=True)
        self.thegroup = QActionGroup(self, exclusive=True)

        for dosmodel in availablDOS:
            a = self.dosgroup.addAction(QAction(dosmodel, checkable=True))
            if dosmodel == 'Couderc_2014':
                a.setChecked(True)
            self.dosmodel.addAction(a)
        for nimodel in availablNI:
            a = self.nigroup.addAction(QAction(nimodel, checkable=True))
            if nimodel == 'Couderc_2014':
                a.setChecked(True)
            self.nimodel.addAction(a)
        for ionmodel in availablIon:
            a = self.iongroup.addAction(QAction(ionmodel, checkable=True))
            if ionmodel == 'Altermatt_2006_table1':
                a.setChecked(True)
            self.ionmodel.addAction(a)
        for themodel in availablthervel:
            a = self.thegroup.addAction(QAction(themodel, checkable=True))
            if themodel == 'Green_1990':
                a.setChecked(True)
            self.themodel.addAction(a)

    def closeEvent(self, event):
        self.mainwindow.show()

    def updateplot(self):
        self.ax1.clear()
        self.ax1.grid()
        for item in self.listWidget_fit.selectedItems():
            for data in self.Rawdat:
                if data.uid == item.data(32):
                    if self.comboBox_plotopt2.currentIndex() == 0:
                        X = self.getXY(nxc=data.nxc, T=data.temp,
                                       Ndop=data.Ndop, doptype=data.doptype)
                        self.ax1.set_xlabel(r'X/Y')
                        self.ax1.set_ylabel('Lifetime [s]')
                        self.ax1.plot(X, data.tau, '.', label=data.name)
                    elif self.comboBox_plotopt2.currentIndex() == 1:
                        self.ax1.semilogx()
                        self.ax1.set_ylabel('Lifetime [s]')
                        self.ax1.set_xlabel(
                            r'Excess carrier density $[cm^{-3}]$')
                        self.ax1.plot(data.nxc, data.tau, '.', label=data.name)
                    elif self.comboBox_plotopt2.currentIndex() == 2:
                        self.ax1.semilogx()
                        self.ax1.set_ylabel('Inverse Lifetime [s-1]')
                        self.ax1.set_xlabel(
                            r'Excess carrier density $[cm^{-3}]$')
                        self.ax1.plot(data.nxc, 1. / data.tau,
                                      '.', label=data.name)

        for item in self.listWidget_fitres.selectedItems():
            for fitres in self.fitreslist:
                if fitres.uid == item.data(32):
                    sigleid = fitres.uid
                    broid = fitres.brotherid
                    praid = fitres.parentid

            for data in self.Rawdat:
                if data.uid == praid:
                    X = self.getXY(nxc=data.nxc, T=data.temp,
                                   Ndop=data.Ndop, doptype=data.doptype)
                    if self.comboBox_plotopt2.currentIndex() == 0:
                        self.ax1.plot(X, data.tau, '.', label=data.name)
                        if broid == None:
                            for fitres2 in self.fitreslist:
                                if fitres2.uid == sigleid:
                                    self.ax1.plot(
                                        X, X * fitres2.m + fitres2.b, label=fitres2.name)
                        else:
                            twolist = []
                            for fitres3 in self.fitreslist:
                                if fitres3.brotherid == broid:
                                    twolist.append(fitres3)
                                    self.ax1.plot(
                                        X, X * fitres3.m + fitres3.b, label=fitres3.name)
                            if len(twolist) == 2:
                                self.ax1.plot(X, 1 / (1. / (X * twolist[0].m + twolist[0].b) + 1 / (
                                    X * twolist[1].m + twolist[1].b)), label='Two defects fit' + data.name)
                        self.ax1.set_xlabel(r'X/Y')
                        self.ax1.set_ylabel('Lifetime [s]')

                    elif self.comboBox_plotopt2.currentIndex() == 1:
                        self.ax1.semilogx()
                        self.ax1.set_ylabel('Lifetime [s]')
                        self.ax1.set_xlabel(
                            r'Excess carrier density $[cm^{-3}]$')
                        self.ax1.plot(data.nxc, data.tau, '.', label=data.name)
                        if broid == None:
                            for fitres2 in self.fitreslist:
                                if fitres2.uid == sigleid:
                                    self.ax1.plot(
                                        data.nxc, X * fitres2.m + fitres2.b, label=fitres2.name)
                        else:
                            twolist = []
                            for fitres3 in self.fitreslist:
                                if fitres3.brotherid == broid:
                                    twolist.append(fitres3)
                                    self.ax1.plot(
                                        data.nxc, X * fitres3.m + fitres3.b, label=fitres3.name)
                            if len(twolist) == 2:
                                self.ax1.plot(data.nxc, 1 / (1. / (X * twolist[0].m + twolist[0].b) + 1 / (
                                    X * twolist[1].m + twolist[1].b)), label='Two defects fit' + data.name)
                    elif self.comboBox_plotopt2.currentIndex() == 2:
                        self.ax1.semilogx()
                        self.ax1.set_ylabel('Inverse Lifetime [s-1]')
                        self.ax1.set_xlabel(
                            r'Excess carrier density $[cm^{-3}]$')
                        self.ax1.plot(data.nxc, 1. / data.tau,
                                      '.', label=data.name)
                        if broid == None:
                            for fitres2 in self.fitreslist:
                                if fitres2.uid == sigleid:
                                    self.ax1.plot(
                                        data.nxc, 1. / (X * fitres2.m + fitres2.b), label=fitres2.name)
                        else:
                            twolist = []
                            for fitres3 in self.fitreslist:
                                if fitres3.brotherid == broid:
                                    twolist.append(fitres3)
                                    self.ax1.plot(
                                        data.nxc, 1. / (X * fitres3.m + fitres3.b), label=fitres3.name)
                            if len(twolist) == 2:
                                self.ax1.plot(data.nxc, (1. / (X * twolist[0].m + twolist[0].b) + 1 / (
                                    X * twolist[1].m + twolist[1].b)), label='Two defects fit' + data.name)

        self.ax1.legend(loc=0)
        self.canvas.draw()
        self.figure.tight_layout()

    def enablewidgts(self):
        if len(self.listWidget_fit.selectedItems()) == 1:
            self.pushButton_set2.setEnabled(True)
        else:
            self.pushButton_set2.setEnabled(False)
        if len(self.listWidget_fit.selectedItems()) == 1 and self.radioButton_ind.isChecked() == True:
            self.pushButton_onefit.setEnabled(True)
            try:
                a = float(self.lineEdit_fitini1.text())
                b = float(self.lineEdit_fitini2.text())
                if a < 100 and a > 0 and b < 100 and b > 0:
                    self.pushButton_twofit.setEnabled(True)
                else:
                    self.pushButton_twofit.setEnabled(False)
            except ValueError:
                self.pushButton_twofit.setEnabled(False)
        else:
            self.pushButton_onefit.setEnabled(False)
            self.pushButton_twofit.setEnabled(False)
        if len(self.listWidget_fitres.selectedItems()) != 0 and self.radioButton_ind.isChecked() == True:
            self.pushButton_delfitres.setEnabled(True)
            self.pushButton_dpss.setEnabled(True)
        else:
            self.pushButton_delfitres.setEnabled(False)
            self.pushButton_dpss.setEnabled(False)

    def opensetparam2(self):
        self.dialog = QtWidgets.QDialog()
        grid = QGridLayout(self.dialog)
        self.ntype = QRadioButton('n-type')
        self.ptype = QRadioButton('p-type')
        self.dop = QLineEdit()
        self.temp = QLineEdit()
        Ldop = QLabel('Ndop (cm-3)')
        Ltemp = QLabel('Temp (K)')
        self.ok = QPushButton('OK')
        self.ok.setEnabled(False)
        grid.addWidget(self.ntype, 0, 0)
        grid.addWidget(self.ptype, 0, 1)
        grid.addWidget(Ldop, 1, 0)
        grid.addWidget(self.dop, 1, 1)
        grid.addWidget(Ltemp, 2, 0)
        grid.addWidget(self.temp, 2, 1)
        grid.addWidget(self.ok, 3, 1)
        for data in self.Rawdat:
            if data.uid == self.listWidget_fit.selectedItems()[0].data(32):
                if data.Ndop is not None:
                    self.dop.setText('{:e}'.format(data.Ndop))
                if data.temp is not None:
                    self.temp.setText(str(data.temp))
                if data.doptype == 'n':
                    self.ntype.setChecked(True)
                if data.doptype == 'p':
                    self.ptype.setChecked(True)
        self.dop.textChanged[str].connect(self.checkparam2)
        self.temp.textChanged[str].connect(self.checkparam2)
        self.ntype.toggled.connect(self.checkparam2)
        self.ok.clicked.connect(self.setparam2)
        self.dialog.exec_()

    def checkparam2(self):
        try:
            float(self.dop.text())
            float(self.temp.text())
            self.ok.setEnabled(True)
        except ValueError:
            self.ok.setEnabled(False)

    def setparam2(self):
        for data in self.Rawdat:
            if data.uid == self.listWidget_fit.selectedItems()[0].data(32):
                data.Ndop = float(self.dop.text())
                data.temp = float(self.temp.text())
                if self.ntype.isChecked():
                    data.doptype = 'n'
                if self.ptype.isChecked():
                    data.doptype = 'p'
        self.dialog.close()

    def itemrightclicked2(self):
        if len(self.listWidget_fit.selectedItems()) == 1:
            menu = QMenu(self.listWidget_fit)
            Rset = menu.addAction('Set')
            Rset.triggered.connect(self.opensetparam2)
            Rrename = menu.addAction('Rename')
            Rrename.triggered.connect(self.changename2)
            menu.popup(QtGui.QCursor.pos())

    def changename2(self):
        if len(self.listWidget_fit.selectedItems()) == 1:
            item = self.listWidget_fit.selectedItems()[0]
            text, ok = QInputDialog.getText(
                self, 'Rename', 'Enter a new name:')
            if ok:
                item.setText(text)
                for data in self.Rawdat:
                    if data.uid == item.data(32):
                        data.name = text
                self.updateplot()

    def getXY(self, nxc, T, Ndop, doptype, **kwarg):
        for action in self.ionmodel.actions():
            if action.isChecked():
                ionauthor = action.text()
        if doptype == 'n':
            Nidop = Ion(temp=T).update_dopant_ionisation(
                N_dop=Ndop, nxc=nxc, impurity='phosphorous', author=ionauthor)
            ne, nh = CF.get_carriers(0, Nidop, nxc, temp=T)
            X = np.divide(nxc, ne)
        elif doptype == 'p':
            Nidop = Ion(temp=T).update_dopant_ionisation(
                N_dop=Ndop, nxc=nxc, impurity='boron', author=ionauthor)
            ne, nh = CF.get_carriers(Nidop, 0, nxc, temp=T)
            X = np.divide(nxc, nh)
        return X

    def twoline(self, x, s1, t1, s2, t2):
        return np.divide(1, np.divide(1, s1 * x + t1) + np.divide(1, s2 * x + t2))

    def twolinenew(self, x, s1, t1, s2, t2):
        return np.divide(1, np.divide(1, (s1 - t1) * x + t1) + np.divide(1, (s2 - t2) * x + t2))

    def twofit(self):
        self.pushButton_acceptfit.setEnabled(True)
        self.updateplot()
        item = self.listWidget_fit.selectedItems()[0]
        for data in self.Rawdat:
            if data.uid == item.data(32):
                data2fit = data
        X = self.getXY(nxc=data2fit.nxc, T=data2fit.temp,
                       Ndop=data2fit.Ndop, doptype=data2fit.doptype)
        index1 = int(X.shape[0] * float(self.lineEdit_fitini1.text()) / 100)
        index2 = int(X.shape[0] * float(self.lineEdit_fitini2.text()) / 100)
        s1, t1 = np.polyfit(X[:index1], data2fit.tau[:index1], 1)
        s2, t2 = np.polyfit(X[index2:], data2fit.tau[index2:], 1)
        try:
            popt, pcov = curve_fit(
                self.twoline, X, data2fit.tau, p0=[s1, t1, s2, t2])
            if popt[3] + popt[2] > popt[1] + popt[0]:
                m = [popt[0], popt[2]]
                b = [popt[1], popt[3]]
            elif popt[3] + popt[2] <= popt[1] + popt[0]:
                m = [popt[2], popt[0]]
                b = [popt[3], popt[1]]
            fittau = self.twoline(X, *popt)

            if b[0] <= 0 or b[1] <= 0 or m[0] + b[0] <= 0 or m[1] + b[1] <= 0:
                t11 = np.amax([t1, 0])
                s11 = np.amax([s1 + t1, 0])
                # s22, t22 = np.polyfit(X[index2:], tau[index2:], 1)
                t22 = np.amax([t2, 0])
                s22 = np.amax([s2 + t2, 0])
                popt, pcov = curve_fit(self.twolinenew, X, data2fit.tau, p0=[s11, t11, s22, t22], bounds=(
                    0, [np.inf, np.inf, np.inf, np.inf]))
                if popt[2] > popt[0]:
                    m = [popt[0] - popt[1], popt[2] - popt[3]]
                    b = [popt[1], popt[3]]
                elif popt[2] <= popt[0]:
                    m = [popt[2] - popt[3], popt[0] - popt[1]]
                    b = [popt[3], popt[1]]
                fittau = self.twolinenew(X, *popt)

            domi = m[0] * X + b[0]
            secon = m[1] * X + b[1]
            residuals = data2fit.tau - fittau
            res = np.sum(residuals**2)
            ss_tot = np.sum((data2fit.tau - np.mean(data2fit.tau))**2)
            r_squared = 1 - (res / ss_tot)
            if self.comboBox_plotopt2.currentIndex() == 0:
                self.ax1.plot(X, domi, label='Dominant')
                self.ax1.plot(X, secon, label='Secondary')
                self.ax1.plot(X, fittau, label='Two defects fit')
            elif self.comboBox_plotopt2.currentIndex() == 1:
                self.ax1.plot(data2fit.nxc, domi, label='Dominant')
                self.ax1.plot(data2fit.nxc, secon, label='Secondary')
                self.ax1.plot(data2fit.nxc, fittau,
                              label='1 defect fit' + data2fit.name)
            elif self.comboBox_plotopt2.currentIndex() == 2:
                self.ax1.plot(data2fit.nxc, 1. / domi, label='Dominant')
                self.ax1.plot(data2fit.nxc, 1. / secon, label='Secondary')
                self.ax1.plot(data2fit.nxc, 1. / fittau,
                              label='1 defect fit' + data2fit.name)
            self.ax1.set_title(
                'Residual={:e}, R2={:f}'.format(res, r_squared), fontsize=10)
            self.ax1.legend(loc=0)
            self.canvas.draw()
            self.figure.tight_layout()
            aa = self.generatebroid()
            self.currentfitres = [fitres(
                name='1st defect_' + data2fit.name, Ndop=data2fit.Ndop, temp=data2fit.temp, doptype=data2fit.doptype, m=m[0], b=b[0], brotherid=aa, parentid=data2fit.uid), fitres(
                    name='2nd defect_' + data2fit.name, Ndop=data2fit.Ndop, temp=data2fit.temp, doptype=data2fit.doptype, m=m[1], b=b[1], brotherid=aa, parentid=data2fit.uid)]

        except RuntimeError:
            self.ax1.set_title(
                'Cannot fit, try to change the initial fit range or go back to data process', fontsize=10)

    def generatebroid(self):
        self.currentvroid += 1
        return (self.currentvroid - 1)

    def onefit(self):
        self.pushButton_acceptfit.setEnabled(True)
        self.updateplot()
        item = self.listWidget_fit.selectedItems()[0]
        for data in self.Rawdat:
            if data.uid == item.data(32):
                data2fit = data
        X = self.getXY(nxc=data2fit.nxc, T=data2fit.temp,
                       Ndop=data2fit.Ndop, doptype=data2fit.doptype)
        p = np.polyfit(X, data2fit.tau, 1)
        m = p[0]
        b = p[1]
        fittau = m * X + b
        residuals = fittau - data2fit.tau
        res = np.sum(residuals**2)
        ss_tot = np.sum((data2fit.tau - np.mean(data2fit.tau))**2)
        r_squared = 1 - (res / ss_tot)
        if self.comboBox_plotopt2.currentIndex() == 0:
            self.ax1.plot(
                X, fittau, label='1 defect fit' + data2fit.name)
        elif self.comboBox_plotopt2.currentIndex() == 1:
            self.ax1.plot(data2fit.nxc, fittau,
                          label='1 defect fit' + data2fit.name)
        elif self.comboBox_plotopt2.currentIndex() == 2:
            self.ax1.plot(data2fit.nxc, 1. / fittau,
                          label='1 defect fit' + data2fit.name)
        self.ax1.set_title(
            'Residual={:e}, R2={:f}'.format(res, r_squared), fontsize=10)
        self.ax1.legend(loc=0)
        self.canvas.draw()
        self.figure.tight_layout()
        self.currentfitres = [fitres(
            name='single defect_' + data2fit.name, Ndop=data2fit.Ndop, temp=data2fit.temp, doptype=data2fit.doptype, m=m, b=b, parentid=data2fit.uid)]

    def acceptfit(self):
        self.pushButton_acceptfit.setEnabled(False)
        self.fitreslist.extend(self.currentfitres)
        self.updateuid()
        self.updatefitreslist()

    def updatefitreslist(self):
        self.listWidget_fitres.clear()
        for data in self.fitreslist:
            dataitem = QListWidgetItem(parent=self.listWidget_fitres)
            dataitem.setText(data.name)
            dataitem.setData(32, data.uid)
            self.listWidget_fitres.addItem(dataitem)

    def updateuid(self):
        for fitres in self.fitreslist:
            fitres.uid = self.fitreslist.index(fitres)

    def onedplot(self):
        pass

    def twodplot(self):
        pass

    def simufit(self):
        pass

    def itemrightclicked(self):
        if len(self.listWidget_fitres.selectedItems()) > 0:
            menu = QMenu(self.listWidget_fitres)
            Rdel = menu.addAction('Delete')
            Rdel.triggered.connect(self.delfit)
            if len(self.listWidget_fitres.selectedItems()) == 1:
                Rrename = menu.addAction('Rename')
                Rrename.triggered.connect(self.changename)
            menu.popup(QtGui.QCursor.pos())

    def delfit(self):
        for item in self.listWidget_fitres.selectedItems():
            for fitres in self.fitreslist:
                if fitres.uid == item.data(32):
                    self.fitreslist.remove(fitres)
        self.updateuid()
        self.updatefitreslist()
        self.updateplot()

    def changename(self):
        if len(self.listWidget_fitres.selectedItems()) == 1:
            item = self.listWidget_fitres.selectedItems()[0]
            text, ok = QInputDialog.getText(
                self, 'Rename', 'Enter a new name:')
            if ok:
                item.setText(text)
                for data in self.fitreslist:
                    if data.uid == item.data(32):
                        data.name = text
        self.updateplot()

    def dpss(self):
        self.dialogdpss = QtWidgets.QDialog()
        self.dialogdpss.setGeometry(400, 100, 600, 800)
        self.dpssvlayout = QVBoxLayout(self.dialogdpss)
        self.NTbuttom = QPushButton('Display/hide Newton method result')
        self.NTbuttom.clicked.connect(self.NT)
        self.dpssvlayout.addWidget(self.NTbuttom)
        self.figure2 = plt.figure()
        self.canvas2 = FigureCanvas(self.figure2)
        self.toolbar2 = NavigationToolbar(self.canvas2, self)
        self.dpssvlayout.addWidget(self.canvas2)
        self.dpssvlayout.addWidget(self.toolbar2)
        self.ax3 = self.figure2.add_subplot(211)
        self.ax4 = self.figure2.add_subplot(212, sharex=self.ax3)
        self.ax4.set_xlabel(r'$E_{t}-E_{i}$ [eV]')
        self.ax4.set_ylabel('k')
        self.ax3.set_ylabel(r'$\tau_{minor}$ [s]')
        self.ax3.semilogy()
        self.ax4.semilogy()
        self.figure2.tight_layout()
        for item in self.listWidget_fitres.selectedItems():
            for fitres in self.fitreslist:
                if fitres.uid == item.data(32):
                    EtRange, k, TauMinor = self.CalDPSS(
                        m=fitres.m, b=fitres.b, T=fitres.temp, Ndop=fitres.Ndop, doptype=fitres.doptype)
                    self.ax3.plot(EtRange, TauMinor, label=fitres.name)
                    self.ax4.plot(EtRange, k, label=fitres.name)
        self.ax4.legend(loc=0)
        self.dialogdpss.exec_()

    def CalDPSS(self, m, b, T, Ndop, doptype, **kwarg):
        kb = const.k / const.e
        for action in self.ionmodel.actions():
            if action.isChecked():
                ionauthor = action.text()
        for action in self.nimodel.actions():
            if action.isChecked():
                niauthor = action.text()
        for action in self.themodel.actions():
            if action.isChecked():
                theauthor = action.text()
        vth_e300, vth_h300 = the_vel().update(temp=300, author=theauthor)
        EtRange = np.linspace(-0.6, 0.6, 401)
        k = np.ones(EtRange.shape[0])
        TauMinor = np.ones(EtRange.shape[0])
        vth_e, vth_h = the_vel().update(temp=T, author=theauthor)
        if doptype == 'n':
            Nidop = Ion(temp=T).update_dopant_ionisation(
                N_dop=Ndop, nxc=0, impurity='phosphorous', author=ionauthor)
            n0, p0 = CF.get_carriers(0, Nidop, 0, temp=T)
        elif doptype == 'p':
            Nidop = Ion(temp=T).update_dopant_ionisation(
                N_dop=Ndop, nxc=0, impurity='boron', author=ionauthor)
            n0, p0 = CF.get_carriers(Nidop, 0, 0, temp=T)
        for j in range(0, EtRange.shape[0]):
            n1 = NI().update(temp=T, author=niauthor) * \
                np.exp(EtRange[j] / kb / T)
            p1 = NI().update(temp=T, author=niauthor) * \
                np.exp(-EtRange[j] / kb / T)
            C = np.divide(m, (m + b))
            if doptype == 'n':
                Q = (1 - p1 / n0 - C) / (C + n1 / n0)
            elif doptype == 'p':
                Q = (C + p1 / p0) / (1 - n1 / p0 - C)
            if Q < 0:
                Q = np.nan
            k[j] = Q * vth_h / vth_e
            if doptype == 'n':
                TauMinor[j] = b / (1 + n1 / n0 + p1 / n0 /
                                   Q) * vth_h / vth_h300
            elif doptype == 'p':
                TauMinor[j] = b / (1 + Q * n1 / p0 + p1 /
                                   p0) * vth_e / vth_e300
        return EtRange, k, TauMinor

    def NT(self):
        pass


if __name__ == '__main__':
    app = QApplication(sys.argv)
    LSana_dpss = LSana_dpss(rawdata=None, uidlist=None, mainwindow=LSana_dpss)
    LSana_dpss.show()
    sys.exit(app.exec_())
