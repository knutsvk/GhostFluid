/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
** Created: Fri Jan 22 16:44:52 2016
**      by: Qt User Interface Compiler version 4.8.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QComboBox>
#include <QtGui/QDoubleSpinBox>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QMainWindow>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QPushButton>
#include <QtGui/QSpinBox>
#include <QtGui/QStatusBar>
#include <QtGui/QToolBar>
#include <QtGui/QWidget>
#include "qcustomplot.h"

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QWidget *centralWidget;
    QPushButton *runButton;
    QCustomPlot *densityPlot;
    QCustomPlot *velocityPlot;
    QCustomPlot *pressurePlot;
    QCustomPlot *energyPlot;
    QComboBox *testCombo;
    QComboBox *schemeCombo;
    QComboBox *limiterCombo;
    QDoubleSpinBox *rhoL;
    QDoubleSpinBox *uL;
    QDoubleSpinBox *pL;
    QDoubleSpinBox *gammaL;
    QDoubleSpinBox *rhoR;
    QDoubleSpinBox *uR;
    QDoubleSpinBox *pR;
    QDoubleSpinBox *gammaR;
    QLabel *label;
    QLabel *label_4;
    QLabel *label_5;
    QSpinBox *Cells;
    QLabel *label_2;
    QLabel *label_3;
    QDoubleSpinBox *CFL;
    QDoubleSpinBox *Time;
    QLabel *label_6;
    QDoubleSpinBox *x0;
    QLabel *label_7;
    QCheckBox *saveResults;
    QLineEdit *fileName;
    QMenuBar *menuBar;
    QMenu *menuEuler_equations_in_1D;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
        MainWindow->resize(898, 665);
        centralWidget = new QWidget(MainWindow);
        centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
        runButton = new QPushButton(centralWidget);
        runButton->setObjectName(QString::fromUtf8("runButton"));
        runButton->setGeometry(QRect(767, 546, 111, 41));
        densityPlot = new QCustomPlot(centralWidget);
        densityPlot->setObjectName(QString::fromUtf8("densityPlot"));
        densityPlot->setGeometry(QRect(30, 20, 411, 211));
        velocityPlot = new QCustomPlot(centralWidget);
        velocityPlot->setObjectName(QString::fromUtf8("velocityPlot"));
        velocityPlot->setGeometry(QRect(460, 20, 411, 211));
        pressurePlot = new QCustomPlot(centralWidget);
        pressurePlot->setObjectName(QString::fromUtf8("pressurePlot"));
        pressurePlot->setGeometry(QRect(30, 250, 411, 211));
        energyPlot = new QCustomPlot(centralWidget);
        energyPlot->setObjectName(QString::fromUtf8("energyPlot"));
        energyPlot->setGeometry(QRect(460, 250, 411, 211));
        testCombo = new QComboBox(centralWidget);
        testCombo->setObjectName(QString::fromUtf8("testCombo"));
        testCombo->setGeometry(QRect(780, 490, 78, 27));
        schemeCombo = new QComboBox(centralWidget);
        schemeCombo->setObjectName(QString::fromUtf8("schemeCombo"));
        schemeCombo->setGeometry(QRect(30, 470, 78, 27));
        limiterCombo = new QComboBox(centralWidget);
        limiterCombo->setObjectName(QString::fromUtf8("limiterCombo"));
        limiterCombo->setGeometry(QRect(30, 500, 78, 27));
        rhoL = new QDoubleSpinBox(centralWidget);
        rhoL->setObjectName(QString::fromUtf8("rhoL"));
        rhoL->setGeometry(QRect(320, 490, 81, 27));
        rhoL->setDecimals(5);
        rhoL->setSingleStep(0.1);
        rhoL->setValue(1);
        uL = new QDoubleSpinBox(centralWidget);
        uL->setObjectName(QString::fromUtf8("uL"));
        uL->setGeometry(QRect(420, 490, 81, 27));
        uL->setDecimals(5);
        uL->setMinimum(-100);
        uL->setMaximum(100);
        uL->setValue(0.75);
        pL = new QDoubleSpinBox(centralWidget);
        pL->setObjectName(QString::fromUtf8("pL"));
        pL->setGeometry(QRect(520, 490, 81, 27));
        pL->setDecimals(5);
        pL->setMaximum(1e+06);
        pL->setSingleStep(0.1);
        pL->setValue(1);
        gammaL = new QDoubleSpinBox(centralWidget);
        gammaL->setObjectName(QString::fromUtf8("gammaL"));
        gammaL->setGeometry(QRect(620, 490, 62, 27));
        gammaL->setMaximum(5);
        gammaL->setSingleStep(0.1);
        gammaL->setValue(1.4);
        rhoR = new QDoubleSpinBox(centralWidget);
        rhoR->setObjectName(QString::fromUtf8("rhoR"));
        rhoR->setGeometry(QRect(320, 520, 81, 27));
        rhoR->setDecimals(5);
        rhoR->setSingleStep(0.1);
        rhoR->setValue(0.125);
        uR = new QDoubleSpinBox(centralWidget);
        uR->setObjectName(QString::fromUtf8("uR"));
        uR->setGeometry(QRect(420, 520, 81, 27));
        uR->setDecimals(5);
        uR->setMinimum(-100);
        uR->setMaximum(100);
        uR->setValue(0);
        pR = new QDoubleSpinBox(centralWidget);
        pR->setObjectName(QString::fromUtf8("pR"));
        pR->setGeometry(QRect(520, 520, 81, 27));
        pR->setDecimals(5);
        pR->setMaximum(1e+06);
        pR->setSingleStep(0.1);
        pR->setValue(0.1);
        gammaR = new QDoubleSpinBox(centralWidget);
        gammaR->setObjectName(QString::fromUtf8("gammaR"));
        gammaR->setGeometry(QRect(620, 520, 62, 27));
        gammaR->setMaximum(5);
        gammaR->setSingleStep(0.1);
        gammaR->setValue(1.4);
        label = new QLabel(centralWidget);
        label->setObjectName(QString::fromUtf8("label"));
        label->setGeometry(QRect(320, 470, 381, 17));
        label_4 = new QLabel(centralWidget);
        label_4->setObjectName(QString::fromUtf8("label_4"));
        label_4->setGeometry(QRect(690, 490, 31, 17));
        label_5 = new QLabel(centralWidget);
        label_5->setObjectName(QString::fromUtf8("label_5"));
        label_5->setGeometry(QRect(690, 520, 41, 20));
        Cells = new QSpinBox(centralWidget);
        Cells->setObjectName(QString::fromUtf8("Cells"));
        Cells->setGeometry(QRect(630, 570, 60, 27));
        Cells->setMinimum(2);
        Cells->setMaximum(2000);
        Cells->setSingleStep(100);
        Cells->setValue(100);
        label_2 = new QLabel(centralWidget);
        label_2->setObjectName(QString::fromUtf8("label_2"));
        label_2->setGeometry(QRect(610, 570, 16, 17));
        label_3 = new QLabel(centralWidget);
        label_3->setObjectName(QString::fromUtf8("label_3"));
        label_3->setGeometry(QRect(500, 570, 16, 17));
        CFL = new QDoubleSpinBox(centralWidget);
        CFL->setObjectName(QString::fromUtf8("CFL"));
        CFL->setGeometry(QRect(520, 570, 62, 27));
        CFL->setMaximum(1);
        CFL->setSingleStep(0.1);
        CFL->setValue(0.9);
        Time = new QDoubleSpinBox(centralWidget);
        Time->setObjectName(QString::fromUtf8("Time"));
        Time->setGeometry(QRect(420, 570, 62, 27));
        Time->setDecimals(3);
        Time->setMaximum(10);
        Time->setSingleStep(0.005);
        Time->setValue(0.2);
        label_6 = new QLabel(centralWidget);
        label_6->setObjectName(QString::fromUtf8("label_6"));
        label_6->setGeometry(QRect(400, 570, 16, 17));
        x0 = new QDoubleSpinBox(centralWidget);
        x0->setObjectName(QString::fromUtf8("x0"));
        x0->setGeometry(QRect(320, 570, 62, 27));
        x0->setDecimals(2);
        x0->setMaximum(1);
        x0->setSingleStep(0.05);
        x0->setValue(0.3);
        label_7 = new QLabel(centralWidget);
        label_7->setObjectName(QString::fromUtf8("label_7"));
        label_7->setGeometry(QRect(300, 570, 16, 17));
        saveResults = new QCheckBox(centralWidget);
        saveResults->setObjectName(QString::fromUtf8("saveResults"));
        saveResults->setGeometry(QRect(30, 550, 121, 22));
        saveResults->setTristate(false);
        fileName = new QLineEdit(centralWidget);
        fileName->setObjectName(QString::fromUtf8("fileName"));
        fileName->setGeometry(QRect(30, 570, 251, 27));
        MainWindow->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(MainWindow);
        menuBar->setObjectName(QString::fromUtf8("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 898, 25));
        menuEuler_equations_in_1D = new QMenu(menuBar);
        menuEuler_equations_in_1D->setObjectName(QString::fromUtf8("menuEuler_equations_in_1D"));
        MainWindow->setMenuBar(menuBar);
        mainToolBar = new QToolBar(MainWindow);
        mainToolBar->setObjectName(QString::fromUtf8("mainToolBar"));
        MainWindow->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(MainWindow);
        statusBar->setObjectName(QString::fromUtf8("statusBar"));
        MainWindow->setStatusBar(statusBar);

        menuBar->addAction(menuEuler_equations_in_1D->menuAction());

        retranslateUi(MainWindow);

        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "MainWindow", 0, QApplication::UnicodeUTF8));
        runButton->setText(QApplication::translate("MainWindow", "RUN", 0, QApplication::UnicodeUTF8));
        testCombo->clear();
        testCombo->insertItems(0, QStringList()
         << QApplication::translate("MainWindow", "Toro 1", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "Toro 2", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "Toro 3", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "Toro 4", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "Toro 5", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "MiniProject 1", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "MiniProject 2", 0, QApplication::UnicodeUTF8)
        );
        schemeCombo->clear();
        schemeCombo->insertItems(0, QStringList()
         << QApplication::translate("MainWindow", "FORCE", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "FLIC", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "SLIC", 0, QApplication::UnicodeUTF8)
        );
        limiterCombo->clear();
        limiterCombo->insertItems(0, QStringList()
         << QApplication::translate("MainWindow", "Minbee", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "Vanleer", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "Superbee", 0, QApplication::UnicodeUTF8)
        );
        label->setText(QApplication::translate("MainWindow", " Density              Velocity              Pressure            Gamma", 0, QApplication::UnicodeUTF8));
        label_4->setText(QApplication::translate("MainWindow", "Left", 0, QApplication::UnicodeUTF8));
        label_5->setText(QApplication::translate("MainWindow", "Right", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("MainWindow", "N", 0, QApplication::UnicodeUTF8));
        label_3->setText(QApplication::translate("MainWindow", "C", 0, QApplication::UnicodeUTF8));
        label_6->setText(QApplication::translate("MainWindow", "T", 0, QApplication::UnicodeUTF8));
        label_7->setText(QApplication::translate("MainWindow", "x0", 0, QApplication::UnicodeUTF8));
        saveResults->setText(QApplication::translate("MainWindow", "Save Results", 0, QApplication::UnicodeUTF8));
        fileName->setText(QApplication::translate("MainWindow", "test.out", 0, QApplication::UnicodeUTF8));
        menuEuler_equations_in_1D->setTitle(QApplication::translate("MainWindow", "Euler equations in 1D", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
