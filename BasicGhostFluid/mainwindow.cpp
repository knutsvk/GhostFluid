#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cstring>
#include <vector>

#include <QDebug>
#include <QString>
#include <QVector>

#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "Functions.h"
#include "Types.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->limiterCombo->setVisible(0);
    ui->fileName->setVisible(0);
    ui->schemeCombo->setCurrentIndex(2);
    ui->limiterCombo->setCurrentIndex(2);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_schemeCombo_currentIndexChanged()
{
    int n = ui->schemeCombo->currentIndex();
    ui->limiterCombo->setVisible(n>0);
}

void MainWindow::on_testCombo_currentIndexChanged()
{
    int n = ui->testCombo->currentIndex();
    double rho_L;       // Left density
    double u_L;         // Left velocity
    double p_L;         // Left pressure
    double gamma_L;     // Left gamma
    double rho_R;       // Right density
    double u_R;         // Right velocity
    double p_R;         // Right pressure
    double gamma_R;     // Right gamma
    double tStop;       // Final time
    double x_0;          
    initiateTestCase(n+1, tStop, x_0, rho_L, u_L, p_L, gamma_L,
            rho_R, u_R, p_R, gamma_R);
    ui->rhoL->setValue(rho_L);
    ui->uL->setValue(u_L);
    ui->pL->setValue(p_L);
    ui->gammaL->setValue(gamma_L); 
    ui->rhoR->setValue(rho_R);
    ui->uR->setValue(u_R);
    ui->pR->setValue(p_R);
    ui->gammaR->setValue(gamma_R);
    ui->Time->setValue(tStop);
    ui->x0->setValue(x_0);
    ui->CFL->setValue(0.9);
    ui->Cells->setValue(100);
}

void MainWindow::on_saveResults_stateChanged(int state)
{
    ui->fileName->setVisible(state!=0);
}

void PlotGraph(QCustomPlot* plot, QVector<double> X, 
        QVector<double> Y, QString label)
{
    plot->addGraph();
    plot->graph(0)->setScatterStyle(QCPScatterStyle(
                QCPScatterStyle::ssCircle, QPen(Qt::red,0.5), 
                QBrush(Qt::red), 3));
    plot->graph(0)->setLineStyle(QCPGraph::lsNone);
    plot->graph(0)->setData(X,Y);
    double min = *std::min_element(Y.begin(),Y.end());
    double max = *std::max_element(Y.begin(),Y.end());

    double ymin = min - 0.05*(max-min);
    double ymax = max + 0.05*(max-min);
    if (ymax-ymin < 1)
    {
        ymin = (ymax+ymin)/2 - 0.5;
        ymax = (ymax+ymin)/2 + 0.5;
    }

    plot->xAxis->setRange(0,1);
    plot->yAxis->setRange(ymin, ymax);
    plot->xAxis->setLabel("x");
    plot->yAxis->setLabel(label);
    plot->replot();
}

void MainWindow::on_runButton_clicked()
{
    int N = ui->Cells->value();
    int schemeChoice = ui->schemeCombo->currentIndex();
    int limitChoice = ui->limiterCombo->currentIndex();

    bool save = ui->saveResults->isChecked(); 

    double rho_L = ui->rhoL->value();   // Left density
    double u_L = ui->uL->value();       // Left velocity
    double p_L = ui->pL->value();       // Left pressure
    double gamma_L = ui->gammaL->value(); 
    double rho_R = ui->rhoR->value();   // Right density
    double u_R = ui->uR->value();       // Right velocity
    double p_R = ui->pR->value();       // Right pressure
    double gamma_R = ui->gammaR->value();
    double a_L;                         // Left sound speed
    double a_R;                         // Right sound speed
    double a;                           // Max sound speed
    double x_0 = ui->x0->value();       // Left/Right boundary 
    double dx = 1.0/N;                  // Width of cells
    double t = 0.0;                     // Time
    double dt;                          // Time step
    double dt_A; 
    double dt_B; 
    double tStop = ui->Time->value();   // Simulation time
    double c = ui->CFL->value();        // CFL number
    double S_max_A;                       // Maximum wave speed
    double S_max_B;                       // Maximum wave speed

    double phi[N+2*NGC];
   
    char scheme[10]=""; // Name of scheme for fluxes
    char limitFunc[10]=""; // Name of limiter function 

    setScheme(scheme, limitFunc, 
            schemeChoice, limitChoice);
    
    // Vectors of primitive vars
    Primitive *W_A = new Primitive[N+2*NGC];
    Primitive *W_B = new Primitive[N+2*NGC];
    // Vectors of conserved vars
    Conserved *U_A = new Conserved[N+2*NGC];
    Conserved *U_A_old = new Conserved[N+2*NGC];
    Conserved *U_B = new Conserved[N+2*NGC];
    Conserved *U_B_old = new Conserved[N+2*NGC];

    QVector<double> Qx(N), Qrho(N), Qu(N), Qp(N), Qe(N);

    // Calculate sound speeds 
    // TODO: Check if update necessary
    // TODO: Check if vacuum generated
    a_L = sqrt(gamma_L*p_L/rho_L); 
    a_R = sqrt(gamma_R*p_R/rho_R);
    a = (fabs(a_L)>fabs(a_R))?fabs(a_L):fabs(a_R);

    // APPROXIMATE SOLVER:
    // Initiate solution vector
    initialConditions(W_A, U_A_old, W_B, U_B_old, phi, x_0, N, dx,
            rho_L, u_L, p_L, gamma_L, rho_R, u_R, p_R, gamma_R);

    // Loop from t=0 to tStop
    while(t < tStop){

        // Calculate ghost BCs for each material
        updateGhostCells(W_A, U_A_old, W_B, U_B_old, phi, 
                gamma_L, gamma_R);

        // Calculate maximum wave speed, update dt
        // TODO: Update for ghost fluid method
        S_max_A = maxWaveSpeed(W_A, a, N);
        S_max_B = maxWaveSpeed(W_B, a, N);
        dt_A = c*dx/S_max_A; 
        dt_B = c*dx/S_max_B; 
        dt = dt_A < dt_B ? dt_A : dt_B;

        // Make sure we don't go past tStop
        if(t+dt > tStop) dt = tStop-t;

        // Apply boundary conditionsn for domain
        boundaryConditions(W_A, U_A_old, gamma_L, N); 
        boundaryConditions(W_B, U_B_old, gamma_R, N); 

        // Advance level set function in time
        advanceLevelSet(phi, W_A, W_B, dt, dx, N);
        int pos = interfacePosition(phi);
        qDebug() << "Interface at: " << (pos-NGC+0.5)*dx; 
        reinitialize(phi, N, dx);
        
        // Advance system in time
        advance(W_A, U_A, U_A_old, dt, dx, N, gamma_L,
                scheme, limitFunc);
        advance(W_B, U_B, U_B_old, dt, dx, N, gamma_R,
                scheme, limitFunc);
        t += dt;
    }

    for(int i=NGC; i<N+NGC; i++){
        Qx[i-NGC] = (i-NGC+0.5)*dx;
        if(i<interfacePosition(phi)){
            Qrho[i-NGC] = W_A[i].rho;
            Qu[i-NGC] = W_A[i].u;
            Qp[i-NGC] = W_A[i].p;
            Qe[i-NGC] = W_A[i].p/(W_A[i].rho*(gamma_L-1));
        }else{
            Qrho[i-NGC] = W_B[i].rho;
            Qu[i-NGC] = W_B[i].u;
            Qp[i-NGC] = W_B[i].p;
            Qe[i-NGC] = W_B[i].p/(W_B[i].rho*(gamma_R-1));
        }
    }

    PlotGraph(ui->densityPlot, Qx, Qrho, "Density");
    PlotGraph(ui->velocityPlot, Qx, Qu, "Velocity");
    PlotGraph(ui->pressurePlot, Qx, Qp, "Pressure");
    PlotGraph(ui->energyPlot, Qx, Qe, "Internal energy");

    if(save){
        std::fstream fs; 
        QString filename = ui->fileName->text(); 
        filename.prepend("../Results/");
        fs.open(filename.toAscii(), std::fstream::out);
        fs << "x\trho\tu\tp\te\n"; 
        for(int i=0; i<N; i++){
            fs << Qx[i] << "\t" << Qrho[i] << "\t" << Qu[i] 
                << "\t" << Qp[i] << "\t" << Qe[i] << "\n";
        }
        fs.close();
        ui->saveResults->setChecked(0);
    }
}
