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
    ui->RPCombo->setVisible(0);
    ui->fileName->setVisible(0);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_schemeCombo_currentIndexChanged()
{
    int n = ui->schemeCombo->currentIndex();
    ui->limiterCombo->setVisible(n==1 || n==2 || n==4);
    ui->RPCombo->setVisible(n==2 || n==3);
}

void MainWindow::on_testCombo_currentIndexChanged()
{
    int n = ui->testCombo->currentIndex();
    double rho_L,       // Left density
           u_L,         // Left velocity
           p_L,         // Left pressure
           rho_R,       // Right density
           u_R,         // Right velocity
           p_R,         // Right pressure
           tStop,       // Final time
           x_0;          
    initiateTestCase(n+1, tStop, x_0, rho_L, u_L, p_L, 
            rho_R, u_R, p_R);
    ui->rhoL->setValue(rho_L);
    ui->uL->setValue(u_L);
    ui->pL->setValue(p_L);
    ui->rhoR->setValue(rho_R);
    ui->uR->setValue(u_R);
    ui->pR->setValue(p_R);
    ui->Time->setValue(tStop);
    ui->x0->setValue(x_0);
    ui->CFL->setValue(0.9);
    ui->Cells->setValue(100);
}

void MainWindow::on_saveResults_stateChanged()
{
    bool save = ui->saveResults->isTristate();
    ui->fileName->setVisible(save);
}

void PlotGraph(QCustomPlot* plot, QVector<double> X, 
        QVector<double> Y, QVector<double> Z, QString label)
{
    plot->addGraph();
    plot->graph(0)->setPen(QPen(Qt::blue));
    plot->graph(0)->setData(X,Z);
    double min0 = *std::min_element(Z.begin(),Z.end());
    double max0 = *std::max_element(Z.begin(),Z.end());

    plot->addGraph();
    plot->graph(1)->setScatterStyle(QCPScatterStyle(
                QCPScatterStyle::ssCircle, QPen(Qt::red,0.5), 
                QBrush(Qt::red), 3));
    plot->graph(1)->setLineStyle(QCPGraph::lsNone);
    plot->graph(1)->setData(X,Y);
    double min1 = *std::min_element(Y.begin(),Y.end());
    double max1 = *std::max_element(Y.begin(),Y.end());

    double MIN = std::min(min0,min1);
    double MAX = std::max(max0,max1);
    double ymin = MIN - 0.05*(MAX-MIN);
    double ymax = MAX + 0.05*(MAX-MIN);
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
    int N = ui->Cells->value(),
        schemeChoice = ui->schemeCombo->currentIndex(),
        limitChoice = ui->limiterCombo->currentIndex(),
        RPChoice = ui->RPCombo->currentIndex();

    bool save = ui->saveResults->isTristate(); 

    double rho_L = ui->rhoL->value(),   // Left density
           u_L = ui->uL->value(),       // Left velocity
           p_L = ui->pL->value(),       // Left pressure
           rho_R = ui->rhoR->value(),   // Right density
           u_R = ui->uR->value(),       // Right velocity
           p_R = ui->pR->value(),       // Right pressure
           a_L,                         // Left sound speed
           a_R,                         // Right sound speed
           a,                           // Max sound speed
           x_0 = ui->x0->value(),       // Left/Right boundary 
           dx = 1.0/N,                  // Width of cells
           t = 0.0,                     // Time
           dt,                          // Time step
           tStop = ui->Time->value(),   // Simulation time
           c = ui->CFL->value(),        // CFL number
           S_max,                       // Maximum wave speed
           // For exact solver only: 
           S,                           // Speed (=x/t)
           p_S,                         // Pressure in star region
           u_S,                         // Velocity in star region
           rho,                         // Density 
           u,                           // Velocity 
           p,                           // Pressure 
           e;                           // Internal energy 
   
    char scheme[10]="", // Name of scheme for fluxes
         limitFunc[10]="", // Name of limiter function 
         RP[10]="";        // Name of Riemann solver

    if(schemeChoice == 0)
        strcpy(scheme, "FORCE");
    else if(schemeChoice == 1){
        strcpy(scheme, "FLIC");
        if(limitChoice == 0)
            strcpy(limitFunc, "minbee");
        else if(limitChoice == 1)
            strcpy(limitFunc, "vanleer");
        else
            strcpy(limitFunc, "superbee");
    }
    else if(schemeChoice == 2){
        strcpy(scheme, "MUSCL");
        if(limitChoice == 0)
            strcpy(limitFunc, "minbee");
        else if(limitChoice == 1)
            strcpy(limitFunc, "vanleer");
        else
            strcpy(limitFunc, "superbee");
        if(RPChoice == 0)
            strcpy(RP, "HLLC");
        else
            strcpy(RP, "Exact");
    }
    else if(schemeChoice == 3){
        strcpy(scheme, "Godunov");
        if(RPChoice == 0)
            strcpy(RP, "HLLC");
        else
            strcpy(RP, "Exact");
    }
    else{
        strcpy(scheme, "SLIC"); 
        if(limitChoice == 0)
            strcpy(limitFunc, "minbee");
        else if(limitChoice == 1)
            strcpy(limitFunc, "vanleer");
        else
            strcpy(limitFunc, "superbee");
    }
    
    // Vectors of primitive vars
    Primitive *W = new Primitive[N+2*NGC];
    // Vectors of conserved vars
    Conserved *U = new Conserved[N+2*NGC],
              *U_old = new Conserved[N+2*NGC];

    QVector<double> Qx(N), Qrho(N), Qu(N), Qp(N), Qe(N),
        QrhoExact(N), QuExact(N), QpExact(N), QeExact(N);

    // Calculate sound speeds to check for vacuum
    a_L = sqrt(g*p_L/rho_L); 
    a_R = sqrt(g*p_R/rho_R);
    if(2*(a_L+a_R)/(g-1) <= u_R-u_L){
       printf("Vacuum generated by data.\n");
       printf("Program stopped.\n");
    }

    // EXACT SOLVER: 
    // Compute pressure and velocity in star region
    starRegionPressureVelocity(p_S, u_S, 
            rho_L, u_L, p_L, a_L, 
            rho_R, u_R, p_R, a_R);

    // Calculate solution for each cell
    for(int i=0; i<N; i++){
        S = ((i+0.5)*dx-x_0)/tStop;
        sample(rho, u, p, p_S, u_S, S, 
                rho_L, u_L, p_L, a_L, 
                rho_R, u_R, p_R, a_R);
        e = p/(rho*(g-1));
        QrhoExact[i] = rho; 
        QuExact[i] = u; 
        QpExact[i] = p; 
        QeExact[i] = e; 
    }

    // APPROXIMATE SOLVER:
    // Initiate solution vector
    initialConditions(W, U_old, x_0, N, dx, 
            rho_L, u_L, p_L, rho_R, u_R, p_R);

    // Loop from t=0 to tStop
    while(t < tStop){

        // Calculate maximum wave speed, update dt
        a = (fabs(a_L)>fabs(a_R))?fabs(a_L):fabs(a_R);
        S_max = maxWaveSpeed(W, a, N);
        dt = c*dx/S_max; 

        // Make sure we don't go past tStop
        if(t+dt > tStop) dt = tStop-t;

        // Advance system in time
        boundaryConditions(W, U_old, N); 
        advance(W, U, U_old, dt, dx, N, 
                scheme, limitFunc, RP);
        t += dt;
    }

    for(int i=NGC; i<N+NGC; i++){
        Qx[i-NGC] = (i-NGC+0.5)*dx;
        Qrho[i-NGC] = W[i].rho;
        Qu[i-NGC] = W[i].u;
        Qp[i-NGC] = W[i].p;
        Qe[i-NGC] = W[i].p/(W[i].rho*(g-1));
    }

    PlotGraph(ui->densityPlot, Qx, Qrho, QrhoExact, "Density");
    PlotGraph(ui->velocityPlot, Qx, Qu, QuExact, "Velocity");
    PlotGraph(ui->pressurePlot, Qx, Qp, QpExact, "Pressure");
    PlotGraph(ui->energyPlot, Qx, Qe, QeExact, "Internal energy");

    if(save){
        std::fstream fs; 
        QString filename = ui->fileName->text(); 
        filename.prepend("../Results/");
        fs.open(filename.toAscii(), std::fstream::out);
        fs << "x\trho\tu\tp\te\trho_E\tu_E\tp_E\te_E\n"; 
        for(int i=0; i<N; i++){
            fs << Qx[i] << "\t" << Qrho[i] << "\t" << Qu[i] 
                << "\t" << Qp[i] << "\t" << Qe[i] << "\t" 
                << QrhoExact[i] << "\t" << QuExact[i] << "\t" 
                << QpExact[i] << "\t" << QeExact[i] << "\n";
        }
        fs.close();
    }
}
