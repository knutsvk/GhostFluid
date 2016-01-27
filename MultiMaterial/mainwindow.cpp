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
    ui->fileName->setVisible(0);
    ui->schemeCombo->setCurrentIndex(2);
    ui->limiterCombo->setCurrentIndex(2);
    ui->isobaricFix->setVisible(0);
    ui->x1->setVisible(0);
    ui->x2->setVisible(0);
    ui->rho2->setVisible(0);
    ui->rho3->setVisible(0);
    ui->u2->setVisible(0);
    ui->u3->setVisible(0);
    ui->p2->setVisible(0);
    ui->p3->setVisible(0);
    ui->mat1->setCurrentIndex(1);
    ui->mat2->setVisible(0);
    ui->mat3->setVisible(0);
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

void MainWindow::on_GFMCombo_currentIndexChanged()
{
    int n = ui->GFMCombo->currentIndex();
    ui->isobaricFix->setVisible(n>0);
}

void MainWindow::on_Interfaces_valueChanged(int i)
{
    ui->x1->setVisible(i>1);
    ui->rho2->setVisible(i>1);
    ui->u2->setVisible(i>1);
    ui->p2->setVisible(i>1);
    ui->mat2->setVisible(i>1);
    ui->x2->setVisible(i>2);
    ui->rho3->setVisible(i>2);
    ui->u3->setVisible(i>2);
    ui->p3->setVisible(i>2);
    ui->mat3->setVisible(i>2);
}

void MainWindow::on_testCombo_currentIndexChanged()
{
    int n = ui->testCombo->currentIndex();
    int nInterfaces; 
    double x[3];          
    double rho[4];       // Left density
    double u[4];         // Left velocity
    double p[4];         // Left pressure
    int mat[4]; 
    double gamma[2];     // Left gamma
    double tStop;       // Final time

    initiateTestCase(n+1, nInterfaces, x, rho, u, p, mat, gamma, 
            tStop);

    ui->Interfaces->setValue(nInterfaces);
    
    ui->x0->setValue(x[0]);
    ui->x1->setValue(x[1]);
    ui->x2->setValue(x[2]);

    ui->rho0->setValue(rho[0]);
    ui->rho1->setValue(rho[1]);
    ui->rho2->setValue(rho[2]);
    ui->rho3->setValue(rho[3]);

    ui->u0->setValue(u[0]);
    ui->u1->setValue(u[1]);
    ui->u2->setValue(u[2]);
    ui->u3->setValue(u[3]);
    
    ui->p0->setValue(p[0]);
    ui->p1->setValue(p[1]);
    ui->p2->setValue(p[2]);
    ui->p3->setValue(p[3]);

    ui->mat0->setCurrentIndex(mat[0]);
    ui->mat1->setCurrentIndex(mat[1]);
    ui->mat2->setCurrentIndex(mat[2]);
    ui->mat3->setCurrentIndex(mat[3]);

    ui->gammaA->setValue(gamma[0]);
    ui->gammaB->setValue(gamma[1]);

    ui->Time->setValue(tStop);
    ui->CFL->setValue(0.9);
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
    int method = ui->GFMCombo->currentIndex();
    int nInterfaces = ui->Interfaces->value();
    int nRegions = nInterfaces+1;
    int nMaterialInterfaces=0;

    bool save = ui->saveResults->isChecked(); 
    bool isobarix = ui->isobaricFix->isChecked();
    if(method==1 && isobarix) method++;

    double x[nInterfaces];
    double rho[nRegions];
    double u[nRegions];
    double p[nRegions];
    double mat[nRegions];
    double gamma[2];
    gamma[0] = ui->gammaA->value();
    gamma[1] = ui->gammaB->value();
    double a = 0.0;   // Maximum sound speed.
    double dx = 1.0/N;          // Width of cells
    double t = 0.0;               // Time
    double dt;      // Time step. 
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

    QVector<double> Qx(N), Qrho(N), Qu(N), Qp(N), Qe(N), Qphi(N);

    // Read info from ui 
    x[0] = ui->x0->value();
    rho[0] = ui->rho0->value();
    u[0] = ui->u0->value();
    p[0] = ui->p0->value();
    mat[0] = ui->mat0->currentIndex();
    rho[1] = ui->rho1->value();
    u[1] = ui->u1->value();
    p[1] = ui->p1->value();
    mat[1] = ui->mat1->currentIndex();
    if(nInterfaces>1){
        x[1] = ui->x1->value();
        rho[2] = ui->rho2->value();
        u[2] = ui->u2->value();
        p[2] = ui->p2->value();
        mat[2] = ui->mat2->currentIndex();
        if(nInterfaces>2){
            x[2] = ui->x2->value();
            rho[3] = ui->rho3->value();
            u[3] = ui->u3->value();
            p[3] = ui->p3->value();
            mat[3] = ui->mat3->currentIndex();
        }
    }

    // Find which interfaces are material interfaces
    for(int i=0; i<nInterfaces; i++)
        if(mat[i]!=mat[i+1]) nMaterialInterfaces++;


    double x_M[nMaterialInterfaces];
    int count = 0;
    for(int i=0; i<nInterfaces; i++){
        if(mat[i]!=mat[i+1]){
            x_M[count]=x[i];
            count++;
        }
    }

    // Calculate sound speeds 
    for(int i=0; i<nRegions; i++){
        int mati = mat[i];
        double a_i = sqrt(gamma[mati]*p[i]/rho[i]);
        if(fabs(a_i)>a) a = fabs(a_i);
    }

    // APPROXIMATE SOLVER:
    // Initiate solution vector
    initialConditions(W_A, U_A_old, W_B, U_B_old, phi, N, dx,
            nInterfaces, x, nMaterialInterfaces, x_M, 
            rho, u, p, mat, gamma);

    // Loop from t=0 to tStop
    while(t < tStop){

        // Calculate ghost BCs for each material interface
        int *pos = interfacePosition(phi, nMaterialInterfaces);
        for(int i=0; i<nMaterialInterfaces; i++){
            if(phi[pos[i]-1]<0){
                updateGhostCells(W_A, U_A_old, W_B, U_B_old, 
                        pos[i], gamma[0], gamma[1], method);
            }else{
                updateGhostCells(W_B, U_B_old, W_A, U_A_old, 
                        pos[i], gamma[1], gamma[0], method);
            }
        }
        delete[] pos;

        // Calculate maximum wave speed, update dt
        S_max_A = maxWaveSpeed(W_A, a, N);
        S_max_B = maxWaveSpeed(W_B, a, N);
        dt = c*dx/S_max_A; 
        if(c*dx/S_max_B < dt) dt = c*dx/S_max_B; 

        // Make sure we don't go past tStop
        if(t+dt > tStop) dt = tStop-t;

        // Apply boundary conditionsn for domain
        boundaryConditions(W_A, U_A_old, gamma[0], N); 
        boundaryConditions(W_B, U_B_old, gamma[1], N); 

        // Advance level set function in time
        advanceLevelSet(phi, W_A, W_B, dt, dx, N);
        reinitialize(phi, N, dx, nMaterialInterfaces);
        
        // Advance system in time
        advance(W_A, U_A, U_A_old, dt, dx, N, gamma[0],
                scheme, limitFunc);
        advance(W_B, U_B, U_B_old, dt, dx, N, gamma[1],
                scheme, limitFunc);
        t += dt;
    }

    for(int i=NGC; i<N+NGC; i++){
        Qx[i-NGC] = (i-NGC+0.5)*dx;
        Qphi[i-NGC] = phi[i];
        if(phi[i]<0.0){
            Qrho[i-NGC] = W_A[i].rho;
            Qu[i-NGC] = W_A[i].u;
            Qp[i-NGC] = W_A[i].p;
            Qe[i-NGC] = W_A[i].p/(W_A[i].rho*(gamma[0]-1));
        }else{
            Qrho[i-NGC] = W_B[i].rho;
            Qu[i-NGC] = W_B[i].u;
            Qp[i-NGC] = W_B[i].p;
            Qe[i-NGC] = W_B[i].p/(W_B[i].rho*(gamma[1]-1));
        }
    }

    // TODO: Include exact when possible
    PlotGraph(ui->densityPlot, Qx, Qrho, "Density");
    PlotGraph(ui->velocityPlot, Qx, Qu, "Velocity");
    PlotGraph(ui->pressurePlot, Qx, Qp, "Pressure");
    PlotGraph(ui->energyPlot, Qx, Qe, "Internal energy");
    PlotGraph(ui->levelsetPlot, Qx, Qphi, "Level set function");

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
