#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    void on_schemeCombo_currentIndexChanged();
    void on_GFMCombo_currentIndexChanged();
    void on_Interfaces_valueChanged(int i);
    void on_testCombo_currentIndexChanged();
    void on_Animate_stateChanged(int state);
    void on_saveResults_stateChanged(int state);
    void on_runButton_clicked();
    
private:
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
