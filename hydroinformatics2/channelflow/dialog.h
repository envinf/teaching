#ifndef DIALOG_H
#define DIALOG_H

#include <QDialog>
#include <vector>
#include <fstream>
using namespace std;
#include <QLineEdit>
#include "plotter.h"

class Dialog : public QDialog
{
    Q_OBJECT

public:
    Dialog(QWidget *parent = 0);
    ~Dialog();

private:
  double dx;
  vector<double> u_new, u_old;
  int n;
  ofstream out_file;
  double Ne;
  double dt;
  QPushButton* pushButtonIC;
  QPushButton* pushButtonBC;
  QPushButton* pushButtonMAT;
  QPushButton* pushButtonRUN;
  QPushButton* pushButtonSHO;
  double* matrix;
  double* vecb;
  double* vecx;
  double* wetted_cross_section;
  double* water_level_elevation;
  double* flow_velocity;
  double* Froude_number;
  double* wetted_perimeter;
  double* hydraulic_radius;
  double* friction_slope;
  double discharge;
  double gravity;
  double friction_law_exponent;
  double error_tolerance;
  double bed_slope;
  double bottom_width;
  double m;
  double friction_coefficient;
  double* x;
  double* bottom_elevation;
  int kn;
  QString sIC;
  QPushButton* pushButtonALL;
  Plotter *plotter;
  Plotter *plotter2;
  int k;
  QLineEdit* lineEditIC;
  QLineEdit* lineEditBCR;
  QLineEdit* lineEditDischarge;
  QLineEdit* lineEditFrictionLawExponent;
  QLineEdit* lineEditNewtonError;
  QLineEdit* lineEditNewtonTolerance;
  QLineEdit* lineEditFrictionCoefficient;
  QLineEdit* lineEditBedSlope;
  QLineEdit* lineEditIterations;
  QString sDummy;
  const QString text;
  QPushButton* pushButtonICChange;
  double ICValue;
  QPushButton* pushButtonBCChange;
  double BCValue;
  QPushButton* pushButtonDischargeChange;
  QPushButton* pushButtonFrictionCoefficientChange;
  QPushButton* pushButtonBedSlopeChange;
  int counter;

private slots:
    void on_pushButtonIC_clicked();
    void on_pushButtonBC_clicked();
    void on_pushButtonMAT_clicked();
    void on_pushButtonRUN_clicked();
    void on_pushButtonALL_clicked();
    double RUN_NewtonStep();
    void on_pushButtonICChange_clicked();
    void on_pushButtonBCChange_clicked();
    void on_pushButtonDischargeChange_clicked();
    void on_pushButtonFrictionLawExponentChange_clicked();
    void on_pushButtonFrictionCoefficientChange_clicked();
    void on_pushButtonBedSlopeChange_clicked();
    void on_pushButtonNewtonToleranceChange_clicked();
};

#endif // DIALOG_H
