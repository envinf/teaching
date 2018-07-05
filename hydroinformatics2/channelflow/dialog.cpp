#include "dialog.h"
#include "plotter.h"

#include <QVBoxLayout>
#include <QPushButton>
#include <QMessageBox>
#include <QLabel>

#include <string>
using namespace std;
#include <cmath>

Dialog::Dialog(QWidget *parent)
    : QDialog(parent)
{
  //elements
  //text labels
  QLabel* labelHeader = new QLabel(tr("Here we can set text ..."));
  QLabel* labelIC = new QLabel(tr("Initial water level:"));
  QLabel* labelBC = new QLabel(tr("Water level boundary:"));
  QLabel* labelST = new QLabel(tr("Channel discharge:"));
  QLabel* labelMAT1 = new QLabel(tr("Friction coefficient 1:"));
  QLabel* labelMAT2 = new QLabel(tr("Friction coefficient 2:"));
  QLabel* labelGEO1 = new QLabel(tr("Bed slope:"));
  QLabel* labelNUM1 = new QLabel(tr("Newton error tolerance:"));
  QLabel* labelNUM2 = new QLabel(tr("Newton error:"));
  QLabel* labelITE = new QLabel(tr("Newton iterations:"));
  //graphic labels
  QLabel *label_ogs = new QLabel();
  label_ogs->setAlignment(Qt::AlignCenter);
  label_ogs->setPixmap(QPixmap("../ogs_teaching_150.png"));
  QLabel *label_exercise = new QLabel();
  label_exercise->setAlignment(Qt::AlignCenter);
  label_exercise->setPixmap(QPixmap("../bhywi-08-08_title.png"));
  QLabel *label_channel = new QLabel();
  label_channel->setAlignment(Qt::AlignCenter);
  label_channel->setPixmap(QPixmap("../gerinne.png"));
  //push buttons
  pushButtonIC = new QPushButton(tr("Initial conditions"));
  pushButtonIC->setEnabled(false);
  pushButtonIC->setToolTip("This is to specify initial conditions for the PDE");
  pushButtonBC = new QPushButton(tr("Boundary conditions"));
  pushButtonBC->setEnabled(false);
  pushButtonMAT = new QPushButton(tr("Material conditions"));
  pushButtonMAT->setEnabled(false);
  pushButtonRUN = new QPushButton(tr("Run Newton step"));
  pushButtonRUN->setEnabled(true);
  pushButtonSHO = new QPushButton(tr("Show results"));
  pushButtonSHO->setEnabled(false);
  pushButtonALL = new QPushButton(tr("All-in-one"));
  pushButtonICChange = new QPushButton(tr("Change IC value"));
  pushButtonBCChange = new QPushButton(tr("Change BC value"));
  pushButtonDischargeChange = new QPushButton(tr("Change discharge value"));
  pushButtonFrictionCoefficientChange = new QPushButton(tr("Change friction value"));
  pushButtonBedSlopeChange = new QPushButton(tr("Change bed slope value"));
  // edits
  lineEditIC = new QLineEdit();
  lineEditBCR = new QLineEdit();
  lineEditFrictionLawExponent = new QLineEdit();
  lineEditFrictionCoefficient = new QLineEdit();
  lineEditBedSlope = new QLineEdit();
  lineEditDischarge = new QLineEdit();
  lineEditNewtonTolerance = new QLineEdit();
  lineEditNewtonError = new QLineEdit();
  lineEditIterations = new QLineEdit();
  //connect
  connect(pushButtonIC,SIGNAL(clicked()),this,SLOT(on_pushButtonIC_clicked()));
  connect(pushButtonBC,SIGNAL(clicked()),this,SLOT(on_pushButtonBC_clicked()));
  connect(pushButtonMAT,SIGNAL(clicked()),this,SLOT(on_pushButtonMAT_clicked()));
  connect(pushButtonRUN,SIGNAL(clicked()),this,SLOT(on_pushButtonRUN_clicked()));
  connect(pushButtonSHO,SIGNAL(clicked()),this,SLOT(on_pushButtonSHO_clicked()));
  connect(pushButtonALL,SIGNAL(clicked()),this,SLOT(on_pushButtonALL_clicked()));
  //connect(lineEditIC,SIGNAL(returnPressed()),this,SLOT(setText(sIC)));
  //connect(lineEditIC,SIGNAL(textChanged(text)),this,SLOT(Dummy()));
  //connect(lineEditIC,SIGNAL(textChanged(text)),this,SLOT(setText(sIC)));
  //connect(lineEditIC,SIGNAL(textEdited(text)),this,SLOT(setText(sIC)));
  //connect(lineEditIC,SIGNAL(textEdited(text)),this,SLOT(Dummy()));
  connect(pushButtonICChange,SIGNAL(clicked()),this,SLOT(on_pushButtonICChange_clicked()));
  connect(pushButtonBCChange,SIGNAL(clicked()),this,SLOT(on_pushButtonBCChange_clicked()));
  connect(pushButtonDischargeChange,SIGNAL(clicked()),this,SLOT(on_pushButtonDischargeChange_clicked()));
  connect(pushButtonFrictionCoefficientChange,SIGNAL(clicked()),this,SLOT(on_pushButtonFrictionCoefficientChange_clicked()));
  connect(pushButtonBedSlopeChange,SIGNAL(clicked()),this,SLOT(on_pushButtonBedSlopeChange_clicked()));
  //layout
  QVBoxLayout *leftLayout = new QVBoxLayout;
  leftLayout->addWidget(label_ogs);
  leftLayout->addWidget(labelHeader);
  leftLayout->addWidget(label_exercise);
  leftLayout->addWidget(pushButtonIC);
  leftLayout->addWidget(pushButtonBC);
  leftLayout->addWidget(pushButtonMAT);
  leftLayout->addWidget(pushButtonRUN);
  leftLayout->addWidget(pushButtonSHO);
  leftLayout->addWidget(pushButtonALL);
  leftLayout->addWidget(label_channel);
  QVBoxLayout *rightLayout = new QVBoxLayout;
  rightLayout->addWidget(labelIC);
  rightLayout->addWidget(lineEditIC);
  rightLayout->addWidget(labelBC);
  rightLayout->addWidget(lineEditBCR);
  rightLayout->addWidget(labelST);
  rightLayout->addWidget(lineEditDischarge);
  rightLayout->addWidget(labelMAT1);
  rightLayout->addWidget(lineEditFrictionLawExponent);
  rightLayout->addWidget(labelMAT2);
  rightLayout->addWidget(lineEditFrictionCoefficient);
  rightLayout->addWidget(labelGEO1);
  rightLayout->addWidget(lineEditBedSlope);
  rightLayout->addWidget(labelNUM1);
  rightLayout->addWidget(lineEditNewtonTolerance);
  rightLayout->addWidget(labelNUM2);
  rightLayout->addWidget(lineEditNewtonError);
  rightLayout->addWidget(labelITE);
  rightLayout->addWidget(lineEditIterations);
  QVBoxLayout *changeLayout = new QVBoxLayout;
  changeLayout->setAlignment(Qt::AlignTop);
  changeLayout->addWidget(pushButtonICChange);
  changeLayout->addWidget(pushButtonBCChange);
  changeLayout->addWidget(pushButtonDischargeChange);
  changeLayout->addWidget(pushButtonFrictionCoefficientChange);
  changeLayout->addWidget(pushButtonBedSlopeChange);
  QHBoxLayout *mainLayout = new QHBoxLayout;
  mainLayout->addLayout(leftLayout);
  mainLayout->addLayout(rightLayout);
  mainLayout->addLayout(changeLayout);
  setLayout(mainLayout);
  setWindowTitle(tr("E9"));
  //mainLayout->addStretch();
  //memory
  out_file.open("out.txt");
  //out_file.setf(ios::scientific);
  out_file.precision(3);
  n = 201;
  u_new.resize(n);
  u_old.resize(n);
  matrix = new double[n*n];
  vecb = new double[n];
  vecx = new double[n];
  wetted_cross_section = new double[n];
  water_level_elevation = new double[n];
  flow_velocity = new double[n];
  Froude_number = new double[n];
  wetted_perimeter = new double[n];
  hydraulic_radius = new double[n];
  friction_slope = new double[n];
  //geometry
  x = new double[n];
  double dx = 100./double(n-1);
  for(int i=0;i<n;i++)
//   x[i] = -100. + i*10.;
    x[i] = -100. + i*dx;
  bottom_elevation = new double[n];
  for(int i=0;i<n;i++)
    bottom_elevation[i] = 0.04 - i*0.004;
  //iteration
  kn = 12;
  k=0;
  sIC = "999";
  ICValue = 0.2;
  BCValue = 0.15;
  bed_slope = 0.0004; // [m/m]
  friction_coefficient = 10.; //
  discharge = 0.1; // Volumenflie�rate [m3/s]
  counter = 0;
  //plotter
  plotter = new Plotter;
  plotter->setWindowTitle(QObject::tr("My Function Plotter"));
  PlotSettings settings;
  settings.minX = -100.0;
  settings.maxX = 0.0;
  settings.minY = 0.0;
  settings.maxY = 1.0;
  plotter->setPlotSettings(settings);
  plotter2 = new Plotter;
  plotter2->setWindowTitle(QObject::tr("My Function Plotter 2"));
  plotter2->setPlotSettings(settings);
  //for bugfixing
  on_pushButtonIC_clicked();
  on_pushButtonBC_clicked();
  on_pushButtonMAT_clicked();
}

Dialog::~Dialog()
{
  delete [] matrix;
  delete [] vecb;
  delete [] vecx;
  delete [] wetted_cross_section;
  delete [] water_level_elevation;
  delete [] flow_velocity;
  delete [] Froude_number;
  delete [] wetted_perimeter;
  delete [] hydraulic_radius;
  delete [] friction_slope;
}

void Dialog::on_pushButtonIC_clicked()
{
  // Anfangsbedingungen setzen
  for(int i=0;i<n;i++)
  {
    u_old[i] = ICValue;
    u_new[i] = ICValue;
  }
  // Daten im Dialog sichtbar machen
  sDummy.setNum(ICValue,'f',5);
  lineEditIC->setText(sDummy);
  // Schnick-Schnack
  pushButtonIC->setStyleSheet("background-color: green");
}

void Dialog::on_pushButtonBC_clicked()
{
  u_old[n-1] = BCValue; // Wasserstand flu�abw�rts [m]
  u_new[n-1] = BCValue; // Wasserstand flu�abw�rts [m]
  // Daten im Dialog sichtbar machen
  sDummy.setNum(BCValue,'f',5);
  lineEditBCR->setText(sDummy);
  pushButtonBC->setStyleSheet("background-color: green");
}

void Dialog::on_pushButtonMAT_clicked()
{
  //ab in den Konstruktor! discharge = 0.05; // Volumenflie�rate [m3/s]
  gravity = 9.81; // [m/s2]
  friction_law_exponent = 0.5; // Chezy, Manning-Strickler [-]
  error_tolerance = 1e-2; // [m]
  bottom_width = 1.; // [m]
  m = 1.; //
  // Daten im Dialog sichtbar machen
  sDummy.setNum(friction_law_exponent,'f',5);
  lineEditFrictionLawExponent->setText(sDummy);
  sDummy.setNum(m,'f',5);
  lineEditFrictionCoefficient->setText(sDummy);
  //out_file << bed_slope << endl;
  sDummy.setNum(bed_slope,'f',5);
  lineEditBedSlope->setText(sDummy);
  sDummy.setNum(discharge,'f',5);
  lineEditDischarge->setText(sDummy);
  sDummy.setNum(error_tolerance,'f',5);
  lineEditNewtonTolerance->setText(sDummy);
  pushButtonMAT->setStyleSheet("background-color: green");
}

void Dialog::on_pushButtonRUN_clicked()
{
  double error;
  error = RUN_NewtonStep();
  QVector<QPointF> points0;
  for(int i=0;i<n-1;i++)
  {
    points0.append(QPointF(x[i],u_new[i]));
  }
  //  points0.append(QPointF(x[10],u_new[10]));
  points0.append(QPointF(x[n-1],u_new[n-1]));
  plotter->setCurveData(k++, points0);
  plotter->show();
  QString sError;
  sError.setNum(error,'f',5);
  lineEditNewtonError->setText(sError);
  counter++;
  sError.setNum(counter,5);
  lineEditIterations->setText(sError);
}

void Dialog::on_pushButtonALL_clicked()
{
  on_pushButtonIC_clicked();
  on_pushButtonBC_clicked();
  on_pushButtonMAT_clicked();
  //on_pushButtonRUN_clicked();
  //
  float error = 1.1*error_tolerance;
  while(error>error_tolerance)
  {
    error = RUN_NewtonStep();
    QVector<QPointF> points0;
    for(int i=0;i<n-1;i++)
    {
      points0.append(QPointF(x[i],u_new[i]));
    }
    points0.append(QPointF(x[n-1],u_new[n-1]));
    plotter->setCurveData(k++, points0);
    plotter2->setCurveData(k, points0);
    //
    sDummy.setNum(error,'f',5);
    lineEditNewtonError->setText(sDummy);
    sDummy.setNum(k,5);
    lineEditIterations->setText(sDummy);
  }
  plotter->show();
  PlotSettings settings;
  settings.minX = -100.0;
  settings.maxX = 0.0;
  settings.minY = 0.0;
  settings.maxY = 1.0;
  plotter2 = new Plotter;
  plotter2->setWindowTitle(QObject::tr("My Function Plotter 2"));
  plotter2->setPlotSettings(settings);
  plotter2->show();
}

double Dialog::RUN_NewtonStep()
{
  //local variables
  double N,N1,N2,N3,D,D1,D2,D21,D22;
  double error = 0;
  {
    //start values
    for(int i=0;i<n;i++)
    {
      wetted_perimeter[i] = bottom_width + 2.*sqrt(1.+m*m)*u_old[i];
      wetted_cross_section[i] = (bottom_width + m*u_old[i])*u_old[i];
      hydraulic_radius[i] = wetted_cross_section[i] / wetted_perimeter[i];
      water_level_elevation[i] = bottom_elevation[i] + u_old[i];
      flow_velocity[i] = discharge/wetted_cross_section[i];
      Froude_number[i] = flow_velocity[i]/(sqrt(gravity*wetted_cross_section[i]\
                       /sqrt(bottom_width*bottom_width+4.*m*wetted_cross_section[i])));
      friction_slope[i] = pow(flow_velocity[i]/(friction_coefficient*pow(hydraulic_radius[i],friction_law_exponent)),2);
    }
    //Newton step
    for(int i=0;i<n-1;i++)
    {
      N1 = pow(discharge,2)/pow(wetted_cross_section[i+1],2) + gravity*u_old[i+1];
      N2 = pow(discharge,2)/pow(wetted_cross_section[i],2) + gravity*u_old[i];
      N3 = gravity*(bed_slope - (friction_slope[i+1]+friction_slope[i])/2.)*(x[i+1]-x[i]);
      N = N1 - N2 - N3;
      D1 = pow(discharge,2)/pow(wetted_cross_section[i],3) * (bottom_width+2.*m*u_old[i]) - gravity;
      D21 = friction_law_exponent*2.*(sqrt(1+m*m))/wetted_perimeter[i];
      D22 = (1.+friction_law_exponent)/wetted_cross_section[i] * (bottom_width+2.*m*u_old[i]);
      D2 = gravity*friction_slope[i]*(D21-D22)*(x[i+1]-x[i]);
      D = D1 + D2;
      u_new[i] = u_old[i] - N/D;
    }
    //calc Newton error
    for(int i=0;i<n-1;i++)
    {
      error += u_old[i] -u_new[i];
    }
//E3
    error = sqrt(error*error);
    //error = abs(error);
    //save Newton step
    for(int i=0;i<n-1;i++)
    {
      u_old[i] = u_new[i];
    }
  }
  return error;
}

void Dialog::on_pushButtonICChange_clicked()
{
  QString sICValue = lineEditIC->text();
  ICValue = lineEditIC->text().toDouble();
  on_pushButtonIC_clicked();
}

void Dialog::on_pushButtonBCChange_clicked()
{
  QString sBCValue = lineEditBCR->text();
  BCValue = lineEditBCR->text().toDouble();
  on_pushButtonBC_clicked();
}

void Dialog::on_pushButtonDischargeChange_clicked()
{
  discharge = lineEditDischarge->text().toDouble();
  on_pushButtonMAT_clicked();
}

void Dialog::on_pushButtonFrictionLawExponentChange_clicked()
{
}

void Dialog::on_pushButtonFrictionCoefficientChange_clicked()
{
  friction_coefficient = lineEditFrictionCoefficient->text().toDouble();
  on_pushButtonMAT_clicked();
}

void Dialog::on_pushButtonBedSlopeChange_clicked()
{
  bed_slope = lineEditBedSlope->text().toDouble();
  on_pushButtonMAT_clicked();
}

void Dialog::on_pushButtonNewtonToleranceChange_clicked()
{
}
