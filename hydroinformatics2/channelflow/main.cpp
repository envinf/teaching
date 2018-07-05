//Hydroinformatics-II
//Channel flow
//by Olaf Kolditz
#include <QApplication>
#include "dialog.h"

int main(int argc, char *argv[])
{
  QApplication a(argc, argv);
  Dialog w;
  w.setWindowTitle("Newton Simulator");
  //w.setFixedWidth(300);
  w.setFixedWidth(450);
  w.show();
  return a.exec();
}
