#include "main_window.h"
#include "ui_main_window.h"
#include "global.h"

MainWindow::MainWindow(QWidget *parent) :
  QMainWindow(parent),
  ui(new Ui::MainWindow)
{
  ui->setupUi(this);
  global::gl = ui->opengl_widget;
}

MainWindow::~MainWindow()
{
  delete ui;
}

