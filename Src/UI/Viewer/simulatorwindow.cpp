// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "simulatorwindow.h"

#include <QAction>
#include <QMenuBar>
#include <QToolBar>
#include <QLabel>
#include <QStatusBar>
#include <QMessageBox>
#include <QApplication>
#include <QHeaderView>
#include <QVBoxLayout>
#include <QPushButton>
#include <QDialogButtonBox>
#include <QTimer>
#include <QComboBox>

#include <QFileDialog>
#include <QSettings>
#include <QImage>
#include <QWheelEvent>
#include <QElapsedTimer>

#include <QBuffer>

#include "version.h"
#include "DB/db.h"

#include "Common/mathalgo.h"
#include "Common/fileutils.h"

#include "UI/GUICommon/Simulator/imageview.h"
#include "UI/GUICommon/Simulator/concentplotwidget.h"
#include "UI/GUICommon/Simulator/simulatorthread.h"

#include "Search/generationindividual.h"
#include "Search/individual.h"
#include "Search/generation.h"
#include "Search/deme.h"
#include "Search/search.h"
#include "Search/searchexperiment.h"
#include "Search/evaluatormorphogens.h"

#include "Model/model.h"

#include "Simulator/simulator.h"
#include "Simulator/simstate.h"
#include "Simulator/simulatorconfig.h"
#include "Simulator/simparams.h"

#include "UI/GUICommon/experimentform.h"
#include "UI/GUICommon/modelform.h"
#include "UI/GUICommon/imagesaver.h"
#include "UI/GUICommon/modelgraphview.h"
#include "Experiment/experiment.h"


namespace LoboLab {

SimulatorWindow::SimulatorWindow(Individual *ind, Search *sea, DB *d, bool autod,
  QWidget *parent)
  : QDialog(parent),
  isAutoDelete_(autod),
  db_(d),
  individual_(ind),
  experiment_(sea->experiment(0)),
  model_(ind->model()),
  search_(sea),
  simParams_(sea->simParams()),
  modelForm_(NULL),
  simulator_(NULL),
  simThread_(NULL),
  comparisonState_(NULL),
  simulating_(false),
  testing_(false),
  internalBlockProdChange_(false),
  movieSaver_(this) {
  setWindowTitle(tr("Simulator"));

  createWidgets();
  readSettings();

  if (isAutoDelete_)
    setAttribute(Qt::WA_DeleteOnClose);

  QTimer::singleShot(0, this, SLOT(loadModel()));
}

SimulatorWindow::~SimulatorWindow() {
  writeSettings();
  delete simThread_;
  delete comparisonState_;
  delete modelForm_;
  delete simulator_;

  if (isAutoDelete_)
    delete individual_;
}

void SimulatorWindow::createWidgets() {
  QPushButton *editExperimentButton = new QPushButton("Exp:");
  editExperimentButton->setMinimumWidth(30);
  editExperimentButton->setMaximumWidth(40);
  connect(editExperimentButton, SIGNAL(clicked()),
    this, SLOT(editExperiment()));

  experimentComboBox_ = new QComboBox();
  int nExp = search_->nExperiments();
  for (int i = 0; i < nExp; ++i)
    experimentComboBox_->addItem(search_->experiment(i)->name());
  experimentComboBox_->setMinimumWidth(50);
  connect(experimentComboBox_, SIGNAL(currentIndexChanged(int)),
    this, SLOT(experimentComboChanged(int)));

  QPushButton *editModelButton = new QPushButton(QString("Ind %1").arg(individual_->id()));
  editModelButton->setMinimumWidth(40);
  connect(editModelButton, SIGNAL(clicked()),
    this, SLOT(editModel()));

  QPushButton *fitButton = new QPushButton("Error");
  fitButton->setMinimumWidth(40);
  connect(fitButton, SIGNAL(clicked()), this, SLOT(calcError()));

  QPushButton *testButton = new QPushButton("Test");
  testButton->setMinimumWidth(40);
  connect(testButton, SIGNAL(clicked()), this, SLOT(testModel()));

  QPushButton *loadButton = new QPushButton("Load");
  loadButton->setMinimumWidth(40);
  connect(loadButton, SIGNAL(clicked()), this, SLOT(loadModel()));

  startStopButton_ = new QPushButton("Start");
  startStopButton_->setMinimumWidth(40);
  connect(startStopButton_, SIGNAL(clicked()), this, SLOT(startStopSim()));

  QPushButton *predButton = new QPushButton("Pred");
  predButton->setMinimumWidth(40);
  connect(predButton, SIGNAL(clicked()), this, SLOT(saveAllPredictions()));

  QPushButton *exportButton = new QPushButton("Export");
  predButton->setMinimumWidth(40);
  connect(exportButton, SIGNAL(clicked()), this, SLOT(exportExperimentData()));

  movieButton_ = new QPushButton("Mov");
  movieButton_->setCheckable(true);
  movieButton_->setMinimumWidth(40);
  connect(movieButton_, SIGNAL(toggled(bool)), this, SLOT(startStopMovie(bool)));

  ModelGraphView *modelGraphView = new ModelGraphView(model_,
    search_->simParams()->nInputMorphogens + search_->simParams()->nTargetMorphogens, true);
  modelGraphView->setStyleSheet("border: 1px solid #C0C0C0");

  cellsImg_ = new ImageView();
  connect(cellsImg_, SIGNAL(rectDrawed(const QRect &, Qt::MouseButton)),
    this, SLOT(cellsRectDrawed(const QRect &, Qt::MouseButton)));

  concentImg_ = new ImageView(true);
  concentImg_->setRectColor(Qt::white);
  connect(concentImg_, SIGNAL(rectDrawed(const QRect &, Qt::MouseButton)),
    this, SLOT(concentRectDrawed(const QRect &, Qt::MouseButton)));
  connect(concentImg_, SIGNAL(mouseWheel(QWheelEvent *)),
    this, SLOT(concentWheeled(QWheelEvent *)));
  connect(concentImg_, SIGNAL(lineHmoved()),
    this, SLOT(concentLineHMoved()));
  connect(concentImg_, SIGNAL(lineVmoved()),
    this, SLOT(concentLineVMoved()));

  plotsWidget_ = new ConcentPlotWidget();

  statusText_ = new QLabel("Individual view");
  statusText_->setMinimumWidth(100);
  speedText_ = new QLabel("Sel. speed: ");
  speedSlider_ = new QSlider(Qt::Horizontal);
  speedSlider_->setMinimum(1);
  speedSlider_->setMaximum(maxSpeedSlider());
  speedSlider_->setSingleStep(1);
  speedSlider_->setPageStep(10);
  connect(speedSlider_, SIGNAL(valueChanged(int)),
    this, SLOT(speedSliderChanged(int)));

  QLabel *modelLabel = new QLabel("Model");
  modelLabel->setAlignment(Qt::AlignCenter);
  QLabel *cellsLabel = new QLabel("Morphogens concentrations");
  cellsLabel->setAlignment(Qt::AlignCenter);
  QLabel *concentLabel = new QLabel("All products concentrations");
  concentLabel->setAlignment(Qt::AlignCenter);

  QHBoxLayout *imgLabelsLay = new QHBoxLayout();
  imgLabelsLay->addWidget(modelLabel);
  imgLabelsLay->addWidget(cellsLabel);
  imgLabelsLay->addWidget(concentLabel);

  QHBoxLayout *simLay = new QHBoxLayout();
  simLay->addWidget(modelGraphView, 1);
  simLay->addWidget(cellsImg_, 1);
  simLay->addWidget(concentImg_, 1);

  QHBoxLayout *buttonsLay = new QHBoxLayout();
  buttonsLay->addWidget(editExperimentButton);
  buttonsLay->addWidget(experimentComboBox_);
  buttonsLay->addWidget(editModelButton);
  buttonsLay->addWidget(fitButton);
  buttonsLay->addWidget(testButton);
  buttonsLay->addWidget(loadButton);
  buttonsLay->addWidget(startStopButton_);
  buttonsLay->addWidget(predButton);
  buttonsLay->addWidget(exportButton);
  buttonsLay->addWidget(movieButton_);

  QHBoxLayout *bottomLay = new QHBoxLayout();
  bottomLay->addWidget(statusText_);
  bottomLay->addSpacing(15);
  bottomLay->addWidget(speedText_);
  bottomLay->addWidget(speedSlider_, 1);

  QVBoxLayout *graphsLay = new QVBoxLayout();
  graphsLay->addLayout(imgLabelsLay);
  graphsLay->addLayout(simLay);
  graphsLay->addWidget(plotsWidget_, 1);

  graphsWidget_ = new QWidget();
  graphsWidget_->setStyleSheet("background-color:white;");
  graphsWidget_->setLayout(graphsLay);

  QVBoxLayout *mainLay = new QVBoxLayout();
  mainLay->addLayout(buttonsLay);
  mainLay->addWidget(graphsWidget_);
  mainLay->addLayout(bottomLay);

  setLayout(mainLay);
}

void SimulatorWindow::readSettings() {
  QSettings settings;
  settings.beginGroup("SimulatorWindow");
  resize(1102, 844);
  move(settings.value("pos", pos()).toPoint());
  if (settings.value("maximized", false).toBool())
    showMaximized();

  speedSlider_->setValue(settings.value("speed", maxSpeedSlider()).toInt());
  speedSliderChanged(speedSlider_->value());
}

void SimulatorWindow::writeSettings() {
  QSettings settings;
  settings.beginGroup("SimulatorWindow");
  if (isMaximized())
    settings.setValue("maximized", isMaximized());
  else {
    settings.setValue("maximized", false);
    settings.setValue("size", size());
    settings.setValue("pos", pos());
  }

  settings.setValue("speed", speedSlider_->value());
}

void SimulatorWindow::loadModel() {
  double concentImageFactor = 1;
  if (simulator_) {
    if (simulating_)
      startStopSim();

    concentImageFactor = simThread_->getConcentImageFactor();
    delete simThread_;
    delete simulator_;
    delete comparisonState_;
  }

  simulator_ = new Simulator(*search_);
  simulator_->loadModel(model_);
  simulator_->loadExperiment(experiment_);
  simulator_->initialize();

  comparisonState_ = new SimState(simParams_->size,
    simParams_->nInputMorphogens, simParams_->nTargetMorphogens, simParams_->distErrorThreshold);
  comparisonState_->loadMorphologyImage(*experiment_->outputMorphology(), simParams_->nTargetMorphogens);

  plotsWidget_->setSimulator(simulator_);
  concentLineHMoved();
  concentLineVMoved();
  plotsWidget_->updatePlots();

  simThread_ = new SimulatorThread(simulator_, concentImageFactor, this);
  cellsImg_->setImage(&simThread_->getCellImage());
  concentImg_->setImage(&simThread_->getConcentImage());

  connect(simThread_, SIGNAL(simCycleDone(double)),
    this, SLOT(simCycleDone(double)));

  testing_ = false;
  simulationStatus_ = "Loaded";
  updateStatusText();
}

void SimulatorWindow::scoreResult() {

}

void SimulatorWindow::calcError() {
  loadModel();

  EvaluatorMorphogens evaluatorMorphogens(*search_);

  QElapsedTimer timer;
  timer.start();
  double fMorph = evaluatorMorphogens.evaluate(*model_, 1000);
  QMessageBox::information(this, "Evaluator", QString("Evaluator error (time=%1s): "
    "products = %2").arg(timer.elapsed() / 1000.0)
    .arg(fMorph));
}

void SimulatorWindow::testModel() {
  loadModel();
  nStepsLeft_ = simParams_->NumSimSteps;
  nJumpSteps_ = 0;
  startStopButton_->setText("Stop");
  simulationStatus_ = "Simulating";
  simulating_ = true;
  testing_ = true;
  simulationTimer_.restart();
  simThread_->start();
}

void SimulatorWindow::startStopSim() {
  if (simulator_) {
    if (simulating_) {
      startStopButton_->setText("Start");
      if (testing_)
        simulationStatus_.prepend("Paused - ");
      else
        simulationStatus_ = "Paused";
      simThread_->stopThread();
      QApplication::processEvents(); // process last simCycleDone
      simulating_ = false;
    }
    else {
      if (testing_)
        simulationStatus_.remove(0, 9);
      else
        simulationStatus_ = "Simulating";
      startStopButton_->setText("Stop");
      simulationTimer_.restart();
      simThread_->start();
      simulating_ = true;
    }
  }
}

// Signaled by the simulator thread
void SimulatorWindow::simCycleDone(double change) {
  int nStepsPerCycle = simThread_->getNumStepsPerCycle();
  nStepsLeft_ -= nStepsPerCycle;

  int elapsed = simulationTimer_.restart();
  double nStepsPerSecond = (nStepsPerCycle * 1000.0) / elapsed;

  if (testing_) {
    if (change < simParams_->minConcChange || 
        nStepsLeft_ == 0) { // sim ended
      simThread_->stopThread();
      simulationStatus_ = "Ended";
      startStopButton_->setText("Start");
      testing_ = false;
      simulating_ = false;
      updateStatusText(nStepsPerSecond, change);
      return;
    }
  }

  updateViews();
  updateStatusText(nStepsPerSecond, change);

  int selectedStepsPerSecond = getSelectedSpeed();
  if (selectedStepsPerSecond == maxSpeedSlider()) {
    // Max speed: simulate the number of steps necessary for a refresh rate 
    // of 30 FPS = 33 msPF
    int newNumStepsPerCycle = (33 * nStepsPerCycle) / elapsed;
    static const int maxNumStepsPerCycle = 50000;
    if (newNumStepsPerCycle < 1)
      simThread_->setNumStepsPerCycle(1);
    else if (testing_ && newNumStepsPerCycle > nStepsLeft_)
      simThread_->setNumStepsPerCycle(nStepsLeft_);
    else if (newNumStepsPerCycle > maxNumStepsPerCycle)
      simThread_->setNumStepsPerCycle(maxNumStepsPerCycle);
    else
      simThread_->setNumStepsPerCycle(newNumStepsPerCycle);

    simNextCycle();
  }
  else {
    // A cycle happens every 33 ms (for 30 FPS)
    double selectedStepsPerCycle = 0.033 * selectedStepsPerSecond;
    if (selectedStepsPerCycle >= 1.0) { // More than one step per cycle
      if (testing_ && selectedStepsPerCycle > nStepsLeft_)
        simThread_->setNumStepsPerCycle(nStepsLeft_);
      else
        simThread_->setNumStepsPerCycle(selectedStepsPerCycle);

      while (simulationTimer_.elapsed() + elapsed < 33);

      simNextCycle();
    }
    else { // A cycle last less that the period between steps
      static int sleepTime = 0;
      simThread_->setNumStepsPerCycle(1);
      sleepTime += (33 / selectedStepsPerCycle) - elapsed;
      if (sleepTime > 0)
        QTimer::singleShot(sleepTime, this, SLOT(simNextCycle()));
      else {
        sleepTime = 0;
        simNextCycle();
      }
    }
  }
}

void SimulatorWindow::simNextCycle() {
  simThread_->simNextCycle();
}

void SimulatorWindow::updateStatusText(double sps, double change) {
  Q_ASSERT(simulator_);
  double dist = simulator_->simulatedState().calcDist(*comparisonState_);
  double distLog = simulator_->simulatedState().negativeLogLikelihood(*comparisonState_);

  statusText_->setText(QString("%1 - Dist=%2 (%3) t=%4 s=%5 SPS=%6 cc=%7")
    .arg(simulationStatus_)
    .arg(dist, 0, 'f', 2)
    .arg(distLog, 0, 'f', 2)
    .arg(simulator_->currStep()*simParams_->dt, 0, 'f', 2)
    .arg(simulator_->currStep())
    .arg(sps, 7, 'f', 2, QChar('0'))
    .arg(change, 0, 'e', 1)
    + QString(" (%1x%2)").arg(graphsWidget_->size().width())
    .arg(graphsWidget_->size().height()));
}

void SimulatorWindow::cellsRectDrawed(const QRect &rect, Qt::MouseButton) {
  if (simulator_) {
    bool wasSimulating;
    if (simulating_) {
      wasSimulating = true;
      startStopSim();
    }
    else
      wasSimulating = false;

    int pos1X, pos1Y, pos2X, pos2Y;
    calcImgPositions(rect, pos1X, pos1Y, pos2X, pos2Y);

    SimState &state = simulator_->simulatedState();

    for (int i = pos1X; i <= pos2X; ++i)
      for (int j = pos1Y; j <= pos2Y; ++j) {
        state.deleteCell(i, j);
      }


    if (wasSimulating)
      startStopSim();
    else {
      simThread_->updateImages();
      updateViews();
    }

    updateStatusText();
  }
}

void SimulatorWindow::updateViews() {
  cellsImg_->update();
  concentImg_->update();
  plotsWidget_->updatePlots();

  if (movieSaver_.isRecording())
    movieSaver_.saveFrame();
}

void SimulatorWindow::concentRectDrawed(const QRect &rect, Qt::MouseButton button) {
  if (simulator_) {
    int pos1X, pos1Y, pos2X, pos2Y;
    calcImgPositions(rect, pos1X, pos1Y, pos2X, pos2Y);

    SimState &state = simulator_->simulatedState();

    int n = simulator_->nProducts();
    for (int k = 0; k < n; ++k) {
      Eigen::MatrixXd &conc = state.product(k);

      if (button == Qt::LeftButton) {
        for (int i = pos1X; i <= pos2X; ++i)
          for (int j = pos1Y; j <= pos2Y; ++j)
            conc(i, j) = 0;
      }
      else {
        double sum;
        if (button == Qt::MidButton)
          sum = 0.25;
        else
          sum = -0.25;

        for (int i = pos1X; i <= pos2X; ++i)
          for (int j = pos1Y; j <= pos2Y; ++j)
            conc(i, j) = MathAlgo::max(0.0, conc(i, j) + sum);
      }
    }

    if (!simulating_) {
      simThread_->updateImages();
      updateViews();
    }

    updateStatusText();
  }
}

void SimulatorWindow::concentWheeled(QWheelEvent *event) {
  if (event->delta()>0)
    simThread_->setConcentImageFactor(
      simThread_->getConcentImageFactor()*1.1);
  else
    simThread_->setConcentImageFactor(
      simThread_->getConcentImageFactor()*0.9);

  if (!simulating_) {
    simThread_->updateImages();
    concentImg_->update();
  }
}

void SimulatorWindow::concentLineHMoved() {
  if (simulator_)
    plotsWidget_->setIndH(MathAlgo::min(
    (int)(concentImg_->lineHpos() * simulator_->domainSize().height()),
      simulator_->domainSize().height() - 1));
}

void SimulatorWindow::concentLineVMoved() {
  if (simulator_)
    plotsWidget_->setIndV(MathAlgo::min(
    (int)(concentImg_->lineVpos() * simulator_->domainSize().width()),
      simulator_->domainSize().width() - 1));
}

void SimulatorWindow::calcImgPositions(const QRect &rect, int &pos1X, int &pos1Y,
  int &pos2X, int &pos2Y) {
  if (simulator_) {
    QRect imgRect = cellsImg_->rect();
    float width = simulator_->domainSize().width();
    float height = simulator_->domainSize().height();
    float propX = width / imgRect.width();
    float propY = height / imgRect.height();

    pos1X = propX * rect.left();
    pos1Y = propY * rect.top();

    pos2X = propX * rect.right();
    pos2Y = propY * rect.bottom();

    if (pos1X > pos2X) {
      int i = pos2X;
      pos2X = pos1X;
      pos1X = i;
    }
    if (pos1Y > pos2Y) {
      int i = pos2Y;
      pos2Y = pos1Y;
      pos1Y = i;
    }

    if (pos1X < 0)
      pos1X = 0;

    if (pos2X >= width)
      pos2X = width - 1;

    if (pos1Y < 0)
      pos1Y = 0;

    if (pos2Y >= height)
      pos2Y = height - 1;
  }
}


void SimulatorWindow::translatorTest() {
}

void SimulatorWindow::manipulationTest() {
}

void SimulatorWindow::experimentComboChanged(int index) {
  if (simulating_)
    startStopSim();

  experiment_ = search_->experiment(index);
  loadModel();
}

void SimulatorWindow::editExperiment() {
  ExperimentForm *form = new ExperimentForm(experiment_, db_, this);
  bool ok = form->exec();
  delete form;

  if (ok)
    experiment_->submit(db_);
}

void SimulatorWindow::editModel() {

  if (!modelForm_)
    modelForm_ = new ModelForm(model_, db_,
      search_->simParams()->nInputMorphogens +
      search_->simParams()->nTargetMorphogens, this);

  modelForm_->show();
  modelForm_->raise();
  modelForm_->activateWindow();
}

void SimulatorWindow::speedSliderChanged(int value) {
  if (value == maxSpeedSlider())
    speedText_->setText("Sel. Speed: MAX");
  else
    speedText_->setText(QString("Sel. Speed: %1").arg(getSelectedSpeed()));
}

// The selected speed is a quadratic function in order to make easier to select
// lower speeds
int SimulatorWindow::getSelectedSpeed() {
  int sliderValue = speedSlider_->value();
  return (1.0 / (1.0 + maxSpeedSlider())) * sliderValue * sliderValue +
    maxSpeedSlider() / (1.0 + maxSpeedSlider());
}

void SimulatorWindow::prodsBlockDifChanged(QStandardItem * item) {

}

void SimulatorWindow::prodsBlockProdChanged(QStandardItem * item) {

}

void SimulatorWindow::applyAllProdBlocks() {

}

void SimulatorWindow::uncheckAllProdBlocks() {
  internalBlockProdChange_ = true;

  internalBlockProdChange_ = false;
}

void SimulatorWindow::loadProdBlocks() {
  internalBlockProdChange_ = true;
  static int firstProduct = 4;

  internalBlockProdChange_ = false;
}

void SimulatorWindow::startStopMovie(bool start) {
  if (start) {
    movieButton_->setStyleSheet("QPushButton {color: red; font: bold}");
    movieSaver_.startMovie(graphsWidget_);
  }
  else {
    if (simulating_)
      startStopSim();

    movieSaver_.endMovie();
    movieButton_->setStyleSheet("QPushButton {color: black}");
  }

}

}