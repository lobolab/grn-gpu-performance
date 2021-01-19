// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include <QDialog>
#include <QLabel>
#include <QPushButton>
#include <QSlider>
#include <QComboBox>
#include <QTime>
#include <QStandardItemModel>

#include "UI/GUICommon/moviesaver.h"

namespace LoboLab {

class Individual;
class Search;
class Simulator;
class SimState;
class Model;
class SimParams;
class SimulatorThread;
class CellsView;
class ImageView;
class ConcentPlotWidget;
class Experiment;
class DataTableView;
class ModelForm;
class DB;

class SimulatorWindow : public QDialog {
  Q_OBJECT

  public:
    SimulatorWindow(Individual *ind, Search *search, DB *db, bool autoDelete,
      QWidget *parent = NULL);
    virtual ~SimulatorWindow();

    private slots:
    void calcError();
    void testModel();
    void loadModel();
    void startStopSim();
    void scoreResult();
    void simCycleDone(double change);
    void cellsRectDrawed(const QRect &rect, Qt::MouseButton button);
    void concentRectDrawed(const QRect &rect, Qt::MouseButton button);
    void concentWheeled(QWheelEvent *event);
    void concentLineHMoved();
    void concentLineVMoved();

    void translatorTest();
    void manipulationTest();

    void experimentComboChanged(int index);
    void editExperiment();
    void editModel();

    void speedSliderChanged(int value);
    void simNextCycle();

    void prodsBlockDifChanged(QStandardItem * item);
    void prodsBlockProdChanged(QStandardItem * item);

    void startStopMovie(bool start);

  private:
    void createWidgets();
    void readSettings();
    void writeSettings();
    void calcImgPositions(const QRect &rect, int &pos1X, int &pos1Y, int &pos2X,
      int &pos2Y);
    void updateViews();
    void updateStatusText(double sps = 0.0, double change = 0.0);
    int getSelectedSpeed();
    static inline int maxSpeedSlider() { return 4000; }

    void applyAllProdBlocks();
    void uncheckAllProdBlocks();
    void loadProdBlocks();

    bool isAutoDelete_; // auto delete the window and individual when closed

    DB *db_;
    Individual *individual_;
    Experiment *experiment_;
    Model *model_;
    Search *search_;
    SimParams *simParams_;

    ModelForm *modelForm_;

    QStandardItemModel *prodsBlockDifModel_;
    QStandardItemModel *prodsBlockProdModel_;

    QComboBox *experimentComboBox_;
    QPushButton *movieButton_;

    ImageView *cellsImg_;
    ImageView *concentImg_;

    ConcentPlotWidget *plotsWidget_;

    QWidget *graphsWidget_;
    QLabel *statusText_;
    QLabel *speedText_;
    QSlider *speedSlider_;

    QPushButton *startStopButton_;

    Simulator *simulator_;
    SimulatorThread *simThread_;
    SimState *comparisonState_;

    int nProducts_;
    QList<int> prodLabels_;
    QHash<int, int> prodsBlockLabelIndex_;
    QStringList prodNames_;


    bool simulating_;
    bool testing_;
    QTime simulationTimer_;
    QString simulationStatus_;
    int nStepsLeft_;
    int nJumpSteps_;
    bool internalBlockProdChange_;

    MovieSaver movieSaver_;
  };

}