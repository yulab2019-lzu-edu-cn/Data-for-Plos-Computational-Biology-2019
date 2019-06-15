#ifndef IZHIKEVICH_SPNET_H
#define IZHIKEVICH_SPNET_H

#include <QMainWindow>
#include <QDir>

#include "qcustomplot.h"

#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <random>
#include <numeric>
#include <string>
#include <cmath>
#include <time.h>
#include <regex>
#include <map>

using std::vector;

namespace Ui {
class izhikevich_SPNET;
}

class izhikevich_SPNET : public QMainWindow {
    Q_OBJECT

public:
    explicit izhikevich_SPNET(QWidget *parent = 0);
    ~izhikevich_SPNET();

    void Control(int it, std::pair<QVector<double>, QVector<double>> temp_Pair);
    void PrintPlot_0_1(QCustomPlot *SPNET_Plot_0_1, std::pair<QVector<double>, QVector<double>> temp_Pair);
    void PrintPlot_0_2(QCustomPlot *SPNET_Plot_0_2, std::pair<QVector<double>, QVector<double>> temp_Pair);
    void PrintPlot_0_3(QCustomPlot *SPNET_Plot_0_3, std::pair<QVector<double>, QVector<double>> temp_Pair);
    void PrintPlot_0_4(QCustomPlot *SPNET_Plot_0_4, std::pair<QVector<double>, QVector<double>> temp_Pair);
    void Display_LCD();

    int Get_Random_int(const int &max);
    int Get_Random_int(const int &min, const int &max);
    double Get_Random_real();
    double Noise();
	vector<vector<double>> Transposition(const vector<vector<double>> &temp_vec);

    void Data();
    void Data_d();
    void Dti();
    void Process_t_Num();
    void Save_t_Num();
    void Save_Post();
    void Save_Weight();
    void Data_Bold();
    void Save_Bold();
    void Branching_Parameter();

    std::pair<QVector<double>, QVector<double>> Data_firings(const int &h); //for plot
    std::pair<QVector<double>, QVector<double>> Data_Num_t(const int &h); //for plot
    std::pair<QVector<double>, QVector<double>> Data_Bold(const int &h); //for plot
    std::pair<QVector<double>, QVector<double>> Data_Power_Law(); //for plot
    void Data_I(); //for plot
    std::pair<QVector<double>, QVector<double>> Data_It(); //for plot

    void Initialize();
    void Initialize_Rest();
    void Simulation();
    void Start();

private slots:
    void on_pushButton_clicked();

    void on_spinBox_editingFinished();

    void on_pushButton_2_clicked();

    void on_doubleSpinBox_2_editingFinished();

    void on_spinBox_2_editingFinished();

private:
    Ui::izhikevich_SPNET *ui;

    const int Neu_num = 90; // number of encephalic regions
    const int Ne = 100; //excitatory neurons
    const int Ni = 25; // inhibitory neurons
    const int N = Ne + Ni; // total number of neurons
    const int M = N / 10; // the number of synapses per neuron
    const int D = 20; // maximal axonal conduction delay
    vector<vector<int>> Post; //[Neu_num * N][M] indeces of postsynaptic neurons
    vector<vector<double>> Weight; //[Neu_num * N][Neu_num * N] matrix of synaptic weights
    vector<vector<double>> a, b, c, d, V, U; //[Neu_num][N] neuronal dynamics parameters and activity variables
    vector<int> N_firings, N_firings_temp; // the number of fired neurons
    const int N_firings_max = 100 * N; // upper limit on the number of fired neurons per sec
    vector<vector<vector<int>>> firings, firings_temp; //[Neu_num][N_frings_max][2] indeces and timings of spikes
    vector<vector<double>> I; // [Neu_num][N] Current
    QVector<double> I_t;
    vector<vector<vector<double>>> Save_I; // [Neu_num][N][N] Current of 1000ms
    vector<vector<double>> Save_tNum, Save_tNum2, Save_tNum_cut; //[Neu_num][] number of firing per 1ms
    long long sec = 0; //Simulation time
    long long T = 1; //Total simulaton time
    double count_sec = 0.0; //Sum of sec
    double Weight_Max = 0.01;  // maximal synaptic strength
    QVector<double> vec_f_x, vec_f_y0; //Save sec and firings
    QVector<double> vec_Nf_x, vec_Nf_y0; //Save sec and N_firings
    vector<int> temp_line;
    vector<vector<int>> vec_Dti; //Dti data
    bool Click = false; //on_spinBox_2_editingFinished()
    vector<std::map<long long, long long>> map_size_count_temp;
    std::map<long long, long long> map_size_count; //Raw data of power law
    const double epsilon = 0.2, kappa = 0.65, gamma = 0.41, tau = 0.98, alpha = 0.32, rho = 0.34, V_0 = 0.02; //The value associated with Bold
    vector<vector<double>> vec_s, vec_f, vec_v, vec_q, vec_y; //The function associated with Bold;
    const int sampling_10 = 10; //ms
    const int sampling_50 = 50; //ms
    const int sampling_100 = 100; //ms
    const int sampling_150 = 150; //ms
    const int sampling_200 = 200; //ms
    const int sampling_500 = 500; //ms
    const int sampling_1000 = 1000; //ms
    const int sampling_2000 = 2000; //ms
    QString dat = ".dat";
};

#endif // IZHIKEVICH_SPNET_H
