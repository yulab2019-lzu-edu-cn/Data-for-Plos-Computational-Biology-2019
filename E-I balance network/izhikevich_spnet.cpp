#include "izhikevich_spnet.h"
#include "ui_izhikevich_spnet.h"

izhikevich_SPNET::izhikevich_SPNET(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::izhikevich_SPNET)
{
    ui->setupUi(this);

    Dti();

    Initialize();

    ui->spinBox->setMinimum(1);
    ui->spinBox->setMaximum(864000);
    ui->doubleSpinBox_2->setMinimum(0.01);
    ui->doubleSpinBox_2->setMaximum(25);
    ui->spinBox_2->setMinimum(1);
    ui->spinBox_2->setMaximum(90);

    ui->lcdNumber_2->setDecMode();
    ui->lcdNumber_2->setDigitCount(6);
    ui->lcdNumber_2->setSegmentStyle(QLCDNumber::Flat);
    ui->lcdNumber->setDecMode();
    ui->lcdNumber->setDigitCount(6);
    ui->lcdNumber->setSegmentStyle(QLCDNumber::Flat);
}

izhikevich_SPNET::~izhikevich_SPNET()
{
    delete ui;
}

void izhikevich_SPNET::Control(int it, std::pair<QVector<double>, QVector<double>> temp_Pair) {
    switch (it) {
    case 1:
        ui->SPNET_Plot_0_1->clearGraphs();
        PrintPlot_0_1(ui->SPNET_Plot_0_1, temp_Pair);
        ui->SPNET_Plot_0_1->replot();
        break;
    case 2:
        ui->SPNET_Plot_0_2->clearGraphs();
        PrintPlot_0_2(ui->SPNET_Plot_0_2, temp_Pair);
        ui->SPNET_Plot_0_2->replot();
        break;
    case 3:
        ui->SPNET_Plot_0_3->clearGraphs();
        PrintPlot_0_3(ui->SPNET_Plot_0_3, temp_Pair);
        ui->SPNET_Plot_0_3->replot();
        break;
    case 4:
        ui->SPNET_Plot_0_4->clearGraphs();
        PrintPlot_0_4(ui->SPNET_Plot_0_4, temp_Pair);
        ui->SPNET_Plot_0_4->replot();
        break;
    default:
        break;
    }
}

void izhikevich_SPNET::PrintPlot_0_1(QCustomPlot *SPNET_Plot_0_1, std::pair<QVector<double>, QVector<double>> temp_Pair) {
    //SPNET_Plot_0_1->setOpenGl(true);
    SPNET_Plot_0_1->addGraph(SPNET_Plot_0_1->xAxis, SPNET_Plot_0_1->yAxis);
    SPNET_Plot_0_1->graph(0)->addData(temp_Pair.first, temp_Pair.second);
    SPNET_Plot_0_1->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, QPen(Qt::black), QBrush(Qt::black), 2));
    SPNET_Plot_0_1->graph(0)->setLineStyle(QCPGraph::lsNone);
    SPNET_Plot_0_1->graph(0)->rescaleAxes();
    SPNET_Plot_0_1->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
}

void izhikevich_SPNET::PrintPlot_0_2(QCustomPlot *SPNET_Plot_0_2, std::pair<QVector<double>, QVector<double>> temp_Pair) {
    //SPNET_Plot_0_2->setOpenGl(true);
    SPNET_Plot_0_2->addGraph(SPNET_Plot_0_2->xAxis, SPNET_Plot_0_2->yAxis);
    SPNET_Plot_0_2->graph(0)->addData(temp_Pair.first, temp_Pair.second);
    SPNET_Plot_0_2->graph(0)->rescaleAxes();
    SPNET_Plot_0_2->addGraph(SPNET_Plot_0_2->xAxis, SPNET_Plot_0_2->yAxis);
    SPNET_Plot_0_2->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
}

void izhikevich_SPNET::PrintPlot_0_3(QCustomPlot *SPNET_Plot_0_3, std::pair<QVector<double>, QVector<double>> temp_Pair) {
    //SPNET_Plot_0_3->setOpenGl(true);
    SPNET_Plot_0_3->addGraph(SPNET_Plot_0_3->xAxis, SPNET_Plot_0_3->yAxis);
    SPNET_Plot_0_3->graph(0)->addData(temp_Pair.first, temp_Pair.second);
    //SPNET_Plot_0_3->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, QPen(Qt::blue), QBrush(Qt::blue), 2));
    //SPNET_Plot_0_3->graph(0)->setLineStyle(QCPGraph::lsNone);
    SPNET_Plot_0_3->graph(0)->rescaleAxes();
    SPNET_Plot_0_3->addGraph(SPNET_Plot_0_3->xAxis, SPNET_Plot_0_3->yAxis);
    for (int i = -20; i <= 400; ++i) {
        SPNET_Plot_0_3->graph(1)->addData(i / 100, -1.5 * (i / 100) + 6);
    }
    SPNET_Plot_0_3->graph(1)->setPen(QPen(Qt::black));
    SPNET_Plot_0_3->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
}

void izhikevich_SPNET::PrintPlot_0_4(QCustomPlot *SPNET_Plot_0_4, std::pair<QVector<double>, QVector<double>> temp_Pair) {
    //SPNET_Plot_0_4->setOpenGl(true);
    SPNET_Plot_0_4->addGraph(SPNET_Plot_0_4->xAxis, SPNET_Plot_0_4->yAxis);
    SPNET_Plot_0_4->graph()->addData(temp_Pair.first, temp_Pair.second);
    SPNET_Plot_0_4->graph()->setPen(QPen(Qt::red));
    SPNET_Plot_0_4->graph()->rescaleAxes();
    SPNET_Plot_0_4->addGraph(SPNET_Plot_0_4->xAxis, SPNET_Plot_0_4->yAxis);
    SPNET_Plot_0_4->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
}

void izhikevich_SPNET::Display_LCD() {
    ui->lcdNumber_2->display(count_sec);
    ui->lcdNumber->display(N_firings[ui->spinBox_2->value() - 1]);
}

int izhikevich_SPNET::Get_Random_int(const int &max) {
//	static std::default_random_engine e(time(nullptr));
//	std::uniform_int_distribution<int> u(0, max - 1);

//	return u(e);

    return Get_Random_int(0, max - 1);
}

int izhikevich_SPNET::Get_Random_int(const int &min, const int &max) {
    static std::default_random_engine e(time(nullptr));
    std::uniform_int_distribution<int> u(min, max);

    return u(e);
}

double izhikevich_SPNET::Get_Random_real() {
    static std::default_random_engine e(time(nullptr));
    std::uniform_real_distribution<double> u(0, 1);

    return u(e);
}

double izhikevich_SPNET::Noise() {
//    static std::default_random_engine e(time(nullptr)); //原内容
//    std::normal_distribution<double> u(0, 1); //Gaussian Distribution

//    return /*5.0 * */u(e);

    static std::default_random_engine e(time(nullptr));
    std::uniform_real_distribution<double> u(0, 1);

    auto a = u(e);
    auto b = u(e);
    const auto pi = 3.14159265359;

    return std::sqrt(-4 * std::pow(10.0, -3) * std::log(a)) * std::cos(2 * pi * b);
}

vector<vector<double>> izhikevich_SPNET::Transposition(const vector<vector<double>> &temp_vec) {
    vector<vector<double>> temp_Trans;
    temp_Trans.resize(temp_vec[0].size());
    for (auto &i : temp_Trans) {
        i.resize(temp_vec.size());
    }

    for (auto i = 0; i != temp_vec.size(); ++i) {
        for (auto j = 0; j != temp_vec[0].size(); ++j) {
            temp_Trans[j][i] = temp_vec[i][j];
        }
    }

    return temp_Trans;
}

void izhikevich_SPNET::Process_t_Num() {  //counting fire pre ms
    for (auto h = 0; h != Neu_num; ++h) {
        std::vector<int> vec_t;
        for (auto i = 0; i != N_firings[h]; ++i) {
            if (firings[h][i][0] >= 0) {
                vec_t.push_back(firings[h][i][0]); //time point
            }
        }

        auto count_temp = 0;
        for (auto i = 0, k = 1; i != 1000; ++i, ++k) {
            auto j = std::count(vec_t.begin(), vec_t.end(), i); //How many time points, the current time how many spiking.

            Save_tNum[h].push_back(j);

            count_temp += j;
            if (k % 2 == 0) {
                Save_tNum2[h].push_back(count_temp);
                count_temp = 0;
            }

            if(j >= 7) {
                Save_tNum_cut[h].push_back(j);
            }else{
                Save_tNum_cut[h].push_back(0);
            }
        }
    }
}

std::pair<QVector<double>, QVector<double>> izhikevich_SPNET::Data_firings(const int &h) {
    vec_f_x.clear();
    vec_f_y0.clear();

    for (auto i = 0; i != N_firings_temp[h]; ++i) {
        if (firings_temp[h][i][0] >= 0) {
            vec_f_x.push_back(firings_temp[h][i][0]);
            vec_f_y0.push_back(firings_temp[h][i][1]);
        }
    }

    if (count_sec == T) {
        std::ofstream out("D:\\Qt_projects\\data9\\Firings\\F" + std::to_string(Weight_Max) + ".dat"), out1; //out,just delet early data.
        out1.open("D:\\Qt_projects\\data9\\Firings\\F" + std::to_string(Weight_Max) + ".dat", std::ios::app);
        std::ofstream out_BP("D:\\Qt_projects\\data9\\B_P\\F" + std::to_string(Weight_Max) + ".dat"), out1_BP; //out,just delet early data.
        out1_BP.open("D:\\Qt_projects\\data9\\B_P\\F" + std::to_string(Weight_Max) + ".dat", std::ios::app);
        for (auto i : vec_f_x) {
            out1 << i << " ";
            out1_BP << i << " ";
        }
        out1 << std::endl;
        out1_BP << std::endl;
        for (auto i : vec_f_y0) {
            out1 << i << " ";
            out1_BP << i << " ";
        }
        out1 << std::endl;
        out1.close();
        out1_BP << std::endl;
        out1_BP.close();
    }

    return std::make_pair(vec_f_x, vec_f_y0);
}

std::pair<QVector<double>, QVector<double>> izhikevich_SPNET::Data_Num_t(const int& h) {
    //std::ofstream out("tNum" + std::to_string(h + 1) + ".dat");
    QVector<double> vec_t, vec_num;
    auto temp_t = 0;

    for (int i = (count_sec  - 1)* 1000; i != Save_tNum[h].size(); ++i) { //current sec
        vec_num.push_back(Save_tNum[h][i]);
        vec_t.push_back(++temp_t);
    }

    return std::make_pair(vec_t, vec_num);
}

std::pair<QVector<double>, QVector<double>> izhikevich_SPNET::Data_Bold(const int& h) {
    //std::ofstream out("bold" + std::to_string(h + 1) + ".dat");
    QVector<double> vec_t, vec_bold;
    auto temp_t = 0;

    for (auto k = 0; k != Neu_num; ++k) {
        vec_s[k].clear();
        vec_f[k].clear();
        vec_v[k].clear();
        vec_q[k].clear();
        vec_y[k].clear();

        vec_s[k].resize(Save_tNum[k].size());
        vec_f[k].resize(Save_tNum[k].size());
        vec_v[k].resize(Save_tNum[k].size());
        vec_q[k].resize(Save_tNum[k].size());
        vec_y[k].resize(Save_tNum[k].size());

        vec_f[k][0] = 1;
        vec_q[k][0] = 1;
        vec_v[k][0] = 1;

        for (auto i = 0; i != Save_tNum[k].size() - 1; ++i) {
            vec_s[k][i + 1] = vec_s[k][i] + 0.05 * (epsilon * Save_tNum[k][i] - kappa * vec_s[k][i] - gamma *(vec_f[k][i] - 1));
            vec_f[k][i + 1] = vec_f[k][i] + 0.05 * vec_s[k][i];
            vec_v[k][i + 1] = vec_v[k][i] + 0.05 * (vec_f[k][i] - (pow(vec_v[k][i], (1 / alpha)))) / tau;
            vec_q[k][i + 1] = vec_q[k][i] + 0.05 * ((vec_f[k][i] * (1 - pow((1 - rho), 1 / vec_f[k][i]))) / rho - (vec_q[k][i] * pow(vec_v[k][i], 1 / alpha)) / vec_v[k][i]) / tau;
            vec_y[k][i + 1] = V_0 * (7 * rho * (1 - vec_q[k][i + 1]) + 2 * (1 - vec_q[k][i + 1] / vec_v[k][i + 1]) + (2 * rho - 0.2) * (1 - vec_v[k][i + 1]));
        }
    }

    auto count_t = 1;
    for (auto i : vec_y[h]) {
        if(/*count_t % sampling_100 == 0 && */count_t >= 100) { //erase first data
            vec_bold.push_back(i);
            vec_t.push_back(++temp_t);
        }
        ++count_t;

//        vec_bold.push_back(i);
//        vec_t.push_back(++temp_t);

        //out << i << " ";
    }

    return std::make_pair(vec_t, vec_bold);
}

void izhikevich_SPNET::Data_Bold() {
    for (auto k = 0; k != Neu_num; ++k) {
        vec_s[k].clear();
        vec_f[k].clear();
        vec_v[k].clear();
        vec_q[k].clear();
        vec_y[k].clear();

        vec_s[k].resize(Save_tNum[k].size());
        vec_f[k].resize(Save_tNum[k].size());
        vec_v[k].resize(Save_tNum[k].size());
        vec_q[k].resize(Save_tNum[k].size());
        vec_y[k].resize(Save_tNum[k].size());

        vec_f[k][0] = 1;
        vec_q[k][0] = 1;
        vec_v[k][0] = 1;

        for (auto i = 0; i != Save_tNum[k].size() - 1; ++i) {
            vec_s[k][i + 1] = vec_s[k][i] + 0.05 * (epsilon * Save_tNum[k][i] - kappa * vec_s[k][i] - gamma *(vec_f[k][i] - 1));
            vec_f[k][i + 1] = vec_f[k][i] + 0.05 * vec_s[k][i];
            vec_v[k][i + 1] = vec_v[k][i] + 0.05 * (vec_f[k][i] - (pow(vec_v[k][i], (1 / alpha)))) / tau;
            vec_q[k][i + 1] = vec_q[k][i] + 0.05 * ((vec_f[k][i] * (1 - pow((1 - rho), 1 / vec_f[k][i]))) / rho - (vec_q[k][i] * pow(vec_v[k][i], 1 / alpha)) / vec_v[k][i]) / tau;
            vec_y[k][i + 1] = V_0 * (7 * rho * (1 - vec_q[k][i + 1]) + 2 * (1 - vec_q[k][i + 1] / vec_v[k][i + 1]) + (2 * rho - 0.2) * (1 - vec_v[k][i + 1]));
        }
    }
}

std::pair<QVector<double>, QVector<double>> izhikevich_SPNET::Data_Power_Law() {
    QVector<double> vec_size, vec_count;

    auto count = 0;
    for (auto h = 0; h != Neu_num; ++h) {
        auto count_t = 0;
        for (auto i : Save_tNum[h]) {
            if (count_t >= 100) {
                if (i >= 1) {
                    count += i;
                }
                if (count > 0 && i == 0) {
                    ++map_size_count_temp[h][count];
                    ++map_size_count[count];
                    count = 0;
                }
            }
            ++count_t;
        }
    }

    for (auto i : map_size_count) {
        if (i.second >= Neu_num) {
            vec_size.push_back(log10(i.first ));
            vec_count.push_back(log10(static_cast<double>(i.second) / Neu_num)); //average
        }
    }

    if (count_sec == T) {
        std::ofstream out("D:\\Qt_projects\\data9\\P_L\\P_L" + std::to_string(Weight_Max) + ".dat"), out1; //out,just delet early data.
        out1.open("D:\\Qt_projects\\data9\\P_L\\P_L" + std::to_string(Weight_Max) + ".dat", std::ios::app);
        std::ofstream out_BP("D:\\Qt_projects\\data9\\B_P\\P_L" + std::to_string(Weight_Max) + ".dat"), out1_BP; //out,just delet early data.
        out1_BP.open("D:\\Qt_projects\\data9\\B_P\\P_L" + std::to_string(Weight_Max) + ".dat", std::ios::app);
        for (auto i : map_size_count) {
            if (i.second >= Neu_num) {
                out1 << i.first << " ";
                out1_BP << i.first << " ";
            }
        }
        out1 << std::endl;
        out1_BP << std::endl;
        for (auto i : map_size_count) {
            if (i.second >= Neu_num) {
                out1 << static_cast<double>(i.second) / Neu_num << " "; //average
                out1_BP << static_cast<double>(i.second) / Neu_num << " "; //average
            }
        }
        out1 << std::endl;
        out1.close();
        out1_BP << std::endl;
        out1_BP.close();

        std::ofstream out2("D:\\Qt_projects\\data9\\P_L_each\\P_L_each" + std::to_string(Weight_Max) + ".dat"), out3; //out2,just delet early data.
        std::ofstream out2_BP("D:\\Qt_projects\\data9\\B_P\\P_L_each" + std::to_string(Weight_Max) + ".dat"), out3_BP; //out2,just delet early data.
        for (auto h = 0; h != Neu_num; ++h) {
            out3.open("D:\\Qt_projects\\data9\\P_L_each\\P_L_each" + std::to_string(Weight_Max) + ".dat", std::ios::app);
            out3_BP.open("D:\\Qt_projects\\data9\\B_P\\P_L_each" + std::to_string(Weight_Max) + ".dat", std::ios::app);
            for (auto i : map_size_count_temp[h]) {
                out3 << i.first << " ";
                out3_BP << i.first << " ";
            }
            out3 << std::endl;
            out3_BP << std::endl;
            for (auto i : map_size_count_temp[h]) {
                out3 << i.second << " ";
                out3_BP << i.second << " ";
            }
            out3 << std::endl;
            out3.close();
            out3_BP << std::endl;
            out3_BP.close();
        }
    }

    return std::make_pair(vec_size, vec_count);
}

void izhikevich_SPNET::Data_I() {
    auto temp_I = 0.0;
    for (auto h = 0; h != Neu_num; ++h) {
        for (auto i : I[h]) {
            temp_I += i;
        }
    }
    I_t.push_back(temp_I);
}

std::pair<QVector<double>, QVector<double>> izhikevich_SPNET::Data_It() {
    if (count_sec == T) {
        std::ofstream out("D:\\Qt_projects\\data9\\I\\I" + std::to_string(Weight_Max) + ".dat");
        I_t.erase(I_t.begin(), I_t.begin() + 200);
        auto sum_I = 0.0;
        for (auto i : I_t) {
            sum_I += i;
        }
        out << sum_I;
    }

    QVector<double> vec_t;
    for (auto i = 0; i != I_t.size(); ++i) {
        vec_t.push_back(i);
    }

    return std::make_pair(vec_t, I_t);
}

void izhikevich_SPNET::Branching_Parameter() {
    if (count_sec == T) {
        auto count = 0;
        vector<double> sigma_pd, sigma_i, sigma_Neu;
        vector<std::pair<double, double>> sigma_ad;
        auto sigma = 0.0, N_descendants = 0.0, N_ancestors = 0.0, Sum_na = 0.0, Sum_nad = 0.0;
        std::ofstream out("D:\\Qt_projects\\data9\\B_P\\B_P" + std::to_string(Weight_Max) + ".dat");

        for (auto h = 0; h != Neu_num; ++h) {
            auto temp_t = 1, count_t = 0;
            sigma_pd.clear();
            sigma_i.clear();            
            for(auto i  = 0; i != Save_tNum_cut[h].size(); ++i) {
                if (Save_tNum_cut[h][i] >= 1) {
                    count = Save_tNum_cut[h][i];

                    if(count_t < 2 && temp_t % 2 != 0) { //a
                        ++count_t;
                        N_ancestors = count;
                        Sum_na += N_ancestors; //Σa

                        ++temp_t;
                        count = 0;

                        if((i + 1) == Save_tNum_cut[h].size() || Save_tNum_cut[h][i + 1] == 0) { //0206...
                            N_descendants = 0;
                            sigma_ad.push_back(std::make_pair(N_ancestors, std::round(N_descendants / N_ancestors))); //d = round(nd/na)
                            ++temp_t;
                        }
                    }

                    if(count > 0 && temp_t % 2 == 0) { //d
                        ++count_t;
                        if(count_t <= 2) {
                            N_descendants = count;
                            Sum_nad += N_ancestors; //Σa|d
                            sigma_ad.push_back(std::make_pair(N_ancestors, std::round(N_descendants / N_ancestors))); //d = round(nd/na)
                            ++temp_t;
                        }
                    }

                }else {
                    count_t = 0;
                }
            }
            auto Find = [&](const int &i) {
                return std::find_if(sigma_ad.begin(), sigma_ad.end(), [=](const std::pair<double, double> &p) {
                    return i == p.second;
                });
            };

            auto Pd = [&](const unsigned long long &i) { //P(d)
                sigma_pd.push_back((Sum_nad / Sum_na) * (N - 1) / (N - sigma_ad[i].first));
                return (Sum_nad / Sum_na) * (N - 1) / (N - sigma_ad[i].first);
            };

            for(auto i = 0; i <= N; ++i) {
                auto at = std::distance(sigma_ad.begin(), Find(i));

                if (at != sigma_ad.size()) {
                    sigma_i.push_back(i * Pd(at));
                }else{
                    sigma_i.push_back(0);
                }
            }
            sigma_Neu.push_back(std::accumulate(sigma_i.begin(), sigma_i.end(), 0.0));
        }
        sigma = std::accumulate(sigma_Neu.begin(), sigma_Neu.end(), 0.0) / sigma_Neu.size(); //average

        out << sigma;
    }
}

void izhikevich_SPNET::Save_t_Num() {
    if (count_sec == T) {
        std::ofstream out("D:\\Qt_projects\\data9\\tNum\\tNum" + std::to_string(Weight_Max) + ".dat");
         std::ofstream out_BP("D:\\Qt_projects\\data9\\B_P\\tNum" + std::to_string(Weight_Max) + ".dat");

        for (auto h = 0; h != Neu_num; ++h) {
            auto count_t = 0;
            for (auto i : Save_tNum[h]) {
                if(count_t >= 60) { //erase first data
                    out << i << " ";
                    if(count_t % sampling_100 == 0) {
                        out_BP << i << " ";
                    }
                }
                ++count_t;
            }
            out << std::endl;
            out_BP << std::endl;
        }
    }
}

void izhikevich_SPNET::Save_Post() {
    std::ofstream out("Post.dat");

    for(auto i : Post) {
        for(auto j : i) {
            out << j << " ";
        }
        out << std::endl;
    }
}

void izhikevich_SPNET::Save_Weight() {
    std::ofstream out("Weight.dat");

    for(auto i : Weight) {
        for(auto j : i) {
            out << j << " ";
        }
        out << std::endl;
    }
}

void izhikevich_SPNET::Save_Bold() {
    Data_Bold();

    if (count_sec == T) {
        std::ofstream out_1("D:\\Qt_projects\\data9\\1_200\\Bold_o" + std::to_string(Weight_Max) + ".dat");
        std::ofstream out_10("D:\\Qt_projects\\data9\\10_200\\Bold" + std::to_string(Weight_Max) + ".dat");
        std::ofstream out_50("D:\\Qt_projects\\data9\\50_200\\Bold" + std::to_string(Weight_Max) + ".dat");
        std::ofstream out_100("D:\\Qt_projects\\data9\\100_200\\Bold" + std::to_string(Weight_Max) + ".dat");
        std::ofstream out_BP("D:\\Qt_projects\\data9\\B_P\\Bold" + std::to_string(Weight_Max) + ".dat");
        std::ofstream out_150("D:\\Qt_projects\\data9\\150_200\\Bold" + std::to_string(Weight_Max) + ".dat");
        std::ofstream out_200("D:\\Qt_projects\\data9\\200_200\\Bold" + std::to_string(Weight_Max) + ".dat");
        std::ofstream out_500("D:\\Qt_projects\\data9\\500_200\\Bold" + std::to_string(Weight_Max) + ".dat");
        std::ofstream out_1000("D:\\Qt_projects\\data9\\1000_200\\Bold" + std::to_string(Weight_Max) + ".dat");
        std::ofstream out_2000("D:\\Qt_projects\\data9\\2000_200\\Bold" + std::to_string(Weight_Max) + ".dat");

        for (auto h = 0; h != Neu_num; ++h) {
            auto count_t = 0;
            for (auto i : vec_y[h]) {
                if(count_t >= 60) { //erase first data
                    out_1 << i << " ";
                    if(count_t % sampling_10 == 0) {
                        out_10 << i << " ";
                    }
                    if(count_t % sampling_50 == 0) {
                        out_50 << i << " ";
                    }
                    if(count_t % sampling_100 == 0) {
                        out_100 << i << " ";
                        out_BP << i << " ";
                    }
                    if(count_t % sampling_150 == 0) {
                        out_150 << i << " ";
                    }
                    if(count_t % sampling_200 == 0) {
                        out_200 << i << " ";
                    }
                    if(count_t % sampling_500 == 0) {
                        out_500 << i << " ";
                    }
                    if(count_t % sampling_1000 == 0) {
                        out_1000 << i << " ";
                    }
                    if(count_t % sampling_2000 == 0) {
                        out_2000 << i << " ";
                    }
                }

                ++count_t;

                //out << i << " ";
            }
            out_1 << std::endl;
            out_10 << std::endl;
            out_50 << std::endl;
            out_100 << std::endl;
            out_BP << std::endl;
            out_150 << std::endl;
            out_200 << std::endl;
            out_500 << std::endl;
            out_1000 << std::endl;
            out_2000 << std::endl;
        }
    }
}

void izhikevich_SPNET::Data() {
    Display_LCD();

    Control(1, Data_firings(ui->spinBox_2->value() - 1));
    Control(2, Data_Num_t(ui->spinBox_2->value() - 1));
    Control(3, Data_Power_Law());
    //Control(4, Data_It());
    Control(4, Data_Bold(ui->spinBox_2->value() - 1));

    Save_Bold();
    Branching_Parameter();
    Save_t_Num();
}

void izhikevich_SPNET::Data_d() {
    Display_LCD();

    Control(1, Data_firings(ui->spinBox_2->value() - 1));
    Control(2, Data_Num_t(ui->spinBox_2->value() - 1));
    Control(4, Data_Bold(ui->spinBox_2->value() - 1));
}

void izhikevich_SPNET::Dti() {
    std::string line;
    std::ifstream in("dti.txt");
    std::regex pat_regex("[[:digit:]]+");

    while (std::getline(in, line)) {
        for (std::sregex_iterator it(line.begin(), line.end(), pat_regex), end_it; it != end_it; ++it) {
            temp_line.push_back(std::stoi(it->str()));
        }
        vec_Dti.push_back(temp_line);
        temp_line.clear();
    }
}

void izhikevich_SPNET::Initialize() {
    Post.resize(Neu_num * N);
    //inPost
    int r;
    bool exists;
    for (auto h = 0; h != Neu_num; ++h) {
        auto pre_i = h * N, post_i = pre_i + N, Ne_i = Ne + pre_i;
        for (auto i = pre_i; i != post_i; ++i) { // create random synapses
            for (auto j = 0; j != M; ++j) {
                do {
                    exists = false; // avoid multiple synapses
//                    if (i < Ne_i) {
//                        r = Get_Random_int(pre_i, post_i - 1);  // exc->exc, exc->inh
//                    }
//                    if (i >= Ne_i) {
//                        r = Get_Random_int(Ne_i, post_i - 1); //  inh->exc only
//                        //r = Get_Random_int(pre_i, post_i - 1);
//                    }

                    r = Get_Random_int(pre_i, post_i - 1);  // exc->exc, exc->inh, inh->inh, inh->exc

                    if (r == i) {
                        exists = true; // no self-synapses
                    }
                    for (auto k = 0; k != j; ++k) {
                        if (Post[i][k] == r) {
                            exists = true; // synapse already exists
                        }
                    }
                } while (exists);
                Post[i].push_back(r);
            }
        }

        //outPost
        for (auto k = 0; k != Neu_num; ++k) {
            auto pre_i_2 = k * N, Ne_i_2 = Ne + k * N;
            if (k != h) {
                for (auto j = 0; j != vec_Dti[h][k]; ++j) {
                    Post[Get_Random_int(pre_i, Ne_i - 1)].push_back(Get_Random_int(pre_i_2, Ne_i_2 - 1));
                }
            }
        }
    }

    Save_Post();
}

void izhikevich_SPNET::Initialize_Rest() {
    vec_s.resize(Neu_num);
    vec_f.resize(Neu_num);
    vec_v.resize(Neu_num);
    vec_q.resize(Neu_num);
    vec_y.resize(Neu_num);
    Weight.clear();
    a.clear();
    b.clear();
    c.clear();
    d.clear();
    V.clear();
    U.clear();
    I.clear();
    I_t.clear();
    Save_I.clear();
    Save_tNum.clear();
    Save_tNum2.clear();
    Save_tNum_cut.clear();
    map_size_count.clear();
    Weight.resize(Neu_num * N);
    map_size_count_temp.clear();
    map_size_count_temp.resize(Neu_num);
    for (auto &i : Weight) { i.resize(Neu_num * N); }
    a.resize(Neu_num);
    for (auto &i : a) { i.resize(N); }
    b.resize(Neu_num);
    for (auto &i : b) { i.resize(N); }
    c.resize(Neu_num);
    for (auto &i : c) { i.resize(N); }
    d.resize(Neu_num);
    for (auto &i : d) { i.resize(N); }
    V.resize(Neu_num);
    for (auto &i : V) { i.resize(N); }
    U.resize(Neu_num);
    for (auto &i : U) { i.resize(N); }
    firings.resize(Neu_num);
    for (auto &i : firings) { i.resize(N_firings_max); }
    for (auto &i : firings) { for (auto &j : i) { j.resize(2); } }
    I.resize(Neu_num);
    for (auto &i : I) { i.resize(N); }
    Save_I.resize(Neu_num);
    for (auto &i : Save_I) { i.resize(1000); }
    Save_tNum.resize(Neu_num);
    Save_tNum2.resize(Neu_num);
    Save_tNum_cut.resize(Neu_num);
    for (auto h = 0; h != Neu_num; ++h) {
        for (auto i = 0; i != Ne; ++i) a[h][i] = 0.02;// RS type
        for (auto i = Ne; i != N; ++i) a[h][i] = 0.02 + 0.08 * Get_Random_real(); // FS type

        for (auto i = 0; i != Ne; ++i) b[h][i] = 0.2;// RS type
        for (auto i = Ne; i != N; ++i) b[h][i] = 0.2 - 0.02 * Get_Random_real(); // FS type

        for (auto i = 0; i != Ne; ++i) c[h][i] = -65.0 + 15.0 * pow(Get_Random_real(), 2);// RS type
        for (auto i = Ne; i != N; ++i) c[h][i] = -65.0; // FS type

        for (auto i = 0; i != Ne; ++i) d[h][i] = 8.0 - 6.0 * pow(Get_Random_real(), 2); // RS type
        for (auto i = Ne; i != N; ++i) d[h][i] = 2.0; // FS type

        for (auto i = 0; i != N; ++i)  V[h][i] = -65.0;	// initial values for V
        for (auto i = 0; i != N; ++i)  U[h][i] = 0.2 * V[h][i];	// initial values for U
    }
//    std::ofstream outa("a.dat");
//    for (auto i : Transposition(a)) {
//        for(auto j : i) {
//            outa << std::setw(10) << std::left << j << " ";
//        }
//        outa << std::endl;
//    }
    N_firings.resize(Neu_num);
    for (auto &i : N_firings) { i = 1; } // spike timings
    for (auto i = 0; i != Neu_num; ++i) {
        firings[i][0][0] = -D;	// put a dummy spike at -D for simulation efficiency
        firings[i][0][1] = 0;	// index of the dummy spike
    }

    //Weight
    for (auto h = 0; h != Neu_num; ++h) {
        auto pre_i = h * N, post_i = pre_i + N, Ne_i = Ne + pre_i;
        for (auto i = pre_i; i != post_i; ++i) { // synaptic weights
            for (auto j = 0; j != Post[i].size(); ++j) {
                if (i < Ne_i && static_cast<int>(Post[i][j] / N) == h && static_cast<int>(Post[i][j] % N) <= (Ne - 1)) {
                    Weight[i][Post[i][j]] = Weight_Max;  // exc->exc, dynamic
                }
                if (i < Ne_i && static_cast<int>(Post[i][j] / N) == h && static_cast<int>(Post[i][j] % N) > (Ne - 1)) {
                    Weight[i][Post[i][j]] = 4.0;  // exc->inh
                }
                if (i < Ne_i && static_cast<int>(Post[i][j] / N) != h && static_cast<int>(Post[i][j] % N) <= (Ne - 1)) {
                    Weight[i][Post[i][j]] = 4.0; // area->area
                }
                if (i >= Ne_i) {
                    Weight[i][Post[i][j]] = -1.0;  // inh->inh, inh->exc
                }
            }
        }
    }

    //Save_Weight();
}

void izhikevich_SPNET::Simulation() {
    for (; sec != T; ++sec) {	// simulation of T sec
        for (auto t = 0; t != 1000; ++t) {			// simulation of 1 sec
            for (auto h = 0; h != Neu_num; ++h) {
                for (auto i = 0; i != N; ++i) I[h][i] = 0.0;	// reset the input
                for (auto i = 0; i != 2; ++i) I[h][Get_Random_int(N)] = 20.0;		// random thalamic input
                for (auto i = 0; i != N; ++i) {
                    if (V[h][i] >= 30) {					// did it fire
                        V[h][i] = c[h][i];					// voltage reset
                        U[h][i] += d[h][i];					// recovery variable reset
                        firings[h][N_firings[h]][0] = t;      // record spiking time
                        firings[h][N_firings[h]][1] = i + h * N;      // record spiking number in that time
                        ++N_firings[h];
                        try {
                            if (N_firings[h] == N_firings_max) {
                                throw std::runtime_error("Two many spikes at t=");
                            }
                        }
                        catch (std::runtime_error &err) {
                            std::cout << err.what() << t << " (ignoring all)";
                            N_firings[h] = 1;
                        }
                    }
                }

                auto k = N_firings[h] - 1;
                auto spiking_t = firings[h][k][0]; // spiking time
                while (firings[h][k][0] == spiking_t && spiking_t >= 0) { // all spiking neuron at that time
                    auto i = firings[h][k][1];  // presynaptic
                    for (auto m = 0; m != Post[i].size(); ++m) {
                        auto j = Post[i][m]; // postsynaptic
                        auto temp_j_0 = j / N;
                        auto temp_j_1 = j % N;

                        I[temp_j_0][temp_j_1] += Weight[i][j];
                    }
                    --k;
                }

                for (auto i = 0; i != N; ++i) {
                    V[h][i] += 0.5*((0.04 * V[h][i] + 5) * V[h][i] + 140 - U[h][i] + I[h][i]) + Noise(); // for numerical stability 改过括号
                    V[h][i] += 0.5*((0.04 * V[h][i] + 5) * V[h][i] + 140 - U[h][i] + I[h][i]) + Noise(); // time step is 0.5 ms 改过括号
                    U[h][i] += a[h][i] * (b[h][i] * V[h][i] - U[h][i]);
                }
            }

            //Data_I();
        }

        ++count_sec; //Sum of sec;

        firings_temp = firings;
        N_firings_temp = N_firings;

        Process_t_Num(); // must befor Data().

        Data(); //must befor rest firings.

        for (auto h = 0; h != Neu_num; ++h) {
            for (auto i = 0; i != N_firings[h]; ++i) { // reset firings
                firings[h][i][0] = -D;
                firings[h][i][1] = 0;
            }
            N_firings[h] = 1;  // reset N_firings
        }
    }
}

void izhikevich_SPNET::Start() {
    for(auto w = 8.0; w <= 11.0; w += 0.1) {
        T = 30;
        Weight_Max = w;
        sec = 0;
        count_sec = 0.0;
        ui->lcdNumber_2->display(count_sec);
        ui->doubleSpinBox_2->setValue(w);
        ui->spinBox->setValue(T);
        ui->spinBox_2->setValue(10);

        Initialize_Rest();
        Simulation();
    }
}

void izhikevich_SPNET::on_pushButton_clicked() {
//    if (count_sec == 0.0 || T < count_sec) {
//        count_sec = 0.0;
//        sec = 0;

//        Initialize_Rest();
//        Simulation();
//    }
//    else {
//        Simulation();
//    }

//    Click = true;

    Start();
}

void izhikevich_SPNET::on_spinBox_editingFinished() {
//    T = ui->spinBox->value();
}

void izhikevich_SPNET::on_pushButton_2_clicked() {
//    ++T;
//    on_pushButton_clicked();
//    ui->spinBox->setValue(T);
//    T = ui->spinBox->value();
}

void izhikevich_SPNET::on_doubleSpinBox_2_editingFinished() {
//    if (Weight_Max != ui->doubleSpinBox_2->value()) {
//        Weight_Max = ui->doubleSpinBox_2->value();
//        sec = 0;
//        count_sec = 0.0;
//        ui->lcdNumber_2->display(count_sec);
//        ui->spinBox->setValue(1);
//        T = 1;

//        Initialize_Rest();
//    }
}

void izhikevich_SPNET::on_spinBox_2_editingFinished() {
//    if (Click == true) {
//        Data_d();
//    }
}
