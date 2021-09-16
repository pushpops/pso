#include <string>
// #include <boost/progress.hpp>
#include <fstream>
#include <Eigen/Core>
#include <Eigen/LU>
#include <vector>
#include <random>
#include <iomanip>
#include <sys/stat.h>
#include <iostream>

class Particle 
{
public:
    Particle();
    void set_x();
    void set_omega(double fg, double threshold);
    void set_p(Eigen::Vector2d x);
    void set_v(Eigen::Vector2d g);
    //各次元でこれまでに到達した位置のxmax、xminを更新し,
    //各粒子の速度の最大値を更新する関数
    void set_vmax();
    void compare_vmax();//vmaxと比較、丸めを行う
    void check_xrange();//探索範囲外に出ないように矯正する
    Eigen::Vector2d get_x();

    Eigen::Vector2d x;//現在の位置ベクトル
    Eigen::Vector2d v = Eigen::VectorXd::Zero(2);//速度ベクトル

    //位置ベクトルではなく単に各パラメータの最大最小のリスト
    Eigen::Vector2d xmax;//各次元でこれまで到達した最大値の位置を格納
    Eigen::Vector2d xmin;//各次元でこれまで到達した最小値の位置を格納
    Eigen::Vector2d vmax;

    double fp;//粒子自身が今までに発見した最良解の評価値
    double c1;//pに依存するパラメータ
    double c2;//gに依存するパラメータ
    double omega;//前の速度に対する依存度
    static std::string r1info;
    static std::string r2info;
    static double rangealpha;
    static Eigen::Vector2d g;//全粒子が今までに発見した最良解
    Eigen::Vector2d r1;//pに依存する乱数パラメータ
    Eigen::Vector2d r2;//gに依存する乱数パラメータ
    Eigen::Vector2d p;//粒子自身が今までに発見した最良解

private:    

};

class Pso
{
public:
  Pso(int pnum, int psomaxsec, int pso_max_loop, bool p_init_flag, bool omega_flag); 
  int pso_max_loop;//ループ回数
  int pnum;//粒子の数
  int psomaxsec;//最大秒数
  int omega_delta;
  int omega_flag; //omegaの関数の種類を表す
  void search();

 private:
  //PSO法でカメラ位置推定を行い，適応度とパラメータ(α，β，γ，t3)を返す関数
  //search((0, αの初期値，βの初期値，γの初期値，t3の初期値))
  

  double calculate_fitness(Eigen::Vector2d p);

  Particle *particles;//粒子
  bool p_init_flag;
  double pso_threshold;
  double hc_threshold;
  Eigen::Vector2d g = Eigen::VectorXd::Zero(2);//全粒子が今までに発見した最良解

};

