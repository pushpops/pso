#include <string>
#include <boost/progress.hpp>
#include <fstream>
#include <Eigen/Core>
#include <Eigen/LU>
#include <vector>
#include <random>
#include <iomanip>
#include <sys/stat.h>

class Particle 
{
public:
    Particle();
    void set_x();
    void set_omega(double fg, double threshold);
    void set_p(Eigen::Matrix<double,4,1> x);
    void set_v(Eigen::Matrix<double,4,1> g);
    //各次元でこれまでに到達した位置のxmax、xminを更新し,
    //各粒子の速度の最大値を更新する関数
    void set_vmax();
    void compare_vmax();//vmaxと比較、丸めを行う
    void check_xrange();//探索範囲外に出ないように矯正する
    Eigen::Matrix<double,4,1> get_x();

    Eigen::Matrix<double,4,1> x;//現在の位置ベクトル
    Eigen::Matrix<double,4,1> v;//速度ベクトル

    //位置ベクトルではなく単に各パラメータの最大最小のリスト
    Eigen::Matrix<double,4,1> xmax;//各次元でこれまで到達した最大値を格納
    Eigen::Matrix<double,4,1> xmin;//各次元でこれまで到達した最小値を格納
    Eigen::Matrix<double,4,1> vmax;

    double fp;//粒子自身が今までに発見した最良解の評価値
    double c1;//pに依存するパラメータ
    double c2;//gに依存するパラメータ
    double omega;//前の速度に対する依存度
    static std::string r1info;
    static std::string r2info;
    static double rangealpha;
    static double rangebeta;
    static double rangegamma;
    static double ranget3;
    Position init_param; //初期値
    static Eigen::Matrix<double,4,1> g;//全粒子が今までに発見した最良解
    Eigen::Matrix<double,4,1> r1;//pに依存する乱数パラメータ
    Eigen::Matrix<double,4,1> r2;//gに依存する乱数パラメータ
    Eigen::Matrix<double,4,1> p;//粒子自身が今までに発見した最良解

private:    

};

class Pso
{
public:
  //コンストラクタにカメラ画像から得られる情報を渡す
  Camera(PixelPosition navel,
         const BorderImage &image,
         RenderSetting param,
         BorderSetting border_param); //腹部塗りつぶし画像(オーバーラップ用)

  //推定を任意の回数行い，
  //もっとも良い推定結果の適応度とパラメータ(α，β，γ，t3)を返す関数
  //search_iterateの中で，searchを呼び出している
  //search_iterate(推定回数，(0, αの初期値，βの初期値，γの初期値，t3の初期値))
  std::tuple<Position, cv::Mat> search_iterate(int iteration, Position parameter, std::string camerafilename);
  int pso_max_loop;//ループ回数
  int pnum;//粒子の数
  int psomaxsec;//最大秒数
  int omega_end_fit; //omegaの関数のための探索の終盤を表すfitness
  int omega_flag; //omegaの関数の種類を表す

 private:
  //PSO法でカメラ位置推定を行い，適応度とパラメータ(α，β，γ，t3)を返す関数
  //search((0, αの初期値，βの初期値，γの初期値，t3の初期値))
  std::tuple<Position, cv::Mat> search(Particle *particles, Position parameter);

  double calculate_fitness();

  Particle *particles;//粒子
  bool p_init_flag;
  double pso_threshold;
  double hc_threshold;

};

inline std::string make_directory()
{
    //現在日時を取得する
    time_t t = time(nullptr);

    //形式を変換する
    const tm* lt = localtime(&t);

    //sに独自フォーマットになるように連結していく
    std::stringstream s;
    s<<"20";
    s<<lt->tm_year-100; //100を引くことで20xxのxxの部分になる
    s<<"-"<<lt->tm_mon+1; //月を0からカウントしているため
    s<<"-"<<lt->tm_mday;
    s<<"-"<<lt->tm_hour;
    s<<"-"<<lt->tm_min;
    s<<"-"<<lt->tm_sec;

    //result = "2015-5-19-20-30"
    std::string foldername = s.str();

    const std::string mkdirpathbuf = "mkdir -p ";
    const std::string pathbuf = "/Users/user/Desktop/Experiment/";
    system((mkdirpathbuf+pathbuf).c_str());                   //mkdir .../Experiment

    const std::string pathbuf1 = pathbuf + foldername;        //Users/.../日付
    system((mkdirpathbuf+pathbuf1).c_str());                  //mkdir .../Experiment/日付

    return pathbuf1;
}

#endif //BPM_BLEND_CAMERA_SEARCH_H
