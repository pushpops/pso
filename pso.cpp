#include "PSO.h"
#include <vector>
#include <fstream>
#include <opencv2/opencv.hpp>
#include <iomanip>
#include <sys/stat.h>

using namespace std;
//staticメンバの初期化, 推定するパラメータの範囲
double Particle::rangealpha(30.0);
double Particle::rangebeta(30.0);
double Particle::rangegamma(30.0);
double Particle::ranget3(100);
std::string Particle::r1info = "αβγ:(0.5–1.0) t3:(0.5–1.0)*(rangegamma/ranget3)";
std::string Particle::r2info = "αβγ:(0.5–1.0) t3:(0.5–1.0)*(rangegamma/ranget3)";

Particle::Particle():
    fp(-1.0),
    c1(1.0),
    c2(2.0),
    omega(0.9)
{
        v << 0,0,0,0;
}

void Particle::set_x(){
    x << x(0,0)+v(0,0), x(1,0)+v(1,0), x(2,0)+v(2,0), x(3,0)+v(3,0);
}

void Particle::set_p(Eigen::Matrix<double,4,1> x){
    p << x(0,0), x(1,0), x(2,0), x(3,0);
}

void Particle::set_v(Eigen::Matrix<double,4,1> g){
    v = omega*v + (r1.array()*c1*(p-x).array() + r2.array()*c2*(g-x).array()).matrix();
}

void Particle::set_omega(double fg, double threshold){
    omega = (-0.5/threshold)*fg+0.9; //0.9-0.4
}

Eigen::Matrix<double,4,1> Particle::get_x(){
    return x;
}

//位置推定に必要な変数の初期化
Pso::Pso(){
    pnum = 20;
    psomaxsec = 15;
    pso_max_loop = 1000;
    //omega_end_fitは基本的にはpso_thresholdと同じにするがthresholdが探索の終わりを意味しない時は変更する
    omega_end_fit = 1000000;
    pso_threshold = 1000000000; 
    p_init_flag = 0; //0のとき正規乱数,1のとき一様関数による粒子の初期化
    omega_flag = 0; //omegaの関数の種類.0オフ,1線形減少
}

std::tuple<Position, cv::Mat> Pso::search(Particle *particles, Position parameter)
{
    std::random_device rnd;     // 非決定的な乱数生成器を生成
    unsigned int seed = rnd();  // saigenjikkennwosurutokihaseednoataiwoshiteisuru
    std::mt19937 mt(seed);      // メルセンヌ・ツイスタの32ビット版、引数は初期シード値
    CameraRender rendering(render_setting_);
    rendering.render(navel, parameter);

    Position pso_position;
    Eigen::Matrix<double,4,1> g;

    const std::string pathbuf1 = "mkdir ./result_camera_search";
    const std::string pathbuf2 = "./result_camera_search";
    const std::string pathbuf3 = "/result.csv";
    const std::string pathbuf4 = "/Image";
    const std::string pathbuf5 = "./result_camera_search/Image/best_of_pso.png";
    system((pathbuf1).c_str());
    system((pathbuf1+pathbuf4).c_str());
    std::ofstream ofs;
    ofs.open(pathbuf2+pathbuf3);
    boost::timer t;

    int pso_t=0;//最大反復回数までのカウンタ
    g << 0,0,0,0;//全ての粒子の最良解

    double fg = 0;//全ての粒子の最良解の適応度
    double f_tmp;//算出した適応度を一時的に保存して比較時に用いる
    int best_t = 0;
    Position best_position;
    double stime = 0;

    //乱数生成の準備
    std::uniform_real_distribution<> unif01(0, 1);     // [0,1] 範囲の一様乱数
    std::uniform_real_distribution<> unif11(-1, 1);     // [-1,1] 範囲の一様乱数

    //引数は範囲，一様乱数の分布関数
    std::uniform_real_distribution<> unifa(-particles[0].rangealpha/2+parameter.GetAlpha(), particles[0].rangealpha/2+parameter.GetAlpha());
    std::uniform_real_distribution<> unifb(-particles[0].rangebeta/2+parameter.GetBeta(), particles[0].rangebeta/2+parameter.GetBeta());
    std::uniform_real_distribution<> unifg(-particles[0].rangegamma/2+parameter.GetGamma(), particles[0].rangegamma/2+parameter.GetGamma());
    std::uniform_real_distribution<> unift(-particles[0].ranget3/2+parameter.GetT3(), particles[0].ranget3/2+parameter.GetT3());

    //https://cpprefjp.github.io/reference/random/normal_distribution.html
    //平均mu、標準偏差sigmaで分布させる正規乱数.片側範囲のため/2,99.7%の確率で探索範囲内の値を生成するため/3よって/6
    std::normal_distribution<> norma(parameter.GetAlpha(), particles[0].rangealpha/6);
    std::normal_distribution<> normb(parameter.GetBeta(), particles[0].rangebeta/6);
    std::normal_distribution<> normg(parameter.GetGamma(), particles[0].rangegamma/6);
    std::normal_distribution<> normt(parameter.GetT3(), particles[0].ranget3/6);

    //実験の設定をcsvファイルに保存する
    std::ofstream ofs2;
    ofs2.open("./result_camera_search/param_info.csv");

    ofs2 << "loop_max,"<< pso_max_loop << "\n"
         << "particle_num," << pnum << "\n"
         << "psomaxtime," << psomaxsec << "\n"
         << "pso_max_loop," << pso_max_loop << "\n"
         << "pso_threshold," << pso_threshold << "\n"
         << "init_type," << p_init_flag << "\n"
         << "omega_type," << omega_flag << "\n"
         << "omegaendfit,"<< omega_end_fit << "\n"
         << "hc_threshold," << hc_threshold << "\n"
         << "set_alpha,"<< parameter.GetAlpha() << "\n"
         << "set_beta," << parameter.GetBeta() << "\n"
         << "set_gamma," << parameter.GetGamma() << "\n"
         << "set_t3," << parameter.GetT3() << "\n"
         << "set_navelx," << navel.GetX() << "\n"
         << "set_navely," << navel.GetY() << "\n"
         << "c1," << particles[0].c1 << "\n"
         << "c2," << particles[0].c2 << "\n"
         << "r1," << particles[0].r1info << "\n"
         << "r2,"<< particles[0].r2info  << std::endl;
    ofs2 << std::fixed << std::setprecision(0) << "seed,"<< seed << std::endl;
    ofs2.close();

    do{ 
        for(int i=0; i<pnum; i++){
            //r1,r2に0-1の一様乱数を格納
            particles[i].r1(0) = unif01(mt);
            particles[i].r1(1) = unif01(mt);
            particles[i].r1(2) = unif01(mt);
            particles[i].r1(3) = (particles[i].rangegamma/particles[i].ranget3)*unif01(mt); //(rangegamma/ranget3)*
            particles[i].r2(0) = unif01(mt);
            particles[i].r2(1) = unif01(mt);
            particles[i].r2(2) = unif01(mt);
            particles[i].r2(3) = (particles[i].rangegamma/particles[i].ranget3)*unif01(mt);
            particles[i].init_param = parameter; //各粒子にUIから与えた初期parameterを保存

            //xの初期値を設定
            if(p_init_flag == 1){//粒子の位置を一様乱数で初期化する場合
                particles[i].x << unifa(mt), unifb(mt), unifg(mt), unift(mt);

            }else if(p_init_flag == 0){ //粒子の位置を正規乱数で初期化.
                particles[i].x << norma(mt),normb(mt),normg(mt),normt(mt);
            }

            //適応度計算のためのセット
            pso_position.SetAlpha(particles[i].get_x()(0,0));
            pso_position.SetBeta(particles[i].get_x()(1,0));
            pso_position.SetGamma(particles[i].get_x()(2,0));
            pso_position.SetT3(particles[i].get_x()(3,0));
            rendering.render(navel, pso_position);//推定過程のパラメータでレンダリング
            f_tmp = calculate_fitness(border_image_, rendering.getGene_img());//現在の適応度
            particles[i].set_p(particles[i].x);
            particles[i].fp = f_tmp;

            //personal bestとglobal bestの更新
//            if(f_tmp > particles[i].fp){    //適応度の比較
//                particles[i].set_p(particles[i].x);

            if(f_tmp > fg){
                g = particles[i].x;
                fg = f_tmp;
            }
        }
        std::cout << "Particles initialized fg=" << fg << ", seed=" << seed << std::endl;
    }while(fg < 0); //loop0でfgが一定値を超えるまで初期化をやり直す

    //粒子No.0はカメラ,メジャーの測定値を初期値とする．
    particles[0].x << particles[0].init_param.GetAlpha(),particles[0].init_param.GetBeta(),
                        particles[0].init_param.GetGamma(),particles[0].init_param.GetT3();
    //particles[0].v << 0, 0, 0, 0; //初期値を持つ粒子No.0が最初に不動にならないようにする
    ofs <<"loop," << "particle," << "alpha," << "beta,"<< "gamma,"<< "t3,"<< "fp," << "fg," << "time"<< std::endl;
    std::cout <<"loop," << "particle," << "alpha," << "beta,"<< "gamma,"<< "t3,"<< "fp," << "fg," << "time" << std::endl;
    while ((pso_t < pso_max_loop && fg < pso_threshold && stime < psomaxsec)){
        for(int i = 0; i < pnum; i++){

            //位置を更新
            particles[i].set_x();   //x(t+1) = x(t)+v(t)

            //適応度計算のためのセット
            pso_position.SetAlpha(particles[i].get_x()(0,0));
            pso_position.SetBeta(particles[i].get_x()(1,0));
            pso_position.SetGamma(particles[i].get_x()(2,0));
            pso_position.SetT3(particles[i].get_x()(3,0));
            rendering.render(navel, pso_position);//推定過程のパラメータでレンダリング
            f_tmp = calculate_fitness(border_image_, rendering.getGene_img());//現在の適応度

            //personal bestとglobal bestの更新
            if(f_tmp > particles[i].fp){    //適応度の比較
                particles[i].set_p(particles[i].x);
                particles[i].fp = f_tmp;   //抜けていた・・・・・・
                if(f_tmp > fg){
                    g = particles[i].x;
                    best_position.SetAlpha(particles[i].get_x()(0,0));
                    best_position.SetBeta(particles[i].get_x()(1,0));
                    best_position.SetGamma(particles[i].get_x()(2,0));
                    best_position.SetT3(particles[i].get_x()(3,0));
                    best_t = pso_t;
                    fg = f_tmp;
                    std::cout << pso_t << "," << i << "," << particles[i].get_x()(0,0) << "," << particles[i].get_x()(1,0) << ","
                        << particles[i].get_x()(2,0) << ","<< particles[i].get_x()(3,0) << "," << particles[i].fp << "," << fg << "," << stime <<  std::endl;
                }
            }

            if(omega_flag == 1){//omegaの線形減少
                particles[i].set_omega(fg,omega_end_fit);
            }
            particles[i].set_v(g);  //v(t)=omegav(t-1)...


            //記録
            stime = t.elapsed();

            ofs << pso_t << "," << i << "," << particles[i].get_x()(0,0) << "," << particles[i].get_x()(1,0) << ","
                << particles[i].get_x()(2,0) << "," << particles[i].get_x()(3,0) << "," << particles[i].fp << "," << fg << "," << stime << std::endl;

        }
        //rendering.render(navel, best_position);
        //cv::imwrite(pathbuf4+pathbuf2+"pso_"+std::to_string(pso_t)+".png", rendering.getGene_img());

        pso_t++;
    }

    double pso_time = t.elapsed();
    pso_position.SetAlpha(g(0,0));
    pso_position.SetBeta(g(1,0));
    pso_position.SetGamma(g(2,0));
    pso_position.SetT3(g(3,0));
    pso_position.SetFitness(fg);
    rendering.render(navel, pso_position);
    cv::imwrite(pathbuf5, rendering.getGene_img());
    std::string overwrappath = "./result_camera_search/Image/overlap_best_of_pso.png";
    make_overwrapping(body_image_, rendering.getGene_img(), overwrappath);

    }
    ofs.close();
    
}

double Pso::calculate_fitness(cv::Mat true_image, cv::Mat rendered_image) {
  ;
}

