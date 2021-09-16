#include "PSO.h"

using namespace std;
//staticメンバの初期化, 推定するパラメータの両側範囲
double Particle::rangealpha(5.0);
std::string Particle::r1info = "αβγ:(0.5–1.0)";
std::string Particle::r2info = "αβγ:(0.5–1.0)";

Particle::Particle():
    fp(20000.0),
    c1(1.0),
    c2(2.0),
    omega(0.9)
{}

void Particle::set_x(){
    x = x+v;
}

void Particle::set_p(Eigen::Vector2d x){
    p = x;
}

void Particle::set_v(Eigen::Vector2d g){
    v = omega*v + (r1.array()*c1*(p-x).array() + r2.array()*c2*(g-x).array()).matrix();
}

void Particle::set_omega(double fg, double threshold){
    omega = 0.9/threshold*fg; //0.9-0.2
}

Eigen::Vector2d Particle::get_x(){
    return x;
}

//位置推定に必要な変数の初期化
Pso::Pso(int pnum, int psomaxsec, int pso_max_loop, bool p_init_flag,bool omega_flag):
    pnum(pnum),
    psomaxsec(psomaxsec),
    pso_max_loop(pso_max_loop),
    omega_delta(20000),  
    p_init_flag(p_init_flag),
    omega_flag(omega_flag) 
{
    particles = new Particle[pnum];
}

void Pso::search()
{
    std::random_device rnd;     // 非決定的な乱数生成器を生成
    unsigned int seed = rnd();  // saigenjikkennwosurutokihaseednoataiwoshiteisuru
    std::mt19937 mt(seed);      // メルセンヌ・ツイスタの32ビット版、引数は初期シード値

    const std::string pathbuf1 = "mkdir -p ./result_search";
    const std::string pathbuf2 = "./result_search";
    const std::string pathbuf3 = "/result.csv";
    system((pathbuf1).c_str());
    std::ofstream ofs;
    ofs.open(pathbuf2+pathbuf3);
    // boost::timer t;

    int pso_t=0;//最大反復回数までのカウンタ
    double fg = 20000;//全ての粒子の最良解の適応度
    double f_tmp;//算出した適応度を一時的に保存して比較時に用いる
    int best_t = 0;
    double stime = 0;

    g << 20000,20000;
  
    //乱数生成の準備
    std::uniform_real_distribution<> unif01(0, 1);     // [0,1] 範囲の一様乱数
    std::uniform_real_distribution<> unif11(-1, 1);     // [-1,1] 範囲の一様乱数

    //引数は範囲，一様乱数の分布関数
    std::uniform_real_distribution<> unifa(-particles[0].rangealpha/2, particles[0].rangealpha/2);

    //https://cpprefjp.github.io/reference/random/normal_distribution.html
    //平均mu、標準偏差sigmaで分布させる正規乱数.片側範囲のため/2,99.7%の確率で探索範囲内の値を生成するため/3よって/6
    std::normal_distribution<> norma(particles[0].rangealpha/2, particles[0].rangealpha/6);

    //実験の設定をcsvファイルに保存する
    std::ofstream ofs2;
    ofs2.open("./result_search/param_info.csv");
    ofs2 << "loop_max,"<< pso_max_loop << "\n"
         << "particle_num," << pnum << "\n"
         << "psomaxtime," << psomaxsec << "\n"
         << "pso_max_loop," << pso_max_loop << "\n"
         << "init_type," << p_init_flag << "\n"
         << "omega_type," << omega_flag << "\n"
         << "omega_delta,"<< omega_delta << "\n"
         << "hc_threshold," << hc_threshold << "\n"
         << "c1," << particles[0].c1 << "\n"
         << "c2," << particles[0].c2 << "\n"
         << "r1," << particles[0].r1info << "\n"
         << "r2,"<< particles[0].r2info  << std::endl;
    ofs2 << std::fixed << std::setprecision(0) << "seed,"<< seed << std::endl;
    ofs2.close();

    cout << "0" << endl;
    //初期化
    do{ 
        for(int i=0; i<pnum; i++){
            //r1,r2に0-1の一様乱数を格納
            particles[i].r1(0) = unif01(mt);
            particles[i].r1(1) = unif01(mt);
            particles[i].r2(0) = unif01(mt);
            particles[i].r2(1) = unif01(mt);

            //xの初期値を設定
            if(p_init_flag == 1){//粒子の位置を一様乱数で初期化する場合
                particles[i].x << unifa(mt),unifa(mt);

            }else if(p_init_flag == 0){ //粒子の位置を正規乱数で初期化.
                particles[i].x << norma(mt),norma(mt);
            }        

            f_tmp = calculate_fitness(particles[i].get_x());//現在の適応度  

            //personal bestとglobal bestの更新
           if(f_tmp < particles[i].fp){    //適応度の比較
               particles[i].set_p(particles[i].x);
               particles[i].fp = f_tmp;
            }

            if(f_tmp < fg){
                g = particles[i].x;
                fg = f_tmp;
            }
        }

        std::cout << "Particles initialized fg=" << fg << ", seed=" << seed << std::endl;
    }while(fg > 20000); //loop0でfgが一定値を下回るまで初期化をやり直す
    
    ofs << "iterate" << "," << "particleNo" << "," << "x1" << "," << "x2" << ","
                << "fp"<< "," << "fg" << "," << "time" << std::endl;

    while ((pso_t < pso_max_loop  && stime < psomaxsec)){
        for(int i = 0; i < pnum; i++){

            //位置を更新
            particles[i].set_x();   //x(t+1) = x(t)+v(t)

            f_tmp = calculate_fitness(particles[i].get_x());//現在の適応度

            //personal bestとglobal bestの更新
            if(f_tmp < particles[i].fp){    //適応度の比較
                particles[i].set_p(particles[i].x);
                particles[i].fp = f_tmp;
                if(f_tmp < fg){
                    g = particles[i].x;                   
                    best_t = pso_t;
                    fg = f_tmp;
                }
            }

            if(omega_flag == 1){//omegaの線形減少
                particles[i].set_omega(fg,omega_delta);
            }
            particles[i].set_v(g);

            //記録
            // stime = t.elapsed();

            ofs << pso_t << "," << i << "," << particles[i].get_x()(0) << "," << particles[i].get_x()(1) << ","
                << particles[i].fp << "," << fg << "," << stime << std::endl;

        }
        pso_t++;
    }
    ofs.close();
}

double Pso::calculate_fitness(Eigen::Vector2d p) {
    //Beale function x1=3, x2=0.5 best=0
  return pow(1.5-p(0)+p(0)*p(1),2)+pow(2.25-p(0)+p(0)*pow(p(1),2),2)+pow(2.625-p(0)+p(0)*pow(p(1),3),2);
}

