#include "Pso.h"

int main(int argc, char *argv[])
{
	int pnum = 20;
	int psomaxsec = 5;
	int pso_max_loop = 100;
	bool p_init_flag = 0; //0のとき正規乱数,1のとき一様関数による粒子の初期化
	bool omega_flag = 0; //omegaの関数の種類.0オフ,1線形減少
	Pso pso(pnum, psomaxsec, pso_max_loop, p_init_flag, omega_flag);
	pso.search();
	return 0;
}