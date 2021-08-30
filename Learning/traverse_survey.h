#pragma once
class traverse_survey
{
private:
	const double pi = 3.1415926;
public:
	traverse_survey();  //导线测量的构造函数
	~traverse_survey(); // 导线测量的析构函数
	void h_d(double* d, double* m, double* s, int N, double* Decimal_degrees);//把观测的度分秒转化了弧度使用
	void transfrom_d_m_s(double Decimal_degrees); // 十进制的度转化为度分秒
	void cal_trigo_function(double Decimal_degrees); // 计算三角函数值
	void position_angle(double x, double y);  //求方位角
	void start_end_position(double A_x, double A_y, double B_x, double B_y, double& start_position);//求起始和末位的坐标方位角。
	void sum_distance(int N, double* distance, double& sum);//求导线的总长
	void Angle_position(double start_positon, int N, double* Decimal_degrees, double* every_position);//求每一站的坐标方位角
	void Correct_Angle_position(int N, double* every_position, double end_position, double* correct_every_position);//改正后的每一站的坐标方位角
	void d_x_y(int N, double* distance, double* correct_every_position, double* d_x, double* d_y);//求d_x_y的值
	void Cf_x(int N, double* d_x, double A_x, double C_x, double& f_x);//求f_x
	void Cf_y(int N, double* d_y, double A_y, double C_y, double& f_y); //求f_y
	void Cf_s(double sum, double f_x, double f_y); //求全长闭合差和全长相对闭合差
	void Correct_d_x_y(int N, double* d_x, double* d_y, double* distance, double sum, double f_x, double f_y, double* correct_d_x, double* correct_d_y);//求改正之后的d_x_y
	void x_y(int N, double A_x, double A_y, double* correct_d_x, double* correct_d_y, long double* lx, long double* ly);//求x,y;

};

