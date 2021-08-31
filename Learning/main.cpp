#include "traverse_survey.h"
#include<math.h>
#include<iostream>
#include<iomanip>
using namespace std;
traverse_survey::traverse_survey()
{
}
traverse_survey::~traverse_survey()
{
	cout << "-----------------------结束---------------------------" << endl;
}
//------------------  把观测的度分秒转化了十进制的度使用------------------------//
//*d----------------输入的度
//*m----------------输入的分
//*s----------------输入的秒
//N-----------------站数
//*Decimal_degrees----------保存的每一站转化的观测角度的--十进制的度。
void traverse_survey::h_d(double* d, double* m, double* s, int N, double* Decimal_degrees)
{
	for (int i = 0; i < N; i++)
	{
		Decimal_degrees[i] = d[i] + m[i] / 60 + s[i] / 3600;
		cout << "第" << i + 1 << "个夹角为：" << Decimal_degrees[i] << endl;
	}
}
//------------------十进制度 转化为度分秒输出------------------------//
//Decimal_degrees--------------十进制的度
void traverse_survey::transfrom_d_m_s(double Decimal_degrees)
{
	double d = 0, m = 0, s = 0;
	d = int(Decimal_degrees);
	m = int((Decimal_degrees - d) * 60);
	s = (Decimal_degrees - d - m / 60) * 3600;
	cout << d << "°" << m << "′" << s << "″" << endl;
}
//------------------求每一站的坐标方位角-和角度闭合差--------------------------//
//start_positon--------------起始方位角
//N--------------------------站数
//*Decimal_degrees-----------每一站的观测角度
//*every_position------------保存每一站的方位角 (十进制的度)
void traverse_survey::Angle_position(double start_positon, int N, double* Decimal_degrees, double* every_position)
{
	for (int i = 0; i < N; i++)
	{
		if (i == 0)
		{
			every_position[0] = start_positon + Decimal_degrees[0] + 180;
			if (every_position[0] < 0)
			{
				every_position[0] += 360;
			}
			else if (every_position[0] > 360)
			{
				every_position[0] -= 360;
			}
		}
		else
		{
			every_position[i] = every_position[i - 1] + Decimal_degrees[i] + 180;
			if (every_position[i] < 0)
			{
				every_position[i] += 360;
			}
			else if (every_position[i] > 360)
			{
				every_position[i] -= 360;
			}
		}
		cout << "第" << i + 1 << "个坐标方位角为：";
		transfrom_d_m_s(every_position[i]);
		//cout << endl;
	}
}
//------------------------------改正后的每一站的坐标方位角------------------------------//
//N------------------------------站数
//*every_position  --------------每一站的坐标方位角
//end_position-------------------最后一站的坐标方位角
// *correct_every_position-------改正后的每一站坐标方位角
void traverse_survey::Correct_Angle_position(int N, double* every_position, double end_position, double* correct_every_position)
{
	double temp = 0;
	cout << "该导线的角度闭合差为：";
	transfrom_d_m_s(every_position[N - 1] - end_position);
	temp = (every_position[N - 1] - end_position) / N;

	for (int i = 0; i < N - 1; i++)
	{
		correct_every_position[i] = every_position[i] - (i + 1) * temp;
		cout << "第" << i + 1 << "个改正后的坐标方位角为：";
		transfrom_d_m_s(correct_every_position[i]);
		//cout << endl;
	}
}
//---------------------求d_x_y的值----------------------------//
//N------------------------------------------站数
//*distance----------------------------------每一站的距离
//*correct_every_position--------------------每一站改正后的坐标方位角
//*d_x---------------------------------------存储每一站的d_x
//*d_y---------------------------------------存储每一站的d_y
void traverse_survey::d_x_y(int N, double* distance, double* correct_every_position, double* d_x, double* d_y)
{
	for (int i = 0; i < N - 1; i++)
	{
		d_x[i] = distance[i] * cos(correct_every_position[i] * (pi / 180));
		d_y[i] = distance[i] * sin(correct_every_position[i] * (pi / 180));
		cout << i + 1 << "号测站的d_x和d_y分别为：" << d_x[i] << "， " << d_y[i] << endl;
	}
}
//-------------------------导线的总长--------------------------//
//N-----------------站数
//*distance---------每一站的距离
//sum --------------导线的总长
void traverse_survey::sum_distance(int N, double* distance, double& sum)
{
	for (int i = 0; i < N - 1; i++)
	{
		sum += distance[i];
	}
	cout << "导线的总长：" << sum << endl;
}

//---------------------------求f_x-----------------------------//
//N---------------------站数
//*d_x------------------每一站的d_x
//A_x-------------------起始点的x坐标
//C_x-------------------终止点的x坐标
//f_x-----------------f_x的值
void traverse_survey::Cf_x(int N, double* d_x, double A_x, double C_x, double& f_x)
{
	double sum = 0;
	double temp = A_x - C_x;
	for (int i = 0; i < N - 1; i++)
	{
		sum += d_x[i];
	}
	f_x = temp + sum;
	cout << "f_x的值为：" << f_x << endl;
}
//---------------------------求f_y-----------------------------//
//N---------------------站数
//*d_y------------------每一站的d_y
//A_y-------------------起始点的y坐标
//C_x-------------------终止点的y坐标
//f_y-----------------f_y的值
void traverse_survey::Cf_y(int N, double* d_y, double A_y, double C_y, double& f_y)
{
	double sum = 0;
	double temp = A_y - C_y;
	for (int i = 0; i < N - 1; i++)
	{
		sum += d_y[i];
	}
	f_y = temp + sum;
	cout << "f_y的值为：" << f_y << endl;
}
//---------------------求全长闭合差和全长相对闭合差---------------------//
//f_x-------------------------导线f_x的值
//f_y-------------------------导线f_y的值
void traverse_survey::Cf_s(double sum, double f_x, double f_y)
{
	double f_s = sqrt(pow(f_x, 2) + pow(f_y, 2));
	cout << "全长闭合差f_s的值为:  " << f_s << endl;
	double K = 1 / (sum / f_s);
	cout << "全长相对闭合差K为：" << K << endl;
	cout << "1/K = " << int(1 / K) << endl;
	/*if (int(1 / K) > 4000)
	{
		cout << "该导线没有超限！" << endl;
	}
	else
	{
		cout << "该导线超限！" << endl;
	}*/

}
//-------------------------求改正之后的d_x_y----------------------//
//N----------------------------站数
//*d_x-------------------------每一站的d_x的值
//*d_y-------------------------每一站的d_y的值
//*distance--------------------每一站的距离
//sum--------------------------导线的总长
//f_x--------------------------f_x的值
//f_y--------------------------f_y的值
void traverse_survey::Correct_d_x_y(int N, double* d_x, double* d_y, double* distance, double sum, double f_x, double f_y, double* correct_d_x, double* correct_d_y)
{
	double* temp_x = new double[N - 1]{ 0 };        
	double* temp_y = new double[N - 1]{ 0 };
	for (int i = 0; i < N - 1; i++)
	{
		temp_x[i] = f_x * (distance[i] / sum);
		temp_y[i] = f_y * (distance[i] / sum);
		cout << "-----------------------------------------------------------------" << endl;
		cout << "第" << i + 1 << "个的" << "改正数d_x和d_y的值为：" << temp_x[i] << ", " << temp_y[i] << endl;
		correct_d_x[i] = d_x[i] - temp_x[i];
		correct_d_y[i] = d_y[i] - temp_y[i];
		cout << "第" << i + 1 << "个的" << "改正之后的dx和dy的值为：" << correct_d_x[i] << ", " << correct_d_y[i] << endl;
	}
}
//-------------------------------------求x,y--------------------------------//;
//N---------------------------站数
//A_x-------------------------起始的x坐标
//A_y-------------------------起始的y坐标
//*correct_d_x----------------改正后的d_x
//*correct_d_y----------------改正后的d_y
//*x--------------------------储存每一站x的坐标
//*y--------------------------储存每一站y的坐标
void traverse_survey::x_y(int N, double A_x, double A_y, double* correct_d_x, double* correct_d_y, long double* x, long double* y)//求x,y;
{
	for (int i = 0; i < N - 1; i++)
	{
		if (i == 0)
		{
			x[0] = correct_d_x[0] + A_x;
			y[0] = correct_d_y[0] + A_y;
		}
		else
		{
			x[i] = x[i - 1] + correct_d_x[i];
			y[i] = y[i - 1] + correct_d_y[i];
		}
		cout << i + 1 << "号测站的坐标为 " << x[i] << ", " << y[i] << endl;

	}
}
//求起始和末位的坐标方位角。
//A_x,A_y,B_x, B_y  起始边
//C_x,C_y,D_x, D_y  终止边
//start_position, 起始坐标方位角
//end_position   终止坐标方位角
void traverse_survey::start_end_position(double A_x, double A_y, double B_x, double B_y, double& start_position)
{
	double start_x = B_x - A_x;
	double start_y = B_y - A_y;
	double temp = atan2(start_y, start_x) * (180 / pi);
	if (start_y > 0)
	{
		start_position = fabs(temp);
	}
	if (start_y < 0)
	{
		start_position = 360 - fabs(temp);
	}
	transfrom_d_m_s(start_position);
}
//计算三角函数值
void traverse_survey::cal_trigo_function(double Decimal_degrees)
{
	cout << "sin = " << sin(Decimal_degrees) << endl;
	cout << "cos = " << cos(Decimal_degrees) << endl;
	cout << "tan = " << tan(Decimal_degrees) << endl;
}

int main() 
{
	int input = 0, N = 0;
	double* d, * m, * s, * Decimal_degrees, * every_position, * correct_every_position, * distance, * d_x, * d_y, * correct_d_x, * correct_d_y,
		x = 0.0, y = 0.0, sum_x = 0.0, sum_y = 0.0,
		A_x = 0.0, B_y = 0.0, A_y = 0.0, B_x = 0.0, start_position = 0.0, end_position = 0.0,
		C_x = 0.0, C_y = 0.0, sum = 0.0, f_x = 0.0, f_y = 0.0;
	long double* lx, * ly;

	traverse_survey ts;

	cout << "输入需要计算的测站数(数据数)" << endl;
	cin >> N;

	//initialize
	distance = new double[N] { 0 };
	Decimal_degrees = new double[N] { 0 };
	every_position = new double[N] { 0 };
	correct_every_position = new double[N] { 0 };
	correct_d_x = new double[N] { 0 };
	correct_d_y = new double[N] { 0 };
	d_x = new double[N] { 0 };
	d_y = new double[N] { 0 };
	lx = new long double[N] { 0 };
	ly = new long double[N] { 0 };
	s = new double[N] { 0 };
	m = new double[N] { 0 };
	d = new double[N] { 0 };

	cout << fixed << setprecision(4);

	while(1)
	{

		cout << "输入选项对应的数字：" << endl;
		cout << "0. 修改测站数(数据数)\n1. 度的六十进制转十进制\n2. \n3. \n"
			"4. 输入起始和末位的坐标\n5. 求导线的总长\n6. 求每一站的坐标方位角\n"
			"7. 改正后的每一站的坐标方位角\n8. 求相对坐标d_x和d_y的值\n9. 求总误差f_x和f_y的值\n"
			"10. 求全长闭合差和全长相对闭合差\n11. 求改正之后的d_x_y\n"
			"12. 求改正后每个测站的x,y值\n13. Debug Mode\n14. 所有结果清零" << endl;
		cin >> input;
		switch (input)
		{
		case 0:
			cout << "请输入新的数" << endl;
			cin >> N;
			distance = new double[N] { 0 };
			Decimal_degrees = new double[N] { 0 };
			every_position = new double[N];
			correct_every_position = new double[N] { 0 };
			correct_d_x = new double[N] { 0 };
			correct_d_y = new double[N] { 0 };
			d_x = new double[N] { 0 };
			d_y = new double[N] { 0 };
			lx = new long double[N] { 0 };
			ly = new long double[N] { 0 };
			s = new double[N] { 0 };
			m = new double[N] { 0 };
			d = new double[N] { 0 };
			break;
		case 1:
			cout << "按顺序输入，并用空格隔开" << endl;
			for (int i = 0; i < N; ++i)
			{
				cin >> d[i] >> m[i] >> s[i];
			}
			ts.h_d(d, m, s, N, Decimal_degrees);
			break;
		case 2:
			cout << "输入角度" << endl;
			for (int i = 0; i < N; ++i)
			{
				cin >> Decimal_degrees[i];
				ts.transfrom_d_m_s(Decimal_degrees[i]);
			}
			break;
		case 3:
			Decimal_degrees = new double[1];
			cin >> Decimal_degrees[0];
			ts.cal_trigo_function(Decimal_degrees[0]);
			break;
		case 4:
			cout << "分别输入前后的x，y值，用空格隔开" << endl;
			cin >> A_x >> A_y >> B_x >> B_y;
			//ts.start_end_position(A_x, A_y, B_x, B_y, start_position);
			break;
		case 5:
			cout << "请按照顺序输入距离" << endl;
			for (int i = 0; i < N; ++i)
			{
				cin >> distance[i];
			}
			ts.sum_distance(N, distance, sum);
			break;
		case 6:
			cout << "请输入起始方位角" << endl;
			cin >>start_position;

			//cout << "请输入起始坐标" << endl;
			//cin >> A_x >> A_y;

			/*cout << "请按顺序输入夹角" << endl;
			for (int i = 0; i < N; ++i)
			{
				cin >> Decimal_degrees[i];
			}*/

			ts.Angle_position(start_position, N, Decimal_degrees, every_position);
			break;
		case 7:
			cout << "输入结束方位角" << endl;
			cin >> end_position;
			//end_position -= start_position;

			ts.Correct_Angle_position(N, every_position, end_position, correct_every_position);
			break;
		case 8:
			ts.d_x_y(N, distance, correct_every_position, d_x, d_y);
			break;
		case 9:
			//cout << "输入结束点的坐标" << endl;
			//cin >> C_x >> C_y;

			ts.Cf_x(N, d_x, A_x, B_x, f_x);
			ts.Cf_y(N, d_y, A_y, B_y, f_y);
			break;
		case 10:
			ts.Cf_s(sum, f_x, f_y);
			break;
		case 11:
			ts.Correct_d_x_y(N, d_x, d_y, distance, sum, f_x, f_y, correct_d_x, correct_d_y);
			break;
		case 12:
			ts.x_y(N, A_x, A_y, correct_d_x, correct_d_y, lx, ly);
			break;
		case 13:
			cout << "按顺序输入各个测站的夹角" << endl;
			for (int i = 0; i < N; ++i)
			{
				cin >> d[i] >> m[i] >> s[i];
			}
			ts.h_d(d, m, s, N, Decimal_degrees);
			
			cout << "*************************************************" << endl;
			cout << "分别输入始末坐标的x，y值，用空格隔开" << endl;
			cin >> A_x >> A_y >> B_x >> B_y;

			cout << "*************************************************" << endl;
			cout << "请按照顺序输入距离" << endl;
			for (int i = 0; i < N; ++i)
			{
				cin >> distance[i];
			}
			ts.sum_distance(N, distance, sum);

			cout << "*************************************************" << endl;
			cout << "请输入起始方位角和终止方位角" << endl;
			cin >> start_position >> end_position;
			ts.Angle_position(start_position, N, Decimal_degrees, every_position);
			ts.Correct_Angle_position(N, every_position, end_position, correct_every_position);

			cout << "*************************************************" << endl;
			ts.d_x_y(N, distance, correct_every_position, d_x, d_y);
			ts.Cf_x(N, d_x, A_x, B_x, f_x);
			ts.Cf_y(N, d_y, A_y, B_y, f_y);
			ts.Cf_s(sum, f_x, f_y);
			ts.Correct_d_x_y(N, d_x, d_y, distance, sum, f_x, f_y, correct_d_x, correct_d_y);

			cout << "*************************************************" << endl;
			ts.x_y(N, A_x, A_y, correct_d_x, correct_d_y, lx, ly);
			break;
		case 14:
			x = 0.0, y = 0.0, sum_x = 0.0, sum_y = 0.0,
				A_x = 0.0, B_y = 0.0, A_y = 0.0, B_x = 0.0, start_position = 0.0, end_position = 0.0,
				C_x = 0.0, C_y = 0.0, sum = 0.0;

			delete[]distance;
			delete[]Decimal_degrees;
			delete[]every_position;
			delete[]correct_d_x;
			delete[]correct_d_y;
			delete[]d_x;
			delete[]d_y;
			delete[]lx;
			delete[]ly;
			delete[]s;
			delete[]m;
			delete[]d;
			break;
		default:
			break;
		}
	};
	
	return 0;
}
