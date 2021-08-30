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
	cout << "-----------------------����---------------------------" << endl;
}
//------------------  �ѹ۲�Ķȷ���ת����ʮ���ƵĶ�ʹ��------------------------//
//*d----------------����Ķ�
//*m----------------����ķ�
//*s----------------�������
//N-----------------վ��
//*Decimal_degrees----------�����ÿһվת���Ĺ۲�Ƕȵ�--ʮ���ƵĶȡ�
void traverse_survey::h_d(double* d, double* m, double* s, int N, double* Decimal_degrees)
{
	for (int i = 0; i < N; i++)
	{
		Decimal_degrees[i] = d[i] + m[i] / 60 + s[i] / 3600;
	}
}
//------------------ʮ���ƶ� ת��Ϊ�ȷ������------------------------//
//Decimal_degrees--------------ʮ���ƵĶ�
void traverse_survey::transfrom_d_m_s(double Decimal_degrees)
{
	double d = 0, m = 0, s = 0;
	d = int(Decimal_degrees);
	m = int((Decimal_degrees - d) * 60);
	s = (Decimal_degrees - d - m / 60) * 3600;
	cout << d << "��" << m << "��" << s << "��" << endl;
}
//------------------��ÿһվ�����귽λ��-�ͽǶȱպϲ�--------------------------//
//start_positon--------------��ʼ��λ��
//N--------------------------վ��
//*Decimal_degrees-----------ÿһվ�Ĺ۲�Ƕ�
//*every_position------------����ÿһվ�ķ�λ�� (ʮ���ƵĶ�)
void traverse_survey::Angle_position(double start_positon, int N, double* Decimal_degrees, double* every_position)
{
	for (int i = 0; i < N; i++)
	{
		if (i == 0)
		{
			every_position[0] = start_positon + Decimal_degrees[0];
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
			every_position[i] = every_position[i - 1] + Decimal_degrees[i] - 360;
			if (every_position[i] < 0)
			{
				every_position[i] += 360;
			}
			else if (every_position[i] > 360)
			{
				every_position[i] -= 360;
			}
		}
		cout << "��" << i + 1 << "�����귽λ��Ϊ��";
		transfrom_d_m_s(every_position[i]);
		//cout << endl;
	}
}
//------------------------------�������ÿһվ�����귽λ��------------------------------//
//N------------------------------վ��
//*every_position  --------------ÿһվ�����귽λ��
//end_position-------------------���һվ�����귽λ��
// *correct_every_position-------�������ÿһվ���귽λ��
void traverse_survey::Correct_Angle_position(int N, double* every_position, double end_position, double* correct_every_position)
{
	double temp = 0;
	cout << "�õ��ߵĽǶȱպϲ�Ϊ��";
	transfrom_d_m_s(every_position[N - 1] - end_position);
	cout << endl;
	temp = (every_position[N - 1] - end_position) / N;

	for (int i = 0; i < N; i++)
	{
		correct_every_position[i] = every_position[i] - (i + 1) * temp;
		cout << "��" << i + 1 << "������������귽λ��Ϊ��";
		transfrom_d_m_s(correct_every_position[i]);
		cout << endl;
	}
}
//---------------------��d_x_y��ֵ----------------------------//
//N------------------------------------------վ��
//*distance----------------------------------ÿһվ�ľ���
//*correct_every_position--------------------ÿһվ����������귽λ��
//*d_x---------------------------------------�洢ÿһվ��d_x
//*d_y---------------------------------------�洢ÿһվ��d_y
void traverse_survey::d_x_y(double N, double* distance, double* correct_every_position, double* d_x, double* d_y)
{
	for (int i = 0; i < N; i++)
	{
		d_x[i] = distance[i] * cos(correct_every_position[i] * (pi / 180));
		d_y[i] = distance[i] * sin(correct_every_position[i] * (pi / 180));
		cout << i + 1 << "�Ų�վ��d_x��d_y�ֱ�Ϊ��" << d_x[i] << "�� " << d_y[i] << endl;
	}
}
//-------------------------���ߵ��ܳ�--------------------------//
//N-----------------վ��
//*distance---------ÿһվ�ľ���
//sum --------------���ߵ��ܳ�
void traverse_survey::sum_distance(int N, double* distance, double& sum)
{
	for (int i = 0; i < N - 1; i++)
	{
		sum += distance[i];
	}
	cout << "���ߵ��ܳ���" << sum << endl;
}

//---------------------------��f_x-----------------------------//
//N---------------------վ��
//*d_x------------------ÿһվ��d_x
//B_x-------------------��ʼ���x����
//C_x-------------------��ֹ���x����
//f_x-----------------f_x��ֵ
void traverse_survey::f_x(double N, double* d_x, double B_x, double C_x, double& f_x)
{
	double sum = 0;
	double temp = B_x - C_x;
	for (int i = 0; i < N - 1; i++)
	{
		sum += d_x[i];
	}
	f_x = temp + sum;
	cout << "f_x��ֵΪ��" << f_x << endl;
}
//---------------------------��f_y-----------------------------//
//N---------------------վ��
//*d_y------------------ÿһվ��d_y
//B_y-------------------��ʼ���y����
//C_x-------------------��ֹ���y����
//f_y-----------------f_y��ֵ
void traverse_survey::f_y(double N, double* d_y, double B_y, double C_y, double& f_y)
{
	double sum = 0;
	double temp = B_y - C_y;
	for (int i = 0; i < N - 1; i++)
	{
		sum += d_y[i];
	}
	f_y = temp + sum;
	cout << "f_x��ֵΪ��" << f_y << endl;
}
//---------------------��ȫ���պϲ��ȫ����Ապϲ�---------------------//
//f_x-------------------------����f_x��ֵ
//f_y-------------------------����f_y��ֵ
void traverse_survey::f_s(double sum, double f_x, double f_y)
{
	double f_s = sqrt(pow(f_x, 2) + pow(f_y, 2));
	cout << "ȫ���պϲ�f_s��ֵΪ:  " << f_s << endl;
	double K = 1 / (sum / f_s);
	cout << "ȫ����Ապϲ�KΪ��" << K << endl;
	cout << "1/K = " << int(1 / K);
	/*if (int(1 / K) > 4000)
	{
		cout << "�õ���û�г��ޣ�" << endl;
	}
	else
	{
		cout << "�õ��߳��ޣ�" << endl;
	}*/

}
//-------------------------�����֮���d_x_y----------------------//
//N----------------------------վ��
//*d_x-------------------------ÿһվ��d_x��ֵ
// *d_y------------------------ÿһվ��d_y��ֵ
//*distance--------------------ÿһվ�ľ���
//sum--------------------------���ߵ��ܳ�
//f_x--------------------------f_x��ֵ
//f_y--------------------------f_y��ֵ
void traverse_survey::Correct_d_x_y(int N, double* d_x, double* d_y, double* distance, double sum, double f_x, double f_y, double* correct_d_x, double* correct_d_y)
{
	double* temp_x = new double[N - 1]{ 0 };        //��̬�Ĵ���һ������
	double* temp_y = new double[N - 1]{ 0 };
	for (int i = 0; i < N - 1; i++)
	{
		temp_x[i] = f_x * (distance[i] / sum);
		temp_y[i] = f_y * (distance[i] / sum);
		cout << "-----------------------------------------------------------------" << endl;
		cout << "��" << i + 1 << "����" << "��������----��d_x��d_y��ֵΪ��" << temp_x[i] << "  ||  " << temp_y[i] << endl;
		correct_d_x[i] = d_x[i] - temp_x[i];
		correct_d_y[i] = d_y[i] - temp_y[i];
		cout << "��" << i + 1 << "����" << "����֮��ġ���d_x��d_y��ֵΪ��" << correct_d_x[i] << "  ||  " << correct_d_y[i] << endl;
	}
}
//-------------------------------------��x,y--------------------------------//;
//N---------------------------վ��
//A_x-------------------------��ʼ��x����
//A_y-------------------------��ʼ��y����
//*correct_d_x----------------�������d_x
//*correct_d_y----------------�������d_y
//*x--------------------------����ÿһվx������
//*y--------------------------����ÿһվy������
void traverse_survey::x_y(int N, double A_x, double A_y, double* correct_d_x, double* correct_d_y, long double* x, long double* y)//��x,y;
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
		printf("��%dվ����%.4lf  ||   %.4lf\n", i + 1, x[i], y[i]);

	}
}
//����ʼ��ĩλ�����귽λ�ǡ�
//A_x,A_y,B_x, B_y  ��ʼ��
//C_x,C_y,D_x, D_y  ��ֹ��
//start_position, ��ʼ���귽λ��
//end_position)   ��ֹ���귽λ��
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
//�������Ǻ���ֵ
void traverse_survey::cal_trigo_function(double Decimal_degrees)
{
	cout << "sin = " << sin(Decimal_degrees) << endl;
	cout << "cos = " << cos(Decimal_degrees) << endl;
	cout << "tan = " << tan(Decimal_degrees) << endl;
}

int main() 
{
	int input, N;
	double *d, *m, *s, *Decimal_degrees, x, y, sum_x, sum_y,
		A_x, B_y, A_y, B_x, start_position, *every_position, end_position,
		correct_every_position, *distance, *d_x, *d_y,
		C_x, C_y, f_x, f_y, sum, correct_d_x, correct_d_y;
	long double lx, ly;

	traverse_survey ts;

	while(1)
	{
		cout << "����ѡ���Ӧ�����֣�" << endl;
		cout << "1. �Ƕ�ת����\n2. �ȵ�ʮ����ת��ʮ����\n3. �����Ǻ���ֵ\n"
			"4. ����ʼ��ĩλ�����귽λ��\n5. ���ߵ��ܳ�\n6. ��ÿһվ�����귽λ��������(��������)\n"
			"7. �������ÿһվ�����귽λ��\n8. ��d_x_y��ֵ\n9. ��f_x\n"
			"10. ��f_y\n11. ��ȫ���պϲ��ȫ����Ապϲ�\n12. �����֮���d_x_y\n"
			"13. ��x,y\n14. ���н������" << endl;

		cin >> input;
		switch (input)
		{
		case 1:
			cout << "����ת������" << endl;
			cin >> N;
			s = new double[N];
			m = new double[N];
			d = new double[N];
			Decimal_degrees = new double[N];
			cout << "��˳�����룬���ÿո����" << endl;
			for (int i = 0; i < N; ++i)
			{
				cin >> d[i] >> m[i] >> s[i];
			}
			ts.h_d(d, m, s, N, Decimal_degrees);
			delete[]s;
			delete[]m;
			delete[]d;
			delete[]Decimal_degrees;
			break;
		case 2:
			cout << "����Ƕ�" << endl;
			Decimal_degrees = new double[1];
			cin >> Decimal_degrees[0];
			ts.transfrom_d_m_s(Decimal_degrees[0]);
			break;
		case 3:
			Decimal_degrees = new double[1];
			cin >> Decimal_degrees[0];
			ts.cal_trigo_function(Decimal_degrees[0]);
			break;
		case 4:
			cout << "�ֱ�����ǰ���x��yֵ���ÿո����" << endl;
			cin >> A_x >> A_y >> B_x >> B_y;
			ts.start_end_position(A_x, A_y, B_x, B_y, start_position);
			break;
		case 5:
			cout << "�����ڴ�" << endl;
			break;
		case 6:
			cout << "��������Ҫ����Ĳ�վ������ʼ��λ��" << endl;
			cin >> N >>start_position;
			Decimal_degrees = new double[N];
			every_position = new double[N];
			distance = new double[N];
			d_x = new double[N];
			d_y = new double[N];

			cout << "��������ʼ����" << endl;
			cin >> A_x >> A_y;

			cout << "�밴˳������н�" << endl;
			for (int i = 0; i < N; ++i)
			{
				cin >> Decimal_degrees[i];
			}

			ts.Angle_position(start_position, N, Decimal_degrees, every_position);
			cout << endl;

			cout << "�밴˳���������" << endl;
			for (int i = 0; i < N; ++i)
			{
				cin >> distance[i];
			}

			ts.d_x_y(N, distance, every_position, d_x, d_y);
			cout << endl;
			sum_x = A_x, sum_y = A_y;
			for (int i = 0; i < N; ++i)
			{
				sum_x += d_x[i];
				sum_y += d_y[i];
				cout << i + 1 << "�Ų�վ������Ϊ " << fixed << setprecision(4) << sum_x <<  ", " << sum_y << endl;
			}

			delete[]Decimal_degrees;
			delete[]every_position;
			delete[]distance;
			delete[]d_x;
			delete[]d_y;
			break;
		case 7:
			cout << "�����ڴ�" << endl;
			break;
		case 8:
			cout << "�����ڴ�" << endl;
			break;
		case 9:
			cout << "�����ڴ�" << endl;
			break;
		case 10:
			cout << "�����ڴ�" << endl;
			break;
		case 11:
			break;
		case 12:
			cout << "�����ڴ�" << endl;
			break;
		case 13:
			cout << "�����ڴ�" << endl;
			break;
		case 14:
			break;
		default:
			break;
		}
	};
	
	return 0;
}